/* ============================================================================
 * scroll_overlap -- driver for the overlap-repair (quilting) module.
 * Loads a UV-carrying OBJ (scroll_ribbon's <id>_uv.obj), runs Quilt_run, and
 * writes the repaired unroll + a decision-colored inspection OBJ + stats.
 * Re-bake the texture with obj_bake_raw on the _repaired_uv.obj.
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../common/arena.h"
#include "../common/pca.h"
#include "../common/ves_platform.h"
#include "../flatten/overlap_quilt.h"

/* ---- minimal OBJ (+vt) reader: verts (z,y,x), vt (u,v) 1:1, faces a/a ---- */
typedef struct { float *v; float *t; int32_t *f; size_t nv, nvt, nf; } Obj;

static void *grow(void *p, size_t need, size_t *cap, size_t el) {
    size_t c = *cap; if (need <= c) return p;
    c = c ? c : 4096; while (c < need) c *= 2;
    p = realloc(p, c*el); if (!p){fprintf(stderr,"OOM\n");exit(1);} *cap=c; return p;
}

static int obj_load(const char *path, Obj *m) {
    FILE *fp = fopen(path, "rb"); char line[1024];
    size_t cv=0, ctv=0, cf=0;
    memset(m, 0, sizeof *m);
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); return -1; }
    while (fgets(line, sizeof line, fp)) {
        if (line[0]=='v' && line[1]==' ') {
            char *s=line+2,*e; double a=strtod(s,&e); s=e; double b=strtod(s,&e); s=e; double c=strtod(s,&e);
            m->v=(float*)grow(m->v,(m->nv+1)*3,&cv,sizeof(float));
            m->v[m->nv*3+0]=(float)a; m->v[m->nv*3+1]=(float)b; m->v[m->nv*3+2]=(float)c; m->nv++;
        } else if (line[0]=='v'&&line[1]=='t'&&line[2]==' ') {
            char *s=line+3,*e; double u=strtod(s,&e); s=e; double v=strtod(s,&e);
            m->t=(float*)grow(m->t,(m->nvt+1)*2,&ctv,sizeof(float));
            m->t[m->nvt*2+0]=(float)u; m->t[m->nvt*2+1]=(float)v; m->nvt++;
        } else if (line[0]=='f'&&(line[1]==' '||line[1]=='\t')) {
            long idx[16]; int ni=0; char *s=line+2;
            while (ni<16){ char*e; while(*s==' '||*s=='\t')s++; if(!*s||*s=='\n'||*s=='\r')break;
                long v=strtol(s,&e,10); if(e==s)break; idx[ni++]=v; s=e;
                while(*s&&*s!=' '&&*s!='\t'&&*s!='\n'&&*s!='\r')s++; }
            for (int k=2;k<ni;k++){ m->f=(int32_t*)grow(m->f,(m->nf+1)*3,&cf,sizeof(int32_t));
                m->f[m->nf*3+0]=(int32_t)(idx[0]-1); m->f[m->nf*3+1]=(int32_t)(idx[k-1]-1);
                m->f[m->nf*3+2]=(int32_t)(idx[k]-1); m->nf++; }
        }
    }
    fclose(fp);
    if (m->nv==0 || m->nf==0 || (m->nvt!=0 && m->nvt!=m->nv)) {
        fprintf(stderr, "%s: nv=%zu nvt=%zu nf=%zu (need faces and either 0 or nv vt)\n",
                path, m->nv, m->nvt, m->nf);
        return -1;
    }
    return 0;
}

/* Retain faces whose 3D triangle bbox intersects the requested crop, then
 * compact vertices. bbox order is zmin,zmax,ymin,ymax,xmin,xmax. */
static int obj_crop(Obj *m, const double bbox[6]) {
    uint8_t *keep = (uint8_t *)calloc(m->nf, 1);
    uint8_t *used = (uint8_t *)calloc(m->nv, 1);
    int32_t *remap = NULL;
    float *new_v = NULL;
    float *new_t = NULL;
    int32_t *new_f = NULL;
    size_t kept_faces = 0;
    size_t kept_verts = 0;

    if (keep == NULL || used == NULL) {
        fprintf(stderr, "scroll_overlap: OOM cropping mesh\n");
        exit(1);
    }
    for (size_t fi = 0; fi < m->nf; fi++) {
        double lo[3] = {1e30, 1e30, 1e30};
        double hi[3] = {-1e30, -1e30, -1e30};
        int intersects = 1;

        for (int k = 0; k < 3; k++) {
            int32_t vi = m->f[fi*3+(size_t)k];
            for (int d = 0; d < 3; d++) {
                double x = m->v[(size_t)vi*3+(size_t)d];
                if (x < lo[d]) lo[d] = x;
                if (x > hi[d]) hi[d] = x;
            }
        }
        for (int d = 0; d < 3; d++) {
            if (hi[d] < bbox[d*2] || lo[d] > bbox[d*2+1]) {
                intersects = 0;
            }
        }
        if (intersects) {
            keep[fi] = 1;
            kept_faces++;
            for (int k = 0; k < 3; k++) {
                used[(size_t)m->f[fi*3+(size_t)k]] = 1;
            }
        }
    }
    if (kept_faces == 0) {
        fprintf(stderr, "scroll_overlap: crop contains no faces\n");
        free(keep);
        free(used);
        return -1;
    }

    remap = (int32_t *)malloc(m->nv * sizeof *remap);
    if (remap == NULL) {
        fprintf(stderr, "scroll_overlap: OOM building crop remap\n");
        exit(1);
    }
    for (size_t vi = 0; vi < m->nv; vi++) {
        remap[vi] = -1;
        if (used[vi]) {
            remap[vi] = (int32_t)kept_verts++;
        }
    }
    new_v = (float *)malloc(kept_verts * 3 * sizeof *new_v);
    if (m->nvt == m->nv) {
        new_t = (float *)malloc(kept_verts * 2 * sizeof *new_t);
    }
    new_f = (int32_t *)malloc(kept_faces * 3 * sizeof *new_f);
    if (new_v == NULL || new_f == NULL ||
        (m->nvt == m->nv && new_t == NULL)) {
        fprintf(stderr, "scroll_overlap: OOM compacting crop\n");
        exit(1);
    }
    for (size_t vi = 0; vi < m->nv; vi++) {
        int32_t dst = remap[vi];
        if (dst < 0) continue;
        memcpy(&new_v[(size_t)dst*3], &m->v[vi*3], 3*sizeof(float));
        if (new_t != NULL) {
            memcpy(&new_t[(size_t)dst*2], &m->t[vi*2], 2*sizeof(float));
        }
    }
    size_t out_face = 0;
    for (size_t fi = 0; fi < m->nf; fi++) {
        if (!keep[fi]) continue;
        for (int k = 0; k < 3; k++) {
            new_f[out_face*3+(size_t)k] = remap[(size_t)m->f[fi*3+(size_t)k]];
        }
        out_face++;
    }

    fprintf(stderr, "  crop: %zu v, %zu f from %zu v, %zu f\n",
            kept_verts, kept_faces, m->nv, m->nf);
    free(m->v);
    free(m->t);
    free(m->f);
    free(keep);
    free(used);
    free(remap);
    m->v = new_v;
    m->t = new_t;
    m->f = new_f;
    m->nv = kept_verts;
    m->nvt = new_t != NULL ? kept_verts : 0;
    m->nf = kept_faces;
    return 0;
}

/* Derive one common tangent-plane projection for a local pair of sheets. */
static int obj_project_pca(Obj *m) {
    float normal[3];
    float centroid[3];
    float u_axis[3];
    float v_axis[3];

    if (PCA_normal(m->v, m->nv, normal, centroid) != 0) {
        fprintf(stderr, "scroll_overlap: PCA projection failed\n");
        return -1;
    }
    PCA_orthonormal_basis(normal, u_axis, v_axis);
    free(m->t);
    m->t = (float *)malloc(m->nv * 2 * sizeof *m->t);
    if (m->t == NULL) {
        fprintf(stderr, "scroll_overlap: OOM projecting mesh\n");
        exit(1);
    }
    for (size_t vi = 0; vi < m->nv; vi++) {
        double d[3];
        d[0] = (double)m->v[vi*3+0] - centroid[0];
        d[1] = (double)m->v[vi*3+1] - centroid[1];
        d[2] = (double)m->v[vi*3+2] - centroid[2];
        m->t[vi*2+0] = (float)(d[0]*u_axis[0] + d[1]*u_axis[1] + d[2]*u_axis[2]);
        m->t[vi*2+1] = (float)(d[0]*v_axis[0] + d[1]*v_axis[1] + d[2]*v_axis[2]);
    }
    m->nvt = m->nv;
    fprintf(stderr,
            "  PCA projection: normal=(%.6f,%.6f,%.6f) "
            "u=(%.6f,%.6f,%.6f) v=(%.6f,%.6f,%.6f)\n",
            (double)normal[0], (double)normal[1], (double)normal[2],
            (double)u_axis[0], (double)u_axis[1], (double)u_axis[2],
            (double)v_axis[0], (double)v_axis[1], (double)v_axis[2]);
    return 0;
}

/* write OBJ with vt, keeping only faces where keep[f] (NULL = all) */
static int obj_write_uv(const char *path, const float *v, size_t nv,
                        const int32_t *f, size_t nf, const float *t,
                        const uint8_t *keep) {
    FILE *fp = fopen(path, "wb"); size_t i;
    if (!fp) return -1;
    for (i=0;i<nv;i++) fprintf(fp,"v %.6f %.6f %.6f\n",(double)v[i*3+0],(double)v[i*3+1],(double)v[i*3+2]);
    for (i=0;i<nv;i++) fprintf(fp,"vt %.6f %.6f\n",(double)t[i*2+0],(double)t[i*2+1]);
    for (i=0;i<nf;i++){ if(keep&&!keep[i])continue; int a=f[i*3+0]+1,b=f[i*3+1]+1,c=f[i*3+2]+1;
        fprintf(fp,"f %d/%d %d/%d %d/%d\n",a,a,b,b,c,c); }
    fclose(fp); return 0;
}

static void obj_write_colored_face(FILE *fp, const float *v,
                                   const int32_t *f, size_t fi,
                                   size_t out_face, const double rgb[3]) {
    size_t base = out_face * 3 + 1;

    for (int k = 0; k < 3; k++) {
        size_t vi = (size_t)f[fi*3+(size_t)k];
        fprintf(fp, "v %.6f %.6f %.6f %.4f %.4f %.4f\n",
                (double)v[vi*3+0], (double)v[vi*3+1], (double)v[vi*3+2],
                rgb[0], rgb[1], rgb[2]);
    }
    fprintf(fp, "f %zu %zu %zu\n", base, base+1, base+2);
}

/* Duplicate vertices per face so layer/decision boundaries retain hard colors. */
static int obj_write_dec(const char *path, const float *v,
                         const int32_t *f, size_t nf, const uint8_t *dec) {
    static const double colors[4][3] = {
        {0.45, 0.45, 0.45},
        {0.10, 0.90, 0.20},
        {0.15, 0.40, 0.95},
        {0.95, 0.12, 0.10}
    };
    FILE *fp = fopen(path, "wb");

    if (fp == NULL) return -1;
    fprintf(fp, "# quilting decision: gray=outside green=selected blue=relocated red=not-selected\n");
    for (size_t fi = 0; fi < nf; fi++) {
        uint8_t d = dec[fi] <= 3 ? dec[fi] : 0;
        obj_write_colored_face(fp, v, f, fi, fi, colors[d]);
    }
    fclose(fp);
    return 0;
}

/* Energy colors: green is the lowest finite layer energy, red the highest;
 * magenta means the layer had no usable shared anchor seam; gray is unscored. */
static int obj_write_energy(const char *path, const float *v,
                            const int32_t *f, size_t nf,
                            const double *energy) {
    FILE *fp = NULL;
    double lo = 1e30;
    double hi = -1e30;

    for (size_t fi = 0; fi < nf; fi++) {
        double e = energy[fi];
        if (isfinite(e) && e >= 0.0 && e < 1e17) {
            if (e < lo) lo = e;
            if (e > hi) hi = e;
        }
    }
    fp = fopen(path, "wb");
    if (fp == NULL) return -1;
    if (hi >= lo) {
        fprintf(fp, "# quilting boundary energy: green=%.9g red=%.9g magenta=no-seam gray=unscored\n",
                lo, hi);
    } else {
        fprintf(fp, "# quilting boundary energy: no finite layer energies; magenta=no-seam gray=unscored\n");
    }
    for (size_t fi = 0; fi < nf; fi++) {
        double e = energy[fi];
        double rgb[3] = {0.45, 0.45, 0.45};

        if (e >= 1e17) {
            rgb[0] = 0.90;
            rgb[1] = 0.10;
            rgb[2] = 0.90;
        } else if (isfinite(e) && e >= 0.0) {
            double t = hi > lo + 1e-12 ? (e-lo)/(hi-lo) : 0.0;
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;
            rgb[0] = fmin(1.0, 2.0*t);
            rgb[1] = fmin(1.0, 2.0*(1.0-t));
            rgb[2] = 0.05;
        }
        obj_write_colored_face(fp, v, f, fi, fi, rgb);
    }
    fclose(fp);
    return 0;
}

typedef struct LayerStat {
    int32_t region;
    int32_t layer;
    size_t faces;
    size_t selected;
    size_t relocated;
    size_t rejected;
    double energy;
} LayerStat;

static int obj_write_layer_csv(const char *path, size_t nf,
                               const QuiltResult *result) {
    LayerStat *stats = NULL;
    size_t nstats = 0;
    size_t cap = 0;
    FILE *fp = NULL;

    for (size_t fi = 0; fi < nf; fi++) {
        int32_t region = result->face_region[fi];
        int32_t layer = result->face_layer[fi];
        size_t si = 0;

        if (region < 0 || layer < 0) continue;
        while (si < nstats &&
               (stats[si].region != region || stats[si].layer != layer)) {
            si++;
        }
        if (si == nstats) {
            stats = (LayerStat *)grow(stats, nstats+1, &cap, sizeof *stats);
            memset(&stats[nstats], 0, sizeof stats[nstats]);
            stats[nstats].region = region;
            stats[nstats].layer = layer;
            stats[nstats].energy = result->face_energy[fi];
            nstats++;
        }
        stats[si].faces++;
        if (result->face_dec[fi] == 1) stats[si].selected++;
        else if (result->face_dec[fi] == 2) stats[si].relocated++;
        else if (result->face_dec[fi] == 3) stats[si].rejected++;
    }

    fp = fopen(path, "wb");
    if (fp == NULL) {
        free(stats);
        return -1;
    }
    fprintf(fp, "region,layer,faces,boundary_energy,status,selected,relocated,not_selected\n");
    for (size_t si = 0; si < nstats; si++) {
        const char *status = "finite";
        if (stats[si].energy < 0.0 || !isfinite(stats[si].energy)) status = "unscored";
        else if (stats[si].energy >= 1e17) status = "no_seam";
        fprintf(fp, "%d,%d,%zu,%.9g,%s,%zu,%zu,%zu\n",
                stats[si].region, stats[si].layer, stats[si].faces,
                stats[si].energy, status, stats[si].selected,
                stats[si].relocated, stats[si].rejected);
    }
    fclose(fp);
    free(stats);
    return 0;
}

static double jget(const char *path, const char *key) {
    /* crude JSON scalar fetch: find "key" then the next number */
    FILE *fp=fopen(path,"rb"); if(!fp) return 0.0; char *buf; long n;
    fseek(fp,0,SEEK_END); n=ftell(fp); fseek(fp,0,SEEK_SET);
    buf=(char*)malloc((size_t)n+1); if(!buf){fclose(fp);return 0.0;}
    if(fread(buf,1,(size_t)n,fp)!=(size_t)n){} buf[n]=0; fclose(fp);
    char *p=strstr(buf,key); double val=0.0;
    if(p){ p+=strlen(key); while(*p&&(*p==' '||*p==':'||*p=='"'))p++; val=strtod(p,NULL); }
    free(buf); return val;
}

int main(int argc, char **argv) {
    const char *in = NULL;
    const char *raw = NULL;
    const char *outdir = NULL;
    const char *id = "overlap";
    const char *statsj = NULL;
    QuiltOpts o;
    Obj m;
    QuiltResult result;
    Arena_T arena = NULL;
    double bbox[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double t0 = 0.0;
    int use_raw = 1;
    int project_pca = 0;
    int have_crop = 0;
    int diagnostic_only = 0;
    char path[4096];

    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) {
        return Quilt_selftest() == 0 ? 0 : 1;
    }
    if (argc < 4) {
        fprintf(stderr,
                "usage: %s <mesh.obj> <raw_cube_dir> <out_dir> [--id S]\n"
                "   [--project-pca] [--crop zmin zmax ymin ymax xmin xmax]\n"
                "   [--diagnostic-only] [--normal-range F] [--normal-samples N]\n"
                "   [--tex-du F] [--tex-dv F] [--radius-gate F] [--cell F]\n"
                "   [--cross-component-only]\n"
                "   [--axis-point z y x] [--axis-dir z y x] [--stats-json F]\n"
                "   [--no-relocate] [--no-raw]\n"
                "       %s --selftest\n",
                argv[0], argv[0]);
        return 1;
    }

    in = argv[1];
    raw = argv[2];
    outdir = argv[3];
    QuiltOpts_default(&o);
    o.raw_dir = raw;
    o.verbose = 1;
    for (int i = 4; i < argc; i++) {
        if (strcmp(argv[i], "--id") == 0 && i+1 < argc) {
            id = argv[++i];
        } else if (strcmp(argv[i], "--stats-json") == 0 && i+1 < argc) {
            statsj = argv[++i];
        } else if (strcmp(argv[i], "--axis-point") == 0 && i+3 < argc) {
            o.axis_point[0] = (float)atof(argv[++i]);
            o.axis_point[1] = (float)atof(argv[++i]);
            o.axis_point[2] = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--axis-dir") == 0 && i+3 < argc) {
            o.axis_dir[0] = (float)atof(argv[++i]);
            o.axis_dir[1] = (float)atof(argv[++i]);
            o.axis_dir[2] = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--radius-gate") == 0 && i+1 < argc) {
            o.radius_gate = atof(argv[++i]);
        } else if (strcmp(argv[i], "--cell") == 0 && i+1 < argc) {
            o.cell = atof(argv[++i]);
        } else if (strcmp(argv[i], "--normal-range") == 0 && i+1 < argc) {
            o.normal_range = atof(argv[++i]);
        } else if (strcmp(argv[i], "--normal-samples") == 0 && i+1 < argc) {
            o.normal_samples = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--tex-du") == 0 && i+1 < argc) {
            o.tex_du = atof(argv[++i]);
        } else if (strcmp(argv[i], "--tex-dv") == 0 && i+1 < argc) {
            o.tex_dv = atof(argv[++i]);
        } else if (strcmp(argv[i], "--crop") == 0 && i+6 < argc) {
            for (int k = 0; k < 6; k++) bbox[k] = atof(argv[++i]);
            have_crop = 1;
        } else if (strcmp(argv[i], "--project-pca") == 0) {
            project_pca = 1;
        } else if (strcmp(argv[i], "--diagnostic-only") == 0) {
            diagnostic_only = 1;
            o.relocate = 0;
        } else if (strcmp(argv[i], "--cross-component-only") == 0) {
            o.cross_component_only = 1;
        } else if (strcmp(argv[i], "--no-relocate") == 0) {
            o.relocate = 0;
        } else if (strcmp(argv[i], "--no-raw") == 0) {
            use_raw = 0;
        } else {
            fprintf(stderr, "unknown or incomplete option %s\n", argv[i]);
            return 1;
        }
    }
    if (!use_raw) o.raw_dir = NULL;
    if (statsj != NULL) {
        double pitch = 0.0;
        o.spiral_a = jget(statsj, "spiral_a");
        o.spiral_b = jget(statsj, "spiral_b");
        pitch = jget(statsj, "pitch_used");
        if (pitch > 0.0) o.pitch = pitch;
    }

    fprintf(stderr, "[scroll_overlap] in=%s raw=%s out=%s id=%s diagnostic=%d\n",
            in, raw, outdir, id, diagnostic_only);
    t0 = ves_clock_sec();
    if (obj_load(in, &m) != 0) return 1;
    fprintf(stderr, "  loaded %zu v, %zu f (%.1fs)\n",
            m.nv, m.nf, ves_clock_sec()-t0);
    if (have_crop && obj_crop(&m, bbox) != 0) {
        free(m.v);
        free(m.t);
        free(m.f);
        return 1;
    }
    if (project_pca && obj_project_pca(&m) != 0) {
        free(m.v);
        free(m.t);
        free(m.f);
        return 1;
    }
    if (m.nvt != m.nv) {
        fprintf(stderr, "scroll_overlap: mesh has no per-vertex vt; use --project-pca\n");
        free(m.v);
        free(m.t);
        free(m.f);
        return 1;
    }

    snprintf(path, sizeof path, "%s/%s_projected_uv.obj", outdir, id);
    ves_ensure_parent_dir(path);
    if (obj_write_uv(path, m.v, m.nv, m.f, m.nf, m.t, NULL) == 0) {
        fprintf(stderr, "  wrote %s\n", path);
    }

    arena = Arena_new();
    t0 = ves_clock_sec();
    if (Quilt_run(arena, m.v, m.nv, m.f, m.nf, m.t, &o, &result) != 0) {
        fprintf(stderr, "Quilt_run failed\n");
        free(m.v);
        free(m.t);
        free(m.f);
        Arena_dispose(&arena);
        return 1;
    }
    fprintf(stderr,
            "  quilt: regions=%zu overlap_faces=%zu layers=%zu max_region=%zu "
            "kept=%zu relocated=%zu dropped=%zu multi_cells %zu->%zu (%.1fs)\n",
            result.n_regions, result.n_overlap_faces, result.n_layers_total,
            result.max_region_faces, result.n_kept, result.n_relocated,
            result.n_dropped, result.multi_cells_before,
            result.multi_cells_after, ves_clock_sec()-t0);

    if (!diagnostic_only) {
        snprintf(path, sizeof path, "%s/%s_repaired_uv.obj", outdir, id);
        if (obj_write_uv(path, m.v, m.nv, m.f, m.nf,
                         result.uv, result.face_keep) == 0) {
            fprintf(stderr, "  wrote %s\n", path);
        }
    }
    snprintf(path, sizeof path, "%s/%s_layers.obj", outdir, id);
    if (obj_write_dec(path, m.v, m.f, m.nf, result.face_dec) == 0) {
        fprintf(stderr, "  wrote %s\n", path);
    }
    snprintf(path, sizeof path, "%s/%s_energy.obj", outdir, id);
    if (obj_write_energy(path, m.v, m.f, m.nf, result.face_energy) == 0) {
        fprintf(stderr, "  wrote %s\n", path);
    }
    snprintf(path, sizeof path, "%s/%s_layer_energy.csv", outdir, id);
    if (obj_write_layer_csv(path, m.nf, &result) == 0) {
        fprintf(stderr, "  wrote %s\n", path);
    }
    snprintf(path, sizeof path, "%s/%s_overlap_stats.json", outdir, id);
    FILE *fp = fopen(path, "wb");
    if (fp != NULL) {
        fprintf(fp,
                "{ \"id\":\"%s\", \"diagnostic_only\":%d, \"project_pca\":%d,\n"
                "  \"regions\":%zu, \"overlap_faces\":%zu, \"layers_total\":%zu,\n"
                "  \"max_region_faces\":%zu, \"kept\":%zu, \"relocated\":%zu, \"dropped\":%zu,\n"
                "  \"multi_cells_before\":%zu, \"multi_cells_after\":%zu,\n"
                "  \"have_texture\":%d, \"energy_mean\":%.9g,\n"
                "  \"normal_range\":%.9g, \"normal_samples\":%d, \"radius_gate\":%.9g,\n"
                "  \"cross_component_only\":%d,\n",
                id, diagnostic_only, project_pca, result.n_regions,
                result.n_overlap_faces, result.n_layers_total,
                result.max_region_faces, result.n_kept, result.n_relocated,
                result.n_dropped, result.multi_cells_before,
                result.multi_cells_after, result.have_texture,
                result.energy_mean, o.normal_range, o.normal_samples,
                o.radius_gate, o.cross_component_only);
        if (have_crop) {
            fprintf(fp,
                    "  \"crop\":[%.9g,%.9g,%.9g,%.9g,%.9g,%.9g] }\n",
                    bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
        } else {
            fprintf(fp, "  \"crop\":null }\n");
        }
        fclose(fp);
        fprintf(stderr, "  wrote %s\n", path);
    }

    free(m.v);
    free(m.t);
    free(m.f);
    Arena_dispose(&arena);
    return 0;
}
