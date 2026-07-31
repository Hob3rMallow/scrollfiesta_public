/*
 * tifxyz_render.c -- render a VC3D tifxyz atlas the way a consumer sees it:
 * read each segment's exported x/y/z.tif POSITIONS, estimate the surface
 * normal from the position field, sample the RAW CT volume (normal-max) at
 * those positions, and composite every segment onto one canvas at its
 * origin_uv. The output is the unrolled papyrus reconstructed ENTIRELY from
 * the exported segment files -- so it validates the export round-trip and the
 * cross-segment co-registration, not just the in-memory strip.
 *
 * usage: tifxyz_render <atlas_root> --raw <cubes_RAW> --out <stitch.png>
 *          [--down N] [--range F] [--nsteps N] [--tiles DIR]
 *        <atlas_root> holds seg/<uuid>/{x,y,z,mask}.tif + meta.json
 */
#include "../common/ves_platform.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
  #include <windows.h>
#else
  #include <dirent.h>
#endif

#include "../common/arena.h"
#include "../common/tiff_io.h"
#include "../common/raw_sample.h"
#include "../common/ves_png.h"

/* crude key scan (same convention as the placed_index readers). */
static const char *jfind(const char *s, const char *key)
{
    char pat[128];
    snprintf(pat, sizeof(pat), "\"%s\"", key);
    const char *p = strstr(s, pat);
    if (p == NULL) return NULL;
    p = strchr(p + strlen(pat), ':');
    return p != NULL ? p + 1 : NULL;
}

typedef struct {
    char name[128]; double ou, ov, du, dv;
    double blo[3], bhi[3];   /* world (x,y,z) bbox from meta */
    int have_bbox;
} SegMeta;

static int read_seg_meta(const char *seg_dir, SegMeta *m)
{
    char path[1200];
    snprintf(path, sizeof(path), "%s/meta.json", seg_dir);
    FILE *f = fopen(path, "rb");
    if (f == NULL) return -1;
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);
    buf[n] = '\0';
    const char *p;
    m->du = m->dv = 1.0; m->ou = m->ov = 0.0; m->have_bbox = 0;
    if ((p = jfind(buf, "du")) != NULL) sscanf(p, " %lf", &m->du);
    if ((p = jfind(buf, "dv")) != NULL) sscanf(p, " %lf", &m->dv);
    if ((p = jfind(buf, "origin_uv")) != NULL) sscanf(p, " [ %lf , %lf", &m->ou, &m->ov);
    if ((p = jfind(buf, "bbox")) != NULL
        && sscanf(p, " [ [ %lf , %lf , %lf ] , [ %lf , %lf , %lf ]",
                  &m->blo[0], &m->blo[1], &m->blo[2],
                  &m->bhi[0], &m->bhi[1], &m->bhi[2]) == 6)
        m->have_bbox = 1;
    return 0;
}

static int cmp_f(const void *x, const void *y)
{
    float fx = *(const float *)x, fy = *(const float *)y;
    return fx < fy ? -1 : (fx > fy ? 1 : 0);
}

/* percentile contrast stretch of the accumulated intensity canvas. */
static void stretch_to_u8(const float *val, const uint8_t *cov, uint8_t *out,
                          size_t n)
{
    size_t m = 0;
    float *tmp = (float *)malloc((n ? n : 1) * sizeof(float));
    for (size_t i = 0; i < n; i++) if (cov[i]) tmp[m++] = val[i];
    if (m < 2) { for (size_t i = 0; i < n; i++) out[i] = 0; free(tmp); return; }
    qsort(tmp, m, sizeof(float), cmp_f);
    float lo = tmp[(size_t)(0.01 * m)], hi = tmp[(size_t)(0.99 * m)];
    if (hi <= lo) hi = lo + 1.0f;
    for (size_t i = 0; i < n; i++) {
        if (!cov[i]) { out[i] = 0; continue; }
        float t = (val[i] - lo) / (hi - lo);
        if (t < 0) t = 0; if (t > 1) t = 1;
        out[i] = (uint8_t)(t * 255.0f + 0.5f);
    }
    free(tmp);
}

/* Derive "<base>_strips.png" from the out path (insert before a trailing .png). */
static void strips_path(const char *out_png, char *dst, size_t dstsz)
{
    const char *dot = strrchr(out_png, '.');
    if (dot != NULL && (strcmp(dot, ".png") == 0 || strcmp(dot, ".PNG") == 0))
        snprintf(dst, dstsz, "%.*s_strips.png", (int)(dot - out_png), out_png);
    else
        snprintf(dst, dstsz, "%s_strips.png", out_png);
}

/* Re-lay the wide CWxCH stitch into `nstrips` stacked row-slices (a scroll ->
 * book-page layout) so an extreme-aspect ribbon is actually viewable. Each
 * slice is CW/nstrips wide and stacks below the previous with a thin separator.
 * nstrips<=0 => auto (target a roughly square result). Writes <out>_strips.png. */
static void write_strips_png(const char *out_png, const uint8_t *img,
                             int CW, int CH, int nstrips)
{
    if (CW < 1 || CH < 1) return;
    if (nstrips <= 0) {
        nstrips = (int)(sqrt((double)CW / (double)CH) + 0.5);
        if (nstrips < 1) nstrips = 1;
        if (nstrips > 64) nstrips = 64;
    }
    int sw = (CW + nstrips - 1) / nstrips;         /* strip width */
    const int gap = 6;                             /* separator rows */
    const uint8_t sep = 40;
    long OW = sw;
    long OH = (long)nstrips * (long)CH + (long)(nstrips - 1) * (long)gap;
    if (OW < 1 || OH < 1 || (double)OW * (double)OH > 500000000.0) {
        fprintf(stderr, "tifxyz_render: strips %ldx%ld too large -- skipped "
                "(raise --down)\n", OW, OH);
        return;
    }
    size_t on = (size_t)OW * (size_t)OH;
    uint8_t *out = (uint8_t *)malloc(on);
    if (out == NULL) return;
    memset(out, 0, on);
    for (int k = 0; k < nstrips; k++) {
        long dr0 = (long)k * (long)(CH + gap);
        int c0 = k * sw;
        int cw = (c0 + sw <= CW) ? sw : (CW - c0);
        for (int r = 0; r < CH; r++) {
            uint8_t *drow = out + (size_t)(dr0 + r) * (size_t)OW;
            const uint8_t *srow = img + (size_t)r * (size_t)CW + (size_t)c0;
            if (cw > 0) memcpy(drow, srow, (size_t)cw);   /* right pad stays 0 */
        }
        if (k < nstrips - 1)
            for (int g = 0; g < gap; g++)
                memset(out + (size_t)(dr0 + CH + g) * (size_t)OW, sep, (size_t)OW);
    }
    char sp[1400];
    strips_path(out_png, sp, sizeof sp);
    if (VesPng_write_gray(sp, out, (int)OW, (int)OH) != 0)
        fprintf(stderr, "tifxyz_render: write %s FAILED\n", sp);
    else
        fprintf(stderr, "tifxyz_render: wrote %s (%ldx%ld, %d strips)\n",
                sp, OW, OH, nstrips);
    free(out);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: tifxyz_render <atlas_root> --raw <dir> "
                "--out <png> [--down N] [--range F] [--nsteps N] [--tiles DIR]\n"
                "         [--strip-rows N]  (also write <out>_strips.png: the\n"
                "         wide stitch re-laid as N stacked row-slices; 0=auto)\n");
        return 1;
    }
    const char *root = argv[1];
    const char *raw_dir = NULL, *out_png = NULL, *tiles_dir = NULL;
    int down = 4, nsteps = 7, strip_rows = 0;
    double range = 3.0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--raw") && i + 1 < argc) raw_dir = argv[++i];
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) out_png = argv[++i];
        else if (!strcmp(argv[i], "--down") && i + 1 < argc) down = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--range") && i + 1 < argc) range = atof(argv[++i]);
        else if (!strcmp(argv[i], "--nsteps") && i + 1 < argc) nsteps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--tiles") && i + 1 < argc) tiles_dir = argv[++i];
        else if (!strcmp(argv[i], "--strip-rows") && i + 1 < argc) strip_rows = atoi(argv[++i]);
        else { fprintf(stderr, "unknown arg %s\n", argv[i]); return 1; }
    }
    if (!raw_dir || !out_png) { fprintf(stderr, "need --raw and --out\n"); return 1; }
    if (down < 1) down = 1;

    char segroot[1024];
    snprintf(segroot, sizeof(segroot), "%s/seg", root);

    /* enumerate segment dirs. */
    char (*names)[128] = NULL; size_t nseg = 0, cap = 0;
#ifdef _WIN32
    char glob[1200]; snprintf(glob, sizeof(glob), "%s\\*", segroot);
    WIN32_FIND_DATAA fd; HANDLE h = FindFirstFileA(glob, &fd);
    if (h != INVALID_HANDLE_VALUE) {
        do {
            if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) continue;
            if (fd.cFileName[0] == '.') continue;
            if (nseg == cap) { cap = cap ? cap * 2 : 256;
                names = (char (*)[128])realloc(names, cap * 128); }
            snprintf(names[nseg++], 128, "%s", fd.cFileName);
        } while (FindNextFileA(h, &fd));
        FindClose(h);
    }
#else
    DIR *d = opendir(segroot);
    if (d) { struct dirent *e;
        while ((e = readdir(d))) {
            if (e->d_name[0] == '.') continue;
            if (nseg == cap) { cap = cap ? cap * 2 : 256;
                names = (char (*)[128])realloc(names, cap * 128); }
            snprintf(names[nseg++], 128, "%s", e->d_name);
        }
        closedir(d);
    }
#endif
    if (nseg == 0) { fprintf(stderr, "no segments under %s\n", segroot); return 1; }
    fprintf(stderr, "tifxyz_render: %zu segments under %s\n", nseg, segroot);

    /* pass 1: read metas, find the global (u,v) extent AND the global world
     * bbox (for one shared CubeTable prewarmed once -- per-segment reload of
     * overlapping cubes across a 1702-cube grid is what made this I/O bound). */
    SegMeta *metas = (SegMeta *)calloc(nseg, sizeof(SegMeta));
    double gu0 = 1e300, gv0 = 1e300, gu1 = -1e300, gv1 = -1e300;
    /* native (z,y,x) bbox accumulated from meta world (x,y,z) bboxes. */
    double bz0 = 1e300, by0 = 1e300, bx0 = 1e300;
    double bz1 = -1e300, by1 = -1e300, bx1 = -1e300;
    for (size_t s = 0; s < nseg; s++) {
        char sd[1200]; snprintf(sd, sizeof(sd), "%s/%s", segroot, names[s]);
        snprintf(metas[s].name, sizeof(metas[s].name), "%s", names[s]);
        if (read_seg_meta(sd, &metas[s]) != 0) continue;
        if (metas[s].ou < gu0) gu0 = metas[s].ou;
        if (metas[s].ov < gv0) gv0 = metas[s].ov;
        if (metas[s].have_bbox) {   /* world (x,y,z) -> native (z,y,x) */
            if (metas[s].blo[2] < bz0) bz0 = metas[s].blo[2];
            if (metas[s].blo[1] < by0) by0 = metas[s].blo[1];
            if (metas[s].blo[0] < bx0) bx0 = metas[s].blo[0];
            if (metas[s].bhi[2] > bz1) bz1 = metas[s].bhi[2];
            if (metas[s].bhi[1] > by1) by1 = metas[s].bhi[1];
            if (metas[s].bhi[0] > bx1) bx1 = metas[s].bhi[0];
        }
    }

    /* one shared CubeTable over the global bbox, prewarmed once. */
    Arena_T ct_arena = Arena_new();
    CubeTable ct;
    int have_ct = 0;
    if (bz1 > bz0) {
        float corners[6] = { (float)bz0,(float)by0,(float)bx0,
                             (float)bz1,(float)by1,(float)bx1 };
        if (cubetable_init(&ct, ct_arena, raw_dir, 128, corners, 2, range + 3.0) == 0) {
            double tp = ves_clock_sec();
            int nl = cubetable_prewarm_all(&ct);
            fprintf(stderr, "tifxyz_render: prewarmed %d cubes (%.1fs)\n",
                    nl, ves_clock_sec() - tp);
            have_ct = 1;
        }
    }
    if (!have_ct) { fprintf(stderr, "tifxyz_render: no bbox/RAW; cannot render\n"); return 1; }

    /* pass 2: load each segment, sample RAW, composite (downsampled). We size
     * the canvas lazily from the running max of origin+W. Two-phase: first
     * load all, record W/H + rendered tiles, then composite. */
    Arena_T ar = Arena_new();
    /* store per-seg rendered intensity + coverage + placement */
    typedef struct { float *val; uint8_t *cov; int W, H; double ou, ov; } Rend;
    Rend *R = (Rend *)calloc(nseg, sizeof(Rend));

    for (size_t s = 0; s < nseg; s++) {
        Arena_free(ar);
        char sd[1200]; snprintf(sd, sizeof(sd), "%s/%s", segroot, names[s]);
        char px[1300], py[1300], pz[1300];
        snprintf(px, sizeof px, "%s/x.tif", sd);
        snprintf(py, sizeof py, "%s/y.tif", sd);
        snprintf(pz, sizeof pz, "%s/z.tif", sd);
        float *X = NULL, *Y = NULL, *Z = NULL; int W = 0, H = 0, W2, H2;
        if (TiffIO_load_float2d(ar, px, &X, &W, &H) != 0) continue;
        if (TiffIO_load_float2d(ar, py, &Y, &W2, &H2) != 0 || W2 != W || H2 != H) continue;
        if (TiffIO_load_float2d(ar, pz, &Z, &W2, &H2) != 0 || W2 != W || H2 != H) continue;

        size_t np = (size_t)W * (size_t)H;
        /* native (z,y,x) verts for the sampler + bbox. */
        float *P = (float *)ARENA_ALLOC(ar, (long)(np * 3 * sizeof(float)));
        uint8_t *valid = (uint8_t *)ARENA_ALLOC(ar, (long)np);
        for (size_t i = 0; i < np; i++) {
            float xx = X[i], yy = Y[i], zz = Z[i];
            valid[i] = (xx > 0.0f && zz > 0.0f && !(xx != xx)) ? 1 : 0;   /* z<=0 dropped */
            P[i*3+0] = zz; P[i*3+1] = yy; P[i*3+2] = xx;
        }

        float *val = (float *)malloc(np * sizeof(float));
        uint8_t *cov = (uint8_t *)malloc(np);
        for (int r = 0; r < H; r++) {
            for (int c = 0; c < W; c++) {
                size_t i = (size_t)r * (size_t)W + (size_t)c;
                if (!valid[i]) { cov[i] = 0; val[i] = 0; continue; }
                /* normal from position gradients (central diff, valid nbrs). */
                double n[3] = {0,0,0}, nn = 0.0;
                if (c > 0 && c + 1 < W && r > 0 && r + 1 < H) {
                    size_t il=i-1, ir=i+1, iu=i-(size_t)W, id=i+(size_t)W;
                    if (valid[il]&&valid[ir]&&valid[iu]&&valid[id]) {
                        double dc[3]={P[ir*3]-P[il*3],P[ir*3+1]-P[il*3+1],P[ir*3+2]-P[il*3+2]};
                        double dr[3]={P[id*3]-P[iu*3],P[id*3+1]-P[iu*3+1],P[id*3+2]-P[iu*3+2]};
                        n[0]=dc[1]*dr[2]-dc[2]*dr[1];
                        n[1]=dc[2]*dr[0]-dc[0]*dr[2];
                        n[2]=dc[0]*dr[1]-dc[1]*dr[0];
                        nn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
                    }
                }
                float nf[3];
                if (nn > 1e-6) { nf[0]=(float)(n[0]/nn); nf[1]=(float)(n[1]/nn); nf[2]=(float)(n[2]/nn); }
                else { nf[0]=nf[1]=nf[2]=0.0f; }
                double v = sample_vertex(&ct, &P[i*3], nf, range, nsteps);
                if (v < 0) { cov[i]=0; val[i]=0; }
                else { cov[i]=1; val[i]=(float)v; }
            }
        }
        R[s].val = val; R[s].cov = cov; R[s].W = W; R[s].H = H;
        R[s].ou = metas[s].ou; R[s].ov = metas[s].ov;
        if (metas[s].ou + W > gu1) gu1 = metas[s].ou + W;
        if (metas[s].ov + H > gv1) gv1 = metas[s].ov + H;

        /* optional full-res per-seg tile. */
        if (tiles_dir) {
            uint8_t *g = (uint8_t *)malloc(np);
            stretch_to_u8(val, cov, g, np);
            char tp[1400]; snprintf(tp, sizeof tp, "%s/%s.png", tiles_dir, names[s]);
            ves_ensure_parent_dir(tp);
            VesPng_write_gray(tp, g, W, H);
            free(g);
        }
        if ((s % 10) == 0) fprintf(stderr, "  rendered %zu/%zu\n", s + 1, nseg);
    }

    /* composite onto the downsampled global canvas (max over covers). */
    int CW = (int)((gu1 - gu0) / down) + 1;
    int CH = (int)((gv1 - gv0) / down) + 1;
    if (CW < 1) CW = 1; if (CH < 1) CH = 1;
    fprintf(stderr, "tifxyz_render: canvas %dx%d (down %d) from u[%.0f,%.0f] v[%.0f,%.0f]\n",
            CW, CH, down, gu0, gu1, gv0, gv1);
    size_t cn = (size_t)CW * (size_t)CH;
    float *cval = (float *)calloc(cn, sizeof(float));
    uint8_t *ccov = (uint8_t *)calloc(cn, 1);
    for (size_t s = 0; s < nseg; s++) {
        if (!R[s].val) continue;
        for (int r = 0; r < R[s].H; r++) {
            for (int c = 0; c < R[s].W; c++) {
                size_t i = (size_t)r * (size_t)R[s].W + (size_t)c;
                if (!R[s].cov[i]) continue;
                int cc = (int)((R[s].ou + c - gu0) / down);
                int cr = (int)((R[s].ov + r - gv0) / down);
                if (cc < 0 || cc >= CW || cr < 0 || cr >= CH) continue;
                size_t ci = (size_t)cr * (size_t)CW + (size_t)cc;
                if (!ccov[ci] || R[s].val[i] > cval[ci]) { cval[ci] = R[s].val[i]; ccov[ci] = 1; }
            }
        }
        free(R[s].val); free(R[s].cov);
    }
    uint8_t *img = (uint8_t *)malloc(cn);
    stretch_to_u8(cval, ccov, img, cn);
    ves_ensure_parent_dir(out_png);
    if (VesPng_write_gray(out_png, img, CW, CH) != 0)
        fprintf(stderr, "tifxyz_render: write %s FAILED\n", out_png);
    else
        fprintf(stderr, "tifxyz_render: wrote %s (%dx%d)\n", out_png, CW, CH);

    /* also emit the strips layout so the wide ribbon is viewable. */
    write_strips_png(out_png, img, CW, CH, strip_rows);

    free(img); free(cval); free(ccov); free(R); free(metas); free(names);
    Arena_dispose(&ar);
    Arena_dispose(&ct_arena);
    return 0;
}
