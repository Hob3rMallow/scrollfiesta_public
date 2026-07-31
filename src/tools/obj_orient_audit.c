/* ============================================================================
 * obj_orient_audit -- per-component orientation & UV-mirror audit for a
 * scroll-surface OBJ (optionally carrying per-vertex vt from an unroll stage).
 *
 * Question it answers: when components "face different directions", is that
 *   (a) stored-winding sign flips (some components present recto, some verso
 *       -- rendering/backface metadata, geometry fine), or
 *   (b) genuinely mirrored charts in the UV unroll (text would read
 *       backwards -- a parameterization bug)?
 *
 * Three independent measurements per connected component:
 *
 *   1. RADIAL sign  s_rad(f) = sign(n_f . r_hat)   (n_f from stored winding,
 *      r_hat = radial direction from the scroll axis at the face centroid).
 *      A wrap fragment consistently wound has one dominant sign: which
 *      physical side of the sheet the stored winding calls "front".
 *
 *   2. same_dir     count of directed edges traversed twice in the SAME
 *      direction = internal winding inconsistencies (0 on a consistently
 *      wound component; >0 means mixed winding INSIDE the component).
 *
 *   3. MIRROR invariant  kappa(f) = sign(uv_area(f)) * s_rad(f).
 *      Flipping a face's stored winding flips BOTH factors, so kappa is
 *      INVARIANT under winding flips -- it measures only the geometry of the
 *      3D->UV map. All non-mirrored charts of one unroll share a single
 *      global kappa; a component with the opposite kappa is mirror-imaged
 *      in UV.
 *
 * Verdict per component vs the largest (main) component:
 *   WINDING-FLIPPED  radial sign opposite to main (agreement gated)
 *   MIRRORED         kappa opposite to main       (agreement gated)
 *   MIXED            internal agreement < 80% (radial) -- fold/core noise or
 *                    genuinely mixed winding (check same_dir)
 *
 * Optional dumps (flat (u,v,0) OBJs, per-vertex color):
 *   <id>_orient_flat.obj  green = radial sign agrees with main, red =
 *                         opposite, gray = ambiguous/edge-on
 *   <id>_mirror_flat.obj  green = kappa agrees with main, magenta = mirrored,
 *                         gray = ambiguous
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"

/* ------------------------------------------------------------------ util */

static void *xmalloc(size_t nbytes)
{
    void *p = malloc(nbytes > 0 ? nbytes : 1);
    if (p == NULL) {
        fprintf(stderr, "ERROR: out of memory (%zu bytes)\n", nbytes);
        exit(1);
    }
    return p;
}

static void *xcalloc(size_t count, size_t size)
{
    void *p = calloc(count > 0 ? count : 1, size);
    if (p == NULL) {
        fprintf(stderr, "ERROR: out of memory (%zu x %zu)\n", count, size);
        exit(1);
    }
    return p;
}

static void *xgrow(void *p, size_t need_elems, size_t *cap_elems, size_t elem)
{
    size_t ncap = *cap_elems;
    if (need_elems <= ncap) return p;
    ncap = (ncap == 0) ? 4096 : ncap;
    while (ncap < need_elems) ncap *= 2;
    p = realloc(p, ncap * elem);
    if (p == NULL) {
        fprintf(stderr, "ERROR: out of memory (%zu bytes)\n", ncap * elem);
        exit(1);
    }
    *cap_elems = ncap;
    return p;
}

/* -------------------------------------------------- OBJ with vt parsing */

typedef struct UvMesh {
    float   *verts;   /* [nv*3] (z,y,x) as stored */
    float   *uv;      /* [nvt*2] */
    int32_t *faces;   /* [nf*3] 0-based */
    size_t   nv, nvt, nf;
} UvMesh;

static int parse_uv_obj(const char *path, UvMesh *m)
{
    FILE *f = NULL;
    char line[1024];
    size_t cap_v = 0, cap_t = 0, cap_f = 0;

    memset(m, 0, sizeof *m);
    f = fopen(path, "rb");
    if (f == NULL) {
        fprintf(stderr, "ERROR: cannot open %s\n", path);
        return -1;
    }
    while (fgets(line, sizeof line, f) != NULL) {
        if (line[0] == 'v' && line[1] == ' ') {
            char *s = line + 2, *end = NULL;
            double a = strtod(s, &end); s = end;
            double b = strtod(s, &end); s = end;
            double c = strtod(s, &end);
            m->verts = (float *)xgrow(m->verts, (m->nv + 1) * 3, &cap_v,
                                      sizeof(float));
            m->verts[m->nv * 3 + 0] = (float)a;
            m->verts[m->nv * 3 + 1] = (float)b;
            m->verts[m->nv * 3 + 2] = (float)c;
            m->nv++;
        } else if (line[0] == 'v' && line[1] == 't' && line[2] == ' ') {
            char *s = line + 3, *end = NULL;
            double u = strtod(s, &end); s = end;
            double v = strtod(s, &end);
            m->uv = (float *)xgrow(m->uv, (m->nvt + 1) * 2, &cap_t,
                                   sizeof(float));
            m->uv[m->nvt * 2 + 0] = (float)u;
            m->uv[m->nvt * 2 + 1] = (float)v;
            m->nvt++;
        } else if (line[0] == 'f' && (line[1] == ' ' || line[1] == '\t')) {
            char *s = line + 2;
            long idx[16];
            int nidx = 0, k = 0;
            while (nidx < 16) {
                char *end = NULL;
                long v = 0;
                while (*s == ' ' || *s == '\t') s++;
                if (*s == '\0' || *s == '\n' || *s == '\r') break;
                v = strtol(s, &end, 10);
                if (end == s) break;
                idx[nidx++] = v;
                s = end;
                while (*s != '\0' && *s != ' ' && *s != '\t'
                       && *s != '\n' && *s != '\r') s++;
            }
            for (k = 2; k < nidx; k++) {
                m->faces = (int32_t *)xgrow(m->faces, (m->nf + 1) * 3, &cap_f,
                                            sizeof(int32_t));
                m->faces[m->nf * 3 + 0] = (int32_t)(idx[0] - 1);
                m->faces[m->nf * 3 + 1] = (int32_t)(idx[k - 1] - 1);
                m->faces[m->nf * 3 + 2] = (int32_t)(idx[k] - 1);
                m->nf++;
            }
        }
    }
    fclose(f);
    if (m->nv == 0) {
        fprintf(stderr, "ERROR: %s has no vertices\n", path);
        return -1;
    }
    {
        size_t t = 0, bad = 0;
        for (t = 0; t < m->nf * 3; t++)
            if (m->faces[t] < 0 || (size_t)m->faces[t] >= m->nv) bad++;
        if (bad > 0) {
            fprintf(stderr, "ERROR: %zu face indices out of range\n", bad);
            return -1;
        }
    }
    return 0;
}

static void uvmesh_free(UvMesh *m)
{
    free(m->verts); free(m->uv); free(m->faces);
    memset(m, 0, sizeof *m);
}

/* ---------------------------------------------------------------- audit */

typedef struct CompStat {
    size_t  nf;
    double  area;               /* total 3D area */
    double  area_out, area_in;  /* radial-clear area by stored-winding sign */
    double  kpos, kneg;         /* kappa votes, area-weighted */
    double  coh;                /* area-weighted mean |cos(n, r_hat)| */
    double  dot111;             /* area-weighted mean n_hat . (1,1,1)/sqrt(3):
                                 * the per-cube MLS anchor's vote direction */
    double  az_y, az_x;         /* area-weighted mean in-plane radial dir
                                 * (unit r_hat summed) -> mean azimuth */
    double  umin, umax, vmin, vmax;
    size_t  same_dir;           /* directed-edge winding violations */
    int32_t root;               /* UF root vertex id */
} CompStat;

typedef struct AuditResult {
    CompStat *comp;      /* sorted by area desc */
    size_t    ncomp;
    int       R0, K0;    /* main component's radial / kappa signs */
    /* per-face data for dumps (indices into sorted comp[]) */
    int32_t  *face_comp;
    int8_t   *face_rad;  /* -1/0/+1 */
    int8_t   *face_kap;  /* -1/0/+1 */
    /* per-vertex comp (sorted idx) */
    int32_t  *vert_comp;
} AuditResult;

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static int cmp_comp_area(const void *a, const void *b)
{
    const CompStat *x = (const CompStat *)a, *y = (const CompStat *)b;
    return (x->area > y->area) ? -1 : (x->area < y->area) ? 1 : 0;
}

/* Core analysis on in-memory arrays (uv may be NULL: radial-only audit).
 * axis_point/axis_dir in the same (z,y,x) component order as verts. */
static int audit_run(Arena_T arena,
                     const float *verts, size_t nv, const float *uv,
                     const int32_t *faces, size_t nf,
                     const float axis_point[3], const float axis_dir[3],
                     double min_cos, AuditResult *R)
{
    UnionFind ufv;
    int32_t *root2dense = NULL;
    int32_t *comp_of_vert = NULL;
    CompStat *cs = NULL;
    size_t ncomp = 0, f = 0, i = 0;
    double ad[3] = { 0.0, 0.0, 0.0 };
    double alen = 0.0;

    memset(R, 0, sizeof *R);
    if (nv == 0 || nf == 0 || nv > (size_t)INT32_MAX) return -1;

    alen = sqrt((double)axis_dir[0] * axis_dir[0]
                + (double)axis_dir[1] * axis_dir[1]
                + (double)axis_dir[2] * axis_dir[2]);
    if (alen < 1e-12) return -1;
    for (i = 0; i < 3; i++) ad[i] = (double)axis_dir[i] / alen;

    /* components over shared vertices */
    ufv = UF_new(arena, (int32_t)nv);
    for (f = 0; f < nf; f++) {
        uf_union(&ufv, faces[f * 3 + 0], faces[f * 3 + 1]);
        uf_union(&ufv, faces[f * 3 + 0], faces[f * 3 + 2]);
    }
    root2dense = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
    comp_of_vert = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
    for (i = 0; i < nv; i++) root2dense[i] = -1;
    /* count components that actually own faces */
    for (f = 0; f < nf; f++) {
        int32_t r = uf_find(&ufv, faces[f * 3 + 0]);
        if (root2dense[r] < 0) root2dense[r] = (int32_t)ncomp++;
    }
    if (ncomp == 0) return -1;
    cs = (CompStat *)ARENA_CALLOC(arena, (long)ncomp, (long)sizeof(CompStat));
    for (i = 0; i < ncomp; i++) {
        cs[i].umin = 1e30; cs[i].umax = -1e30;
        cs[i].vmin = 1e30; cs[i].vmax = -1e30;
    }
    for (i = 0; i < nv; i++) {
        int32_t r = uf_find(&ufv, (int32_t)i);
        comp_of_vert[i] = root2dense[r];   /* -1 for isolated verts */
        if (root2dense[r] >= 0) cs[root2dense[r]].root = r;
    }

    R->face_comp = (int32_t *)ARENA_ALLOC(arena, (long)(nf * sizeof(int32_t)));
    R->face_rad = (int8_t *)ARENA_CALLOC(arena, (long)nf, 1);
    R->face_kap = (int8_t *)ARENA_CALLOC(arena, (long)nf, 1);

    /* per-face radial + kappa votes */
    for (f = 0; f < nf; f++) {
        size_t a = (size_t)faces[f * 3 + 0];
        size_t b = (size_t)faces[f * 3 + 1];
        size_t c = (size_t)faces[f * 3 + 2];
        int32_t ci = comp_of_vert[a];
        CompStat *S = &cs[ci];
        double e1[3], e2[3], gn[3], cen[3], rv[3];
        double area = 0.0, dax = 0.0, rlen = 0.0, cosr = 0.0;
        int k = 0;

        R->face_comp[f] = ci;
        S->nf++;
        for (k = 0; k < 3; k++) {
            e1[k] = (double)verts[b * 3 + k] - (double)verts[a * 3 + k];
            e2[k] = (double)verts[c * 3 + k] - (double)verts[a * 3 + k];
            cen[k] = ((double)verts[a * 3 + k] + (double)verts[b * 3 + k]
                      + (double)verts[c * 3 + k]) / 3.0
                     - (double)axis_point[k];
        }
        gn[0] = e1[1] * e2[2] - e1[2] * e2[1];
        gn[1] = e1[2] * e2[0] - e1[0] * e2[2];
        gn[2] = e1[0] * e2[1] - e1[1] * e2[0];
        area = 0.5 * sqrt(gn[0] * gn[0] + gn[1] * gn[1] + gn[2] * gn[2]);
        S->area += area;
        if (area < 1e-12) continue;

        dax = cen[0] * ad[0] + cen[1] * ad[1] + cen[2] * ad[2];
        for (k = 0; k < 3; k++) rv[k] = cen[k] - dax * ad[k];
        rlen = sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]);
        if (rlen < 1.0) continue;   /* on-axis: radial dir undefined */
        cosr = (gn[0] * rv[0] + gn[1] * rv[1] + gn[2] * rv[2])
               / (2.0 * area * rlen);
        S->coh += area * fabs(cosr);
        S->dot111 += area * (gn[0] + gn[1] + gn[2])
                     / (2.0 * area * 1.7320508075688772);
        S->az_y += area * rv[1] / rlen;
        S->az_x += area * rv[2] / rlen;
        if (cosr >= min_cos) { S->area_out += area; R->face_rad[f] = 1; }
        else if (cosr <= -min_cos) { S->area_in += area; R->face_rad[f] = -1; }

        if (uv != NULL && R->face_rad[f] != 0) {
            double ua = (double)uv[a * 2 + 0], va = (double)uv[a * 2 + 1];
            double ub = (double)uv[b * 2 + 0], vb = (double)uv[b * 2 + 1];
            double uc = (double)uv[c * 2 + 0], vc = (double)uv[c * 2 + 1];
            double A2 = (ub - ua) * (vc - va) - (vb - va) * (uc - ua);
            if (fabs(A2) > 1e-9) {
                int kap = ((A2 > 0.0) == (cosr > 0.0)) ? 1 : -1;
                R->face_kap[f] = (int8_t)kap;
                if (kap > 0) S->kpos += area; else S->kneg += area;
            }
        }
    }

    /* per-comp u/v ranges */
    if (uv != NULL) {
        for (i = 0; i < nv; i++) {
            int32_t ci = comp_of_vert[i];
            double u = 0.0, v = 0.0;
            if (ci < 0) continue;
            u = (double)uv[i * 2 + 0]; v = (double)uv[i * 2 + 1];
            if (u < cs[ci].umin) cs[ci].umin = u;
            if (u > cs[ci].umax) cs[ci].umax = u;
            if (v < cs[ci].vmin) cs[ci].vmin = v;
            if (v > cs[ci].vmax) cs[ci].vmax = v;
        }
    }

    /* same-direction directed-edge violations, attributed per comp */
    {
        uint64_t *keys = (uint64_t *)xmalloc(nf * 3 * sizeof(uint64_t));
        size_t n = 0;
        for (f = 0; f < nf; f++) {
            int32_t v3[3] = { faces[f * 3 + 0], faces[f * 3 + 1],
                              faces[f * 3 + 2] };
            int e = 0;
            for (e = 0; e < 3; e++) {
                uint64_t aa = (uint64_t)(uint32_t)v3[e];
                uint64_t bb = (uint64_t)(uint32_t)v3[(e + 1) % 3];
                keys[n++] = (aa << 32) | bb;
            }
        }
        qsort(keys, n, sizeof(uint64_t), cmp_u64);
        for (i = 0; i + 1 < n; i++) {
            if (keys[i] == keys[i + 1]) {
                int32_t a = (int32_t)(keys[i] >> 32);
                if (comp_of_vert[a] >= 0) cs[comp_of_vert[a]].same_dir++;
                while (i + 1 < n && keys[i] == keys[i + 1]) i++;
            }
        }
        free(keys);
    }

    /* sort by area desc; remap face_comp / vert_comp to sorted order */
    {
        int32_t *old2new = (int32_t *)ARENA_ALLOC(arena,
                                                  (long)(ncomp * sizeof(int32_t)));
        CompStat *sorted = (CompStat *)ARENA_ALLOC(arena,
                                                   (long)(ncomp * sizeof(CompStat)));
        memcpy(sorted, cs, ncomp * sizeof(CompStat));
        qsort(sorted, ncomp, sizeof(CompStat), cmp_comp_area);
        /* old->new by matching root ids (unique per comp); O(n^2) is fine
         * for the tens-to-hundreds of comps this tool sees */
        for (i = 0; i < ncomp; i++) old2new[i] = -1;
        for (i = 0; i < ncomp; i++) {
            size_t j = 0;
            for (j = 0; j < ncomp; j++)
                if (sorted[j].root == cs[i].root) { old2new[i] = (int32_t)j; break; }
        }
        for (f = 0; f < nf; f++) R->face_comp[f] = old2new[R->face_comp[f]];
        R->vert_comp = (int32_t *)ARENA_ALLOC(arena,
                                              (long)(nv * sizeof(int32_t)));
        for (i = 0; i < nv; i++)
            R->vert_comp[i] = (comp_of_vert[i] >= 0)
                              ? old2new[comp_of_vert[i]] : -1;
        R->comp = sorted;
        R->ncomp = ncomp;
    }

    R->R0 = (R->comp[0].area_out >= R->comp[0].area_in) ? 1 : -1;
    R->K0 = (R->comp[0].kpos >= R->comp[0].kneg) ? 1 : -1;
    return 0;
}

/* verdict helpers */
static int comp_rad_sign(const CompStat *c)
{ return (c->area_out >= c->area_in) ? 1 : -1; }
static double comp_rad_agr(const CompStat *c)
{
    double s = c->area_out + c->area_in;
    return (s > 0.0) ? ((c->area_out > c->area_in ? c->area_out : c->area_in) / s)
                     : 0.0;
}
static int comp_kap_sign(const CompStat *c)
{ return (c->kpos >= c->kneg) ? 1 : -1; }
static double comp_kap_agr(const CompStat *c)
{
    double s = c->kpos + c->kneg;
    return (s > 0.0) ? ((c->kpos > c->kneg ? c->kpos : c->kneg) / s) : 0.0;
}

/* ------------------------------------------------------------- dumps */

/* Flat (u,v,0) OBJ colored by per-vertex area-weighted vote. mode 0 =
 * radial-vs-main (green agree / red opposite / gray unclear), mode 1 =
 * kappa-vs-main (green / magenta / gray). */
static int write_vote_flat(const char *path,
                           const float *uv, size_t nv,
                           const int32_t *faces, size_t nf,
                           const float *verts3d,
                           const AuditResult *R, int mode)
{
    float *acc = (float *)xcalloc(nv, sizeof(float));
    float *wsum = (float *)xcalloc(nv, sizeof(float));
    FILE *fp = NULL;
    size_t f = 0, i = 0;
    int ref = (mode == 0) ? R->R0 : R->K0;

    for (f = 0; f < nf; f++) {
        int8_t s = (mode == 0) ? R->face_rad[f] : R->face_kap[f];
        size_t a = (size_t)faces[f * 3 + 0], b = (size_t)faces[f * 3 + 1],
               c = (size_t)faces[f * 3 + 2];
        double e1[3], e2[3], gn[3];
        float area = 0.0f;
        int k = 0;
        for (k = 0; k < 3; k++) {
            e1[k] = (double)verts3d[b * 3 + k] - (double)verts3d[a * 3 + k];
            e2[k] = (double)verts3d[c * 3 + k] - (double)verts3d[a * 3 + k];
        }
        gn[0] = e1[1] * e2[2] - e1[2] * e2[1];
        gn[1] = e1[2] * e2[0] - e1[0] * e2[2];
        gn[2] = e1[0] * e2[1] - e1[1] * e2[0];
        area = (float)(0.5 * sqrt(gn[0] * gn[0] + gn[1] * gn[1]
                                  + gn[2] * gn[2]));
        acc[a] += (float)(s * ref) * area;
        acc[b] += (float)(s * ref) * area;
        acc[c] += (float)(s * ref) * area;
        wsum[a] += area; wsum[b] += area; wsum[c] += area;
    }
    fp = fopen(path, "wb");
    if (fp == NULL) { free(acc); free(wsum); return -1; }
    for (i = 0; i < nv; i++) {
        double t = (wsum[i] > 1e-12f) ? (double)acc[i] / (double)wsum[i] : 0.0;
        double r = 0.55, g = 0.55, bl = 0.55;   /* gray = unclear */
        if (t > 0.3) { r = 0.15; g = 0.80; bl = 0.25; }          /* agree */
        else if (t < -0.3) {
            if (mode == 0) { r = 0.90; g = 0.12; bl = 0.12; }    /* flipped */
            else           { r = 0.85; g = 0.10; bl = 0.85; }    /* mirrored */
        }
        fprintf(fp, "v %.4f %.4f 0 %.3f %.3f %.3f\n",
                (double)uv[i * 2 + 0], (double)uv[i * 2 + 1], r, g, bl);
    }
    for (f = 0; f < nf; f++)
        fprintf(fp, "f %d %d %d\n", faces[f * 3 + 0] + 1,
                faces[f * 3 + 1] + 1, faces[f * 3 + 2] + 1);
    fclose(fp);
    free(acc); free(wsum);
    return 0;
}

/* ------------------------------------------------------------ reporting */

static void print_report(const AuditResult *R, const float *uv, size_t max_rows)
{
    size_t i = 0, n_flip = 0, n_mirror = 0, n_mixed = 0, n_samedir = 0;
    double tot_area = 0.0;

    for (i = 0; i < R->ncomp; i++) tot_area += R->comp[i].area;
    fprintf(stderr,
        "comp   faces      area%%  sameDir  radial       kappa        dot111  azim    u-range          verdict\n");
    for (i = 0; i < R->ncomp && i < max_rows; i++) {
        const CompStat *c = &R->comp[i];
        int rs = comp_rad_sign(c), ks = comp_kap_sign(c);
        double ra = comp_rad_agr(c), ka = comp_kap_agr(c);
        double d111 = (c->area > 0.0) ? c->dot111 / c->area : 0.0;
        double azim = atan2(c->az_x, c->az_y) * 57.29577951308232;
        char verdict[80] = "ok";
        if (ra < 0.8) snprintf(verdict, sizeof verdict, "MIXED(%.0f%%)", ra * 100.0);
        else if (rs != R->R0) snprintf(verdict, sizeof verdict, "WINDING-FLIPPED");
        if (uv != NULL && ka >= 0.8 && ks != R->K0) {
            size_t len = strlen(verdict);
            snprintf(verdict + len, sizeof verdict - len, "%sMIRRORED",
                     (strcmp(verdict, "ok") == 0) ? "" : "+");
            if (strncmp(verdict, "ok", 2) == 0)
                memmove(verdict, verdict + 2, strlen(verdict + 2) + 1);
        }
        fprintf(stderr,
            "%4zu %8zu  %8.3f  %7zu  %s %5.1f%%   %s %5.1f%%   %+6.2f  %+4.0f  [%7.0f,%7.0f]  %s\n",
            i, c->nf, 100.0 * c->area / (tot_area > 0 ? tot_area : 1.0),
            c->same_dir,
            rs > 0 ? "out" : "in ", ra * 100.0,
            uv ? (ks > 0 ? "+" : "-") : "?", uv ? ka * 100.0 : 0.0,
            d111, azim,
            uv ? c->umin : 0.0, uv ? c->umax : 0.0, verdict);
    }
    if (R->ncomp > max_rows)
        fprintf(stderr, "  ... %zu more components\n", R->ncomp - max_rows);
    for (i = 0; i < R->ncomp; i++) {
        const CompStat *c = &R->comp[i];
        double ra = comp_rad_agr(c), ka = comp_kap_agr(c);
        if (ra < 0.8) n_mixed++;
        else if (comp_rad_sign(c) != R->R0) n_flip++;
        if (uv != NULL && ka >= 0.8 && comp_kap_sign(c) != R->K0) n_mirror++;
        if (c->same_dir > 0) n_samedir++;
    }
    fprintf(stderr,
        "[orient_audit] %zu comps: main radial=%s | %zu WINDING-FLIPPED, "
        "%zu MIXED(<80%%), %zu with same_dir>0%s",
        R->ncomp, R->R0 > 0 ? "out" : "in", n_flip, n_mixed, n_samedir,
        uv != NULL ? ", " : "\n");
    if (uv != NULL)
        fprintf(stderr, "%zu MIRRORED (kappa %s)\n",
                n_mirror, R->K0 > 0 ? "+" : "-");
}

/* -------------------------------------------------------------- selftest */

#define CHECK(cond, msg) \
    do { if (!(cond)) { fprintf(stderr, "  FAIL: %s\n", (msg)); nfail++; } \
         else { fprintf(stderr, "  ok: %s\n", (msg)); } } while (0)

/* Cylinder-arc sheet around axis z through (y,x)=(0,0): radius rad,
 * theta in [th0,th1] (nth+1 cols), z in [0,nz] (nz+1 rows).
 * verts (z,y,x) = (j, rad cos th, rad sin th); uv = (rad*th, j).
 * flip != 0 reverses face winding. mirror != 0 negates u.
 * Appends into caller arrays (offset = existing vert count). */
static void build_arc(float rad, double th0, double th1, int nth, int nz,
                      int flip, int mirror,
                      float **verts, float **uv, int32_t **faces,
                      size_t *nv, size_t *nf,
                      size_t *cap_v, size_t *cap_t, size_t *cap_f)
{
    int i = 0, j = 0;
    size_t base = *nv;
    for (j = 0; j <= nz; j++) {
        for (i = 0; i <= nth; i++) {
            double th = th0 + (th1 - th0) * (double)i / (double)nth;
            double u = (double)rad * th * (mirror ? -1.0 : 1.0);
            *verts = (float *)xgrow(*verts, (*nv + 1) * 3, cap_v, sizeof(float));
            *uv = (float *)xgrow(*uv, (*nv + 1) * 2, cap_t, sizeof(float));
            (*verts)[*nv * 3 + 0] = (float)j;
            (*verts)[*nv * 3 + 1] = rad * (float)cos(th);
            (*verts)[*nv * 3 + 2] = rad * (float)sin(th);
            (*uv)[*nv * 2 + 0] = (float)u;
            (*uv)[*nv * 2 + 1] = (float)j;
            (*nv)++;
        }
    }
    for (j = 0; j < nz; j++) {
        for (i = 0; i < nth; i++) {
            int32_t a = (int32_t)(base + (size_t)j * (size_t)(nth + 1) + (size_t)i);
            int32_t b = a + 1;
            int32_t c = a + nth + 1;
            int32_t d = c + 1;
            int32_t t1[3], t2[3];
            /* (a,b,c): gn ~ d_theta x d_z = radial OUT (see derivation) */
            t1[0] = a; t1[1] = b; t1[2] = c;
            t2[0] = b; t2[1] = d; t2[2] = c;
            if (flip) {
                int32_t t = t1[1]; t1[1] = t1[2]; t1[2] = t;
                t = t2[1]; t2[1] = t2[2]; t2[2] = t;
            }
            *faces = (int32_t *)xgrow(*faces, (*nf + 2) * 3, cap_f,
                                      sizeof(int32_t));
            memcpy(&(*faces)[*nf * 3], t1, sizeof t1); (*nf)++;
            memcpy(&(*faces)[*nf * 3], t2, sizeof t2); (*nf)++;
        }
    }
}

static int selftest(void)
{
    int nfail = 0;
    const float AXP[3] = { 0.0f, 0.0f, 0.0f };
    const float AXD[3] = { 1.0f, 0.0f, 0.0f };
    Arena_T arena = Arena_new();
    int K1 = 0;

    fprintf(stderr, "[selftest] obj_orient_audit\n");

    /* t1: outward-wound arc, standard uv */
    {
        float *v = NULL, *t = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, ct = 0, cf = 0;
        AuditResult R;
        build_arc(50.0f, 0.2, 2.2, 40, 10, 0, 0, &v, &t, &fc, &nv, &nf,
                  &cv, &ct, &cf);
        CHECK(audit_run(arena, v, nv, t, fc, nf, AXP, AXD, 0.3, &R) == 0,
              "t1 audit runs");
        CHECK(R.ncomp == 1, "t1 one component");
        CHECK(comp_rad_sign(&R.comp[0]) == 1 && comp_rad_agr(&R.comp[0]) > 0.99,
              "t1 radial OUT, coherent");
        CHECK(R.comp[0].same_dir == 0, "t1 same_dir == 0");
        CHECK(comp_kap_agr(&R.comp[0]) > 0.99, "t1 kappa coherent");
        K1 = comp_kap_sign(&R.comp[0]);
        fprintf(stderr, "  (t1 kappa reference = %+d)\n", K1);
        free(v); free(t); free(fc);
    }

    /* t2: all faces flipped -> radial IN, kappa UNCHANGED (invariance) */
    {
        float *v = NULL, *t = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, ct = 0, cf = 0;
        AuditResult R;
        build_arc(50.0f, 0.2, 2.2, 40, 10, 1, 0, &v, &t, &fc, &nv, &nf,
                  &cv, &ct, &cf);
        audit_run(arena, v, nv, t, fc, nf, AXP, AXD, 0.3, &R);
        CHECK(comp_rad_sign(&R.comp[0]) == -1 && comp_rad_agr(&R.comp[0]) > 0.99,
              "t2 flipped winding -> radial IN");
        CHECK(comp_kap_sign(&R.comp[0]) == K1 && comp_kap_agr(&R.comp[0]) > 0.99,
              "t2 kappa INVARIANT under winding flip");
        CHECK(R.comp[0].same_dir == 0, "t2 still consistently wound");
        free(v); free(t); free(fc);
    }

    /* t3: mirrored uv (u -> -u), original winding -> kappa flips */
    {
        float *v = NULL, *t = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, ct = 0, cf = 0;
        AuditResult R;
        build_arc(50.0f, 0.2, 2.2, 40, 10, 0, 1, &v, &t, &fc, &nv, &nf,
                  &cv, &ct, &cf);
        audit_run(arena, v, nv, t, fc, nf, AXP, AXD, 0.3, &R);
        CHECK(comp_rad_sign(&R.comp[0]) == 1, "t3 winding untouched");
        CHECK(comp_kap_sign(&R.comp[0]) == -K1 && comp_kap_agr(&R.comp[0]) > 0.99,
              "t3 mirrored uv -> kappa flipped");
        free(v); free(t); free(fc);
    }

    /* t4: two components, second winding-flipped at another azimuth ->
     * WINDING-FLIPPED detected, NOT mirrored */
    {
        float *v = NULL, *t = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, ct = 0, cf = 0;
        AuditResult R;
        build_arc(50.0f, 0.2, 2.2, 40, 10, 0, 0, &v, &t, &fc, &nv, &nf,
                  &cv, &ct, &cf);                       /* main (bigger) */
        build_arc(60.0f, 3.6, 4.6, 12, 6, 1, 0, &v, &t, &fc, &nv, &nf,
                  &cv, &ct, &cf);                       /* far side, flipped */
        audit_run(arena, v, nv, t, fc, nf, AXP, AXD, 0.3, &R);
        CHECK(R.ncomp == 2, "t4 two components");
        CHECK(R.comp[0].area > R.comp[1].area, "t4 main is largest");
        CHECK(R.R0 == 1, "t4 main radial OUT");
        CHECK(comp_rad_sign(&R.comp[1]) == -1
              && comp_rad_agr(&R.comp[1]) > 0.99,
              "t4 second comp WINDING-FLIPPED");
        CHECK(comp_kap_sign(&R.comp[1]) == R.K0, "t4 second comp NOT mirrored");
        free(v); free(t); free(fc);
    }

    /* t5: mixed winding inside one comp -> same_dir > 0, agreement ~50% */
    {
        float *v = NULL, *t = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, ct = 0, cf = 0, f = 0;
        AuditResult R;
        build_arc(50.0f, 0.2, 2.2, 40, 10, 0, 0, &v, &t, &fc, &nv, &nf,
                  &cv, &ct, &cf);
        for (f = 0; f < nf; f += 2) {   /* flip every other face */
            int32_t tmp = fc[f * 3 + 1];
            fc[f * 3 + 1] = fc[f * 3 + 2];
            fc[f * 3 + 2] = tmp;
        }
        audit_run(arena, v, nv, t, fc, nf, AXP, AXD, 0.3, &R);
        CHECK(R.comp[0].same_dir > 0, "t5 same_dir > 0 on mixed winding");
        CHECK(comp_rad_agr(&R.comp[0]) < 0.8, "t5 agreement < 80% (MIXED)");
        free(v); free(t); free(fc);
    }

    /* t6: no uv -> radial-only path; degenerate input -> error */
    {
        float *v = NULL, *t = NULL; int32_t *fc = NULL;
        size_t nv = 0, nf = 0, cv = 0, ct = 0, cf = 0;
        AuditResult R;
        build_arc(50.0f, 0.2, 2.2, 8, 3, 0, 0, &v, &t, &fc, &nv, &nf,
                  &cv, &ct, &cf);
        CHECK(audit_run(arena, v, nv, NULL, fc, nf, AXP, AXD, 0.3, &R) == 0
              && comp_rad_sign(&R.comp[0]) == 1,
              "t6 radial-only audit ok");
        CHECK(audit_run(arena, v, 0, NULL, fc, 0, AXP, AXD, 0.3, &R) != 0,
              "t6 empty input -> error");
        free(v); free(t); free(fc);
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] %s (%d failures)\n",
            nfail == 0 ? "ALL PASS" : "FAILURES", nfail);
    return nfail;
}

/* ------------------------------------------------------------------ main */

static void usage(void)
{
    fprintf(stderr,
        "usage: obj_orient_audit <mesh.obj> [options]\n"
        "       obj_orient_audit --selftest\n"
        "  Per-component stored-winding orientation + UV-mirror audit.\n"
        "  mesh.obj: source-space (z,y,x) verts, optional per-vertex vt.\n"
        "options:\n"
        "  --axis-point Z Y X  scroll axis point (default 0 3405 2878)\n"
        "  --axis-dir Z Y X    scroll axis direction (default 1 0 0)\n"
        "  --min-cos F         radial-clarity gate |cos| (default 0.3)\n"
        "  --out-dir D --id S  write colored flat dumps (needs vt):\n"
        "                      <id>_orient_flat.obj, <id>_mirror_flat.obj\n"
        "  --max-rows N        table rows to print (default 40)\n");
}

int main(int argc, char **argv)
{
    const char *in_path = NULL, *out_dir = NULL, *id = "audit";
    float axp[3] = { 0.0f, 3405.0f, 2878.0f };
    float axd[3] = { 1.0f, 0.0f, 0.0f };
    double min_cos = 0.3;
    size_t max_rows = 40;
    int i = 0;
    UvMesh m;
    AuditResult R;
    Arena_T arena = NULL;
    double t0 = 0.0;

    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return selftest() == 0 ? 0 : 1;
    if (argc < 2) { usage(); return 1; }
    in_path = argv[1];
    for (i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--axis-point") == 0 && i + 3 < argc) {
            axp[0] = (float)strtod(argv[++i], NULL);
            axp[1] = (float)strtod(argv[++i], NULL);
            axp[2] = (float)strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--axis-dir") == 0 && i + 3 < argc) {
            axd[0] = (float)strtod(argv[++i], NULL);
            axd[1] = (float)strtod(argv[++i], NULL);
            axd[2] = (float)strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--min-cos") == 0 && i + 1 < argc) {
            min_cos = strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (strcmp(argv[i], "--id") == 0 && i + 1 < argc) {
            id = argv[++i];
        } else if (strcmp(argv[i], "--max-rows") == 0 && i + 1 < argc) {
            max_rows = (size_t)strtol(argv[++i], NULL, 10);
        } else {
            fprintf(stderr, "unknown option %s\n", argv[i]);
            usage();
            return 1;
        }
    }

    fprintf(stderr, "[obj_orient_audit] in=%s axis=(%.0f,%.0f,%.0f)+"
            "t(%.0f,%.0f,%.0f)\n", in_path,
            (double)axp[0], (double)axp[1], (double)axp[2],
            (double)axd[0], (double)axd[1], (double)axd[2]);
    t0 = ves_clock_sec();
    if (parse_uv_obj(in_path, &m) != 0) return 1;
    if (m.nvt != 0 && m.nvt != m.nv) {
        fprintf(stderr, "  WARN: nvt (%zu) != nv (%zu) -- ignoring vt\n",
                m.nvt, m.nv);
        free(m.uv);
        m.uv = NULL;
        m.nvt = 0;
    }
    fprintf(stderr, "  parsed: %zu verts, %zu vt, %zu faces (%.1fs)\n",
            m.nv, m.nvt, m.nf, ves_clock_sec() - t0);

    arena = Arena_new();
    if (audit_run(arena, m.verts, m.nv, m.uv, m.faces, m.nf, axp, axd,
                  min_cos, &R) != 0) {
        fprintf(stderr, "ERROR: audit failed (degenerate input?)\n");
        return 1;
    }
    print_report(&R, m.uv, max_rows);

    if (out_dir != NULL && m.uv != NULL) {
        char path[2600];
        snprintf(path, sizeof path, "%s/%s_orient_flat.obj", out_dir, id);
        ves_ensure_parent_dir(path);
        if (write_vote_flat(path, m.uv, m.nv, m.faces, m.nf, m.verts, &R, 0) == 0)
            fprintf(stderr, "  wrote %s\n", path);
        snprintf(path, sizeof path, "%s/%s_mirror_flat.obj", out_dir, id);
        if (write_vote_flat(path, m.uv, m.nv, m.faces, m.nf, m.verts, &R, 1) == 0)
            fprintf(stderr, "  wrote %s\n", path);
    }

    uvmesh_free(&m);
    Arena_dispose(&arena);
    return 0;
}
