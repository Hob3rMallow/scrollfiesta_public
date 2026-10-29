/* asm_clean.c -- one cube mesh -> cleaned charts with intrinsic flattenings. */
#include "asm_clean.h"
#include "../common/kdtree.h"
#include "../common/closest_tri.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/csr.h"
#include "../common/mesh_manifold.h"
#include "../common/mesh_normals.h"
#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"
#include "../flatten/mesh_topo.h"
#include "asm_flatten.h"

void AsmClean_default_opts(AsmCleanOpts *o)
{
    o->blob_area_per_face = ASM_BLOB_AREA_PER_FACE;
    o->blob_double_sided = ASM_BLOB_DOUBLE_SIDED;
    o->min_chart_area = ASM_MIN_CHART_AREA;
    o->stress_band = ASM_STRESS_BAND;
    o->flatten_iters = ASM_FLATTEN_ITERS;
}

static double ac_tri_area(const float *v, int32_t a, int32_t b, int32_t c)
{
    double ax = v[(size_t)a*3], ay = v[(size_t)a*3+1], az = v[(size_t)a*3+2];
    double e1[3] = { v[(size_t)b*3]-ax, v[(size_t)b*3+1]-ay, v[(size_t)b*3+2]-az };
    double e2[3] = { v[(size_t)c*3]-ax, v[(size_t)c*3+1]-ay, v[(size_t)c*3+2]-az };
    double cx = e1[1]*e2[2] - e1[2]*e2[1];
    double cy = e1[2]*e2[0] - e1[0]*e2[2];
    double cz = e1[0]*e2[1] - e1[1]*e2[0];
    return 0.5 * sqrt(cx*cx + cy*cy + cz*cz);
}

void AsmClean_boundary_mask(Arena_T scratch, size_t nv,
                            const int32_t *faces, size_t nf, uint8_t *mark)
{
    Arena_Mark m = Arena_save(scratch);
    memset(mark, 0, nv);
    CSR_T adj = CSR_from_faces(scratch, faces, nf, nv);
    const int32_t *off = CSR_offset(adj);
    const int32_t *col = CSR_target(adj);
    size_t nnz = (size_t)CSR_nnz(adj);
    int32_t *cnt = ARENA_CALLOC(scratch, nnz, sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f*3 + (size_t)e], b = faces[f*3 + (size_t)((e+1)%3)];
            for (int32_t k = off[a]; k < off[a+1]; k++) if (col[k] == b) { cnt[k]++; break; }
            for (int32_t k = off[b]; k < off[b+1]; k++) if (col[k] == a) { cnt[k]++; break; }
        }
    }
    for (size_t i = 0; i < nv; i++)
        for (int32_t k = off[i]; k < off[i+1]; k++)
            if (cnt[k] == 1) { mark[i] = 1; break; }
    Arena_restore(scratch, m);
}

/* DOUBLE-SIDED test: the fraction of sampled vertices that find the chart's OWN
 * surface between 1.5 and 8 vox along their normal (either way).  A fused
 * prediction blob is meshed as a closed ball or a thick slab whose far side sits
 * a few voxels behind every point; a papyrus sheet is one face thick and its
 * neighbours are other charts a layer away.  Independent of mesh density. */
static double ac_double_sided(Arena_T scratch, const AsmChart *ch)
{
    size_t nf = ch->nf, nv = ch->nv;
    if (nf < 8 || nv < 8) return 0.0;
    Arena_Mark mark = Arena_save(scratch);
    float *fc = ARENA_ALLOC(scratch, nf * 3 * sizeof(float));
    for (size_t f = 0; f < nf; f++) {
        const int32_t *t = &ch->faces[f*3];
        for (int d = 0; d < 3; d++)
            fc[f*3 + (size_t)d] = (ch->xyz[(size_t)t[0]*3 + (size_t)d] + ch->xyz[(size_t)t[1]*3 + (size_t)d] + ch->xyz[(size_t)t[2]*3 + (size_t)d]) / 3.0f;
    }
    KDTree_T tree = KDTree_new(scratch, fc, nf);
    if (tree == NULL) { Arena_restore(scratch, mark); return 0.0; }
    size_t stride = nv / 1500; if (stride < 1) stride = 1;
    size_t sampled = 0, hit = 0;
    int32_t ball[64];
    for (size_t v = 0; v < nv; v += stride) {
        const float *p = &ch->xyz[v*3], *n = &ch->nrm[v*3];
        double nl = sqrt((double)n[0]*n[0] + (double)n[1]*n[1] + (double)n[2]*n[2]);
        if (nl < 1e-6) continue;
        sampled++;
        int found = 0;
        for (int side = -1; side <= 1 && !found; side += 2) {
            for (double s = 2.0; s <= 7.0 && !found; s += 1.0) {
                float q[3] = { (float)(p[0] + side * s * n[0] / nl), (float)(p[1] + side * s * n[1] / nl), (float)(p[2] + side * s * n[2] / nl) };
                size_t m = KDTree_ball_query(tree, q, 9.0f, ball, 64);
                for (size_t b = 0; b < m; b++) {
                    const int32_t *t = &ch->faces[(size_t)ball[b]*3];
                    double a[3], bb[3], cc[3], qd[3] = { q[0], q[1], q[2] }, cp[3], cu, cv;
                    for (int d = 0; d < 3; d++) { a[d] = ch->xyz[(size_t)t[0]*3 + (size_t)d]; bb[d] = ch->xyz[(size_t)t[1]*3 + (size_t)d]; cc[d] = ch->xyz[(size_t)t[2]*3 + (size_t)d]; }
                    ClosestTri_point(qd, a, bb, cc, cp, &cu, &cv);
                    double dq = sqrt((cp[0]-qd[0])*(cp[0]-qd[0]) + (cp[1]-qd[1])*(cp[1]-qd[1]) + (cp[2]-qd[2])*(cp[2]-qd[2]));
                    if (dq > 1.0) continue;
                    /* the closest point must be a second surface: 1.5..8 vox from p along the normal */
                    double along = ((cp[0]-p[0])*n[0] + (cp[1]-p[1])*n[1] + (cp[2]-p[2])*n[2]) / nl;
                    if (fabs(along) >= 1.5 && fabs(along) <= 8.0) { found = 1; break; }
                }
            }
        }
        hit += (size_t)found;
    }
    Arena_restore(scratch, mark);
    return sampled ? (double)hit / (double)sampled : 0.0;
}

/* bbox and vertex mean of a chart */
static void ac_summaries(AsmChart *ch)
{
    ch->bbox_lo[0] = ch->bbox_lo[1] = ch->bbox_lo[2] = 1e30f;
    ch->bbox_hi[0] = ch->bbox_hi[1] = ch->bbox_hi[2] = -1e30f;
    double sum[3] = { 0.0, 0.0, 0.0 };
    for (size_t l = 0; l < ch->nv; l++)
        for (int d = 0; d < 3; d++) {
            float x = ch->xyz[l*3 + (size_t)d];
            if (x < ch->bbox_lo[d]) ch->bbox_lo[d] = x;
            if (x > ch->bbox_hi[d]) ch->bbox_hi[d] = x;
            sum[d] += x;
        }
    for (int d = 0; d < 3; d++) ch->centroid[d] = ch->nv ? sum[d] / (double)ch->nv : 0.0;
}

/* The faces of `src` selected by `take` (1 = take) as a chart with its own vertex arrays, vertices
 * renumbered in their relative order (the order AsmFlatten_drop_degenerate uses). */
static void ac_subchart(Arena_T scratch, Arena_T persist, const AsmChart *src, const uint8_t *take, AsmChart *out)
{
    Arena_Mark mark = Arena_save(scratch);
    int32_t *vm = ARENA_ALLOC(scratch, src->nv * sizeof(int32_t));
    for (size_t i = 0; i < src->nv; i++) vm[i] = -1;
    size_t nf = 0, nv = 0;
    for (size_t f = 0; f < src->nf; f++) if (take[f]) {
        nf++;
        for (int e = 0; e < 3; e++) vm[src->faces[3*f + (size_t)e]] = 0;
    }
    for (size_t i = 0; i < src->nv; i++) if (vm[i] >= 0) vm[i] = (int32_t)nv++;
    out->nv = nv; out->nf = nf;
    out->faces = ARENA_ALLOC(persist, (nf ? nf : 1) * 3 * sizeof(int32_t));
    out->vid = ARENA_ALLOC(persist, (nv ? nv : 1) * sizeof(int32_t));
    out->xyz = ARENA_ALLOC(persist, (nv ? nv : 1) * 3 * sizeof(float));
    out->nrm = ARENA_ALLOC(persist, (nv ? nv : 1) * 3 * sizeof(float));
    out->boundary = ARENA_CALLOC(persist, nv ? nv : 1, sizeof(uint8_t));
    size_t k = 0;
    for (size_t f = 0; f < src->nf; f++) if (take[f]) {
        for (int e = 0; e < 3; e++) out->faces[3*k + (size_t)e] = vm[src->faces[3*f + (size_t)e]];
        k++;
    }
    for (size_t i = 0; i < src->nv; i++) if (vm[i] >= 0) {
        out->vid[vm[i]] = src->vid[i];
        memcpy(out->xyz + 3*(size_t)vm[i], src->xyz + 3*i, 3 * sizeof(float));
        memcpy(out->nrm + 3*(size_t)vm[i], src->nrm + 3*i, 3 * sizeof(float));
    }
    Arena_restore(scratch, mark);
    if (nv && nf) AsmClean_boundary_mask(scratch, nv, out->faces, nf, out->boundary);
    ac_summaries(out);
}

/* The chart failed to flatten only because of zero-area faces (the mesher closes zero-width slits
 * with pairs of collinear triangles).  Flatten it without them; on success the chart keeps its other
 * faces with that map and `rest` becomes an excluded TINY chart of the zero-area faces, so every
 * source face is still owned by exactly one chart.  On failure nothing changes.  Returns 0 on
 * success, with the map in uv and its statistics in fs. */
static int ac_split_degenerate(Arena_T scratch, Arena_T persist, AsmChart *ch, float *uv,
                               const AsmCleanOpts *opts, AsmChart *rest, AsmFlattenStats *fs, AsmCleanStats *st)
{
    Arena_Mark mark = Arena_save(scratch);
    size_t nv = ch->nv, nf = ch->nf;
    int32_t *kf = ARENA_ALLOC(scratch, nf * 3 * sizeof(int32_t));
    int32_t *vm = ARENA_ALLOC(scratch, nv * sizeof(int32_t));
    size_t knv = 0, knf = AsmFlatten_drop_degenerate(ch->xyz, nv, ch->faces, nf, kf, vm, &knv);
    if (knf == 0 || knf == nf) { Arena_restore(scratch, mark); return -1; }
    float *kx = ARENA_ALLOC(scratch, knv * 3 * sizeof(float));
    for (size_t i = 0; i < nv; i++) if (vm[i] >= 0) memcpy(kx + 3*(size_t)vm[i], ch->xyz + 3*i, 3 * sizeof(float));
    AsmFlattenStats ks;
    if (AsmFlatten_chart_certified(scratch, kx, knv, kf, knf, opts->stress_band, opts->flatten_iters, uv, &ks) != 0) {
        Arena_restore(scratch, mark);
        return -1;
    }
    uint8_t *keep = ARENA_ALLOC(scratch, nf), *drop = ARENA_ALLOC(scratch, nf);
    for (size_t f = 0; f < nf; f++) {
        int degenerate = AsmFlatten_face_degenerate(ch->xyz, ch->faces + 3*f);
        keep[f] = (uint8_t)!degenerate; drop[f] = (uint8_t)degenerate;
    }
    memset(rest, 0, sizeof *rest);
    rest->id = -1; rest->cube = ch->cube; rest->comp = ch->comp; rest->parent = -1;
    rest->component = -1; rest->placed = 0; rest->flags = ASM_CHART_TINY;
    ac_subchart(scratch, persist, ch, drop, rest);
    for (size_t f = 0; f < rest->nf; f++)
        rest->area3d += ac_tri_area(rest->xyz, rest->faces[3*f], rest->faces[3*f+1], rest->faces[3*f+2]);
    AsmChart main = *ch;
    ac_subchart(scratch, persist, ch, keep, &main);
    /* the flattened arrays and the kept chart number their vertices the same way */
    if (main.nv != knv || main.nf != knf || memcmp(main.faces, kf, knf * 3 * sizeof(int32_t))) {
        Arena_restore(scratch, mark);
        return -1;
    }
    *ch = main;
    *fs = ks;
    st->degenerate_split++;
    st->degenerate_faces += nf - knf;
    Arena_restore(scratch, mark);
    return 0;
}

/* Charts split off a component, appended after the components.  The list lives in `persist`: a
 * scratch mark restored inside a split would otherwise free it. */
typedef struct AcExtras { AsmChart *chart; size_t n, cap; } AcExtras;

static void ac_push_extra(Arena_T persist, AcExtras *x, const AsmChart *c)
{
    if (x->n == x->cap) {
        size_t cap = x->cap ? 2 * x->cap : 8;
        AsmChart *grown = ARENA_ALLOC(persist, cap * sizeof(AsmChart));
        if (x->n) memcpy(grown, x->chart, x->n * sizeof(AsmChart));
        x->chart = grown; x->cap = cap;
    }
    x->chart[x->n++] = *c;
}

/* area, singular-value range and stress of a chart's map, as AsmClean_cube records them */
static void ac_map_stats(Arena_T scratch, AsmChart *ch, double band)
{
    Arena_Mark mark = Arena_save(scratch);
    double *lo = ARENA_ALLOC(scratch, ch->nf * sizeof(double)), *hi = ARENA_ALLOC(scratch, ch->nf * sizeof(double));
    int8_t *sign = ARENA_ALLOC(scratch, ch->nf);
    AsmFlatten_face_singular(ch->xyz, ch->faces, ch->nf, ch->uv, lo, hi, sign);
    ch->area3d = 0.0; ch->area_uv = 0.0; ch->sigma_lo = 1e300; ch->sigma_hi = 0.0; ch->n_flipped = 0;
    size_t stressed = 0;
    for (size_t f = 0; f < ch->nf; f++) {
        double a = ac_tri_area(ch->xyz, ch->faces[3*f], ch->faces[3*f+1], ch->faces[3*f+2]);
        ch->area3d += a;
        ch->area_uv += lo[f] * hi[f] * a;
        if (lo[f] < ch->sigma_lo) ch->sigma_lo = lo[f];
        if (hi[f] > ch->sigma_hi) ch->sigma_hi = hi[f];
        if (sign[f] < 0) ch->n_flipped++;
        if (lo[f] < 1.0 - band || hi[f] > 1.0 + band || sign[f] <= 0) stressed++;
    }
    ch->stress_frac = ch->nf ? (double)stressed / (double)ch->nf : 0.0;
    Arena_restore(scratch, mark);
}

/* No seed flattens the chart: excise its crumpled regions (AsmFlatten_rescue_crumpled).  On success
 * the chart becomes its largest certified piece (map in uv, stats in fs), every other certified
 * piece is appended as a chart in the layout, and all remaining faces (excised, zero-area,
 * uncertified or small) are appended as one excluded FLAT_FAILED chart, so every source face is
 * still owned by exactly one chart.  On failure nothing changes.  Returns 0 on success. */
static int ac_split_crumpled(Arena_T scratch, Arena_T persist, AsmChart *ch, float *uv,
                             const AsmCleanOpts *opts, AcExtras *extras, AsmFlattenStats *fs, AsmCleanStats *st)
{
    Arena_Mark mark = Arena_save(scratch);
    size_t nv = ch->nv, nf = ch->nf;
    /* zero-area faces take no part: they stay with the remainder */
    int32_t *kf = ARENA_ALLOC(scratch, nf * 3 * sizeof(int32_t)), *kidx = ARENA_ALLOC(scratch, nf * sizeof(int32_t));
    size_t knf = 0;
    for (size_t f = 0; f < nf; f++) if (!AsmFlatten_face_degenerate(ch->xyz, ch->faces + 3*f)) {
        memcpy(kf + 3*knf, ch->faces + 3*f, 3 * sizeof(int32_t));
        kidx[knf++] = (int32_t)f;
    }
    int32_t *piece = ARENA_ALLOC(scratch, (knf ? knf : 1) * sizeof(int32_t));
    float *puv = ARENA_ALLOC(scratch, nv * 2 * sizeof(float));
    size_t np = 0;
    AsmFlattenStats ps;
    if (!knf || AsmFlatten_rescue_crumpled(scratch, ch->xyz, nv, kf, knf, opts->stress_band, opts->flatten_iters,
                                           opts->min_chart_area, piece, puv, &np, &ps) != 0 || !np) {
        Arena_restore(scratch, mark);
        return -1;
    }
    int32_t *label = ARENA_ALLOC(scratch, nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) label[f] = -1;
    for (size_t k = 0; k < knf; k++) label[kidx[k]] = piece[k];
    uint8_t *take = ARENA_ALLOC(scratch, nf);
    int32_t *vm = ARENA_ALLOC(scratch, nv * sizeof(int32_t));
    AsmChart source = *ch, main_piece;
    memset(&main_piece, 0, sizeof main_piece);
    for (size_t p = 0; p <= np; p++) {
        int32_t want = p < np ? (int32_t)p : -1;             /* the pieces, then the remainder */
        size_t count = 0;
        for (size_t f = 0; f < nf; f++) { take[f] = (uint8_t)(label[f] == want); count += take[f]; }
        if (!count) continue;
        AsmChart c = source;
        c.flags = 0; c.uv = NULL; c.placed_uv = NULL; c.repaired = NULL;
        c.area_uv = 0.0; c.sigma_lo = c.sigma_hi = c.stress_frac = 0.0; c.n_flipped = 0;
        ac_subchart(scratch, persist, &source, take, &c);
        if (want < 0) {
            c.flags = ASM_CHART_FLAT_FAILED;
            c.area3d = 0.0;
            for (size_t f = 0; f < c.nf; f++) c.area3d += ac_tri_area(c.xyz, c.faces[3*f], c.faces[3*f+1], c.faces[3*f+2]);
            st->crumple_faces += c.nf;
            st->crumple_area_excluded += c.area3d;
            ac_push_extra(persist, extras, &c);
            continue;
        }
        /* the piece's vertices in the order ac_subchart numbers them */
        for (size_t v = 0; v < nv; v++) vm[v] = -1;
        for (size_t f = 0; f < nf; f++) if (take[f]) for (int e = 0; e < 3; e++) vm[source.faces[3*f + (size_t)e]] = 0;
        size_t n = 0;
        for (size_t v = 0; v < nv; v++) if (vm[v] >= 0) vm[v] = (int32_t)n++;
        float *cuv = p == 0 ? uv : ARENA_ALLOC(persist, c.nv * 2 * sizeof(float));
        for (size_t v = 0; v < nv; v++) if (vm[v] >= 0) memcpy(cuv + 2*(size_t)vm[v], puv + 2*v, 2 * sizeof(float));
        c.uv = cuv;
        ac_map_stats(scratch, &c, opts->stress_band);
        if (c.stress_frac > 0.0 || c.n_flipped > 0) { c.flags |= ASM_CHART_SUSPECT; if (p) st->suspect++; }
        if (p == 0) { c.uv = NULL; main_piece = c; continue; }   /* the caller sets its uv and counts it */
        st->kept++; st->verts_kept += c.nv; st->faces_kept += c.nf; st->area_kept += c.area3d;
        ac_push_extra(persist, extras, &c);
    }
    *ch = main_piece;
    *fs = ps;
    st->crumple_split++;
    st->crumple_pieces += np;
    Arena_restore(scratch, mark);
    return 0;
}

int AsmClean_cube(Arena_T scratch, Arena_T persist,
                  const MeshBinData *mesh, int32_t cube,
                  const AsmCleanOpts *opts,
                  AsmChart **out_charts, size_t *out_n,
                  AsmCleanStats *st)
{
    memset(st, 0, sizeof *st);
    *out_charts = NULL; *out_n = 0;
    size_t nv = mesh->nv, nf = mesh->nf;
    st->verts_in = nv; st->faces_in = nf;
    if (nv == 0 || nf == 0) return 0;
    Arena_Mark mark = Arena_save(scratch);

    /* components over vertices through faces */
    UnionFind uf = UF_new(scratch, (int32_t)nv);
    for (size_t f = 0; f < nf; f++) {
        uf_union(&uf, mesh->faces[f*3], mesh->faces[f*3+1]);
        uf_union(&uf, mesh->faces[f*3], mesh->faces[f*3+2]);
    }
    int32_t *root = ARENA_ALLOC(scratch, nv * sizeof(int32_t));
    int32_t *comp_of_root = ARENA_ALLOC(scratch, nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) { root[i] = uf_find(&uf, (int32_t)i); comp_of_root[i] = -1; }
    /* order components by first face occurrence so ids are deterministic */
    size_t ncomp = 0;
    int32_t *comp_of_face = ARENA_ALLOC(scratch, nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) {
        int32_t r = root[mesh->faces[f*3]];
        if (comp_of_root[r] < 0) comp_of_root[r] = (int32_t)ncomp++;
        comp_of_face[f] = comp_of_root[r];
    }
    st->components = ncomp;
    size_t *c_nf = ARENA_CALLOC(scratch, ncomp, sizeof(size_t));
    double *c_area = ARENA_CALLOC(scratch, ncomp, sizeof(double));
    for (size_t f = 0; f < nf; f++) {
        double a = ac_tri_area(mesh->verts, mesh->faces[f*3], mesh->faces[f*3+1], mesh->faces[f*3+2]);
        c_nf[(size_t)comp_of_face[f]]++;
        c_area[(size_t)comp_of_face[f]] += a;
        st->area_in += a;
    }
    /* per-vertex normals for the whole cube (winding-consistent per chart) */
    float *nrm_all = MeshNormals_compute(mesh->verts, nv, mesh->faces, nf);

    AsmChart *charts = ARENA_CALLOC(persist, ncomp, sizeof(AsmChart));
    AcExtras extras = { NULL, 0, 0 };                 /* charts split off a component, after the components */
    int32_t *local = ARENA_ALLOC(scratch, nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) local[i] = -1;
    /* face lists per component: counting sort */
    size_t *c_off = ARENA_CALLOC(scratch, ncomp + 1, sizeof(size_t));
    for (size_t c = 0; c < ncomp; c++) c_off[c+1] = c_off[c] + c_nf[c];
    size_t *c_fill = ARENA_CALLOC(scratch, ncomp, sizeof(size_t));
    int32_t *face_order = ARENA_ALLOC(scratch, nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) {
        size_t c = (size_t)comp_of_face[f];
        face_order[c_off[c] + c_fill[c]++] = (int32_t)f;
    }

    double t_flat = 0.0;
    for (size_t c = 0; c < ncomp; c++) {
        AsmChart *ch = &charts[c];
        ch->id = -1; ch->cube = cube; ch->comp = (int32_t)c; ch->parent = -1;
        ch->component = -1; ch->placed = 0;
        ch->nf = c_nf[c];
        ch->area3d = c_area[c];
        ch->faces = ARENA_ALLOC(persist, ch->nf * 3 * sizeof(int32_t));
        /* local vertex numbering in order of first appearance */
        size_t lnv = 0;
        for (size_t k = 0; k < ch->nf; k++) {
            int32_t f = face_order[c_off[c] + k];
            for (int e = 0; e < 3; e++) {
                int32_t g = mesh->faces[(size_t)f*3 + (size_t)e];
                if (local[g] < 0) local[g] = (int32_t)lnv++;
                ch->faces[k*3 + (size_t)e] = local[g];
            }
        }
        ch->nv = lnv;
        ch->vid = ARENA_ALLOC(persist, lnv * sizeof(int32_t));
        ch->xyz = ARENA_ALLOC(persist, lnv * 3 * sizeof(float));
        ch->nrm = ARENA_ALLOC(persist, lnv * 3 * sizeof(float));
        ch->boundary = ARENA_CALLOC(persist, lnv, sizeof(uint8_t));
        ch->repaired = NULL;
        ch->uv = NULL;
        for (size_t k = 0; k < ch->nf; k++) {
            int32_t f = face_order[c_off[c] + k];
            for (int e = 0; e < 3; e++) {
                int32_t g = mesh->faces[(size_t)f*3 + (size_t)e];
                int32_t l = local[g];
                ch->vid[l] = g;
                memcpy(&ch->xyz[(size_t)l*3], &mesh->verts[(size_t)g*3], 3 * sizeof(float));
                memcpy(&ch->nrm[(size_t)l*3], &nrm_all[(size_t)g*3], 3 * sizeof(float));
            }
        }
        /* reset local map for the next component */
        for (size_t l = 0; l < lnv; l++) local[ch->vid[l]] = -1;

        /* geometry summaries */
        ch->bbox_lo[0] = ch->bbox_lo[1] = ch->bbox_lo[2] = 1e30f;
        ch->bbox_hi[0] = ch->bbox_hi[1] = ch->bbox_hi[2] = -1e30f;
        double cz = 0.0, cy = 0.0, cx = 0.0;
        for (size_t l = 0; l < lnv; l++) {
            for (int d = 0; d < 3; d++) {
                float x = ch->xyz[l*3 + (size_t)d];
                if (x < ch->bbox_lo[d]) ch->bbox_lo[d] = x;
                if (x > ch->bbox_hi[d]) ch->bbox_hi[d] = x;
            }
            cz += ch->xyz[l*3]; cy += ch->xyz[l*3+1]; cx += ch->xyz[l*3+2];
        }
        ch->centroid[0] = cz / (double)lnv;
        ch->centroid[1] = cy / (double)lnv;
        ch->centroid[2] = cx / (double)lnv;
        AsmClean_boundary_mask(scratch, lnv, ch->faces, ch->nf, ch->boundary);

        /* classification */
        /* Triangulation density is not evidence of a prediction blob.
         * Measure the original surface even when faces are very small. */
        double dsided = ac_double_sided(scratch, ch);
        if (dsided >= opts->blob_double_sided) {
            ch->flags |= ASM_CHART_BLOB;
            st->blobs++; st->verts_blob += lnv; st->faces_blob += ch->nf; st->area_blob += ch->area3d;
            continue;
        }
        if (ch->area3d < opts->min_chart_area) {
            ch->flags |= ASM_CHART_TINY;
            st->tiny++;
            continue;
        }
        {
            Arena_Mark m2 = Arena_save(scratch);
            MeshManifoldStats ms = MeshManifold_audit(scratch, lnv, ch->faces, ch->nf);
            Arena_restore(scratch, m2);
            if (!MeshManifold_ok(&ms)) {
                ch->flags |= ASM_CHART_NONMANIFOLD;
                st->nonmanifold++;
                continue;
            }
            MeshTopoInfo ti;
            m2 = Arena_save(scratch);
            int trc = MeshTopo_analyze(scratch, ch->xyz, lnv, ch->faces, ch->nf, &ti);
            Arena_restore(scratch, m2);
            if (trc == 0 && ti.genus > 0.5) {
                ch->flags |= ASM_CHART_HANDLE | ASM_CHART_EXTRAS;
                st->handles++;
                continue;
            }
        }
        /* flatten */
        double t0 = ves_clock_sec();
        float *uv = ARENA_ALLOC(persist, lnv * 2 * sizeof(float));
        AsmFlattenStats fs;
        int frc = AsmFlatten_chart(scratch, ch->xyz, lnv, ch->faces, ch->nf,
                                   opts->stress_band, opts->flatten_iters, uv, &fs);
        if (frc != 0 && fs.fail_reason == ASM_FLAT_DEGENERATE) {
            AsmChart rest;
            if (ac_split_degenerate(scratch, persist, ch, uv, opts, &rest, &fs, st) == 0) {
                ac_push_extra(persist, &extras, &rest);
                lnv = ch->nv;
                frc = 0;
            }
        }
        if (frc != 0 && ac_split_crumpled(scratch, persist, ch, uv, opts, &extras, &fs, st) == 0) {
            lnv = ch->nv;
            frc = 0;
        }
        t_flat += ves_clock_sec() - t0;
        if (frc != 0) {
            fprintf(stderr,"[assemble clean] cube %d chart %zu faces %zu flatten failed: reason %d, iterations %d, flipped %zu\n",
                cube,c,ch->nf,fs.fail_reason,fs.arap_iters,fs.n_flipped);
            ch->flags |= ASM_CHART_FLAT_FAILED;
            st->flat_failed++;
            continue;
        }
        ch->uv = uv;
        ch->area_uv = fs.area_uv;
        ch->sigma_lo = fs.sigma_lo;
        ch->sigma_hi = fs.sigma_hi;
        ch->stress_frac = fs.stress_frac;
        ch->n_flipped = fs.n_flipped;
        if (fs.stress_frac > 0.0 || fs.n_flipped > 0) { ch->flags |= ASM_CHART_SUSPECT; st->suspect++; }
        st->kept++;
        st->verts_kept += lnv; st->faces_kept += ch->nf; st->area_kept += ch->area3d;
    }
    st->flatten_sec = t_flat;
    free(nrm_all);
    Arena_restore(scratch, mark);
    if (extras.n) {
        AsmChart *all = ARENA_ALLOC(persist, (ncomp + extras.n) * sizeof(AsmChart));
        memcpy(all, charts, ncomp * sizeof(AsmChart));
        memcpy(all + ncomp, extras.chart, extras.n * sizeof(AsmChart));
        charts = all;
    }
    *out_charts = charts;
    *out_n = ncomp + extras.n;
    return 0;
}

/* ---- selftest ------------------------------------------------------------- */

/* Refining the same physical sheet must not turn it into a prediction blob.
 * A closed four-voxel slab remains a geometric blob at both densities. */
static int ac_density_controls(void)
{
    int fails = 0;
    const int sizes[] = {8, 32};
    for (int slab = 0; slab < 2; slab++) for (int resolution = 0; resolution < 2; resolution++) {
        int n = sizes[resolution];
        size_t plane_nv = (size_t)(n+1)*(size_t)(n+1);
        size_t nv = plane_nv*(size_t)(slab ? 2 : 1);
        size_t nf = (size_t)n*(size_t)n*(size_t)(slab ? 4 : 2)+(slab ? (size_t)8*n : 0);
        Arena_T scratch = Arena_new(), persist = Arena_new();
        float *v = ARENA_ALLOC(scratch, nv*3*sizeof(float));
        int32_t *f = ARENA_ALLOC(scratch, nf*3*sizeof(int32_t));
        for (int side = 0; side <= slab; side++)
            for (int j = 0; j <= n; j++) for (int i = 0; i <= n; i++) {
                size_t k = (size_t)side*plane_nv+(size_t)j*(size_t)(n+1)+(size_t)i;
                v[k*3] = 4.0f*(float)side;
                v[k*3+1] = 16.0f*(float)j/(float)n;
                v[k*3+2] = 16.0f*(float)i/(float)n;
            }
        size_t k = 0;
        for (int side = 0; side <= slab; side++)
            for (int j = 0; j < n; j++) for (int i = 0; i < n; i++) {
                int32_t a = (int32_t)((size_t)side*plane_nv)+j*(n+1)+i;
                int32_t b = a+1, c = a+n+1, d = c+1;
                f[k*3] = a; f[k*3+1] = side ? d : b; f[k*3+2] = side ? b : d; k++;
                f[k*3] = a; f[k*3+1] = side ? c : d; f[k*3+2] = side ? d : c; k++;
            }
        if (slab) for (int wall = 0; wall < 4; wall++) for (int i = 0; i < n; i++) {
            int32_t a, b;
            if (wall == 0) { a = i; b = i+1; }
            else if (wall == 1) { a = i*(n+1)+n; b = a+n+1; }
            else if (wall == 2) { a = n*(n+1)+n-i; b = a-1; }
            else { a = (n-i)*(n+1); b = a-n-1; }
            int32_t aa = a+(int32_t)plane_nv, bb = b+(int32_t)plane_nv;
            f[k*3] = b; f[k*3+1] = a; f[k*3+2] = aa; k++;
            f[k*3] = b; f[k*3+1] = aa; f[k*3+2] = bb; k++;
        }
        MeshBinData mesh = {v, NULL, f, nv, nf};
        AsmCleanOpts options; AsmClean_default_opts(&options);
        AsmChart *charts = NULL; size_t count = 0; AsmCleanStats stats;
        int rc = AsmClean_cube(scratch, persist, &mesh, 0, &options, &charts, &count, &stats);
        int ok = rc == 0 && k == nf && count == 1 && stats.blobs == (size_t)slab
            && stats.kept == (size_t)!slab && stats.tiny == 0;
        if (count == 1) {
            AsmChart *c = charts;
            ok = ok && c->nv == nv && c->nf == nf && fabs(c->area3d-(slab ? 768.0 : 256.0)) < 1e-8;
            for (size_t i = 0; i < c->nv; i++) for (int d = 0; d < 3; d++)
                ok = ok && c->xyz[i*3+(size_t)d] == v[(size_t)c->vid[i]*3+(size_t)d];
            for (size_t i = 0; i < c->nf*3; i++) ok = ok && c->vid[c->faces[i]] == f[i];
            if (!slab) ok = ok && c->uv != NULL && c->sigma_lo > .99 && c->sigma_hi < 1.01;
        }
        fprintf(stderr, "  asm_clean density control: %s %dx%d, area %.9g, faces %zu, kept %zu, blobs %zu: %s\n",
                slab ? "closed slab" : "planar sheet", n, n, stats.area_in, nf, stats.kept, stats.blobs,
                ok ? "passed" : "FAILED");
        fails += !ok;
        Arena_dispose(&scratch); Arena_dispose(&persist);
    }
    return fails;
}

/* A sheet whose only defect is the mesher's closed slit (two collinear triangles over a split
 * vertex) is flattened without them; the zero-area faces become an excluded remainder chart, and
 * every source face is owned by exactly one chart. */
static int ac_zero_area_controls(void)
{
    int fails = 0;
    Arena_T scratch = Arena_new(), persist = Arena_new();
    const int nx = 16, ny = 10;
    size_t row = (size_t)nx + 1, gnv = row * (size_t)(ny + 1), gnf = (size_t)nx * (size_t)ny * 2;
    size_t nv = gnv + 1, nf = gnf + 2;
    float *verts = ARENA_ALLOC(scratch, nv * 3 * sizeof(float));
    int32_t *faces = ARENA_ALLOC(scratch, nf * 3 * sizeof(int32_t));
    for (int j = 0; j <= ny; j++) for (int i = 0; i <= nx; i++) {
        size_t k = (size_t)j * row + (size_t)i;
        verts[k*3] = 0.0f; verts[k*3+1] = 2.0f * (float)j; verts[k*3+2] = 2.0f * (float)i;
    }
    /* the split vertex (8,5): right-hand faces of rows 4..5 take its copy */
    int32_t a = (int32_t)(4 * row + 8), m = (int32_t)(5 * row + 8), b = (int32_t)(6 * row + 8), mc = (int32_t)gnv;
    memcpy(verts + 3*(size_t)mc, verts + 3*(size_t)m, 3 * sizeof(float));
    size_t f = 0;
    for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) {
        int32_t p = (int32_t)((size_t)j * row + (size_t)i), q = p + 1, r = p + (int32_t)row, s = r + 1;
        int32_t tri[6] = { p, q, s, p, s, r };
        if (i == 8 && (j == 4 || j == 5)) for (int k = 0; k < 6; k++) if (tri[k] == m) tri[k] = mc;
        memcpy(faces + 3*f, tri, sizeof tri);
        f += 2;
    }
    /* the zero-area closure, wound against the slit's half-edges: the left faces hold m -> a */
    int32_t fill[6] = { a, m, b, a, b, mc };
    for (size_t g = 0; g < gnf; g++) for (int e = 0; e < 3; e++)
        if (faces[3*g + (size_t)e] == a && faces[3*g + (size_t)((e + 1) % 3)] == m) {
            fill[1] = b; fill[2] = m; fill[4] = mc; fill[5] = b;
        }
    memcpy(faces + 3*gnf, fill, sizeof fill);
    MeshBinData mesh = { verts, NULL, faces, nv, nf, NULL };
    AsmCleanOpts o; AsmClean_default_opts(&o);
    o.min_chart_area = 10.0;
    AsmChart *charts = NULL; size_t n = 0; AsmCleanStats st;
    AsmClean_cube(scratch, persist, &mesh, 3, &o, &charts, &n, &st);
    int ok = n == 2 && st.components == 1 && st.kept == 1 && st.flat_failed == 0 && st.degenerate_split == 1 && st.degenerate_faces == 2;
    if (!ok) { fprintf(stderr, "  asm_clean selftest FAIL: zero-area split: charts %zu kept %zu failed %zu split %zu faces %zu\n",
                       n, st.kept, st.flat_failed, st.degenerate_split, st.degenerate_faces); fails++; }
    if (ok) {
        const AsmChart *main_chart = &charts[0], *rest = &charts[1];
        if (!AsmChart_in_layout(main_chart) || main_chart->nf != gnf || main_chart->nv != nv || !(main_chart->sigma_lo > 0.75 && main_chart->sigma_hi < 1.25)) {
            fprintf(stderr, "  asm_clean selftest FAIL: zero-area split main chart\n"); fails++;
        }
        if (AsmChart_in_layout(rest) || !(rest->flags & ASM_CHART_TINY) || rest->nf != 2 || rest->cube != 3 || rest->comp != main_chart->comp) {
            fprintf(stderr, "  asm_clean selftest FAIL: zero-area remainder chart\n"); fails++;
        }
        /* ownership: every source face exactly once, by its vertex triple */
        uint8_t *owned = ARENA_CALLOC(scratch, nf, 1);
        size_t bad = 0;
        for (size_t c = 0; c < n; c++) for (size_t t = 0; t < charts[c].nf; t++) {
            int32_t v[3];
            for (int e = 0; e < 3; e++) v[e] = charts[c].vid[charts[c].faces[3*t + (size_t)e]];
            size_t hit = nf;
            for (size_t g = 0; g < nf && hit == nf; g++)
                if (faces[3*g] == v[0] && faces[3*g+1] == v[1] && faces[3*g+2] == v[2]) hit = g;
            if (hit == nf || owned[hit]++) bad++;
            for (int e = 0; e < 3; e++) {
                int32_t lv = charts[c].faces[3*t + (size_t)e];
                if (memcmp(charts[c].xyz + 3*(size_t)lv, verts + 3*(size_t)v[e], 3 * sizeof(float))) bad++;
            }
        }
        for (size_t g = 0; g < nf; g++) if (!owned[g]) bad++;
        if (bad) { fprintf(stderr, "  asm_clean selftest FAIL: zero-area split ownership (%zu)\n", bad); fails++; }
    }
    Arena_dispose(&scratch);
    Arena_dispose(&persist);
    return fails;
}

/* The crumpled split on a cone disk (apex defect 1.5 rad): certified pieces in the layout, the
 * excised faces in one excluded FLAT_FAILED chart, every source face owned exactly once. */
static int ac_crumple_controls(void)
{
    int fails = 0;
    Arena_T scratch = Arena_new(), persist = Arena_new();
    const int rings = 12, sectors = 48;
    size_t ns = (size_t)sectors, nv = 1 + (size_t)rings * ns, nf = ns + 2 * (size_t)(rings - 1) * ns;
    float *verts = ARENA_ALLOC(scratch, nv * 3 * sizeof(float));
    int32_t *faces = ARENA_ALLOC(scratch, nf * 3 * sizeof(int32_t));
    double k = (2.0 * 3.14159265358979323846 - 1.5) / (2.0 * 3.14159265358979323846);
    verts[0] = verts[1] = verts[2] = 0.0f;
    for (int i = 1; i <= rings; i++) for (int j = 0; j < sectors; j++) {
        double phi = 2.0 * 3.14159265358979323846 * j / sectors;
        float *p = verts + 3 * (1 + (size_t)(i - 1) * ns + (size_t)j);
        p[0] = (float)(2.0 * i * sqrt(1.0 - k * k)); p[1] = (float)(2.0 * k * i * sin(phi)); p[2] = (float)(2.0 * k * i * cos(phi));
    }
    size_t f = 0;
    for (int j = 0; j < sectors; j++) { faces[3*f] = 0; faces[3*f+1] = 1 + j; faces[3*f+2] = 1 + (j + 1) % sectors; f++; }
    for (int i = 1; i < rings; i++) for (int j = 0; j < sectors; j++) {
        int32_t p = (int32_t)(1 + (size_t)(i - 1) * ns + (size_t)j), q = (int32_t)(1 + (size_t)(i - 1) * ns + (size_t)((j + 1) % sectors));
        int32_t tri[6] = { p, p + sectors, q + sectors, p, q + sectors, q };
        memcpy(faces + 3*f, tri, sizeof tri);
        f += 2;
    }
    MeshBinData mesh = { verts, NULL, faces, nv, nf, NULL };
    AsmCleanOpts o; AsmClean_default_opts(&o);
    o.min_chart_area = 10.0;
    AsmChart *charts = NULL; size_t n = 0; AsmCleanStats st;
    AsmClean_cube(scratch, persist, &mesh, 5, &o, &charts, &n, &st);
    if (n < 1 || !charts[0].vid) { fprintf(stderr, "  asm_clean selftest FAIL: cone disk chart\n"); fails++; }
    else {
        AsmChart ch = charts[0];
        AcExtras extras = { NULL, 0, 0 };
        AsmFlattenStats fs;
        AsmCleanStats cs; memset(&cs, 0, sizeof cs);
        float *uv = ARENA_ALLOC(persist, ch.nv * 2 * sizeof(float));
        int rc = ac_split_crumpled(scratch, persist, &ch, uv, &o, &extras, &fs, &cs);
        ch.uv = uv;
        size_t remainder = 0, pieces = 1;
        for (size_t e = 0; e < extras.n; e++) { if (extras.chart[e].flags & ASM_CHART_FLAT_FAILED) remainder++; else pieces++; }
        if (rc != 0 || remainder != 1 || cs.crumple_split != 1 || cs.crumple_pieces != pieces || !(ch.sigma_lo >= 0.75 - 1e-9 && ch.sigma_hi <= 1.25 + 1e-9)) {
            fprintf(stderr, "  asm_clean selftest FAIL: crumpled split rc %d remainder %zu pieces %zu sigma [%g, %g]\n",
                    rc, remainder, pieces, ch.sigma_lo, ch.sigma_hi);
            fails++;
        }
        uint8_t *owned = ARENA_CALLOC(scratch, nf, 1);
        size_t bad = 0;
        for (size_t c = 0; c <= extras.n; c++) {
            const AsmChart *x = c ? &extras.chart[c - 1] : &ch;
            for (size_t t = 0; t < x->nf; t++) {
                int32_t v[3];
                for (int e = 0; e < 3; e++) v[e] = x->vid[x->faces[3*t + (size_t)e]];
                size_t hit = nf;
                for (size_t g = 0; g < nf && hit == nf; g++)
                    if (faces[3*g] == v[0] && faces[3*g+1] == v[1] && faces[3*g+2] == v[2]) hit = g;
                if (hit == nf || owned[hit]++) bad++;
            }
        }
        for (size_t g = 0; g < nf; g++) if (!owned[g]) bad++;
        if (bad) { fprintf(stderr, "  asm_clean selftest FAIL: crumpled split ownership (%zu)\n", bad); fails++; }
    }
    Arena_dispose(&scratch);
    Arena_dispose(&persist);
    return fails;
}

int AsmClean_selftest(void)
{
    int fails = ac_zero_area_controls() + ac_crumple_controls();
    Arena_T scratch = Arena_new();
    Arena_T persist = Arena_new();
    /* Two planar components at different scales. The small component is
     * a tiny sheet fragment; its fine triangulation does not make it a blob. */
    int nx = 6, ny = 6;
    size_t gnv = (size_t)(nx+1)*(size_t)(ny+1), gnf = (size_t)nx*(size_t)ny*2;
    size_t bnv = 16, bnf = 9*2;
    size_t nv = gnv + bnv, nf = gnf + bnf;
    float *verts = ARENA_ALLOC(scratch, nv * 3 * sizeof(float));
    int32_t *faces = ARENA_ALLOC(scratch, nf * 3 * sizeof(int32_t));
    for (int j = 0; j <= ny; j++) for (int i = 0; i <= nx; i++) {
        size_t k = (size_t)j*(size_t)(nx+1)+(size_t)i;
        verts[k*3] = 0.0f; verts[k*3+1] = 7.0f*(float)j; verts[k*3+2] = 7.0f*(float)i;
    }
    size_t f = 0;
    for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) {
        int32_t a = (int32_t)(j*(nx+1)+i), b = a+1, c = a+nx+1, d = c+1;
        faces[f*3] = a; faces[f*3+1] = b; faces[f*3+2] = d; f++;
        faces[f*3] = a; faces[f*3+1] = d; faces[f*3+2] = c; f++;
    }
    for (int j = 0; j < 4; j++) for (int i = 0; i < 4; i++) {
        size_t k = gnv + (size_t)(j*4+i);
        verts[k*3] = 50.0f; verts[k*3+1] = 0.3f*(float)j; verts[k*3+2] = 0.3f*(float)i;
    }
    for (int j = 0; j < 3; j++) for (int i = 0; i < 3; i++) {
        int32_t a = (int32_t)(gnv + (size_t)(j*4+i)), b = a+1, c = a+4, d = c+1;
        faces[f*3] = a; faces[f*3+1] = b; faces[f*3+2] = d; f++;
        faces[f*3] = a; faces[f*3+1] = d; faces[f*3+2] = c; f++;
    }
    MeshBinData mesh = { verts, NULL, faces, nv, nf };
    AsmCleanOpts o; AsmClean_default_opts(&o);
    o.min_chart_area = 10.0;
    AsmChart *charts = NULL; size_t n = 0; AsmCleanStats st;
    AsmClean_cube(scratch, persist, &mesh, 7, &o, &charts, &n, &st);
    if (n != 2) { fprintf(stderr, "  asm_clean selftest FAIL: components %zu\n", n); fails++; }
    if (st.blobs != 0 || st.tiny != 1 || st.kept != 1) { fprintf(stderr, "  asm_clean selftest FAIL: blobs %zu tiny %zu kept %zu\n", st.blobs, st.tiny, st.kept); fails++; }
    if (n == 2) {
        const AsmChart *g = &charts[0];
        if (g->cube != 7 || g->nv != gnv || g->nf != gnf || g->uv == NULL) { fprintf(stderr, "  asm_clean selftest FAIL: grid chart\n"); fails++; }
        size_t nb = 0; for (size_t i = 0; i < g->nv; i++) nb += g->boundary[i];
        if (nb != (size_t)(2*nx + 2*ny)) { fprintf(stderr, "  asm_clean selftest FAIL: boundary count %zu\n", nb); fails++; }
        if (!(g->sigma_lo > 0.99 && g->sigma_hi < 1.01)) { fprintf(stderr, "  asm_clean selftest FAIL: grid sigma\n"); fails++; }
        if (fabs(g->nrm[0]) < 0.99f) { fprintf(stderr, "  asm_clean selftest FAIL: normal\n"); fails++; }
        if (!(charts[1].flags & ASM_CHART_TINY) || (charts[1].flags & ASM_CHART_BLOB)) { fprintf(stderr, "  asm_clean selftest FAIL: tiny sheet classification\n"); fails++; }
    }
    Arena_dispose(&scratch);
    Arena_dispose(&persist);
    fails += ac_density_controls();
    if (fails == 0) fprintf(stderr, "  asm_clean selftest: all passed\n");
    return fails;
}
