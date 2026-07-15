/*
 * mesh_resplit.c -- re-LOP-from-original after a split. See mesh_resplit.h.
 *
 * MeshResplit_remesh_pieces is the reusable core: given a parent's original
 * point cloud and the sub-meshes it was split into, it labels each original
 * point to its nearest piece and re-LOP+re-BPA each piece from ITS original
 * points in isolation. MeshResplit_run applies it after the connectivity split
 * (step 2); pipeline_cube applies it again after bridge-cut/overlap (step 4).
 */
#include "mesh_resplit.h"

#include "../common/union_find.h"
#include "../common/kdtree.h"
#include "../common/mls_project.h"
#include "../common/vert_weld.h"
#include "../common/pca.h"
#include "../common/dump_obj.h"
#include "../common/pipeline_constants.h"
#include "../remesh/ball_pivot.h"
#include "../remesh/normal_orient.h"
#include "../remesh/orient_mesh.h"

#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifndef RESPLIT_MIN_COMP_VERTS
#define RESPLIT_MIN_COMP_VERTS 200   /* drop connectivity-components smaller than this */
#endif

/* ---- Seam-band pin (cross-cube convergence) -------------------------------
 * The extract MLS-projects owned verts within MLS_PROJECT_RADIUS_VOX of the cube
 * boundary using HALO (two-sided) support, so adjacent cubes converge to identical
 * seam geometry there. The re-LOP passes below re-MLS WITHOUT the halo (owned-only
 * support) and drift the seam independently per cube -- the welded "offset" / high-
 * lambda good-join. Pin verts within SEAM_PIN_BAND_VOX of the owned-box boundary
 * [0, RESPLIT_CUBE_D) at their extract positions through the re-LOP so the seam
 * stays where the halo put it. See SEAM_PIN_BAND_VOX in pipeline_constants.h. */
#define RESPLIT_CUBE_D 128.0   /* owned cube dimension (chunk_size, fixed) */
static double seam_pin_band(void)
{
    if (getenv("VES_SEAM_PIN_OFF")) return 0.0;
    const char *e = getenv("SEAM_PIN_BAND_VOX");
    if (e) { double v = atof(e); if (v >= 0.0) return v; }
    return (double)SEAM_PIN_BAND_VOX;
}
static int in_seam_band(const float *v, double band)
{
    return v[0] < band || v[0] > RESPLIT_CUBE_D - band ||
           v[1] < band || v[1] > RESPLIT_CUBE_D - band ||
           v[2] < band || v[2] > RESPLIT_CUBE_D - band;
}
/* Tapered pin weight: 1 at the cube face (full freeze), smoothstep down to 0 at
 * the band edge. A HARD freeze (1 inside the band, 0 outside) left a non-
 * developable CREASE at the band boundary -- 54% of the post-pin high-lambda verts
 * sat exactly there. Tapering blends the frozen seam into the re-LOP'd interior,
 * killing the crease while still holding the actual seam edge (weight ~1). */
static float seam_pin_weight(const float *v, double band)
{
    double dz = v[0] < RESPLIT_CUBE_D - v[0] ? v[0] : RESPLIT_CUBE_D - v[0];
    double dy = v[1] < RESPLIT_CUBE_D - v[1] ? v[1] : RESPLIT_CUBE_D - v[1];
    double dx = v[2] < RESPLIT_CUBE_D - v[2] ? v[2] : RESPLIT_CUBE_D - v[2];
    double d = dz; if (dy < d) d = dy; if (dx < d) d = dx;
    if (d >= band) return 0.0f;
    if (d < 0.0)   d = 0.0;
    double t = d / band;                       /* 0 at face, 1 at band edge */
    return (float)(1.0 - t*t*(3.0 - 2.0*t));   /* smoothstep: 1 -> 0 */
}

/* Diagnostic dump: materialize the stage dir "<dir>/<cube_id>_<stage>" and
 * return it in `out` (created on demand). Returns 0 if dumping is enabled. */
static int resplit_stage_dir(const MeshResplitDump *d, const char *stage,
                             char *out, size_t out_sz)
{
    if (!d || !d->dir || !d->cube_id || !stage) return -1;
    snprintf(out, out_sz, "%s/%s_%s", d->dir, d->cube_id, stage);
    DumpObj_ensure_dir(d->dir);
    DumpObj_ensure_dir(out);
    return 0;
}

/* ---- small mesh-cleanup helpers (mirror the statics in mesh_extract.c) ---- */

/* Drop connectivity sub-components with < MIN_FRAGMENT_FACES faces; compact
 * faces in place; return the kept face count. */
static size_t drop_small_fragments(Arena_T arena, size_t nv,
                                   int32_t *faces, size_t nf)
{
    if (nf == 0) return 0;
    UnionFind uf = UF_new(arena, (int32_t)nv);
    for (size_t i = 0; i < nf; i++) {
        uf_union(&uf, faces[i*3+0], faces[i*3+1]);
        uf_union(&uf, faces[i*3+1], faces[i*3+2]);
        uf_union(&uf, faces[i*3+0], faces[i*3+2]);
    }
    int32_t *fc = (int32_t *)ARENA_CALLOC(arena, (long)nv, (long)sizeof(int32_t));
    for (size_t i = 0; i < nf; i++) fc[uf_find(&uf, faces[i*3])]++;
    size_t kept = 0;
    for (size_t i = 0; i < nf; i++) {
        if (fc[uf_find(&uf, faces[i*3])] >= MIN_FRAGMENT_FACES) {
            faces[kept*3+0] = faces[i*3+0];
            faces[kept*3+1] = faces[i*3+1];
            faces[kept*3+2] = faces[i*3+2];
            kept++;
        }
    }
    return kept;
}

/* Remove unreferenced verts; remap faces in place; new compact verts out. */
static void compact_unreferenced(Arena_T arena,
                                 const float *old_v, size_t old_nv,
                                 int32_t *faces, size_t nf,
                                 float **new_v, size_t *new_nv)
{
    int32_t *remap = (int32_t *)ARENA_CALLOC(arena, (long)old_nv, (long)sizeof(int32_t));
    for (size_t i = 0; i < nf*3; i++) remap[faces[i]] = 1;
    int32_t count = 0;
    for (size_t i = 0; i < old_nv; i++) remap[i] = remap[i] ? count++ : -1;
    float *v = (float *)ARENA_ALLOC(arena, (long)((size_t)count*3*sizeof(float)));
    for (size_t i = 0; i < old_nv; i++) if (remap[i] >= 0) {
        v[remap[i]*3+0] = old_v[i*3+0];
        v[remap[i]*3+1] = old_v[i*3+1];
        v[remap[i]*3+2] = old_v[i*3+2];
    }
    for (size_t i = 0; i < nf*3; i++) faces[i] = remap[faces[i]];
    *new_v = v; *new_nv = (size_t)count;
}

static void finalize_mesh(ComponentMesh *m, float *verts, size_t nv,
                          int32_t *faces, size_t nf)
{
    float pca[3] = {0,0,1}, cen[3] = {0,0,0};
    PCA_normal(verts, nv, pca, cen);
    memset(m, 0, sizeof(*m));
    m->verts = verts; m->faces = faces; m->nv = nv; m->nf = nf;
    m->pca_normal[0]=pca[0]; m->pca_normal[1]=pca[1]; m->pca_normal[2]=pca[2];
    m->centroid[0]=cen[0]; m->centroid[1]=cen[1]; m->centroid[2]=cen[2];
    m->self = m;
}

/* Build a ComponentMesh from the faces of `src` whose vertices belong to UF root
 * `root` (the original geometry of one connectivity-component). */
static int extract_submesh(Arena_T arena, const ComponentMesh *src,
                           UnionFind *uf, int32_t root, ComponentMesh *out)
{
    size_t nf = 0;
    for (size_t f = 0; f < src->nf; f++)
        if (uf_find(uf, src->faces[f*3]) == root) nf++;
    if (nf == 0) return -1;
    int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(nf*3*sizeof(int32_t)));
    size_t k = 0;
    for (size_t f = 0; f < src->nf; f++)
        if (uf_find(uf, src->faces[f*3]) == root) {
            faces[k*3+0] = src->faces[f*3+0];
            faces[k*3+1] = src->faces[f*3+1];
            faces[k*3+2] = src->faces[f*3+2];
            k++;
        }
    float *fv = NULL; size_t fnv = 0;
    compact_unreferenced(arena, src->verts, src->nv, faces, nf, &fv, &fnv);
    if (fnv < 3) return -1;
    finalize_mesh(out, fv, fnv, faces, nf);
    return 0;
}

/* Re-LOP (from original points, in isolation) + weld + Hoppe + BPA + winding +
 * cleanup -> a clean ComponentMesh. Mirrors extract's point->mesh chain.
 * Returns 0 on success, -1 (out untouched) on any failure. */
static int relop_rebpa(Arena_T arena, const float *pts, size_t n,
                       const float cell_origin[3],
                       const MeshResplitDump *dump, const char *tag,
                       int lop_iters, const float *pin_w, ComponentMesh *out)
{
    if (n < 3) return -1;
    float *pv = (float *)ARENA_ALLOC(arena, (long)(n*3*sizeof(float)));
    float *pn = (float *)ARENA_ALLOC(arena, (long)(n*3*sizeof(float)));
    float *sn = (float *)ARENA_ALLOC(arena, (long)(n*3*sizeof(float)));
    float *sv = (lop_iters >= 2)
        ? (float *)ARENA_ALLOC(arena, (long)(n*3*sizeof(float))) : NULL;
    const float *src = pts;
    float *dst = pv;
    for (int it = 0; it < lop_iters; it++) {
        float *ndst = (it == lop_iters - 1) ? pn : sn;
        MLS_project_verts(arena, src, n, MLS_PROJECT_RADIUS_VOX,
                          cell_origin, dst, ndst);
        /* Blend pinned (seam-band) points back toward their extract positions --
         * they were MLS'd with halo support this re-LOP lacks, so re-projecting
         * them freely drifts the cross-cube seam. The weight tapers from 1 at the
         * cube face to 0 at the band edge (no hard crease). They still serve as MLS
         * support for the (free) interior points being projected. */
        if (pin_w) {
            for (size_t i = 0; i < n; i++) {
                float w = pin_w[i];
                if (w > 0.0f) {
                    float iw = 1.0f - w;
                    dst[i*3+0] = w*pts[i*3+0] + iw*dst[i*3+0];
                    dst[i*3+1] = w*pts[i*3+1] + iw*dst[i*3+1];
                    dst[i*3+2] = w*pts[i*3+2] + iw*dst[i*3+2];
                }
            }
        }
        src = dst;
        dst = (dst == pv) ? sv : pv;
    }
    if (src != pv) memcpy(pv, src, n*3*sizeof(float));

    /* Point-only weld to a fixed point (collapse LOP duplicates), as extract. */
    int32_t dummy_face = 0; size_t dummy_nf = 0;
    float *wv = pv, *wn = pn; size_t wnv = n;
    {
        const float *sv2 = pv; const float *sn2 = pn; size_t snv = n;
        for (int iter = 0; iter < 5; iter++) {
            size_t prev = wnv;
            Weld_verts(arena, sv2, snv, sn2, &dummy_face, 0, &dummy_nf,
                       MLS_WELD_EPS_VOX, /*guard_orient=*/false, &wv, &wnv, &wn);
            if (wnv == prev && iter > 0) break;
            sv2 = wv; sn2 = wn; snv = wnv;
        }
    }
    if (wnv < 3) return -1;

    if (dump && dump->mls_stage) {
        char sd[1024];
        if (resplit_stage_dir(dump, dump->mls_stage, sd, sizeof sd) == 0) {
            char path[1280];
            snprintf(path, sizeof path, "%s/%s_%s_%s.obj",
                     sd, dump->cube_id, dump->mls_stage, tag ? tag : "x");
            DumpObj_write_points_world(arena, path, dump->cube_id, wv, NULL, wnv);
            /* pre-MLS: the raw input points fed to this re-LOP (pre-projection,
             * pre-weld), as a sibling _pre file. */
            snprintf(path, sizeof path, "%s/%s_%s_%s_pre.obj",
                     sd, dump->cube_id, dump->mls_stage, tag ? tag : "x");
            DumpObj_write_points_world(arena, path, dump->cube_id, pts, NULL, n);
        }
    }

    /* Hoppe consistent normal sign (BPA's per-vertex normal gate needs it). */
    { size_t hf = 0; NormalOrient_consistent(arena, wv, wn, wnv, 1.2f, &hf); }

    int32_t *faces = NULL; size_t nf = 0;
    if (BallPivot_reconstruct(arena, wv, wn, wnv, BPA_RHO_VOX, &faces, &nf) != 0
        || nf == 0)
        return -1;

    OrientMesh_consistent(arena, wv, wnv, wn, faces, nf, NULL, NULL, NULL);

    nf = drop_small_fragments(arena, wnv, faces, nf);
    if (nf == 0) return -1;
    float *fv = NULL; size_t fnv = 0;
    compact_unreferenced(arena, wv, wnv, faces, nf, &fv, &fnv);
    if (fnv < 3) return -1;

    finalize_mesh(out, fv, fnv, faces, nf);
    return 0;
}

/* Pass `pieces` through to *out unchanged (no re-LOP possible / no cloud). */
static int pass_pieces_through(Arena_T arena,
                               const ComponentMesh *pieces, size_t n_pieces,
                               ComponentMesh **out, size_t *n_out,
                               MeshResplitCloud **out_clouds)
{
    ComponentMesh *o = (ComponentMesh *)ARENA_ALLOC(arena,
                           (long)(n_pieces * sizeof(ComponentMesh)));
    MeshResplitCloud *oc = out_clouds
        ? (MeshResplitCloud *)ARENA_ALLOC(arena, (long)(n_pieces * sizeof(MeshResplitCloud)))
        : NULL;
    for (size_t p = 0; p < n_pieces; p++) {
        o[p] = pieces[p]; o[p].self = &o[p]; o[p].comp_id = (int)(p + 1);
        if (oc) memset(&oc[p], 0, sizeof(oc[p]));
    }
    *out = o; *n_out = n_pieces;
    if (out_clouds) *out_clouds = oc;
    return 0;
}

/* Re-LOP + re-BPA a mesh from its OWN vertices (no cloud relabel, so no
 * vacuuming/fold) -- the clean re-triangulation QEM needs. See mesh_resplit.h. */
int MeshResplit_resurface_own(Arena_T arena, ComponentMesh *m,
                              const float cell_origin[3], int iters)
{
    assert(arena && m && m->self == m);
    if (m->nv < 3 || iters <= 0) return 0;
    float cori[3] = { 0.0f, 0.0f, 0.0f };
    if (cell_origin) { cori[0] = cell_origin[0]; cori[1] = cell_origin[1]; cori[2] = cell_origin[2]; }

    /* Pin the seam band: these owned verts were placed with halo two-sided support
     * at extract, and re-LOPing them from this piece's owned-only verts would drift
     * the cross-cube seam. Keep them fixed; the interior still re-LOPs cleanly. */
    float *pinw = NULL;
    double band = seam_pin_band();
    if (band > 0.0) {
        pinw = (float *)ARENA_CALLOC(arena, (long)m->nv, (long)sizeof(float));
        for (size_t i = 0; i < m->nv; i++)
            pinw[i] = seam_pin_weight(&m->verts[i*3], band);
    }

    /* re-LOP + weld + Hoppe + Ball-Pivot + orient on the piece's OWN verts */
    ComponentMesh nm;
    if (relop_rebpa(arena, m->verts, m->nv, cori, NULL, NULL, iters, pinw, &nm) != 0)
        return 0;                       /* keep m unchanged on failure */

    int cid = m->comp_id;
    nm.comp_id      = cid;
    nm.pin_mask     = NULL;             /* fresh vertex set -> old per-vert masks invalid */
    nm.vert_normals = NULL;
    nm.nv_pre_fill  = 0;
    *m = nm;
    m->self = m;
    return 0;
}

int MeshResplit_remesh_pieces(Arena_T arena,
                              const MeshResplitCloud *cloud,
                              const ComponentMesh    *pieces, size_t n_pieces,
                              const char             *dump_tag,
                              const MeshResplitDump  *dump,
                              ComponentMesh         **out, size_t *n_out,
                              MeshResplitCloud      **out_clouds)
{
    assert(arena && out && n_out);
    *out = NULL; *n_out = 0;
    if (out_clouds) *out_clouds = NULL;
    if (n_pieces == 0) return 0;

    size_t totv = 0;
    for (size_t p = 0; p < n_pieces; p++) totv += pieces[p].nv;
    if (!cloud || cloud->n == 0 || !cloud->orig_pts || !cloud->lop_pts || totv == 0)
        return pass_pieces_through(arena, pieces, n_pieces, out, n_out, out_clouds);

    /* One KD-tree per piece. Assign each parent original point to the piece its
     * LOP image is nearest to, but ONLY when that piece is unambiguously closest
     * -- the next-nearest piece must be at least MLS_RESPLIT_ASSIGN_MARGIN_VOX
     * farther. Points that straddle the split seam (or sit on an adjacent
     * close-wrap ~2-3 vox away) are DROPPED instead of vacuumed into the piece;
     * vacuuming them is what grew the piece ~60% and folded its re-BPA'd sheet
     * (the step7_cc_bpa fold). margin = 0 reverts to greedy nearest-vertex. */
    KDTree_T *ptree = (KDTree_T *)ARENA_ALLOC(arena, (long)(n_pieces*sizeof(KDTree_T)));
    for (size_t p = 0; p < n_pieces; p++)
        ptree[p] = KDTree_new(arena, pieces[p].verts, pieces[p].nv);

    const float margin = MLS_RESPLIT_ASSIGN_MARGIN_VOX;
    size_t  *bsz    = (size_t *)ARENA_CALLOC(arena, (long)n_pieces, (long)sizeof(size_t));
    int32_t *plabel = (int32_t *)ARENA_ALLOC(arena, (long)(cloud->n*sizeof(int32_t)));
    size_t n_drop = 0;
    for (size_t j = 0; j < cloud->n; j++) {
        float best = FLT_MAX, second = FLT_MAX;
        int32_t bp = -1;
        for (size_t p = 0; p < n_pieces; p++) {
            float d2 = FLT_MAX;
            (void)KDTree_nearest(ptree[p], &cloud->lop_pts[j*3], &d2);
            float d = sqrtf(d2);
            if (d < best)        { second = best; best = d; bp = (int32_t)p; }
            else if (d < second) { second = d; }
        }
        if (bp >= 0 && (second - best) >= margin) {
            plabel[j] = bp; bsz[(size_t)bp]++;
        } else {
            plabel[j] = -1; n_drop++;   /* seam-ambiguous: drop, don't vacuum */
        }
    }

    /* gather each piece's original + LOP point subsets (dropped points skipped) */
    float **porig = (float **)ARENA_ALLOC(arena, (long)(n_pieces*sizeof(float*)));
    float **plop  = (float **)ARENA_ALLOC(arena, (long)(n_pieces*sizeof(float*)));
    size_t *bpos  = (size_t *)ARENA_CALLOC(arena, (long)n_pieces, (long)sizeof(size_t));
    for (size_t p = 0; p < n_pieces; p++) {
        size_t m = bsz[p] ? bsz[p] : 1;
        porig[p] = (float *)ARENA_ALLOC(arena, (long)(m*3*sizeof(float)));
        plop[p]  = (float *)ARENA_ALLOC(arena, (long)(m*3*sizeof(float)));
    }
    for (size_t j = 0; j < cloud->n; j++) {
        if (plabel[j] < 0) continue;
        size_t p = (size_t)plabel[j];
        memcpy(&porig[p][bpos[p]*3], &cloud->orig_pts[j*3], 3*sizeof(float));
        memcpy(&plop[p][bpos[p]*3],  &cloud->lop_pts[j*3],  3*sizeof(float));
        bpos[p]++;
    }
    if (n_drop)
        fprintf(stderr,
                "    re-LOP assign: %zu pieces, dropped %zu/%zu seam-ambiguous "
                "points (margin %.2f vox)\n",
                n_pieces, n_drop, cloud->n, (double)margin);

    /* diagnostic: the connectivity/piece split that feeds the re-LOP */
    if (dump && dump->cc_stage) {
        char sd[1024];
        if (resplit_stage_dir(dump, dump->cc_stage, sd, sizeof sd) == 0) {
            static const float pal[6][3] =
                {{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1}};
            for (size_t p = 0; p < n_pieces; p++) {
                char path[1280];
                snprintf(path, sizeof path, "%s/%s_%s_%s_p%02zu.obj", sd,
                         dump->cube_id, dump->cc_stage, dump_tag ? dump_tag : "x", p);
                DumpObj_write_one_world(arena, path, dump->cube_id,
                                        pieces[p].verts, pieces[p].nv,
                                        pieces[p].faces, pieces[p].nf, pal[p % 6]);
            }
        }
    }

    ComponentMesh   *o  = (ComponentMesh *)ARENA_ALLOC(arena,
                              (long)(n_pieces*sizeof(ComponentMesh)));
    MeshResplitCloud *oc = out_clouds
        ? (MeshResplitCloud *)ARENA_ALLOC(arena, (long)(n_pieces*sizeof(MeshResplitCloud)))
        : NULL;
    /* Only a parent that actually SPLIT (yielded >1 pieces) gets the higher
     * re-LOP pass count -- a pass-through single piece keeps the lighter extract
     * count so it is not over-smoothed (which would flatten real curvature). */
    int lop_iters = (n_pieces > 1) ? MLS_RESPLIT_ITERS : MLS_PROJECT_ITERS;
    double band = seam_pin_band();
    size_t cnt = 0;
    for (size_t p = 0; p < n_pieces; p++) {
        char tag[96];
        snprintf(tag, sizeof(tag), "%s_p%02zu", dump_tag ? dump_tag : "x", p);

        /* Re-LOP input = the labelled cloud points PLUS this piece's extract
         * seam-band verts, PINNED. The seam-ambiguous-point drop above shed the
         * halo support exactly at the seam; re-adding the extract seam verts as
         * fixed anchors keeps the cross-cube seam at its halo-converged positions
         * while the (non-pinned) cloud points re-LOP the interior. */
        const float *rpts = porig[p];
        size_t rn = bsz[p];
        const float *rpin = NULL;
        if (band > 0.0 && bsz[p] > 0) {
            size_t ns = 0;
            for (size_t v = 0; v < pieces[p].nv; v++)
                if (in_seam_band(&pieces[p].verts[v*3], band)) ns++;
            if (ns > 0) {
                size_t ncomb = ns + bsz[p];
                float *comb = (float *)ARENA_ALLOC(arena, (long)(ncomb*3*sizeof(float)));
                float *pinc = (float *)ARENA_CALLOC(arena, (long)ncomb, (long)sizeof(float));
                size_t w = 0;
                for (size_t v = 0; v < pieces[p].nv; v++)
                    if (in_seam_band(&pieces[p].verts[v*3], band)) {
                        comb[w*3+0] = pieces[p].verts[v*3+0];
                        comb[w*3+1] = pieces[p].verts[v*3+1];
                        comb[w*3+2] = pieces[p].verts[v*3+2];
                        pinc[w] = seam_pin_weight(&pieces[p].verts[v*3], band);
                        w++;
                    }
                /* re-LOP only the INTERIOR cloud points; the seam band is owned by
                 * the pinned extract verts above. Including band cloud points too
                 * would lay a second, drifting layer over the pinned seam (a
                 * through-thickness double sheet that raises lambda), so drop them. */
                size_t w2 = ns;
                for (size_t j = 0; j < bsz[p]; j++)
                    if (!in_seam_band(&porig[p][j*3], band)) {
                        comb[w2*3+0] = porig[p][j*3+0];
                        comb[w2*3+1] = porig[p][j*3+1];
                        comb[w2*3+2] = porig[p][j*3+2];
                        w2++;
                    }
                rpts = comb; rn = w2; rpin = pinc;
            }
        }

        ComponentMesh m;
        if (relop_rebpa(arena, rpts, rn, cloud->cell_origin,
                        dump, tag, lop_iters, rpin, &m) != 0) {
            m = pieces[p];                 /* fall back to the piece unchanged */
            float fpca[3] = {0,0,1}, fcen[3] = {0,0,0};   /* keep centroid/pca valid */
            PCA_normal(m.verts, m.nv, fpca, fcen);
            m.pca_normal[0]=fpca[0]; m.pca_normal[1]=fpca[1]; m.pca_normal[2]=fpca[2];
            m.centroid[0]=fcen[0]; m.centroid[1]=fcen[1]; m.centroid[2]=fcen[2];
        }
        o[cnt] = m; o[cnt].self = &o[cnt]; o[cnt].comp_id = (int)(cnt + 1);
        if (oc) {
            oc[cnt].orig_pts = porig[p];
            oc[cnt].lop_pts  = plop[p];
            oc[cnt].n        = bsz[p];
            oc[cnt].cell_origin[0] = cloud->cell_origin[0];
            oc[cnt].cell_origin[1] = cloud->cell_origin[1];
            oc[cnt].cell_origin[2] = cloud->cell_origin[2];
        }
        cnt++;
    }
    *out = o; *n_out = cnt;
    if (out_clouds) *out_clouds = oc;
    return 0;
}

int MeshResplit_run(Arena_T arena,
                    const ComponentMesh   *in, size_t n_in,
                    const MeshResplitCloud *clouds,
                    int                    n_threads,
                    const MeshResplitDump *dump,
                    ComponentMesh        **out, size_t *n_out,
                    MeshResplitCloud     **out_clouds)
{
    (void)n_threads;
    assert(arena && out && n_out);
    *out = NULL; *n_out = 0;
    if (out_clouds) *out_clouds = NULL;
    if (n_in == 0) return 0;

    size_t cap = n_in + 8, cnt = 0;
    ComponentMesh    *acc  = (ComponentMesh *)ARENA_ALLOC(arena,
                                 (long)(cap * sizeof(ComponentMesh)));
    MeshResplitCloud *accc = out_clouds
        ? (MeshResplitCloud *)ARENA_ALLOC(arena, (long)(cap * sizeof(MeshResplitCloud)))
        : NULL;
    MeshResplitCloud zerocl; memset(&zerocl, 0, sizeof zerocl);

    #define ENSURE_CAP() do { if (cnt >= cap) { size_t nc = cap * 2;            \
        ComponentMesh *na = (ComponentMesh *)ARENA_ALLOC(arena,                 \
            (long)(nc * sizeof(ComponentMesh)));                               \
        memcpy(na, acc, cnt * sizeof(ComponentMesh));                          \
        for (size_t _i = 0; _i < cnt; _i++) na[_i].self = &na[_i]; acc = na;    \
        if (accc) { MeshResplitCloud *nq = (MeshResplitCloud *)ARENA_ALLOC(     \
                arena, (long)(nc * sizeof(MeshResplitCloud)));                  \
            memcpy(nq, accc, cnt * sizeof(MeshResplitCloud)); accc = nq; }      \
        cap = nc; } } while (0)
    #define EMIT(MSH, CLD) do { ENSURE_CAP();                                   \
        acc[cnt] = (MSH); acc[cnt].self = &acc[cnt];                            \
        acc[cnt].comp_id = (int)(cnt + 1);                                      \
        if (accc) accc[cnt] = (CLD); cnt++; } while (0)

    for (size_t i = 0; i < n_in; i++) {
        const ComponentMesh    *cm = &in[i];
        const MeshResplitCloud *cl = clouds ? &clouds[i] : NULL;

        /* connectivity-split this mesh's vertices over its faces */
        UnionFind uf = UF_new(arena, (int32_t)cm->nv);
        for (size_t f = 0; f < cm->nf; f++) {
            uf_union(&uf, cm->faces[f*3+0], cm->faces[f*3+1]);
            uf_union(&uf, cm->faces[f*3+1], cm->faces[f*3+2]);
            uf_union(&uf, cm->faces[f*3+0], cm->faces[f*3+2]);
        }
        int32_t *vcount = (int32_t *)ARENA_CALLOC(arena, (long)cm->nv,
                                                  (long)sizeof(int32_t));
        for (size_t v = 0; v < cm->nv; v++) vcount[uf_find(&uf, (int32_t)v)]++;
        int32_t *compid = (int32_t *)ARENA_ALLOC(arena, (long)(cm->nv*sizeof(int32_t)));
        for (size_t v = 0; v < cm->nv; v++) compid[v] = -1;
        int32_t n_real = 0;
        for (size_t v = 0; v < cm->nv; v++) {
            int32_t r = uf_find(&uf, (int32_t)v);
            if ((size_t)r == v && vcount[r] >= RESPLIT_MIN_COMP_VERTS)
                compid[r] = n_real++;
        }

        if (n_real <= 1 || !cl || cl->n == 0 || !cl->orig_pts || !cl->lop_pts) {
            EMIT(*cm, cl ? *cl : zerocl);          /* pass through with its cloud */
            continue;
        }

        /* build connectivity sub-meshes for the real components */
        int32_t *root_of = (int32_t *)ARENA_ALLOC(arena, (long)(n_real*sizeof(int32_t)));
        for (size_t v = 0; v < cm->nv; v++)
            if (compid[v] >= 0) root_of[compid[v]] = (int32_t)v;
        ComponentMesh *pieces = (ComponentMesh *)ARENA_ALLOC(arena,
                                    (long)(n_real*sizeof(ComponentMesh)));
        size_t np = 0;
        for (int32_t c = 0; c < n_real; c++)
            if (extract_submesh(arena, cm, &uf, root_of[c], &pieces[np]) == 0) np++;
        if (np == 0) { EMIT(*cm, cl ? *cl : zerocl); continue; }

        /* re-LOP + re-BPA each connectivity-component from its original points */
        char tag[32]; snprintf(tag, sizeof(tag), "in%02zu", i);
        ComponentMesh   *po  = NULL; size_t npo = 0; MeshResplitCloud *poc = NULL;
        MeshResplit_remesh_pieces(arena, cl, pieces, np, tag, dump,
                                  &po, &npo, out_clouds ? &poc : NULL);
        for (size_t k = 0; k < npo; k++)
            EMIT(po[k], (poc ? poc[k] : zerocl));
    }

    #undef EMIT
    #undef ENSURE_CAP
    *out = acc; *n_out = cnt;
    if (out_clouds) *out_clouds = accc;
    return 0;
}
