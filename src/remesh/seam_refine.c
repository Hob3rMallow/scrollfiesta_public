/*
 * seam_refine.c -- weld-time seam-band refinement (see seam_refine.h).
 *
 * Split mechanics mirror isotropic_remesh.c::split_round (independent-set
 * face-locked midpoint splits, winding-preserving) with the boundary-edge
 * (single-face) case added -- split_round structurally never splits boundary
 * edges, but the seam rim IS boundary, so the rim case is the whole point
 * here. Flip relief between rounds reuses WeldCleanup_flip_rounds (boundary-
 * frozen, min-angle, normal-guarded).
 *
 * Scratch is arena-allocated under the caller's arena; grown vert/face arrays
 * are fresh arena allocations each round (the previous round's arrays become
 * arena garbage, reclaimed when the caller disposes the weld arena). No
 * Arena_save/restore across the growth allocations (the twice-hit aliasing
 * rule: outputs must never sit above a restored mark).
 */
#include "seam_refine.h"
#include "weld_cleanup.h"
#include "../common/pipeline_constants.h"
#include "../common/u64_radix.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include "../common/ves_omp.h"
#endif

static double edge_len3(const float *V, int32_t a, int32_t b)
{
    double dx=(double)V[(size_t)a*3+0]-V[(size_t)b*3+0];
    double dy=(double)V[(size_t)a*3+1]-V[(size_t)b*3+1];
    double dz=(double)V[(size_t)a*3+2]-V[(size_t)b*3+2];
    return sqrt(dx*dx+dy*dy+dz*dz);
}

static double tri_area3(const float *V, int32_t a, int32_t b, int32_t c)
{
    double e1x=(double)V[(size_t)b*3+0]-V[(size_t)a*3+0];
    double e1y=(double)V[(size_t)b*3+1]-V[(size_t)a*3+1];
    double e1z=(double)V[(size_t)b*3+2]-V[(size_t)a*3+2];
    double e2x=(double)V[(size_t)c*3+0]-V[(size_t)a*3+0];
    double e2y=(double)V[(size_t)c*3+1]-V[(size_t)a*3+1];
    double e2z=(double)V[(size_t)c*3+2]-V[(size_t)a*3+2];
    double cx=e1y*e2z-e1z*e2y, cy=e1z*e2x-e1x*e2z, cz=e1x*e2y-e1y*e2x;
    return 0.5*sqrt(cx*cx+cy*cy+cz*cz);
}

/* min altitude = 2*area / longest edge (deliberate local copy, see
 * weld_cleanup.c header note on self-containment). */
static double tri_min_alt3(const float *V, int32_t a, int32_t b, int32_t c)
{
    double l0=edge_len3(V,a,b), l1=edge_len3(V,b,c), l2=edge_len3(V,c,a);
    double lmax = l0>l1?(l0>l2?l0:l2):(l1>l2?l1:l2);
    if (lmax < 1e-12) return 0.0;
    return 2.0*tri_area3(V,a,b,c)/lmax;
}

/* Key first so the generic record radix can move the complete half-edge while
 * sorting on (min_vertex,max_vertex).  face_dir packs the direction bit into
 * the low bit; grid_weld's face budget is far below the remaining 31 bits. */
typedef VesU64Record16 RHE;

static uint64_t rhe_key(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)(a < b ? a : b);
    uint32_t hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static int32_t rhe_v0(const RHE *h) { return (int32_t)(h->key >> 32); }
static int32_t rhe_v1(const RHE *h) { return (int32_t)(uint32_t)h->key; }
static int32_t rhe_face(const RHE *h)
{ return (int32_t)((uint32_t)h->data >> 1); }
static int rhe_fwd(const RHE *h) { return (int)((uint32_t)h->data & 1u); }

static int rhe_cmp(const void *pa, const void *pb)
{
    const RHE *a=(const RHE*)pa, *b=(const RHE*)pb;
    return a->key < b->key ? -1 : a->key > b->key ? 1 : 0;
}

static int seam_work_threads(size_t n, size_t grain)
{
#ifdef _OPENMP
    int nt = omp_get_max_threads();
    size_t useful = (n + grain - 1) / grain;
    if (useful < 1) useful = 1;
    if ((size_t)nt > useful) nt = (int)useful;
    return nt > 0 ? nt : 1;
#else
    (void)n; (void)grain;
    return 1;
#endif
}

/* One split request: interior (fd >= 0) or boundary (fd == -1). fc traverses
 * the edge a->b in face order; fd (if any) traverses b->a. */
typedef struct { int32_t a, b, fc, fd; } RSplit;

void SeamRefine_default_params(SeamRefineParams *p)
{
    assert(p);
    p->target_len      = SEAM_REFINE_TARGET_VOX;
    p->band            = 6.0f;                 /* the bridge's seam band */
    p->min_parent_alt  = SEAM_SLIVER_MIN_ALT;
    p->max_rounds      = SEAM_REFINE_MAX_ROUNDS;
    p->flip_max_rounds = WELD_CLEANUP_FLIP_ROUNDS;
}

static int seam_refine_process_impl(Arena_T arena,
                       const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       const SeamPlane *planes, size_t np,
                       const SeamRefineParams *params,
                       const uint8_t *freeze_source_faces,
                       size_t freeze_source_nf,
                       float **out_verts, size_t *out_nv,
                       int32_t **out_faces, size_t *out_nf,
                       int32_t **out_new_vert_parent0,
                       int32_t **out_new_vert_parent1,
                       int32_t **out_face_source,
                       size_t *out_n_new,
                       SeamRefineStats *st)
{
    SeamRefineParams p;
    float   *V = NULL;
    int32_t *F = NULL;
    int32_t *src0 = NULL, *src1 = NULL; /* per NEW vert: edge endpoints */
    int32_t *face_source = NULL;         /* final face -> input source face */
    uint8_t *face_frozen = NULL;         /* propagated source-face freeze */
    size_t   cnv = nv, cnf = nf, nsrc = 0, src_cap = 0;
    int      round = 0;

    assert(arena);
    assert(out_verts && out_nv && out_faces && out_nf);
    assert(out_new_vert_parent0 && out_new_vert_parent1 && out_n_new);

    if (params) p = *params; else SeamRefine_default_params(&p);
    if (st) memset(st, 0, sizeof *st);
    if (st) st->planes = np;

    /* No-op path: hand back the inputs by reference. */
    *out_verts = (float *)verts;   *out_nv = nv;
    *out_faces = (int32_t *)faces; *out_nf = nf;
    *out_new_vert_parent0 = NULL;
    *out_new_vert_parent1 = NULL;  *out_n_new = 0;
    if (out_face_source) *out_face_source = NULL;
    if (nf == 0 || nv == 0 || np == 0 || planes == NULL) return 0;
    if (p.target_len <= 0.0f || p.max_rounds <= 0) return 0;
    if (freeze_source_faces && freeze_source_nf != nf) return -1;

    for (round = 0; round < p.max_rounds; round++) {
        const float   *cV = V ? V : verts;
        const int32_t *cF = F ? F : faces;
        size_t n_he = 0, i = 0, nreq = 0;
        size_t bnd_this = 0, int_this = 0;
        int nt = seam_work_threads(cnf > cnv ? cnf : cnv, 65536u);
        uint8_t *vband = (uint8_t *)malloc(cnv ? cnv : 1);
        size_t *edge_off = (size_t *)calloc((size_t)nt + 1, sizeof(*edge_off));
        RHE *he = NULL;
        uint8_t *fdone = NULL;
        RSplit *req = NULL;
        if (vband == NULL || edge_off == NULL ||
            cnf > (size_t)UINT32_MAX / 2u) {
            free(vband); free(edge_off);
            return -1;
        }

        /* The old loop evaluated distance-to-any-seam twice for every unique
         * edge.  Cache it once per vertex, then omit half-edges whose endpoints
         * both fail the exact same predicate.  Every incident copy of a kept
         * edge has the same endpoints, so multiplicity and boundary detection
         * are unchanged for every edge the refiner can act on. */
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            size_t lo = cnv * (size_t)tid / (size_t)nt;
            size_t hi = cnv * (size_t)(tid + 1) / (size_t)nt;
            size_t v;
            for (v = lo; v < hi; v++)
                vband[v] = (uint8_t)
                    (SeamPlanes_vert_dist(cV, (int32_t)v, planes, np)
                        <= (double)p.band);
        }
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
        {
            int tid = 0;
            size_t count = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            size_t lo = cnf * (size_t)tid / (size_t)nt;
            size_t hi = cnf * (size_t)(tid + 1) / (size_t)nt;
            size_t f;
            for (f = lo; f < hi; f++) {
                int32_t v[3] = {
                    cF[f*3+0], cF[f*3+1], cF[f*3+2]
                };
                int e;
                for (e = 0; e < 3; e++) {
                    int32_t a = v[e], b = v[(e+1)%3];
                    count += (size_t)(vband[(size_t)a] ||
                                      vband[(size_t)b]);
                }
            }
            edge_off[(size_t)tid + 1] = count;
        }
        {
            int t;
            for (t = 0; t < nt; t++)
                edge_off[(size_t)t + 1] += edge_off[(size_t)t];
        }
        n_he = edge_off[(size_t)nt];
        he = (RHE *)malloc((n_he ? n_he : 1) * sizeof(*he));
        fdone = (uint8_t *)calloc(cnf ? cnf : 1, 1);
        req = (RSplit *)malloc((n_he/2 + 2) * sizeof(*req));
        if (he == NULL || fdone == NULL || req == NULL) {
            free(vband); free(edge_off);
            free(he); free(fdone); free(req);
            return -1;
        }
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
#endif
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            size_t lo = cnf * (size_t)tid / (size_t)nt;
            size_t hi = cnf * (size_t)(tid + 1) / (size_t)nt;
            size_t pos = edge_off[(size_t)tid];
            size_t f;
            for (f = lo; f < hi; f++) {
                int32_t v[3] = {
                    cF[f*3+0], cF[f*3+1], cF[f*3+2]
                };
                int e;
                for (e = 0; e < 3; e++) {
                    int32_t a = v[e], b = v[(e+1)%3];
                    if (!vband[(size_t)a] && !vband[(size_t)b]) continue;
                    he[pos].key = rhe_key(a, b);
                    he[pos].data =
                        (uint64_t)(((uint32_t)f << 1) |
                                   (uint32_t)(a < b));
                    pos++;
                }
            }
            assert(pos == edge_off[(size_t)tid + 1]);
        }
        free(vband);
        free(edge_off);
        if (ves_u64_record16_sort(he, n_he) != 0)
            qsort(he, n_he, sizeof(*he), rhe_cmp);

        while (i < n_he) {
            size_t j = i + 1;
            while (j < n_he && he[j].key == he[i].key) j++;
            size_t run = j - i;
            int32_t a = rhe_v0(&he[i]), b = rhe_v1(&he[i]);
            do {
                if (run > 2) break;                       /* non-manifold: leave alone */
                if (edge_len3(cV, a, b) <= (double)p.target_len) break;
                if (run == 2) {
                    /* interior edge: needs oppositely-wound faces (same_dir
                     * pairs are left for the orient passes, as split_round) */
                    int32_t fc, fd;
                    if (rhe_fwd(&he[i]) == rhe_fwd(&he[i+1])) break;
                    if (rhe_fwd(&he[i])) {
                        fc = rhe_face(&he[i]);
                        fd = rhe_face(&he[i+1]);
                    } else {
                        fc = rhe_face(&he[i+1]);
                        fd = rhe_face(&he[i]);
                    }
                    if ((face_frozen &&
                         (face_frozen[(size_t)fc] ||
                          face_frozen[(size_t)fd])) ||
                        (!face_frozen && freeze_source_faces &&
                         (freeze_source_faces[(size_t)fc] ||
                          freeze_source_faces[(size_t)fd])))
                        break;
                    if (fdone[fc] || fdone[fd]) break;
                    if (tri_min_alt3(cV, cF[fc*3+0], cF[fc*3+1], cF[fc*3+2])
                            < (double)p.min_parent_alt) break;
                    if (tri_min_alt3(cV, cF[fd*3+0], cF[fd*3+1], cF[fd*3+2])
                            < (double)p.min_parent_alt) break;
                    req[nreq].a = a; req[nreq].b = b;
                    req[nreq].fc = fc; req[nreq].fd = fd;
                    nreq++; fdone[fc] = 1; fdone[fd] = 1; int_this++;
                } else {
                    /* boundary edge (run 1): split its single face. Recover the
                     * DIRECTED orientation: fwd means the face traverses a->b
                     * (v0->v1); else b->a. Normalize so fc traverses a->b. */
                    int32_t fc = rhe_face(&he[i]);
                    int32_t da = rhe_fwd(&he[i]) ? a : b;
                    int32_t db = rhe_fwd(&he[i]) ? b : a;
                    if ((face_frozen && face_frozen[(size_t)fc]) ||
                        (!face_frozen && freeze_source_faces &&
                         freeze_source_faces[(size_t)fc]))
                        break;
                    if (fdone[fc]) break;
                    if (tri_min_alt3(cV, cF[fc*3+0], cF[fc*3+1], cF[fc*3+2])
                            < (double)p.min_parent_alt) break;
                    req[nreq].a = da; req[nreq].b = db;
                    req[nreq].fc = fc; req[nreq].fd = -1;
                    nreq++; fdone[fc] = 1; bnd_this++;
                }
            } while (0);
            i = j;
        }

        if (nreq == 0) { free(he); free(fdone); free(req); break; }

        /* Grow: +1 vert per split; +2 faces per interior, +1 per boundary. */
        {
            size_t add_f = 2*int_this + bnd_this;
            size_t new_nv = cnv + nreq, new_nf = cnf + add_f;
            float   *nV = (float *)ARENA_ALLOC(arena,
                              (new_nv*3*sizeof(float)));
            int32_t *nF = (int32_t *)ARENA_ALLOC(arena,
                              (new_nf*3*sizeof(int32_t)));
            int32_t *nFaceSource = out_face_source
                ? (int32_t *)ARENA_ALLOC(
                      arena, (new_nf*sizeof(int32_t)))
                : NULL;
            uint8_t *nFaceFrozen = freeze_source_faces
                ? (uint8_t *)ARENA_ALLOC(arena, new_nf)
                : NULL;
            memcpy(nV, cV, cnv*3*sizeof(float));
            memcpy(nF, cF, cnf*3*sizeof(int32_t));
            if (nFaceSource) {
                if (face_source)
                    memcpy(nFaceSource, face_source,
                           cnf*sizeof(int32_t));
                else
                    for (size_t sf = 0; sf < cnf; sf++)
                        nFaceSource[sf] = (int32_t)sf;
            }
            if (nFaceFrozen) {
                if (face_frozen)
                    memcpy(nFaceFrozen, face_frozen,
                           cnf*sizeof(uint8_t));
                else
                    memcpy(nFaceFrozen, freeze_source_faces,
                           cnf*sizeof(uint8_t));
            }
            if (nsrc + nreq > src_cap) {
                size_t ncap = src_cap ? src_cap : 1024;
                while (ncap < nsrc + nreq) ncap <<= 1;
                int32_t *ns0 = (int32_t *)ARENA_ALLOC(arena,
                                  (ncap*sizeof(int32_t)));
                int32_t *ns1 = (int32_t *)ARENA_ALLOC(arena,
                                  (ncap*sizeof(int32_t)));
                if (src0) memcpy(ns0, src0, nsrc*sizeof(int32_t));
                if (src1) memcpy(ns1, src1, nsrc*sizeof(int32_t));
                src0 = ns0; src1 = ns1; src_cap = ncap;
            }

            size_t wf = cnf;
            for (size_t r = 0; r < nreq; r++) {
                int32_t a = req[r].a, b = req[r].b;
                int32_t fc = req[r].fc, fd = req[r].fd;
                int32_t m = (int32_t)(cnv + r);
                int k;
                /*
                 * Accumulate the midpoint in double.  At PHerc world
                 * coordinates (~4500), adding two floats in float precision
                 * loses a bit before the divide; repeated subdivision can then
                 * drift a child plane by one ULP and turn a sub-voxel-clear
                 * neighbouring chart into a real stab.  The result is still a
                 * float OBJ vertex, but it receives the correctly-rounded chord
                 * midpoint in one rounding step.
                 */
                nV[(size_t)m*3+0] = (float)(0.5*
                    ((double)nV[(size_t)a*3+0]+
                     (double)nV[(size_t)b*3+0]));
                nV[(size_t)m*3+1] = (float)(0.5*
                    ((double)nV[(size_t)a*3+1]+
                     (double)nV[(size_t)b*3+1]));
                nV[(size_t)m*3+2] = (float)(0.5*
                    ((double)nV[(size_t)a*3+2]+
                     (double)nV[(size_t)b*3+2]));
                src0[nsrc + r] = a;
                src1[nsrc + r] = b;
                /* fc traverses a->b: copy a->m (b-side sub-tri), in place b->m */
                {
                    int32_t Fc2[3];
                    for (k=0;k<3;k++) Fc2[k] = (nF[(size_t)fc*3+(size_t)k]==a)
                                                   ? m : nF[(size_t)fc*3+(size_t)k];
                    for (k=0;k<3;k++) if (nF[(size_t)fc*3+(size_t)k]==b)
                                          nF[(size_t)fc*3+(size_t)k] = m;
                    for (k=0;k<3;k++) nF[wf*3+(size_t)k] = Fc2[k];
                    if (nFaceSource)
                        nFaceSource[wf] = nFaceSource[(size_t)fc];
                    if (nFaceFrozen)
                        nFaceFrozen[wf] = nFaceFrozen[(size_t)fc];
                    wf++;
                }
                if (fd >= 0) {
                    /* fd traverses b->a: copy b->m (a-side), in place a->m */
                    int32_t Fd2[3];
                    for (k=0;k<3;k++) Fd2[k] = (nF[(size_t)fd*3+(size_t)k]==b)
                                                   ? m : nF[(size_t)fd*3+(size_t)k];
                    for (k=0;k<3;k++) if (nF[(size_t)fd*3+(size_t)k]==a)
                                          nF[(size_t)fd*3+(size_t)k] = m;
                    for (k=0;k<3;k++) nF[wf*3+(size_t)k] = Fd2[k];
                    if (nFaceSource)
                        nFaceSource[wf] = nFaceSource[(size_t)fd];
                    if (nFaceFrozen)
                        nFaceFrozen[wf] = nFaceFrozen[(size_t)fd];
                    wf++;
                }
            }
            assert(wf == new_nf);
            V = nV; F = nF;
            face_source = nFaceSource;
            face_frozen = nFaceFrozen;
            cnv = new_nv; cnf = new_nf;
            nsrc += nreq;
        }
        free(he); free(fdone); free(req);

        if (st) {
            st->rounds++;
            st->bnd_splits += bnd_this;
            st->int_splits += int_this;
        }

        /* Flip relief so the next round splits well-shaped triangles (and the
         * final geometry is near-Delaunay, circumradius ~ edge/sqrt(3)).
         * This is a seam-band operation: the previous global pass repeatedly
         * sorted every edge in the scroll and could alter remote chart
         * geometry.  Recompute the exact vertex predicate after appending the
         * midpoints and expose only that band to the guarded flip machinery. */
        if (p.flip_max_rounds > 0) {
            uint8_t *flip_active =
                (uint8_t *)malloc(cnv ? cnv : 1);
            int flip_nt = seam_work_threads(cnv, 65536u);
            size_t flips;
            if (flip_active == NULL) return -1;
#ifdef _OPENMP
#pragma omp parallel num_threads(flip_nt)
#endif
            {
                int tid = 0;
#ifdef _OPENMP
                tid = omp_get_thread_num();
#endif
                size_t lo = cnv * (size_t)tid / (size_t)flip_nt;
                size_t hi = cnv * (size_t)(tid + 1) / (size_t)flip_nt;
                size_t v;
                for (v = lo; v < hi; v++)
                    flip_active[v] = (uint8_t)
                        (SeamPlanes_vert_dist(
                            V, (int32_t)v, planes, np) <= (double)p.band);
            }
            flips = WeldCleanup_flip_rounds_active_masked(
                arena, V, cnv, F, cnf, p.flip_max_rounds,
                face_frozen, flip_active);
            free(flip_active);
            if (st) st->flips += flips;
        }
    }

    if (V == NULL) return 0;       /* nothing split anywhere: no-op outputs stand */

    *out_verts = V;   *out_nv = cnv;
    *out_faces = F;   *out_nf = cnf;
    *out_new_vert_parent0 = src0;
    *out_new_vert_parent1 = src1;
    if (out_face_source) *out_face_source = face_source;
    *out_n_new = nsrc;
    if (st) {
        st->verts_added = cnv - nv;
        st->faces_added = cnf - nf;
    }
    return 0;
}

int SeamRefine_process_masked_with_roots(
                       Arena_T arena,
                       const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       const SeamPlane *planes, size_t np,
                       const SeamRefineParams *params,
                       const uint8_t *freeze_source_faces,
                       size_t freeze_source_nf,
                       float **out_verts, size_t *out_nv,
                       int32_t **out_faces, size_t *out_nf,
                       int32_t **out_new_vert_parent0,
                       int32_t **out_new_vert_parent1,
                       int32_t **out_face_source,
                       size_t *out_n_new,
                       SeamRefineStats *st)
{
    return seam_refine_process_impl(
        arena, verts, nv, faces, nf, planes, np, params,
        freeze_source_faces, freeze_source_nf,
        out_verts, out_nv, out_faces, out_nf,
        out_new_vert_parent0, out_new_vert_parent1, out_face_source,
        out_n_new, st);
}

int SeamRefine_process_with_parents(
                       Arena_T arena,
                       const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       const SeamPlane *planes, size_t np,
                       const SeamRefineParams *params,
                       float **out_verts, size_t *out_nv,
                       int32_t **out_faces, size_t *out_nf,
                       int32_t **out_new_vert_parent0,
                       int32_t **out_new_vert_parent1,
                       size_t *out_n_new,
                       SeamRefineStats *st)
{
    return seam_refine_process_impl(
        arena, verts, nv, faces, nf, planes, np, params,
        NULL, 0,
        out_verts, out_nv, out_faces, out_nf,
        out_new_vert_parent0, out_new_vert_parent1, NULL, out_n_new, st);
}

int SeamRefine_process(Arena_T arena,
                       const float *verts, size_t nv,
                       const int32_t *faces, size_t nf,
                       const SeamPlane *planes, size_t np,
                       const SeamRefineParams *params,
                       float **out_verts, size_t *out_nv,
                       int32_t **out_faces, size_t *out_nf,
                       int32_t **out_new_vert_src, size_t *out_n_new,
                       SeamRefineStats *st)
{
    int32_t *unused_parent1 = NULL;
    return seam_refine_process_impl(
        arena, verts, nv, faces, nf, planes, np, params,
        NULL, 0,
        out_verts, out_nv, out_faces, out_nf,
        out_new_vert_src, &unused_parent1, NULL, out_n_new, st);
}
