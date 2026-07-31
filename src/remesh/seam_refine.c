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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

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

typedef struct { int32_t v0,v1,face,opposite; int8_t fwd; } RHE;
static int rhe_cmp(const void *pa, const void *pb)
{
    const RHE *a=(const RHE*)pa, *b=(const RHE*)pb;
    if (a->v0 != b->v0) return a->v0 < b->v0 ? -1 : 1;
    if (a->v1 != b->v1) return a->v1 < b->v1 ? -1 : 1;
    return 0;
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
                       float **out_verts, size_t *out_nv,
                       int32_t **out_faces, size_t *out_nf,
                       int32_t **out_new_vert_parent0,
                       int32_t **out_new_vert_parent1,
                       size_t *out_n_new,
                       SeamRefineStats *st)
{
    SeamRefineParams p;
    float   *V = NULL;
    int32_t *F = NULL;
    int32_t *src0 = NULL, *src1 = NULL; /* per NEW vert: edge endpoints */
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
    if (nf == 0 || nv == 0 || np == 0 || planes == NULL) return 0;
    if (p.target_len <= 0.0f || p.max_rounds <= 0) return 0;

    for (round = 0; round < p.max_rounds; round++) {
        const float   *cV = V ? V : verts;
        const int32_t *cF = F ? F : faces;
        size_t n_he = cnf * 3, i = 0, nreq = 0;
        size_t bnd_this = 0, int_this = 0;

        RHE *he = (RHE *)malloc((n_he ? n_he : 1) * sizeof(RHE));
        uint8_t *fdone = (uint8_t *)calloc(cnf ? cnf : 1, 1);
        RSplit *req = (RSplit *)malloc((n_he/2 + 2) * sizeof(RSplit));
        if (!he || !fdone || !req) { free(he); free(fdone); free(req); return -1; }

        for (size_t f = 0; f < cnf; f++) {
            int32_t v[3] = { cF[f*3+0], cF[f*3+1], cF[f*3+2] };
            for (int e = 0; e < 3; e++) {
                int32_t a = v[e], b = v[(e+1)%3], opp = v[(e+2)%3];
                size_t idx = f*3 + (size_t)e;
                he[idx].v0 = (a<b)?a:b; he[idx].v1 = (a<b)?b:a;
                he[idx].face = (int32_t)f; he[idx].opposite = opp;
                he[idx].fwd = (int8_t)((a<b)?1:0);
            }
        }
        qsort(he, n_he, sizeof(RHE), rhe_cmp);

        while (i < n_he) {
            size_t j = i + 1;
            while (j < n_he && he[j].v0 == he[i].v0 && he[j].v1 == he[i].v1) j++;
            size_t run = j - i;
            int32_t a = he[i].v0, b = he[i].v1;
            do {
                if (run > 2) break;                       /* non-manifold: leave alone */
                if (edge_len3(cV, a, b) <= (double)p.target_len) break;
                /* band test: either endpoint within band of a seam plane */
                if (SeamPlanes_vert_dist(cV, a, planes, np) > (double)p.band &&
                    SeamPlanes_vert_dist(cV, b, planes, np) > (double)p.band)
                    break;
                if (run == 2) {
                    /* interior edge: needs oppositely-wound faces (same_dir
                     * pairs are left for the orient passes, as split_round) */
                    int32_t fc, fd;
                    if (he[i].fwd == he[i+1].fwd) break;
                    if (he[i].fwd) { fc = he[i].face; fd = he[i+1].face; }
                    else           { fc = he[i+1].face; fd = he[i].face; }
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
                    int32_t fc = he[i].face;
                    int32_t da = he[i].fwd ? a : b;
                    int32_t db = he[i].fwd ? b : a;
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
                              (long)(new_nv*3*sizeof(float)));
            int32_t *nF = (int32_t *)ARENA_ALLOC(arena,
                              (long)(new_nf*3*sizeof(int32_t)));
            memcpy(nV, cV, cnv*3*sizeof(float));
            memcpy(nF, cF, cnf*3*sizeof(int32_t));
            if (nsrc + nreq > src_cap) {
                size_t ncap = src_cap ? src_cap : 1024;
                while (ncap < nsrc + nreq) ncap <<= 1;
                int32_t *ns0 = (int32_t *)ARENA_ALLOC(arena,
                                  (long)(ncap*sizeof(int32_t)));
                int32_t *ns1 = (int32_t *)ARENA_ALLOC(arena,
                                  (long)(ncap*sizeof(int32_t)));
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
                nV[(size_t)m*3+0] = 0.5f*(nV[(size_t)a*3+0]+nV[(size_t)b*3+0]);
                nV[(size_t)m*3+1] = 0.5f*(nV[(size_t)a*3+1]+nV[(size_t)b*3+1]);
                nV[(size_t)m*3+2] = 0.5f*(nV[(size_t)a*3+2]+nV[(size_t)b*3+2]);
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
                    wf++;
                }
            }
            assert(wf == new_nf);
            V = nV; F = nF;
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
         * final geometry is near-Delaunay, circumradius ~ edge/sqrt(3)). */
        if (st) st->flips += WeldCleanup_flip_rounds(arena, V, cnv, F, cnf,
                                                     p.flip_max_rounds);
        else    (void)WeldCleanup_flip_rounds(arena, V, cnv, F, cnf,
                                              p.flip_max_rounds);
    }

    if (V == NULL) return 0;       /* nothing split anywhere: no-op outputs stand */

    *out_verts = V;   *out_nv = cnv;
    *out_faces = F;   *out_nf = cnf;
    *out_new_vert_parent0 = src0;
    *out_new_vert_parent1 = src1;
    *out_n_new = nsrc;
    if (st) {
        st->verts_added = cnv - nv;
        st->faces_added = cnf - nf;
    }
    return 0;
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
        out_verts, out_nv, out_faces, out_nf,
        out_new_vert_parent0, out_new_vert_parent1, out_n_new, st);
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
        out_verts, out_nv, out_faces, out_nf,
        out_new_vert_src, &unused_parent1, out_n_new, st);
}
