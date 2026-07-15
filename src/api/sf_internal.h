/*
 * sf_internal.h — shared internals of the public-API implementation
 * (src/api/*.c only; never installed, never included by modules).
 */
#ifndef SF_INTERNAL_H_INCLUDED
#define SF_INTERNAL_H_INCLUDED

#include "scrollfiesta.h"

#include "../common/arena.h"
#include "../common/except.h"
#include "../common/mesh_types.h"
#include "../common/run_ctx.h"
#include "../common/ves_platform.h"

/* One ScrollFiesta operation at a time per process (file-scope module config
 * globals + the Triangle crash guard are not concurrency-safe yet). */
extern ves_mutex_t g_sf_op_mutex;

/*
 * Per-operation guard. Usage pattern in every wrapper:
 *
 *     SfOpGuard g;
 *     volatile sf_status rc = sf_op_begin(&g, &cfg->common);
 *     if (rc != SF_OK) return rc;
 *     TRY {
 *         ... arena work; module calls; copy-out LAST (copy-out never raises)
 *     }
 *     EXCEPT(Arena_Failed) rc = SF_ERROR_OOM;
 *     EXCEPT(IO_Failed)    rc = SF_ERROR_IO;
 *     EXCEPT(Timeout)      rc = SF_ERROR_TIMEOUT;
 *     EXCEPT(Sf_Cancelled) rc = g.ctx.timed_out ? SF_ERROR_TIMEOUT : SF_CANCELLED;
 *     END_TRY;
 *     if (rc != SF_OK) { free any partially-filled outputs }
 *     sf_op_end(&g);
 *     return rc;
 *
 * No longjmp ever escapes the wrapper: the TRY frame lives here, and every
 * internal RAISE unwinds to it; the arena makes the unwind leak-free.
 */
typedef struct SfOpGuard {
    Arena_T  arena;
    SfRunCtx ctx;
    int      locked;
} SfOpGuard;

sf_status sf_op_begin(SfOpGuard *g, const sf_common_opts *opts);
void      sf_op_end(SfOpGuard *g);

/* Shared tail of every wrapper's EXCEPT chain: map the exception that
 * unwound to the boundary into a status. Must appear between the TRY body's
 * closing brace and END_TRY. */
#define SF_EXCEPT_TAIL(g_, rc_)                                              \
    EXCEPT(Arena_Failed) (rc_) = SF_ERROR_OOM;                               \
    EXCEPT(IO_Failed)    (rc_) = SF_ERROR_IO;                                \
    EXCEPT(Timeout)      (rc_) = SF_ERROR_TIMEOUT;                           \
    EXCEPT(Sf_Cancelled) (rc_) = (g_).ctx.timed_out ? SF_ERROR_TIMEOUT       \
                                                    : SF_CANCELLED;

/* ── conversions (public xyz <-> internal zyx) ─────────────────────────────
 * The axis permutation has determinant -1, so faces swap indices [1]<->[2]
 * in the SAME step; applying both directions round-trips to the identity and
 * keeps winding/normal agreement intact on the inside. All conversion
 * allocations are arena-backed (may RAISE Arena_Failed). */
float   *sf_xyz_to_zyx(Arena_T a, const float *xyz, size_t n);
int32_t *sf_faces_swap(Arena_T a, const int32_t *faces, size_t nf);

/* Deep-copy an sf_mesh into an internal ComponentMesh (arena-backed, zyx,
 * pca_normal/centroid computed, self sentinel set). Returns 0 on success. */
int sf_cm_from_mesh(Arena_T a, const sf_mesh *in, ComponentMesh *cm);

/* ── copy-out (malloc, xyz) ────────────────────────────────────────────────
 * These never RAISE; they return SF_ERROR_OOM and free their own partial
 * work, so they are safe as the final step inside a TRY body. */
sf_status sf_mesh_out(sf_mesh *out,
                      const float *verts_zyx, size_t nv,
                      const int32_t *faces_internal, size_t nf,
                      const float *normals_zyx,
                      const uint8_t *pins,
                      const int32_t *vmap);
sf_status sf_mesh_out_cm(sf_mesh *out, const ComponentMesh *cm, const int32_t *vmap);
/* Raw variant: buffers copied verbatim (already public-convention). */
sf_status sf_mesh_out_raw(sf_mesh *out,
                          const float *verts_xyz, size_t nv,
                          const int32_t *faces, size_t nf,
                          const float *normals_xyz,
                          const uint8_t *pins,
                          const int32_t *vmap);

/* ── provenance ────────────────────────────────────────────────────────────
 * Exact (bit-identical) position join: for each vertex of `out_verts`, the
 * index of a matching vertex in `in_verts`, or -1. Works because every
 * ScrollFiesta split/cleanup op copies surviving positions verbatim.
 * Duplicate positions are consumed in ascending input order; extra claims
 * re-use the first match. Both arrays must be in the SAME axis convention. */
int32_t *sf_vmap_by_position(Arena_T a, const float *in_verts, size_t in_nv,
                             const float *out_verts, size_t out_nv);
int32_t *sf_vmap_identity(Arena_T a, size_t nv);

/* ── misc mesh utilities ───────────────────────────────────────────────────*/
long sf_count_components(Arena_T a, size_t nv, const int32_t *faces, size_t nf);
int  sf_concat_cms(Arena_T a, const ComponentMesh *ms, size_t n, ComponentMesh *out);

#endif /* SF_INTERNAL_H_INCLUDED */
