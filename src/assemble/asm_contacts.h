#ifndef ASM_CONTACTS_INCLUDED
#define ASM_CONTACTS_INCLUDED
#include "../common/arena.h"
#include <stddef.h>
#include <stdint.h>

/* Complete triangle-domain query. No raster cells, chart-size filters,
 * same-chart exclusions or pair-count truncation. Refit retains the tree's
 * topology but recomputes every box from the actual trial coordinates. */
typedef struct AsmContacts AsmContacts;
/* Fixed material policy shared by fields, workers and audits. An overlap
 * is tolerated only while every point is within half this width of an
 * original material boundary. Internal triangulation edges do not count. */
double AsmContacts_material_tolerance(void);
const char *AsmContacts_material_policy(void);
/* leaf_pairs = leaf-leaf node pairs the walk visited; complete = 0 when the
 * walk stopped at its budget (then pairs/area are NOT an answer). */
typedef struct AsmContactStats { size_t pairs; double area; size_t leaf_pairs; int complete; } AsmContactStats;
typedef void (*AsmContactVisit)(void *context, size_t a, size_t b, double area);
AsmContacts *AsmContacts_new(Arena_T arena, size_t nv, const int32_t *faces,
                             size_t nf, const double *uv, double tolerance);
/* Repartition after a large placement change, reusing the allocated storage.
 * This changes search order only; all faces and query masks are retained. */
int AsmContacts_rebuild(AsmContacts *query, const double *uv);
/* As AsmContacts_new, with the worker count set BEFORE the first build so
 * that build runs in parallel (the tree is identical to the serial one). */
AsmContacts *AsmContacts_new_threads(Arena_T arena, size_t nv, const int32_t *faces,
                                     size_t nf, const double *uv, double tolerance, int threads);
/* UV aliases of validated source cuts retain their original material
 * adjacency. material_vertex is an idempotent map into [0,nv), or NULL.
 * Collision geometry still uses every actual UV vertex and original face. */
AsmContacts *AsmContacts_new_source(Arena_T arena, size_t nv, const int32_t *faces,
                                    size_t nf, const double *uv, double tolerance, int threads,
                                    const int32_t *material_vertex);
/* Returns 0 with the complete answer, -1 on a non-finite input, or -2 when the
 * walk exceeded its leaf-pair budget (out->complete = 0, out->leaf_pairs = the
 * count at the stop): a field whose leaves overlap that much is a stack or a
 * tree built for other coordinates, and no caller may read pairs/area then. */
int AsmContacts_measure(AsmContacts *query, const double *uv, double *gradient,
                        AsmContactVisit visit, void *context, AsmContactStats *out);
/* Same material contacts/callbacks, plus a continuous full-intersection-area
 * objective and its gradient, accumulated in the same walk. Allowed boundary
 * pairs contribute to this soft objective, never to the contact certificate.
 * Dropping them from the derivative erases the force on fine-mesh boundaries.
 * Optimization uses this API; final acceptance still uses out->pairs. */
int AsmContacts_objective(AsmContacts *query,const double *uv,double *gradient,
    AsmContactVisit visit,void *context,AsmContactStats *out,double *raw_area);
/* Weight only the soft area/gradient; the complete material statistics and
 * visit callback remain unfiltered. Called in serial reduction order, even
 * for tolerated boundary overlap. A weight must be finite and nonnegative.
 * A collecting callback may return zero while observing every raw area. */
typedef double (*AsmContactWeight)(void *context,size_t a,size_t b,double raw_area);
int AsmContacts_weighted_objective(AsmContacts *query,const double *uv,double *gradient,
    AsmContactVisit visit,void *context,AsmContactWeight weight,void *weight_context,
    AsmContactStats *out,double *raw_area);
/* WALK BUDGET (ASM_CONTACT_MAX_LEAF_PAIRS_PER_LEAF x leaves, capped at
 * ASM_CONTACT_MAX_LEAF_PAIRS, set at every rebuild).  The serial walk, the
 * threaded enumeration and the touching walk count the same leaf-pair sequence,
 * so they stop at the same pair whatever the thread count.  An explicit budget
 * overrides the default until the query is destroyed; 0 = unbounded. */
void AsmContacts_leaf_pair_budget(AsmContacts *query, size_t budget);
/* Bring every box up to date for `uv` (the full refit, ignoring any active
 * mask).  0 on success, -1 on a non-finite coordinate. */
int AsmContacts_refit(AsmContacts *query, const double *uv);
/* Refit only leaves containing changed faces and their ancestors. Every face
 * incident to a changed vertex must be marked; all other coordinates must
 * match the last successful refit. The mask is not retained and does not
 * filter subsequent probes. NULL requests a full refit. */
int AsmContacts_refit_changed(AsmContacts *query,const double *uv,const uint8_t *changed_faces);
/* EXTERNAL TRIANGLES against this tree: `tri` holds n_ext triangles of 6
 * doubles (u,v per corner) in the tree's own frame, and every face not marked
 * in exclude_face[nf] (NULL = none) is tested against each of them with the
 * same positive-area predicate as a measure.  This is how a LOCAL field checks
 * its candidate against the placed material that its own field leaves out; the
 * active mask and the frame labels are not consulted, because an external
 * triangle belongs to no chart of this tree.  `visit` receives (context,
 * external index, tree face, area).  0, -1, or -2 over the walk budget. */
int AsmContacts_probe(AsmContacts *query, const double *uv, size_t n_ext, const double *tri,
                      const uint8_t *exclude_face, AsmContactVisit visit, void *context, AsmContactStats *out);
/* Reuse the boxes and borrowed UVs from the last successful refit/probe.
 * The caller must keep those UVs byte-unchanged between queries. This is for
 * repeated candidate checks against an immutable incumbent, not moved fields. */
/* Optional branch pruning for repeated prepared probes. The borrowed mask
 * must stay byte-identical until cleared with NULL or another refit. */
void AsmContacts_prepare_exclusion(AsmContacts *query,const uint8_t *exclude_face);
int AsmContacts_probe_prepared(AsmContacts *query, size_t n_ext, const double *tri,
                               const uint8_t *exclude_face, AsmContactVisit visit, void *context, AsmContactStats *out);
/* Mesh-aware external queries preserve the source patch's original boundary,
 * including edges outside the selected faces. source_face maps each external
 * index to a source face (NULL = identity). Callbacks retain external indices.
 * Only the source boundary index is refitted; a prepared incumbent's UVs and
 * triangle boxes stay unchanged. A query cannot probe itself at different UVs.
 * Raw triangle probes above treat each external triangle as isolated material;
 * use these mesh variants for repair, admission and material audits. */
int AsmContacts_probe_mesh(AsmContacts *query,const double *uv,AsmContacts *source,const double *source_uv,
    size_t n_ext,const size_t *source_face,const uint8_t *exclude_face,AsmContactVisit visit,void *context,AsmContactStats *out);
int AsmContacts_probe_mesh_prepared(AsmContacts *query,AsmContacts *source,const double *source_uv,
    size_t n_ext,const size_t *source_face,const uint8_t *exclude_face,AsmContactVisit visit,void *context,AsmContactStats *out);
size_t AsmContacts_budget(const AsmContacts *query);
size_t AsmContacts_leaves(const AsmContacts *query);
/* Initializer-only complete SAT query, including borders within 1e-10 voxels.
 * This supplies contact regions for temporary chart translations. It does
 * not change the positive-area predicate used by repair and final audits. */
int AsmContacts_touching(AsmContacts *query, const double *uv,
                         AsmContactVisit visit, void *context, AsmContactStats *out);
/* Only for an explicitly labelled collection of independent coordinate frames
 * (e.g. unregistered extras). Never use material IDs to filter a global sheet.
 * Each frame is queried completely; unrelated frame coordinates are not
 * comparable. The caller retains frame_per_face[nf] for the query lifetime. */
void AsmContacts_separate_frames(AsmContacts *query, const int32_t *frame_per_face);
/* Solver-only affected-pair query: retain every pair with at least one moving
 * face, including collisions against the entire fixed field. NULL restores
 * the complete query used by final audits. The caller owns the mask. */
void AsmContacts_active_faces(AsmContacts *query, const uint8_t *moving_faces);
/* Walk memo tallies since the query was built: masked walks measured afresh
 * and masked requests replayed from an identical stored walk (a replay is
 * the fresh walk's exact committed sequence; see asm_contacts.c). */
void AsmContacts_memo_counts(const AsmContacts *query, size_t *measured, size_t *replayed);

/* Measure the field's contacts on `threads` workers.  <= 1 keeps the exact
 * serial walk.  Enumeration and the reduction stay serial and in the walk's
 * own order, so the answer does not depend on the thread count. */
void AsmContacts_threads(AsmContacts *query, int threads);
/* Independent polygon clipping for audit; boundary integration for gradient.
 * Coordinates are recentered before area products. */
double AsmContacts_pair(const double a[6], const double b[6], double tolerance,
                        double ga[6], double gb[6]);
/* Whole-material predicate for one original face pair, ignoring query masks.
 * UVs matching the last refit use the boundary index and must be unchanged;
 * another array is scanned without mutating the index (before/trial guards).
 * Returns NAN on incomplete geometry/work. The low-level pair above retains
 * numerical triangle SAT semantics and cannot certify material width. */
double AsmContacts_pair_faces(const AsmContacts *query,const double *uv,size_t a,size_t b,double ga[6],double gb[6]);
int AsmContacts_selftest(void);
#endif
