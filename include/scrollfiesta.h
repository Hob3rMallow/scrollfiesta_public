/*
 * scrollfiesta.h — the public C API of ScrollFiesta.
 *
 * This is the ONLY header consumers include. It exposes ScrollFiesta's mesh
 * operations (cleanup, detangling/splitting, hole filling, decimation,
 * developability fairing, topology surgery, BPA reconstruction, and the full
 * per-cube volume->mesh pipeline) over plain C types.
 *
 * ── Binding models ─────────────────────────────────────────────────────────
 * 1. Direct linking: link the `scrollfiesta` library (static or shared) and
 *    call the sf_* functions declared below.
 * 2. Runtime loading (the canonical embedding, e.g. VC3D): dlopen()/
 *    LoadLibrary() the shared library, resolve the single symbol
 *    "sf_get_api", and call everything through the returned function-pointer
 *    table:
 *
 *        sf_get_api_fn get_api = (sf_get_api_fn)dlsym(h, "sf_get_api");
 *        const sf_api *sf = get_api(SCROLLFIESTA_ABI_VERSION);
 *        if (!sf) { incompatible library; }
 *        ... sf->topology_audit(...); sf->mesh_free(...);
 *
 *    The table is append-only within an ABI version: fields are never
 *    reordered or removed, so a newer library keeps working with an older
 *    header. Any breaking change bumps SCROLLFIESTA_ABI_VERSION and
 *    sf_get_api returns NULL for mismatched requests. A table entry may be
 *    NULL when the library was built without that capability — check before
 *    calling entries documented as optional.
 *
 * ── Conventions ────────────────────────────────────────────────────────────
 * - Coordinates at this boundary are (x, y, z) in voxel units. (Internally
 *   ScrollFiesta is (z,y,x); the conversion, including the face-index swap
 *   that keeps winding/normal semantics intact under the axis permutation,
 *   happens inside the library.) NB: OBJ files written by the ScrollFiesta
 *   CLI tools use the internal "v z y x" column order — those are NOT
 *   interchangeable with the true-xyz OBJs read/written by this API.
 * - Ownership: structs FILLED BY the library own their buffers; release them
 *   with sf_mesh_free / sf_mesh_list_free / sf_handle_loops_free / sf_free.
 *   Structs filled by the CALLER (op inputs) are never modified or retained.
 *   All library allocations cross the boundary via sf_malloc/sf_free, so a
 *   library built with a different compiler/CRT than the host is safe.
 * - Errors: every operation returns sf_status; no exception, signal, or
 *   longjmp ever crosses this boundary.
 * - Threads: any single operation may be called from any thread, but at most
 *   ONE ScrollFiesta operation runs at a time per process (an internal mutex
 *   serializes them); operations parallelize internally per
 *   sf_common_opts.n_threads.
 * - Progress/cancel: pass a callback in sf_common_opts; return nonzero from
 *   it to request cancellation (the call returns SF_CANCELLED, outputs
 *   untouched).
 */
#ifndef SCROLLFIESTA_H_INCLUDED
#define SCROLLFIESTA_H_INCLUDED

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ══════════════════════════════════════════════════════════════════════════
 * Version / export
 * ════════════════════════════════════════════════════════════════════════ */

#define SCROLLFIESTA_VERSION_MAJOR 0
#define SCROLLFIESTA_VERSION_MINOR 10
#define SCROLLFIESTA_VERSION_PATCH 0

/* Bumped on ANY breaking change to a struct, enum, or signature in this
 * header. sf_get_api(N) returns NULL unless the library implements ABI N. */
#define SCROLLFIESTA_ABI_VERSION   1u

#if defined(_WIN32) && defined(SCROLLFIESTA_SHARED)
  #if defined(SCROLLFIESTA_BUILD)
    #define SF_API __declspec(dllexport)
  #else
    #define SF_API __declspec(dllimport)
  #endif
#else
  #define SF_API
#endif

/* ══════════════════════════════════════════════════════════════════════════
 * Status
 * ════════════════════════════════════════════════════════════════════════ */

typedef enum sf_status {
    SF_OK              = 0,
    SF_ERROR           = -1,  /* operation-specific internal failure        */
    SF_ERROR_BAD_ARG   = -2,  /* NULL/inconsistent input, failed validation */
    SF_ERROR_OOM       = -3,  /* out of memory                              */
    SF_ERROR_IO        = -4,  /* file I/O failure                           */
    SF_ERROR_TIMEOUT   = -5,  /* sf_common_opts.timeout_sec exceeded        */
    SF_CANCELLED       = -6,  /* progress callback requested cancellation   */
    SF_ERROR_NO_TIFF   = -7,  /* built without TIFF support                 */
    SF_ERROR_ABI       = -8,  /* header/library ABI mismatch                */
    SF_ERROR_UNSUPPORTED = -9 /* capability not available in this build     */
} sf_status;

/* ══════════════════════════════════════════════════════════════════════════
 * Mesh
 * ════════════════════════════════════════════════════════════════════════ */

/* Indexed triangle mesh. Vertices/normals are (x,y,z) voxel-space float32.
 * Faces are 0-based index triples; winding is consistent with the normals
 * (counter-clockwise seen from the normal side).
 *
 * vmap (library-filled outputs only): per-vertex provenance back to the
 * operation's INPUT mesh — vmap[i] is the input-vertex index that output
 * vertex i continues, or -1 for newly created vertices (e.g. hole-fill
 * Steiner points). NULL when the operation cannot provide provenance
 * (decimation, BPA/pipeline reconstruction). On input meshes vmap is
 * ignored. */
typedef struct sf_mesh {
    float   *vertices;        /* [n_vertices*3]                     */
    int32_t *faces;           /* [n_faces*3]                        */
    float   *vertex_normals;  /* [n_vertices*3] or NULL             */
    uint8_t *pin_mask;        /* [n_vertices], 1 = do not move; or NULL */
    int32_t *vmap;            /* [n_vertices] or NULL (outputs only) */
    size_t   n_vertices;
    size_t   n_faces;
} sf_mesh;

typedef struct sf_mesh_list {
    sf_mesh *items;
    size_t   count;
} sf_mesh_list;

/* ══════════════════════════════════════════════════════════════════════════
 * Common options (embedded as the first field of every op config)
 * ════════════════════════════════════════════════════════════════════════ */

/* Return 0 to continue, nonzero to request cancellation. `stage` is a short
 * stable identifier ("cleanup", "bpa", "holefill", ...); `fraction` is a
 * best-effort monotone 0..1 within the operation. Called on the operation's
 * calling thread. Must be fast and must not call back into ScrollFiesta. */
typedef int  (*sf_progress_fn)(void *user, const char *stage, double fraction);
typedef void (*sf_log_fn)(void *user, const char *line);

typedef struct sf_kv { const char *key; const char *value; } sf_kv;

typedef struct sf_common_opts {
    int            n_threads;     /* <=0 = library default (cpu count)      */
    double         timeout_sec;   /* <=0 = none; SF_ERROR_TIMEOUT on expiry */
    sf_progress_fn progress;      /* NULL = no progress, not cancellable    */
    void          *progress_user;
    sf_log_fn      log;           /* NULL = library logs to stderr          */
    void          *log_user;
    /* Expert tuning knobs by their documented env-var names, applied for
     * this call only (e.g. {"BPA_WALL_GUARD_COS","0.5"}). */
    const sf_kv   *tuning;
    size_t         n_tuning;
} sf_common_opts;

/* ══════════════════════════════════════════════════════════════════════════
 * Diagnostics
 * ════════════════════════════════════════════════════════════════════════ */

typedef struct sf_topology_report {
    long   n_verts, n_faces, n_edges;
    long   n_boundary_edges;      /* incident to exactly 1 triangle        */
    long   n_manifold_edges;      /* exactly 2                             */
    long   n_nonmanifold_edges;   /* 3 or more                             */
    long   n_boundary_loops;
    long   euler;                 /* V - E + F                             */
    double genus;                 /* (2 - loops - euler)/2 when manifold   */
    long   n_unref_verts;
    int    is_disk;               /* manifold, 1 boundary loop, genus 0    */
    long   n_components;          /* vertex-connectivity components        */
} sf_topology_report;

typedef struct sf_developability_report {
    double nondev_fraction;  /* interior verts with lambda > eps (0..1) */
    double mean_lambda;      /* mean interior Crane energy              */
} sf_developability_report;

/* ══════════════════════════════════════════════════════════════════════════
 * Cleanup family
 * ════════════════════════════════════════════════════════════════════════ */

typedef struct sf_cleanup_config {
    sf_common_opts common;
    int    manifold_repair;     /* resolve >2-face edges + split pinch verts (default 1) */
    int    reorient;            /* re-assert consistent winding after repair (default 1) */
    int    fill_pinholes;       /* close exact single-triangle holes (default 1)         */
    int    sliver_cleanup;      /* edge-flip + guarded collapse pass (default 0)         */
    double sliver_min_alt;      /* sliver threshold, vox (default 0.10)                  */
    double degen_area;          /* degenerate-face area, vox^2 (default 1e-3)            */
    double max_collapse_len;    /* longest collapsible edge, vox (default 5.0)           */
    float  cull_min_area_frac;  /* drop components under this fraction of total area;
                                 * <=0 = no cull (default 0; pipeline uses 0.02)         */
    int    respect_pins;        /* honor in->pin_mask (default 1)                        */
} sf_cleanup_config;

typedef struct sf_cleanup_report {
    size_t nm_edges_resolved, faces_deleted, pinch_splits, orient_flips;
    size_t pinholes_filled, tris_added, loops_skipped;
    size_t n_flips, n_collapses;
    size_t components_culled;
} sf_cleanup_report;

typedef struct sf_holefill_config {
    sf_common_opts common;
    int interior_only;   /* 1 = fill only geometrically-interior holes (default 1) */
} sf_holefill_config;

typedef struct sf_holefill_report {
    size_t n_loops, n_interior, n_filled;
} sf_holefill_report;

typedef struct sf_decimate_config {
    sf_common_opts common;
    float  target_ratio;   /* keep this fraction of faces (default 0.075)   */
    size_t target_faces;   /* nonzero overrides target_ratio                */
    int    respect_pins;   /* honor in->pin_mask (default 1)                */
} sf_decimate_config;

typedef struct sf_decimate_report {
    size_t faces_in, faces_out, verts_out;
} sf_decimate_report;

typedef struct sf_fair_config {
    sf_common_opts common;   /* Stein/Grinspun/Crane developability descent  */
    int    max_iters;        /* <=0 -> 100                                   */
    double tol_grad;         /* <=0 -> 1e-6                                  */
    double tol_e;            /* <=0 -> 1e-6                                  */
    double max_disp_vox;     /* per-vertex displacement cap (<=0 -> 1.0)     */
    int    pin_boundary;     /* hold the boundary fixed (default 1)          */
} sf_fair_config;

typedef struct sf_fair_report {
    double e_initial, e_final, max_disp;
    int    iters, converged;
} sf_fair_report;

typedef struct sf_orient_config {
    sf_common_opts common;
    int   align_components;  /* also align component signs spatially (default 1) */
    float align_radius;      /* neighbour ball for the cross-component vote (vox, default 3.0) */
} sf_orient_config;

typedef struct sf_orient_report {
    size_t flips, components, residual_same_dir, components_flipped;
} sf_orient_report;

/* ══════════════════════════════════════════════════════════════════════════
 * Detangle / split family (multi-mesh outputs)
 * ════════════════════════════════════════════════════════════════════════ */

typedef struct sf_split_report {
    size_t n_pieces;          /* count == 1 means self-gated: no split      */
    int    self_gated;
    double nondev_before;     /* developability split only, else 0          */
    double nondev_after_max;
} sf_split_report;

typedef struct sf_split_peel_config {
    sf_common_opts common;
    double min_gap;           /* depth jump across an edge, vox (default 1.5) */
    double gap_depth;         /* seam exclusion radius, vox (default 1.0)     */
    size_t min_comp_verts;    /* smallest surviving layer (default 200)       */
} sf_split_peel_config;

typedef struct sf_split_dev_config {
    sf_common_opts common;
    double eps;               /* Crane-energy seam threshold (default 0.05)   */
    double gap_depth;         /* seam exclusion radius, vox (default 3.5)     */
    size_t min_comp_verts;    /* smallest surviving piece (default 200)       */
    int    apply_gate;        /* accept only if pieces are measurably more
                               * developable than the parent (default 1)      */
    double gate_margin;       /* required improvement (default 0.02)          */
    double gate_abs_ceil;     /* absolute piece ceiling (default 0.05)        */
} sf_split_dev_config;

typedef struct sf_split_bridge_config {
    sf_common_opts common;
    int    oracle_grid_size;      /* sheet-count raycast grid (default 320)   */
    double component_timeout_sec; /* per-component wall clock (default 30;
                                   * <=0 = unlimited)                         */
} sf_split_bridge_config;

typedef struct sf_split_overlap_config {
    sf_common_opts common;
    int    sheet_count;           /* expected sheets; 0 = auto (default 0)    */
    double component_timeout_sec; /* per-component wall clock (default 30)    */
} sf_split_overlap_config;

/* Umbrella detangle: depth-peel -> developability cut (gated) -> bridge cut
 * -> overlap separation, applied per piece. Mirrors the per-cube pipeline's
 * split cascade EXCEPT that pieces are not re-projected from their source
 * point clouds (a mesh-only input has none), so the developability gate is
 * evaluated on the raw pieces. */
typedef struct sf_detangle_config {
    sf_common_opts common;
    int    enable_peel, enable_dev, enable_bridge, enable_overlap;
    size_t min_comp_verts;        /* default 200 */
    double component_timeout_sec; /* per-piece wall clock for the bridge/
                                   * overlap stages (default 30)              */
} sf_detangle_config;

typedef struct sf_detangle_report {
    size_t n_pieces;
    size_t peel_splits, dev_splits, bridge_splits, overlap_splits;
} sf_detangle_report;

/* ══════════════════════════════════════════════════════════════════════════
 * Topology surgery
 * ════════════════════════════════════════════════════════════════════════ */

typedef struct sf_seamcut_config {
    sf_common_opts common;    /* Seamster (Sheffer & Hart 2002)              */
    double accept_frac;       /* fraction of distortion left uncut (default 0) */
    int    r_max;             /* max region-radius level (default 5)          */
} sf_seamcut_config;

typedef struct sf_seamcut_report {
    long genus_before, genus_after;
    long loops_before, loops_after;
    long n_terminals, seam_edges;
    int  is_disk_after;
    long n_components;        /* components cut (each becomes a disk)        */
} sf_seamcut_report;

typedef struct sf_handle_loop {
    int32_t *vertex_indices;  /* indices into the INPUT mesh                 */
    size_t   n;
    double   length;          /* summed edge length (vox)                    */
} sf_handle_loop;

/* ══════════════════════════════════════════════════════════════════════════
 * Reconstruction
 * ════════════════════════════════════════════════════════════════════════ */

typedef struct sf_bpa_config {
    sf_common_opts common;
    float rho;                /* pivot ball radius, vox (<=0 -> 1.2)         */
} sf_bpa_config;

typedef struct sf_bpa_report {
    size_t n_faces;
} sf_bpa_report;

/* ══════════════════════════════════════════════════════════════════════════
 * Volume pipeline
 * ════════════════════════════════════════════════════════════════════════ */

/* In-memory binary prediction volume: nonzero = surface voxel.
 * data[z*(ny*nx) + y*nx + x]; nx==ny==nz==cube_size + 2*halo_voxels. */
typedef struct sf_volume {
    const uint8_t *data;
    int nx, ny, nz;
} sf_volume;

typedef struct sf_pipeline_config {
    sf_common_opts common;
    int     cube_size;          /* owned region edge (default 128)           */
    int     halo_voxels;        /* neighbour padding (default 13)            */
    int64_t origin_xyz[3];      /* owned-cube world origin, (x,y,z)          */
    float   qem_target_ratio;   /* <=0 -> default 0.075                      */
    int     skip_qem;
    int     enable_depth_peel;  /* default 1                                 */
    int     enable_dev_cut;     /* default 1                                 */
    int     enable_sever;       /* short-handle sever (default 1)            */
    const char *dump_dir;       /* per-stage OBJ debug dumps; NULL = none    */
} sf_pipeline_config;

typedef struct sf_pipeline_report {
    double t_extract, t_qem, t_trim, t_dump;
    size_t n_bad_sheets;        /* components flagged by the dev gate        */
} sf_pipeline_report;

/* Gradual multi-winding growth gate (optional).
 *
 * A purely local pairwise winding test cannot see a Ball-Pivoting front that
 * walks gradually up a through-thickness wall: every individual edge looks
 * innocent while one growth front accumulates several complete turns, fusing
 * wraps. Arming this gate integrates the branch-cut-free local phase step
 * dw = dr/pitch - dtheta/(2*pi) from each seed and refuses a pivot once the
 * accumulated phase leaves [-tol,+tol]; re-seeding then reconstructs the next
 * wrap as a distinct component.
 *
 * It has no dedicated config field — like the weld's rho_max it is reached
 * through documented knobs, so pass them in sf_common_opts.tuning (per call,
 * no process environment involved). ALL FOUR are required to arm it; any
 * missing or out-of-range value leaves the historical ungated behaviour:
 *
 *   BPA_GROW_UMBILICUS_Y  scroll axis y, source-space voxels
 *   BPA_GROW_UMBILICUS_X  scroll axis x, source-space voxels
 *   BPA_GROW_WRAP_PITCH   radial spacing per turn, voxels (> 0)
 *   BPA_GROW_WIND_TOL     max accumulated seed phase, turns (0 < tol < 0.5)
 *
 * PHerc0139, for example: {"3405","2878","9.5","0.45"}. The origin the phase
 * is measured against is sf_pipeline_config.origin_xyz, so no extra plumbing
 * is needed. The grid_pipeline CLI arms the same gate from
 * --umb-y/--umb-x/--wrap-pitch/--grow-wind-tol. */

/* ══════════════════════════════════════════════════════════════════════════
 * Weld (experimental — probe the table entry for NULL before calling; a
 * build may ship without it)
 * ════════════════════════════════════════════════════════════════════════ */

typedef struct sf_weld_config {
    sf_common_opts common;
    float cube_size;        /* seam-plane spacing (default 128)              */
    float rho, rho_max;     /* bridge ball radius floor/cap (1.5 / 3.0)      */
    float weld_eps;         /* duplicate-vertex merge tolerance (1e-4)       */
    int   fill_holes;       /* interior hole fill after bridging (default 1) */
    int   cleanup;          /* sliver cleanup + manifold guard (default 1)   */
    float band;             /* seam-edge selection half-width, vox (6.0)     */
} sf_weld_config;

typedef struct sf_weld_report {
    size_t verts_welded, bridge_faces, holes_filled;
} sf_weld_report;

/* ══════════════════════════════════════════════════════════════════════════
 * Functions (direct-link surface; the sf_api table mirrors all of these)
 * ════════════════════════════════════════════════════════════════════════ */

SF_API void        sf_version(int *major, int *minor, int *patch);
SF_API const char *sf_version_string(void);
SF_API const char *sf_status_str(sf_status s);

SF_API void *sf_malloc(size_t n);
SF_API void  sf_free(void *p);
SF_API void  sf_mesh_free(sf_mesh *m);            /* frees buffers, zeroes struct */
SF_API void  sf_mesh_list_free(sf_mesh_list *l);
SF_API void  sf_handle_loops_free(sf_handle_loop *loops, size_t n);
SF_API sf_status sf_mesh_validate(const sf_mesh *m, char *err, size_t err_len);

/* True-xyz OBJ convenience I/O (NOT the internal zyx tool format). */
SF_API sf_status sf_mesh_load_obj(const char *path, sf_mesh *out);
SF_API sf_status sf_mesh_save_obj(const char *path, const sf_mesh *m);

/* Defaults */
SF_API sf_common_opts          sf_common_opts_default(void);
SF_API sf_cleanup_config       sf_cleanup_config_default(void);
SF_API sf_holefill_config      sf_holefill_config_default(void);
SF_API sf_decimate_config      sf_decimate_config_default(void);
SF_API sf_fair_config          sf_fair_config_default(void);
SF_API sf_orient_config        sf_orient_config_default(void);
SF_API sf_split_peel_config    sf_split_peel_config_default(void);
SF_API sf_split_dev_config     sf_split_dev_config_default(void);
SF_API sf_split_bridge_config  sf_split_bridge_config_default(void);
SF_API sf_split_overlap_config sf_split_overlap_config_default(void);
SF_API sf_detangle_config      sf_detangle_config_default(void);
SF_API sf_seamcut_config       sf_seamcut_config_default(void);
SF_API sf_bpa_config           sf_bpa_config_default(void);
SF_API sf_pipeline_config      sf_pipeline_config_default(void);
SF_API sf_weld_config          sf_weld_config_default(void);

/* Diagnostics */
SF_API sf_status sf_topology_audit(const sf_mesh *in, const sf_common_opts *opts,
                                   sf_topology_report *out);
SF_API sf_status sf_developability_energy(const sf_mesh *in, double eps,
                                          const sf_common_opts *opts,
                                          double **out_lambda /* [nv], sf_free */,
                                          sf_developability_report *rep);
SF_API sf_status sf_count_sheets(const sf_mesh *in, int grid_size /* <=0 -> 320 */,
                                 const sf_common_opts *opts, int *out_sheets);

/* Cleanup family */
SF_API sf_status sf_cleanup(const sf_mesh *in, const sf_cleanup_config *cfg,
                            sf_mesh *out, sf_cleanup_report *rep);
SF_API sf_status sf_fill_holes(const sf_mesh *in, const sf_holefill_config *cfg,
                               sf_mesh *out, sf_holefill_report *rep);
SF_API sf_status sf_decimate(const sf_mesh *in, const sf_decimate_config *cfg,
                             sf_mesh *out, sf_decimate_report *rep);
SF_API sf_status sf_fair_developable(const sf_mesh *in, const sf_fair_config *cfg,
                                     sf_mesh *out, sf_fair_report *rep);
SF_API sf_status sf_orient(const sf_mesh *in, const sf_orient_config *cfg,
                           sf_mesh *out, sf_orient_report *rep);

/* Detangle / split family */
SF_API sf_status sf_split_depth_peel(const sf_mesh *in, const sf_split_peel_config *cfg,
                                     sf_mesh_list *out, sf_split_report *rep);
SF_API sf_status sf_split_developability(const sf_mesh *in, const sf_split_dev_config *cfg,
                                         sf_mesh_list *out, sf_split_report *rep);
SF_API sf_status sf_split_bridges(const sf_mesh *in, const sf_split_bridge_config *cfg,
                                  sf_mesh_list *out, sf_split_report *rep);
SF_API sf_status sf_split_overlap(const sf_mesh *in, const sf_split_overlap_config *cfg,
                                  sf_mesh_list *out, sf_split_report *rep);
SF_API sf_status sf_detangle(const sf_mesh *in, const sf_detangle_config *cfg,
                             sf_mesh_list *out, sf_detangle_report *rep);

/* Topology surgery */
SF_API sf_status sf_cut_to_disk(const sf_mesh *in, const sf_seamcut_config *cfg,
                                sf_mesh *out, sf_seamcut_report *rep);
SF_API sf_status sf_sever_short_handles(const sf_mesh *in,
                                        double max_loop_len /* <=0 -> cut all */,
                                        const sf_common_opts *opts,
                                        sf_mesh *out, long *out_loops_cut);
SF_API sf_status sf_find_handle_loops(const sf_mesh *in, const sf_common_opts *opts,
                                      sf_handle_loop **out_loops, size_t *out_n);

/* Reconstruction */
SF_API sf_status sf_reconstruct_bpa(const float *points_xyz,
                                    const float *normals_xyz /* required */,
                                    size_t n_points, const sf_bpa_config *cfg,
                                    sf_mesh *out, sf_bpa_report *rep);
SF_API sf_status sf_mls_project(const float *points_xyz, size_t n_points,
                                float radius_vox, const sf_common_opts *opts,
                                float **out_points_xyz /* sf_free */,
                                float **out_normals_xyz /* sf_free, may be NULL */);

/* Volume pipeline */
SF_API sf_status sf_pipeline_run(const sf_volume *vol, const sf_pipeline_config *cfg,
                                 sf_mesh_list *out_full,
                                 sf_mesh_list *out_trimmed,
                                 sf_pipeline_report *rep);

/* Weld adjacent per-cube meshes (world coordinates, each consistently wound
 * — pipeline outputs are) into one mesh: duplicate-vertex merge + face dedup
 * across the halo overlap, BPA seam bridging across cube-boundary planes,
 * pinhole + interior hole fill, cross-component orientation, sliver cleanup,
 * manifold guard. Composes the library weld modules; the grid_weld TOOL has
 * additional env-gated experimental stages this call does not run. */
SF_API sf_status sf_weld(const sf_mesh *meshes, size_t n_meshes,
                         const sf_weld_config *cfg,
                         sf_mesh *out, sf_weld_report *rep);

/* Runs the built-in module self-tests; 0 = all pass. */
SF_API int sf_selftest(void);

/* ══════════════════════════════════════════════════════════════════════════
 * The API table — the dlopen surface.
 * Append-only within an ABI version; never reorder or remove fields.
 * ════════════════════════════════════════════════════════════════════════ */

typedef struct sf_api {
    uint32_t abi_version;   /* == SCROLLFIESTA_ABI_VERSION of the library */
    uint32_t struct_size;   /* sizeof(sf_api) the library was built with  */

    void        (*version)(int *major, int *minor, int *patch);
    const char *(*version_string)(void);
    const char *(*status_str)(sf_status s);

    void      *(*malloc_)(size_t n);
    void       (*free_)(void *p);
    void       (*mesh_free)(sf_mesh *m);
    void       (*mesh_list_free)(sf_mesh_list *l);
    void       (*handle_loops_free)(sf_handle_loop *loops, size_t n);
    sf_status  (*mesh_validate)(const sf_mesh *m, char *err, size_t err_len);

    sf_status  (*mesh_load_obj)(const char *path, sf_mesh *out);
    sf_status  (*mesh_save_obj)(const char *path, const sf_mesh *m);

    sf_common_opts          (*common_opts_default)(void);
    sf_cleanup_config       (*cleanup_config_default)(void);
    sf_holefill_config      (*holefill_config_default)(void);
    sf_decimate_config      (*decimate_config_default)(void);
    sf_fair_config          (*fair_config_default)(void);
    sf_orient_config        (*orient_config_default)(void);
    sf_split_peel_config    (*split_peel_config_default)(void);
    sf_split_dev_config     (*split_dev_config_default)(void);
    sf_split_bridge_config  (*split_bridge_config_default)(void);
    sf_split_overlap_config (*split_overlap_config_default)(void);
    sf_detangle_config      (*detangle_config_default)(void);
    sf_seamcut_config       (*seamcut_config_default)(void);
    sf_bpa_config           (*bpa_config_default)(void);
    sf_pipeline_config      (*pipeline_config_default)(void);
    sf_weld_config          (*weld_config_default)(void);

    sf_status  (*topology_audit)(const sf_mesh *in, const sf_common_opts *opts,
                                 sf_topology_report *out);
    sf_status  (*developability_energy)(const sf_mesh *in, double eps,
                                        const sf_common_opts *opts,
                                        double **out_lambda,
                                        sf_developability_report *rep);
    sf_status  (*count_sheets)(const sf_mesh *in, int grid_size,
                               const sf_common_opts *opts, int *out_sheets);

    sf_status  (*cleanup)(const sf_mesh *in, const sf_cleanup_config *cfg,
                          sf_mesh *out, sf_cleanup_report *rep);
    sf_status  (*fill_holes)(const sf_mesh *in, const sf_holefill_config *cfg,
                             sf_mesh *out, sf_holefill_report *rep);
    sf_status  (*decimate)(const sf_mesh *in, const sf_decimate_config *cfg,
                           sf_mesh *out, sf_decimate_report *rep);
    sf_status  (*fair_developable)(const sf_mesh *in, const sf_fair_config *cfg,
                                   sf_mesh *out, sf_fair_report *rep);
    sf_status  (*orient)(const sf_mesh *in, const sf_orient_config *cfg,
                         sf_mesh *out, sf_orient_report *rep);

    sf_status  (*split_depth_peel)(const sf_mesh *in, const sf_split_peel_config *cfg,
                                   sf_mesh_list *out, sf_split_report *rep);
    sf_status  (*split_developability)(const sf_mesh *in, const sf_split_dev_config *cfg,
                                       sf_mesh_list *out, sf_split_report *rep);
    sf_status  (*split_bridges)(const sf_mesh *in, const sf_split_bridge_config *cfg,
                                sf_mesh_list *out, sf_split_report *rep);
    sf_status  (*split_overlap)(const sf_mesh *in, const sf_split_overlap_config *cfg,
                                sf_mesh_list *out, sf_split_report *rep);
    sf_status  (*detangle)(const sf_mesh *in, const sf_detangle_config *cfg,
                           sf_mesh_list *out, sf_detangle_report *rep);

    sf_status  (*cut_to_disk)(const sf_mesh *in, const sf_seamcut_config *cfg,
                              sf_mesh *out, sf_seamcut_report *rep);
    sf_status  (*sever_short_handles)(const sf_mesh *in, double max_loop_len,
                                      const sf_common_opts *opts,
                                      sf_mesh *out, long *out_loops_cut);
    sf_status  (*find_handle_loops)(const sf_mesh *in, const sf_common_opts *opts,
                                    sf_handle_loop **out_loops, size_t *out_n);

    sf_status  (*reconstruct_bpa)(const float *points_xyz, const float *normals_xyz,
                                  size_t n_points, const sf_bpa_config *cfg,
                                  sf_mesh *out, sf_bpa_report *rep);
    sf_status  (*mls_project)(const float *points_xyz, size_t n_points,
                              float radius_vox, const sf_common_opts *opts,
                              float **out_points_xyz, float **out_normals_xyz);

    sf_status  (*pipeline_run)(const sf_volume *vol, const sf_pipeline_config *cfg,
                               sf_mesh_list *out_full, sf_mesh_list *out_trimmed,
                               sf_pipeline_report *rep);

    /* Experimental — NULL when unavailable in this library build/version. */
    sf_status  (*weld)(const sf_mesh *meshes, size_t n_meshes,
                       const sf_weld_config *cfg,
                       sf_mesh *out, sf_weld_report *rep);

    int        (*selftest)(void);
} sf_api;

/* The single dlopen entry point. Returns the library's function table, or
 * NULL if the library does not implement `requested_abi`. The returned
 * pointer is static — valid for the lifetime of the loaded library. */
SF_API const sf_api *sf_get_api(uint32_t requested_abi);
typedef const sf_api *(*sf_get_api_fn)(uint32_t requested_abi);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SCROLLFIESTA_H_INCLUDED */
