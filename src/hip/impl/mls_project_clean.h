#ifndef MLS_PROJECT_CLEAN_H
#define MLS_PROJECT_CLEAN_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Required drop-in symbol. This clean implementation accepts only an arena
 * returned by mls_clean_arena_create; an upstream arena is not compatible. */
void MLS_project_verts(void* arena,
                       const float* verts, size_t nv,
                       float radius_vox,
                       const float* cell_origin,
                       float* out_verts,
                       float* out_normals);

/* Clean-side primitive with an explicit fixed-support update count. Compatibility
 * callers must use one update per call and rebuild support between LOP passes. */
void mls_clean_project_verts_iterations(void* arena,
                                        const float* verts, size_t nv,
                                        float radius_vox,
                                        const float* cell_origin,
                                        uint32_t iterations,
                                        float* out_verts,
                                        float* out_normals);

typedef struct mls_clean_arena mls_clean_arena;

typedef enum mls_clean_status {
    MLS_CLEAN_STATUS_OK = 0,
    MLS_CLEAN_STATUS_INVALID_ARGUMENT = 1,
    MLS_CLEAN_STATUS_NONFINITE_INPUT = 2,
    MLS_CLEAN_STATUS_ORIGIN_MISMATCH = 3,
    MLS_CLEAN_STATUS_SIZE_OVERFLOW = 4,
    MLS_CLEAN_STATUS_OUT_OF_MEMORY = 5,
    MLS_CLEAN_STATUS_HIP_ERROR = 6,
    MLS_CLEAN_STATUS_DEGENERATE_NEIGHBORHOOD = 7,
    MLS_CLEAN_STATUS_NUMERICAL_FAILURE = 8,
    MLS_CLEAN_STATUS_WRONG_ARENA = 9
} mls_clean_status;

typedef struct mls_clean_arena_info {
    size_t struct_size;
    size_t point_count;
    size_t query_capacity;
    size_t last_problem_query_count;
    float cell_size_vox;
    float index_origin[3];
    int device_id;
    mls_clean_status last_status;
} mls_clean_arena_info;

/* Construct a clean-side arena from point_count packed x,y,z FP32 support
 * points. index_origin defines only the uniform-grid coordinate frame.
 * cell_origin passed to MLS_project_verts must numerically match it. */
mls_clean_status mls_clean_arena_create(const float* support_points_xyz,
                                        size_t point_count,
                                        float cell_size_vox,
                                        const float index_origin[3],
                                        mls_clean_arena** out_arena);

void mls_clean_arena_destroy(mls_clean_arena* arena);

/* Status for the most recent clean-side call on the current host thread. */
mls_clean_status mls_clean_last_status(void);

/* Snapshot clean-side state without altering the required void ABI. */
mls_clean_status mls_clean_arena_get_info(const mls_clean_arena* arena,
                                          mls_clean_arena_info* out_info);

const char* mls_clean_status_string(mls_clean_status status);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MLS_PROJECT_CLEAN_H */
