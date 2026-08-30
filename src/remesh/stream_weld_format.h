#ifndef STREAM_WELD_FORMAT_INCLUDED
#define STREAM_WELD_FORMAT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "chart_zipper.h"

/* Binary, independently readable full-resolution weld shards.  A source shard
 * owns vertices/faces for exactly one cube.  A patch shard owns only zipper
 * triangles and addresses vertices in two source shards by their local index.
 * No global vertex numbering is serialized, so adding/removing another cube
 * never rewrites an unchanged artifact. */

typedef struct {
    char cube_id[48];
    size_t nv, nf, n_charts;
    float *verts;                 /* [nv*3] */
    int32_t *faces;               /* [nf*3] */
    int32_t *vertex_chart;        /* [nv]   */
    int32_t *parent0, *parent1;   /* [nv], GWLIN2 semantics */
} StreamWeldShard;

enum {
    STREAM_WELD_CONFLICT_GEOMETRY = 1u
};

/* Pairwise incompatibility between two transactions in one patch shard.
 * Shared source-boundary ports are reconstructed globally by the manifest
 * planner; this record carries exact geometry conflicts found by the bounded
 * two-cube worker. */
typedef struct {
    uint32_t transaction_a;
    uint32_t transaction_b;
    uint32_t reason_mask;
} StreamWeldConflict;

int StreamWeld_write_shard(const char *path, const char *cube_id,
                           const float *verts, size_t nv,
                           const int32_t *faces, size_t nf,
                           const int32_t *vertex_chart, size_t n_charts,
                           const int32_t *parent0, const int32_t *parent1);

int StreamWeld_read_shard(Arena_T arena, const char *path,
                          StreamWeldShard *out);

/* faces contains the complete two-shard source prefix followed by the zipper
 * strips. Each transaction's face_first/face_count names its strip and its
 * orient_a/orient_b pair records the source-chart winding constraint. */
int StreamWeld_write_patch(const char *path,
                           const StreamWeldShard *a,
                           const StreamWeldShard *b,
                           const int32_t *faces, size_t nf,
                           const ChartZipperTransaction *tx, size_t ntx,
                           const StreamWeldConflict *conflicts,
                           size_t nconflicts);

int StreamWeldFormat_selftest(void);

#endif
