#ifndef MESH_CHART_INDEX_INCLUDED
#define MESH_CHART_INDEX_INCLUDED

/* mesh_chart_index -- page the charts of an assembled mesh in, instead of reading all of it.
 *
 * The pipeline's assembled meshes (stage6_repaired.vmesh + stage6_repaired_chart.i32, a close-up's _surface.vmesh +
 * _surface_chart.i32) store each chart's vertices contiguously, and its faces contiguously, referencing only its own
 * vertices.  A sidecar <vmesh>.chartidx records per chart run its vertex and face ranges and its 3-D and UV bounding
 * boxes, so a reader that needs the faces near a box (a close-up: 4.7M of a 21^3 sheet's 700M) seeks to the charts
 * whose boxes meet it -- conservative boxes, so the reader's own exact face tests give the same answer as on the
 * whole mesh (2026-10-13, user: "can we chunk the mesh somehow or make it so we can efficiently page in?").
 *
 * A record is a maximal run of vertices with one chart id.  A chart id >= 0 owns exactly one run; runs of -1
 * (unassigned vertices) may repeat.  A mesh that breaks either property, or a face that spans two runs, is not
 * indexable (MeshChartIndex_build returns 1) and its readers keep the whole-mesh read.
 *
 * The sidecar carries the mesh's header numbers (nv, nf, flags, size) and the chart file's size; a sidecar that no
 * longer matches its mesh is refused by MeshChartIndex_load (stale), never trusted.
 */

#include <stddef.h>
#include <stdint.h>

#include "mesh_bin.h"

typedef struct MeshChartRec {
    int32_t chart;                  /* the run's chart id (-1 = unassigned) */
    uint64_t v0, v1;                /* vertex range [v0, v1) */
    uint64_t f0, f1;                /* face range [f0, f1); f0 == f1 when the run has no faces */
    double lo[3], hi[3];            /* 3-D bounding box of the run's vertices (z, y, x as stored) */
    double uvlo[2], uvhi[2];        /* UV bounding box; lo > hi when the mesh has no UV */
} MeshChartRec;

typedef struct MeshChartIndex {
    MeshChartRec *rec;              /* [n], in vertex order */
    size_t n;
    MeshBinHeader mesh;             /* the indexed mesh's header */
    uint64_t chart_bytes;           /* the chart file's size when indexed */
} MeshChartIndex;

/* Build <out_idx> for <vmesh> with its per-vertex chart file: 0 written, 1 the mesh is not chart-contiguous (nothing
 * written), -1 error.  Streams the mesh in blocks: memory is the records plus one block. */
int MeshChartIndex_build(const char *vmesh, const char *chart_i32, const char *out_idx);

/* The identity of a mesh and its chart file (header numbers, size, mtime), taken BEFORE a whole read of them. */
typedef struct MeshChartStamp {
    MeshBinHeader mesh;
    uint64_t chart_bytes;
    int64_t mtime;
} MeshChartStamp;
int MeshChartIndex_stamp(const char *vmesh, const char *chart_i32, MeshChartStamp *out);   /* 0, -1 */

/* The same index from a mesh already read whole (MeshBin_read_precise_malloc) with its RAW per-vertex charts: a reader
 * that had no index writes it while it holds the data, so the next reader pages.  Nothing is written when the files no
 * longer match `stamp` (changed during the read).  0 written, 1 not chart-contiguous, -1 error / changed. */
int MeshChartIndex_build_mem(const MeshChartStamp *stamp, const char *vmesh, const char *chart_i32,
                             const MeshBinData *m, const int32_t *vchart, const char *out_idx);

/* Load <idx> for <vmesh>: 0 ok, -1 absent, unreadable or stale (the mesh changed since it was indexed). */
int MeshChartIndex_load(const char *idx, const char *vmesh, MeshChartIndex *out);
void MeshChartIndex_free(MeshChartIndex *ix);

/* The selected records (indices into ix->rec, ascending) as one compact mesh, in file order: verts, faces renumbered
 * to the compact vertices, uv64 when the mesh has UV (doubles; floats widened exactly).  *vchart (malloc'd, [nv]) gets
 * each compact vertex's chart id. Reads each file section separately, joining adjacent selected ranges without
 * loading omitted charts; float-UV conversion uses a bounded buffer. MeshBin_dispose frees *out. 0 ok, -1 error. */
int MeshChartIndex_read(const MeshChartIndex *ix, const char *vmesh, const size_t *sel, size_t nsel,
                        MeshBinData *out, int32_t **vchart);

/* A per-vertex side file of `es`-byte elements (e.g. _attr.u8) for the same selection, in the same order (malloc'd);
 * the file must hold exactly ix->mesh.nv elements.  0 ok, -1 error. */
int MeshChartIndex_read_side(const MeshChartIndex *ix, const char *path, size_t es, const size_t *sel, size_t nsel,
                             void **out);

/* <vmesh>.chartidx */
int MeshChartIndex_path(const char *vmesh, char *out, size_t cap);

int MeshChartIndex_selftest(void);

#endif
