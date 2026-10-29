#ifndef ASM_EMIT_INCLUDED
#define ASM_EMIT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "asm_types.h"

/* ============================================================================
 * asm_emit.h -- the deliverable: the mesh itself carrying (u, v).
 *
 * sheet.vmesh (+ sheet.obj companion) holds every chart of the globally
 * registered global field with its source triangles unchanged and per-vertex uv
 * from the chart pose; sidecars carry provenance (_component.i32,
 * _chart.i32, _cube.i32).  sheet_extras.vmesh holds every other in-layout
 * unregistered chart in its own component frame. The repair and complete source
 * audit keep that missing global coverage explicit. The RAW bake is produced by spawning
 * obj_bake_raw (it reads the companion vmesh), so the texture is sampled
 * fresh from the CT at every pixel.
 * ==========================================================================*/

typedef struct AsmEmitOpts {
    const char *raw_source;      /* complete RAW cube dir or zarr; NULL skips the bake */
    const char *bake_exe;        /* obj_bake_raw.exe path */
    double      normal_range;    /* bake: max-sample reach along the normal, vox */
    double      raster_du, raster_dv;
    int         level_v_to_z;    /* rotate the layout so +v follows world z */
} AsmEmitOpts;

typedef struct AsmEmitStats {
    size_t sheet_charts, sheet_verts, sheet_faces;
    size_t extras_charts, extras_verts, extras_faces;
    size_t metric_checked_charts, metric_regressions;
    double sheet_area, extras_area;
    double u_span, v_span;
    double gauge_theta;
    int    bake_rc;
    double sec;
} AsmEmitStats;

void AsmEmit_default_opts(AsmEmitOpts *o);

/* Zero only when mesh emission and every requested bake succeed. A failed bake
 * leaves the emitted source mesh available for diagnosis and sets st->bake_rc. */
int AsmEmit_run(AsmRun *run, const AsmEmitOpts *o, const char *out_dir, AsmEmitStats *st);
int AsmEmit_selftest(void);

#endif
