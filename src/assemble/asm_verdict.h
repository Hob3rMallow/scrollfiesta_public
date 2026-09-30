#ifndef ASM_VERDICT_INCLUDED
#define ASM_VERDICT_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "asm_axis.h"
#include "asm_types.h"
#include "asm_continuity.h"

/* ============================================================================
 * asm_verdict.h -- the lattice adapters: the assembled ribbon (rows = z,
 * the --solidify contract) and the verdict lattice (rows = the sheet's v,
 * every contested claimant its own vertex, material = the lineage of
 * measured joins, faces wherever the assembly claims continuity), which is
 * what ribbon_verdict measures.  See asm_verdict.c.
 * ==========================================================================*/

typedef struct AsmVerdictOpts {
    double      du, dv;         /* lattice steps, vox (2, 2) */
    double      umb_y, umb_x;   /* measurement-only axis for the cross-wrap gate */
    double      pitch;
    double      core_radius;
    const char *verdict_exe;    /* ribbon_verdict.exe */
    int         include_blobs_in_source;
    const AsmAxis *axis;          /* the axis polyline: the verdict lattice and its source are measured in
                                   * its frame (s, e1, e2) per vertex, the axis straight along +z through the
                                   * origin; NULL = world frame with umb_y/umb_x */
} AsmVerdictOpts;

typedef struct AsmVerdictStats {
    AsmContinuityStats continuity;
    /* the VERDICT lattice (rows = v): what ribbon_verdict measured */
    size_t lattice_verts, lattice_faces;
    size_t lattice_rows, lattice_cols;
    size_t cells_contested;     /* cells with a claimant a layer apart from the primary claimant */
    size_t stacked_verts;       /* the extra vertices those claimants become (each fails single-cover) */
    size_t lineages;            /* material lineages among the sheet's charts */
    size_t faces_by_relation;   /* faces across two charts joined by an accepted seam relation */
    size_t faces_by_geometry;   /* faces across unrelated charts within the wrap gate */
    size_t quads_refused;       /* quads with an edge between unrelated charts beyond the gate */
    size_t quads_refused_lineage;   /* quads with an edge between two lineages (never joined by the raster) */
    size_t long_join_edges, long_join_pairs;   /* relation edges over the wrap gate, and the joins carrying them */
    size_t long_join_edges_tail;   /* of those, within 8 vox: the per-cube trim gap (~4.4 vox) plus a lattice step */
    size_t cells_bridged_uv;       /* generated cells closing gaps between joined charts in the verdict lattice */
    /* GAPS: in-row gaps of the verdict lattice between occupied cells, split by the radius of the two
     * bounding cells: within one wrap (a hole or a mis-hop: a defect) or between wraps (the next page of
     * the strip: box physics).  Unmeasured (n = 0) without an axis. */
    size_t gap_within_n, gap_between_n;
    double gap_within_vox, gap_between_vox;          /* totals */
    double gap_within_p50, gap_within_p90, gap_between_p50, gap_between_p90;
    int    gap_measured;
    /* the ASSEMBLED ribbon (rows = z): the --solidify contract */
    size_t solid_verts, solid_faces, solid_rows, solid_cols;
    size_t solid_contested;     /* z-row cells with a claimant a layer apart (folded into the primary there) */
    size_t cells_bridged;       /* generated cells closing cube-seam rows between joined charts */
    size_t source_verts, source_faces;
    int    axis_frame;          /* 1 = measured in the local axis frame of the configured table */
    int    frame_measured;      /* 0 = no axis and no umbilicus: the verdict's cross-wrap and depth-seam gates are ABSENT */
    int    verdict_rc;
    double sec;
} AsmVerdictStats;

void AsmVerdict_default_opts(AsmVerdictOpts *o);

/* Writes <out>/assembled_ribbon.vmesh (+ _lane.i32, _material_identity.i32, _support.u8,
 * _provenance.u8, _phase.f32, _stats.json: the solid stage's input contract),
 * <out>/source.vmesh, then spawns ribbon_verdict with --report
 * <out>/verdict.json (log in <out>/verdict.log).  Returns 0 when the
 * lattice was written; the verdict's own exit code is in stats. */
int AsmVerdict_run(AsmRun *run, const AsmVerdictOpts *o, const char *out_dir, AsmVerdictStats *st);

/* The word after a gate's name in verdict.log (PASS / FAIL / ABSENT): returns the text after
 * the word (the gate's detail) or NULL when the line carries none; *pass and *measured as the
 * word says (ABSENT = not measured, counted as a failed gate by ribbon_verdict). */
const char *AsmVerdict_gate_word(const char *line, int *pass, int *measured);

int AsmVerdict_selftest(void);

#endif
