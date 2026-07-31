#ifndef ATLAS_SOLUTION_INCLUDED
#define ATLAS_SOLUTION_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "../unroll/piece_set.h"
#include "atlas_overlap_fix.h"

/*
 * Versioned, authoritative handoff from the discrete/continuous atlas solve
 * to downstream ribbon fitting.
 *
 * Visualization OBJ files deliberately are not used as state: they omit chart
 * ranks, killed-face identity, and collision provenance.  The checkpoint is
 * tied to one exact PieceSet ordering by a 64-bit content fingerprint.
 */

typedef struct {
    int32_t chart;
    int32_t group;
    int32_t winding;   /* round-local tabu correction k -- diagnostic only;
                        * all zero when the alternation loop converged */
    uint64_t rank;
    double qc_max;
    double wind;       /* absolute continuous wind in turns (sense*phi
                        * convention), phase-lift structure anchored per
                        * gauge island; consumers snap per sample */
} AtlasSolutionChart;

typedef struct {
    int32_t chart0;
    int32_t chart1;
    int32_t face0;
    int32_t face1;
    int32_t allowed;
    int32_t reason;
    double min_xyz;
    double radius_delta;
} AtlasSolutionResidual;

typedef struct {
    size_t nvertices;
    size_t nfaces;
    size_t ncharts;
    size_t nresiduals;
    uint64_t piece_fingerprint;

    double *u;                  /* [nvertices] */
    double *v;                  /* [nvertices], exact axial coordinate */
    uint8_t *face_keep;         /* [nfaces] */
    int32_t *vertex_chart;      /* [nvertices], -1 if unused */
    AtlasSolutionChart *chart;  /* [ncharts] */
    AtlasSolutionResidual *residual; /* [nresiduals] */
} AtlasSolution;

uint64_t AtlasSolution_piece_fingerprint(const PieceSet *ps);

int AtlasSolution_write(const char *path,
                        const PieceSet *ps,
                        const double *u,
                        const double *v,
                        const uint8_t *face_keep,
                        const AtlasOverlapFixTabuReport *report);

int AtlasSolution_read(Arena_T arena,
                       const char *path,
                       const PieceSet *ps,
                       AtlasSolution *out);

/* Re-rank charts by (wind, producer rank) so downstream order follows the
 * physical winding order rather than the u layout, which the tabu solve
 * leaves globally scrambled.  No-op (returns 0) when the wind field spans
 * half a turn or less; returns 1 when ranks were rewritten, -1 on error. */
int AtlasSolution_rank_by_wind(Arena_T arena, AtlasSolution *solution);

/* Deterministic checkpoint round-trip fixture. */
int AtlasSolution_selftest(Arena_T arena, const char *path);

#endif
