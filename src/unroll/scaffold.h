#ifndef SCAFFOLD_INCLUDED
#define SCAFFOLD_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include "../common/arena.h"
#include "piece_set.h"

/* ============================================================================
 * scaffold.h -- the L1 SCAFFOLD: structure (not parameterization) over a
 * scroll_whole placed dir. Per-vertex (wrap k, phase phi, z), derived from the
 * registered winding phase, plus the pinned spiral calibration. The scaffold
 * is what local repair and the export atlas gate on; it is NOT itself a
 * parameterization.
 *
 * u(phi)   = a*phi + b*phi^2/(4pi)      (arc length along the spiral)
 * r(phi)   = a + b*phi/(2pi)            (radius from the umbilicus)
 * k        = floor(sense*phi / 2pi)     (integer wrap index, monotone outward)
 *
 * phi is the REGISTERED winding angle (third channel of <id>_uvphi.f32), so it
 * is continuous and multi-turn; k is a pure function of it (derived, never
 * stored -- E0 verifies this is safe). `sense` flips the winding direction so
 * k increases outward regardless of the scroll's handedness.
 * ==========================================================================*/

#ifndef SCAFFOLD_2PI
#define SCAFFOLD_2PI 6.283185307179586476925286766559
#endif

/* Pinned spiral calibration + axis, read off placed_index.json's calibration
 * block. One calibration is shared across every cube (and every slab) -- it is
 * fit on the seed cube and pinned grid-wide. */
typedef struct {
    double spiral_a, spiral_b;   /* u(phi) = a*phi + b*phi^2/(4pi) */
    double pitch;                /* vox/turn (== |spiral_b|) */
    int    sense;                /* winding sense (+/-1); k = floor(sense*phi/2pi) */
    float  axis_point[3];        /* (z,y,x) umbilicus */
    float  axis_dir[3];          /* (z,y,x) scroll axis (unit) */
} ScaffoldCalib;

/* 0139 defaults (umbilicus (0,3405,2878), pitch 9.5, sense -1). */
void Scaffold_calib_default(ScaffoldCalib *c);

/* Fill from <placed_dir>/placed_index.json (missing keys keep defaults, incl.
 * the calibration{ spiral_a, spiral_b, sense } sub-block). Returns 0 if the
 * file was read, -1 otherwise. */
int Scaffold_read_calib(const char *placed_dir, ScaffoldCalib *c);

/* Wrap index k = floor(sense*phi / 2pi). */
static inline int Scaffold_wrap_index(const ScaffoldCalib *c, double phi)
{
    return (int)floor((double)c->sense * phi / SCAFFOLD_2PI);
}

/* Unwound turn coordinate (radians, increasing outward): w = sense*phi. */
static inline double Scaffold_unwound(const ScaffoldCalib *c, double phi)
{
    return (double)c->sense * phi;
}

/* ---- E0: derived-k coverage + per-face winding coherence ----
 * Answers "is per-face k well-defined?" (piece selection depends on it) and
 * "how many turns does the assembly span?" (sanity vs the gold turn count). */
typedef struct {
    size_t nv_used;           /* verts referenced by kept faces */
    size_t nf_used;           /* kept faces */
    int    k_min, k_max;      /* over nv_used verts */
    int    k_extent;          /* k_max - k_min + 1 (turns spanned) */
    size_t face_spread_ge_pi; /* faces whose 3 verts span >= pi of unwound phi */
    double face_spread_frac;  /* face_spread_ge_pi / nf_used */
    double p50_face_spread;   /* median per-face unwound-phi spread (rad) */
    double p99_face_spread;   /* 99th pct (rad) */
    size_t n_groups;          /* distinct (cube, gid>=0) groups */
    size_t group_spread_ge_pi;/* groups spanning >= pi of unwound phi (informative;
                                 core wraps legitimately exceed) */
} ScaffoldKStats;

/* Compute E0 stats over the piece set. If hist_path != NULL, also writes a
 * k x cube-z coverage histogram (JSON) there. Returns 0, -1 on empty input. */
int Scaffold_kstats(Arena_T arena, const PieceSet *ps, const ScaffoldCalib *c,
                    ScaffoldKStats *out, const char *hist_path);

int Scaffold_selftest(void);

#endif
