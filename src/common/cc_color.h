#ifndef CC_COLOR_INCLUDED
#define CC_COLOR_INCLUDED

#include <stddef.h>
#include <stdint.h>

/*
 * Connected-component vertex coloring.
 *
 * Union-find over face connectivity groups vertices into connected components,
 * then each component gets a distinct golden-angle hue by DESCENDING face-count
 * rank (rank 0 = largest); components with fewer than `min_faces` faces are
 * painted dim grey ("dust"). This is the coloring we render with mesh_render
 * --vcolor to see wrap shatter / cross-wrap fusion at a glance.
 *
 * The algorithm needs only face topology + the vertex count -- positions are
 * irrelevant, so it is agnostic to the pipeline's (z,y,x) vertex order.
 */

typedef struct {
    double sat;         /* HSV saturation (0..1); default 0.62                 */
    size_t min_faces;   /* components with < this many faces -> grey; default 0 */
} CCColorOpts;

typedef struct {
    size_t ncomp;                    /* number of connected components          */
    size_t nv, nf;                   /* echoed input sizes                      */
    size_t largest_faces;            /* faces in the biggest component          */
    size_t cover50, cover90, cover99;/* #comps to cover 50/90/99% of faces      */
    size_t dust_comps, dust_faces;   /* components (and faces) below min_faces   */
} CCColorStats;

/* Fill p with defaults (sat 0.62, min_faces 0). */
void CCColor_default_opts(CCColorOpts *o);

/* Compute per-vertex RGB (each in [0,1]) into out_rgb[nv*3]. `faces` are 0-based
 * [nf*3]. `opts` may be NULL (defaults used). `stats` may be NULL. Returns the
 * component count. On allocation failure fills out_rgb with mid-grey and returns 0. */
size_t CCColor_compute(size_t nv, const int32_t *faces, size_t nf,
                       const CCColorOpts *opts, float *out_rgb, CCColorStats *stats);

/* Unit test: two disjoint triangles plus a triangle sharing one vertex with the
 * first (a bowtie) -> exactly 2 components. Returns 0 on pass, else failure count. */
int CCColor_selftest(void);

#endif /* CC_COLOR_INCLUDED */
