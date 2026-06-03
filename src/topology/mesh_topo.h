#ifndef MESH_TOPO_INCLUDED
#define MESH_TOPO_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"

/* ============================================================================
 * mesh_topo.h -- triangle-mesh topology analysis for the flatten/cut stage.
 *
 * Edge classification (by #incident triangles), boundary-loop count, Euler
 * characteristic and genus, and non-manifold / unreferenced-vertex detection.
 * Intended to be run per connected component. Used to diagnose why a sheet
 * won't flatten (LSCM needs a manifold topological disk) and to drive the
 * Seamster cut.
 * ==========================================================================*/

typedef struct {
    long   n_verts, n_faces, n_edges;
    long   n_boundary_edges;     /* incident to exactly 1 triangle */
    long   n_manifold_edges;     /* exactly 2 */
    long   n_nonmanifold_edges;  /* 3 or more */
    long   n_boundary_loops;     /* connected components of the boundary graph */
    long   euler;                /* V - E + F */
    double genus;                /* (2 - boundary_loops - euler)/2 (manifold) */
    long   n_unref_verts;        /* vertices in no triangle */
    int    is_disk;              /* 1 iff manifold, 1 boundary loop, genus 0 */
} MeshTopoInfo;

/* Analyze a (preferably single-component) mesh. verts may be NULL (topology
 * only uses faces). Returns 0 on success, -1 on empty input. */
int MeshTopo_analyze(Arena_T arena, const float *verts, size_t nv,
                     const int32_t *faces, size_t nf, MeshTopoInfo *out);

/* In-process unit tests (disk / annulus / hole). Returns 0 if all pass. */
int MeshTopo_selftest(void);

#endif
