#ifndef SNAP_CG_INCLUDED
#define SNAP_CG_INCLUDED

#include <stddef.h>
#include <stdint.h>
#include "arena.h"
#include "csr.h"

/* ILU(0)-preconditioned CG solve for Laplacian snap-back.
 *
 * Solves  (L_graph + diag(w)) v = diag(w) target  per channel (z, y, x).
 * L_graph = D - A  is the combinatorial graph Laplacian from adjacency adj.
 * w[i] = per-vertex data weight.  w[i]=0 -> pure Laplacian smoothing.
 *
 * adj must be unweighted adjacency from CSR_from_faces (sorted cols,
 * no self-loops).  verts[nv*3] is both initial guess and output. */
int SnapCG_solve(Arena_T arena, const CSR_T adj,
                 const float *w, const float *target,
                 float *verts, size_t nv,
                 int max_iter, double tol);

#endif
