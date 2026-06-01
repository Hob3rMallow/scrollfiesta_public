#ifndef CSR_INCLUDED
#define CSR_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "arena.h"

typedef struct CSR_T *CSR_T;

/* Build from mesh faces - undirected adjacency, deduped.
 * Two passes: degree count -> prefix sum -> fill.
 * O(V + E), no sorting. All arrays arena-allocated. */
CSR_T CSR_from_faces(Arena_T arena, const int32_t *faces, size_t nf,
                     size_t nv);

/* Build from COO (coordinate) triples - for weighted matrices.
 * Expects pre-sorted or will sort internally. */
CSR_T CSR_from_coo(Arena_T arena, const int32_t *rows, const int32_t *cols,
                   const float *vals, size_t nnz, size_t nrows);

/* Access */
int32_t  CSR_nrows(const CSR_T csr);
int32_t  CSR_nnz(const CSR_T csr);
const int32_t *CSR_offset(const CSR_T csr);   /* [nrows + 1] */
const int32_t *CSR_target(const CSR_T csr);    /* [nnz] */
const float   *CSR_weight(const CSR_T csr);    /* [nnz] or NULL */

/* Build row-stochastic uniform Laplacian from mesh faces.
 * L[i,j] = 1/degree_i for each neighbor j.
 * CSR_matvec3(L, v, out) computes the neighbor-mean (smoothing). */
CSR_T CSR_uniform_laplacian(Arena_T arena, const int32_t *faces, size_t nf,
                            size_t nv);

/* Sparse mat-vec: y = A @ x (3-channel: x and y are [N*3]) */
void CSR_matvec3(const CSR_T csr, const float *x, float *y);

/* Sparse mat-vec: y = A @ x (single channel) */
void CSR_matvec(const CSR_T csr, const float *x, float *y);

/* Validation (debug builds) */
void CSR_validate(const CSR_T csr);

#endif
