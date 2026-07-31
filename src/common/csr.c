#include "csr.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>

struct CSR_T {
    int32_t *offset;    /* [nrows + 1] */
    int32_t *target;    /* [nnz] */
    float   *weight;    /* [nnz] or NULL */
    int32_t  nrows;
    int32_t  nnz;
    void    *self;      /* validation sentinel */
};

/* Comparison function for sorting int32_t */
static int cmp_int32(const void *a, const void *b)
{
    int32_t va = *(const int32_t *)a;
    int32_t vb = *(const int32_t *)b;
    if (va < vb) { return -1; }
    if (va > vb) { return  1; }
    return 0;
}

CSR_T CSR_from_faces(Arena_T arena, const int32_t *faces, size_t nf,
                     size_t nv)
{
    assert(arena);

    CSR_T csr;
    ARENA_NEW(arena, csr);
    csr->nrows = (int32_t)nv;
    csr->weight = NULL;
    csr->self = csr;

    if (nf == 0 || nv == 0) {
        csr->offset = (int32_t *)ARENA_CALLOC(arena,
                                               (long)(nv + 1),
                                               (long)sizeof(int32_t));
        csr->target = NULL;
        csr->nnz = 0;
        return csr;
    }

    assert(faces);

    /* Pass 1: count degree (overcount: 2 neighbors per face vertex) */
    int32_t *degree = (int32_t *)ARENA_CALLOC(arena, (size_t)nv,
                                               (long)sizeof(int32_t));

    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f * 3 + 0];
        int32_t b = faces[f * 3 + 1];
        int32_t c = faces[f * 3 + 2];
        assert(a >= 0 && a < (int32_t)nv);
        assert(b >= 0 && b < (int32_t)nv);
        assert(c >= 0 && c < (int32_t)nv);
        degree[a] += 2;
        degree[b] += 2;
        degree[c] += 2;
    }

    /* Prefix sum -> offset */
    csr->offset = (int32_t *)ARENA_ALLOC(arena,
                                          (long)(nv + 1) * (long)sizeof(int32_t));
    csr->offset[0] = 0;
    for (size_t v = 0; v < nv; v++) {
        csr->offset[v + 1] = csr->offset[v] + degree[v];
    }

    int32_t total_entries = csr->offset[nv];
    int32_t *raw_target = (int32_t *)ARENA_ALLOC(arena,
                                                   (long)total_entries * (long)sizeof(int32_t));

    /* Pass 2: fill adjacency using cursor array */
    int32_t *cursor = (int32_t *)ARENA_ALLOC(arena,
                                              (long)nv * (long)sizeof(int32_t));
    memcpy(cursor, csr->offset, nv * sizeof(int32_t));

    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f * 3 + 0];
        int32_t b = faces[f * 3 + 1];
        int32_t c = faces[f * 3 + 2];

        raw_target[cursor[a]++] = b;
        raw_target[cursor[a]++] = c;
        raw_target[cursor[b]++] = a;
        raw_target[cursor[b]++] = c;
        raw_target[cursor[c]++] = a;
        raw_target[cursor[c]++] = b;
    }

    /* Dedup: sort each row and remove duplicates, compact in-place */
    int32_t write_pos = 0;
    for (size_t v = 0; v < nv; v++) {
        int32_t start = csr->offset[v];
        int32_t end   = csr->offset[v + 1];
        int32_t row_len = end - start;

        if (row_len > 1) {
            qsort(raw_target + start, (size_t)row_len, sizeof(int32_t),
                  cmp_int32);
        }

        int32_t new_start = write_pos;
        for (int32_t j = start; j < end; j++) {
            /* Skip self-loops and duplicates */
            if (raw_target[j] == (int32_t)v) {
                continue;
            }
            if (j > start && raw_target[j] == raw_target[j - 1]) {
                continue;
            }
            raw_target[write_pos++] = raw_target[j];
        }
        csr->offset[v] = new_start;
    }
    csr->offset[nv] = write_pos;
    csr->nnz = write_pos;

    /* Copy compacted target to final arena-allocated array */
    csr->target = (int32_t *)ARENA_ALLOC(arena,
                                          (long)write_pos * (long)sizeof(int32_t));
    memcpy(csr->target, raw_target, (size_t)write_pos * sizeof(int32_t));

    return csr;
}

/* COO entry for sorting */
typedef struct {
    int32_t row;
    int32_t col;
    float   val;
} COOEntry;

static int cmp_coo(const void *a, const void *b)
{
    const COOEntry *ea = (const COOEntry *)a;
    const COOEntry *eb = (const COOEntry *)b;
    if (ea->row != eb->row) {
        return (ea->row < eb->row) ? -1 : 1;
    }
    return (ea->col < eb->col) ? -1 : (ea->col > eb->col) ? 1 : 0;
}

CSR_T CSR_from_coo(Arena_T arena, const int32_t *rows, const int32_t *cols,
                   const float *vals, size_t nnz, size_t nrows)
{
    assert(arena);

    CSR_T csr;
    ARENA_NEW(arena, csr);
    csr->nrows = (int32_t)nrows;
    csr->self = csr;

    if (nnz == 0) {
        csr->offset = (int32_t *)ARENA_CALLOC(arena,
                                               (long)(nrows + 1),
                                               (long)sizeof(int32_t));
        csr->target = NULL;
        csr->weight = NULL;
        csr->nnz = 0;
        return csr;
    }

    assert(rows && cols);

    /* Sort COO entries by (row, col) */
    COOEntry *entries = (COOEntry *)ARENA_ALLOC(arena,
                                                 (long)nnz * (long)sizeof(COOEntry));
    for (size_t i = 0; i < nnz; i++) {
        entries[i].row = rows[i];
        entries[i].col = cols[i];
        entries[i].val = vals ? vals[i] : 1.0f;
    }
    qsort(entries, nnz, sizeof(COOEntry), cmp_coo);

    /* Count rows for offset */
    csr->offset = (int32_t *)ARENA_CALLOC(arena,
                                           (long)(nrows + 1),
                                           (long)sizeof(int32_t));
    for (size_t i = 0; i < nnz; i++) {
        assert(entries[i].row >= 0 && entries[i].row < (int32_t)nrows);
        csr->offset[entries[i].row + 1]++;
    }
    for (size_t i = 1; i <= nrows; i++) {
        csr->offset[i] += csr->offset[i - 1];
    }

    csr->nnz = (int32_t)nnz;
    csr->target = (int32_t *)ARENA_ALLOC(arena,
                                          (long)nnz * (long)sizeof(int32_t));
    csr->weight = (float *)ARENA_ALLOC(arena,
                                        (long)nnz * (long)sizeof(float));

    for (size_t i = 0; i < nnz; i++) {
        csr->target[i] = entries[i].col;
        csr->weight[i] = entries[i].val;
    }

    return csr;
}

CSR_T CSR_uniform_laplacian(Arena_T arena, const int32_t *faces, size_t nf,
                            size_t nv)
{
    assert(arena);
    if (nf == 0 || nv == 0) {
        /* Return empty CSR */
        return CSR_from_faces(arena, faces, nf, nv);
    }

    /* Get adjacency structure */
    CSR_T adj = CSR_from_faces(arena, faces, nf, nv);

    /* Build COO triples with weight = 1/degree */
    int32_t total_nnz = CSR_nnz(adj);
    const int32_t *off = CSR_offset(adj);
    const int32_t *tgt = CSR_target(adj);

    int32_t *coo_rows = (int32_t *)ARENA_ALLOC(arena,
                                                 (long)total_nnz * (long)sizeof(int32_t));
    int32_t *coo_cols = (int32_t *)ARENA_ALLOC(arena,
                                                 (long)total_nnz * (long)sizeof(int32_t));
    float *coo_vals = (float *)ARENA_ALLOC(arena,
                                            (long)total_nnz * (long)sizeof(float));

    for (int32_t i = 0; i < (int32_t)nv; i++) {
        int32_t degree = off[i + 1] - off[i];
        float w = (degree > 0) ? 1.0f / (float)degree : 0.0f;
        for (int32_t j = off[i]; j < off[i + 1]; j++) {
            coo_rows[j] = i;
            coo_cols[j] = tgt[j];
            coo_vals[j] = w;
        }
    }

    return CSR_from_coo(arena, coo_rows, coo_cols, coo_vals,
                        (size_t)total_nnz, nv);
}

int32_t CSR_nrows(const CSR_T csr)
{
    assert(csr && csr->self == csr);
    return csr->nrows;
}

int32_t CSR_nnz(const CSR_T csr)
{
    assert(csr && csr->self == csr);
    return csr->nnz;
}

const int32_t *CSR_offset(const CSR_T csr)
{
    assert(csr && csr->self == csr);
    return csr->offset;
}

const int32_t *CSR_target(const CSR_T csr)
{
    assert(csr && csr->self == csr);
    return csr->target;
}

const float *CSR_weight(const CSR_T csr)
{
    assert(csr && csr->self == csr);
    return csr->weight;
}

void CSR_matvec(const CSR_T csr, const float *x, float *y)
{
    assert(csr && csr->self == csr);
    assert(x && y);

    const int32_t *off = csr->offset;
    const int32_t *tgt = csr->target;
    const float   *w   = csr->weight;
    int32_t nr = csr->nrows;

    if (w != NULL) {
        for (int32_t i = 0; i < nr; i++) {
            float sum = 0.0f;
            for (int32_t j = off[i]; j < off[i + 1]; j++) {
                sum += w[j] * x[tgt[j]];
            }
            y[i] = sum;
        }
    } else {
        /* Unweighted: effectively weight = 1 */
        for (int32_t i = 0; i < nr; i++) {
            float sum = 0.0f;
            for (int32_t j = off[i]; j < off[i + 1]; j++) {
                sum += x[tgt[j]];
            }
            y[i] = sum;
        }
    }
}

void CSR_matvec3(const CSR_T csr, const float *x, float *y)
{
    assert(csr && csr->self == csr);
    assert(x && y);

    const int32_t *off = csr->offset;
    const int32_t *tgt = csr->target;
    const float   *w   = csr->weight;
    int32_t nr = csr->nrows;

    if (w != NULL) {
        for (int32_t i = 0; i < nr; i++) {
            float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f;
            for (int32_t j = off[i]; j < off[i + 1]; j++) {
                int32_t k = tgt[j];
                float wj = w[j];
                s0 += wj * x[k * 3 + 0];
                s1 += wj * x[k * 3 + 1];
                s2 += wj * x[k * 3 + 2];
            }
            y[i * 3 + 0] = s0;
            y[i * 3 + 1] = s1;
            y[i * 3 + 2] = s2;
        }
    } else {
        for (int32_t i = 0; i < nr; i++) {
            float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f;
            for (int32_t j = off[i]; j < off[i + 1]; j++) {
                int32_t k = tgt[j];
                s0 += x[k * 3 + 0];
                s1 += x[k * 3 + 1];
                s2 += x[k * 3 + 2];
            }
            y[i * 3 + 0] = s0;
            y[i * 3 + 1] = s1;
            y[i * 3 + 2] = s2;
        }
    }
}

void CSR_validate(const CSR_T csr)
{
    assert(csr && csr->self == csr);
    assert(csr->offset[0] == 0);

    for (int32_t i = 0; i < csr->nrows; i++) {
        assert(csr->offset[i] <= csr->offset[i + 1]);
    }
    assert(csr->offset[csr->nrows] == csr->nnz);

    for (int32_t i = 0; i < csr->nnz; i++) {
        assert(csr->target[i] >= 0 && csr->target[i] < csr->nrows);
    }
}
