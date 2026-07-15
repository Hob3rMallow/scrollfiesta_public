/*
 * sf_mesh_convert.c — memory model, validation, and the xyz<->zyx boundary
 * conversions of the public C API (see sf_internal.h).
 */
#include "sf_internal.h"

#include "../common/pca.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ══════════════════════════════════════════════════════════════════════════
 * Public memory model: everything the library hands out is sf_malloc'd and
 * released through these exported functions, so a host built with a
 * different CRT never frees across the boundary.
 * ════════════════════════════════════════════════════════════════════════ */

SF_API void *sf_malloc(size_t n)
{
    return malloc(n ? n : 1);
}

SF_API void sf_free(void *p)
{
    free(p);
}

SF_API void sf_mesh_free(sf_mesh *m)
{
    if (!m)
        return;
    free(m->vertices);
    free(m->faces);
    free(m->vertex_normals);
    free(m->pin_mask);
    free(m->vmap);
    memset(m, 0, sizeof *m);
}

SF_API void sf_mesh_list_free(sf_mesh_list *l)
{
    if (!l)
        return;
    for (size_t i = 0; i < l->count; i++)
        sf_mesh_free(&l->items[i]);
    free(l->items);
    memset(l, 0, sizeof *l);
}

SF_API void sf_handle_loops_free(sf_handle_loop *loops, size_t n)
{
    if (!loops)
        return;
    for (size_t i = 0; i < n; i++)
        free(loops[i].vertex_indices);
    free(loops);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Validation
 * ════════════════════════════════════════════════════════════════════════ */

static sf_status validate_fail(char *err, size_t err_len, const char *msg)
{
    if (err && err_len)
        snprintf(err, err_len, "%s", msg);
    return SF_ERROR_BAD_ARG;
}

SF_API sf_status sf_mesh_validate(const sf_mesh *m, char *err, size_t err_len)
{
    if (!m)
        return validate_fail(err, err_len, "mesh is NULL");
    if (!m->vertices || m->n_vertices == 0)
        return validate_fail(err, err_len, "mesh has no vertices");
    if (!m->faces || m->n_faces == 0)
        return validate_fail(err, err_len, "mesh has no faces");

    for (size_t i = 0; i < m->n_vertices * 3; i++) {
        if (!isfinite(m->vertices[i]))
            return validate_fail(err, err_len, "mesh has non-finite vertex coordinates");
    }
    for (size_t i = 0; i < m->n_faces * 3; i++) {
        int32_t v = m->faces[i];
        if (v < 0 || (size_t)v >= m->n_vertices)
            return validate_fail(err, err_len, "face index out of range");
    }
    if (err && err_len)
        err[0] = '\0';
    return SF_OK;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Boundary conversions (arena-backed; may RAISE Arena_Failed)
 * ════════════════════════════════════════════════════════════════════════ */

float *sf_xyz_to_zyx(Arena_T a, const float *xyz, size_t n)
{
    float *out = ARENA_ALLOC(a, n * 3 * sizeof *out);
    for (size_t i = 0; i < n; i++) {
        out[i * 3 + 0] = xyz[i * 3 + 2];
        out[i * 3 + 1] = xyz[i * 3 + 1];
        out[i * 3 + 2] = xyz[i * 3 + 0];
    }
    return out;
}

int32_t *sf_faces_swap(Arena_T a, const int32_t *faces, size_t nf)
{
    int32_t *out = ARENA_ALLOC(a, nf * 3 * sizeof *out);
    for (size_t i = 0; i < nf; i++) {
        out[i * 3 + 0] = faces[i * 3 + 0];
        out[i * 3 + 1] = faces[i * 3 + 2];
        out[i * 3 + 2] = faces[i * 3 + 1];
    }
    return out;
}

int sf_cm_from_mesh(Arena_T a, const sf_mesh *in, ComponentMesh *cm)
{
    memset(cm, 0, sizeof *cm);
    cm->verts = sf_xyz_to_zyx(a, in->vertices, in->n_vertices);
    cm->faces = sf_faces_swap(a, in->faces, in->n_faces);
    cm->nv    = in->n_vertices;
    cm->nf    = in->n_faces;
    if (in->vertex_normals)
        cm->vert_normals = sf_xyz_to_zyx(a, in->vertex_normals, in->n_vertices);
    if (in->pin_mask) {
        cm->pin_mask = ARENA_ALLOC(a, in->n_vertices);
        memcpy(cm->pin_mask, in->pin_mask, in->n_vertices);
    }
    cm->comp_id = 1;
    if (PCA_normal(cm->verts, cm->nv, cm->pca_normal, cm->centroid) != 0) {
        cm->pca_normal[0] = 1.0f;
        cm->pca_normal[1] = 0.0f;
        cm->pca_normal[2] = 0.0f;
        cm->centroid[0] = cm->centroid[1] = cm->centroid[2] = 0.0f;
    }
    cm->self = cm;
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Copy-out (malloc; never raises; frees partial work on failure)
 * ════════════════════════════════════════════════════════════════════════ */

sf_status sf_mesh_out(sf_mesh *out,
                      const float *verts_zyx, size_t nv,
                      const int32_t *faces_internal, size_t nf,
                      const float *normals_zyx,
                      const uint8_t *pins,
                      const int32_t *vmap)
{
    memset(out, 0, sizeof *out);
    out->vertices = sf_malloc(nv * 3 * sizeof(float));
    out->faces    = sf_malloc(nf * 3 * sizeof(int32_t));
    if (!out->vertices || !out->faces)
        goto oom;
    for (size_t i = 0; i < nv; i++) {
        out->vertices[i * 3 + 0] = verts_zyx[i * 3 + 2];
        out->vertices[i * 3 + 1] = verts_zyx[i * 3 + 1];
        out->vertices[i * 3 + 2] = verts_zyx[i * 3 + 0];
    }
    for (size_t i = 0; i < nf; i++) {
        out->faces[i * 3 + 0] = faces_internal[i * 3 + 0];
        out->faces[i * 3 + 1] = faces_internal[i * 3 + 2];
        out->faces[i * 3 + 2] = faces_internal[i * 3 + 1];
    }
    if (normals_zyx) {
        out->vertex_normals = sf_malloc(nv * 3 * sizeof(float));
        if (!out->vertex_normals)
            goto oom;
        for (size_t i = 0; i < nv; i++) {
            out->vertex_normals[i * 3 + 0] = normals_zyx[i * 3 + 2];
            out->vertex_normals[i * 3 + 1] = normals_zyx[i * 3 + 1];
            out->vertex_normals[i * 3 + 2] = normals_zyx[i * 3 + 0];
        }
    }
    if (pins) {
        out->pin_mask = sf_malloc(nv);
        if (!out->pin_mask)
            goto oom;
        memcpy(out->pin_mask, pins, nv);
    }
    if (vmap) {
        out->vmap = sf_malloc(nv * sizeof(int32_t));
        if (!out->vmap)
            goto oom;
        memcpy(out->vmap, vmap, nv * sizeof(int32_t));
    }
    out->n_vertices = nv;
    out->n_faces    = nf;
    return SF_OK;

oom:
    sf_mesh_free(out);
    return SF_ERROR_OOM;
}

sf_status sf_mesh_out_cm(sf_mesh *out, const ComponentMesh *cm, const int32_t *vmap)
{
    return sf_mesh_out(out, cm->verts, cm->nv, cm->faces, cm->nf,
                       cm->vert_normals, cm->pin_mask, vmap);
}

sf_status sf_mesh_out_raw(sf_mesh *out,
                          const float *verts_xyz, size_t nv,
                          const int32_t *faces, size_t nf,
                          const float *normals_xyz,
                          const uint8_t *pins,
                          const int32_t *vmap)
{
    memset(out, 0, sizeof *out);
    out->vertices = sf_malloc(nv * 3 * sizeof(float));
    out->faces    = sf_malloc(nf * 3 * sizeof(int32_t));
    if (!out->vertices || !out->faces)
        goto oom;
    memcpy(out->vertices, verts_xyz, nv * 3 * sizeof(float));
    memcpy(out->faces, faces, nf * 3 * sizeof(int32_t));
    if (normals_xyz) {
        out->vertex_normals = sf_malloc(nv * 3 * sizeof(float));
        if (!out->vertex_normals)
            goto oom;
        memcpy(out->vertex_normals, normals_xyz, nv * 3 * sizeof(float));
    }
    if (pins) {
        out->pin_mask = sf_malloc(nv);
        if (!out->pin_mask)
            goto oom;
        memcpy(out->pin_mask, pins, nv);
    }
    if (vmap) {
        out->vmap = sf_malloc(nv * sizeof(int32_t));
        if (!out->vmap)
            goto oom;
        memcpy(out->vmap, vmap, nv * sizeof(int32_t));
    }
    out->n_vertices = nv;
    out->n_faces    = nf;
    return SF_OK;

oom:
    sf_mesh_free(out);
    return SF_ERROR_OOM;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Provenance by exact position (bit-identical float match)
 * ════════════════════════════════════════════════════════════════════════ */

typedef struct {
    uint32_t k[3];   /* the three floats, bit pattern */
    int32_t  idx;
} SfPosEntry;

static int pos_entry_cmp(const void *pa, const void *pb)
{
    const SfPosEntry *a = (const SfPosEntry *)pa;
    const SfPosEntry *b = (const SfPosEntry *)pb;
    for (int i = 0; i < 3; i++) {
        if (a->k[i] < b->k[i]) return -1;
        if (a->k[i] > b->k[i]) return 1;
    }
    if (a->idx < b->idx) return -1;
    if (a->idx > b->idx) return 1;
    return 0;
}

static void pos_key(uint32_t k[3], const float *p)
{
    memcpy(k, p, 3 * sizeof(uint32_t));
}

static int pos_key_cmp(const uint32_t a[3], const uint32_t b[3])
{
    for (int i = 0; i < 3; i++) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

int32_t *sf_vmap_by_position(Arena_T a, const float *in_verts, size_t in_nv,
                             const float *out_verts, size_t out_nv)
{
    Arena_Mark  mark    = Arena_save(a);
    SfPosEntry *entries = ARENA_ALLOC(a, in_nv * sizeof *entries);
    uint8_t    *used    = ARENA_CALLOC(a, in_nv, 1);
    int32_t    *vmap    = ARENA_ALLOC(a, out_nv * sizeof *vmap);

    for (size_t i = 0; i < in_nv; i++) {
        pos_key(entries[i].k, &in_verts[i * 3]);
        entries[i].idx = (int32_t)i;
    }
    qsort(entries, in_nv, sizeof *entries, pos_entry_cmp);

    for (size_t i = 0; i < out_nv; i++) {
        uint32_t key[3];
        pos_key(key, &out_verts[i * 3]);

        /* lower_bound over entries by key */
        size_t lo = 0, hi = in_nv;
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (pos_key_cmp(entries[mid].k, key) < 0)
                lo = mid + 1;
            else
                hi = mid;
        }

        int32_t match = -1;
        if (lo < in_nv && pos_key_cmp(entries[lo].k, key) == 0) {
            /* prefer the first unconsumed duplicate; else re-use the first */
            size_t j = lo;
            match = entries[lo].idx;
            while (j < in_nv && pos_key_cmp(entries[j].k, key) == 0) {
                if (!used[entries[j].idx]) {
                    match = entries[j].idx;
                    used[entries[j].idx] = 1;
                    break;
                }
                j++;
            }
        }
        vmap[i] = match;
    }

    /* keep only the vmap allocation live */
    int32_t *result = vmap;
    (void)mark; /* entries/used stay allocated until op arena disposal — fine
                 * for per-op arenas; do NOT restore (would free vmap too). */
    return result;
}

int32_t *sf_vmap_identity(Arena_T a, size_t nv)
{
    int32_t *vmap = ARENA_ALLOC(a, nv * sizeof *vmap);
    for (size_t i = 0; i < nv; i++)
        vmap[i] = (int32_t)i;
    return vmap;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Misc mesh utilities
 * ════════════════════════════════════════════════════════════════════════ */

static int32_t uf_find(int32_t *parent, int32_t x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

long sf_count_components(Arena_T a, size_t nv, const int32_t *faces, size_t nf)
{
    Arena_Mark mark    = Arena_save(a);
    int32_t   *parent  = ARENA_ALLOC(a, nv * sizeof *parent);
    uint8_t   *touched = ARENA_CALLOC(a, nv, 1);
    long       count   = 0;

    for (size_t i = 0; i < nv; i++)
        parent[i] = (int32_t)i;
    for (size_t f = 0; f < nf; f++) {
        int32_t v0 = faces[f * 3 + 0];
        int32_t v1 = faces[f * 3 + 1];
        int32_t v2 = faces[f * 3 + 2];
        touched[v0] = touched[v1] = touched[v2] = 1;
        int32_t r0 = uf_find(parent, v0);
        int32_t r1 = uf_find(parent, v1);
        int32_t r2 = uf_find(parent, v2);
        if (r1 != r0) parent[r1] = r0;
        r2 = uf_find(parent, v2);
        r0 = uf_find(parent, v0);
        if (r2 != r0) parent[r2] = r0;
    }
    for (size_t i = 0; i < nv; i++) {
        if (touched[i] && uf_find(parent, (int32_t)i) == (int32_t)i)
            count++;
    }
    Arena_restore(a, mark);
    return count;
}

int sf_concat_cms(Arena_T a, const ComponentMesh *ms, size_t n, ComponentMesh *out)
{
    size_t total_nv = 0, total_nf = 0;
    int    all_normals = 1, any_pins = 0;

    memset(out, 0, sizeof *out);
    for (size_t i = 0; i < n; i++) {
        total_nv += ms[i].nv;
        total_nf += ms[i].nf;
        if (!ms[i].vert_normals)
            all_normals = 0;
        if (ms[i].pin_mask)
            any_pins = 1;
    }
    if (total_nv == 0)
        return -1;

    out->verts = ARENA_ALLOC(a, total_nv * 3 * sizeof(float));
    out->faces = ARENA_ALLOC(a, total_nf * 3 * sizeof(int32_t));
    if (all_normals)
        out->vert_normals = ARENA_ALLOC(a, total_nv * 3 * sizeof(float));
    if (any_pins)
        out->pin_mask = ARENA_CALLOC(a, total_nv, 1);

    size_t voff = 0, foff = 0;
    for (size_t i = 0; i < n; i++) {
        memcpy(out->verts + voff * 3, ms[i].verts, ms[i].nv * 3 * sizeof(float));
        if (all_normals)
            memcpy(out->vert_normals + voff * 3, ms[i].vert_normals,
                   ms[i].nv * 3 * sizeof(float));
        if (any_pins && ms[i].pin_mask)
            memcpy(out->pin_mask + voff, ms[i].pin_mask, ms[i].nv);
        for (size_t f = 0; f < ms[i].nf * 3; f++)
            out->faces[foff * 3 + f] = ms[i].faces[f] + (int32_t)voff;
        voff += ms[i].nv;
        foff += ms[i].nf;
    }
    out->nv = total_nv;
    out->nf = total_nf;
    out->comp_id = 1;
    if (PCA_normal(out->verts, out->nv, out->pca_normal, out->centroid) != 0) {
        out->pca_normal[0] = 1.0f;
        out->pca_normal[1] = 0.0f;
        out->pca_normal[2] = 0.0f;
    }
    out->self = out;
    return 0;
}
