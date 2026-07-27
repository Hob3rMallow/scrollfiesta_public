#include "mesh_trim.h"

#include <assert.h>
#include <string.h>

int Mesh_trim_to_owned_box(Arena_T arena,
                           const float *verts_in, size_t nv_in,
                           const int32_t *faces_in, size_t nf_in,
                           const uint8_t *pin_mask_in,
                           float owned_lo, float owned_hi,
                           float **out_verts, size_t *out_nv,
                           int32_t **out_faces, size_t *out_nf,
                           uint8_t **out_pin_mask)
{
    assert(arena);
    assert(verts_in && faces_in);
    assert(out_verts && out_nv && out_faces && out_nf);

    *out_verts = NULL;
    *out_nv = 0;
    *out_faces = NULL;
    *out_nf = 0;
    if (out_pin_mask) *out_pin_mask = NULL;
    if (nf_in == 0 || nv_in == 0) return 0;

    /* Pass 1: keep faces with ALL THREE vertices inside [owned_lo, owned_hi).
     * Strict (all-vertex), not centroid: a centroid test let a vertex poke
     * ~edge/2 past the boundary, which shrank the inter-cube gap at an insetted
     * seam back into z-fight range (~0.5 vox). All-vertex containment makes the
     * inset exact, so adjacent cubes end up a clean 2*inset apart. */
    int32_t *kept_faces = (int32_t *)ARENA_ALLOC(arena,
                              (long)nf_in * 3L * (long)sizeof(int32_t));
    size_t nf_kept = 0;
    for (size_t f = 0; f < nf_in; f++) {
        int32_t a = faces_in[f * 3 + 0];
        int32_t b = faces_in[f * 3 + 1];
        int32_t c = faces_in[f * 3 + 2];
        if ((size_t)a >= nv_in || (size_t)b >= nv_in || (size_t)c >= nv_in) {
            continue;  /* bogus index, skip defensively */
        }
        int32_t tri[3] = { a, b, c };
        int inside = 1;
        for (int t = 0; t < 3; t++) {
            float vz = verts_in[tri[t]*3+0];
            float vy = verts_in[tri[t]*3+1];
            float vx = verts_in[tri[t]*3+2];
            if (vz < owned_lo || vz >= owned_hi ||
                vy < owned_lo || vy >= owned_hi ||
                vx < owned_lo || vx >= owned_hi) { inside = 0; break; }
        }
        if (!inside) continue;
        kept_faces[nf_kept * 3 + 0] = a;
        kept_faces[nf_kept * 3 + 1] = b;
        kept_faces[nf_kept * 3 + 2] = c;
        nf_kept++;
    }
    if (nf_kept == 0) return 0;

    /* Pass 2: mark referenced vertices. */
    uint8_t *used = (uint8_t *)ARENA_CALLOC(arena, (long)nv_in, 1L);
    for (size_t f = 0; f < nf_kept; f++) {
        used[kept_faces[f * 3 + 0]] = 1;
        used[kept_faces[f * 3 + 1]] = 1;
        used[kept_faces[f * 3 + 2]] = 1;
    }

    /* Pass 3: build remap and compact verts. */
    int32_t *remap = (int32_t *)ARENA_ALLOC(arena,
                         (long)nv_in * (long)sizeof(int32_t));
    size_t nv_kept = 0;
    for (size_t v = 0; v < nv_in; v++) {
        if (used[v]) {
            remap[v] = (int32_t)nv_kept;
            nv_kept++;
        } else {
            remap[v] = -1;
        }
    }

    float *new_verts = (float *)ARENA_ALLOC(arena,
                          (long)nv_kept * 3L * (long)sizeof(float));
    size_t wi = 0;
    for (size_t v = 0; v < nv_in; v++) {
        if (used[v]) {
            new_verts[wi * 3 + 0] = verts_in[v * 3 + 0];
            new_verts[wi * 3 + 1] = verts_in[v * 3 + 1];
            new_verts[wi * 3 + 2] = verts_in[v * 3 + 2];
            wi++;
        }
    }

    int32_t *new_faces = (int32_t *)ARENA_ALLOC(arena,
                            (long)nf_kept * 3L * (long)sizeof(int32_t));
    for (size_t f = 0; f < nf_kept; f++) {
        new_faces[f * 3 + 0] = remap[kept_faces[f * 3 + 0]];
        new_faces[f * 3 + 1] = remap[kept_faces[f * 3 + 1]];
        new_faces[f * 3 + 2] = remap[kept_faces[f * 3 + 2]];
    }

    *out_verts = new_verts;
    *out_nv = nv_kept;
    *out_faces = new_faces;
    *out_nf = nf_kept;

    /* Pin mask: pure index remap (no vert position changes). */
    if (pin_mask_in && out_pin_mask) {
        uint8_t *new_pins = (uint8_t *)ARENA_CALLOC(arena,
                                (long)nv_kept, 1L);
        for (size_t v = 0; v < nv_in; v++) {
            if (used[v] && pin_mask_in[v]) {
                new_pins[remap[v]] = 1;
            }
        }
        *out_pin_mask = new_pins;
    }

    return 0;
}
