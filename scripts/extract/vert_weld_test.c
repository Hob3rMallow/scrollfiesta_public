/*
 * vert_weld_test.c -- unit tests for the orientation guard in Weld_verts
 * (src/common/vert_weld.c). The micro-weld in grid_weld merges exactly-
 * coincident vertex pairs; without a guard it can fuse two opposite-facing
 * (recto vs verso) surfaces, creating a non-orientable `same_dir` edge that a
 * winding-repair BFS cannot fix. `guard_orient` rejects those merges.
 *
 * Gates (synthetic meshes; verts stored (z,y,x), all at z=0 so winding
 * normals lie along +-z; same_dir measured with MeshManifold_audit):
 *   W1  recto-verso must NOT merge: two coincident triangles sharing edge
 *       p->q wound OPPOSITELY (normals +z / -z). guard on -> not merged
 *       (out_nv=6, same_dir=0); guard off -> merged (out_nv=4, same_dir=1,
 *       reproduces the pathology the fix targets).
 *   W2  no over-blocking: same coincident pair but ANTI-PARALLEL shared edge
 *       (both normals +z, consistent). guard on -> still merges (out_nv=4,
 *       same_dir=0) -- legitimate slit zips are preserved.
 *   W3  point-dedup callers unaffected: nf=0 with OPPOSING in_normals and
 *       guard off -> still merges (out_nv=1). Proves the mesh_extract /
 *       mesh_resplit double-envelope collapse path is untouched.
 *
 * Build: vert_weld_test.vcxproj (links vert_weld.c + mesh_manifold.c +
 *        arena.c + except.c). Run: exit 0 = all pass, 1 = a gate failed.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/arena.h"
#include "common/vert_weld.h"
#include "common/mesh_manifold.h"

static int fails = 0;
#define CHECK(cond, ...) do { \
        if (!(cond)) { printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); fails++; } \
        else        { printf("  ok:   "); printf(__VA_ARGS__); printf("\n"); } \
    } while (0)

static size_t count_same_dir(Arena_T a, size_t nv, const int32_t *faces, size_t nf)
{
    MeshManifoldStats s = MeshManifold_audit(a, nv, faces, nf);
    return s.same_dir_edges;
}

/* Two triangles sharing coincident verts 0==3 and 1==4; verts 2,5 differ so
 * the pair is a "slit". T2's winding is chosen by the caller (recto/verso). */
static const float V6[6*3] = {
    0,0,0,   /* 0  p  */
    0,1,0,   /* 1  q  */
    0,0,1,   /* 2  r  */
    0,0,0,   /* 3  p' (coincident with 0) */
    0,1,0,   /* 4  q' (coincident with 1) */
    0,0,-1,  /* 5  s  */
};

int main(void)
{
    printf("=== vert_weld_test ===\n");
    Arena_T arena = Arena_new();

    /* ---- W1: recto-verso, T2 wound SAME direction on the shared edge -> opposing normals ---- */
    {
        float v[6*3]; memcpy(v, V6, sizeof V6);
        int32_t f_off[6] = { 0,1,2,  3,4,5 };   /* T2 edge 3->4 same dir as T1 0->1 */
        int32_t f_on[6];  memcpy(f_on, f_off, sizeof f_off);
        float *ov = NULL; size_t onv = 0, onf = 0;

        Arena_Mark m = Arena_save(arena);
        Weld_verts(arena, v, 6, NULL, f_off, 2, &onf, 1e-3f, /*guard=*/false, &ov, &onv, NULL);
        size_t sd_off = count_same_dir(arena, onv, f_off, onf);
        CHECK(onv == 4 && sd_off == 1,
              "W1 guard OFF: recto-verso fuses -> out_nv=%zu (exp 4), same_dir=%zu (exp 1)", onv, sd_off);
        Arena_restore(arena, m);

        m = Arena_save(arena);
        Weld_verts(arena, v, 6, NULL, f_on, 2, &onf, 1e-3f, /*guard=*/true, &ov, &onv, NULL);
        size_t sd_on = count_same_dir(arena, onv, f_on, onf);
        CHECK(onv == 6 && sd_on == 0,
              "W1 guard ON:  recto-verso refused -> out_nv=%zu (exp 6), same_dir=%zu (exp 0)", onv, sd_on);
        Arena_restore(arena, m);
    }

    /* ---- W2: legitimate slit, T2 wound ANTI-PARALLEL (edge 4->3) -> agreeing normals ---- */
    {
        float v[6*3]; memcpy(v, V6, sizeof V6);
        int32_t f[6] = { 0,1,2,  4,3,5 };   /* T2 edge 4->3 opposite dir to T1 0->1 */
        float *ov = NULL; size_t onv = 0, onf = 0;

        Arena_Mark m = Arena_save(arena);
        Weld_verts(arena, v, 6, NULL, f, 2, &onf, 1e-3f, /*guard=*/true, &ov, &onv, NULL);
        size_t sd = count_same_dir(arena, onv, f, onf);
        CHECK(onv == 4 && sd == 0,
              "W2 guard ON:  legit slit still zips -> out_nv=%zu (exp 4), same_dir=%zu (exp 0)", onv, sd);
        Arena_restore(arena, m);
    }

    /* ---- W3: point-dedup path (nf=0) with opposing in_normals must still merge ---- */
    {
        float v[2*3]  = { 0,0,0,  0,0,0 };          /* coincident */
        float n[2*3]  = { 1,0,0,  -1,0,0 };         /* opposing normals (recto/verso) */
        int32_t dummy = 0; size_t onf = 0;
        float *ov = NULL, *on = NULL; size_t onv = 0;

        Arena_Mark m = Arena_save(arena);
        Weld_verts(arena, v, 2, n, &dummy, 0, &onf, 1e-3f, /*guard=*/false, &ov, &onv, &on);
        CHECK(onv == 1,
              "W3 nf=0 dedup: opposing-normal pair still merges -> out_nv=%zu (exp 1)", onv);
        Arena_restore(arena, m);
    }

    Arena_dispose(&arena);
    printf("=== vert_weld_test %s (%d failure%s) ===\n",
           fails ? "FAILED" : "PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 1 : 0;
}
