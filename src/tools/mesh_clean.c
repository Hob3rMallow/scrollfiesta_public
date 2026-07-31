/* mesh_clean — make an OBJ 2-manifold via ManifoldGuard (resolve non-manifold
 * edges, split bowtie/pinch vertices; optional re-orient). Used to mop up the
 * few CVT-dual bowties left by a standalone per-component RVD remesh (the
 * per-cube pipeline runs the same guard automatically).
 *
 *   mesh_clean <in.obj> <out.obj> [--reorient]
 */
#include "../remesh/manifold_guard.h"
#include "../common/mesh_types.h"
#include "../common/obj_io.h"
#include "../common/arena.h"

#include <stdio.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <in.obj> <out.obj> [--reorient]\n", argv[0]);
        return 2;
    }
    int reorient = 0;
    for (int i = 3; i < argc; i++)
        if (!strcmp(argv[i], "--reorient")) reorient = 1;

    Arena_T a = Arena_new();
    float *V = NULL; size_t nv = 0; int32_t *F = NULL; size_t nf = 0;
    if (ObjIO_read(a, argv[1], &V, &nv, &F, &nf) != 0) {
        fprintf(stderr, "error: cannot read %s\n", argv[1]);
        Arena_dispose(&a); return 1;
    }

    ComponentMesh cm;
    memset(&cm, 0, sizeof(cm));
    cm.verts = V; cm.nv = nv;
    cm.faces = F; cm.nf = nf;
    cm.comp_id = 1;
    cm.self = &cm;

    ManifoldGuardStats mg;
    memset(&mg, 0, sizeof(mg));
    if (ManifoldGuard_process(a, &cm, 1, reorient, &mg) != 0) {
        fprintf(stderr, "error: ManifoldGuard_process failed\n");
        Arena_dispose(&a); return 1;
    }
    fprintf(stderr, "mesh_clean: %zu->%zu v  %zu->%zu f  "
            "(nm_edges=%zu faces_del=%zu pinch_splits=%zu orient_flips=%zu)\n",
            nv, cm.nv, nf, cm.nf,
            mg.nm_edges_resolved, mg.faces_deleted, mg.pinch_splits, mg.orient_flips);

    if (ObjIO_write(argv[2], cm.verts, cm.nv, cm.faces, cm.nf) != 0) {
        fprintf(stderr, "error: cannot write %s\n", argv[2]);
        Arena_dispose(&a); return 1;
    }
    Arena_dispose(&a);
    return 0;
}
