#include "extract/mesh_extract.h"
#include "common/arena.h"
#include "common/obj_io.h"
#include "common/obj_colors.h"
#include "common/mesh_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef _MSC_VER
#include <sys/stat.h>
#include <sys/types.h>
#endif

static double timespec_diff(struct timespec *a, struct timespec *b)
{
    return (double)(b->tv_sec - a->tv_sec) +
           (double)(b->tv_nsec - a->tv_nsec) / 1e9;
}

#ifdef TEST_HARNESS
int step0_test_main(void) { return 0; /* needs file args */ }
#else
int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input.tiff> <output_dir>\n", argv[0]);
        return 1;
    }

    const char *tiff_path = argv[0 + 1];
    const char *out_dir   = argv[0 + 2];

    /* Create output directory (ignore error if exists) */
    mkdir(out_dir, 0755);

    Arena_T arena = Arena_new();

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    ComponentMesh *meshes = NULL;
    size_t n_meshes = 0;
    int rc = MeshExtract_run(arena, tiff_path, NULL, 0,
                             0, 0, 0, 1, NULL, NULL, 0,
                             &meshes, &n_meshes);

    clock_gettime(CLOCK_MONOTONIC, &t1);

    if (rc != 0) {
        fprintf(stderr, "MeshExtract_run failed: %d\n", rc);
        Arena_dispose(&arena);
        return 1;
    }

    fprintf(stdout, "Step 0 complete: %.3f s, %zu component(s)\n",
            timespec_diff(&t0, &t1), n_meshes);

    for (size_t i = 0; i < n_meshes; i++) {
        ComponentMesh *cm = &meshes[i];
        fprintf(stdout, "  comp %d: %zu verts, %zu faces, "
                "normal=(%.3f, %.3f, %.3f), "
                "centroid=(%.1f, %.1f, %.1f)\n",
                cm->comp_id, cm->nv, cm->nf,
                (double)cm->pca_normal[0], (double)cm->pca_normal[1],
                (double)cm->pca_normal[2],
                (double)cm->centroid[0], (double)cm->centroid[1],
                (double)cm->centroid[2]);

        /* Write OBJ with per-component color */
        char path[512];
        snprintf(path, sizeof(path), "%s/comp_%03d.obj", out_dir, cm->comp_id);
        if (ObjIO_write_colored(path, cm->verts, cm->nv, cm->faces, cm->nf,
                                OBJ_COLORS[i % OBJ_NUM_COLORS]) != 0) {
            fprintf(stderr, "Warning: failed to write %s\n", path);
        } else {
            fprintf(stdout, "    -> %s\n", path);
        }
    }

    Arena_dispose(&arena);
    return 0;
}
#endif /* !TEST_HARNESS */
