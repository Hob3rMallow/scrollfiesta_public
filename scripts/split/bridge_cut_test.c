#define _POSIX_C_SOURCE 200809L

#include "split/bridge_cut.h"
#include "split/oracle.h"
#include "common/arena.h"
#include "common/obj_io.h"
#include "common/obj_colors.h"
#include "common/pca.h"
#include "common/pipeline_constants.h"
#include "common/mesh_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef _MSC_VER
#include <dirent.h>
#endif

static double tdiff(struct timespec *a, struct timespec *b)
{
    return (double)(b->tv_sec - a->tv_sec) +
           (double)(b->tv_nsec - a->tv_nsec) / 1e9;
}

static int cmp_str(const void *a, const void *b)
{
    return strcmp(*(const char **)a, *(const char **)b);
}

#ifdef TEST_HARNESS
int bridge_cut_test_main(void) { return 0; /* needs file args */ }
#else
int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr,
                "Usage: %s <obj_dir> <sheet_counts_file> <out_dir>\n",
                argv[0]);
        return 1;
    }

    const char *obj_dir = argv[1];
    const char *sheet_file = argv[2];
    const char *out_dir = argv[3];

    /* Read sheet counts (one int per line) */
    FILE *sf = fopen(sheet_file, "r");
    if (!sf) {
        fprintf(stderr, "Cannot open sheet file: %s\n", sheet_file);
        return 1;
    }
    int sheet_counts[64];
    int n_sheets_read = 0;
    while (n_sheets_read < 64 &&
           fscanf(sf, "%d", &sheet_counts[n_sheets_read]) == 1) {
        n_sheets_read++;
    }
    fclose(sf);

    /* Collect and sort OBJ filenames */
    DIR *dir = opendir(obj_dir);
    if (!dir) {
        fprintf(stderr, "Cannot open directory: %s\n", obj_dir);
        return 1;
    }

    char *names[64];
    int n_files = 0;
    struct dirent *ent;
    while ((ent = readdir(dir)) != NULL && n_files < 64) {
        size_t len = strlen(ent->d_name);
        if (len < 4 || strcmp(ent->d_name + len - 4, ".obj") != 0)
            continue;
        names[n_files] = strdup(ent->d_name);
        n_files++;
    }
    closedir(dir);

    qsort(names, (size_t)n_files, sizeof(char *), cmp_str);

    Arena_T arena = Arena_new();

    struct timespec t_all0, t_all1;
    clock_gettime(CLOCK_MONOTONIC, &t_all0);

    int total_pieces = 0;
    int next_color = 0;

    for (int i = 0; i < n_files; i++) {
        char path[1024];
        snprintf(path, sizeof(path), "%s/%s", obj_dir, names[i]);

        Arena_Mark file_mark = Arena_save(arena);

        float   *verts = NULL;
        int32_t *faces = NULL;
        size_t   nv = 0, nf = 0;

        if (ObjIO_read(arena, path, &verts, &nv, &faces, &nf) != 0) {
            fprintf(stderr, "Failed to read: %s\n", path);
            Arena_restore(arena, file_mark);
            continue;
        }

        /* Build ComponentMesh */
        ComponentMesh cm;
        memset(&cm, 0, sizeof(cm));
        cm.verts = verts;
        cm.faces = faces;
        cm.nv    = nv;
        cm.nf    = nf;
        PCA_normal(verts, nv, cm.pca_normal, cm.centroid);
        cm.self = &cm;

        int sheets = (i < n_sheets_read) ? sheet_counts[i] : 1;

        if (sheets <= 1) {
            /* Single sheet: pass through */
            char outpath[1024];
            snprintf(outpath, sizeof(outpath), "%s/%s", out_dir, names[i]);
            ObjIO_write_colored(outpath, verts, nv, faces, nf,
                                OBJ_COLORS[next_color % OBJ_NUM_COLORS]);
            next_color++;
            printf("  %s: %zu verts, %zu faces, %d sheet(s) "
                   "-> pass through\n",
                   names[i], nv, nf, sheets);
            total_pieces++;
        } else {
            /* Multi-sheet: run bridge cut */
            int parent_color = next_color;
            next_color++;

            struct timespec t0, t1;
            clock_gettime(CLOCK_MONOTONIC, &t0);

            ComponentMesh *out_meshes = NULL;
            size_t out_count = 0;

            BridgeCut_process(arena, &cm, ORACLE_UV_GRID_SIZE,
                              30.0, &out_meshes, &out_count);

            clock_gettime(CLOCK_MONOTONIC, &t1);

            printf("  %s: %zu verts, %zu faces, %d sheet(s) "
                   "-> %zu piece(s)  (%.3f s)\n",
                   names[i], nv, nf, sheets,
                   out_count, tdiff(&t0, &t1));

            /* Write output pieces */
            for (size_t p = 0; p < out_count; p++) {
                /* Name: comp_XXX_pYYY.obj */
                char outpath[1024];
                snprintf(outpath, sizeof(outpath), "%s/%.*s_p%03d.obj",
                         out_dir,
                         (int)(strlen(names[i]) - 4), names[i],
                         (int)p);
                int cidx = (p == 0) ? parent_color
                                    : next_color++;
                ObjIO_write_colored(outpath, out_meshes[p].verts,
                                    out_meshes[p].nv,
                                    out_meshes[p].faces, out_meshes[p].nf,
                                    OBJ_COLORS[cidx % OBJ_NUM_COLORS]);
                printf("    piece %d: %zu verts, %zu faces\n",
                       (int)p, out_meshes[p].nv, out_meshes[p].nf);
            }
            total_pieces += (int)out_count;
        }

        Arena_restore(arena, file_mark);
    }

    clock_gettime(CLOCK_MONOTONIC, &t_all1);
    printf("\nTotal: %d file(s), %d piece(s), %.3f s\n",
           n_files, total_pieces, tdiff(&t_all0, &t_all1));

    /* Cleanup */
    for (int i = 0; i < n_files; i++) {
        free(names[i]);
    }
    Arena_dispose(&arena);
    return 0;
}
#endif /* !TEST_HARNESS */
