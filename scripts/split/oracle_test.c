/*
 * oracle_test.c — Test harness for Step 1: Oracle sheet counting.
 *
 * Usage: ./oracle_test <obj_dir>
 *
 * Reads all *.obj files from obj_dir (sorted alphabetically),
 * runs Oracle_count_sheets on each, prints results, and writes
 * sheets.txt to obj_dir (one sheet count per line, matching sort order
 * for consumption by bridge_cut_test).
 */
#define _POSIX_C_SOURCE 200809L

#include "split/oracle.h"
#include "common/arena.h"
#include "common/obj_io.h"
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

#define MAX_COMPS 256

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
int oracle_test_main(void) { return 0; /* needs file args */ }
#else
int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <obj_dir>\n", argv[0]);
        return 1;
    }

    const char *obj_dir = argv[1];
    Arena_T arena = Arena_new();

    /* Collect and sort OBJ filenames */
    DIR *dir = opendir(obj_dir);
    if (!dir) {
        fprintf(stderr, "Cannot open directory: %s\n", obj_dir);
        Arena_dispose(&arena);
        return 1;
    }

    char *names[MAX_COMPS];
    int n_files = 0;
    struct dirent *ent;
    while ((ent = readdir(dir)) != NULL && n_files < MAX_COMPS) {
        size_t len = strlen(ent->d_name);
        if (len < 4 || strcmp(ent->d_name + len - 4, ".obj") != 0)
            continue;
        names[n_files] = strdup(ent->d_name);
        n_files++;
    }
    closedir(dir);

    qsort(names, (size_t)n_files, sizeof(char *), cmp_str);

    /* Process each component */
    int sheet_counts[MAX_COMPS];
    struct timespec t_all0, t_all1;
    clock_gettime(CLOCK_MONOTONIC, &t_all0);

    for (int i = 0; i < n_files; i++) {
        char path[1024];
        snprintf(path, sizeof(path), "%s/%s", obj_dir, names[i]);

        Arena_Mark file_mark = Arena_save(arena);

        float   *verts = NULL;
        int32_t *faces = NULL;
        size_t   nv = 0, nf = 0;

        if (ObjIO_read(arena, path, &verts, &nv, &faces, &nf) != 0) {
            fprintf(stderr, "Failed to read: %s\n", path);
            sheet_counts[i] = 1;
            Arena_restore(arena, file_mark);
            continue;
        }

        ComponentMesh cm;
        memset(&cm, 0, sizeof(cm));
        cm.verts = verts;
        cm.faces = faces;
        cm.nv    = nv;
        cm.nf    = nf;
        cm.self  = &cm;
        PCA_normal(verts, nv, cm.pca_normal, cm.centroid);

        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);

        int sheets = Oracle_count_sheets(arena, &cm, ORACLE_UV_GRID_SIZE);

        clock_gettime(CLOCK_MONOTONIC, &t1);

        sheet_counts[i] = sheets;

        printf("  %s: %zu verts, %zu faces -> %d sheet(s)  (%.3f s)\n",
               names[i], nv, nf, sheets, tdiff(&t0, &t1));
        printf("    normal=[%.4f, %.4f, %.4f]  centroid=[%.1f, %.1f, %.1f]\n",
               (double)cm.pca_normal[0], (double)cm.pca_normal[1],
               (double)cm.pca_normal[2],
               (double)cm.centroid[0], (double)cm.centroid[1],
               (double)cm.centroid[2]);

        Arena_restore(arena, file_mark);
    }

    clock_gettime(CLOCK_MONOTONIC, &t_all1);
    printf("\nTotal: %d component(s), %.3f s\n", n_files, tdiff(&t_all0, &t_all1));

    /* Write sheets.txt to obj_dir */
    char sheets_path[1024];
    snprintf(sheets_path, sizeof(sheets_path), "%s/sheets.txt", obj_dir);

    FILE *fp = fopen(sheets_path, "w");
    if (fp) {
        for (int i = 0; i < n_files; i++) {
            fprintf(fp, "%d\n", sheet_counts[i]);
        }
        fclose(fp);
        printf("Wrote %s\n", sheets_path);
    } else {
        fprintf(stderr, "WARNING: cannot write %s\n", sheets_path);
    }

    /* Cleanup */
    for (int i = 0; i < n_files; i++)
        free(names[i]);
    Arena_dispose(&arena);
    return 0;
}
#endif /* !TEST_HARNESS */
