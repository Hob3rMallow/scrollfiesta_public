/*
 * pred_reject.c -- standalone garbage-prediction-cube auditor.
 *
 *   pred_reject <cube.tif>                         print metrics + verdict
 *   pred_reject --scan <dir> [--out f] [--csv f]   scan every *.tif in dir
 *   pred_reject --selftest                         run the detector self-test
 *
 * The --scan mode loads each cube's OWN 128^3 prediction (no halo -- the
 * verdict is per-cube), prints a one-line verdict, and optionally writes a
 * reject list (one cube_id per line) and a full per-cube stats CSV. This finds
 * the garbage slabs in a grid before any meshing is spawned.
 *
 * Exit: 0 ok, 1 IO error, 2 usage error, 3 self-test failure.
 */
#include "../common/ves_platform.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
  #include <windows.h>
#else
  #include <dirent.h>
#endif

#include "../common/arena.h"
#include "../common/tiff_io.h"
#include "../extract/pred_reject.h"

typedef struct {
    char id[128];
    char path[1024];
} Entry;

/* basename without directory or ".tif" extension. */
static void basename_noext(const char *path, char *out, size_t out_sz)
{
    const char *base = path;
    for (const char *p = path; *p; p++) {
        if (*p == '/' || *p == '\\') base = p + 1;
    }
    size_t len = strlen(base);
    if (len >= 4 && strcmp(base + len - 4, ".tif") == 0) len -= 4;
    if (len >= out_sz) len = out_sz - 1;
    memcpy(out, base, len);
    out[len] = '\0';
}

static void print_stats(const char *id, const PredRejectStats *s)
{
    printf("%-30s fg=%-7zu cc=%-7zu fill=%.4f interior=%.3f(%zu) thick=%d "
           "rect[ax%d]=%.2f run=%-3d  %s\n",
           id, s->fg_count, s->largest_cc, s->fill_frac, s->interior_frac,
           s->interior_count, s->max_thickness, s->rect_axis, s->rect_frac,
           s->max_rect_run, s->reason);
}

static int cmp_entry(const void *a, const void *b)
{
    return strcmp(((const Entry *)a)->id, ((const Entry *)b)->id);
}

/* Collect every <dir>/*.tif into a malloc'd, caller-freed Entry array. */
static int collect_tifs(const char *dir, Entry **out, size_t *n_out)
{
    size_t cap = 256, n = 0;
    Entry *e = (Entry *)calloc(cap, sizeof(Entry));
    if (!e) return -1;

#ifdef _WIN32
    char glob[1024];
    snprintf(glob, sizeof(glob), "%s/*.tif", dir);
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA(glob, &fd);
    if (h == INVALID_HANDLE_VALUE) { free(e); return -1; }
    do {
        if (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) continue;
        const char *name = fd.cFileName;
        size_t nlen = strlen(name);
        if (nlen < 5 || strcmp(name + nlen - 4, ".tif") != 0) continue;
        if (n >= cap) {
            cap *= 2;
            Entry *g = (Entry *)realloc(e, cap * sizeof(Entry));
            if (!g) { FindClose(h); free(e); return -1; }
            e = g;
        }
        memset(&e[n], 0, sizeof(Entry));
        basename_noext(name, e[n].id, sizeof(e[n].id));
        snprintf(e[n].path, sizeof(e[n].path), "%s/%s", dir, name);
        n++;
    } while (FindNextFileA(h, &fd));
    FindClose(h);
#else
    DIR *d = opendir(dir);
    if (!d) { free(e); return -1; }
    struct dirent *de = NULL;
    while ((de = readdir(d)) != NULL) {
        size_t nlen = strlen(de->d_name);
        if (nlen < 5 || strcmp(de->d_name + nlen - 4, ".tif") != 0) continue;
        if (n >= cap) {
            cap *= 2;
            Entry *g = (Entry *)realloc(e, cap * sizeof(Entry));
            if (!g) { closedir(d); free(e); return -1; }
            e = g;
        }
        memset(&e[n], 0, sizeof(Entry));
        basename_noext(de->d_name, e[n].id, sizeof(e[n].id));
        snprintf(e[n].path, sizeof(e[n].path), "%s/%s", dir, de->d_name);
        n++;
    }
    closedir(d);
#endif

    qsort(e, n, sizeof(Entry), cmp_entry);
    *out = e;
    *n_out = n;
    return 0;
}

static int run_one(Arena_T arena, const char *path, PredRejectStats *out)
{
    Arena_Mark mark = Arena_save(arena);
    uint8_t *vol = NULL;
    int D = 0, H = 0, W = 0;
    if (TiffIO_load(arena, path, &vol, &D, &H, &W) != 0) {
        fprintf(stderr, "  load failed: %s\n", path);
        Arena_restore(arena, mark);
        return -1;
    }
    PredReject_is_garbage(arena, vol, D, H, W, out);
    Arena_restore(arena, mark);
    return out->verdict;
}

static int scan_dir(Arena_T arena, const char *dir,
                    const char *out_path, const char *csv_path)
{
    Entry *e = NULL;
    size_t n = 0;
    if (collect_tifs(dir, &e, &n) != 0) {
        fprintf(stderr, "scan: cannot read dir %s\n", dir);
        return 1;
    }

    FILE *fout = out_path ? fopen(out_path, "w") : NULL;
    FILE *fcsv = csv_path ? fopen(csv_path, "w") : NULL;
    if (out_path && !fout) { fprintf(stderr, "scan: cannot write %s\n", out_path); }
    if (csv_path && !fcsv) { fprintf(stderr, "scan: cannot write %s\n", csv_path); }
    if (fcsv) {
        fprintf(fcsv, "cube_id,verdict,fill_frac,interior_frac,interior_count,"
                      "max_thickness,rect_axis,rect_frac,max_rect_run,reason\n");
    }

    size_t n_rej = 0, n_keep = 0, n_err = 0;
    for (size_t i = 0; i < n; i++) {
        PredRejectStats s;
        int v = run_one(arena, e[i].path, &s);
        if (v < 0) { n_err++; continue; }
        print_stats(e[i].id, &s);
        if (s.verdict) {
            n_rej++;
            if (fout) { fprintf(fout, "%s\n", e[i].id); fflush(fout); }
        } else {
            n_keep++;
        }
        if (fcsv) {
            fprintf(fcsv, "%s,%d,%.4f,%.4f,%zu,%d,%d,%.4f,%d,%s\n",
                    e[i].id, s.verdict, s.fill_frac, s.interior_frac,
                    s.interior_count, s.max_thickness, s.rect_axis,
                    s.rect_frac, s.max_rect_run, s.reason);
        }
    }

    if (fout) fclose(fout);
    if (fcsv) fclose(fcsv);
    free(e);

    fprintf(stderr, "scan: %zu cubes -> %zu REJECT, %zu keep, %zu err\n",
            n, n_rej, n_keep, n_err);
    if (out_path) fprintf(stderr, "  reject list -> %s\n", out_path);
    if (csv_path) fprintf(stderr, "  stats csv   -> %s\n", csv_path);
    return 0;
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s <cube.tif>                          # single-cube verdict\n"
        "       %s --scan <dir> [--out f] [--csv f]    # scan every *.tif\n"
        "       %s --selftest\n", prog, prog, prog);
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return PredReject_selftest();

    if (argc < 2) { usage(argv[0]); return 2; }

    if (strcmp(argv[1], "--scan") == 0) {
        if (argc < 3) { usage(argv[0]); return 2; }
        const char *dir = argv[2];
        const char *out = NULL, *csv = NULL;
        for (int i = 3; i < argc; i++) {
            if (!strcmp(argv[i], "--out") && i + 1 < argc) out = argv[++i];
            else if (!strcmp(argv[i], "--csv") && i + 1 < argc) csv = argv[++i];
            else { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 2; }
        }
        Arena_T arena = Arena_new();
        int rc = scan_dir(arena, dir, out, csv);
        Arena_dispose(&arena);
        return rc;
    }

    /* single cube */
    const char *path = argv[1];
    char id[128];
    basename_noext(path, id, sizeof(id));
    Arena_T arena = Arena_new();
    PredRejectStats s;
    int v = run_one(arena, path, &s);
    if (v >= 0) print_stats(id, &s);
    Arena_dispose(&arena);
    return (v < 0) ? 1 : 0;
}
