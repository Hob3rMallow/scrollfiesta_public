#include "mesh_pile.h"
#include "ves_platform.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

static void mp_join(char *out, size_t cap, const char *a, const char *b)
{
    size_t la = strlen(a);
    if (la > 0 && (a[la - 1] == '/' || a[la - 1] == '\\'))
        snprintf(out, cap, "%s%s", a, b);
    else
        snprintf(out, cap, "%s/%s", a, b);
}

int MeshPile_parse_cube_id(const char *name, char id[24],
                           long *oz, long *oy, long *ox)
{
    long z = 0, y = 0, x = 0;
    if (name == NULL || strlen(name) < 20) return -1;
    if (name[0] != 'z' || name[6] != '_' || name[7] != 'y' ||
        name[13] != '_' || name[14] != 'x')
        return -1;
    if (sscanf(name, "z%5ld_y%5ld_x%5ld", &z, &y, &x) != 3) return -1;
    memcpy(id, name, 20);
    id[20] = '\0';
    *oz = z; *oy = y; *ox = x;
    return 0;
}

static int mp_push(MeshPileEntry *pile, size_t cap, size_t *n,
                   const char *dir, const char *name, int depth,
                   MeshPile_Filter filter, void *ctx)
{
    size_t len = strlen(name);
    char id[24] = "";
    long oz = 0, oy = 0, ox = 0;
    int has_id = 0;
    MeshPileEntry *e = NULL;
    if (len < 7 || strcmp(name + len - 6, ".vmesh") != 0) return 0;
    has_id = MeshPile_parse_cube_id(name, id, &oz, &oy, &ox) == 0;
    /* per-cube final dumps carry the id AND the _final_all marker; a plain
     * world-frame mesh (top level only) has neither restriction */
    if (has_id && strstr(name, "_final_all.vmesh") == NULL) return 0;
    if (!has_id && depth > 0) return 0;
    if (has_id && filter != NULL && !filter(ctx, oz, oy, ox)) return 0;
    if (has_id) {
        for (size_t i = 0; i < *n; i++)
            if (strcmp(pile[i].cube_id, id) == 0) return 0;   /* dedup */
    }
    if (*n >= cap) return -1;
    e = &pile[(*n)++];
    mp_join(e->path, sizeof e->path, dir, name);
    memcpy(e->cube_id, id, sizeof id);
    e->oz = oz; e->oy = oy; e->ox = ox;
    e->has_id = has_id;
    return 0;
}

static int mp_scan(const char *dir, int depth, MeshPileEntry *pile,
                   size_t cap, size_t *n, MeshPile_Filter filter, void *ctx)
{
    if (depth > 6) return 0;
#ifdef _WIN32
    {
        char glob[MESH_PILE_MAX_PATH];
        WIN32_FIND_DATAA fd;
        HANDLE h;
        snprintf(glob, sizeof glob, "%s/*", dir);
        h = FindFirstFileA(glob, &fd);
        if (h == INVALID_HANDLE_VALUE) return -1;
        do {
            if (strcmp(fd.cFileName, ".") == 0 ||
                strcmp(fd.cFileName, "..") == 0)
                continue;
            if (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
                char sub[MESH_PILE_MAX_PATH];
                mp_join(sub, sizeof sub, dir, fd.cFileName);
                if (mp_scan(sub, depth + 1, pile, cap, n, filter, ctx) != 0) {
                    FindClose(h);
                    return -1;
                }
            } else if (mp_push(pile, cap, n, dir, fd.cFileName, depth,
                               filter, ctx) != 0) {
                FindClose(h);
                return -1;
            }
        } while (FindNextFileA(h, &fd));
        FindClose(h);
    }
#else
    {
        DIR *d = opendir(dir);
        struct dirent *de = NULL;
        if (d == NULL) return -1;
        while ((de = readdir(d)) != NULL) {
            char sub[MESH_PILE_MAX_PATH];
            struct stat st;
            if (strcmp(de->d_name, ".") == 0 || strcmp(de->d_name, "..") == 0)
                continue;
            mp_join(sub, sizeof sub, dir, de->d_name);
            if (stat(sub, &st) != 0) continue;
            if (S_ISDIR(st.st_mode)) {
                if (mp_scan(sub, depth + 1, pile, cap, n, filter, ctx) != 0) {
                    closedir(d);
                    return -1;
                }
            } else if (mp_push(pile, cap, n, dir, de->d_name, depth,
                               filter, ctx) != 0) {
                closedir(d);
                return -1;
            }
        }
        closedir(d);
    }
#endif
    return 0;
}

static int mp_cmp(const void *a, const void *b)
{
    const MeshPileEntry *p = (const MeshPileEntry *)a;
    const MeshPileEntry *q = (const MeshPileEntry *)b;
    if (p->oz != q->oz) return p->oz < q->oz ? -1 : 1;
    if (p->oy != q->oy) return p->oy < q->oy ? -1 : 1;
    if (p->ox != q->ox) return p->ox < q->ox ? -1 : 1;
    return strcmp(p->path, q->path);
}

void MeshPile_sort(MeshPileEntry *pile, size_t n)
{
    if (pile != NULL && n > 1) qsort(pile, n, sizeof *pile, mp_cmp);
}

int MeshPile_scan(const char *dir, MeshPileEntry *pile, size_t cap,
                  size_t *n, MeshPile_Filter filter, void *ctx)
{
    int rc = 0;
    if (dir == NULL || pile == NULL || n == NULL) return -1;
    *n = 0;
    rc = mp_scan(dir, 0, pile, cap, n, filter, ctx);
    MeshPile_sort(pile, *n);
    return rc;
}

/* ---- selftest -------------------------------------------------------------*/

static int mp_test_filter(void *ctx, long oz, long oy, long ox)
{
    (void)ctx; (void)oy; (void)ox;
    return oz >= 4352;
}

int MeshPile_selftest(void)
{
    int fails = 0;
    char id[24] = "";
    long oz = 0, oy = 0, ox = 0;
    MeshPileEntry pile[8];
    size_t n = 0;
    if (MeshPile_parse_cube_id("z04352_y03200_x02688_step12_final_all.vmesh",
                               id, &oz, &oy, &ox) != 0 ||
        strcmp(id, "z04352_y03200_x02688") != 0 ||
        oz != 4352 || oy != 3200 || ox != 2688)
        fails++;
    if (MeshPile_parse_cube_id("welded.vmesh", id, &oz, &oy, &ox) == 0)
        fails++;
    if (MeshPile_parse_cube_id("z4352_y3200_x2688.vmesh", id, &oz, &oy, &ox)
        == 0)
        fails++;
    /* in-memory push semantics: dedup, filter, marker, depth */
    n = 0;
    if (mp_push(pile, 8, &n, "d", "z04352_y03200_x02688_step12_final_all.vmesh",
                2, mp_test_filter, NULL) != 0 || n != 1) fails++;
    if (mp_push(pile, 8, &n, "d2", "z04352_y03200_x02688_step12_final_all.vmesh",
                2, mp_test_filter, NULL) != 0 || n != 1) fails++;   /* dedup */
    if (mp_push(pile, 8, &n, "d", "z04224_y03200_x02688_step12_final_all.vmesh",
                2, mp_test_filter, NULL) != 0 || n != 1) fails++;   /* filtered */
    if (mp_push(pile, 8, &n, "d", "z04480_y03200_x02688_step11.vmesh",
                2, mp_test_filter, NULL) != 0 || n != 1) fails++;   /* no marker */
    if (mp_push(pile, 8, &n, "d", "plain.vmesh", 1, NULL, NULL) != 0 || n != 1)
        fails++;                                                    /* deep plain */
    if (mp_push(pile, 8, &n, "d", "plain.vmesh", 0, NULL, NULL) != 0 || n != 2)
        fails++;                                                    /* top plain */
    if (mp_push(pile, 8, &n, "d", "z04224_y03200_x02688_step12_final_all.vmesh",
                2, NULL, NULL) != 0 || n != 3) fails++;
    MeshPile_sort(pile, n);
    if (!(pile[0].oz == 0 && !pile[0].has_id && pile[1].oz == 4224 &&
          pile[2].oz == 4352))
        fails++;
    fprintf(stderr, "[selftest] mesh_pile %s (%d failures)\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails;
}
