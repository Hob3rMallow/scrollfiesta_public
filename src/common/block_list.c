/* block_list.c -- see block_list.h. */

#include "block_list.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ves_platform.h"

int BlockList_read(const char *path, int *n_out, char ***blocks_out)
{
    FILE *f = NULL;
    char *text = NULL, *store = NULL, *line = NULL, **blocks = NULL;
    size_t cap = 0, len = 0, got = 0, i = 0;
    int n = 0, lines = 0;
    if (n_out == NULL || blocks_out == NULL) return -1;
    *n_out = 0; *blocks_out = NULL;
    if (path == NULL || (f = fopen(path, "rb")) == NULL) {
        fprintf(stderr, "block_list: cannot read the block list %s\n", path ? path : "(null)");
        return -1;
    }
    do {                                                   /* the whole file, then the lines in place */
        if (len + 65536 + 1 > cap) {
            char *grown = (char *)realloc(text, cap = cap ? cap * 2 : 262144);
            if (grown == NULL) { free(text); fclose(f); return -1; }
            text = grown;
        }
        got = fread(text + len, 1, 65536, f);
        len += got;
    } while (got > 0);
    fclose(f);
    text[len] = '\0';
    for (i = 0; i < len; i++) if (text[i] == '\n') lines++;
    blocks = (char **)malloc(((size_t)lines + 1) * 2 * sizeof *blocks + len + 1);
    if (blocks == NULL) { free(text); return -1; }
    store = (char *)(blocks + ((size_t)lines + 1) * 2);    /* the text lives behind the pointer table: one free */
    memcpy(store, text, len + 1);
    free(text);
    line = store;
    while (*line != '\0') {
        char *end = strchr(line, '\n'), *sp = NULL, *next = NULL;
        size_t L = 0;
        next = end != NULL ? end + 1 : line + strlen(line);
        if (end != NULL) *end = '\0';
        L = strlen(line);
        while (L > 0 && (line[L - 1] == '\r' || line[L - 1] == ' ' || line[L - 1] == '\t')) line[--L] = '\0';
        while (*line == ' ' || *line == '\t') line++;
        if (*line != '\0' && *line != '#') {
            sp = strchr(line, ' ');
            if (sp == NULL || sp == line) {
                fprintf(stderr, "block_list: %s: bad block line '%s' (want: <col0> <raster>)\n", path, line);
                free(blocks);
                return -1;
            }
            *sp = '\0';
            while (sp[1] == ' ' || sp[1] == '\t') sp++;
            blocks[(size_t)n * 2] = line;
            blocks[(size_t)n * 2 + 1] = sp + 1;
            n++;
        }
        line = next;
    }
    if (n == 0) { fprintf(stderr, "block_list: %s: no blocks\n", path); free(blocks); return -1; }
    *n_out = n; *blocks_out = blocks;
    return 0;
}

const char *BlockList_arg(const char *arg)
{
    return arg != NULL && arg[0] == '@' && arg[1] != '\0' ? arg + 1 : NULL;
}

int BlockList_selftest(void)
{
    int fails = 0, n = 0;
    char dir[1024], list[1200], bad[1200];
    char **blocks = NULL;
    FILE *f = NULL;
    if (ves_selftest_dir(dir, sizeof dir, "block_list") != 0) { fprintf(stderr, "[selftest] block list: no scratch dir\n"); return 1; }
    snprintf(list, sizeof list, "%s/blocks.txt", dir);
    snprintf(bad, sizeof bad, "%s/bad.txt", dir);
    f = fopen(list, "wb");
    if (f == NULL || fputs("# the blocks\r\n\r\n0 C:/a dir/b00/x_sheet_mask.tif\r\n  37   D:\\b01\\mask.tif  \n", f) < 0 || fclose(f) != 0) fails++;
    else if (BlockList_read(list, &n, &blocks) != 0 || n != 2) { fprintf(stderr, "  FAIL: list parse (%d blocks)\n", n); fails++; }
    else {
        fails += strcmp(blocks[0], "0") != 0 || strcmp(blocks[1], "C:/a dir/b00/x_sheet_mask.tif") != 0 ||
                 strcmp(blocks[2], "37") != 0 || strcmp(blocks[3], "D:\\b01\\mask.tif") != 0;
        if (fails) fprintf(stderr, "  FAIL: parsed '%s' '%s' '%s' '%s'\n", blocks[0], blocks[1], blocks[2], blocks[3]);
    }
    free(blocks); blocks = NULL;
    f = fopen(bad, "wb");
    if (f == NULL || fputs("12\n", f) < 0 || fclose(f) != 0) fails++;
    else if (BlockList_read(bad, &n, &blocks) == 0) { fprintf(stderr, "  FAIL: a line without a path accepted\n"); fails++; free(blocks); blocks = NULL; }
    fails += BlockList_arg("@x.txt") == NULL || strcmp(BlockList_arg("@x.txt"), "x.txt") != 0 || BlockList_arg("x.txt") != NULL ||
             BlockList_arg("@") != NULL || BlockList_arg(NULL) != NULL;
    remove(list); remove(bad);
    if (fails == 0) ves_selftest_clean(dir);
    fprintf(stderr, "[selftest] block list %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}
