#ifndef BLOCK_LIST_INCLUDED
#define BLOCK_LIST_INCLUDED

/* ============================================================================
 * block_list.h -- a sheet's column blocks named in a FILE instead of on the
 * command line (2026-10-25: the PHerc0358 32^3 sheet's 569 blocks were a
 * 60 k-character command line, past Windows' 32 k limit -- WinError 206 at
 * the end of a 2.3 h prep).  One block per line, "<col0> <raster path>", the
 * path the rest of the line (spaces allowed); blank lines and '#' comments
 * skipped.  Tools take it as "@<file>" in place of the "<col0> <raster>"
 * argument pairs (sheet_strip --plan/--repack, sheet_layers --sheet-compose).
 * ==========================================================================*/

/* *blocks_out: 2n strings in argv order (col0, path, col0, path, ...), ONE
 * malloc the caller frees.  0 = ok; -1 = unreadable, a line without a path, or
 * no blocks (a message on stderr). */
int BlockList_read(const char *path, int *n_out, char ***blocks_out);

/* "@<file>" -> <file>, else NULL. */
const char *BlockList_arg(const char *arg);

/* 0 = pass: comments, blank lines, CRLF, a path with spaces, a bad line refused. Scratch under <repo>/selftest/. */
int BlockList_selftest(void);

#endif
