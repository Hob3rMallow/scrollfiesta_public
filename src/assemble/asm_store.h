#ifndef ASM_STORE_INCLUDED
#define ASM_STORE_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"
#include "asm_types.h"

/* ============================================================================
 * asm_store.h -- per-cube chart store on disk.
 *
 *   <dir>/<cube_id>.asc     fixed-width little-endian binary, one file per
 *                           cube, published by atomic rename.
 *
 * The store makes every stage resumable and keeps the run's memory bounded
 * by what a stage needs: stage 1 writes charts as it finishes each cube,
 * later stages read them back.  Extras are the same records with flags.
 *
 * IDENTITY.  A record is only reusable for the SAME input and the SAME
 * cleaning: every file carries the source mesh's path, size and mtime and a
 * fingerprint of the cleaning options + ASM_STORE_CODE_VERSION, and a read
 * with a different identity reports the record stale (-2) instead of
 * reusing it (reviewer, 2026-09-09: resume validated nothing).  Bump
 * ASM_STORE_CODE_VERSION whenever asm_clean / asm_flatten change what they
 * write.
 * ==========================================================================*/

#define ASM_STORE_CODE_VERSION "asm-clean-2026-09-14-chart-area-recovery"

/* FNV-1a of a text: the caller composes the cleaning options into it */
uint64_t AsmStore_fingerprint(const char *text);

int AsmStore_write_cube(const char *dir, const char *cube_id,
                        const char *src_path, uint64_t fingerprint,
                        const AsmChart *charts, size_t n);

/* Reads the charts of one cube into `arena`.  Returns 0; -1 when the file is
 * absent or malformed; -2 when it was written for another source file (path,
 * size or mtime) or another fingerprint (out_n is then 0). */
int AsmStore_read_cube(Arena_T arena, const char *dir, const char *cube_id,
                       const char *src_path, uint64_t fingerprint,
                       AsmChart **out_charts, size_t *out_n);

int AsmStore_exists(const char *dir, const char *cube_id);

int AsmStore_selftest(void);

/* Complete native placement/repair checkpoint. Preserves global float64 UV,
 * omitted charts, source frames and the entire proposal ledger. A versioned
 * ABI header and payload SHA reject incompatible or damaged snapshots.
 * v2 also retains opaque, uncommitted repair state under the same checksum;
 * v1 remains readable and is written when no such state exists. */
int AsmStore_write_run(const char *path, const AsmRun *run);
int AsmStore_read_run(const char *path, AsmRun *run);

#endif
