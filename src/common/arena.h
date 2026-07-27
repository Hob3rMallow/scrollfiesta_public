#ifndef ARENA_INCLUDED
#define ARENA_INCLUDED

#include <stddef.h>

#include "except.h"

typedef struct Arena_T *Arena_T;

/* Forward-declare chunk for Arena_Mark */
struct Arena_chunk;

/* Save/restore checkpoints for scratch memory */
typedef struct {
    char               *avail;
    struct Arena_chunk  *chunk;
} Arena_Mark;

/* Lifecycle */
Arena_T Arena_new(void);
void    Arena_dispose(Arena_T *ap);
void    Arena_free(Arena_T arena);

/* Allocation — never returns NULL (raises Arena_Failed on OOM).
 * Sizes are size_t (not long): a single allocation may exceed 2 GB on LLP64
 * (Win64), where `long` is 32-bit and would silently truncate. */
void   *Arena_alloc(Arena_T arena, size_t nbytes, const char *file, int line);
void   *Arena_calloc(Arena_T arena, size_t count, size_t size,
                     const char *file, int line);

/* Save/restore */
Arena_Mark Arena_save(Arena_T arena);
void       Arena_restore(Arena_T arena, Arena_Mark mark);

/* Convenience macros — inject __FILE__, __LINE__ for diagnostics.
 * The explicit (size_t) casts keep the ~1000 existing call sites (many of which
 * still pass a legacy `(long)(...)` size) warning-clean under gcc -Wconversion
 * without editing them; new code passing a plain size_t gets the full >2 GB
 * range. Do NOT wrap the size in (long) in new code — that would re-truncate. */
#define ARENA_ALLOC(a, nbytes) \
    Arena_alloc((a), (size_t)(nbytes), __FILE__, __LINE__)
#define ARENA_CALLOC(a, count, size) \
    Arena_calloc((a), (size_t)(count), (size_t)(size), __FILE__, __LINE__)
#define ARENA_NEW(a, p) \
    ((p) = Arena_alloc((a), sizeof *(p), __FILE__, __LINE__))

#endif
