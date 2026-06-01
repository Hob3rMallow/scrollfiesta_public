#ifndef ARENA_INCLUDED
#define ARENA_INCLUDED

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

/* Allocation — never returns NULL (raises Arena_Failed on OOM) */
void   *Arena_alloc(Arena_T arena, long nbytes, const char *file, int line);
void   *Arena_calloc(Arena_T arena, long count, long size,
                     const char *file, int line);

/* Save/restore */
Arena_Mark Arena_save(Arena_T arena);
void       Arena_restore(Arena_T arena, Arena_Mark mark);

/* Convenience macros — inject __FILE__, __LINE__ for diagnostics */
#define ARENA_ALLOC(a, nbytes) \
    Arena_alloc((a), (nbytes), __FILE__, __LINE__)
#define ARENA_CALLOC(a, count, size) \
    Arena_calloc((a), (count), (size), __FILE__, __LINE__)
#define ARENA_NEW(a, p) \
    ((p) = Arena_alloc((a), (long)sizeof *(p), __FILE__, __LINE__))

#endif
