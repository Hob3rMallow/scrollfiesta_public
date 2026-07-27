#include "arena.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

/* Alignment union — guarantees 16-byte alignment on x86_64 [Hanson Ch.6] */
union align {
    double      d;
    long        l;
    void       *p;
    long double ld;
};

#define ALIGNMENT (sizeof(union align))

/* 64 MB default chunk size [common.md §1] */
#define CHUNK_SIZE ((size_t)(64 * 1024 * 1024))

/* Max free list length */
#define MAX_FREE_CHUNKS 10

struct Arena_chunk {
    struct Arena_chunk *prev;
    char               *avail;
    char               *limit;
    /* data follows immediately after this header */
};

struct Arena_T {
    struct Arena_chunk *current;
    struct Arena_chunk *freelist;
    int                 free_count;
};

/* Round a byte count up to ALIGNMENT. size_t throughout so >2 GB requests do
 * not truncate on LLP64 (Win64). */
static size_t roundup(size_t nbytes)
{
    size_t rem = nbytes % ALIGNMENT;
    if (rem != 0) {
        nbytes += ALIGNMENT - rem;
    }
    return nbytes;
}

static struct Arena_chunk *chunk_new(size_t min_bytes)
{
    size_t header_size = roundup(sizeof(struct Arena_chunk));
    size_t data_size = min_bytes > CHUNK_SIZE ? min_bytes : CHUNK_SIZE;
    size_t total = header_size + data_size;
    struct Arena_chunk *c = malloc(total);
    if (c == NULL) {
        return NULL;
    }
    c->prev  = NULL;
    c->avail = (char *)c + header_size;
    c->limit = (char *)c + total;
    return c;
}

Arena_T Arena_new(void)
{
    Arena_T arena = malloc(sizeof(*arena));
    if (arena == NULL) {
        RAISE(Arena_Failed);
        return NULL; /* unreachable */
    }
    arena->current    = NULL;
    arena->freelist   = NULL;
    arena->free_count = 0;
    return arena;
}

void Arena_dispose(Arena_T *ap)
{
    assert(ap && *ap);
    Arena_free(*ap);

    /* Free the free list */
    struct Arena_chunk *c = (*ap)->freelist;
    while (c != NULL) {
        struct Arena_chunk *prev = c->prev;
        free(c);
        c = prev;
    }

    free(*ap);
    *ap = NULL;
}

void Arena_free(Arena_T arena)
{
    assert(arena);
    /* Move all active chunks to the free list (up to MAX_FREE_CHUNKS) */
    struct Arena_chunk *c = arena->current;
    while (c != NULL) {
        struct Arena_chunk *prev = c->prev;
        if (arena->free_count < MAX_FREE_CHUNKS) {
            /* Reset chunk for reuse */
            size_t header_size = roundup(sizeof(struct Arena_chunk));
            c->avail = (char *)c + header_size;
            c->prev  = arena->freelist;
            arena->freelist = c;
            arena->free_count++;
        } else {
            free(c);
        }
        c = prev;
    }
    arena->current = NULL;
}

void *Arena_alloc(Arena_T arena, size_t nbytes, const char *file, int line)
{
    assert(arena);
    /* Empty/degenerate input (e.g. a CSR with 0 adjacency entries, an isolated-
     * face mesh) legitimately asks for 0 bytes; hand back a valid, unique,
     * minimally-sized pointer instead of aborting -- a 0-length array is never
     * dereferenced. Rule 18: handle empty input gracefully, never crash. */
    if (nbytes == 0) { nbytes = 1; }

    nbytes = roundup(nbytes);

    /* Try current chunk (avail <= limit always, so the diff is non-negative) */
    struct Arena_chunk *c = arena->current;
    if (c != NULL && (size_t)(c->limit - c->avail) >= nbytes) {
        void *ptr = c->avail;
        c->avail += nbytes;
        return ptr;
    }

    /* Try free list */
    struct Arena_chunk **pp = &arena->freelist;
    while (*pp != NULL) {
        struct Arena_chunk *fc = *pp;
        if ((size_t)(fc->limit - fc->avail) >= nbytes) {
            /* Unlink from free list */
            *pp = fc->prev;
            arena->free_count--;
            /* Push onto active list */
            fc->prev = arena->current;
            arena->current = fc;
            void *ptr = fc->avail;
            fc->avail += nbytes;
            return ptr;
        }
        pp = &(*pp)->prev;
    }

    /* Allocate new chunk */
    struct Arena_chunk *nc = chunk_new(nbytes);
    if (nc == NULL) {
        if (file) {
            fprintf(stderr, "Arena_alloc: OOM at %s:%d (%zu bytes)\n",
                    file, line, nbytes);
        }
        RAISE(Arena_Failed);
        return NULL; /* unreachable */
    }
    nc->prev = arena->current;
    arena->current = nc;
    void *ptr = nc->avail;
    nc->avail += nbytes;
    return ptr;
}

void *Arena_calloc(Arena_T arena, size_t count, size_t size,
                   const char *file, int line)
{
    /* A zero count or size flows through to Arena_alloc (which returns a valid
     * pointer for a 0-byte request); memset of 0 bytes is a no-op. Guard
     * count*size overflow BEFORE it wraps (size_t is unsigned). */
    if (size != 0 && count > SIZE_MAX / size) {
        if (file) {
            fprintf(stderr, "Arena_calloc: size overflow at %s:%d "
                    "(%zu * %zu)\n", file, line, count, size);
        }
        RAISE(Arena_Failed);
        return NULL; /* unreachable */
    }
    size_t nbytes = count * size;
    void *ptr = Arena_alloc(arena, nbytes, file, line);
    memset(ptr, 0, nbytes);
    return ptr;
}

Arena_Mark Arena_save(Arena_T arena)
{
    Arena_Mark mark;
    assert(arena);
    if (arena->current == NULL) {
        /* Force a chunk so we have something to save */
        struct Arena_chunk *nc = chunk_new(0);
        if (nc == NULL) {
            RAISE(Arena_Failed);
        }
        nc->prev = NULL;
        arena->current = nc;
    }
    mark.avail = arena->current->avail;
    mark.chunk = arena->current;
    return mark;
}

void Arena_restore(Arena_T arena, Arena_Mark mark)
{
    assert(arena);
    /* Free any chunks allocated after the mark */
    while (arena->current != mark.chunk) {
        struct Arena_chunk *c = arena->current;
        assert(c != NULL);
        arena->current = c->prev;
        /* Move to free list or free */
        if (arena->free_count < MAX_FREE_CHUNKS) {
            size_t header_size = roundup(sizeof(struct Arena_chunk));
            c->avail = (char *)c + header_size;
            c->prev  = arena->freelist;
            arena->freelist = c;
            arena->free_count++;
        } else {
            free(c);
        }
    }
    arena->current->avail = mark.avail;
}
