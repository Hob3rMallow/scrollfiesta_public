#ifndef __BUF_H_
#define __BUF_H_

// ============================================================================
// STRETCHY BUFFERS
// ============================================================================
//
// C stretchy buffers, a la Sean Barrett; heavily based on an implementation
// attributed to Per Vognsen. Supports custom allocators.
//
// STANDALONE USAGE:
//   Define BUF_STANDALONE before including to use without allocators.h.
//   In standalone mode, uses standard malloc/realloc/free.
//
//   #define BUF_STANDALONE
//   #include "buf.h"
//
// WITH CUSTOM ALLOCATORS:
//   Include allocators.h before buf.h (or don't define BUF_STANDALONE).
//   This allows stretchy buffers to use custom memory allocators.
//
// USAGE:
//   int* arr = NULL;
//   buf_push(arr, 10);
//   buf_push(arr, 20);
//   for (int i = 0; i < buf_len(arr); i++) { ... }
//   buf_free(arr);
//
// ============================================================================

#include <string.h>    // for memmove, memcpy
#include <stdlib.h>    // for malloc/realloc/free (standalone mode)

#ifndef BUF_STANDALONE
#include "../core/allocators.h"
#endif

#ifdef _WIN32
#pragma warning(disable:4200)  // nonstandard extension: zero-sized array in struct/union
#endif

// ============================================================================
// Standalone allocator shim
// ============================================================================

#ifdef BUF_STANDALONE
struct Allocator
{
	void*	(*alloc)(size_t s);
	void*	(*realloc)(void* ptr, size_t new_size);
	void	(*free)(void* ptr);
};

static inline void* buf__std_alloc(size_t s) { return malloc(s); }
static inline void* buf__std_realloc(void* p, size_t s) { return realloc(p, s); }
static inline void  buf__std_free(void* p) { free(p); }

static inline Allocator* buf__get_default_allocator(void)
{
	static Allocator a = { buf__std_alloc, buf__std_realloc, buf__std_free };
	return &a;
}
#define allocator_get_cpp_allocator() buf__get_default_allocator()
#endif

// ============================================================================
// Buffer internals
// ============================================================================

#ifndef buf_offsetof
#define buf_offsetof(TYPE, MEMBER) ((size_t)&((TYPE *)0)->MEMBER)
#endif

#ifndef buf_containerof
#define buf_containerof(PTR, TYPE, MEMBER) \
	((TYPE *)((char *)(PTR) - buf_offsetof(TYPE, MEMBER)))
#endif

struct BufHeader
{
	size_t length;
	size_t capacity;
	struct Allocator* allocator;
	char data[];
};

// Returns the length of a buffer; 0 if NULL
#define buf_len(a)		((a) ? buf_get_internals(a)->length : 0)

// Returns the current storage capacity of a buffer; 0 if NULL
#define buf_cap(a)		((a) ? buf_get_internals(a)->capacity : 0)

// Frees a buffer using the allocator that was originally used to allocate it, and nulls the pointer
#define buf_free(a)		do { if (a) { buf_get_internals(a)->allocator->free(buf_get_internals(a)); a = nullptr; } } while(0)
// Determines the amount of capacity of a buffer; expands it to at least 16 elements, otherwise, 3/2 of its old capacity.
// `old_cap` MUST be parenthesised in the body -- callers pass compound expressions
// like `buf_cap(a) + n`, and without the parens C precedence parses
// `3 * buf_cap(a) + n / 2` as `(3 * buf_cap(a)) + (n / 2)` rather than
// `3 * (buf_cap(a) + n) / 2`. That under-allocates and the next memcpy walks
// off the end of the new buffer (caught by buf_push_n on growth).
#define buf__new_cap(old_cap) ((3 * (old_cap) / 2) < 16 ? 16 : (3 * (old_cap) / 2))

// Returns TRUE if this will overflow the buffer; FALSE otherwise
#define buf__needsgrow(a,n) ((a) == nullptr || buf_len(a) + (n) >= buf_cap(a))

// Maybe expands a buffer by n elements if the capacity isn't where it needs to be.
#define buf__maybegrow(a,n) (buf__needsgrow(a,(n)) ? buf__sbgrow(a,buf__new_cap(buf_cap(a)+(n))) : 0)
#define buf__sbgrow(a,n)      (*((void **)&(a)) = buf_do_realloc_((a), (n), sizeof(*(a))))

// ============================================================================
// Function declarations
// ============================================================================

#ifndef BUF_STANDALONE
void* buf_alloc_ex(struct Allocator* a);
void* buf_alloc(void);
#endif
/* static inline to match the BUF_IMPLEMENTATION definition below: g++
 * rejects extern-then-static (MSVC permissively accepts it). The sole
 * consumer is the implementing TU (uinteger.cpp). */
static inline void* buf_do_realloc_(void *a, size_t nr, size_t sz);

// Backing function for buf_clone(): copies the `elem_size`-byte elements of the
// stretchy buffer `src` into a fresh buffer. Returns NULL for a NULL/empty src.
void* buf_clone_(const void* src, size_t elem_size);

// Append a null-terminated string `s` (without its trailing '\0') to a
// stretchy `char*` buffer. The buffer is NOT auto null-terminated. NULL `s`
// is a no-op.
void buf_push_string(char** buf, const char* s);

// Append a printf-style formatted string to a stretchy `char*` buffer.
// Like buf_push_string, does not null-terminate. Internally uses a small
// scratch then a heap fallback if the result is large.
void buf_printf(char** buf, const char* fmt, ...);

// ============================================================================
// Macros
// ============================================================================

// Iterate over every element of `a`. `a` is bound to a local so it's
// evaluated exactly once per loop -- `buf_for_each(it, get_some_buf())`
// won't call `get_some_buf` per iteration.
#define buf_for_each(it, a)						\
	for (auto __bfe_buf = (a), it = __bfe_buf; (it) < __bfe_buf + buf_len(__bfe_buf); (it)++)

// Pre: buf_len(a) >= 1. NULL or empty buffers are UB (length-1 underflows
// size_t and (a)[SIZE_MAX] reads off the heap).
#define buf_last(a)		((a)[buf_len(a) - 1])

#define buf_push(a,v) (buf__maybegrow(a,1), (a)[buf_get_internals(a)->length++] = (v))

// Pre: buf_len(a) >= 1. NULL or empty buffers are UB (size_t length
// underflows past 0). Return value is the popped element.
#define buf_pop(a) (buf_get_internals(a)->length--, (a)[buf_get_internals(a)->length])

#define buf_erase(a,i)      buf_erasen(a,i,1)
#define buf_erasen(a,i,n)   (memmove(&(a)[i], &(a)[(i)+(n)], sizeof *(a) * (buf_get_internals(a)->length-(n)-(i))), buf_get_internals(a)->length -= (n))

// Set both length and capacity to n. New elements (when growing) are
// uninitialised; pre-existing elements past n are dropped on shrink.
// Use buf_resize_zeroed if you need new elements zeroed.
#define buf_resize(a,n) (buf__sbgrow(a,(n)),buf_get_internals(a)->length=(n))

// Resize and zero any newly-added slots.
#define buf_resize_zeroed(a,n) do { \
	size_t __brz_old = buf_len(a); \
	buf_resize((a), (n)); \
	if ((size_t)(n) > __brz_old) \
		memset(&(a)[__brz_old], 0, sizeof(*(a)) * ((size_t)(n) - __brz_old)); \
} while(0)

// Set capacity to exactly n. Does NOT touch length -- if n < length, the
// buffer ends up with cap < len and subsequent reads/writes past index n
// are out of bounds. Prefer buf_reserve_at_least unless you specifically
// want to lock the capacity to n.
#define buf_reserve(a,n) (*((void **)&(a)) = buf_do_realloc_((a), (n), sizeof(*(a))))

// Grow capacity to at least n; never shrinks. No-op if cap is already >= n.
#define buf_reserve_at_least(a,n) do { \
	if (buf_cap(a) < (size_t)(n)) buf_reserve((a),(n)); \
} while(0)

#define buf_clear(a) ((a) ? (buf_get_internals(a)->length = 0) : 0)

// Set the buffer's logical length without touching capacity. Used for
// in-place compaction (a write-pointer loop fills [0, writeIdx), then
// buf_truncate(a, writeIdx) discards the tail). Caller MUST ensure
// n <= buf_cap(a); growing past capacity via this macro is UB. Use
// buf_resize when you need to grow.
#define buf_truncate(a,n) ((a) ? (buf_get_internals(a)->length = (size_t)(n)) : (size_t)0)

// O(1) unordered remove: swap element i with the last, then pop.
// Pre: buf_len(a) >= 1 and i < buf_len(a).
#define buf_swap_erase(a,i) do { (a)[i] = (a)[buf_get_internals(a)->length - 1]; buf_get_internals(a)->length--; } while(0)

// Insert element v at index i, shifting subsequent elements right.
// Pre: i <= buf_len(a). i == buf_len(a) appends.
#define buf_insert(a,i,v) do { buf__maybegrow(a,1); memmove(&(a)[(i)+1], &(a)[i], sizeof(*(a)) * (buf_get_internals(a)->length - (i))); (a)[i] = (v); buf_get_internals(a)->length++; } while(0)

// Bulk append n elements from pointer ptr. ptr must NOT alias a
// (memcpy semantics; reallocation may move a out from under ptr).
#define buf_push_n(a,ptr,n) do { buf__maybegrow(a,n); memcpy(&(a)[buf_get_internals(a)->length], (ptr), sizeof(*(a)) * (n)); buf_get_internals(a)->length += (n); } while(0)

// Append every element of src onto dst. Both buffers must be the same
// element type. NULL/empty src is a no-op.
#define buf_concat(dst, src) do { \
	size_t __bc_n = buf_len(src); \
	if (__bc_n > 0) buf_push_n((dst), (src), __bc_n); \
} while(0)

/* Private */
#define buf_get_internals(a)							\
	((a) ? ((struct BufHeader *)(buf_containerof(a, struct BufHeader, data))) : (struct BufHeader *)nullptr)

// Deep-copy a stretchy buffer. Returns a fresh buffer with the same
// elements as src; src and the returned buffer are independent.
// NULL/empty src returns NULL (matches the convention that an empty
// stretchy buf is the NULL pointer).
#define buf_clone(src) ((decltype(src))buf_clone_((src), sizeof(*(src))))

// ============================================================================
// Inline implementation (for standalone single-header usage)
// ============================================================================

#ifdef BUF_IMPLEMENTATION

static inline void* buf_do_realloc_(void *a, size_t nr, size_t sz)
{
	struct BufHeader* b;
	Allocator* allocator = a ? buf_get_internals(a)->allocator : allocator_get_cpp_allocator();

	b = (struct BufHeader*)allocator->realloc(buf_get_internals(a),
		sizeof(struct BufHeader) + nr * sz);
	if (!b) { abort(); }
	b->capacity = nr;
	b->allocator = allocator;
	if (!a) b->length = 0;
	return (void *)&(b->data);
}

#endif // BUF_IMPLEMENTATION

#endif // __BUF_H_
