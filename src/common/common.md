# COMMON.md — Shared Modules Reference

**Location**: `common/`
**Role**: Foundation layer. Every pipeline step links against these modules.
Build them first. Test them independently. No step-specific logic lives here.

All modules follow the Hanson pattern (§2 of c-style-guide): one header, one
`.c`, one prefix. Opaque types where appropriate. Every public function asserts
its preconditions. Everything file-private is `static`.

**Compile flags** (all modules):
```
-std=c99 -Wall -Wextra -Wshadow -Wformat=2 -Werror -Wconversion
Debug:   -fsanitize=address,undefined -DVESUVIUS_DEBUG
Release: -O2 -DNDEBUG -fopenmp
```
[c-style-guide §10.5; CLAUDE.md rule 10]

---

## Build Order (dependency DAG)

```
arena.c          (no deps — build first)
    |
except.c         (depends on arena for Except_T defs only)
    |
    +--- mesh_types.h   (header-only, depends on nothing)
    |
    +--- union_find.c   (depends on arena)
    |
    +--- pca.c          (depends on arena, libm)
    |
    +--- tiff_io.c      (depends on arena, libtiff, pthreads)
    |
    +--- csr.c          (depends on arena)
    |
    +--- kdtree.c       (depends on arena)
    |
    +--- bfs.c          (depends on arena, csr)
    |
    |
    +--- obj_io.c       (depends on arena)
```

Link line for a step:
```
gcc -o step0 step0.o arena.o tiff_io.o union_find.o pca.o csr.o \
    -ltiff -lm -fopenmp -lpthread
```

---

## 1. arena.h — Region Allocator

**Prefix**: `Arena_`
**Source**: `arena.c`, `arena.h`
**Dependencies**: `<stdlib.h>`, `<string.h>`, `<assert.h>`, `<pthread.h>` (optional, for thread-safe variant)
**Book refs**: Hanson Ch.6 (arena), Ch.16 (save/restore), Ch.20 (threading); Sanglard Ch.2 (DOOM zone allocator)

### Purpose

One arena per volume. Create at the start of `run_pipeline_for_cube()`. Pass
as explicit `Arena_T arena` parameter to every function that allocates. Dispose
at end. No individual `free()`. No leaks. No use-after-free. This is the single
most important pattern in the pipeline.
[Hanson rule 9; CLAUDE.md rule 1]

### Interface

```c
#ifndef ARENA_INCLUDED
#define ARENA_INCLUDED

typedef struct Arena_T *Arena_T;

/* Lifecycle */
Arena_T Arena_new(void);
void    Arena_dispose(Arena_T *ap);          /* free all + NULL the pointer */
void    Arena_free(Arena_T arena);           /* free all memory, keep arena alive */

/* Allocation — never returns NULL (raises Arena_Failed on OOM) */
void   *Arena_alloc(Arena_T arena, long nbytes, const char *file, int line);
void   *Arena_calloc(Arena_T arena, long count, long size,
                     const char *file, int line);

/* Save/restore checkpoints for scratch memory */
typedef struct { char *avail; struct Arena_chunk *chunk; } Arena_Mark;
Arena_Mark Arena_save(Arena_T arena);
void       Arena_restore(Arena_T arena, Arena_Mark mark);

/* Convenience macros — inject __FILE__, __LINE__ for diagnostics */
#define ARENA_ALLOC(a, nbytes) \
    Arena_alloc((a), (nbytes), __FILE__, __LINE__)
#define ARENA_CALLOC(a, count, size) \
    Arena_calloc((a), (count), (size), __FILE__, __LINE__)
#define ARENA_NEW(a, p) \
    ((p) = Arena_alloc((a), (long)sizeof *(p), __FILE__, __LINE__))

extern const Except_T Arena_Failed;  /* raised on OOM */

#endif
```

### Internals

Linked list of chunks. Each chunk: contiguous block from `malloc`. Arena struct
holds `{avail, limit, chunks}`. Allocation rounds up to 16-byte alignment
(union of `double`, `long`, `void *`, `long double`). If request > remaining,
allocate new chunk of `max(request, 64MB)`. Free list capped at ~10 chunks for
reuse across volumes.
[Hanson rules 11–16; c-style-guide §1.1–1.2]

### Key Rules

1. **Arena is an explicit parameter — never a global.** Every function that
   allocates takes `Arena_T arena` as its first argument. No module-level
   `static Arena_T`.
   [CLAUDE.md rule 1; Hanson rule 9]

2. **Use save/restore for scratch within a step.** Pattern:
   ```c
   Arena_Mark mark = Arena_save(arena);
   /* allocate scratch: flow network, BFS queue, etc. */
   /* ... do work, copy results to persistent storage ... */
   Arena_restore(arena, mark);  /* reclaim scratch */
   ```
   Never call `Arena_restore` after `Arena_free` — the mark is invalid.
   [Hanson rules 62–64; CLAUDE.md rule 2]

3. **No malloc in hot paths.** All vertex arrays, face arrays, BFS queues,
   CSR arrays — arena-allocated before the loop. Zero allocations inside
   marching cubes, Edmonds-Karp, or KD-tree ball query loops.
   [Bentley Ch.10; c-style-guide §1.4.2]

4. **Per-thread arenas for OpenMP.** Shared read-only data (KD-tree, mesh)
   needs no synchronization. Per-thread scratch (neighbor buffer, covariance
   matrix) uses a per-thread arena or stack allocation. Never serialize
   bump-allocate behind a mutex in hot loops.
   [Hanson rules 55, 65, 76; c-style-guide §5.6]

5. **Size initial chunks at 64 MB.** Typical volume working set: 50–200 MB.
   64 MB chunks minimize chunk-list traversal. Arena's `+ request_size` term
   handles rare oversized allocations.
   [Hanson rules 12, 77]

6. **External library allocations bypass the arena.** libtiff, PoissonRecon,
   Triangle, Clipper2 manage their own memory. The arena is for our data.
   [Hanson rule 14]

### Thread Safety

`Arena_T` is **not thread-safe**. Use one per thread. The pattern for OpenMP:
```c
Arena_T thread_arenas[MAX_THREADS];
for (int t = 0; t < n_threads; t++)
    thread_arenas[t] = Arena_new();

#pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
for (size_t i = 0; i < n_components; i++) {
    Arena_T my_arena = thread_arenas[omp_get_thread_num()];
    process_component(my_arena, components[i]);
}
```
[Hanson rule 55; c-style-guide §5.6.5]

---

## 2. except.h — Structured Exception Handling

**Prefix**: `Except_`, `TRY`, `EXCEPT`, `FINALLY`, `RAISE`, `RETURN`
**Source**: `except.c`, `except.h`
**Dependencies**: `<setjmp.h>`, `<stdio.h>`
**Book refs**: Hanson Ch.4 (exceptions)

### Purpose

Wraps `setjmp`/`longjmp` into structured exception handling. Three use cases
in the pipeline:
- `Arena_Failed`: OOM during allocation
- `IO_Failed`: TIFF read/write failure, PLY parse failure
- `Timeout`: per-volume time budget exceeded

### Interface

```c
typedef struct { const char *reason; } Except_T;

extern const Except_T Arena_Failed;
extern const Except_T IO_Failed;
extern const Except_T Timeout;

/* Usage pattern — per-volume dispatch: */
Arena_T arena = Arena_new();
TRY
    run_pipeline(arena, cube_id, out_dir);
EXCEPT(Timeout)
    fprintf(stderr, "cube %s: timeout\n", cube_id);
    copy_raw_prediction(cube_id, out_dir);
EXCEPT(Arena_Failed)
    fprintf(stderr, "cube %s: OOM\n", cube_id);
    copy_raw_prediction(cube_id, out_dir);
FINALLY
    Arena_dispose(&arena);
END_TRY;
```

### Key Rules

1. **Use RETURN, not bare `return`, inside TRY blocks.** A bare `return`
   skips the exception stack pop → corrupts the handler chain. The `RETURN`
   macro pops the stack first.
   [Hanson rule 20]

2. **Exceptions are for rare conditions.** If you write TRY in inner loops,
   the design is wrong. Exceptions appear at: per-volume dispatch, I/O calls,
   external process calls. Algorithmic control flow (no bridge found, oracle
   returns 1) uses return codes, not exceptions.
   [Hanson rules 21–22]

3. **Every exception fires a fallback, never silent corruption.** If bridge
   cut hits OOM or timeout, return the original mesh unchanged. The outer
   loop copies the raw prediction. A wrong topology is worse than no
   post-processing.
   [Hanson rule 25; CLAUDE.md rule 18]

---

## 3. mesh_types.h — Canonical ComponentMesh

**Source**: `mesh_types.h` (header-only)
**Dependencies**: `<stdint.h>`, `<stddef.h>`
**Book refs**: Hanson Ch.14 (self-pointer validation); Elements rule 14

### Purpose

Single struct definition shared by ALL stages. Defined once here. Do not
redefine per-stage. Do not add step-specific fields — use a separate context
struct if a step needs extra per-component data.

### Definition

```c
#ifndef MESH_TYPES_INCLUDED
#define MESH_TYPES_INCLUDED

#include <stdint.h>
#include <stddef.h>

typedef struct ComponentMesh {
    float   *verts;          /* [nv * 3], float32, (z,y,x) order */
    int32_t *faces;          /* [nf * 3], int32, 0-based indices */
    size_t   nv;             /* vertex count */
    size_t   nf;             /* face (triangle) count */
    int      comp_id;        /* 1-based, assigned by Step 4 */
    float    pca_normal[3];  /* unit normal from Step 0 PCA */
    float    centroid[3];    /* mean vertex position from Step 0 */
    void    *self;           /* validation sentinel: self == &this */
} ComponentMesh;

#endif
```

### Contracts

- `verts` and `faces` are arena-allocated flat arrays. Never `free()` them
  individually.
- Coordinate convention: **(z, y, x)** voxel coordinates everywhere.
  Document this in every function that touches vertices.
  [STAGE_0 §6.1; c-style-guide §6.1.1]
- `float32` for vertices (not float64). Sub-microvoxel accuracy with ~7
  decimal digits. Halves cache footprint vs double.
  [CLAUDE.md rule 4; c-style-guide §5.4.6–5.4.7]
- `int32_t` for face indices (not int64). Max 2.1B vertices — a single
  320³ component never exceeds ~2M vertices.
  [CLAUDE.md rule 4]
- `pca_normal` and `centroid` are set by Step 0 and reused by Steps 1–3.
  Sign convention: largest-magnitude component is positive.
  [STAGE_0 Phase 4; PIPELINE_REFERENCE §4 step 6d]
- `self` validation sentinel: set `cm->self = cm` after init. Assert
  `cm->self == cm` on entry to every public function that takes a
  `ComponentMesh *`. Set `cm->self = NULL` on invalidation.
  [Hanson rules 61, 74; c-style-guide §10.3]

### Validation Function

```c
static inline int ComponentMesh_valid(const ComponentMesh *cm) {
    return cm && cm->self == (void *)cm && cm->nv > 0 && cm->verts
           && cm->faces && cm->nf > 0;
}
```

Call this with `assert(ComponentMesh_valid(cm))` at the top of every step
entry point. In debug builds, also verify: all face indices in
`[0, nv - 1]`, no NaN in verts, no degenerate faces.
[CLAUDE.md rule 17; c-style-guide §10.3]

---

## 4. csr.h — Compressed Sparse Row Graph / Matrix

**Prefix**: `CSR_`
**Source**: `csr.c`, `csr.h`
**Dependencies**: `arena.h`
**Book refs**: CLRS Ch.22.1 (adjacency lists); Skiena Ch.7.2, 15.4 (packed arrays); c-style-guide §7.1; Nocedal Ch.5.1 (mat-vec)

### Purpose

Single representation for all graph adjacency and sparse matrices in the
pipeline. Used for: mesh adjacency (Steps 0–5), flow network (Step 2),
Laplacian and cotangent Laplacian (Steps 5–6), UV grid cell lists (Step 3).

### Interface

```c
#ifndef CSR_INCLUDED
#define CSR_INCLUDED

typedef struct CSR_T *CSR_T;

/* Build from mesh faces — undirected adjacency, deduped.
 * Two passes: degree count → prefix sum → fill.
 * O(V + E), no sorting. All arrays arena-allocated. */
CSR_T CSR_from_faces(Arena_T arena, const int32_t *faces, size_t nf,
                     size_t nv);

/* Build from COO (coordinate) triples — for weighted matrices.
 * Expects pre-sorted or will sort internally. */
CSR_T CSR_from_coo(Arena_T arena, const int32_t *rows, const int32_t *cols,
                   const float *vals, size_t nnz, size_t nrows);

/* Access: neighbors of vertex u are target[offset[u] .. offset[u+1]).
 * Direct field access is permitted for performance-critical loops. */
int32_t  CSR_nrows(const CSR_T csr);
int32_t  CSR_nnz(const CSR_T csr);
const int32_t *CSR_offset(const CSR_T csr);   /* [nrows + 1] */
const int32_t *CSR_target(const CSR_T csr);    /* [nnz] */
const float   *CSR_weight(const CSR_T csr);    /* [nnz] or NULL */

/* Sparse mat-vec: y = A @ x (3-channel: x and y are [N*3]) */
void CSR_matvec3(const CSR_T csr, const float *x, float *y);

/* Sparse mat-vec: y = A @ x (single channel: x and y are [N]) */
void CSR_matvec(const CSR_T csr, const float *x, float *y);

/* Validation (debug builds) */
void CSR_validate(const CSR_T csr);

#endif
```

### Construction (Two-Pass)

```
Pass 1: For each face (a,b,c), increment degree[a], degree[b], degree[c]
        by 2 each (two neighbors per face vertex). This overcounts duplicates.
        Compute prefix sum: offset[0] = 0, offset[v] = offset[v-1] + degree[v-1].

Pass 2: Reset a cursor array = copy of offset. For each face (a,b,c):
        Insert b,c into a's list; a,c into b's list; a,b into c's list.
        Use cursor[v]++ as write position.

Dedup:  Sort each row (offset[v]..offset[v+1]) and remove duplicates.
        Compact target[] in-place, update offset[].
```
[CLRS Ch.22.1; c-style-guide §7.1.3–7.1.4]

### Key Rules

1. **CSR is immutable after construction.** Never insert/delete edges. If
   mesh topology changes (bridge cut, split), build a new CSR.
   [CLRS Ch.22.1; c-style-guide §7.1.2]

2. **Arena-allocate all three arrays in one shot.** `offset`, `target`,
   `weight` are allocated contiguously for cache locality.
   [c-style-guide §7.1.3]

3. **For the flow network, pre-allocate forward + reverse edges.** Each edge
   stores the index of its reverse via `reverse[E]` for O(1) updates during
   Edmonds-Karp augmentation.
   [Skiena Ch.8.5; CLRS Ch.26.2; c-style-guide §7.1.5]

4. **Uniform Laplacian = CSR with weight = 1/degree.** Build adjacency CSR,
   then normalize: each row's weights are `1.0f / (float)(offset[v+1] - offset[v])`.
   The mat-vec `L @ verts` computes the average of each vertex's neighbors.
   [STAGE_5 §9.2; c-style-guide §8.1.3]

5. **CSR_validate checks**: `offset[0] == 0`, `offset` is monotonically
   non-decreasing, `offset[nrows] == nnz`, all `target[i]` in `[0, nrows-1]`.
   Run in debug builds after every construction.
   [CLAUDE.md rule 17]

### Users

| Step | CSR usage |
|------|-----------|
| 0 | Not used (Union-Find for mesh CC instead) |
| 2 | Mesh adjacency for k-ring BFS; flow network CSR |
| 3 | Mesh adjacency for BFS flood fill; UV grid cell lists |
| 5 | Laplacian (snap-back mat-vec); cotangent Laplacian (CG solve) |
| 6 | Not used |

---

## 5. kdtree.h — Flat-Array 3D KD-Tree

**Prefix**: `KDTree_`
**Source**: `kdtree.c`, `kdtree.h`
**Dependencies**: `arena.h`
**Book refs**: Skiena Ch.15.6 (KD-tree); Sanglard Ch.5 (DOOM blockmap alternative); c-style-guide §6.4; Abrash Ch.15–17 (cache)

### Purpose

All spatial proximity queries: ball query (Step 3 local PCA, r=10), nearest-
neighbor (Step 2 gap exclusion, Step 5 snap-back), point location (optional).
Build once, query many times. This is the performance-critical data structure
for Step 3.

### Interface

```c
#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

typedef struct KDTree_T *KDTree_T;

/* Build balanced KD-tree from N points in R^3.
 * Points are copied into internal storage (sorted by tree order).
 * Construction: O(N log N). Arena-allocated. */
KDTree_T KDTree_new(Arena_T arena, const float *points, size_t n);

/* 1-NN query. Returns index of nearest point.
 * *out_dist_sq receives squared distance. */
size_t KDTree_nearest(const KDTree_T tree, const float query[3],
                      float *out_dist_sq);

/* Ball query. Returns count of points within squared radius.
 * Writes up to max_results indices into out_indices.
 * If more than max_results points fall within radius, only
 * the first max_results found are returned (no ordering guarantee). */
size_t KDTree_ball_query(const KDTree_T tree, const float center[3],
                         float radius_sq,
                         int32_t *out_indices, size_t max_results);

#endif
```

### Internal Layout

Flat array, heap layout: children of node `i` at `2*i+1` and `2*i+2`. No
pointers. N nodes for N points. Each node stores:

```c
struct KDNode {
    float split_value;     /* coordinate value at split plane */
    int32_t point_index;   /* index into original point array */
    uint8_t split_axis;    /* 0=z, 1=y, 2=x */
    uint8_t is_leaf;       /* 1 if no children */
    /* pad to 12 or 16 bytes for alignment */
};
/* Alternatively, for minimal size: */
struct KDNode {
    float   point[3];      /* copied point coordinates */
    int32_t orig_index;    /* original index for result mapping */
    /* split_axis = depth % 3, or stored as 2-bit field */
};  /* 16 bytes → 4 nodes per 64-byte cacheline */
```

Construction: median-of-three partitioning. Cycle split axes (z→y→x) or
cut along largest-extent dimension. The point array is reordered by
construction order — the original indices are preserved in `orig_index`.
[Skiena rule 56; Sanglard rule 28; c-style-guide §6.4.2–6.4.4]

### Key Rules

1. **Compare squared distances. Never compute `sqrt()` in queries.**
   `dx*dx + dy*dy + dz*dz < radius_sq`. One fewer transcendental per
   comparison.
   [Skiena Ch.15.6; Bentley Ch.9; c-style-guide §6.4.5]

2. **Prune subtrees by bounding-box vs. query-sphere test.** If the node's
   bounding box is entirely outside the query sphere, skip the subtree.
   For ball query with fixed r=10, `radius_sq = 100.0f`.
   [Skiena rule 54; c-style-guide §6.4.5]

3. **Build once, query N times.** Construction O(N log N). Each query
   O(√N) expected for 3D. Total for Step 3: O(N^{3/2}).
   [Skiena rule 55; c-style-guide §6.4.6]

4. **KD-tree is read-only after construction.** Safe for concurrent queries
   from multiple OpenMP threads without synchronization. Per-thread scratch
   (result buffer) is thread-local.
   [Hanson rule 65; c-style-guide §5.6]

5. **Keep node struct small.** 16 bytes per node → 4 nodes per 64-byte
   cacheline. For N=200K: 3.2 MB, fits L3 on both platforms.
   [Sanglard rule 29; c-style-guide §6.4.3]

6. **Consider the grid alternative.** For roughly uniform vertex density, a
   3D uniform grid (cell size = query radius) converts each ball query into
   visiting 27 cells. Implementation as CSR: `cell_offsets[ncells+1]` +
   `vertex_indices[N]`. May be 8× faster than KD-tree for uniform data.
   Profile both. Start with KD-tree (proven correct in Python reference),
   swap to grid if profiling shows the ball query is slower than estimated.
   [Sanglard rules 33–35; c-style-guide §6.4.7–6.4.8]

### Users

| Step | Query type | Typical N |
|------|------------|-----------|
| 2 | 1-NN (gap exclusion distance) | ~1K boundary verts |
| 3 | Ball query r=10 (local PCA normals) | 50K–200K |
| 5 | 1-NN (snap-back target, 20 iterations) | 50K |

---

## 6. union_find.h — Path-Splitting Union-Find

**Prefix**: `UF_` or `uf_`
**Source**: `union_find.c`, `union_find.h`
**Dependencies**: `arena.h`
**Book refs**: CLRS Ch.21.3; Skiena Ch.8.1.3; c-style-guide §7.3

### Purpose

Mesh connected components when edges arrive incrementally (iterate faces,
union edge endpoints). Used in Step 0 (mesh cleanup after backface cull),
Step 3 (small-component cleanup after label assignment), Step 5 (fragment
removal after Poisson).

### Interface

```c
#ifndef UNION_FIND_INCLUDED
#define UNION_FIND_INCLUDED

typedef struct {
    int32_t *parent;   /* parent[i] = parent of i, or i if root */
    int32_t *rank;     /* rank[i] = upper bound on height */
    int32_t  count;    /* current number of distinct components */
    int32_t  n;        /* total elements */
} UnionFind;

/* Create UF for n elements. Arena-allocated. Each element starts as its
 * own component (parent[i] = i, rank[i] = 0, count = n). */
UnionFind UF_new(Arena_T arena, int32_t n);

/* Find with path splitting (iterative, no recursion). */
int32_t uf_find(UnionFind *uf, int32_t x);

/* Union by rank. Decrements uf->count on merge. */
void uf_union(UnionFind *uf, int32_t a, int32_t b);

#endif
```

### Implementation

```c
int32_t uf_find(UnionFind *uf, int32_t x) {
    assert(x >= 0 && x < uf->n);
    while (uf->parent[x] != x) {
        uf->parent[x] = uf->parent[uf->parent[x]];  /* path splitting */
        x = uf->parent[x];
    }
    return x;
}

void uf_union(UnionFind *uf, int32_t a, int32_t b) {
    a = uf_find(uf, a);
    b = uf_find(uf, b);
    if (a == b) return;
    if (uf->rank[a] < uf->rank[b]) { int32_t t = a; a = b; b = t; }
    uf->parent[b] = a;
    if (uf->rank[a] == uf->rank[b]) uf->rank[a]++;
    uf->count--;
}
```
Path splitting (two-pointer) instead of recursive path compression to avoid
stack depth on large meshes.
[CLRS Ch.21.3; c-style-guide §7.3.2]

### Key Rules

1. **Amortized O(α(N)) per op** — effectively constant. No need for
   anything fancier.
   [CLRS Ch.21.3; Skiena Ch.8.1.3]

2. **Use for mesh CC when edges arrive incrementally.** Iterate faces:
   for each face `(a,b,c)`, call `uf_union(a,b); uf_union(b,c);`. After
   all faces, `uf_find(v)` gives the component root. Group faces by root.
   [CLRS Ch.23.2 Kruskal pattern; c-style-guide §7.3.3]

3. **Size the parent/rank arrays to the original vertex count.** After
   backface cull, vertex indices may have gaps. That's fine — unused indices
   just remain as singleton components. Compact vertices AFTER Union-Find
   grouping, not before.
   [STAGE_0 Pitfalls: "Union-Find with Non-Contiguous Vertex Indices"]

---

## 7. bfs.h — BFS Utilities

**Prefix**: `BFS_`
**Source**: `bfs.c`, `bfs.h`
**Dependencies**: `arena.h`, `csr.h`
**Book refs**: CLRS Ch.22.2; Skiena Ch.7.6–7.8; c-style-guide §7.2

### Purpose

Core graph traversal. Three query patterns share the same BFS engine with
different termination/action callbacks:

1. **k-ring expansion** (Step 2): expand from seed vertex k hops → return
   all vertices within k geodesic hops.
2. **Augmenting path** (Step 2): BFS on residual flow graph, stop when sink
   is reached, trace path via `parent_edge[]`.
3. **Flood fill** (Step 3): multi-source BFS on mesh adjacency restricted to
   same-label edges.

### Interface

```c
#ifndef BFS_INCLUDED
#define BFS_INCLUDED

/* k-ring BFS: expand from seed by k hops on CSR graph.
 * Writes visited vertex indices into out_vertices (arena-allocated).
 * Returns count of visited vertices. */
size_t BFS_kring(const CSR_T graph, int32_t seed, int32_t k,
                 int32_t *out_vertices, size_t max_out,
                 uint32_t *dist);  /* [nrows], caller-provided scratch */

/* Multi-source BFS: expand from all seeds simultaneously.
 * seeds[0..n_seeds-1] are starting vertices. dist[seed] = 0.
 * Fills dist[] for all reachable vertices. */
void BFS_multi_source(const CSR_T graph, const int32_t *seeds,
                      size_t n_seeds, uint32_t *dist);

/* Edmonds-Karp BFS: find shortest path from source to sink in residual
 * graph. Returns 1 if path found, 0 if sink unreachable.
 * parent_edge[v] = edge index used to reach v (-1 = unvisited).
 * The caller traces the path from sink to source via parent_edge[]. */
int BFS_augmenting_path(const int32_t *offset, const int32_t *target,
                        const int32_t *cap, const int32_t *flow,
                        int32_t n_nodes, int32_t source, int32_t sink,
                        int32_t *parent_edge,  /* [n_nodes] scratch */
                        int32_t *queue);       /* [n_nodes] scratch */

#endif
```

### Key Rules

1. **Flat queue with head/tail indices.** Never a linked list. Each vertex
   enqueued at most once → queue never exceeds V entries. Zero allocation
   during execution. Arena-allocate once, reuse across BFS invocations.
   [CLRS Ch.22.2; c-style-guide §7.2.1]

2. **Eliminate the color array.** `dist[v] == UINT32_MAX` means unvisited.
   Saves one array.
   [CLRS Ch.22.2; c-style-guide §7.2.2]

3. **Reset with `memset(dist, 0xFF, n * sizeof(uint32_t))`** to set all
   entries to `UINT32_MAX` (0xFFFFFFFF). This is correct because `0xFF`
   fills every byte.

4. **Customize via early termination, not callbacks.** The BFS core is
   written once. k-ring BFS stops when `dist[v] > k`. Edmonds-Karp BFS
   stops when `dist[sink] != UINT32_MAX`. Flood fill checks label equality.
   These are inlined if-checks, not function-pointer callbacks.
   [Skiena Ch.7.6; c-style-guide §7.2.4]

---

## 8. pca.h — 3×3 PCA Normal Computation

**Prefix**: `PCA_`
**Source**: `pca.c`, `pca.h`
**Dependencies**: `arena.h`, `<math.h>`
**Book refs**: Nocedal & Wright (numerical stability); PIPELINE_REFERENCE §5–7; c-style-guide §9.4

### Purpose

Compute the sheet-normal direction of a vertex cloud. The smallest eigenvector
of the 3×3 covariance matrix gives the direction perpendicular to the sheet
surface. Used by Steps 0, 1, 2, 3 (global PCA) and Step 3 (per-vertex local
PCA on neighbor sets).

### Interface

```c
#ifndef PCA_INCLUDED
#define PCA_INCLUDED

/* Compute PCA normal of N points.
 * out_normal receives the unit eigenvector of the smallest eigenvalue.
 * out_centroid receives the mean position (may be NULL).
 * Sign convention: largest-magnitude component of out_normal is positive.
 *
 * Uses double internally for the covariance accumulation, converts
 * result to float. This is critical for numerical stability when
 * N > 100K vertices.
 *
 * Returns 0 on success, -1 if N < 3 (degenerate). */
int PCA_normal(const float *points, size_t n,
               float out_normal[3], float out_centroid[3]);

/* Project N points onto a direction vector. out_proj[i] = dot(point[i], dir).
 * Returns (min_proj, max_proj) via out pointers. */
void PCA_project(const float *points, size_t n, const float dir[3],
                 float *out_proj, float *out_min, float *out_max);

/* Build orthonormal UV basis perpendicular to a given normal.
 * Uses Hughes-Möller technique for robustness near axis-aligned normals. */
void PCA_orthonormal_basis(const float normal[3],
                           float out_u[3], float out_v[3]);

#endif
```

### Key Rules

1. **Use `double` for covariance accumulation.** The 3×3 covariance matrix
   sums `(v[i] - mean) * (v[j] - mean)` over 100K+ vertices. Float32
   accumulation introduces catastrophic cancellation. Accumulate in double,
   convert the final eigenvector to float.
   [STAGE_0 Pitfalls: "Floating-Point PCA on Large Vertex Sets"; c-style-guide §9.4]

2. **3×3 symmetric eigensolver: Jacobi rotation or closed-form cubic.**
   Do not link LAPACK for a 3×3 problem. Jacobi: 3–5 iterations for full
   convergence on a symmetric matrix. Closed-form: Cardano's formula for
   the characteristic polynomial. Either works. Jacobi is more numerically
   stable for near-degenerate cases.
   [STAGE_0 Phase 4; Nocedal — numerical hygiene]

3. **Sign convention: make the largest-magnitude component positive.**
   The eigenvector's sign is arbitrary. After computing it, find the
   component with the largest absolute value, and negate the whole vector
   if that component is negative. This must be consistent across ALL steps.
   [PIPELINE_REFERENCE §4 step 6d; STAGE_0 Pitfalls: "PCA Sign Convention"]

4. **Reuse PCA results across steps.** Step 0 computes `pca_normal` and
   `centroid` into `ComponentMesh`. Steps 1, 2, 3 read them. Do not
   recompute unless the mesh has been modified (bridge cut, split).
   [c-style-guide §6.6.3; Abrash Ch.50]

---

## 9. tiff_io.h — libtiff Wrapper

**Prefix**: `TiffIO_`
**Source**: `tiff_io.c`, `tiff_io.h`
**Dependencies**: `arena.h`, `<tiff.h>`, `<tiffio.h>`, `<pthread.h>`
**Book refs**: Hanson rule 67; c-style-guide §11.2; STAGE_0 §11.1

### Purpose

Load/save multi-page TIFF volumes. Wraps libtiff behind a mutex (libtiff is
NOT thread-safe). Single entry/exit points with error propagation.

### Interface

```c
#ifndef TIFF_IO_INCLUDED
#define TIFF_IO_INCLUDED

/* Load multi-page TIFF as flat uint8 volume.
 * Arena-allocates *out_vol of size D*H*W.
 * Returns 0 on success, -1 on failure.
 * Page 0 = z=0, page 1 = z=1, etc. */
int TiffIO_load(Arena_T arena, const char *path,
                uint8_t **out_vol, int *out_D, int *out_H, int *out_W);

/* Save flat uint8 volume as multi-page TIFF.
 * Returns 0 on success, -1 on failure. */
int TiffIO_save(const char *path,
                const uint8_t *vol, int D, int H, int W);

#endif
```

### Key Rules

1. **Mutex-serialize all libtiff calls.**
   ```c
   static pthread_mutex_t tiff_lock = PTHREAD_MUTEX_INITIALIZER;
   ```
   Lock before `TIFFOpen`, unlock after `TIFFClose`. Even if only one thread
   calls TIFF I/O, the mutex protects against future refactoring.
   [Hanson rule 67; STAGE_0 §11.2]

2. **Check every return value.** `TIFFOpen` returns NULL on failure.
   `TIFFReadEncodedStrip` returns -1 on failure. Propagate errors to the
   outer loop for fallback to raw prediction.
   [c-style-guide §2.6; CLAUDE.md rule 18]

3. **Verify page ordering.** The pipeline assumes page 0 = z=0. Some TIFF
   writers produce pages in reverse order. Check `TIFFCurrentDirectory()`
   after each page read.
   [STAGE_0 Pitfalls: "libtiff Page Ordering"]

4. **Use strip-based API.** `TIFFReadEncodedStrip()` is more efficient than
   `TIFFReadScanline()` for bulk reads. Not a bottleneck (~0.1–0.3s for
   32.8 MB), but no reason to be slow.
   [STAGE_0 §11.1]

5. **Zero-mesh case.** If zero meshes survive processing, save an all-zero
   volume. The TIFF must still be valid (all-zero pages). Must not crash.
   [STAGE_6 §9.3.2]

---

## 11. obj_io.h — OBJ I/O (Debug Only)

**Prefix**: `ObjIO_`
**Source**: `obj_io.c`, `obj_io.h`
**Dependencies**: none (uses `<stdio.h>` directly)

### Purpose

Write OBJ files for debugging / visual inspection in MeshLab. Read OBJ files
for test fixtures. Not performance-critical — used only in debug/test paths.

### Interface

```c
/* Write mesh as text OBJ.
 * Vertices as "v z y x\n" (note: OBJ convention may need reorder).
 * Faces as "f a+1 b+1 c+1\n" (OBJ is 1-indexed). */
int ObjIO_write(const char *path, const float *verts, size_t nv,
                const int32_t *faces, size_t nf);

/* Read mesh from text OBJ. Arena-allocates output arrays.
 * Returns 0 on success. */
int ObjIO_read(Arena_T arena, const char *path,
               float **out_verts, size_t *out_nv,
               int32_t **out_faces, size_t *out_nf);
```

---

## 12. Pipeline Constants

Shared constants referenced by multiple stages. Define in a single header
`pipeline_constants.h`:

```c
#ifndef PIPELINE_CONSTANTS_INCLUDED
#define PIPELINE_CONSTANTS_INCLUDED

/* Step 0 */
#define MIN_CC_SIZE          500    /* voxels — discard smaller components */
#define MAX_COMPONENTS        20    /* keep top N by size */
#define CLEANUP_MICRO_HOLE_MAX 6   /* max boundary loop verts for fill */
#define MIN_FRAGMENT_FACES   100    /* discard sub-components smaller */

/* Step 1 */
#define ORACLE_UV_GRID_SIZE  320    /* UV grid resolution for raycasting */
#define ORACLE_MIN_GAP        10    /* voxel gap to count as sheet boundary */

/* Step 2 */
#define SEED_RING             20    /* BFS expansion hops for source/sink */
#define MAX_FLOW_LIMIT       100    /* Edmonds-Karp early exit */
#define CUT_GAP_DEPTH       7.0f   /* exclusion zone around cut (voxels) */
#define BRIDGE_MAX_DEPTH      10    /* max recursion depth */
#define FLOW_INF       (1 << 30)   /* ~1 billion — never INT_MAX */

/* Step 3 */
#define BALL_QUERY_RADIUS  10.0f   /* KD-tree ball query radius */
#define BALL_QUERY_R_SQ   100.0f   /* squared, for comparison */
#define GAP_THRESHOLD      8.0f    /* w-gap between sheets */
#define CUT_EXPANSION_HOPS     6   /* BFS expansion at label boundaries */
#define MIN_FRAGMENT_VERTS_S3 30   /* percent threshold for small-comp merge */

/* Step 5 */
#define MIN_POISSON_VERTS   1000   /* skip Poisson below this */
#define SNAP_ITERATIONS       20   /* Laplacian snap-back iterations */
#define SNAP_ALPHA          0.3f   /* blend factor (target weight) */
#define DENSITY_TRIM_PCT    0.05f  /* 5th percentile density trim */
#define MIN_FRAGMENT_VERTS_S5  50  /* post-Poisson fragment removal */

/* Step 6 */
#define MIN_BARY_SUBDIV        6   /* minimum barycentric samples per edge */

/* Threading */
#define PARALLEL_VERTEX_CUTOFF 1000  /* don't parallelize below this */
#define PARALLEL_FACE_CUTOFF   2000

#endif
```

---

## 13. Cross-Cutting Concerns

### Volume Indexing Convention

**Everywhere**: `vol[z * H * W + y * W + x]`. Always `size_t` for the
product `(size_t)z * H * W + (size_t)y * W + (size_t)x`. Never compute
`z * H * W` as `int` — 320×320×320 = 32.8M fits int32 but leaves no headroom.
[CLAUDE.md rule 19; c-style-guide §6.1.1; STAGE_0 §9.5]

### Neighbor Offset Tables

Static arrays for 6-connectivity and 26-connectivity:
```c
static const int NBR6[6][3] = {
    {-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}
};
static const int NBR26[26][3] = { /* all (dz,dy,dx) except (0,0,0) */ };
```
Bounds-check every neighbor access: `if (nz >= 0 && nz < D && ...)`.
[Elements Ch.2; c-style-guide §6.1.4]

### Float Discipline

- Never compare floats for exact equality. Use `fabsf(a - b) < eps`.
- Never use float loop counters. Integer for counting/indexing, float for
  geometry.
- Guard divisions: `if (denom > 1e-12f)` before `x / denom`.
- Use `double` for: PCA covariance accumulation, CG/PCG solver internals,
  any summation over 10K+ terms.
[CLAUDE.md rule 20; c-style-guide §9.4]

### Error Handling Hierarchy

```
assert()           — contract violations (bugs). Abort immediately.
                     Face index out of range, NULL arena, self != self.

return -1          — runtime failures. Caller handles.
                     Zero components, degenerate mesh, can't open file.

RAISE(Except_T)    — resource exhaustion, timeouts. Outer loop catches.
                     OOM, libtiff failure, PoissonRecon timeout.
```
[Hanson rules 17–22; CLAUDE.md rules 7, 18; c-style-guide §2.4–2.6]

### Validation Functions (Debug Builds)

Every module provides a `_validate` function. Call after construction in
debug builds (`#ifndef NDEBUG`). Strip in release (`-DNDEBUG`).

| Module | Validation |
|--------|------------|
| `ComponentMesh` | `self == self`, face indices in range, no NaN |
| `CSR` | offset monotonic, offset[0]==0, target in range |
| `KDTree` | all points present, tree balanced within 1 level |
| `UnionFind` | parent[i] in range, root is self-referential |
| `FlowNet` | flow conservation, capacity constraint, reverse consistency |

[CLAUDE.md rule 17; c-style-guide §10.3]

---

## 14. Memory Budget Summary

| Data structure | Size formula | Typical (N=100K verts, 320³ vol) |
|----------------|-------------|----------------------------------|
| Volume (uint8) | D×H×W | 32.8 MB |
| Labels (int32) | D×H×W×4 | 131 MB |
| Ownership (uint16) | D×H×W×2 | 65.5 MB |
| Vertices (float×3) | N×12 | 1.2 MB |
| Faces (int32×3) | M×12 | 1.2 MB (M≈N) |
| CSR adjacency | ~(N + 6N)×4 | 2.8 MB |
| KD-tree nodes | N×16 | 1.6 MB |
| Flow network | ~82 MB (see STAGE_2) | 82 MB (N=200K) |
| BFS queue+dist | 2×D×H×W×4 or 2×V×4 | 263 MB (voxel) or 1.6 MB (mesh) |

**Kaggle RAM budget**: ~13 GB usable. Worst case: volume arrays (~230 MB) +
flow network (~82 MB) + misc (~50 MB) ≈ 362 MB per volume. Comfortable.

---

## 15. Makefile Sketch

```makefile
CC       = gcc
CFLAGS   = -std=c99 -Wall -Wextra -Wshadow -Wformat=2 -Wconversion
LDFLAGS  = -lm -ltiff -fopenmp -lpthread

# Debug
DEBUG_FLAGS = -g -O0 -fsanitize=address,undefined -DVESUVIUS_DEBUG
# Release
RELEASE_FLAGS = -O2 -DNDEBUG -fopenmp

COMMON_SRC = common/arena.c common/except.c common/csr.c common/kdtree.c \
             common/union_find.c common/bfs.c common/pca.c \
             common/tiff_io.c common/obj_io.c
COMMON_OBJ = $(COMMON_SRC:.c=.o)

# External deps
EXT_INC  = -Ideps/include
EXT_LIB  = -Ldeps/lib -ltriangle -lClipper2 -lClipper2Z -lstdc++

all: pipeline

common/%.o: common/%.c common/%.h
	$(CC) $(CFLAGS) $(RELEASE_FLAGS) -c $< -o $@

pipeline: main.c $(COMMON_OBJ) step0.o step1.o step2.o step3.o step4.o step5.o step6.o
	$(CC) $(CFLAGS) $(RELEASE_FLAGS) $^ -o $@ $(LDFLAGS) $(EXT_LIB)

test: test_common test_step0 test_step2 test_step3

test_common: tests/test_common.c $(COMMON_OBJ)
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $^ -o $@ $(LDFLAGS)
	./$@

clean:
	rm -f common/*.o *.o pipeline test_common
```
