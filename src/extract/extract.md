# STAGE 0 — Mesh Extraction

**Module prefix**: `MeshExtract_`
**Source files**: `mesh_extract.c`, `mesh_extract.h`
**Dependencies**: `arena.c`, `tiff_io.c`, `union_find.c`, `marching_cubes.c`

---

## Purpose

Step 0 converts the raw nnUNet binary prediction TIFF into a list of triangle
meshes — one per connected component. This is the foundation of the entire
pipeline: every subsequent step operates on meshes, not voxels. If Step 0
produces bad meshes (disconnected fragments, missing faces, wrong normals,
merged components), every downstream step inherits the error.

**Input**: Multi-page TIFF (`uint8 [D×H×W]`, values 0 or 1).
**Output**: Array of `(component_id, vertices[N×3], faces[M×3])` tuples, plus
per-component PCA normal vectors. Optionally saved as OBJ files for debugging.

**Why it matters for the metric**:
- The connected component labeling determines which voxels become which mesh.
  A labeling error that merges two wraps into one component means the bridge
  cutter (Step 2) must fix what should have been separate from the start.
- The backface cull removes the "wrong side" of each surface. If PCA sign
  convention is inconsistent, the cull keeps the wrong half and downstream
  steps see an inverted surface.
- Micro-hole fill ensures clean manifold boundaries for Poisson reconstruction
  (Step 5). Unfilled holes cause PoissonRecon to produce bloated, incorrect
  surfaces.
- VOI is computed on 26-connected components of the final output. If mesh
  extraction merges fragments that should be separate, VOI score degrades.

---

## Algorithm

### Phase 1: Load and Threshold
1. Load multi-page TIFF via `TiffIO_load()` into a flat `uint8 *vol` of size
   `D * H * W`, arena-allocated. [PIPELINE_REFERENCE §4, step 1–2]
2. Threshold: `binary[i] = (vol[i] > 0) ? 1 : 0`. In practice, nnUNet output
   is already binary; the threshold guards against intermediate values.

### Phase 2: 3D Connected Components (6-connectivity)
3. Allocate `int32_t *labels` of size `D*H*W` from the arena, zero-initialized
   via `Arena_calloc`.
4. BFS flood fill: iterate all voxels in linear order. For each unvisited
   foreground voxel, start BFS with a flat arena-allocated queue. Mark each
   visited voxel with the current component ID. Use 6-connectivity
   (face-adjacent only: ±z, ±y, ±x). Count voxels per component.
   [CLRS Ch.22.2; Skiena Ch.7.7; c-style-guide §6.2.1]
5. Filter: discard components with fewer than `MIN_CC_SIZE` (500) voxels.
   [PIPELINE_REFERENCE §4, step 4]
6. Sort remaining components by voxel count descending. Keep top
   `MAX_COMPONENTS` (20). [PIPELINE_REFERENCE §4, step 5]

### Phase 3: Per-Component Marching Cubes
For each surviving component (up to 20):

7. **Extract bounding box** of the component from the label volume. Pad by 1
   voxel in each direction (clamp to volume bounds). Allocate a local
   `uint8 *padded` subvolume of size `(dz+2) × (dy+2) × (dx+2)` from the
   arena. Copy component voxels into it. The padding ensures marching cubes
   can produce triangles at the boundary.
   [PIPELINE_REFERENCE §4, step 6a–6b; c-style-guide §6.3.2]

8. **Run marching cubes** on the padded subvolume with `level=0.5`. Use the
   Lorensen & Cline lookup table (256 entries, 16 cases with rotational
   symmetry). For each cube of 8 voxels, look up the triangle configuration,
   interpolate edge vertices, emit triangles. Output: `float *verts` (N×3)
   and `int32_t *faces` (M×3), 0-indexed.
   [PIPELINE_REFERENCE §4, step 6b; c-style-guide §6.3.1]

9. **Undo padding offset**: subtract 1 from each vertex coordinate, then add
   the bounding-box origin to get back to global voxel coordinates (z, y, x).
   [PIPELINE_REFERENCE §4, step 6c]

### Phase 4: Backface Cull
10. **Compute PCA normal** of the component's vertex cloud. Build the 3×3
    covariance matrix of vertex positions (subtract mean first). Find the
    eigenvector corresponding to the smallest eigenvalue — this is the sheet
    normal direction. Use a direct 3×3 symmetric eigensolver (Jacobi rotation
    or the closed-form cubic formula), not a general eigenvalue library.
    **Sign convention**: flip the eigenvector so its largest-magnitude component
    is positive. [PIPELINE_REFERENCE §4, step 6d]

11. **Per-face normal**: for each triangle `(a, b, c)`, compute
    `n = (verts[b] - verts[a]) × (verts[c] - verts[a])`. No normalization
    needed — only the sign of `dot(n, pca_normal)` matters.
    [Abrash Ch.51 — compute minimum for a binary decision]

12. **Keep front faces**: retain face `i` iff `dot(face_normal[i], pca_normal) > 0`.
    This removes the "back side" of the sheet, approximately halving the mesh.
    Compact the face array in-place. [PIPELINE_REFERENCE §4, step 6d]

### Phase 5: Mesh Cleanup (connected components + micro-hole fill)
13. **Mesh connected components via Union-Find**: iterate faces; for each
    face's three edges `(a,b)`, `(b,c)`, `(a,c)`, call `uf_union(a, b)` etc.
    After all faces, `uf_find(v)` gives the component representative for each
    vertex. Group faces by their vertices' root. Discard sub-components with
    fewer than 100 faces.
    [CLRS Ch.21.3; Skiena Ch.8.1.3; c-style-guide §7.3.3]

    **Why Union-Find instead of BFS**: the face list gives us edges
    incrementally. Union-Find naturally groups vertices into components during
    a single pass over the faces — no need to build an explicit adjacency graph
    first. This replaces Python's `trimesh.split()` — see PIPELINE_REFERENCE
    §4 Cleanup for the Python approach.
    [Skiena §18.1; CLRS Ch.23.2 (Kruskal pattern)]

14. **Clean unreferenced vertices**: after face removal, some vertices may be
    orphaned. Build a `used[N]` boolean array, mark vertices referenced by
    surviving faces, compact the vertex array, and remap face indices.
    This is a standard reindex pass.

15. **Micro-hole fill**: find boundary edges (edges appearing in exactly 1
    triangle). Chain into closed loops. For loops with ≤ `CLEANUP_MICRO_HOLE_MAX`
    (6) vertices:
    - Loops with ≤ 4 verts: fan-triangulate from the first vertex.
    - Loops with 5–6 verts: insert centroid, connect each boundary edge to centroid.
    - Match orientation to the adjacent face's half-edge direction.
    [PIPELINE_REFERENCE §4, Micro-hole filling algorithm]

    **Boundary edge detection**: build a hash table or sorted edge list. For
    each face `(a,b,c)`, insert edges `(min(a,b), max(a,b))` etc. Count
    occurrences. Edges with count == 1 are boundary. Count == 2 means interior.
    Count > 2 means non-manifold — log a warning but continue.
    [PIPELINE_REFERENCE §14 — boundary detection pattern]

16. **Output**: for each component, store `(comp_id, float *verts, int32_t *faces,
    size_t n_verts, size_t n_faces, float pca_normal[3])` in an arena-allocated
    array of component descriptors. Optionally write OBJ files for debugging.

---

## Book References

### Connected Components (Phase 2)
- **CLRS Ch.22.2** (BFS): the BFS algorithm for graph traversal, O(V+E).
  Adapt for voxel grid: vertices are voxels, edges are 6-connectivity neighbors.
  Use `dist == UINT32_MAX` as unvisited sentinel to eliminate the color array.
- **Skiena Ch.7.7.1**: connected components via repeated BFS from unvisited
  vertices. Clean implementation pattern: loop over all vertices, BFS from
  each undiscovered one, label all discovered with current component number.
- **Skiena Ch.7.6.1**: customize traversal via callbacks — the same BFS core
  serves connected components, k-ring expansion, and Edmonds-Karp.
- **CLRS Ch.21.3**: Union-Find with union-by-rank and path compression,
  O(α(N)) amortized. Used for mesh-level connected components in Phase 5.
- **Skiena Ch.8.1.3**: Union-Find for incremental component tracking. Iterate
  faces, UNION edge endpoints. Avoids building an explicit adjacency graph.

### Marching Cubes (Phase 3)
- **Bentley Ch.7**: feasibility check arithmetic. 320³ voxels × ~50 ns/voxel
  = ~1.6s on Kaggle. Fits easily. No need for a faster isosurface algorithm.
- **Lorensen & Cline (1987)**: the original marching cubes paper. The lookup
  table (256 entries) is public domain and widely available. scikit-image's
  implementation is a reference.

### Backface Cull (Phase 4)
- **Abrash Ch.51**: compute the minimum quantity needed for a binary decision.
  For backface test, only the sign of the Z-component of the cross product
  matters (in our case, the dot product with the PCA normal). Compute one
  dot product per face, not a full normalized normal vector.
- **Abrash Ch.51**: establish conventions (consistent winding order) that make
  per-element decisions O(1). Define winding once at marching cubes time;
  every subsequent step inherits it.

### PCA / Eigensolver (Phase 4)
- **Nocedal & Wright Ch.5.1**: context for iterative methods, though the 3×3
  eigensolver is direct (not iterative). The PCA is a direct eigenvalue
  decomposition of a 3×3 symmetric matrix — use the closed-form cubic formula
  or 3–5 Jacobi rotations.
- **Skiena Ch.16.1** (numerical stability): test with known covariance matrices.
  Degenerate cases: all vertices coplanar (one eigenvalue ≈ 0), all vertices
  collinear (two eigenvalues ≈ 0).

### Micro-Hole Fill (Phase 5)
- **PIPELINE_REFERENCE §4**: boundary edge chaining, fan triangulation for
  small loops, centroid insertion for 5–6 vertex loops.
- **PIPELINE_REFERENCE §14**: boundary detection via edge occurrence counting.
  The same sparse-matrix pattern used throughout the pipeline.

### Memory Management (all phases)
- **Hanson Ch.6**: arena allocation. One arena per volume. All Phase 2–5
  allocations come from the arena. Scratch buffers (BFS queue, padded
  subvolume) use `Arena_save`/`Arena_restore` marks if memory pressure is
  a concern, otherwise let the arena hold everything until volume disposal.
- **Sanglard Ch.2**: DOOM's zone allocator — the historical precedent for
  arena allocation with purge tags. The pipeline's arena marks serve the
  same role as DOOM's `PU_LEVEL` / `PU_CACHE` tags.
- **Bentley Ch.10**: malloc overhead is 36–47 bytes per allocation. Arena
  eliminates this. Never use per-vertex or per-face malloc.

### Performance Framework
- **Abrash Ch.1–3**: measure first, then optimize. Profile Step 0 in isolation
  before optimizing. The Python bottleneck (trimesh split, 15–45s) is Python
  overhead that vanishes in C. The C bottleneck may be marching cubes or I/O.
- **Bentley Ch.6**: Appel's design levels — attack in order: problem definition,
  algorithm, data structure, code tuning, system, hardware.
- **Abrash Ch.16**: the "levels of optimization" framework. Step 0 is mostly
  Level 1 (straightforward port from Python to C). Level 2 opportunities:
  fuse connected-component labeling with marching cubes to avoid a separate
  pass. Level 3: unlikely needed given the time budget.

### Module Design
- **Hanson Ch.4–6**: opaque types, arena parameters, `sizeof *p` allocation.
- **Elements of Style Ch.2–3**: one module does one thing. Don't combine TIFF
  loading with connected components with marching cubes in one function. Split
  into composable operations.
- **TPOP Ch.4**: pass everything explicitly. The arena, the volume, the
  dimensions — all explicit parameters, no globals.

---

## Performance Budget

**Total pipeline budget**: ~180s per volume on Kaggle (2 cores @ 2.2 GHz). This 180s is the C preprocessing share of the ~270s/volume total (the remaining ~90s is nnUNet inference).
**Step 0 budget**: ~10–15s (combining mesh extraction 5–8s + cleanup 15–45s
from the Python profile, but the Python cleanup is dominated by trimesh overhead
that vanishes in C).

### Back-of-Envelope for Kaggle (2 cores @ 2.2 GHz)

| Sub-step | Work | Est. Time |
|---|---|---|
| TIFF load (320³ uint8) | 32.8 MB sequential read | ~0.1–0.3s |
| Threshold | 32.8M comparisons | ~0.05s (memory-bound) |
| BFS flood fill (6-conn) | 32.8M voxels, each visited once | ~0.5–1.0s |
| Component filtering + sort | ≤20 components, O(N) scan | negligible |
| Marching cubes (per comp) | ~50 ns/voxel × up to 32.8M total | ~1.6s total |
| Backface cull | M faces × 1 dot product | ~0.01s per component |
| Union-Find mesh CC | 3 unions per face × M faces | ~0.01s per component |
| Vertex reindex | O(N) scan + remap | ~0.01s per component |
| Micro-hole fill | O(boundary edges), typically small | ~0.01s per component |
| OBJ write (optional) | Text I/O, buffered | ~0.1–0.5s per component |
| **Total** | | **~3–5s** |

**Verdict**: Step 0 fits easily within budget. The Python profile shows 5–8s
for extraction + 15–45s for cleanup, but the cleanup cost is almost entirely
Python/trimesh overhead. In C with Union-Find instead of trimesh.split(), the
cleanup becomes O(M) — negligible. The dominant cost is marching cubes (~1.6s)
and BFS flood fill (~0.5–1.0s).

### Scaling: 7950X → Kaggle
- 7950X single-core is ~2.6× faster than Kaggle (5.7/2.2 GHz, plus IPC and
  cache advantages).
- Step 0 is not parallelizable at the voxel level (BFS is inherently serial),
  but marching cubes per-component can be parallelized across components.
- With 2 threads on Kaggle processing 2 components simultaneously, the per-
  component work overlaps. But with ≤20 components and the largest dominating,
  expect ~1.3–1.5× speedup from 2 threads.
- **Kaggle estimate**: 3–5s. Well within the 10–15s allocation.

### Does it fit in memory?
- Input volume: 32.8 MB (uint8, 320³)
- Label volume: 131 MB (int32, 320³) — **this is the largest single allocation
  in Step 0**. Consider using `int16_t` if component count ≤ 32767 (we keep
  at most 20), saving 65 MB. Or use a `uint16_t` label volume (65 MB).
- Padded subvolume per component: variable, but ≤ 35 MB for a full-volume
  component
- Vertex array (100K verts × 3 × float32): 1.2 MB
- Face array (200K faces × 3 × int32): 2.4 MB
- BFS queue: up to 32.8M entries × 4 bytes = 131 MB in worst case (all
  voxels foreground). **Optimization**: use a `uint32_t` packed index
  (z*H*W + y*W + x fits in uint32 for volumes up to ~1625³).
- Total peak: ~200–300 MB. Fits in Kaggle's ~13 GB.

---

## Key C Rules

Cite by section of `c-style-guide.md`:

**§1.1 Arena Architecture**: Allocate the entire Step 0 working set from the
per-volume arena. The BFS queue, label volume, padded subvolumes, vertex arrays,
face arrays — all arena-allocated. No `malloc` calls except inside libtiff
(which manages its own memory). The label volume and padded subvolumes are
scratch that can use `Arena_save`/`Arena_restore` to reclaim memory after each
component's marching cubes run. Vertex and face arrays for surviving components
persist in the arena for downstream steps.

**§1.3 Allocation Macros**: Use `ARENA_CALLOC(arena, D*H*W, sizeof(int32_t))`
for the label volume (needs zero-init). Use `ARENA_ALLOC(arena, n_verts * 3 * sizeof(float))`
for vertex arrays (filled immediately by marching cubes).

**§1.6 Stack Hazards**: Never put the volume or label array on the stack.
`uint8_t vol[320][320][320]` is 32 MB — instant stack overflow.

**§2.1 Module Structure**: `mesh_extract.h` exposes one public function and the
component descriptor type. All internal functions (BFS, marching cubes wrapper,
backface cull, Union-Find mesh split, micro-hole fill) are `static` in
`mesh_extract.c`. Or factor marching cubes into its own `marching_cubes.c`
module if reuse is anticipated.

**§4.6.1 Dense Arrays**: Vertices as `float verts[N*3]`, faces as
`int32_t faces[M*3]`. Not arrays of structs with pointers. Not linked lists.
Sequential access for cache-friendly processing.

**§5.1 Feasibility Checks**: The back-of-envelope above confirms Step 0 fits.
Always re-verify if volume dimensions change (400³ volumes would quadruple BFS
time but still fit).

**§5.4 Cache Locality**: The BFS flood fill accesses 6-neighbors per voxel.
In row-major layout, ±x neighbors are adjacent in memory (stride 1). ±y
neighbors have stride W (320). ±z neighbors have stride H*W (102,400). The z-
neighbor access is cache-unfriendly but unavoidable for 3D BFS. The working set
(vol + labels = ~165 MB with int32 labels, ~98 MB with uint16 labels) exceeds
L3. Expect DRAM-bound BFS. Reducing label size from int32 to uint16 saves a
cache level transition for the label array.

**§5.6 Threading**: Parallelize per-component marching cubes with OpenMP:
```c
#pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
for (size_t i = 0; i < n_components; i++) {
    extract_one_component(arena, vol, labels, comp_info[i], &out_meshes[i]);
}
```
Give each thread its own scratch arena. The shared volume and label arrays are
read-only during this phase. Use cutoff: don't parallelize if `n_components < 2`.

**§6.1 Volume Layout**: `vol[z*H*W + y*W + x]`. Coordinate convention `(z,y,x)`.
Document everywhere.

**§6.2 Connected Components**: BFS with flat arena-allocated queue. 6-connectivity.
Use `UINT32_MAX` as unvisited sentinel in distance array, or use the label
array itself (`labels[i] == 0` means unvisited foreground).

**§9.3 Defensive Programming**: Every function handles the empty case. Zero
foreground voxels → return empty component list. Zero faces from marching cubes
→ skip component. All face indices validated: `assert(faces[i] >= 0 &&
faces[i] < n_verts)`.

**§9.4 Float Discipline**: PCA eigenvector computation uses `double` internally
for numerical stability (3×3 covariance accumulation over 100K+ vertices).
Convert to `float` for the final normal vector. Vertex positions are `float`
throughout.

**§9.5 UB Avoidance**: Volume index arithmetic uses `size_t`:
`size_t idx = (size_t)z * H * W + (size_t)y * W + (size_t)x`. Never compute
`z*H*W` as `int` — 320×320×320 = 32.8M, fits int32 but leaves no headroom.

**§11.2 libtiff**: Wrap in `TiffIO_load()`. Serialize behind a mutex (libtiff
is not thread-safe). Check every return value. Propagate errors to the outer
loop for fallback to raw prediction.

---

## Data Structures

### From `common/` (shared modules)

| Module | Type / Function | Used For |
|---|---|---|
| `arena.h` | `Arena_T`, `Arena_save`, `Arena_restore` | All allocation |
| `tiff_io.h` | `TiffIO_load()` | Load input TIFF |
| `union_find.h` | `UnionFind`, `uf_find()`, `uf_union()` | Mesh CC in Phase 5 |

### Created internally by Step 0

```c
/* Component descriptor — output of Step 0.
 * Canonical definition: see CLAUDE.md "Canonical ComponentMesh".
 * All STAGE docs share this single struct layout. */
typedef struct ComponentMesh {
    float   *verts;          /* [nv * 3], float32, (z,y,x) */
    int32_t *faces;          /* [nf * 3], int32, 0-based */
    size_t   nv;
    size_t   nf;
    int      comp_id;        /* 1-based component number */
    float    pca_normal[3];  /* unit normal, sign-corrected (set by Step 0) */
    float    centroid[3];    /* mean vertex position (set by Step 0) */
    void    *self;           /* validation sentinel: self == &this_struct */
} ComponentMesh;
```

**Scratch data** (lives within Step 0, reclaimed via arena marks):
- `uint8_t *padded_subvol` — per-component padded binary subvolume
- `int32_t *labels` — voxel-level component labels (or `uint16_t` for memory)
- `uint32_t *bfs_queue` — flat queue, max size = D×H×W (packed linear index)
- `int32_t *comp_sizes` — voxel count per component (temporary for sorting)

**Vertex storage note**: The Python pipeline uses `float64` for vertices. The
C port uses `float32`. The SurfaceDice tolerance is τ=2.0 voxels. Float32 has
~7 decimal digits of precision, giving sub-microvoxel accuracy. Float64 is
wasted precision that halves cache utilization. See c-style-guide §5.4.6–5.4.7.

**Face index storage note**: Python uses `int64`. The C port uses `int32_t`.
A mesh with 2M vertices fits int32 (max 2.1B). A single 320³ component
will never exceed ~2M vertices. If somehow a 512³ volume produces a mesh
exceeding 2B vertices, we have bigger problems.

---

## Interfaces

### Input
```c
/*
 * MeshExtract_run — Step 0 entry point.
 *
 * Loads the TIFF at `tiff_path`, extracts connected components,
 * runs marching cubes + backface cull + cleanup on each.
 *
 * Returns 0 on success, nonzero on failure.
 * On success, *out_meshes points to an arena-allocated array of
 * ComponentMesh structs, and *out_n_meshes is the count.
 *
 * On failure (TIFF load error, zero components, OOM), the caller
 * should fall back to the raw prediction.
 */
int MeshExtract_run(Arena_T          arena,
                    const char      *tiff_path,
                    size_t           cube_D,
                    size_t           cube_H,
                    size_t           cube_W,
                    int              n_threads,
                    ComponentMesh  **out_meshes,
                    size_t          *out_n_meshes);
```

### Output

The `ComponentMesh` array is arena-allocated and lives until `Arena_dispose`.
Downstream steps (Oracle, Bridge Cut, etc.) receive individual `ComponentMesh`
pointers.

### Connection to Step 1 (Oracle)

Step 1 iterates over the output array:
```c
for (size_t i = 0; i < n_meshes; i++) {
    int sheet_count = Oracle_count_sheets(arena, &meshes[i], grid_size);
    if (sheet_count > 1) {
        /* mark for bridge cut in Step 2 */
    }
}
```

The PCA normal computed in Step 0 is reused by the Oracle (Step 1) for
raycast direction, by Bridge Cut (Step 2) for source/sink seed projection,
and by Raycast Separator (Step 3) for global normal sign consistency.

---

## External Dependencies

### libtiff (Step 0 input)

- **Used for**: reading the multi-page TIFF prediction volume.
- **Call pattern**: `TIFFOpen(path, "r")` → loop over pages with
  `TIFFReadScanline()` or `TIFFReadEncodedStrip()` → `TIFFClose()`.
- **Wrap in**: `tiff_io.c` with mutex. Return arena-allocated `uint8_t *vol`
  and dimensions `(D, H, W)`.
- **Error handling**: if `TIFFOpen` returns NULL or any read fails, return
  error code. Caller falls back to raw prediction.
- **Performance**: sequential read of 32.8 MB. Not a bottleneck (~0.1–0.3s).
  Use strip-based API for efficiency.
- **Thread safety**: libtiff is NOT thread-safe. All calls serialized behind
  `static pthread_mutex_t tiff_lock`.

### No other external dependencies in Step 0

Marching cubes, Union-Find, PCA, and micro-hole fill are all hand-rolled C.
No PoissonRecon, no Triangle, no Clipper2 in this step.

---

## Test Strategy

### Comparison with Python Reference
The Python pipeline's `extract_meshes_from_tiff()` is the reference. For each
test volume, run both Python and C, then compare:

1. **Component count**: must match exactly. If C finds 12 components and Python
   finds 13, there's a bug in BFS connectivity or size filtering.

2. **Per-component vertex count**: should match within ~1% (marching cubes
   implementations may differ in edge interpolation or degenerate handling,
   producing slightly different vertex counts). If counts differ by >5%,
   investigate.

3. **Per-component face count after backface cull**: should be approximately
   half the pre-cull count. Exact match with Python is unlikely (PCA sign
   convention must match exactly; floating-point tiebreakers in eigenvector
   computation may flip a few faces). Verify that the cull removes ~45–55% of
   faces consistently.

4. **PCA normal direction**: compare `dot(c_normal, python_normal)`. Should be
   >0.99 (nearly identical direction). A negative dot product means the sign
   convention is flipped — catastrophic for downstream steps.

5. **Mesh validity invariants** (check programmatically after every run):
   - All face indices in `[0, n_verts - 1]`.
   - No degenerate triangles (three identical vertex indices).
   - No unreferenced vertices (every vertex appears in at least one face).
   - Boundary edge count after micro-hole fill ≤ boundary count before
     (holes filled, not created).

### Boundary-Value Tests
[c-style-guide §10.1.4; Elements Ch.4; TPOP Ch.6]

- **Empty volume** (all zeros): should return 0 components, no crash.
- **Single voxel** (1×1×1, value=1): marching cubes on a 3×3×3 padded volume.
  May produce 0 faces (single voxel below MC threshold). Handle gracefully.
- **All-ones volume** (320³ all foreground): one giant component. Marching cubes
  produces a box mesh. Backface cull keeps three faces of the box. Verify.
- **Two components separated by 1 voxel gap**: must label as separate. Verify
  6-connectivity does NOT merge them (26-connectivity would).
- **Component with exactly MIN_CC_SIZE (500) voxels**: must survive filtering.
  Component with 499 voxels: must be discarded.
- **More than MAX_COMPONENTS (20) components**: keep top 20 by size, discard rest.
- **Volume dimensions not 320³**: test 256³, 400³, 128×256×320 (non-cubic).
  All index arithmetic must derive from actual dimensions.
- **Volume with 2 label regions for which marching cubes produces 0 faces**:
  handle gracefully (skip component, log warning).

### Stress Tests
[TPOP Ch.6 — machine-generated input]

- Random binary volumes of sizes 1³ through 512³.
- Volumes with pathological topology: a thin bridge of 1 voxel connecting
  two large regions (should be one component in 6-connectivity).
- Checkerboard patterns (maximizes component count).

### Regression Suite
- One test per bug found, stored permanently. The suite grows monotonically.
- `make test_step0` runs all Step 0 tests in <30s on the 7950X.
- Compare with `taskset -c 0,1` periodically to validate Kaggle timing.

---

## Pitfalls

### Off-by-One in Padding
The padded subvolume is `(dz+2) × (dy+2) × (dx+2)`. Vertex positions from
marching cubes are relative to the padded volume. You must subtract 1 (undo
padding) AND add the bounding-box origin. Getting either offset wrong shifts
all vertices by 1 voxel — invisible unless you compare against Python output.
[Elements Ch.4 — off-by-one; c-style-guide §9.3.3]

### PCA Sign Convention
The eigenvector of the smallest eigenvalue has arbitrary sign. If C and Python
use different eigensolver implementations, the sign may flip. The backface cull
keeps faces whose normal dots **positively** with the PCA normal — a sign flip
keeps the wrong half. Always apply the sign convention: "make the largest-
magnitude component positive." Test with a known planar mesh whose normal you
can verify by hand. [PIPELINE_REFERENCE §4, step 6d]

### 6-Connectivity vs 26-Connectivity
Step 0 uses **6-connectivity** for voxel-level connected components. The final
VOI metric uses **26-connectivity**. Two voxel regions that are 6-disconnected
but 26-connected will be separate components in Step 0 but could be merged in
the final evaluation. This is by design — the pipeline operates on 6-connected
components and the revoxelization (Step 6) with bridge erasure handles the
26-connectivity requirement. Do NOT use 26-connectivity in Step 0.
[PIPELINE_REFERENCE §4]

### Integer Overflow in Index Arithmetic
`z * H * W` for z=319, H=320, W=320: 319 × 320 × 320 = 32,665,600. This fits
`int32_t` but `z * H * W + y * W + x` for the maximum voxel is 32,767,999 —
still fits int32 but barely. For a 512³ volume: 511 × 512 × 512 = 134,086,656.
Fits int32. But for safety, always use `size_t` for index computation.
[c-style-guide §9.4.4; Pitfalls — semantic]

### Union-Find with Non-Contiguous Vertex Indices
After backface cull removes faces, vertex indices may have gaps (e.g., vertex 17
is unreferenced). The Union-Find array is sized `n_verts` (the original count).
All indices are valid. After Union-Find groups components, the vertex reindex
pass compacts the array. Do NOT run Union-Find on already-reindexed vertices
unless you rebuild it. Order matters: (1) backface cull faces, (2) Union-Find
on full vertex set, (3) filter small components, (4) reindex vertices.

### Degenerate Marching Cubes Output
Some cube configurations produce zero-area triangles (three collinear vertices)
or duplicate triangles. The Lorensen & Cline table doesn't guarantee non-
degenerate output. After marching cubes, optionally scan for degenerate faces
(cross product magnitude < ε) and remove them. This prevents NaN in normal
computation downstream. [Elements Ch.4 — validate output]

### BFS Queue Sizing
In the worst case (entire volume is foreground), the BFS queue needs D×H×W
entries. For a 320³ volume, that's 32.8M × 4 bytes = 131 MB. This is large
but fits in the arena. Do NOT assume a smaller queue. If memory is tight,
use a two-pass approach: first pass counts component sizes, second pass
extracts each component separately.

### libtiff Page Ordering
libtiff reads pages in order, and the pipeline assumes page 0 = z=0,
page 1 = z=1, etc. Verify that the TIFF files in the test set follow this
convention. Some TIFF writers produce pages in reverse order. Check
`TIFFCurrentDirectory()` after each page read.

### Floating-Point PCA on Large Vertex Sets
Computing the covariance matrix as `sum((v - mean) * (v - mean)^T)` over 100K+
vertices in float32 accumulates significant rounding error. Use `double` for the
accumulator (6 unique elements of the 3×3 symmetric matrix) and convert the
resulting eigenvector to `float`. The mean should also be computed in `double`.
[c-style-guide §9.4; Nocedal — numerical hygiene]
