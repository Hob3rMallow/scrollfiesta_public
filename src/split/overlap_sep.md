# STAGE 3 — Raycast Sheet Separator

**Module prefix**: `RaycastSep_`
**Source files**: `raycast_sep.c`, `raycast_sep.h`
**Dependencies**: `arena.c`, `kdtree.c`, `pca.c`, `csr.c`, `bfs.c`, `mesh.c` (from `common/`)

---

## Purpose

The raycast sheet separator handles multi-sheet components that bridge cut
(Step 2) could not resolve. Bridge cut severs thin necks; this step separates
**broadly overlapping sheets** where the connection is not a thin bridge but
a wide zone of shared faces. These components passed through bridge cut with
`is_bridge == false` despite the oracle reporting ≥ 2 sheets.

The algorithm assigns each vertex a **sheet label** (integer ≥ 1) by projecting
vertices into a UV coordinate system derived from local surface normals, sorting
along the stacking direction (w), and detecting sheet boundaries via gaps. It
then cleans these labels with BFS flood fill and cuts the mesh at label
boundaries, producing per-sheet sub-meshes.

**Input**: A `ComponentMesh` (vertices, faces) that the oracle says has
sheet_count ≥ 2, and that bridge cut returned unchanged (no thin bridge found).

**Output**: A list of single-sheet sub-meshes (or the original mesh unchanged
if separation fails or isn't needed).

**Why it matters for the metric**:

- **VOI (35% of score)**: Any inter-sheet merger that survives to the output
  degrades VOI. The raycast separator is the *last chance* to separate sheets
  that bridge cut missed. If this step fails, overlapping wraps merge in the
  final revoxelization, producing catastrophic VOI penalties.
- **TopoScore (30% of score)**: Merged sheets create wrong Betti numbers.
  Successful separation restores correct topology.
- **SurfaceDice (35%)**: Cutting too aggressively (removing real surface)
  degrades surface proximity. The 6-hop gap expansion at cut boundaries is
  calibrated to produce clean cuts without excessive surface loss.
- **Failure mode**: If this step crashes, the outer loop falls back to the
  raw prediction. Design so that every degenerate case returns the mesh
  unchanged rather than crashing. A wrong cut (merging labels that should be
  distinct) is worse than no cut at all.

---

## Algorithm

Four phases, executed sequentially within a single component.

### Phase 1: Compute Per-Vertex Local Normals

1. **Global PCA normal.** Compute the PCA normal of all vertices (smallest
   eigenvector of the 3×3 covariance matrix). Sign correction: make the
   largest absolute component positive. This establishes a consistent
   "up" direction for orientation flipping.
   [PIPELINE_REFERENCE §7, Phase 1 step 1; reuse PCA from Step 0/Step 2 if
   already computed — c-style-guide §6.6.3]

2. **KD-tree construction.** Build a single KD-tree over all N vertices.
   This tree will be queried N times (ball query, radius = 10.0).
   [c-style-guide §6.4.1–6.4.6; Skiena Ch.15.6; Sanglard Ch.5]

3. **Per-vertex ball query + local PCA.** For each vertex `v`:
   - Query the KD-tree: find all vertices within Euclidean distance 10.0.
   - Compute the 3×3 covariance matrix of the neighbor positions.
   - Eigendecompose → smallest eigenvector = local normal.
   - If `dot(local_normal, global_normal) < 0`, flip local normal.
   - Store the local normal in a flat `float local_normals[N*3]` array.
   [PIPELINE_REFERENCE §7, Phase 1 steps 2–4]

4. **Parallelization.** This loop is embarrassingly parallel — no shared
   mutable state (the KD-tree is read-only, each vertex writes to its own
   slot in the output array). Use OpenMP:
   ```c
   #pragma omp parallel for schedule(dynamic, 64) num_threads(n_threads)
   for (size_t i = 0; i < n_verts; i++) {
       /* ball query + local PCA for vertex i */
   }
   ```
   Per-thread scratch (neighbor index buffer, 3×3 covariance matrix) should
   be thread-local stack variables or per-thread arena marks.
   [c-style-guide §5.6.1–5.6.6; Hanson Ch.20; Sanglard Ch.7]

### Phase 2: Ray Casting to Assign Sheet Labels

1. **Mean normal.** Compute the mean of all local normals → `mean_normal`.
   Normalize to unit length.

2. **UV basis.** Build an orthonormal basis `(u_axis, v_axis, mean_normal)`
   perpendicular to `mean_normal`. Use Gram-Schmidt or the Hughes-Möller
   technique for robustness when `mean_normal` is near-axis-aligned.

3. **Project vertices to (u, v, w).** For each vertex:
   - `u = dot(vertex, u_axis)`
   - `v = dot(vertex, v_axis)`
   - `w = dot(vertex, mean_normal)`

4. **Direction correction.** Compute the component's centroid. If the centroid
   is below the cube center along `mean_normal`, negate all w-values (flip
   the stacking direction so w increases toward the umbilicus).
   [PIPELINE_REFERENCE §7, Phase 2 step 4]

5. **Adaptive UV grid.** Grid size = `(int)sqrt(N_verts / 50.0)`, clamped
   to a reasonable range (e.g., [4, 200]). Compute the UV bounding box, divide
   into `grid_size × grid_size` cells.

6. **Bin vertices by (u, v) cell.** For each vertex, compute its cell index
   via integer division. Store as a CSR-like structure: `cell_offsets[]` and
   `cell_vertex_indices[]`, with vertices sorted by cell.
   [c-style-guide §6.4.8; Sanglard Ch.5 — CSR cell lists]

7. **Per-cell gap labeling.** For each non-empty cell:
   - Sort the cell's vertices by w-coordinate (use `qsort` or an in-place
     sort on the cell's index range).
   - Walk sorted vertices: assign label = 1. When the w-gap between
     consecutive vertices exceeds `GAP_THRESHOLD = 8.0`, increment label.
   - Write labels to a flat `int32_t labels[N]` array indexed by vertex.
   [PIPELINE_REFERENCE §7, Phase 2 step 7]

   Vertices not assigned to any cell (shouldn't happen if the grid covers
   the bounding box, but handle defensively) get label = 0 (unlabeled).

### Phase 3: Boundary Cleanup

The raw gap-based labels are noisy near mesh boundaries. This phase cleans
them up with a three-stage process.

1. **Find mesh boundary vertices.** Scan all edges (from the face array).
   Boundary edges appear in exactly one triangle. Vertices on boundary edges
   are boundary vertices. Implementation: build a half-edge count array
   (or use the CSR adjacency to check edge valence).
   [PIPELINE_REFERENCE §7, Phase 3 step 1]

2. **Clear labels near boundary.** BFS from all boundary vertices, expanding
   25 hops along mesh edges (using the CSR adjacency graph). Set
   `labels[v] = 0` for all vertices within 25 hops of any boundary vertex.
   This removes unreliable labels that were assigned near the ragged edge.
   [PIPELINE_REFERENCE §7, Phase 3 step 2]

3. **Flood fill from interior.** BFS from all vertices that still have a
   nonzero label. For each visited unlabeled neighbor: assign the propagating
   label. This fills the cleared boundary zone from the reliable interior
   labels outward. Standard multi-source BFS on the CSR adjacency graph.
   [PIPELINE_REFERENCE §7, Phase 3 step 3; CLRS Ch.22.2]

4. **Label disconnected unlabeled components.** Any vertices still unlabeled
   after flood fill are in disconnected regions. For each such vertex, find
   its nearest labeled vertex using a 3D KD-tree nearest-neighbor query.
   Assign the nearest label. (Reuse the same KD-tree from Phase 1, but
   restrict the search to labeled vertices — or build a second smaller tree
   over labeled vertices only.)
   [PIPELINE_REFERENCE §7, Phase 3 step 4; Skiena Ch.15.6]

5. **Small-component cleanup.** Compute connected components of same-label
   regions (Union-Find on the CSR graph restricted to same-label edges).
   For each label-component: if its size is < 30% of the largest adjacent
   label's component, merge it into that adjacent label. Use **frozen
   original sizes** to prevent cascading merges — compute all sizes first,
   then apply all merges.
   [PIPELINE_REFERENCE §7, Phase 3 step 5]

### Phase 4: Cut at Label Boundaries

Function: `cut_at_label_boundaries(vertices, faces, labels, gap_hops=6)`

1. **Find label-boundary vertices.** For each vertex, check its CSR
   neighbors. If any neighbor has a different label, the vertex is on a
   label boundary.

2. **Expand boundary by 6 hops.** BFS from all boundary vertices, 6 hops
   deep. Mark these vertices as "in the cut zone."
   [PIPELINE_REFERENCE §7, Phase 4 step 2]

3. **Remove cut-zone faces.** For each face (a, b, c): if ANY vertex is in
   the cut zone, discard the face.

4. **Extract per-label sub-meshes.** For each unique label:
   - Collect all faces whose three vertices share that label and are NOT in
     the cut zone.
   - Build a vertex reindex map (old → new).
   - Emit reindexed faces and the corresponding vertex subset.
   - Run mesh connected-components (Union-Find) on the sub-mesh. Discard
     fragments with < `MIN_FRAGMENT_FACES` (e.g., 50) faces.

5. **Return sub-meshes.** If only one label was found (separation failed),
   return the original mesh unchanged.

---

## Book References

### KD-Tree Construction and Ball Query
- **Skiena Ch.15.6**: KD-tree construction (O(N log N)), ball query
  algorithm (prune subtrees by bounding-box/sphere intersection), 3D
  optimality, flat-array vs linked representation.
  [c-style-guide §6.4.1–6.4.8; Skiena rules 51–56]
- **Sanglard Ch.5 (DOOM blockmap)**: uniform grid as alternative to
  KD-tree for roughly uniform vertex distributions. CSR storage for
  variable-length per-cell lists. Profile grid vs KD-tree.
  [c-style-guide §6.4.7–6.4.8; Sanglard rules 33–35]
- **Bentley Ch.9**: squared-distance comparison to avoid sqrt in ball
  query inner loop. Feasibility check methodology.
  [c-style-guide §6.4.5]

### PCA and Eigendecomposition
- **PIPELINE_REFERENCE §5, §6, §7**: PCA normal computation (covariance
  matrix → smallest eigenvector). Sign correction convention. Reuse
  across steps.
- **Nocedal & Wright** (if distilled): numerical stability of 3×3
  eigendecomposition. For the pipeline's small 3×3 case, a closed-form
  solution or Jacobi iteration is appropriate — no need for LAPACK.

### BFS and Flood Fill
- **CLRS Ch.22.2**: BFS algorithm, correctness proof, O(V+E) time.
  Multi-source BFS for the flood-fill step (initialize queue with all
  labeled vertices simultaneously).
  [c-style-guide §7.2]
- **Skiena Ch.7.7**: connected components via BFS. Union-Find alternative
  for the small-component cleanup step.
  [c-style-guide §6.2.1]

### CSR Graph Representation
- **CLRS Ch.22.1**: adjacency list representation. CSR is the
  array-packed version.
- **Skiena Ch.15.4**: packed arrays for static graphs — 4× speedup
  over linked lists.
  [c-style-guide §7.1.1–7.1.5]

### Optimization Framework
- **Abrash Ch.3–7 (c-style-guide §5.2)**: profile before optimizing.
  The KD-tree ball query loop is the hot spot — focus there. Use
  `perf record` to confirm.
- **Abrash Ch.8–10**: the three levels of optimization. Level 1: tighter
  KD-tree traversal (prune earlier, smaller nodes). Level 2: replace
  KD-tree with uniform grid (eliminate tree traversal entirely). Level 3:
  compute local normals from mesh adjacency (k-ring CSR walk) instead
  of spatial proximity (ball query), eliminating the KD-tree entirely.
  [Abrash rules 8–10; Bentley rule 10]
- **Abrash Ch.17–18**: sparse processing — if most vertices are on a
  single sheet, the gap detection rarely triggers. Design the inner loop
  for the common case (no gap).
  [Abrash rules 35–38]
- **Bentley Ch.1**: define the actual problem before choosing an
  algorithm. The problem is "compute local surface normal." The KD-tree
  ball query is one solution. A k-ring walk on mesh adjacency is another.
  A grid-binned neighbor lookup is a third. Evaluate all three.
  [Bentley rules 9–10]
- **Bentley Ch.7**: back-of-envelope estimation. See §Performance Budget.
  [c-style-guide §5.1.1–5.1.5]

### Memory Management
- **Hanson Ch.6**: arena allocation for all Phase 1–4 data. Save/restore
  for scratch within phases.
  [c-style-guide §1.1–1.2]
- **Hanson Ch.16**: save/restore checkpoints. Use one mark per phase
  if intermediate data can be discarded.
  [Hanson rules 62–64]
- **Hanson Ch.20**: threading patterns. Per-thread scratch arenas for the
  OpenMP ball-query loop. Read-only shared KD-tree requires no mutex.
  [Hanson rules 53–55, 65]
- **Sanglard Ch.2**: DOOM zone allocator precedent. Tag-by-lifetime
  pattern maps directly to arena marks.

### Data Layout and Cache
- **Abrash Ch.15–17 (c-style-guide §5.4)**: working-set analysis. Vertex
  array + KD-tree + local normals must fit in L3 for acceptable
  performance. See §Performance Budget for size estimates.
- **Sanglard Ch.3**: cacheline-oriented layout. Vertices as flat
  `float[N*3]` for stride-1 access during projection.
  [c-style-guide §5.4.2–5.4.3]

### Mesh Operations
- **Abrash Ch.38–39 (c-style-guide §6.5)**: separate "what to process"
  from "do the processing." Phase 2 (labeling) is separated from Phase 4
  (cutting). The intermediate label array is the "horizontal line list"
  equivalent.
  [Abrash rules 60–61]
- **Abrash Ch.38 (c-style-guide §6.5)**: ownership rule for vertices
  at cut boundaries — each vertex belongs to exactly one label. The
  cut zone removes ambiguous vertices entirely.
  [Abrash rules 62–63]

---

## Performance Budget

**Total pipeline budget**: ~180s per volume on Kaggle (2 cores @ 2.2 GHz). This 180s is the C preprocessing share of the ~270s/volume total (the remaining ~90s is nnUNet inference).
**Step 3 budget**: 30–110s in Python. Target for C: **1–5s** per component
on Kaggle. For a volume with ≤ 3 components reaching Step 3: **2–10s total.**

### Back-of-Envelope: Single Component (N=100K verts, M=200K faces)

| Sub-step | Work | Est. Time (Kaggle) |
|---|---|---|
| Global PCA | O(N) = 100K × ~10 ops | ~0.2 ms |
| KD-tree build | O(N log N) = 100K × 17 | ~20 ms |
| N ball queries (r=10) | 100K × ~5 µs/query | **~0.5 s** |
| N local PCAs (3×3 eigen, ~12 neighbors avg) | 100K × ~200 ops | ~10 ms |
| Mean normal + UV projection | O(N) | ~1 ms |
| UV grid binning (CSR build) | O(N) | ~2 ms |
| Per-cell sort (N verts, ~50 cells avg) | O(N log(N/G²)) ≈ 100K × ~7 | ~5 ms |
| Gap labeling | O(N) | ~1 ms |
| Boundary detection | O(M) = 200K | ~3 ms |
| 25-hop BFS clear | O(boundary × 25 × degree) ≈ ~50K | ~1 ms |
| Flood fill BFS | O(N + E) ≈ 700K | ~7 ms |
| Nearest-label KD-tree queries | O(K × √N), K ≈ small | ~1 ms |
| Small-comp cleanup (Union-Find) | O(N α(N)) ≈ N | ~2 ms |
| Cut zone expansion (6-hop BFS) | O(boundary × 6 × degree) | ~1 ms |
| Face filtering + sub-mesh extraction | O(M) | ~5 ms |
| **Total** | | **~0.6 s** |

**The KD-tree ball query dominates** at ~80% of runtime. Everything else
combined is ~100ms. This matches the Python profile where the per-vertex
KD-tree query is the bottleneck.

### Feasibility Cross-Check (Bentley §7 method)

Second estimate via memory bandwidth: N ball queries, each touching ~√N
KD-tree nodes × 32 bytes/node = 100K × 316 × 32B = ~1 GB of reads.
At Kaggle's memory bandwidth (~20 GB/s effective with L3): ~50 ms if fully
cached. But random access defeats prefetch → effective bandwidth ~2 GB/s
for pointer-chasing → ~500 ms. Consistent with the 0.5s estimate above.

### Working Set Analysis

| Data structure | Size (N=100K) | Cache level |
|---|---|---|
| Vertices float[N×3] | 1.2 MB | L2 |
| KD-tree nodes (N × 32B) | 3.2 MB | L3 |
| Local normals float[N×3] | 1.2 MB | L2 |
| Labels int32[N] | 400 KB | L2 |
| CSR adjacency (from faces) | ~4.8 MB | L3 |
| UV grid CSR | ~1.2 MB | L2 |
| BFS queue + visited | ~500 KB | L2 |
| **Total working set** | **~12.5 MB** | **Fits L3** |

On Kaggle (Xeon with ≥ 20 MB L3): fits comfortably. On 7950X (64 MB L3):
fits trivially. No DRAM spill expected. This is much better than Step 2's
flow network (~42 MB).

For larger meshes (N=200K): doubles to ~25 MB. Still fits L3 on both
platforms.

### Scaling: 7950X → Kaggle

- **CPU speed**: 7950X single-core ~2.6× faster than Kaggle @ 2.2 GHz.
  A 0.6s run on Kaggle ≈ 0.23s on 7950X.
- **Threading**: Phase 1 (ball queries) parallelizes perfectly. With 2
  Kaggle threads: ~0.3s for queries. With 16 dev threads: ~0.03s. The
  serial phases (~100ms) don't benefit from threading.
- **Total per component on Kaggle (2 threads)**: ~0.4s. For 3 components:
  ~1.2s. Well within budget.
- **Worst case (N=200K, 5 components)**: ~4s. Still fits.

### The Grid Alternative (Abrash Level 2)

If the KD-tree ball query proves slower than estimated (e.g., non-uniform
vertex density in some volumes causes 100+ neighbors per query):

A **3D uniform grid** with cell size = 10.0 (matching the ball query radius)
converts each ball query into visiting the 3×3×3 = 27 neighboring cells.
With CSR storage, this is 27 array lookups + distance checks on the
cell contents. For roughly uniform vertex density:

- Grid dimensions: ~32×32×32 = 32K cells for a 320³ volume.
- Grid CSR: cell_offsets[32K+1] + vertex_indices[100K] ≈ 530 KB.
- Per query: 27 cell lookups × ~12 verts/cell × distance check = ~324 ops.
- Total: 100K × 324 × ~2 ns = ~65 ms. **8× faster than KD-tree.**

The grid is simpler to implement, uses less memory, and has fully
predictable access patterns (no tree traversal). **Profile both. Start
with KD-tree (proven correct), fall back to grid if needed for speed.**
[c-style-guide §6.4.7; Sanglard Ch.5; Bentley rule 10]

---

## Key C Rules

### §1 Memory Management
Use arena for everything in this step. Allocate on entry, dispose is handled
by the per-volume arena at the end of `run_pipeline_for_cube()`.

**§1.2.2 Tag by lifetime.** Phase 1's KD-tree and local normals persist
through Phase 4. Phase 2's UV grid and cell lists are scratch — use
`Arena_save`/`Arena_restore` if memory pressure is a concern, but given
the ~12 MB working set, this is optional. The CSR adjacency graph persists
across all phases.
[Hanson Ch.6, Ch.16; c-style-guide §1.2.2]

**§1.4.2 No malloc in hot paths.** The per-vertex ball query returns a
variable-length neighbor list. Do NOT malloc per query. Pre-allocate a
fixed-size buffer (e.g., `int32_t neighbors[MAX_NEIGHBORS]` where
`MAX_NEIGHBORS = 512`) on the stack or in per-thread scratch. If a query
exceeds MAX_NEIGHBORS, truncate (document the truncation). For N=100K
in a 320³ volume, expected neighbors ≈ 12.5; 512 provides >40× headroom.
[Hanson Ch.6; Bentley Ch.10; c-style-guide §1.4.2]

### §2 Module Design
**§2.1 Opaque type.** The raycast separator's internal state (KD-tree,
normals, UV grid, labels) should not leak into the header. The header
exposes only `RaycastSep_process()`.
[Hanson Ch.4; c-style-guide §2.1]

**§2.4 Return int status.** `RaycastSep_process()` returns 0 on success,
-1 on failure (timeout, degenerate input). On failure, `*out_meshes`
is set to the input mesh and `*out_count = 1`.
[c-style-guide §2.4–2.5]

### §4 Type Safety
**§4.4.1 Signed/unsigned.** Loop counters for vertices are `size_t`.
Labels are `int32_t` (can be 0 for unlabeled). The gap threshold (8.0)
is `float`. The hop counts (25, 6) are `int`. Never compare `size_t`
with a negative `int`.
[K&R Ch.2; Pitfalls — semantic; c-style-guide §4.4.1–4.4.3]

**§4.4.4 Float comparison.** The gap detection uses `if (w_diff > 8.0f)`.
This is a strict inequality — not an epsilon comparison. Document this.
If the Python reference uses `> 8.0` with float64, the float32 port may
produce slightly different label assignments due to rounding. Test on all
validation volumes and accept minor differences if VOI/TopoScore are
unaffected.
[Elements Ch.5; c-style-guide §9.4.1–9.4.2]

### §5 Performance
**§5.2 Profile first.** The Phase 1 ball query loop is the predicted
hot spot. Confirm with `perf record` before optimizing anything else.
If Phase 1 is <50% of runtime, something else is wrong.
[Abrash Ch.3; c-style-guide §5.2.1]

**§5.4 Cache locality.** Vertices as `float[N*3]` gives stride-1 access
during UV projection (Phase 2). The KD-tree is read-only during queries —
its flat-array layout ensures spatial locality for tree traversal.
[Abrash Ch.15; Sanglard Ch.3; c-style-guide §5.4.1–5.4.3]

**§5.6 Threading.** Phase 1 is the only phase worth parallelizing.
Phases 2–4 are O(N) serial operations totaling ~100ms — threading
overhead exceeds the benefit. Use `#pragma omp parallel for` with
`schedule(dynamic, 64)` for the ball query loop. Chunk size 64 balances
thread overhead against load imbalance.
[Sanglard Ch.7; Hanson Ch.20; c-style-guide §5.6.1–5.6.5]

### §6 Algorithm Implementations
**§6.4.2 Flat-array KD-tree.** Heap layout: children of node `i` at
`2*i+1` and `2*i+2`. No pointer chasing. 32 bytes per node fits 2 nodes
per 64-byte cacheline.
[Skiena Ch.15.6; Sanglard Ch.5; c-style-guide §6.4.2]

**§6.4.5 Squared distance.** Never compute `sqrt()` in the ball query.
Compare `dx*dx + dy*dy + dz*dz < r_sq` where `r_sq = 100.0f`.
[Skiena Ch.15.6; Bentley Ch.9; c-style-guide §6.4.5]

**§6.6.1 Build all structures in one pass.** When entering Step 3, build
the KD-tree, CSR adjacency, and global PCA normal in a single pass through
the vertex/face data while it's hot in cache. Don't read the vertex array
three separate times.
[Sanglard Ch.2; c-style-guide §6.6.1]

**§6.6.3 Reuse the KD-tree.** The same KD-tree built in Phase 1 serves
Phase 3 (nearest-labeled-vertex queries) and Phase 4 if needed. Don't
build a second tree.
[Abrash Ch.15; c-style-guide §6.6.3]

---

## Data Structures

### From `common/` (shared with other steps)

```c
/* Arena — region allocator */
typedef struct Arena_T *Arena_T;

/* KD-tree — flat-array spatial index */
typedef struct KDTree_T *KDTree_T;
KDTree_T KDTree_new(Arena_T arena, const float *points, size_t n, int dim);
size_t   KDTree_ball_query(const KDTree_T tree, const float *center,
                           float radius_sq, int32_t *out_indices,
                           size_t max_results);
size_t   KDTree_nearest(const KDTree_T tree, const float *query,
                        float *out_dist_sq);

/* CSR — compressed sparse row graph */
typedef struct CSR_T *CSR_T;
CSR_T  CSR_from_faces(Arena_T arena, const int32_t *faces, size_t nf,
                      size_t nv);
/* Access: neighbors of v are target[offset[v] .. offset[v+1]) */

/* Mesh — flat vertex/face arrays */
typedef struct ComponentMesh {
    float   *verts;          /* [nv * 3], float32, (z,y,x) */
    int32_t *faces;          /* [nf * 3], int32, 0-based */
    size_t   nv;
    size_t   nf;
    int      comp_id;
    float    pca_normal[3];  /* unit normal from Step 0 PCA */
    float    centroid[3];    /* mean vertex position from Step 0 */
    void    *self;           /* validation sentinel: self == &this_struct */
} ComponentMesh;

/* PCA — 3×3 covariance eigenvectors */
void PCA_normal(const float *points, size_t n, float out_normal[3]);
```

### Internal to `raycast_sep.c`

> **NOTE**: The variables below illustrate memory layout only. In the actual
> implementation, these MUST be local variables or fields of a context struct
> passed explicitly to functions — NOT module-level statics. Module-level
> `static` pointers create hidden state, prevent concurrent processing of
> multiple components, and violate the project convention that arena pointers
> are always explicit parameters (CLAUDE.md rule 4, c-style-guide §2,
> TPOP Ch.4). The `static` keyword here means file-scoped visibility only.

```c
/* Per-vertex local normal array — flat float[N*3] */
float *local_normals;  /* arena-allocated */

/* Per-vertex sheet labels — flat int32[N] */
int32_t *labels;       /* arena-allocated, 0 = unlabeled */

/* UV grid binning — CSR layout */
int32_t *cell_offsets;     /* [grid_size * grid_size + 1] */
int32_t *cell_vert_ids;    /* [N], vertices sorted by cell */

/* Per-vertex UV projections — flat float[N*3] for (u, v, w) */
float *uvw;               /* arena-allocated */
```

All arrays are flat 1D, arena-allocated. No VLAs, no large stack arrays.
[c-style-guide §CLAUDE.md rule 3]

### Memory Layout (N=100K)

```
Arena allocation sequence (single contiguous region):
  local_normals:   float[100K × 3]  =  1.2 MB
  labels:          int32[100K]       =  400 KB
  uvw:             float[100K × 3]  =  1.2 MB
  cell_offsets:    int32[G² + 1]     =  ~16 KB  (G≈45)
  cell_vert_ids:   int32[100K]       =  400 KB
  CSR offsets:     int32[100K + 1]   =  400 KB
  CSR targets:     int32[600K]       =  2.4 MB
  KD-tree nodes:   32B × 100K       =  3.2 MB
  BFS queue:       int32[100K]       =  400 KB
  BFS visited:     uint8[100K]       =  100 KB
  -----------------------------------------------
  Total:                              ~9.7 MB
```

---

## Interfaces

### Primary Entry Point

```c
/*
 * RaycastSep_process — Step 3 entry point.
 *
 * Separate overlapping sheets in a multi-sheet component using
 * ray-based label assignment and mesh cutting.
 *
 * Parameters:
 *   arena       — arena for all allocations
 *   mesh        — input component mesh (oracle sheet_count >= 2,
 *                 bridge cut returned it unchanged)
 *   n_threads   — thread count for Phase 1 parallelization
 *   timeout_sec — max wall time for this component (e.g., 30.0)
 *   out_meshes  — output: array of resulting sub-meshes
 *   out_count   — output: number of resulting sub-meshes
 *
 * Returns 0 on success. On failure (timeout, degenerate mesh,
 * single-label result), returns -1 and sets *out_meshes pointing
 * to the original mesh, *out_count = 1.
 */
int RaycastSep_process(Arena_T              arena,
                       const ComponentMesh  *mesh,
                       int                   n_threads,
                       double                timeout_sec,
                       ComponentMesh       **out_meshes,
                       size_t               *out_count);
```

### Internal Functions (all `static`)

```c
/* Phase 1: compute per-vertex local normals via KD-tree ball query */
static int compute_local_normals(Arena_T arena, const ComponentMesh *mesh,
                                 KDTree_T tree, const float global_normal[3],
                                 int n_threads, float *out_normals);

/* Phase 2: project to UV, bin, gap-label */
static int assign_sheet_labels(Arena_T arena, const ComponentMesh *mesh,
                               const float *local_normals,
                               int32_t *out_labels);

/* Phase 3: boundary cleanup (clear + flood fill + nearest-label + merge) */
static int cleanup_labels(Arena_T arena, const ComponentMesh *mesh,
                          CSR_T adj, KDTree_T tree,
                          int32_t *labels);

/* Phase 4: cut at label boundaries with gap expansion */
static int cut_at_labels(Arena_T arena, const ComponentMesh *mesh,
                         CSR_T adj, const int32_t *labels,
                         int gap_hops,
                         ComponentMesh **out_meshes, size_t *out_count);
```

### Connection to Step 2 (upstream)

After bridge cut, each resulting sub-mesh is re-checked by the Oracle. If
any sub-mesh still has `sheet_count ≥ 2`, it enters `RaycastSep_process()`.
Bridge cut handles thin bridges; Step 3 handles broad overlapping sheets
that don't have a thin neck.

The pipeline driver's dispatch logic:
```
for each component from Step 0:
    sheet_count = Oracle_count(component)
    if sheet_count >= 2:
        pieces = BridgeCut_process(component)
        for each piece in pieces:
            piece_sheets = Oracle_count(piece)
            if piece_sheets >= 2:
                sub_pieces = RaycastSep_process(piece)  // <- Step 3
                add sub_pieces to final list
            else:
                add piece to final list
    else:
        add component to final list
```

### Connection to Step 4 (downstream)

Step 4 (Rebuild Component List) collects all sub-meshes from bridge cut
and raycast separation, assigns new sequential component IDs, and passes
them to Step 5. The sub-meshes from `RaycastSep_process` are
`ComponentMesh` structs with valid vertex/face arrays, ready for Poisson
reconstruction.

---

## External Dependencies

**None.** The raycast separator is pure C. No external libraries are
required. All operations (KD-tree, PCA, BFS, CSR, Union-Find, sorting)
are implemented in `common/`.

The only indirect dependency is the Oracle (Step 1), which is called by
the pipeline driver before entering Step 3 — but the Oracle is not called
from within `raycast_sep.c`.

---

## Test Strategy

### Comparison with Python Reference

The algorithm is **not fully deterministic** due to:
- Floating-point rounding differences (float32 in C vs float64 in Python)
  may shift some vertices across gap boundaries, changing their labels.
- Sort order of vertices with identical w-values may differ.
- BFS visit order depends on CSR neighbor ordering.

**Primary comparison**: for each multi-sheet component that reaches Step 3,
compare:

1. **Number of output sub-meshes**: must match Python exactly. If the
   counts differ, the label assignment or cut logic has a significant bug.

2. **Sub-mesh sizes (vertex count)**: sort by size descending, compare.
   Allow ±5% tolerance on individual piece sizes due to rounding-induced
   label differences at sheet boundaries.

3. **Label distribution**: dump per-vertex labels from both Python and C.
   Compute the label confusion matrix. Diagonal should dominate (>95%
   agreement). Off-diagonal entries indicate systematic labeling errors.

4. **Visual spot-check**: for 5 representative volumes, export labeled
   meshes as colored OBJ (color by label). Visually verify that sheets
   are correctly separated and no sheet is split into fragments.

### Unit Tests

- **Single-sheet component (trivial)**: Oracle says 1 sheet → Step 3
  should never be called. But if it is, it must return the mesh unchanged
  (all vertices get label 1, no cut).

- **Two parallel planes (canonical)**: Two flat planes separated by a gap
  > 8.0 in the w-direction. Expected: 2 labels, 2 sub-meshes, clean cut.

- **Two planes with partial contact**: Planes that overlap in some UV cells
  (gap < 8.0 locally) but separate in most. Tests that flood fill
  propagates labels correctly from the clean region into the ambiguous zone.

- **Three sheets**: A mesh with 3 distinct planes. Expected: 3 labels.
  Tests the gap-labeling logic with gap_threshold = 8.0 producing
  labels 1, 2, 3.

- **Boundary-heavy component**: A mesh where >50% of vertices are within
  25 hops of the boundary. Tests that the boundary-clear + flood-fill
  pipeline recovers labels rather than leaving most vertices unlabeled.

- **Degenerate: < 10 vertices**: Return unchanged without crashing.

- **Degenerate: all vertices on a single plane (gap detection finds no
  gaps)**: Return unchanged (single label → no cut).

- **Degenerate: UV grid size computes to 0 or 1**: Clamp grid_size to
  minimum 2. Return unchanged if labeling fails.

- **Timeout**: Set a 0.001s timeout. Verify that the function returns the
  original mesh unchanged with status -1.

### Regression Tests

Maintain a set of 5–10 multi-sheet components (extracted from test volumes)
where the Python reference produces known label assignments. After every
code change, run the C version on these components and compare sub-mesh
counts and sizes. `make test` catches regressions automatically.

---

## Pitfalls

### From Koenig (C Traps and Pitfalls)

1. **Signed/unsigned mismatch in vertex loops.** Loop counters iterating
   over `nv` (size_t) must be `size_t`. If a signed `int i` is used and
   compared with `nv`, the comparison is unsigned-promoted — a negative `i`
   becomes a huge number. This is the #1 C bug for Python porters.
   [Pitfalls — semantic; c-style-guide §4.4.1]

2. **Off-by-one in BFS hop count.** The 25-hop and 6-hop BFS expansions
   must count hops correctly. An off-by-one (24 or 26 hops) changes which
   vertices are cleared/cut and can cause subtle label differences vs
   Python. Use `depth < max_hops` not `depth <= max_hops` — the seed
   vertex is at depth 0.
   [Elements Ch.4; Pitfalls — semantic; c-style-guide §9.3.3]

3. **Uninitialized labels.** The `labels[N]` array must be zero-initialized
   (unlabeled). Use `Arena_calloc`. An uninitialized label that happens to
   be nonzero will be treated as a valid sheet assignment, causing wrong
   cuts.
   [Elements Ch.4; c-style-guide §1.3.3]

4. **KD-tree neighbor buffer overflow.** If a ball query returns more
   neighbors than `MAX_NEIGHBORS`, and the code doesn't handle this,
   it writes past the buffer. Pre-check the buffer size; truncate and
   document.
   [K&R Ch.5; Pitfalls — semantic; c-style-guide §4.5.1]

5. **Float32 vs float64 gap threshold.** Python uses `gap > 8.0` with
   float64 vertices. C uses `gap > 8.0f` with float32 vertices. Near the
   threshold, rounding may produce different labels. This is acceptable
   if sub-mesh counts match, but document it and test boundary cases.
   [Elements Ch.5; c-style-guide §9.4.1–9.4.2]

6. **Integer overflow in grid index.** `cell_idx = row * grid_size + col`.
   If `grid_size` is large (>46340 on 32-bit int), the product overflows.
   Use `size_t` for the multiplication. In practice, grid_size ≤ 200,
   so this is safe with `int32_t`, but document the bound.
   [Pitfalls — semantic; c-style-guide §9.4.4]

7. **CSR neighbor ordering.** If the CSR is built from faces in a different
   order than Python's adjacency list, BFS visit order changes, flood fill
   produces different label assignments, and sub-mesh vertex counts differ.
   This is the most likely source of "close but not identical" outputs.
   Build the CSR deterministically (process faces in index order, sort
   neighbors per vertex after construction) to match Python's ordering.
   [CLRS Ch.22.1; c-style-guide §7.1.3]

8. **Small-component merge cascading.** If component sizes are not frozen
   before merging, a cascade can occur: component A merges into B, which
   makes B large enough to absorb C, etc. The Python code uses frozen
   sizes. Replicate exactly.
   [PIPELINE_REFERENCE §7, Phase 3 step 5]

### From Linden (Expert C Programming)

9. **sizeof(pointer) vs sizeof(array).** `sizeof(verts)` in a function
   receiving `float *verts` returns 8 (pointer size), not the array size.
   Always pass `nv` explicitly.
   [Linden Ch.4; c-style-guide §4.5.2]

10. **Structure assignment aliases pointers.** `ComponentMesh a = *mesh;`
    makes `a.verts` point to the same memory as `mesh->verts`. If you
    later overwrite `a.verts` (e.g., during sub-mesh extraction), the
    original is NOT affected — but if you modify the *contents* through
    either pointer, both see the change. With arena allocation, this is
    safe (no double-free) but can cause logic bugs.
    [Linden Ch.5; K&R Ch.6; c-style-guide §4.3.4]

### General Step 3 Pitfalls

11. **PCA instability on near-planar or near-linear point clouds.** If all
    vertices lie nearly on a line, the covariance matrix is rank-1 and the
    "smallest eigenvector" is arbitrary. Guard against this: if the smallest
    eigenvalue is < 1e-6 × the largest, fall back to using the global PCA
    normal for that vertex.
    [Abrash rule 87 — precision errors are the enemy of elegant algorithms]

12. **Empty UV cells produce no labels.** If the UV grid is too coarse,
    some cells have no vertices and are skipped. This is correct behavior.
    But if the grid is too fine, most cells have ≤ 1 vertex and gap
    detection becomes meaningless. The adaptive `sqrt(N/50)` formula
    prevents this, but clamp to [4, 200] defensively.

13. **The KD-tree ball query must use squared radius.** Passing `r=10.0`
    instead of `r_sq=100.0` silently returns too few neighbors (only
    vertices within distance √10 ≈ 3.16). This halves the effective
    neighborhood and produces noisy, unstable local normals. Test by
    verifying that the average neighbor count matches the Python reference
    (expected: ~12.5 for N=100K in a 320³ volume).

14. **Boundary vertex detection must count directed half-edges, not
    undirected edges.** An edge (a,b) in face (a,b,c) is a half-edge.
    A boundary edge has exactly 1 half-edge (not 2). If you count
    undirected edges (ignoring direction), non-manifold edges (shared by
    3+ faces) may be misclassified. Use the half-edge count approach:
    for each face (a,b,c), increment counts for (a,b), (b,c), (c,a).
    Edges with count == 1 are boundary.
