# STAGE 1 — Oracle (Sheet Counter)

**Module prefix**: `Oracle_`
**Source files**: `oracle.c`, `oracle.h`
**Dependencies**: `arena.c` (from `common/`)

---

## Purpose

The Oracle determines how many distinct papyrus sheets each mesh component
contains. This integer — the **sheet count** — is the single most important
routing decision in the pipeline. A component with sheet count = 1 skips
bridge cut (Step 2) and raycast separation (Step 3) entirely. A component
with sheet count ≥ 2 enters the expensive topological repair path.

The Oracle is the pipeline's **high-level pruning gate**. Its cost is trivial
(1–2s in Python, sub-second in C) but its effect is enormous: it prevents
the pipeline from spending 15–75s on Edmonds-Karp max-flow or 30–110s on
per-vertex KD-tree queries for components that don't need it. Skiena §9.5's
chessboard covering war story is the exact analog — pruning at a higher
abstraction level (whole-component sheet count) eliminates vastly more
downstream work than any per-vertex or per-edge check could.

**Input**: a single `ComponentMesh` from Step 0 (vertices, faces, PCA normal).
**Output**: integer sheet count (1, 2, 3, …) for that component.

**Why it matters for the metric**:
- A false negative (oracle says 1 sheet when there are 2) means the bridge
  cutter never runs. The inter-wrap bridge survives to the final output.
  VOI and TopoScore are catastrophically degraded.
- A false positive (oracle says 2 sheets when there is 1) wastes time on
  bridge cut but is not catastrophic — bridge cut will find no bridge
  (flow is high, split ratio is bad) and return the mesh unchanged. The
  cost is wasted cycles, not wrong topology.
- Therefore: **bias toward sensitivity over specificity**. It is far better
  to send a single-sheet component into bridge cut (which correctly does
  nothing) than to let a multi-sheet component skip it.

---

## Algorithm

The Oracle works by rasterizing the mesh into a temporary 3D volume, then
casting virtual rays through it along the PCA normal direction. Where a ray
passes through multiple separated surfaces (detected by gaps in the voxel
occupancy), the component contains multiple sheets.

### Phase 1: PCA Normal and Coordinate Setup
1. **Reuse PCA normal** from Step 0's `ComponentMesh.pca_normal[3]`. Do NOT
   recompute — the normal is already computed and sign-corrected.
   [c-style-guide §6.6.3 — store intermediate forms for reuse;
   Abrash Ch.50 — avoid redundant recomputation]

2. **Build orthonormal UV basis** perpendicular to the PCA normal. Use the
   Gram-Schmidt-like construction: pick the axis (x, y, or z) least aligned
   with the normal, cross-product to get U, cross again to get V. The triple
   `(U, V, normal)` forms a right-handed orthonormal frame.

3. **Shift vertices** so all coordinates are non-negative. Compute per-axis
   minimum, subtract from all vertices. This ensures voxel indices are ≥ 0
   when rounding. [PIPELINE_REFERENCE §5, step 2]

### Phase 2: Mesh Rasterization to Temporary Volume
4. **Compute rasterization volume bounds**: from the shifted vertices, find
   per-axis max. The temporary volume dimensions are
   `ceil(max_z)+1 × ceil(max_y)+1 × ceil(max_x)+1`. Allocate from arena
   as `uint8_t *rast_vol`, zero-initialized via `Arena_calloc`.
   [PIPELINE_REFERENCE §5, step 3]

5. **Rasterize each triangle** using barycentric sampling:
   ```
   For each face (a, b, c):
       Compute edge lengths: e0 = |b-a|, e1 = |c-a|, e2 = |c-b|
       max_edge = max(e0, e1, e2)
       subdiv = max(6, ceil(max_edge))    // sample density
       step = 1.0 / subdiv
       For ba = 0, step, 2*step, ... up to 1.0:
           For bb = 0, step, 2*step, ... up to (1.0 - ba):
               bc = 1.0 - ba - bb
               pos = ba * verts[a] + bb * verts[b] + bc * verts[c]
               iz = round(pos.z), iy = round(pos.y), ix = round(pos.x)
               if in bounds: rast_vol[iz*H*W + iy*W + ix] = 1
   ```
   This is a tight arithmetic loop — the Python version uses Numba JIT,
   the C version is a direct nested loop with no function calls in the
   inner body. [PIPELINE_REFERENCE §5, step 3; c-style-guide §6.5.1]

### Phase 3: Ray Casting Along PCA Normal
6. **Project all occupied voxels** to `(u, v, w)` coordinates:
   ```
   For each voxel (z, y, x) where rast_vol[...] == 1:
       pos = (z, y, x)   // float
       u = dot(pos, U_axis)
       v = dot(pos, V_axis)
       w = dot(pos, normal)
   ```
   Store the `(u, v, w)` triples in a temporary array.
   [PIPELINE_REFERENCE §5, step 4]

7. **Bin by (u, v)** on a `GRID_SIZE × GRID_SIZE` grid (GRID_SIZE = 320):
   - Compute `u_min, u_max, v_min, v_max` from the projected points.
   - Cell `(cu, cv)` for a point:
     `cu = floor((u - u_min) / (u_max - u_min) * GRID_SIZE)`, clamped to
     `[0, GRID_SIZE-1]`. Same for `cv`.
   - Use a flat index `cu * GRID_SIZE + cv` to bin.
   [PIPELINE_REFERENCE §5, step 4]

8. **Per-cell gap counting**: for each non-empty cell, sort the w-values.
   Scan sorted w-values; count gaps where `w[i+1] - w[i] > MIN_GAP` (10).
   Sheet count for this cell = 1 + number of gaps.
   [PIPELINE_REFERENCE §5, step 4; Bentley Ch.8 — scanning algorithm with
   running state]

### Phase 4: Aggregation
9. **Consensus voting**: build a histogram of per-cell sheet counts. Find the
   highest sheet count where ≥ `MIN_CELL_RATIO` (10%) of non-empty cells agree.
   Return that as the component's sheet count. If no count reaches 10%, return 1
   (single sheet). [PIPELINE_REFERENCE §5, step 5]

### Refinements from Distillations

- **Prefix-sum optimization** (Bentley Ch.2, c-style-guide §6.6.5): instead
  of sorting w-values per cell, consider precomputing a cumulative histogram
  of w-values over a fixed w-grid (e.g., 1-voxel resolution along the normal).
  Then counting gaps per UV cell is an O(1) lookup per w-bin instead of an
  O(K log K) sort per cell. Profile both approaches — for the typical case
  (few voxels per cell), the sort is fast enough. The prefix-sum optimization
  is a Level 2 (Abrash) improvement worth trying only if profiling shows the
  per-cell sort dominates.

- **Resolution tuning** (sanglard_rules_part2 section N "Console Port Lessons", rule 41): the 320×320 UV grid may be
  overkill. A 160×160 grid has 4× fewer cells and may detect the same
  multi-sheet components. Profile on training volumes with 160, 240, 320 grids
  and compare sheet-count accuracy. If 160 is as accurate, use it — it halves
  the sort work and reduces the temporary allocation.

---

## Book References

### Core Algorithm
- **PIPELINE_REFERENCE §5**: the complete Oracle algorithm specification
  including all constants (GRID_SIZE=320, MIN_GAP=10, MIN_CELL_RATIO=0.10).
- **Bentley Ch.8** (scanning algorithms): the per-cell gap counting is a
  classic scanning algorithm — sort once, scan linearly maintaining a running
  gap state. This is the O(N) pattern from Bentley's maximum-subarray analysis:
  maintain running state, never rescan.
- **Bentley Ch.2** (prefix sums): the cumulative-histogram optimization for
  gap counting is a direct application of Bentley's prefix-sum pattern —
  precompute a summary, answer queries in O(1).

### Performance and Pruning Role
- **Skiena Ch.9.5** (chessboard covering war story): the Oracle is the "weak
  attack" pruning at a higher abstraction level. It tests the whole component
  for sheet count before any per-vertex/per-edge work is done. This is
  orders of magnitude cheaper than the per-edge flow computation it gates.
- **Skiena Ch.9.4** (Sudoku): look-ahead pruning + most-constrained-first
  together gave >1000× speedup. The Oracle is the look-ahead: check the
  child before recursing. Combine with process-largest-first ordering for
  multiplicative gain.
- **Abrash Ch.51** (backface removal): "instant rejection of half the work
  is the best optimization you'll never notice." The Oracle rejects single-
  sheet components from the expensive path — potentially 50%+ of components.
- **Abrash Ch.86** (match sophistication to scale): small or simple components
  don't need the full topological repair pipeline. The Oracle is the dispatch
  that routes each component through minimum necessary processing.

### Rasterization Loop
- **Abrash Ch.38–39** (polygon fill): the rasterization phase is a direct
  analog of Abrash's scan-conversion — iterate over barycentric coordinates,
  compute 3D position, round to grid. The inner loop is pure arithmetic with
  no branches (just bounds-check and write). This is Level 1 optimization:
  same algorithm as Python, just without interpreter overhead.
- **c-style-guide §6.5.1**: barycentric subdivision pattern for mesh
  rasterization. Determine sample density from max edge length.
- **Abrash Ch.19** (branch prediction): the rasterization inner loop has a
  predictable branch (bounds check — almost always true). Don't unroll. The
  branch predictor handles it. Do NOT add complexity to eliminate this branch.

### PCA and Coordinate Transform
- **Abrash Ch.61** (frame of reference): the UV basis construction is a
  coordinate frame change. The dot product is the universal workhorse — three
  dot products transform a 3D point into the (u, v, w) frame. If the PCA
  normal looks wrong, the frame is wrong, and all raycasting results are
  garbage. Verify the frame visually on a known test case.
- **Abrash Ch.74**: "the dot product is the universal workhorse of 3D
  geometry." Every operation in the Oracle — basis construction, projection,
  gap measurement — reduces to dot products.

### Memory Management
- **Hanson Ch.6**: the temporary rasterization volume is scratch memory.
  Allocate from the arena with a save/restore mark. When the Oracle finishes,
  restore the mark to reclaim the volume. This is the DOOM-style
  `PU_LEVEL`/`PU_CACHE` pattern from Sanglard Ch.2.
- **Hanson Ch.16** (save/restore): `Arena_save` before allocating the rasterization
  volume, `Arena_restore` after extracting the sheet count. Clean, zero-leak.

### Sorting
- **Bentley Ch.8** (Column 8): the per-cell sort of w-values is on small arrays
  (typically 1–50 elements). Use insertion sort for N ≤ ~20, `qsort` for
  larger. Or since the per-cell count is almost always small, just use
  insertion sort unconditionally — it's optimal for small N and has no
  function-call overhead.
- **Skiena Ch.3.6** (bounded-height PQ): since sheet counts are small integers
  (typically 1–5), the aggregation histogram is a tiny array indexed by sheet
  count. No need for any sorting structure — direct array indexing.

---

## Performance Budget

**Total pipeline budget**: ~180s per volume on Kaggle (2 cores @ 2.2 GHz). This 180s is the C preprocessing share of the ~270s/volume total (the remaining ~90s is nnUNet inference).
**Step 1 budget**: ~2s (from the Python profile: 1–2s, bottleneck is
Numba-JIT rasterization). In C, expect significant speedup over the JIT —
the rasterization loop is pure arithmetic with perfect branch prediction.

### Back-of-Envelope for Kaggle (2 cores @ 2.2 GHz)

The Oracle runs once per component, up to 20 components. The cost is
dominated by the largest component.

| Sub-step | Work | Est. Time |
|---|---|---|
| UV basis construction | O(1) — 2 cross products | negligible |
| Vertex shift | N vertices × 3 subtracts | ~0.001s |
| Rasterize (per component) | M faces × subdiv² samples × 6 ops | see below |
| Voxel projection | V_occupied × 3 dot products | see below |
| UV binning | V_occupied × 2 divides + array write | see below |
| Per-cell sort | Σ(K_cell × log K_cell) over all cells | see below |
| Aggregation | O(n_cells) scan | negligible |

**Rasterization detail**: a typical component has M ≈ 100K faces (after
backface cull). Average `subdiv` ≈ 8 (for faces with max edge length ~8
voxels). Barycentric samples per face ≈ subdiv²/2 ≈ 32. Total samples:
100K × 32 = 3.2M. Per sample: 3 multiplies + 3 adds (barycentric interp)
+ 3 rounds + 1 bounds check + 1 array write ≈ 15 ops. Total: 3.2M × 15
= 48M ops × ~0.5 ns/op = ~24 ms. **Rasterization is cheap.**

**Projection detail**: occupied voxels V_occ ≈ 100K–500K (depending on
component size and surface density). Per voxel: 3 dot products (9 muls +
6 adds) + 2 divides (for UV binning) ≈ 20 ops. Total: 500K × 20 = 10M ops
× ~0.5 ns = ~5 ms. **Projection is cheap.**

**Sorting detail**: 320² = 102,400 cells. Average occupancy ≈ 5 voxels/cell
(500K voxels / 100K non-empty cells). Insertion sort on 5 elements: ~10
comparisons. Total: 100K cells × 10 = 1M ops × ~1 ns = ~1 ms. Even with
`qsort` overhead, this is negligible.

**Total per component**: ~30–50 ms. **For 20 components: ~0.6–1.0s.**

**Memory**: the temporary rasterization volume for a 320³ bounding box is
32.8 MB (uint8). The (u,v,w) projection array is up to 500K × 12 bytes =
6 MB. The UV grid bin structure (CSR or array-of-arrays) is ~1 MB. Total
scratch: ~40 MB. Fits easily in the arena; reclaimed after Oracle returns.

### Scaling: 7950X → Kaggle
- The Oracle is pure computation (no I/O, no external calls).
- Single-core 7950X is ~2.6× faster than Kaggle.
- Expect ~0.2–0.4s on 7950X, ~0.6–1.0s on Kaggle.
- No parallelism needed within the Oracle — it's already fast enough.
  If you want, parallelize across components (run Oracle on 2 components
  simultaneously), but the per-component cost is so low it's not worth
  the thread-management overhead on 2 cores.

**Verdict**: Step 1 fits trivially within budget. This is the cheapest
computational step in the pipeline. Spend zero optimization effort here.
The Oracle's value is not in its own speed but in the downstream work it
prevents.

---

## Key C Rules

**§1.2 Arena Lifetime** (c-style-guide): the rasterization volume, projection
arrays, and UV grid bins are all scratch data. Use `Arena_save` before
Oracle entry, allocate everything, extract the integer sheet count, then
`Arena_restore` to reclaim. The only output that survives is the integer.
[Hanson Ch.6; Sanglard Ch.2]

**§1.6 Stack Hazards**: the rasterization volume for a 320³ component could
be ~35 MB. Never allocate on the stack. Always arena.

**§4.6.1 Dense Arrays**: vertices as `float verts[N*3]`, the rasterization
volume as `uint8_t rast_vol[D*H*W]`, the UV bins as flat arrays. No linked
structures, no per-cell malloc.

**§5.1 Feasibility Check**: done above. Step 1 fits with >100× headroom.
No optimization needed. Write the simplest correct code.

**§5.2.1 Correctness First**: "The first principle of optimization is don't."
The Oracle is already cheap. Prioritize correctness (matching Python output)
over speed. A wrong sheet count that passes a single-sheet bridge into
Step 2–3 is harmless (wasted time). A wrong sheet count that prevents a
multi-sheet bridge from entering Step 2 is catastrophic (wrong topology).

**§5.4.6 Reducing Space**: use `uint8_t` for the rasterization volume (binary:
0 or 1). Not `int32_t`. The 32.8 MB volume fits in L3 at uint8; it would be
131 MB at int32 and spill to DRAM.

**§5.5.1 No Loop Unrolling**: the rasterization inner loop is already
branch-predicted perfectly. Don't unroll. Don't add SIMD. Write clean scalar C.

**§6.1.1 Volume Layout**: `rast_vol[z*H*W + y*W + x]`. Same convention as
everywhere else. Document it.

**§6.6.3 Store Intermediate Forms**: the PCA normal from Step 0 is reused.
The UV basis vectors computed here may also be useful for Step 3 (raycast
separator) if the same component enters it. Consider storing them in the
`ComponentMesh` struct for reuse. [Abrash Ch.50]

**§9.3.2 Graceful Empty Cases**: zero vertices → return sheet_count=1. Zero
occupied voxels after rasterization → return 1. Zero non-empty UV cells →
return 1. Never crash on degenerate input.

**§9.4.3 Integer Arithmetic for Counting**: sheet counts, cell counts, gap
counts — all integers. The `w` values and UV coordinates are float. The gap
comparison `w[i+1] - w[i] > MIN_GAP` uses float comparison with integer
threshold (10.0). No exact-equality comparison needed.

**§9.5 UB Avoidance**: the rasterization volume index
`(size_t)iz * H * W + (size_t)iy * W + (size_t)ix` must use `size_t` to
avoid signed overflow on intermediate products. Bounds-check before write:
`if (iz >= 0 && iz < (int)D && iy >= 0 && iy < (int)H && ix >= 0 && ix < (int)W)`.

---

## Data Structures

### From `common/` (shared modules)

| Module | Type / Function | Used For |
|---|---|---|
| `arena.h` | `Arena_T`, `Arena_save`, `Arena_restore` | Scratch allocation |
| `mesh_extract.h` | `ComponentMesh` (from Step 0) | Input mesh + PCA normal |

### Created internally by Step 1

**Temporary rasterization volume**:
```c
uint8_t *rast_vol;   /* [D_rast * H_rast * W_rast], arena-allocated, zero-init */
size_t   D_rast, H_rast, W_rast;  /* dimensions of the temporary volume */
```

**Projected voxel list** (one of two approaches):

*Approach A* — collect then bin:
```c
typedef struct {
    float u, v, w;
} UVW_Point;

UVW_Point *proj_points;   /* [n_occupied], arena-allocated */
size_t     n_occupied;
```

*Approach B* — bin directly during projection (avoids the intermediate array):
scan the rasterization volume linearly. For each occupied voxel, compute (u,v),
determine cell index, append w to that cell's list. Store cells as a CSR-like
structure:
```c
int32_t *cell_offsets;  /* [GRID_SIZE * GRID_SIZE + 1] */
float   *cell_w_values; /* [n_occupied], sorted per cell */
```
Build in two passes: pass 1 counts per-cell occupancy → prefix sum for offsets.
Pass 2 writes w-values using a write cursor per cell. Then sort w-values within
each cell. This is the standard CSR construction pattern.
[c-style-guide §7.1.3; Sanglard Ch.5 §6.4.8]

**Approach B is preferred** — it avoids the intermediate UVW array and produces
the binned structure directly. The two-pass CSR construction is O(V_occupied)
with excellent cache behavior (sequential scan of the volume, sequential writes
to the cell arrays).

**Aggregation histogram**:
```c
int32_t sheet_hist[MAX_SHEET_COUNT + 1];  /* on the stack — MAX_SHEET_COUNT ≈ 10 */
```
This is tiny (≤ 44 bytes). Stack allocation is fine.

---

## Interfaces

### Input
```c
/*
 * Oracle_count_sheets — Step 1 entry point.
 *
 * Given a component mesh (with PCA normal already computed by Step 0),
 * rasterize it into a temporary volume, cast rays along the PCA normal
 * on a UV grid, and count sheets by gap detection.
 *
 * Returns the sheet count (≥ 1). Returns 1 for degenerate inputs.
 *
 * All scratch memory is arena-allocated and reclaimed internally
 * via Arena_save/Arena_restore. The arena state is unchanged on return.
 */
int Oracle_count_sheets(Arena_T              arena,
                        const ComponentMesh *mesh,
                        int                  grid_size);  /* typically 320 */
```

### Output

A single `int` — the sheet count. No heap allocation survives this call.

### Connection to Step 0 (upstream)

The Oracle receives `ComponentMesh` structs from Step 0. It reads:
- `mesh->verts`, `mesh->n_verts` — vertex positions
- `mesh->faces`, `mesh->n_faces` — triangle indices
- `mesh->pca_normal` — the PCA sheet normal (precomputed, sign-corrected)

### Connection to Step 2 (downstream)

The caller (the pipeline driver) uses the Oracle's result to decide routing:
```c
for (size_t i = 0; i < n_meshes; i++) {
    int sheets = Oracle_count_sheets(arena, &meshes[i], GRID_SIZE);
    if (sheets >= 2) {
        /* Enter bridge cut (Step 2) */
        BridgeCut_process(arena, &meshes[i], ...);
    }
    /* else: single sheet, skip to Step 4 (rebuild) → Step 5 (Poisson) */
}
```

The sheet count is also used as the recursive bridge-cut termination criterion:
after each split, re-run the Oracle on each piece. If Oracle returns 1, stop
recursing. [Skiena Ch.9.4 — look-ahead pruning]

---

## External Dependencies

**None.** The Oracle is pure C with no external library calls. No libtiff,
no PoissonRecon, no Triangle, no Clipper2. It uses only the arena allocator
and basic floating-point arithmetic.

This is by design — the Oracle must be fast, simple, and self-contained.
It's the cheapest step in the pipeline and should have zero dependency
complexity.

---

## Test Strategy

### Comparison with Python Reference

The Python Oracle is deterministic for a given mesh. The C Oracle should
produce **identical sheet counts** for every component on every test volume.
Sheet count is an integer — there is no tolerance. If they disagree, one of
them is wrong.

However, there are subtle sources of divergence:

1. **Rasterization sampling order**: the Python Numba kernel and the C loop
   may iterate barycentric coordinates in slightly different order, producing
   slightly different sets of occupied voxels (due to float rounding of
   `round(pos.z)` etc.). This is acceptable as long as the gap detection
   produces the same sheet count. If it doesn't, investigate the specific
   divergent cell.

2. **Sort stability**: the per-cell sort of w-values. If two w-values are
   identical (two voxels at the same depth along the normal), their relative
   order doesn't matter for gap counting. But if a sort is unstable, the gap
   at the boundary between identical values might be counted differently.
   Use a stable sort or (better) deduplicate identical w-values within each
   cell before gap counting.

3. **PCA normal sign**: if the C PCA normal differs from Python's by a sign
   flip, the w-projection is negated, the UV basis is different, and the
   entire ray casting grid is rotated. Sheet counts may still match (the gap
   structure is symmetric) but if they don't, the PCA sign is the first thing
   to check. **Verify `dot(c_normal, python_normal) > 0.99` for every
   component before comparing sheet counts.**

### Unit Tests

- **Known single-sheet mesh**: a flat disk or a single-sheet plane segment.
  Oracle must return 1.
- **Known two-sheet mesh**: two parallel planes separated by 15 voxels
  (> MIN_GAP=10). Oracle must return 2.
- **Known three-sheet mesh**: three parallel planes. Oracle must return 3.
- **Borderline gap**: two planes separated by exactly 10 voxels.
  `w[i+1] - w[i] > MIN_GAP` with MIN_GAP=10 means a gap of exactly 10 is
  NOT a gap (> 10, not ≥ 10). Verify this boundary. A gap of 11 IS a gap.
  [Elements Ch.4 — branch the right way on equality]
- **Empty mesh**: 0 vertices, 0 faces. Must return 1, not crash.
- **Single-face mesh**: 1 triangle. Must return 1.
- **All voxels in one UV cell**: a very small mesh. Must still count correctly.
- **No cells reaching 10% threshold**: should return 1.

### Integration Tests

Run the full pipeline on 5–10 representative training volumes. For each
component, compare Oracle sheet counts between Python and C. Log any
disagreements with the component ID, vertex count, and the actual vs expected
sheet count. Zero disagreements is the target.

### Regression Suite

- `make test_step1` runs all Oracle tests in <5s on the 7950X.
- Add a new test for every bug found (growing monotonically).
- Include the training volumes that have the most components (stress the
  per-component loop) and the largest components (stress the rasterization).

---

## Pitfalls

### The > vs ≥ Gap Threshold
The gap detection uses `w[i+1] - w[i] > MIN_GAP`, which is strict inequality.
A gap of exactly `MIN_GAP` (10.0) voxels is NOT counted as a sheet separator.
If you write `>=` instead of `>`, you'll count more gaps and inflate sheet
counts, sending single-sheet components into unnecessary (but harmless) bridge
cut. If you write `>` when Python uses `>=` (or vice versa), sheet counts
diverge. **Match the Python implementation exactly.** Write a unit test for
the boundary case: gap = 9.99 (no), gap = 10.0 (depends on > vs >=),
gap = 10.01 (yes).
[Elements Ch.4 — "take care to branch the right way on equality";
c-style-guide §9.3.4]

### Rasterization Volume Dimensions
The temporary volume dimensions come from the shifted vertex bounding box:
`D_rast = ceil(max_z) + 1`, etc. If you forget the `+1`, vertices exactly at
`max_z` round to an out-of-bounds index. If you use `floor` instead of `ceil`,
the volume is too small. If you don't shift vertices first, negative coordinates
produce negative indices — array underflow, silent memory corruption.
[c-style-guide §9.3.3 — off-by-one; Pitfalls — semantic]

### Barycentric Sampling Step Size
The subdivision level `subdiv = max(6, ceil(max_edge))` determines sample
density. If `subdiv` is too low, large triangles have gaps in the
rasterization volume, creating false gaps in the ray cast. If `subdiv` is
too high, you waste cycles. The Python code uses `ceil(max_edge)` — match
this exactly. Edge length is in voxel units, so a 10-voxel-long edge needs
at least 10 samples to avoid gaps. The `max(6, ...)` floor ensures small
triangles still get adequate coverage.

### UV Grid Clamp
When computing the UV cell index:
`cu = floor((u - u_min) / (u_range) * GRID_SIZE)`, the result can be
exactly `GRID_SIZE` when `u == u_max` (due to floating point). Clamp to
`[0, GRID_SIZE-1]`. Without the clamp, you get an out-of-bounds array write.
This is a classic off-by-one at the upper boundary.
[c-style-guide §9.3.3; Pitfalls — semantic]

### Division by Zero in UV Range
If all vertices project to the same u-value (degenerate: mesh is a line
perpendicular to U axis), then `u_range = u_max - u_min = 0.0`, and the UV
binning division produces NaN or infinity. Guard against this:
`if (u_range < epsilon) u_range = 1.0;`. The sheet count for such a
degenerate component is 1 (it can't have multiple sheets if it has no extent
in one UV direction).
[c-style-guide §9.4.2 — never compare float for exact equality;
guard divisions against zero]

### Rasterization Volume Memory for Large Components
A component whose bounding box spans the full 320³ volume requires a 32.8 MB
temporary volume. This is fine. But if a volume is 512³ and a component spans
the full extent, the temporary volume is 134 MB. This still fits in the arena
but is a large allocation. The arena must be sized for the worst case. If
memory is tight, consider rasterizing at half resolution (every other voxel)
and doubling MIN_GAP to 20 — same detection quality, 8× less memory.
[c-style-guide §1.1.5 — arena sizing; sanglard_rules_part2 section N "Console Port Lessons", rule 41 — reduce resolution
before removing functionality]

### Integer Overflow in Volume Index
The rasterization volume index `iz * H_rast * W_rast + iy * W_rast + ix`
can overflow `int32_t` for large volumes. For a 400³ volume:
400 × 400 × 400 = 64M — fits int32 (max 2.1B) but compute the intermediate
product `iz * H_rast * W_rast` in `size_t` to be safe. Use `size_t` for all
index arithmetic as mandated by c-style-guide §9.4.4.

### The 10% Consensus Threshold
`MIN_CELL_RATIO = 0.10` means at least 10% of non-empty cells must agree on a
sheet count for it to be returned. If the mesh is very irregular (surface folds,
non-planar), many cells may disagree. The aggregation takes the **highest**
sheet count meeting the 10% threshold, not the most common count. This biases
toward detecting multi-sheet components (the safe direction). Make sure the
aggregation logic matches Python: iterate from highest possible sheet count
downward, return the first one meeting the threshold.

### Coordinate Convention Consistency
Step 0 produces vertices in `(z, y, x)` voxel coordinates. The Oracle's PCA
normal is also in `(z, y, x)`. The UV basis must be in the same coordinate
system. If you accidentally treat coordinates as `(x, y, z)` anywhere in the
Oracle, the projections are wrong and sheet counts are garbage. The symptom:
correct sheet counts for roughly cubic components, wrong counts for
elongated ones.
[PIPELINE_REFERENCE §16 — coordinate convention]
