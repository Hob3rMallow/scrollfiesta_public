# STAGE 2 — Bridge Cut

**Module prefix**: `BridgeCut_`
**Source files**: `bridge_cut.c`, `bridge_cut.h`
**Dependencies**: `arena.c`, `csr.c`, `kdtree.c`, `bfs.c`, `pca.c` (from `common/`)

---

## Purpose

Bridge cut is the pipeline's core topological repair step. It detects and
severs thin "bridges" — spurious connections between adjacent scroll wraps
produced by the nnUNet segmentation. These bridges merge distinct papyrus
sheets into a single connected component, catastrophically degrading VOI
(which penalizes mergers) and TopoScore (which penalizes wrong Betti numbers).

The Oracle (Step 1) identifies components with ≥ 2 sheets. Bridge cut
attempts to separate them by finding the minimum vertex cut — the thinnest
neck connecting two sides of the mesh. If the cut is thin relative to the
component's extent (the `is_bridge` heuristic), the mesh is split at the cut
with a geometric exclusion zone to prevent ragged edges.

The process is **recursive**: after splitting, each piece is re-evaluated
(oracle + bridge detection) and split again if needed, up to depth 10.

**Input**: a multi-sheet `ComponentMesh` from Step 0 (oracle sheet count ≥ 2).
**Output**: a list of single-sheet (or unsplittable) sub-meshes.

**Why it matters for the metric**:
- **VOI (35% of score)**: every inter-wrap merger that survives to the output
  adds to VOI. Bridge cut is the *only* mechanism for separating bridges that
  are too narrow for the raycast separator (Step 3) to detect. If bridge cut
  fails silently, VOI degrades catastrophically.
- **TopoScore (30% of score)**: bridges create spurious handles (k=1 Betti
  defects). Cutting them restores correct topology.
- **SurfaceDice (35%)**: cutting too aggressively (removing real surface) or
  leaving ragged edges degrades surface proximity. The CUT_GAP_DEPTH
  exclusion zone mitigates this.
- **Failure mode**: a crash or assertion failure in bridge cut wastes the
  volume. The outer loop falls back to the raw prediction. Design so that
  every degenerate case returns the mesh unchanged rather than crashing.

---

## Algorithm

The algorithm has three layers: bridge detection (max-flow), mesh splitting
(KD-tree gap exclusion), and recursive decomposition (backtracking with oracle
pruning).

### Layer 1: Bridge Detection (`detect_bridge`)

Given a component mesh (vertices, faces):

1. **PCA normal and height.** Compute the PCA normal (smallest-eigenvalue
   eigenvector of the vertex covariance matrix). Project all vertices onto
   this normal → scalar value per vertex. `height` = max projection − min
   projection. This is the extent of the mesh along the sheet-stacking
   direction.
   [PIPELINE_REFERENCE §6; reuse PCA from Step 0 if available]

2. **Source/sink seed selection.** Source seed = vertex with minimum
   projection (one extreme of the mesh along the normal). Sink seed = vertex
   with maximum projection (other extreme).
   [PIPELINE_REFERENCE §6, step 3]

3. **k-ring BFS expansion.** Expand source seed to `SEED_RING=20` geodesic
   hops using BFS on the mesh adjacency graph. Same for sink seed. This
   creates two "anchor" regions on opposite sides of the component. Remove
   any vertices that appear in both sets (overlap). If either set is empty
   after overlap removal, the mesh is too small or too connected — return
   NO BRIDGE.
   [PIPELINE_REFERENCE §6, step 5–6; CLRS §22.2 for BFS]

4. **Build vertex-split flow network.** The key transformation: convert the
   vertex connectivity problem into an edge connectivity problem solvable by
   standard max-flow.

   For each mesh vertex v (N total):
   - Create `v_in` (node index `2*v`) and `v_out` (node index `2*v + 1`).
   - Internal edge: `(v_in → v_out)` with capacity 1.
     Exception: source/sink seed vertices get capacity INF.
   - Reverse edge: `(v_out → v_in)` with capacity 0.

   For each mesh edge (u, v) (undirected, from face adjacency):
   - Forward: `(u_out → v_in)` with capacity INF. Reverse: `(v_in → u_out)` cap 0.
   - Forward: `(v_out → u_in)` with capacity INF. Reverse: `(u_in → v_out)` cap 0.

   Super-source S (node `2*N`): connects to all source-seed `v_in` with capacity INF.
   Super-sink T (node `2*N + 1`): connects from all sink-seed `v_out` with capacity INF.

   Total nodes: `2*N + 2`. Total directed edges (with reverses):
   `2*(N + 4*E_undirected + |src_seeds| + |sink_seeds|)`. For N=100K, E=300K:
   ~2.6M edges.
   [CLRS §26 Problem 26-1; Skiena §18.8; c-style-guide §7.5.1–7.5.2]

5. **Edmonds-Karp max-flow.** Standard BFS-based augmenting-path algorithm.
   ```
   total_flow = 0
   while total_flow < MAX_FLOW_LIMIT (100):
       BFS from S in residual graph (edges with cap[e] - flow[e] > 0)
       if sink T not reached: break   // max-flow found
       trace path S→T via parent_edge[], find bottleneck
       augment: for each edge e on path:
           flow[e] += bottleneck
           flow[reverse[e]] -= bottleneck
       total_flow += bottleneck
   ```
   The `MAX_FLOW_LIMIT = 100` early exit is **essential**. Without it,
   worst-case Edmonds-Karp is O(VE²) — intractable for 400K nodes × 2.6M
   edges. With the cap, we do at most 100 BFS passes through a 2.6M-edge
   graph, giving ~260M operations — feasible.
   [CLRS §26.2 Theorem 26.8; c-style-guide §7.4.1–7.4.6]

6. **Min-cut extraction.** After Edmonds-Karp terminates, the last failed
   BFS's visited set IS the S-side of the min-cut. No separate BFS needed —
   save the visited flags from the final iteration. Mark all visited nodes
   as reachable, all others as unreachable. Map back to original vertex
   indices: vertex v is reachable iff `v_in` (node `2*v`) is reachable.
   [CLRS §26.2 Theorem 26.6; c-style-guide §7.6.1]

7. **Compute heuristic metrics.**
   - `flow` = total_flow (the min-cut size).
   - `split_ratio` = min(|reachable|, |unreachable|) / N.
   - `flow_height_ratio` = flow / height.

8. **Bridge decision (`is_bridge`).**
   ```
   if split_ratio < 0.20:          → NO BRIDGE (too asymmetric — one side
                                      is tiny, this is a fragment not a bridge)
   if flow_height_ratio ≤ 0.36:    → BRIDGE (thin neck relative to extent)
   if flow_height_ratio ≤ 1.3
      AND split_ratio ≥ 0.40
      AND height ≥ 20:             → BRIDGE (gray zone — moderately thin,
                                      but balanced and tall enough)
   else:                           → NO BRIDGE
   ```
   [PIPELINE_REFERENCE §6, is_bridge rule]

### Layer 2: Mesh Splitting (`split_mesh`)

When a bridge is detected:

1. **Find boundary vertices.** Iterate mesh edges; a vertex v is a boundary
   vertex if `reachable[v] != reachable[neighbor]` for any neighbor. Collect
   boundary vertex positions.
   [PIPELINE_REFERENCE §6, split_mesh step 1]

2. **Build KD-tree on boundary vertices.** The boundary set is typically
   small (~100–2K vertices). A KD-tree on this set enables fast nearest-
   boundary-distance queries. Alternatively, brute-force is acceptable for
   |boundary| < 2000 (2000² = 4M comparisons at ~1ns = 4ms).
   [PIPELINE_REFERENCE §6, split_mesh step 2; Skiena §15.6]

3. **Compute exclusion zone.** For every mesh vertex, query its distance to
   the nearest boundary vertex. If `dist < CUT_GAP_DEPTH` (7.0 voxel units),
   mark the vertex as excluded. This creates a geometric gap around the cut
   so that the two resulting sub-meshes don't share ragged edges.
   [PIPELINE_REFERENCE §6, split_mesh step 3–4]

4. **Partition faces.** For each face (a, b, c):
   - If any vertex is excluded: skip the face (it falls in the gap).
   - If all three vertices are reachable: assign to source side.
   - If all three vertices are unreachable: assign to sink side.
   - Mixed (some reachable, some unreachable, none excluded): assign to the
     side that has the majority of vertices. (This shouldn't happen if the
     exclusion zone is wide enough, but handle it defensively.)
   [PIPELINE_REFERENCE §6, split_mesh step 5]

5. **Extract sub-meshes.** For each side, collect the relevant faces, build
   a vertex reindex map (old index → new index), emit reindexed faces and
   the corresponding vertex subset. Run mesh connected-components (Union-Find)
   on each side; discard fragments with < 10 faces.
   [PIPELINE_REFERENCE §6, split_mesh step 6]

### Layer 3: Recursive Decomposition (`process_piece_recursive`)

```
process_piece_recursive(mesh, depth, timeout_at):
    if depth > MAX_RECURSION_DEPTH (10): return [mesh]
    if mesh.n_faces < 20: return [mesh]
    if clock() > timeout_at: return [mesh]  // timeout safety

    sheet_count = Oracle_count_sheets(mesh)
    if sheet_count <= 1: return [mesh]      // look-ahead prune

    (flow, reachable, split_ratio, height) = detect_bridge(mesh)
    if not is_bridge(flow, height, split_ratio): return [mesh]

    (side_A, side_B) = split_mesh(mesh, reachable)

    if side_A.n_faces < 10 or side_B.n_faces < 10:
        return [mesh]  // degenerate split, keep original

    results = []
    // Process larger piece first (most-constrained-first heuristic)
    if side_A.n_faces >= side_B.n_faces:
        results += process_piece_recursive(side_A, depth+1, timeout_at)
        results += process_piece_recursive(side_B, depth+1, timeout_at)
    else:
        results += process_piece_recursive(side_B, depth+1, timeout_at)
        results += process_piece_recursive(side_A, depth+1, timeout_at)

    return results
```
[Skiena §9.1 backtracking template; Skiena §9.4 look-ahead pruning;
Skiena §9.4 most-constrained-first; Bentley Ch.8 divide-and-conquer]

### Refinements from Distillations

- **Articulation vertex pre-check** (c-style-guide §7.8.1; Skiena §7.9.2):
  before running the expensive max-flow, check for articulation vertices
  with a linear-time DFS O(V+E). If one exists, cutting there is free — no
  flow computation needed. This is a cheap pre-filter that can save ~1.8s
  per component in the best case. Implementation: Tarjan's algorithm
  (single DFS with `low[]` array).

- **Oracle caching** (Skiena §9.3): if a sub-mesh was already tested by the
  oracle at a previous recursion level and found to be single-sheet, cache
  that result (keyed by vertex count + hash of face indices). Don't re-run
  the oracle on it. The oracle is cheap (~30ms), but over 10 recursion
  levels × 2 pieces × 20 components, the calls add up.

- **Timeout discipline**: use `clock_gettime(CLOCK_MONOTONIC)` for the
  timeout check. Check at the top of each recursive call, not inside
  Edmonds-Karp (which would add overhead to the hot inner loop). If timeout
  fires, return the mesh unchanged — partial results are worse than the
  original unsplit mesh.

---

## Book References

### Max-Flow Algorithm
- **CLRS Ch.26 §26.1–26.2**: flow network definitions, Ford-Fulkerson method,
  Edmonds-Karp algorithm (BFS augmenting paths), O(VE²) bound (Theorem 26.8),
  max-flow min-cut theorem (Theorem 26.6). The full pseudocode in the
  distillation (clrs_rules, rule 28) maps directly to C.
- **CLRS Ch.26 Problem 26-1**: vertex splitting — the construction that
  transforms vertex capacity constraints into edge capacity constraints. This
  is the theoretical foundation for our vertex-split flow graph.
- **CLRS Ch.26.3 Theorem 26.10**: integrality theorem — since all capacities
  are integer, the max-flow is integer. Each augmenting path adds ≥ 1 unit.
  Number of BFS passes = flow value = min vertex cut size.
- **Skiena Ch.8.5 §8.5.2**: practical Edmonds-Karp implementation. The
  `augment_path` code, `valid_edge` predicate, and residual graph management.
  The #1 implementation bug: forgetting to update the reverse edge.
- **Skiena Ch.18.8**: vertex connectivity requires vertex splitting to reduce
  to edge connectivity / max-flow. Confirms the construction. Also:
  biconnected components (articulation vertices) can be found in linear time
  as a cheap pre-filter.
- **Skiena Ch.18.9**: practical advice on flow algorithms. Unit-capacity
  networks have specialized faster algorithms. Preflow-push is theoretically
  faster but Edmonds-Karp with flow cap is simpler and sufficient here.

### Graph Representation
- **CLRS Ch.22.1**: adjacency-list representation. CSR is the cache-friendly
  analog for static graphs.
- **Skiena Ch.7.2**: "adjacency lists are the right data structure for most
  applications of graphs." Never use adjacency matrices for our sparse meshes.
- **Skiena Ch.15.4**: packed arrays for static graphs — exactly CSR.
  4× speedup over LEDA's general graph type just from compactness.
- **c-style-guide §7.1**: CSR construction (two-pass: degree count → prefix
  sum → fill). Pre-allocate forward + reverse edges together.

### BFS
- **CLRS Ch.22.2**: BFS algorithm, O(V+E), shortest-path guarantee. Flat
  queue with head/tail indices. Eliminate the color array — `dist[v] == UINT32_MAX`
  means unvisited.
- **Skiena Ch.7.6**: customize BFS via callbacks/early-termination, not by
  rewriting. Same BFS core for k-ring expansion and Edmonds-Karp path finding.

### Recursive Decomposition
- **Skiena Ch.9.1 §9.1–9.3**: the five-callback backtracking template. Maps
  directly onto recursive bridge cut: `is_a_solution` = oracle says 1 sheet;
  `construct_candidates` = Edmonds-Karp min-cut; `make_move` = split mesh.
  DFS recursion uses O(depth) space vs O(width) for BFS — critical since
  depth ≤ 10 but width could be 2^10 = 1024.
- **Skiena Ch.9.4**: look-ahead pruning (oracle check before recursing) +
  most-constrained-first (process larger piece first) together gave >1000×
  speedup in the Sudoku experiment.
- **Skiena Ch.9.5**: "chessboard covering" war story — pruning at a higher
  abstraction level (whole-component oracle check) eliminates vastly more
  work than per-vertex/per-edge checks. The oracle is the "weak attack."
- **Bentley Ch.8**: divide-and-conquer yields O(N log K) total work when
  subproblems don't overlap and depth is bounded.

### KD-Tree for Gap Exclusion
- **Skiena Ch.15.6**: KD-tree for nearest-neighbor queries. Build on boundary
  vertices (~1K points), query for all mesh vertices. Construction O(B log B),
  each query O(√B). Or brute-force O(B) per query if B < 2000.
- **c-style-guide §6.4.1–6.4.6**: flat-array KD-tree, heap layout, compare
  squared distances, prune subtrees by bounding box.
- **Bentley Ch.9**: replace sqrt with squared-distance comparison.

### Performance and Cache
- **Abrash Ch.19 (c-style-guide §5.5.1–5.5.2)**: the BFS inner loop has
  unpredictable `visited[v]` branches. Consider branchless techniques or
  restructuring. Don't unroll loops — branch prediction handles the backwards
  branch perfectly.
- **Abrash Ch.15–17 (c-style-guide §5.4)**: working-set analysis. The flow
  network is ~85 MB for a 200K-vertex mesh — fits in Kaggle's L3 but is
  tight. If it spills to DRAM, BFS performance degrades catastrophically.
  Monitor cache miss rates with `perf stat`.
- **Bentley Ch.7 (c-style-guide §5.1)**: always do the arithmetic before
  writing the code. The feasibility estimates in §Performance Budget below
  are Bentley-style calculations.
- **Sanglard Ch.7 (c-style-guide §5.6.1)**: prefer a serial pipeline with
  parallelism inside each stage. Bridge cut is recursive — don't try to
  parallelize across recursion levels.

### Arena and Memory Management
- **Hanson Ch.6 (c-style-guide §1.1–1.2)**: arena allocation for all flow
  graph memory. Use `Arena_save`/`Arena_restore` marks to reclaim scratch at
  each recursion level. The save/restore pattern from Hanson Ch.16 (Text
  module) is the exact fit: snapshot before allocating the flow network,
  restore after extracting the min-cut result.
- **Sanglard Ch.2**: DOOM zone allocator precedent. The per-recursion-level
  save/restore maps directly to DOOM's PU_LEVEL purge tag pattern.

### Error Handling
- **Hanson Ch.4 (c-style-guide §2.6)**: Exceptions (`TRY`/`EXCEPT`) for
  timeout and OOM. Checked runtime errors (`assert`) for contract violations.
  Never silently produce wrong output — assert and fail the volume.
- **TPOP Ch.5**: detect errors at low level, handle at high level. The bridge
  cut module detects inconsistencies; the per-volume dispatch handles fallback.

---

## Performance Budget

**Total pipeline budget**: ~180s per volume on Kaggle (2 cores @ 2.2 GHz). This 180s is the C preprocessing share of the ~270s/volume total (the remaining ~90s is nnUNet inference).
**Step 2 budget**: 15–75s in Python. Target for C: **2–10s** per component
on Kaggle, depending on mesh size and recursion depth. For a volume with ≤ 5
multi-sheet components: **5–30s total.**

### Back-of-Envelope: Single Bridge Detection (N=100K verts, E=300K edges)

| Sub-step | Work | Est. Time (Kaggle) |
|---|---|---|
| PCA + projection | O(N) = 100K ops | ~0.1 ms |
| k-ring BFS (SEED_RING=20) | O(ball) ≈ 5K verts | ~0.05 ms |
| Build vertex-split CSR | O(N + E) ≈ 400K ops | ~5 ms |
| Edmonds-Karp (100 BFS passes) | 100 × O(V'+E') = 100 × 2.6M | ~2.6 s |
| Path trace + augment (100×) | 100 × ~20 hops | ~0.2 ms |
| Min-cut extraction | O(V') = 200K | ~2 ms |
| is_bridge heuristic | O(1) | negligible |
| Build KD-tree (boundary, ~1K pts) | O(B log B) = 10K | ~0.1 ms |
| Nearest-boundary query (100K pts) | 100K × O(√B) ≈ 100K × 30 | ~3 ms |
| Face partitioning | O(M) = 200K | ~2 ms |
| Vertex reindex + sub-mesh extract | O(N + M) | ~5 ms |
| **Total per detection+split** | | **~2.6 s** |

**Edmonds-Karp BFS dominates.** The 2.6s estimate assumes each BFS access
costs ~10ns (L3 cache hit for the ~85 MB working set). If the flow network
spills to DRAM (likely on Kaggle's smaller L3), cost rises to ~20ns/access
→ ~5.2s per detection. This is the critical uncertainty.

### Recursive Depth Analysis

Typical case: 2–3 recursion levels per component (multi-sheet components
rarely need > 3 cuts). Each level operates on a smaller sub-mesh (roughly
halved per split). Total work is geometric: level 0 on N verts, level 1 on
~N/2 verts, level 2 on ~N/4 verts. Total ≈ 2N — approximately 2× the cost
of a single detection.

Worst case: 10 recursion levels, each splitting roughly evenly. Total ≈ 10×
cost of level 0 (since sub-meshes shrink). For N=100K: ~26s. This is tight
but fits the 30s budget.

Pathological case: a component that produces highly asymmetric splits
(split_ratio ~0.20), so each level barely shrinks the larger piece. This
converges slowly. The depth limit (10) and the timeout prevent runaway.

### Scaling: 7950X → Kaggle

- **CPU speed**: 7950X single-core is ~2.6× faster than Kaggle @ 2.2 GHz.
  A 2.6s detection on Kaggle ≈ 1.0s on 7950X.
- **Cache**: 7950X L3 is 64 MB (shared) with 1 MB L2 per core. Kaggle's
  Xeon L3 may be smaller. The flow network (~85 MB for 100K verts) may
  spill on Kaggle but fit in 7950X L3. This is a significant variable —
  profile with `perf stat` on both.
- **Threading**: Bridge cut is recursive and operates on one component at a
  time. No parallelism within bridge detection (BFS is inherently serial).
  Parallelism across components is possible but limited by Kaggle's 2 cores.
  Process 2 components simultaneously with OpenMP tasks or pthreads, but
  ensure each thread gets its own arena mark.

### Memory Budget

For a single bridge detection on N=100K vertices:
```
CSR flow graph:
  offset[200K+3]:   ~800 KB
  target[2.6M]:     ~10.4 MB
  cap[2.6M]:        ~10.4 MB
  flow[2.6M]:       ~10.4 MB
  reverse[2.6M]:    ~10.4 MB
  Subtotal:         ~42.4 MB

BFS workspace:
  parent_edge[200K]: ~800 KB
  queue[200K]:       ~800 KB
  Subtotal:          ~1.6 MB

Min-cut:
  reachable[200K]:   ~200 KB

KD-tree on boundary:
  nodes[~1K]:        ~32 KB

Total scratch per recursion level: ~44 MB
```

With arena save/restore, memory from level N is reclaimed before level N+1
(DFS order). Peak memory = one level's scratch + the accumulated sub-meshes.
For 10 levels: ~44 MB scratch + ~12 MB sub-meshes = ~56 MB. Fits easily
in 13 GB.

**For larger meshes (N=200K)**: scratch doubles to ~85 MB. Still fits.

---

## Key C Rules

**§1.2 Arena Save/Restore** (c-style-guide): use `Arena_save` before
allocating the flow network at each recursion level. After extracting the
min-cut result (reachable array) into persistent storage, `Arena_restore`
to reclaim the flow graph, BFS arrays, and KD-tree. This is the exact
pattern from Hanson Ch.16 (Text save/restore) and Sanglard's DOOM zone
allocator. Each recursion level is a save/restore bracket:
```c
Arena_Mark mark = Arena_save(arena);
/* allocate flow network, BFS arrays, boundary KD-tree */
/* ... run Edmonds-Karp, extract min-cut ... */
/* copy sub-mesh vertex/face arrays to persistent region BEFORE restore */
Arena_restore(arena, mark);
```
[Hanson Ch.6, Ch.16; Sanglard Ch.2; c-style-guide §1.2.2]

**§5.1 Feasibility Check**: done above. Single detection ~2.6s on Kaggle
with L3-resident working set. Must monitor for DRAM spill with `perf stat`.
The MAX_FLOW_LIMIT = 100 is the essential safety valve that caps BFS passes.
Without it, worst-case is intractable.
[Bentley Ch.7; c-style-guide §7.4.2]

**§5.4 Working Set**: the flow network for N=100K is ~42 MB. On 7950X
(64 MB L3) it fits. On Kaggle's Xeon, it may not. If profiling shows DRAM
spills, consider: (a) using `int16_t` for capacity (all caps are 1 or INF=32767),
saving ~5 MB, (b) compressing the CSR by not storing zero-capacity reverse
edges initially, (c) shrinking the mesh before bridge detection (remove
clearly-interior vertices that can't be in the cut).
[Abrash Ch.19; c-style-guide §5.4.1–5.4.7]

**§5.5 Branch Prediction**: the BFS inner loop's `if (residual > 0 && parent_edge[v] == -1)`
is an unpredictable branch on early iterations (roughly 50/50 whether a
neighbor is visited). On later iterations, most neighbors are visited →
the branch becomes predictable. Don't try to optimize this away — the
BFS is memory-bound, not branch-bound. Profile first.
[Abrash Ch.19; c-style-guide §5.5.2]

**§7.1 CSR Construction**: build the vertex-split flow graph in two passes.
Pass 1: count the out-degree of each node (iterate all edges). Pass 2: prefix
sum for offsets, then fill `target[]`, `cap[]`, `reverse[]` using a write
cursor per node. Pre-allocate forward + reverse edges together.
[CLRS Ch.22.1; c-style-guide §7.1.3–7.1.5]

**§7.4.4 Reverse Edge Bug**: the #1 max-flow implementation bug. When
augmenting flow along a path, BOTH the forward edge and its reverse must be
updated: `flow[e] += bottleneck; flow[rev[e]] -= bottleneck`. Test this with
a 2-vertex, 1-edge graph where the answer is flow=1.
[CLRS §26.2; Skiena §8.5; c-style-guide §7.4.4]

**§9.3.2 Graceful Empty Cases**: handle every degenerate case:
- Zero-face mesh → return unchanged.
- All vertices in one seed set (overlap removes everything) → return NO BRIDGE.
- Flow = 0 (source disconnected from sink in residual graph) → degenerate,
  return NO BRIDGE.
- Either sub-mesh has < 10 faces → degenerate split, return original.
- Depth > 10 → return unchanged.
[Elements Ch.4; c-style-guide §9.3.2]

**§9.4.4 Index Overflow**: the vertex-split graph has nodes indexed 0 to
2N+1. For N=200K, max node index = 400,001. Products like
`node * max_degree_estimate` must use `size_t` for intermediate computation.
CSR offset array should be `int32_t` (max value ~2.6M for 100K verts, fits
int32 easily). But target array indices reference nodes in [0, 400K], also
fine for int32. Use `int32_t` for all graph arrays.
[c-style-guide §9.4.4; Pitfalls — semantic]

**§2.6 Error Handling**: wrap the per-component bridge cut in a TRY-EXCEPT
at the caller level. If bridge cut raises `Timeout`, return the original mesh
unchanged. If it raises `Arena_Failed` (OOM on flow graph allocation), return
unchanged. Never crash. Never produce partial results (half-split mesh).
[Hanson Ch.4; c-style-guide §2.6.1–2.6.6]

---

## Data Structures

### From `common/` (shared modules)

| Module | Type / Function | Used For |
|---|---|---|
| `arena.h` | `Arena_T`, `Arena_save`, `Arena_restore` | All allocation |
| `csr.h` | `CSR_T`, `CSR_new_from_coo` | Mesh adjacency (for k-ring BFS) |
| `kdtree.h` | `KDTree_T`, `KDTree_nearest` | Gap exclusion distance queries |
| `bfs.h` | `BFS_kring`, `BFS_find_path` | Seed expansion, Edmonds-Karp |
| `pca.h` | `PCA_normal`, `PCA_project` | PCA normal + vertex projection |
| `mesh_types.h` | `ComponentMesh` | Input/output mesh struct |
| `oracle.h` | `Oracle_count_sheets` | Look-ahead pruning in recursion |

### Created internally by Step 2

**Flow network** (vertex-split CSR graph):
```c
typedef struct FlowNet {
    int32_t  n_nodes;          /* 2*N + 2 */
    int32_t  n_edges;          /* total directed edges including reverses */
    int32_t *offset;           /* [n_nodes + 1] — CSR row offsets */
    int32_t *target;           /* [n_edges] — edge target node */
    int32_t *cap;              /* [n_edges] — edge capacity (1 or INF) */
    int32_t *flow;             /* [n_edges] — current flow (init 0) */
    int32_t *reverse;          /* [n_edges] — index of reverse edge */
    int32_t  super_source;     /* node index 2*N */
    int32_t  super_sink;       /* node index 2*N + 1 */
} FlowNet;
```
All arrays arena-allocated in a single save/restore bracket.

**BFS workspace** (reused across augmenting path iterations):
```c
int32_t *parent_edge;   /* [n_nodes] — edge used to reach each node, -1 = unvisited */
int32_t *queue;         /* [n_nodes] — flat BFS queue */
```
`parent_edge` is `memset` to -1 at the start of each BFS. `parent_edge[source] = -2`
(visited sentinel, no parent). Reuse the same arrays for all 100 BFS passes —
no reallocation.

**Reachable array** (min-cut result):
```c
uint8_t *reachable;     /* [N] — per original vertex, 0 or 1 */
```
Derived from the last failed BFS: `reachable[v] = (parent_edge[2*v] != -1)`.
This is the persistent output that must be copied before `Arena_restore`.

**Seed sets** (temporary):
```c
uint8_t *is_source_seed; /* [N] — bitmap */
uint8_t *is_sink_seed;   /* [N] — bitmap */
```
Populated by k-ring BFS, consumed during graph construction.

**Boundary vertex set** (for split_mesh):
```c
float   *boundary_pts;   /* [B * 3] — positions of boundary vertices */
int32_t *boundary_ids;   /* [B] — original vertex indices */
size_t   n_boundary;
```
Built by scanning edges where `reachable[u] != reachable[v]`.

---

## Interfaces

### Primary Entry Point

```c
/*
 * BridgeCut_process — Step 2 entry point.
 *
 * Given a multi-sheet component mesh, recursively detect bridges
 * and split until all pieces are single-sheet (or unsplittable).
 *
 * Parameters:
 *   arena       — arena for all allocations (uses save/restore internally)
 *   mesh        — input component mesh (oracle sheet_count >= 2)
 *   grid_size   — oracle UV grid size (typically 320)
 *   timeout_sec — max wall time for this component (e.g. 30.0)
 *   n_threads   — thread count (unused within bridge cut, reserved)
 *   out_meshes  — output: array of resulting sub-meshes
 *   out_count   — output: number of resulting sub-meshes
 *
 * Returns 0 on success. On failure (timeout, OOM, degenerate mesh),
 * returns -1 and sets *out_meshes = &mesh, *out_count = 1 (original
 * mesh returned unchanged).
 */
int BridgeCut_process(Arena_T               arena,
                      const ComponentMesh   *mesh,
                      int                    grid_size,
                      double                 timeout_sec,
                      int                    n_threads,
                      ComponentMesh        **out_meshes,
                      size_t                *out_count);
```

### Internal Functions (all `static`)

```c
/* Detect bridge: run max-flow, return flow value and reachable array */
static int detect_bridge(Arena_T arena, const ComponentMesh *mesh,
                         int32_t *out_flow, float *out_height,
                         float *out_split_ratio, uint8_t **out_reachable);

/* Evaluate is_bridge heuristic */
static bool is_bridge(int32_t flow, float height, float split_ratio);

/* Split mesh at the min-cut boundary with gap exclusion */
static int split_mesh(Arena_T arena, const ComponentMesh *mesh,
                      const uint8_t *reachable,
                      ComponentMesh *out_side_a, ComponentMesh *out_side_b);

/* Recursive decomposition */
static int process_recursive(Arena_T arena, const ComponentMesh *mesh,
                             int depth, int grid_size,
                             struct timespec *timeout_at,
                             ComponentMesh **results, size_t *n_results);
```

### Connection to Step 1 (upstream)

The pipeline driver runs the Oracle on each component from Step 0. Components
with sheet_count ≥ 2 are routed to `BridgeCut_process`. Single-sheet
components skip to Step 4.

### Connection to Step 3 (downstream)

After bridge cut, each resulting sub-mesh is re-checked by the Oracle. If
any sub-mesh still has sheet_count ≥ 2, it enters the raycast separator
(Step 3). Bridge cut handles thin bridges; Step 3 handles broad overlapping
sheets that don't have a thin neck.

### Connection to Step 4 (bookkeeping)

Step 4 collects all sub-meshes (from original single-sheet components,
bridge-cut pieces, and raycast-separated pieces) and assigns new sequential
component IDs.

---

## External Dependencies

**None for the core algorithm.** Bridge cut is pure C: CSR graph construction,
BFS, flow augmentation, KD-tree (or brute-force) distance queries. No external
libraries.

The only indirect dependency is the Oracle (Step 1), which is called within
the recursive loop for look-ahead pruning. The Oracle itself has no external
dependencies.

---

## Test Strategy

### Comparison with Python Reference

Bridge cut's output is deterministic for a given mesh (Edmonds-Karp with BFS
is deterministic, PCA is deterministic after sign correction, k-ring BFS is
deterministic with consistent neighbor ordering).

**Primary comparison**: for each multi-sheet component on each test volume,
compare the number and sizes of resulting sub-meshes between Python and C.
An exact match of sub-mesh vertex counts is the target. If counts differ,
investigate:

1. **Flow value mismatch**: compare `total_flow` values. They must be
   identical — Edmonds-Karp is deterministic given the same graph.
2. **Reachable set mismatch**: compare the reachable array vertex by vertex.
   If the flow values match but reachable sets differ, there may be a
   tie-breaking difference in BFS neighbor ordering. Ensure CSR neighbor
   ordering matches Python's adjacency list ordering.
3. **is_bridge heuristic**: compare `flow`, `height`, `split_ratio`,
   `flow_height_ratio` values. Differences here indicate PCA or projection
   bugs.
4. **Sub-mesh face count**: after split_mesh, compare per-side face counts.
   If they differ, the exclusion zone or face partitioning logic diverges.

### Unit Tests

- **Minimal bridge**: two triangles connected by a single edge. Flow = 1.
  is_bridge should return true (if height conditions are met). Split should
  produce two single-triangle meshes.
- **No bridge**: a flat planar mesh (single sheet, well-connected). Flow
  should hit MAX_FLOW_LIMIT. is_bridge returns false. Output = input.
- **Two parallel planes connected by a thin tube**: the canonical bridge
  test case. The tube vertices should be the min-cut. Verify the exclusion
  zone removes the tube, producing two clean planes.
- **Asymmetric split (split_ratio < 0.20)**: a mesh with a tiny dangling
  fragment. is_bridge returns false despite low flow, because the split is
  too asymmetric. Output = input unchanged.
- **Gray zone**: a mesh with flow_height_ratio between 0.36 and 1.3, with
  split_ratio ≥ 0.40 and height ≥ 20. Verify the gray-zone rule correctly
  identifies it as BRIDGE.
- **Degenerate: zero faces**: returns empty, doesn't crash.
- **Degenerate: < 10 faces on one side after split**: returns original unchanged.
- **Recursion depth limit**: a mesh that triggers splits at every level
  (synthetically constructed). Verify recursion stops at depth 10.
- **Timeout**: set a 0.1s timeout on a large mesh. Verify the function returns
  the original mesh within ~0.2s (timeout + cleanup overhead).

### Flow Correctness Checks (Debug Mode)

In debug builds (`#ifndef NDEBUG`), run these after every Edmonds-Karp
invocation:

1. **Flow conservation**: for every node except S and T, verify
   Σ(flow on incoming edges) = Σ(flow on outgoing edges).
   [CLRS §26.1]
2. **Capacity constraint**: for all edges, verify `0 ≤ flow[e] ≤ cap[e]`.
   [CLRS §26.1]
3. **Max-flow = min-cut**: after extracting the reachable set, compute the
   capacity of the cut (sum of `cap[e]` for edges crossing from S to T)
   and verify it equals `total_flow`.
   [CLRS §26.2 Theorem 26.6]
4. **Reverse edge consistency**: for all edges e, verify
   `flow[e] + flow[reverse[e]] == 0`.

These checks are O(V+E) each — negligible compared to the BFS passes.
Leave them in until the flow implementation is proven correct on ≥ 20
test volumes.

### Regression Suite

- `make test_step2` runs all bridge cut tests in < 60s on the 7950X.
- One test per bug, growing monotonically.
- Include a stress test with synthetic meshes of N=200K vertices (matching
  the 200K upper bound from the performance estimates).
- Include a `taskset -c 0,1` run on the 7950X to simulate Kaggle.

---

## Pitfalls

### 1. The Reverse-Edge Bug (CLRS §26.2; Skiena §8.5; c-style-guide §7.4.4)

The #1 most common max-flow implementation bug: forgetting to update the
reverse edge during augmentation. Must update BOTH:
```c
flow[e] += bottleneck;
flow[reverse_edge[e]] -= bottleneck;
```
If only the forward edge is updated, the residual graph is inconsistent,
BFS finds invalid paths, and the flow value is wrong (typically too high).
The symptom: flow exceeds the true min-cut capacity. The debug-mode
check (max-flow = min-cut capacity) will catch this.

### 2. Flow Initialization (CLRS §26.2 rule 30)

Initialize `flow[]` to 0, not to `cap[]`. Residual capacity of a fresh
edge is `cap[e] - flow[e]` = `cap[e] - 0` = `cap[e]`. If you initialize
`flow[]` to `cap[]`, all edges start saturated and BFS finds no path.
Symptom: flow = 0 on every component.

### 3. Node Indexing: v_in = 2v, v_out = 2v+1

Getting the vertex-split indexing wrong silently corrupts the graph. The
convention must be consistent across all four contexts:
- Graph construction (internal edges, cross edges, super-source/sink edges)
- BFS traversal
- Min-cut extraction (mapping node reachability back to vertex reachability)
- Capacity assignment (internal edges get cap 1, others get INF)

A single misplaced `2*v` vs `2*v+1` in any of these contexts produces a
subtly wrong flow value. Write helper macros:
```c
#define V_IN(v)  ((v) * 2)
#define V_OUT(v) ((v) * 2 + 1)
```
Use them everywhere. Never compute the indices inline.
[c-style-guide §3.3 — hide confusing expressions in named macros]

### 4. Seed Overlap (PIPELINE_REFERENCE §6 step 6)

After expanding source and sink seeds by SEED_RING=20 hops, the two sets
may overlap (if the mesh is small or highly connected). Overlap vertices
must be removed from BOTH sets. If they're not, the super-source and
super-sink are connected through the overlapping vertices with INF capacity
— the flow is immediately INF (or MAX_FLOW_LIMIT), and every component
looks like "no bridge." Symptom: bridge cut never detects any bridges.

### 5. PCA Sign Convention (same pitfall as Steps 0 and 1)

The PCA normal has arbitrary sign. The "largest absolute component positive"
convention must be applied consistently. If Step 2's PCA sign differs from
Step 1's, the projection direction flips, source/sink swap, and the
is_bridge heuristic may evaluate differently (since `height` is the same
but which side is "source" and which is "sink" swaps). For the flow value
this doesn't matter (the graph is undirected), but for the split_ratio
(which depends on which side is "reachable"), it can change results.
Mitigate by reusing Step 0's PCA normal when possible.

### 6. The ≤ vs < in is_bridge (Elements Ch.4 rule 25)

`flow_height_ratio ≤ 0.36` means BRIDGE. If you write `< 0.36`, exact
boundary cases (ratio exactly 0.36) are classified as NO BRIDGE. Write a
unit test for the boundary: flow=36, height=100 → ratio=0.36 → should be
BRIDGE. Same for `split_ratio < 0.20` (strict — exactly 0.20 IS a bridge).
The gray-zone rule uses `≤ 1.3` AND `≥ 0.40` — verify all four boundary
conditions independently.

### 7. INF Capacity Value

Define `INF` as a specific integer (not `INT_MAX`, which risks overflow
during `bottleneck = min(cap[e] - flow[e], ...)` if flow is negative).
Use a large but not maximal value: `#define FLOW_INF (1 << 30)` (≈ 1 billion).
This is safely larger than MAX_FLOW_LIMIT (100) and won't overflow in
any addition or subtraction. If you use `INT_MAX` and `flow[e]` is -1
(from a reverse-edge update), `cap[e] - flow[e]` = `INT_MAX + 1` = overflow.
[Pitfalls — semantic; c-style-guide §9.5.1]

### 8. Sub-Mesh Extraction: Vertex Reindexing

When extracting sub-meshes after splitting, vertices must be reindexed
to [0, n_sub-1]. The reindex map must be built before remapping faces.
Common bug: reindex map is off-by-one because excluded vertices create
gaps. Use a two-pass approach: pass 1 assigns new indices (skip excluded
and wrong-side vertices); pass 2 remaps face indices through the map.
If any face index maps to -1 (vertex was excluded), the face should have
been removed in step 4.
[c-style-guide §9.3.3 — asymmetric bounds]

### 9. Arena Restore Before Sub-Mesh Copy

The sub-mesh vertex and face arrays are allocated in the current arena
save bracket. If you `Arena_restore` BEFORE copying them to persistent
storage, the data is reclaimed and you read freed memory. The discipline:
(a) allocate sub-mesh arrays, (b) populate them, (c) copy/transfer them
to the output array (which must be allocated in the OUTER save bracket or
in persistent arena space), (d) THEN restore.

Alternatively, allocate the output sub-mesh arrays from a *separate* arena
mark (saved before the recursive call), so they survive the inner restore.
[Hanson Ch.6; Hanson Ch.16 rule 64]

### 10. Timeout Race Condition

The timeout check (`clock_gettime(CLOCK_MONOTONIC)`) is at the top of
each recursive call. If a single Edmonds-Karp invocation takes 30s (large
mesh, many BFS passes), the timeout fires only after that invocation
completes. The worst-case overshoot is one full Edmonds-Karp run. Mitigate
by also checking the timeout inside the BFS loop (every 10th iteration),
accepting the small overhead. Use a global volatile flag set by the timeout
check, not a per-iteration clock call:
```c
if (iter % 10 == 0 && *timeout_flag) break;
```
[c-style-guide §2.6.2 — timeout exception]

### 11. CSR Neighbor Ordering Affects BFS Tie-Breaking

Edmonds-Karp's BFS finds the shortest augmenting path, but when multiple
shortest paths exist, the tie is broken by the order in which neighbors
are visited — which depends on CSR construction order. If the C CSR builds
neighbors in a different order than Python's adjacency list, different
augmenting paths are chosen, potentially leading to different min-cuts (of
the same value, but with different reachable sets). This means sub-mesh
sizes may differ even though the flow value matches.

This is not a bug — different min-cuts of the same value are all correct.
But it complicates Python-vs-C comparison. For debugging, sort CSR neighbors
by index during construction to make the ordering deterministic and
comparable.
[CLRS §26.2 — multiple optimal solutions]
