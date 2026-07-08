#ifndef PIPELINE_CONSTANTS_INCLUDED
#define PIPELINE_CONSTANTS_INCLUDED

/* Pre-Step 0 — Voxel denoise (tunnel removal) */
#define TUNNEL_GAP_FILL_THRESH   4   /* min 6-conn FG neighbors to fill BG voxel */
#define TUNNEL_GAP_FILL_ITERS   10   /* max gap-fill passes */
#define TUNNEL_BG_CC_MAX        27   /* max BG CC size to fill (26-conn) */
#define MORPH_CLOSE_ITERS        3   /* multi-pass morph close iterations */
#define SIMPLE_FILL_MIN_FG       4   /* min FG 6-conn neighbors to be candidate */
#define SIMPLE_FILL_MAX_ITERS   20   /* max iterations for simple-point fill */
#define BRIDGE_FILL_REACH       10   /* max distance to check for opposing FG */
#define GENUS0_BBOX_PAD          2   /* voxel padding around component bbox */
#define GENUS0_MAX_ITERS       200   /* max thinning passes for genus-0 */
#define GENUS0_MAX_FILLABLE  15000000  /* skip genus0 if >15M fillable voxels */

/* Pre-Step 0 — Garbage prediction rejection (input validation).
 * Reject cubes whose prediction is a big SOLID slab/rectangle (nnUNet failure)
 * rather than thin recto-surface sheets. Two AND-combined signals; tuned
 * conservatively (only reject unambiguous slabs). See src/extract/pred_reject.c.
 * Calibrated 2026-06-03 on the PHerc0139-4x21x21 grid (4 labelled garbage cubes
 * vs the real dense tangle z04736_y04608_x03968). */
#define GARBAGE_ERODE_R          2     /* 3D 6-conn erosion passes (thickness probe) */
#define GARBAGE_ERODE_MAXPASS   16     /* cap passes when measuring max_thickness */
#define GARBAGE_INTERIOR_FRAC  0.50f   /* reject needs >= this frac of FG surviving erosion.
                                        * Calibration found a clean EMPTY gap on the grid:
                                        * real cubes <=0.32, garbage slabs >=0.68. 0.50 sits
                                        * dead-centre for maximal margin both ways. */
#define GARBAGE_INTERIOR_MIN  2000     /* ...and >= this many surviving voxels (abs floor) */
#define GARBAGE_RECT_AREA_FRAC 0.10f   /* slice's largest comp >= this frac of slice area */
#define GARBAGE_RECT_FILL      0.75f   /* ...AND fills >= this frac of its bbox (solid rect) */
#define GARBAGE_RECT_FRAC      0.20f   /* reject needs >= this frac of slices be rect frames */
#define GARBAGE_RECT_RUN        12     /* ...AND a run of >= this many consecutive rect frames */

/* Step 0 */
#define MIN_CC_SIZE          500    /* voxels - discard smaller components */
#define MAX_COMPONENTS        20    /* keep top N by size */
#define CLEANUP_MICRO_HOLE_MAX 6   /* max boundary loop verts for fill */
#define MIN_FRAGMENT_FACES   100    /* discard sub-components smaller */
#define RESPLIT_MIN_COMP_VERTS 200  /* mesh_resplit: re-mesh a connectivity-component
                                     * only if it has >= this many verts (smaller
                                     * disconnected fragments are dropped) */
#define KIBBLE_AREA_FRAC     0.02f  /* component_cull: drop connectivity-components
                                     * whose surface area is < this fraction of the
                                     * total meshed cube area (post hole-fill/QEM) */

/* Step 0 — MLS-midpoint projection (LOP family, μ=0). Collapses the
 * MC double-envelope to its single-sided centerline. Halo-deterministic:
 * neighbour set at each vertex depends only on the halo-fixed point cloud
 * within MLS_PROJECT_RADIUS_VOX. Constants must satisfy
 *   MLS_PROJECT_RADIUS_VOX + 1 <= halo_voxels
 * so a vertex inside the cube never queries a neighbour outside the halo.
 * Satisfied since 2026-05-29: the halo default was raised to 13 (>= 12+1) so
 * near-seam boundary verts get full two-sided support and adjacent cubes
 * converge to the same height (no seam walls).
 * See plan: the-main-problem-with-composed-badger.md */
#define MLS_PROJECT_RADIUS_VOX  12.0f /* Wendland kernel radius (voxels).
                                       * Must be ~2x layer thickness so the
                                       * Wendland tail still has meaningful
                                       * weight on the far face. For ~5 vox
                                       * layers, R=12 gives weight ~0.05 on
                                       * the far envelope face — enough to
                                       * actually collapse it. */
/* Iterative LOP with tangent-plane projection (V' = V − ((V−C)·N) N).
 * Each iter only moves vertices perpendicular to the local sheet normal,
 * so in-plane geometry is preserved across iterations regardless of
 * how the cloud's center of mass drifts. Multi-iter LOP converges the
 * point cloud to near-zero through-thickness (~0.01 vox at ITERS=100)
 * without the catastrophic in-plane shrinkage that centroid replacement
 * suffers at high iter counts. Step0 ping-pongs between mls_verts and
 * a scratch buffer so the read/write aliasing doesn't corrupt the
 * per-vertex result; normals are computed every iter (not just the
 * last) because tangent projection needs them. */
#define MLS_PROJECT_ITERS       5    /* 2026-05-29: cut 30->5. Overlap-split +
                                      * per-sheet re-LOP now carries flatness;
                                      * a single global pass no longer has to
                                      * over-flatten, and 5 tangent-plane iters
                                      * suffice for BPA -- ~6x faster Step 0. */
#define MLS_RESPLIT_ITERS      20    /* 2026-06-03: the re-LOP-from-original after a
                                      * split (mesh_resplit) carries the FINAL per-sheet
                                      * flatness, so it gets more passes than the
                                      * speed-sensitive extract LOP. Raised 5->10->20 to
                                      * smooth residual through-thickness "lumps" left on
                                      * post-split components (sheets 1/2 of
                                      * z04736_y04224_x01920 still bumpy at 10). */

/* Seam-band pin (cross-cube convergence). Extract MLS-projects owned verts within
 * MLS_PROJECT_RADIUS_VOX of the cube boundary with HALO (two-sided) support, so
 * adjacent cubes agree on the seam geometry. The post-split re-LOP passes
 * (resurface_own, remesh_pieces) re-MLS the seam WITHOUT the halo and drift it
 * independently per cube, breaking that agreement (the welded "offset" / high-
 * lambda good-join). Pin verts within this band of the owned-box boundary at their
 * extract positions through every re-LOP. Env VES_SEAM_PIN_OFF disables;
 * SEAM_PIN_BAND_VOX overrides the width. */
#define SEAM_PIN_BAND_VOX  MLS_PROJECT_RADIUS_VOX

/* Winding-gate default tolerance, in TURNS about the umbilicus. The seam-bridge
 * winding gate rejects a bridge whose winding phase w = r/pitch - theta/2pi
 * differs by more than this between front edge and candidate. At seam-bridge
 * scale (1-3 vox chords) this is effectively a radial gate of tol*pitch vox.
 * 2026-07-08 A/B on the core-containing 4x5x5 (true pitch ~9.5 vox): with the
 * seam pin ON, ungated = 27 cross-wrap handles; effective 2.0-2.5 vox (this
 * tol at the swept pitches) = 4-5 handles, +0.3-1.3k unpaired; looser 3.5 vox
 * = 9. On non-core dumps 2.0 vox zeroes all handles. Residual core handles are
 * sub-gate wrap contacts -> need a fusion-line cutter, not a tighter radius
 * (tightening to 1.5 vox only fragments). Env SEAM_WIND_TOL overrides. */
#define SEAM_WIND_TOL_DEFAULT_TURNS 0.25
#define MLS_RESPLIT_ASSIGN_MARGIN_VOX 1.5f
                                     /* 2026-06-03: re-LOP point->piece assignment
                                      * margin. A parent original point is assigned
                                      * to its NEAREST split piece only when the
                                      * next-nearest piece is at least this much
                                      * farther; otherwise the point straddles the
                                      * split seam (or belongs to an adjacent
                                      * close-wrap ~2-3 vox away) and is DROPPED,
                                      * not vacuumed in. Stops the re-LOP+re-BPA
                                      * from grabbing the other wrap and folding the
                                      * sheet (the step7_cc_bpa_003 fold). 0 = old
                                      * greedy nearest-vertex (no margin). */
#define MLS_WELD_EPS_VOX        0.25f /* merge verts within this distance
                                       * after MLS projection. raw-snap will
                                       * later pull the merged vert onto the
                                       * NN prediction, so coarse weld is OK. */
#define MLS_MIN_NEIGHBOURS       4    /* below this, leave vert in place
                                       * (rare — only at extreme corners). */

/* Step 0 — Ball-Pivoting Algorithm (Bernardini 1999). Triangulates the
 * LOP-collapsed oriented point cloud into a manifold mesh. Replaces
 * the marching-cubes + backface-cull pair: voxel-center sampling gives
 * a single-layer point cloud (no double envelope to undo), LOP smooths
 * it, BPA stitches it. Pivot radius needs to be > 1 vox (the voxel-
 * center grid spacing) so the ball sees enough neighbours to pivot
 * smoothly, but small enough that the empty-ball test still discrim-
 * inates between the surface and itself across thin sheets.
 *
 * rho is a COVERAGE/ANTI-MERGE knob, NOT a quality knob. A sweep on the
 * baseline cube (z04480_y03584_x02816, rho ∈ {0.85, 1.0, 1.5}) left mean-min-
 * angle (30.4 deg) and mean-edge (0.771 vox) BIT-IDENTICAL: triangle shape is
 * fixed by the point ARRANGEMENT (a cubic-lattice surface shell projected onto
 * tilted planes — anisotropic), which the ball cannot move. What rho does
 * change is the band edges: the voxel-derived spacing is ~0.77 vox, the safe
 * band is rho ≈ 1.3–2x that. Below ~1.1x (rho 0.85) sparse sheets fall through
 * — comp007 went genus 0 -> -18, 86 boundary loops. Above ~2x the empty ball
 * crowds (more cocircular pinhole intruders) for no quality gain. 1.2 vox
 * (~1.55x spacing) sits central in the band: clear of the fall-through floor,
 * with more anti-merge and seam-bridge headroom than the old 1.5 (~2x, top of
 * band). Per-cube quality is unchanged; downstream PinholeFill/HoleFill close
 * the boundary loops either way. */
#define BPA_RHO_VOX             1.2f

/* Step 0 — BPA wall-guard (anti-fusion). Ball-Pivoting at rho=1.2 with no
 * inter-wrap clearance can bridge two papyrus wraps that sit only ~2-3 vox apart,
 * fusing a stack of near-parallel sheets into ONE component (the downstream
 * stacked-wrap "monster"). A fusion bridge is a "wall" triangle: it stands nearly
 * PERPENDICULAR to the local surface (MLS) normal -- |dot(face_n, mean_vertex_n)|
 * ~ 0 -- and reaches across the gap with a LONG edge. A true surface triangle has
 * its face normal ~PARALLEL to the vertex normals (|dot| ~ 1) and short edges.
 * Reject a candidate that is BOTH steep (|dot| < COS) AND long (spanning edge >
 * MIN_EDGE); the long-edge conjunct spares sharp folds/rims (short edges) so this
 * fires only on through-thickness bridges. Distinct from the face-coherence guard
 * (which compares to a neighbour FACE and is fooled by a tilted neighbour on a
 * curved stack). Env: BPA_WALL_GUARD_COS / BPA_WALL_MIN_EDGE_VOX sweep,
 * BPA_NO_WALL_GUARD disables; disabled in the seam bridge (it legitimately spans
 * gaps). The downstream depth-peel splitter separates any stack this misses. */
#define BPA_WALL_GUARD_COS      0.34f   /* reject faces tilted > ~70 deg from the surface normal */
#define BPA_WALL_MIN_EDGE_VOX   2.0f    /* ...but only if the spanning edge is at least this long (vox) */

/* Step 0 — pre-BPA owned-region cloud trim (halo mode). The READ halo + LOP give
 * near-boundary owned verts full two-sided support (denoised, cross-cube-stable
 * positions), but we then drop POINTS within INSET voxels of every cube face
 * before triangulating, so BPA's open boundary lands strictly INSIDE the face.
 * Insetting (not the old PAD=0 split exactly at the shared plane) is what keeps
 * adjacent cubes' charts from TOUCHING: a wrap GRAZING the seam plane otherwise
 * lands in BOTH cubes' owned regions (interleaved samples a fraction of a voxel
 * apart) and doubles — a z-fighting seam. With the inset each cube ends INSET vox
 * short of the face, leaving a 2*INSET gap the cross-cube seam bridge spans. Keep
 * 2*INSET + the ~1-vox natural row gap below 2*bridge-rho_max (~4.8 vox) so the
 * bridge can still reach across. Applies to all 6 faces (the outer-grid faces
 * lose a negligible INSET-vox band of real surface). */
#define BPA_OWNED_TRIM_INSET    1.0f

/* Seam weld -- pre-bridge sliver cull. A boundary triangle whose open edge lies
 * in a seam plane is a bridge-front candidate; if it is a near-degenerate sliver
 * (min altitude = 2*area/longest-edge, in voxels, below this) the bridge primes
 * off a needle and emits holes / flipped-normal faces. Drop such triangles before
 * bridging (their other edges become fresh boundary); a post-bridge remesh can
 * refine later. Env override: SEAM_SLIVER_MIN_ALT. 0 disables. */
#define SEAM_SLIVER_MIN_ALT     0.3f

/* Cross-cube seam bridge (no eat-back). The bridge ball radius adapts to the
 * post-QEM boundary spacing: base = BRIDGE_RHO_K * median near-seam boundary-
 * edge length, clamped to [MIN, MAX], with the caller's SEAM_RHO as a floor.
 * The front escalates rho -> 1.5rho -> 2rho capped at MAX so sparse boundaries
 * still close, while 2*MAX stays below the inter-wrap clearance (NO inter-wrap
 * mergers). A bridge triangle whose longest edge exceeds 2*MAX is rejected.
 * Override at runtime with env SEAM_RHO (floor) / SEAM_RHO_MAX (cap). */
#define BRIDGE_RHO_MIN          1.5f
#define BRIDGE_RHO_MAX          3.0f
#define BRIDGE_RHO_K            0.75f

/* Cross-cube seam weld bridges the cubes' clean trimmed boundaries directly (no
 * peelback): the directed-glue BPA front rolls across the ~1-vox seam gap at an
 * escalating radius capped so 2*BRIDGE_RHO_MAX (=6 vox) < ~CUT_GAP_DEPTH (7 vox),
 * which guarantees NO inter-wrap mergers. An earlier ring-peel was removed --
 * peeling receded the boundary and starved the bridge, leaving MORE open seam on
 * real dense data. See src/remesh/seam_weld.c. */

/* Step 0 — Spatial-pairing cull. Replaces the directional dot-test cull
 * (backface_cull_per_vert) which intrinsically tears curved sheets at any
 * fixed sign-convention boundary. After LOP collapses MC's double envelope,
 * each front face has a co-located back face within ~residual-thickness.
 * Spatial pairing finds these pairs by centroid hash + opposite-normal
 * test and drops one of each pair by a deterministic position tie-break.
 * Halo-deterministic by construction; no global sign reference. */
#define PAIRING_RADIUS_VOX     3.0f /* max centroid distance to call a pair.
                                     * Must exceed post-LOP residual through-
                                     * thickness (~1-2 vox for ~5-vox layers
                                     * at R=12), but small enough that genuine
                                     * single-sided geometry doesn't pair.
                                     * Connected components THICKER than this
                                     * stay as closed double envelopes — that
                                     * is the correct behaviour for things
                                     * that genuinely aren't thin sheets. */
#define PAIRING_DOT_THRESHOLD -0.5f /* normals dot < this counts as opposite-
                                     * winding pair. -0.5 ≈ 120° tolerance,
                                     * loose enough for noisy MC normals on
                                     * curved sheets. */

/* Step 0 — exported halo width. The READ halo (--halo CLI arg) governs
 * how much neighbour data is loaded for MC + LOP; the EXPORT halo is the
 * narrower band that's actually written to per-cube OBJ. We export only
 * what grid_weld needs to find duplicate verts (1 vox is plenty given
 * grid_weld's 1e-4 vox hash precision and MC's half-vox positioning).
 * The wider read halo is still in effect for MC and LOP correctness. */
#define EXPORT_HALO_VOX          1

/* Step 1 */
#define ORACLE_UV_GRID_SIZE  320    /* UV grid resolution for raycasting */
#define ORACLE_MIN_GAP        10    /* voxel gap to count as sheet boundary */

/* Flatten -- scroll UV unwrap (depth + winding). See src/flatten/. */
#define FLATTEN_AXIAL_BINS        256   /* centerline slices along the scroll axis */
#define FLATTEN_CENTERLINE_SMOOTH   5   /* moving-average half-window over slice centers */
#define FLATTEN_PX_PER_VOX        1.0f  /* tifxyz cells per voxel, both axes. UVs are
                                         * length-like (vox), matching vc_obj2tifxyz's
                                         * metric mode (~1 cell per UV unit). */
#define FLATTEN_MAX_GRID_DIM     8192   /* cap on either tifxyz grid dimension */
#define FLATTEN_SPIRAL_MIN_PTS     64   /* min component verts to trust a spiral fit */

/* Step 1c -- Developability cut (Stein/Grinspun/Crane 2018, sec 4.4). The FIRST
 * splitter: compute per-vertex Crane energy lambda_v = lambda_min(A_v) on the raw
 * BPA mesh (no optimization), mark seam vertices (lambda > EPS), remove the seam
 * band (+ geometric dilation) and keep the connected pieces -- a cut THROUGH the
 * non-developable ridge, i.e. exactly where the papyrus stops being one
 * developable sheet. A clean wrap has no spanning seam, so nothing disconnects
 * (self-gates: NO intra-sheet split). The piece set is accepted only if it is
 * measurably MORE developable than the parent (lambda-energy gate); else the
 * component falls through to bridge-cut/overlap unchanged. All env-overridable
 * (read in pipeline_cube.c) so they sweep without a rebuild; calibrate EPS on a
 * known merged-wrap vs clean cube. See src/split/dev_cut.c. */
#define DEV_CUT_EPS            0.05  /* lambda_min above which an interior vertex is
                                      * a seam (rad-weighted, scale/tessellation-
                                      * invariant). PROVISIONAL -- calibrate in the
                                      * bimodal gap between body and seam. Env:
                                      * DEV_CUT_EPS. */
#define DEV_CUT_GAP_DEPTH     3.5f   /* exclusion zone around seam verts (voxels),
                                      * half of CUT_GAP_DEPTH so the dev-cut barrier
                                      * is narrower than the bridge cut. Env:
                                      * DEV_CUT_GAP_DEPTH. */
#define DEV_CUT_MIN_COMP_VERTS 200   /* a survivor piece needs >= this many verts to
                                      * count as a real sheet (== RESPLIT_MIN_COMP_
                                      * VERTS). Smaller crumbs are ignored. */
#define DEV_CUT_GATE_MARGIN   0.02   /* accept the split only if every piece's
                                      * non-developable interior fraction is at
                                      * least this much BELOW the parent's. Env:
                                      * DEV_CUT_GATE_MARGIN. */
#define DEV_CUT_GATE_ABS_CEIL 0.05   /* ...AND below this absolute fraction, so a
                                      * still-tangled "piece" can't be accepted just
                                      * because it is marginally better than a very
                                      * bad parent. Env: DEV_CUT_GATE_ABS_CEIL. */

/* Step 1b-peel — stacked-wrap separator (src/split/depth_peel.c). BPA (rho=1.2,
 * no inter-wrap clearance) can fuse papyrus wraps ~2-3 vox apart into ONE
 * component; depth_peel projects verts onto the component PCA normal and cuts mesh
 * edges whose endpoints JUMP in depth by more than MIN_GAP -- the inter-wrap bridge
 * edges -- then returns the connected layers. A single sheet (however folded) has
 * only short edges, so no large jump: it self-gates to one piece (the "NO
 * intra-sheet split" guarantee). Runs ahead of dev-cut/bridge/overlap; whatever it
 * can't cleanly separate (a strongly curved stack) falls through to overlap-sep.
 * Env: DEPTH_PEEL_MIN_GAP_VOX / DEPTH_PEEL_GAP_DEPTH; VES_DEPTHPEEL_OFF disables. */
#define DEPTH_PEEL_MIN_GAP_VOX  1.5  /* depth jump (vox) across a triangle edge that
                                      * marks an inter-wrap bridge: above the intra-
                                      * sheet edge length (~0.6-1.0 vox), below the
                                      * 2-3 vox inter-wrap spacing. THE key knob;
                                      * calibrate on a stacked vs a clean cube. */
#define DEPTH_PEEL_GAP_DEPTH    1.0f /* exclusion zone (vox) dilated around the cut
                                      * (KD-tree, like CUT_GAP_DEPTH) so a slightly
                                      * ragged bridge band still severs cleanly.
                                      * Env: DEPTH_PEEL_GAP_DEPTH. */
#define DEPTH_PEEL_MIN_COMP_VERTS 200 /* a peeled layer needs >= this many verts to be
                                      * kept: a fold flap is small and dropped, a real
                                      * wrap is large (== RESPLIT_MIN_COMP_VERTS). */

/* Step 2 */
#define SEED_RING             20    /* BFS expansion hops for source/sink */
#define MAX_FLOW_LIMIT       100    /* Edmonds-Karp early exit */
#define CUT_GAP_DEPTH       7.0f   /* exclusion zone around cut (voxels) */
#define BRIDGE_MAX_DEPTH      10    /* max recursion depth */
#define FLOW_INF       (1 << 30)   /* ~1 billion - never INT_MAX */

/* Short-handle sever (post hole-fill). BPA can roll a thin self-bridge that
 * fuses a sheet into a genus>0 tangle (seen on ~4/100 4x5x5 cubes, up to
 * genus 10, every loop < 35 vox). Open every non-separating loop shorter than
 * this (the thin bridges) and keep longer ones -- a per-cube sheet patch has no
 * legitimate scroll-scale handle, but the length bound is a safety margin. Set
 * env VES_SEVER_OFF=1 to disable (for A/B). */
#define SEVER_MAX_LOOP_VOX   60.0   /* double: max handle-loop length to sever */

/* Step 3 — Overlap separator */
#define MIN_FRAGMENT_VERTS_S3   30     /* percent threshold for small-comp merge */
#define OVERLAP_GRID_SCALE       2.0    /* cell_size = scale * sqrt(median_tri_area) */
#define OVERLAP_NEG_WEIGHT       1e6    /* negative weight multiplier for overlaps */
#define OVERLAP_MIN_FACES        50     /* minimum faces to attempt separation */
#define OVERLAP_EDGE_EPS         1e-6   /* epsilon for 2D intersection tests (double) */

/* Step 5 — Cotangent Laplacian constants */
#define COT_WEIGHT_MIN          1e-8f  /* minimum cotangent weight (clamp) */
#define COT_SIN_GUARD           1e-12f /* guard for degenerate triangles */

/* Snap-back CG solver */
#define SNAP_CG_TOL          1e-6   /* relative residual tolerance */

/* Step 6 */
#define MIN_BARY_SUBDIV        6   /* minimum barycentric samples per edge */

/* Step 0 — QEM mesh simplification */
#define QEM_TARGET_RATIO          0.075f  /* keep 7.5% of faces. History: 0.10
                                            * -> 0.15 (2026-05, interior was
                                            * oversimplified) -> 0.075 (2026-05-30,
                                            * ~50% more aggressive; the post-QEM
                                            * edge-flip pass cleans up slivers the
                                            * heavier decimation leaves). */
#define QEM_BOUNDARY_WEIGHT       10.0    /* penalty for boundary edge planes */
#define QEM_MIN_FACES_FOR_SIMPLIFY 400    /* skip simplification below this */
#define QEM_DET_THRESHOLD         1e-12   /* 3x3 solve singularity guard */
#define QEM_NORMAL_DOT_THRESHOLD  0.0f    /* reject collapse if dot(n_v0, n_v1) < this */
#define QEM_FLIP_COS_THRESHOLD    0.0f    /* fold guard angle as a cosine. 0.0 = the original
                                           * >90deg reversal-only behaviour (reverted from the
                                           * 0.5/60deg experiment, which was the wrong lever --
                                           * the real fault is upstream in the multicut). */
#define QEM_DISPLACEMENT_CLAMP    1.0f    /* max opt_pos distance = clamp * edge_len (was 2.0) */
#define QEM_SAFE_RADIUS_FACTOR    0.5f    /* proximity radius = factor * median_edge_len */
#define QEM_KDTREE_REBUILD_INTERVAL 500   /* rebuild spatial index every N collapses */
#define QEM_PROB_SIGMA_FACTOR     0.05    /* σ_n = factor * median_edge_len */
#define QEM_MAINTENANCE_INTERVAL  10000   /* collapses between maintenance passes */
#define QEM_MAINT_SMOOTH_LAMBDA   0.1f    /* gentle smoothing factor (vs 0.5 in tri_quality) */
#define QEM_FLIP_MAX_ROUNDS       8       /* post-QEM Surazhsky-Gotsman edge-flip rounds:
                                           * iterate to convergence. One locked pass only
                                           * flips non-adjacent edges, so the heavier 7.5%
                                           * decimation needs several rounds to clear the
                                           * slivers it leaves (flips never move verts). */
#define QEM_LME_MAX_ROUNDS        50      /* safety limit on outer LME round count
                                             (Ozaki & Kanai 2015 algorithm) */
#define QEM_LME_MIN_PROGRESS_FRAC 0.001   /* stop LME-rounds when a single round
                                             collapses fewer than 0.1% of the remaining
                                             faces — converged short of target_nf */

/* Post-weld cleanup (grid_weld terminal step — WeldCleanup_process).
 * The BPA seam bridge leaves slivers / zero-area faces / T-junctions the
 * per-cube QEM never sees (QEM runs before the weld). Flip-first, then guarded
 * collapse. Thresholds mirror seam_audit's --degen census so what the audit
 * flags is what the cleanup targets. */
#define WELD_CLEANUP_SLIVER_MIN_ALT   0.10  /* triangle min-altitude < this (vox) => sliver  */
#define WELD_CLEANUP_DEGEN_AREA       1e-3  /* triangle area < this (vox^2)       => degenerate */
#define WELD_CLEANUP_MAX_COLLAPSE_LEN 5.0   /* never collapse an edge longer than this (vox);
                                             * < the ~7 vox inter-wrap clearance, so a collapse
                                             * can never fuse two scroll wraps */
#define WELD_CLEANUP_FLIP_ROUNDS      8     /* Surazhsky-Gotsman flip rounds (to convergence) */
#define WELD_CLEANUP_COLLAPSE_ROUNDS  6     /* guarded short-edge collapse rounds             */

/* Threading */
#define PARALLEL_VERTEX_CUTOFF 1000  /* don't parallelize below this */
#define PARALLEL_FACE_CUTOFF   2000

/* Halo-overlap grid stitching */
#define HALO_PIN_EPS 1.5f  /* voxels: pin verts within this distance of an
                            * integer-cube_size boundary. This is the FULL
                            * pin radius — verts further into the halo are
                            * LOP-projected and QEM-simplifiable. Wider
                            * halo (halo_voxels) is still loaded for MC
                            * determinism; only the pin set is tight.
                            *
                            * Was halo_voxels + 0.5 (=8.5 at halo=8) until
                            * 2026-05 — that pinned the entire 8-voxel
                            * halo band, leaving a dense ring of
                            * un-collapsible triangles around each seam.
                            * Narrowing to 1.5 vox keeps only the verts
                            * that will actually weld (seam plane ±0.5
                            * MC half-voxel slack on either side). */

/* Zipper — cross-cube T&L chain-stitch (no welding, bridge triangles only) */
#define ZIP_WALL_CAP_TOL          1.0f  /* voxels: cap detected by median Y within tol of wall */
#define ZIP_WALL_CAP_NORMAL_MIN   0.9f  /* |n·ŷ| threshold for a face to count as a Y cap */
#define ZIP_WALL_BAND             6.0f  /* voxels from wall plane to qualify as wall boundary */
#define ZIP_CHAIN_MATCH_DIST_MAX  60.0f /* voxels: reject pair if (X,Z) centroid dist exceeds.
                                           Tolerant because long open-path chains can have
                                           centroids near the cube center, far from any
                                           individual closed loop's centroid in the other cube. */
#define ZIP_CHAIN_MATCH_NDOT_MIN  0.8f  /* reject pair if |avg(n_A) · avg(n_B)| below */
#define ZIP_CHAIN_LEN_WEIGHT      0.1f  /* score weight for |len_A - len_B| */
#define ZIP_CHAIN_NORMAL_WEIGHT   2.0f  /* score weight for (1 - |n_A · n_B|) */
#define ZIP_FOLDOVER_NEAR_VOX     4.0f  /* bbox padding for triangle-tri intersection diagnostic */

#endif
