#ifndef PIPELINE_CONSTANTS_INCLUDED
#define PIPELINE_CONSTANTS_INCLUDED

/* Pre-Step 0 — Voxel denoise (tunnel removal) */
#define TUNNEL_GAP_FILL_THRESH   4   /* min 6-conn FG neighbors to fill BG voxel */
#define TUNNEL_GAP_FILL_ITERS   10   /* max gap-fill passes */
#define TUNNEL_BG_CC_MAX        27   /* max BG CC size to fill (26-conn) */
#define MORPH_CLOSE_ITERS        3   /* multi-pass morph close iterations */
#define SIMPLE_FILL_MIN_FG       4   /* min FG 6-conn neighbors to be candidate */
/* CT-weighted LOP/MLS (experiment, armed by VES_CT_WEIGHTED_LOP=1): the
 * voxel-centre cloud is weighted by raw-CT brightness inside the predicted
 * shell so the projected surface sits on the papyrus centre instead of the
 * shell's own midsurface (measured 2026-09-03: 57% of ribbon vertices sit in
 * the dark gap 2.5 vox inside the sheet, only 20-28% within 1 vox of the
 * sheet centre).  w = FLOOR + (1-FLOOR) * clamp((raw-LO)/(HI-LO), 0, 1)^GAMMA. */
#define EXTRACT_CT_LOP_LO       47.0f
#define EXTRACT_CT_LOP_HI      193.0f
#define EXTRACT_CT_LOP_GAMMA     2.0f
#define EXTRACT_CT_LOP_FLOOR     0.05f
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
#define MIN_CC_SIZE          500    /* voxels - discard smaller components. Also
                                     * the "empty garbage" floor: a cube whose
                                     * LARGEST 6-conn component is below this can
                                     * never mesh, so pred_reject rejects it
                                     * pre-spawn (sibling of the solid-slab reject)
                                     * rather than mesh to a guaranteed empty FAIL. */
#define MAX_COMPONENTS        20    /* keep top N by size */
#define CLEANUP_MICRO_HOLE_MAX 6   /* max boundary loop verts for fill */
#define MIN_FRAGMENT_FACES   100    /* discard sub-components smaller */
#define RESPLIT_MIN_COMP_VERTS 200  /* mesh_resplit: re-mesh a connectivity-component
                                     * only if it has >= this many verts (smaller
                                     * disconnected fragments are dropped) */
#define KIBBLE_AREA_FRAC     0.02f  /* component_cull: drop connectivity-components
                                     * whose surface area is < this fraction of the
                                     * total meshed cube area (post hole-fill/CVT) */

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
 * (tightening to 1.5 vox only fragments). Env SEAM_WIND_TOL overrides.
 *
 * 2026-07-22 raised 0.25 -> 0.35 for the CVT-coarse pipeline. At CVT edge
 * lengths (~8-12 vox, refined to 3.5 at the seam) the within-wrap radial WOBBLE
 * of a single sheet across a seam reaches ~2.85 vox -- above 0.25*9.5=2.4 vox --
 * so the bridge-face FILTER (ball_pivot.c) dropped the seam faces needed to weld
 * two charts of the SAME wrap, leaving intra-sheet splits (wind_audit: the
 * dominant #0<->#2 split, dw=-0.12, unbridged). 0.35 (=3.3 vox) clears the
 * wobble while staying under the 0.40 hard cap that guards the ~0.50 delamination
 * membrane. Faithful grid_weld A/B (8-cube L1 node, full recoarsen+band-CVT):
 * split_comp_pairs 10->5, FULL-TURN fusions 4599->4620 (flat), MANIFOLD; core
 * node (r~99): fusions 572->556 (down), splits unchanged. No fusion cost at mid
 * OR core. See [[project_wind_audit_topology]]. */
#define SEAM_WIND_TOL_DEFAULT_TURNS 0.35

/* Winding-gate HARD cap, in turns: reject a seam bridge whose phase gap
 * exceeds this REGARDLESS of chord direction. The radial-dominance conjunct
 * exempts lateral chords so grazing-seam closures can weld, but that lane
 * also re-admitted lateral CROSS-LAYER stitches: at a delaminated sheet
 * (same sheet, two predicted surfaces ~half a pitch apart -- observed
 * |dw| = 0.50 at the z4608/y3584/x3072 corner, "red-to-pink puddle") the
 * bridge laid a flat membrane between the two layers. Same-wrap grazing
 * offsets stay under ~0.3 turn (seam disagreement <= ~3 vox radial; in-sheet
 * motion in a grazing zone is radial-free), folds sit at ~0, so 0.4 turn
 * separates every legitimate closure from the cross-layer/next-wrap step.
 * Env SEAM_WIND_HARD_TOL overrides. */
#define SEAM_WIND_HARD_TOL_DEFAULT_TURNS 0.40

/* Phase-2 sheet-correspondence re-weld (src/remesh/sheet_reweld.c). After the
 * conservative phase-1 seam bridge, recover legit same-sheet closures it gated
 * (divots) by confirming a 1:1 sheet correspondence across each seam and
 * permissively re-welding each confirmed pair on a cloud restricted to those
 * two sheets. Correspondence evidence = phase-1 bridge VOTES (each phase-1
 * bridge face voting for the sheet-pair it joins) weighted vs geometric near-
 * seam boundary OVERLAP; a pair is confirmed only when it is the mutual-best
 * match with the runner-up below MARGIN_FRAC (one sheet and ONLY one). The pair
 * weld runs with the winding gate OFF (cloud restriction is the safety) and a
 * wider ball (RHO_MAX) to span the divot gaps. */
#define SHEET_REWELD_OVERLAP_R    3.0   /* boundary-overlap match radius (vox) */
#define SHEET_REWELD_VOTE_WEIGHT  2.0   /* one phase-1 vote vs one overlap match */
#define SHEET_REWELD_MIN_SCORE    8.0   /* min combined evidence to confirm */
#define SHEET_REWELD_MARGIN_FRAC  0.5   /* runner-up < this * best, both sides */
#define SHEET_REWELD_RHO_MAX      6.0f  /* wider ball for the permissive weld */
#define SHEET_REWELD_BAND         6.0f  /* seam band (match phase-1) */
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

/* Cut-at-plane trim (Mesh_trim_cut_to_owned_box, the halo-mode default when
 * inset > 0). Vertices/crossings within this of an owned-box plane are clamped
 * onto it: far above the double->float ulp at 128-vox coords (~1e-5, the CVT
 * boundary-site jitter that cost whole coarse rings under the drop-only rule),
 * far below any real feature. Max positional error = this. Also the minimum
 * output edge length the cut can create (crossings nearer an endpoint reuse
 * the endpoint). Env kill switch for the whole cut path: VES_TRIM_CUT=0. */
#define TRIM_CUT_SNAP_EPS_VOX   0.05f

/* Seam weld -- pre-bridge sliver cull. A boundary triangle whose open edge lies
 * in a seam plane is a bridge-front candidate; if it is a near-degenerate sliver
 * (min altitude = 2*area/longest-edge, in voxels, below this) the bridge primes
 * off a needle and emits holes / flipped-normal faces. Drop such triangles before
 * bridging (their other edges become fresh boundary); a post-bridge remesh can
 * refine later. Env override: SEAM_SLIVER_MIN_ALT. 0 disables. */
#define SEAM_SLIVER_MIN_ALT     0.3f

/* Cross-cube seam bridge (no eat-back). The bridge ball radius adapts to the
 * post-CVT boundary spacing: base = BRIDGE_RHO_K * median near-seam boundary-
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

/* Flatten -- quadribbon pipeline (src/tools/quadribbon_pipeline.c) and the
 * solidify stage (src/flatten/quadribbon_solidify.c).  The fit lattice numbers
 * are the ones the 2026-09-01 4x5x5 champion (step11_peel) was produced with;
 * the pipeline must reproduce that certificate and fit byte-identically, so
 * they are constants, not config. */
#define QP_FIT_GRID_U        2.0    /* vox per lattice column (scroll_ribbon --grid-u) */
#define QP_FIT_SLICE_H       2.0    /* vox per lattice row / slice (--slice-h) */
#define QP_FIT_SAMPLE_H      2.0    /* chain sampling step (--sample-h) */
#define QP_FIT_ITERS         1      /* --iters */
#define QP_FIT_FINAL_ITERS   1      /* --final-iters */
#define QP_FIT_METRIC_ITERS  0      /* --metric-iters */
#define QP_FIT_PEEL_LAYERS   4      /* collision depth.  At 10x, 2/3/4 layers
                                     * gave exact primary+extras coverage of
                                     * 96.787/99.575/99.881%; four reaches the
                                     * accepted 4x coverage regime. */
#define QP_FIT_PROMOTE_PEEL_MIN_SHARE 0.005
    /* Collision-peel U runs this large relative to the primary fitted vertex
     * count are real missing sheet-scale coverage, not fragment debris.  Pack
     * them beside the primary atlas; retain smaller runs in *_extras. */
#define QP_FIT_PROMOTE_PEEL_MIN_VERTICES 8192
    /* The relative rule alone scales its cutoff with total block surface area
     * (about 9k vertices at 4x, but about 70k at 10x) and stranded 34 coherent
     * 10x runs / 960,780 recovered vertices in extras.  A physical-size OR
     * criterion keeps the cutoff stable across rungs; tiny runs remain extras. */
#define QP_VERDICT_CORE_RADIUS 64.0 /* ribbon_verdict --core-radius: the shipped
                                     * view excludes the umbilicus core */

/* Solidify: the solid ribbon is FITTED where there is an observation and
 * harmonic (cylindrical residual) where there is not; fills are CAPPED so a
 * large hole never becomes a chord through the adjacent wrap (the 10x marble
 * lesson).  Regime-bound: revalidate per rung. */
#define QS_GOOD_COLUMN_FRAC    0.04  /* a lattice column is GOOD when at least this
                                      * fraction of its rows is observed; the fill
                                      * domain is every row of every column between
                                      * the first and the last good column of a run
                                      * (the V domain is the box's z extent and every
                                      * wrap spans it -- user rule 2026-09-01) */
#define QS_FILL_MAX_U_GAP      64    /* a stretch of non-good columns between two
                                      * good ones wider than this stays empty */
#define QS_SKIP_DR_PITCH       1.5   /* across such a gap -- and across every row or
                                      * column bridge up to the gap reach -- the
                                      * radius extrapolated from each side's local
                                      * slope must land within this many pitches of
                                      * the other side, or the bridge is a wrap
                                      * skip and is refused */
#define QS_ORDER_TOL_FRAC      0.30  /* u is arc length, so a row bridge's angle
                                      * change must be its column span times the
                                      * run's median angle step per column, within
                                      * this fraction (ellipticity) plus
                                      * QS_ORDER_TOL_RAD; a bridge between fragments
                                      * placed out of angular order fails and is
                                      * refused at any length */
#define QS_ORDER_TOL_RAD       0.02
#define QS_SKIP_MIN_ROWS_FRAC  0.10  /* a gap between good columns must share at
                                      * least this fraction of the rows (>= 4) on
                                      * both sides to be verified; unverifiable
                                      * gaps are refused */
#define QS_RUN_MIN_OBSERVED    1024  /* a column run with fewer observed cells is
                                      * confetti the certificate could not place;
                                      * it is not emitted (counted, painted) */
#define QS_RUN_MIN_DENSITY     0.30  /* ... and so is a SMALL run (fewer than
                                      * QS_RUN_CONFETTI_MAX_OBS observed cells)
                                      * whose observed cells cover less than this
                                      * fraction of its rectangle over the rows it
                                      * occupies: 54-86 column slivers at 4-16%
                                      * baked as full-height black bars (user
                                      * review 2026-09-02).  Measured on the 10x3x3
                                      * tube: a density over ALL rows dropped real
                                      * 12-16k-cell arc fragments (coverage 99.1 ->
                                      * 88.3), hence the occupied-row density and
                                      * the size ceiling. */
#define QS_RUN_CONFETTI_MAX_OBS 8192  /* a run with this many observed cells is a
                                      * real fragment whatever its density */
#define QS_RUN_GAP_COLS        2     /* LAYOUT: emitted runs are placed consecutively
                                      * in u with this separator; the atlas gaps
                                      * between runs (unoccupied columns, confetti)
                                      * never reach the sheet.  Runs stay unbridged;
                                      * the original column range is in the report. */
#define QS_HOLE_VETO_FRAC      0.50  /* a hole (connected non-observed region inside
                                      * the fill domain) is refused AS A WHOLE when at
                                      * least this fraction of its cells is vetoed by
                                      * a bridge gate, and filled as a whole otherwise
                                      * (per-bridge vetoes drew comb teeth at hole
                                      * rims: alternating vetoed / filled rows) */
#define QS_HOLE_FUSION_FRAC    0.25  /* ... or when at least this fraction is vetoed
                                      * by a fusion-class gate (material, wrap-skip,
                                      * angular order): topology outranks cosmetics */
#define QS_RUN_MERGE_GAP_COLS  64    /* LAYOUT: a run spans unoccupied column gaps up
                                      * to this wide (= QS_FILL_MAX_U_GAP, the fill
                                      * rule then bridges them); measured 2026-09-02:
                                      * a single empty column split a run and drew a
                                      * black crack through a continuous fill (user
                                      * boxes 6, 7) */
#define QS_FOLD_STEP_FRAC      0.25  /* FOLD: an observed u-step shorter than this
                                      * fraction of du is compressed (the lattice
                                      * folded back; the bake smears) ... */
#define QS_FOLD_COLUMN_FRAC    0.50  /* (unused since the per-cell rule: the right
                                      * cell of a compressed pair is demoted to a
                                      * hole and re-interpolated; user box 1 had 29%
                                      * of its u-steps under 0.6 vox spread across
                                      * columns, so no column ever reached 50%) */
#define QS_DARK_CELL_SNAP_ENABLE 0   /* DARK observed cells (no CT ridge within the
                                      * bake's normal reach: rows that bake as black
                                      * dashes, user boxes 2/4/5).  MEASURED 2026-09-02:
                                      * moving them onto the ridge (<= 4 vox) cleared
                                      * dashes 19 -> 13 but cost DISTORTION 0.018 ->
                                      * 0.147% (4x5x5) and 0.36 -> 0.60% = a failed
                                      * gate (10x3x3): a row moved 3 vox against
                                      * unmoved neighbours is anisotropy by
                                      * definition.  The bake's normal reach carries
                                      * the fix instead (bake.normal_reach_vox); the
                                      * probe stays as the DARK gate. */
#define QS_DARK_U8             40    /* no papyrus at all along the normal (60 = the ridge
                                      * level demoted 22% of observed cells: dim papyrus is
                                      * still papyrus; dim rows are a display matter) */
#define QS_DARK_PROBE_REACH    4.0   /* the gate's probe = the bake's normal reach */
#define QS_DARK_SNAP_REACH     4.0
#define QS_MERGE_PHASE_AGREE   0.75  /* LAYOUT: a same-lineage column gap is bridged
                                      * only when the lifted phase is CONTINUOUS
                                      * across it (the angular-order test, local
                                      * radius) on at least this fraction of the
                                      * rows present on both sides; measured
                                      * 2026-09-02: a lineage-only merge bridged
                                      * runs whose phases jump by a turn and the
                                      * merged holes were refused and baked black */
#define QS_MIDLINE_REACH       6.0   /* MIDLINE test (user rule 2026-09-02: weight
                                      * darkness by the distance to the papyrus
                                      * midline, judged against the neighbours):
                                      * every observed cell probes the CT along its
                                      * lattice normal within this reach for the
                                      * ridge (the midline) ... */
#define QS_MIDLINE_WINDOW      4     /* ... the median ridge offset over the
                                      * (2w+1)^2 observed neighbours is the local
                                      * reference (a systematic face offset is not
                                      * an error) ... */
#define QS_MIDLINE_OUTLIER_VOX 2.0   /* ... a cell is SUSPECT when it bakes DARK (no
                                      * CT >= QS_DARK_U8 within +-1 vox) AND the
                                      * nearest papyrus face along its normal (first
                                      * sample >= QS_FILL_SNAP_MIN_U8) is further
                                      * than this, or when no face is within reach:
                                      * darkness weighted by the distance to the
                                      * papyrus (user rule).  Suspects are demoted
                                      * and refit through their neighbours.  A first
                                      * version compared the ridge MAXIMUM offset
                                      * with the neighbours' median and flagged 33%
                                      * of the sheet: the maximum flips between the
                                      * two bright faces (p50 offset 3 vox). */
#define QS_MIDLINE_ISOLATED_MAX 3    /* a suspect is REFIT (demoted, interpolated from
                                      * its neighbours) only when at most this many of
                                      * its 8 neighbours are suspects too: a dash or a
                                      * stray row.  Dense suspect regions are the
                                      * fit's systematic depth error (measured
                                      * 2026-09-02: 20% of observed cells sit dark
                                      * and > 2 vox from papyrus; demoting them all
                                      * shattered coherence 96 -> 78%): reported,
                                      * left in place, owned by the shell-dedup
                                      * campaign */
#define QS_DARK_DEMOTE         0     /* MEASURED 2026-09-02: refitting even the ISOLATED
                                      * midline suspects through their neighbours raised
                                      * depth seams 0.31 -> 1.57% (4x5x5) and 0.64 -> 2.83%
                                      * (10x3x3) with more dashes, because the neighbours
                                      * themselves sit ~3 vox off the papyrus (dense suspects
                                      * 17%): the MIDLINE test stays a GATE + report, the
                                      * refit is off until the fit's depth error is fixed
                                      * upstream.  1 = refit isolated suspects. */
#define QS_DARK_DEMOTE_DOC     1     /* an observed cell with no CT ridge within
                                      * QS_DARK_PROBE_REACH along its normal (it
                                      * would bake as a black dash) is demoted to a
                                      * hole and interpolated from its neighbours,
                                      * which do sit on papyrus; provenance FILL */
/* --- RIDGE TRACK (src/flatten/ridge_track.c): put the ribbon ON the papyrus ---
 * MEASURED 2026-09-03: the source mesh itself sits p50 2.5 vox from the brightest
 * papyrus within 4 vox of it (63% of its vertices read >= 70), and the solid ribbon
 * inherits that.  The reach is kept at or under the verdict's COVERAGE radius (3.0
 * vox) so a tracked cell stays within it of the mesh it came from.
 * (python/scripts/ribbon_wrap_tracker.py): the solid ribbon
 * sits p50 3-6 vox from the nearest CT wrap ridge, in a stack whose local ridge
 * spacing is 15-17 vox, and 0% of the samples were FUSED (every one had a dark
 * gap either side).  So the surface runs on a wrap's shoulder or in the gap and
 * drifts onto the next wrap over a run -- the reviewer's swirl.  The tracker
 * chooses a per-cell offset along the cell normal that lands on a ridge, stays
 * near the fit, and varies smoothly over the lattice (so neighbours pick the
 * SAME ridge).  An independent per-cell snap was measured worse in 2026-09-02
 * (`QS_FILL_SNAP_ENABLE`); the smoothness term is the difference. */
#define QS_RIDGE_SNAP_ENABLE   1     /* PROMOTED 2026-09-03 after a measured A/B on both
                                      * rungs, same fit, same exe.  It puts the ribbon ON
                                      * the papyrus and is GATE-NEUTRAL:
                                      *   4x5x5   on-papyrus 57.2 -> 81.9%, CT at the cell
                                      *           p50 92 -> 128; verdict 1/8 both ways;
                                      *           coverage 90.5 -> 90.3, coherence 96.2 both
                                      *   10x3x3  on-papyrus 46.3 -> 65.7%, CT 53 -> 83;
                                      *           verdict 1/8 both ways; coverage 99.2 ->
                                      *           99.1, distortion 0.436 -> 0.407
                                      * The only cost is TEXTURE/depth-seams (0.35 -> 0.48%
                                      * and 0.66 -> 0.82%), the gate that already failed.
                                      * By eye the 1:1 bake crops of the reviewer's boxes
                                      * are brighter and sharper (ridge_ab_final.png).
                                      * Evidence: output/pherc0139_*_ridge_final_20260903. */
#define QS_RIDGE_EDGE_LIMIT    5.7   /* no displaced lattice edge may exceed this (the
                                      * verdict's TANGLE/long-edges gate is 6.0 and
                                      * tolerates none).  Offsets on a violating edge
                                      * are shrunk until it fits. */
#define QS_RIDGE_REACH         2.5   /* half-width of the offset search, vox: under
                                      * half the measured 15-17 vox ridge spacing,
                                      * so a cell cannot reach the next wrap */
#define QS_RIDGE_STEP          0.25  /* offset discretisation.  With MAX_STEP 1 this is also
                                      * the largest offset difference between neighbours, so
                                      * it sets how much depth step the tracker can add: at
                                      * 0.5 it cost TEXTURE/depth-seams 0.35 -> 0.48%. */
#define QS_RIDGE_PRIOR         0.030 /* cost per vox of |offset| (stay near the fit) */
#define QS_RIDGE_SMOOTH_U      0.250 /* cost per vox of offset change along u.  At 0.060
                                      * (measured 2026-09-03) neighbours chose ridges up to
                                      * 12 vox apart: on-papyrus rose 37 -> 74% but the
                                      * verdict went 1 -> 4 failed gates (TANGLE/long-edges,
                                      * TANGLE/cross-wrap).  The data term spans [0,1], so a
                                      * weight of 0.25/vox makes a 4-vox disagreement cost
                                      * more than any data gain. */
#define QS_RIDGE_SMOOTH_V      0.250 /* ... and across v */
#define QS_RIDGE_MIN_U8        70.0  /* CT value that counts as papyrus */
#define QS_RIDGE_PROFILE_WIN   0     /* cells (odd; 0 or 1 = off) over which a cell's CT
                                      * profile is AVERAGED before the tracker chooses its
                                      * offset.  MEASURED 2026-09-03: the brightest voxel
                                      * along a normal flips between a wrap's two fibre
                                      * layers from cell to cell, so the per-cell ridge
                                      * field steps over 0.75 vox on 31-35% of neighbour
                                      * pairs while the mesh does so on 0.24-0.34% -- the
                                      * CT ridge is ~100x rougher than the surface it is
                                      * meant to define.  Averaging the profile over a few
                                      * cells first gives a ridge that is a SURFACE fitted
                                      * to a neighbourhood rather than a per-cell maximum,
                                      * so the tracker can follow the papyrus without
                                      * fighting the smoothness terms.
                                      * MEASURED 2026-09-03 and NEUTRAL: windows 5 and 11
                                      * give seams 0.188% and 0.190% against 0.189% with
                                      * none.  The row DP's own L1 smoothness already
                                      * picks a coherent path through the flicker, so
                                      * pre-averaging the profile changes nothing.  Left
                                      * off; the finding it came from (the CT ridge is
                                      * ~100x rougher than the mesh) still stands and is
                                      * why re-meshing is NOT the lever. */
#define QS_RIDGE_FLATTEN       0.05  /* FLATTEN term: cost per voxel of the cell's depth
                                      * DEVIATION FROM ITS NEIGHBOURHOOD, measured exactly
                                      * as the verdict's TEXTURE/depth-seams gate measures
                                      * it -- radius minus a 9-cell local median.  That
                                      * gate exists ("BORN FAILING by design") to be driven
                                      * down by a claim-depth fix: blocks of cells sitting
                                      * 1-2 vox proud of their neighbours because the claim
                                      * DP hopped between the front and back face of the
                                      * prediction shell.  With this term the tracker
                                      * chooses, among the offsets that keep a cell on
                                      * papyrus, the one that also puts it level with its
                                      * neighbourhood -- so it repairs the plateaus instead
                                      * of inheriting them.  0 disables.
                                      * SWEPT 2026-09-03 on the 4x5x5 (depth-seams /
                                      * on-papyrus): 0 -> 0.400% / 79.8%, 0.05 -> 0.368% /
                                      * 75.5%, 0.15 -> 0.341% / 64.8%, 0.30 -> 0.336% /
                                      * 60.0%, 0.60 -> 0.333% / 58.1%; the champion without
                                      * the tracker is 0.348% / 57.2%.  Seams asymptote at
                                      * 0.333% -- the gate wants 0.100%, so the rest is the
                                      * FIT's claim-depth, not reachable from this stage.
                                      * 0.05 keeps the texture win at champion-level seams.
                                      * Relaxing QS_RIDGE_MAX_STEP to 2/6/10 does NOT help
                                      * (0.350 / 0.359 / 0.362): the smoothness weight, not
                                      * the bound, is what holds a plateau together. */
#define QS_RIDGE_FLATTEN_WIN   9     /* local median window, cells (matches the gate) */
#define QS_RIDGE_ITERS         3     /* row DP rounds (round 2 on sees the v pull) */
#define QS_RIDGE_MAX_STEP      1     /* HARD bound: adjacent cells' offsets may differ by
                                      * at most this many QS_RIDGE_STEPs.  An unbounded L1
                                      * smoothness (measured 2026-09-03, even at 0.25/vox)
                                      * still let two neighbours pick offsets 8 vox apart
                                      * where the data was emphatic, which broke
                                      * TANGLE/long-edges and TANGLE/cross-wrap.  With this
                                      * bound an emitted edge can grow by at most
                                      * 2 * STEP * MAX_STEP voxels, so the tangle gates are
                                      * safe by construction, and the field can still drift
                                      * a whole wrap over 20 cells. */

#define QS_RIDGE_FILL_REACH    6.0   /* DARK FILLS ARE RECOVERABLE (user, 2026-09-08): a fill
                                      * cell has no observation to stay near, so its search
                                      * half-width along the normal is this, not QS_RIDGE_REACH.
                                      * Under half the 9.5-17 vox layer gap so a fill cannot reach
                                      * the next wrap; the row DP's MAX_STEP bound (0.25 vox per
                                      * cell) still makes the move coherent: a fill 14 cells from
                                      * its rim can move the full width, the rim itself cannot.
                                      * MEASURED 2026-09-08 (same fit, same exe, reach 2.5 -> 6):
                                      *   assembled F0 4x5x5 ribbon: fills on papyrus 39.2 ->
                                      *   48.8%, dark bake fills 68.6k -> 67.0k of 488k, seams
                                      *   0.166 -> 0.169%, coverage 89.5 both, tangle gates hold
                                      *   (max edge 6.00);
                                      *   champion F0 fit ribbon: fills 36.2 -> 43.4%, seams
                                      *   0.438 -> 0.449%, coverage 98.9 both.
                                      * Of the dark fills 70% have no source vertex within 10
                                      * vox (a hole in the prediction, not a missed recto):
                                      * those are not recoverable by any local search.
                                      * Evidence: output/assemble_f0_4x5x5_20260907/solid_fillreach,
                                      * output/quadribbon_f0_4x5x5_20260908/stage4_solid{,_prefillreach}. */
#define QS_RIDGE_FILL_PRIOR    0.0   /* fills pay no cost for |offset| (observed cells: QS_RIDGE_PRIOR) */

/* --- EXACT LOCAL FRAME -------------------------------------------------------
 * The lane slices at z = const and measures radius in that plane.  MEASURED
 * 2026-09-03 (python/scripts/axis_frame_compare.py): the papyrus is
 * perpendicular to the LOCAL AXIS TANGENT, not to +Z -- area-weighted
 * |cos(mesh normal, tangent)| 0.121 on the 10x3x3 and 0.123 on the 4x5x5,
 * against |cos(normal, +Z)| 0.201 and 0.181; the local frame is the better wrap
 * frame over 76-79% of mesh area.  A z-plane therefore cuts a tilted wrap in an
 * ellipse and its apparent radius is wrong by r*(1/cos t - 1): measured
 * |r_plane - r_perp| p90 6.1 vox on the 4x5x5 and 15.8 vox on the 10x3x3
 * (1.66 pitches, where wraps are 9.5 apart).  Straightening does NOT fix this --
 * it removes the translation, not the tilt.
 * With this armed the solid stage measures radius PERPENDICULAR to the local
 * tangent, so the wrap-skip, radial-order, core and RAMP tests all see the
 * true radius.  Needs the axis curve (geometry.axis_table); no-op without it. */
#define QS_COLUMN_JUMP_PITCH   1.0   /* COLUMN-CONTINUITY repair: a lattice column is one
                                      * angular position (MEASURED 2026-09-03: a column
                                      * holds its angle to p50 0.000 and p90 0.017 vox of
                                      * arc per vox of z), so its radius must vary
                                      * smoothly down the rows.  A cell whose radius
                                      * differs from its column's local median by more
                                      * than this many pitches is ON ANOTHER WRAP -- the
                                      * certificate's lift changed between adjacent z --
                                      * and is demoted so the fill bridges its rows
                                      * instead.  This is what refuses 284,773 row
                                      * neighbours on the 10x10x10 (they sit at the same
                                      * angle, 2 vox apart in z, 27 vox apart in radius)
                                      * and shatters its largest run into 34,908
                                      * components.  0 disables. */
#define QS_COLUMN_JUMP_WIN     15    /* rows in that local median (odd) */
#define QS_COLUMN_RELOCATE     0     /* instead of DEMOTING a wrong-wrap cell, move it
                                      * radially onto its column's local median -- the
                                      * radius its neighbours agree on.  Demoting costs
                                      * coverage (the cell stops representing its source
                                      * vertex: 75.7% -> 72.7% on the 10x), and the peel
                                      * band cannot give the coverage back (measured: only
                                      * 0.5% of demoted cells have a correct-wrap twin
                                      * there).  But the correct-wrap MATERIAL is present
                                      * -- the wraps are separated, measured 0% fused -- so
                                      * a cell moved to the column median lands near real
                                      * material of the right wrap and may keep the
                                      * coverage that demotion loses.  0 = demote.
                                      * MEASURED 2026-09-03 and NOT EFFECTIVE AS WRITTEN:
                                      * it writes `field[]`, but the build reads observed
                                      * positions through `sel[]` into the input mesh, so
                                      * the cell does not actually move.  The 10x run with
                                      * it on reads coverage 75.8%, 34,967 components,
                                      * seams 2.902% -- the same as NO repair (34,908
                                      * components), which is the tell.  Relocating for
                                      * real needs an override array for observed
                                      * positions, and it is worth questioning first:
                                      * placing a cell at a computed radius INVENTS an
                                      * observation, which is what the fill already does
                                      * honestly after a demotion. */
#define QS_OBS_SMOOTH_VOX      0.0   /* DEPTH SMOOTHING of the emitted lattice: bounded
                                      * Jacobi smoothing of the per-cell radius, no cell
                                      * moving more than this many voxels in total.
                                      * MEASURED 2026-09-03: after the ridge tracker, the
                                      * column repair and the row-trend smoothing, 84-89%
                                      * of the remaining depth seams join TWO REAL
                                      * OBSERVATIONS -- the extracted surface's own
                                      * roughness at the 2-cell scale, which is exactly
                                      * what the gate's 9-cell high-pass measures.  The
                                      * ridge tracker cannot remove it because its CT term
                                      * defends every local ridge; a small unconditional
                                      * smoothing can, and a bound well inside the
                                      * verdict's 3 vox COVERAGE radius keeps every cell
                                      * on its source vertex.  0 disables.
                                      * MEASURED 2026-09-03 and REJECTED: bounds 0.30 and
                                      * 0.75 vox give seams 0.278% and 0.829% against
                                      * 0.189% with no smoothing, and the verdict goes
                                      * 1/8 -> 2/8.  The BOUND is the problem: a cell that
                                      * wants to move further is clamped, and the clamp
                                      * boundary is a new step.  That is the third repair
                                      * to fail this way (fill cap, depth-outlier
                                      * demotion, this) -- every threshold-based repair
                                      * creates steps at its own threshold.  What works is
                                      * continuous: the ridge tracker's flatten term and
                                      * the row-trend median. */
#define QS_OBS_SMOOTH_ITERS    8     /* Jacobi passes */
#define QS_OBS_SMOOTH_LAMBDA   0.5   /* per-pass step toward the neighbour mean */
#define QS_ROW_TREND_SMOOTH    17     /* rows in a median smoothing of the fill's per-row RADIUS
                                      * baseline (odd; 0 or 1 = off).  That baseline is the
                                      * median radius of each row's own observations, so two
                                      * adjacent rows observing different columns can differ by
                                      * more than the seam gate's 0.75 vox, and the fills that
                                      * ride on it inherit the jump.  MEASURED 2026-09-03 on the
                                      * 4x5x5: fill-fill pairs carry 29% of the depth seams
                                      * along v (3,707 of 12,724) even though a harmonic solve
                                      * should be smooth.  A wrap's radius varies smoothly with
                                      * z, so the baseline may be smoothed; the harmonic
                                      * residual still reproduces the observations exactly at
                                      * the cells that have them. */
#define QS_COLUMN_JUMP_2D      1     /* also test across COLUMNS (a separable median down the
                                      * rows and then along the columns).  A wrong-wrap cell is
                                      * a pitch off its neighbourhood in BOTH directions; the
                                      * column-only test catches it only where the column has
                                      * enough rows.  0 = column only. */

#define QS_ARCLEN_V            1     /* the SHEET's v is arc length along the umbilicus rather
                                      * than world z.  The lattice rows are perpendicular to the
                                      * local axis, so a v-edge spans sqrt(1+slope^2) voxels of
                                      * surface per voxel of z: MEASURED 2026-09-03, the ribbon's
                                      * v-edge stretch tracks that band for band (1.0123 at slope
                                      * < 0.15 up to 1.0674 at 0.35-0.45), so one sheet is
                                      * stretched 1% at its upright end and 7% at its steepest --
                                      * a 5.5-point differential that NO gate sees, because none
                                      * measures v scale.  Arc length is a function of z alone,
                                      * so this is a 1-D monotone remap of the row coordinate,
                                      * applied ONLY where the sheet is composited
                                      * (`quadribbon_strips --axis-table`): the lattice, the
                                      * verdict and every sidecar keep v = z. */
#define QS_PERP_RADIUS         1     /* measure radius perpendicular to the axis tangent */

#define QS_DEPTH_OUTLIER_VOX   0.0   /* DEPTH-OUTLIER demotion: an observed cell whose
                                      * radius differs from its 9-cell local median by
                                      * more than this is demoted to a hole, so the
                                      * harmonic fill replaces it with the surface
                                      * through its neighbours.  This is the same
                                      * high-pass the verdict's TEXTURE/depth-seams gate
                                      * uses, applied as a repair instead of a report:
                                      * a block of cells the claim DP put a voxel or two
                                      * proud of its surroundings (the front/back faces
                                      * of the prediction shell) becomes fill and the
                                      * plateau goes away.  0 disables.
                                      * MEASURED 2026-09-03 and REJECTED: thresholds 0.75
                                      * and 0.50 vox give depth-seams 0.464% and 0.484%
                                      * against 0.368% with no demotion.  Replacing an
                                      * outlier with fill does not remove the step, it
                                      * moves it to the patch BOUNDARY, and the demoted
                                      * cells shrink the edge count they are measured
                                      * against.  The repair has to pick a different
                                      * CLAIM, not delete the cell. */
#define QS_DEPTH_OUTLIER_WIN   9     /* its local median window, cells */
#define QS_GATE_RAMP_PITCHES   3.0   /* RAMP gate: within one run, u is arc length
                                      * and the phase is the geometric angle, so
                                      * the ribbon must gain ONE pitch of radius
                                      * per turn.  A run that gains more than this
                                      * many pitches per turn is walking ACROSS the
                                      * wrap stack -- the reviewer sees it as a
                                      * swirl / smear.  Measured 2026-09-03 on the
                                      * 10x3x3 tube: the big material islands gain
                                      * 1.0-1.5 pitches per turn (correct), while
                                      * seven islands at r 235-475 gain 19-49
                                      * (they follow mesh where the prediction has
                                      * FUSED the wraps: the source radius
                                      * histogram there is solid, comb period 4 vox
                                      * instead of 9.5). */
#define QS_GATE_RAMP_MAX       0.10  /* ... and at most this fraction of observed
                                      * cells may sit in such runs */
#define QS_RAMP_MIN_TURNS      0.15  /* a run must span at least this many turns
                                      * before its ramp is measured.  Below that
                                      * the least-squares slope is noise (measured
                                      * 2026-09-03: runs under 0.15 turns report
                                      * 33-117 pitches/turn on BOTH rungs), and the
                                      * ramp is confounded by the wraps' own
                                      * non-circularity -- PHerc0139's wraps deviate
                                      * by ~10% of r, which at r 450 is 4.7 pitches
                                      * over a third of a turn.  RAMP is therefore a
                                      * SCREEN: a flagged run is confirmed by the
                                      * source-mesh radius comb (fused wraps show a
                                      * solid histogram, period 4 vox instead of the
                                      * 9.5 pitch) and by a cross-section. */
#define QS_GATE_SCALE_LOG2     0.4   /* SCALE gate: a lattice triangle whose 3-D
                                      * area differs from its (u,v) area by more
                                      * than this many octaves (0.4 = 1.32x) is
                                      * SMEARED: the bake reads a stretched or
                                      * crumpled patch of CT into one pixel.
                                      * The verdict's DISTORTION gate measures
                                      * anisotropy (sigma1/sigma2) and is BLIND to
                                      * uniform scale, which is why the reviewer's
                                      * 2026-09-03 smear boxes passed it.
                                      * Measured on the reviewer's boxes: inside
                                      * p90 1.7, outside p90 1.18; observed cells
                                      * p90 0.27-0.29 on both rungs. */
#define QS_GATE_SCALE_MAX      0.08  /* ... and at most this fraction of emitted
                                      * triangles may exceed it (4x5x5 measured
                                      * 4.5%, 10x3x3 tube ~15%) */
#define QS_FILL_SUPPORT_CELLS  64    /* FILL-SUPPORT gate: a fill cell this many
                                      * lattice steps from the nearest observation
                                      * is unsupported invention.  Measured
                                      * 2026-09-03: 4x5x5 fills p50 35 p90 126,
                                      * 10x3x3 tube p50 66 p90 297 max 642 (a
                                      * column with QS_GOOD_COLUMN_FRAC of its
                                      * rows observed is filled over EVERY row).
                                      * Reported, not capped: the cap is a
                                      * separate measured decision. */
#define QS_FILL_MAX_DIST_CELLS 0     /* CAP: a fill cell further than this many
                                      * lattice steps from the nearest observation
                                      * is not emitted (0 = no cap, today's
                                      * behaviour).  The fill domain rule fills
                                      * every row of a GOOD column, and a column is
                                      * good at QS_GOOD_COLUMN_FRAC (4%) of its
                                      * rows, so a column observed over 26 of 639
                                      * rows is invented over the other 613.  Those
                                      * cells land on papyrus as often as the
                                      * observations do (measured), but they carry
                                      * 91% of the SCALE violations on the 10x3x3.
                                      * MEASURED 2026-09-03 at 64 and REJECTED:
                                      * 4x5x5 fills 420,680 -> 297,453 (-29%) with
                                      * SCALE 2.47% -> 2.52% (no better); 10x3x3
                                      * fills 2,097,910 -> 1,049,948 (-50%) with
                                      * SCALE 9.89% -> 8.33% (still FAIL) -- the
                                      * scale violations are NOT concentrated in the
                                      * far fills, so the cap only trades marbling
                                      * for black.  The reviewer said the same by
                                      * eye. Evidence: output/pherc0139_*_fillcap64_20260903. */
#define QS_GATE_FILL_FAR_MAX   0.35  /* ... and at most this fraction of fill cells
                                      * may be further than that */
#define QS_GATE_FOLD_MAX       0.01  /* solid gates (fractions of observed cells /
                                      * of hole cells): fold-demoted, dark and not
                                      * recovered, refused hole area */
#define QS_GATE_DARK_MAX       0.05  /* MIDLINE gate: fraction of observed cells that are
                                      * ISOLATED suspects (dashes; measured 2.9% on the
                                      * 4x5x5, 3.7% on the 10x3x3 -- they are refit);
                                      * the dense suspect fraction is reported, not gated */
#define QS_GATE_REFUSED_HOLE_MAX 0.20  /* measured 12-18% on the accepted state */
#define QS_FILL_SNAP_ENABLE    0     /* MEASURED 2026-09-02 on the 4x5x5 + 10x3x3
                                      * (smooth axis, layout): the snap moved 97% of
                                      * fills by p50 1.0 vox but the wide black fill
                                      * bands stayed black (CT max along the normal
                                      * ~59 within 6 vox: dim or absent papyrus) and
                                      * depth seams rose 0.37 -> 1.69% / 0.56 -> 5.5%;
                                      * no visible gain -> OFF (champion rule). 1
                                      * re-arms it for an A/B. */
#define QS_FILL_SNAP_REACH     6.0   /* CT-RIDGE SNAP of FILL cells (measured
                                      * 2026-09-02: fills baked black lay 2-6 vox off
                                      * the papyrus with bright CT within reach along
                                      * the lattice normal for 61-90% of them): a fill
                                      * moves along its normal to the NEAREST voxel
                                      * at or above QS_FILL_SNAP_MIN_U8 within this
                                      * reach (the neighbouring wrap's face is >= 8
                                      * vox away), refined to the local ridge maximum
                                      * within QS_FILL_SNAP_RIDGE_HALF; observed cells
                                      * never move */
#define QS_FILL_SNAP_MIN_U8    60
#define QS_FILL_SNAP_RIDGE_HALF 1.5
#define QS_FILL_SNAP_SMOOTH_ITERS 3  /* Jacobi passes on the snap displacement over
                                      * fill cells (observed neighbours pin 0) */
#define QS_WRAP_GATE_VOX       6.0   /* mirrors RIB_WRAP_GATE (ribbon.c) and
                                      * RV_WRAP_GATE (ribbon_verdict.c): no emitted
                                      * edge may span the inter-wrap gap.  ribbon.c
                                      * is not linkable into the pipeline. */
#define QS_CORE_RADIUS_VOX     64.0  /* no fill inside the umbilicus core */
#define QS_BRIDGE_MAX_CELLS    8     /* a row bridge between two observations at
                                      * most this many cells apart is checked against
                                      * the fit's chord gates (RIB_UFILL_MAX reach) */
#define QS_FILL_STRETCH        1.25  /* ... its rim-to-rim chord must be at most this
                                      * times the lattice arclength it replaces
                                      * (RIB_FILL_STRETCH): a chord to the adjacent
                                      * wrap fails it and VETOES the row's cells */
#define QS_FILL_CHORD_ABS      24.0  /* ... and at most this long in voxels
                                      * (RIB_FILL_CHORD_ABS) */
#define QS_BBOX_MARGIN_VOX     32.0  /* a fill may extend this far beyond the input
                                      * ribbon's bounding box: the RAW is read
                                      * anywhere and the fit ribbon stops short of
                                      * the mesh box faces (6 vox refused 70,912
                                      * bridged-gap fills on the 4x5x5, 2026-09-02) */
#define QS_BAND1_ADOPT_DR_PITCH 0.35 /* a peel-band observation is the SAME wrap as
                                      * its band-0 neighbourhood when its radius is
                                      * within this many pitches (same-sheet wobble
                                      * stays <= 0.3 turn; a ply sits >= 0.5) */
#define QS_BAND1_NEIGHBOUR_CELLS 8   /* search window (cells) for that band-0
                                      * reference radius */
#define QS_LATTICE_EPS         1e-3  /* |u/du - round| tolerance when recovering
                                      * the lattice from UV */
#define QS_PHASE_MISMATCH_MAX  1.5707963 /* sidecar-vs-atan2 lifted-phase mismatch
                                      * above which an observation is reported as
                                      * phase-inconsistent (pi/2) */
#define QS_OPTIMIZE_ROUNDS     2     /* untangle -> refit rounds in stage 5 */
#define QS_OPTIMIZE_MIN_GAIN   0.05  /* stop when long conflicts drop by less */
#define QS_REPARAM_FINAL       0     /* 1 = scroll_ribbon --preserve-input-topology
                                      * on the optimized ribbon (metric re-solve;
                                      * the installed certificate gauge is retained
                                      * by default -- measured lever) */

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

/* Step 0 — CVT/RVD variational remesher. The CWF energy converges
 * by ~12-15 iterations with the decaying lambda_CVT, well short of the standalone
 * default of 50; the pipeline trades the tail for throughput. */
#define CVT_PIPELINE_ITERS          12    /* Lloyd/CWF iterations per component       */
#define CVT_WIND_TOL_DEFAULT      0.30    /* max branch-free turn span on any edge     */
#define CVT_TARGET_SPACING_PITCH  0.44    /* initial uniform edge spacing / wrap pitch */
#define CVT_WIND_MAX_RETRIES         4    /* winding shortcut: density may cure it     */
#define CVT_TOPOLOGY_MAX_RETRIES     1    /* one 2x attempt, then retain source chart */
#define CVT_MIN_FACES_FOR_REMESH    400    /* retain already-small dense charts        */

/* GRADED density (the pipeline default): the CVT target edge length varies with
 * depth into the owned box -- fine (weldable) near the cube faces, coarse in the
 * interior. The seam edge length is DERIVED from the seam-weld bridge cap so seams
 * close by construction: a bridge triangle survives only if its longest edge
 * <= 2*BRIDGE_RHO_MAX (=6 vox, pinned by the 7-vox inter-wrap clearance), so the
 * seam band targets well under that. This replaces hand-tuning one global ratio to
 * straddle the cap. See src/remesh/sizing_field.h. */
#define CVT_SEAM_EDGE     (1.2f * BRIDGE_RHO_MAX)  /* ~3.6 vox: seam-band target edge, DERIVED;
                                                    * invariant CVT_SEAM_EDGE < 2*BRIDGE_RHO_MAX. */
#define CVT_INTERIOR_EDGE   13.0f  /* vox: deep-interior target edge (the one coarseness knob) */
#define CVT_SEAM_BAND        2.0f  /* vox: flat full-density pad inward from each owned-box plane.
                                    * THIN RIM (was 6.0): just past the 1-vox trim inset, so only
                                    * the ring that survives to become the seam is dense. */
#define CVT_SEAM_RAMP       10.0f  /* vox: smoothstep ramp width seam->interior (was 30.0; the
                                    * 6+30 band covered most of a core-dense cube -> 1.32M faces).
                                    * |grad h| ~0.94 at 10 vox -- steeper than the documented 0.47
                                    * comfort point, so watch ramp-ring min-angle in the A/B;
                                    * ~8-12 is the floor before ramp-ring quality degrades. */

#define CVT_TARGET_RATIO         0.0025f  /* Legacy/diagnostic uniform-density ratio.
                                           * Scroll runs with an axis and pitch derive their
                                           * default site count from physical area and
                                           * CVT_TARGET_SPACING_PITCH instead.  An explicit
                                           * --cvt-ratio or VES_CVT_RATIO overrides that default
                                           * for density ablations.  Axis-free callers retain
                                           * this ratio as a compatibility fallback. */
#define CVT_MIN_SITES               50    /* floor: a component needs at least this many
                                           * generators to form a valid RVD dual; small sheets
                                           * never simplify below it. */
/* Post-weld cleanup (grid_weld terminal step — WeldCleanup_process).
 * The BPA seam bridge leaves slivers / zero-area faces / T-junctions the
 * per-cube CVT never sees (CVT runs before the weld). Flip-first, then guarded
 * collapse. Thresholds mirror seam_audit's --degen census so what the audit
 * flags is what the cleanup targets. */
#define WELD_CLEANUP_SLIVER_MIN_ALT   0.10  /* triangle min-altitude < this (vox) => sliver  */
#define WELD_CLEANUP_DEGEN_AREA       1e-3  /* triangle area < this (vox^2)       => degenerate */
#define WELD_CLEANUP_MAX_COLLAPSE_LEN 5.0   /* never collapse an edge longer than this (vox);
                                             * < the ~7 vox inter-wrap clearance, so a collapse
                                             * can never fuse two scroll wraps */
#define WELD_CLEANUP_FLIP_ROUNDS      8     /* Surazhsky-Gotsman flip rounds (to convergence) */

/* Weld-time seam-band refinement (SeamRefine_process, grid_weld pre-bridge).
 * Makes uniform-coarse CVT cubes bridgeable without per-cube dense rims: the
 * band next to each detected seam plane is subdivided (+flipped) into
 * well-shaped ~target triangles. The target is pinned by the bridge's priming
 * geometry: a hinge ball reconstructs only when the front triangle's
 * CIRCUMRADIUS <= rho <= BRIDGE_RHO_MAX (ball_pivot.c sphere_center), and a
 * near-equilateral triangle of edge e has circumradius e/sqrt(3), so
 * 3.5/1.73 ~ 2.0 <= the adaptive rho (~2.6 at a 3.5-vox front median). Also
 * matches CVT_SEAM_EDGE (~3.6), the rim size the graded runs PROVED weldable.
 * NOTE plain boundary-midpoint splitting alone cannot work: fan children keep
 * the far apex and their circumradius stays ~5 > rho at every depth -- the
 * interior band splits + flips are what make the band well-shaped. */
#define SEAM_REFINE_TARGET_VOX      3.5f  /* band edge-length target (vox) */
#define SEAM_REFINE_MAX_ROUNDS      8     /* split-round cap: 13 -> 3.5 needs 2
                                           * halvings, but the independent-set
                                           * face locks throttle a face to one
                                           * split/round; 8 converges with margin */

/* Post-recoarsen seam-band CVT beautification (SeamBandCvt_process). The
 * guarded collapse leaves the band collapse-scarred and anisotropic next to the
 * blue-noise CVT interior; this re-meshes each band patch with the SAME
 * CVT/RVD engine, pinned-boundary (junction ring + hole rims immutable,
 * bit-exact -> conforming stitch, open boundaries unchanged) and FAIL-CLOSED
 * (per-patch conformity + Euler/manifold/connectivity gates; a rejected patch
 * keeps its recoarsen-quality geometry). Runs AFTER recoarsen: the coarse band
 * is a far cheaper RVD input, and the collapse pass doubles as the fallback. */
#define SEAM_BANDCVT_BAND_VOX      14.0   /* face eligibility: any vert within this of a
                                           * detected plane (== recoarsen band -> the CVT
                                           * sees exactly the collapse-touched region) */
#define SEAM_BANDCVT_H_VOX         13.0   /* band target edge == CVT_INTERIOR_EDGE, so the
                                           * band matches the per-cube interior density and
                                           * the weld disappears visually */
#define SEAM_BANDCVT_ITERS         10     /* Lloyd/CWF iterations per patch (12 is the
                                           * per-cube default; bands are simpler strips) */
#define SEAM_BANDCVT_MIN_FACES     64     /* leave dust patches alone */
#define SEAM_BANDCVT_TILE_VOX      96.0   /* in-plane tile size: the fail-closed blast
                                           * radius. Untiled, the band lattice is ONE
                                           * grid-spanning patch and any single gate trip
                                           * rejects the whole band (measured). ~Cube-scale
                                           * tiles keep patches cheap and rejections local
                                           * (a core fold rejects its tile only). */
#define SEAM_BANDCVT_FINE_FRAC     0.10   /* process a patch only if >= this fraction of
                                           * its edges is < 0.6*h: pass B (half-tile
                                           * offset) then re-meshes only pass A's border
                                           * strips + rejects instead of churning accepted
                                           * CVT interiors. */

/* Post-bridge seam-band recoarsening (WeldCleanup_recoarsen_seam). The weldable
 * fine seam band (thin graded CVT rim now; weld-time refinement later) is only
 * needed WHILE bridging + filling; once the seam is closed the band is collapsed
 * back toward the coarse budget. Length-based, band-restricted guarded collapse:
 * same link-condition / normal-flip / boundary-preserving guards as the sliver
 * cleanup, so open boundary loops (grid edges, the next hierarchical level's
 * seams) exit BIT-IDENTICAL -- that boundary preservation is the hierarchical
 * composability guarantee. Fusion safety does NOT come from the band (which
 * only SELECTS candidates spatially): it comes from collapse contracting
 * EXISTING edges only, plus the WELD_CLEANUP_MAX_COLLAPSE_LEN cap (5 < 7-vox
 * inter-wrap clearance). The band settles at ~5-10 vox edges, not the 13-vox
 * interior -- do NOT raise the cap to chase the last factor. ORDERING is
 * load-bearing: recoarsen runs AFTER every hole closer (fills gate on extent;
 * a pre-fill recoarsen coarsened slit rims past the gates: 9 -> 188 open seam
 * loops on the 4x5x5). */
#define WELD_RECOARSEN_BAND_VOX      14.0   /* verts within this of a detected seam plane are
                                             * candidates. Covers the graded rim + ramp
                                             * (CVT_SEAM_BAND 2 + CVT_SEAM_RAMP 10 + margin);
                                             * harmless when wider -- interior edges (~13 vox)
                                             * are never below the length threshold. */
#define WELD_RECOARSEN_BELOW_VOX      5.0   /* collapse band edges shorter than this (vox);
                                             * == WELD_CLEANUP_MAX_COLLAPSE_LEN so the ~3.6-vox
                                             * rim edges are candidates (3.0 missed them: the
                                             * rim merged nothing and only bridge tris shrank) */
#define WELD_RECOARSEN_MAX_ROUNDS     8     /* collapse-round cap (1-ring locking clears ~half
                                             * a chain per round; 2.0 -> ~5 vox needs ~3) */

/* Isotropic incremental remeshing (Remesh_isotropic quality pass).
 * Botsch-Kobbelt / Surazhsky-Gotsman: split long edges, collapse short edges,
 * flip for max-min-angle, tangential relax. Interior-only (boundary frozen),
 * same target face count, fail-closed. Runs AFTER decimation to repair the
 * slivers / coarse-center-vs-fine-rim anisotropy left by collapse-only meshes. */
#define REMESH_ITERS               5      /* outer split/collapse/flip/relax iterations */
/* Split/collapse band deliberately WIDER than the textbook 4/3-4/5: the input is
 * an already-decimated mesh whose flip+smooth maintenance has converged,
 * so a tight band just churns good geometry into slivers with no net gain. Only
 * genuinely-coarse edges (>1.5 L) are split and genuine needles (<0.5 L) removed;
 * flip + tangential relax carry the rest. */
#define REMESH_SPLIT_RATIO         1.5f   /* split   edges > ratio*L (coarse only) */
#define REMESH_COLLAPSE_RATIO      0.5f   /* collapse edges < ratio*L (needles only) */
#define REMESH_RELAX_LAMBDA        0.1f    /* tangential Laplacian step             */
#define REMESH_SPLIT_MAX_ROUNDS    4       /* inner split rounds per iteration      */
#define REMESH_COLLAPSE_MAX_ROUNDS 6       /* inner collapse rounds per iteration   */
#define REMESH_FLIP_MAX_ROUNDS     8       /* inner flip rounds per iteration       */
#define REMESH_FACE_COUNT_TOL      0.10f   /* accept only if |nf_out - F| <= tol*F  */
#define REMESH_MAX_COLLAPSE_LEN    6.0f    /* never collapse an edge longer than this (vox);
                                            * < the ~7 vox inter-wrap clearance, so a collapse
                                            * can never fuse two scroll wraps (belt-and-braces
                                            * atop the interior-only + link/fold guards) */
#define WELD_CLEANUP_COLLAPSE_ROUNDS  6     /* guarded short-edge collapse rounds             */

/* Threading */
#define PARALLEL_VERTEX_CUTOFF 1000  /* don't parallelize below this */
#define PARALLEL_FACE_CUTOFF   2000

/* Halo-overlap grid stitching */
#define HALO_PIN_EPS 1.5f  /* voxels: pin verts within this distance of an
                            * integer-cube_size boundary. This is the FULL
                            * pin radius — verts further into the halo are
                            * LOP-projected and CVT-remeshable. Wider
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

/* Step 0 — patch repair (src/remesh/patch_repair.c): local tangent-plane
 * Delaunay re-triangulation of BPA "lightning-bolt" tear patches. Per
 * connected sheet only; every gate failure leaves the patch untouched.
 * Opt-in via env PATCH_REPAIR=1 while under validation. ball_pivot-style
 * local mirrors of these live in patch_repair.c. */
#define PR_TEAR_MAX_EDGES     64    /* boundary cluster larger than this is not a tear */
#define PR_PATCH_RINGS         2    /* face rings grown around a tear                   */
#define PR_PATCH_MAX_FACES   256    /* patch bigger than this is skipped                */
#define PR_PATCH_MAX_VERTS   512
#define PR_PLANAR_MAX_THICK 0.75f   /* sqrt(smallest PCA eigenvalue) cap, vox           */


/* ---- scroll axis tracking (scroll_axis_track, src/whole/axis_track.c) ------
 * Regime: PHerc0139 at 9.36 um (pitch ~9.5 vox, drift up to ~0.24 vox/vox).
 * Revalidate the windows and caps when the input scale changes. */
#define AXIS_TRACE_SLAB_VOX        64.0   /* z-slab thickness per table row */
#define AXIS_TRACE_LOCAL_R        250.0   /* samples within this radius of the running centre feed the primary fit */
#define AXIS_TRACE_TRIM_FRAC        0.60  /* IRLS keeps the best 60% by normal-line residual */
#define AXIS_TRACE_IRLS_ROUNDS      5
#define AXIS_TRACE_MIN_SAMPLES   1500     /* primary fit needs this many in-window face samples */
#define AXIS_TRACE_RESERVOIR_R    400.0   /* vox: with a prior table only samples this close to the prior curve are kept (pass 2) */
#define AXIS_TRACE_MIN_COND         0.10  /* lambda_min / lambda_max of the in-window normal scatter */
#define AXIS_TRACE_MIN_SECTORS      9     /* of 12 30-degree sectors around the centre with kept samples */
#define AXIS_TRACE_VOID_R          50.0   /* vox: the umbilicus void test disc */
#define AXIS_TRACE_VOID_MAX         0.5   /* seed candidates need inner/annulus wrap density below this */
#define AXIS_TRACE_MAX_RESID_P50   40.0   /* vox: normal-line residual median above this abstains */
#define AXIS_TRACE_MAX_JUMP       150.0   /* vox from the trend prediction: else abstain (never hold) */
#define AXIS_TRACE_FWD_BWD_AGREE   25.0   /* vox: forward/backward pass disagreement abstains both */
#define AXIS_TRACE_SEED_GRID      128.0   /* vox: spacing of the seed start grid on the bottom slab */
#define AXIS_TRACE_RIDGE_MIN_ARC  270.0   /* deg: a polar ridge must span this much theta to vote */
#define AXIS_TRACE_POLAR_DTH        2.0   /* deg per theta bin of the polar raster */
#define AXIS_TRACE_POLAR_DR         1.0   /* vox per r bin of the polar raster */
#define AXIS_TRACE_POLAR_RMAX    1200.0   /* vox: polar raster radius */
#define AXIS_TRACE_RESERVOIR_FACES 400000 /* face samples kept per slab (reservoir; uniform over the slab in pass 1, local to the prior curve in pass 2) */
#define AXIS_TRACE_RESERVOIR_VERTS 500000 /* vertex positions kept per slab (reservoir) */
#define AXIS_TRACE_RASTER_VOX       4.0   /* occupancy raster cell for the contact sheet */
#define AXIS_TABLE_MAX_SLOPE        0.5   /* |d centre / dz| above this abstains the row */
#define AXIS_TABLE_MAX_GAP_ROWS     3     /* abstained gaps up to this many slabs are interpolated */
#define AXIS_TABLE_EXTRAP_VOX     128.0   /* explicit status=extrap rows emitted past the last lock */
#define AXIS_TABLE_MANIFEST_RETIRE 20.0   /* vox: a manifest point further than this from the derived seed is retired */

/* Fixed-source multiresolution scroll model. New experimental path; the
 * existing production lane is unchanged until the rung ladder qualifies it. */
#define SCROLL_SOURCE_VERSION 2u
#define SCROLL_MODEL_VERSION 3u
#define SCROLL_MODEL_EVIDENCE_REVISION 5u
#define SCROLL_MODEL_MAX_CUBES 32768u
#define SCROLL_MODEL_REGION_MAX_VERTICES 8000000u
#define SCROLL_MODEL_MAX_FACTORS 16000000u
#define SCROLL_MODEL_PHASE_BINS 64
#define SCROLL_MODEL_MAX_PROFILE_BINS 1048576
#define SCROLL_MODEL_MRF_ROUNDS 8
#define SCROLL_MODEL_MRF_WINDOW 4

/* Evidence units are source material intervals, never output triangles. */
#define MATERIAL_EVIDENCE_MIN_UNITS 4u
#define MATERIAL_EVIDENCE_MIN_SUPPORT 3.0
#define MATERIAL_EVIDENCE_CONTAMINATION_BUDGET 1.0
#define MATERIAL_FRONT_GAP 6.0
#define MATERIAL_FRONT_TANGENT_COS 0.9
#define MATERIAL_FRONT_PHASE_TOL 0.1
/* Full PHerc0139 count: 8,076,891 intervals / 904,611,792 bytes. Retain all
 * intervals or refuse, never a prefix. This bounds records, not evidence mass. */
#define MATERIAL_FRONT_INTERVAL_BYTES UINT64_C(1073741824)
#define MATERIAL_FRONT_MAX_BLOCK_SPLITS 4000000u
#define MATERIAL_FRONT_LENGTH_SCALE 8.0
/* Neighboring original triangles are correlated evidence. A boundary interval
 * of this physical length shares at most one MRF unit, across ALL partners.
 * This is an explicit robustness scale, not a calibrated independence length. */
#define MATERIAL_FRONT_CORRELATION_LENGTH 16.0

/* Nine-array observation grid payload only; source/scratch/export are extra.
 * Keep the established 48 GB ceiling. Serialize independent projective peels
 * separately instead of raising it or changing sampling/claimant selection. */
#define RIBBON_OBSERVATION_GRID_BYTES UINT64_C(48000000000)

/* Bounded initial constrained Q1 fitter. Refuse a too-large FREE direct solve;
 * never silently decimate observations or allocate a dense UV rectangle. */
#define QUAD_FIELD_MAX_CELLS 1000000u
#define QUAD_FIELD_MAX_UNKNOWNS 500000u
/* Node records plus packed hanging stencils, NOT total pipeline memory.
 * 400 MiB is below the old 424,000,000-byte million-node array alone. Cell
 * and free-factor limits remain separate; do not mistake known data for DOFs. */
#define QUAD_FIELD_MAX_NODE_BYTES (400u * 1024u * 1024u)
#define QUAD_FIELD_NULLSPACE_PIN 1e-11
#define SCROLL_QUAD_ROOT_CELLS 256
#define SCROLL_QUAD_MIN_STEP 2.0
#define SCROLL_QUAD_FIT_ERROR 0.5

/* Supported global fit: the old claimant-separated ribbon is the observation
 * surface, not a frozen fitted solution. All Q1 interiors solve together.
 * Finest elements retain the original PL diagonal and source samples.
 * Four checkpoints: 4, 2, 1 and 0.5 voxel continuous/mesh bounds. */
#define QUAD_RIBBON_GRID_STEP 2.0
#define QUAD_RIBBON_ROOT_CELLS 32
#define QUAD_RIBBON_INITIAL_ERROR 4.0
#define QUAD_RIBBON_LEVELS 4
/* Work scheduling only: never split a connected source component or change
 * its fitting tolerance. Oversized components remain explicit and still face
 * the constructor/free-factor/storage guards. These are not LOD thresholds. */
#define QUAD_RIBBON_BATCH_FACES 1000000u
#define QUAD_RIBBON_BATCH_VERTICES 500000u
/* Independent whole-component jobs; each owns its arenas and output path.
 * Bound simultaneous copies/fields, not just each local solve. TAUCS remains
 * serialized because its legacy ordering routines have process-global state. */
#define QUAD_RIBBON_BATCH_WORKERS 4u
/* Complete 21^3 preflight: the intact 19,255,933-face primary component
 * needs 425,408,160 aligned copy bytes; every other member needs <22 MB
 * per part. At most BATCH_WORKERS copies are resident. This is only input
 * storage: independent cell, free-factor and node-storage guards remain. */
#define QUAD_RIBBON_BATCH_INPUT_BYTES (512u * 1024u * 1024u)

/* Exact raster overlap diagnostics: sparse encountered component pairs, never
 * a C*C table. Excluding isolated points, 21^3 has 47,560 primary surface
 * components (18,095,628,800 B dense) and 83,905 in extras (56,320,392,200 B).
 * Bounds both hash tables while rehashing; other raster/arena storage is extra.
 * Refuse the bake on exhaustion, never omit counts or alter pixel ownership. */
#define RAWTEX_CONFLICT_HASH_BYTES (64u * 1024u * 1024u)

/* Ink-model surface volumes (sheet_layers, 2026-09-29).  The 9 um ink family
 * (ink_9um, reader-v2) reads the centred 17 layers of an aligned 21-layer stack
 * (villa aligned21 corpus, labels at depth//2), one native voxel per layer,
 * layer index increasing INWARD toward the axis (docs/ink9um.md, W043 control:
 * outward order scored AUC 0.389, inward 0.768).  Depth zero is the mesh, i.e.
 * the F0 recto-surface midline; the model's +-2 layer training jitter and
 * reader-v2's best-quality window mean absorb a small reference offset.  The
 * grid is the review bake's (1 px per UV unit), so the ink map overlays
 * sheet_rawtex.tif; the face gates and missing-RAW reach are the bake's
 * (obj_bake_raw defaults; review --normal-range 4 + 2 trilinear pad). */
#define SHEET_LAYERS_DEPTH 21
#define SHEET_LAYERS_STEP_VOX 1.0
#define SHEET_LAYERS_OFFSET_VOX 0.0
/* --v8in: the 24-layer contract of the ink-8um v8-in models (YoussefMoNader/ink-8um-v8in and its
 * PHerc1447 fine-tune): offsets -11.5..+11.5 vox at SHEET_LAYERS_STEP_VOX, layer index inward as
 * ever (2026-09-29; the default 21-layer volume is reader-v2's). */
#define SHEET_LAYERS_V8IN_DEPTH 24
#define SHEET_LAYERS_NORMAL_SMOOTH 2      /* one-ring passes over oriented normals */
#define SHEET_LAYERS_TILE_COLS 1024       /* work tile: ~70 MB per thread at H = 3,000 */
#define SHEET_LAYERS_CHUNK 128
#define SHEET_LAYERS_BAKE_REACH_VOX 6.0
#define SHEET_LAYERS_STRETCH_RATIO 4.0
#define SHEET_LAYERS_STRETCH_FLOOR 25.0
/* Seam fill between abutting charts: vmesh_tifxyz v2 (First Letters, 2026-09-27) closes gaps of
 * at most 2 px whose two sides lie within 4.5 vox in 3D; on the PHerc0139 10^3 sheet 92% of the
 * review bake's extra pixels are such gaps.  Passes let seam crossings close. */
#define SHEET_LAYERS_SEAM_REACH_PX 2
#define SHEET_LAYERS_SEAM_TOL_VOX 4.5
#define SHEET_LAYERS_SEAM_PASSES 3
/* Context apron (sheet_layers --apron, 2026-09-29): reader-v2 on our PHerc0139 10^3 render loses ~0.11 AUC
 * within 16 px of a piece edge (vesuvius-c-f8 held-out strata); edge pixels get context from the sheet's own
 * quadric continuation up to APRON_PX out (half the 128 px model patch), fitted over APRON_FIT_PX with one fit
 * per APRON_CELL_PX cell of border pixels, sampling every APRON_STEP_PX.  The quadric keeps a radius-100 turn
 * within 0.1 vox 9 px out and 0.3 vox 17 px out (selftest). */
#define SHEET_LAYERS_APRON_PX 64
#define SHEET_LAYERS_APRON_FIT_PX 48
#define SHEET_LAYERS_APRON_CELL_PX 8
#define SHEET_LAYERS_APRON_STEP_PX 4
/* Surface regularisation (sheet_layers --smooth, 2026-09-29): on held-out PHerc0139 reader-v2 reads our
 * render better than the team's away from edges and chart steps (0.891 vs 0.884) EXCEPT where our mesh
 * wobbles: normal > 8 deg off the team's smooth surface (-0.07..-0.09 AUC) or depth jitter >= 0.4 vox over
 * 9 px (to -0.17; output/ink_sf_improve_20260929/geomdiag).  A masked Gaussian of SIGMA_PX over the sheet's
 * points and normals; pixels move along the smoothed normal only, clamped to CAP_VOX.  A radius-100 turn
 * stays within 0.12 vox at sigma 4 (selftest). */
#define SHEET_LAYERS_SMOOTH_SIGMA_PX 4.0
#define SHEET_LAYERS_SMOOTH_CAP_VOX 1.5
/* Step blending (sheet_layers --stepblend, 2026-09-29): within 32 px of a 3-D step between charts (6% of the
 * PHerc0139 10^3 yardstick, 18% of the PHerc0343_c text band) reader-v2 loses ~0.07 AUC and normal-only
 * smoothing does not recover it (the steps are tangential tears, p50 2.6 vox).  Near steps each point moves
 * toward a Gaussian-weighted local linear fit (SIGMA_PX) of its neighbours, weighted by the step strength
 * spread over SPREAD_PX; REACH_VOX pads the RAW table for the moved points. */
#define SHEET_LAYERS_STEP_SIGMA_PX 4.0
#define SHEET_LAYERS_STEP_SPREAD_PX 8.0
#define SHEET_LAYERS_STEP_REACH_VOX 4.0
/* sheet_layers --window: faces within this many UV units (px) of the requested columns are kept, so the
 * seam fill, the apron (64 px + a 48 px fit + an 8 px cell) and every tile halo see the same neighbourhood
 * as in a whole-sheet render. */
#define SHEET_LAYERS_WINDOW_PAD_PX 512.0

/* Parallelize only independent clearance projections. Candidate traversal and
 * minimum fold retain serial order. 65,536 pairs use 1 MiB of chunk arrays,
 * independent of scroll size; this is not a cap on considered pairs. */
#define INTERSECTION_CLEARANCE_CHUNK ((size_t)1 << 16)
#define INTERSECTION_CLEARANCE_PARALLEL_MIN 128u

/* Inspection textures: bounded images at native bake density, with room for
 * whole-face overlap gutters. A too-long UV edge fails, never downsamples. */
#define TEXOBJ_TILE_PIXELS 4096
#define TEXOBJ_MAX_TILE_PIXELS 8192
#define TEXOBJ_GUTTER_PIXELS 2

/* ============================================================================
 * Chart assembly (src/assemble, sheet_assemble.exe) -- the 2026-09-07 lane.
 * The units are the per-cube mesh charts kept whole; no axis, no pitch
 * constant, no integer winding.  Every value below is REGIME-BOUND (rule 26):
 * re-derive per rung and per scroll, and record the derivation.
 * ==========================================================================*/
/* Charts with less 3-D area per face than this are prediction blobs meshed at
 * voxel pitch (measured 2026-09-05: blobs 0.15-0.5 vox^2/face, papyrus 30-65,
 * fixtures 1.9): 78% of the 4x5x5 pile's vertices and 1.9% of its area. */
#define ASM_BLOB_DOUBLE_SIDED  0.5   /* BLOB = a surface that has ITSELF within 1.5-8 vox along its
                                      * own normal on half its vertices (a fused prediction meshed
                                      * as a closed ball or a thick slab).  A thin sheet never does:
                                      * the next wrap is another chart 9+ vox away.  MEASURED
                                      * 2026-09-08 on the F0 centre crop: the old rule (area/face
                                      * < 1.0) dropped 18.5% of the pile's area as "blobs" because
                                      * the CVT winding-certificate retries remesh good sheets at
                                      * 4-8x density (area/face quantized at 4.49/2.23/1.11/0.55),
                                      * so the sheet 'died' where the prediction was solid. */
#define ASM_BLOB_AREA_PER_FACE 0.2   /* Legacy config/cache compatibility only. A face-area
                                      * threshold cannot identify physical material. */
/* Confetti below this 3-D area (vox^2) never enters the layout. */
#define ASM_MIN_CHART_AREA 100.0
/* A chart whose flattening has a singular value outside [1-band, 1+band]
 * on any face is SUSPECT (lower prior weight, first cut candidate). 0.05 is
 * the unroller_v2 acceptance band; papyrus is developable to ~1%. */
#define ASM_STRESS_BAND 0.05
#define ASM_FLATTEN_ITERS 30
/* Free-boundary conformal (LSCM) seed offered to a chart whose Tutte-seeded
 * map stays outside the stress band (2026-09-25). A crushed scroll's layers
 * turn back inside single cubes; those charts are developable, but a circular
 * Tutte seed compresses them so far that guarded ARAP stops strained, and the
 * relate gates then refuse physically continuous turn seams as cross-wrap
 * (PHerc1447 eye slab: 26 of 32 refused turn seams pass on near-isometric
 * maps). Accepted only through the certified recovery gate; charts inside the
 * band are bitwise unchanged. 0 disables. */
#define ASM_FLATTEN_CONFORMAL_RECOVERY 1
/* Local original-metric (band) recovery budget, 2026-09-28. Over 3,054
 * recoveries in 20 First Letters logs the success rate falls with the number
 * of strict failures it starts from: 95.6% at 1, 52% at 5-16, 20.8% at
 * 65-256, 7.4% at 257-1,024 and 2.0% above, while every failure runs the full
 * 1,200 L-BFGS steps. Successes need a p90 of 24-466 steps. On PHerc0343_c
 * crumpled-chart rescues spent 998 of 1,800 s here (196 charts ran to
 * the 60 s cap with nothing certified). A recovery starting from more than
 * MAX_BAD strict failures is not attempted; the rest stop after MAX_STEPS
 * (keeps 97% of the measured successes, saves ~62% of the steps). */
#define ASM_FLATTEN_BAND_RECOVER_MAX_BAD 256
#define ASM_FLATTEN_BAND_RECOVER_MAX_STEPS 400
/* Seam correspondences across the cube gap: trim inset 1 vox per side plus
 * BPA erosion ~1.2 gives ~4.4-vox gaps; 3.5 vox is the in-tree pair gate
 * (group_graph.h) and stays under half the 7-vox inter-wrap clearance. */
#define ASM_SEAM_MAX_VOX 3.5
/* Continuation runs follow original mesh boundary arcs on both sides. These
 * are the existing run budgets, measured in source length instead of UV/PCA
 * extent. Short observations remain obligations once a full run is present. */
#define ASM_CONT_RUN_MIN_SAMPLES 6u
#define ASM_CONT_RUN_MIN_LENGTH 8.0
#define ASM_CONT_RUN_MAX_STEP 16.0
#define ASM_CONT_RUN_MAX_LENGTH_ERROR 4.0
#define ASM_SEAM_MIN_LEN 30.0
#define ASM_SEAM_NORMAL_DOT 0.8
/* Tangency gate: a pair whose displacement runs more than this along the
 * normal (per-pair median) is a face-to-face CONTACT of two wraps. */
#define ASM_TANGENCY_MAX_VOX 1.75
/* Tangency rejections are ledgered as ASM_REL_CONTACT relations without
 * correspondences (2026-09-17: the flag was documented as that ledger but was
 * never set in production, so the layout-join veto ASM_JOIN_VETO_FLAGS never
 * fired).  0 = the pair leaves no relation behind. */
#define ASM_RELATE_LEDGER_CONTACT 1
/* THE AXIS DERIVED FROM THE MATERIAL (asm_axis_derive.c, stage `axis`, 2026-09-17).  Face
 * samples of the cleaned charts, at most this many per chart (area-weighted), projected
 * onto the current polyline and binned by arc length into slabs in the LOCAL tangent frame;
 * per slab the point closest to the in-plane normal lines by Cauchy IRLS, the scale
 * annealed from the first value to the second over the rounds; a slab locks with at least
 * the sample count, a conditioning (lambda_min / lambda_max of the weighted normal scatter,
 * ~35 degrees of aperture at 0.03) and a residual median under the limit; the locked
 * centres go through the LOESS of a loaded table and the whole is iterated.  Without a
 * prior the seed is found on z-slabs from a grid of starts: the void-qualified attractor
 * with the lowest residual (the tracker's doctrine: never rank by support). */
#define ASM_AXIS_DERIVE_SAMPLES_PER_CHART 256
#define ASM_AXIS_DERIVE_SLAB_VOX 64.0
#define ASM_AXIS_DERIVE_ITERS 3
#define ASM_AXIS_DERIVE_MAX_AXIAL 0.7        /* |n . t| above this: a face carries no centre information */
#define ASM_AXIS_DERIVE_LOCAL_R 250.0        /* the window about the running centre that feeds a slab's fit (the first two rounds use twice it);
                                              * measured 2026-09-17: without it the whole-slab fit lands 39 vox rms from the 4x5x5's trusted table
                                              * (b 1.02 -> 0.90) and 130-360 vox off on the big rungs -- the evolute cusps of flattened wraps */
#define ASM_AXIS_DERIVE_WINDOW_ROUNDS 8
#define ASM_AXIS_DERIVE_CAUCHY_VOX0 48.0
#define ASM_AXIS_DERIVE_CAUCHY_VOX 12.0
#define ASM_AXIS_DERIVE_IRLS_ROUNDS 6
#define ASM_AXIS_DERIVE_MIN_SAMPLES 500
#define ASM_AXIS_DERIVE_MIN_COND 0.03
#define ASM_AXIS_DERIVE_MAX_RESID 80.0       /* the tracker abstains at 40; pherc343's sections read 47-69 on every z-slab (2026-09-17), so the
                                              * bar admits them and the pose / placement judge the frame (A/B: 40 locks nothing there) */
#define ASM_AXIS_DERIVE_SEED_SLAB 128.0
#define ASM_AXIS_DERIVE_SEED_GRID 128.0
#define ASM_AXIS_DERIVE_SEED_SAMPLES 20000
#define ASM_AXIS_DERIVE_VOID_R 50.0
#define ASM_AXIS_DERIVE_VOID_MAX 0.5
/* A derived curve replaces the configured frame only when at least this fraction of its slabs
 * locked and its rows sit within this rms of the smoothed curve; otherwise it is reported and the
 * configured frame stands (2026-09-17: the 21x5x5, a column beside its umbilicus, locks 20 of
 * 51 slabs at 182 vox rms -- evolute attractors inside the column, not an axis). */
#define ASM_AXIS_DERIVE_ACCEPT_LOCK 0.5
#define ASM_AXIS_DERIVE_ACCEPT_RMS 48.0      /* the 4x5x5's own tracker table sits 14.6 vox rms off its curve, the CT prior 21.5 */

/* Unimodality a seam must reach to be the FIRST join between two otherwise
 * disconnected blocks (0 = off).  Such a seam is uncorroborated and
 * load-bearing: on the familiar 4x5x5, 14 relations in the bottom tail
 * (unimodality 0.80-0.99 against an accepted median of 0.991) chained old
 * components of 239+59+14+12+11 charts into one 358-chart primary, whose
 * layout then blocked the correct seam placement of a 72-chart, 754k vox^2
 * component.  A relation that only corroborates an existing connection keeps
 * the ordinary ASM_SEAM_UNIMODAL bar. */
#define ASM_SEAM_BRIDGE_UNIMODAL 0.95
/* The premium's pass unions in descending unimodality and asks "is this a
 * bridge SO FAR", so two independent marginal seams joining the same pair of
 * blocks demote each other: the first is uncorroborated when it is tried, the
 * second because the first was demoted.  With this on, a seam is only
 * load-bearing when it is a CUT EDGE of the whole gate-accepted graph.
 * Replayed from the 09-17 ledgers (analysis only): 384 of the 21x5x5's 969
 * demotions and 421 of pherc343's 1,888 are corroborated by another seam, and
 * the 4x5x5's 67 include 19 -- components after relate 2,325 -> 2,136,
 * 6,313 -> 6,085 and 131 -> 122.  0 until the A/B through placement. */
#define ASM_SEAM_BRIDGE_CORROBORATED 0
/* A block below this area (vox^2) is a leaf, not a weld: a seam attaching it
 * keeps the ordinary ASM_SEAM_UNIMODAL bar instead of the bridge premium
 * (0 = off).  The premium guards joins between components; at the 21x21x21
 * it also strands single charts whose only seam is marginal, and region
 * workers grow along required-source seams only, so such a chart is never a
 * candidate.  The 2026-09-24 masked-ends run left 100.5 M vox^2 (11.3% of the
 * source) in required-source components below the 20,000 vox^2 region minimum;
 * 5,423 of those charts (30.4 M) hang off placed material by seams relate had
 * ACCEPTED (median unimodality 0.897, rms 0.91) and the premium demoted: 23.3 M
 * of it true leaves (one such seam), 8.7 M corroborated by a second seam
 * (which ASM_SEAM_BRIDGE_CORROBORATED would restore).  The value is the
 * repair's region minimum.  A/B pending (needs a new inventory). */
#define ASM_SEAM_BRIDGE_LEAF_AREA 0.0
/* MEASURED DEAD 2026-09-17, and the measurement is the reason the 21x5x5 and
 * pherc343 sheets looked the way they did.  The idea (09-16, never committed
 * but left in the working tree, so every big-rung number of that day carries
 * it): let the unimodality statistic reject a seam only when its RAW rigid fit
 * ALSO misses this rms, on the argument that correspondences fitting one rigid
 * transform within the audit's seam bar cannot contain a sheet switch.  On the
 * 4x5x5 that promotes 44 of 81 unimodality rejections back to joins.  It is
 * wrong: the statistic is the cross-wrap detector, and a fused pair of wraps at
 * the crumpled core fits a single rigid transform very well.  Three piles, one
 * exe apart, same caches:
 *   4x5x5     gates 6/8 -> 5/8, cross-wrap 0 -> 2, long edges 7 -> 18,
 *             overlaps 11,339 -> 12,185, area represented 97.44% -> 97.21%
 *   21x5x5    coverage 69.6% -> 52.0%, overlaps 724,508 -> 1,541,975, seams
 *             closed 4,345 -> 3,208, stacked (u,v) 8,158 -> 15,140; and before
 *             placement it merges 52% of the area into one component whose
 *             layout carries 41,403 contradicting chart pairs against 1,259
 *   pherc343  coverage 51.0% -> 15.5%, area placed 49.1% -> 20.7%
 * The constant stays as the record; nothing reads it. */
#define ASM_SEAM_UNIMODAL_RMS 2.0
/* The unimodality gate's verdict is recorded as ASM_REL_CROSSWRAP and measured
 * for EVERY pair that gets a rigid fit (2026-09-17).  Two blind spots made it
 * partial: an isometry rejection set WEAK and jumped straight past the
 * unimodality block, so the cross-wrap detector never ran on 1,971 of
 * pherc343's pairs and 356 of the 21x5x5's; and a seam the detector DID reject
 * was handed to placement as evidence anyway (pherc343 3,041 of them, 21x5x5
 * 1,220, the 4x5x5 68) because WEAK means both "rejected" and "placement
 * anchor".  The statistic costs ~200 random correspondence pairs, so measuring
 * it always is free.  A SHORT seam (12-30 vox) still cannot be measured: the
 * test needs 8 sampled spans of 8+ vox. */
#define ASM_WRAP_MEASURE_WEAK 1
/* Negative evidence: a seam carrying ASM_REL_CROSSWRAP is refused as a placement anchor, as a
 * layout join and as a discovery hypothesis, instead of being demoted to placement evidence.
 *
 * MEASURED 2026-09-18, all three together, one exe apart, both arms sharing the machine and
 * the 4,725-s budget.  It is the 21x5x5's first 3/8: FRAGMENTATION/coherence 49.7% FAIL ->
 * 82.9% PASS, coverage 72.7% -> 74.4%, stacked (u,v) 8,258 -> 6,414, (u,z) 15,250 -> 11,039,
 * post-placement contradiction 5.19e4 -> 4.51e4.  No gate regressed.  But the 4x5x5, which
 * carries only 70 such relations, places 72 of its 153 components instead of 105 on 8,475
 * anchors instead of 12,045: on a pile whose axis frame is armed those 70 seams are
 * load-bearing evidence the audits were already policing.
 *
 * So the master ships 0 and the veto is split three ways, because the three do different work:
 * the anchor veto is what cleans the 21x5x5's hop poses, the join veto is what costs it
 * lineages (39 layout joins contact-vetoed, 588 -> 738), and the discovery veto is untested
 * alone.  Set ONE of them to measure it. */
#define ASM_WRAP_NEGATIVE 0
/* apl_placement_seam: 0 = still an anchor, 1 = never an anchor, 2 = never an anchor on a pile
 * whose AXIS FRAME IS DISARMED.  2 is the reading the A/B points at: negative evidence should
 * bind where there is nothing better to check a pose against.  The 4x5x5 has 73,971 anchors
 * from the radius, a v pin and a turn model measured over 568 real turns, so its audits police
 * a bad cross-wrap pose and its 70 flagged seams are worth keeping; the 21x5x5 and pherc343
 * have no frame at all, and there the flagged seam simply wins. */
#define ASM_WRAP_VETO_ANCHOR   ASM_WRAP_NEGATIVE
#define ASM_WRAP_VETO_JOIN     ASM_WRAP_NEGATIVE   /* ASM_JOIN_VETO_FLAGS: never a layout join */
#define ASM_WRAP_VETO_DISCOVER ASM_WRAP_NEGATIVE   /* ad_usable: never a discovery hypothesis */
/* Isometry check: the MEDIAN relative discrepancy between the two charts'
 * seam parameterizations (|d_a - d_b| / d_a over random correspondence
 * pairs).  Two independent flattenings differ by up to twice the stress band
 * on a seam; a sheet switch inside the seam is a step of a wrap spacing. */
#define ASM_ISO_TOL 0.12
/* Position noise of one boundary correspondence, vox: the two charts' boundary vertices are
 * the mesher's own samples on either side of the 2-vox cube gap, so a pair is only known to
 * about a vertex spacing.  The isometry and unimodality tests compare pairs of
 * correspondences 8+ vox apart; without this budget a 1-vox error over 8 vox is a 12%
 * length discrepancy and a 0.125 rad angle, and both gates rejected genuine seams by the
 * hundreds on the centre-region piles (measured 2026-09-08: unimodal rejections 17% of
 * chart pairs on the 21x5x5 against 7% on the familiar, isometry 11% on the coarse
 * no-winding pile against 3%; the rejected seams fit rigidly at rms 1.9). */
/* MEASURED 2026-09-08: at 1.0 the budget raises connectivity everywhere (familiar 929 -> 965
 * joins, 21x5x5 largest component 510 -> 858 charts) and 23 of 24 checkable new joins are
 * genuine, but the merged components then carry the rare wrong join that no co-location audit
 * can see (the two wraps do not overlap in the component's own frame) and placement refuses
 * them whole: familiar 88.5% -> 50.8% placed.  OFF until placement can cut a component at the
 * join between its disagreeing anchor clusters; `relate.seam_noise_vox` in the config for A/Bs. */
#define ASM_SEAM_NOISE_VOX 0.0
/* clean2's single-drop pass (every flagged relation tried alone, kept when the cleaning score
 * rises).  MEASURED DEAD 2026-09-08: the score has no connectivity term, so a drop that splits a
 * component "resolves" its contradiction for free; the pass dropped 22-28 relations per pile,
 * the familiar's primary shrank 2.55M -> 0.76M vox^2 and placement fell 88.5% -> 43.4%.  Kept
 * in-tree, off; a split penalty in the score is the precondition for turning it on. */
#define ASM_CLEAN_SINGLE_DROPS 0
/* finish clean2 on its best-scoring audited state (else on the last).  OFF: the best-scoring
 * state is not the best for placement (the familiar's round-6 state scored 7.92e6 with a
 * non-sheet-like big component; the end of the extras pass scored 8.07e6 and placed 88.5%). */
#define ASM_CLEAN_KEEP_BEST 0
/* Layer-neighbour probe reach along the normal (vox) before a local layer
 * distance has been measured. */
#define ASM_LAYER_PROBE_VOX 25.0
/* LAYER CUT (asm_layer_cut, between relate and pose): a join path that connects two charts
 * relate found stacked one or more wraps apart WITHOUT winding around the axis crosses a
 * bridge between wraps; each such path nominates its most used join, and the nominated joins
 * are dropped until no witness path is left.  Measured 2026-09-22 on the source inventories:
 * the 10x10x10's largest source region holds 1.0 radial layers per (unwrapped angle, z) bin
 * and its 1,076 closing cycles all wind 0; the 21x21x21's region of 40,436 charts (60% of the
 * area) spans 3.1 turns with 20.5 layers per bin, and ~39k of its layer pairs are joined within
 * 16 joins where one turn at r 600-1500 is 30-80 joins.  Charts nearer the axis than
 * MIN_RADIUS have no usable angle and are neither traversed nor cut.  1 = on. */
#define ASM_LAYER_CUT 1
#define ASM_LAYER_CUT_MAX_HOPS 16
#define ASM_LAYER_CUT_MIN_RADIUS 128.0
#define ASM_LAYER_CUT_MIN_HITS 3
#define ASM_LAYER_CUT_MAX_PASSES 512
/* Contradiction raster cell (vox) and the continuous contradiction weight
 * smoothstep bounds as fractions of the LOCAL layer distance. */
#define ASM_CELL_VOX 8.0
#define ASM_CONTRA_LO 0.25
#define ASM_CONTRA_HI 0.60
#define ASM_CLEAN_ROUNDS 8
/* clean2 continues past ASM_CLEAN_ROUNDS while its last round was KEPT and gained at least
 * EXTEND_GAIN of the score (measured progress), up to MAX_ROUNDS and EXTEND_SECONDS since
 * clean2 began.  2026-09-22: the 10x10x10's 8th round gained 0.06% (the rule leaves it at 8
 * rounds); the layer-cut 21x21x21's 8th gained 1.5% with contradicted pairs still falling
 * (28,383 -> 17,768 over 8 rounds).  0 gain = the fixed round count. */
#define ASM_CLEAN_EXTEND_GAIN 1.0e-3
#define ASM_CLEAN_MAX_ROUNDS 48
#define ASM_CLEAN_EXTEND_SECONDS 900.0
/* SPLIT GUARD: a cleaning batch that separates material from its component (a dropped bridge
 * join, or a switch the drop triggers) must gain at least this fraction of the separated area in
 * score.  Splits were free (the audit compares claimants within one component, so parting two
 * blocks removes every contradiction between them while both still count as covered): on the
 * familiar the straightened poses of 2026-09-09 let clean2 part the 2.57M primary into 1.11M +
 * 0.70M + 0.47M vox^2 for gains of 0.8k-52k, and the 0.47M part then failed to place.  A refused
 * bridge is protected for the rest of the loop.  0 = off (the pre-2026-09-09 behaviour). */
#define ASM_CLEAN_SPLIT_MIN_GAIN 0.20   /* familiar A/B 2026-09-09: 0.02 -> gates 4/8, primary 1.87M, depth seams 0.101%; 0.20 -> 5/8, 2.16M, 0.056% (+273 stacked (u,z) cells, spine p90 8 -> 12) */
/* Relations whose robust weight falls below this after the pose solve are
 * switched off. */
#define ASM_SWITCH_OFF 0.2
/* Placement of disconnected components by measured turn vectors. */
#define ASM_PLACE_MIN_ANCHORS 12
#define ASM_PLACE_MAX_RMS 12.0
/* THE TURN LENGTH'S TWO MEASURED DEFECTS (2026-09-19, diagnosed offline on the 10x10x10's
 * stage5_anchors.csv against the 4x5x5 champion's).  Both estimators are in the same family and
 * both err in proportion to how far the pile sits from 2*pi*d, which is exactly the
 * champion/10x asymmetry -- the champion's layer coordinate IS a wrap index (dr/dlidx 16.3
 * against d 15.0), the 10x10x10's is not (7-35 against d 13).
 *
 * D1, GROWTH IS IN THE WRONG UNITS.  apl_measure_growth returns the median of |T_out| - |T_in|
 * per LIDX STEP -- 196.3 vox on the 10x10x10, 129.1 on the champion -- but apl_turn_for and
 * apl_turn_k consume it per PHYSICAL WRAP (lp->k).  Regressing the measured one-wrap turn on
 * the axis radius gives dT/dr = 0.96 x 2*pi (10x) and 0.95 x 2*pi (champion), so the offset-
 * curve rule 2*pi*d = 81.7 / 94.2 is the right number: the shipped value is 2.40x / 1.37x too
 * large, on 53.2% / 57.2% of all anchors.
 * 0 = the measured median (shipped).  1 = 2*pi*d, which is already the n < 3 fallback inside
 * apl_measure_growth and the loop's own initial value, so 1 is a REMOVAL, not new physics.
 * 2 = the fitted slope of T on r carried to d, clamped to [0.8, 1.25] x 2*pi*d. */
#define ASM_PLACE_GROWTH_PER_WRAP 0
/* D2, THE RADIUS TURN HAS NO INTERCEPT.  apl_turn_radius returns 2*pi*k*r -- a circular
 * section.  The measured turns fit T(r) = 2*pi*r + c with c = 498..692 vox on the 10x10x10 and
 * 361..536 on the champion: the section's perimeter excess, 1.6-2.7x the hop gate.  Dropping it
 * makes the radius family the WORST-FITTING anchors on both piles (kept |residual| p50 60.1 vox
 * against 19.9 for measured on the 10x; 30.1 against 5.3 on the champion), and of the candidates
 * holding both families T_measured/T_radius reads p50 1.27 / p90 2.09 with the two landing in
 * disjoint clusters 32% of the time.  A candidate one intercept short of its wrap lands on the
 * previous wrap, which is what 2,017 of 2,071 `contradiction mass` audit rejections are.
 * 0 = shipped (no intercept).  1 = add k * c from apl_fit_turn_radius's robust fit. */
#define ASM_PLACE_TURN_INTERCEPT 0
#define ASM_PLACE_TURN_C_MIN_N 24      /* charts needed before an intercept is adopted at all */
/* HOW WELL THE TURN LENGTH IS KNOWN.  When no axis table is armed, the distance in u
 * between two adjacent wraps comes from the layer stretch: T = 2*pi*d/(s-1), a ratio of
 * small differences.  On the 21x5x5 `s` sits within 4% of 1, so a 1% error in `s` is a 25%
 * error in T, and the per-component turns the estimator emits run 645, 1,290, 1,361, …,
 * 3,791, 8,506, 9,105 vox where the truth varies maybe 1.5x (2026-09-17).  The slope's
 * standard error is free from the same accumulators, inflated by the design effect
 * sqrt(hits/pairs) because one layer pair's hits are one look at the geometry, and
 * se(s-1) = se(s).  0 keeps the shipped fixed bar (|s-1| < 0.02 falls back to the curvature
 * turn); positive replaces it with this bound on se(s-1)/|s-1|, which is the same statement
 * about T's own relative error.  Reported either way in the place log and
 * stage5_components.csv.  A/B pending. */
#define ASM_STRETCH_MAX_REL_ERR 0.0
/* Layout-confirmed joins (after placement).  Two charts of different lineages
 * whose adjacent or co-located audit cells AGREE in 3-D (the same layer:
 * normal-direction separation under CONTRA_LO of the local layer distance;
 * normals aligned; the 3-D step under the cell step plus the wrap gate) over
 * at least ASM_JOIN_MIN_CELLS cells, with at most ASM_JOIN_MAX_DISAGREE of
 * their contact cells contradicting, are JOINED by a relation taken from
 * their placed poses (ASM_REL_LAYOUT).  Geometric continuity the layout shows
 * thereby becomes a decision of the stitcher, ledgered; the verdict lattice
 * joins nothing else (reviewer, 2026-09-09: raster coincidences were joining
 * lineages).  4 cells x 8 vox = the seam minimum length. */
#define ASM_JOIN_MIN_CELLS 4
/* A contact with NO contradicting cell needs only this many agreeing cells: on the
 * familiar 4x5x5 (dev48) 14 of the 23 "short" refusals had 2-3 agreeing cells, none
 * disagreeing, normal separation 0.00-0.06 vox (A/B 2026-09-09; 4 = off). */
#define ASM_JOIN_MIN_CELLS_CLEAN 2
#define ASM_JOIN_MAX_DISAGREE 0.10
#define ASM_JOIN_NORMAL_DOT 0.8
#define ASM_JOIN_STEP_SLACK 6.0
/* A claim pair AGREES up to this many slacks of in-plane step beyond the cell step (the
 * normal direction and the normals decide the layer; the plane offset is what the
 * alignment below corrects), and the alignment shift is capped at the same reach.  On the
 * familiar (dev50) pair 141-144 agreed on 8 cells and "contradicted" on 6 only because it
 * sat 12 vox off in v (A/B 2026-09-09; 1.0 = the old rule). */
#define ASM_JOIN_AGREE_REACH 2.0
/* A pair that already has a WEAK seam relation (rigidly fitting, gate-rejected) needs only
 * this many clean agreeing cells. */
#define ASM_JOIN_MIN_CELLS_WEAK 1
/* Contact evidence is also taken across a gap of up to this many empty audit cells (8 vox
 * each) to the right and below: the expected 3-D step grows with the gap, so a block placed
 * short of its neighbour measures its own offset and is pulled in by the alignment. */
#define ASM_JOIN_GAP_CELLS 3   /* 8 (64 vox) measured 2026-09-09: 52 joins instead of 36 on the familiar but edges over the gate 846 -> 1,095 and no rotation recovered */
/* The verdict lattice closes a gap of up to this many cells (du = 2 vox each) between two
 * cells of charts of ONE LINEAGE (joined directly or through other charts), in u and in v,
 * with generated cells (support 0) when the 3-D step is a sheet's: a measured join whose
 * charts stop short of touching at lattice resolution is one component, as the evidence
 * says; a gap inside one chart stays a hole.  0 = off. */
#define ASM_VERDICT_BRIDGE_CELLS 3
/* POST-PLACEMENT SEAM RE-SOLVE (2026-09-10).  A piece placed by the turn model lands tens of vox off
 * its seam neighbours, and nothing rejoined it (the layout join reaches 24 vox and aligns 12); its
 * seams are still in the run: accepted seams the cleaning dropped or the solve switched off, and the
 * WEAK / SHORT placement-only seams.  After the layout joins, every such seam between placed charts
 * whose pair now sits within ASM_READMIT_VOX and ASM_READMIT_RAD of the pose it implies is READMITTED
 * as a join, and the placed charts outside the largest lineage are re-solved over the joins with
 * that lineage fixed (AsmPose_refine).  Kept only if the contradiction after does not exceed the
 * contradiction before by more than ASM_REFINE_MAX_CONTRA_GROWTH plus ASM_REFINE_CONTRA_FLOOR vox^2
 * (the 10x's contra after placement is ~8e4 vox^2: 10% alone is 128 cells); round 2 un-readmits the
 * seams whose charts' contradiction grew; else the whole step is reverted.  MEASURED 2026-09-10 (10x10x10):
 * with the deferral on, lineages 133 -> 116 (dev84) / 122 (dev86, rotation held, 50-vox gap) while long
 * lattice edges rose 14.1k -> 22.2k and cross-wrap edges 558 -> 877 (the readmitted seams carry 42 long
 * edges per join against 24 for accepted seams); without the deferral lineages 543 -> 515 (dev85).  The
 * layout joins already merge what the deferral leaves; OFF by the champion rule, kept in-tree. */
#define ASM_REFINE_SEAMS 0
/* AsmContinuity_confirm promotes a gate-rejected relation to a join when its
 * charts share a component.  That promotion must clear the same seam bar the
 * audit applies (AsmField_seam), or the isometry and unimodality
 * gates are silently undone downstream -- and it bypassed ASM_REFINE_SEAMS,
 * which the champion rule had already switched off for exactly this effect. */
#define ASM_CONTINUITY_CONFIRM_RMS 2.0

/* The continuity measure's pass gate, in vox.  It used to be fmax(8, 4*rms)
 * with rms the relation's own relate-stage residual -- a tolerance set by the
 * error it was supposed to judge.  While relate's rms was laundered (~0.2) the
 * 4*rms term never bound and every measurement in the record was made against
 * a flat 8; once rms was honest again the term handed the WORST seams the
 * widest gates: on the familiar 4x5x5 five seams passed at 8.6-27.7 vox
 * through gates of 35-67, and 25 of 1,035 passing seams sat above the audit's
 * own 2-vox seam bar.  A flat 8 is the gate the accepted results were actually
 * measured under. */
#define ASM_CONTINUITY_GATE_VOX 8.0
/* SCALE multiplies the relation's own relate residual into the gate: 4.0 is the
 * self-referential law measured above, 0 the flat floor.  Both stay in-tree for
 * the A/B; ship what the sheet measures better. */
#define ASM_CONTINUITY_GATE_SCALE 0.0
#define ASM_CONTINUITY_GATE(rms) fmax(ASM_CONTINUITY_GATE_VOX, ASM_CONTINUITY_GATE_SCALE*(rms))

/* Restoration QP solver.  1 = Goldfarb-Idnani dual active set on one TAUCS
 * factor of H: the exact optimum of the linearised model, cost scaling with
 * the BINDING rows (src/flatten/active_set_qp.c).  0 = OSQP ADMM on the full
 * (n+m) KKT, set up and factored on every call, up to 4000 iterations, its
 * iteration-budget iterate accepted as the step -- the 2026-09-15 profile
 * showed QDLDL triangular solves at 0.98 cores for minutes per iteration on
 * 480k-coordinate levels.  A refused active-set solve falls back to OSQP. */
#define ASM_RESTORATION_ACTIVE_SET 1
/* Source-repair trials: 1 = the queued seams whose causal core is a large
 * component are restored TOGETHER, one ladder per component, and each seam
 * is then judged on its own; a batch that regresses a passing seam or moves
 * fixed context is reverted and its seams take the sequential path.  0 = one
 * whole-context ladder per seam (2026-09-16: 11 of the first 23 trials on
 * the 4x5x5 each restored the entire 973k-face primary, ~1 h apiece at three
 * closure iterations, ~4 h at twelve). */
#define ASM_SOURCE_TRIALS_BATCHED 1
#define ASM_SOURCE_TRIALS_BATCH_MIN_CHARTS 32
/* ...and a MAXIMUM on the batch's moving charts, which nothing bounded before 2026-09-19.  The
 * moving set is the union of the members' causal cores, so the more coherent the placement the
 * larger the batch: the 21x5x5, whose components are rubble, batches 70-150 charts, while the
 * 10x10x10, whose primary holds 63% of its sheet, produced a FIRST batch of 1,975 moving charts
 * and sat in one solve for over an hour.  Better placement must not mean a slower repair.  The
 * overflow simply forms the next batch, and a lone seam whose own core exceeds the cap still
 * runs by itself, so no seam becomes unreachable.  0 = no cap (the old behaviour). */
#define ASM_SOURCE_TRIALS_BATCH_MAX_CHARTS 512
#define ASM_SOURCE_TRIALS_BATCH_ROUNDS 2
/* A queued seam whose initial RMS is beyond this many vox is a placement
 * error, not a source-repair target: on the 4x5x5 (2026-09-16, 35 trials)
 * every trial that started above 59 vox failed (59, 86, 1,301, 2,035,
 * 2,037, 2,037, 2,757, 2,768) and the 1,301-vox one held the whole primary
 * for three hours; every accepted trial started below 18. */
#define ASM_SOURCE_TRIAL_MAX_RMS 32.0
/* A trial moves only the charts within this many relation hops of the
 * seam's two endpoints (0 = the whole causal core).  The causal core of a
 * seam in the primary is the primary's whole cyclic block: 278 moving
 * charts / 1.3M coordinates on the 4x5x5, and the profile of 2026-09-16
 * (36 stack samples, one busy thread) put 28% of the repair in TAUCS's
 * symbolic analysis of that matrix and 36% in the constrained step over
 * it.  A seam within reach (RMS <= ASM_SOURCE_TRIAL_MAX_RMS) closes by
 * deforming its neighbourhood inside the metric bands; the rest of the
 * sheet is collision context.  Measured A/B: see the changelog. */
#define ASM_SOURCE_TRIAL_CONTEXT_HOPS 0   /* A/B 2026-09-16 (same checkpoint, closure 3): hops 2 closed 1,042 of 1,106 seams with 14,555 overlaps, hops 0 closed 1,046 with 12,185; the audit decides, the time comes from the bounded admission */
/* The sequential path retries a failed trial with the complete source
 * component moving; only for components up to this many charts (the
 * whole primary retried is the hours the hop cap removes). */
#define ASM_SOURCE_TRIAL_EXPAND_MAX_CHARTS 64
/* apr_solve repartitions the contact tree only when at least this fraction
 * of the field's faces can move.  Every contact measurement refits the
 * boxes, so a small trial (13 charts of 700 on the 4x5x5) that skips the
 * rebuild measures the same contacts; the rebuild was 19% of the hop-capped
 * trial phase (90 samples, 2026-09-16). */
#define ASM_REPAIR_REBUILD_MIN_MOVING 0.25
/* The clearance QP (chart translations against contact slacks, OSQP) is
 * bounded: one solve may take at most this many seconds and one clearance
 * step at most the budget over its cutting-plane rounds.  2026-09-16: a
 * final-admission candidate of five charts on top of 24,738 contact pairs
 * (another wrap) held the 4x5x5 repair in ADMM for over 20 minutes; a
 * solve that cannot finish in this time is refused and the step is not
 * taken, which is what a hopeless clearance deserves. */
#define ASM_CLEARANCE_QP_TIME_LIMIT 30.0
#define ASM_CLEARANCE_QP_BUDGET 90.0
/* Compact work is bounded by elapsed time and local workspace bytes. These
 * counts now protect signed indices only, not the density of original mesh. */
#define ASM_REPAIR_PATCH_FACES 2147483647u
#define ASM_REPAIR_PATCH_FIELD_FACES 2147483647u
#define ASM_REPAIR_PATCH_COORDINATES 2147483647u
/* Per-transaction storage, excluding the reusable immutable whole sheet.
 * Includes conservative field/assembly preflight and sampled process growth. */
#define ASM_REPAIR_PATCH_BYTES ((size_t)2*1024*1024*1024)
/* Dimensionless singular-value comparison allowance. Packing subtracts large
 * encoded coordinates; a few last-place bits must not drop a whole face from
 * an area percentile. This is far below source float precision. */
#define ASM_METRIC_ROUNDOFF 1e-8
/* Fixed boundary-alignment policy. Both RMS and supported sample distances
 * use this bound in admission, repair constraints and the complete audit. */
#define ASM_SEAM_TOLERANCE_VOX 4.0
/* Automatic repair of missing neighborhoods in represented material.
 * Reserve the remainder of normal repair for present seam/contact defects. */
#define ASM_REPAIR_NEIGHBORHOOD_WORK_FRACTION .85
#define ASM_REGION_REFINE_ITERATIONS 64
#define ASM_REGION_REFINE_CONTACT_THREADS 4
/* The refinement's interior metric reserve against rigid packing round-off.
 * 0 since ASM_REPAIR_STRICT_BAND keeps changed faces inside the nominal band
 * (2026-09-23: the 1e-8 reserve refused candidates that passed every other
 * certificate -- 554 refinement targets ended on packing_metric_margin in the
 * 5 h run -- while the strict band holds the audit with 0 strict failures). */
#define ASM_REGION_REFINE_METRIC_MARGIN 0.0
#define ASM_REGION_REFINE_ATTEMPT_SEC 8.0
#define ASM_REGION_REFINE_FIRST_SEC 4.0

/* Read-only physical navigation sampling; no geometry admission policy. */
#define ASM_READING_SAMPLE_AREA 128.0
#define ASM_READING_MIN_SAMPLES 8u
#define ASM_READING_MAX_SAMPLES 256u
#define ASM_READING_AXIAL_BIN 8.0
#define ASM_READING_PHASE_BINS 1024
#define ASM_READING_CORE_RADIUS 64.0
#define ASM_READING_PHASE_RESIDUAL 0.20
#define ASM_READING_SENSE_MIN_TURNS 1.0
#define ASM_READING_SENSE_MIN_PITCH 0.25
#define ASM_READING_SENSE_MAX_PITCH 4.0
/* Presentation spacing only. Does not authorize a material continuation. */
#define ASM_READING_RIBBON_GAP 32.0
/* Give a contact-only residual an early opportunity without consuming the
 * fitting fallback's whole quantum. These bound work, not material size. */
#define ASM_REPAIR_PRE_FIT_CLEARANCE_SEC 2.0
#define ASM_REPAIR_PRE_FIT_CLEARANCE_STEPS 8
#define ASM_REPAIR_PRE_FIT_CLEARANCE_SHARE 0.5
#define ASM_REGION_REFINE_ESTIMATE_FLOOR_SEC 1.0
#define ASM_REGION_REFINE_GROWTH_SEC 5.0
/* A refinement target stops at the repair patch budget (45 s) unless its
 * last joint quantum still cut the measured violation below ..._RATIO of its
 * starting value; such a converging target may continue to ..._SEC in total.
 * Stalled targets stop exactly as before; 0 disables the extension.
 * 2026-09-23, 5 h 21x21x21 run: of the 2,560 frontier charts left at
 * target_budget_spent, 712 were still falling at their last quantum (violation
 * at a median 0.40 of the first joint value; the pass's own forecast of the
 * time still needed p25/p50/p75 12/45/105 s) and 1,847 were flat or rising.
 * Measured neutral at a 1 h budget (180 s: 60.666% against 60.753%, 2,033
 * refinement admissions against 2,092 -- the extra quanta displace other
 * targets when there is no spare time); off, untested at long budgets. */
#define ASM_REGION_REFINE_CONVERGING_SEC 0.0
/* Metric steps keep changed faces inside the nominal band.  The step guard
 * compares singular values with the audit's comparison allowance
 * (ASM_METRIC_ROUNDOFF), so a face inside [0.75, 1.25] may step to
 * 1.25 + 1e-8 and rigid packing round-off (~1e-10) then crosses the strict
 * audit bound -- the reason for the refinement's 1e-8 interior reserve, which a
 * solve with an active metric constraint cannot meet.  1: a face inside the
 * nominal band before a step may not leave it by more than its own arithmetic
 * noise (128 ulp, or the rigid-motion round-off bound); faces already outside
 * keep the non-regression rule.  2026-09-23, 21x21x21 1 h: dropping the reserve
 * gained 0.40 points (61.155% vs 60.753%) with one strict failure, a face at
 * sigma_max 1.2500000100716; with this rule and the reserve at 0: 61.145%,
 * 62,485 passing seams (61,118), 0 strict failures, 0 overlaps.  Default on,
 * with ASM_REGION_REFINE_METRIC_MARGIN 0. */
#define ASM_REPAIR_STRICT_BAND 1
/* Holes before edges, half the time.  The refinement queue surveys every
 * unsurveyed frontier target before any joint continuation, and in a large
 * region each admission exposes new edge targets, so an enclosed target (at
 * least ..._HOLE_SIDES placed required-source partners: a hole inside the
 * piece) waits forever for its joint rounds.  1: every other pick takes the
 * best surveyed, still-requeued enclosed target when one exists.
 * 2026-09-23, 5 h 21x21x21 run, holes of the five largest pieces: region 1180
 * left all 250 surveyed-and-waiting (160 never got a joint attempt), 1986
 * 140 of 180 and 1972 99 of 131; chart 78304 (six placed partners) cut its
 * violation 16.8 -> 13.5 vox in its only joint quantum and was never picked
 * again.  Measured 2026-09-23 (1 h, pherc0139_21x21x21_holes_fresh_3600):
 * 60.653% against 61.145% (-0.49 points, 61,811 passing seams against
 * 62,485).  It does fill holes -- 78304 converged 16.8 -> 0.17 vox and was
 * admitted in its 9th quantum; the five largest pieces kept 658 holes against
 * 696 -- but with the same 54,000 refinement seconds it admitted 1,905 charts
 * against 2,594: a joint quantum makes 0-3 solver iterations, so half the
 * picks spent on holes cost more edge admissions than they gain.  Off. */
#define ASM_REGION_REFINE_HOLE_SHARE 0
#define ASM_REGION_REFINE_HOLE_SIDES 3
/* Gateways before edges, half the time.  The refinement ranks targets by
 * survey-first, then their own area per forecast second; the missing area
 * BEHIND a target (its unplaced branch) only breaks ties.  2026-09-23, 5 h
 * 21x21x21 run on the new defaults: 148.0 M vox^2 (16.6% of the source) lies in
 * 23 missing groups behind <= 5 refused boundary charts; region 1180's two
 * gateways (74312, 74334) hold 6,379 charts / 76.9 M vox^2 behind them, got one
 * 4 s joint quantum each (violation cut to 0.035 and 0.028; contact down to
 * one pair of < 0.05 vox^2) and were never picked again.  1: every other pick
 * takes the target, surveyed or not, with the most missing branch area per
 * forecast second, when something lies behind it (a single-chart hole is not
 * a gateway; ASM_REGION_REFINE_HOLE_SHARE chose those and lost).  Measured
 * 2026-09-24 with the requeued-only turn (pherc0139_21x21x21_gateway_fresh_3600):
 * 62.590% against 63.291% (-0.70 points; 1,506 refinement admissions against
 * 2,324).  Region 1180's gateway 74312 did get five joint quanta and stalled at
 * a residual contact violation of 0.026 (0.081, 0.052, 0.034, 0.029, 0.026,
 * 0.026) until its budget was spent: more time does not open these walls, so the
 * turn only takes time from edge admissions.  The surveyed-or-not variant
 * above was not run.  Off. */
#define ASM_REGION_REFINE_GATEWAY_SHARE 0
#define ASM_REGION_REFINE_CONVERGING_RATIO 0.95
/* Attached holes do not create standalone pieces. Keep this independent of
 * the caller's minimum output-region area (formerly both were 20,000). */
#define ASM_REGION_REFINE_MIN_BRANCH_AREA 1.0
/* In-memory diagnostic cache; full capacity stops insertion, never evicts
 * geometry evidence or changes acceptance. Each entry is under 100 bytes. */
#define ASM_ADMISSION_FAILURE_CACHE_ENTRIES 8192
/* An admission candidate whose initial contacts with the placed sheet
 * exceed this many pairs per candidate face lies on another wrap; it is
 * refused before the whole-field solve (the same 3-D truth the discovery
 * pass's same-surface guard uses, measured in UV contacts). */
#define ASM_ADMISSION_MAX_CONTACT_RATIO 1.0
/* The same pre-solve shortcut decided on overlap AREA and the 3-D truth
 * instead of a triangle-pair count, which scales with both charts' triangle
 * sizes (a coarse CVT candidate over a 0.6-vox source-retained chart meets
 * ~40 placed triangles per face).  Each blocking pair's UV overlap area is
 * split by the 3-D gap between the two materials at the same uv spot (each
 * triangle's own uv->xyz map, so triangle size does not enter): within
 * ASM_DISCOVER_SAME_SURFACE_VOX it is the same place on the sheet (a start
 * misfit a bounded solve may clear); beyond it the placed material lies
 * elsewhere -- another wrap or a misplacement no local solve moves.  The
 * candidate is refused before the solve when that far overlap exceeds
 * ..._MAX_FAR_OVERLAP of its area (the discovery pass's measured 15%) or any
 * overlap exceeds ..._MAX_PLACED_OVERLAP.  External pairs, which carry no
 * geometry here, count as far.  0 = the pair-count ratio above.
 * 2026-09-23, 21x21x21 seam-gate run: the ratio refused 885 charts for good
 * (1.1% of the source area; start overlap p10/p50/p90 16/32/61% of the
 * candidate), and in 3-D 44% of those attempts collided with material more
 * than 8 vox away while 21% touched only their neighbours.  Measured neutral
 * 2026-09-23 (1 h budget: 58.764% against 58.770%, +22 passing seams; it moves
 * refusals between reasons without admitting more), so it stays off. */
#define ASM_ADMISSION_AREA_3D_TEST 0
#define ASM_ADMISSION_MAX_FAR_OVERLAP 0.15
#define ASM_ADMISSION_MAX_PLACED_OVERLAP 0.5
/* The seam twin of the contact ratio: with only the candidate moving, a required incident
 * seam that starts more than this many vox off (5x the 4-vox seam bound) is refused before
 * the solve, with the refusal the solve would return.  Measured 2026-09-22 on the 21x21x21's
 * largest layer-cut region (1212): of 558 solves, the 118 that started > 20 vox off
 * accepted 0 and took 504 s; 10-20 vox accepted 17 of 107; 4-6 vox 70 of 82.  0 = off. */
#define ASM_ADMISSION_MAX_START_SEAM_RMS 20.0
/* Masked ends of a solve with fixed charts.  apr_solve opens and closes with a
 * complete-field objective (energy/overlap/contacts before and after) that also
 * gates the solve on the whole field being finite and within the walk budget.
 * In a growth repair only the candidate moves, so those two walks re-measure
 * the fixed context against itself: 39% of a sampled growth worker's time
 * (30 of 76 samples, 2026-09-23, 5 h run on the memo defaults).  Their results
 * reach only reports, and no chart of any 2026-09-23 run was finally refused on
 * a walk budget.  1: with fixed charts, both use the solve's own mask (the
 * closing one is usually a memo replay) and the fixed faces are checked for
 * validity face by face; the admission certificates that follow are unchanged.
 * Measured 2026-09-24, 1 h 21x21x21 arm on the closure-first defaults
 * (pherc0139_21x21x21_masked_ends_fresh_3600): 63.666% coverage against
 * 63.291% (+0.38 points), 48,565 placed charts (+233), 66,854 passing seams
 * (+488), 2,156 pieces, zero overlap, zero strict metric failures; faster
 * growth solves left refinement 61,432 s instead of 54,204 s (2,510 refinement
 * admissions against 2,324).  Default on. */
#define ASM_REPAIR_MASKED_ENDS 1
/* Incremental replacement index of a joint refinement.  The pass keeps the
 * current versions of every chart it admitted or moved in a delta index beside
 * the immutable incumbent, and rebuilt that whole index after each admission:
 * O(changed faces) per admission, O(n^2) over a region.  Region 1180 of the
 * 21x21x21 5 h run (worker 0, 2026-09-23) spent 10,505 of its 15,623 refine
 * seconds rebuilding it (1,957 admissions, 16.7 M faces at the end) and ended
 * with 444 targets still waiting, its two wall gateways among them.
 * 1: binary-carry blocks as in the growth pass; each chart's current version
 * lives in exactly one block and queries skip the stale copies left in older
 * blocks, which are purged when a block merges or is more than half stale.
 * Query answers are those of the full rebuild (apr_overlay_incremental_control).
 * Measured 2026-09-24: worker 0 rerun alone under both builds with one shared
 * deadline (region1180_ab/): identical growth (1,696 charts, 21,569,723.873
 * vox^2), then after ~2.15 h of refinement 1,972 admissions against 1,295
 * (3,973 attempts against 2,450) -- more than the whole 5 h run's 1,957.
 * It was off for a night because that worker committed 37.9 GB against 18.0
 * GB (a 32-worker run with it reached the machine's 287.8 GB commit limit).
 * The cause was the arena, not the index: every arena committed a 64 MB chunk,
 * and each waiting continuation and index block had its own.  With
 * Arena_new_sized for those (2026-09-24) the same worker holds 11.3 GB after
 * 2.2 h of refinement with 2,020 admissions (region1180_ab/inc_sized).
 * Default on. */
#define ASM_REPAIR_OVERLAY_INCREMENTAL 1
/* Rigid residual clearance at the admission endpoint.  A refinement quantum
 * that ends within one vox^2 of admission is often never cleared: in the 5 h
 * 21x21x21 run (2026-09-23) 2,420 targets outside the three deadline-bound
 * regions ended a quantum under 1 vox^2 of blocking contact; 1,138 were
 * admitted later (median one more quantum) and 1,281 never were.  The solve's
 * contact separation is a finite penalty balanced by metric and seam forces,
 * and the single-vertex polish follows an area gradient that vanishes when one
 * triangle contains the other (zero_gradient_visits).  1: when the complete
 * certificate refuses on at most APA_RESIDUAL_PAIRS new-material pairs, the
 * single new chart is translated rigidly away from its colliders in steps of
 * 1/16 to 2 vox; each trial re-checks metric, every incident seam and the
 * complete contact certificate, and a failed search restores the UVs.
 * 2: also growth's rigid and repair attempts (frontier 2 and 3).  Few apply:
 * in that run 2,392 of the rigid growth refusals on new-material contact came
 * from the exterior probe (no local pair; the nudge sees only the local field)
 * and only 55 rigid and 112 repair growth refusals ended on 1-8 local pairs.
 * Measured 2026-09-24, 1 h 21x21x21 arm with 2 on the masked-ends defaults
 * (pherc0139_21x21x21_nudge_growth_fresh_3600): NEUTRAL, 63.665% against
 * 63.666%, 48,547 placed charts against 48,565, 66,769 passing seams against
 * 66,854, zero overlap; 58 of 770 nudges qualified (new-material refusals
 * -28, seam refusals +16).  Off. */
#define ASM_REPAIR_RESIDUAL_NUDGE 0
/* The nudge's share of a quantum: a joint solve always spends its deadline,
 * so when the nudge applies the solve stops this much earlier (only in a
 * quantum longer than twice it) and the endpoint searches within the rest;
 * past the deadline the endpoint changes nothing.  At most six certificates. */
#define ASM_REPAIR_RESIDUAL_NUDGE_SEC 1.0
/* Each admission pass attempts its candidates in the order of their fixed
 * source anchors (then area) and stops after this many seconds; the rest
 * are logged as admission_budget_exhausted.  2026-09-16: the final pass on
 * the 4x5x5 took 89 min for 3 admitted charts even with the bounded
 * clearance -- every candidate builds and solves the whole field. */
#define ASM_ADMISSION_BUDGET_SEC 600.0
/* CONTACT WALK BUDGET (2026-09-17).  The 21x5x5 repair died asking the arena
 * for 256 GiB of leaf pairs: an admission's contact tree had been built with
 * the candidate at the uv origin and refitted after the seam fit moved it
 * thousands of vox, so every leaf that mixed candidate and origin material
 * spanned the move and the walk enumerated ~2^34 leaf pairs.  A sound single
 * cover visits ~8 overlapping leaves per leaf (closed-interval boxes), a
 * double cover ~30.  The walk stops at the budget (the serial, threaded and
 * touching walks count the same sequence, so they stop at the same pair),
 * AsmContacts_measure returns -2, and the callers refuse the trial, mark the
 * audit's overlaps unmeasured or end the phase: the sheet still exports. */
/* Fixed material-contact policy: the requested 4-voxel triangle-pair SAT band
 * applies to placement, repair, workers and audits without runtime switches. */
#define ASM_MATERIAL_CONTACT_TOLERANCE_VOX 4.0
#define ASM_CONTACT_MAX_LEAF_PAIRS_PER_LEAF 64
/* REGIME-BOUND (revalidated 2026-09-22 for the 21x21x21): a full audit walks
 * ~462,000 leaf pairs per audited Mface (4x5x5 403,790; 10x10x10 462,290),
 * so the 21x21x21's ~650 audited Mfaces need ~300 M pairs, past 2^27 =
 * 134 M, and an audit stopped at its budget can never qualify. 2^29 gives
 * ~1.8x headroom. The list grows by doubling only as a walk needs it, and
 * small fields stay governed by the per-leaf term. AcHit.lp is int32 and
 * Windows long is 32-bit, so this must stay below 2^31. */
#define ASM_CONTACT_MAX_LEAF_PAIRS ((size_t)1 << 29)
/* Contact walk memo (asm_contacts.c, 2026-09-23): a masked solver walk
 * requested again with byte-identical coordinates, mask, partition, mode and
 * gradient flag replays its committed hit sequence instead of re-measuring
 * the geometry (the sampled refinement worker spent ~74% of its solver time
 * in walks, most of them at already-measured coordinates).  1 = on every
 * replay also re-measure the geometry without committing and compare the hit
 * list bit for bit; a mismatch is reported and the fresh walk is used.  A
 * diagnostic build setting, 0 in production.  Measured 2026-09-23: the verify
 * build re-measured 4,814 real replays with 0 mismatches; the 1 h 21x21x21
 * arm with the memo replays 51% of masked walks, makes 52% more joint solver
 * iterations per attempt-second, and reaches 61.841% against 61.145% (+0.70
 * points, 63,696 passing seams against 62,485, zero overlap). */
#define ASM_CONTACT_MEMO_VERIFY 0
/* The admission repartitions the trial's contact tree AFTER the seam fit has
 * transported the candidate (the cause above); 0 = the refit-only path. */
#define ASM_ADMISSION_REBUILD_AFTER_FIT 1
/* LOCAL TRIAL FIELDS (2026-09-17).  An admission candidate used to build its trial field over
 * ALL registered charts -- 21-33M vertices of metric, witnesses and BVH per candidate -- so the
 * initial pass attempted 1 of 7,070 candidates on pherc343 and 4 of 1,098 on the 21x5x5 inside
 * its 600-s budget.  The field is built over the candidate's CONTEXT instead: its own charts,
 * every partner of an incident source relation, and the registered charts whose uv box comes
 * within (the candidate's own extent + ASM_ADMISSION_CONTEXT_VOX) of those partners' box.  The
 * material the field leaves out is checked by probing the candidate's faces against the incumbent
 * whole field's tree and against this pass's earlier admissions, with the same predicate, so
 * `pairs == 0` still means what it meant.  A candidate with no registered partner to localise
 * around (a layer seed) keeps the whole field.  0 = the whole field always. */
#define ASM_ADMISSION_LOCAL_FIELD 1
#define ASM_ADMISSION_CONTEXT_VOX 256.0   /* two cubes: apr_rigid caps a step at 4 vox and 0.05 rad, twelve rounds plus clearance */
/* A candidate holding a chart above this many faces whose OWN metric fails is refused instead
 * of cold re-flattening it: measured 2026-09-17, the 21x5x5 and pherc343 spend 100-190 s per
 * candidate on 4,000-step recoveries of 45-136k-face charts that end "proposal refused" and
 * the candidate is refused anyway; the 4x5x5 recovers nothing in its admissions. */
#define ASM_ADMISSION_RECOVER_MAX_FACES 20000
/* A placed layout whose 8-vox cells carry two charts a layer apart on more
 * than this fraction of its cells is a stack, not a sheet with seams to
 * close: the repair runs its bounded admissions and skips the source-repair
 * trials and the regions (whole-field solves over a field that cannot be
 * made sound).  The pre-check logs every pile's value (stage6_repair.json:
 * stacked_fraction). */
#define ASM_REPAIR_MAX_OVERLAP_FRACTION 0.05
/* DISCOVERY, the second pass over the placement (2026-09-16, the user's
 * two-pass rule): a leftover whose placed seam partners agree on its pose is
 * placed, even onto material the ribbon already carries when that material
 * is the same surface in 3-D.  Measured on the accepted 09-13 4x5x5 layout
 * with the Python instrument that preceded this module: partner-implied
 * poses of the PLACED charts land within 0.4 vox p50 of their actual poses
 * (98% within 10 vox); the 7 charts with 2-4 agreeing sides are exactly the
 * ones the Python postfix had admitted; 7 charts whose seam pose lands on
 * material 8-100 vox away in 3-D are other wraps that an overlap-tolerant
 * rule without the 3-D guard would have stacked; 18 single-sided charts sit
 * at the crop faces (the other sides do not exist). */
#define ASM_DISCOVER 1                    /* 0 = report only */
#define ASM_DISCOVER_MIN_SIDES 2          /* agreeing quadrants an interior chart needs; a crop-face chart needs 1 */
#define ASM_DISCOVER_SOFT_MIN_SIDES 2     /* ... and when every supporting seam is gate-rejected (switched off, weak, short): unanimous from this many sides */
#define ASM_DISCOVER_ROUNDS 8             /* admitted charts become partners in the next round */
#define ASM_DISCOVER_AGREE_VOX 10.0       /* two implied poses agree when no vertex moves more than this between them */
#define ASM_DISCOVER_OVERLAP_VOX 4.0      /* a vertex overlaps placed material within this uv radius */
#define ASM_DISCOVER_MAX_OVERLAP 0.15     /* more overlap than this needs the same-surface test */
#define ASM_DISCOVER_SAME_SURFACE_VOX 6.0 /* overlapped material within this in 3-D is the same surface (pitch 9.5) */
#define ASM_DISCOVER_RIM_VOX 12.0         /* boundary vertices: placed material within this in uv must be within it in 3-D */
#define ASM_DISCOVER_MIN_RIM 0.8          /* fraction of the tested rim that must be consistent */
#define ASM_DISCOVER_EDGE_MARGIN_VOX 16.0 /* a chart this close to the pile's extent touches a crop face */
/* A restoration level ends after an accepted step the line search had to
 * shrink below this alpha.  Real steps run 1e-4..0.5; the ones below are
 * no-ops (alpha 1e-10 left the seam violation unchanged to seven digits,
 * 2026-09-15) that each cost a refused QP and its OSQP fallback, 200 s at
 * 211k coordinates, for the rest of the level's twelve iterations.  The
 * level keeps its archive either way. */
#define ASM_RESTORATION_STALL_ALPHA 1e-5
/* ... and after an accepted step that cut the seam violation by less than
 * this fraction while the energy did not fall.  Measured 2026-09-16 on a
 * 219k-coordinate level: iterations 7-12 bound 2,200-2,300 rows, took
 * 300-480 s each, and moved the violation by 0.05% at alpha 5e-4 while the
 * energy rose; healthy iterations move it by 2-5%. */
#define ASM_RESTORATION_STALL_PROGRESS 1e-3
/* A closure step whose linearised model, at its unconstrained minimiser,
 * lies more than this many times the current energy below zero is not a
 * model of anything: the PSD-projected Hessian has lost its curvature
 * there and the ridge is all that bounds the step.  Measured 2026-09-15:
 * sane steps sit at 0.1x the energy; the four pathological solves (one
 * outgrew a 2,048-row working set in 227 s, three were infeasible) had
 * OSQP values of 1e21..1e24 against an energy of 1e7, and the steps they
 * yielded were a 12x energy blow-up and three no-ops.  The level ends and
 * keeps its archive. */
#define ASM_RESTORATION_MODEL_LIMIT 1e3

#define ASM_READMIT_VOX 100.0
#define ASM_READMIT_RAD 0.25
#define ASM_REFINE_READMIT_DROPPED 1     /* accepted seams the cleaning dropped are candidates too */
#define ASM_REFINE_MAX_CONTRA_GROWTH 0.10
#define ASM_REFINE_CONTRA_FLOOR 2.0e3
#define ASM_REFINE_ROUNDS 2
#define ASM_JOIN_AFTER_REFINE 1          /* the layout joins run once more on the refined layout */
/* Layout joins also MEASURE the in-plane offset of the joined piece: at every agreeing
 * cell pair the 3-D offset between the two claims is compared with what chart a's own
 * map predicts for the lattice step, and the difference, brought back to the sheet plane
 * through the face's Jacobian, is the shift that aligns b's drawing with a's
 * continuation.  The smaller lineage of every join is shifted by the mean over its joins
 * (capped at the step slack) BEFORE the relations are taken from the poses, so the join
 * closes the seam instead of freezing its misalignment (the familiar's layout joins
 * carried 77-101 lattice edges over the wrap gate each, 2026-09-09).  0 = off (A/B). */
#define ASM_JOIN_ALIGN 1
/* The alignment is a rigid 2-D fit per moving lineage (weighted Procrustes over its agreeing
 * cell pairs and seam offsets); the rotation is capped here (radians) and the translation at
 * the agree reach. */
#define ASM_JOIN_ALIGN_ROT_MAX 0.10
/* A rotation is fitted only when the evidence spreads over at least this radius (vox) about
 * its centroid; a short contact fixes the translation alone. */
#define ASM_JOIN_ALIGN_ROT_SPREAD 60.0
/* A pair with an existing relation carrying any of these flags is never re-joined by
 * the layout.  CONTACT only: the tangency gate's verdict stands, but a join the
 * cleaning DROPPED on a witness path is re-joined when the final layout shows the
 * two charts continuous (measured on the familiar 4x5x5, 2026-09-09: pairs
 * 285-312 / 285-480 / 493-498 / 296-328 with 39 / 16 / 11 / 10 agreeing cells and
 * none disagreeing were dropped joins). */
#define ASM_JOIN_VETO_FLAGS (ASM_REL_CONTACT | ASM_REL_DROPPED | \
                             (ASM_WRAP_VETO_JOIN ? ASM_REL_CROSSWRAP : 0u))
/* Output bands: every raster writer caps (TifXYZ at FLATTEN_MAX_GRID_DIM,
 * Rawtex at 2^28 px), so sheets are emitted in fixed-width u-bands. */
#define ASM_BAND_COLS 8192

/* Strip pages of a winding-ribbon bake (strip_layout.c, sheet_strip.exe; user
 * 2026-09-29: sheets "stripped", in spiral order, with the kibble that
 * survived the 40,000 vox^2 minimum piece moved to the end).  Presentation
 * only: pixels are moved, never changed, and no geometry is touched.
 * KIBBLE: a piece under 250,000 vox^2 of covered area (~500 x 500 vox, ~4.7 mm
 * square at 9.4 um) goes to the tray.  On the 12 boxes published 2026-09-29
 * such pieces are 76% of all pieces but hold 16% of the area: the piece count
 * falls steadily from a pile of chips just above the 40,000 minimum and levels
 * off at ~0.5-2 Mvox^2, where full-height sheets carry the area (no sharper
 * valley exists in the pooled histogram).  It trays 191 of 261 pieces (23% of
 * the area) on PHerc0813 10^3 and 48 of 101 (5.6%) on PHerc0139 10^3.
 * Measured by eye on both: 100,000 left the spiral lines busy with chips,
 * 500,000 trayed medium sheets (36% of 0813's area). */
#define STRIP_KIBBLE_AREA_VOX2 250000.0
#define STRIP_PAGE_ASPECT 1.6           /* page width / height: a landscape screen */
#define STRIP_LINE_GAP_VOX 96.0         /* blank rows between lines */
#define STRIP_PIECE_GAP_VOX 48.0        /* blank columns between pieces on a line */
#define STRIP_MARGIN_VOX 96.0           /* blank border */
#define STRIP_TRAY_BREAK_LINES 3.0      /* extra line gaps before the tray */
/* Line flow: a main piece narrower than this share of a line moves whole to
 * the next line instead of being split; a split falls on the piece's thinnest
 * column within the last STRIP_SOFT_BREAK_FRACTION of the line; with less than
 * STRIP_MIN_ROOM_FRACTION of the line left, the next piece starts a new line. */
#define STRIP_WRAP_WHOLE_FRACTION 0.12
#define STRIP_SOFT_BREAK_FRACTION 0.06
#define STRIP_MIN_ROOM_FRACTION 0.05
/* Balance: no lone piece on the last spiral line before the tray.  Among line
 * widths within +-STRIP_BALANCE_SPAN of the best aspect fit (and within
 * STRIP_BALANCE_ASPECT_SLACK of its |log aspect error|), the page takes the
 * best aspect whose last spiral line is at least STRIP_BALANCE_GOOD_FILL
 * full, else the fullest.  A widening of ~5-7% pulls a widow up into the line
 * above; 0.06 slack refused it on PHerc0813 (last line 18% full) and
 * PHerc0175B (7%), 2026-09-29. */
#define STRIP_BALANCE_SPAN 0.10
#define STRIP_BALANCE_STEPS 64
#define STRIP_BALANCE_ASPECT_SLACK 0.15
#define STRIP_BALANCE_GOOD_FILL 0.5
/* sheet_strip's quick-look PNG: longest side at most this many pixels. */
#define STRIP_PREVIEW_MAX_PX 4096

#endif
