# Quadribbon big-sheet pipeline — canonical order of operations

Status: 2026-09-01, fifth revision. The lane is ONE tool in the order the
problem demands: mesh → initial winding/unwrap → parameterize → fit the
solid quad ribbon → optimize it (refitting as needed) → big sheet. The
2026-09-01 4x5x5 champion (`output/pherc0139_4x5x5_step11_peel_20260901`,
7/7 gates) was produced by a hand-driven chain; the tool now reproduces its
certificate and fit byte-identically (`output/pherc0139_4x5x5_solid_lane_
20260902`, stages 2-3, every sidecar) and continues from there.

1. **mesh** — the meshed region: a welded world-frame VMESH
   (`meshes/PHerc0139-<rung>/head_*/weld/welded.vmesh`; the measured 4x5x5
   champion path — the 08-30 "never weld" rule was overturned by the 09-01
   measurement, pile certificate coherence 17.2% vs 95%) or the pre-weld
   per-cube pile (`meshes/PHerc0139-<rung>/head_*/mesh/dump/`, concatenated
   without welding; the 21x route once relation error reaches ~0).
   Artifact: `stage1_mesh/mesh.vmesh` (+ `mesh_source.json` or
   `mesh_cubes.tsv`).
2. **unwrap** — `scroll_ribbon --winding-only`: initial winding field,
   winding registration + MRF, gauge sync, the arc-length U map; written as
   the winding certificate `stage2_unwrap/unwrap_winding.vmesh` (+ the four
   `_winding_*.f32` sidecars, topology preserved).
3. **fit** (parameterize) — `scroll_ribbon --scaffold-solve
   --component-global` on the certificate: the claims-mode lattice ribbon
   (`QP_FIT_*` constants: grid_u 2, slice_h 2, iters 1/1/0; peel bands;
   wrap-safe emission). These are the OBSERVATIONS. `stage3_fit/fit_ribbon.
   vmesh` + `_support.u8 _phase.f32 _material_identity.i32
   _reconstruction_component.i32 _stats.json`, `verdict.json`
   (ribbon_verdict, `--source` = the certificate, core radius 64), root
   `stage3_fit_bake.png`.
4. **solid** — `src/flatten/quadribbon_solidify.c`: ONE solid quad ribbon
   per atlas column run (runs are separate wraps, never bridged). Observed
   cells exact; holes filled by quad_strip's cylindrical residual harmonic
   solve on TAUCS (factor once per run, all channels in one solve; per-row
   trend, so fragments out of angular order cannot flip a one-sided fill).
   Fill domain (user rule 2026-09-01): every row of every column between the
   first and last GOOD column of a run (>= `QS_GOOD_COLUMN_FRAC` of rows
   observed) -- interpolation along u, continuation along v to the box faces
   -- inside the observed box + one wrap gate and outside the umbilicus core.
   Non-good stretches wider than `QS_FILL_MAX_U_GAP` stay empty; a gap or a
   row/column bridge is refused when it fails the wrap-skip test (radial
   extrapolation from each side, `QS_SKIP_DR_PITCH`), the angular-order test
   (u is arc length: angle change = span x the run's angle per column,
   `QS_ORDER_TOL_*`), the fit's stretch gate on short bridges, or joins two
   material lineages; runs below `QS_RUN_MIN_OBSERVED` cells are confetti and
   are not emitted. Every vertex labeled (`_provenance.u8`: 0 observed-fixed, 1 observed-generated,
   2 adopted peel-band, 3 fill). Peel band 0 wins a cell; a band-1 cell is
   adopted only when its radius is within `QS_BAND1_ADOPT_DR_PITCH` of its
   band-0 neighbourhood, otherwise it is another wrap: counted in
   `solid_ribbon_report.json` (duplicate cells + |dr|/pitch quantiles), never
   planted. Emission re-applies the 6-vox wrap gate and material purity per
   triangle. `stage4_solid/solid_ribbon.vmesh` + sidecars,
   `solid_ribbon_provenance.png`, `verdict.json`, root `stage4_solid_bake.png`.
   4x5x5: 2.3 s for 25 runs; 6/7 gates — coverage 81.9% vs the fit's 95.9%
   because 215,118 peel-band cells are a different wrap at the same (u,v);
   142k fills, distortion beyond 4x 0.0031%.
5. **optimize** — rounds (`QS_OPTIMIZE_ROUNDS`) of untangle
   (`src/flatten/quadribbon_untangle.c`: turn-order repair + elastic shell
   + settle, radial motion only, all vertices movable) followed by a REFIT:
   the same solidify call on the moved ribbon, whose provenance sidecar says
   which cells are observations (kept where the shell put them) and which
   were fills (re-solved). Stops at a conflict-free state, the round cap, or
   when long conflicts drop by less than `QS_OPTIMIZE_MIN_GAIN`. Optional
   final metric re-solve (`QS_REPARAM_FINAL`, default off: the certificate
   gauge is retained). `stage5_optimize/round_N/{untangled_ribbon.vmesh,
   untangled_ribbon_dr.png, untangle.json, refit_ribbon.*, verdict.json}`,
   `optimized_ribbon.*`, `optimize_report.json`, root `stage5_optimize_bake.png`.
6. **sheet** — layer strips (hard/soft conflict split), RAW bakes, C
   first-cover composite with provenance and the quality-masked shipped
   view. "pre" = the solid ribbon, "post" = the optimized one.

Phase 2 (gated on the phase-1 bakes): the soft-anchor consistency fit
(`QuadStrip_metric_arap` on each run rectangle before compaction, observations
as robust Huber anchors, `QS_ARAP_*`). Phase 3: CT-ridge snap of fills
(`QuadStripSnap_run`, `QS_SNAP_ENABLE`).

## The one command

```
build\Release\quadribbon.exe <mesh_dir|welded.vmesh> <out_dir> [config.json] [--stop-after <stage>]
```

- Policy comes from `configs/default.json` (axis, wrap
  pitch, RAW zarr, bake window/dark, `compute.threads` — 8 matches the
  original reference's thread count) + `pipeline_constants.h` (`QP_FIT_*`, `QS_*`).
  Set the CT and executable paths for the local build. The shared release
  configuration includes an axis table; historical byte-identity claims above
  describe the original September profile, not a validation of today's defaults.
  No tuning flags.
- Stages resume: a stage whose primary artifact exists is skipped (the solid
  and optimize stages also require every sidecar at the right length).
- `--stop-after mesh|unwrap|fit|solid|optimize|sheet` (or 1..6) finishes the
  named stage, verdict and bake included, then exits.
- Outputs: `stage1_mesh/`, `stage2_unwrap/`, `stage3_fit/`, `stage4_solid/`,
  `stage5_optimize/`, `stage6_sheet_{pre,post}/`, root
  `stage{3,4,5}_*_bake.png`, `big_sheet_{pre,post,provenance,pre_vs_post}.png`
  (+ `_masked`). Every spawned stage logs to `logs/<stage>.log`;
  `logs/pipeline.log` is the elapsed-stamped timeline.
- `quadribbon.exe --selftest` runs the composite, untangle, TAUCS wrapper,
  quad_strip and solidify module tests; `scroll_ribbon --selftest` covers the
  unwrap/fit machinery.

## Diagnostics

`python/scripts/mesh_cross_section_compare.py` is the ONE cross-section
viewer (exact plane intersections on RAW zarr slices, U-colored). It now
also takes `--fitted-labels file.i32` (label-colored panel) and
`--uv-box U0 U1 V0 V1` (ring a big-sheet UV box's vertices, print its 3-D
locus, auto-pick z slices from the box) — the artifact-chase workflow.

## Retired from this lane (2026-08-30 cleanup)

Deleted; recover from git history at the cleanup commit if ever needed:

- `scripts/flatten/rebuild_quad_ribbon_big_sheet.ps1`,
  `iterate_quad_strip_bake_refit.ps1`, `run_quad_ribbon_tracks.ps1`,
  `run_10x_quad_ribbon_when_ready.ps1`, `bake_gmg_level_big_sheets*.ps1`,
  `make_quad_strip_{compare,preview}.ps1` — the driver/runner sprawl the
  one-command tool replaces (the v92-proxy solid lane they drove is not part
  of this workflow).
- `scripts/flatten/dump_ribbon_cross_sections.py` — one-off viewer, folded
  into `mesh_cross_section_compare.py`.
- `src/flatten/disk_flatten.{c,h}` — dead (no build target, no callers).
- `python/scripts/composite_layer_sheets.py` remains for ad-hoc use but the
  deliverable chain composites in C.

`scroll_ribbon.exe`, `quadribbon_strips.exe`, `obj_bake_raw.exe` remain as
stage-local debug entry points over the same modules.

## Key findings this lane is built on (2026-08-29/30)

- The consistency solve is TAUCS IRLS over assembled normal equations
  (supernodal multifrontal Cholesky, 2M unknowns in ~8 s; 512 MB stack
  required). The matrix-free CG + ad-hoc compensations were deleted.
- Crack ties: horizontal fit cracks had NO tie mechanism; the contact scan
  now probes slice+1 gated by per-side face-tie coverage. 28,597 ties at 4x
  carried the visible ±20–25 px slips; sheet seam shift median 3 px after.
- Ply gates: perp ≤ 3.5 vox (a true ply sits ≥ 7 away), tangent alignment,
  authoritative-phase Δφ; contacts get a wider Huber (c=6) than verticals
  (c=1).
- Island placement runs TWICE: the carried-u phi→u map is atlas-scrambled at
  the inner turns (PAVA flattens it; charts collapse onto one u plateau);
  the second round maps against the solved u.
- The untangler's gate in solid_quad_ribbon (`--cyl-axis` only) had silently
  disarmed it in the production solid step (`--snap-axis`); the extraction
  makes it an unconditional pipeline stage instead.

## Open items

- Phase 2 consistency fit and phase 3 CT-ridge snap (see stage 4/5 above),
  each gated on bakes.
- The single-ribbon coverage trade (peel-band duplicates are other wraps at
  the same certificate u) is the winding-certificate campaign's problem:
  island-scale relation aggregation, then the 10x and the pile lane.
- 4x21x21 / 21x21x21 need stage-1 per-cube meshing (running off the frozen
  fleet) before they can run.
- Remaining TAUCS migrations (legacy Stage-C strip solve incl. the parked
  RibAMG/GMG, quad_strip metric ARAP multigrid, snap_cg, SLIM, sheet_probe,
  adc ADMM) — see the changelog's solver inventory. The quad_strip HARMONIC
  fill is migrated (TAUCS factor-once; legacy PCG kept under
  `harmonic_solver = 0` / `QUAD_STRIP_NO_TAUCS`).
