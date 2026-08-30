# Quadribbon big-sheet pipeline — canonical order of operations

Status: 2026-08-30, fourth revision. The lane was consolidated into ONE tool
after the metric-projection campaign kept running on inputs in unknown states
(a fit that never saw the collision untangler; a v144 solid ribbon from a
different lineage). The order of operations is now fixed and enforced by the
tool itself:

1. **mesh the region** — the PER-CUBE pipeline only (`cube_mesh` via
   `grid_pipeline`). **No welding.** Welding is slow and introduces
   artifacts; avoiding it is one of the reasons the quadribbon approach
   exists. Canonical per-cube piles live at
   `meshes/PHerc0139-<rung>/head_*/mesh/dump/`.
2. **fit quadribbon + reoptimize UV** — claims-mode fit (winding
   registration, branch lanes, phase-continuation island JOIN) followed by
   the metric projection (row slicing, F(u_ref), island placement by
   authoritative phase, TAUCS IRLS consistency solve with crack/contact
   ties, a second placement round against the solved u, global re-zero).
3. **untangle** — the collision shell (`src/flatten/quadribbon_untangle.c`:
   turn-order repair + elastic-shell contact ledger + settle-toward-rest,
   extracted from solid_quad_ribbon's lattice version and generalized to UV
   buckets + mesh adjacency), then **re-optimize UV** on the moved geometry.
4. **big sheet** — layer strips (hard/soft conflict split), RAW bakes, C
   first-cover composite (`src/flatten/sheet_composite.c`) with provenance.

## The one command

```
build\Release\quadribbon.exe <mesh_dir> <out_dir> [config.json]
```

- `<mesh_dir>`: the pre-weld per-cube pile (recursive `*_final_all.vmesh`;
  plain world-frame `*.vmesh` at the top level also accepted). Cube ids are
  parsed for offsets when a mesh is cube-local; the canonical dumps are
  already world-frame.
- Policy comes from `configs/flatten/pherc0139.production.json` (axis, wrap
  pitch, RAW zarr, bake window/dark) + `pipeline_constants.h`. No tuning
  flags.
- Stages resume: a stage whose primary artifact exists is skipped.
- Outputs: `stage1_soup.vmesh`, `stage2_fit/`, `stage3_reopt/`,
  `stage4_untangle/`, `stage5_reopt/`, `stage6_sheet_{pre,post}/`, and at
  the root `big_sheet_{pre,post,provenance,pre_vs_post}.png` ("pre" =
  before untangle, "post" = final). Every spawned stage logs to
  `logs/<stage>.log`.
- `quadribbon.exe --selftest` runs the composite + untangle module tests;
  `scroll_ribbon --selftest` covers the fit/reopt machinery (incl. the
  phase-continuation join and the two-round placement).

Cross-cube continuity is the FIT's job (claims ownership + continuation
gates + the phase-interval join in `src/flatten/unwrap.c`). Measured on the
4x5x5 pile: 100 per-cube meshes, 783 raw components → ~25 material
candidates with continuation saturation 1.000 — cube seams close without
any geometric welding.

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

- 10x10x10 rung from its pre-weld pile (replaces the unknown-state 0163
  input).
- 4x21x21 / 21x21x21 need stage-1 per-cube meshing before they can run.
- Remaining TAUCS migrations (legacy Stage-C strip solve incl. the parked
  RibAMG/GMG, quad_strip MG, snap_cg, SLIM, sheet_probe, adc ADMM) — see the
  changelog's solver inventory.
