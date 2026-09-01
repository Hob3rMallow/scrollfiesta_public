# Supported pipeline

This document expands the quick start in the
[repository README](../README.md). It describes the supported native
whole-scroll workflow and the artifacts expected from a reviewable run.

```text
                         +-> scroll_unroll
grid_pipeline -> scroll_whole
                         +-> atlas_overlap_fix -> atlas_ribbon_fit
```

The unroll and atlas/ribbon paths consume the same registered placement. They
can be run independently after `scroll_whole` completes.

## 1. Mesh the grid

```powershell
build\Release\grid_pipeline.exe GRID output\run --halo 13
```

The current default uses CVT/RVD simplification and a zero trim inset. Important
outputs include:

- `output/run/dump/`: resumable per-cube geometry;
- `output/run/logs/`: per-cube logs;
- `output/run/pipeline_summary.csv`: aggregate status and metrics;
- `output/run/welded.obj`: optional direct weld when requested and practical.

Use `--skip-existing` to resume an interrupted grid run.

The zero inset is important for whole-scroll registration: registration pairs
skins across cube seams, so charts must reach the cube faces. `--trim-inset 1` is an
explicit legacy choice for weld-focused experiments where the seam bridge
expects a gap. On the PHerc0139 4×5×5 fixture, the zero
inset produced 6,441 seam pairs and 2.98% whole-turn error; a one-voxel inset
reduced the evidence to 2,713 pairs and raised the error to 15.78%.

The per-cube multi-winding BPA growth gate is optional and requires scroll
geometry. For PHerc0139, the current calibration is:

```text
--umb-y 3405 --umb-x 2878 --wrap-pitch 9.5
```

Those parameters also enable the winding/phase gates in `grid_weld`. Keep them
off for uncalibrated datasets until the scroll geometry is known.

## 2. Register and place cubes

```powershell
build\Release\scroll_whole.exe output\run\dump output\run_placed
build\Release\scroll_whole.exe output\run_placed --reregister --audit
```

The placed directory contains one `*_placed.obj` plus UV, skin, and group
records for each usable cube. It also contains `placed_index.json`,
`audit.json`, and logs. Treat a failing registration or seam audit as a failed
run; do not promote only the visually plausible artifacts.

## 3. Unroll, snap, and bake textures

```powershell
build\Release\scroll_unroll.exe output\run_placed output\run_unroll `
  --raw GRID\cubes_RAW --steps 12345 --id run
```

The default stage sequence is:

| Stage | Operation |
|---|---|
| `1` | Base texture bake. |
| `2` | Seam join and relaxation. |
| `3` | Overlap ownership. |
| `4` | CT surface snap. |
| `5` | Final light relaxation. |

For each retained stage, the unroll directory can contain:

| Pattern | Meaning |
|---|---|
| `*_rawtex.tif` | Full-resolution grayscale CT texture. |
| `*_rawtex_strip.tif` / `*.png` | Texture tiled into a readable strip. |
| `*_rawtex_preview.png` | Downsampled review preview. |
| `*_diagclass.tif` / `*.png` | Validity and ownership diagnostics. |
| `*_xyzmap*.png` | World-coordinate diagnostics. |
| `*_seamzoom.png` | Focused seam review. |
| `*_stats.json` | Stage metrics. |

The directory also contains final scroll/leftover OBJs and
`pipeline_stats.json`. Run `scroll_unroll --help` for optional TIFXYZ and atlas
exports.

## 4. Resolve atlas overlaps

Registration decides where a chart sits within a winding. This stage decides
which winding owns it:

```powershell
build\Release\atlas_overlap_fix.exe output\run_placed output\run_atlas `
  --rounds 7
```

Each round alternates a tabu search over integer winding depths with a
continuous per-chart relayout, feeds the field into the next round, and stops
early when no chart moves.

Expected outputs include before/after UV OBJs, per-round winding-layer OBJs,
`charts.csv`, `groups.csv`, `diag_before.png`, `diag_after.png`,
`overlap_fix_stats.json`, and `atlas_solution.bin`.

## 5. Fit the ribbon

The ribbon stage slices the registered atlas into constant-v rows and fits one
collision-free ribbon through them:

```powershell
build\Release\atlas_ribbon_fit.exe output\run_placed `
  output\run_atlas\atlas_solution.bin output\run_ribbon --mode ribbon
```

Expected outputs include `ribbon.obj`, `ribbon_rows_world.obj`,
`observations_world.obj`, `u_overlap_audit.json`, the `u_*` collision and
registration ledgers, and `ribbon_fit_stats.json`.

`--mode` selects intermediate ablation stages (`f0` observations through `f6`
collision). `--raw DIR` additionally audits adjacent-row RAW texture phase.

## Canonical-run contract

[`scripts/run_canonical_grid.ps1`](../scripts/run_canonical_grid.ps1) runs the
five-stage unroll review alongside the atlas/ribbon branch. A canonical run is
incomplete unless its `snap_relax/` directory includes the step 4 snapped and
step 5 final-relaxed:

- full-resolution TIFF;
- readable TIFF/PNG strip;
- preview PNG;
- metrics JSON.

The script creates a `CANONICAL_OUTPUTS.md` inside that generated run, records
the artifact paths in `logs/SUMMARY.txt`, and fails when a required artifact is
missing or empty. `CANONICAL_OUTPUTS.md` is a generated per-run index, not a
repository-root document.

When promoting a run into the tracked result bundle, preserve that evidence
and update the checksums described in
[`output/canonical_best/README.md`](../output/canonical_best/README.md). Large
OBJ, BIN, TIFF, and PNG files in that bundle are stored with Git LFS.

## Output conventions

- CLI-written OBJ vertices use internal `(z,y,x)` column order.
- The public C API uses ordinary `(x,y,z)` coordinates.
- Texture rasters use unroll `(u,v)` coordinates.
- JSON metrics record the canvas origin, dimensions, fill, multi-coverage,
  darkness, and seam measurements.
- Probe, ablation, worker, and ordinary generated outputs belong under
  `output/` and should remain untracked.
