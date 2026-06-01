# Best Results — final meshes

Curated showcase of the final output meshes from the two best pipeline runs
(both 2026-05-31). Only the **final-stage** mesh is included for each — no
intermediate stage dumps (step0 / pre_simplify / post_lop / holefill / bridge_cut).

## `grid_weld/` — grid-weld cube test (final welded grid)

Source: `output/grid_4x5x5_20260531/` (latest full run, 16:44).
The end product of the grid path: 100 per-cube charts welded across the grid
(`grid_weld` → `SeamWeld_bridge` → interior hole-fill → sliver/T-junction cleanup).

| File | What |
|------|------|
| `welded.obj` | Final welded surface mesh, 100 cubes stitched (165 MB) |
| `welded.obj.weld_report.json` | Weld + manifold audit |
| `pipeline_summary.csv` | Per-cube pipeline timing/component table |

Weld audit (`welded.obj.weld_report.json`):
- cubes processed: **100**
- unique verts / faces: **1,628,178 / 3,022,017**
- manifold pairs: **4,487,120** · non-manifold edges: **133** · same-dir pairs: **76**
- boundary (unpaired) edges: **91,260** (expected — grid-edge open boundary)

## `kaggle_cubes/` — final per-cube meshes

Source: `output/train_individual_20260531/` (7 canonical Kaggle eval cubes,
nnunet-preds, run individually with `--halo 0`). Final stage = QEM (7.5% faces);
`<id>_qem_all.obj` is the all-components mesh for that cube.

| Cube | Verts | Faces | Comps |
|------|------:|------:|------:|
| 86701140    | 178,329 | 338,059 | 14 |
| 118632705   | 149,865 | 284,102 | 11 |
| 1006462223  | 338,400 | 639,787 | 22 |
| 1013184726  | 243,669 | 457,612 | 21 |
| 3290306825  | 193,229 | 353,202 | 23 |
| 3294954456  | 177,339 | 333,421 | 18 |
| 3394433588  | 140,036 | 265,323 | 15 |

All 7 cubes completed cleanly (exit 0). See `run_summary.txt` for full timings.
