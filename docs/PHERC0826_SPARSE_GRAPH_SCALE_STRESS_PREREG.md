# PHerc0826 sparse-graph scale and complete-RAW stress test

Status: **frozen by the public commit containing this document, before any
voxel outside the already published 25-cube PHerc0826 box was accessed for
this experiment**. The fixed output root did not exist when this protocol was
written.

This is a post-result **scale and visual stress test**, not a second
independent holdout. The earlier preregistered 1 x 5 x 5 PHerc0826 result is
already known: the candidate reduced whole-turn seam errors from 44/1,076 to
7/1,076. This protocol mechanically expands that known center, separates all
new-only seams from already observed seams, supplies a complete one-cube RAW
halo, and runs the supported full-strip consumer. The known center and its
boundary can never make the scale-efficacy decision pass.

No box, arm, binary, axis, threshold, visual range, stage, retry setting, or
decision rule may change after the first expanded-box voxel read. Positive,
null, inconclusive, and negative results will be published.

## Question and non-duplication check

The tested implementation remains candidate
`a47a21c3a9d445c90d28e218b8f0fcbecd49c079` against its exact parent
`e0cf51cab03f6fe84e43e8c7ff5be9051d951b1d`. It admits sparse observed
cross-cube group-pair constraints after the existing geometry and phase gates,
then prevents graph-supported relations from being overwritten by local
polish. It does not generate winding annotations or replace the global
unroller.

A fresh check covered upstream ScrollFiesta, open Villa PRs, indexed GitHub
code, `abundantjoe/winding-sync`, and the supplied 15-channel Discord archive
(3,523 messages through 2026-08-08). No matching implementation or comparable
parent/candidate multi-layer test was found. `winding-sync` is complementary:
it generates relative constraints from CT and performs global L1 integer
synchronization, whereas this change consumes ScrollFiesta's existing mesh-
seam correspondences. The Discord archive also emphasizes that local pitch
varies sharply, verified constraints matter, and a poor spiral can look
visually convincing. Exact pair metrics and complete visuals are therefore
both mandatory here. This search is bounded: private work, unindexed code,
August 9-10 Discord messages, and the dedicated ScrollFiesta channel are not
present in the archive.

The experiment addresses the two disclosed limitations of the first result:

1. registration was evaluated on only one 25-cube z slab; and
2. the downstream consumer found only 25 RAW cubes and missed 122 requested
   neighbors, leaving the complete texture visually sparse.

It does not claim a new surface model, full-scroll tracing, ink recovery,
legible letters, or statistical independence from the first PHerc0826 result.

## Frozen data and mechanical expansion

- RAW level-0 CT:
  `s3://vesuvius-challenge-open-data/PHerc0826/volumes/20250821151701-9.362um-1.2m-113keV-masked.zarr`
- m7 level-0 surface prediction:
  `s3://vesuvius-challenge-open-data/PHerc0826/representations/predictions/surfaces/20250821151701-surface-20260413222639-surface-m7-L0-th0.2.zarr`
- Both arrays: shape `[16920,8169,8169]`, dtype `uint8`.
- Fixed chunk/cube edge: 128 level-0 voxels.
- Fixed axis point z/y/x: `[0,4829,4118]`.
- Fixed axis direction z/y/x: `[1,0,0]`.

The already observed center is exactly the earlier frozen box:

- known box z/y/x:
  `[12800,12928) x [4480,5120) x [3840,4480)`;
- known lattice: `1 x 5 x 5 = 25` cubes.

Expand that box by exactly one cube on every face to define the registered
core:

- core box z/y/x:
  `[12672,13056) x [4352,5248) x [3712,4608)`;
- core lattice: `3 x 7 x 7 = 147` possible prediction cubes.

Expand the core by exactly one further cube on every face to define the RAW
halo:

- halo box z/y/x:
  `[12544,13184) x [4224,5376) x [3584,4736)`;
- halo lattice: `5 x 9 x 9 = 405` RAW cubes.

The full halo is carved once for both sources. Prediction TIFFs inside the
core are then hard-linked into a separate meshing grid; prediction cubes in
the outer RAW shell are never registered. The same constant axis and all
registration defaults are deliberately retained. Drift or failure under this
mechanical expansion is a result, not permission to refit the axis.

## Frozen code, binaries, and verification tools

Candidate registration executable:

- source behavior is unchanged from candidate commit
  `a47a21c3a9d445c90d28e218b8f0fcbecd49c079` in the current documentation and
  evaluation branch;
- MSVC Release/x64 `scroll_whole.exe`: 273,408 bytes, SHA-256
  `713500073e0bd48ff0a0611b6c68ccca117bee9537d61260132242982f6aa4b8`;
- complete in-process selftest: PASS, including sparse-edge and graph-polish
  guard cases.

Exact parent comparator:

- commit `e0cf51cab03f6fe84e43e8c7ff5be9051d951b1d`, tree
  `fcf7f285eaead55589918437007d5362c06b8b20`;
- `scroll_whole.exe`: 270,336 bytes, SHA-256
  `eb1445504d69375f36891081f39ee6aa7d65a8af22b46b164cc75bbabce2db8a`.

Shared meshing and downstream executables:

- `grid_pipeline.exe`: 414,720 bytes, SHA-256
  `98b9e8e6a6de9083d47c3527b79d77f6338afa212f40bbb37ab0c0c929876eb6`;
- `cube_mesh.exe`: 1,173,504 bytes, SHA-256
  `c08437db47a73158f5d1ca7abded140ab7333e9a304ca24d19ea5818829ea72e`;
- `scroll_unroll.exe`: 768,000 bytes, SHA-256
  `848969a62c9765e07629fe1bb1ad217606ba6ecab6e02dd963ea9928f8335844`;
- complete `scroll_unroll --selftest`: PASS.

Frozen Python tools:

- carver `python/scripts/carve_grid_tifs.py`: 6,049 bytes, SHA-256
  `19e98e83d3b54259f23ac22bbf8286fe7cfdab85fbbb442528182870bf394811`;
- input validator/core linker
  `python/scripts/prepare_winding_scale_grid.py`: 9,717 bytes, SHA-256
  `9e7ae57ced49cc5d5f937daa9f0a092a948a1881363649e9c175d26827099018`;
- scale comparator/three-view renderer
  `python/scripts/compare_winding_scale_stress.py`: 24,351 bytes, SHA-256
  `5294ae7253edaef02084af1f9837d8a177532e270f7d21bfc0295bc45db2ee97`;
- downstream comparator `python/scripts/compare_unroll_outputs.py`: 15,659
  bytes, SHA-256
  `9893c56edf0d2d967ab7cc578540e065c41880a1872a8b0caea9c482b9472bdb`.

Nine focused comparator, mismatch-rejection, repeat, renderer, and grid-
validation tests pass. The new test files are themselves frozen:

- `test_compare_winding_scale_stress.py`: 6,113 bytes, SHA-256
  `f961cfe18c52f6de91b3aef65fd5e71d9549d337ce0e8a6e748662f810b38b77`;
- `test_prepare_winding_scale_grid.py`: 2,820 bytes, SHA-256
  `bfb5ba5a3392550e2bcb8e68456dd8e75df2c9ce0a45243bf357c5fd926dfebb`.

## Frozen registration execution

Run from the focused PR checkout. The fixed root is
`output/experiments/pherc0826_sparse_graph_scale_stress_20260810`; it must not
exist at start.

```powershell
$root = 'output\experiments\pherc0826_sparse_graph_scale_stress_20260810'
$pred = 's3://vesuvius-challenge-open-data/PHerc0826/representations/predictions/surfaces/20250821151701-surface-20260413222639-surface-m7-L0-th0.2.zarr'
$raw = 's3://vesuvius-challenge-open-data/PHerc0826/volumes/20250821151701-9.362um-1.2m-113keV-masked.zarr'
$py = '..\_sf_sparse\python\.venv\Scripts\python.exe'

& $py python\scripts\carve_grid_tifs.py `
  --pred-zarr $pred --raw-zarr $raw `
  --bbox 12544 13184 4224 5376 3584 4736 `
  --umbilicus 4829 4118 --out "$root\halo_grid"

& $py python\scripts\prepare_winding_scale_grid.py `
  --halo-grid "$root\halo_grid" --core-grid "$root\core_grid" `
  --halo-bbox 12544 13184 4224 5376 3584 4736 `
  --core-bbox 12672 13056 4352 5248 3712 4608 `
  --required-pred-bbox 12800 12928 4480 5120 3840 4480 `
  --pred-zarr $pred --raw-zarr $raw --umbilicus 4829 4118 `
  --min-core-pred 100
```

The validator must find exactly 405 RAW TIFFs, at least 100 nonempty core
prediction TIFFs, and all 25 known-center prediction TIFFs. Every file must be
uint8 `128^3`; prediction values must be only 0/255 and nonempty. The manifest,
inventory, per-file nonzero count, bytes, and SHA-256 are retained. Globally
zero RAW is a hard failure; individual all-zero RAW cubes are disclosed but
allowed because empty physical space is valid.

Mesh exactly the hard-linked core inventory once:

```powershell
& '..\_sf_calibration\build\Release\grid_pipeline.exe' `
  "$root\core_grid" "$root\mesh" `
  --halo 13 --threads-per-cube 1 --max-concurrent 8 `
  --simplify cvt --trim-inset 0 --skip-weld `
  --exe '..\_sf_calibration\build\Release\cube_mesh.exe'
```

Every present core prediction cube must exit zero and emit a nonempty complete
`step12_final` OBJ. The reject file must remain empty. No retry with a different
simplifier, concurrency, halo, reject override, cube subset, or box is allowed.

Run the exact parent, candidate, and one fresh candidate repeat on the same
mesh dump:

```powershell
& '..\_sf_baseline_e0cf\build\Release\scroll_whole.exe' `
  "$root\mesh\dump" "$root\baseline" `
  --axis-point 0 4829 4118 --axis-dir 1 0 0 `
  --pair-gate 3.5 --skin 4.0 --max-concurrent 8
& '..\_sf_baseline_e0cf\build\Release\scroll_whole.exe' "$root\baseline" --reregister
& '..\_sf_baseline_e0cf\build\Release\scroll_whole.exe' "$root\baseline" --audit

& build\Release\scroll_whole.exe `
  "$root\mesh\dump" "$root\candidate" `
  --axis-point 0 4829 4118 --axis-dir 1 0 0 `
  --pair-gate 3.5 --skin 4.0 --max-concurrent 8
& build\Release\scroll_whole.exe "$root\candidate" --reregister
& build\Release\scroll_whole.exe "$root\candidate" --audit

& build\Release\scroll_whole.exe `
  "$root\mesh\dump" "$root\candidate_repeat" `
  --axis-point 0 4829 4118 --axis-dir 1 0 0 `
  --pair-gate 3.5 --skin 4.0 --max-concurrent 8
& build\Release\scroll_whole.exe "$root\candidate_repeat" --reregister
& build\Release\scroll_whole.exe "$root\candidate_repeat" --audit
```

Audit exit 2 is a scientific gate failure, not permission to change or rerun
the arm. Run the frozen comparator once:

```powershell
& $py python\scripts\compare_winding_scale_stress.py `
  --baseline-audit "$root\baseline\audit.json" `
  --candidate-audit "$root\candidate\audit.json" `
  --candidate-repeat-audit "$root\candidate_repeat\audit.json" `
  --baseline-index "$root\baseline\placed_index.json" `
  --candidate-index "$root\candidate\placed_index.json" `
  --candidate-repeat-index "$root\candidate_repeat\placed_index.json" `
  --candidate-dir "$root\candidate" `
  --candidate-repeat-dir "$root\candidate_repeat" `
  --core-bbox 12672 13056 4352 5248 3712 4608 `
  --known-bbox 12800 12928 4480 5120 3840 4480 `
  --out-json "$root\paired_scale_result.json" `
  --figure-prefix "$root\complete_scale_seam_map"
```

## Frozen strata, gates, and decision

Every audited face-adjacent pair set is assigned mechanically:

- `known_center`: both endpoint cubes are in the prior 25-cube box;
- `bridge`: exactly one endpoint is in the prior box;
- `new_only`: neither endpoint is in the prior box.

All strata are reported and rendered. Only `new_only` can satisfy the primary
effect gate. The minimum inventory is 100 placed cubes overall, 1,000 audited
correspondences overall, and at least 20 pair sets / 500 correspondences in
`new_only`.

Comparability requires identical initial calibration, cube IDs and geometry,
pair keys, and per-pair correspondence counts; all placed inputs OK; zero
skipped and low-confidence inputs; all metrics finite; and the complete known
25-cube box present. Candidate repeat `audit.json`, `placed_index.json`, and
every registered `_skin.f32` must be byte-identical, with the skin count equal
to the placed cube count.

`SCALE_STRESS_PASS` requires every structural gate plus all of:

1. built-in candidate audit PASS and overall candidate turn-off `<=5%`;
2. candidate `new_only` turn-off `<=5%`;
3. `new_only` turn-off at least **1.00 absolute percentage point** below the
   parent on identical pairs;
4. overall and `new_only` turn-off may not worsen by more than 0.50 point;
5. overall candidate `|du|<2` may be at most 0.50 point below parent; and
6. overall and `new_only` collision counts may not increase.

Passing structural and safety gates but missing only the 1.00-point effect is
`SCALE_STRESS_COMPATIBILITY_ONLY`. A turn-off safety, `|du|<2`, or collision
failure after valid structural comparison is `SCALE_STRESS_REGRESSION`.
Missing inventory, comparability, repeat, minimum-pair, or candidate audit
requirements is `SCALE_STRESS_INCONCLUSIVE`. No bridge, known-center, pitch,
sign, energy, visual, `|du|<6`, or downstream metric can rescue a failed
primary gate.

The complete seam figure contains every pair set in x-y, x-z, and y-z
projections for parent, candidate, and delta, using fixed 0-20% absolute and
-20 to +20 percentage-point delta color ranges. The prior box is outlined.
No favorable edge, layer, crop, or orientation may replace this figure.

## Frozen complete-RAW downstream execution

If input, mesh, and pair comparability are valid, run this downstream section
regardless of whether registration improves, is null, or regresses. All arms
use the complete 405-cube halo directory and the exact parameters from the
earlier downstream protocol.

```powershell
$env:OMP_NUM_THREADS = '1'
$down = "$root\downstream"
$common = @(
  '--raw', "$root\halo_grid\cubes_RAW",
  '--steps', '12345', '--du', '1', '--dv', '1',
  '--band-cols', '32768', '--raw-chunk', '128',
  '--range', '2', '--nsteps', '5', '--threads', '1',
  '--snap-recto-iters', '0', '--snap-recto-range', '3',
  '--stretch-ratio', '4', '--stretch-floor', '25',
  '--max-edge', '0', '--synth-max-edge', '6',
  '--strip-w', '4096', '--dark-thresh', '40',
  '--seam-top', '0', '--no-xyzmap'
)

& build\Release\scroll_unroll.exe "$root\baseline" "$down\parent" `
  --id parent @common
& build\Release\scroll_unroll.exe "$root\candidate" "$down\candidate" `
  --id candidate @common
& build\Release\scroll_unroll.exe "$root\candidate" "$down\candidate_repeat" `
  --id candidate_repeat @common

& $py python\scripts\compare_unroll_outputs.py `
  --parent-dir "$down\parent" --parent-id parent `
  --candidate-dir "$down\candidate" --candidate-id candidate `
  --candidate-repeat-dir "$down\candidate_repeat" `
  --candidate-repeat-id candidate_repeat `
  --out-json "$down\paired_downstream_result.json" `
  --figure-prefix "$down\complete_unroll_strips"
```

Every arm must report zero missing RAW cubes during prewarm and stage 1. If
the one-cube halo is unexpectedly insufficient, that is disclosed as a failed
complete-RAW gate; the halo is not enlarged after seeing the log.

The existing downstream gates remain unchanged: exact arm comparability and
repeat artifacts; at most -0.20 point stage-1/final fill loss; at most +0.20
point final multi-cover; at most +0.05 final seam excess; stage-1/final canvas
area ratio in `[0.95,1.05]`; and efficacy only when stage-1 seam excess falls
by at least 0.05 without increased raw seam-column discontinuity. Its decision
is reported separately and cannot change the registration decision.

The combined complete stage-1/final figure and every original full-resolution
readable-strip PNG are retained and published. No ink, text, or legibility
claim is permitted without separate ground truth.
