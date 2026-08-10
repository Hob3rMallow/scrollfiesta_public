# PHerc0211 held-out placement calibration validation

Status: **frozen before access to any voxel, cube, mesh, or placement outcome
from the selected slab** at 2026-08-10T17:44:36+01:00
(Africa/Casablanca).

This protocol tests the multi-seed winding calibration and consensus component
gauge in commit `19152e49acf855d2bbe215329050ab0cc0085b5e`. It was written
after development on PHerc0139, PHerc0332, and PHerc1203, but before opening
the PHerc0211 slab fixed below. No validation threshold, axis, box, binary, or
algorithm may change after selected-slab access.

## Prior evidence and fixed question

The complete development evidence is disclosed rather than treated as held
out. With the candidate defaults, the audited turn-off fractions changed as
follows on already observed meshes:

| Development set | Frozen former default | Candidate default |
|---|---:|---:|
| PHerc0139 original 100 cubes | 2.903% | 2.593% |
| PHerc0139 adjacent 25 cubes | 3.418% | 1.172% |
| PHerc0332 available 81 cubes | 5.123% | 3.796% |
| PHerc1203 physical slab, 25 cubes | 9.858% | 3.312% |

Those comparisons diagnosed two order-sensitive behaviors: the former
calibrator accepted the first supported seed, and raw component gauging left
locally unanimous integer-radius corrections unapplied. They are development
results, not independent evidence.

The fixed validation question is: on one previously uninspected PHerc0211
placement slab, does the candidate default reduce cross-cube turn-off errors
against the frozen former default, while preserving join completeness and
using the identical cube meshes and audit pairs?

## Data provenance and metadata-only slab selection

- RAW level-0 CT:
  `s3://vesuvius-challenge-open-data/PHerc0211/volumes/20250821151803-9.362um-1.2m-113keV-masked.zarr`
- m7 level-0 surface prediction:
  `s3://vesuvius-challenge-open-data/PHerc0211/representations/predictions/surfaces/20250821151803-surface-20260413222639-surface-m7-L0-th0.2.zarr`
- Both level-0 arrays have shape `[19416,7948,7948]`, dtype `uint8`, C
  order, and `/` dimension separators. RAW chunks are `[128,128,128]` and
  uncompressed. Prediction chunks are `[192,192,192]` and Blosc/Zstd
  compressed.
- Exact remote `.zarray` byte hashes are RAW
  `c754c949661243c8e60a4b3212d924882908f480c8249d8f9557eb83c84856c4`
  (239 bytes) and prediction
  `e94ae5d7a63dfbf5923a0030f3d2fe969383995acf02ffb60f73e72959e6319f`
  (402 bytes).
- The public manual umbilicus has 87 control points. The exact original-byte
  local copy is
  `automatic_umbilicus/inputs/PHerc0211_umbilicus.original-bytes.json`,
  7,551 bytes, SHA-256
  `aa07c0cba4458def8092fb54e85299558a038ced408851366463f3a38655595c`.

The selection rule uses only array metadata and that manual curve. Among
interior control points with `4096 <= z <= 16000`, compute centered transverse
drift

`128 * hypot(x[i+1]-x[i-1], y[i+1]-y[i-1]) / (z[i+1]-z[i-1])`.

Require positive centered `dz`; choose the minimum drift, breaking an exact
tie by lower `z` and then original file order. The winner is curve index 30,
`(x,y,z)=(3527,3905,6007)`, centered `dz=475`, drift
`4.2179303323 voxels per 128 z voxels`.

The z origin is the 128-aligned interval containing the winning point. For y
and x, choose the 128-aligned five-cube start whose central cube center is
nearest the winning coordinate, with a lower-coordinate tie break. This fixes:

- box z/y/x: `[5888,6016) x [3584,4224) x [3200,3840)`;
- grid: `1 x 5 x 5 = 25` cubes of `128^3` voxels;
- axis point z/y/x: `[0,3905,3527]`;
- axis direction z/y/x: `[1,0,0]`.

This is a held-out **slab**, not a claim that no PHerc0211 voxel has ever been
seen. An earlier model-artifact reproduction used box
`[9708,9964) x [3974,4230) x [3974,4230)`, which is disjoint from this box in
z and x. Before this freeze the workspace contained no local PHerc0211 Zarr,
no `z05888_y*.tif`, and no PHerc0211 mesh or placement result. A Lasagna object
inventory and acquisition plan existed, but its guarded runner records that
no PHerc0211 dataset was materialized. If later evidence contradicts this
access audit, the result must be relabeled as a cross-scroll replication, not
an untouched-slab validation.

## Frozen source, build, and tests

Candidate source:

- commit: `19152e49acf855d2bbe215329050ab0cc0085b5e`;
- tree: `168ee535c6353ddcab6924c2d3a52ac46e5ba5fd`;
- parent: `b45504664df1a027dd2d7aeb9d61c60dc04c19bf`;
- `src/tools/scroll_whole.c`: 96,574 bytes,
  `266039d3ef5ce8828177416e9479d8e7d9af2d273fd2b6130a72c03f6e9f0996`;
- `src/whole/group_graph.c`: 54,100 bytes,
  `f98383f2bb28fefe87671c052cb03c9563e410739ca0b27500ba5d5e243971ce`;
- `src/whole/group_graph.h`: 8,517 bytes,
  `489ca48e689cd140442653ff5545cab5777a126bffde9e274d813185b5eb9f0a`.

Frozen former-default comparator:

- source commit `d9c67f137aeb16dd6e3d4b2ef6149079ffcdfe42`, tree
  `1265f49b0602febb2b958e4088fc811d508a358b`;
- `scroll_whole.exe`: 264,192 bytes,
  `aec8b814d6c0bfa59d29dccefea5389ce562172fa3569e014629335dd1be5a0f`.

Release/x64 candidate binaries rebuilt after the candidate commit:

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| `scroll_whole.exe` | 270336 | `ed78f26f8220f7199f35cd8a29fe8bb1618f72f176392dfb638cc0f110eb3715` |
| `grid_pipeline.exe` | 414720 | `98b9e8e6a6de9083d47c3527b79d77f6338afa212f40bbb37ab0c0c929876eb6` |
| `cube_mesh.exe` | 1173504 | `c08437db47a73158f5d1ca7abded140ab7333e9a304ca24d19ea5818829ea72e` |

The unchanged carver is 6,049 bytes, SHA-256
`19e98e83d3b54259f23ac22bbf8286fe7cfdab85fbbb442528182870bf394811`.
The frozen Python environment is CPython 3.13.5 with NumPy 2.4.6,
tifffile 2026.6.1, Zarr 3.2.1, and fsspec 2026.4.0.

Post-build `scroll_whole --selftest` and `grid_pipeline --selftest` both pass
with zero failures. The former includes synthetic multi-seed selection and
unanimous/disputed component-gauge cases. The last full native aggregate run
was 8/9 harness groups: its sole failure is the pre-existing
`pipeline_cube_smoke_test` hard-coded external TIFF path; the focused common
suite is 24/24. This is disclosed and is not called a fully green aggregate
suite.

## Frozen execution

Use output root
`output/experiments/pherc0211_multiseed_holdout_20260810`. The output root
must not exist at start.

1. Run the unmodified carver once with this exact semantic command:

   ```powershell
   & '..\_sf_sparse\python\.venv\Scripts\python.exe' python\scripts\carve_grid_tifs.py `
     --pred-zarr 's3://vesuvius-challenge-open-data/PHerc0211/representations/predictions/surfaces/20250821151803-surface-20260413222639-surface-m7-L0-th0.2.zarr' `
     --raw-zarr 's3://vesuvius-challenge-open-data/PHerc0211/volumes/20250821151803-9.362um-1.2m-113keV-masked.zarr' `
     --bbox 5888 6016 3584 4224 3200 3840 `
     --umbilicus 3905 3527 `
     --out output\experiments\pherc0211_multiseed_holdout_20260810\grid
   ```

2. Validate and hash every carved file before meshing. RAW must contain
   exactly 25 uint8 `128^3` TIFFs and may not be globally all zero. Prediction
   must contain at least 20 nonempty uint8 `128^3` TIFFs, and every value must
   be exactly 0 or 255. Missing prediction files are permitted only when the
   carver reports the corresponding cube as all zero. Source shape/dtype and
   the bbox/umbilicus manifest must match this document.
3. Mesh the fixed present prediction inventory once:

   ```powershell
   & build\Release\grid_pipeline.exe `
     output\experiments\pherc0211_multiseed_holdout_20260810\grid `
     output\experiments\pherc0211_multiseed_holdout_20260810\mesh `
     --halo 13 --threads-per-cube 1 --max-concurrent 8 `
     --simplify cvt --trim-inset 0 --skip-weld `
     --exe build\Release\cube_mesh.exe
   ```

   Every present prediction cube must exit 0 and produce a complete
   `step12_final` OBJ. No outcome-dependent retry, alternative simplifier,
   changed concurrency, reject override, or alternate slab is allowed.
4. Feed the identical dump to the frozen former binary and then the candidate
   binary, in that order, using no pitch or gauge override:

   ```powershell
   & '..\_sf_crossscroll1203\build\Release\scroll_whole.exe' `
     output\experiments\pherc0211_multiseed_holdout_20260810\mesh\dump `
     output\experiments\pherc0211_multiseed_holdout_20260810\baseline `
     --axis-point 0 3905 3527 --axis-dir 1 0 0 `
     --pair-gate 3.5 --skin 4.0 --max-concurrent 8
   & '..\_sf_crossscroll1203\build\Release\scroll_whole.exe' `
     output\experiments\pherc0211_multiseed_holdout_20260810\baseline --reregister
   & '..\_sf_crossscroll1203\build\Release\scroll_whole.exe' `
     output\experiments\pherc0211_multiseed_holdout_20260810\baseline --audit

   & build\Release\scroll_whole.exe `
     output\experiments\pherc0211_multiseed_holdout_20260810\mesh\dump `
     output\experiments\pherc0211_multiseed_holdout_20260810\candidate `
     --axis-point 0 3905 3527 --axis-dir 1 0 0 `
     --pair-gate 3.5 --skin 4.0 --max-concurrent 8
   & build\Release\scroll_whole.exe `
     output\experiments\pherc0211_multiseed_holdout_20260810\candidate --reregister
   & build\Release\scroll_whole.exe `
     output\experiments\pherc0211_multiseed_holdout_20260810\candidate --audit
   ```

All stdout/stderr, `placed_index.json`, `audit.json`, the input inventory, and
SHA-256 manifests are retained. Audit exit 2 is a scientific gate failure,
not a reason to rerun.

## Frozen gates and interpretation

Input/mesh/comparability gates all must pass:

- at least 20 prediction cubes, exactly 25 RAW cubes, and all carved-type and
  value checks above;
- every present prediction cube meshes with exit 0; no rejected, failed, or
  truncated final mesh;
- baseline and candidate see exactly the same cube IDs, vertex/face counts,
  adjacency pair-sets, total audit-pair count, and at least 100 audit pairs;
- both placement passes finish with zero failed/skipped cubes and, after
  reregistration, zero low-confidence accepted inputs;
- candidate calibration has at least 10 supported seed estimates, at least
  60% of supported estimates on its selected winding sign, and selects a real
  evaluated candidate nearest the within-sign median; and
- all required JSON files parse and all reported metrics are finite.

The candidate is an **independent efficacy pass** only if all comparability
gates pass and all of the following hold:

1. candidate built-in audit status is PASS: turn-off <=5%, `|du|<2` >=50%,
   and at least 100 pairs;
2. candidate turn-off fraction is at least **1.00 absolute percentage point**
   lower than baseline on the identical pairs; and
3. candidate `|du|<2` completeness is no more than 2.00 percentage points
   below baseline.

If the candidate passes its audit but improvement is below 1.00 point, the
result is only held-out compatibility/no-regression evidence. If the candidate
turn-off fraction is more than 0.50 point worse, drops `|du|<2` by more than
2.00 points, or fails while baseline passes, record a regression. Every other
outcome is inconclusive. Component radius/raw gauge counts, pitch dispersion,
collisions, and `|du|<6` are mandatory descriptive outputs but cannot rescue a
failed primary gate.

If any input or mesh gate fails, stop without opening another slab. If a code,
binary, axis, or threshold defect is found after access, publish the failure
and write a new protocol for a different future scroll; do not amend this run
into a pass. Publish the complete positive, negative, or inconclusive result
and clearly distinguish the four development sets from this held-out slab.
