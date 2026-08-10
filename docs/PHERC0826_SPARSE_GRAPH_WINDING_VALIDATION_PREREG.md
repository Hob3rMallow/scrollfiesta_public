# PHerc0826 sparse-graph winding validation

Status: **frozen before access to any voxel, cube, mesh, or placement outcome
from the selected slab** at 2026-08-10T20:19:23+01:00
(Africa/Casablanca).

This protocol tests only the sparse graph-edge admission and graph-supported
polish guard in candidate commit
`a47a21c3a9d445c90d28e218b8f0fcbecd49c079`. Its exact parent
`e0cf51cab03f6fe84e43e8c7ff5be9051d951b1d` is the comparator. Both arms
therefore share the already developed multi-seed pitch selection and consensus
component gauge; the held-out question is whether preserving sparse, observed
cross-cube winding constraints improves registration on a separate scroll.

No axis, box, binary, threshold, visual range, or decision rule may change
after selected-slab voxel access.

## Development evidence and fixed mechanism

The earlier PHerc0211 held-out run improved under the parent implementation
but failed its frozen ceiling: 10.012% former-default turn-off became 6.876%,
still above 5%. That result was published as inconclusive before follow-up.
PHerc0211 then became development data.

Failure attribution found that all 54 residual parent errors lay in 38
observed group-pair buckets with only one or two correspondences. The old
three-pair edge floor discarded every one. Admitting one-pair observed edges
removed those errors, while unrestricted Gauss-Seidel polish could reintroduce
them. The candidate therefore:

1. admits an observed group-pair bucket from one correspondence;
2. records every group supported by an admitted graph edge; and
3. restores both integer winding and matching residual shift after polish for
   graph-supported groups, while leaving polish active for graph-isolated
   groups.

Across five now-development blocks (15,437 identical audit pairs), the parent
to candidate comparison was 501/15437 (3.245%) to 271/15437 (1.756%), a 45.9%
relative reduction; collisions fell 1250 to 721. This is mechanism-development
evidence, not the independent result registered here.

## Data and metadata-only selection

- RAW level-0 CT:
  `s3://vesuvius-challenge-open-data/PHerc0826/volumes/20250821151701-9.362um-1.2m-113keV-masked.zarr`
- m7 level-0 surface prediction:
  `s3://vesuvius-challenge-open-data/PHerc0826/representations/predictions/surfaces/20250821151701-surface-20260413222639-surface-m7-L0-th0.2.zarr`
- Both arrays: shape `[16920,8169,8169]`, dtype `uint8`, C order, `/`
  dimension separator. RAW chunks are `128^3`, uncompressed; prediction chunks
  are `192^3`, Blosc/Zstd.
- Exact remote `.zarray` bytes: RAW 239 bytes,
  SHA-256 `93860ad88aee4fd67cf98b819dfb6e9afffaaae116d07123ad095e1ef8559c96`;
  prediction 402 bytes,
  SHA-256 `a7fd8efb718d57000790fc393414cc51761a74d8e9dd5a180ab96b6e4686e549`.
- Public manual curve: 49 points; exact reversible cache 5,917 bytes,
  SHA-256 `b82b854dac4990b09dc980d0d0fe227c2cf6766541db985e173f4000a6770307`;
  decoded source 4,379 bytes,
  SHA-256 `ddf2ffa2ab91270b4ccc443d22f10587090f1c5b5561d34dfc4c838ed7a451f3`.

The complete conservative pre-access reconstruction is frozen in
`docs/results/pherc0826_prior_access_manifest.json`. It includes the earlier
regional reproduction, all 30 deterministic 64-cube diagnostic requests, and
all 72 possible deterministic 96-cube attempts from a conditional diagnostic.
Catalog, ROI arithmetic, Lasagna inventory, and curve reads were metadata-only.

Selection uses only that access manifest, array shape, and the manual curve.
For each interior curve point (`4096 <= z <= 16000`), form its containing
128-aligned slab and require a full 128-voxel z guard from every prior request.
Among eligible points, minimize centered transverse drift

`128 * hypot(x[i+1]-x[i-1], y[i+1]-y[i-1]) / (z[i+1]-z[i-1])`,

breaking ties by lower z then file order. Exactly one point is eligible:
zero-based index 35, `(x,y,z)=(4118,4829,12846)`, centered `dz=1044`, drift
`27.3082634844` transverse voxels per 128 z voxels. The nearest prior z request
remains 348 voxels away. Centering a 5 by 5 aligned transverse grid fixes:

- box z/y/x: `[12800,12928) x [4480,5120) x [3840,4480)`;
- grid: `1 x 5 x 5 = 25` cubes of `128^3` voxels;
- axis point z/y/x: `[0,4829,4118]`;
- axis direction z/y/x: `[1,0,0]`.

This is a held-out **slab**, not a claim that PHerc0826 has never been used.
The prior regional reproduction used
`[8460,8716) x [4084,4340) x [4084,4340)`. Two old collection-wide
diagnostics accessed sparse random blocks, which is why their full request
sequences are disclosed rather than omitted.

Jeff Chen's public `jeff-j-chen/vesuvius` commit
`3cb86a6bb1781354fa31c4395d277b3da6f4e2dd` documents a PHerc0826
VC3D-grown merged patch (`475 x 227` tifxyz) used for ink-model tests. Its raw
coordinates are not published, so external spatial non-overlap cannot be
proved. The repository has no declared license; neither code nor artifacts are
copied. The claim here is bounded algorithmic cross-cube registration, not an
untouched-scroll or first-surface claim.

## Frozen source, binaries, and overlap check

Candidate:

- commit `a47a21c3a9d445c90d28e218b8f0fcbecd49c079`, tree
  `0f2ed0cc9638c4d7021a79440c25e100684e4bb8`;
- `src/tools/scroll_whole.c`: 101,338 bytes,
  `4df29f7cc0b6b3d5a688f6dd2a4ff7d49eee29b84728f0a87fcc1b1b874c57a5`;
- `src/whole/group_graph.c`: 55,472 bytes,
  `a0ba27ba842407a73ad7e18ad7fa53ea2b1ac106498a09476f8a6863b29d4b55`;
- `src/whole/group_graph.h`: 8,516 bytes,
  `856628e174d5117dd8767676fa4b5d4eb4eb3c8d77787b6d4890c6d7a5475d00`;
- MSVC Release/x64 `scroll_whole.exe`: 273,408 bytes,
  `6624f7fe3ca062d10c03e7a799b1a32e275b8c2c61695a7f1f8360bc8539818e`.

Parent comparator:

- commit `e0cf51cab03f6fe84e43e8c7ff5be9051d951b1d`, tree
  `fcf7f285eaead55589918437007d5362c06b8b20`;
- clean detached build `scroll_whole.exe`: 270,336 bytes,
  `eb1445504d69375f36891081f39ee6aa7d65a8af22b46b164cc75bbabce2db8a`.

Shared meshing inputs:

- `grid_pipeline.exe`: 414,720 bytes,
  `98b9e8e6a6de9083d47c3527b79d77f6338afa212f40bbb37ab0c0c929876eb6`;
- `cube_mesh.exe`: 1,173,504 bytes,
  `c08437db47a73158f5d1ca7abded140ab7333e9a304ca24d19ea5818829ea72e`;
- carver: 6,049 bytes,
  `19e98e83d3b54259f23ac22bbf8286fe7cfdab85fbbb442528182870bf394811`.

Both exact `scroll_whole` binaries pass their complete in-process selftests.
The candidate includes the new sparse-edge and polish-guard cases. The last
full native aggregate remains 8/9 groups because of the pre-existing
hard-coded external TIFF path in `pipeline_cube_smoke_test`; the focused common
suite is 24/24. This is not described as a fully green aggregate suite.

The frozen paired verifier/renderer is
`python/scripts/compare_winding_audits.py`, 18,815 bytes,
SHA-256 `e8c9d492e4eb914d005614f9c50bc7f45089ba3cef45ccf01ba70d5f9b578def`.
Its three focused tests pass, including mismatch rejection and output hashing;
two complete renders of the same real development comparison were byte
identical in PNG and SVG.

At freeze time, upstream ScrollFiesta main was
`4f43cfc242d2d7f90147ceb5a38790fb36517dcd`; its only open PRs were our
unrelated #11 and #12. Current open Villa PRs and GitHub code search found no
second implementation of this sparse cross-cube constraint guard. GitHub
search is bounded and cannot see private or unindexed work.

The updated 15-channel Discord archive (3,523 messages through 2026-08-08)
likewise contains no matching implementation. It does contain the directly
relevant public observations that correct same/relative winding descriptors
are crucial and that a poor fit can look visually solid. This protocol
therefore requires both exact pair metrics and a complete, non-selected seam
visual; appearance alone cannot pass it. The dedicated ScrollFiesta channel is
not in that archive, so the overlap statement remains bounded.

## Frozen execution

Output root:
`output/experiments/pherc0826_sparse_graph_holdout_20260810`. It must not exist
at start.

1. Carve the fixed box once:

   ```powershell
   & '..\_sf_sparse\python\.venv\Scripts\python.exe' python\scripts\carve_grid_tifs.py `
     --pred-zarr 's3://vesuvius-challenge-open-data/PHerc0826/representations/predictions/surfaces/20250821151701-surface-20260413222639-surface-m7-L0-th0.2.zarr' `
     --raw-zarr 's3://vesuvius-challenge-open-data/PHerc0826/volumes/20250821151701-9.362um-1.2m-113keV-masked.zarr' `
     --bbox 12800 12928 4480 5120 3840 4480 `
     --umbilicus 4829 4118 `
     --out output\experiments\pherc0826_sparse_graph_holdout_20260810\grid
   ```

2. Hash and validate all carved files before meshing. RAW must contain exactly
   25 uint8 `128^3` TIFFs and may not be globally zero. Prediction must contain
   at least 20 nonempty uint8 `128^3` TIFFs with values only 0/255. A missing
   prediction tile is permitted only when the carver records that tile as all
   zero. Source metadata, box, and umbilicus must match this document.
3. Mesh the fixed present inventory once:

   ```powershell
   & build\Release\grid_pipeline.exe `
     output\experiments\pherc0826_sparse_graph_holdout_20260810\grid `
     output\experiments\pherc0826_sparse_graph_holdout_20260810\mesh `
     --halo 13 --threads-per-cube 1 --max-concurrent 8 `
     --simplify cvt --trim-inset 0 --skip-weld `
     --exe build\Release\cube_mesh.exe
   ```

   Every present prediction cube must exit 0 and produce a complete
   `step12_final` OBJ. No retry with changed simplifier, concurrency, reject
   override, or slab is allowed.
4. Run parent, then candidate, from the identical dump with no pitch, gauge, or
   graph override:

   ```powershell
   & '..\_sf_baseline_e0cf\build\Release\scroll_whole.exe' `
     output\experiments\pherc0826_sparse_graph_holdout_20260810\mesh\dump `
     output\experiments\pherc0826_sparse_graph_holdout_20260810\baseline `
     --axis-point 0 4829 4118 --axis-dir 1 0 0 `
     --pair-gate 3.5 --skin 4.0 --max-concurrent 8
   & '..\_sf_baseline_e0cf\build\Release\scroll_whole.exe' `
     output\experiments\pherc0826_sparse_graph_holdout_20260810\baseline --reregister
   & '..\_sf_baseline_e0cf\build\Release\scroll_whole.exe' `
     output\experiments\pherc0826_sparse_graph_holdout_20260810\baseline --audit

   & build\Release\scroll_whole.exe `
     output\experiments\pherc0826_sparse_graph_holdout_20260810\mesh\dump `
     output\experiments\pherc0826_sparse_graph_holdout_20260810\candidate `
     --axis-point 0 4829 4118 --axis-dir 1 0 0 `
     --pair-gate 3.5 --skin 4.0 --max-concurrent 8
   & build\Release\scroll_whole.exe `
     output\experiments\pherc0826_sparse_graph_holdout_20260810\candidate --reregister
   & build\Release\scroll_whole.exe `
     output\experiments\pherc0826_sparse_graph_holdout_20260810\candidate --audit
   ```

5. Run one fresh `candidate_repeat` from the same dump and command, then require
   byte-identical `audit.json`, `placed_index.json`, and every `_skin.f32` file.
6. Run the frozen comparator with its defaults and retain both complete seam
   maps:

   ```powershell
   & '..\_sf_sparse\python\.venv\Scripts\python.exe' `
     python\scripts\compare_winding_audits.py `
     --baseline-audit output\experiments\pherc0826_sparse_graph_holdout_20260810\baseline\audit.json `
     --candidate-audit output\experiments\pherc0826_sparse_graph_holdout_20260810\candidate\audit.json `
     --baseline-index output\experiments\pherc0826_sparse_graph_holdout_20260810\baseline\placed_index.json `
     --candidate-index output\experiments\pherc0826_sparse_graph_holdout_20260810\candidate\placed_index.json `
     --out-json output\experiments\pherc0826_sparse_graph_holdout_20260810\paired_result.json `
     --figure-prefix output\experiments\pherc0826_sparse_graph_holdout_20260810\all_pair_seam_map
   ```

All logs, source/input hashes, manifests, indices, audits, repeat comparison,
machine result, PNG, and SVG are retained. Audit exit 2 is a scientific gate
failure, not a reason to rerun.

## Frozen gates and interpretation

Input, mesh, and comparability gates all must pass:

- exactly 25 RAW cubes, at least 20 prediction cubes, and all fixed type/value
  checks above;
- every present prediction cube meshes successfully with no rejected, failed,
  or truncated final mesh;
- both arms share exact initial calibration fields, cube IDs, origins,
  vertex/face/skin/group counts, statuses, pair keys, and per-key pair counts;
- at least 100 audit pairs; both arms end with all inputs OK, zero skipped, and
  zero low-confidence inputs; and
- every parsed metric is finite.

The candidate is an **independent efficacy pass** only if all comparability
gates and all of these primary gates pass:

1. built-in candidate audit PASS: turn-off `<=5%`, `|du|<2 >=50%`, at least
   100 pairs;
2. candidate turn-off is at least **1.00 absolute percentage point** below the
   parent on identical pairs;
3. candidate `|du|<2` is no more than **0.50 percentage point** below parent;
4. candidate collision count does not increase;
5. the fresh candidate repeat is byte-identical for audit, index, and all skin
   sidecars; and
6. both complete seam figures render legibly, contain every audit pair set,
   and reconcile with the machine JSON. No selected-edge figure can satisfy
   this gate.

Pitch candidates, sign vote, graph node/edge/component counts, moves, energy,
per-edge changes, `|du|<6`, and collision delta are mandatory descriptive
outputs. Pitch/sign are not a between-arm efficacy gate because the exact
parent and candidate share initial calibration; the change under test occurs
only during reregistration.

If candidate passes audit and safety but improves by less than 1.00 point, the
result is held-out compatibility only. Candidate turn-off more than 0.50 point
worse, `|du|<2` below the 0.50-point safety floor, or increased collisions is a
regression. Other failures are inconclusive. No descriptive metric rescues a
failed primary gate.

Only a full independent efficacy pass authorizes a focused upstream
default-change PR. Any other result is still published exactly, with no
alternate slab, threshold, axis, or post-outcome gate change.
