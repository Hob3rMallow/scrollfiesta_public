# PHerc0139 adjacent-slab replication of ScrollFiesta PR #11

Status: **frozen before voxel/chunk access** at
2026-08-10T14:30:50+01:00 (Africa/Casablanca).

At freeze time the workspace contained no file whose cube id begins with
`z04864`; only Zarr metadata and the already published block ending at z=4864
had been inspected. This is a spatial replication on the same scroll, not a
cross-scroll claim.

## Fixed question

Does the PR #11 repair-only default (`recto_iters=0`) retain its lower
CT-ridge residual than the prior behavior (`recto_iters=4`, range 3) on the
first untouched 128-voxel slab immediately adjacent to the published
PHerc0139 block?

## Frozen data and selection rule

- Prediction, level 0:
  `s3://vesuvius-challenge-open-data/PHerc0139/representations/predictions/surfaces/20250728140407-surface-20260413222639-surface-m7-L0-th0.2.zarr`
- Raw CT, level 0:
  `s3://vesuvius-challenge-open-data/PHerc0139/volumes/20250728140407-9.362um-1.2m-113keV-masked.zarr`
- Bounding box, z/y/x:
  `z [4864,4992), y [3072,3712), x [2560,3200)`.
- Grid: `1 x 5 x 5 = 25` cubes of `128^3` voxels.
- The rule is the next 128-aligned z slab after the published
  `z [4352,4864)` block, with its unchanged y/x footprint. No voxel value or
  outcome was used to select it.
- Axis calibration remains the published PHerc0139 value:
  `axis_point_zyx=[0,3405,2878]`, `axis_dir_zyx=[1,0,0]`.

## Frozen implementation

- ScrollFiesta PR #11 exact source head:
  `b8803b901846eb49d5ecfa962735c3ea6a7214ce`.
- The exact-head binaries will be rebuilt from the clean `_sf_sparse`
  worktree before data access and their SHA-256 values recorded.
- Unmodified public carver SHA-256:
  `19E98E83D3B54259F23AC22BBF8286FE7CFDAB85FBBB442528182870BF394811`.
- The already tested grid-parameterized scorer used for the independent
  PHerc0332 protocol has SHA-256:
  `2F7AABAB3E83A53C06B9B3C9F2D7CD5654DF99A91C9C87A607E2D3BDB277B796`.
  Its numerical method is unchanged: sigma-1 CT, fixed pre-snap normals,
  +/-4 voxel profile, 0.25 step, parabolic peak refinement, cube as the
  inferential unit, seed 20260808 and 10,000 bootstrap resamples.

## Frozen execution

1. Carve the fixed RAW and thresholded prediction grid. RAW must contain
   exactly 25 uint8 `128^3` TIFFs. Record every source and output hash.
2. Run the exact-head cube pipeline with halo 13 and documented defaults.
   Do not tune from this slab.
3. Run `scroll_whole`, then `--reregister --audit`, using the fixed axis.
4. Input gate: at least 20 cubes must place; failed/skipped/low-confidence
   counts must be zero among accepted inputs; turn-off pairs must be <=5%;
   `|du|<2` join completeness must be >=75%. Failure makes the result
   inconclusive and no snap arm is scored.
5. From one common placed input export the shared-index pre-snap mesh and
   exactly two stage-4 arms:
   - candidate: `--steps 124 --snap-recto-iters 0 --snap-recto-range 3`
   - production comparator: `--steps 124 --snap-recto-iters 4
     --snap-recto-range 3`
6. Both arms must preserve identical vertex count, face indices, cube ranges,
   UV records and finite coordinates. Score only vertices eligible in every
   arm, with the pre-snap normal fixed for both.
7. Re-run the candidate with no `--snap-recto-iters` flag. It must be
   byte-identical to explicit iteration 0 before making a default-equivalence
   claim.

## Frozen decision gates

The replication passes only if all of the following hold:

1. Paired median cube improvement
   `abs_residual(iter4) - abs_residual(iter0) >= 0.20 voxel`.
2. The fixed 95% bootstrap interval for that paired median has lower bound
   greater than zero.
3. At least 80% of scored cubes have positive paired improvement.
4. All input, shared-index, topology/UV, finiteness and default-equivalence
   gates above pass.

If it passes, publish the complete result and a fixed midpoint cross-section
at global z=4928 as an additional spatial-replication amendment to PR #11.
If any gate fails, record the failure or inconclusive result without changing
the region, metric, thresholds or existing PR #11 claims. No new snap method
may be developed on this slab.

## Post-freeze pre-access build receipt

Before any new Zarr chunk was read, exact head `b8803b9` was rebuilt from the
clean `_sf_sparse` source with the repository's vendored dependency builder
and Visual Studio 2022 Release/x64 solution. The solution completed with zero
errors. `scroll_unroll --selftest` passed with zero failures. The generalized
carver/scorer focused suite passed 4/4 in the locked project environment.

Built-binary SHA-256 values:

- `scroll_unroll.exe`: `6AFDCDE2B8E3347EC708BA247DB8027BFAE8B5BA4C6B0E6BF8D04398331568CA`
- `scroll_whole.exe`: `8111CD6F51DBF0FD8B769D2FFCABFFF479DC31D7953BEA3E13D0BDBC656556B6`
- `grid_pipeline.exe`: `5725BD5F4FB7BBB1C388A79F9680BE3D4C6106654878783DA88FAFA35831981C`
- `cube_mesh.exe`: `66DA665E5927B092659F0230664C96099A47CFF154B9EF708299404E89F31EC8`
- `grid_weld.exe`: `443BE94EF371464F9A0452CF61A434D15234EA88E321CD2BAAF7F32F3A10644B`
