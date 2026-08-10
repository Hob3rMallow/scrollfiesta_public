# PHerc1203 physical-truth snap replication preregistration

Status: **frozen before any PHerc1203 voxel, label chunk, mesh, or local
outcome was read** at 2026-08-10T15:38:13+01:00 (Africa/Casablanca).

Pre-access amendment 1: **frozen at 2026-08-10T15:54:31+01:00, still before
any PHerc1203 voxel, label chunk, mesh, or local outcome was read**. A
skeptical review identified five material weaknesses in the first scorer:
the lack of the reference evaluator's shifted null, arm-dependent attrition,
an IID bootstrap on one contiguous slab, an effect floor below registration
error, and no separation of `boundary_poor` tissue. The amended metric below
supersedes those parts of the first freeze. Git history preserves the original
protocol and scorer; no outcome was available when this correction was made.

Before this freeze, only public catalog/Zarr metadata, the physical-audit
README and source code, and its already published whole-scroll aggregate
numbers had been inspected. The `labels1203_L1.tar` asset had not been
downloaded or opened. No result localized to the box below was known.

This is an independent-scroll, independent-physical-reference test of the
default decision in ScrollFiesta PR #11. It explicitly extends, and does not
duplicate or claim authorship of, the PHerc1203 physical labels and binary
volume evaluator published by `7jycwjmbfn-eng/pherc0139-physical-audit` and
proposed in Villa PR #1382. The new question is narrower: which of two
same-topology ScrollFiesta mesh arms is closer to the published physical
recto reference?

## Fixed question and interpretation

On one untouched, metadata-selected PHerc1203 slab, is the PR #11
repair-only default (`snap_recto_iters=0`) closer to an independently derived
physical recto band than the former four-iteration behavior
(`snap_recto_iters=4`, range 3)?

The physical score is the primary decision metric. The already published
CT-ridge score is secondary and mechanistic. If they disagree, the physical
reference controls any practical default claim, and the disagreement must be
reported rather than averaged away.

The physical labels are still a registered, algorithmically derived
reference rather than infallible manual truth. Their published held-out
registration median is 2.38 um. Claims must retain that limitation.

## Frozen code and data

- Candidate implementation: ScrollFiesta PR #11 exact public head
  `4777630e5c81111684d156ef4ffdb5964cedc57d`.
- The original pre-access protocol is public at `7726035`; the first scorer is
  public at `b9c495a`. This amendment is the commit containing this paragraph;
  its exact SHA and the combined-build hashes will be recorded in the build
  receipt before any data access.
- Reliability patch used for this run: ScrollFiesta PR #12 exact public head
  `f0d9d2e54823e7ba2460725e81290eead8ed6e5e`. The experiment build is the
  PR #12 patch cherry-picked onto PR #11; its resulting commit/tree and
  binary hashes will be recorded before data access. No geometry or scoring
  behavior may otherwise change.
- RAW CT, level 0:
  `s3://vesuvius-challenge-open-data/PHerc1203/volumes/20250820131727-9.362um-1.2m-113keV-masked.zarr`
- m7 surface prediction, level 0:
  `s3://vesuvius-challenge-open-data/PHerc1203/representations/predictions/surfaces/20250820131727-surface-20260413222639-surface-m7-L0-th0.2.zarr`
- Both level-0 arrays have metadata shape `[18977,6844,6844]`; RAW chunks
  are `[128,128,128]` and prediction chunks are `[192,192,192]`.
- Physical reference: release `v1.0` asset `labels1203_L1.tar`, 515,379,200
  bytes, SHA-256
  `32a09f6081342b0f015b258ec577d0296ff23a55892af9785689d8a55bff344c`.
- The release tag commit predates the attached PHerc1203 asset and is not its
  code provenance. The PHerc1203 label generator first appears at exact
  commit `5a86b43743adbd816cb53d115be105c6b2b81e5e`; the immutable asset digest
  above is the data authority.
- The contained uint8 label array has L1 origin `[3936,0,0]`, shape
  `[2016,3456,3456]`, and flags `valid=1`, `material=2`, `centerline=4`,
  `recto_band=8`, `boundary_poor=16`.

The fixed L0 box is:

`z [9856,9984), y [3072,3712), x [3072,3712)`

This is `1 x 5 x 5 = 25` cubes of `128^3` voxels. The selection rule uses
metadata only: choose the 128-aligned one-slab z interval nearest the center
of the physical-label z overlap, and the 128-aligned 640-by-640 y/x footprint
nearest the center of the level-0 array. Ties choose the lower coordinate.
The fixed straight axis is the level-0 array center,
`axis_point_zyx=[0,3421.5,3421.5]`, `axis_dir_zyx=[1,0,0]`. The box and axis
may not move after any data value is read.

## Frozen pipeline and arms

1. Build the exact combined source above in Release/x64, run
   `scroll_unroll --selftest` and the relevant common/Python tests, and record
   source, scorer, test, and binary hashes before data access.
2. Carve exactly the fixed RAW and nonzero-to-255 prediction cubes with the
   unmodified public carver. Record source metadata and output manifests.
3. Run `grid_pipeline` with documented defaults, halo 13, one thread per
   cube, and no outcome-dependent retry settings.
4. Run `scroll_whole`, then `--reregister --audit`, with the fixed axis. No
   axis, registration, audit, or placement parameter may be tuned on this
   box.
5. From one common placed input export a pre-snap mesh with `--steps 12` and
   exactly two stage-4 arms, all other settings fixed:
   - candidate: `--steps 124 --snap-recto-iters 0 --snap-recto-range 3`
   - former behavior: `--steps 124 --snap-recto-iters 4 --snap-recto-range 3`
6. Rebuild once and export the candidate without a `--snap-recto-iters`
   flag. It must be byte-identical to explicit iteration 0.

## Amended frozen primary physical metric

The scorer must pass synthetic coordinate, topology, eligibility,
shifted-null, area-coverage, registration-sensitivity, spatial-bootstrap, and
end-to-end decision/JSON tests before the label asset or fixed-box chunks are
opened.

- OBJ coordinates are z/y/x level-0 index coordinates. Convert them to the
  physical-label local L1 coordinates by
  `p_label = (p_L0 - 0.5) / 2 - [3936,0,0]`. This is the voxel-center mapping
  consistent with the reference evaluator's pairwise L0-to-L1 max pooling.
- Read only the fixed label-local window `z [976,1072), y [1408,2048),
  x [1408,2048)`. The expanded y/x halo keeps both the real and shifted-null
  samples at least 64 L1 voxels from a crop edge. Record a SHA-256 of the
  loaded window in addition to the archive and Zarr metadata hashes.
- For every label z plane, compute the two-dimensional Euclidean distance
  transform to bit 8 (`recto_band`). Sample the resulting distance stack
  trilinearly at each mesh vertex. Per-plane 2-D distances match the
  reference's axial side instrument and prevent a 3-D shortcut to a surface
  in another z plane.
- Use the physical evaluator's fixed null: shift every mesh sample by +64 L1
  voxels in y. For each cube, compute the real effect and the shifted-null
  effect as `distance(iter4) - distance(iter0)`. The primary effect is their
  difference, `real_effect - null_effect`; positive favors iteration 0.
- Define eligibility from the pre-snap mesh only. A pre-reference vertex must
  be face-used, have a finite positive pre-snap area weight, lie at least four
  L0 voxels inside the fixed box, and have valid real and +64-y null label
  support. A target vertex must additionally remain within 3.0 L1 voxels of a
  real recto band under every frozen registration offset below.
- `boundary_poor` is also frozen from pre only: if any of the eight pre-snap
  label corners carries bit 16 under any registration offset, the vertex is
  excluded from the primary resolved-boundary score and reported in a
  separate descriptive stratum. Its predeclared, scored, and arm-support-
  failure counts are reported explicitly. It cannot be used to select the box
  or tune a gate.
- An arm may not remove its own failures. Candidate and former-arm box,
  validity, finiteness, null-support, and crop-halo checks are safety gates on
  the complete predeclared primary set. Any failure kills the metric decision;
  it never drops the affected vertex and continues as if it had not existed.
- Give each vertex one third of every incident pre-snap triangle area. For
  each provenance cube and arm, report fixed-area-weighted median real, null,
  and null-corrected effects, plus unweighted sensitivity results.
- Registration sensitivity uses the published 6.1 um p95 error, or
  `6.1 / 18.724 = 0.326` L1 voxel. Evaluate the zero offset and all 26
  nonzero directions in `{-1,0,1}^3`, normalizing every nonzero direction to
  radius 0.326. The eligible set stays fixed; no offset is selected from the
  result.
- This is one contiguous 5x5 slab, not 25 IID samples. Report every cube and
  all five y-row and five x-column medians. Use 10,000 resamples of the five
  row medians (seed `20260810`) and separately the five column medians (seed
  `20260811`). These spatial-cluster intervals are consistency checks, not a
  population-level cross-scroll confidence claim.

## Amended frozen gates and decision rule

Input, coverage, and safety gates:

- exactly 25 nonzero uint8 RAW cubes and 25 prediction cube files, with no
  grid holes;
- at least 20 cubes placed, with zero failed, skipped, or low-confidence
  accepted inputs; audit turn-off pairs <=5% and `|du|<2` completeness >=75%;
- identical vertex count, face indices, UVs, cube ranges and finite
  coordinates across pre and both arms;
- zero arm-support failures on the complete predeclared primary vertex set,
  and every sampled real/null EDT distance <64 L1 voxels;
- at least 20 scored cubes, each with at least 100 primary vertices and at
  least 20% of its pre-reference surface area retained; at least 40% retained
  area globally; and at least three scored cubes in every y row and x column;
- explicit iteration 0 and the no-flag default are byte-identical; and
- two scorer invocations emit byte-identical JSON.

The candidate's **physical metric gate passes, with external pipeline gates
still required**, only if every item above passes and all of these hold:

1. central median null-corrected cube effect >=0.33 L1 voxel (6.18 um,
   above the published 6.1 um registration p95);
2. the central median uncorrected real effect is positive;
3. both the row-cluster and column-cluster bootstrap 95% lower bounds are >0;
4. at least 80% of scored cubes have positive null-corrected effect, every
   y-row and x-column median is positive, and no more than 10% of cubes favor
   iter4 by more than 0.25 L1 voxel; and
5. the overall median real and null-corrected effects remain positive under
   all 27 frozen registration offsets.

The former four-iteration metric gate passes only if the same five gates pass
after reversing every sign and interval bound. Otherwise the result is a
tie/inconclusive default decision. A full claim is permitted only after the
external placement/audit/default-equivalence gates also pass. This is
independent-scroll, one-slab evidence, not a population-wide cross-scroll
inference. No threshold, box, eligibility rule, offset, or arm may change
after access.

Run the published sigma-1, fixed-pre-normal CT-ridge scorer as a secondary
mechanistic metric on the same 25 cubes. It cannot override the null-corrected
physical decision.

## Frozen visual and publication rule

The only qualitative plane is global L0 z=9920, fixed as the midpoint of the
slab. Assemble all 25 RAW tiles and overlay the physical recto band plus the
pre, iter0, and iter4 mesh intersections. Also render the complete fixed 5x5
per-cube null-corrected resolved-boundary effect heatmap. No slice or crop
search is allowed.

Publish the complete positive, negative, or inconclusive result with credit
to the physical-audit/#1382 authors, exact commands, hashes, per-cube data,
and limitations. If the physical result favors iter4, amend the PR #11
default claim rather than hiding the contradiction. If an input/safety gate
fails, publish only a clearly labeled feasibility failure; it supports no
physical arm claim.

## Pre-access build and test receipt

Frozen at 2026-08-10T16:11:23+01:00, still with
`labels1203_L1.tar` absent and without reading any PHerc1203 voxel, label
chunk, mesh, or local outcome.

- Public amended protocol/scorer commit:
  `120690c5b3d34bc9786505811684b092fdc18d79`.
- Exact combined experiment commit:
  `d9c67f137aeb16dd6e3d4b2ef6149079ffcdfe42`, tree
  `1265f49b0602febb2b958e4088fc811d508a358b`. Its PR #12 cherry-pick
  `2b4de7632303d0baca6df943ad071a304a6151b8` has the same stable patch ID
  (`98c41e386847549c05da95d2045fb88b33df907b`) as public PR #12 head
  `f0d9d2e54823e7ba2460725e81290eead8ed6e5e`.
- MSBuild `17.10.4.21802`, Release/x64, full `scrollfiesta.sln` build: exit 0.
- `scroll_unroll --selftest`: PASS, zero failures.
- Full Python suite: 58 passed, 2 skipped. The focused native common suite is
  24/24, including `clipper2_union_2d_abi`.
- The inherited aggregate native runner is 8/9: its only failure occurs
  before pipeline execution because upstream `main` hard-codes another
  developer's `C:/Users/mordr/...` TIFF path in
  `scripts/extract/pipeline_cube_smoke_test.c:47`. This is a disclosed test-
  harness portability failure, not a passing pipeline test and not evidence
  for either arm.

Frozen SHA-256 and byte sizes:

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| `scroll_unroll.exe` | 772096 | `bd85a0bb7d7e2481b535e9293cc8e026cff3e371a7e4269435ee2ce3c5ec37b4` |
| `scroll_whole.exe` | 264192 | `aec8b814d6c0bfa59d29dccefea5389ce562172fa3569e014629335dd1be5a0f` |
| `grid_pipeline.exe` | 414720 | `670dcabe05e854e93b30a0ce97ceeea0e0bd4b4370aac7a69016d4c78e8819a3` |
| `cube_mesh.exe` | 1174016 | `53a8ac8ed7980d27d39a4ca1aabf76c8fefe47230bd2115ae117c614ba5066e7` |
| `grid_weld.exe` | 489472 | `8bae8ce7888696e3bdc8b32275467e1682f1370c3386c7c8835115623760c2bc` |
| `all_tests.exe` | 1331712 | `c43056e39ce93a6e8f9c01918f390ece11c8b731427d76af2d70c777bd79356e` |
| `Clipper2.lib` | 1165150 | `b2d6f07c03237f770ab72cebfbdd1036210a9f5345dddb0c086f17c1418fd818` |
| physical scorer | 36639 | `3c8b7f4ec724bfe9e64cd55d2cbfa31dfd297dfcbaf75d651a417fae6786b64e` |
| physical scorer tests | 10257 | `3b163cd07552f87911801c773c367ee17c109d88c1d58e42a965ae689b0e247e` |
| preregistration before this receipt | 12364 | `3327ce96a721ae643699dfa4ccd92ed7c24db9a236afbc1449ce164d0f8bfb10` |

## Outcome

The fixed run stopped at the external placement gate: 125/1268 turn-off
pairs (9.858%) exceeded the frozen 5% maximum, although all 25 cubes completed
and `|du|<2` completeness was 79.18%. No snap arm or physical arm score was
run. See
[`PHERC1203_PHYSICAL_SNAP_FEASIBILITY_RESULT.md`](PHERC1203_PHYSICAL_SNAP_FEASIBILITY_RESULT.md)
for the complete feasibility result and hashes.
