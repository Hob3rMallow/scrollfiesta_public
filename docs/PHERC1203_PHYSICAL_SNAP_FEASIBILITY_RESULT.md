# PHerc1203 physical snap replication: feasibility result

Status: **stopped at the preregistered external placement gate; no snap arm
was run and no physical arm claim is supported**.

The public protocol was frozen before data access in
[`PHERC1203_PHYSICAL_SNAP_REPLICATION_PREREG.md`](PHERC1203_PHYSICAL_SNAP_REPLICATION_PREREG.md).
The fixed 1x5x5 PHerc1203 slab was carved and meshed successfully, but its
post-reregistration placement audit exceeded the frozen turn-off limit:

| Gate | Frozen requirement | Observed | Result |
|---|---:|---:|---|
| cube pipeline | at least 20; zero failed/skipped | 25/25 exit 0; 0 failed/skipped | pass |
| turn-off pairs | at most 5% | 125/1268 = 9.858% | **fail** |
| join completeness | at least 75% with `|du|<2` | 79.18% | pass |

The initial placement completed 25/25 cubes and flagged five low-confidence
cubes; the required reregistration rewrote the final index with 25/25 status
0 entries and no failed finalization. The turn-off failure is independently
sufficient to stop the experiment. The worst concentration was the central
`y=3200`/`y=3328, x=3328` seam, with 58 turn-off pairs among 122.

This is useful negative feasibility evidence, not a comparison of the snap
methods. It shows that on this fixed compressed-scroll slab, the current
placement stage is not reliable enough for the physical labels to adjudicate
the PR #11 default. It does not weaken the two PHerc0139 results and it does
not establish cross-scroll generalization.

## Fixed inputs and provenance

- Box: L0 `z [9856,9984), y [3072,3712), x [3072,3712)`; 25 cubes.
- Axis: `point=[0,3421.5,3421.5]`, `dir=[1,0,0]`.
- Combined experiment source: `d9c67f137aeb16dd6e3d4b2ef6149079ffcdfe42`,
  tree `1265f49b0602febb2b958e4088fc811d508a358b`.
- Physical label archive: 515379200 bytes, SHA-256
  `32a09f6081342b0f015b258ec577d0296ff23a55892af9785689d8a55bff344c`.
- Fixed label-window SHA-256:
  `60986ca12c6282096a9a8e35dda9543b8155ab00fd453d2c5e1cfac5efaaf62c`.
- Public carver SHA-256:
  `19e98e83d3b54259f23ac22bbf8286fe7cfdab85fbbb442528182870bf394811`.
- Prediction TIFF manifest SHA-256:
  `560120bb7633b449c4b5ed0028c87ed00b5f58a3911dcd875572fd5bbf5be2d6`.
- RAW TIFF manifest SHA-256:
  `a2fbf03c8f5ab3ff815bba94bb7ad9f6a0131a203a8d76f15c062be39a54e737`.
  Each manifest is SHA-256 over sorted UTF-8 `filename sha256` lines with a
  final newline.

All 25 prediction TIFFs and all 25 RAW TIFFs are uint8 `128^3`; no cube is
all zero. Prediction values are strictly `{0,255}`.

## Exact commands

```powershell
python python/scripts/carve_grid_tifs.py `
  --pred-zarr s3://vesuvius-challenge-open-data/PHerc1203/representations/predictions/surfaces/20250820131727-surface-20260413222639-surface-m7-L0-th0.2.zarr `
  --raw-zarr s3://vesuvius-challenge-open-data/PHerc1203/volumes/20250820131727-9.362um-1.2m-113keV-masked.zarr `
  --bbox 9856 9984 3072 3712 3072 3712 `
  --umbilicus 3421.5 3421.5 --out data/pherc1203_physical_snap_fixed25 `
  --s3-anon yes

build/Release/grid_pipeline.exe data/pherc1203_physical_snap_fixed25 `
  output/experiments/pherc1203_physical_fixed25/run `
  --halo 13 --threads-per-cube 1 --max-concurrent 8

build/Release/scroll_whole.exe `
  output/experiments/pherc1203_physical_fixed25/run/dump `
  output/experiments/pherc1203_physical_fixed25/placed `
  --axis-point 0 3421.5 3421.5 --axis-dir 1 0 0

build/Release/scroll_whole.exe `
  output/experiments/pherc1203_physical_fixed25/placed --reregister --audit
```

The final command exits 1 because the audit gate fails; its output artifact is
valid and records `quality_gate.status="FAIL"`.

## Output hashes

- Pipeline summary: `e732cae0640db7cfdb4c1e9bce111e130cc711a2745d7554c9b059eac4c1e842`.
- Welded OBJ: `534c22fcb61fbba75ac8ffbd8a5fbfe0133dd5e6d5845458c7247e35a3908123`.
- Weld report: `0ad847bae16be6a09227cf2c271c6dd09a4b17c8f73fd42b075ef32ba5fbac95`.
- Final audit JSON: `da994f6db3357bc58d8fb12d5fa7733d2d33c5f28d3a678692b11d66738d1705`.
- Final placed index: `0f014e9cb10bad6fbfdebf662c202d304091ee94a0df2ebe2eefe1fe2b945792`.
- Pre-reregistration placed index: `90d69012853d032a907d99bea12a6efccb4b6dc121bf1a81c3e4c6ac36cab344`.

Machine-readable summary and every cube exit/timing row are in
[`results/pherc1203_physical_snap_feasibility.json`](results/pherc1203_physical_snap_feasibility.json)
and
[`results/pherc1203_physical_snap_pipeline_summary.csv`](results/pherc1203_physical_snap_pipeline_summary.csv).

The physical labels and evaluator are credited to
`7jycwjmbfn-eng/pherc0139-physical-audit` / Villa PR #1382. No pre-snap,
iter0, iter4, or no-flag mesh was exported; `score_snap_physical.py` was not
run on arm outcomes.
