# PHerc0826 complete downstream texture validation

Status: **frozen before any downstream `scroll_unroll` output was generated**
at 2026-08-10T21:51:55+01:00 (Africa/Casablanca).

The sparse-winding registration result is already known: candidate
`a47a21c3a9d445c90d28e218b8f0fcbecd49c079` reduced held-out whole-turn
errors from 44/1,076 to 7/1,076 on the fixed PHerc0826 slab. This follow-up is
therefore not a second independent test of that effect. It prospectively asks
whether the frozen registration difference survives the supported downstream
unroll and produces a measurable, completely rendered texture result.

Before this freeze, no PHerc0826 `scroll_unroll` raw texture, pipeline-stats
file, or downstream result was present anywhere in the workspace. The fixed
output root below did not exist. The parent and candidate registration audits,
placed UVs, and complete seam map had already been inspected and published.

## Why this test and why it is not duplicate work

The August prize criteria favor quantitative improvement to a core virtual-
unwrapping stage on real data, with visually verifiable output. PR #13 already
measures its registration stage directly, but its complete seam map is a
diagnostic rather than a flattened CT texture. This protocol closes that
specific evidence gap using ScrollFiesta's existing supported downstream
consumer; it does not invent another unroller, winding solver, or analytic
dashboard.

A fresh check of open ScrollFiesta and Villa PRs, indexed GitHub code, the
Discord archive through August 8, and the adjacent `abundantjoe/winding-sync`
repository found no paired parent/candidate full-strip validation of this
sparse cross-cube guard. `winding-sync` generates and globally reconciles CT-
derived winding constraints at a different layer. The earlier PHerc1203
physical-snap protocol already failed its frozen placement gate and will not be
repeated or presented as downstream evidence here.

## Frozen code and executable

- Focused upstream PR branch before this protocol: commit
  `94dc19337449422ebd4cf732b9aa0e9f50de827c`.
- Exact MSVC Release/x64 `scroll_unroll.exe`: 768,000 bytes, SHA-256
  `848969a62c9765e07629fe1bb1ad217606ba6ecab6e02dd963ea9928f8335844`.
- The exact binary's complete `--selftest` passed with zero failures.
- Frozen comparator `python/scripts/compare_unroll_outputs.py`: 15,659 bytes,
  SHA-256
  `9893c56edf0d2d967ab7cc578540e065c41880a1872a8b0caea9c482b9472bdb`.
- Comparator tests `python/tests/test_compare_unroll_outputs.py`: 6,322 bytes,
  SHA-256
  `4f97660a804ac9b4593f89a819d9fbd547540f48131a49a9971db5c3a6d4224f`.
- The focused non-network Python suite passed: 50 passed, one skipped, four
  network tests deselected.

The executable was built from the focused PR tree. Its unchanged TIFF/zlib
build inputs were copied byte-for-byte from the already audited Release/x64
dependency build; no ScrollFiesta source or generated placement artifact was
copied from another branch.

## Frozen real inputs

Shared source experiment:
`output/experiments/pherc0826_sparse_graph_holdout_20260810` under the public
development checkout.

- Exact parent audit SHA-256:
  `450965d08d1a0c7bc1bfb25011265f75291685be351fd4efa92f1c061b17f159`.
- Exact parent index SHA-256:
  `22bab257ff3cde175aa29c9004dddc5e9b0d3243f70afa1588964a361d748f06`.
- Candidate audit SHA-256:
  `c3900e6dbb501ab51665a0b05a754fd47c5bf5b673e74923be225a1d5e41653a`.
- Candidate index SHA-256:
  `405000fc1b99264c93a99413a7433fb1731fa8b3aef38b4dcebb0f28b64ca90f`.
- RAW/input validation SHA-256:
  `1453720bd2c78f94d7316cb3818539b9afacb7123197bfa6c93f629eaaba5d2a`.
- Parent complete placed-input manifest: 125 files; canonical
  `name<TAB>bytes<TAB>sha256` manifest SHA-256
  `95575fa8195d5d46810b99816af57a80b06d6d92446a5677765875a7c08295a2`.
- Candidate complete placed-input manifest: 125 files; same schema, SHA-256
  `973cde2138fd09561201508144924151e418074b438736479da1be2c3aa54e89`.

Both arms contain the identical 25 physical meshes and differ only in their
registered winding/UV sidecars. The fixed RAW inventory is the same 25-cube
PHerc0826 slab used by the preregistered registration comparison.

## Frozen execution

Run from the focused PR checkout. The environment and CLI both use one worker
thread. Each arm uses the same explicit parameters; the recto-boundary pass is
disabled in both arms so the known PR #11 target distinction cannot confound
this registration test.

```powershell
$env:OMP_NUM_THREADS = '1'
$source = '..\_sf_calibration\output\experiments\pherc0826_sparse_graph_holdout_20260810'
$out = '..\_sf_calibration\output\experiments\pherc0826_sparse_graph_downstream_20260810'
$common = @(
  '--raw', "$source\grid\cubes_RAW",
  '--steps', '12345', '--du', '1', '--dv', '1',
  '--band-cols', '32768', '--raw-chunk', '128',
  '--range', '2', '--nsteps', '5', '--threads', '1',
  '--snap-recto-iters', '0', '--snap-recto-range', '3',
  '--stretch-ratio', '4', '--stretch-floor', '25',
  '--max-edge', '0', '--synth-max-edge', '6',
  '--strip-w', '4096', '--dark-thresh', '40',
  '--seam-top', '0', '--no-xyzmap'
)

& build\Release\scroll_unroll.exe "$source\baseline" "$out\parent" `
  --id parent @common
& build\Release\scroll_unroll.exe "$source\candidate" "$out\candidate" `
  --id candidate @common
& build\Release\scroll_unroll.exe "$source\candidate" "$out\candidate_repeat" `
  --id candidate_repeat @common
```

No parameter, stage, arm order, output root, or retry setting may change after
the first downstream process starts. A nonzero exit or missing output is a
result, not permission to rerun with a different setting.

Run the frozen comparator once with its checked-in defaults:

```powershell
& '..\_sf_sparse\python\.venv\Scripts\python.exe' `
  python\scripts\compare_unroll_outputs.py `
  --parent-dir "$out\parent" --parent-id parent `
  --candidate-dir "$out\candidate" --candidate-id candidate `
  --candidate-repeat-dir "$out\candidate_repeat" `
  --candidate-repeat-id candidate_repeat `
  --out-json "$out\paired_downstream_result.json" `
  --figure-prefix "$out\complete_unroll_strips"
```

## Frozen outputs, gates, and interpretation

Every arm must emit exactly five stage rows and, at every stage, a nonempty
uint8 RAW texture TIF, diagnostic-class TIF, complete readable-strip PNG, and
full preview PNG whose dimensions reconcile with the pipeline JSON. Parent,
candidate, and repeat must have identical cube, vertex, and face counts.

The fresh candidate repeat must have identical non-timing metrics and
byte-identical required TIF/PNG artifacts. Failure of comparability or repeat
determinism is `DOWNSTREAM_INCONCLUSIVE`.

Safety gates, fixed before the output exists:

1. stage-1 and final-stage fill may each fall by at most 0.20 percentage point;
2. final multi-cover fraction may rise by at most 0.20 percentage point;
3. final seam excess `max(0, seam_ratio - 1)` may rise by at most 0.05; and
4. candidate/parent canvas area ratio must remain in `[0.95, 1.05]` at stage 1
   and stage 5.

A safety failure is `DOWNSTREAM_REGRESSION`. Otherwise the downstream efficacy
gate passes only when stage-1 seam excess falls by at least 0.05 and raw
stage-1 seam-column discontinuity does not increase. Passing all gates is
`DOWNSTREAM_EFFICACY_PASS`; passing safety without the efficacy effect is
`DOWNSTREAM_COMPATIBILITY_ONLY`.

The visual contains the complete stage-1 and stage-5 readable-strip products
for both arms. It may not select, crop, hide, or reorder any seam or texture
window. The output can support a fiber-continuity/downstream-texture claim, but
never an ink, legibility, full-scroll, or new-surface-model claim. Positive,
null, inconclusive, or negative outcomes will be published without changing
these rules.
