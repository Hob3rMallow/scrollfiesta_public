# Sparse graph winding constraints: post-hoc development ablation

Status: **development evidence, not a new held-out claim**  
Frozen candidate: `a47a21c3a9d445c90d28e218b8f0fcbecd49c079`  
Public freeze: <https://github.com/TAUIL-Abd-Elilah/scrollfiesta_public/commit/a47a21c3a9d445c90d28e218b8f0fcbecd49c079>

The PHerc0211 preregistered calibration validation failed its fixed decision
rule at commit `e0cf51c`: the candidate improved turn-off from 10.012% to
6.876%, but remained above the 5% ceiling. The result was published as
inconclusive before any follow-up arm was run. Everything below is therefore
post-hoc development, and PHerc0211 is no longer held out for this change.

## Failure attribution

On the first-seed/radius diagnostic arm, all **54/54** remaining whole-turn
errors were in 38 `(cube A group, cube B group)` buckets supported by only one
or two nearest-skin correspondences. There were **zero** errors in buckets
meeting the old `min_edge_pairs=3` floor. These sparse buckets had already
passed the 3-D distance, radial-distance, and near-integer phase gates, but the
floor discarded their graph relations and left the groups on a noisier
component gauge.

Two controlled arms separated this from the later coordinate polish:

| PHerc0211 arm | Turn-off | `|du|<2` | Result |
|---|---:|---:|---|
| three-pair graph, no polish | 83/829 (10.012%) | 81.79% | FAIL |
| one-pair graph, no polish | 5/829 (0.603%) | 91.07% | PASS |
| one-pair graph, unconstrained polish | 35/829 (4.222%) | 87.58% | PASS |
| one-pair graph, graph-supported groups protected | 5/829 (0.603%) | 91.07% | PASS |

The unconstrained `CubeReg_solve` polish pooled all neighbor groups by only the
current group id. It could therefore overwrite explicit `(gid A, gid B)` graph
relations. The candidate now:

1. admits one- and two-correspondence buckets after the existing geometric and
   phase gates; and
2. preserves both integer turn and matching residual-u for every group with an
   admitted graph edge, while retaining polish for graph-isolated groups.

This is not a threshold-only rescue: disabling polish outright failed PHerc1203
at 64/1268 (5.047%). The support-aware guard retains the isolated-group recovery
and finishes at 23/1268 (1.814%).

## Five-block development matrix

Each before/after comparison reuses the exact same raw skins and audit pair
inventory. All 5 pair-set counts and all 15,437 individual pair counts are
unchanged.

| Dataset | Before turn-off | Candidate turn-off | Reduction | `|du|<2` before → after | Collisions before → after |
|---|---:|---:|---:|---:|---:|
| PHerc0139 original, 100 cubes | 167/6441 (2.593%) | 126/6441 (1.956%) | 24.6% | 83.05% → 83.43% | 416 → 257 |
| PHerc0139 adjacent, 25 cubes | 12/1024 (1.172%) | 4/1024 (0.391%) | 66.7% | 89.75% → 90.43% | 88 → 74 |
| PHerc0332, 81 cubes | 223/5875 (3.796%) | 113/5875 (1.923%) | 49.3% | 78.49% → 80.36% | 531 → 299 |
| PHerc1203, 25 cubes | 42/1268 (3.312%) | 23/1268 (1.814%) | 45.2% | 80.21% → 82.18% | 96 → 65 |
| PHerc0211, 25 cubes | 57/829 (6.876%) | 5/829 (0.603%) | 91.2% | 84.80% → 91.07% | 119 → 26 |
| **Weighted total** | **501/15437 (3.245%)** | **271/15437 (1.756%)** | **45.9%** | — | **1250 → 721** |

All five candidate audits pass. The machine-readable result, including audit
hashes, is in
`docs/results/sparse_graph_winding_ablation.json`.

## Reproduction and checks

For each existing placed directory, the candidate was run with its default
settings and no geometry rewrite:

```powershell
& build\Release\scroll_whole.exe <placed_dir> --reregister --no-obj
& build\Release\scroll_whole.exe <placed_dir> --audit
```

- `scroll_whole --selftest`: PASS, including a one-pair graph-edge contract and
  a graph-supported polish-guard contract.
- MSVC Release build: PASS. Only the pre-existing alignment/unreachable-code
  warnings were emitted.
- PHerc0211 repeat determinism: two complete reregister/audit runs produced
  byte-identical `audit.json`, `placed_index.json`, and all 25 registered skin
  sidecars.
- PHerc0211 audit SHA-256:
  `938c1d9484aa08824754fe5e00e8ebef0bd16d17ba0a80a9fe0c40140f30e230`.
- Built `scroll_whole.exe` SHA-256:
  `6624f7fe3ca062d10c03e7a799b1a32e275b8c2c61695a7f1f8360bc8539818e`.

## Scope and next gate

This matrix supports a mechanism and checks regressions, but it is not an
independent efficacy test because every listed block is now development data.
The code was pushed publicly before selecting or reading the next validation
slab. A separate, provenance-audited scroll/slab must be preregistered and must
pass the frozen quantitative and visual gates before an upstream success claim
or default-change PR is opened.

As of the freeze, searches of the open ScrollPrize/villa and upstream
ScrollFiesta PRs found no implementation of sparse group-edge admission plus a
graph-supported polish guard. The Discord export through 2026-08-08 emphasizes
high-confidence same/relative winding relations and warns that attractive
unrolls can still be wrong; it does not describe this fix. The work therefore
complements, rather than duplicates, the public track/patch graph work.
