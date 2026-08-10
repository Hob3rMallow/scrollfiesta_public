# PHerc0211 held-out placement validation: improvement, but no pass

Result: **inconclusive under the frozen protocol; no held-out efficacy or
production-default claim**.

The preregistration was published at
[`0b5493c`](https://github.com/TAUIL-Abd-Elilah/scrollfiesta_public/commit/0b5493cfbfc6f287a4ecbb0813db07fb4bc33195)
before the selected PHerc0211 slab was opened. The algorithm under test was
the already committed `19152e49acf855d2bbe215329050ab0cc0085b5e`.

## Outcome

Both arms consumed the same 25 successfully meshed cubes and the same 34
adjacent pair-sets / 829 audit pairs.

| Frozen metric | Former default | Candidate default | Change |
|---|---:|---:|---:|
| Turn-off pairs | 83/829 (10.012%) | 57/829 (6.876%) | -26, **-3.136 pp** |
| Join completeness `|du|<2` | 82.27% | 84.80% | **+2.53 pp** |
| Join completeness `|du|<6` | 89.51% | 92.16% | **+2.65 pp** |
| Collisions below 5 voxels | 171 | 119 | **-52** |
| Built-in audit | FAIL | FAIL | candidate remains above 5% |

The turn-off reduction is 31.3% relative. It was not produced by sacrificing
another pair-set: 13 of 34 pair-sets improved, none worsened, and 21 were
unchanged. Cube IDs and every cube's vertex count, face count, skin count, and
status were identical across arms. Both ended with 25 OK, zero skipped, and
zero low-confidence inputs after reregistration.

The candidate used pitch 8.3661 vox/turn versus the former first-seed estimate
of 8.2617. All 25 candidate estimates were supported, but their sign vote was
only 13 positive to 12 negative. That is 52%, below the preregistered 60%
stability floor. Absolute candidate pitches ranged from 4.4350 to 14.8074
vox/turn, with quartiles 7.3431 / 8.3661 / 9.2877. The consensus gauge used
radius shifts for 74 graph components and retained the raw chart for 116.

## Gate decision

The input, mesh, topology, pair-identity, minimum-pair, and completeness-safety
gates passed. The 3.136-point improvement exceeded the fixed 1.00-point effect
floor, and `|du|<2` improved rather than declining. Two mandatory gates failed:

1. candidate turn-off remained 6.876%, above the fixed 5% audit ceiling; and
2. selected winding-sign support was 52%, below the fixed 60% floor.

The result is therefore neither a pass nor a regression. It is evidence that
multi-seed/consensus placement removes a broad subset of former errors on a
fourth scroll, but it does not establish that the current implementation is a
safe production default. No alternate pitch, sign, gauge, slab, or threshold
was tried after observing this outcome.

## Reproducibility anchors

- selected box: z/y/x `[5888,6016) x [3584,4224) x [3200,3840)`;
- input manifest SHA-256:
  `ce034435a2ec98e67cec656383c67827fb77777fbf9e1cae5dadf157e2fd307e`;
- input validation SHA-256:
  `2af7224c09c13b50b6dee3dec228b03e4835b265854e4fa4a1cf1f9046dcdd13`;
- pipeline summary SHA-256:
  `9a5c962f7f6f4da5cfbd90e2367633e94a122cf8b847454062e53ca426f40477`;
- former audit JSON SHA-256:
  `1cd070508123448d6a8a0489e99bb09b000214d5e687f5c5597ebbf24d6d610c`;
- candidate audit JSON SHA-256:
  `51ea67d956dc47cd5f807935509d2a18de1d9dab570968d28d37cd45256352bf`;
- structured result: `docs/results/pherc0211_placement_validation.json`.

The complete development results on PHerc0139, PHerc0332, and PHerc1203 are
listed in the preregistration and are not relabeled as held-out evidence.
