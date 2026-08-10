# PHerc0826 scale comparator schema erratum

Status: **published before the corrected comparator produced any stratum
aggregate, decision JSON, or figure**.

The first execution of the frozen PHerc0826 scale comparator stopped before
classification or rendering with:

```text
ValueError: baseline
('z12672_y04608_x04352', 'z12672_y04608_x04480'): invalid counts
```

No `paired_scale_result.json` or `complete_scale_seam_map` existed after that
execution. The data, registrations, repeat, thresholds, boxes, strata, and
decision rules are unchanged.

## Cause

The inherited audit validator required `0 <= collide <= n`. That relationship
is not part of the audit schema. In `scroll_whole.c`, `collide` counts nearest
matches inside the fixed 5-voxel `AUDIT_COLLIDE_GATE`, while `n` and
`turn_off` count the subset inside the narrower configured pair gate (3.5
voxels here). The source increments `collide` before rejecting matches outside
the pair gate. Therefore `collide > n` is valid by construction.

The larger baseline contains eight such pair sets and the candidate contains
five. In both arms the per-pair sums exactly equal the reported aggregates:

- parent: 10,285 pair-gated correspondences and 1,188 collision-gated
  matches;
- candidate: 10,285 pair-gated correspondences and 900 collision-gated
  matches.

All `n`, `turn_off`, and `collide` fields are nonnegative integers,
`n > 0`, and `turn_off <= n`. This is a validator assumption error, not a
malformed native audit.

## Minimal correction

Only this validation predicate changes:

```diff
- n <= 0 or not (0 <= turn <= n) or not (0 <= collide <= n)
+ n <= 0 or not (0 <= turn <= n) or collide < 0
```

No metric is capped, discarded, recomputed, or reweighted. Collision totals
and the frozen nonincrease gates remain exact. A regression test now supplies
a valid wider-gate collision count above `n` and requires the original totals
to survive comparison.

- preregistered inherited validator: 18,815 bytes, SHA-256
  `e8c9d492e4eb914d005614f9c50bc7f45089ba3cef45ccf01ba70d5f9b578def`;
- corrected `compare_winding_audits.py`: 19,554 bytes, SHA-256
  `2e6b0e41ff3cfab59186b6dddeaf1b8399435396dfb51771ac9f4e7809ab02c5`;
- corrected tests: 5,530 bytes, SHA-256
  `7b8557984bff5232025aeafb21622226f7fec56b1876593cb4b7682c614085af`;
- eight focused base/scale comparator tests: PASS.

The scale script itself remains byte-identical to the preregistration at
SHA-256
`5294ae7253edaef02084af1f9837d8a177532e270f7d21bfc0295bc45db2ee97`.

## Frozen result inputs

These files already existed before the failed parser call and are fixed for
the corrected call:

| input | SHA-256 |
|---|---|
| parent audit | `85efdd6b916dd65589add0956a8a793650d88d6add77b1e4e371ea9faa066d59` |
| candidate audit | `7a06b908f45b20537a8eaafc828b29f0f1e0b3c0007d9ed1b2a4d81b05611f2f` |
| candidate-repeat audit | `7a06b908f45b20537a8eaafc828b29f0f1e0b3c0007d9ed1b2a4d81b05611f2f` |
| parent index | `4f0793f9f1a05a94bea13bb1dc921c4b80a8b3c9b3a8b44afe456f0e58255972` |
| candidate index | `eb6bd4ffb47905e0e1899b946378a52c44af5bbf94cd81fd4aedd7568d51f7d` |
| candidate-repeat index | `eb6bd4ffb47905e0e1899b946378a52c44af5bbf94cd81fd4aedd7568d51f7d` |

The corrected comparator will be run once on exactly these inputs after this
erratum and patch are public. The result will explicitly disclose this
protocol-preserving implementation correction; it will not be described as
an execution of the original validator hash.
