#!/usr/bin/env python3
"""Validate and render a paired ``scroll_whole --audit`` comparison.

The comparison is deliberately pair-preserving: both arms must report the
same cube geometry, the same initial calibration, and the same adjacent seam
keys with the same correspondence counts.  The generated figure shows every
audited seam, so it is a complete visual rendering of the quantitative result
rather than a selected best/worst crop.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm


CUBE_RE = re.compile(r"^z(?P<z>\d+)_y(?P<y>\d+)_x(?P<x>\d+)$")
INDEX_SHARED_FIELDS = (
    "leaf_stage",
    "chunk",
    "axis_point_zyx",
    "axis_dir_zyx",
    "pitch",
    "pitch_mode",
    "pair_gate",
    "skin_dist",
    "calibration",
)
CUBE_SHARED_FIELDS = (
    "id",
    "origin",
    "status",
    "nv",
    "nf",
    "skin",
    "flood_comp",
    "n_groups",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return value


def cube_origin(cube_id: str) -> tuple[int, int, int]:
    match = CUBE_RE.fullmatch(cube_id)
    if match is None:
        raise ValueError(f"invalid cube id: {cube_id!r}")
    return tuple(int(match.group(axis)) for axis in ("z", "y", "x"))


def _finite_number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label}: expected a number")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{label}: expected a finite number")
    return number


def validate_audit(audit: dict[str, Any], label: str) -> dict[tuple[str, str], dict[str, Any]]:
    pairs = audit.get("pairs")
    if not isinstance(pairs, list) or not pairs:
        raise ValueError(f"{label}: pairs must be a nonempty list")

    indexed: dict[tuple[str, str], dict[str, Any]] = {}
    total_n = total_turn = total_collide = 0
    for position, pair in enumerate(pairs):
        if not isinstance(pair, dict):
            raise ValueError(f"{label}.pairs[{position}]: expected an object")
        a, b = pair.get("a"), pair.get("b")
        if not isinstance(a, str) or not isinstance(b, str):
            raise ValueError(f"{label}.pairs[{position}]: invalid cube ids")
        cube_origin(a)
        cube_origin(b)
        key = (a, b)
        if key in indexed:
            raise ValueError(f"{label}: duplicate pair key {key}")
        n = pair.get("n")
        turn = pair.get("turn_off")
        collide = pair.get("collide")
        if not all(isinstance(v, int) and not isinstance(v, bool) for v in (n, turn, collide)):
            raise ValueError(f"{label} {key}: n/turn_off/collide must be integers")
        # ``collide`` uses AUDIT_COLLIDE_GATE (5 voxels), while ``n`` and
        # ``turn_off`` use the narrower configured pair gate (3.5 voxels by
        # default).  Collision matches can therefore legitimately outnumber
        # the pair-gated correspondences reported in ``n``.
        if n <= 0 or not (0 <= turn <= n) or collide < 0:
            raise ValueError(f"{label} {key}: invalid counts")
        for field in ("dphi_med", "du_med", "du_mad", "dv_absmed"):
            _finite_number(pair.get(field), f"{label} {key}.{field}")
        total_n += n
        total_turn += turn
        total_collide += collide
        indexed[key] = pair

    expected = {
        "n_pair_sets": len(indexed),
        "n_pairs": total_n,
        "turn_off_pairs": total_turn,
        "collisions": total_collide,
    }
    for field, value in expected.items():
        if audit.get(field) != value:
            raise ValueError(f"{label}.{field}: reported {audit.get(field)!r}, recomputed {value}")
    join = audit.get("join_completeness")
    if not isinstance(join, dict) or join.get("n") != total_n:
        raise ValueError(f"{label}: invalid join_completeness count")
    for field in ("lt2", "lt6"):
        fraction = _finite_number(join.get(field), f"{label}.join_completeness.{field}")
        if not 0.0 <= fraction <= 1.0:
            raise ValueError(f"{label}.join_completeness.{field}: outside [0,1]")
    quality = audit.get("quality_gate")
    if not isinstance(quality, dict) or quality.get("status") not in {"PASS", "FAIL"}:
        raise ValueError(f"{label}: invalid quality_gate status")
    return indexed


def validate_indexes(baseline: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    mismatches = [
        field for field in INDEX_SHARED_FIELDS if baseline.get(field) != candidate.get(field)
    ]
    if mismatches:
        raise ValueError(f"placed indexes differ before reregistration: {mismatches}")

    baseline_cubes = baseline.get("cubes")
    candidate_cubes = candidate.get("cubes")
    if not isinstance(baseline_cubes, list) or not isinstance(candidate_cubes, list):
        raise ValueError("placed indexes must contain cube lists")
    if len(baseline_cubes) != len(candidate_cubes) or not baseline_cubes:
        raise ValueError("placed indexes contain different cube counts")

    cube_differences: list[str] = []
    for left, right in zip(baseline_cubes, candidate_cubes, strict=True):
        if not isinstance(left, dict) or not isinstance(right, dict):
            raise ValueError("placed-index cube record is not an object")
        changed = [field for field in CUBE_SHARED_FIELDS if left.get(field) != right.get(field)]
        if changed:
            cube_differences.append(f"{left.get('id')}:{','.join(changed)}")
    if cube_differences:
        raise ValueError(f"cube geometry differs: {cube_differences}")

    for label, index in (("baseline", baseline), ("candidate", candidate)):
        count = len(baseline_cubes)
        if index.get("n_cubes") != count or index.get("n_ok") != count:
            raise ValueError(f"{label}: not every cube is OK")
        if index.get("n_skipped") != 0 or index.get("n_low_conf") != 0:
            raise ValueError(f"{label}: skipped or low-confidence cubes remain")
        if index.get("reregistered") != 1:
            raise ValueError(f"{label}: index is not reregistered")

    return {
        "pass": True,
        "cube_count": len(baseline_cubes),
        "shared_initial_calibration": True,
        "shared_cube_geometry": True,
        "baseline_ok_skipped_low_conf": [
            baseline["n_ok"],
            baseline["n_skipped"],
            baseline["n_low_conf"],
        ],
        "candidate_ok_skipped_low_conf": [
            candidate["n_ok"],
            candidate["n_skipped"],
            candidate["n_low_conf"],
        ],
    }


def compare(
    baseline_audit: dict[str, Any],
    candidate_audit: dict[str, Any],
    baseline_index: dict[str, Any],
    candidate_index: dict[str, Any],
    *,
    min_pairs: int = 100,
    max_candidate_turnoff: float = 0.05,
    min_reduction_pp: float = 1.0,
    max_lt2_drop_pp: float = 0.5,
) -> dict[str, Any]:
    baseline_pairs = validate_audit(baseline_audit, "baseline")
    candidate_pairs = validate_audit(candidate_audit, "candidate")
    if baseline_pairs.keys() != candidate_pairs.keys():
        raise ValueError("audit arms have different pair keys")
    for key in baseline_pairs:
        if baseline_pairs[key]["n"] != candidate_pairs[key]["n"]:
            raise ValueError(f"audit arms have different counts for {key}")

    index_gate = validate_indexes(baseline_index, candidate_index)
    n_pairs = int(baseline_audit["n_pairs"])
    baseline_fraction = baseline_audit["turn_off_pairs"] / n_pairs
    candidate_fraction = candidate_audit["turn_off_pairs"] / n_pairs
    reduction_pp = 100.0 * (baseline_fraction - candidate_fraction)
    baseline_lt2 = float(baseline_audit["join_completeness"]["lt2"])
    candidate_lt2 = float(candidate_audit["join_completeness"]["lt2"])
    lt2_change_pp = 100.0 * (candidate_lt2 - baseline_lt2)

    improved = worsened = equal = 0
    per_pair: list[dict[str, Any]] = []
    for key in sorted(baseline_pairs):
        left, right = baseline_pairs[key], candidate_pairs[key]
        delta = right["turn_off"] - left["turn_off"]
        improved += delta < 0
        worsened += delta > 0
        equal += delta == 0
        per_pair.append(
            {
                "a": key[0],
                "b": key[1],
                "n": left["n"],
                "baseline_turn_off": left["turn_off"],
                "candidate_turn_off": right["turn_off"],
                "delta_turn_off": delta,
                "baseline_fraction": left["turn_off"] / left["n"],
                "candidate_fraction": right["turn_off"] / right["n"],
            }
        )

    gates = {
        "index_and_pair_comparability": index_gate["pass"],
        "at_least_min_pairs": n_pairs >= min_pairs,
        "candidate_builtin_audit_pass": candidate_audit["quality_gate"]["status"] == "PASS",
        "candidate_turnoff_at_most_ceiling": candidate_fraction <= max_candidate_turnoff,
        "turnoff_reduction_at_least_floor": reduction_pp >= min_reduction_pp,
        "join_lt2_drop_within_floor": lt2_change_pp >= -max_lt2_drop_pp,
        "collisions_nonincreasing": candidate_audit["collisions"] <= baseline_audit["collisions"],
    }
    if all(gates.values()):
        decision = "INDEPENDENT_EFFICACY_PASS"
    elif (
        candidate_fraction - baseline_fraction > 0.005
        or lt2_change_pp < -max_lt2_drop_pp
        or candidate_audit["collisions"] > baseline_audit["collisions"]
    ):
        decision = "REGRESSION"
    elif (
        gates["index_and_pair_comparability"]
        and gates["at_least_min_pairs"]
        and gates["candidate_builtin_audit_pass"]
        and gates["candidate_turnoff_at_most_ceiling"]
        and gates["join_lt2_drop_within_floor"]
        and gates["collisions_nonincreasing"]
    ):
        decision = "HELDOUT_COMPATIBILITY_ONLY"
    else:
        decision = "INCONCLUSIVE"

    return {
        "format": "scrollfiesta-paired-winding-audit-v1",
        "decision": decision,
        "comparability": index_gate
        | {
            "pair_keys_and_counts_equal": True,
            "pair_sets": len(baseline_pairs),
            "audit_pairs": n_pairs,
        },
        "baseline": {
            "audit_status": baseline_audit["quality_gate"]["status"],
            "turn_off_pairs": baseline_audit["turn_off_pairs"],
            "turn_off_fraction": baseline_fraction,
            "join_lt2_fraction": baseline_lt2,
            "join_lt6_fraction": baseline_audit["join_completeness"]["lt6"],
            "collisions": baseline_audit["collisions"],
        },
        "candidate": {
            "audit_status": candidate_audit["quality_gate"]["status"],
            "turn_off_pairs": candidate_audit["turn_off_pairs"],
            "turn_off_fraction": candidate_fraction,
            "join_lt2_fraction": candidate_lt2,
            "join_lt6_fraction": candidate_audit["join_completeness"]["lt6"],
            "collisions": candidate_audit["collisions"],
        },
        "paired_effect": {
            "turn_off_pairs_removed": baseline_audit["turn_off_pairs"]
            - candidate_audit["turn_off_pairs"],
            "turn_off_absolute_reduction_pp": reduction_pp,
            "turn_off_relative_reduction": (
                (baseline_fraction - candidate_fraction) / baseline_fraction
                if baseline_fraction > 0.0
                else None
            ),
            "join_lt2_change_pp": lt2_change_pp,
            "join_lt6_change_pp": 100.0
            * (
                candidate_audit["join_completeness"]["lt6"]
                - baseline_audit["join_completeness"]["lt6"]
            ),
            "collisions_removed": baseline_audit["collisions"]
            - candidate_audit["collisions"],
            "pair_sets_improved_worse_equal": [improved, worsened, equal],
        },
        "frozen_thresholds": {
            "min_pairs": min_pairs,
            "max_candidate_turnoff_fraction": max_candidate_turnoff,
            "min_turnoff_reduction_pp": min_reduction_pp,
            "max_join_lt2_drop_pp": max_lt2_drop_pp,
            "collisions_must_not_increase": True,
        },
        "gates": gates,
        "per_pair": per_pair,
    }


def render(result: dict[str, Any], output_prefix: Path) -> tuple[Path, Path]:
    records = result["per_pair"]
    origins: dict[str, tuple[int, int, int]] = {}
    for record in records:
        origins[record["a"]] = cube_origin(record["a"])
        origins[record["b"]] = cube_origin(record["b"])
    z_values = {origin[0] for origin in origins.values()}
    if len(z_values) != 1:
        raise ValueError("the frozen visual renderer requires a one-z-layer slab")

    max_n = max(record["n"] for record in records)
    segments = [
        [(origins[r["a"]][2], origins[r["a"]][1]), (origins[r["b"]][2], origins[r["b"]][1])]
        for r in records
    ]
    widths = [1.5 + 5.0 * math.sqrt(r["n"] / max_n) for r in records]
    baseline_values = [r["baseline_fraction"] for r in records]
    candidate_values = [r["candidate_fraction"] for r in records]
    delta_values = [r["candidate_fraction"] - r["baseline_fraction"] for r in records]

    matplotlib.rcParams["svg.hashsalt"] = "scrollfiesta-winding-audit-v1"
    figure, axes = plt.subplots(1, 3, figsize=(18, 6.4), constrained_layout=True)
    absolute_norm = Normalize(vmin=0.0, vmax=0.20, clip=True)
    delta_norm = TwoSlopeNorm(vmin=-0.20, vcenter=0.0, vmax=0.20)
    panels = (
        ("Parent baseline", baseline_values, "magma", absolute_norm),
        ("Sparse-graph candidate", candidate_values, "magma", absolute_norm),
        ("Candidate - baseline", delta_values, "RdBu_r", delta_norm),
    )
    xs = [origin[2] for origin in origins.values()]
    ys = [origin[1] for origin in origins.values()]
    pad = 80

    for axis, (title, values, cmap, norm) in zip(axes, panels, strict=True):
        collection = LineCollection(
            segments,
            array=values,
            cmap=cmap,
            norm=norm,
            linewidths=widths,
            capstyle="round",
            zorder=2,
        )
        axis.add_collection(collection)
        axis.scatter(xs, ys, s=48, facecolor="white", edgecolor="#202020", linewidth=1.0, zorder=3)
        for cube_id, (_, y, x) in sorted(origins.items()):
            axis.annotate(
                f"{y},{x}",
                (x, y),
                xytext=(0, 7),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=6.5,
                color="#202020",
            )
        axis.set_xlim(min(xs) - pad, max(xs) + pad)
        axis.set_ylim(min(ys) - pad, max(ys) + pad)
        axis.set_aspect("equal")
        axis.set_title(title, fontsize=13, weight="bold")
        axis.set_xlabel("cube x origin (voxels)")
        axis.set_ylabel("cube y origin (voxels)")
        scalar = ScalarMappable(norm=norm, cmap=cmap)
        scalar.set_array([])
        colorbar = figure.colorbar(scalar, ax=axis, fraction=0.046, pad=0.035)
        colorbar.set_label("turn-off fraction" if title != "Candidate - baseline" else "fraction change")

    baseline = result["baseline"]
    candidate = result["candidate"]
    effect = result["paired_effect"]
    figure.suptitle(
        "Complete held-out seam map (all pair sets; no selected edges)\n"
        f"turn-off {baseline['turn_off_pairs']}/{result['comparability']['audit_pairs']} "
        f"({100*baseline['turn_off_fraction']:.2f}%) -> "
        f"{candidate['turn_off_pairs']}/{result['comparability']['audit_pairs']} "
        f"({100*candidate['turn_off_fraction']:.2f}%); "
        f"change {-effect['turn_off_absolute_reduction_pp']:+.2f} pp; "
        f"decision {result['decision']}",
        fontsize=13,
    )
    figure.text(
        0.5,
        0.005,
        "Line width encodes sqrt(correspondence count); fixed color ranges: absolute 0-20%, change -20 to +20 pp.",
        ha="center",
        fontsize=9,
    )

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    png = output_prefix.with_suffix(".png")
    svg = output_prefix.with_suffix(".svg")
    figure.savefig(png, dpi=180, metadata={"Software": "compare_winding_audits.py v1"})
    figure.savefig(svg, metadata={"Date": None, "Creator": "compare_winding_audits.py v1"})
    plt.close(figure)
    return png, svg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-audit", type=Path, required=True)
    parser.add_argument("--candidate-audit", type=Path, required=True)
    parser.add_argument("--baseline-index", type=Path, required=True)
    parser.add_argument("--candidate-index", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--figure-prefix", type=Path, required=True)
    parser.add_argument("--min-pairs", type=int, default=100)
    parser.add_argument("--max-candidate-turnoff", type=float, default=0.05)
    parser.add_argument("--min-reduction-pp", type=float, default=1.0)
    parser.add_argument("--max-lt2-drop-pp", type=float, default=0.5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inputs = {
        "baseline_audit": args.baseline_audit,
        "candidate_audit": args.candidate_audit,
        "baseline_index": args.baseline_index,
        "candidate_index": args.candidate_index,
    }
    result = compare(
        load_json(args.baseline_audit),
        load_json(args.candidate_audit),
        load_json(args.baseline_index),
        load_json(args.candidate_index),
        min_pairs=args.min_pairs,
        max_candidate_turnoff=args.max_candidate_turnoff,
        min_reduction_pp=args.min_reduction_pp,
        max_lt2_drop_pp=args.max_lt2_drop_pp,
    )
    png, svg = render(result, args.figure_prefix)
    result["evidence"] = {
        "input_sha256": {name: sha256(path) for name, path in inputs.items()},
        "figure_png": {"path": str(png), "sha256": sha256(png)},
        "figure_svg": {"path": str(svg), "sha256": sha256(svg)},
        "visual_scope": "Every audited pair set is rendered; no edge is selected or hidden.",
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"decision": result["decision"], "out": str(args.out_json)}, sort_keys=True))


if __name__ == "__main__":
    main()
