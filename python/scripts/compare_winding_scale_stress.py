#!/usr/bin/env python3
"""Compare a multi-layer winding-registration scale stress test.

The primary stratum contains only seams whose two endpoint cubes are outside
the previously evaluated box.  Seams internal to that known box and seams
crossing its boundary remain visible and are reported, but cannot make the
scale-stress decision pass.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_winding_audits as base  # noqa: E402


CUBE = 128


def validate_bbox(values: tuple[int, ...], label: str) -> tuple[int, ...]:
    if len(values) != 6:
        raise ValueError(f"{label}: expected six z0 z1 y0 y1 x0 x1 values")
    z0, z1, y0, y1, x0, x1 = values
    if not (z0 < z1 and y0 < y1 and x0 < x1):
        raise ValueError(f"{label}: each lower bound must be below its upper bound")
    if any(value % CUBE for value in values):
        raise ValueError(f"{label}: every bound must be a multiple of {CUBE}")
    return values


def contains(bbox: tuple[int, ...], origin: tuple[int, int, int]) -> bool:
    z0, z1, y0, y1, x0, x1 = bbox
    z, y, x = origin
    return z0 <= z < z1 and y0 <= y < y1 and x0 <= x < x1


def lattice_ids(bbox: tuple[int, ...]) -> set[str]:
    z0, z1, y0, y1, x0, x1 = bbox
    return {
        f"z{z:05d}_y{y:05d}_x{x:05d}"
        for z in range(z0, z1, CUBE)
        for y in range(y0, y1, CUBE)
        for x in range(x0, x1, CUBE)
    }


def _manifest(root: Path) -> dict[str, Any]:
    files = sorted(root.glob("*_skin.f32"))
    rows = []
    for path in files:
        rows.append((path.name, path.stat().st_size, base.sha256(path)))
    canonical = "".join(f"{name}\t{size}\t{digest}\n" for name, size, digest in rows)
    return {
        "count": len(rows),
        "canonical_sha256": hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
        "files": {name: {"bytes": size, "sha256": digest} for name, size, digest in rows},
    }


def verify_repeat(
    candidate_audit_path: Path,
    candidate_index_path: Path,
    repeat_audit_path: Path,
    repeat_index_path: Path,
    candidate_dir: Path,
    repeat_dir: Path,
) -> dict[str, Any]:
    candidate_index = base.load_json(candidate_index_path)
    expected_skins = candidate_index.get("n_cubes")
    candidate_manifest = _manifest(candidate_dir)
    repeat_manifest = _manifest(repeat_dir)
    audit_equal = candidate_audit_path.read_bytes() == repeat_audit_path.read_bytes()
    index_equal = candidate_index_path.read_bytes() == repeat_index_path.read_bytes()
    skin_inventory_equal = candidate_manifest["files"] == repeat_manifest["files"]
    skin_count_equal_index = (
        isinstance(expected_skins, int)
        and candidate_manifest["count"] == expected_skins
        and repeat_manifest["count"] == expected_skins
    )
    return {
        "pass": audit_equal and index_equal and skin_inventory_equal and skin_count_equal_index,
        "audit_byte_identical": audit_equal,
        "index_byte_identical": index_equal,
        "skin_inventory_byte_identical": skin_inventory_equal,
        "skin_count_equals_index_cube_count": skin_count_equal_index,
        "candidate_skin_count": candidate_manifest["count"],
        "repeat_skin_count": repeat_manifest["count"],
        "candidate_skin_manifest_sha256": candidate_manifest["canonical_sha256"],
        "repeat_skin_manifest_sha256": repeat_manifest["canonical_sha256"],
    }


def _summarize(records: list[dict[str, Any]]) -> dict[str, Any]:
    n = sum(record["n"] for record in records)
    baseline_turn = sum(record["baseline_turn_off"] for record in records)
    candidate_turn = sum(record["candidate_turn_off"] for record in records)
    baseline_collisions = sum(record["baseline_collisions"] for record in records)
    candidate_collisions = sum(record["candidate_collisions"] for record in records)
    baseline_fraction = baseline_turn / n if n else None
    candidate_fraction = candidate_turn / n if n else None
    reduction_pp = (
        100.0 * (baseline_fraction - candidate_fraction)
        if baseline_fraction is not None and candidate_fraction is not None
        else None
    )
    improved = sum(record["delta_turn_off"] < 0 for record in records)
    worsened = sum(record["delta_turn_off"] > 0 for record in records)
    equal = len(records) - improved - worsened
    return {
        "pair_sets": len(records),
        "audit_pairs": n,
        "baseline_turn_off_pairs": baseline_turn,
        "candidate_turn_off_pairs": candidate_turn,
        "baseline_turn_off_fraction": baseline_fraction,
        "candidate_turn_off_fraction": candidate_fraction,
        "turn_off_pairs_removed": baseline_turn - candidate_turn,
        "turn_off_absolute_reduction_pp": reduction_pp,
        "turn_off_relative_reduction": (
            (baseline_fraction - candidate_fraction) / baseline_fraction
            if baseline_fraction not in (None, 0.0) and candidate_fraction is not None
            else None
        ),
        "baseline_collisions": baseline_collisions,
        "candidate_collisions": candidate_collisions,
        "collisions_removed": baseline_collisions - candidate_collisions,
        "pair_sets_improved_worse_equal": [improved, worsened, equal],
    }


def compare_scale(
    baseline_audit: dict[str, Any],
    candidate_audit: dict[str, Any],
    baseline_index: dict[str, Any],
    candidate_index: dict[str, Any],
    repeat: dict[str, Any],
    *,
    core_bbox: tuple[int, ...],
    known_bbox: tuple[int, ...],
    min_core_cubes: int = 100,
    min_total_pairs: int = 1000,
    min_new_pair_sets: int = 20,
    min_new_pairs: int = 500,
    max_candidate_turnoff: float = 0.05,
    min_new_reduction_pp: float = 1.0,
    max_turnoff_worsening_pp: float = 0.5,
    max_lt2_drop_pp: float = 0.5,
    require_complete_known_box: bool = True,
) -> dict[str, Any]:
    core_bbox = validate_bbox(core_bbox, "core_bbox")
    known_bbox = validate_bbox(known_bbox, "known_bbox")
    if not (
        core_bbox[0] <= known_bbox[0] < known_bbox[1] <= core_bbox[1]
        and core_bbox[2] <= known_bbox[2] < known_bbox[3] <= core_bbox[3]
        and core_bbox[4] <= known_bbox[4] < known_bbox[5] <= core_bbox[5]
    ):
        raise ValueError("known_bbox must be fully contained in core_bbox")

    baseline_pairs = base.validate_audit(baseline_audit, "baseline")
    candidate_pairs = base.validate_audit(candidate_audit, "candidate")
    if baseline_pairs.keys() != candidate_pairs.keys():
        raise ValueError("audit arms have different pair keys")
    for key in baseline_pairs:
        if baseline_pairs[key]["n"] != candidate_pairs[key]["n"]:
            raise ValueError(f"audit arms have different counts for {key}")
    index_gate = base.validate_indexes(baseline_index, candidate_index)

    cubes = baseline_index["cubes"]
    cube_ids: set[str] = set()
    for cube in cubes:
        cube_id = cube.get("id")
        if not isinstance(cube_id, str):
            raise ValueError("placed-index cube has no valid id")
        origin = base.cube_origin(cube_id)
        if tuple(cube.get("origin", ())) != origin:
            raise ValueError(f"placed-index origin disagrees with id: {cube_id}")
        if not contains(core_bbox, origin):
            raise ValueError(f"cube lies outside frozen core bbox: {cube_id}")
        if cube_id in cube_ids:
            raise ValueError(f"duplicate placed-index cube: {cube_id}")
        cube_ids.add(cube_id)

    known_expected = lattice_ids(known_bbox)
    known_present = cube_ids & known_expected
    known_complete = known_present == known_expected
    if require_complete_known_box and not known_complete:
        missing = sorted(known_expected - known_present)
        raise ValueError(f"known box is incomplete; missing {missing}")

    records: list[dict[str, Any]] = []
    for key in sorted(baseline_pairs):
        if key[0] not in cube_ids or key[1] not in cube_ids:
            raise ValueError(f"audit pair endpoint is absent from placed index: {key}")
        origin_a, origin_b = base.cube_origin(key[0]), base.cube_origin(key[1])
        delta = tuple(abs(left - right) for left, right in zip(origin_a, origin_b, strict=True))
        if sorted(delta) != [0, 0, CUBE]:
            raise ValueError(f"audit pair is not face-adjacent: {key}")
        a_known = key[0] in known_expected
        b_known = key[1] in known_expected
        stratum = "known_center" if a_known and b_known else "bridge" if a_known or b_known else "new_only"
        left, right = baseline_pairs[key], candidate_pairs[key]
        records.append(
            {
                "a": key[0],
                "b": key[1],
                "n": left["n"],
                "stratum": stratum,
                "baseline_turn_off": left["turn_off"],
                "candidate_turn_off": right["turn_off"],
                "delta_turn_off": right["turn_off"] - left["turn_off"],
                "baseline_fraction": left["turn_off"] / left["n"],
                "candidate_fraction": right["turn_off"] / right["n"],
                "baseline_collisions": left["collide"],
                "candidate_collisions": right["collide"],
            }
        )

    strata = {
        name: _summarize([record for record in records if record["stratum"] == name])
        for name in ("new_only", "bridge", "known_center")
    }
    new = strata["new_only"]
    total_pairs = int(baseline_audit["n_pairs"])
    baseline_fraction = baseline_audit["turn_off_pairs"] / total_pairs
    candidate_fraction = candidate_audit["turn_off_pairs"] / total_pairs
    overall_reduction_pp = 100.0 * (baseline_fraction - candidate_fraction)
    baseline_lt2 = float(baseline_audit["join_completeness"]["lt2"])
    candidate_lt2 = float(candidate_audit["join_completeness"]["lt2"])
    lt2_change_pp = 100.0 * (candidate_lt2 - baseline_lt2)

    gates = {
        "index_and_pair_comparability": index_gate["pass"],
        "at_least_min_core_cubes": len(cube_ids) >= min_core_cubes,
        "known_center_complete": known_complete,
        "at_least_min_total_pairs": total_pairs >= min_total_pairs,
        "at_least_min_new_only_pair_sets": new["pair_sets"] >= min_new_pair_sets,
        "at_least_min_new_only_pairs": new["audit_pairs"] >= min_new_pairs,
        "candidate_builtin_audit_pass": candidate_audit["quality_gate"]["status"] == "PASS",
        "candidate_overall_turnoff_at_most_ceiling": candidate_fraction <= max_candidate_turnoff,
        "candidate_new_only_turnoff_at_most_ceiling": (
            new["candidate_turn_off_fraction"] is not None
            and new["candidate_turn_off_fraction"] <= max_candidate_turnoff
        ),
        "new_only_reduction_at_least_floor": (
            new["turn_off_absolute_reduction_pp"] is not None
            and new["turn_off_absolute_reduction_pp"] >= min_new_reduction_pp
        ),
        "overall_turnoff_worsening_within_floor": (
            100.0 * (candidate_fraction - baseline_fraction) <= max_turnoff_worsening_pp
        ),
        "new_only_turnoff_worsening_within_floor": (
            new["candidate_turn_off_fraction"] is not None
            and new["baseline_turn_off_fraction"] is not None
            and 100.0
            * (new["candidate_turn_off_fraction"] - new["baseline_turn_off_fraction"])
            <= max_turnoff_worsening_pp
        ),
        "join_lt2_drop_within_floor": lt2_change_pp >= -max_lt2_drop_pp,
        "overall_collisions_nonincreasing": (
            candidate_audit["collisions"] <= baseline_audit["collisions"]
        ),
        "new_only_collisions_nonincreasing": (
            new["candidate_collisions"] <= new["baseline_collisions"]
        ),
        "candidate_repeat_byte_identical": bool(repeat.get("pass")),
    }

    safety_names = (
        "overall_turnoff_worsening_within_floor",
        "new_only_turnoff_worsening_within_floor",
        "join_lt2_drop_within_floor",
        "overall_collisions_nonincreasing",
        "new_only_collisions_nonincreasing",
    )
    structural_names = (
        "index_and_pair_comparability",
        "at_least_min_core_cubes",
        "known_center_complete",
        "at_least_min_total_pairs",
        "at_least_min_new_only_pair_sets",
        "at_least_min_new_only_pairs",
        "candidate_repeat_byte_identical",
    )
    efficacy_names = (
        "candidate_builtin_audit_pass",
        "candidate_overall_turnoff_at_most_ceiling",
        "candidate_new_only_turnoff_at_most_ceiling",
        "new_only_reduction_at_least_floor",
    )
    if all(gates.values()):
        decision = "SCALE_STRESS_PASS"
    elif all(gates[name] for name in structural_names) and not all(
        gates[name] for name in safety_names
    ):
        decision = "SCALE_STRESS_REGRESSION"
    elif all(gates[name] for name in structural_names + safety_names) and all(
        gates[name] for name in efficacy_names[:-1]
    ):
        decision = "SCALE_STRESS_COMPATIBILITY_ONLY"
    else:
        decision = "SCALE_STRESS_INCONCLUSIVE"

    return {
        "format": "scrollfiesta-winding-scale-stress-v1",
        "decision": decision,
        "scope": (
            "Post-result scale/visual stress test. Only new_only seams can satisfy the "
            "efficacy floor; this is not a second independent holdout."
        ),
        "boxes_l0_zyx": {"core": list(core_bbox), "known_center": list(known_bbox)},
        "comparability": index_gate
        | {
            "pair_keys_and_counts_equal": True,
            "core_cube_count": len(cube_ids),
            "known_cube_count": len(known_present),
            "new_cube_count": len(cube_ids - known_expected),
            "pair_sets": len(records),
            "audit_pairs": total_pairs,
        },
        "overall": {
            "baseline_turn_off_pairs": baseline_audit["turn_off_pairs"],
            "candidate_turn_off_pairs": candidate_audit["turn_off_pairs"],
            "baseline_turn_off_fraction": baseline_fraction,
            "candidate_turn_off_fraction": candidate_fraction,
            "turn_off_absolute_reduction_pp": overall_reduction_pp,
            "baseline_join_lt2_fraction": baseline_lt2,
            "candidate_join_lt2_fraction": candidate_lt2,
            "join_lt2_change_pp": lt2_change_pp,
            "baseline_collisions": baseline_audit["collisions"],
            "candidate_collisions": candidate_audit["collisions"],
        },
        "strata": strata,
        "repeat": repeat,
        "frozen_thresholds": {
            "min_core_cubes": min_core_cubes,
            "min_total_pairs": min_total_pairs,
            "min_new_only_pair_sets": min_new_pair_sets,
            "min_new_only_pairs": min_new_pairs,
            "max_candidate_turnoff_fraction": max_candidate_turnoff,
            "min_new_only_turnoff_reduction_pp": min_new_reduction_pp,
            "max_turnoff_worsening_pp": max_turnoff_worsening_pp,
            "max_join_lt2_drop_pp": max_lt2_drop_pp,
            "overall_and_new_only_collisions_must_not_increase": True,
            "candidate_repeat_must_be_byte_identical": True,
        },
        "gates": gates,
        "per_pair": records,
    }


def render(result: dict[str, Any], output_prefix: Path) -> tuple[Path, Path]:
    records = result["per_pair"]
    if not records:
        raise ValueError("cannot render an empty pair inventory")
    origins: dict[str, tuple[int, int, int]] = {}
    for record in records:
        origins[record["a"]] = base.cube_origin(record["a"])
        origins[record["b"]] = base.cube_origin(record["b"])

    views = (
        ("x-y projection", 2, 1, "x origin", "y origin"),
        ("x-z projection", 2, 0, "x origin", "z origin"),
        ("y-z projection", 1, 0, "y origin", "z origin"),
    )
    metrics = (
        ("Parent baseline", [r["baseline_fraction"] for r in records], "magma", Normalize(0, 0.20, clip=True)),
        ("Sparse-graph candidate", [r["candidate_fraction"] for r in records], "magma", Normalize(0, 0.20, clip=True)),
        (
            "Candidate - parent",
            [r["candidate_fraction"] - r["baseline_fraction"] for r in records],
            "RdBu_r",
            TwoSlopeNorm(vmin=-0.20, vcenter=0.0, vmax=0.20),
        ),
    )
    max_n = max(record["n"] for record in records)
    widths = [0.8 + 3.2 * math.sqrt(record["n"] / max_n) for record in records]
    known_bbox = tuple(result["boxes_l0_zyx"]["known_center"])

    matplotlib.rcParams["svg.hashsalt"] = "scrollfiesta-winding-scale-stress-v1"
    figure, axes = plt.subplots(3, 3, figsize=(18, 16), constrained_layout=True)
    for row, (view_name, horizontal, vertical, xlabel, ylabel) in enumerate(views):
        segments = [
            [
                (origins[r["a"]][horizontal], origins[r["a"]][vertical]),
                (origins[r["b"]][horizontal], origins[r["b"]][vertical]),
            ]
            for r in records
        ]
        h_values = [origin[horizontal] for origin in origins.values()]
        v_values = [origin[vertical] for origin in origins.values()]
        for col, (metric_name, values, cmap, norm) in enumerate(metrics):
            axis = axes[row, col]
            collection = LineCollection(
                segments,
                array=values,
                cmap=cmap,
                norm=norm,
                linewidths=widths,
                capstyle="round",
                alpha=0.88,
                zorder=2,
            )
            axis.add_collection(collection)
            axis.scatter(h_values, v_values, s=7, color="#444444", alpha=0.45, zorder=3)
            horizontal_bounds = (known_bbox[2 * horizontal], known_bbox[2 * horizontal + 1])
            vertical_bounds = (known_bbox[2 * vertical], known_bbox[2 * vertical + 1])
            axis.add_patch(
                Rectangle(
                    (horizontal_bounds[0], vertical_bounds[0]),
                    horizontal_bounds[1] - horizontal_bounds[0],
                    vertical_bounds[1] - vertical_bounds[0],
                    fill=False,
                    edgecolor="#00a6d6",
                    linewidth=1.4,
                    linestyle="--",
                    zorder=4,
                )
            )
            axis.set_xlim(min(h_values) - 80, max(h_values) + 80)
            axis.set_ylim(min(v_values) - 80, max(v_values) + 80)
            axis.set_aspect("equal")
            axis.set_title(f"{metric_name} — {view_name}", fontsize=11, weight="bold")
            axis.set_xlabel(f"cube {xlabel} (voxels)")
            axis.set_ylabel(f"cube {ylabel} (voxels)")

    for col, (_, _, cmap, norm) in enumerate(metrics):
        scalar = ScalarMappable(norm=norm, cmap=cmap)
        scalar.set_array([])
        colorbar = figure.colorbar(scalar, ax=axes[:, col], fraction=0.02, pad=0.015)
        colorbar.set_label("turn-off fraction" if col < 2 else "fraction change")

    new = result["strata"]["new_only"]
    if new["audit_pairs"]:
        new_summary = (
            f"new-only primary stratum {new['baseline_turn_off_pairs']}/{new['audit_pairs']} "
            f"({100 * new['baseline_turn_off_fraction']:.2f}%) -> "
            f"{new['candidate_turn_off_pairs']}/{new['audit_pairs']} "
            f"({100 * new['candidate_turn_off_fraction']:.2f}%); "
            f"reduction {new['turn_off_absolute_reduction_pp']:.2f} pp"
        )
    else:
        new_summary = "new-only primary stratum contains no audited pairs"
    figure.suptitle(
        "Complete multi-layer seam map: all pair sets in three orthogonal projections\n"
        f"{new_summary}; decision {result['decision']}",
        fontsize=14,
    )
    figure.text(
        0.5,
        0.002,
        "Every audited seam is supplied to every projection; dashed cyan box is the previously tested 25-cube center. "
        "Line width encodes sqrt(correspondence count); fixed color ranges are 0-20% and -20 to +20 pp.",
        ha="center",
        fontsize=9,
    )
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    png = output_prefix.with_suffix(".png")
    svg = output_prefix.with_suffix(".svg")
    figure.savefig(png, dpi=180, metadata={"Software": "compare_winding_scale_stress.py v1"})
    figure.savefig(svg, metadata={"Date": None, "Creator": "compare_winding_scale_stress.py v1"})
    plt.close(figure)
    return png, svg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-audit", type=Path, required=True)
    parser.add_argument("--candidate-audit", type=Path, required=True)
    parser.add_argument("--candidate-repeat-audit", type=Path, required=True)
    parser.add_argument("--baseline-index", type=Path, required=True)
    parser.add_argument("--candidate-index", type=Path, required=True)
    parser.add_argument("--candidate-repeat-index", type=Path, required=True)
    parser.add_argument("--candidate-dir", type=Path, required=True)
    parser.add_argument("--candidate-repeat-dir", type=Path, required=True)
    parser.add_argument("--core-bbox", nargs=6, type=int, required=True)
    parser.add_argument("--known-bbox", nargs=6, type=int, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--figure-prefix", type=Path, required=True)
    parser.add_argument("--min-core-cubes", type=int, default=100)
    parser.add_argument("--min-total-pairs", type=int, default=1000)
    parser.add_argument("--min-new-pair-sets", type=int, default=20)
    parser.add_argument("--min-new-pairs", type=int, default=500)
    parser.add_argument("--max-candidate-turnoff", type=float, default=0.05)
    parser.add_argument("--min-new-reduction-pp", type=float, default=1.0)
    parser.add_argument("--max-turnoff-worsening-pp", type=float, default=0.5)
    parser.add_argument("--max-lt2-drop-pp", type=float, default=0.5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repeat = verify_repeat(
        args.candidate_audit,
        args.candidate_index,
        args.candidate_repeat_audit,
        args.candidate_repeat_index,
        args.candidate_dir,
        args.candidate_repeat_dir,
    )
    inputs = {
        "baseline_audit": args.baseline_audit,
        "candidate_audit": args.candidate_audit,
        "candidate_repeat_audit": args.candidate_repeat_audit,
        "baseline_index": args.baseline_index,
        "candidate_index": args.candidate_index,
        "candidate_repeat_index": args.candidate_repeat_index,
    }
    result = compare_scale(
        base.load_json(args.baseline_audit),
        base.load_json(args.candidate_audit),
        base.load_json(args.baseline_index),
        base.load_json(args.candidate_index),
        repeat,
        core_bbox=tuple(args.core_bbox),
        known_bbox=tuple(args.known_bbox),
        min_core_cubes=args.min_core_cubes,
        min_total_pairs=args.min_total_pairs,
        min_new_pair_sets=args.min_new_pair_sets,
        min_new_pairs=args.min_new_pairs,
        max_candidate_turnoff=args.max_candidate_turnoff,
        min_new_reduction_pp=args.min_new_reduction_pp,
        max_turnoff_worsening_pp=args.max_turnoff_worsening_pp,
        max_lt2_drop_pp=args.max_lt2_drop_pp,
    )
    png, svg = render(result, args.figure_prefix)
    result["evidence"] = {
        "input_sha256": {name: base.sha256(path) for name, path in inputs.items()},
        "figure_png": {"path": str(png), "sha256": base.sha256(png)},
        "figure_svg": {"path": str(svg), "sha256": base.sha256(svg)},
        "visual_scope": (
            "Every audited pair set is rendered in three orthogonal projections; "
            "no seam is selected or hidden."
        ),
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"decision": result["decision"], "out": str(args.out_json)}, sort_keys=True))


if __name__ == "__main__":
    main()
