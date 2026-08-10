#!/usr/bin/env python3
"""Validate and render a paired ``scroll_unroll`` downstream comparison.

The comparison consumes complete outputs from one parent placement, one
candidate placement, and a fresh repeat of the candidate.  It validates every
stage inventory and renders the complete readable-strip products from stage 1
and stage 5.  It never selects a seam, crop, or favorable window.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import tifffile


EXPECTED_STAGES = (
    "step1_ribbon",
    "step2_join",
    "step3_overlap",
    "step4_snap",
    "step5_relax",
)
DISPLAY_STAGES = ("step1_ribbon", "step5_relax")
REQUIRED_ARTIFACTS = (
    "rawtex.tif",
    "diagclass.tif",
    "rawtex_strip.png",
    "rawtex_preview.png",
)
FINITE_STAGE_FIELDS = (
    "fill",
    "multi_frac",
    "dark_frac",
    "seam_dcol_mean",
    "base_dcol_mean",
    "seam_ratio",
    "vseam_gap_fill",
    "vseam_ratio",
    "synth_frac",
    "fill_total",
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


def _stage_map(stats: dict[str, Any], label: str) -> dict[str, dict[str, Any]]:
    if stats.get("tool") != "scroll_unroll":
        raise ValueError(f"{label}: unexpected tool {stats.get('tool')!r}")
    stages = stats.get("stages")
    if not isinstance(stages, list):
        raise ValueError(f"{label}: missing stages list")
    names = tuple(stage.get("name") for stage in stages if isinstance(stage, dict))
    if names != EXPECTED_STAGES:
        raise ValueError(f"{label}: expected stages {EXPECTED_STAGES}, got {names}")

    mapped: dict[str, dict[str, Any]] = {}
    for stage in stages:
        name = stage["name"]
        for key in ("W", "H", "bands", *FINITE_STAGE_FIELDS):
            if key not in stage:
                raise ValueError(f"{label}/{name}: missing {key}")
        if int(stage["W"]) <= 0 or int(stage["H"]) <= 0:
            raise ValueError(f"{label}/{name}: invalid canvas")
        for key in FINITE_STAGE_FIELDS:
            if not math.isfinite(float(stage[key])):
                raise ValueError(f"{label}/{name}: nonfinite {key}")
        mapped[name] = stage
    return mapped


def _artifact_path(directory: Path, run_id: str, stage: str, suffix: str) -> Path:
    return directory / f"{run_id}_{stage}_{suffix}"


def inspect_artifacts(
    directory: Path,
    run_id: str,
    stages: dict[str, dict[str, Any]],
    label: str,
) -> dict[str, dict[str, dict[str, Any]]]:
    inventory: dict[str, dict[str, dict[str, Any]]] = {}
    for stage_name, stage in stages.items():
        stage_inventory: dict[str, dict[str, Any]] = {}
        expected_shape = (int(stage["H"]), int(stage["W"]))
        for suffix in REQUIRED_ARTIFACTS:
            path = _artifact_path(directory, run_id, stage_name, suffix)
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"{label}/{stage_name}: missing {suffix}")
            row: dict[str, Any] = {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            if suffix.endswith(".tif"):
                with tifffile.TiffFile(path) as tif:
                    series = tif.series[0]
                    shape = tuple(int(value) for value in series.shape)
                    dtype = str(series.dtype)
                if shape != expected_shape:
                    raise ValueError(
                        f"{label}/{stage_name}/{suffix}: shape {shape}, expected {expected_shape}"
                    )
                if dtype != "uint8":
                    raise ValueError(f"{label}/{stage_name}/{suffix}: expected uint8, got {dtype}")
                row.update({"shape": list(shape), "dtype": dtype})
            stage_inventory[suffix] = row
        inventory[stage_name] = stage_inventory
    return inventory


def _metric_signature(stats: dict[str, Any]) -> dict[str, Any]:
    stages = []
    for stage in stats["stages"]:
        stages.append({key: value for key, value in stage.items() if key != "seconds"})
    return {
        key: stats[key]
        for key in ("tool", "n_cubes", "nv", "nf", "n_seam_cols")
    } | {"stages": stages}


def _artifact_hash_signature(
    inventory: dict[str, dict[str, dict[str, Any]]]
) -> dict[str, dict[str, str]]:
    return {
        stage: {suffix: row["sha256"] for suffix, row in artifacts.items()}
        for stage, artifacts in inventory.items()
    }


def _stage_effect(parent: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
    parent_excess = max(0.0, float(parent["seam_ratio"]) - 1.0)
    candidate_excess = max(0.0, float(candidate["seam_ratio"]) - 1.0)
    parent_area = int(parent["W"]) * int(parent["H"])
    candidate_area = int(candidate["W"]) * int(candidate["H"])
    return {
        "parent": {
            key: parent[key]
            for key in (
                "W",
                "H",
                "fill",
                "multi_frac",
                "seam_dcol_mean",
                "base_dcol_mean",
                "seam_ratio",
            )
        },
        "candidate": {
            key: candidate[key]
            for key in (
                "W",
                "H",
                "fill",
                "multi_frac",
                "seam_dcol_mean",
                "base_dcol_mean",
                "seam_ratio",
            )
        },
        "change": {
            "fill_pp": 100.0 * (float(candidate["fill"]) - float(parent["fill"])),
            "multi_pp": 100.0
            * (float(candidate["multi_frac"]) - float(parent["multi_frac"])),
            "seam_dcol": float(candidate["seam_dcol_mean"])
            - float(parent["seam_dcol_mean"]),
            "seam_ratio": float(candidate["seam_ratio"]) - float(parent["seam_ratio"]),
            "seam_excess_reduction": parent_excess - candidate_excess,
            "canvas_area_ratio": candidate_area / parent_area,
        },
    }


def compare(
    parent_stats: dict[str, Any],
    candidate_stats: dict[str, Any],
    repeat_stats: dict[str, Any],
    parent_inventory: dict[str, dict[str, dict[str, Any]]],
    candidate_inventory: dict[str, dict[str, dict[str, Any]]],
    repeat_inventory: dict[str, dict[str, dict[str, Any]]],
    *,
    max_fill_drop_pp: float = 0.20,
    max_multi_increase_pp: float = 0.20,
    min_seam_excess_reduction: float = 0.05,
    max_final_seam_excess_increase: float = 0.05,
    min_canvas_ratio: float = 0.95,
    max_canvas_ratio: float = 1.05,
) -> dict[str, Any]:
    parent_stages = _stage_map(parent_stats, "parent")
    candidate_stages = _stage_map(candidate_stats, "candidate")
    repeat_stages = _stage_map(repeat_stats, "candidate_repeat")

    topology_fields = ("n_cubes", "nv", "nf")
    comparability = {
        key: parent_stats.get(key) == candidate_stats.get(key) == repeat_stats.get(key)
        for key in topology_fields
    }
    comparability["stage_names"] = (
        tuple(parent_stages) == tuple(candidate_stages) == tuple(repeat_stages)
    )
    comparability_pass = all(comparability.values())

    repeat_metrics_equal = _metric_signature(candidate_stats) == _metric_signature(repeat_stats)
    repeat_artifacts_equal = _artifact_hash_signature(candidate_inventory) == _artifact_hash_signature(
        repeat_inventory
    )
    repeat_pass = repeat_metrics_equal and repeat_artifacts_equal

    stage_effects = {
        stage: _stage_effect(parent_stages[stage], candidate_stages[stage])
        for stage in EXPECTED_STAGES
    }
    first = stage_effects["step1_ribbon"]
    final = stage_effects["step5_relax"]

    safety_gates = {
        "stage1_fill_drop_within_pp": first["change"]["fill_pp"] >= -max_fill_drop_pp,
        "final_fill_drop_within_pp": final["change"]["fill_pp"] >= -max_fill_drop_pp,
        "final_multi_increase_within_pp": final["change"]["multi_pp"]
        <= max_multi_increase_pp,
        "final_seam_excess_not_worse": final["change"]["seam_excess_reduction"]
        >= -max_final_seam_excess_increase,
        "stage1_canvas_bounded": min_canvas_ratio
        <= first["change"]["canvas_area_ratio"]
        <= max_canvas_ratio,
        "final_canvas_bounded": min_canvas_ratio
        <= final["change"]["canvas_area_ratio"]
        <= max_canvas_ratio,
    }
    safety_pass = all(safety_gates.values())
    efficacy_gates = {
        "stage1_seam_excess_reduction": first["change"]["seam_excess_reduction"]
        >= min_seam_excess_reduction,
        "stage1_raw_seam_nonincrease": first["change"]["seam_dcol"] <= 0.0,
    }
    efficacy_pass = all(efficacy_gates.values())

    if not comparability_pass or not repeat_pass:
        decision = "DOWNSTREAM_INCONCLUSIVE"
    elif not safety_pass:
        decision = "DOWNSTREAM_REGRESSION"
    elif efficacy_pass:
        decision = "DOWNSTREAM_EFFICACY_PASS"
    else:
        decision = "DOWNSTREAM_COMPATIBILITY_ONLY"

    return {
        "schema_version": 1,
        "decision": decision,
        "comparability": {"pass": comparability_pass, **comparability},
        "candidate_repeat": {
            "pass": repeat_pass,
            "metrics_equal_excluding_timing_and_paths": repeat_metrics_equal,
            "required_artifacts_byte_identical": repeat_artifacts_equal,
        },
        "thresholds": {
            "max_fill_drop_pp": max_fill_drop_pp,
            "max_multi_increase_pp": max_multi_increase_pp,
            "min_seam_excess_reduction": min_seam_excess_reduction,
            "max_final_seam_excess_increase": max_final_seam_excess_increase,
            "canvas_area_ratio": [min_canvas_ratio, max_canvas_ratio],
        },
        "safety_gates": {"pass": safety_pass, **safety_gates},
        "efficacy_gates": {"pass": efficacy_pass, **efficacy_gates},
        "stage_effects": stage_effects,
        "scope": (
            "Complete stage-1 and stage-5 readable strips; no seam, crop, or favorable "
            "window is selected. This is downstream texture evidence, not ink recovery."
        ),
    }


def render(
    result: dict[str, Any],
    parent_dir: Path,
    candidate_dir: Path,
    parent_id: str,
    candidate_id: str,
    output_prefix: Path,
) -> tuple[Path, Path]:
    figure, axes = plt.subplots(2, 2, figsize=(18, 14), constrained_layout=True)
    arms = (("Exact parent", parent_dir, parent_id), ("PR candidate", candidate_dir, candidate_id))
    for row_index, stage in enumerate(DISPLAY_STAGES):
        for column_index, (arm_label, directory, run_id) in enumerate(arms):
            path = _artifact_path(directory, run_id, stage, "rawtex_strip.png")
            image = mpimg.imread(path)
            axis = axes[row_index, column_index]
            axis.imshow(image, interpolation="nearest")
            axis.set_axis_off()
            metrics = result["stage_effects"][stage][
                "parent" if column_index == 0 else "candidate"
            ]
            axis.set_title(
                f"{arm_label} - {stage}\n"
                f"fill={100.0 * float(metrics['fill']):.2f}% | "
                f"seam ratio={float(metrics['seam_ratio']):.3f}",
                fontsize=12,
            )
    figure.suptitle(
        "PHerc0826 complete held-out unroll strips\n"
        "Every emitted strip tile is shown; no seam or crop selected | "
        f"decision: {result['decision']}",
        fontsize=16,
    )
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    png = output_prefix.with_suffix(".png")
    svg = output_prefix.with_suffix(".svg")
    figure.savefig(png, dpi=160, facecolor="white")
    figure.savefig(svg, facecolor="white")
    plt.close(figure)
    return png, svg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for arm in ("parent", "candidate", "candidate-repeat"):
        name = arm.replace("-", "_")
        parser.add_argument(f"--{arm}-dir", dest=f"{name}_dir", type=Path, required=True)
        parser.add_argument(f"--{arm}-id", dest=f"{name}_id", required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--figure-prefix", type=Path, required=True)
    parser.add_argument("--max-fill-drop-pp", type=float, default=0.20)
    parser.add_argument("--max-multi-increase-pp", type=float, default=0.20)
    parser.add_argument("--min-seam-excess-reduction", type=float, default=0.05)
    parser.add_argument("--max-final-seam-excess-increase", type=float, default=0.05)
    parser.add_argument("--min-canvas-ratio", type=float, default=0.95)
    parser.add_argument("--max-canvas-ratio", type=float, default=1.05)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    arms = {}
    for name in ("parent", "candidate", "candidate_repeat"):
        directory: Path = getattr(args, f"{name}_dir")
        run_id: str = getattr(args, f"{name}_id")
        stats_path = directory / f"{run_id}_pipeline_stats.json"
        stats = load_json(stats_path)
        stages = _stage_map(stats, name)
        inventory = inspect_artifacts(directory, run_id, stages, name)
        arms[name] = {
            "directory": directory,
            "run_id": run_id,
            "stats_path": stats_path,
            "stats": stats,
            "inventory": inventory,
        }

    result = compare(
        arms["parent"]["stats"],
        arms["candidate"]["stats"],
        arms["candidate_repeat"]["stats"],
        arms["parent"]["inventory"],
        arms["candidate"]["inventory"],
        arms["candidate_repeat"]["inventory"],
        max_fill_drop_pp=args.max_fill_drop_pp,
        max_multi_increase_pp=args.max_multi_increase_pp,
        min_seam_excess_reduction=args.min_seam_excess_reduction,
        max_final_seam_excess_increase=args.max_final_seam_excess_increase,
        min_canvas_ratio=args.min_canvas_ratio,
        max_canvas_ratio=args.max_canvas_ratio,
    )
    png, svg = render(
        result,
        arms["parent"]["directory"],
        arms["candidate"]["directory"],
        arms["parent"]["run_id"],
        arms["candidate"]["run_id"],
        args.figure_prefix,
    )
    result["evidence"] = {
        "pipeline_stats": {
            name: {
                "path": str(arm["stats_path"]),
                "sha256": sha256(arm["stats_path"]),
            }
            for name, arm in arms.items()
        },
        "required_artifacts": {
            name: arm["inventory"] for name, arm in arms.items()
        },
        "figure_png": {"path": str(png), "sha256": sha256(png)},
        "figure_svg": {"path": str(svg), "sha256": sha256(svg)},
        "visual_scope": (
            "Complete stage-1 and stage-5 readable-strip products are rendered; "
            "no seam, crop, or favorable window is selected."
        ),
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"decision": result["decision"], "out": str(args.out_json)}, sort_keys=True))


if __name__ == "__main__":
    main()
