"""Semantic audit and candidate selection for adaptive BPA rerolls.

The ordinary winding audit catches a one-edge jump between turns.  Its local
phase-span mode also catches a broad BPA fusion wall assembled from many short
edges.  This module runs both checks at the BPA-bearing cube stages and applies
the deadline-safe selection policy:

* try rho 1.20 first, then 1.10 and 1.00 only when the prior result is fused;
* never accept a candidate that loses more than five percent of baseline points;
* never accept one that grows the initial-BPA boundary by more than 25 percent;
* prefer the highest-radius candidate that is clean at every audited stage;
* keep the baseline if no candidate passes, and record the failure explicitly.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


CRITICAL_STAGES = (
    ("step1_bpa", "initial_bpa"),
    ("step4_bpa", "post_split_bpa"),
    ("step7_cc_bpa", "post_component_bpa"),
    ("step12_final", "final_cube"),
)


@dataclass(frozen=True)
class AdaptiveBpaConfig:
    wind_audit_bin: Path
    umb_y: float
    umb_x: float
    pitch: float = 9.5
    rhos: tuple[float, ...] = (1.20, 1.10, 1.00)
    local_span_tol: float = 1.5
    min_faces: int = 64
    point_tolerance: float = 1.5
    min_point_coverage: float = 0.95
    max_boundary_ratio: float = 1.25
    keep_candidates: bool = False


def _stage_obj(cube_dir: Path, token: str) -> Path | None:
    """Return the aggregate OBJ for one named dump stage, if it exists."""
    dirs = sorted(cube_dir.glob(f"*_{token}"))
    for stage_dir in dirs:
        exact = stage_dir / f"{stage_dir.name}_all.obj"
        if exact.is_file():
            return exact
        all_objs = sorted(stage_dir.glob("*_all.obj"))
        if all_objs:
            return all_objs[0]
    return None


def final_obj(cube_dir: Path) -> Path | None:
    return _stage_obj(Path(cube_dir), "step12_final")


def initial_bpa_obj(cube_dir: Path) -> Path | None:
    return _stage_obj(Path(cube_dir), "step1_bpa")


def _audit_one(obj: Path, out_json: Path, cfg: AdaptiveBpaConfig) -> dict:
    out_json.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(cfg.wind_audit_bin), str(obj),
        "--umb-y", str(cfg.umb_y), "--umb-x", str(cfg.umb_x),
        "--pitch", str(cfg.pitch),
        "--local-span-tol", str(cfg.local_span_tol),
        "--min-faces", str(cfg.min_faces),
        "--top", "8", "--json", str(out_json),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180.0)
    if proc.returncode != 0 or not out_json.is_file():
        raise RuntimeError(
            f"wind_audit failed (rc={proc.returncode}) for {obj}\n"
            f"stderr: {proc.stderr[-1000:]}")
    result = json.loads(out_json.read_text())
    result["obj"] = str(obj)
    return result


def audit_candidate(cube_dir: Path, cfg: AdaptiveBpaConfig) -> dict:
    """Audit every topology-critical stage produced for one cube run."""
    cube_dir = Path(cube_dir)
    stages: list[dict] = []
    missing: list[str] = []
    for token, classification in CRITICAL_STAGES:
        obj = _stage_obj(cube_dir, token)
        if obj is None:
            missing.append(token)
            continue
        metrics = _audit_one(obj, cube_dir / "wind_audit" / f"{token}.json", cfg)
        metrics["stage"] = token
        metrics["classification"] = classification
        metrics["fused"] = bool(
            metrics.get("merge_ft", 0)
            or metrics.get("broad_fusion_components", 0)
        )
        stages.append(metrics)

    by_stage = {s["stage"]: s for s in stages}
    required_present = "step1_bpa" in by_stage and "step12_final" in by_stage
    first_failure = next((s["classification"] for s in stages if s["fused"]), None)
    return {
        "stages": stages,
        "missing_stages": missing,
        "first_failure": first_failure,
        "topology_clean": required_present and not any(s["fused"] for s in stages),
    }


def load_obj_vertices(path: Path) -> np.ndarray:
    """Load only OBJ vertex positions, avoiding a full trimesh parse."""
    rows: list[tuple[float, float, float]] = []
    with open(path, "rb") as fh:
        for raw in fh:
            if not raw.startswith(b"v "):
                continue
            fields = raw.split()
            if len(fields) >= 4:
                rows.append((float(fields[1]), float(fields[2]), float(fields[3])))
    return np.asarray(rows, dtype=np.float32).reshape((-1, 3))


def point_coverage(reference_obj: Path, candidate_obj: Path, tolerance: float) -> float:
    """Fraction of baseline BPA vertices represented within ``tolerance``."""
    reference = load_obj_vertices(reference_obj)
    candidate = load_obj_vertices(candidate_obj)
    if len(reference) == 0 or len(candidate) == 0:
        return 0.0
    tree = cKDTree(candidate)
    distances, _ = tree.query(reference, k=1, distance_upper_bound=tolerance,
                              workers=-1)
    return float(np.count_nonzero(np.isfinite(distances)) / len(reference))


def stage_metrics(evaluation: dict, stage: str) -> dict | None:
    return next((s for s in evaluation.get("stages", [])
                 if s.get("stage") == stage), None)


def apply_acceptance(evaluation: dict, baseline: dict | None,
                     cfg: AdaptiveBpaConfig) -> dict:
    """Attach coverage/boundary acceptance metrics to an audited candidate."""
    current_initial = stage_metrics(evaluation, "step1_bpa")
    if baseline is None:
        coverage = 1.0
        boundary_ratio = 1.0
    else:
        base_initial = stage_metrics(baseline, "step1_bpa")
        if not base_initial or not current_initial:
            coverage = 0.0
            boundary_ratio = None
        else:
            coverage = point_coverage(
                Path(base_initial["obj"]), Path(current_initial["obj"]),
                cfg.point_tolerance,
            )
            base_boundary = int(base_initial.get("boundary_edges", 0))
            cur_boundary = int(current_initial.get("boundary_edges", 0))
            boundary_ratio = (cur_boundary / base_boundary
                              if base_boundary > 0 else (1.0 if cur_boundary == 0 else None))
    evaluation["point_coverage"] = coverage
    evaluation["boundary_ratio"] = boundary_ratio
    evaluation["coverage_ok"] = coverage >= cfg.min_point_coverage
    evaluation["boundary_ok"] = (boundary_ratio is not None
                                 and boundary_ratio <= cfg.max_boundary_ratio)
    evaluation["process_ok"] = int(evaluation.get("returncode", 0)) == 0
    evaluation["passes"] = bool(
        evaluation["process_ok"]
        and evaluation.get("topology_clean")
        and evaluation["coverage_ok"]
        and evaluation["boundary_ok"]
    )
    return evaluation


def choose_candidate(candidates: list[dict]) -> tuple[dict, bool]:
    """Return first passing (therefore highest-rho) candidate, else baseline."""
    for candidate in candidates:
        if candidate.get("passes"):
            return candidate, True
    available = [c for c in candidates if c.get("produced_obj")]
    if not available:
        raise RuntimeError("adaptive BPA produced no candidate OBJ")
    return available[0], False


def diagnose_failure(candidates: list[dict]) -> str | None:
    """Classify why an adaptive radius sweep did not produce a clean mesh.

    A stable broad-fusion signature across two or more successful radii is
    evidence that simply shrinking BPA again is unlikely to help.  We call it
    ``radius_insensitive_fusion`` rather than claiming the source cloud is at
    fault: ``step1_bpa`` is still post-BPA, so the audit cannot prove where the
    connection originated.
    """
    if any(candidate.get("passes") for candidate in candidates):
        return None
    if not any(candidate.get("produced_obj") for candidate in candidates):
        return "no_candidate_output"
    if any(int(candidate.get("returncode", 0)) != 0 for candidate in candidates):
        return "candidate_execution_failure"

    initial: list[dict] = []
    for candidate in candidates:
        stage = stage_metrics(candidate, "step1_bpa")
        if stage is not None:
            initial.append(stage)
    if len(initial) >= 2 and all(stage.get("fused") for stage in initial):
        broad_counts = [int(stage.get("broad_fusion_components", 0))
                        for stage in initial]
        spans = [float(stage.get("broad_fusion_max_span", 0.0))
                 for stage in initial]
        span_tolerance = max(0.25, 0.05 * max(spans, default=0.0))
        stable_broad_signature = (
            broad_counts[0] > 0
            and len(set(broad_counts)) == 1
            and max(spans) - min(spans) <= span_tolerance
        )
        if stable_broad_signature:
            return "radius_insensitive_fusion"
    return "no_clean_radius"
