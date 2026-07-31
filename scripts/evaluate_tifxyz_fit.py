#!/usr/bin/env python3
"""Measure deterministic geometric agreement between TIFXYZ surfaces.

The comparison is deliberately independent of Spiral's stochastic training
loss.  Each candidate and the reference are compared in both directions with
nearest-neighbour distances in voxel coordinates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import tifffile
from scipy.spatial import cKDTree


def _surface_digest(directory: Path) -> str:
    digest = hashlib.sha256()
    for name in ("meta.json", "x.tif", "y.tif", "z.tif"):
        path = directory / name
        if not path.is_file():
            continue
        digest.update(name.encode("utf-8"))
        with path.open("rb") as stream:
            while block := stream.read(1024 * 1024):
                digest.update(block)
    return digest.hexdigest()


def _load_points(directory: Path) -> tuple[np.ndarray, dict]:
    arrays = [
        np.asarray(tifffile.imread(directory / f"{axis}.tif"), dtype=np.float64)
        for axis in "xyz"
    ]
    if arrays[0].shape != arrays[1].shape or arrays[0].shape != arrays[2].shape:
        raise ValueError(f"coordinate shapes differ in {directory}")
    valid = np.logical_and.reduce(
        [np.isfinite(array) & (array >= 0.0) for array in arrays]
    )
    points = np.column_stack([array[valid] for array in arrays])
    if not len(points):
        raise ValueError(f"no valid coordinates in {directory}")
    info = {
        "path": str(directory.resolve()),
        "sha256": _surface_digest(directory),
        "grid_shape": list(arrays[0].shape),
        "valid_points": int(len(points)),
        "valid_fraction": float(valid.mean()),
        "bbox_xyz": [
            [float(value) for value in points.min(axis=0)],
            [float(value) for value in points.max(axis=0)],
        ],
    }
    return points, info


def _sample(points: np.ndarray, maximum: int, seed: int) -> np.ndarray:
    if len(points) <= maximum:
        return points
    rng = np.random.default_rng(seed)
    indices = np.sort(rng.choice(len(points), size=maximum, replace=False))
    return points[indices]


def _distance_summary(distances: np.ndarray) -> dict:
    percentiles = np.percentile(distances, [50, 90, 95, 99])
    return {
        "samples": int(len(distances)),
        "mean_vx": float(distances.mean()),
        "rms_vx": float(np.sqrt(np.mean(np.square(distances)))),
        "p50_vx": float(percentiles[0]),
        "p90_vx": float(percentiles[1]),
        "p95_vx": float(percentiles[2]),
        "p99_vx": float(percentiles[3]),
        "max_vx": float(distances.max()),
        "within_1vx_fraction": float(np.mean(distances <= 1.0)),
        "within_2vx_fraction": float(np.mean(distances <= 2.0)),
        "within_5vx_fraction": float(np.mean(distances <= 5.0)),
        "within_10vx_fraction": float(np.mean(distances <= 10.0)),
        "within_20vx_fraction": float(np.mean(distances <= 20.0)),
    }


def _directed(source: np.ndarray, target: np.ndarray) -> dict:
    distances, _ = cKDTree(target).query(source, workers=-1)
    return _distance_summary(np.asarray(distances))


def _candidate(value: str) -> tuple[str, Path]:
    label, separator, path = value.partition("=")
    if not separator or not label or not path:
        raise argparse.ArgumentTypeError("candidate must be LABEL=PATH")
    return label, Path(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument(
        "--candidate",
        action="append",
        required=True,
        type=_candidate,
        metavar="LABEL=PATH",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-points", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=137)
    args = parser.parse_args()
    if args.max_points <= 0:
        parser.error("--max-points must be positive")

    reference, reference_info = _load_points(args.reference)
    sampled_reference = _sample(reference, args.max_points, args.seed)
    report = {
        "schema": "tifxyz-geometric-fit-v1",
        "units": "voxels",
        "seed": args.seed,
        "max_points_per_surface": args.max_points,
        "reference": reference_info,
        "candidates": [],
    }
    for index, (label, path) in enumerate(args.candidate):
        candidate, candidate_info = _load_points(path)
        sampled_candidate = _sample(
            candidate, args.max_points, args.seed + index + 1
        )
        candidate_to_reference = _directed(sampled_candidate, sampled_reference)
        reference_to_candidate = _directed(sampled_reference, sampled_candidate)
        symmetric_rms = np.sqrt(
            (
                candidate_to_reference["rms_vx"] ** 2
                + reference_to_candidate["rms_vx"] ** 2
            )
            / 2.0
        )
        report["candidates"].append(
            {
                "label": label,
                "surface": candidate_info,
                "candidate_to_reference": candidate_to_reference,
                "reference_to_candidate": reference_to_candidate,
                "symmetric_rms_vx": float(symmetric_rms),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
