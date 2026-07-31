#!/usr/bin/env python3
"""Replay a captured ScrollFiesta MLS cloud and compare it to a C result.

The step0 diagnostic OBJs preserve vertex order, so this converts their world
ZYX coordinates back to cube-local coordinates, runs the vendored Rust/CubeCL
benchmark, and records both timing and position parity as JSON.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_BENCH = (REPO / "deps" / "scrollfiesta-mls-cubecl" / "target"
                 / "release" / "deps" / "mls_real_bench.exe")


def obj_vertices(path: Path) -> np.ndarray:
    rows: list[tuple[float, float, float]] = []
    with path.open("rb") as handle:
        for line in handle:
            if line.startswith(b"v "):
                fields = line.split()
                rows.append((float(fields[1]), float(fields[2]), float(fields[3])))
    result = np.asarray(rows, dtype=np.float32).reshape((-1, 3))
    if not len(result):
        raise ValueError(f"no OBJ vertices in {path}")
    return result


def infer_origin(path: Path) -> tuple[float, float, float] | None:
    match = re.search(r"z(\d+)_y(\d+)_x(\d+)", str(path))
    if not match:
        return None
    return tuple(float(value) for value in match.groups())


def native_env(repetitions: int) -> dict[str, str]:
    env = dict(os.environ)
    env["MLS_BENCH_REPS"] = str(repetitions)
    local_cuda = REPO / ".toolchains" / "cuda-13.2-wheels" / "nvidia" / "cu13"
    configured = env.get("CUDA_PATH")
    cuda_root = Path(configured) if configured else local_cuda
    if (cuda_root / "include" / "cccl").exists():
        env["CUDA_PATH"] = str(cuda_root)
        old_path = ""
        for key in list(env):
            if key.upper() == "PATH":
                old_path = old_path or env[key]
                del env[key]
        bins = [cuda_root / "bin" / "x86_64", cuda_root / "bin"]
        env["PATH"] = os.pathsep.join(
            [str(path) for path in bins if path.exists()]
            + ([old_path] if old_path else []))
    return env


def scalar_output(stdout: str) -> dict[str, object]:
    values: dict[str, object] = {}
    for line in stdout.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        value = value.strip()
        if value.lower() in {"true", "false"}:
            values[key] = value.lower() == "true"
            continue
        try:
            values[key] = float(value) if any(ch in value for ch in ".eE") else int(value)
        except ValueError:
            values[key] = value
    return values


def comparison(reference: np.ndarray, actual: np.ndarray) -> dict[str, object]:
    if reference.shape != actual.shape:
        raise ValueError(f"vertex-count mismatch: reference {reference.shape}, "
                         f"candidate {actual.shape}")
    distances = np.linalg.norm(
        reference.astype(np.float64) - actual.astype(np.float64), axis=1)
    return {
        "count": int(len(distances)),
        "max_position_error_vox": float(np.max(distances)),
        "rms_position_error_vox": float(np.sqrt(np.mean(distances ** 2))),
        "mean_position_error_vox": float(np.mean(distances)),
        "p99_position_error_vox": float(np.quantile(distances, 0.99)),
        "outliers_gt_0p00124": int(np.count_nonzero(distances > 0.00124)),
        "outliers_gt_0p25": int(np.count_nonzero(distances > 0.25)),
        "strict_position_parity_pass": bool(np.max(distances) <= 0.00124),
        "weld_safety_pass": bool(np.max(distances) < 0.25),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pre-obj", required=True, type=Path)
    parser.add_argument("--reference-post-obj", type=Path)
    parser.add_argument("--bench", type=Path, default=DEFAULT_BENCH)
    parser.add_argument("--backend", default="cubecl-cuda",
                        choices=("rust-cpu", "cubecl-cpu", "cubecl-hip",
                                 "cubecl-cuda"))
    parser.add_argument("--passes", type=int, default=5)
    parser.add_argument("--radius", type=float, default=12.0)
    parser.add_argument("--origin", type=float, nargs=3,
                        metavar=("Z", "Y", "X"))
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--out-prefix", required=True, type=Path)
    args = parser.parse_args()

    if not args.bench.is_file():
        parser.error(f"benchmark binary not found: {args.bench}")
    origin_tuple = tuple(args.origin) if args.origin else infer_origin(args.pre_obj)
    if origin_tuple is None:
        parser.error("could not infer z/y/x cube origin; pass --origin")
    if args.passes <= 0 or args.repetitions <= 0 or args.radius <= 0:
        parser.error("passes, repetitions, and radius must be positive")

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    origin = np.asarray(origin_tuple, dtype=np.float32)
    pre = obj_vertices(args.pre_obj) - origin
    input_f32 = Path(f"{args.out_prefix}.pre.positions_zyx.f32")
    pre.astype("<f4", copy=False).tofile(input_f32)

    command = [
        str(args.bench), args.backend, str(input_f32), str(args.out_prefix),
        str(args.passes), str(args.radius),
        *(str(value) for value in origin_tuple),
    ]
    proc = subprocess.run(command, capture_output=True, text=True,
                          env=native_env(args.repetitions))
    if proc.returncode != 0:
        print(proc.stdout, end="")
        print(proc.stderr, file=sys.stderr, end="")
        return proc.returncode

    generated = Path(f"{args.out_prefix}.cubecl.positions_zyx.f32")
    actual = np.fromfile(generated, dtype="<f4").reshape((-1, 3))
    report: dict[str, object] = {
        "schema": "scrollfiesta.mls-capture-verification.v1",
        "backend": args.backend,
        "passes": args.passes,
        "radius_vox": args.radius,
        "origin_zyx": list(origin_tuple),
        "pre_obj": str(args.pre_obj.resolve()),
        "reference_post_obj": (str(args.reference_post_obj.resolve())
                               if args.reference_post_obj else None),
        "benchmark": scalar_output(proc.stdout),
    }
    exit_code = 0
    if args.reference_post_obj:
        reference = obj_vertices(args.reference_post_obj) - origin
        report["c_reference_vs_backend"] = comparison(reference, actual)
        exit_code = 0 if report["c_reference_vs_backend"][
            "strict_position_parity_pass"] else 2

    report_path = Path(f"{args.out_prefix}.verification.json")
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))
    print(f"report={report_path}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
