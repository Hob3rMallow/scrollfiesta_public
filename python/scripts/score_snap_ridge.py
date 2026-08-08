"""Controlled same-vertex CT-ridge benchmark for ScrollFiesta snap arms.

The benchmark fixes the sampling normal to the pre-snap mesh, uses cube-level
summaries, and exposes only one spatial split per invocation. OBJ inputs must
come from ``scroll_unroll --export-mesh`` and retain identical shared indices.

Examples
--------
Prepare the sigma-1 CT cache once::

  python python/scripts/score_snap_ridge.py prepare --raw-dir GRID/cubes_RAW

Score development arms without touching the held-out z band::

  python python/scripts/score_snap_ridge.py score --split dev --pre pre.obj \
    --arm iter0=iter0.obj --arm iter4=iter4.obj --production iter4 --out dev.json

Reveal the held-out band only after naming the selected candidate::

  python python/scripts/score_snap_ridge.py score --split test --pre pre.obj \
    --arm candidate=candidate.obj --arm production=iter4.obj \
    --production production --candidate candidate --out test.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import tifffile
from scipy.ndimage import gaussian_filter, map_coordinates


ORIGIN = np.array([4352.0, 3072.0, 2560.0], dtype=np.float64)
SHAPE = (512, 640, 640)
CHUNK = 128
SIGMA = 1.0
HALF = 4.0
STEP = 0.25
MARGIN = HALF + 2.0
BOOTSTRAP_SEED = 20260808
BOOTSTRAP_N = 10_000
CUBE_RE = re.compile(
    r"^# cube (z\d{5}_y\d{5}_x\d{5}) vertex_range_0based \[(\d+),(\d+)\)$"
)
TIFF_RE = re.compile(r"^z(-?\d+)_y(-?\d+)_x(-?\d+)\.tif$")


def sha256(path: Path, block: int = 8 << 20) -> str:
    """Return the lowercase SHA-256 digest of *path*."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while data := handle.read(block):
            digest.update(data)
    return digest.hexdigest()


def prepare_cache(raw_dir: Path, cache: Path) -> dict:
    """Assemble and sigma-1 smooth the frozen 4x5x5 PHerc0139 RAW grid."""
    files = sorted(raw_dir.glob("z*_y*_x*.tif"))
    if len(files) != 100:
        raise ValueError(f"expected 100 RAW cubes, found {len(files)} in {raw_dir}")
    volume = np.empty(SHAPE, dtype=np.uint8)
    occupied = np.zeros((4, 5, 5), dtype=bool)
    names: list[str] = []
    for path in files:
        match = TIFF_RE.match(path.name)
        if not match:
            raise ValueError(f"bad cube filename: {path.name}")
        world = np.array([int(value) for value in match.groups()], dtype=np.int64)
        index = ((world - ORIGIN.astype(np.int64)) // CHUNK).astype(np.int64)
        if np.any(index < 0) or np.any(index >= np.array(occupied.shape)):
            raise ValueError(f"cube outside registered bbox: {path.name}")
        if occupied[tuple(index)]:
            raise ValueError(f"duplicate cube index: {path.name}")
        cube = np.asarray(tifffile.imread(path))
        if cube.shape != (CHUNK,) * 3 or cube.dtype != np.uint8:
            raise ValueError(
                f"bad cube shape/dtype: {path.name} {cube.shape} {cube.dtype}"
            )
        slices = tuple(
            slice(int(i * CHUNK), int((i + 1) * CHUNK)) for i in index
        )
        volume[slices] = cube
        occupied[tuple(index)] = True
        names.append(path.name)
    if not occupied.all():
        raise ValueError("RAW cube grid has holes")

    cache.parent.mkdir(parents=True, exist_ok=True)
    smooth = np.lib.format.open_memmap(
        cache, mode="w+", dtype=np.float32, shape=SHAPE
    )
    gaussian_filter(volume, SIGMA, output=smooth, mode="nearest")
    smooth.flush()
    del smooth
    manifest = {
        "raw_dir": raw_dir.as_posix(),
        "origin_zyx": ORIGIN.astype(int).tolist(),
        "shape_zyx": list(SHAPE),
        "chunk": CHUNK,
        "sigma": SIGMA,
        "n_cubes": len(files),
        "cube_names_sha256": hashlib.sha256(
            "\n".join(names).encode()
        ).hexdigest(),
        "cache": cache.as_posix(),
        "cache_sha256": sha256(cache),
    }
    metadata_path = cache.with_suffix(cache.suffix + ".json")
    metadata_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def load_obj(path: Path) -> dict:
    """Load and validate a shared-index OBJ emitted by ``--export-mesh``."""
    verts: list[list[float]] = []
    uv: list[list[float]] = []
    normals: list[list[float]] = []
    faces: list[list[int]] = []
    cubes: list[tuple[str, int, int]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("v "):
                verts.append([float(value) for value in line.split()[1:4]])
            elif line.startswith("vt "):
                uv.append([float(value) for value in line.split()[1:3]])
            elif line.startswith("vn "):
                normals.append([float(value) for value in line.split()[1:4]])
            elif line.startswith("f "):
                faces.append(
                    [int(value.split("/", 1)[0]) - 1 for value in line.split()[1:4]]
                )
            elif line.startswith("# cube "):
                match = CUBE_RE.match(line.rstrip("\n"))
                if not match:
                    raise ValueError(
                        f"malformed cube range in {path}: {line.rstrip()}"
                    )
                cubes.append(
                    (match.group(1), int(match.group(2)), int(match.group(3)))
                )
    out = {
        "path": path.as_posix(),
        "sha256": sha256(path),
        "verts": np.asarray(verts, dtype=np.float32),
        "uv": np.asarray(uv, dtype=np.float32),
        "normals": np.asarray(normals, dtype=np.float32),
        "faces": np.asarray(faces, dtype=np.int32),
        "cubes": cubes,
    }
    nv = len(out["verts"])
    if out["verts"].shape != (nv, 3) or out["uv"].shape != (nv, 2):
        raise ValueError(f"incomplete v/vt records in {path}")
    if out["normals"].shape != (nv, 3) or out["faces"].ndim != 2:
        raise ValueError(f"incomplete vn/f records in {path}")
    if out["faces"].shape[1:] != (3,):
        raise ValueError(f"non-triangular face records in {path}")
    if len(out["faces"]) and (
        out["faces"].min() < 0 or out["faces"].max() >= nv
    ):
        raise ValueError(f"face index outside vertex range in {path}")
    if not cubes or cubes[0][1] != 0 or cubes[-1][2] != nv:
        raise ValueError(f"missing/non-covering cube ranges in {path}")
    for (name, lo, hi), next_cube in zip(
        cubes, cubes[1:] + [("END", nv, nv)]
    ):
        if not (0 <= lo < hi <= nv) or hi != next_cube[1]:
            raise ValueError(f"non-contiguous cube ranges at {name} in {path}")
    return out


def parse_arm(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("arm must be NAME=OBJ")
    name, path = value.split("=", 1)
    if not name or not path:
        raise argparse.ArgumentTypeError("arm must be NAME=OBJ")
    return name, Path(path)


def ridge_offsets(
    smooth: np.ndarray, pts_world: np.ndarray, normals: np.ndarray
) -> np.ndarray:
    """Return parabolically refined intensity-ridge offsets along fixed normals."""
    ts = np.arange(-HALF, HALF + 1e-8, STEP, dtype=np.float32)
    points = pts_world.astype(np.float32) - ORIGIN.astype(np.float32)
    coords = (
        points[:, :, None]
        + normals[:, :, None].astype(np.float32) * ts[None, None, :]
    )
    profile = map_coordinates(
        smooth,
        coords.transpose(1, 0, 2).reshape(3, -1),
        order=1,
        mode="constant",
        cval=np.nan,
        prefilter=False,
    ).reshape(len(points), len(ts))
    finite = np.isfinite(profile).all(axis=1)
    safe = np.where(np.isfinite(profile), profile, -np.inf)
    peak = np.argmax(safe, axis=1)
    out = ts[peak].astype(np.float64)
    valid = finite & (peak > 0) & (peak + 1 < len(ts))
    rows = np.flatnonzero(valid)
    if len(rows):
        k = peak[rows]
        y0 = profile[rows, k - 1]
        y1 = profile[rows, k]
        y2 = profile[rows, k + 1]
        denominator = y0 - 2.0 * y1 + y2
        shift = np.zeros(len(rows), dtype=np.float64)
        nonzero = np.abs(denominator) > 1e-9
        shift[nonzero] = (
            0.5 * (y0[nonzero] - y2[nonzero]) / denominator[nonzero]
        )
        out[rows] = ts[k] + np.clip(shift, -1.0, 1.0) * STEP
    out[~valid] = np.nan
    return out


def percentiles(values: np.ndarray) -> dict:
    if not len(values):
        return {"p50": None, "p95": None, "p99": None, "max": None}
    quantiles = np.percentile(values, [50, 95, 99, 100])
    return {
        "p50": float(quantiles[0]),
        "p95": float(quantiles[1]),
        "p99": float(quantiles[2]),
        "max": float(quantiles[3]),
    }


def bootstrap_median_ci(values: np.ndarray) -> list[float]:
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    n_values = len(values)
    if n_values == 0:
        return [float("nan"), float("nan")]
    draws = rng.integers(0, n_values, size=(BOOTSTRAP_N, n_values))
    medians = np.median(values[draws], axis=1)
    return [float(value) for value in np.percentile(medians, [2.5, 97.5])]


def _validated_cache_metadata(cache: Path) -> dict:
    metadata_path = cache.with_suffix(cache.suffix + ".json")
    if not cache.exists() or not metadata_path.exists():
        raise FileNotFoundError("run the prepare command first")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    expected = {
        "origin_zyx": ORIGIN.astype(int).tolist(),
        "shape_zyx": list(SHAPE),
        "chunk": CHUNK,
        "sigma": SIGMA,
        "n_cubes": 100,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(f"cache metadata {key!r} does not match frozen protocol")
    actual_sha256 = sha256(cache)
    if metadata.get("cache_sha256") != actual_sha256:
        raise ValueError("cache SHA-256 does not match its metadata")
    return {
        key: metadata[key]
        for key in (
            "origin_zyx",
            "shape_zyx",
            "chunk",
            "sigma",
            "n_cubes",
            "cube_names_sha256",
            "cache_sha256",
        )
    }


def score(args: argparse.Namespace) -> dict:
    """Score one locked split and write the complete JSON result."""
    if args.split == "test" and not args.candidate:
        raise ValueError(
            "test split is locked: pass --candidate NAME after development selection"
        )
    cache = Path(args.cache)
    cache_metadata = _validated_cache_metadata(cache)
    smooth = np.load(cache, mmap_mode="r")

    pre = load_obj(Path(args.pre))
    arm_pairs = list(args.arm)
    arm_paths = dict(arm_pairs)
    if len(arm_paths) != len(arm_pairs):
        raise ValueError("duplicate arm name")
    if args.production not in arm_paths:
        raise ValueError(f"production arm {args.production!r} was not supplied")
    if args.candidate and args.candidate not in arm_paths:
        raise ValueError(f"candidate arm {args.candidate!r} was not supplied")
    arms = {name: load_obj(path) for name, path in arm_paths.items()}
    for name, arm in arms.items():
        if arm["verts"].shape != pre["verts"].shape:
            raise ValueError(f"{name}: vertex count differs from pre")
        if not np.array_equal(arm["faces"], pre["faces"]):
            raise ValueError(f"{name}: topology/face indices differ from pre")
        if arm["cubes"] != pre["cubes"]:
            raise ValueError(f"{name}: cube vertex ranges differ from pre")
        if not np.array_equal(arm["uv"], pre["uv"]):
            raise ValueError(f"{name}: registered UV changed")

    used = np.zeros(len(pre["verts"]), dtype=bool)
    used[pre["faces"].ravel()] = True
    normals = pre["normals"].astype(np.float64)
    normal_length = np.linalg.norm(normals, axis=1)
    eligible = (
        used
        & np.isfinite(normals).all(axis=1)
        & (normal_length > 0.99)
        & (normal_length < 1.01)
    )
    all_positions = [pre["verts"]] + [arm["verts"] for arm in arms.values()]
    low = ORIGIN + MARGIN
    high = ORIGIN + np.array(SHAPE) - 1.0 - MARGIN
    for positions in all_positions:
        eligible &= np.all((positions >= low) & (positions <= high), axis=1)

    split_vertex = np.zeros(len(pre["verts"]), dtype=bool)
    cube_id = np.full(len(pre["verts"]), -1, dtype=np.int32)
    selected_cubes: list[tuple[int, str, int, int]] = []
    for cube_index, (name, lo, hi) in enumerate(pre["cubes"]):
        z_origin = int(name[1:6])
        wanted = (
            args.split == "dev" and z_origin in (4352, 4480, 4608)
        ) or (args.split == "test" and z_origin == 4736)
        cube_id[lo:hi] = cube_index
        if wanted:
            split_vertex[lo:hi] = True
            selected_cubes.append((cube_index, name, lo, hi))
    eligible &= split_vertex

    indices = np.flatnonzero(eligible)
    fixed_normals = normals[indices].astype(np.float32)
    offsets = {
        "pre": ridge_offsets(smooth, pre["verts"][indices], fixed_normals)
    }
    for name, arm in arms.items():
        offsets[name] = ridge_offsets(
            smooth, arm["verts"][indices], fixed_normals
        )
    common = np.ones(len(indices), dtype=bool)
    for values in offsets.values():
        common &= np.isfinite(values)
    common_indices = indices[common]

    per_cube: list[dict] = []
    for cube_index, name, _lo, _hi in selected_cubes:
        mask = cube_id[common_indices] == cube_index
        row: dict[str, object] = {"cube": name, "n": int(mask.sum())}
        for arm_name, values in offsets.items():
            row[arm_name] = (
                float(np.median(np.abs(values[common][mask])))
                if mask.any()
                else None
            )
        per_cube.append(row)
    if any(row["n"] == 0 for row in per_cube):
        raise ValueError("one or more selected cubes has zero common valid vertices")

    summaries: dict[str, dict] = {}
    production = np.array(
        [row[args.production] for row in per_cube], dtype=float
    )
    for name in ["pre", *arms.keys()]:
        cube_scores = np.array([row[name] for row in per_cube], dtype=float)
        item: dict[str, object] = {
            "median_cube_abs_ridge": float(np.median(cube_scores)),
            "cube_p25_p75": [
                float(value) for value in np.percentile(cube_scores, [25, 75])
            ],
        }
        if name != "pre":
            improvement = production - cube_scores
            item["vs_production_median_improvement"] = float(
                np.median(improvement)
            )
            item["vs_production_bootstrap95"] = bootstrap_median_ci(improvement)
            displacement = np.linalg.norm(
                arms[name]["verts"] - pre["verts"], axis=1
            )
            item["displacement_vs_pre"] = percentiles(displacement[used])
        summaries[name] = item

    result = {
        "tool": "score_snap_ridge",
        "protocol": {
            "split": args.split,
            "dev_z_origins": [4352, 4480, 4608],
            "test_z_origin": 4736,
            "sigma": SIGMA,
            "half": HALF,
            "step": STEP,
            "outer_margin": MARGIN,
            "normal_source": "pre-snap OBJ, fixed for all arms",
            "inferential_unit": "cube",
            "bootstrap_seed": BOOTSTRAP_SEED,
            "bootstrap_resamples": BOOTSTRAP_N,
            "production": args.production,
            "candidate": args.candidate,
        },
        "cache_metadata": cache_metadata,
        "inputs": {
            "pre": {"path": pre["path"], "sha256": pre["sha256"]},
            "arms": {
                name: {"path": arm["path"], "sha256": arm["sha256"]}
                for name, arm in arms.items()
            },
        },
        "invariants": {
            "vertices": len(pre["verts"]),
            "faces": len(pre["faces"]),
            "cubes_total": len(pre["cubes"]),
            "cubes_scored": len(per_cube),
            "eligible_before_bracketing": len(indices),
            "common_bracketed": int(common.sum()),
            "topology_uv_cube_ranges_identical": True,
        },
        "summary": summaries,
        "per_cube": per_cube,
    }
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--raw-dir", required=True, type=Path)
    prepare.add_argument(
        "--cache",
        type=Path,
        default=Path("data/PHerc0139-4x5x5/ct_sigma1.npy"),
    )
    run = commands.add_parser("score")
    run.add_argument("--split", choices=("dev", "test"), required=True)
    run.add_argument("--pre", required=True, type=Path)
    run.add_argument("--arm", action="append", required=True, type=parse_arm)
    run.add_argument("--production", required=True)
    run.add_argument("--candidate")
    run.add_argument("--cache", default="data/PHerc0139-4x5x5/ct_sigma1.npy")
    run.add_argument("--out", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        result = prepare_cache(args.raw_dir, args.cache)
        print(json.dumps(result, indent=2))
    else:
        result = score(args)
        print(
            json.dumps(
                {
                    "split": args.split,
                    "invariants": result["invariants"],
                    "summary": result["summary"],
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
