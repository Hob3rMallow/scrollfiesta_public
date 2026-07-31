#!/usr/bin/env python3
"""Bake local chunked CT intensity through a TIFXYZ coordinate map.

The default sampler matches the native manual-unroll tools: trilinear CT
sampling at five points over +/-2 voxels along a normal derived from the
TIFXYZ position field, keeping the brightest sample. That small normal search
is important for papyrus surfaces because a fitted surface can sit just to one
side of the bright sheet core.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image


_CUBE_RE = re.compile(r"z(-?\d+)_y(-?\d+)_x(-?\d+)\.tif$", re.IGNORECASE)


def _load_map(directory: Path) -> tuple[list[np.ndarray], np.ndarray]:
    arrays = [
        np.asarray(tifffile.imread(directory / f"{axis}.tif"), dtype=np.float32)
        for axis in "xyz"
    ]
    if arrays[0].shape != arrays[1].shape or arrays[0].shape != arrays[2].shape:
        raise ValueError(f"coordinate shapes differ in {directory}")
    valid = np.logical_and.reduce(
        [
            np.isfinite(arrays[0]) & (arrays[0] > 0.0),
            np.isfinite(arrays[1]),
            np.isfinite(arrays[2]) & (arrays[2] > 0.0),
        ]
    )
    return arrays, valid


def _discover_cubes(directory: Path) -> dict[tuple[int, int, int], Path]:
    cubes = {}
    for path in sorted(directory.glob("*.tif")):
        match = _CUBE_RE.fullmatch(path.name)
        if match:
            key = tuple(int(value) for value in match.groups())
            cubes[key] = path
    if not cubes:
        raise ValueError(f"no zNNN_yNNN_xNNN.tif cubes found in {directory}")
    return cubes


def _position_normals(
    coordinates: list[np.ndarray], valid: np.ndarray
) -> tuple[np.ndarray, int]:
    """Estimate native-compatible normals from central row/column differences.

    Coordinates are converted from TIFXYZ (x,y,z) to the sampler's (z,y,x).
    As in ``tifxyz_render``, a normal is emitted only when all four central
    neighbours are valid. Degenerate/border pixels fall back to point sampling.
    """

    x, y, z = coordinates
    positions = np.stack((z, y, x), axis=-1)
    normals = np.zeros(positions.shape, dtype=np.float32)
    if positions.shape[0] < 3 or positions.shape[1] < 3:
        return normals, 0

    neighbours_valid = (
        valid[1:-1, 1:-1]
        & valid[1:-1, :-2]
        & valid[1:-1, 2:]
        & valid[:-2, 1:-1]
        & valid[2:, 1:-1]
    )
    dc = positions[1:-1, 2:] - positions[1:-1, :-2]
    dr = positions[2:, 1:-1] - positions[:-2, 1:-1]
    cross = np.cross(dc, dr)
    lengths = np.linalg.norm(cross, axis=-1)
    good = neighbours_valid & np.isfinite(lengths) & (lengths > 1.0e-6)
    interior = normals[1:-1, 1:-1]
    interior[good] = cross[good] / lengths[good, None]
    return normals, int(good.sum())


def _cube_origins_for_bbox(
    cubes: dict[tuple[int, int, int], Path],
    positions: np.ndarray,
    normal_range: float,
    chunk: int = 128,
) -> list[tuple[int, int, int]]:
    """Return discovered cubes intersecting the padded position bbox."""

    pad = max(normal_range, 0.0) + 1.0
    lower = np.floor(positions.min(axis=0) - pad)
    upper = np.ceil(positions.max(axis=0) + pad)
    selected = []
    for key in cubes:
        origin = np.asarray(key, dtype=np.float64)
        if np.all(origin <= upper) and np.all(origin + chunk > lower):
            selected.append(key)
    if not selected:
        raise ValueError("no RAW cubes overlap the coordinate-map bounding box")
    return sorted(selected)


def _load_volume(
    cubes: dict[tuple[int, int, int], Path],
    selected: list[tuple[int, int, int]],
    chunk: int = 128,
) -> tuple[np.ndarray, np.ndarray]:
    """Load a contiguous selected cube brick and return it plus (z,y,x) origin."""

    origins = np.asarray(selected, dtype=np.int64)
    lower = origins.min(axis=0)
    upper = origins.max(axis=0)
    axis_origins = [
        list(range(int(lower[axis]), int(upper[axis]) + chunk, chunk))
        for axis in range(3)
    ]
    required = {
        (z, y, x)
        for z in axis_origins[0]
        for y in axis_origins[1]
        for x in axis_origins[2]
    }
    missing = sorted(required.difference(cubes))
    if missing:
        raise ValueError(
            "RAW cubes overlapping the surface do not form a contiguous brick; "
            f"missing {len(missing)} cube(s), first={missing[0]}"
        )

    shape = tuple(int(upper[axis] - lower[axis] + chunk) for axis in range(3))
    volume = np.empty(shape, dtype=np.uint8)
    for key in sorted(required):
        cube = np.asarray(tifffile.imread(cubes[key]))
        if cube.shape != (chunk, chunk, chunk):
            raise ValueError(
                f"{cubes[key]} has shape {cube.shape}, expected "
                f"{(chunk, chunk, chunk)}"
            )
        offsets = tuple(int(key[axis] - lower[axis]) for axis in range(3))
        slices = tuple(
            slice(offsets[axis], offsets[axis] + chunk) for axis in range(3)
        )
        volume[slices] = cube
    return volume, lower


def _sample_trilinear(
    volume: np.ndarray, local_positions: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Sample like raw_sample.c, including missing-corner renormalization."""

    base = np.floor(local_positions).astype(np.int64)
    fraction = local_positions - base
    values = np.zeros(len(local_positions), dtype=np.float32)
    weights = np.zeros(len(local_positions), dtype=np.float32)
    shape = np.asarray(volume.shape, dtype=np.int64)
    for oz in (0, 1):
        wz = fraction[:, 0] if oz else 1.0 - fraction[:, 0]
        for oy in (0, 1):
            wy = fraction[:, 1] if oy else 1.0 - fraction[:, 1]
            for ox in (0, 1):
                wx = fraction[:, 2] if ox else 1.0 - fraction[:, 2]
                index = base + np.asarray((oz, oy, ox), dtype=np.int64)
                inside = np.logical_and.reduce(
                    [
                        index[:, 0] >= 0,
                        index[:, 1] >= 0,
                        index[:, 2] >= 0,
                        index[:, 0] < shape[0],
                        index[:, 1] < shape[1],
                        index[:, 2] < shape[2],
                    ]
                )
                if not inside.any():
                    continue
                weight = (wz * wy * wx)[inside].astype(np.float32, copy=False)
                ii = index[inside]
                values[inside] += (
                    volume[ii[:, 0], ii[:, 1], ii[:, 2]].astype(np.float32)
                    * weight
                )
                weights[inside] += weight
    sampled = weights > 1.0e-9
    values[sampled] /= weights[sampled]
    values[~sampled] = -1.0
    return values, sampled


def _sample_normal_max(
    volume: np.ndarray,
    volume_origin: np.ndarray,
    positions: np.ndarray,
    normals: np.ndarray,
    normal_range: float,
    normal_samples: int,
    sample_chunk: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Trilinearly point-sample or max-sample symmetrically along normals."""

    count = len(positions)
    output = np.full(count, -1.0, dtype=np.float32)
    sampled = np.zeros(count, dtype=bool)
    if normal_range <= 0.0 or normal_samples < 2:
        taps = np.asarray((0.0,), dtype=np.float32)
    else:
        taps = np.linspace(
            -normal_range, normal_range, normal_samples, dtype=np.float32
        )

    for start in range(0, count, sample_chunk):
        stop = min(start + sample_chunk, count)
        points = positions[start:stop]
        directions = normals[start:stop]
        best = np.full(stop - start, -1.0, dtype=np.float32)
        any_sampled = np.zeros(stop - start, dtype=bool)
        for tap in taps:
            local = points + tap * directions - volume_origin
            values, tap_sampled = _sample_trilinear(volume, local)
            improve = tap_sampled & ((~any_sampled) | (values > best))
            best[improve] = values[improve]
            any_sampled |= tap_sampled
        output[start:stop] = best
        sampled[start:stop] = any_sampled
    return output, sampled


def _contrast_preview(
    texture: np.ndarray, valid: np.ndarray
) -> tuple[np.ndarray, dict]:
    values = texture[valid]
    if not len(values):
        raise ValueError("the coordinate map did not sample any source voxels")
    low, high = np.percentile(values, [1.0, 99.0])
    if high <= low:
        low = float(values.min())
        high = float(values.max())
    scale = 255.0 / max(float(high - low), 1.0)
    preview = np.zeros(texture.shape, dtype=np.uint8)
    preview[valid] = np.clip((values - low) * scale, 0.0, 255.0).astype(np.uint8)
    return preview, {
        "minimum": float(values.min()),
        "maximum": float(values.max()),
        "mean": float(values.mean()),
        "contrast_p01": float(low),
        "contrast_p99": float(high),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface", required=True, type=Path)
    parser.add_argument("--raw-cubes", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--preview-scale", type=int, default=8)
    parser.add_argument(
        "--normal-range",
        type=float,
        default=2.0,
        help="sample +/- this many voxels along the position-field normal (default: 2)",
    )
    parser.add_argument(
        "--normal-samples",
        type=int,
        default=5,
        help="number of evenly spaced normal samples (default: 5)",
    )
    parser.add_argument(
        "--sample-chunk",
        type=int,
        default=200000,
        help="coordinate points processed per vectorized sampling chunk",
    )
    args = parser.parse_args()
    if args.preview_scale <= 0:
        parser.error("--preview-scale must be positive")
    if args.normal_range < 0:
        parser.error("--normal-range cannot be negative")
    if args.normal_samples <= 0:
        parser.error("--normal-samples must be positive")
    if args.normal_range > 0 and args.normal_samples < 2:
        parser.error("--normal-samples must be at least 2 when --normal-range > 0")
    if args.sample_chunk <= 0:
        parser.error("--sample-chunk must be positive")

    coordinates, map_valid = _load_map(args.surface)
    cube_paths = _discover_cubes(args.raw_cubes)
    x, y, z = coordinates
    normals, normal_points = _position_normals(coordinates, map_valid)
    positions = np.stack((z[map_valid], y[map_valid], x[map_valid]), axis=-1)
    selected = _cube_origins_for_bbox(cube_paths, positions, args.normal_range)
    volume, volume_origin = _load_volume(cube_paths, selected)
    values, valid_samples = _sample_normal_max(
        volume,
        volume_origin.astype(np.float32),
        positions,
        normals[map_valid],
        args.normal_range,
        args.normal_samples,
        args.sample_chunk,
    )

    texture_float = np.zeros(x.shape, dtype=np.float32)
    sampled = np.zeros(x.shape, dtype=bool)
    texture_float[map_valid] = np.where(valid_samples, values, 0.0)
    sampled[map_valid] = valid_samples
    texture = np.clip(np.rint(texture_float), 0.0, 255.0).astype(np.uint8)
    preview, intensity_stats = _contrast_preview(texture_float, sampled)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    raw_path = args.output_prefix.with_name(args.output_prefix.name + "_raw.tif")
    preview_path = args.output_prefix.with_name(args.output_prefix.name + "_preview.png")
    mask_path = args.output_prefix.with_name(args.output_prefix.name + "_mask.png")
    report_path = args.output_prefix.with_name(args.output_prefix.name + "_bake.json")
    tifffile.imwrite(raw_path, texture)
    target_size = (
        preview.shape[1] * args.preview_scale,
        preview.shape[0] * args.preview_scale,
    )
    Image.fromarray(preview).resize(
        target_size, Image.Resampling.NEAREST
    ).save(preview_path)
    Image.fromarray((sampled * 255).astype(np.uint8)).resize(
        target_size, Image.Resampling.NEAREST
    ).save(mask_path)

    report = {
        "schema": "tifxyz-ct-bake-v2",
        "surface": str(args.surface.resolve()),
        "raw_cubes": str(args.raw_cubes.resolve()),
        "grid_shape": list(texture.shape),
        "coordinate_valid_points": int(map_valid.sum()),
        "sampled_points": int(sampled.sum()),
        "sampled_fraction_of_coordinates": float(sampled.sum() / map_valid.sum()),
        "sampling": {
            "method": "trilinear-normal-max"
            if args.normal_range > 0
            else "trilinear-point",
            "normal_range_voxels": args.normal_range,
            "normal_samples": args.normal_samples if args.normal_range > 0 else 1,
            "normal_source": "central differences of TIFXYZ position field",
            "points_with_usable_normals": normal_points,
            "points_falling_back_to_point_sample": int(map_valid.sum() - normal_points),
        },
        "loaded_cubes": len(selected),
        "loaded_cube_origins_zyx": [list(key) for key in selected],
        "loaded_volume_origin_zyx": volume_origin.tolist(),
        "loaded_volume_shape_zyx": list(volume.shape),
        "missing_cubes": [],
        "intensity": intensity_stats,
        "outputs": {
            "raw_tif": raw_path.name,
            "preview_png": preview_path.name,
            "mask_png": mask_path.name,
        },
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if bool(sampled[map_valid].all()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
