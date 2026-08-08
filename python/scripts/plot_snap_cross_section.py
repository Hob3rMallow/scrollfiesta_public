#!/usr/bin/env python3
"""Render an identical PHerc0139 CT cross-section with two snap arms.

The OBJ files produced by ``scroll_unroll --export-mesh`` store source-space
coordinates as z, y, x.  This script intersects their shared triangle topology
with one fixed z plane and overlays the resulting y/x contours on the matching
RAW CT slice.  It is deliberately visualization-only: it does not calculate or
claim a sheet-jump metric.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from matplotlib.collections import LineCollection


CUBE_RE = re.compile(r"^z(?P<z>\d+)_y(?P<y>\d+)_x(?P<x>\d+)\.tif$")


@dataclass(frozen=True)
class ObjMesh:
    vertices: np.ndarray
    faces: np.ndarray


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_obj(path: Path) -> ObjMesh:
    vertices: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    with path.open("r", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            if line.startswith("v "):
                fields = line.split()
                if len(fields) != 4:
                    raise ValueError(f"{path}:{line_number}: malformed vertex")
                vertices.append(tuple(map(float, fields[1:4])))
            elif line.startswith("f "):
                fields = line.split()
                if len(fields) != 4:
                    raise ValueError(f"{path}:{line_number}: non-triangle face")
                face = tuple(int(field.split("/", 1)[0]) - 1 for field in fields[1:4])
                faces.append(face)

    vertex_array = np.asarray(vertices, dtype=np.float64)
    face_array = np.asarray(faces, dtype=np.int64)
    if vertex_array.ndim != 2 or vertex_array.shape[1] != 3:
        raise ValueError(f"{path}: no valid vertices")
    if face_array.ndim != 2 or face_array.shape[1] != 3:
        raise ValueError(f"{path}: no valid triangle faces")
    if face_array.min() < 0 or face_array.max() >= len(vertex_array):
        raise ValueError(f"{path}: face index outside vertex array")
    return ObjMesh(vertex_array, face_array)


def load_raw_slice(raw_dir: Path, z_plane: float) -> tuple[np.ndarray, tuple[int, int, int, int], int]:
    records: list[tuple[int, int, int, Path]] = []
    for path in sorted(raw_dir.glob("*.tif")):
        match = CUBE_RE.fullmatch(path.name)
        if match is None:
            continue
        z0, y0, x0 = (int(match.group(axis)) for axis in ("z", "y", "x"))
        if z0 <= z_plane < z0 + 128:
            records.append((z0, y0, x0, path))
    if not records:
        raise ValueError(f"no RAW cube in {raw_dir} contains z={z_plane}")

    z_origins = {record[0] for record in records}
    if len(z_origins) != 1:
        raise ValueError(f"z={z_plane} unexpectedly maps to multiple z origins: {z_origins}")
    z0 = records[0][0]
    y_min = min(record[1] for record in records)
    x_min = min(record[2] for record in records)
    y_max = max(record[1] + 128 for record in records)
    x_max = max(record[2] + 128 for record in records)
    expected = ((y_max - y_min) // 128) * ((x_max - x_min) // 128)
    if len(records) != expected:
        raise ValueError(f"incomplete RAW tile grid: found {len(records)}, expected {expected}")

    mosaic = np.zeros((y_max - y_min, x_max - x_min), dtype=np.uint8)
    occupancy = np.zeros_like(mosaic, dtype=bool)
    local_z = int(round(z_plane)) - z0
    if abs(z_plane - round(z_plane)) > 1e-9:
        raise ValueError("the CT background requires an integer voxel-centre z plane")
    for _, y0, x0, path in records:
        cube = tifffile.imread(path)
        if cube.shape != (128, 128, 128) or cube.dtype != np.uint8:
            raise ValueError(f"{path}: expected uint8 (128,128,128), got {cube.dtype} {cube.shape}")
        ys = slice(y0 - y_min, y0 - y_min + 128)
        xs = slice(x0 - x_min, x0 - x_min + 128)
        if occupancy[ys, xs].any():
            raise ValueError(f"{path}: overlapping tile")
        mosaic[ys, xs] = cube[local_z]
        occupancy[ys, xs] = True
    if not occupancy.all():
        raise ValueError("RAW mosaic has uncovered pixels")
    return mosaic, (y_min, y_max, x_min, x_max), len(records)


def plane_segments(mesh: ObjMesh, z_plane: float) -> np.ndarray:
    """Return triangle/plane intersections as [segment, endpoint, (x,y)]."""
    triangles = mesh.vertices[mesh.faces]
    z_values = triangles[:, :, 0]
    candidates = triangles[(z_values.min(axis=1) <= z_plane) & (z_values.max(axis=1) >= z_plane)]
    segments: list[np.ndarray] = []
    eps = 1e-9

    for triangle in candidates:
        delta = triangle[:, 0] - z_plane
        if np.all(np.abs(delta) <= eps):
            continue
        points: list[np.ndarray] = []
        for first, second in ((0, 1), (1, 2), (2, 0)):
            d0, d1 = delta[first], delta[second]
            if abs(d0) <= eps:
                points.append(triangle[first])
            if d0 * d1 < 0.0:
                fraction = -d0 / (d1 - d0)
                points.append(triangle[first] + fraction * (triangle[second] - triangle[first]))

        unique: list[np.ndarray] = []
        for point in points:
            if not any(np.linalg.norm(point - prior) <= 1e-7 for prior in unique):
                unique.append(point)
        if len(unique) < 2:
            continue
        if len(unique) > 2:
            pairs = [
                (np.linalg.norm(unique[a] - unique[b]), a, b)
                for a in range(len(unique))
                for b in range(a + 1, len(unique))
            ]
            _, first, second = max(pairs)
        else:
            first, second = 0, 1
        # OBJ is z,y,x; LineCollection expects x,y.
        segment = np.asarray(
            [[unique[first][2], unique[first][1]], [unique[second][2], unique[second][1]]],
            dtype=np.float64,
        )
        if np.linalg.norm(segment[1] - segment[0]) > 1e-7:
            segments.append(segment)

    if not segments:
        raise ValueError(f"mesh has no non-degenerate intersection at z={z_plane}")
    return np.stack(segments)


def add_scale_bar(axis: plt.Axes, bounds: tuple[int, int, int, int], voxel_um: float) -> None:
    y_min, y_max, x_min, x_max = bounds
    length_voxels = 1000.0 / voxel_um
    x0 = x_min + 0.055 * (x_max - x_min)
    y0 = y_min + 0.06 * (y_max - y_min)
    axis.plot([x0, x0 + length_voxels], [y0, y0], color="white", linewidth=4, solid_capstyle="butt")
    axis.text(
        x0 + length_voxels / 2,
        y0 + 0.018 * (y_max - y_min),
        "1 mm",
        color="white",
        fontsize=9,
        ha="center",
        va="bottom",
        weight="bold",
    )


def render(args: argparse.Namespace) -> dict[str, object]:
    repair = read_obj(args.repair_obj)
    recto = read_obj(args.recto_obj)
    if repair.vertices.shape != recto.vertices.shape:
        raise ValueError("the two arms have different vertex-array shapes")
    if repair.faces.shape != recto.faces.shape or not np.array_equal(repair.faces, recto.faces):
        raise ValueError("the two arms do not have identical face topology")

    raw, bounds, tile_count = load_raw_slice(args.raw_dir, args.z)
    repair_segments = plane_segments(repair, args.z)
    recto_segments = plane_segments(recto, args.z)
    y_min, y_max, x_min, x_max = bounds
    vmin, vmax = np.percentile(raw, [1.0, 99.0])
    if not vmin < vmax:
        raise ValueError("RAW slice has no display contrast")

    figure, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    common = dict(
        origin="lower",
        cmap="gray",
        vmin=vmin,
        vmax=vmax,
        extent=(x_min, x_max, y_min, y_max),
        interpolation="nearest",
    )
    panels = (
        ("Repair only (recto iterations 0)", [(repair_segments, "#00e676", "solid", 0.9)]),
        ("Oriented boundary (recto iterations 4)", [(recto_segments, "#ff3d71", "solid", 0.9)]),
        (
            "Identical CT slice, both contours",
            [
                (repair_segments, "#00e676", "solid", 1.0),
                (recto_segments, "#ff3d71", "dashed", 0.9),
            ],
        ),
    )
    for axis, (title, overlays) in zip(axes, panels, strict=True):
        axis.imshow(raw, **common)
        for segments, color, style, alpha in overlays:
            axis.add_collection(
                LineCollection(segments, colors=color, linewidths=0.75, linestyles=style, alpha=alpha)
            )
        axis.set_xlim(x_min, x_max)
        axis.set_ylim(y_min, y_max)
        axis.set_aspect("equal")
        axis.set_title(title, fontsize=13)
        axis.set_xlabel("world x (voxels)")
        axis.set_ylabel("world y (voxels)")
        add_scale_bar(axis, bounds, args.voxel_um)

    axes[2].plot([], [], color="#00e676", linewidth=2, label="repair only")
    axes[2].plot([], [], color="#ff3d71", linewidth=2, linestyle="--", label="4 recto iterations")
    axes[2].legend(loc="upper right", framealpha=0.85, fontsize=9)
    figure.suptitle(
        f"PHerc0139 held-out midpoint cross-section z={args.z:.0f} | "
        f"RAW CT p1–p99 display, {args.voxel_um:g} µm voxels\n"
        "Contours are triangle/plane intersections; visualization is not a sheet-jump metric",
        fontsize=15,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out, dpi=180, facecolor="white")
    plt.close(figure)
    return {
        "output": str(args.out),
        "output_sha256": sha256(args.out),
        "z": args.z,
        "raw_tiles": tile_count,
        "raw_bounds_yx": [y_min, y_max, x_min, x_max],
        "raw_display_p1_p99": [float(vmin), float(vmax)],
        "vertices": int(repair.vertices.shape[0]),
        "faces": int(repair.faces.shape[0]),
        "repair_segments": int(len(repair_segments)),
        "recto_segments": int(len(recto_segments)),
        "repair_obj_sha256": sha256(args.repair_obj),
        "recto_obj_sha256": sha256(args.recto_obj),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--repair-obj", type=Path, required=True)
    parser.add_argument("--recto-obj", type=Path, required=True)
    parser.add_argument("--z", type=float, default=4800.0)
    parser.add_argument("--voxel-um", type=float, default=9.362)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    import json

    print(json.dumps(render(parse_args()), indent=2, sort_keys=True))
