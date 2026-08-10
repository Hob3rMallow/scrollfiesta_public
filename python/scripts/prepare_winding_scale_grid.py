#!/usr/bin/env python3
"""Validate a carved halo grid and hard-link its frozen meshing core.

``carve_grid_tifs.py`` writes RAW and prediction cubes for one box.  A scale
stress test needs the complete outer RAW halo downstream, but must mesh only a
smaller registered core.  This helper validates every carved TIFF, records
hashes, and creates a core ``cubes_PRED`` inventory with hard links so no
prediction cube is silently recomputed or copied differently.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import tifffile


CUBE = 128
CUBE_RE = re.compile(r"^z\d+_y\d+_x\d+$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_bbox(values: tuple[int, ...], label: str) -> tuple[int, ...]:
    if len(values) != 6:
        raise ValueError(f"{label}: expected six values")
    z0, z1, y0, y1, x0, x1 = values
    if not (z0 < z1 and y0 < y1 and x0 < x1):
        raise ValueError(f"{label}: invalid bounds")
    if any(value % CUBE for value in values):
        raise ValueError(f"{label}: bounds must be multiples of {CUBE}")
    return values


def contains(outer: tuple[int, ...], inner: tuple[int, ...]) -> bool:
    return (
        outer[0] <= inner[0] < inner[1] <= outer[1]
        and outer[2] <= inner[2] < inner[3] <= outer[3]
        and outer[4] <= inner[4] < inner[5] <= outer[5]
    )


def lattice_ids(bbox: tuple[int, ...]) -> set[str]:
    z0, z1, y0, y1, x0, x1 = bbox
    return {
        f"z{z:05d}_y{y:05d}_x{x:05d}"
        for z in range(z0, z1, CUBE)
        for y in range(y0, y1, CUBE)
        for x in range(x0, x1, CUBE)
    }


def _tiff_row(path: Path, kind: str) -> dict[str, Any]:
    cube_id = path.stem
    if CUBE_RE.fullmatch(cube_id) is None:
        raise ValueError(f"invalid cube filename: {path.name}")
    with tifffile.TiffFile(path) as tif:
        array = tif.asarray()
    if array.shape != (CUBE, CUBE, CUBE):
        raise ValueError(f"{path}: shape {array.shape}, expected {(CUBE,) * 3}")
    if array.dtype != np.uint8:
        raise ValueError(f"{path}: dtype {array.dtype}, expected uint8")
    nonzero = int(np.count_nonzero(array))
    if kind == "PRED":
        unique = np.unique(array)
        if nonzero == 0:
            raise ValueError(f"{path}: carved PRED TIFF is all zero")
        if not set(map(int, unique)).issubset({0, 255}):
            raise ValueError(f"{path}: PRED values are not binary 0/255")
    return {
        "id": cube_id,
        "kind": kind,
        "bytes": path.stat().st_size,
        "nonzero_voxels": nonzero,
        "sha256": sha256(path),
    }


def prepare(
    halo_grid: Path,
    core_grid: Path,
    *,
    halo_bbox: tuple[int, ...],
    core_bbox: tuple[int, ...],
    required_pred_bbox: tuple[int, ...],
    pred_zarr: str,
    raw_zarr: str,
    umbilicus_yx: tuple[float, float],
    min_core_pred: int = 100,
) -> dict[str, Any]:
    halo_bbox = validate_bbox(halo_bbox, "halo_bbox")
    core_bbox = validate_bbox(core_bbox, "core_bbox")
    required_pred_bbox = validate_bbox(required_pred_bbox, "required_pred_bbox")
    if not contains(halo_bbox, core_bbox):
        raise ValueError("core_bbox must be contained in halo_bbox")
    if not contains(core_bbox, required_pred_bbox):
        raise ValueError("required_pred_bbox must be contained in core_bbox")
    if core_grid.exists():
        raise ValueError(f"core output already exists: {core_grid}")

    manifest_path = halo_grid / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_chunks = [
        (halo_bbox[1] - halo_bbox[0]) // CUBE,
        (halo_bbox[3] - halo_bbox[2]) // CUBE,
        (halo_bbox[5] - halo_bbox[4]) // CUBE,
    ]
    expected_manifest = {
        "chunk_size": CUBE,
        "bbox_l0_zyx": list(halo_bbox),
        "n_chunks": expected_chunks,
        "sources": {"pred": pred_zarr, "raw": raw_zarr},
        "umbilicus_yx": list(umbilicus_yx),
        "created_by": "carve_grid_tifs.py",
    }
    mismatches = [key for key, value in expected_manifest.items() if manifest.get(key) != value]
    if mismatches:
        raise ValueError(f"halo manifest differs in fields: {mismatches}")

    expected_halo = lattice_ids(halo_bbox)
    expected_core = lattice_ids(core_bbox)
    required_pred = lattice_ids(required_pred_bbox)
    raw_paths = sorted((halo_grid / "cubes_RAW").glob("*.tif"))
    pred_paths = sorted((halo_grid / "cubes_PRED").glob("*.tif"))
    raw_ids = {path.stem for path in raw_paths}
    pred_ids = {path.stem for path in pred_paths}
    if raw_ids != expected_halo:
        raise ValueError(
            f"RAW inventory mismatch: missing={sorted(expected_halo - raw_ids)}, "
            f"extra={sorted(raw_ids - expected_halo)}"
        )
    if not pred_ids <= expected_halo:
        raise ValueError(f"PRED inventory has out-of-box ids: {sorted(pred_ids - expected_halo)}")
    present_path = halo_grid / "cubes_PRED" / "present.json"
    present = json.loads(present_path.read_text(encoding="utf-8"))
    if present != sorted(pred_ids):
        raise ValueError("PRED present.json does not exactly match the TIFF inventory")

    core_pred_ids = pred_ids & expected_core
    if len(core_pred_ids) < min_core_pred:
        raise ValueError(
            f"only {len(core_pred_ids)} nonempty core PRED cubes; required {min_core_pred}"
        )
    if not required_pred <= core_pred_ids:
        raise ValueError(f"required prediction box is incomplete: {sorted(required_pred - core_pred_ids)}")

    rows = [_tiff_row(path, "PRED") for path in pred_paths]
    rows.extend(_tiff_row(path, "RAW") for path in raw_paths)
    if not any(row["kind"] == "RAW" and row["nonzero_voxels"] > 0 for row in rows):
        raise ValueError("every RAW cube is all zero")

    core_pred_dir = core_grid / "cubes_PRED"
    core_pred_dir.mkdir(parents=True)
    for cube_id in sorted(core_pred_ids):
        os.link(
            halo_grid / "cubes_PRED" / f"{cube_id}.tif",
            core_pred_dir / f"{cube_id}.tif",
        )
    (core_pred_dir / "present.json").write_text(
        json.dumps(sorted(core_pred_ids), indent=0), encoding="utf-8"
    )
    core_manifest = {
        "chunk_size": CUBE,
        "bbox_l0_zyx": list(core_bbox),
        "n_chunks": [
            (core_bbox[1] - core_bbox[0]) // CUBE,
            (core_bbox[3] - core_bbox[2]) // CUBE,
            (core_bbox[5] - core_bbox[4]) // CUBE,
        ],
        "sources": {"pred": pred_zarr, "raw": raw_zarr},
        "umbilicus_yx": list(umbilicus_yx),
        "created_by": "prepare_winding_scale_grid.py (hard links from halo grid)",
        "halo_grid": str(halo_grid),
        "halo_bbox_l0_zyx": list(halo_bbox),
    }
    (core_grid / "manifest.json").write_text(
        json.dumps(core_manifest, indent=1) + "\n", encoding="utf-8"
    )

    canonical = "".join(
        f"{row['kind']}\t{row['id']}\t{row['bytes']}\t{row['nonzero_voxels']}\t{row['sha256']}\n"
        for row in sorted(rows, key=lambda row: (row["kind"], row["id"]))
    )
    result = {
        "format": "scrollfiesta-winding-scale-input-v1",
        "status": "PASS",
        "halo_bbox_l0_zyx": list(halo_bbox),
        "core_bbox_l0_zyx": list(core_bbox),
        "required_pred_bbox_l0_zyx": list(required_pred_bbox),
        "halo_lattice_cubes": len(expected_halo),
        "raw_cubes": len(raw_ids),
        "pred_cubes_in_halo": len(pred_ids),
        "pred_cubes_in_core": len(core_pred_ids),
        "required_pred_cubes": len(required_pred),
        "raw_all_zero_cubes": sum(
            row["kind"] == "RAW" and row["nonzero_voxels"] == 0 for row in rows
        ),
        "manifest_sha256": sha256(manifest_path),
        "file_rows_canonical_sha256": hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
        "files": sorted(rows, key=lambda row: (row["kind"], row["id"])),
    }
    (core_grid / "input_validation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--halo-grid", type=Path, required=True)
    parser.add_argument("--core-grid", type=Path, required=True)
    parser.add_argument("--halo-bbox", nargs=6, type=int, required=True)
    parser.add_argument("--core-bbox", nargs=6, type=int, required=True)
    parser.add_argument("--required-pred-bbox", nargs=6, type=int, required=True)
    parser.add_argument("--pred-zarr", required=True)
    parser.add_argument("--raw-zarr", required=True)
    parser.add_argument("--umbilicus", nargs=2, type=float, required=True)
    parser.add_argument("--min-core-pred", type=int, default=100)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = prepare(
        args.halo_grid,
        args.core_grid,
        halo_bbox=tuple(args.halo_bbox),
        core_bbox=tuple(args.core_bbox),
        required_pred_bbox=tuple(args.required_pred_bbox),
        pred_zarr=args.pred_zarr,
        raw_zarr=args.raw_zarr,
        umbilicus_yx=tuple(args.umbilicus),
        min_core_pred=args.min_core_pred,
    )
    print(
        json.dumps(
            {
                "status": result["status"],
                "raw_cubes": result["raw_cubes"],
                "pred_cubes_in_core": result["pred_cubes_in_core"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
