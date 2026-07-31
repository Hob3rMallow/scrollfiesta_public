"""Build and validate the file contract consumed by Villa's spiral fitter."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import tifffile


FORMAT = "scrollfiesta-villa-dataset-v1"
PCL_FILES = ("abs_winding.json", "patch-overlap-pcls.json",
             "relative_windings.json", "same_windings.json")


def _empty_point_collection() -> dict:
    return {"vc_pointcollections_json_version": "1", "collections": {}}


def init_dataset(root: Path) -> Path:
    root = Path(root)
    (root / "verified_patches").mkdir(parents=True, exist_ok=True)
    (root / "unverified_patches").mkdir(parents=True, exist_ok=True)
    for name in PCL_FILES:
        path = root / name
        if not path.exists():
            path.write_text(json.dumps(_empty_point_collection(), indent=2))
    manifest = root / "scrollfiesta_manifest.json"
    if not manifest.exists():
        manifest.write_text(json.dumps({
            "format": FORMAT,
            "coordinate_order": "xyz_level0_voxels",
            "patches": [],
            "point_collections": list(PCL_FILES),
        }, indent=2))
    config = root / "villa_dataset.json"
    config.write_text(json.dumps({
        "format": FORMAT,
        "verified_patches": "verified_patches",
        "unverified_patches": "unverified_patches",
        "point_collections": list(PCL_FILES),
        "manifest": "scrollfiesta_manifest.json",
    }, indent=2))
    return root


def validate_patch(path: Path, *, require_mask: bool = True) -> dict:
    path = Path(path)
    required = ["x.tif", "y.tif", "z.tif", "meta.json"]
    if require_mask:
        required.append("mask.tif")
    missing = [name for name in required if not (path / name).is_file()]
    if missing:
        raise ValueError(f"{path}: missing Villa tifxyz files: {missing}")

    xyz = [tifffile.imread(path / f"{axis}.tif") for axis in "xyz"]
    if any(a.ndim != 2 for a in xyz) or len({a.shape for a in xyz}) != 1:
        raise ValueError(f"{path}: x/y/z channels must be equal-size 2D rasters")
    if any(a.dtype.kind != "f" for a in xyz):
        raise ValueError(f"{path}: x/y/z channels must be floating point")
    finite = np.isfinite(xyz[0]) & np.isfinite(xyz[1]) & np.isfinite(xyz[2])
    valid = finite & (xyz[0] >= 0) & (xyz[1] >= 0) & (xyz[2] >= 0)
    if require_mask:
        mask = tifffile.imread(path / "mask.tif")
        if mask.shape != xyz[0].shape:
            raise ValueError(f"{path}: mask shape {mask.shape} != xyz {xyz[0].shape}")
        valid &= mask != 0
        # A contested/invalid cell must never remain visible to a consumer that
        # ignores mask.tif. This duplicates the safety represented by the mask.
        hidden_but_positive = ((mask == 0) & finite & (xyz[0] >= 0)
                               & (xyz[1] >= 0) & (xyz[2] >= 0))
        if np.any(hidden_but_positive):
            raise ValueError(
                f"{path}: {int(hidden_but_positive.sum())} masked cells retain "
                "positive coordinates")
    quad = (valid[:-1, :-1] & valid[1:, :-1]
            & valid[:-1, 1:] & valid[1:, 1:])
    if not np.any(quad):
        raise ValueError(f"{path}: no valid tifxyz quad for Villa")
    meta = json.loads((path / "meta.json").read_text())
    if meta.get("format") != "tifxyz":
        raise ValueError(f"{path}: meta.format is not tifxyz")
    scale = meta.get("scale")
    if not (isinstance(scale, list) and len(scale) == 2
            and all(float(v) > 0 for v in scale)):
        raise ValueError(f"{path}: invalid meta.scale {scale!r}")
    return {
        "shape": list(xyz[0].shape),
        "valid_cells": int(valid.sum()),
        "valid_quads": int(quad.sum()),
        "scale": [float(scale[0]), float(scale[1])],
        "uuid": meta.get("uuid"),
    }


def materialize_mask(path: Path) -> None:
    """Make the mask contract explicit after a Villa transform.

    Villa flatteners preserve invalid (-1) vertices but do not all emit a mask.
    This writes mask.tif and forces every invalid coordinate triplet to -1, so
    both mask-aware and legacy consumers see the same topology.
    """
    path = Path(path)
    xyz = [tifffile.imread(path / f"{axis}.tif").astype(np.float32)
           for axis in "xyz"]
    valid = (np.isfinite(xyz[0]) & np.isfinite(xyz[1]) & np.isfinite(xyz[2])
             & (xyz[0] >= 0) & (xyz[1] >= 0) & (xyz[2] >= 0))
    old_mask = path / "mask.tif"
    if old_mask.is_file():
        valid &= tifffile.imread(old_mask) != 0
    for axis, values in zip("xyz", xyz):
        values[~valid] = -1.0
        tifffile.imwrite(path / f"{axis}.tif", values)
    tifffile.imwrite(old_mask, np.where(valid, 255, 0).astype(np.uint8))


def register_patch(root: Path, patch_id: str, classification: str,
                   patch_dir: Path, *, component_report: dict) -> dict:
    root = init_dataset(root)
    if classification not in {"verified", "unverified"}:
        raise ValueError(f"bad patch classification {classification!r}")
    expected = root / f"{classification}_patches" / patch_id
    if Path(patch_dir).resolve() != expected.resolve():
        raise ValueError(f"patch path {patch_dir} is not dataset path {expected}")
    validation = validate_patch(expected, require_mask=True)
    manifest_path = root / "scrollfiesta_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    entry = {
        "id": patch_id,
        "classification": classification,
        "path": str(Path(f"{classification}_patches") / patch_id).replace("\\", "/"),
        "validation": validation,
        "component": component_report,
    }
    manifest["patches"] = [p for p in manifest.get("patches", [])
                           if p.get("id") != patch_id]
    manifest["patches"].append(entry)
    manifest["patches"].sort(key=lambda p: p["id"])
    tmp = manifest_path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2))
    tmp.replace(manifest_path)
    return entry


def import_scroll_hint(root: Path, source: Path, patch_id: str) -> dict:
    """Import a ScrollFiesta tifxyz output as a low-trust spiral-fit hint.

    Automated ScrollFiesta geometry is useful evidence for Villa's scroll
    diffeomorphism, but it is not a verified surface. Keep that trust boundary
    structural: this helper always writes to ``unverified_patches``.
    """
    root = init_dataset(root)
    source = Path(source).resolve()
    if not source.is_dir():
        raise ValueError(f"ScrollFiesta hint directory not found: {source}")
    validation = validate_patch(source, require_mask=True)
    destination = root / "unverified_patches" / patch_id
    if destination.exists():
        raise ValueError(f"refusing to overwrite existing hint: {destination}")
    shutil.copytree(source, destination)
    return register_patch(
        root, patch_id, "unverified", destination,
        component_report={
            "source_type": "scrollfiesta_tifxyz_hint",
            "source": str(source),
            "trust": "unverified",
            "source_validation": validation,
        },
    )


def classify_patch(uv_report: dict, mode: str = "auto") -> str:
    if mode in {"verified", "unverified"}:
        return mode
    uv = uv_report.get("uv", {})
    clean = (uv.get("flipped_uv_triangles", 1) == 0
             and uv.get("degenerate_uv_triangles", 1) == 0)
    return "verified" if clean else "unverified"
