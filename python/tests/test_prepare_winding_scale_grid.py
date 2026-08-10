import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
import tifffile


def _load_script():
    path = Path(__file__).resolve().parents[1] / "scripts" / "prepare_winding_scale_grid.py"
    spec = importlib.util.spec_from_file_location("prepare_winding_scale_grid_tested", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


BOX = (0, 128, 0, 128, 0, 128)
PRED = "s3://example/pred.zarr"
RAW = "s3://example/raw.zarr"
ID = "z00000_y00000_x00000"


def _halo(tmp_path):
    halo = tmp_path / "halo"
    pred_dir = halo / "cubes_PRED"
    raw_dir = halo / "cubes_RAW"
    pred_dir.mkdir(parents=True)
    raw_dir.mkdir()
    pred = np.zeros((128, 128, 128), dtype=np.uint8)
    raw = np.zeros_like(pred)
    pred[1, 2, 3] = 255
    raw[4, 5, 6] = 17
    tifffile.imwrite(pred_dir / f"{ID}.tif", pred, compression=None, rowsperstrip=128)
    tifffile.imwrite(raw_dir / f"{ID}.tif", raw, compression=None, rowsperstrip=128)
    (pred_dir / "present.json").write_text(json.dumps([ID], indent=0), encoding="utf-8")
    manifest = {
        "chunk_size": 128,
        "bbox_l0_zyx": list(BOX),
        "n_chunks": [1, 1, 1],
        "sources": {"pred": PRED, "raw": RAW},
        "umbilicus_yx": [10.0, 20.0],
        "created_by": "carve_grid_tifs.py",
    }
    (halo / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return halo


def test_validates_and_hardlinks_exact_core(tmp_path):
    module = _load_script()
    halo = _halo(tmp_path)
    core = tmp_path / "core"
    result = module.prepare(
        halo,
        core,
        halo_bbox=BOX,
        core_bbox=BOX,
        required_pred_bbox=BOX,
        pred_zarr=PRED,
        raw_zarr=RAW,
        umbilicus_yx=(10.0, 20.0),
        min_core_pred=1,
    )
    assert result["status"] == "PASS"
    assert result["raw_cubes"] == result["pred_cubes_in_core"] == 1
    assert (core / "cubes_PRED" / f"{ID}.tif").read_bytes() == (
        halo / "cubes_PRED" / f"{ID}.tif"
    ).read_bytes()
    stored = json.loads((core / "input_validation.json").read_text(encoding="utf-8"))
    assert stored["file_rows_canonical_sha256"] == result["file_rows_canonical_sha256"]


def test_rejects_present_inventory_mismatch(tmp_path):
    module = _load_script()
    halo = _halo(tmp_path)
    (halo / "cubes_PRED" / "present.json").write_text("[]", encoding="utf-8")
    with pytest.raises(ValueError, match="present.json"):
        module.prepare(
            halo,
            tmp_path / "core",
            halo_bbox=BOX,
            core_bbox=BOX,
            required_pred_bbox=BOX,
            pred_zarr=PRED,
            raw_zarr=RAW,
            umbilicus_yx=(10.0, 20.0),
            min_core_pred=1,
        )
