import importlib.util
from pathlib import Path

import numpy as np
import pytest
import tifffile


def _load_script():
    path = Path(__file__).resolve().parents[1] / "scripts" / "score_snap_ridge.py"
    spec = importlib.util.spec_from_file_location("score_snap_ridge_tested", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_prepare_and_validate_explicit_grid(tmp_path):
    score = _load_script()
    raw = tmp_path / "raw"
    raw.mkdir()
    cube = np.zeros((128, 128, 128), dtype=np.uint8)
    cube[63:66] = np.uint8(200)
    tifffile.imwrite(
        raw / "z00384_y00512_x00640.tif",
        cube,
        photometric="minisblack",
        compression=None,
        rowsperstrip=128,
    )
    origin = np.array([384.0, 512.0, 640.0])
    shape = (128, 128, 128)
    cache = tmp_path / "sigma1.npy"
    metadata = score.prepare_cache(raw, cache, origin=origin, shape=shape)
    assert metadata["origin_zyx"] == [384, 512, 640]
    assert metadata["shape_zyx"] == [128, 128, 128]
    assert metadata["n_cubes"] == 1
    assert np.load(cache, mmap_mode="r").shape == shape
    validated = score._validated_cache_metadata(cache, origin, shape)
    assert validated["cache_sha256"] == metadata["cache_sha256"]


def test_prepare_rejects_unaligned_grid(tmp_path):
    score = _load_script()
    with pytest.raises(ValueError, match="origin must contain multiples"):
        score.prepare_cache(
            tmp_path,
            tmp_path / "cache.npy",
            origin=np.array([1.0, 0.0, 0.0]),
            shape=(128, 128, 128),
        )


def test_ridge_offsets_are_translation_invariant_with_explicit_origin():
    score = _load_script()
    smooth = np.zeros((32, 32, 32), dtype=np.float32)
    smooth[16] = np.float32(100.0)
    local = np.array([[16.0, 16.0, 16.0]], dtype=np.float32)
    normals = np.array([[1.0, 0.0, 0.0]], dtype=np.float32)
    origin = np.array([384.0, 512.0, 640.0], dtype=np.float32)

    baseline = score.ridge_offsets(
        smooth, local, normals, origin=np.zeros(3, dtype=np.float32)
    )
    translated = score.ridge_offsets(
        smooth, local + origin, normals, origin=origin
    )
    np.testing.assert_allclose(translated, baseline)
