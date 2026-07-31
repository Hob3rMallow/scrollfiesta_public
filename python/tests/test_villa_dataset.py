import json

import numpy as np
import pytest
import tifffile

from scrollunwrap.villa_dataset import (classify_patch, import_scroll_hint,
                                        init_dataset, materialize_mask,
                                        register_patch, validate_patch)


def _write_patch(path, *, masked_positive=False):
    path.mkdir(parents=True)
    y, x = np.mgrid[:4, :4]
    tifffile.imwrite(path / "x.tif", (x + 10).astype(np.float32))
    tifffile.imwrite(path / "y.tif", (y + 20).astype(np.float32))
    tifffile.imwrite(path / "z.tif", np.full((4, 4), 30, np.float32))
    mask = np.full((4, 4), 255, np.uint8)
    if masked_positive:
        mask[0, 0] = 0
    tifffile.imwrite(path / "mask.tif", mask)
    (path / "meta.json").write_text(json.dumps(
        {"format": "tifxyz", "uuid": path.name, "scale": [1.0, 1.0]}))


def test_dataset_contract_and_registration(tmp_path):
    root = init_dataset(tmp_path / "dataset")
    patch = root / "verified_patches" / "p0"
    _write_patch(patch)
    stats = validate_patch(patch)
    assert stats["valid_quads"] == 9
    register_patch(root, "p0", "verified", patch,
                   component_report={"faces": 2})
    manifest = json.loads((root / "scrollfiesta_manifest.json").read_text())
    assert manifest["format"] == "scrollfiesta-villa-dataset-v1"
    assert manifest["patches"][0]["path"] == "verified_patches/p0"
    for name in manifest["point_collections"]:
        pcl = json.loads((root / name).read_text())
        assert pcl["vc_pointcollections_json_version"] == "1"


def test_masked_cells_must_also_be_coordinate_invalid(tmp_path):
    patch = tmp_path / "p"
    _write_patch(patch, masked_positive=True)
    with pytest.raises(ValueError, match="masked cells retain"):
        validate_patch(patch)


def test_auto_classification_is_conservative():
    assert classify_patch({"uv": {"flipped_uv_triangles": 0,
                                  "degenerate_uv_triangles": 0}}) == "verified"
    assert classify_patch({"uv": {"flipped_uv_triangles": 1,
                                  "degenerate_uv_triangles": 0}}) == "unverified"


def test_materialize_mask_hides_invalid_coordinates(tmp_path):
    patch = tmp_path / "p"
    _write_patch(patch)
    (patch / "mask.tif").unlink()
    x = tifffile.imread(patch / "x.tif")
    x[0, 0] = -1
    tifffile.imwrite(patch / "x.tif", x)
    materialize_mask(patch)
    mask = tifffile.imread(patch / "mask.tif")
    assert mask[0, 0] == 0
    for axis in "xyz":
        assert tifffile.imread(patch / f"{axis}.tif")[0, 0] == -1


def test_scroll_hint_is_always_imported_as_unverified(tmp_path):
    source = tmp_path / "scrollfiesta_hint"
    _write_patch(source)
    root = tmp_path / "dataset"

    entry = import_scroll_hint(root, source, "hint0")

    assert entry["classification"] == "unverified"
    assert entry["path"] == "unverified_patches/hint0"
    assert entry["component"]["trust"] == "unverified"
    assert (root / entry["path"] / "x.tif").is_file()
    assert not (root / "verified_patches" / "hint0").exists()


def test_scroll_hint_import_refuses_overwrite(tmp_path):
    source = tmp_path / "scrollfiesta_hint"
    _write_patch(source)
    root = tmp_path / "dataset"
    import_scroll_hint(root, source, "hint0")

    with pytest.raises(ValueError, match="refusing to overwrite"):
        import_scroll_hint(root, source, "hint0")
