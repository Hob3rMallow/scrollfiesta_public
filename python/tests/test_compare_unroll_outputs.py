import copy
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
import tifffile


def _load_script():
    path = Path(__file__).resolve().parents[1] / "scripts" / "compare_unroll_outputs.py"
    spec = importlib.util.spec_from_file_location("compare_unroll_outputs_tested", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _stats(*, seam1, seam5, dcol1, fill=0.40, multi=0.02):
    names = (
        "step1_ribbon",
        "step2_join",
        "step3_overlap",
        "step4_snap",
        "step5_relax",
    )
    stages = []
    for index, name in enumerate(names):
        ratio = seam1 if index == 0 else seam5 if index == 4 else (seam1 + seam5) / 2
        dcol = dcol1 if index == 0 else dcol1 * 0.9
        stages.append(
            {
                "name": name,
                "W": 12,
                "H": 8,
                "bands": 1,
                "window": [10.0, 240.0],
                "fill": fill,
                "filled_px": 38,
                "multi_px": 2,
                "skip_uv_faces": 0,
                "skip_3d_faces": 0,
                "skip_own_faces": 0,
                "multi_frac": multi,
                "dark_frac": 0.01,
                "seam_dcol_mean": dcol,
                "base_dcol_mean": dcol / ratio,
                "seam_ratio": ratio,
                "vseam_gap_fill": 0.0,
                "vseam_ratio": 0.0,
                "n_vseam_rows": 0,
                "synth_px": 0,
                "synth_frac": 0.0,
                "fill_total": fill,
                "seconds": 1.0 + index,
            }
        )
    return {
        "tool": "scroll_unroll",
        "placed_dir": "placed",
        "raw_dir": "raw",
        "n_cubes": 25,
        "nv": 100,
        "nf": 150,
        "n_seam_cols": 20,
        "stages": stages,
        "total_seconds": 9.0,
    }


def _write_artifacts(module, directory, run_id, stats, *, delta=0):
    directory.mkdir(parents=True)
    stages = module._stage_map(stats, run_id)
    for stage_index, (stage_name, stage) in enumerate(stages.items()):
        image = np.full((stage["H"], stage["W"]), 80 + stage_index + delta, dtype=np.uint8)
        diag = np.ones_like(image)
        tifffile.imwrite(directory / f"{run_id}_{stage_name}_rawtex.tif", image)
        tifffile.imwrite(directory / f"{run_id}_{stage_name}_diagclass.tif", diag)
        plt.imsave(
            directory / f"{run_id}_{stage_name}_rawtex_strip.png", image, cmap="gray", vmin=0, vmax=255
        )
        plt.imsave(
            directory / f"{run_id}_{stage_name}_rawtex_preview.png", image, cmap="gray", vmin=0, vmax=255
        )
    return module.inspect_artifacts(directory, run_id, stages, run_id)


def test_downstream_pass_and_complete_render(tmp_path):
    module = _load_script()
    parent = _stats(seam1=1.80, seam5=1.30, dcol1=10.0)
    candidate = _stats(seam1=1.40, seam5=1.28, dcol1=8.0)
    repeat = copy.deepcopy(candidate)
    repeat["total_seconds"] = 12.0
    for stage in repeat["stages"]:
        stage["seconds"] += 0.5

    parent_inventory = _write_artifacts(module, tmp_path / "parent", "parent", parent, delta=2)
    candidate_inventory = _write_artifacts(module, tmp_path / "candidate", "candidate", candidate)
    repeat_inventory = _write_artifacts(module, tmp_path / "repeat", "repeat", repeat)
    result = module.compare(
        parent,
        candidate,
        repeat,
        parent_inventory,
        candidate_inventory,
        repeat_inventory,
    )
    assert result["decision"] == "DOWNSTREAM_EFFICACY_PASS"
    assert result["candidate_repeat"]["pass"] is True
    assert result["stage_effects"]["step1_ribbon"]["change"]["seam_excess_reduction"] == pytest.approx(
        0.4
    )
    png, svg = module.render(
        result,
        tmp_path / "parent",
        tmp_path / "candidate",
        "parent",
        "candidate",
        tmp_path / "complete_strips",
    )
    assert png.is_file() and png.stat().st_size > 1000
    assert svg.is_file() and svg.stat().st_size > 1000


def test_repeat_pixel_change_is_inconclusive(tmp_path):
    module = _load_script()
    parent = _stats(seam1=1.80, seam5=1.30, dcol1=10.0)
    candidate = _stats(seam1=1.40, seam5=1.28, dcol1=8.0)
    repeat = copy.deepcopy(candidate)
    parent_inventory = _write_artifacts(module, tmp_path / "parent", "parent", parent, delta=2)
    candidate_inventory = _write_artifacts(module, tmp_path / "candidate", "candidate", candidate)
    repeat_inventory = _write_artifacts(module, tmp_path / "repeat", "repeat", repeat, delta=1)
    result = module.compare(
        parent,
        candidate,
        repeat,
        parent_inventory,
        candidate_inventory,
        repeat_inventory,
    )
    assert result["decision"] == "DOWNSTREAM_INCONCLUSIVE"
    assert result["candidate_repeat"]["required_artifacts_byte_identical"] is False


def test_rejects_changed_topology(tmp_path):
    module = _load_script()
    parent = _stats(seam1=1.80, seam5=1.30, dcol1=10.0)
    candidate = _stats(seam1=1.40, seam5=1.28, dcol1=8.0)
    repeat = copy.deepcopy(candidate)
    candidate["nv"] += 1
    parent_inventory = _write_artifacts(module, tmp_path / "parent", "parent", parent, delta=2)
    candidate_inventory = _write_artifacts(module, tmp_path / "candidate", "candidate", candidate)
    repeat_inventory = _write_artifacts(module, tmp_path / "repeat", "repeat", repeat)
    result = module.compare(
        parent,
        candidate,
        repeat,
        parent_inventory,
        candidate_inventory,
        repeat_inventory,
    )
    assert result["decision"] == "DOWNSTREAM_INCONCLUSIVE"
    assert result["comparability"]["nv"] is False


def test_artifact_shape_must_match_stats(tmp_path):
    module = _load_script()
    stats = _stats(seam1=1.4, seam5=1.2, dcol1=8.0)
    directory = tmp_path / "bad"
    inventory = _write_artifacts(module, directory, "bad", stats)
    del inventory
    tifffile.imwrite(
        directory / "bad_step1_ribbon_rawtex.tif", np.zeros((7, 12), dtype=np.uint8)
    )
    with pytest.raises(ValueError, match="shape"):
        module.inspect_artifacts(directory, "bad", module._stage_map(stats, "bad"), "bad")
