import importlib.util
import json
from pathlib import Path

import pytest


def _load_script():
    path = Path(__file__).resolve().parents[1] / "scripts" / "compare_winding_scale_stress.py"
    spec = importlib.util.spec_from_file_location("compare_winding_scale_stress_tested", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


IDS = (
    "z12800_y04480_x03840",  # the one-cube known box used by these tests
    "z12672_y04480_x03840",
    "z12672_y04608_x03840",
    "z12672_y04608_x03968",
    "z12800_y04608_x03840",
)
CORE = (12672, 12928, 4480, 4736, 3840, 4096)
KNOWN = (12800, 12928, 4480, 4608, 3840, 3968)


def _audit(turns, collisions, *, status, lt2):
    definitions = (
        (IDS[0], IDS[1], 100),
        (IDS[1], IDS[2], 300),
        (IDS[2], IDS[3], 300),
        (IDS[0], IDS[4], 100),
    )
    pairs = []
    for position, (a, b, n) in enumerate(definitions):
        pairs.append(
            {
                "a": a,
                "b": b,
                "n": n,
                "dphi_med": 0.01,
                "du_med": 0.1,
                "du_mad": 0.2,
                "dv_absmed": 1.0,
                "turn_off": turns[position],
                "collide": collisions[position],
            }
        )
    total = sum(pair["n"] for pair in pairs)
    return {
        "pairs": pairs,
        "n_pair_sets": len(pairs),
        "n_pairs": total,
        "turn_off_pairs": sum(turns),
        "collisions": sum(collisions),
        "join_completeness": {"lt2": lt2, "lt6": 0.9, "n": total},
        "quality_gate": {"status": status},
    }


def _index(ids=IDS):
    cubes = []
    for cube_id in ids:
        z, y, x = (int(part[1:]) for part in cube_id.split("_"))
        cubes.append(
            {
                "id": cube_id,
                "origin": [z, y, x],
                "status": 0,
                "nv": 10,
                "nf": 12,
                "skin": 5,
                "flood_comp": 0,
                "n_groups": 2,
            }
        )
    return {
        "leaf_stage": "step12_final",
        "chunk": 128,
        "axis_point_zyx": [0.0, 4829.0, 4118.0],
        "axis_dir_zyx": [1.0, 0.0, 0.0],
        "pitch": 8.2,
        "pitch_mode": "auto",
        "pair_gate": 3.5,
        "skin_dist": 4.0,
        "calibration": {"seed_id": cubes[0]["id"], "spiral_a": 1.0, "spiral_b": 8.2, "sense": 1},
        "reregistered": 1,
        "n_cubes": len(cubes),
        "n_ok": len(cubes),
        "n_skipped": 0,
        "n_low_conf": 0,
        "cubes": cubes,
    }


def _compare(module, baseline=None, candidate=None, repeat=None):
    baseline = baseline or _audit((4, 30, 30, 4), (1, 4, 4, 1), status="FAIL", lt2=0.80)
    candidate = candidate or _audit((2, 6, 6, 2), (1, 2, 2, 1), status="PASS", lt2=0.82)
    return module.compare_scale(
        baseline,
        candidate,
        _index(),
        _index(),
        repeat or {"pass": True},
        core_bbox=CORE,
        known_bbox=KNOWN,
        min_core_cubes=5,
        min_total_pairs=500,
        min_new_pair_sets=2,
        min_new_pairs=500,
    )


def test_new_only_stratum_controls_pass_and_complete_3d_render(tmp_path):
    module = _load_script()
    result = _compare(module)
    assert result["decision"] == "SCALE_STRESS_PASS"
    assert result["strata"]["new_only"]["audit_pairs"] == 600
    assert result["strata"]["new_only"]["turn_off_absolute_reduction_pp"] == pytest.approx(8.0)
    assert result["strata"]["bridge"]["audit_pairs"] == 200
    assert result["strata"]["known_center"]["audit_pairs"] == 0
    png, svg = module.render(result, tmp_path / "all_pairs_3d")
    assert png.is_file() and png.stat().st_size > 1000
    assert svg.is_file() and svg.stat().st_size > 1000


def test_new_only_worsening_is_regression_even_if_bridge_improves():
    module = _load_script()
    candidate = _audit((0, 45, 45, 0), (0, 2, 2, 0), status="PASS", lt2=0.82)
    result = _compare(module, candidate=candidate)
    assert result["decision"] == "SCALE_STRESS_REGRESSION"
    assert not result["gates"]["new_only_turnoff_worsening_within_floor"]


def test_incomplete_known_box_is_rejected():
    module = _load_script()
    index = _index(IDS[1:])
    with pytest.raises(ValueError, match="known box is incomplete"):
        module.compare_scale(
            _audit((4, 30, 30, 4), (1, 4, 4, 1), status="FAIL", lt2=0.80),
            _audit((2, 6, 6, 2), (1, 2, 2, 1), status="PASS", lt2=0.82),
            index,
            index,
            {"pass": True},
            core_bbox=CORE,
            known_bbox=KNOWN,
            min_core_cubes=4,
            min_total_pairs=500,
            min_new_pair_sets=2,
            min_new_pairs=500,
        )


def test_repeat_verifier_requires_exact_skin_inventory(tmp_path):
    module = _load_script()
    candidate = tmp_path / "candidate"
    repeat = tmp_path / "repeat"
    candidate.mkdir()
    repeat.mkdir()
    audit_path = candidate / "audit.json"
    repeat_audit_path = repeat / "audit.json"
    index_path = candidate / "placed_index.json"
    repeat_index_path = repeat / "placed_index.json"
    audit_bytes = json.dumps(_audit((4, 30, 30, 4), (1, 4, 4, 1), status="PASS", lt2=0.8)).encode()
    index_bytes = json.dumps(_index()).encode()
    audit_path.write_bytes(audit_bytes)
    repeat_audit_path.write_bytes(audit_bytes)
    index_path.write_bytes(index_bytes)
    repeat_index_path.write_bytes(index_bytes)
    for cube_id in IDS:
        (candidate / f"{cube_id}_skin.f32").write_bytes(b"same")
        (repeat / f"{cube_id}_skin.f32").write_bytes(b"same")
    result = module.verify_repeat(
        audit_path,
        index_path,
        repeat_audit_path,
        repeat_index_path,
        candidate,
        repeat,
    )
    assert result["pass"]
    (repeat / f"{IDS[-1]}_skin.f32").write_bytes(b"changed")
    assert not module.verify_repeat(
        audit_path,
        index_path,
        repeat_audit_path,
        repeat_index_path,
        candidate,
        repeat,
    )["pass"]
