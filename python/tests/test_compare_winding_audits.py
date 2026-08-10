import importlib.util
import json
from pathlib import Path

import pytest


def _load_script():
    path = Path(__file__).resolve().parents[1] / "scripts" / "compare_winding_audits.py"
    spec = importlib.util.spec_from_file_location("compare_winding_audits_tested", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _audit(turns, collisions, *, status, lt2):
    ids = (
        "z007936_y003968_x002944",
        "z007936_y004096_x002944",
        "z007936_y004096_x003072",
    )
    pairs = [
        {
            "a": ids[0],
            "b": ids[1],
            "n": 50,
            "dphi_med": 0.01,
            "du_med": 0.1,
            "du_mad": 0.2,
            "dv_absmed": 1.0,
            "turn_off": turns[0],
            "collide": collisions[0],
        },
        {
            "a": ids[1],
            "b": ids[2],
            "n": 50,
            "dphi_med": -0.01,
            "du_med": -0.1,
            "du_mad": 0.3,
            "dv_absmed": 1.1,
            "turn_off": turns[1],
            "collide": collisions[1],
        },
    ]
    return {
        "pairs": pairs,
        "n_pair_sets": 2,
        "n_pairs": 100,
        "turn_off_pairs": sum(turns),
        "collisions": sum(collisions),
        "join_completeness": {"lt2": lt2, "lt6": 0.9, "n": 100},
        "quality_gate": {"status": status},
    }


def _index():
    cubes = []
    for cube_id in (
        "z007936_y003968_x002944",
        "z007936_y004096_x002944",
        "z007936_y004096_x003072",
    ):
        z, y, x = map(int, (cube_id[1:7], cube_id[9:15], cube_id[17:23]))
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
        "axis_point_zyx": [0.0, 4296.0, 3324.0],
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


def test_pass_and_complete_render(tmp_path):
    module = _load_script()
    baseline = _audit((10, 10), (5, 4), status="FAIL", lt2=0.80)
    candidate = _audit((2, 3), (2, 1), status="PASS", lt2=0.82)
    result = module.compare(baseline, candidate, _index(), _index())
    assert result["decision"] == "INDEPENDENT_EFFICACY_PASS"
    assert result["paired_effect"]["turn_off_absolute_reduction_pp"] == pytest.approx(15.0)
    assert result["paired_effect"]["pair_sets_improved_worse_equal"] == [2, 0, 0]
    png, svg = module.render(result, tmp_path / "all_pairs")
    assert png.is_file() and png.stat().st_size > 1000
    assert svg.is_file() and svg.stat().st_size > 1000


def test_rejects_changed_pair_inventory():
    module = _load_script()
    baseline = _audit((10, 10), (5, 4), status="FAIL", lt2=0.80)
    candidate = _audit((2, 3), (2, 1), status="PASS", lt2=0.82)
    candidate["pairs"][0]["n"] = 49
    candidate["n_pairs"] = 99
    candidate["join_completeness"]["n"] = 99
    with pytest.raises(ValueError, match="different counts"):
        module.compare(baseline, candidate, _index(), _index())


def test_collision_gate_count_may_exceed_pair_gate_count():
    module = _load_script()
    baseline = _audit((10, 10), (75, 4), status="FAIL", lt2=0.80)
    candidate = _audit((2, 3), (71, 1), status="PASS", lt2=0.82)
    result = module.compare(baseline, candidate, _index(), _index())
    assert result["baseline"]["collisions"] == 79
    assert result["candidate"]["collisions"] == 72
    assert result["gates"]["collisions_nonincreasing"]


def test_cli_writes_hashed_machine_result(tmp_path, monkeypatch):
    module = _load_script()
    paths = {}
    for name, value in {
        "baseline_audit": _audit((10, 10), (5, 4), status="FAIL", lt2=0.80),
        "candidate_audit": _audit((2, 3), (2, 1), status="PASS", lt2=0.82),
        "baseline_index": _index(),
        "candidate_index": _index(),
    }.items():
        path = tmp_path / f"{name}.json"
        path.write_text(json.dumps(value), encoding="utf-8")
        paths[name] = path
    out_json = tmp_path / "result.json"
    argv = [
        "compare_winding_audits.py",
        "--baseline-audit",
        str(paths["baseline_audit"]),
        "--candidate-audit",
        str(paths["candidate_audit"]),
        "--baseline-index",
        str(paths["baseline_index"]),
        "--candidate-index",
        str(paths["candidate_index"]),
        "--out-json",
        str(out_json),
        "--figure-prefix",
        str(tmp_path / "figure"),
    ]
    monkeypatch.setattr("sys.argv", argv)
    module.main()
    result = json.loads(out_json.read_text(encoding="utf-8"))
    assert result["decision"] == "INDEPENDENT_EFFICACY_PASS"
    assert set(result["evidence"]["input_sha256"]) == set(paths)
    assert result["evidence"]["visual_scope"].startswith("Every audited pair")
