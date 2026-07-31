from pathlib import Path

from scrollunwrap.adaptive_bpa import (AdaptiveBpaConfig, apply_acceptance,
                                       choose_candidate, diagnose_failure,
                                       point_coverage)


def _obj(path: Path, vertices):
    path.write_text("".join(f"v {z} {y} {x}\n" for z, y, x in vertices))
    return path


def _evaluation(obj: Path, *, boundary=100, clean=True):
    return {
        "topology_clean": clean,
        "stages": [{"stage": "step1_bpa", "obj": str(obj),
                    "boundary_edges": boundary}],
    }


def test_point_coverage_uses_physical_tolerance(tmp_path):
    reference = _obj(tmp_path / "ref.obj", [(0, 0, 0), (0, 1, 0), (0, 2, 0)])
    candidate = _obj(tmp_path / "candidate.obj", [(0, 0.2, 0), (0, 1.2, 0)])
    assert point_coverage(reference, candidate, 0.3) == 2 / 3


def test_acceptance_enforces_coverage_and_boundary(tmp_path):
    ref = _obj(tmp_path / "ref.obj", [(0, i, 0) for i in range(20)])
    same = _obj(tmp_path / "same.obj", [(0, i, 0) for i in range(20)])
    cfg = AdaptiveBpaConfig(Path("wind_audit"), 0, 0)
    baseline = apply_acceptance(_evaluation(ref), None, cfg)
    good = apply_acceptance(_evaluation(same, boundary=124), baseline, cfg)
    assert good["point_coverage"] == 1.0
    assert good["boundary_ratio"] == 1.24
    assert good["passes"]

    too_open = apply_acceptance(_evaluation(same, boundary=126), baseline, cfg)
    assert not too_open["passes"]


def test_candidate_choice_prefers_first_clean_and_falls_back_to_baseline():
    baseline = {"rho": 1.2, "produced_obj": True, "passes": False}
    clean = {"rho": 1.1, "produced_obj": True, "passes": True}
    chosen, accepted = choose_candidate([baseline, clean])
    assert accepted and chosen["rho"] == 1.1

    chosen, accepted = choose_candidate([baseline, {"rho": 1.0,
                                                     "produced_obj": True,
                                                     "passes": False}])
    assert not accepted and chosen is baseline


def test_failure_diagnosis_separates_radius_insensitive_fusion():
    def candidate(rho, count, span):
        return {
            "rho": rho, "returncode": 0, "produced_obj": True,
            "passes": False,
            "stages": [{
                "stage": "step1_bpa", "fused": True,
                "broad_fusion_components": count,
                "broad_fusion_max_span": span,
            }],
        }

    stable = [candidate(1.2, 9, 18.70), candidate(1.1, 9, 18.68),
              candidate(1.0, 9, 18.71)]
    assert diagnose_failure(stable) == "radius_insensitive_fusion"

    improving = [candidate(1.2, 9, 18.7), candidate(1.1, 4, 8.0)]
    assert diagnose_failure(improving) == "no_clean_radius"

    stable[-1]["passes"] = True
    assert diagnose_failure(stable) is None
