"""Drive the scrollfiesta mesher directly from streamed zarr cubes.

For each owned cube we stream a padded block from the zarr (single chunked
slice, halo included) and pipe its bytes to ``cube_mesh --stdin-raw`` — no TIFF
cubes on disk. cube_mesh writes its per-cube OBJ dump; ``grid_weld`` then
stitches those dumps into ``welded.obj``.

stdin-raw contract (must match the C patch in main.c / mesh_extract.c):
    cube_mesh --stdin-raw <p_size> <oz> <oy> <ox> --halo <h> --dump-obj <dir> [--no-qem]
    stdin = p_size**3 bytes of uint8, C-order (z,y,x); buffer index (0,0,0) maps
    to world voxel (oz-h, oy-h, ox-h); cube owned origin is (oz,oy,ox).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from tqdm import tqdm

from ._paths import default_env
from .adaptive_bpa import (AdaptiveBpaConfig, apply_acceptance, audit_candidate,
                           choose_candidate, diagnose_failure, final_obj)
from .cube_planner import Cube
from .zarr_source import ZarrVolume, apply_threshold, read_padded_cube


@dataclass
class CubeResult:
    cube_id: str
    returncode: int
    produced_obj: bool
    stderr_tail: str
    adaptive_report: dict | None = None


@dataclass
class MesherResult:
    welded_obj: Path
    weld_report: dict | None
    dump_dir: Path
    n_cubes: int
    n_ok: int
    n_obj: int


def _mesh_one_cube(vol: ZarrVolume, cube: Cube, dump_dir: Path, *, cube_mesh_bin,
                   halo: int, threshold: str | None, skip_qem: bool, size: int,
                   env: dict, timeout: float,
                   adaptive_bpa: AdaptiveBpaConfig | None = None) -> CubeResult:
    p = size + 2 * halo
    try:
        buf = read_padded_cube(vol, cube.oz, cube.oy, cube.ox, size=size, halo=halo)
        mask = apply_threshold(buf, threshold)          # uint8 {0,255}, C-contiguous
    except subprocess.TimeoutExpired:
        return CubeResult(cube.cube_id, -9, False, f"TIMEOUT after {timeout}s")
    except Exception as e:  # e.g. an S3/credential error (expired STS token)
        return CubeResult(cube.cube_id, -1, False,
                          f"{type(e).__name__}: {str(e)[:300]}")

    payload = mask.tobytes()

    def run_once(one_dump: Path, rho: float | None):
        one_dump.mkdir(parents=True, exist_ok=True)
        cmd = [str(cube_mesh_bin), "--stdin-raw", str(p),
               str(cube.oz), str(cube.oy), str(cube.ox),
               "--halo", str(halo), "--dump-obj", str(one_dump)]
        if skip_qem:
            cmd.append("--no-qem")
        run_env = dict(env)
        if rho is not None:
            run_env["BPA_RHO"] = f"{rho:.6g}"
        return subprocess.run(cmd, input=payload, capture_output=True,
                              env=run_env, timeout=timeout)

    if adaptive_bpa is None:
        exp = (dump_dir / cube.cube_id / f"{cube.cube_id}_step12_final"
               / f"{cube.cube_id}_step12_final_all.obj")
        try:
            proc = run_once(dump_dir, None)
        except subprocess.TimeoutExpired:
            return CubeResult(cube.cube_id, -9, exp.exists(),
                              f"TIMEOUT after {timeout}s")
        tail = proc.stderr.decode("utf-8", "replace")[-400:] if proc.stderr else ""
        return CubeResult(cube.cube_id, proc.returncode, exp.exists(), tail)

    # Adaptive mode runs candidates outside the canonical dump tree. Only the
    # selected cube directory is copied into dump/, so grid_weld sees exactly
    # one internally consistent run.
    candidate_root = dump_dir.parent / "adaptive_bpa_candidates" / cube.cube_id
    candidates: list[dict] = []
    baseline_eval: dict | None = None
    last_tail = ""
    for candidate_index, rho in enumerate(adaptive_bpa.rhos):
        tag = f"rho_{rho:.2f}".replace(".", "p")
        candidate_dump = candidate_root / tag / "dump"
        candidate_cube = candidate_dump / cube.cube_id
        try:
            proc = run_once(candidate_dump, rho)
            last_tail = (proc.stderr.decode("utf-8", "replace")[-400:]
                         if proc.stderr else "")
        except subprocess.TimeoutExpired:
            candidates.append({"rho": rho, "returncode": -9,
                               "produced_obj": False,
                               "error": f"TIMEOUT after {timeout}s"})
            continue

        produced = final_obj(candidate_cube) is not None
        entry: dict = {"rho": rho, "returncode": proc.returncode,
                       "produced_obj": produced,
                       "candidate_cube_dir": str(candidate_cube)}
        if produced:
            try:
                entry.update(audit_candidate(candidate_cube, adaptive_bpa))
                reference = None if candidate_index == 0 else baseline_eval
                apply_acceptance(entry, reference, adaptive_bpa)
                # A lower-radius run cannot pass the preservation gates if the
                # requested baseline did not produce a usable audit.
                if candidate_index > 0 and baseline_eval is None:
                    entry.update({"point_coverage": 0.0,
                                  "coverage_ok": False,
                                  "boundary_ratio": None,
                                  "boundary_ok": False,
                                  "passes": False,
                                  "acceptance_error": "baseline audit unavailable"})
            except Exception as e:
                entry.update({"topology_clean": False, "passes": False,
                              "audit_error": f"{type(e).__name__}: {e}"})
            if candidate_index == 0:
                # Keep even a failed baseline audit as the explicit reference;
                # missing metrics make later candidates fail conservatively.
                baseline_eval = entry
        else:
            entry["passes"] = False
        candidates.append(entry)

        # Highest-rho passing candidate is preferred. A clean baseline avoids
        # any lower-radius reruns; a flagged baseline triggers 1.10 then 1.00.
        if entry.get("passes"):
            break

    try:
        selected, accepted = choose_candidate(candidates)
    except Exception as e:
        if candidate_root.exists() and not adaptive_bpa.keep_candidates:
            shutil.rmtree(candidate_root)
        return CubeResult(cube.cube_id, -1, False,
                          f"adaptive BPA failed: {e}")

    selected_cube = Path(selected["candidate_cube_dir"])
    canonical_cube = dump_dir / cube.cube_id
    if canonical_cube.exists():
        shutil.rmtree(canonical_cube)
    canonical_cube.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(selected_cube, canonical_cube)

    # Candidate paths are ephemeral; retain stage filenames and all metrics in
    # the durable report without leaving misleading dead absolute paths.
    for candidate in candidates:
        candidate.pop("candidate_cube_dir", None)
        for stage in candidate.get("stages", []):
            if "obj" in stage:
                stage["obj"] = Path(stage["obj"]).name
    fallback_to_baseline = bool(
        not accepted and candidates and selected is candidates[0])
    report = {
        "cube_id": cube.cube_id,
        "policy": {
            "rhos": list(adaptive_bpa.rhos),
            "local_span_tol": adaptive_bpa.local_span_tol,
            "min_point_coverage": adaptive_bpa.min_point_coverage,
            "point_tolerance_vox": adaptive_bpa.point_tolerance,
            "max_boundary_ratio": adaptive_bpa.max_boundary_ratio,
        },
        "selected_rho": selected["rho"],
        "accepted_clean_candidate": accepted,
        "no_passing_candidate": not accepted,
        "failure_mode": diagnose_failure(candidates),
        "fallback_to_baseline": fallback_to_baseline,
        "candidates": candidates,
    }
    (canonical_cube / "adaptive_bpa.json").write_text(
        json.dumps(report, indent=2, allow_nan=False))
    if candidate_root.exists() and not adaptive_bpa.keep_candidates:
        shutil.rmtree(candidate_root)

    produced = final_obj(canonical_cube) is not None
    selected_rc = int(selected.get("returncode", 0))
    return CubeResult(cube.cube_id, selected_rc, produced, last_tail, report)


def run_grid_weld(dump_dir, welded_obj, *, grid_weld_bin, env, timeout=1800.0):
    cmd = [str(grid_weld_bin), str(dump_dir), str(welded_obj)]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout)
    welded_obj = Path(welded_obj)
    has_verts = welded_obj.exists() and any(
        ln.startswith("v ") for ln in open(welded_obj))
    if not has_verts:
        raise RuntimeError(
            f"grid_weld produced no welded.obj (rc={proc.returncode})\n"
            f"cmd: {' '.join(cmd)}\nstdout:\n{proc.stdout[-2000:]}\n"
            f"stderr:\n{proc.stderr[-2000:]}")
    if proc.returncode != 0:
        # grid_weld exits nonzero when its manifold audit flags residual
        # non-manifold edges, but it still writes a usable welded.obj; meshprep
        # cleans those edges downstream. Treat the audit as a warning.
        print(f"  [warn] grid_weld audit rc={proc.returncode} (using welded.obj anyway)")
    rep = Path(f"{welded_obj}.weld_report.json")
    return json.loads(rep.read_text()) if rep.exists() else None


def run_mesher(vol: ZarrVolume, cubes: list[Cube], mesh_dir, *, cube_mesh_bin,
               grid_weld_bin, halo: int = 13, threshold: str | None = None,
               skip_qem: bool = False, size: int = 128, max_concurrent: int | None = None,
               threads_per_cube: int = 1, cube_timeout: float = 1200.0,
               adaptive_bpa: AdaptiveBpaConfig | None = None,
               mls_backend: str | None = None) -> MesherResult:
    mesh_dir = Path(mesh_dir)
    dump_dir = mesh_dir / "dump"
    dump_dir.mkdir(parents=True, exist_ok=True)
    env = default_env(threads_per_cube)
    if mls_backend:
        env["MLS_BACKEND"] = mls_backend
    if max_concurrent is None:
        max_concurrent = max(1, (os.cpu_count() or 4) // max(1, threads_per_cube))

    results: list[CubeResult] = []
    with ThreadPoolExecutor(max_workers=max_concurrent) as ex:
        futs = [ex.submit(_mesh_one_cube, vol, c, dump_dir, cube_mesh_bin=cube_mesh_bin,
                          halo=halo, threshold=threshold, skip_qem=skip_qem, size=size,
                          env=env, timeout=cube_timeout,
                          adaptive_bpa=adaptive_bpa) for c in cubes]
        for f in tqdm(as_completed(futs), total=len(futs), desc="mesh cubes"):
            results.append(f.result())

    n_ok = sum(r.returncode == 0 for r in results)
    n_obj = sum(r.produced_obj for r in results)
    if n_obj == 0:
        fails = [r for r in results if r.returncode != 0][:3]
        raise RuntimeError(
            f"no cube produced a mesh ({n_ok}/{len(cubes)} exited 0). "
            f"Region may be empty, or the binary failed. Sample stderr: "
            + " | ".join(r.stderr_tail for r in fails))

    if adaptive_bpa is not None:
        reports = [r.adaptive_report for r in results if r.adaptive_report]
        summary = {
            "cubes": len(results),
            "audited": len(reports),
            "clean_candidate_selected": sum(
                bool(r.get("accepted_clean_candidate")) for r in reports),
            "baseline_fallbacks": sum(bool(r.get("fallback_to_baseline")) for r in reports),
            "no_passing_candidate": sum(bool(r.get("no_passing_candidate")) for r in reports),
            "failure_modes": {
                mode: sum(r.get("failure_mode") == mode for r in reports)
                for mode in sorted({r.get("failure_mode") for r in reports
                                    if r.get("failure_mode")})
            },
            "reports": reports,
        }
        (mesh_dir / "adaptive_bpa_summary.json").write_text(
            json.dumps(summary, indent=2, allow_nan=False))

    welded = mesh_dir / "welded.obj"
    report = run_grid_weld(dump_dir, welded, grid_weld_bin=grid_weld_bin, env=env)
    return MesherResult(welded, report, dump_dir, len(cubes), n_ok, n_obj)
