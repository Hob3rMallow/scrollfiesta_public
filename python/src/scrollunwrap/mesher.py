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
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from tqdm import tqdm

from ._paths import default_env
from .cube_planner import Cube
from .zarr_source import ZarrVolume, apply_threshold, read_padded_cube


@dataclass
class CubeResult:
    cube_id: str
    returncode: int
    produced_obj: bool
    stderr_tail: str


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
                   env: dict, timeout: float) -> CubeResult:
    p = size + 2 * halo
    cmd = [str(cube_mesh_bin), "--stdin-raw", str(p),
           str(cube.oz), str(cube.oy), str(cube.ox),
           "--halo", str(halo), "--dump-obj", str(dump_dir)]
    if skip_qem:
        cmd.append("--no-qem")
    exp = (dump_dir / cube.cube_id / f"{cube.cube_id}_raw_snap"
           / f"{cube.cube_id}_raw_snap_all.obj")
    try:
        buf = read_padded_cube(vol, cube.oz, cube.oy, cube.ox, size=size, halo=halo)
        mask = apply_threshold(buf, threshold)          # uint8 {0,255}, C-contiguous
        proc = subprocess.run(cmd, input=mask.tobytes(), capture_output=True,
                              env=env, timeout=timeout)
    except subprocess.TimeoutExpired:
        # A pathologically dense cube can stall the mesher; skip it rather than
        # hang the whole run (try a coarser --level or smaller ROI).
        return CubeResult(cube.cube_id, -9, exp.exists(), f"TIMEOUT after {timeout}s")
    except Exception as e:  # e.g. an S3/credential error (expired STS token)
        # Skip this cube so the run welds whatever else succeeded instead of
        # crashing the whole ROI. The message surfaces in run_mesher's report.
        return CubeResult(cube.cube_id, -1, exp.exists(),
                          f"{type(e).__name__}: {str(e)[:300]}")
    tail = proc.stderr.decode("utf-8", "replace")[-400:] if proc.stderr else ""
    return CubeResult(cube.cube_id, proc.returncode, exp.exists(), tail)


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
               threads_per_cube: int = 1, cube_timeout: float = 1200.0) -> MesherResult:
    mesh_dir = Path(mesh_dir)
    dump_dir = mesh_dir / "dump"
    dump_dir.mkdir(parents=True, exist_ok=True)
    env = default_env(threads_per_cube)
    if max_concurrent is None:
        max_concurrent = max(1, (os.cpu_count() or 4) // max(1, threads_per_cube))

    results: list[CubeResult] = []
    with ThreadPoolExecutor(max_workers=max_concurrent) as ex:
        futs = [ex.submit(_mesh_one_cube, vol, c, dump_dir, cube_mesh_bin=cube_mesh_bin,
                          halo=halo, threshold=threshold, skip_qem=skip_qem, size=size,
                          env=env, timeout=cube_timeout) for c in cubes]
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

    welded = mesh_dir / "welded.obj"
    report = run_grid_weld(dump_dir, welded, grid_weld_bin=grid_weld_bin, env=env)
    return MesherResult(welded, report, dump_dir, len(cubes), n_ok, n_obj)
