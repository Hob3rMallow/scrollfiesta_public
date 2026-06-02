"""Flatten a prepared OBJ with villa's C++ `flatboi` (SLIM).

`flatboi <mesh.obj> [iters] [energy] [--tol=x]` writes `<stem>_flatboi.obj`
(with per-corner UVs) plus UV heatmap PNGs next to the input. With the harmonic
fallback patch it computes the initial UV itself, so the input needs no UVs.
"""

from __future__ import annotations

import subprocess
from pathlib import Path


def flatboi_output_path(prepped_obj) -> Path:
    p = Path(prepped_obj)
    return p.parent / f"{p.stem}_flatboi.obj"


def run_flatboi(prepped_obj, *, binary, iters: int = 50,
                energy: str = "symmetric_dirichlet", tol: float | None = None,
                env=None, timeout: float = 3600.0) -> Path:
    prepped_obj = Path(prepped_obj)
    cmd = [str(binary), str(prepped_obj), str(int(iters)), energy]
    if tol:
        cmd.append(f"--tol={tol}")
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout)
    out = flatboi_output_path(prepped_obj)
    if proc.returncode != 0 or not out.exists():
        raise RuntimeError(
            f"flatboi failed (rc={proc.returncode}, output_exists={out.exists()})\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout[-2000:]}\nstderr:\n{proc.stderr[-2000:]}")
    # sanity: the legacy converter requires `vt` lines downstream
    with open(out, "r") as fh:
        if not any(line.startswith("vt ") for line in fh):
            raise RuntimeError(f"flatboi output {out} has no UV (vt) lines")
    return out
