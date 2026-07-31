"""stdin_tiff_parity.py -- prove cube_mesh --stdin-raw == TIFF+halo mode.

Assembles the (128+2*halo)^3 padded volume for one cube exactly like
HaloLoader_load (center TIFF + 26 neighbor TIFFs from the same cubes_PRED
dir, missing neighbors -> 0), pipes it to `cube_mesh --stdin-raw`, runs the
same cube through the normal TIFF+--halo path, and compares the
step12_final dumps. Run single-threaded (VESUVIUS_THREADS=1) so both paths
are deterministic.

Usage:
  python scripts/extract/stdin_tiff_parity.py \
      [--grid PHerc0139-4x5x5] [--cube z04480_y03328_x02816] [--halo 13] \
      [--bin build/Release/cube_mesh.exe] [--out output/stdin_tiff_parity]

Requires numpy + tifffile (e.g. the scrollunwrap uv env:
  py -m uv run --project python python scripts/extract/stdin_tiff_parity.py).
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import tifffile

CUBE = 128


def parse_id(cube_id: str) -> tuple[int, int, int]:
    m = re.fullmatch(r"z(\d{5})_y(\d{5})_x(\d{5})", cube_id)
    if not m:
        sys.exit(f"bad cube id: {cube_id}")
    return tuple(int(g) for g in m.groups())


def build_padded(pred_dir: Path, cube_id: str, halo: int) -> np.ndarray:
    """Padded (CUBE+2*halo)^3 uint8 volume, index (0,0,0) = world origin-halo."""
    oz, oy, ox = parse_id(cube_id)
    p = CUBE + 2 * halo
    vol = np.zeros((p, p, p), np.uint8)
    lo = np.array([oz - halo, oy - halo, ox - halo])   # world coord of padded (0,0,0)
    hi = lo + p
    n_used = 0
    for dz in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                nz, ny, nx = oz + dz * CUBE, oy + dy * CUBE, ox + dx * CUBE
                if nz < 0 or ny < 0 or nx < 0:
                    continue
                tif = pred_dir / f"z{nz:05d}_y{ny:05d}_x{nx:05d}.tif"
                if not tif.exists():
                    continue
                nlo = np.array([nz, ny, nx])
                a = np.maximum(lo, nlo)                # world intersection
                b = np.minimum(hi, nlo + CUBE)
                if np.any(a >= b):
                    continue
                cube = tifffile.imread(tif)
                assert cube.shape == (CUBE, CUBE, CUBE), (tif, cube.shape)
                s_dst = tuple(slice(a[i] - lo[i], b[i] - lo[i]) for i in range(3))
                s_src = tuple(slice(a[i] - nlo[i], b[i] - nlo[i]) for i in range(3))
                vol[s_dst] = cube[s_src]
                n_used += 1
    print(f"  padded volume: {p}^3, {n_used}/27 cubes contributed, "
          f"{int((vol > 0).sum())} fg voxels")
    return vol


def final_obj(dump_root: Path, cube_id: str) -> Path:
    return (dump_root / cube_id / f"{cube_id}_step12_final"
            / f"{cube_id}_step12_final_all.obj")


def obj_stats(path: Path) -> tuple[int, int]:
    nv = nf = 0
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b"v "):
                nv += 1
            elif line.startswith(b"f "):
                nf += 1
    return nv, nf


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="PHerc0139-4x5x5")
    ap.add_argument("--cube", default="z04480_y03328_x02816")
    ap.add_argument("--halo", type=int, default=13)
    ap.add_argument("--bin", default="build/Release/cube_mesh.exe")
    ap.add_argument("--out", default="output/stdin_tiff_parity")
    args = ap.parse_args()

    args.bin = str(Path(args.bin).resolve())   # CreateProcess dislikes rel/fwd-slash
    if not Path(args.bin).exists():
        sys.exit(f"missing binary: {args.bin}")
    pred_dir = Path(args.grid) / "cubes_PRED"
    tif = pred_dir / f"{args.cube}.tif"
    if not tif.exists():
        sys.exit(f"missing {tif}")
    out = Path(args.out)
    (out / "tiff").mkdir(parents=True, exist_ok=True)
    (out / "stdin").mkdir(parents=True, exist_ok=True)

    env = dict(os.environ, VESUVIUS_THREADS="1")
    oz, oy, ox = parse_id(args.cube)
    p = CUBE + 2 * args.halo

    print(f"[1/3] TIFF+halo run ({args.cube}, halo={args.halo}, 1 thread)")
    r = subprocess.run([args.bin, str(tif), str(out / "tiff" / "unused.tif"),
                        "--halo", str(args.halo), "--dump-obj", str(out / "tiff")],
                       env=env, capture_output=True)
    if r.returncode != 0:
        sys.exit("TIFF run failed:\n" + r.stderr.decode(errors="replace")[-2000:])

    print(f"[2/3] stdin-raw run (p={p})")
    vol = build_padded(pred_dir, args.cube, args.halo)
    r = subprocess.run([args.bin, "--stdin-raw", str(p), str(oz), str(oy), str(ox),
                        "--halo", str(args.halo), "--dump-obj", str(out / "stdin")],
                       input=vol.tobytes(), env=env, capture_output=True)
    if r.returncode != 0:
        sys.exit("stdin run failed:\n" + r.stderr.decode(errors="replace")[-2000:])

    print("[3/3] compare step12_final")
    fa, fb = final_obj(out / "tiff", args.cube), final_obj(out / "stdin", args.cube)
    for f in (fa, fb):
        if not f.exists():
            sys.exit(f"missing dump {f}")
    va, ta = obj_stats(fa)
    vb, tb = obj_stats(fb)
    same_bytes = fa.read_bytes() == fb.read_bytes()
    print(f"  tiff : {va} verts {ta} faces  {fa}")
    print(f"  stdin: {vb} verts {tb} faces  {fb}")
    print(f"  byte-identical: {same_bytes}")
    if (va, ta) != (vb, tb):
        print("PARITY FAIL: vert/face counts differ")
        return 1
    print("PARITY OK" + ("" if same_bytes else " (counts equal, bytes differ)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
