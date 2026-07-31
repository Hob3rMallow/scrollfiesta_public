"""Carve a cube-grid directory (cubes_PRED/ + cubes_RAW/) from OME-Zarr stores.

Bridges the zarr world to the C pipeline's TIFF-grid input format: for a
level-0 voxel bbox (multiples of 128), stream the prediction (and optionally
the raw masked volume) in z-bands and write per-cube 128^3 multipage TIFFs
(128 pages of 128x128, uint8, uncompressed, little-endian) named
``z#####_y#####_x#####.tif`` — the exact layout ``grid_pipeline`` /
``HaloLoader`` / ``scroll_unroll --raw`` consume.

  * all-zero PRED cubes are skipped (grid_pipeline scans *.tif; absent cube =
    not meshed; HaloLoader treats a missing neighbour as unpadded), unless
    --keep-empty-pred;
  * RAW cubes are always written (the unroll bake samples free space too);
  * writes cubes_PRED/present.json + a minimal manifest.json (chunk_size,
    bbox, sources, optional umbilicus) for provenance.

Usage::

    py -m uv run --project python python/scripts/carve_grid_tifs.py \
        --pred-zarr s3://... --raw-zarr s3://... \
        --bbox 11904 12416 2304 4992 2432 5120 \
        --umbilicus 3624 3815 --out PHerc1447-4x21x21
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import tifffile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from scrollunwrap.zarr_source import open_volume, read_block  # noqa: E402

CUBE = 128


def carve_one(store: str, out_dir: Path, bbox, *, skip_empty: bool,
              threshold_binary: bool, anon) -> tuple[int, int, list[str], int]:
    """Stream bbox in z-bands of CUBE, write per-cube TIFFs.

    Returns (written, skipped_empty, ids, all_zero_cubes)."""
    z0, z1, y0, y1, x0, x1 = bbox
    vol = open_volume(store, level=0, anon=anon)
    D, H, W = vol.shape
    assert 0 <= z0 < z1 <= D and 0 <= y0 < y1 <= H and 0 <= x0 < x1 <= W, \
        f"bbox {bbox} outside volume {vol.shape}"
    out_dir.mkdir(parents=True, exist_ok=True)
    written = skipped = all_zero = 0
    ids: list[str] = []
    for bz in range(z0, z1, CUBE):
        t0 = time.perf_counter()
        band = read_block(vol, bz, bz + CUBE, y0, y1, x0, x1)
        t_read = time.perf_counter() - t0
        if threshold_binary:
            band = (band != 0).astype(np.uint8) * np.uint8(255)
        n_band = 0
        for by in range(0, y1 - y0, CUBE):
            for bx in range(0, x1 - x0, CUBE):
                cube = np.ascontiguousarray(
                    band[:, by:by + CUBE, bx:bx + CUBE])
                cid = f"z{bz:05d}_y{y0 + by:05d}_x{x0 + bx:05d}"
                is_empty = not cube.any()
                if is_empty:
                    all_zero += 1
                if skip_empty and is_empty:
                    skipped += 1
                    continue
                tifffile.imwrite(out_dir / f"{cid}.tif", cube,
                                 photometric="minisblack",
                                 compression=None, rowsperstrip=CUBE)
                written += 1
                n_band += 1
                ids.append(cid)
        print(f"  band z{bz:05d}: read {t_read:.1f}s, wrote {n_band} cubes "
              f"(total {written}, empty-skipped {skipped})", flush=True)
    return written, skipped, ids, all_zero


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-zarr", required=True)
    ap.add_argument("--raw-zarr", default=None)
    ap.add_argument("--bbox", nargs=6, type=int, required=True,
                    metavar=("Z0", "Z1", "Y0", "Y1", "X0", "X1"))
    ap.add_argument("--umbilicus", nargs=2, type=float, default=None,
                    metavar=("Y", "X"))
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--keep-empty-pred", action="store_true")
    ap.add_argument("--s3-anon", choices=["auto", "yes", "no"], default="yes")
    args = ap.parse_args(argv)

    bbox = args.bbox
    for v in bbox:
        if v % CUBE:
            ap.error(f"bbox coordinate {v} is not a multiple of {CUBE}")
    anon = {"auto": None, "yes": True, "no": False}[args.s3_anon]

    print(f"PRED <- {args.pred_zarr}")
    t0 = time.perf_counter()
    nw, ns, ids, _ = carve_one(args.pred_zarr, args.out / "cubes_PRED", bbox,
                            skip_empty=not args.keep_empty_pred,
                            threshold_binary=True, anon=anon)
    print(f"PRED: {nw} cubes written, {ns} empty skipped "
          f"[{time.perf_counter() - t0:.0f}s]")
    (args.out / "cubes_PRED" / "present.json").write_text(
        json.dumps(sorted(ids), indent=0), encoding="utf-8")

    if args.raw_zarr:
        print(f"RAW  <- {args.raw_zarr}")
        t0 = time.perf_counter()
        nw_r, _, _, nz_r = carve_one(args.raw_zarr, args.out / "cubes_RAW", bbox,
                               skip_empty=False, threshold_binary=False,
                               anon=anon)
        print(f"RAW : {nw_r} cubes written [{time.perf_counter() - t0:.0f}s]")
        if nw_r > 0 and nz_r == nw_r:
            print("ERROR: every RAW cube is all zero; the source store may "
                  "be incomplete and returning its fill value", file=sys.stderr)
            return 1
        if nz_r > 0:
            print(f"WARNING: {nz_r}/{nw_r} RAW cubes are all zero",
                  file=sys.stderr)

    z0, z1, y0, y1, x0, x1 = bbox
    manifest = {
        "chunk_size": CUBE,
        "bbox_l0_zyx": bbox,
        "n_chunks": [(z1 - z0) // CUBE, (y1 - y0) // CUBE, (x1 - x0) // CUBE],
        "sources": {"pred": args.pred_zarr, "raw": args.raw_zarr},
        "umbilicus_yx": args.umbilicus,
        "created_by": "carve_grid_tifs.py",
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=1),
                                            encoding="utf-8")
    print(f"manifest -> {args.out / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
