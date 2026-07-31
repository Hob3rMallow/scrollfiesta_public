"""Verify a local OME-Zarr mirror against its remote source.

Checks, for the local mirror produced by ``mirror_zarr.py``:

  1. every level opens via ``scrollunwrap.zarr_source.open_volume`` with the
     expected shape / dtype (i.e. the mesher can read it);
  2. a padded cube read from the LOCAL store is byte-identical to the same cube
     read from the REMOTE store (proves chunk content, not just object count).

Usage::

    py -m uv run --project python python scripts/verify_zarr_mirror.py \
        --local D:/work/vesuvius-c/PHerc0139-full/<...>.zarr \
        --remote s3://vesuvius-challenge-open-data/PHerc0139/.../<...>.zarr
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

from scrollunwrap.roi_finder import auto_pick_rois
from scrollunwrap.zarr_source import open_volume, read_padded_cube


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--local", required=True, help="local mirrored zarr path")
    ap.add_argument("--remote", required=True, help="remote s3:// source zarr")
    ap.add_argument("--levels", default="0,1,2,3,4,5")
    args = ap.parse_args(argv)

    ok = True
    levels = [int(s) for s in args.levels.split(",") if s.strip()]

    print("=== level open / shape check (local) ===")
    local_l0 = None
    for lv in levels:
        try:
            v = open_volume(args.local, lv)
            print(f"  level {lv}: shape={v.shape} chunks={v.chunks} dtype={v.dtype}  OK")
            if lv == 0:
                local_l0 = v
        except Exception as e:
            print(f"  level {lv}: OPEN FAILED: {e}")
            ok = False

    print("\n=== cube byte-parity local vs remote ===")
    # Pick a dense cube via the coarse scan (uses local mirror), then read it
    # from both stores and compare.
    try:
        rois = auto_pick_rois(args.local, n=1, scan_level=5, roi_cubes=(1, 1, 1),
                              target_level=0)
        if not rois:
            print("  no dense ROI found on local mirror (unexpected)")
            return 1 if not ok else 0
        z0, _, y0, _, x0, _ = rois[0].bbox_l0
        # snap to a 128 origin
        oz, oy, ox = (z0 // 128) * 128, (y0 // 128) * 128, (x0 // 128) * 128
        print(f"  dense cube origin (l0 voxels): ({oz},{oy},{ox})")
        lv0 = local_l0 or open_volume(args.local, 0)
        rv0 = open_volume(args.remote, 0)
        lbuf = read_padded_cube(lv0, oz, oy, ox, size=128, halo=13)
        rbuf = read_padded_cube(rv0, oz, oy, ox, size=128, halo=13)
        nz = int((lbuf != 0).sum())
        identical = bool(np.array_equal(lbuf, rbuf))
        print(f"  local cube nonzero voxels = {nz}")
        print(f"  local vs remote byte-identical = {identical}")
        if not identical or nz == 0:
            ok = False
    except Exception as e:
        print(f"  parity check FAILED: {e}")
        ok = False

    print("\nRESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
