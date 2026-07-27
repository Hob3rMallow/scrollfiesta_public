"""Benchmark chunked padded-cube reads from the remote OME-Zarr.

Proves the streaming property: wall-time + throughput per 154^3 padded cube, with
flat peak memory regardless of how many cubes are read (we never hold the whole
volume). Run manually (hits S3):

    uv run python bench/bench_zarr_read.py [N_cubes] [level]
"""

import resource
import sys
import time

from scrollunwrap.roi_finder import auto_pick_rois
from scrollunwrap.zarr_source import open_volume, read_padded_cube

ZARR = ("s3://vesuvius-challenge-open-data/PHerc0139/representations/predictions/"
        "surfaces/20250728140407-surface-20260413222639-surface-m7-L0-th0.2.zarr")


def _rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    level = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    rois = auto_pick_rois(ZARR, n=1, scan_level=5, roi_cubes=(4, 4, 4),
                          target_level=level)
    if not rois:
        print("no ROI found")
        return
    z0, _, y0, _, x0, _ = rois[0].bbox_l0
    vol = open_volume(ZARR, level)
    s = vol.scale
    oz, oy, ox = z0 // s, y0 // s, x0 // s          # level-L origin
    bytes_per = (154 ** 3)
    print(f"level={level} scale={s} cube_origin(lvl)={(oz, oy, ox)}  reading {n} cubes")
    t0 = time.perf_counter()
    nonzero = 0
    for i in range(n):
        buf = read_padded_cube(vol, oz, oy, ox + i * 128, size=128, halo=13)
        nonzero += int((buf != 0).sum())
    dt = time.perf_counter() - t0
    mb = n * bytes_per / 1e6
    print(f"  {n} cubes in {dt:.2f}s  = {dt/n:.2f}s/cube  "
          f"{mb/dt:.1f} MB/s (uncompressed)  peakRSS={_rss_mb():.0f} MB  "
          f"surface_voxels={nonzero}")


if __name__ == "__main__":
    main()
