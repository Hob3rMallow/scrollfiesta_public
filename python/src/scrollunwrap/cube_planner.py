"""Plan the 128-cube grid covering a region of interest.

The ROI is given in **level-0** voxel coordinates. Meshing happens at level
``L`` (downsample ``2**L``), so cube origins are emitted in **level-L** voxel
coordinates, snapped to multiples of ``cube`` (128) so that each cube id encodes
its world origin exactly (``halo_loader.c`` / ``dump_obj.c`` parse the filename
number as the voxel origin).
"""

from __future__ import annotations

import math
from dataclasses import dataclass


def cube_id_str(oz: int, oy: int, ox: int) -> str:
    """scrollfiesta cube id, e.g. ``z00128_y00256_x00384`` (matches the
    ``z%05d_y%05d_x%05d`` format ``halo_loader.c`` expects)."""
    return f"z{oz:05d}_y{oy:05d}_x{ox:05d}"


@dataclass(frozen=True)
class Cube:
    oz: int
    oy: int
    ox: int

    @property
    def cube_id(self) -> str:
        return cube_id_str(self.oz, self.oy, self.ox)


def snap_bbox_to_level_grid(bbox_l0: tuple[int, int, int, int, int, int],
                            level: int, cube: int = 128
                            ) -> tuple[int, int, int, int, int, int]:
    """Convert a level-0 bbox ``(z0,z1,y0,y1,x0,x1)`` to level-``L`` voxels and
    snap lo down / hi up to multiples of ``cube``. Returns the snapped level-L
    bbox (origins are the lo corners)."""
    z0, z1, y0, y1, x0, x1 = bbox_l0
    s = 2 ** int(level)
    out = []
    for lo0, hi0 in ((z0, z1), (y0, y1), (x0, x1)):
        lo = (lo0 // s)
        hi = math.ceil(hi0 / s)
        lo_snap = (lo // cube) * cube
        hi_snap = math.ceil(hi / cube) * cube
        if hi_snap <= lo_snap:  # empty/degenerate -> at least one cube
            hi_snap = lo_snap + cube
        out.append((max(0, lo_snap), hi_snap))
    (zlo, zhi), (ylo, yhi), (xlo, xhi) = out
    return (zlo, zhi, ylo, yhi, xlo, xhi)


def plan_cubes(bbox_l0: tuple[int, int, int, int, int, int], level: int = 0,
               cube: int = 128) -> list[Cube]:
    """Enumerate owned-cube origins (level-L voxels, multiples of ``cube``)
    covering the snapped ROI. No halo ring is materialized — the halo for each
    cube is read inline from the zarr at mesh time."""
    zlo, zhi, ylo, yhi, xlo, xhi = snap_bbox_to_level_grid(bbox_l0, level, cube)
    cubes: list[Cube] = []
    for oz in range(zlo, zhi, cube):
        for oy in range(ylo, yhi, cube):
            for ox in range(xlo, xhi, cube):
                cubes.append(Cube(oz, oy, ox))
    return cubes
