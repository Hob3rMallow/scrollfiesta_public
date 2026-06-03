"""Auto-pick regions of interest by scanning a coarse resolution level.

Reads a downsampled level (default 5, ~32x) in z-slabs, aggregates surface
(nonzero) counts over candidate ROI footprints, and returns the top-N
non-overlapping, spatially-spread densest blocks as **level-0** voxel bboxes.
Never holds the whole coarse volume — one z-slab at a time.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .zarr_source import ZarrVolume, open_volume, read_block


@dataclass(frozen=True)
class Roi:
    bbox_l0: tuple[int, int, int, int, int, int]  # z0,z1,y0,y1,x0,x1 (level-0)
    surface_voxels: int       # nonzero coarse voxels in the block
    density: float            # fraction of the block that is surface
    scan_cell: tuple[int, int, int]


def _footprint_scan_voxels(roi_cubes, cube, target_level, scan_level) -> tuple[int, int, int]:
    s = (2 ** target_level) * cube
    f = 2 ** scan_level
    return tuple(max(1, (rc * s) // f) for rc in roi_cubes)  # type: ignore[return-value]


def auto_pick_rois(zarr_uri: str, n: int = 3, *, scan_level: int = 5,
                   roi_cubes: tuple[int, int, int] = (2, 2, 2),
                   target_level: int = 0, cube: int = 128,
                   min_separation_cells: int = 2, density_min: float = 0.02,
                   density_max: float = 0.30, anon: bool | None = None,
                   storage_options: dict | None = None) -> list[Roi]:
    """Return up to ``n`` ROIs (level-0 bboxes) with sheet-like surface density.

    We rank by surface count but only accept cells whose density falls in
    [density_min, density_max]: below that is mostly-empty space, above it is a
    dense/noisy blob (many overlapping wraps) that stalls the mesher. Thin
    scroll sheets occupy only a small fraction of a volume, so a moderate band
    selects clean, meshable regions rather than the absolute densest one."""
    coarse = open_volume(zarr_uri, scan_level, anon=anon, storage_options=storage_options)
    fz, fy, fx = _footprint_scan_voxels(roi_cubes, cube, target_level, scan_level)
    Dc, Hc, Wc = coarse.shape
    nz, ny, nx = Dc // fz, Hc // fy, Wc // fx
    if nz == 0 or ny == 0 or nx == 0:
        return []

    dens = np.zeros((nz, ny, nx), dtype=np.int64)
    yext, xext = ny * fy, nx * fx
    # Read chunk-aligned z-blocks (a multiple of fz) so we don't re-fetch the
    # same S3 chunks for every density row.
    cz = int(coarse.chunks[0]) if coarse.chunks else fz
    read_rows = max(fz, (cz // fz) * fz)
    k = 0
    while k < nz:
        nb = min(read_rows // fz, nz - k)        # density cells along z this block
        if nb <= 0:
            break
        z0 = k * fz
        block = read_block(coarse, z0, z0 + nb * fz, 0, yext, 0, xext)
        if block.size:
            dens[k:k + nb] = (block != 0).reshape(
                nb, fz, ny, fy, nx, fx).sum(axis=(1, 3, 5))
        k += nb

    block_voxels = fz * fy * fx
    s0 = 2 ** scan_level                       # level-0 voxels per coarse voxel
    roi_l0 = tuple(rc * cube * (2 ** target_level) for rc in roi_cubes)

    order = np.argsort(dens, axis=None)[::-1]
    picked: list[Roi] = []
    chosen_cells: list[tuple[int, int, int]] = []
    for flat in order:
        cnt = int(dens.flat[flat])
        if cnt == 0:
            break
        dval = cnt / block_voxels
        if dval < density_min or dval > density_max:
            continue                       # skip empty space and dense blobs
        cz, cy, cx = (int(v) for v in np.unravel_index(flat, dens.shape))
        if any(max(abs(cz - pz), abs(cy - py), abs(cx - px)) < min_separation_cells
               for pz, py, px in chosen_cells):
            continue
        z0 = cz * fz * s0
        y0 = cy * fy * s0
        x0 = cx * fx * s0
        bbox = (z0, z0 + roi_l0[0], y0, y0 + roi_l0[1], x0, x0 + roi_l0[2])
        picked.append(Roi(bbox_l0=bbox, surface_voxels=cnt,
                          density=cnt / block_voxels, scan_cell=(cz, cy, cx)))
        chosen_cells.append((cz, cy, cx))
        if len(picked) >= n:
            break
    return picked
