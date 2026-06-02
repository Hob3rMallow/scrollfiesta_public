"""Open a remote OME-Zarr surface volume over s3:// and read it in chunks.

Mirrors the proven open pattern from dinovol
(``dinovol_2/dataset/ssl_zarr_dataset.py``): anonymous s3 access via
``zarr.open(uri, path=str(level), storage_options={"anon": True})``, region
reads by ZYX slicing. We never materialize the whole array — the largest read
is a single padded cube (``(128+2*halo)**3`` bytes).
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import zarr


@dataclass(frozen=True)
class ZarrVolume:
    """A single OME-Zarr resolution level (ZYX), with cached shape/chunks."""

    array: "zarr.Array"
    uri: str
    level: int
    scale: int  # voxels of this level per... i.e. level-0 == level * 2**level
    shape: tuple[int, int, int]
    chunks: tuple[int, int, int]
    dtype: np.dtype


def open_volume(zarr_uri: str, level: int = 0, *, anon: bool = True) -> ZarrVolume:
    """Open one resolution ``level`` of an OME-Zarr store.

    Works for both ``s3://`` (anonymous by default) and local filesystem URIs,
    and for both zarr v2-format and v3-format stores.
    """
    uri = str(zarr_uri).rstrip("/")
    if uri.startswith("s3://"):
        arr = zarr.open(uri, path=str(level), mode="r",
                        storage_options={"anon": bool(anon)})
    else:
        arr = zarr.open(uri, path=str(level), mode="r")
    shape = tuple(int(s) for s in arr.shape)
    if len(shape) != 3:
        raise ValueError(f"expected a 3D (z,y,x) array, got shape {shape}")
    return ZarrVolume(
        array=arr,
        uri=uri,
        level=int(level),
        scale=2 ** int(level),
        shape=shape,  # type: ignore[arg-type]
        chunks=tuple(int(c) for c in arr.chunks),  # type: ignore[arg-type]
        dtype=np.dtype(arr.dtype),
    )


def read_block(vol: ZarrVolume, z0: int, z1: int, y0: int, y1: int,
               x0: int, x1: int) -> np.ndarray:
    """Read ``array[z0:z1, y0:y1, x0:x1]`` (clipped to bounds) as a numpy array.

    NEVER reads the whole array. Returns an empty array if the (clipped) box is
    degenerate.
    """
    D, H, W = vol.shape
    z0, z1 = max(0, int(z0)), min(D, int(z1))
    y0, y1 = max(0, int(y0)), min(H, int(y1))
    x0, x1 = max(0, int(x0)), min(W, int(x1))
    if z0 >= z1 or y0 >= y1 or x0 >= x1:
        return np.empty((0, 0, 0), dtype=vol.dtype)
    return np.asarray(vol.array[z0:z1, y0:y1, x0:x1])


def read_padded_cube(vol: ZarrVolume, oz: int, oy: int, ox: int,
                     size: int = 128, halo: int = 13) -> np.ndarray:
    """Read one cube at voxel origin ``(oz,oy,ox)`` padded by ``halo`` on every
    side, into a C-contiguous ``(p,p,p)`` uint8 buffer (``p = size + 2*halo``).

    Reproduces ``HaloLoader_load``: buffer index ``(0,0,0)`` corresponds to world
    voxel ``(oz-halo, oy-halo, ox-halo)``; out-of-volume halo stays zero. This is
    a single chunk-streamed slice — it replaces both the per-cube TIFF loader and
    the 26-neighbour halo loader.
    """
    p = int(size) + 2 * int(halo)
    buf = np.zeros((p, p, p), dtype=np.uint8)
    D, H, W = vol.shape
    # Desired window (may extend outside the array on the low/high side).
    z0d, y0d, x0d = oz - halo, oy - halo, ox - halo
    z1d, y1d, x1d = oz + size + halo, oy + size + halo, ox + size + halo
    z0, z1 = max(0, z0d), min(D, z1d)
    y0, y1 = max(0, y0d), min(H, y1d)
    x0, x1 = max(0, x0d), min(W, x1d)
    if z0 < z1 and y0 < y1 and x0 < x1:
        data = np.asarray(vol.array[z0:z1, y0:y1, x0:x1])
        buf[z0 - z0d:z1 - z0d, y0 - y0d:y1 - y0d, x0 - x0d:x1 - x0d] = data
    return buf


_THRESH_RE = re.compile(r"^\s*(>=|<=|==|!=|>|<)\s*([0-9]+)\s*$")


def apply_threshold(arr: np.ndarray, spec: str | None) -> np.ndarray:
    """Binarize a uint8 mask to {0,255} for the mesher (nonzero = surface).

    ``spec`` None -> pass-through (any nonzero -> 255). Otherwise a comparison
    like ``">=1"``, ``">128"``, ``"==255"``.
    """
    if spec is None:
        mask = arr != 0
    else:
        m = _THRESH_RE.match(spec)
        if not m:
            raise ValueError(f"bad --threshold {spec!r}; use e.g. '>=1', '>128', '==255'")
        op, val = m.group(1), int(m.group(2))
        cmp = {
            ">=": arr >= val, "<=": arr <= val, "==": arr == val,
            "!=": arr != val, ">": arr > val, "<": arr < val,
        }[op]
        mask = cmp
    return (mask.astype(np.uint8)) * np.uint8(255)
