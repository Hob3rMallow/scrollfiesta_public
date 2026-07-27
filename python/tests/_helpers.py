"""Shared test helpers (synthetic local zarr + meshes). No network, no binaries."""

from __future__ import annotations

import numpy as np
import trimesh
import zarr


def make_local_zarr(path, shape=(40, 40, 40), chunks=(16, 16, 16),
                    blob=(slice(8, 24), slice(8, 24), slice(8, 24)), value=200,
                    level_name="0"):
    """Write a zarr v2 store with one array at ``level_name`` and a dense blob."""
    g = zarr.open_group(store=str(path), mode="w", zarr_format=2)
    a = g.create_array(name=level_name, shape=shape, chunks=chunks, dtype="uint8")
    vol = np.zeros(shape, np.uint8)
    if blob is not None:
        vol[blob] = value
    a[:] = vol
    return str(path), vol


def grid_disk(n=6, hole=False, translate=(0.0, 0.0, 0.0)):
    """A flat triangulated n×n grid in the z=0 plane: a topological disk.
    With ``hole=True`` one central quad is removed (an interior hole)."""
    xs, ys = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    V = np.stack([xs.ravel(), ys.ravel(), np.zeros(n * n)], axis=1).astype(float)
    V += np.asarray(translate, dtype=float)
    F = []

    def vid(i, j):
        return i * n + j

    for i in range(n - 1):
        for j in range(n - 1):
            if hole and i == n // 2 - 1 and j == n // 2 - 1:
                continue
            F.append([vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)])
            F.append([vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)])
    return trimesh.Trimesh(np.asarray(V), np.asarray(F, dtype=np.int64), process=False)
