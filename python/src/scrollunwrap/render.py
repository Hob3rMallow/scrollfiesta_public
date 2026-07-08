"""Offscreen PNG renders for vision validation (matplotlib Agg, no GPU)."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def _subsample(n: int, cap: int) -> np.ndarray:
    if n <= cap:
        return np.arange(n)
    return np.linspace(0, n - 1, cap).astype(np.int64)


def render_mesh_views(mesh, out_png, *, point_cap: int = 60000) -> Path:
    """Three orthographic vertex-scatter views (ZY, ZX, YX), coloured by the
    third axis — a fast, robust sense of the sheet's shape."""
    V = np.asarray(mesh.vertices, dtype=np.float64)
    idx = _subsample(len(V), point_cap)
    P = V[idx]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    combos = [((0, 1), 2, "z vs y (col by x)"), ((0, 2), 1, "z vs x (col by y)"),
              ((1, 2), 0, "y vs x (col by z)")]
    for ax, ((i, j), c, title) in zip(axes, combos):
        sc = ax.scatter(P[:, i], P[:, j], c=P[:, c], s=1, cmap="viridis")
        ax.set_title(title)
        ax.set_aspect("equal")
        fig.colorbar(sc, ax=ax, shrink=0.7)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
    return Path(out_png)


def render_uv_layout(flat_obj, out_png, *, face_cap: int = 120000) -> Path:
    """UV triangulation; flipped (orientation-reversed) triangles tinted red."""
    from .geometry_report import load_obj_uv, _signed_areas
    _, _, face_uv = load_obj_uv(flat_obj)
    fig, ax = plt.subplots(figsize=(8, 8))
    if len(face_uv):
        s = _signed_areas(face_uv)
        majority_pos = (s > 0).sum() >= (s < 0).sum()
        flipped = (s < 0) if majority_pos else (s > 0)
        sel = _subsample(len(face_uv), face_cap)
        from matplotlib.collections import PolyCollection
        tris = face_uv[sel]
        cols = np.where(flipped[sel][:, None], np.array([1.0, 0.0, 0.0, 0.6]),
                        np.array([0.2, 0.4, 0.8, 0.25]))
        pc = PolyCollection(tris, facecolors=cols, edgecolors=(0, 0, 0, 0.15),
                            linewidths=0.2)
        ax.add_collection(pc)
        ax.autoscale_view()
        ax.set_title(f"UV map  (flipped triangles: {int(flipped.sum())}/{len(s)})")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
    return Path(out_png)


def render_tifxyz_channels(tifxyz_dir, out_png) -> Path:
    """x/y/z channels + validity mask (z>0), as heatmaps. Discontinuities in the
    gradients flag flattening/coordinate problems."""
    import tifffile
    d = Path(tifxyz_dir)
    x = tifffile.imread(d / "x.tif").astype(np.float64)
    y = tifffile.imread(d / "y.tif").astype(np.float64)
    z = tifffile.imread(d / "z.tif").astype(np.float64)
    valid = z > 0
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    for ax, (img, name) in zip(axes[:3], ((x, "x.tif"), (y, "y.tif"), (z, "z.tif"))):
        m = np.where(valid, img, np.nan)
        im = ax.imshow(m, cmap="turbo", origin="upper")
        ax.set_title(name)
        fig.colorbar(im, ax=ax, shrink=0.7)
    axes[3].imshow(valid, cmap="gray", origin="upper")
    axes[3].set_title(f"valid (z>0): {int(valid.sum())}/{valid.size}")
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
    return Path(out_png)


def collect_flatboi_heatmaps(flat_obj) -> list[Path]:
    """flatboi writes <stem>_heat_*.png next to its output; return any found."""
    p = Path(flat_obj)
    return sorted(p.parent.glob(f"{p.stem}_heat_*.png")) + \
        sorted(p.parent.glob(f"{p.stem}*heat*.png"))
