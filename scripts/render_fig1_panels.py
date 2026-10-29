"""Panels for Figure 1, the end-to-end pipeline strip.

Every panel comes from the same PHerc0139 4x5x5 block, so the figure reads as
one region carried through the whole system rather than five unrelated runs.
All five share the dark BACKGROUND below, which is what makes the strip read
as a single object instead of five loose thumbnails.

Panels 1-3 come from other tools; the exact invocations are:

  1  pipe_ct_pred.png -- raw CT with the surface prediction in red, exported
     from the grid the block was meshed from.  scrollslice insists on a mesh,
     so pass a single triangle placed inside the volume and let --full drive
     the crop; it contributes no visible contour:
       scrollslice.exe --raw PHerc0139-4x5x5/raw.zarr \
                       --pred PHerc0139-4x5x5/pred.zarr \
                       --mesh one_triangle.obj --mesh-axes zyx \
                       --origin 4352,3072,2560 --full \
                       --export z:4416 --png figures/pipe_ct_pred.png

  2  pipe_cube_extract.png -- one 128^3 cube after extraction and cleanup, each
     sheet a separate component:
       render_mesh.py --components --min-comp=400 --bg=0a0c14 \
           output/figstage_dumps_20260731/z04352_y03072_x02560/\
z04352_y03072_x02560_step8_holefill/\
z04352_y03072_x02560_step8_holefill_all.obj \
           raw_cube.png 1000 900 "0.60,-0.55,-0.50" 1.0 3.5
       render_fig1_panels.py crop raw_cube.png figures/pipe_cube_extract.png

  3  pipe_weld_axial.png -- the welded block down the scroll axis:
       render_mesh.py --components --min-comp=200 --bg=0a0c14 \
           output/canonical_best/pherc0139_4x5x5/weld/welded.obj \
           raw_weld.png 900 780 "0.92,0.26,0.29" 1.0 3.0
       render_fig1_panels.py crop raw_weld.png figures/pipe_weld_axial.png

This script renders panels 4 and 5, deliberately over the *same* window of the
atlas domain so the reader can see the charts beside the CT they carry:

  4  pipe_atlas_charts.png  -- flat atlas, filled with its group colours
  5  pipe_readback_page.png -- atlas_rawtex.png over the identical columns

The atlas-to-raster map is exact and was read off atlas_bake.obj's vt values:
raster_x = u + 39124.4062, raster_y = w - 4353.  ATLAS_U_OFFSET below is that
constant; verify_offset() re-derives it so a regenerated run cannot silently
inherit a stale value.

Panel heights in submission.tex assume the aspect ratios these crops produce.
If a window or camera changes, re-derive \rowah and \rowbh there so each row
still fills \textwidth exactly.

Usage: render_fig1_panels.py <canonical_best/pherc0139_4x5x5> <figures-dir>
                             [x0] [width]
       render_fig1_panels.py crop <render.png> <panel.png>
"""
import os
import sys

import numpy as np
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

ATLAS_U_OFFSET = 39124.4062
ATLAS_W_ORIGIN = 4353.0
DEFAULT_X0, DEFAULT_WIDTH = 20800, 1400
BACKGROUND = (10, 12, 20)          # matches the dark atlas review renders
SUPERSAMPLE = 3


def load_flat_atlas(path):
    """after_uv.obj -> (Nx2 uv, Mx3 faces, Nx3 uint8 vertex colours)."""
    vs, fs = [], []
    with open(path, 'rb') as fh:
        for line in fh:
            if line.startswith(b'v '):
                vs.append(line[2:].split())
            elif line.startswith(b'f '):
                fs.append(line[2:].replace(b'/', b' ').split())
    varr = np.array(vs, dtype=np.float64)
    if varr.shape[1] < 6:
        raise SystemExit('after_uv.obj has no vertex colours')
    uv = varr[:, :2]
    col = np.clip(varr[:, 3:6], 0, 1)
    step = len(fs[0]) // 3
    faces = np.array([[f[0], f[step], f[2 * step]] for f in fs],
                     dtype=np.int64) - 1
    return uv, faces, (col * 255.0 + 0.5).astype(np.uint8)


def verify_offset(bake_path, uv, tol=1e-2):
    """Re-derive raster_x - u from atlas_bake.obj's vt block and check it."""
    vts = []
    with open(bake_path, 'rb') as fh:
        for line in fh:
            if line.startswith(b'vt '):
                vts.append(line[3:].split()[:2])
                if len(vts) >= 5000:
                    break
    vt = np.array(vts, dtype=np.float64)
    n = min(len(vt), len(uv))
    dx = vt[:n, 0] - uv[:n, 0]
    dy = vt[:n, 1] - uv[:n, 1]
    if dx.std() > tol or dy.std() > tol:
        raise SystemExit('atlas uv -> raster map is not a pure translation')
    if abs(dx.mean() - ATLAS_U_OFFSET) > tol or \
       abs(dy.mean() + ATLAS_W_ORIGIN) > tol:
        raise SystemExit(
            f'stale constants: this run maps u+{dx.mean():.4f}, '
            f'w{dy.mean():.4f}; update ATLAS_U_OFFSET / ATLAS_W_ORIGIN')
    return float(dx.mean()), float(dy.mean())


def rasterize(uv, faces, colors, x0, width, height, ss=SUPERSAMPLE):
    """Flat-fill the atlas triangles into a window of the raster domain.

    Barycentric fill per triangle over its integer bounding box; the charts do
    not overlap in the atlas, so draw order does not matter and no depth test
    is needed.  Gaps between charts stay background, which is what makes the
    quad-cell structure visible.
    """
    W, H = width * ss, height * ss
    img = np.empty((H, W, 3), np.uint8)
    img[:] = BACKGROUND

    px = (uv[:, 0] + ATLAS_U_OFFSET - x0) * ss
    py = (uv[:, 1] - ATLAS_W_ORIGIN) * ss
    tx, ty = px[faces], py[faces]

    keep = ((tx.max(1) >= 0) & (tx.min(1) < W) &
            (ty.max(1) >= 0) & (ty.min(1) < H))
    tx, ty, tf = tx[keep], ty[keep], faces[keep]
    print(f'  {len(tf)} triangles touch the window')

    lox = np.clip(np.floor(tx.min(1)).astype(int), 0, W - 1)
    hix = np.clip(np.ceil(tx.max(1)).astype(int), 0, W - 1)
    loy = np.clip(np.floor(ty.min(1)).astype(int), 0, H - 1)
    hiy = np.clip(np.ceil(ty.max(1)).astype(int), 0, H - 1)

    x0f, y0f = tx[:, 0], ty[:, 0]
    e1x, e1y = tx[:, 1] - x0f, ty[:, 1] - y0f
    e2x, e2y = tx[:, 2] - x0f, ty[:, 2] - y0f
    det = e1x * e2y - e2x * e1y
    tcol = colors[tf].mean(1)

    for i in range(len(tf)):
        if det[i] == 0.0:
            continue
        ys = np.arange(loy[i], hiy[i] + 1)
        xs = np.arange(lox[i], hix[i] + 1)
        if not len(ys) or not len(xs):
            continue
        rx = xs[None, :] - x0f[i]
        ry = ys[:, None] - y0f[i]
        a = (rx * e2y[i] - ry * e2x[i]) / det[i]
        b = (ry * e1x[i] - rx * e1y[i]) / det[i]
        inside = (a >= 0) & (b >= 0) & (a + b <= 1)
        if inside.any():
            img[loy[i]:hiy[i] + 1, lox[i]:hix[i] + 1][inside] = tcol[i]

    out = Image.fromarray(img)
    if ss > 1:
        out = out.resize((width, height), Image.LANCZOS)
    return out


def tight_crop(src, dst, frac=0.05, pad=6, tol=24):
    """Trim a render down to the rows and columns that actually carry ink.

    A plain bounding box would keep the frame wide open for a handful of stray
    outer fragments, so drop any row or column below `frac` coverage; the
    panels then sit at a comparable scale to each other in the strip.
    """
    a = np.asarray(Image.open(src).convert('RGB')).astype(np.int32)
    ink = np.abs(a - np.asarray(BACKGROUND, np.int32)).sum(2) > tol
    ys = np.where(ink.mean(1) >= frac)[0]
    xs = np.where(ink.mean(0) >= frac)[0]
    if not len(ys) or not len(xs):
        raise SystemExit(f'{src}: nothing above the ink threshold')
    y0, y1 = max(0, ys[0] - pad), min(a.shape[0], ys[-1] + 1 + pad)
    x0, x1 = max(0, xs[0] - pad), min(a.shape[1], xs[-1] + 1 + pad)
    out = Image.fromarray(a[y0:y1, x0:x1].astype(np.uint8))
    out.save(dst)
    print(f'{os.path.basename(dst)}  {a.shape[1]}x{a.shape[0]} -> '
          f'{out.size[0]}x{out.size[1]}  aspect {out.size[0] / out.size[1]:.4f}')


def main():
    if len(sys.argv) >= 4 and sys.argv[1] == 'crop':
        tight_crop(sys.argv[2], sys.argv[3])
        return
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    root, figdir = sys.argv[1], sys.argv[2]
    x0 = int(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_X0
    width = int(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_WIDTH

    uv_path = os.path.join(root, 'atlas', 'after_uv.obj')
    bake_path = os.path.join(root, 'atlas', 'atlas_bake.obj')
    tex_path = os.path.join(root, 'texture', 'atlas_rawtex.png')

    print(f'atlas window: raster x [{x0}, {x0 + width})')
    uv, faces, colors = load_flat_atlas(uv_path)
    print(f'  {len(uv)} vertices, {len(faces)} faces, '
          f'{len(np.unique(colors, axis=0))} group colours')
    dx, dy = verify_offset(bake_path, uv)
    print(f'  uv -> raster verified: x = u + {dx:.4f}, y = w {dy:.4f}')

    raw = np.asarray(Image.open(tex_path).convert('L'))
    height = raw.shape[0]
    if x0 + width > raw.shape[1]:
        raise SystemExit(f'window exceeds the {raw.shape[1]}px raster')

    readback = raw[:, x0:x0 + width]
    cov = float((readback > 0).mean())
    out_e = os.path.join(figdir, 'pipe_readback_page.png')
    Image.fromarray(readback).save(out_e)
    print(f'  (e) {out_e}  {width}x{height}  coverage {cov:.1%}')

    out_d = os.path.join(figdir, 'pipe_atlas_charts.png')
    rasterize(uv, faces, colors, x0, width, height).save(out_d)
    print(f'  (d) {out_d}  {width}x{height}')


if __name__ == '__main__':
    main()
