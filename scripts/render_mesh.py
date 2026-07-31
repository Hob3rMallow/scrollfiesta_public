"""Offscreen OBJ renderer: orthographic z-buffered surface-sample splatting.

Loads one OBJ or a directory of *_placed.obj cube meshes, shades by surface
normal with a two-sided key light, and writes a white-background PNG.
numpy + PIL only (no pandas / matplotlib).

Coordinates in the scrollfiesta OBJs are (z, y, x); we treat them abstractly.

Usage: render_mesh.py <obj-or-placed-dir> <out.png> [W H "dz,dy,dx" zcut spp]
  W H        output size (default 1600x1200; rendered 2x supersampled)
  dz,dy,dx   orthographic view direction (default 0.55,0.5,0.67)
  zcut       keep triangles with centroid z below lo+zcut*extent (default 1.0)
  spp        surface samples per projected pixel (default 3.0)

Flags: --components colours each connected component from a fixed palette
(instead of normal shading); --min-comp=N drops components below N faces;
--box=z,y,x,size crops faces to a cube-aligned box before rendering.

The submission figures were produced with:
  render_mesh.py output/repro_10x_preserved_production/atlas_bake.obj \
      submission_update/figures/mesh_10x10x10.png 1500 1300 "0.85,0.38,0.4" 0.55 3.0
  render_mesh.py output/pherc0139_4x21_atlas_proxy/placed \
      submission_update/figures/mesh_4x21x21.png 2200 1500 "0.55,0.5,0.67" 1.0 2.5
  render_mesh.py --components --min-comp=100 --box=4480,3328,2688,128 \
      output/pherc0139_4x5x5/1_mesh/welded.obj \
      submission_update/figures/pipe_cube_sheets.png 1100 980 "0.60,-0.55,-0.50" 1.0 4.0
"""
import sys, glob, math, time
import numpy as np
from PIL import Image


def load_obj(fp):
    with open(fp, 'rb') as fh:
        data = fh.read()
    vlines, flines = [], []
    for line in data.splitlines():
        if line.startswith(b'v '):
            vlines.append(line[2:])
        elif line.startswith(b'f '):
            flines.append(line[2:].replace(b'/', b' '))
    if not vlines or not flines:
        return None, None
    vcol = len(vlines[0].split())          # 3, or 6 with RGB vertex colours
    verts = np.array(b' '.join(vlines).split(), dtype=np.float64)
    verts = verts.reshape(-1, vcol)[:, :3]
    ncol = len(flines[0].split())
    ftok = np.array(b' '.join(flines).split(), dtype=np.int64)
    faces = ftok.reshape(-1, ncol)
    step = ncol // 3
    faces = faces[:, [0, step, 2 * step]] - 1
    return verts.astype(np.float32), faces.astype(np.int64)


PALETTE = np.array([
    (0.86, 0.30, 0.30), (0.36, 0.65, 0.32), (0.31, 0.48, 0.79),
    (0.62, 0.42, 0.72), (0.93, 0.61, 0.25), (0.89, 0.79, 0.30),
    (0.35, 0.72, 0.68), (0.83, 0.45, 0.65), (0.55, 0.36, 0.28),
    (0.55, 0.62, 0.79), (0.72, 0.72, 0.40), (0.44, 0.30, 0.55),
], np.float32)


def face_components(faces, nv):
    """Union-find over shared vertices; returns per-face component ids,
    relabelled densest-first so PALETTE order follows component size."""
    parent = np.arange(nv, dtype=np.int64)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for f in faces:
        a, b, c = find(f[0]), find(f[1]), find(f[2])
        if b != a:
            parent[b] = a
        if c != a:
            parent[find(c)] = a
    roots = np.array([find(f[0]) for f in faces])
    uniq, inv, counts = np.unique(roots, return_inverse=True, return_counts=True)
    order = np.argsort(-counts)
    remap = np.empty(len(uniq), np.int64)
    remap[order] = np.arange(len(uniq))
    return remap[inv]


class Renderer:
    def __init__(self, W, H, ss=2):
        self.W, self.H, self.ss = W * ss, H * ss, ss
        self.zbuf = np.full(self.W * self.H, np.inf, np.float32)
        self.rgb = np.ones((self.W * self.H, 3), np.float32)

    def set_camera(self, cam_dir, up_hint, center, half_w, half_h):
        d = np.asarray(cam_dir, np.float64)
        d = d / np.linalg.norm(d)
        u0 = np.asarray(up_hint, np.float64)
        r = np.cross(u0, d)
        r /= np.linalg.norm(r)
        u = np.cross(d, r)
        self.basis = np.stack([r, u, d]).astype(np.float64)   # rows: right,up,dir
        self.center = np.asarray(center, np.float64)
        self.sx = (self.W / 2) / half_w
        self.sy = (self.H / 2) / half_h

    def splat(self, tris, light=(0.35, 0.45, 0.85), spp=3.0, chunk=400_000,
              face_colors=None):
        """tris: (M,3,3) float32 world triangles; face_colors: optional (M,3)
        flat colours that replace the normal-map base."""
        M = len(tris)
        for lo in range(0, M, chunk):
            t = tris[lo:lo + chunk].astype(np.float64)
            fc = face_colors[lo:lo + chunk] if face_colors is not None else None
            e1 = t[:, 1] - t[:, 0]
            e2 = t[:, 2] - t[:, 0]
            n = np.cross(e1, e2)
            nl = np.linalg.norm(n, axis=1)
            keep = nl > 1e-12
            t, e1, e2, n, nl = t[keep], e1[keep], e2[keep], n[keep], nl[keep]
            if fc is not None:
                fc = fc[keep]
            if not len(t):
                continue
            n /= nl[:, None]
            area = 0.5 * nl
            # projected pixel area ~ area * sx * sy (upper bound; fine)
            ns = np.maximum(1, np.ceil(area * self.sx * self.sy * spp).astype(np.int64))
            ns = np.minimum(ns, 100_000)
            tot = int(ns.sum())
            fid = np.repeat(np.arange(len(t)), ns)
            a = np.random.rand(tot); b = np.random.rand(tot)
            flip = a + b > 1.0
            a[flip] = 1.0 - a[flip]; b[flip] = 1.0 - b[flip]
            P = t[fid, 0] + a[:, None] * e1[fid] + b[:, None] * e2[fid]
            self._shade_splat(P, n[fid], light,
                              fc[fid] if fc is not None else None)

    def _shade_splat(self, P, N, light, C=None):
        Pc = (P - self.center) @ self.basis.T
        x = Pc[:, 0] * self.sx + self.W / 2
        y = self.H / 2 - Pc[:, 1] * self.sy
        z = Pc[:, 2]
        xi = x.astype(np.int64); yi = y.astype(np.int64)
        ok = (xi >= 0) & (xi < self.W) & (yi >= 0) & (yi < self.H)
        if not ok.any():
            return
        xi, yi, z, N = xi[ok], yi[ok], z[ok].astype(np.float32), N[ok]
        if C is not None:
            C = C[ok]
        Nc = N @ self.basis.T
        Nc[Nc[:, 2] < 0] *= -1.0                            # two-sided: face the camera
        L = np.asarray(light, np.float64); L /= np.linalg.norm(L)
        lam = np.abs(Nc @ L)
        if C is not None:                                   # flat component colours
            shade = (0.55 + 0.45 * lam).astype(np.float32)
            col = np.clip(C * shade[:, None], 0, 1)
        else:
            shade = (0.60 + 0.40 * lam).astype(np.float32)
            base = (0.58 + 0.42 * Nc).astype(np.float32)    # camera-space normal -> pastel
            col = np.clip(base * shade[:, None], 0, 1) ** 0.82  # mild gamma lift
        pix = yi * self.W + xi
        order = np.lexsort((z, pix))
        pix, z, col = pix[order], z[order], col[order]
        first = np.ones(len(pix), bool)
        first[1:] = pix[1:] != pix[:-1]
        pix, z, col = pix[first], z[first], col[first]
        upd = z < self.zbuf[pix]
        self.zbuf[pix[upd]] = z[upd]
        self.rgb[pix[upd]] = col[upd]

    def image(self):
        img = self.rgb.reshape(self.H, self.W, 3)
        s = self.ss
        img = img.reshape(self.H // s, s, self.W // s, s, 3).mean(axis=(1, 3))
        return Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8))


def bbox_of_files(files):
    """Fast path: cube-name corners (z####_y####_x####, 128-vox cubes)."""
    import re
    lo = np.full(3, np.inf); hi = np.full(3, -np.inf)
    pat = re.compile(r'z(\d+)_y(\d+)_x(\d+)')
    for fp in files:
        m = pat.search(fp)
        if not m:
            v, f = load_obj(fp)
            if v is None:
                continue
            lo = np.minimum(lo, v.min(0)); hi = np.maximum(hi, v.max(0))
            continue
        c = np.array([float(m.group(1)), float(m.group(2)), float(m.group(3))])
        lo = np.minimum(lo, c); hi = np.maximum(hi, c + 128.0)
    return lo, hi


def fit_camera(rn, corners, cam_dir, up_hint, margin=1.04):
    d = np.asarray(cam_dir, np.float64); d /= np.linalg.norm(d)
    u0 = np.asarray(up_hint, np.float64)
    r = np.cross(u0, d); r /= np.linalg.norm(r)
    u = np.cross(d, r)
    c = corners.mean(0)
    pr = (corners - c) @ np.stack([r, u]).T
    half_w = np.abs(pr[:, 0]).max() * margin
    half_h = np.abs(pr[:, 1]).max() * margin
    return c, half_w, half_h


def corners(lo, hi):
    return np.array([[lo[0], lo[1], lo[2]], [lo[0], lo[1], hi[2]],
                     [lo[0], hi[1], lo[2]], [lo[0], hi[1], hi[2]],
                     [hi[0], lo[1], lo[2]], [hi[0], lo[1], hi[2]],
                     [hi[0], hi[1], lo[2]], [hi[0], hi[1], hi[2]]])


def main():
    np.random.seed(7)
    by_component = '--components' in sys.argv
    if by_component:
        sys.argv.remove('--components')
    min_comp_faces = 0
    box = None
    for arg in list(sys.argv):
        if arg.startswith('--min-comp='):
            min_comp_faces = int(arg.split('=')[1])
            sys.argv.remove(arg)
        elif arg.startswith('--box='):
            box = [float(x) for x in arg.split('=')[1].split(',')]  # z,y,x,size
            sys.argv.remove(arg)
    mode = sys.argv[1]
    out = sys.argv[2]
    W = int(sys.argv[3]) if len(sys.argv) > 3 else 1600
    H = int(sys.argv[4]) if len(sys.argv) > 4 else 1200
    t0 = time.time()

    if mode.endswith('.obj'):
        files = [mode]
    else:
        files = sorted(glob.glob(mode + '/*_placed.obj'))
    print(len(files), 'files')

    # camera params from argv: dir z,y,x + optional z-cut fraction
    cd = [float(v) for v in sys.argv[5].split(',')] if len(sys.argv) > 5 else [0.55, 0.5, 0.67]
    zcut = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
    spp = float(sys.argv[7]) if len(sys.argv) > 7 else 3.0

    lo, hi = bbox_of_files(files) if len(files) > 1 else (None, None)
    if lo is None:
        v, f = load_obj(files[0])
        lo, hi = v.min(0), v.max(0)
    if box is not None:
        lo = np.array(box[:3]); hi = lo + box[3]
    zmax = lo[0] + zcut * (hi[0] - lo[0])
    hi_c = hi.copy(); hi_c[0] = min(hi[0], zmax)
    rn = Renderer(W, H, ss=2)
    c, hw, hh = fit_camera(rn, corners(lo, hi_c), cd, (1.0, 0, 0))
    rn.set_camera(cd, (1.0, 0, 0), c, hw, hh)

    for i, fp in enumerate(files):
        v, f = load_obj(fp)
        if v is None:
            continue
        if box is not None:
            cen = v[f].mean(1)
            b0 = np.array(box[:3]); sz = box[3]
            inside = np.all((cen >= b0) & (cen <= b0 + sz), axis=1)
            f = f[inside]
            print(f'box crop: {len(f)} faces')
        cols = None
        if by_component:
            comp = face_components(f, len(v))
            counts = np.bincount(comp)
            keep_comp = counts >= min_comp_faces
            sel = keep_comp[comp]
            f = f[sel]
            comp = comp[sel]
            print(f'{counts.size} components, {int(keep_comp.sum())} kept '
                  f'(>= {min_comp_faces} faces)')
            cols = PALETTE[comp % len(PALETTE)]
        tris = v[f]                                     # (M,3,3)
        if zcut < 1.0:
            cz = tris[:, :, 0].mean(1)
            keep = cz <= zmax
            tris = tris[keep]
            if cols is not None:
                cols = cols[keep]
        if len(tris):
            rn.splat(tris, spp=spp, face_colors=cols)
        if (i + 1) % 200 == 0:
            print(f'{i+1}/{len(files)} {time.time()-t0:.0f}s', flush=True)

    rn.image().save(out)
    print('wrote', out, f'{time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
