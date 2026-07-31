"""Render the four Fig-3 (fig:cube_stages) panels in the MeshLab-screenshot style.

July 2026 figure source: cube z04352_y03456_x02944 comp 005 (cloud step0_mls_003,
cvt step10_qem_005) from a no-halo cube_mesh run on PHerc0139-4x5x5 cubes_PRED;
see render_fig3_drive.py for the auto-mapping driver.

Colors sampled from the existing overview0*_crop.png:
  mesh fill  (255,139,205)  wireframe (64,64,64)  points (82,76,26)
  hole-fill patches (27,123,255)  background white

All four panels share one orthographic camera derived from the hole-filled
mesh's PCA frame, so the stages align like the original figure.

Usage:
  python render_panels.py <points.obj> <bpa.obj> <fill.obj> <cvt.obj> <outdir>
         [--rot DEG] [--mirror] [--size N]
"""
import sys, os
import numpy as np
from PIL import Image, ImageDraw

PINK = np.array([255, 139, 205], float)
EDGE = (64, 64, 64)
OLIVE = (82, 76, 26)
BLUE = np.array([27, 123, 255], float)
SS = 1  # no supersampling: originals are crisp aliased MeshLab screenshots

def load_obj(path):
    vs, fs = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                p = line.split()
                vs.append((float(p[1]), float(p[2]), float(p[3])))
            elif line.startswith("f "):
                idx = [int(t.split("/")[0]) - 1 for t in line.split()[1:]]
                for k in range(1, len(idx) - 1):
                    fs.append((idx[0], idx[k], idx[k + 1]))
    return np.asarray(vs, np.float64), np.asarray(fs, np.int32)

def face_keys(vs, fs, tol=1e-3):
    q = np.round(vs / tol).astype(np.int64)
    keys = set()
    for f in fs:
        keys.add(frozenset((tuple(q[f[0]]), tuple(q[f[1]]), tuple(q[f[2]]))))
    return keys

def frame_from(vs, rot_deg=0.0, mirror=False):
    c = vs.mean(0)
    _, _, Vt = np.linalg.svd(vs - c, full_matrices=False)
    up, right, normal = Vt[0], Vt[1], Vt[2]
    # right-handed: right x up = toward viewer
    if np.dot(np.cross(right, up), normal) < 0:
        normal = -normal
    if mirror:
        right = -right
        normal = -normal
    # auto-upright: rotate in-plane so the world axis with the largest
    # in-plane footprint maps to screen vertical
    best, ba = 0.0, 0.0
    for ax in np.eye(3):
        ip = np.array([np.dot(ax, right), np.dot(ax, up)])
        l = np.hypot(*ip)
        if l > best:
            best, ba = l, np.arctan2(ip[0], ip[1])  # angle from +up
    a = ba + np.deg2rad(rot_deg)
    ca, sa = np.cos(a), np.sin(a)
    r2 = ca * right - sa * up
    u2 = sa * right + ca * up
    return c, r2, u2, normal

def project(vs, frame):
    c, r, u, n = frame
    d = vs - c
    return np.stack([d @ r, d @ u, d @ n], 1)

class Canvas:
    def __init__(self, bbox, height=1200, margin=0.04):
        (x0, y0), (x1, y1) = bbox
        w, h = x1 - x0, y1 - y0
        pad = margin * max(w, h)
        x0 -= pad; y0 -= pad; x1 += pad; y1 += pad
        self.scale = (height * SS) / (y1 - y0)
        self.W = int(round((x1 - x0) * self.scale))
        self.H = height * SS
        self.x0, self.y1 = x0, y1

    def to_px(self, p2):
        x = (p2[:, 0] - self.x0) * self.scale
        y = (self.y1 - p2[:, 1]) * self.scale
        return np.stack([x, y], 1)

def render_mesh(p3, fs, canvas, out, blue_keys=None, keys=None, shade=True):
    img = Image.new("RGB", (canvas.W, canvas.H), (255, 255, 255))
    dr = ImageDraw.Draw(img)
    px = canvas.to_px(p3[:, :2])
    z = p3[:, 2]
    # face normals for shading
    e1 = p3[fs[:, 1]] - p3[fs[:, 0]]
    e2 = p3[fs[:, 2]] - p3[fs[:, 0]]
    fn = np.cross(e1, e2)
    nl = np.linalg.norm(fn, axis=1); nl[nl == 0] = 1
    ndv = np.abs(fn[:, 2] / nl)          # |n . view|
    lum = 0.88 + 0.16 * ndv if shade else np.ones(len(fs))
    order = np.argsort(z[fs].mean(1))     # back to front
    for i in order:
        a, b, c = fs[i]
        col = PINK * lum[i]
        if blue_keys is not None:
            k = keys[i]
            if k in blue_keys:
                col = BLUE * (0.85 + 0.15 * ndv[i])
        col = tuple(int(min(255, v)) for v in col)
        poly = [tuple(px[a]), tuple(px[b]), tuple(px[c])]
        dr.polygon(poly, fill=col)
        dr.line(poly + [poly[0]], fill=EDGE, width=SS)
    if SS > 1:
        img = img.resize((canvas.W // SS, canvas.H // SS), Image.LANCZOS)
    img.save(out)
    return img

def render_points(p3, canvas, out, r=None):
    img = Image.new("RGB", (canvas.W, canvas.H), (255, 255, 255))
    dr = ImageDraw.Draw(img)
    px = canvas.to_px(p3[:, :2])
    if r is None:
        r = max(1, int(round(0.0012 * canvas.H)))
    for x, y in px:
        dr.rectangle([x - r, y - r, x + r, y + r], fill=OLIVE)
    if SS > 1:
        img = img.resize((canvas.W // SS, canvas.H // SS), Image.LANCZOS)
    img.save(out)
    return img

def crop_white(img, pad=10):
    a = np.array(img.convert("L"))
    ys, xs = np.where(a < 250)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, a.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, a.shape[1])
    return img.crop((x0, y0, x1, y1))

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    pts_p, bpa_p, fill_p, cvt_p, outdir = args
    rot = 0.0; mirror = False; height = 1200
    for i, a in enumerate(sys.argv):
        if a == "--rot":
            rot = float(sys.argv[i + 1])
        elif a == "--mirror":
            mirror = True
        elif a == "--size":
            height = int(sys.argv[i + 1])
    os.makedirs(outdir, exist_ok=True)

    pts, _ = load_obj(pts_p)
    bv, bf = load_obj(bpa_p)
    fv, ff = load_obj(fill_p)
    cv, cf = load_obj(cvt_p)

    frame = frame_from(fv, rot, mirror)
    P = {k: project(v, frame) for k, v in
         dict(pts=pts, bpa=bv, fill=fv, cvt=cv).items()}
    allp = np.concatenate([P["pts"], P["bpa"], P["fill"], P["cvt"]])
    bbox = (allp[:, :2].min(0), allp[:, :2].max(0))
    canvas = Canvas(bbox, height)

    b_keys = [frozenset(map(tuple, np.round(fv[f] / 1e-3).astype(np.int64)))
              for f in ff]
    bpa_keys = face_keys(bv, bf)
    blue = {k for k in b_keys if k not in bpa_keys}
    print(f"fill-only faces: {len(blue)}")

    imgs = [
        render_points(P["pts"], canvas, os.path.join(outdir, "overview00.png")),
        render_mesh(P["bpa"], bf, canvas, os.path.join(outdir, "overview01.png")),
        render_mesh(P["fill"], ff, canvas, os.path.join(outdir, "overview02.png"),
                    blue_keys=blue, keys=b_keys),
        render_mesh(P["cvt"], cf, canvas, os.path.join(outdir, "overview03.png")),
    ]
    for i, im in enumerate(imgs):
        crop_white(im).save(os.path.join(outdir, f"overview0{i}_crop.png"))
    print("wrote 8 images ->", outdir)

if __name__ == "__main__":
    main()
