"""Weld z-slab tifxyz strips into one full-height atlas.

Each 1x21x21 strip was exported independently from the SAME global solve
(scroll_unroll --z-range), so same-wrap segments across strips share the u
lattice and tile in v (world z). For each wrap k we stack its (<=4) strip
segments into one full-height tifxyz, and emit a unified atlas.json.

DEMO harness (tifxyz raster tiling). Productionizing as a C tool
(tifxyz_weld + --selftest) is the follow-up for the full-scroll pipeline.
"""
import numpy as np, tifffile, json, os, re, glob, sys

base = r"D:\work\vesuvius-c\output\run_4x21x21"
outdir = os.path.join(base, "welded")
outseg = os.path.join(outdir, "seg")
os.makedirs(outseg, exist_ok=True)
DU = DV = 1.0

# ---- collect strip segments, grouped by wrap k ----
segs = {}   # k -> [(seg_dir, origin_u, origin_v, W, H)]
strip_valid = 0
strip_dirs = sorted(glob.glob(os.path.join(base, "strip_z*", "seg")))
print("strips:", [os.path.basename(os.path.dirname(d)) for d in strip_dirs])
for d in strip_dirs:
    for name in sorted(os.listdir(d)):
        sd = os.path.join(d, name)
        mp = os.path.join(sd, "meta.json")
        if not os.path.isfile(mp):
            continue
        m = json.load(open(mp))["scrollfiesta"]
        ou, ov = m["origin_uv"]
        assert m.get("flip_u", 0) == 0 and m.get("flip_v", 0) == 0, "flip unsupported"
        mk = re.search(r"_w(\d+)_", name)
        if not mk:
            continue
        k = int(mk.group(1))
        mask = tifffile.imread(os.path.join(sd, "mask.tif"))
        H, W = mask.shape
        strip_valid += int((mask >= 255).sum())
        segs.setdefault(k, []).append((sd, float(ou), float(ov), W, H))

# ---- weld each wrap into one full-height segment ----
pieces = []
weld_valid = 0
kept = 0
for k in sorted(segs):
    mem = segs[k]
    u0 = min(s[1] for s in mem)
    u1 = max(s[1] + s[3] * DU for s in mem)
    v0 = min(s[2] for s in mem)
    v1 = max(s[2] + s[4] * DV for s in mem)
    Wk = int(round((u1 - u0) / DU))
    Hk = int(round((v1 - v0) / DV))
    if Wk < 1 or Hk < 1 or Wk * Hk > 950_000_000:
        print(f"  w{k:03d}: skip ({Wk}x{Hk})")
        continue
    X = np.full((Hk, Wk), -1.0, np.float32); Y = X.copy(); Z = X.copy()
    M = np.zeros((Hk, Wk), np.uint8); P = np.zeros((Hk, Wk), np.uint8)
    for (sd, ou, ov, W, H) in mem:
        sx = tifffile.imread(os.path.join(sd, "x.tif"))
        sy = tifffile.imread(os.path.join(sd, "y.tif"))
        sz = tifffile.imread(os.path.join(sd, "z.tif"))
        sm = tifffile.imread(os.path.join(sd, "mask.tif"))
        pp = os.path.join(sd, "provenance.tif")
        sp = tifffile.imread(pp) if os.path.isfile(pp) else np.zeros_like(sm)
        col0 = int(round((ou - u0) / DU)); row0 = int(round((ov - v0) / DV))
        rr, cc = np.nonzero(sm >= 255)
        dr = rr + row0; dc = cc + col0
        empty = M[dr, dc] == 0                     # first strip wins any halo overlap
        di = (dr[empty], dc[empty]); si = (rr[empty], cc[empty])
        X[di] = sx[si]; Y[di] = sy[si]; Z[di] = sz[si]
        M[di] = 255; P[di] = sp[si]
    vm = M >= 255
    nv = int(vm.sum()); weld_valid += nv; kept += 1
    if nv:
        bbox = [[float(X[vm].min()), float(Y[vm].min()), float(Z[vm].min())],
                [float(X[vm].max()), float(Y[vm].max()), float(Z[vm].max())]]
    else:
        bbox = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    nm = f"welded_w{k:03d}"
    od = os.path.join(outseg, nm); os.makedirs(od, exist_ok=True)
    tifffile.imwrite(os.path.join(od, "x.tif"), X)
    tifffile.imwrite(os.path.join(od, "y.tif"), Y)
    tifffile.imwrite(os.path.join(od, "z.tif"), Z)
    tifffile.imwrite(os.path.join(od, "mask.tif"), M)
    tifffile.imwrite(os.path.join(od, "provenance.tif"), P)
    contested = int((P == 3).sum())
    meta = {"scale": [1.0/DU, 1.0/DV], "type": "seg", "format": "tifxyz",
            "uuid": nm, "source": "scrollfiesta-welded",
            "scrollfiesta": {"du": DU, "dv": DV, "flip_u": 0, "flip_v": 0,
                             "origin_uv": [u0, v0], "wrap_k": k,
                             "n_slabs": len(mem), "valid_px": nv,
                             "contested_px": contested,
                             "provenance": "provenance.tif: 0=empty 1=single 2=overlap 3=contested"},
            "bbox": bbox}
    json.dump(meta, open(os.path.join(od, "meta.json"), "w"))   # compact: robust sscanf
    pieces.append({"uuid": nm, "k": k, "n_slabs": len(mem), "W": Wk, "H": Hk,
                   "origin_uv": [u0, v0], "valid_px": nv, "contested_px": contested})

atlas = {"format": "scrollfiesta-atlas-welded", "version": 1, "prefix": "welded",
         "du": DU, "dv": DV, "n_pieces": kept,
         "v_full": [min(p["origin_uv"][1] for p in pieces),
                    max(p["origin_uv"][1] + p["H"] for p in pieces)] if pieces else [0, 0],
         "total_valid_px": weld_valid, "pieces": pieces}
json.dump(atlas, open(os.path.join(outdir, "atlas.json"), "w"), indent=2)
print(f"\nWELD: {kept} full-height wraps written -> {outseg}")
print(f"  strip valid px (sum of 4 strips) = {strip_valid}")
print(f"  welded valid px                  = {weld_valid}  "
      f"({100.0*weld_valid/max(strip_valid,1):.2f}% of strip sum; "
      f"loss = halo double-cover dropped)")
if pieces:
    print(f"  full v-range = {atlas['v_full']}  (~{atlas['v_full'][1]-atlas['v_full'][0]:.0f} vox tall)")
    big = sorted(pieces, key=lambda p: -p["valid_px"])[:5]
    print("  biggest welded wraps:", [(p['uuid'], p['W'], p['H'], p['n_slabs']) for p in big])
