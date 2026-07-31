"""Driver: given cube+comp, auto-map the step0 cloud and step10 comp, then
invoke render_panels. Usage: drive_panels.py <root> <cube> <comp> <outdir> [extra args]
"""
import sys, os, glob, subprocess
import numpy as np
from render_fig3_panels import load_obj

root, cube, comp, outdir = sys.argv[1:5]
extra = sys.argv[5:]
d = os.path.join(root, cube)
fill_p = os.path.join(d, f"{cube}_step8_holefill", f"{cube}_step8_holefill_{comp}.obj")
bpa_p = os.path.join(d, f"{cube}_step7_cc_bpa", f"{cube}_step7_cc_bpa_{comp}.obj")
fv, _ = load_obj(fill_p)
lo, hi = fv.min(0), fv.max(0)
c8 = fv.mean(0)

best, bp = -1, None
for f in sorted(glob.glob(os.path.join(d, f"{cube}_step0_mls", "*_0*.obj"))):
    if "pre" in os.path.basename(f):
        continue
    pv, _ = load_obj(f)
    plo, phi = pv.min(0), pv.max(0)
    inter = np.maximum(0, np.minimum(hi, phi) - np.maximum(lo, plo)).prod()
    frac = inter / max((hi - lo).prod(), 1e-9)
    # prefer full coverage of the sheet with the least extra volume
    score = frac - 0.1 * ((phi - plo).prod() / max((hi - lo).prod(), 1e-9))
    if score > best:
        best, bp = score, f
print("cloud:", os.path.basename(bp))

bestd, cp = 1e9, None
for f in sorted(glob.glob(os.path.join(d, f"{cube}_step10_qem", "*_0*.obj"))):
    cv, _ = load_obj(f)
    dd = float(np.linalg.norm(cv.mean(0) - c8))
    if dd < bestd:
        bestd, cp = dd, f
print(f"cvt: {os.path.basename(cp)} (centroid dist {bestd:.1f})")

subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "render_fig3_panels.py"),
                bp, bpa_p, fill_p, cp, outdir] + extra, check=True)
