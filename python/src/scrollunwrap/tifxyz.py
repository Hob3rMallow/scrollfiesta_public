"""Convert a UV-mapped OBJ to a tifxyz directory via vc_obj2tifxyz_legacy.

The OBJ fed here already has vertices in (x,y,z) level-0 order (done in meshprep)
and per-corner UVs (added by flatboi), so the converter — which does no axis
swap — writes x.tif/y.tif/z.tif = true scroll x/y/z.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path


def run_obj2tifxyz(uv_obj, out_dir, *, binary, step_size: int = 20,
                   timeout: float = 1800.0, require_mask: bool = False
                   ) -> tuple[Path, dict, str]:
    """Run the legacy converter and validate its outputs.

    Returns (out_dir, meta_json, stdout). Raises on missing outputs.
    """
    out_dir = Path(out_dir)
    # The converter creates the output dir itself and refuses if it exists.
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    cmd = [str(binary), str(uv_obj), str(out_dir), str(int(step_size))]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    missing = [fn for fn in ("x.tif", "y.tif", "z.tif", "meta.json")
               if not (out_dir / fn).exists()]
    if require_mask and not (out_dir / "mask.tif").exists():
        missing.append("mask.tif")
    if missing or proc.returncode != 0:
        raise RuntimeError(
            f"vc_obj2tifxyz_legacy failed (rc={proc.returncode}, missing={missing})\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout[-2000:]}\nstderr:\n{proc.stderr[-2000:]}")
    meta = json.loads((out_dir / "meta.json").read_text())
    if meta.get("format") != "tifxyz":
        raise RuntimeError(f"unexpected meta.format={meta.get('format')!r}")
    return out_dir, meta, proc.stdout
