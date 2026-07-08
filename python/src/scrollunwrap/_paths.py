"""Default locations for the native binaries this pipeline drives."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

# .../scrollfiesta_public/python/src/scrollunwrap/_paths.py -> scrollfiesta_public
REPO = Path(__file__).resolve().parents[3]
PROJECTS = REPO.parent
VILLA_BIN = PROJECTS / "villa" / "volume-cartographer" / "build-macos" / "bin"


def _first_existing(*cands: Path) -> Path | None:
    for c in cands:
        if c and Path(c).exists():
            return Path(c)
    return None


def default_cube_mesh() -> Path | None:
    return _first_existing(REPO / "src" / "cube_mesh")


def default_grid_weld() -> Path | None:
    return _first_existing(REPO / "src" / "grid_weld")


def default_flatboi() -> Path | None:
    return _first_existing(VILLA_BIN / "flatboi", Path(shutil.which("flatboi") or ""))


def default_obj2tifxyz() -> Path | None:
    return _first_existing(VILLA_BIN / "vc_obj2tifxyz_legacy",
                           Path(shutil.which("vc_obj2tifxyz_legacy") or ""))


def default_env(threads_per_cube: int | None = None) -> dict:
    """Environment for the native binaries: make Homebrew dylibs findable and
    cap per-process threads."""
    env = dict(os.environ)
    extra = ["/opt/homebrew/lib", "/opt/homebrew/opt/libomp/lib",
             "/opt/homebrew/opt/libtiff/lib"]
    prev = env.get("DYLD_FALLBACK_LIBRARY_PATH", "")
    env["DYLD_FALLBACK_LIBRARY_PATH"] = ":".join([p for p in extra if Path(p).exists()] +
                                                 ([prev] if prev else []))
    if threads_per_cube:
        env["VESUVIUS_THREADS"] = str(threads_per_cube)
        env["OMP_NUM_THREADS"] = str(threads_per_cube)
    return env
