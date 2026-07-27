"""Default locations for the native binaries this pipeline drives."""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

# .../scrollfiesta_public/python/src/scrollunwrap/_paths.py -> scrollfiesta_public
REPO = Path(__file__).resolve().parents[3]
PROJECTS = REPO.parent
VILLA_BIN = PROJECTS / "villa" / "volume-cartographer" / "build-macos" / "bin"

# The Linux/macOS Makefile leaves binaries in src/; Windows MSVC
# (scrollfiesta.sln) drops them in build/{Release,Debug}/ with .exe.
_EXE = ".exe" if sys.platform == "win32" else ""


def _first_existing(*cands: Path) -> Path | None:
    for c in cands:
        if c and Path(c).exists():
            return Path(c)
    return None


def default_cube_mesh() -> Path | None:
    return _first_existing(REPO / "src" / "cube_mesh",
                           REPO / "build" / "Release" / f"cube_mesh{_EXE}",
                           REPO / "build" / "Debug" / f"cube_mesh{_EXE}")


def default_grid_weld() -> Path | None:
    return _first_existing(REPO / "src" / "grid_weld",
                           REPO / "build" / "Release" / f"grid_weld{_EXE}",
                           REPO / "build" / "Debug" / f"grid_weld{_EXE}")


def _which(name: str) -> Path | None:
    """shutil.which as Path | None. (Path(which(...) or "") is Path('.'),
    which exists — it would defeat the not-installed skip logic.)"""
    w = shutil.which(name)
    return Path(w) if w else None


def default_flatboi() -> Path | None:
    return _first_existing(VILLA_BIN / "flatboi", _which("flatboi"))


def default_obj2tifxyz() -> Path | None:
    return _first_existing(VILLA_BIN / "vc_obj2tifxyz_legacy",
                           _which("vc_obj2tifxyz_legacy"))


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
