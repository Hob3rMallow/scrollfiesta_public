"""Geometry / topology validation and report.

Checks connected components, manifoldness, Euler characteristic + genus, and
boundary loops on meshes; flipped-triangle count and bijectivity proxy on the
UV map. Results go to report.json and are asserted in tests.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import trimesh

from .meshprep import euler_genus


def analyze_welded(mesh: trimesh.Trimesh) -> dict:
    comps = mesh.split(only_watertight=False)
    comps = list(comps) if len(comps) else [mesh]
    sizes = sorted((int(len(c.faces)) for c in comps), reverse=True)
    chi, genus, b = euler_genus(mesh)
    return {
        "n_vertices": int(len(mesh.vertices)),
        "n_faces": int(len(mesh.faces)),
        "n_components": len(comps),
        "component_face_counts": sizes[:10],
        "is_winding_consistent": bool(mesh.is_winding_consistent),
        "euler_characteristic": int(chi),
        "genus": int(genus),
        "boundary_loops": int(b),
    }


def load_obj_uv(path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse an OBJ with per-corner UVs. Returns (V, faces, face_uv) where
    face_uv has shape (T,3,2) — the UV of each triangle corner."""
    V: list[list[float]] = []
    VT: list[list[float]] = []
    faces: list[tuple[int, int, int]] = []
    fuv: list[tuple[int, int, int]] = []
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("v "):
                p = line.split()
                V.append([float(p[1]), float(p[2]), float(p[3])])
            elif line.startswith("vt "):
                p = line.split()
                VT.append([float(p[1]), float(p[2])])
            elif line.startswith("f "):
                toks = line.split()[1:]
                vi, ti = [], []
                for t in toks:
                    parts = t.split("/")
                    vi.append(int(parts[0]) - 1)
                    ti.append(int(parts[1]) - 1 if len(parts) > 1 and parts[1] else -1)
                for i in range(1, len(vi) - 1):
                    faces.append((vi[0], vi[i], vi[i + 1]))
                    fuv.append((ti[0], ti[i], ti[i + 1]))
    Va = np.asarray(V, dtype=np.float64)
    VTa = np.asarray(VT, dtype=np.float64) if VT else np.zeros((0, 2))
    fuv_a = np.asarray(fuv, dtype=np.int64)
    if len(VTa) and (fuv_a >= 0).all():
        face_uv = VTa[fuv_a]                       # (T,3,2)
    else:
        face_uv = np.zeros((len(faces), 3, 2))
    return Va, np.asarray(faces, dtype=np.int64), face_uv


def _signed_areas(face_uv: np.ndarray) -> np.ndarray:
    a = face_uv[:, 1] - face_uv[:, 0]
    b = face_uv[:, 2] - face_uv[:, 0]
    return 0.5 * (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])


def analyze_uv(face_uv: np.ndarray) -> dict:
    """Flipped-triangle count (orientation reversals vs the majority sign) and
    the UV bbox. A valid flattening has all triangle area-signs equal."""
    if len(face_uv) == 0:
        return {"n_faces": 0, "flipped_uv_triangles": 0, "flipped_fraction": 0.0}
    s = _signed_areas(face_uv)
    pos = int((s > 0).sum())
    neg = int((s < 0).sum())
    flipped = min(pos, neg)
    uv = face_uv.reshape(-1, 2)
    return {
        "n_faces": int(len(s)),
        "flipped_uv_triangles": int(flipped),
        "flipped_fraction": float(flipped / len(s)),
        "degenerate_uv_triangles": int((s == 0).sum()),
        "uv_bbox": [float(uv[:, 0].min()), float(uv[:, 1].min()),
                    float(uv[:, 0].max()), float(uv[:, 1].max())],
    }


def build_report(*, welded=None, prepped_info=None, flat_obj=None,
                 weld_report=None, tifxyz_meta=None, extra=None) -> dict:
    report: dict = {}
    if welded is not None:
        report["welded"] = analyze_welded(welded)
    if prepped_info is not None:
        report["prepped"] = prepped_info
    if flat_obj is not None and Path(flat_obj).exists():
        _, _, face_uv = load_obj_uv(flat_obj)
        report["uv"] = analyze_uv(face_uv)
    if weld_report is not None:
        report["weld_report"] = weld_report
    if tifxyz_meta is not None:
        report["tifxyz_meta"] = tifxyz_meta
    if extra:
        report.update(extra)
    return report


def write_report(path, report: dict) -> None:
    Path(path).write_text(json.dumps(report, indent=2))
