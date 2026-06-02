"""Prepare a welded scrollfiesta mesh for SLIM flattening.

The welded mesh is multi-component and may not be a topological disk. We split
connected components, keep the largest, clean it to a 2-manifold, fill small
interior holes, and gate on disk topology (one boundary loop, genus 0) — what
``flatboi``'s harmonic + SLIM flattening requires.

We also reorder the vertices from scrollfiesta's **(z,y,x)** convention to the
**(x,y,z)** order ``vc_obj2tifxyz_legacy`` expects, and scale by ``2**level`` so
the coordinates are level-0 scroll voxels (see the plan's coordinate chain).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import trimesh


class PrepError(RuntimeError):
    """Raised when a component cannot be made into a topological disk."""


# --------------------------------------------------------------------------- #
# OBJ load (matches the C reader: first 3 floats of `v`, vertex index of `f`).
# --------------------------------------------------------------------------- #
def load_welded_obj(path) -> trimesh.Trimesh:
    verts: list[tuple[float, float, float]] = []
    faces: list[tuple[int, int, int]] = []
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("v "):
                p = line.split()
                verts.append((float(p[1]), float(p[2]), float(p[3])))
            elif line.startswith("f "):
                idx = [int(tok.split("/")[0]) for tok in line.split()[1:]]
                for i in range(1, len(idx) - 1):           # fan-triangulate n-gons
                    faces.append((idx[0] - 1, idx[i] - 1, idx[i + 1] - 1))
    if not verts or not faces:
        raise PrepError(f"{path}: no vertices/faces parsed")
    return trimesh.Trimesh(vertices=np.asarray(verts, dtype=np.float64),
                           faces=np.asarray(faces, dtype=np.int64), process=False)


def write_obj_vf(path, V: np.ndarray, F: np.ndarray) -> None:
    """Write a plain v/f OBJ (no UV) for flatboi to consume."""
    with open(path, "w") as fh:
        for x, y, z in V:
            fh.write(f"v {x:.6f} {y:.6f} {z:.6f}\n")
        for a, b, c in F + 1:                              # OBJ is 1-indexed
            fh.write(f"f {a} {b} {c}\n")


# --------------------------------------------------------------------------- #
# Topology helpers
# --------------------------------------------------------------------------- #
def split_components(mesh: trimesh.Trimesh) -> list[trimesh.Trimesh]:
    comps = mesh.split(only_watertight=False)
    return list(comps) if len(comps) else [mesh]


def keep_largest(comps: list[trimesh.Trimesh]) -> trimesh.Trimesh:
    return max(comps, key=lambda m: len(m.faces))


def make_manifold(mesh: trimesh.Trimesh) -> trimesh.Trimesh:
    m = mesh.copy()
    m.merge_vertices()
    m.update_faces(m.nondegenerate_faces())
    m.update_faces(m.unique_faces())
    m.remove_unreferenced_vertices()
    return m


def boundary_loops(mesh: trimesh.Trimesh) -> list[list[int]]:
    """Ordered vertex-index loops of the mesh boundary (edges used by 1 face)."""
    es = mesh.edges_sorted
    once = trimesh.grouping.group_rows(es, require_count=1)
    if len(once) == 0:
        return []
    bedges = es[once]
    adj: dict[int, list[int]] = defaultdict(list)
    for a, b in bedges:
        adj[int(a)].append(int(b))
        adj[int(b)].append(int(a))
    seen: set[frozenset] = set()
    loops: list[list[int]] = []
    for start in list(adj):
        for first in adj[start]:
            e0 = frozenset((start, first))
            if e0 in seen:
                continue
            seen.add(e0)
            loop = [start]
            prev, cur = start, first
            while cur != start:
                loop.append(cur)
                nxt = next((x for x in adj[cur]
                            if frozenset((cur, x)) not in seen), None)
                if nxt is None:
                    break
                seen.add(frozenset((cur, nxt)))
                prev, cur = cur, nxt
            loops.append(loop)
    return loops


def euler_genus(mesh: trimesh.Trimesh) -> tuple[int, int, int]:
    """Return (euler_characteristic, genus, n_boundary_loops) for an orientable
    surface: chi = V - E + F, genus = (2 - chi - b) / 2."""
    V = len(mesh.vertices)
    E = len(mesh.edges_unique)
    F = len(mesh.faces)
    chi = V - E + F
    b = len(boundary_loops(mesh))
    genus = (2 - chi - b) // 2
    return chi, genus, b


def is_disk(mesh: trimesh.Trimesh) -> bool:
    chi, genus, b = euler_genus(mesh)
    one_component = mesh.body_count == 1
    return one_component and b == 1 and genus == 0 and chi == 1


def fill_small_interior_holes(mesh: trimesh.Trimesh, max_edges: int) -> trimesh.Trimesh:
    """Fan-fill interior boundary loops with <= max_edges edges, keeping the
    largest loop (the outer boundary) open. Winding is fixed afterward."""
    if max_edges <= 0:
        return mesh
    loops = boundary_loops(mesh)
    if len(loops) <= 1:
        return mesh
    outer = max(range(len(loops)), key=lambda i: len(loops[i]))
    V = mesh.vertices.tolist()
    F = mesh.faces.tolist()
    for i, loop in enumerate(loops):
        if i == outer or len(loop) > max_edges:
            continue
        c = len(V)
        V.append(np.asarray([mesh.vertices[j] for j in loop]).mean(axis=0).tolist())
        for k in range(len(loop)):
            F.append([loop[k], loop[(k + 1) % len(loop)], c])
    out = trimesh.Trimesh(vertices=np.asarray(V), faces=np.asarray(F), process=False)
    out = make_manifold(out)
    trimesh.repair.fix_normals(out)                        # consistent winding
    return out


# --------------------------------------------------------------------------- #
# Coordinate reorder + level scaling (the crux — tested by test_tifxyz_coords)
# --------------------------------------------------------------------------- #
def to_level0_xyz(V: np.ndarray, level: int) -> np.ndarray:
    """scrollfiesta (z,y,x) at level L  ->  (x,y,z) in level-0 voxels.

    Reverse the triplet and multiply by 2**level. After flatboi preserves these
    coordinates, vc_obj2tifxyz_legacy writes x.tif/y.tif/z.tif = true scroll x/y/z.
    """
    return np.asarray(V, dtype=np.float64)[:, ::-1] * float(2 ** int(level))


# --------------------------------------------------------------------------- #
@dataclass
class PreppedPatch:
    mesh: trimesh.Trimesh           # single disk, vertices in (x,y,z) level-0
    info: dict


def select_components(welded: trimesh.Trimesh, *, min_faces: int = 1,
                      max_components: int = 0, mode: str = "all"
                      ) -> list[trimesh.Trimesh]:
    """Connected components of the welded mesh, dropped below ``min_faces`` and
    ranked largest-first. ``mode='largest'`` keeps only the biggest; a positive
    ``max_components`` caps the count. Each surviving component is flattened
    independently (one sheet -> one tifxyz)."""
    comps = [c for c in split_components(welded) if len(c.faces) >= min_faces]
    comps.sort(key=lambda c: len(c.faces), reverse=True)
    if mode == "largest":
        comps = comps[:1]
    elif max_components and max_components > 0:
        comps = comps[:max_components]
    return comps


def prep_component(component: trimesh.Trimesh, level: int = 0, *,
                   fill_max_edges: int = 30) -> PreppedPatch:
    """Clean ONE connected component to a topological disk and reorder/scale to
    level-0 (x,y,z). Raises PrepError if it cannot be made a disk."""
    patch = make_manifold(component)
    patch = fill_small_interior_holes(patch, fill_max_edges)
    patch = make_manifold(patch)
    chi, genus, b = euler_genus(patch)
    info = {
        "kept_vertices": int(len(patch.vertices)),
        "kept_faces": int(len(patch.faces)),
        "euler_characteristic": int(chi),
        "genus": int(genus),
        "boundary_loops": int(b),
    }
    if not (patch.body_count == 1 and b == 1 and genus == 0):
        raise PrepError(
            f"component is not a topological disk "
            f"(boundary_loops={b}, genus={genus}, chi={chi}); "
            f"raise --fill-max-edges or skip. info={info}")
    patch = patch.copy()
    patch.vertices = to_level0_xyz(patch.vertices, level)
    return PreppedPatch(mesh=patch, info=info)
