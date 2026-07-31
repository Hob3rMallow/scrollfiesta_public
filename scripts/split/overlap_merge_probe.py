#!/usr/bin/env python3
"""Inspection-only A/B probe for split versus local overlap healing.

The production overlap separator keeps the lifted-multicut (split) result.  This
tool constructs the competing merge hypothesis without writing it back:

1. detect the close contact between two already-split sheets;
2. grow a small face neighbourhood on both sheets;
3. project that neighbourhood to a common PCA plane;
4. coalesce only near-duplicate, patch-interior samples and Delaunay the union;
5. keep the exterior fixed and fit the healed topology to the nnU-Net surface;
6. dump every stage plus topology, fit, quality, and RAW seam diagnostics.

It is deliberately conservative.  A result can be called ``eligible`` but this
calibration tool never replaces a production component.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage, sparse
from scipy.sparse.linalg import spsolve
from scipy.spatial import Delaunay, cKDTree


def cross2(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def load_obj(path: Path):
    verts, faces, colors = [], [], []
    with path.open("r", encoding="utf-8") as src:
        for line in src:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "v" and len(fields) >= 4:
                verts.append([float(x) for x in fields[1:4]])
                if len(fields) >= 7:
                    colors.append([float(x) for x in fields[4:7]])
                else:
                    colors.append([math.nan, math.nan, math.nan])
            elif fields[0] == "f" and len(fields) >= 4:
                poly = [int(x.split("/")[0]) - 1 for x in fields[1:]]
                for i in range(1, len(poly) - 1):
                    faces.append([poly[0], poly[i], poly[i + 1]])
    if not verts or not faces:
        raise ValueError(f"empty or unreadable OBJ: {path}")
    return (np.asarray(verts, dtype=np.float64),
            np.asarray(faces, dtype=np.int32),
            np.asarray(colors, dtype=np.float64))


def write_obj(path: Path, verts, faces, vertex_colors=None, comment=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as dst:
        if comment:
            dst.write(f"# {comment}\n")
        dst.write(f"# {len(verts)} vertices, {len(faces)} faces\n")
        if vertex_colors is None:
            for p in verts:
                dst.write(f"v {p[0]:.9g} {p[1]:.9g} {p[2]:.9g}\n")
        else:
            for p, c in zip(verts, vertex_colors):
                dst.write("v %.9g %.9g %.9g %.5g %.5g %.5g\n" %
                          (p[0], p[1], p[2], c[0], c[1], c[2]))
        for tri in faces:
            dst.write(f"f {tri[0] + 1} {tri[1] + 1} {tri[2] + 1}\n")


def write_face_colors(path: Path, verts, faces, colors, comment=""):
    """Write exact per-face colors by intentionally duplicating face vertices."""
    soup_v = verts[faces].reshape((-1, 3))
    soup_f = np.arange(len(soup_v), dtype=np.int32).reshape((-1, 3))
    soup_c = np.repeat(np.asarray(colors, dtype=np.float64), 3, axis=0)
    write_obj(path, soup_v, soup_f, soup_c, comment)


class UnionFind:
    def __init__(self, n):
        self.parent = np.arange(n, dtype=np.int64)

    def find(self, x):
        p = self.parent
        while p[x] != x:
            p[x] = p[p[x]]
            x = int(p[x])
        return x

    def join(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[b] = a


def edge_incidence(faces):
    edges = defaultdict(list)
    for fi, tri in enumerate(faces):
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            edges[(min(int(a), int(b)), max(int(a), int(b)))].append(fi)
    return edges


def face_adjacency(faces):
    adjacency = [[] for _ in range(len(faces))]
    for incident in edge_incidence(faces).values():
        if len(incident) == 2:
            a, b = incident
            adjacency[a].append(b)
            adjacency[b].append(a)
    return adjacency


def grow_faces(seed, adjacency, rings):
    selected = np.array(seed, dtype=bool, copy=True)
    for _ in range(rings):
        grown = selected.copy()
        for fi in np.flatnonzero(selected):
            grown[adjacency[fi]] = True
        selected = grown
    return selected


def common_frame(verts, faces, seed):
    ids = np.unique(faces[seed])
    if len(ids) < 6:
        raise ValueError("too few contact vertices for a common plane")
    points = verts[ids]
    center = points.mean(axis=0)
    cov = np.cov((points - center).T)
    values, vectors = np.linalg.eigh(cov)
    normal = vectors[:, 0]
    axis_u = vectors[:, 2]
    axis_v = np.cross(normal, axis_u)
    axis_v /= np.linalg.norm(axis_v)
    if normal[np.argmax(np.abs(normal))] < 0:
        normal = -normal
        axis_v = -axis_v
    projected = np.column_stack(((verts - center) @ axis_u,
                                 (verts - center) @ axis_v))
    thickness = float(np.sqrt(max(values[0], 0.0)))
    return center, axis_u, axis_v, normal, projected, thickness


def auto_contact_seed(va, fa, vb, fb, distance, min_vertices):
    da, _ = cKDTree(vb).query(va, workers=-1)
    db, _ = cKDTree(va).query(vb, workers=-1)
    sa = np.count_nonzero(da[fa] <= distance, axis=1) >= min_vertices
    sb = np.count_nonzero(db[fb] <= distance, axis=1) >= min_vertices
    if np.count_nonzero(sa) < 8 or np.count_nonzero(sb) < 8:
        sa = np.count_nonzero(da[fa] <= distance, axis=1) >= 2
        sb = np.count_nonzero(db[fb] <= distance, axis=1) >= 2
    return sa, sb


def highlighted_seed(path, faces, expected_nv, expected_nf):
    verts, hfaces, colors = load_obj(path)
    if len(verts) != expected_nv or len(hfaces) != expected_nf:
        raise ValueError("highlight OBJ does not match the combined sheet pair")
    red = ((colors[:, 0] > 0.9) & (colors[:, 1] < 0.2) &
           (colors[:, 2] < 0.2))
    return np.all(red[faces], axis=1)


class TriangleUnion:
    """Point-in-union queries over projected source triangles."""
    def __init__(self, projected, triangles, cell=0.5):
        self.p = projected
        self.triangles = triangles
        self.cell = cell
        self.grid = defaultdict(list)
        for ti, tri in enumerate(triangles):
            xy = projected[tri]
            lo = np.floor(xy.min(axis=0) / cell).astype(int)
            hi = np.floor(xy.max(axis=0) / cell).astype(int)
            for y in range(lo[1], hi[1] + 1):
                for x in range(lo[0], hi[0] + 1):
                    self.grid[(x, y)].append(ti)

    @staticmethod
    def _inside(point, tri):
        a, b, c = tri
        cross = np.array((cross2(b - a, point - a),
                          cross2(c - b, point - b),
                          cross2(a - c, point - c)))
        return bool(np.all(cross >= -1e-7) or np.all(cross <= 1e-7))

    def contains(self, points):
        answer = np.zeros(len(points), dtype=bool)
        for pi, point in enumerate(points):
            key = tuple(np.floor(point / self.cell).astype(int))
            for ti in self.grid.get(key, ()):
                if self._inside(point, self.p[self.triangles[ti]]):
                    answer[pi] = True
                    break
        return answer


def cluster_points(ids, projected, distance):
    if len(ids) == 0:
        return []
    uf = UnionFind(len(ids))
    if distance > 0:
        for a, b in cKDTree(projected[ids]).query_pairs(distance):
            uf.join(a, b)
    groups = defaultdict(list)
    for i, vertex in enumerate(ids):
        groups[uf.find(i)].append(int(vertex))
    return list(groups.values())


def add_steiner(points, domain, spacing):
    if spacing <= 0:
        return points
    lo, hi = points.min(axis=0), points.max(axis=0)
    rows = []
    y = lo[1]
    row = 0
    while y <= hi[1]:
        x = lo[0] + (0.5 * spacing if row & 1 else 0.0)
        while x <= hi[0]:
            rows.append((x, y))
            x += spacing
        y += spacing * math.sqrt(3.0) / 2.0
        row += 1
    if not rows:
        return points
    extra = np.asarray(rows, dtype=np.float64)
    extra = extra[domain.contains(extra)]
    if len(extra):
        distance, _ = cKDTree(points).query(extra)
        extra = extra[distance >= 0.35 * spacing]
    return np.vstack((points, extra)) if len(extra) else points


def build_candidate(verts, faces, selected, projected, center, axis_u, axis_v,
                    merge_distance, steiner_spacing):
    kept_faces = faces[~selected]
    source_patch = faces[selected]
    patch_vertices = np.unique(source_patch)
    kept_vertices = set(int(x) for x in np.unique(kept_faces))
    pinned = np.asarray([x for x in patch_vertices if int(x) in kept_vertices],
                        dtype=np.int32)
    interior = np.asarray([x for x in patch_vertices if int(x) not in kept_vertices],
                          dtype=np.int32)
    groups = cluster_points(interior, projected, merge_distance)
    interior_2d = np.asarray([projected[g].mean(axis=0) for g in groups],
                             dtype=np.float64)
    points_2d = np.vstack((projected[pinned], interior_2d))
    n_without_steiner = len(points_2d)
    domain = TriangleUnion(projected, source_patch)
    points_2d = add_steiner(points_2d, domain, steiner_spacing)

    added_3d = np.asarray([center + p[0] * axis_u + p[1] * axis_v
                           for p in points_2d[len(pinned):]], dtype=np.float64)
    candidate_verts = np.vstack((verts, added_3d))
    point_global = np.concatenate((
        pinned.astype(np.int64),
        np.arange(len(verts), len(verts) + len(added_3d), dtype=np.int64)))

    triangulation = Delaunay(points_2d, qhull_options="Qbb Qc Qz Q12")
    simplex = triangulation.simplices
    probes = np.vstack((
        points_2d[simplex].mean(axis=1),
        0.5 * (points_2d[simplex[:, 0]] + points_2d[simplex[:, 1]]),
        0.5 * (points_2d[simplex[:, 1]] + points_2d[simplex[:, 2]]),
        0.5 * (points_2d[simplex[:, 2]] + points_2d[simplex[:, 0]])))
    inside = domain.contains(probes).reshape((4, len(simplex)))
    simplex = simplex[np.all(inside, axis=0)]
    patch_faces = point_global[simplex].astype(np.int32)

    a = points_2d[simplex[:, 0]]
    b = points_2d[simplex[:, 1]]
    c = points_2d[simplex[:, 2]]
    reverse = cross2(b - a, c - a) < 0
    patch_faces[reverse] = patch_faces[reverse][:, [0, 2, 1]]

    candidate_faces = np.vstack((kept_faces, patch_faces))
    patch_mask = np.zeros(len(candidate_faces), dtype=bool)
    patch_mask[len(kept_faces):] = True

    def orientation_counts(all_faces):
        same = opposite = 0
        for edge, incident in edge_incidence(all_faces).items():
            if len(incident) != 2 or patch_mask[incident[0]] == patch_mask[incident[1]]:
                continue
            directions = []
            for fi in incident:
                tri = all_faces[fi]
                direction = 0
                for j in range(3):
                    a, b = int(tri[j]), int(tri[(j + 1) % 3])
                    if (a, b) == edge:
                        direction = 1
                    elif (b, a) == edge:
                        direction = -1
                directions.append(direction)
            if directions[0] == directions[1]:
                same += 1
            else:
                opposite += 1
        return same, opposite

    same, opposite = orientation_counts(candidate_faces)
    if same > opposite:
        candidate_faces[patch_mask] = candidate_faces[patch_mask][:, [0, 2, 1]]
        same, opposite = orientation_counts(candidate_faces)
    pin_mask = np.zeros(len(candidate_verts), dtype=bool)
    pin_mask[np.unique(kept_faces)] = True

    used = np.unique(candidate_faces)
    remap = np.full(len(candidate_verts), -1, dtype=np.int64)
    remap[used] = np.arange(len(used), dtype=np.int64)
    result = {
        "verts": candidate_verts[used],
        "faces": remap[candidate_faces].astype(np.int32),
        "patch_mask": patch_mask,
        "pin_mask": pin_mask[used],
        "points_input": int(n_without_steiner),
        "points_steiner": int(len(points_2d) - n_without_steiner),
        "groups": int(len(groups)),
        "coalesced": int(len(interior) - len(groups)),
        "pinned_patch_vertices": int(len(pinned)),
        "source_patch_faces": int(len(source_patch)),
        "new_patch_faces": int(len(patch_faces)),
        "interface_edges_same_direction": int(same),
        "interface_edges_opposite_direction": int(opposite),
        "points_2d": points_2d,
        "patch_faces_2d": simplex.astype(np.int32),
    }
    return result


def ordered_boundary_cycles(faces):
    graph = defaultdict(list)
    for edge, incident in edge_incidence(faces).items():
        if len(incident) == 1:
            a, b = edge
            graph[a].append(b)
            graph[b].append(a)
    if any(len(nbr) != 2 for nbr in graph.values()):
        return []
    remaining = set(graph)
    cycles = []
    while remaining:
        start = min(remaining)
        previous, current = -1, start
        cycle = []
        while True:
            cycle.append(current)
            remaining.discard(current)
            nbr = graph[current]
            following = nbr[0] if nbr[0] != previous else nbr[1]
            previous, current = current, following
            if current == start:
                break
            if len(cycle) > len(graph):
                return []
        cycles.append(cycle)
    return cycles


def point_in_triangle_2d(point, tri):
    a, b, c = tri
    signs = np.array((cross2(b - a, point - a),
                      cross2(c - b, point - b),
                      cross2(a - c, point - c)))
    return bool(np.all(signs >= -1e-9) or np.all(signs <= 1e-9))


def ear_clip_loop(loop, uv):
    loop = list(loop)
    polygon = uv[loop]
    area2 = float(np.sum(cross2(polygon, np.roll(polygon, -1, axis=0))))
    if area2 < 0:
        loop.reverse()
    remaining = list(range(len(loop)))
    triangles = []
    guard = 0
    while len(remaining) > 3 and guard < len(loop) * len(loop):
        guard += 1
        found = False
        for j, current in enumerate(remaining):
            previous = remaining[j - 1]
            following = remaining[(j + 1) % len(remaining)]
            tri = uv[[loop[previous], loop[current], loop[following]]]
            if cross2(tri[1] - tri[0], tri[2] - tri[0]) <= 1e-10:
                continue
            blocked = False
            for other in remaining:
                if other in (previous, current, following):
                    continue
                if point_in_triangle_2d(uv[loop[other]], tri):
                    blocked = True
                    break
            if blocked:
                continue
            triangles.append([loop[previous], loop[current], loop[following]])
            del remaining[j]
            found = True
            break
        if not found:
            return None
    if len(remaining) == 3:
        triangles.append([loop[x] for x in remaining])
    return np.asarray(triangles, dtype=np.int32)


def fill_small_boundary_loops(candidate, center, axis_u, axis_v,
                              max_edges, max_area):
    verts = candidate["verts"]
    faces = candidate["faces"]
    uv = np.column_stack(((verts - center) @ axis_u,
                          (verts - center) @ axis_v))
    old_edges = edge_incidence(faces)
    additions = []
    filled = []
    skipped = []
    for cycle in ordered_boundary_cycles(faces):
        polygon = uv[cycle]
        area = 0.5 * abs(float(np.sum(cross2(
            polygon, np.roll(polygon, -1, axis=0)))))
        if len(cycle) > max_edges or area > max_area:
            skipped.append({"edges": len(cycle), "projected_area": area})
            continue
        triangles = ear_clip_loop(cycle, uv)
        if triangles is None or len(triangles) != len(cycle) - 2:
            skipped.append({"edges": len(cycle), "projected_area": area,
                            "reason": "triangulation_failed"})
            continue
        same = opposite = 0
        for edge in ((min(int(a), int(b)), max(int(a), int(b)))
                     for tri in triangles
                     for a, b in ((tri[0], tri[1]), (tri[1], tri[2]),
                                  (tri[2], tri[0]))):
            incident = old_edges.get(edge, ())
            if len(incident) != 1:
                continue
            old_face = faces[incident[0]]
            def direction(tri, target):
                for k in range(3):
                    a, b = int(tri[k]), int(tri[(k + 1) % 3])
                    if (a, b) == target:
                        return 1
                    if (b, a) == target:
                        return -1
                return 0
            old_direction = direction(old_face, edge)
            new_direction = 0
            for tri in triangles:
                value = direction(tri, edge)
                if value:
                    new_direction = value
                    break
            if old_direction == new_direction:
                same += 1
            else:
                opposite += 1
        if same > opposite:
            triangles = triangles[:, [0, 2, 1]]
        additions.append(triangles)
        filled.append({"edges": len(cycle), "projected_area": area,
                       "faces_added": len(triangles)})
    if additions:
        added = np.vstack(additions)
        candidate["faces"] = np.vstack((faces, added))
        candidate["patch_mask"] = np.concatenate((
            candidate["patch_mask"], np.ones(len(added), dtype=bool)))
    candidate["small_holes_filled"] = filled
    candidate["boundary_loops_left_unfilled"] = skipped
    candidate["holefill_faces_added"] = int(sum(x["faces_added"] for x in filled))
    return candidate


def mesh_topology(verts, faces):
    edges = edge_incidence(faces)
    boundary = [edge for edge, incident in edges.items() if len(incident) == 1]
    nonmanifold = sum(len(incident) > 2 for incident in edges.values())
    same_direction = 0
    for edge, incident in edges.items():
        if len(incident) != 2:
            continue
        directions = []
        for fi in incident:
            tri = faces[fi]
            direction = 0
            for k in range(3):
                a, b = int(tri[k]), int(tri[(k + 1) % 3])
                if (a, b) == edge:
                    direction = 1
                elif (b, a) == edge:
                    direction = -1
            directions.append(direction)
        if directions[0] == directions[1]:
            same_direction += 1
    boundary_graph = defaultdict(list)
    for a, b in boundary:
        boundary_graph[a].append(b)
        boundary_graph[b].append(a)
    bad_boundary_vertices = sum(len(nbr) != 2 for nbr in boundary_graph.values())
    unseen = set(boundary_graph)
    loop_sizes = []
    while unseen:
        start = unseen.pop()
        todo = [start]
        count = 0
        while todo:
            vertex = todo.pop()
            count += 1
            for nbr in boundary_graph[vertex]:
                if nbr in unseen:
                    unseen.remove(nbr)
                    todo.append(nbr)
        loop_sizes.append(count)

    adjacency = [[] for _ in range(len(faces))]
    for incident in edges.values():
        if len(incident) == 2:
            a, b = incident
            adjacency[a].append(b)
            adjacency[b].append(a)
    unseen_faces = set(range(len(faces)))
    component_sizes = []
    while unseen_faces:
        start = unseen_faces.pop()
        todo = [start]
        count = 0
        while todo:
            face = todo.pop()
            count += 1
            for nbr in adjacency[face]:
                if nbr in unseen_faces:
                    unseen_faces.remove(nbr)
                    todo.append(nbr)
        component_sizes.append(count)

    chi = len(verts) - len(edges) + len(faces)
    components = len(component_sizes)
    genus = 0.5 * (2 * components - len(loop_sizes) - chi)
    return {
        "vertices": int(len(verts)),
        "edges": int(len(edges)),
        "faces": int(len(faces)),
        "euler_characteristic": int(chi),
        "face_components": int(components),
        "face_component_sizes": sorted((int(x) for x in component_sizes),
                                       reverse=True),
        "boundary_edges": int(len(boundary)),
        "boundary_loops": int(len(loop_sizes)),
        "boundary_loop_sizes": sorted((int(x) for x in loop_sizes), reverse=True),
        "bad_boundary_vertices": int(bad_boundary_vertices),
        "nonmanifold_edges": int(nonmanifold),
        "same_direction_interior_edges": int(same_direction),
        "orientable_genus_if_consistent": float(genus),
    }


def face_geometry(verts, faces):
    tri = verts[faces]
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    double_area = np.linalg.norm(cross, axis=1)
    normals = np.zeros_like(cross)
    good = double_area > 1e-12
    normals[good] = cross[good] / double_area[good, None]
    e0 = np.linalg.norm(tri[:, 1] - tri[:, 0], axis=1)
    e1 = np.linalg.norm(tri[:, 2] - tri[:, 1], axis=1)
    e2 = np.linalg.norm(tri[:, 0] - tri[:, 2], axis=1)
    denom = e0 * e0 + e1 * e1 + e2 * e2
    quality = np.zeros(len(faces), dtype=np.float64)
    quality[denom > 0] = 2.0 * math.sqrt(3.0) * double_area[denom > 0] / denom[denom > 0]
    return 0.5 * double_area, normals, quality


def distribution(values, weights=None):
    values = np.asarray(values, dtype=np.float64)
    if len(values) == 0:
        return {"count": 0}
    if weights is None:
        weights = np.ones(len(values), dtype=np.float64)
    order = np.argsort(values)
    v, w = values[order], np.asarray(weights, dtype=np.float64)[order]
    cumulative = np.cumsum(w)
    total = cumulative[-1]
    def quantile(q):
        return float(v[min(np.searchsorted(cumulative, q * total), len(v) - 1)])
    return {
        "count": int(len(v)),
        "mean": float(np.average(v, weights=w)),
        "p01": quantile(0.01),
        "p05": quantile(0.05),
        "p50": quantile(0.50),
        "p90": quantile(0.90),
        "p95": quantile(0.95),
        "p99": quantile(0.99),
        "min": float(v[0]),
        "max": float(v[-1]),
    }


def quality_metrics(verts, faces, mask=None):
    area, _, quality = face_geometry(verts, faces)
    if mask is not None:
        area, quality = area[mask], quality[mask]
    result = distribution(quality, area)
    result["degenerate_faces"] = int(np.count_nonzero(area <= 1e-10))
    result["area"] = float(area.sum())
    return result


def cube_origin_from_path(path):
    match = re.search(r"z(\d+)_y(\d+)_x(\d+)", Path(path).name)
    if not match:
        raise ValueError("cannot infer cube origin; pass --cube-origin z y x")
    return np.asarray([float(x) for x in match.groups()], dtype=np.float64)


def prediction_surface(volume):
    foreground = volume > 0
    structure = ndimage.generate_binary_structure(3, 1)
    interior = ndimage.binary_erosion(foreground, structure=structure,
                                      border_value=0)
    return np.argwhere(foreground & ~interior).astype(np.float64)


def fit_to_prediction(verts, faces, pin_mask, surface, origin, alpha,
                      distance_max, outer_iterations):
    fitted = np.array(verts, dtype=np.float64, copy=True)
    local = fitted - origin
    tree = cKDTree(surface)
    edges = np.asarray(list(edge_incidence(faces)), dtype=np.int64)
    neighbours = [[] for _ in range(len(fitted))]
    for a, b in edges:
        neighbours[a].append(int(b))
        neighbours[b].append(int(a))
    unknown = np.flatnonzero(~pin_mask)
    unknown_index = np.full(len(fitted), -1, dtype=np.int64)
    unknown_index[unknown] = np.arange(len(unknown), dtype=np.int64)
    ratio = alpha / max(1.0 - alpha, 1e-9)
    last_distance = np.zeros(len(unknown), dtype=np.float64)
    matched = np.zeros(len(unknown), dtype=bool)

    for _ in range(max(outer_iterations, 1)):
        last_distance, nearest = tree.query(local[unknown], workers=-1)
        targets = surface[nearest]
        matched = last_distance <= distance_max
        rows, cols, data = [], [], []
        rhs = np.zeros((len(unknown), 3), dtype=np.float64)
        for row, vertex in enumerate(unknown):
            degree = len(neighbours[vertex])
            weight = ratio * degree if matched[row] else 0.0
            rows.append(row)
            cols.append(row)
            data.append(degree + weight)
            rhs[row] += weight * targets[row]
            for nbr in neighbours[vertex]:
                col = unknown_index[nbr]
                if col >= 0:
                    rows.append(row)
                    cols.append(int(col))
                    data.append(-1.0)
                else:
                    rhs[row] += local[nbr]
        matrix = sparse.csr_matrix((data, (rows, cols)),
                                   shape=(len(unknown), len(unknown)))
        for coordinate in range(3):
            local[unknown, coordinate] = spsolve(matrix, rhs[:, coordinate])

    fitted = local + origin
    final_distance, _ = tree.query(local[unknown], workers=-1)
    return fitted, {
        "unknown_vertices": int(len(unknown)),
        "matched_vertices": int(np.count_nonzero(matched)),
        "match_distance_before": distribution(last_distance),
        "match_distance_after": distribution(final_distance),
        "alpha": float(alpha),
        "distance_max": float(distance_max),
        "outer_iterations": int(max(outer_iterations, 1)),
    }


def prediction_fit(verts, faces, mask, surface_tree, origin):
    area, _, _ = face_geometry(verts, faces)
    centroids = verts[faces].mean(axis=1) - origin
    distance, _ = surface_tree.query(centroids[mask], workers=-1)
    return distribution(distance, area[mask])


def vertex_normals(verts, faces):
    area, face_normals, _ = face_geometry(verts, faces)
    normals = np.zeros_like(verts)
    for corner in range(3):
        np.add.at(normals, faces[:, corner], face_normals * area[:, None])
    length = np.linalg.norm(normals, axis=1)
    good = length > 1e-12
    normals[good] /= length[good, None]
    return normals


def trilinear(volume, point):
    z, y, x = point
    d, h, w = volume.shape
    if z < 0 or y < 0 or x < 0 or z > d - 1 or y > h - 1 or x > w - 1:
        return math.nan
    z0, y0, x0 = int(math.floor(z)), int(math.floor(y)), int(math.floor(x))
    z1, y1, x1 = min(z0 + 1, d - 1), min(y0 + 1, h - 1), min(x0 + 1, w - 1)
    fz, fy, fx = z - z0, y - y0, x - x0
    result = 0.0
    for iz, wz in ((z0, 1.0 - fz), (z1, fz)):
        for iy, wy in ((y0, 1.0 - fy), (y1, fy)):
            for ix, wx in ((x0, 1.0 - fx), (x1, fx)):
                result += wz * wy * wx * float(volume[iz, iy, ix])
    return result


def normal_max_sample(volume, point, normal, reach=2.0, samples=5):
    values = []
    for step in np.linspace(-reach, reach, samples):
        value = trilinear(volume, point + step * normal)
        if math.isfinite(value):
            values.append(value)
    return max(values) if values else math.nan


def interface_diagnostics(verts, faces, patch_mask, raw, origin,
                          normal_reach=2.0, normal_samples=5):
    area, face_normals, _ = face_geometry(verts, faces)
    centroids = verts[faces].mean(axis=1) - origin
    pairs = []
    for _, incident in edge_incidence(faces).items():
        if len(incident) != 2 or patch_mask[incident[0]] == patch_mask[incident[1]]:
            continue
        patch_face = incident[0] if patch_mask[incident[0]] else incident[1]
        anchor_face = incident[1] if patch_mask[incident[0]] else incident[0]
        patch_value = normal_max_sample(raw, centroids[patch_face],
                                        face_normals[patch_face], normal_reach,
                                        normal_samples)
        anchor_value = normal_max_sample(raw, centroids[anchor_face],
                                         face_normals[anchor_face], normal_reach,
                                         normal_samples)
        if math.isfinite(patch_value) and math.isfinite(anchor_value):
            cosine = float(np.clip(np.dot(face_normals[patch_face],
                                          face_normals[anchor_face]), -1.0, 1.0))
            angle = math.degrees(math.acos(cosine))
            pairs.append((patch_face, anchor_face, patch_value, anchor_value,
                          angle, 0.5 * (area[patch_face] + area[anchor_face])))
    face_energy = np.full(len(faces), math.nan, dtype=np.float64)
    if not pairs:
        return {"pairs": 0}, face_energy
    offset = float(np.mean([p[2] - p[3] for p in pairs]))
    errors = np.asarray([abs((p[2] - offset) - p[3]) for p in pairs])
    angles = np.asarray([p[4] for p in pairs])
    weights = np.asarray([p[5] for p in pairs])
    accum = defaultdict(list)
    for pair, error in zip(pairs, errors):
        accum[pair[0]].append(float(error))
        accum[pair[1]].append(float(error))
    for face, values in accum.items():
        face_energy[face] = float(np.mean(values))
    return {
        "pairs": int(len(pairs)),
        "brightness_offset": offset,
        "quilt_value_discontinuity": distribution(errors, weights),
        "normal_angle_degrees": distribution(angles, weights),
    }, face_energy


def projected_patch_metrics(projected, faces, selected, candidate_points,
                            candidate_faces_2d):
    tri_a = projected[faces[selected]]
    area_sum = 0.5 * np.abs(cross2(tri_a[:, 1] - tri_a[:, 0],
                                  tri_a[:, 2] - tri_a[:, 0])).sum()
    tri_c = candidate_points[candidate_faces_2d]
    union_area = 0.5 * np.abs(cross2(tri_c[:, 1] - tri_c[:, 0],
                                    tri_c[:, 2] - tri_c[:, 0])).sum()
    overlap = max(float(area_sum - union_area), 0.0)
    return {
        "split_projected_area_sum": float(area_sum),
        "healed_projected_union_area": float(union_area),
        "duplicate_projected_area": overlap,
        "duplicate_fraction_of_union": overlap / max(float(union_area), 1e-12),
    }


def heat_colors(values):
    colors = np.full((len(values), 3), (0.65, 0.65, 0.65), dtype=np.float64)
    finite = np.isfinite(values)
    if not np.any(finite):
        return colors
    low, high = np.percentile(values[finite], (5, 95))
    scale = max(high - low, 1e-9)
    t = np.clip((values[finite] - low) / scale, 0.0, 1.0)
    colors[finite, 0] = t
    colors[finite, 1] = 1.0 - 0.55 * t
    colors[finite, 2] = 0.08
    return colors


def write_side_by_side(path, split_verts, split_faces, split_face_colors,
                       merge_verts, merge_faces, merge_face_colors):
    width = max(float(np.ptp(split_verts[:, 2])), float(np.ptp(merge_verts[:, 2])))
    separation = width + 24.0
    left = split_verts[split_faces].reshape((-1, 3)).copy()
    right = merge_verts[merge_faces].reshape((-1, 3)).copy()
    left[:, 2] -= 0.5 * separation
    right[:, 2] += 0.5 * separation
    verts = np.vstack((left, right))
    faces_left = np.arange(len(left), dtype=np.int32).reshape((-1, 3))
    faces_right = (np.arange(len(right), dtype=np.int32) + len(left)).reshape((-1, 3))
    faces = np.vstack((faces_left, faces_right))
    colors = np.vstack((np.repeat(split_face_colors, 3, axis=0),
                        np.repeat(merge_face_colors, 3, axis=0)))
    write_obj(path, verts, faces, colors,
              "left=current split; right=healed and prediction-fitted merge")


def render_readme(path, metrics):
    verdict = metrics["decision"]["verdict"]
    hard = metrics["decision"]["hard_gates"]
    failed = [name for name, passed in hard.items() if not passed]
    pred_split = metrics["prediction_fit"]["split_patch"]
    pred_merge = metrics["prediction_fit"]["merge_patch"]
    quilt_split = metrics["raw_interface"]["split"]["quilt_value_discontinuity"]
    quilt_merge = metrics["raw_interface"]["merge"]["quilt_value_discontinuity"]
    text = f"""# Split versus merge probe

This is an inspection-only hypothesis test. It did not modify a cube component,
the audit grid, or any weld input.

## Provisional result

- Verdict: `{verdict}`
- Failed hard gates: `{', '.join(failed) if failed else 'none'}`
- Split prediction distance: mean `{pred_split['mean']:.4f}`, p90 `{pred_split['p90']:.4f}` vox
- Merge prediction distance: mean `{pred_merge['mean']:.4f}`, p90 `{pred_merge['p90']:.4f}` vox
- Split RAW boundary energy: `{quilt_split['mean']:.4f}`
- Merge RAW boundary energy: `{quilt_merge['mean']:.4f}`

`eligible` means only that the candidate cleared the geometric safety gates.
During calibration, even a provisionally preferred merge still requires visual
review before production write-back.

## Debug stages

- `00_split_pair.obj`: current two-component split hypothesis.
- `01_contact_seed.obj`: detected close-contact faces in red.
- `02_overlap_plus_rings.obj`: the exact local region replaced by the probe.
- `03_merge_planar_unfilled.obj`: common-plane Delaunay before micro-hole fill.
- `04_merge_planar_holefilled.obj`: bounded 4–8-edge hole fill applied.
- `05_merge_pred_fitted.obj`: healed topology after local nnU-Net snap-back.
- `06_merge_patch_only.obj`: fitted healed patch alone.
- `07_merge_quilt_energy.obj`: RAW seam-energy heat at the patch boundary.
- `08_merge_pred_displacement.obj`: plane-to-pred displacement heat.
- `A_split_output.obj`, `B_merge_output.obj`: the two direct alternatives.
- `C_split_vs_merge_side_by_side.obj`: left/right colored comparison.
- Files ending in `_topology.obj` preserve shared vertex indices for audits;
  the exact per-face color stages intentionally use face-soup duplication.
- `metrics.json`: full gates and unrounded measurements.

## Acceptance rule used here

Hard gates precede scoring: one intended face component, edge manifoldness,
simple boundary loops, expected boundary-loop count, consistently oriented
interface edges, no fitted patch fold, adequate low-tail triangle quality, and
prediction p95 no more than the configured margin worse than split. Evidence is
then Pareto-like: prediction mean/p90, RAW seam energy, and interface bend may
each worsen only by a small hysteresis margin. A tie remains `review`, never an
automatic merge.
"""
    path.write_text(text, encoding="utf-8", newline="\n")


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sheet_a", type=Path)
    parser.add_argument("sheet_b", type=Path)
    parser.add_argument("pred", type=Path, help="original nnU-Net prediction TIFF")
    parser.add_argument("raw", type=Path, help="original volumetric RAW TIFF")
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--seed-highlight", type=Path,
                        help="combined OBJ whose red vertices mark contact faces")
    parser.add_argument("--contact-distance", type=float, default=0.9)
    parser.add_argument("--seed-vertices", type=int, choices=(2, 3), default=3)
    parser.add_argument("--rings", type=int, default=3)
    parser.add_argument("--merge-distance", type=float, default=0.10)
    parser.add_argument("--steiner-spacing", type=float, default=0.0)
    parser.add_argument("--holefill-max-edges", type=int, default=8)
    parser.add_argument("--holefill-max-area", type=float, default=2.0)
    parser.add_argument("--cube-origin", type=float, nargs=3,
                        metavar=("Z", "Y", "X"))
    parser.add_argument("--snap-alpha", type=float, default=0.30)
    parser.add_argument("--snap-distance", type=float, default=3.0)
    parser.add_argument("--snap-outer", type=int, default=1)
    parser.add_argument("--pred-p95-margin", type=float, default=0.35)
    parser.add_argument("--pred-score-margin", type=float, default=0.15)
    parser.add_argument("--quilt-margin", type=float, default=2.0)
    parser.add_argument("--bend-p90-max", type=float, default=45.0)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    va, fa, _ = load_obj(args.sheet_a)
    vb, fb, _ = load_obj(args.sheet_b)
    verts = np.vstack((va, vb))
    faces = np.vstack((fa, fb + len(va)))
    face_cut = len(fa)

    if args.seed_highlight:
        seed = highlighted_seed(args.seed_highlight, faces, len(verts), len(faces))
        seed_a, seed_b = seed[:face_cut], seed[face_cut:]
        seed_mode = "highlight"
    else:
        seed_a, seed_b = auto_contact_seed(va, fa, vb, fb,
                                           args.contact_distance,
                                           args.seed_vertices)
        seed = np.concatenate((seed_a, seed_b))
        seed_mode = "nearest_surface"
    if np.count_nonzero(seed_a) == 0 or np.count_nonzero(seed_b) == 0:
        raise RuntimeError("no two-sided close contact found")

    center, axis_u, axis_v, normal, projected, thickness = common_frame(
        verts, faces, seed)
    selected_a = grow_faces(seed_a, face_adjacency(fa), args.rings)
    selected_b = grow_faces(seed_b, face_adjacency(fb), args.rings)
    selected = np.concatenate((selected_a, selected_b))
    candidate = build_candidate(
        verts, faces, selected, projected, center, axis_u, axis_v,
        args.merge_distance, args.steiner_spacing)
    unfilled_verts = candidate["verts"].copy()
    unfilled_faces = candidate["faces"].copy()
    unfilled_patch_mask = candidate["patch_mask"].copy()
    candidate = fill_small_boundary_loops(
        candidate, center, axis_u, axis_v,
        args.holefill_max_edges, args.holefill_max_area)

    origin = (np.asarray(args.cube_origin, dtype=np.float64)
              if args.cube_origin else cube_origin_from_path(args.pred))
    pred_volume = tifffile.imread(args.pred)
    raw_volume = tifffile.imread(args.raw)
    if pred_volume.ndim != 3 or raw_volume.shape != pred_volume.shape:
        raise ValueError("PRED and RAW must be equal-shaped 3D TIFF volumes")
    surface = prediction_surface(pred_volume)
    surface_tree = cKDTree(surface)
    fitted, snap_stats = fit_to_prediction(
        candidate["verts"], candidate["faces"], candidate["pin_mask"],
        surface, origin, args.snap_alpha, args.snap_distance, args.snap_outer)

    split_topology = mesh_topology(verts, faces)
    merge_topology = mesh_topology(fitted, candidate["faces"])
    split_quality = quality_metrics(verts, faces, selected)
    merge_quality = quality_metrics(fitted, candidate["faces"],
                                    candidate["patch_mask"])
    pred_split = prediction_fit(verts, faces, selected, surface_tree, origin)
    pred_merge = prediction_fit(fitted, candidate["faces"],
                                candidate["patch_mask"], surface_tree, origin)
    split_raw, _ = interface_diagnostics(verts, faces, selected, raw_volume, origin)
    merge_raw, merge_energy = interface_diagnostics(
        fitted, candidate["faces"], candidate["patch_mask"], raw_volume, origin)
    if split_raw.get("pairs", 0) == 0 or merge_raw.get("pairs", 0) == 0:
        raise RuntimeError("could not form a RAW seam-energy comparison")

    planar_area = projected_patch_metrics(
        projected, faces, selected, candidate["points_2d"],
        candidate["patch_faces_2d"])
    initial_area, initial_normals, _ = face_geometry(candidate["verts"],
                                                     candidate["faces"])
    _, fitted_normals, _ = face_geometry(fitted, candidate["faces"])
    patch = candidate["patch_mask"]
    reference_sign = np.sign(np.median(initial_normals[patch] @ normal))
    if reference_sign == 0:
        reference_sign = 1.0
    fitted_folds = int(np.count_nonzero(reference_sign *
                                        (fitted_normals[patch] @ normal) <= 0))
    expected_loops = max(1, split_topology["boundary_loops"] - 1)

    hard_gates = {
        "one_intended_component": merge_topology["face_components"] == 1,
        "edge_manifold": merge_topology["nonmanifold_edges"] == 0,
        "consistently_oriented":
            merge_topology["same_direction_interior_edges"] == 0,
        "simple_boundary": merge_topology["bad_boundary_vertices"] == 0,
        "expected_boundary_loops": merge_topology["boundary_loops"] == expected_loops,
        "consistent_interface_orientation":
            candidate["interface_edges_same_direction"] == 0,
        "no_fitted_patch_folds": fitted_folds == 0,
        "triangle_quality_floor":
            merge_quality["p01"] >= max(0.10, 0.5 * split_quality["p01"]),
        "prediction_p95_safety":
            pred_merge["p95"] <= pred_split["p95"] + args.pred_p95_margin,
    }
    split_quilt = split_raw["quilt_value_discontinuity"]
    merge_quilt = merge_raw["quilt_value_discontinuity"]
    evidence = {
        "prediction_mean_no_worse":
            pred_merge["mean"] <= pred_split["mean"] + args.pred_score_margin,
        "prediction_p90_no_worse":
            pred_merge["p90"] <= pred_split["p90"] + args.pred_score_margin,
        "raw_quilt_no_worse":
            merge_quilt["mean"] <= split_quilt["mean"] + args.quilt_margin,
        "interface_bend_sane":
            merge_raw["normal_angle_degrees"]["p90"] <= args.bend_p90_max,
    }
    eligible = all(hard_gates.values())
    if not eligible:
        verdict = "split_preferred_failed_merge_gates"
    elif all(evidence.values()):
        verdict = "merge_provisionally_preferred_visual_review_required"
    else:
        verdict = "review_evidence_ambiguous"

    construction = {key: value for key, value in candidate.items()
                    if key not in ("verts", "faces", "patch_mask", "pin_mask",
                                   "points_2d", "patch_faces_2d")}
    metrics = {
        "inputs": {
            "sheet_a": str(args.sheet_a), "sheet_b": str(args.sheet_b),
            "pred": str(args.pred), "raw": str(args.raw),
            "cube_origin_zyx": origin.tolist(),
        },
        "parameters": {
            "seed_mode": seed_mode, "contact_distance": args.contact_distance,
            "seed_vertices": args.seed_vertices, "rings": args.rings,
            "merge_distance": args.merge_distance,
            "steiner_spacing": args.steiner_spacing,
            "holefill_max_edges": args.holefill_max_edges,
            "holefill_max_area": args.holefill_max_area,
        },
        "contact": {
            "seed_faces_a": int(np.count_nonzero(seed_a)),
            "seed_faces_b": int(np.count_nonzero(seed_b)),
            "selected_faces_a": int(np.count_nonzero(selected_a)),
            "selected_faces_b": int(np.count_nonzero(selected_b)),
            "common_plane_center": center.tolist(),
            "common_plane_u": axis_u.tolist(),
            "common_plane_v": axis_v.tolist(),
            "common_plane_normal": normal.tolist(),
            "contact_plane_rms_thickness": thickness,
        },
        "construction": construction,
        "projected_area": planar_area,
        "topology": {"split": split_topology, "merge": merge_topology},
        "quality": {"split_patch": split_quality, "merge_patch": merge_quality},
        "prediction_snap": snap_stats,
        "prediction_fit": {"split_patch": pred_split, "merge_patch": pred_merge},
        "raw_interface": {"split": split_raw, "merge": merge_raw},
        "fitted_patch_folds": fitted_folds,
        "decision": {
            "eligible": eligible, "hard_gates": hard_gates,
            "evidence_gates": evidence, "verdict": verdict,
            "calibration_only": True,
        },
    }

    red, purple = np.asarray((0.95, 0.18, 0.18)), np.asarray((0.58, 0.22, 0.88))
    vertex_colors = np.vstack((np.tile(red, (len(va), 1)),
                               np.tile(purple, (len(vb), 1))))
    write_obj(args.out_dir / "00_split_pair.obj", verts, faces, vertex_colors,
              "current split hypothesis: red=sheet A, purple=sheet B")
    seed_colors = np.tile((0.65, 0.65, 0.65), (len(faces), 1))
    seed_colors[seed] = (1.0, 0.05, 0.05)
    write_face_colors(args.out_dir / "01_contact_seed.obj", verts, faces,
                      seed_colors, "red=detected close contact")
    selected_colors = np.vstack((np.tile(red, (len(fa), 1)),
                                 np.tile(purple, (len(fb), 1))))[selected]
    write_face_colors(args.out_dir / "02_overlap_plus_rings.obj", verts,
                      faces[selected], selected_colors,
                      "overlap/contact plus face rings replaced by merge probe")
    unfilled_colors = np.tile((0.62, 0.62, 0.62), (len(unfilled_faces), 1))
    unfilled_colors[unfilled_patch_mask] = (0.10, 0.90, 0.30)
    write_face_colors(args.out_dir / "03_merge_planar_unfilled.obj", unfilled_verts,
                      unfilled_faces, unfilled_colors,
                      "green=common-plane Delaunay before bounded hole fill")
    candidate_colors = np.tile((0.62, 0.62, 0.62), (len(candidate["faces"]), 1))
    candidate_colors[patch] = (0.10, 0.90, 0.30)
    candidate_vertex_colors = np.tile((0.62, 0.62, 0.62),
                                      (len(candidate["verts"]), 1))
    candidate_vertex_colors[np.unique(candidate["faces"][patch])] = (0.10, 0.90, 0.30)
    write_face_colors(args.out_dir / "04_merge_planar_holefilled.obj", candidate["verts"],
                      candidate["faces"], candidate_colors,
                      "green=common-plane Delaunay after bounded micro-hole fill")
    write_obj(args.out_dir / "04_merge_planar_holefilled_topology.obj",
              candidate["verts"], candidate["faces"], candidate_vertex_colors,
              "topology-preserving common-plane merge; green=healed patch")
    write_face_colors(args.out_dir / "05_merge_pred_fitted.obj", fitted,
                      candidate["faces"], candidate_colors,
                      "green=healed patch fitted to nnU-Net prediction")
    write_obj(args.out_dir / "05_merge_pred_fitted_topology.obj", fitted,
              candidate["faces"], candidate_vertex_colors,
              "topology-preserving fitted merge; green=healed patch")
    write_obj(args.out_dir / "06_merge_patch_only.obj", fitted,
              candidate["faces"][patch], comment="fitted healed patch only")
    write_face_colors(args.out_dir / "07_merge_quilt_energy.obj", fitted,
                      candidate["faces"], heat_colors(merge_energy),
                      "green-to-red=low-to-high RAW seam discontinuity; gray=unscored")
    displacement = np.linalg.norm(fitted - candidate["verts"], axis=1)
    disp_face = displacement[candidate["faces"]].mean(axis=1)
    disp_face[~patch] = math.nan
    write_face_colors(args.out_dir / "08_merge_pred_displacement.obj", fitted,
                      candidate["faces"], heat_colors(disp_face),
                      "green-to-red=plane-to-pred displacement; gray=fixed exterior")

    split_face_colors = np.vstack((np.tile(red, (len(fa), 1)),
                                   np.tile(purple, (len(fb), 1))))
    write_obj(args.out_dir / "A_split_output.obj", verts, faces, vertex_colors,
              "current split hypothesis")
    write_obj(args.out_dir / "B_merge_output.obj", fitted, candidate["faces"],
              candidate_vertex_colors,
              "healed, hole-filled, prediction-fitted merge hypothesis")
    write_side_by_side(args.out_dir / "C_split_vs_merge_side_by_side.obj",
                       verts, faces, split_face_colors,
                       fitted, candidate["faces"], candidate_colors)

    with (args.out_dir / "metrics.json").open("w", encoding="utf-8", newline="\n") as dst:
        json.dump(metrics, dst, indent=2, sort_keys=True)
        dst.write("\n")
    render_readme(args.out_dir / "README.md", metrics)
    print(json.dumps({"verdict": verdict, "eligible": eligible,
                      "out_dir": str(args.out_dir),
                      "failed_hard_gates": [k for k, v in hard_gates.items() if not v]},
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
