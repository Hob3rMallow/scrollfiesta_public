#!/usr/bin/env python3
"""Compare face normal orientations between two OBJ files.

Usage:
    python check_face_orientation.py pre.obj post.obj
    python check_face_orientation.py --batch <pre_dir> <post_dir> <prefix_pre> <prefix_post> <num_components>

Diagnoses whether QEM mesh simplification produces inverted faces.
"""
import sys
import math
import os


def read_obj(path):
    """Read an OBJ file, return (verts, faces). Faces are 0-based index triples."""
    verts = []
    faces = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('v '):
                parts = line.split()
                verts.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith('f '):
                parts = line.split()
                # OBJ faces can be "f v" or "f v/vt" or "f v/vt/vn" or "f v//vn"
                indices = []
                for p in parts[1:]:
                    idx = int(p.split('/')[0]) - 1  # OBJ is 1-based
                    indices.append(idx)
                if len(indices) >= 3:
                    # Triangulate if polygon (fan triangulation)
                    for i in range(1, len(indices) - 1):
                        faces.append((indices[0], indices[i], indices[i + 1]))
    return verts, faces


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def length(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def compute_face_normals(verts, faces):
    """Compute per-face normals (unnormalized = area-weighted)."""
    normals = []
    for f in faces:
        v0 = verts[f[0]]
        v1 = verts[f[1]]
        v2 = verts[f[2]]
        e1 = sub(v1, v0)
        e2 = sub(v2, v0)
        n = cross(e1, e2)
        normals.append(n)
    return normals


def analyze_mesh(path):
    """Analyze a single OBJ file. Returns dict with stats."""
    verts, faces = read_obj(path)
    nv = len(verts)
    nf = len(faces)

    if nf == 0:
        return {
            'path': path,
            'nv': nv,
            'nf': nf,
            'mean_normal': (0, 0, 0),
            'inverted_count': 0,
            'inverted_pct': 0.0,
            'degenerate_count': 0,
            'total_area': 0.0,
        }

    normals = compute_face_normals(verts, faces)

    # Area-weighted mean normal
    mean_n = [0.0, 0.0, 0.0]
    total_area = 0.0
    degenerate = 0
    for n in normals:
        area = 0.5 * length(n)
        if area < 1e-12:
            degenerate += 1
            continue
        total_area += area
        mean_n[0] += n[0]
        mean_n[1] += n[1]
        mean_n[2] += n[2]

    mean_len = length(mean_n)
    if mean_len > 1e-12:
        mean_n = (mean_n[0] / mean_len, mean_n[1] / mean_len, mean_n[2] / mean_len)
    else:
        mean_n = (0, 0, 0)

    # Count inverted faces (dot with mean normal < 0)
    inverted = 0
    for n in normals:
        nl = length(n)
        if nl < 1e-12:
            continue
        unit_n = (n[0] / nl, n[1] / nl, n[2] / nl)
        d = dot(unit_n, mean_n)
        if d < 0:
            inverted += 1

    non_degen = nf - degenerate
    inverted_pct = (inverted / non_degen * 100.0) if non_degen > 0 else 0.0

    return {
        'path': path,
        'nv': nv,
        'nf': nf,
        'mean_normal': mean_n,
        'inverted_count': inverted,
        'inverted_pct': inverted_pct,
        'degenerate_count': degenerate,
        'total_area': total_area,
    }


def print_analysis(label, stats):
    print(f"\n{'=' * 60}")
    print(f"  {label}")
    print(f"  {stats['path']}")
    print(f"{'=' * 60}")
    print(f"  Vertices:         {stats['nv']:>10,}")
    print(f"  Faces:            {stats['nf']:>10,}")
    print(f"  Degenerate faces: {stats['degenerate_count']:>10,}")
    mn = stats['mean_normal']
    print(f"  Mean normal:      ({mn[0]:+.4f}, {mn[1]:+.4f}, {mn[2]:+.4f})")
    print(f"  Total area:       {stats['total_area']:>14,.2f}")
    print(f"  Inverted faces:   {stats['inverted_count']:>10,}  ({stats['inverted_pct']:.2f}%)")


def compare_two(path_pre, path_post):
    stats_pre = analyze_mesh(path_pre)
    stats_post = analyze_mesh(path_post)

    print_analysis("PRE-SIMPLIFY", stats_pre)
    print_analysis("POST-SIMPLIFY", stats_post)

    # Check if mean normals are consistent
    d = dot(stats_pre['mean_normal'], stats_post['mean_normal'])
    print(f"\n  Mean normal dot product (pre vs post): {d:+.4f}")
    if d < 0:
        print("  WARNING: Mean normals are OPPOSITE - global orientation flip detected!")
    elif d < 0.5:
        print("  WARNING: Mean normals diverge significantly.")
    else:
        print("  OK: Mean normals are consistent.")

    # Summary
    delta_inv = stats_post['inverted_pct'] - stats_pre['inverted_pct']
    print(f"\n  Inversion change: {delta_inv:+.2f}% points")
    if delta_inv > 5.0:
        print("  ALERT: QEM is significantly increasing face inversions!")
    elif delta_inv > 1.0:
        print("  NOTE: QEM is moderately increasing face inversions.")
    else:
        print("  OK: QEM inversion impact is small.")
    print()


def batch_mode(pre_dir, post_dir, prefix_pre, prefix_post, num_components):
    """Run on all components and produce a summary table."""
    results = []
    for i in range(1, num_components + 1):
        comp_str = f"comp{i:03d}"
        pre_path = os.path.join(pre_dir, f"{prefix_pre}_{comp_str}.obj")
        post_path = os.path.join(post_dir, f"{prefix_post}_{comp_str}.obj")

        if not os.path.exists(pre_path):
            continue
        if not os.path.exists(post_path):
            continue

        stats_pre = analyze_mesh(pre_path)
        stats_post = analyze_mesh(post_path)
        d = dot(stats_pre['mean_normal'], stats_post['mean_normal'])

        results.append({
            'comp': i,
            'pre_nv': stats_pre['nv'],
            'pre_nf': stats_pre['nf'],
            'post_nv': stats_post['nv'],
            'post_nf': stats_post['nf'],
            'pre_inv_pct': stats_pre['inverted_pct'],
            'post_inv_pct': stats_post['inverted_pct'],
            'normal_dot': d,
        })

    # Print summary table
    print()
    print(f"{'Comp':>5}  {'Pre NV':>10}  {'Pre NF':>10}  {'Post NV':>10}  {'Post NF':>10}  "
          f"{'Pre Inv%':>9}  {'Post Inv%':>10}  {'Delta':>7}  {'NormDot':>8}  {'Status':>8}")
    print("-" * 110)

    for r in results:
        delta = r['post_inv_pct'] - r['pre_inv_pct']
        if r['normal_dot'] < 0:
            status = "FLIP!"
        elif delta > 5.0:
            status = "BAD"
        elif delta > 1.0:
            status = "WARN"
        else:
            status = "OK"

        print(f"{r['comp']:>5}  {r['pre_nv']:>10,}  {r['pre_nf']:>10,}  {r['post_nv']:>10,}  "
              f"{r['post_nf']:>10,}  {r['pre_inv_pct']:>8.2f}%  {r['post_inv_pct']:>9.2f}%  "
              f"{delta:>+6.2f}%  {r['normal_dot']:>+7.4f}  {status:>8}")

    print("-" * 110)

    # Summary stats
    if results:
        avg_pre = sum(r['pre_inv_pct'] for r in results) / len(results)
        avg_post = sum(r['post_inv_pct'] for r in results) / len(results)
        flips = sum(1 for r in results if r['normal_dot'] < 0)
        print(f"\n  Components analyzed: {len(results)}")
        print(f"  Average pre-simplify inversion:  {avg_pre:.2f}%")
        print(f"  Average post-simplify inversion: {avg_post:.2f}%")
        print(f"  Global orientation flips:        {flips}")
    print()


def main():
    if len(sys.argv) == 3:
        compare_two(sys.argv[1], sys.argv[2])
    elif len(sys.argv) >= 3 and sys.argv[1] == '--batch':
        if len(sys.argv) != 7:
            print("Usage: python check_face_orientation.py --batch <pre_dir> <post_dir> <prefix_pre> <prefix_post> <num_components>")
            sys.exit(1)
        batch_mode(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6]))
    else:
        print("Usage:")
        print("  python check_face_orientation.py <pre.obj> <post.obj>")
        print("  python check_face_orientation.py --batch <pre_dir> <post_dir> <prefix_pre> <prefix_post> <num_components>")
        sys.exit(1)


if __name__ == '__main__':
    main()
