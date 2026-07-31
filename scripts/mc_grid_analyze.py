#!/usr/bin/env python3
"""
mc_grid_analyze.py -- post-process an mc_grid_dump OBJ and report how
well adjacent cubes' MC output coincides at the seams.

For EVERY vertex (not just halo), quantize the source-space position
to 1e-4 voxel and group by (qz, qy, qx). The cluster size tells you
how many cubes contributed a vert at that exact source position.

Expected if MC is perfectly deterministic across cubes:
  - Verts deep inside ONE cube's owned region: cluster size 1.
  - Verts at a 2-cube face seam: cluster size 2 (both cubes' MC produces
    a vert at the shared source coord).
  - Verts at a 4-cube edge seam: cluster size 4.
  - Verts at an 8-cube corner: cluster size 8.

If MC is non-deterministic, most overlap verts stay at cluster size 1,
revealing that adjacent cubes produce DIFFERENT vert sets in regions
where they should agree by construction.

Usage: python mc_grid_analyze.py <combined.obj> [--break-by-halo]
"""
import sys
from collections import Counter


def main():
    if len(sys.argv) < 2:
        print("Usage: mc_grid_analyze.py <combined.obj> [--break-by-halo]",
              file=sys.stderr)
        sys.exit(2)

    path = sys.argv[1]
    break_by_halo = len(sys.argv) > 2 and sys.argv[2] == '--break-by-halo'

    positions = []
    halo_flags = []  # parallel: True if vert was colored as halo (red)
    n_total = 0
    n_halo = 0
    with open(path, 'r') as f:
        for line in f:
            if not line.startswith('v '):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            n_total += 1
            z, y, x = float(parts[1]), float(parts[2]), float(parts[3])
            qz = int(round(z * 10000))
            qy = int(round(y * 10000))
            qx = int(round(x * 10000))
            positions.append((qz, qy, qx))
            if len(parts) >= 7:
                r, g, b = float(parts[4]), float(parts[5]), float(parts[6])
                is_halo = r > 0.9 and g < 0.3 and b < 0.3
            else:
                is_halo = False
            halo_flags.append(is_halo)
            if is_halo:
                n_halo += 1

    print(f"total verts: {n_total}")
    print(f"halo verts:  {n_halo} ({100.0 * n_halo / n_total:.2f}%)")

    if not positions:
        return

    # Bucket by source position; count cube contributions.
    pos_counts = Counter(positions)
    unique_positions = len(pos_counts)
    print(f"\nunique source positions: {unique_positions}")
    print(f"of which appear in >=2 cubes: "
          f"{sum(1 for c in pos_counts.values() if c >= 2)} "
          f"({100.0 * sum(1 for c in pos_counts.values() if c >= 2) / unique_positions:.2f}%)")

    cluster_hist = Counter(pos_counts.values())
    print("\ncluster size histogram (ALL verts):")
    print(f"  size  unique_positions  contributing_verts  %of_unique  %of_verts")
    for size in sorted(cluster_hist.keys()):
        n_pos = cluster_hist[size]
        n_verts = size * n_pos
        print(f"  {size:>4}  {n_pos:>16}  {n_verts:>18}  "
              f"{100.0 * n_pos / unique_positions:>9.2f}%  "
              f"{100.0 * n_verts / n_total:>8.2f}%")

    if break_by_halo:
        # Same but breaking down by halo status: was the position contributed
        # by at least one halo-marked vert?
        pos_has_halo = {}
        for p, h in zip(positions, halo_flags):
            pos_has_halo[p] = pos_has_halo.get(p, False) or h

        halo_pos = [p for p in pos_counts if pos_has_halo[p]]
        owned_only_pos = [p for p in pos_counts if not pos_has_halo[p]]
        print(f"\nbreakdown by halo touch:")
        print(f"  positions touched by halo vert (cube boundary region): "
              f"{len(halo_pos)}")
        print(f"  positions never touched by halo vert (cube interior):  "
              f"{len(owned_only_pos)}")

        halo_cluster_hist = Counter(pos_counts[p] for p in halo_pos)
        print(f"\n  halo-touched cluster size histogram:")
        for size in sorted(halo_cluster_hist.keys()):
            n_pos = halo_cluster_hist[size]
            print(f"    {size:>4}: {n_pos:>10} positions "
                  f"({100.0 * n_pos / len(halo_pos):.2f}%)")


if __name__ == '__main__':
    main()
