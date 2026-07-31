#!/usr/bin/env python3
"""
viz_voi_loss.py — VOI component analysis.

Shows WHERE VOI split/merge happens by analyzing the contingency table between
GT components and prediction components. Pipeline component labels come from
the uint16 output TIF directly (no re-labeling needed).

Usage: python scripts/viz_voi_loss.py <cube_id> [--base-dir DIR]

Outputs:
  output/{cube_id}/{cube_id}_viz_voi_components.png — Component coloring MIPs
  output/{cube_id}/{cube_id}_viz_voi_splits.png     — Split severity heatmap
  output/{cube_id}/{cube_id}_viz_voi_merges.png     — Merge severity heatmap
  output/{cube_id}/{cube_id}_viz_voi_diff.png       — Voxel difference pipeline vs pred-nn
  output/{cube_id}/{cube_id}_viz_voi_bars.png       — Per-component VOI contribution bars
"""
import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import ndimage
from skimage.metrics import variation_of_information
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from viz_common import (
    load_volumes, mip_any, mip_max, make_rgb_mip_row, save_png, DEFAULT_BASE_DIR
)


def cc3d_26_scipy(binary):
    """26-connected component labeling using scipy (matching compare_metrics.py)."""
    struct = ndimage.generate_binary_structure(3, 3)
    labels, n = ndimage.label(binary.astype(np.uint8), structure=struct)
    return labels.astype(np.int32), n


def random_colors(n, seed=42):
    """Generate n distinct random colors for component visualization."""
    rng = np.random.RandomState(seed)
    colors = np.zeros((n + 1, 3), dtype=np.uint8)  # index 0 = background = black
    for i in range(1, n + 1):
        # HSV-like distribution for good contrast
        hue = (i * 137.508) % 360  # golden angle
        sat = 0.7 + rng.random() * 0.3
        val = 0.7 + rng.random() * 0.3
        colors[i] = hsv_to_rgb(hue / 360, sat, val)
    return colors


def hsv_to_rgb(h, s, v):
    """Convert HSV to RGB uint8."""
    import colorsys
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    return (int(r * 255), int(g * 255), int(b * 255))


def label_to_rgb_mip(labels, colors, axis=0):
    """Create RGB MIP of labeled volume using max-label projection.

    Uses max-projection of label values, then maps to colors.
    """
    proj = np.max(labels, axis=axis)
    h, w = proj.shape
    rgb = np.zeros((h, w, 3), dtype=np.uint8)
    for label_id in range(1, len(colors)):
        mask = (proj == label_id)
        if mask.any():
            rgb[mask] = colors[label_id]
    return rgb


def compute_contingency(gt_labels, pred_labels):
    """Compute contingency table between GT and pred labels.

    Returns:
        table: dict of (gt_id, pred_id) -> voxel_count (FG only)
        gt_sizes: dict of gt_id -> total voxel count
        pred_sizes: dict of pred_id -> total voxel count
    """
    # Only consider foreground voxels
    fg_mask = (gt_labels > 0) | (pred_labels > 0)

    gt_fg = gt_labels[fg_mask]
    pred_fg = pred_labels[fg_mask]

    table = defaultdict(int)
    gt_sizes = defaultdict(int)
    pred_sizes = defaultdict(int)

    for g, p in zip(gt_fg.ravel(), pred_fg.ravel()):
        table[(int(g), int(p))] += 1
        gt_sizes[int(g)] += 1
        pred_sizes[int(p)] += 1

    return dict(table), dict(gt_sizes), dict(pred_sizes)


def compute_contingency_fast(gt_labels, pred_labels):
    """Fast contingency table using numpy unique pairs."""
    fg_mask = (gt_labels > 0) | (pred_labels > 0)
    gt_fg = gt_labels[fg_mask].ravel()
    pred_fg = pred_labels[fg_mask].ravel()

    # Encode pair as single int
    max_pred = int(pred_fg.max()) + 1
    encoded = gt_fg.astype(np.int64) * max_pred + pred_fg.astype(np.int64)
    unique_pairs, counts = np.unique(encoded, return_counts=True)

    table = {}
    gt_sizes = defaultdict(int)
    pred_sizes = defaultdict(int)

    for enc, cnt in zip(unique_pairs, counts):
        g = int(enc // max_pred)
        p = int(enc % max_pred)
        table[(g, p)] = int(cnt)
        gt_sizes[g] += int(cnt)
        pred_sizes[p] += int(cnt)

    return table, dict(gt_sizes), dict(pred_sizes)


def analyze_splits_merges(table, gt_sizes, pred_sizes):
    """Analyze split and merge counts from contingency table.

    Returns:
        gt_split_count: dict gt_id -> number of pred components overlapping it
        pred_merge_count: dict pred_id -> number of GT components it straddles
        gt_dominant_pred: dict gt_id -> pred_id with most overlap
        pred_dominant_gt: dict pred_id -> gt_id with most overlap
    """
    # For each GT component, count how many pred components overlap it
    gt_to_preds = defaultdict(dict)  # gt_id -> {pred_id: count}
    pred_to_gts = defaultdict(dict)  # pred_id -> {gt_id: count}

    for (g, p), cnt in table.items():
        if g > 0 and p > 0:
            gt_to_preds[g][p] = cnt
            pred_to_gts[p][g] = cnt

    gt_split_count = {}
    gt_dominant_pred = {}
    for g, preds in gt_to_preds.items():
        gt_split_count[g] = len(preds)
        gt_dominant_pred[g] = max(preds, key=preds.get)

    pred_merge_count = {}
    pred_dominant_gt = {}
    for p, gts in pred_to_gts.items():
        pred_merge_count[p] = len(gts)
        pred_dominant_gt[p] = max(gts, key=gts.get)

    # Also track GT components that are completely missed (no pred overlap)
    all_gt_ids = set(g for g in gt_sizes if g > 0)
    for g in all_gt_ids:
        if g not in gt_split_count:
            gt_split_count[g] = 0
            gt_dominant_pred[g] = 0

    all_pred_ids = set(p for p in pred_sizes if p > 0)
    for p in all_pred_ids:
        if p not in pred_merge_count:
            pred_merge_count[p] = 0
            pred_dominant_gt[p] = 0

    return gt_split_count, pred_merge_count, gt_dominant_pred, pred_dominant_gt, gt_to_preds, pred_to_gts


def severity_heatmap_mip(labels, severity_map, axis=0):
    """Create RGB MIP where each component is colored by its severity.

    severity_map: dict label_id -> severity count (1, 2, 3+)
    Colors: 1=green, 2=yellow, 3+=red
    """
    proj = np.max(labels, axis=axis)
    h, w = proj.shape
    rgb = np.zeros((h, w, 3), dtype=np.uint8)

    for label_id, sev in severity_map.items():
        mask = (proj == label_id)
        if not mask.any():
            continue
        if sev <= 1:
            color = (0, 200, 0)      # green = clean (1:1 match)
        elif sev == 2:
            color = (220, 220, 0)    # yellow = moderate (split/merge into 2)
        else:
            color = (220, 0, 0)      # red = severe (3+ way split/merge)
        rgb[mask] = color

    return rgb


def compute_voi_cropped(gt_bin, pred_bin, alpha=0.3):
    """Compute VOI with bbox crop matching compare_metrics.py."""
    union = gt_bin | pred_bin
    if not union.any():
        return 0.0, 0.0, 0.0, 1.0

    crop_ratio = 0.7
    coords = np.argwhere(union)
    bb_min = coords.min(axis=0)
    bb_max = coords.max(axis=0) + 1
    bbox_vol = int(np.prod(bb_max - bb_min))
    full_vol = int(np.prod(gt_bin.shape))

    if bbox_vol <= crop_ratio * full_vol:
        slc = tuple(slice(lo, hi) for lo, hi in zip(bb_min, bb_max))
        gt_use = gt_bin[slc]
        pr_use = pred_bin[slc]
    else:
        gt_use = gt_bin
        pr_use = pred_bin

    gt_lab, _ = cc3d_26_scipy(gt_use)
    pr_lab, _ = cc3d_26_scipy(pr_use)

    m = (gt_lab > 0) | (pr_lab > 0)
    a = gt_lab[m]
    b = pr_lab[m]

    voi_split, voi_merge = variation_of_information(a, b)
    voi_total = float(voi_split + voi_merge)
    voi_score = 1.0 / (1.0 + alpha * voi_total)

    return float(voi_split), float(voi_merge), voi_total, voi_score


def main():
    parser = argparse.ArgumentParser(description="VOI component analysis visualization")
    parser.add_argument("cube_id", type=str)
    parser.add_argument("--base-dir", default=DEFAULT_BASE_DIR)
    args = parser.parse_args()

    cube_id = args.cube_id
    base_dir = args.base_dir

    gt_bin, nn_bin, pipe_raw, pipe_bin, shape = load_volumes(cube_id, base_dir)

    # ─── Label connected components ───
    print("\nLabeling GT components (26-conn)...")
    gt_labels, gt_n = cc3d_26_scipy(gt_bin)
    print(f"  GT components: {gt_n}")

    print("Labeling pred-nn components (26-conn)...")
    nn_labels, nn_n = cc3d_26_scipy(nn_bin)
    print(f"  Pred-NN components: {nn_n}")

    # If pipeline output has multi-valued component labels, use them directly;
    # otherwise fall back to 26-conn CC labeling (matching the competition metric)
    pipe_unique = np.unique(pipe_raw[pipe_raw > 0])
    if len(pipe_unique) > 1 and pipe_unique.max() > 1:
        print(f"  Pipeline has {len(pipe_unique)} component labels in TIF, using directly")
        pipe_labels = pipe_raw.astype(np.int32)
        pipe_n = int(pipe_labels.max())
    else:
        print("  Pipeline output is binary, labeling with 26-conn CC...")
        pipe_labels, pipe_n = cc3d_26_scipy(pipe_bin)

    print(f"  Pipeline components: {pipe_n}")

    # ─── Compute VOI ───
    print("\nComputing VOI (GT vs pred-nn)...")
    nn_split, nn_merge, nn_voi, nn_vscore = compute_voi_cropped(gt_bin, nn_bin)
    print(f"  VOI: {nn_voi:.4f} (split={nn_split:.4f}, merge={nn_merge:.4f}), score={nn_vscore:.4f}")

    print("Computing VOI (GT vs pipeline)...")
    p_split, p_merge, p_voi, p_vscore = compute_voi_cropped(gt_bin, pipe_bin)
    print(f"  VOI: {p_voi:.4f} (split={p_split:.4f}, merge={p_merge:.4f}), score={p_vscore:.4f}")

    # ─── Contingency tables ───
    print("\nComputing contingency tables...")
    nn_table, nn_gt_sz, nn_pred_sz = compute_contingency_fast(gt_labels, nn_labels)
    p_table, p_gt_sz, p_pred_sz = compute_contingency_fast(gt_labels, pipe_labels)

    nn_gt_sc, nn_pred_mc, nn_gt_dp, nn_pred_dg, nn_g2p, nn_p2g = \
        analyze_splits_merges(nn_table, nn_gt_sz, nn_pred_sz)

    p_gt_sc, p_pred_mc, p_gt_dp, p_pred_dg, p_g2p, p_p2g = \
        analyze_splits_merges(p_table, p_gt_sz, p_pred_sz)

    # ─── Console report ───
    print(f"\n{'='*60}")
    print(f"VOI REPORT — {cube_id}")
    print(f"{'='*60}")
    print(f"             {'GT-vs-NN':>14s}  {'GT-vs-Pipe':>14s}  {'Delta':>14s}")
    print(f"  VOI total: {nn_voi:>14.4f}  {p_voi:>14.4f}  {p_voi - nn_voi:>+14.4f}")
    print(f"  VOI split: {nn_split:>14.4f}  {p_split:>14.4f}  {p_split - nn_split:>+14.4f}")
    print(f"  VOI merge: {nn_merge:>14.4f}  {p_merge:>14.4f}  {p_merge - nn_merge:>+14.4f}")
    print(f"  VOI score: {nn_vscore:>14.4f}  {p_vscore:>14.4f}  {p_vscore - nn_vscore:>+14.4f}")
    print(f"  #GT comps: {gt_n:>14d}")
    print(f"  #Pred:     {nn_n:>14d}  {pipe_n:>14d}")
    print(f"{'='*60}")

    # Top split offenders
    print(f"\n  GT components — split severity (how many pred comps overlap each GT comp):")
    print(f"  {'GT ID':>6s}  {'Size':>10s}  {'NN splits':>10s}  {'Pipe splits':>12s}")
    for g in sorted(nn_gt_sc.keys()):
        nn_s = nn_gt_sc.get(g, 0)
        p_s = p_gt_sc.get(g, 0)
        sz = nn_gt_sz.get(g, 0)
        if nn_s > 1 or p_s > 1:
            print(f"  {g:>6d}  {sz:>10,}  {nn_s:>10d}  {p_s:>12d}")

    # Top merge offenders
    print(f"\n  Pred components — merge severity (how many GT comps each pred straddles):")
    print(f"  {'--- Pred-NN ---':}")
    for p in sorted(nn_pred_mc.keys()):
        mc = nn_pred_mc[p]
        if mc > 1:
            sz = nn_pred_sz.get(p, 0)
            gts = nn_p2g.get(p, {})
            gt_str = ", ".join(f"GT{g}:{c}" for g, c in sorted(gts.items(), key=lambda x: -x[1])[:5])
            print(f"    Pred {p:>4d}: size={sz:>8,}, merges={mc}, overlaps: {gt_str}")

    print(f"  {'--- Pipeline ---':}")
    for p in sorted(p_pred_mc.keys()):
        mc = p_pred_mc[p]
        if mc > 1:
            sz = p_pred_sz.get(p, 0)
            gts = p_p2g.get(p, {})
            gt_str = ", ".join(f"GT{g}:{c}" for g, c in sorted(gts.items(), key=lambda x: -x[1])[:5])
            print(f"    Pipe {p:>4d}: size={sz:>8,}, merges={mc}, overlaps: {gt_str}")

    print(f"{'='*60}")

    # ─── Generate colors ───
    gt_colors = random_colors(gt_n, seed=42)
    nn_colors = random_colors(nn_n, seed=123)
    pipe_colors_raw = random_colors(pipe_n, seed=456)

    # Color pred components by their dominant GT match
    nn_gt_colors = np.zeros((nn_n + 1, 3), dtype=np.uint8)
    for p, g in nn_pred_dg.items():
        if 0 < p <= nn_n and 0 < g <= gt_n:
            nn_gt_colors[p] = gt_colors[g]
        elif 0 < p <= nn_n:
            nn_gt_colors[p] = (100, 100, 100)  # no GT match -> gray

    pipe_gt_colors = np.zeros((pipe_n + 1, 3), dtype=np.uint8)
    for p, g in p_pred_dg.items():
        if 0 < p <= pipe_n and 0 < g <= gt_n:
            pipe_gt_colors[p] = gt_colors[g]
        elif 0 < p <= pipe_n:
            pipe_gt_colors[p] = (100, 100, 100)

    # ─── Figure 1: Component coloring MIPs ───
    fig1, axes = plt.subplots(3, 1, figsize=(16, 14), facecolor="black",
                               gridspec_kw={"height_ratios": [1, 1, 1]})

    # Row 1: GT components
    gt_mips = [label_to_rgb_mip(gt_labels, gt_colors, axis=a) for a in [0, 1, 2]]
    gt_row = make_rgb_mip_row(gt_mips[0], gt_mips[1], gt_mips[2])

    # Row 2: Pred-NN components colored by GT match
    nn_mips = [label_to_rgb_mip(nn_labels, nn_gt_colors, axis=a) for a in [0, 1, 2]]
    nn_row = make_rgb_mip_row(nn_mips[0], nn_mips[1], nn_mips[2])

    # Row 3: Pipeline components colored by GT match
    pipe_mips = [label_to_rgb_mip(pipe_labels, pipe_gt_colors, axis=a) for a in [0, 1, 2]]
    pipe_row = make_rgb_mip_row(pipe_mips[0], pipe_mips[1], pipe_mips[2])

    for ax, row, title in zip(axes,
                               [gt_row, nn_row, pipe_row],
                               [f"GT ({gt_n} components)",
                                f"Pred-NN ({nn_n} comps, colored by GT match)",
                                f"Pipeline ({pipe_n} comps, colored by GT match)"]):
        ax.imshow(row)
        ax.set_axis_off()
        ax.set_title(title, color="white", fontsize=11, pad=5)

    plt.tight_layout()
    save_png(fig1, cube_id, "viz_voi_components", base_dir)
    plt.close(fig1)

    # ─── Figure 2: Split severity heatmap ───
    fig2, axes = plt.subplots(2, 1, figsize=(16, 10), facecolor="black")

    nn_split_mips = [severity_heatmap_mip(gt_labels, nn_gt_sc, axis=a) for a in [0, 1, 2]]
    nn_split_row = make_rgb_mip_row(nn_split_mips[0], nn_split_mips[1], nn_split_mips[2])

    p_split_mips = [severity_heatmap_mip(gt_labels, p_gt_sc, axis=a) for a in [0, 1, 2]]
    p_split_row = make_rgb_mip_row(p_split_mips[0], p_split_mips[1], p_split_mips[2])

    for ax, row, title in zip(axes,
                               [nn_split_row, p_split_row],
                               ["GT comps colored by split count (Pred-NN)",
                                "GT comps colored by split count (Pipeline)"]):
        ax.imshow(row)
        ax.set_axis_off()
        ax.set_title(title, color="white", fontsize=11, pad=5)

    legend_sev = [
        Patch(facecolor=(0, 200/255, 0), label="1 (clean match)"),
        Patch(facecolor=(220/255, 220/255, 0), label="2 (split into 2)"),
        Patch(facecolor=(220/255, 0, 0), label="3+ (severe split)"),
    ]
    axes[1].legend(handles=legend_sev, loc="lower right", fontsize=8,
                   facecolor="black", edgecolor="gray", labelcolor="white")

    plt.tight_layout()
    save_png(fig2, cube_id, "viz_voi_splits", base_dir)
    plt.close(fig2)

    # ─── Figure 3: Merge severity heatmap ───
    fig3, axes = plt.subplots(2, 1, figsize=(16, 10), facecolor="black")

    nn_merge_mips = [severity_heatmap_mip(nn_labels, nn_pred_mc, axis=a) for a in [0, 1, 2]]
    nn_merge_row = make_rgb_mip_row(nn_merge_mips[0], nn_merge_mips[1], nn_merge_mips[2])

    p_merge_mips = [severity_heatmap_mip(pipe_labels, p_pred_mc, axis=a) for a in [0, 1, 2]]
    p_merge_row = make_rgb_mip_row(p_merge_mips[0], p_merge_mips[1], p_merge_mips[2])

    for ax, row, title in zip(axes,
                               [nn_merge_row, p_merge_row],
                               ["Pred comps colored by merge count (Pred-NN)",
                                "Pred comps colored by merge count (Pipeline)"]):
        ax.imshow(row)
        ax.set_axis_off()
        ax.set_title(title, color="white", fontsize=11, pad=5)

    axes[1].legend(handles=legend_sev, loc="lower right", fontsize=8,
                   facecolor="black", edgecolor="gray", labelcolor="white")

    plt.tight_layout()
    save_png(fig3, cube_id, "viz_voi_merges", base_dir)
    plt.close(fig3)

    # ─── Figure 4: Pipeline vs pred-nn voxel difference ───
    pipe_only = pipe_bin & ~nn_bin   # cyan
    nn_only = nn_bin & ~pipe_bin     # magenta
    both = nn_bin & pipe_bin          # white

    diff_mips = []
    for axis in [0, 1, 2]:
        po_proj = np.any(pipe_only, axis=axis)
        no_proj = np.any(nn_only, axis=axis)
        bo_proj = np.any(both, axis=axis)
        rgb = np.zeros((*po_proj.shape, 3), dtype=np.uint8)
        rgb[bo_proj] = (120, 120, 120)   # gray = both
        rgb[po_proj] = (0, 200, 200)     # cyan = pipeline only
        rgb[no_proj] = (200, 0, 200)     # magenta = pred-nn only
        diff_mips.append(rgb)

    diff_row = make_rgb_mip_row(diff_mips[0], diff_mips[1], diff_mips[2])

    fig4, ax = plt.subplots(figsize=(16, 5), facecolor="black")
    ax.imshow(diff_row)
    ax.set_axis_off()
    ax.set_title(f"{cube_id} — Pipeline vs Pred-NN Voxel Difference", color="white", fontsize=11, pad=10)

    legend4 = [
        Patch(facecolor=(120/255, 120/255, 120/255), label="Both"),
        Patch(facecolor=(0, 200/255, 200/255), label="Pipeline only"),
        Patch(facecolor=(200/255, 0, 200/255), label="Pred-NN only"),
    ]
    ax.legend(handles=legend4, loc="lower right", fontsize=8,
              facecolor="black", edgecolor="gray", labelcolor="white")

    save_png(fig4, cube_id, "viz_voi_diff", base_dir)
    plt.close(fig4)

    # ─── Figure 5: Per-component VOI contribution bars ───
    # Show how each GT component contributes to the VOI difference
    gt_ids = sorted([g for g in nn_gt_sz if g > 0])

    fig5, axes = plt.subplots(1, 2, figsize=(16, 6), facecolor="black")

    # Split counts bar chart
    ax_s = axes[0]
    ax_s.set_facecolor("black")
    x = np.arange(len(gt_ids))
    nn_splits = [nn_gt_sc.get(g, 0) for g in gt_ids]
    p_splits = [p_gt_sc.get(g, 0) for g in gt_ids]

    width = 0.35
    ax_s.bar(x - width/2, nn_splits, width, color="lime", alpha=0.7, label="Pred-NN")
    ax_s.bar(x + width/2, p_splits, width, color="cyan", alpha=0.7, label="Pipeline")
    ax_s.set_xticks(x)
    ax_s.set_xticklabels([f"GT{g}" for g in gt_ids], rotation=45, fontsize=8, color="white")
    ax_s.set_ylabel("# pred components overlapping", color="white")
    ax_s.set_title("Split count per GT component", color="white", fontsize=11)
    ax_s.axhline(1, color="gray", linewidth=0.5, linestyle="--")
    ax_s.legend(fontsize=8, facecolor="black", edgecolor="gray", labelcolor="white")
    ax_s.tick_params(colors="white")
    ax_s.spines["bottom"].set_color("gray")
    ax_s.spines["left"].set_color("gray")
    ax_s.spines["top"].set_visible(False)
    ax_s.spines["right"].set_visible(False)

    # GT component sizes
    ax_sz = axes[1]
    ax_sz.set_facecolor("black")
    gt_szs = [nn_gt_sz.get(g, 0) for g in gt_ids]
    bar_colors = []
    for g in gt_ids:
        c = gt_colors[g] if g <= gt_n else (128, 128, 128)
        bar_colors.append(tuple(v/255 for v in c))

    ax_sz.bar(x, gt_szs, color=bar_colors, alpha=0.8)
    ax_sz.set_xticks(x)
    ax_sz.set_xticklabels([f"GT{g}" for g in gt_ids], rotation=45, fontsize=8, color="white")
    ax_sz.set_ylabel("Voxel count", color="white")
    ax_sz.set_title("GT component sizes", color="white", fontsize=11)
    ax_sz.tick_params(colors="white")
    ax_sz.spines["bottom"].set_color("gray")
    ax_sz.spines["left"].set_color("gray")
    ax_sz.spines["top"].set_visible(False)
    ax_sz.spines["right"].set_visible(False)

    plt.tight_layout()
    save_png(fig5, cube_id, "viz_voi_bars", base_dir)
    plt.close(fig5)

    print(f"\nDone. Output in output/{cube_id}/")


if __name__ == "__main__":
    main()
