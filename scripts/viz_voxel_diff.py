#!/usr/bin/env python3
"""
viz_voxel_diff.py — TP/FP/FN voxel comparison between pred-nn and pipeline.

Shows WHERE material was lost/gained by the pipeline relative to raw nnUNet prediction.

Usage: python scripts/viz_voxel_diff.py <cube_id> [--base-dir DIR]

Outputs:
  output/{cube_id}/{cube_id}_viz_diff_mip.png    — MIP classification grids
  output/{cube_id}/{cube_id}_viz_diff_profile.png — Per-Z-slice FG count profiles
"""
import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(__file__))
from viz_common import (
    load_volumes, mip_any, classify_mips, make_rgb_mip_row,
    bool_to_rgb, save_png, DEFAULT_BASE_DIR
)


def main():
    parser = argparse.ArgumentParser(description="Voxel diff visualization")
    parser.add_argument("cube_id", type=str)
    parser.add_argument("--base-dir", default=DEFAULT_BASE_DIR)
    args = parser.parse_args()

    cube_id = args.cube_id
    base_dir = args.base_dir

    # Load volumes
    gt_bin, nn_bin, pipe_raw, pipe_bin, shape = load_volumes(cube_id, base_dir)

    # ─── Compute TP/FP/FN for both predictions ───
    nn_tp = int((gt_bin & nn_bin).sum())
    nn_fp = int((~gt_bin & nn_bin).sum())
    nn_fn = int((gt_bin & ~nn_bin).sum())
    nn_prec = nn_tp / max(nn_tp + nn_fp, 1)
    nn_rec = nn_tp / max(nn_tp + nn_fn, 1)

    pipe_tp = int((gt_bin & pipe_bin).sum())
    pipe_fp = int((~gt_bin & pipe_bin).sum())
    pipe_fn = int((gt_bin & ~pipe_bin).sum())
    pipe_prec = pipe_tp / max(pipe_tp + pipe_fp, 1)
    pipe_rec = pipe_tp / max(pipe_tp + pipe_fn, 1)

    # Differential analysis
    lost_tp = gt_bin & nn_bin & ~pipe_bin      # was TP, now FN
    gained_tp = gt_bin & ~nn_bin & pipe_bin    # was FN, now TP
    removed_fp = ~gt_bin & nn_bin & ~pipe_bin  # was FP, now TN
    added_fp = ~gt_bin & ~nn_bin & pipe_bin    # was TN, now FP

    n_lost_tp = int(lost_tp.sum())
    n_gained_tp = int(gained_tp.sum())
    n_removed_fp = int(removed_fp.sum())
    n_added_fp = int(added_fp.sum())

    # Console report
    print(f"\n{'='*60}")
    print(f"VOXEL DIFF REPORT — {cube_id}")
    print(f"{'='*60}")
    print(f"             {'Pred-NN':>12s}  {'Pipeline':>12s}  {'Delta':>12s}")
    print(f"  TP:        {nn_tp:>12,}  {pipe_tp:>12,}  {pipe_tp - nn_tp:>+12,}")
    print(f"  FP:        {nn_fp:>12,}  {pipe_fp:>12,}  {pipe_fp - nn_fp:>+12,}")
    print(f"  FN:        {nn_fn:>12,}  {pipe_fn:>12,}  {pipe_fn - nn_fn:>+12,}")
    print(f"  Precision: {nn_prec:>12.4f}  {pipe_prec:>12.4f}  {pipe_prec - nn_prec:>+12.4f}")
    print(f"  Recall:    {nn_rec:>12.4f}  {pipe_rec:>12.4f}  {pipe_rec - nn_rec:>+12.4f}")
    print(f"  FG total:  {int(nn_bin.sum()):>12,}  {int(pipe_bin.sum()):>12,}  {int(pipe_bin.sum()) - int(nn_bin.sum()):>+12,}")
    print(f"{'='*60}")
    print(f"  Differential (pipeline vs pred-nn):")
    print(f"    Lost TP (TP->FN):      {n_lost_tp:>12,}")
    print(f"    Gained TP (FN->TP):    {n_gained_tp:>12,}")
    print(f"    Removed FP (FP->TN):   {n_removed_fp:>12,}")
    print(f"    Added FP (TN->FP):     {n_added_fp:>12,}")
    print(f"    Net TP change:         {n_gained_tp - n_lost_tp:>+12,}")
    print(f"    Net FP change:         {n_added_fp - n_removed_fp:>+12,}")
    print(f"{'='*60}")

    # ─── Figure 1: MIP classification grid ───
    # Row 1: pred-nn classification (TP/FP/FN)
    nn_mip_z, nn_mip_y, nn_mip_x = classify_mips(gt_bin, nn_bin)
    nn_row = make_rgb_mip_row(nn_mip_z, nn_mip_y, nn_mip_x)

    # Row 2: pipeline classification (TP/FP/FN)
    pipe_mip_z, pipe_mip_y, pipe_mip_x = classify_mips(gt_bin, pipe_bin)
    pipe_row = make_rgb_mip_row(pipe_mip_z, pipe_mip_y, pipe_mip_x)

    # Row 3: differential (red=lost TP, green=gained TP, cyan=removed FP, magenta=added FP)
    lost_z, lost_y, lost_x = mip_any(lost_tp)
    gain_z, gain_y, gain_x = mip_any(gained_tp)
    rmfp_z, rmfp_y, rmfp_x = mip_any(removed_fp)
    adfp_z, adfp_y, adfp_x = mip_any(added_fp)

    diff_mips = []
    for lm, gm, rm, am in [(lost_z, gain_z, rmfp_z, adfp_z),
                            (lost_y, gain_y, rmfp_y, adfp_y),
                            (lost_x, gain_x, rmfp_x, adfp_x)]:
        rgb = np.zeros((*lm.shape, 3), dtype=np.uint8)
        rgb[rm] = (0, 150, 150)      # cyan = removed FP (good)
        rgb[am] = (180, 0, 180)      # magenta = added FP (bad)
        rgb[gm] = (0, 220, 0)        # green = gained TP (good)
        rgb[lm] = (220, 0, 0)        # red = lost TP (bad)
        diff_mips.append(rgb)

    diff_row = make_rgb_mip_row(diff_mips[0], diff_mips[1], diff_mips[2])

    # Compose figure
    gap_h = 8
    total_h = nn_row.shape[0] + pipe_row.shape[0] + diff_row.shape[0] + 2 * gap_h
    max_w = max(nn_row.shape[1], pipe_row.shape[1], diff_row.shape[1])

    canvas = np.zeros((total_h, max_w, 3), dtype=np.uint8)
    y = 0
    canvas[y:y + nn_row.shape[0], :nn_row.shape[1]] = nn_row
    y += nn_row.shape[0] + gap_h
    canvas[y:y + pipe_row.shape[0], :pipe_row.shape[1]] = pipe_row
    y += pipe_row.shape[0] + gap_h
    canvas[y:y + diff_row.shape[0], :diff_row.shape[1]] = diff_row

    fig, ax = plt.subplots(figsize=(16, 12), facecolor="black")
    ax.imshow(canvas)
    ax.set_axis_off()
    ax.set_title(f"{cube_id} — Voxel Classification MIP\n"
                 f"Row 1: Pred-NN | Row 2: Pipeline | Row 3: Differential\n"
                 f"Each: Z-proj | Y-proj | X-proj",
                 color="white", fontsize=11, pad=10)

    legend_elements = [
        Patch(facecolor=(0, 200/255, 0), label="TP"),
        Patch(facecolor=(200/255, 0, 0), label="FP"),
        Patch(facecolor=(50/255, 50/255, 200/255), label="FN"),
        Patch(facecolor=(220/255, 0, 0), label="Lost TP (bad)"),
        Patch(facecolor=(0, 220/255, 0), label="Gained TP (good)"),
        Patch(facecolor=(0, 150/255, 150/255), label="Removed FP (good)"),
        Patch(facecolor=(180/255, 0, 180/255), label="Added FP (bad)"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8,
              facecolor="black", edgecolor="gray", labelcolor="white")

    save_png(fig, cube_id, "viz_diff_mip", base_dir)
    plt.close(fig)

    # ─── Figure 2: Per-Z-slice profile ───
    D = shape[0]
    z_slices = np.arange(D)

    gt_per_z = np.array([gt_bin[z].sum() for z in range(D)])
    nn_per_z = np.array([nn_bin[z].sum() for z in range(D)])
    pipe_per_z = np.array([pipe_bin[z].sum() for z in range(D)])

    lost_per_z = np.array([lost_tp[z].sum() for z in range(D)])
    gained_per_z = np.array([gained_tp[z].sum() for z in range(D)])

    fig2, axes = plt.subplots(2, 1, figsize=(14, 8), facecolor="black",
                               gridspec_kw={"height_ratios": [2, 1]})

    # Top: FG counts
    ax1 = axes[0]
    ax1.set_facecolor("black")
    ax1.plot(z_slices, gt_per_z, color="white", linewidth=1.5, label="GT", alpha=0.9)
    ax1.plot(z_slices, nn_per_z, color="lime", linewidth=1.0, label="Pred-NN", alpha=0.8)
    ax1.plot(z_slices, pipe_per_z, color="cyan", linewidth=1.0, label="Pipeline", alpha=0.8)
    ax1.fill_between(z_slices, pipe_per_z, nn_per_z,
                     where=(nn_per_z > pipe_per_z),
                     color="red", alpha=0.3, label="Lost FG")
    ax1.fill_between(z_slices, pipe_per_z, nn_per_z,
                     where=(pipe_per_z > nn_per_z),
                     color="green", alpha=0.3, label="Gained FG")
    ax1.set_ylabel("Foreground voxels per slice", color="white")
    ax1.set_title(f"{cube_id} — Per-Z-slice FG counts", color="white", fontsize=12)
    ax1.legend(fontsize=8, facecolor="black", edgecolor="gray", labelcolor="white")
    ax1.tick_params(colors="white")
    ax1.spines["bottom"].set_color("gray")
    ax1.spines["left"].set_color("gray")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Bottom: Lost/gained TP per slice
    ax2 = axes[1]
    ax2.set_facecolor("black")
    ax2.fill_between(z_slices, 0, -lost_per_z, color="red", alpha=0.6, label="Lost TP")
    ax2.fill_between(z_slices, 0, gained_per_z, color="green", alpha=0.6, label="Gained TP")
    ax2.axhline(0, color="gray", linewidth=0.5)
    ax2.set_xlabel("Z slice", color="white")
    ax2.set_ylabel("Voxel count", color="white")
    ax2.set_title("TP change per slice (positive=gained, negative=lost)", color="white", fontsize=10)
    ax2.legend(fontsize=8, facecolor="black", edgecolor="gray", labelcolor="white")
    ax2.tick_params(colors="white")
    ax2.spines["bottom"].set_color("gray")
    ax2.spines["left"].set_color("gray")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    save_png(fig2, cube_id, "viz_diff_profile", base_dir)
    plt.close(fig2)

    print(f"\nDone. Output in output/{cube_id}/")


if __name__ == "__main__":
    main()
