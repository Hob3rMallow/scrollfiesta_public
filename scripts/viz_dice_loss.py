#!/usr/bin/env python3
"""
viz_dice_loss.py — Surface Dice spatial analysis.

Shows WHERE surface dice is lost by the pipeline vs raw nnUNet prediction.
Uses the exact same border extraction as the surface_distance package
(2x2x2 neighbourhood encoding, bbox crop+pad, EDT tolerance matching).

Usage: python scripts/viz_dice_loss.py <cube_id> [--base-dir DIR] [--tau 2.0]

Outputs:
  output/{cube_id}/{cube_id}_viz_dice_mip.png    — GT/pred-nn/pipeline surface match MIPs
  output/{cube_id}/{cube_id}_viz_dice_diff.png    — GT surface colored by loss/gain
  output/{cube_id}/{cube_id}_viz_dice_slices.png  — 8 equidistant Z-slices with overlay
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
from surface_distance import lookup_tables

sys.path.insert(0, os.path.dirname(__file__))
from viz_common import (
    load_volumes, extract_borders_pair, compute_surface_dice,
    make_rgb_mip_row, save_png, DEFAULT_BASE_DIR
)


def extract_all_borders(gt_bin, nn_bin, pipe_bin):
    """Extract borders for GT-vs-NN and GT-vs-Pipeline pairs.

    Each pair shares its own bounding box (matching the surface_distance package).
    Returns border masks and area maps in the cropped+padded coordinate space.
    """
    kernel = lookup_tables.ENCODE_NEIGHBOURHOOD_3D_KERNEL
    area_table = lookup_tables.create_table_neighbour_code_to_surface_area((1.0, 1.0, 1.0))

    results = {}
    for name, vol_a, vol_b in [("nn", gt_bin, nn_bin), ("pipe", gt_bin, pipe_bin)]:
        ba, aa, bb, ab, bbox_min, crop_shape = extract_borders_pair(vol_a, vol_b)
        results[name] = {
            "border_gt": ba, "area_gt": aa,
            "border_pred": bb, "area_pred": ab,
            "bbox_min": bbox_min, "crop_shape": crop_shape,
        }
    return results


def classify_surface_voxels(border_gt, area_gt, border_pred, area_pred, tau=2.0):
    """Classify each GT border voxel as matched/unmatched by pred, and vice versa.

    Returns:
        gt_matched:   bool array — GT border voxels within tau of pred border
        gt_unmatched: bool array — GT border voxels farther than tau
        pred_matched: bool array
        pred_unmatched: bool array
        gt_dists: float array — distance of each GT border voxel to nearest pred border
        pred_dists: float array — distance of each pred border voxel to nearest GT border
    """
    distmap_pred = ndimage.distance_transform_edt(~border_pred, sampling=(1.0, 1.0, 1.0))
    distmap_gt = ndimage.distance_transform_edt(~border_gt, sampling=(1.0, 1.0, 1.0))

    gt_dists_full = np.full(border_gt.shape, np.inf)
    gt_dists_full[border_gt] = distmap_pred[border_gt]

    pred_dists_full = np.full(border_pred.shape, np.inf)
    pred_dists_full[border_pred] = distmap_gt[border_pred]

    gt_matched = border_gt & (gt_dists_full <= tau)
    gt_unmatched = border_gt & (gt_dists_full > tau)
    pred_matched = border_pred & (pred_dists_full <= tau)
    pred_unmatched = border_pred & (pred_dists_full > tau)

    return gt_matched, gt_unmatched, pred_matched, pred_unmatched, gt_dists_full, pred_dists_full


def surface_mip_rgb(border, matched, unmatched, axis=0, matched_color=(0, 200, 0), unmatched_color=(200, 0, 0)):
    """Create RGB MIP of surface voxels colored by match status.

    Projections collapse along the given axis. Matched (green) and unmatched (red).
    """
    m_proj = np.any(matched, axis=axis)
    u_proj = np.any(unmatched, axis=axis)

    rgb = np.zeros((*m_proj.shape, 3), dtype=np.uint8)
    rgb[u_proj] = unmatched_color
    rgb[m_proj] = matched_color  # matched overwrites unmatched in projection
    return rgb


def main():
    parser = argparse.ArgumentParser(description="Surface Dice loss visualization")
    parser.add_argument("cube_id", type=str)
    parser.add_argument("--base-dir", default=DEFAULT_BASE_DIR)
    parser.add_argument("--tau", type=float, default=2.0, help="Surface dice tolerance")
    args = parser.parse_args()

    cube_id = args.cube_id
    base_dir = args.base_dir
    tau = args.tau

    gt_bin, nn_bin, pipe_raw, pipe_bin, shape = load_volumes(cube_id, base_dir)

    # ─── Extract borders for both pairs ───
    print("\nExtracting borders (GT vs pred-nn)...")
    nn_bgt, nn_agt, nn_bpred, nn_apred, nn_bbox, nn_cs = extract_borders_pair(gt_bin, nn_bin)

    print("Extracting borders (GT vs pipeline)...")
    p_bgt, p_agt, p_bpred, p_apred, p_bbox, p_cs = extract_borders_pair(gt_bin, pipe_bin)

    # ─── Compute surface dice ───
    print(f"\nComputing Surface Dice (tau={tau})...")
    nn_sd, nn_ov_gt, nn_ov_pred, nn_total = compute_surface_dice(nn_bgt, nn_agt, nn_bpred, nn_apred, tau)
    p_sd, p_ov_gt, p_ov_pred, p_total = compute_surface_dice(p_bgt, p_agt, p_bpred, p_apred, tau)

    # ─── Classify border voxels ───
    print("Classifying GT-vs-NN surface matches...")
    nn_gt_m, nn_gt_u, nn_pred_m, nn_pred_u, nn_gt_d, nn_pred_d = \
        classify_surface_voxels(nn_bgt, nn_agt, nn_bpred, nn_apred, tau)

    print("Classifying GT-vs-Pipeline surface matches...")
    p_gt_m, p_gt_u, p_pred_m, p_pred_u, p_gt_d, p_pred_d = \
        classify_surface_voxels(p_bgt, p_agt, p_bpred, p_apred, tau)

    # ─── Console report ───
    print(f"\n{'='*60}")
    print(f"SURFACE DICE REPORT — {cube_id} (tau={tau})")
    print(f"{'='*60}")
    print(f"             {'GT-vs-NN':>14s}  {'GT-vs-Pipe':>14s}")
    print(f"  GT border voxels:  {int(nn_bgt.sum()):>14,}  {int(p_bgt.sum()):>14,}")
    print(f"  Pred border voxels:{int(nn_bpred.sum()):>14,}  {int(p_bpred.sum()):>14,}")
    print(f"  GT area:           {nn_agt[nn_bgt].sum():>14.1f}  {p_agt[p_bgt].sum():>14.1f}")
    print(f"  Pred area:         {nn_apred[nn_bpred].sum():>14.1f}  {p_apred[p_bpred].sum():>14.1f}")
    print(f"  GT overlap area:   {nn_ov_gt:>14.1f}  {p_ov_gt:>14.1f}")
    print(f"  Pred overlap area: {nn_ov_pred:>14.1f}  {p_ov_pred:>14.1f}")
    print(f"  Total area:        {nn_total:>14.1f}  {p_total:>14.1f}")
    print(f"  Surface Dice:      {nn_sd:>14.6f}  {p_sd:>14.6f}")
    print(f"  Delta:             {'':>14s}  {p_sd - nn_sd:>+14.6f}")
    print(f"{'='*60}")

    # GT surface loss analysis (in GT-vs-NN space vs GT-vs-Pipe space)
    # These are in DIFFERENT coordinate spaces, so we analyze them separately
    nn_gt_matched_count = int(nn_gt_m.sum())
    nn_gt_unmatched_count = int(nn_gt_u.sum())
    p_gt_matched_count = int(p_gt_m.sum())
    p_gt_unmatched_count = int(p_gt_u.sum())

    print(f"  GT matched by NN:       {nn_gt_matched_count:>10,}  ({100*nn_gt_matched_count/max(1,nn_gt_matched_count+nn_gt_unmatched_count):.1f}%)")
    print(f"  GT unmatched by NN:     {nn_gt_unmatched_count:>10,}")
    print(f"  GT matched by Pipeline: {p_gt_matched_count:>10,}  ({100*p_gt_matched_count/max(1,p_gt_matched_count+p_gt_unmatched_count):.1f}%)")
    print(f"  GT unmatched by Pipe:   {p_gt_unmatched_count:>10,}")
    print(f"{'='*60}")

    # ─── Figure 1: Surface match MIPs (3 rows: GT surface, NN surface, Pipeline surface) ───
    fig1, axes = plt.subplots(3, 1, figsize=(16, 14), facecolor="black",
                               gridspec_kw={"height_ratios": [1, 1, 1]})

    # Row 1: GT surface (white) — using nn pair's coordinate space
    gt_mips = []
    for axis in [0, 1, 2]:
        proj = np.any(nn_bgt, axis=axis)
        rgb = np.zeros((*proj.shape, 3), dtype=np.uint8)
        rgb[proj] = (200, 200, 200)
        gt_mips.append(rgb)
    gt_row = make_rgb_mip_row(gt_mips[0], gt_mips[1], gt_mips[2])

    # Row 2: NN prediction surface (green=matched, red=unmatched)
    nn_mips = [surface_mip_rgb(nn_bpred, nn_pred_m, nn_pred_u, axis=a) for a in [0, 1, 2]]
    nn_row = make_rgb_mip_row(nn_mips[0], nn_mips[1], nn_mips[2])

    # Row 3: Pipeline prediction surface (green=matched, red=unmatched)
    pipe_mips = [surface_mip_rgb(p_bpred, p_pred_m, p_pred_u, axis=a) for a in [0, 1, 2]]
    pipe_row = make_rgb_mip_row(pipe_mips[0], pipe_mips[1], pipe_mips[2])

    for ax, row, title in zip(axes,
                               [gt_row, nn_row, pipe_row],
                               ["GT surface", f"Pred-NN surface (SDice={nn_sd:.4f})",
                                f"Pipeline surface (SDice={p_sd:.4f})"]):
        ax.imshow(row)
        ax.set_axis_off()
        ax.set_title(title, color="white", fontsize=11, pad=5)

    legend1 = [
        Patch(facecolor=(200/255, 200/255, 200/255), label="GT surface"),
        Patch(facecolor=(0, 200/255, 0), label=f"Matched (d <= {tau})"),
        Patch(facecolor=(200/255, 0, 0), label=f"Unmatched (d > {tau})"),
    ]
    axes[2].legend(handles=legend1, loc="lower right", fontsize=8,
                   facecolor="black", edgecolor="gray", labelcolor="white")

    plt.tight_layout()
    save_png(fig1, cube_id, "viz_dice_mip", base_dir)
    plt.close(fig1)

    # ─── Figure 2: GT surface colored by NN-vs-Pipeline loss ───
    # We need to map GT border voxels to original volume coords to compare
    # NN and Pipeline match status for the SAME GT voxels.
    # Since the two pairs may have different bboxes, we work in original volume coords.

    print("\nComputing GT-centric loss map in original coordinates...")

    # Recompute borders on the full GT vs NN and GT vs Pipe, but use a shared bbox
    # For the differential, we need GT border positions in the same coordinate space
    # Use the full volume for the GT border, then EDT from each pred
    kernel = lookup_tables.ENCODE_NEIGHBOURHOOD_3D_KERNEL
    area_table = lookup_tables.create_table_neighbour_code_to_surface_area((1.0, 1.0, 1.0))

    # Compute GT borders on the full volume (with union bbox = GT itself)
    # Actually, let's do it more carefully: union of all three
    all_union = (gt_bin | nn_bin | pipe_bin).astype(bool)

    proj_0 = np.amax(all_union, axis=(1, 2))
    idx_0 = np.nonzero(proj_0)[0]
    proj_1 = np.amax(all_union, axis=(0, 2))
    idx_1 = np.nonzero(proj_1)[0]
    proj_2 = np.amax(all_union, axis=(0, 1))
    idx_2 = np.nonzero(proj_2)[0]

    bbox_min = np.array([idx_0.min(), idx_1.min(), idx_2.min()])
    bbox_max = np.array([idx_0.max(), idx_1.max(), idx_2.max()])
    padded_shape = (bbox_max - bbox_min) + 2

    def crop_pad(vol):
        cm = np.zeros(padded_shape, np.uint8)
        cm[0:-1, 0:-1, 0:-1] = vol[bbox_min[0]:bbox_max[0]+1,
                                     bbox_min[1]:bbox_max[1]+1,
                                     bbox_min[2]:bbox_max[2]+1]
        return cm

    cm_gt = crop_pad(gt_bin.astype(np.uint8))
    cm_nn = crop_pad(nn_bin.astype(np.uint8))
    cm_pipe = crop_pad(pipe_bin.astype(np.uint8))

    code_gt = ndimage.correlate(cm_gt, kernel, mode="constant", cval=0)
    code_nn = ndimage.correlate(cm_nn, kernel, mode="constant", cval=0)
    code_pipe = ndimage.correlate(cm_pipe, kernel, mode="constant", cval=0)

    border_gt = (code_gt != 0) & (code_gt != 255)
    border_nn = (code_nn != 0) & (code_nn != 255)
    border_pipe = (code_pipe != 0) & (code_pipe != 255)

    # EDT from each prediction's border to GT border
    distmap_nn = ndimage.distance_transform_edt(~border_nn, sampling=(1.0, 1.0, 1.0))
    distmap_pipe = ndimage.distance_transform_edt(~border_pipe, sampling=(1.0, 1.0, 1.0))

    gt_dist_to_nn = np.full(border_gt.shape, np.inf)
    gt_dist_to_nn[border_gt] = distmap_nn[border_gt]

    gt_dist_to_pipe = np.full(border_gt.shape, np.inf)
    gt_dist_to_pipe[border_gt] = distmap_pipe[border_gt]

    # Classify each GT border voxel
    gt_matched_nn = border_gt & (gt_dist_to_nn <= tau)
    gt_matched_pipe = border_gt & (gt_dist_to_pipe <= tau)

    gt_lost = gt_matched_nn & ~gt_matched_pipe      # was matched by NN, now unmatched
    gt_gained = ~gt_matched_nn & gt_matched_pipe     # was unmatched by NN, now matched
    gt_unchanged_m = gt_matched_nn & gt_matched_pipe # still matched
    gt_unchanged_u = ~gt_matched_nn & ~gt_matched_pipe  # still unmatched

    n_gt_lost = int(gt_lost.sum())
    n_gt_gained = int(gt_gained.sum())
    n_gt_unch_m = int(gt_unchanged_m.sum())
    n_gt_unch_u = int(gt_unchanged_u.sum())

    print(f"  GT surfels lost by pipeline:    {n_gt_lost:,}")
    print(f"  GT surfels gained by pipeline:  {n_gt_gained:,}")
    print(f"  GT surfels unchanged (matched): {n_gt_unch_m:,}")
    print(f"  GT surfels unchanged (unmatch): {n_gt_unch_u:,}")

    # MIP of the differential
    diff_mips = []
    for axis in [0, 1, 2]:
        lost_proj = np.any(gt_lost, axis=axis)
        gain_proj = np.any(gt_gained, axis=axis)
        unch_m_proj = np.any(gt_unchanged_m, axis=axis)
        unch_u_proj = np.any(gt_unchanged_u & border_gt, axis=axis)

        rgb = np.zeros((*lost_proj.shape, 3), dtype=np.uint8)
        rgb[unch_u_proj] = (80, 80, 80)    # dark gray = still unmatched
        rgb[unch_m_proj] = (140, 140, 140)  # light gray = still matched
        rgb[gain_proj] = (0, 220, 0)        # green = gained by pipeline
        rgb[lost_proj] = (220, 0, 0)        # red = lost by pipeline
        diff_mips.append(rgb)

    diff_row = make_rgb_mip_row(diff_mips[0], diff_mips[1], diff_mips[2])

    fig2, ax = plt.subplots(figsize=(16, 5), facecolor="black")
    ax.imshow(diff_row)
    ax.set_axis_off()
    ax.set_title(f"{cube_id} — GT Surface: Pipeline vs Pred-NN Dice Differential\n"
                 f"Lost: {n_gt_lost:,} | Gained: {n_gt_gained:,} | "
                 f"Still matched: {n_gt_unch_m:,} | Still unmatched: {n_gt_unch_u:,}",
                 color="white", fontsize=11, pad=10)

    legend2 = [
        Patch(facecolor=(220/255, 0, 0), label="Lost by pipeline (was matched)"),
        Patch(facecolor=(0, 220/255, 0), label="Gained by pipeline"),
        Patch(facecolor=(140/255, 140/255, 140/255), label="Still matched"),
        Patch(facecolor=(80/255, 80/255, 80/255), label="Still unmatched"),
    ]
    ax.legend(handles=legend2, loc="lower right", fontsize=8,
              facecolor="black", edgecolor="gray", labelcolor="white")

    save_png(fig2, cube_id, "viz_dice_diff", base_dir)
    plt.close(fig2)

    # ─── Figure 3: Z-slices with overlay ───
    D_crop = padded_shape[0]
    n_slices = 8
    slice_indices = np.linspace(0, D_crop - 1, n_slices, dtype=int)

    fig3, axes = plt.subplots(2, 4, figsize=(20, 10), facecolor="black")
    axes_flat = axes.flatten()

    for i, zi in enumerate(slice_indices):
        ax = axes_flat[i]
        ax.set_facecolor("black")

        # Create RGB overlay of this Z-slice
        h, w = border_gt.shape[1], border_gt.shape[2]
        rgb = np.zeros((h, w, 3), dtype=np.uint8)

        # Background: GT border in gray
        gt_sl = border_gt[zi]
        rgb[gt_sl] = (100, 100, 100)

        # Overlay classification
        rgb[gt_unchanged_m[zi]] = (140, 140, 140)
        rgb[(gt_unchanged_u & border_gt)[zi]] = (60, 60, 60)
        rgb[gt_gained[zi]] = (0, 220, 0)
        rgb[gt_lost[zi]] = (220, 0, 0)

        # Also show pred borders faintly
        rgb[border_nn[zi] & ~border_gt[zi]] = (0, 80, 0)     # NN-only border in dark green
        rgb[border_pipe[zi] & ~border_gt[zi]] = (0, 0, 120)  # Pipe-only border in dark blue

        ax.imshow(rgb)
        orig_z = zi + bbox_min[0] if zi < D_crop - 1 else -1
        ax.set_title(f"Z={orig_z}", color="white", fontsize=9)
        ax.set_axis_off()

    plt.suptitle(f"{cube_id} — Surface Dice Loss: 8 Z-slices\n"
                 f"Red=lost, Green=gained, Gray=unchanged",
                 color="white", fontsize=12)
    plt.tight_layout()
    save_png(fig3, cube_id, "viz_dice_slices", base_dir)
    plt.close(fig3)

    print(f"\nDone. Output in output/{cube_id}/")


if __name__ == "__main__":
    main()
