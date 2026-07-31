#!/usr/bin/env python3
"""
viz_common.py — Shared utilities for metric loss visualization scripts.

Provides volume loading, MIP projection, border extraction (matching the
surface_distance package exactly), and figure saving.
"""
import os
import numpy as np
import tifffile
from scipy import ndimage
from surface_distance import lookup_tables

IGNORE_LABEL = 2
DEFAULT_BASE_DIR = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))


def load_volumes(cube_id, base_dir=None):
    """Load GT, pred-nn, and pipeline TIFs; apply ignore mask; return binarized volumes.

    Returns:
        gt_bin:   bool array — ground truth (ignore-masked)
        nn_bin:   bool array — nnUNet prediction (ignore-masked)
        pipe_raw: uint16 array — pipeline output (component labels, ignore-masked)
        pipe_bin: bool array — pipeline binarized (ignore-masked)
        shape:    tuple — volume shape (D, H, W)
    """
    if base_dir is None:
        base_dir = DEFAULT_BASE_DIR

    gt_path = os.path.join(base_dir, "gt", f"{cube_id}.tif")
    nn_path = os.path.join(base_dir, "nnunet-preds", f"{cube_id}.tif")
    pipe_path = os.path.join(base_dir, "output", str(cube_id), f"{cube_id}_output.tif")

    print(f"Loading GT:       {gt_path}")
    print(f"Loading pred-nn:  {nn_path}")
    print(f"Loading pipeline: {pipe_path}")

    gt_raw = tifffile.imread(gt_path)
    nn_raw = tifffile.imread(nn_path)
    pipe_raw = tifffile.imread(pipe_path).astype(np.uint16)

    assert gt_raw.shape == nn_raw.shape == pipe_raw.shape, \
        f"Shape mismatch: GT={gt_raw.shape}, NN={nn_raw.shape}, Pipe={pipe_raw.shape}"

    shape = gt_raw.shape
    print(f"Volume shape: {shape}")

    # Apply ignore mask (label=2 in GT)
    ign = (gt_raw == IGNORE_LABEL)
    n_ign = int(ign.sum())
    print(f"Ignored voxels (label={IGNORE_LABEL}): {n_ign}")

    gt_eval = np.where(ign, 0, gt_raw)
    nn_eval = np.where(ign, 0, nn_raw)
    pipe_masked = np.where(ign, 0, pipe_raw).astype(np.uint16)

    gt_bin = (gt_eval != 0).astype(bool)
    nn_bin = (nn_eval != 0).astype(bool)
    pipe_bin = (pipe_masked != 0).astype(bool)

    print(f"GT FG:       {gt_bin.sum():,}")
    print(f"Pred-nn FG:  {nn_bin.sum():,}")
    print(f"Pipeline FG: {pipe_bin.sum():,}")

    return gt_bin, nn_bin, pipe_masked, pipe_bin, shape


def mip_any(bool_vol):
    """Boolean OR projection along Z, Y, X axes.

    Returns:
        mip_z: (H, W) — looking from top
        mip_y: (D, W) — looking from front
        mip_x: (D, H) — looking from side
    """
    mip_z = np.any(bool_vol, axis=0)
    mip_y = np.any(bool_vol, axis=1)
    mip_x = np.any(bool_vol, axis=2)
    return mip_z, mip_y, mip_x


def mip_max(float_vol):
    """Max-intensity projection along Z, Y, X axes.

    Returns:
        mip_z: (H, W)
        mip_y: (D, W)
        mip_x: (D, H)
    """
    mip_z = np.max(float_vol, axis=0)
    mip_y = np.max(float_vol, axis=1)
    mip_x = np.max(float_vol, axis=2)
    return mip_z, mip_y, mip_x


def mip_sum(float_vol):
    """Sum projection along Z, Y, X axes.

    Returns:
        mip_z: (H, W)
        mip_y: (D, W)
        mip_x: (D, H)
    """
    mip_z = np.sum(float_vol, axis=0)
    mip_y = np.sum(float_vol, axis=1)
    mip_x = np.sum(float_vol, axis=2)
    return mip_z, mip_y, mip_x


def extract_borders(bin_vol_a, bin_vol_b=None):
    """Extract border voxels + surfel areas matching surface_distance package exactly.

    Uses the 2x2x2 neighbourhood encoding: correlate with ENCODE_NEIGHBOURHOOD_3D_KERNEL,
    border voxels have code != 0 and code != 255 (i.e., mixed occupancy in their 2x2x2 cube).

    If bin_vol_b is provided, does the bounding-box crop + pad on the union (matching the
    package's _crop_to_bounding_box behavior).

    Args:
        bin_vol_a: First binary volume
        bin_vol_b: Optional second binary volume (for union bbox crop)

    Returns:
        border_a: bool array (same shape as cropped+padded volume)
        area_map_a: float array of per-voxel surfel areas
        code_a: int array of neighbourhood codes
        bbox_min: array of crop origin (for mapping back to original coords)
        crop_shape: shape of the cropped volume before padding
    """
    kernel = lookup_tables.ENCODE_NEIGHBOURHOOD_3D_KERNEL
    area_table = lookup_tables.create_table_neighbour_code_to_surface_area((1.0, 1.0, 1.0))

    vol_a = bin_vol_a.astype(np.uint8)
    if bin_vol_b is not None:
        vol_b = bin_vol_b.astype(np.uint8)
        union = (vol_a | vol_b).astype(bool)
    else:
        union = vol_a.astype(bool)

    # Bounding box on union (matching package's _compute_bounding_box)
    proj_0 = np.amax(union, axis=(1, 2))
    idx_0 = np.nonzero(proj_0)[0]
    proj_1 = np.amax(union, axis=(0, 2))
    idx_1 = np.nonzero(proj_1)[0]
    proj_2 = np.amax(union, axis=(0, 1))
    idx_2 = np.nonzero(proj_2)[0]

    if len(idx_0) == 0:
        # Empty volume
        bbox_min = np.array([0, 0, 0])
        empty = np.zeros((2, 2, 2), dtype=bool)
        return empty, np.zeros_like(empty, dtype=float), np.zeros_like(empty, dtype=int), bbox_min, (0, 0, 0)

    bbox_min = np.array([idx_0.min(), idx_1.min(), idx_2.min()])
    bbox_max = np.array([idx_0.max(), idx_1.max(), idx_2.max()])
    crop_shape = tuple(bbox_max - bbox_min + 1)

    # Crop + pad (matching package's _crop_to_bounding_box — note +2 padding)
    cropmask_a = np.zeros((bbox_max - bbox_min) + 2, np.uint8)
    cropmask_a[0:-1, 0:-1, 0:-1] = vol_a[bbox_min[0]:bbox_max[0]+1,
                                           bbox_min[1]:bbox_max[1]+1,
                                           bbox_min[2]:bbox_max[2]+1]

    # Correlate
    code_a = ndimage.correlate(cropmask_a, kernel, mode="constant", cval=0)
    border_a = (code_a != 0) & (code_a != 255)
    area_map_a = area_table[code_a]

    return border_a, area_map_a, code_a, bbox_min, crop_shape


def extract_borders_pair(bin_vol_a, bin_vol_b):
    """Extract borders for a pair of volumes using shared bounding box.

    Returns:
        border_a, area_map_a, border_b, area_map_b, bbox_min, crop_shape
    """
    kernel = lookup_tables.ENCODE_NEIGHBOURHOOD_3D_KERNEL
    area_table = lookup_tables.create_table_neighbour_code_to_surface_area((1.0, 1.0, 1.0))

    vol_a = bin_vol_a.astype(np.uint8)
    vol_b = bin_vol_b.astype(np.uint8)
    union = (vol_a | vol_b).astype(bool)

    # Bounding box on union
    proj_0 = np.amax(union, axis=(1, 2))
    idx_0 = np.nonzero(proj_0)[0]
    proj_1 = np.amax(union, axis=(0, 2))
    idx_1 = np.nonzero(proj_1)[0]
    proj_2 = np.amax(union, axis=(0, 1))
    idx_2 = np.nonzero(proj_2)[0]

    if len(idx_0) == 0:
        bbox_min = np.array([0, 0, 0])
        empty = np.zeros((2, 2, 2), dtype=bool)
        z = np.zeros_like(empty, dtype=float)
        return empty, z, empty, z.copy(), bbox_min, (0, 0, 0)

    bbox_min = np.array([idx_0.min(), idx_1.min(), idx_2.min()])
    bbox_max = np.array([idx_0.max(), idx_1.max(), idx_2.max()])
    crop_shape = tuple(bbox_max - bbox_min + 1)

    padded_shape = (bbox_max - bbox_min) + 2

    # Crop + pad both volumes
    cropmask_a = np.zeros(padded_shape, np.uint8)
    cropmask_a[0:-1, 0:-1, 0:-1] = vol_a[bbox_min[0]:bbox_max[0]+1,
                                           bbox_min[1]:bbox_max[1]+1,
                                           bbox_min[2]:bbox_max[2]+1]

    cropmask_b = np.zeros(padded_shape, np.uint8)
    cropmask_b[0:-1, 0:-1, 0:-1] = vol_b[bbox_min[0]:bbox_max[0]+1,
                                           bbox_min[1]:bbox_max[1]+1,
                                           bbox_min[2]:bbox_max[2]+1]

    code_a = ndimage.correlate(cropmask_a, kernel, mode="constant", cval=0)
    code_b = ndimage.correlate(cropmask_b, kernel, mode="constant", cval=0)

    border_a = (code_a != 0) & (code_a != 255)
    border_b = (code_b != 0) & (code_b != 255)

    area_map_a = area_table[code_a]
    area_map_b = area_table[code_b]

    return border_a, area_map_a, border_b, area_map_b, bbox_min, crop_shape


def compute_surface_dice(border_a, area_map_a, border_b, area_map_b, tolerance=2.0):
    """Compute Surface Dice from pre-extracted borders.

    Returns:
        sd_val: float — surface dice value
        overlap_a: float — area of A surfels matched by B
        overlap_b: float — area of B surfels matched by A
        total_area: float — total surfel area
    """
    if not border_a.any() and not border_b.any():
        return 1.0, 0.0, 0.0, 0.0

    distmap_b = ndimage.distance_transform_edt(~border_b, sampling=(1.0, 1.0, 1.0))
    distmap_a = ndimage.distance_transform_edt(~border_a, sampling=(1.0, 1.0, 1.0))

    a_dists = distmap_b[border_a]
    b_dists = distmap_a[border_b]

    overlap_a = float(area_map_a[border_a][a_dists <= tolerance].sum())
    overlap_b = float(area_map_b[border_b][b_dists <= tolerance].sum())
    total_area = float(area_map_a[border_a].sum() + area_map_b[border_b].sum())

    if total_area == 0:
        return 1.0, 0.0, 0.0, 0.0

    sd_val = (overlap_a + overlap_b) / total_area
    return sd_val, overlap_a, overlap_b, total_area


def save_png(fig, cube_id, name, base_dir=None):
    """Save figure to output/{cube_id}/{cube_id}_{name}.png"""
    if base_dir is None:
        base_dir = DEFAULT_BASE_DIR
    out_dir = os.path.join(base_dir, "output", str(cube_id))
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{cube_id}_{name}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="black")
    print(f"Saved: {path}")
    return path


def make_mip_row(vol_z, vol_y, vol_x, gap=4, bg=0):
    """Concatenate Z/Y/X MIPs into a single row image with gaps.

    All inputs should be 2D arrays with the same dtype. Y and X MIPs are
    resized to match the Z MIP's height.

    Returns: 2D array (H, W_total)
    """
    h = vol_z.shape[0]

    # Pad Y and X to match height
    def pad_to_height(img, target_h):
        if img.shape[0] == target_h:
            return img
        if img.shape[0] < target_h:
            pad = target_h - img.shape[0]
            top = pad // 2
            bot = pad - top
            return np.pad(img, ((top, bot), (0, 0)), mode="constant", constant_values=bg)
        # Crop if taller
        start = (img.shape[0] - target_h) // 2
        return img[start:start + target_h]

    vy = pad_to_height(vol_y, h)
    vx = pad_to_height(vol_x, h)

    gap_col = np.full((h, gap), bg, dtype=vol_z.dtype) if gap > 0 else np.empty((h, 0), dtype=vol_z.dtype)
    return np.concatenate([vol_z, gap_col, vy, gap_col, vx], axis=1)


def make_rgb_mip_row(mips_z, mips_y, mips_x, gap=4):
    """Concatenate RGB MIP images (H,W,3) into a row with black gaps."""
    h = mips_z.shape[0]

    def pad_to_height(img, target_h):
        if img.shape[0] == target_h:
            return img
        if img.shape[0] < target_h:
            pad = target_h - img.shape[0]
            top = pad // 2
            bot = pad - top
            return np.pad(img, ((top, bot), (0, 0), (0, 0)), mode="constant", constant_values=0)
        start = (img.shape[0] - target_h) // 2
        return img[start:start + target_h]

    vy = pad_to_height(mips_y, h)
    vx = pad_to_height(mips_x, h)

    gap_col = np.zeros((h, gap, 3), dtype=mips_z.dtype) if gap > 0 else np.empty((h, 0, 3), dtype=mips_z.dtype)
    return np.concatenate([mips_z, gap_col, vy, gap_col, vx], axis=1)


def bool_to_rgb(mip_z, mip_y, mip_x, color=(255, 255, 255)):
    """Convert boolean MIPs to RGB images with the given color."""
    results = []
    for mip in [mip_z, mip_y, mip_x]:
        rgb = np.zeros((*mip.shape, 3), dtype=np.uint8)
        rgb[mip] = color
        results.append(rgb)
    return results


def classify_mips(gt_bin, pred_bin):
    """Compute TP/FP/FN boolean volumes and return RGB MIPs.

    Colors: green=TP, red=FP, blue=FN
    Returns: (mip_z, mip_y, mip_x) each (H,W,3) uint8
    """
    tp = gt_bin & pred_bin
    fp = ~gt_bin & pred_bin
    fn = gt_bin & ~pred_bin

    tp_z, tp_y, tp_x = mip_any(tp)
    fp_z, fp_y, fp_x = mip_any(fp)
    fn_z, fn_y, fn_x = mip_any(fn)

    results = []
    for tp_m, fp_m, fn_m in [(tp_z, fp_z, fn_z), (tp_y, fp_y, fn_y), (tp_x, fp_x, fn_x)]:
        rgb = np.zeros((*tp_m.shape, 3), dtype=np.uint8)
        rgb[tp_m] = (0, 200, 0)     # green = TP
        rgb[fp_m] = (200, 0, 0)     # red = FP
        rgb[fn_m] = (50, 50, 200)   # blue = FN
        # Where both TP and FP overlap in projection, TP wins
        results.append(rgb)

    return results[0], results[1], results[2]
