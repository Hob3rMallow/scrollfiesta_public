#!/usr/bin/env python3
"""Select and render a fixed high-structure crop from ScrollFiesta textures.

Crop selection uses only the common pre-snap (stage-2) RAW texture, never either
candidate arm.  Every full-height window on a fixed stride is scored by valid
coverage, local structure-tensor coherence and robust gradient energy.  The
winning coordinates are then applied unchanged to both stage-4 outputs.

The score is a transparent visualization aid, not a validated fiber-quality or
legibility metric.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from scipy.ndimage import binary_erosion, gaussian_filter


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_texture(path: Path) -> np.ndarray:
    image = tifffile.imread(path)
    if image.ndim != 2 or image.dtype != np.uint8:
        raise ValueError(f"{path}: expected one uint8 image, got {image.dtype} {image.shape}")
    return image


def score_window(crop: np.ndarray, min_coverage: float) -> dict[str, float] | None:
    values = crop.astype(np.float32)
    valid = values > 0
    coverage = float(valid.mean())
    if coverage < min_coverage:
        return None
    interior = binary_erosion(valid, iterations=5)
    if interior.sum() < crop.size // 2:
        return None

    smooth = gaussian_filter(values, 1.0)
    grad_y, grad_x = np.gradient(smooth)
    j_xx = gaussian_filter(grad_x * grad_x, 3.0)
    j_yy = gaussian_filter(grad_y * grad_y, 3.0)
    j_xy = gaussian_filter(grad_x * grad_y, 3.0)
    coherence = np.sqrt((j_xx - j_yy) ** 2 + 4.0 * j_xy**2) / (j_xx + j_yy + 1e-6)
    gradient = np.sqrt(grad_x**2 + grad_y**2)
    grad_p50, grad_p90, grad_p99 = np.percentile(gradient[interior], [50.0, 90.0, 99.0])
    eligible = interior & (gradient <= grad_p99)
    median_coherence = float(np.median(coherence[eligible]))
    score = float(coverage * median_coherence * np.log1p(grad_p90))
    return {
        "score": score,
        "coverage": coverage,
        "median_local_coherence": median_coherence,
        "gradient_p50": float(grad_p50),
        "gradient_p90": float(grad_p90),
        "gradient_p99": float(grad_p99),
    }


def select_crop(image: np.ndarray, width: int, stride: int, min_coverage: float) -> dict[str, float | int]:
    if width <= 0 or width > image.shape[1]:
        raise ValueError("window width must fit the texture")
    if stride <= 0:
        raise ValueError("stride must be positive")
    candidates: list[dict[str, float | int]] = []
    for x0 in range(0, image.shape[1] - width + 1, stride):
        metrics = score_window(image[:, x0 : x0 + width], min_coverage)
        if metrics is not None:
            candidates.append({"x0": x0, "x1": x0 + width, **metrics})
    if not candidates:
        raise ValueError("no window passed the minimum-coverage gate")
    # Stable and explicit tie-break: highest score, then smallest x0.
    candidates.sort(key=lambda item: (-float(item["score"]), int(item["x0"])))
    winner = dict(candidates[0])
    winner["eligible_windows"] = len(candidates)
    return winner


def render(args: argparse.Namespace) -> dict[str, object]:
    common = read_texture(args.common_pre)
    repair = read_texture(args.repair)
    recto = read_texture(args.recto)
    if common.shape != repair.shape or common.shape != recto.shape:
        raise ValueError("common, repair and recto textures must have the same shape")

    selection = select_crop(common, args.width, args.stride, args.min_coverage)
    x0, x1 = int(selection["x0"]), int(selection["x1"])
    repair_crop = repair[:, x0:x1]
    recto_crop = recto[:, x0:x1]
    common_crop = common[:, x0:x1]
    valid = common_crop > 0
    display_values = common_crop[valid]
    vmin, vmax = np.percentile(display_values, [1.0, 99.0])
    difference = repair_crop.astype(np.int16) - recto_crop.astype(np.int16)
    common_valid = (repair_crop > 0) & (recto_crop > 0)
    if not common_valid.any():
        raise ValueError("selected stage-4 crops have no common valid pixels")

    abs_difference = np.abs(difference[common_valid])
    limit = max(1.0, float(np.percentile(abs_difference, 99.0)))
    figure, axes = plt.subplots(1, 3, figsize=(18, 5.5), constrained_layout=True)
    extent = (x0, x1, 0, common.shape[0])
    image_options = dict(
        cmap="gray",
        vmin=vmin,
        vmax=vmax,
        origin="lower",
        extent=extent,
        interpolation="nearest",
        aspect="equal",
    )
    axes[0].imshow(repair_crop, **image_options)
    axes[1].imshow(recto_crop, **image_options)
    diff_artist = axes[2].imshow(
        difference,
        cmap="coolwarm",
        vmin=-limit,
        vmax=limit,
        origin="lower",
        extent=extent,
        interpolation="nearest",
        aspect="equal",
    )
    axes[0].set_title("Repair only (recto iterations 0)")
    axes[1].set_title("Oriented boundary (recto iterations 4)")
    axes[2].set_title("Pixel difference: repair only − 4 iterations")
    for axis in axes:
        axis.set_xlabel("unrolled u (pixels)")
        axis.set_ylabel("unrolled v (pixels)")
    colorbar = figure.colorbar(diff_artist, ax=axes[2], fraction=0.047, pad=0.03)
    colorbar.set_label(f"gray levels (symmetric p99 = {limit:.0f})")
    figure.suptitle(
        "PHerc0139 deterministic high-structure texture crop\n"
        f"Selected only on common pre-snap texture: width={args.width}, stride={args.stride}, "
        f"coverage≥{args.min_coverage:.2f}; identical crop and p1–p99 grayscale",
        fontsize=14,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out, dpi=180, facecolor="white")
    plt.close(figure)
    return {
        "output": str(args.out),
        "output_sha256": sha256(args.out),
        "texture_shape": list(common.shape),
        "selection": selection,
        "display_p1_p99": [float(vmin), float(vmax)],
        "difference_symmetric_p99": limit,
        "common_stage4_pixels": int(common_valid.sum()),
        "mean_abs_difference": float(abs_difference.mean()),
        "median_abs_difference": float(np.median(abs_difference)),
        "common_pre_sha256": sha256(args.common_pre),
        "repair_sha256": sha256(args.repair),
        "recto_sha256": sha256(args.recto),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--common-pre", type=Path, required=True)
    parser.add_argument("--repair", type=Path, required=True)
    parser.add_argument("--recto", type=Path, required=True)
    parser.add_argument("--width", type=int, default=768)
    parser.add_argument("--stride", type=int, default=256)
    parser.add_argument("--min-coverage", type=float, default=0.90)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(render(parse_args()), indent=2, sort_keys=True))
