#!/usr/bin/env python3
"""Calibrate an original-RAW continuity score for atlas candidate bundles."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile


def load_bundle_module(path: Path):
    spec = importlib.util.spec_from_file_location("bundle_audit", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RawGrid:
    def __init__(self, directory: Path, chunk: int = 128) -> None:
        self.directory = directory
        self.chunk = chunk
        self.cache: dict[tuple[int, int, int], np.ndarray | None] = {}

    def cube(self, cz: int, cy: int, cx: int) -> np.ndarray | None:
        key = (cz, cy, cx)
        if key not in self.cache:
            path = self.directory / f"z{cz:05d}_y{cy:05d}_x{cx:05d}.tif"
            self.cache[key] = tifffile.imread(path) if path.exists() else None
        return self.cache[key]

    def fetch(self, iz: int, iy: int, ix: int) -> float:
        c = self.chunk
        cz = math.floor(iz / c) * c
        cy = math.floor(iy / c) * c
        cx = math.floor(ix / c) * c
        volume = self.cube(cz, cy, cx)
        if volume is None:
            return math.nan
        return float(volume[iz - cz, iy - cy, ix - cx])

    def sample(self, p: np.ndarray) -> float:
        lo = np.floor(p).astype(np.int64)
        d = p - lo
        total = 0.0
        weight = 0.0
        for mask in range(8):
            oz = mask & 1
            oy = (mask >> 1) & 1
            ox = (mask >> 2) & 1
            value = self.fetch(
                int(lo[0] + oz), int(lo[1] + oy), int(lo[2] + ox)
            )
            if not math.isfinite(value):
                continue
            w = (
                (d[0] if oz else 1.0 - d[0])
                * (d[1] if oy else 1.0 - d[1])
                * (d[2] if ox else 1.0 - d[2])
            )
            total += w * value
            weight += w
        return total / weight if weight > 1.0e-9 else math.nan


def normalize(v: np.ndarray) -> np.ndarray | None:
    length = float(np.linalg.norm(v))
    return v / length if length > 1.0e-9 else None


def sample_tensor(
    raw: RawGrid,
    p: np.ndarray,
    tu_in: np.ndarray,
    tv_in: np.ndarray,
    radius: float = 2.0,
) -> dict[str, float] | None:
    """Python transcription of src/common/raw_sample.c:sample_tangent_tensor."""
    tu = normalize(tu_in)
    if tu is None:
        return None
    tv = tv_in - float(np.dot(tv_in, tu)) * tu
    tv = normalize(tv)
    if tv is None:
        return None
    step = 0.5 * radius
    image = np.empty((5, 5), dtype=np.float64)
    for j in range(-2, 3):
        for i in range(-2, 3):
            q = p + step * (i * tu + j * tv)
            image[j + 2, i + 2] = raw.sample(q)
    if not np.isfinite(image).all():
        return None
    sxx = sxy = syy = wsum = 0.0
    for j in range(1, 4):
        for i in range(1, 4):
            gx = (image[j, i + 1] - image[j, i - 1]) / (2.0 * step)
            gy = (image[j + 1, i] - image[j - 1, i]) / (2.0 * step)
            w = (2.0 if i == 2 else 1.0) * (2.0 if j == 2 else 1.0)
            sxx += w * gx * gx
            sxy += w * gx * gy
            syy += w * gy * gy
            wsum += w
    sxx /= wsum
    sxy /= wsum
    syy /= wsum
    energy = sxx + syy
    a = sxx - syy
    b = 2.0 * sxy
    magnitude = math.hypot(a, b)
    coherence = magnitude / energy if energy > 1.0e-12 else 0.0
    rms = math.sqrt(max(0.0, energy))
    strength = rms / (rms + 8.0)
    axis = (
        0.5 * (1.0 + (a * a - b * b) / (magnitude * magnitude))
        if magnitude > 1.0e-12
        else 0.5
    )
    normal = np.cross(tu, tv)
    minus = raw.sample(p - 0.75 * normal)
    plus = raw.sample(p + 0.75 * normal)
    if not math.isfinite(minus) or not math.isfinite(plus):
        return None
    asymmetry = abs(plus - minus)
    ridge = 1.0 / (1.0 + (asymmetry / 16.0) ** 2)
    quality = ridge * strength * (0.25 + 0.75 * coherence * axis)
    return {
        "quality": min(1.0, max(0.0, quality)),
        "ridge": ridge,
        "strength": strength,
        "coherence": coherence,
        "axis": axis,
        "intensity": raw.sample(p),
    }


def sample_tangents(samples: list[dict]) -> list[np.ndarray]:
    by_stroke: dict[int, list[int]] = defaultdict(list)
    for sample in samples:
        by_stroke[sample["stroke"]].append(sample["sample"])
    tangent = [np.zeros(3, dtype=np.float64) for _ in samples]
    for ids in by_stroke.values():
        ids.sort(key=lambda sid: samples[sid]["ordinal"])
        for j, sid in enumerate(ids):
            a = ids[max(0, j - 1)]
            b = ids[min(len(ids) - 1, j + 1)]
            pa = np.array(
                [samples[a]["z"], samples[a]["y"], samples[a]["x"]],
                dtype=np.float64,
            )
            pb = np.array(
                [samples[b]["z"], samples[b]["y"], samples[b]["x"]],
                dtype=np.float64,
            )
            direction = normalize(pb - pa)
            tangent[sid] = (
                direction if direction is not None else np.array([0.0, 1.0, 0.0])
            )
    return tangent


def target_tangent(link: dict, samples: list[dict], tangents: list[np.ndarray]):
    # target_s is enough to interpolate tangent from the target stroke samples.
    ids = [
        sample["sample"]
        for sample in samples
        if sample["stroke"] == link["target_stroke"]
    ]
    ids.sort(key=lambda sid: samples[sid]["s"])
    if not ids:
        return np.array([0.0, 1.0, 0.0])
    nearest = min(ids, key=lambda sid: abs(samples[sid]["s"] - link["target_s"]))
    return tangents[nearest]


def score_link(
    raw: RawGrid,
    link: dict,
    samples: list[dict],
    tangents: list[np.ndarray],
    target_tangent_by_stroke: dict[int, list[tuple[float, np.ndarray]]],
) -> dict[str, float] | None:
    p0 = np.array(
        [link["source_z"], link["source_y"], link["source_x"]], dtype=np.float64
    )
    p1 = np.array(
        [link["target_z"], link["target_y"], link["target_x"]], dtype=np.float64
    )
    tv = normalize(p1 - p0)
    if tv is None:
        return None
    tu0 = tangents[link["source_sample"]]
    table = target_tangent_by_stroke[link["target_stroke"]]
    tu1 = min(table, key=lambda item: abs(item[0] - link["target_s"]))[1]
    values: dict[str, list[float]] = defaultdict(list)
    for t in (0.20, 0.35, 0.50, 0.65, 0.80):
        p = (1.0 - t) * p0 + t * p1
        tu = normalize((1.0 - t) * tu0 + t * tu1)
        if tu is None:
            continue
        tensor = sample_tensor(raw, p, tu, tv)
        if tensor is None:
            continue
        for name, value in tensor.items():
            values[name].append(value)
    if len(values["quality"]) < 4:
        return None
    result: dict[str, float] = {}
    for name, xs in values.items():
        result[f"{name}_median"] = statistics.median(xs)
        result[f"{name}_min"] = min(xs)
    result["quality_p20"] = float(np.quantile(values["quality"], 0.20))
    result["ridge_p20"] = float(np.quantile(values["ridge"], 0.20))
    result["intensity_p20"] = float(np.quantile(values["intensity"], 0.20))
    return result


def auc(positive: list[float], negative: list[float]) -> float:
    """Probability that a random positive has a larger score than a negative."""
    tagged = [(x, 1) for x in positive] + [(x, 0) for x in negative]
    tagged.sort(key=lambda item: item[0])
    rank_sum = 0.0
    i = 0
    while i < len(tagged):
        j = i + 1
        while j < len(tagged) and tagged[j][0] == tagged[i][0]:
            j += 1
        average_rank = 0.5 * ((i + 1) + j)
        rank_sum += average_rank * sum(label for _, label in tagged[i:j])
        i = j
    np_ = len(positive)
    nn = len(negative)
    return (rank_sum - np_ * (np_ + 1) / 2.0) / (np_ * nn)


def summarize(scores: list[dict[str, float]], metric: str) -> dict[str, float]:
    xs = [score[metric] for score in scores]
    return {
        "n": len(xs),
        "p10": float(np.quantile(xs, 0.10)),
        "p25": float(np.quantile(xs, 0.25)),
        "p50": float(np.quantile(xs, 0.50)),
        "p75": float(np.quantile(xs, 0.75)),
        "p90": float(np.quantile(xs, 0.90)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--atlas-dir",
        type=Path,
        default=Path("output/atlas_strip_4x5x5_fem_gauge"),
    )
    parser.add_argument(
        "--raw-dir", type=Path, default=Path("PHerc0139-4x5x5/cubes_RAW")
    )
    parser.add_argument("--limit", type=int, default=2500)
    parser.add_argument("--seed", type=int, default=13)
    args = parser.parse_args()
    atlas_dir: Path = args.atlas_dir
    audit = load_bundle_module(Path("tmp/analyze_atlas_bundle_hypotheses.py"))
    components, _, _ = audit.parse_mesh_components(atlas_dir / "atlas_world_uv.obj")
    samples = audit.parse_samples(atlas_dir / "samples.csv", components)
    links = audit.load_links(atlas_dir / "cross_sections.csv", samples)
    tangents = sample_tangents(samples)
    target_tangent_by_stroke: dict[int, list[tuple[float, np.ndarray]]] = defaultdict(list)
    for sample in samples:
        target_tangent_by_stroke[sample["stroke"]].append(
            (sample["s"], tangents[sample["sample"]])
        )

    positive = [
        link
        for link in links
        if link["topology_certified"] and link["prior_rank"] == 1
    ]
    negative = [
        link
        for link in links
        if not link["topology_certified"]
        and abs(link["du_initial"]) >= audit.GROSS_DU
        and link["has_topology_alternative"]
    ]
    bounded = [
        link
        for link in links
        if link["source_component"] == 206 and link["target_component"] == 205
    ]
    rng = random.Random(args.seed)
    rng.shuffle(positive)
    rng.shuffle(negative)
    positive = positive[: args.limit]
    negative = negative[: args.limit]

    raw = RawGrid(args.raw_dir)

    def score_many(rows: list[dict]) -> list[dict[str, float]]:
        out = []
        for link in rows:
            score = score_link(
                raw, link, samples, tangents, target_tangent_by_stroke
            )
            if score is not None:
                out.append(score)
        return out

    positive_score = score_many(positive)
    negative_score = score_many(negative)
    bounded_score = score_many(bounded)
    metrics = [
        "quality_p20",
        "quality_median",
        "ridge_p20",
        "ridge_median",
        "intensity_p20",
        "intensity_median",
    ]
    result = {
        "definitions": {
            "positive": "top-prior candidate on the same retained mesh component",
            "negative": (
                "gross cross-component candidate admitted beside a retained-"
                "topology continuation"
            ),
            "bounded": (
                "component pair 206->205, the sole coherent gross bundle with "
                "a bounded close front in the topology audit"
            ),
            "volume": "original cubes_RAW, not nnU-Net prediction",
        },
        "requested": {
            "positive": len(positive),
            "negative": len(negative),
            "bounded": len(bounded),
        },
        "scored": {
            "positive": len(positive_score),
            "negative": len(negative_score),
            "bounded": len(bounded_score),
        },
        "metrics": {},
        "loaded_raw_cubes": sum(v is not None for v in raw.cache.values()),
        "missing_raw_cubes": sum(v is None for v in raw.cache.values()),
    }
    for metric in metrics:
        p = [score[metric] for score in positive_score]
        n = [score[metric] for score in negative_score]
        b = [score[metric] for score in bounded_score]
        result["metrics"][metric] = {
            "auc_positive_over_negative": auc(p, n),
            "positive": summarize(positive_score, metric),
            "negative": summarize(negative_score, metric),
            "bounded": summarize(bounded_score, metric),
            "bounded_percentile_in_positive": (
                sum(x <= statistics.median(b) for x in p) / len(p)
            ),
            "bounded_percentile_in_negative": (
                sum(x <= statistics.median(b) for x in n) / len(n)
            ),
        }
    out_dir = atlas_dir / "bundle_hypothesis_audit"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "raw_continuity_summary.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
