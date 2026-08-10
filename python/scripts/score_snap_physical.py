"""Score paired PHerc1203 snap meshes against the physical recto labels.

This is the frozen scorer for ``PHERC1203_PHYSICAL_SNAP_REPLICATION_PREREG``.
It compares the same vertices and topology in the PR #11 repair-only arm and
the former four-iteration arm.  The reference labels and their provenance are
from 7jycwjmbfn-eng/pherc0139-physical-audit / Villa PR #1382.

Distances are computed independently in every axial L1 label plane, then
sampled trilinearly.  Eligibility and surface-area weights are fixed from the
common pre-snap mesh before either arm is summarized.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
from pathlib import Path

import numpy as np
import zarr
from scipy.ndimage import distance_transform_edt, map_coordinates

from score_snap_ridge import load_obj, sha256


LABEL_ORIGIN_L1 = np.array([3936.0, 0.0, 0.0], dtype=np.float64)
LABEL_SHAPE = (2016, 3456, 3456)
LABEL_WINDOW = (
    slice(976, 1072),
    slice(1408, 2048),
    slice(1408, 2048),
)
WINDOW_ORIGIN = np.array(
    [part.start for part in LABEL_WINDOW], dtype=np.float64
)
BOX_LOW_L0 = np.array([9856.0, 3072.0, 3072.0], dtype=np.float64)
BOX_HIGH_L0 = np.array([9984.0, 3712.0, 3712.0], dtype=np.float64)
PRE_BOX_MARGIN_L0 = 4.0
ARM_BOX_MARGIN_L0 = 1.0
PRE_DISTANCE_MAX_L1 = 3.0
NULL_SHIFT_L1 = 64.0
NULL_SHIFT_VECTOR = np.array([0.0, NULL_SHIFT_L1, 0.0], dtype=np.float64)
L1_VOXEL_UM = 18.724
REGISTRATION_P95_UM = 6.1
REGISTRATION_RADIUS_L1 = REGISTRATION_P95_UM / L1_VOXEL_UM
CROP_SAFE_DISTANCE_L1 = 64.0
MIN_VERTICES_PER_CUBE = 100
MIN_CUBES = 20
MIN_CUBES_PER_ROW_COLUMN = 3
MIN_CUBE_AREA_FRACTION = 0.20
MIN_GLOBAL_AREA_FRACTION = 0.40
EFFECT_MIN_L1 = 0.33
LARGE_WORSENING_L1 = 0.25
ROW_BOOTSTRAP_SEED = 20260810
COLUMN_BOOTSTRAP_SEED = 20260811
BOOTSTRAP_N = 10_000
EXPECTED_Z = 9856
EXPECTED_YX = tuple(range(3072, 3712, 128))
CUBE_NAME_RE = re.compile(r"^z(\d{5})_y(\d{5})_x(\d{5})$")


def l0_to_window_coordinates(points_l0: np.ndarray) -> np.ndarray:
    """Map z/y/x L0 index coordinates to the frozen label-window indices."""
    points = np.asarray(points_l0, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n, 3) in z/y/x order")
    return (points - 0.5) / 2.0 - LABEL_ORIGIN_L1 - WINDOW_ORIGIN


def all_trilinear_corners_true(mask: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Return whether all eight floor/ceil corners are true for each point."""
    mask = np.asarray(mask, dtype=bool)
    coords = np.asarray(coords, dtype=np.float64)
    if mask.ndim != 3 or coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("mask must be 3-D and coords must have shape (n, 3)")
    finite = np.isfinite(coords).all(axis=1)
    safe = np.where(np.isfinite(coords), coords, -1.0)
    base = np.floor(safe).astype(np.int64)
    inside = finite & np.all(base >= 0, axis=1)
    inside &= np.all(base + 1 < np.asarray(mask.shape), axis=1)
    result = inside.copy()
    rows = np.flatnonzero(inside)
    if not len(rows):
        return result
    b = base[rows]
    corner_ok = np.ones(len(rows), dtype=bool)
    for dz in (0, 1):
        for dy in (0, 1):
            for dx in (0, 1):
                corner_ok &= mask[b[:, 0] + dz, b[:, 1] + dy, b[:, 2] + dx]
    result[rows] = corner_ok
    return result


def any_trilinear_corner_true(mask: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Return whether any floor/ceil corner is true for each in-bounds point."""
    mask = np.asarray(mask, dtype=bool)
    coords = np.asarray(coords, dtype=np.float64)
    if mask.ndim != 3 or coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("mask must be 3-D and coords must have shape (n, 3)")
    finite = np.isfinite(coords).all(axis=1)
    safe = np.where(np.isfinite(coords), coords, -1.0)
    base = np.floor(safe).astype(np.int64)
    inside = finite & np.all(base >= 0, axis=1)
    inside &= np.all(base + 1 < np.asarray(mask.shape), axis=1)
    result = np.zeros(len(coords), dtype=bool)
    rows = np.flatnonzero(inside)
    if not len(rows):
        return result
    b = base[rows]
    corner_hit = np.zeros(len(rows), dtype=bool)
    for dz in (0, 1):
        for dy in (0, 1):
            for dx in (0, 1):
                corner_hit |= mask[
                    b[:, 0] + dz, b[:, 1] + dy, b[:, 2] + dx
                ]
    result[rows] = corner_hit
    return result


def registration_offsets() -> np.ndarray:
    """Return zero plus 26 fixed directions on the registration-p95 sphere."""
    offsets = [np.zeros(3, dtype=np.float64)]
    for values in itertools.product((-1.0, 0.0, 1.0), repeat=3):
        vector = np.asarray(values, dtype=np.float64)
        length = float(np.linalg.norm(vector))
        if length == 0:
            continue
        offsets.append(vector / length * REGISTRATION_RADIUS_L1)
    return np.asarray(offsets, dtype=np.float64)


def fixed_vertex_area_weights(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Assign one third of every pre-snap triangle area to each incident vertex."""
    vertices = np.asarray(vertices, dtype=np.float64)
    faces = np.asarray(faces, dtype=np.int64)
    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError("vertices must have shape (n, 3)")
    if faces.ndim != 2 or faces.shape[1] != 3:
        raise ValueError("faces must have shape (m, 3)")
    if len(faces) and (faces.min() < 0 or faces.max() >= len(vertices)):
        raise ValueError("face index outside vertex array")
    tri = vertices[faces]
    area = 0.5 * np.linalg.norm(
        np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1
    )
    weights = np.zeros(len(vertices), dtype=np.float64)
    for corner in range(3):
        np.add.at(weights, faces[:, corner], area / 3.0)
    return weights


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    """Return a deterministic lower weighted median of finite positive weight."""
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if values.shape != weights.shape or values.ndim != 1:
        raise ValueError("values and weights must be equal-length 1-D arrays")
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        raise ValueError("weighted median has no finite positive-weight values")
    v = values[valid]
    w = weights[valid]
    order = np.argsort(v, kind="stable")
    v, w = v[order], w[order]
    cutoff = 0.5 * float(w.sum())
    index = int(np.searchsorted(np.cumsum(w), cutoff, side="left"))
    return float(v[min(index, len(v) - 1)])


def bootstrap_median_ci(
    values: np.ndarray, *, seed: int = ROW_BOOTSTRAP_SEED
) -> list[float | None]:
    """Frozen percentile interval for five spatial-cluster medians."""
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1:
        raise ValueError("bootstrap values must be one-dimensional")
    if not len(values):
        return [None, None]
    if not np.isfinite(values).all():
        raise ValueError("bootstrap values must all be finite")
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(values), size=(BOOTSTRAP_N, len(values)))
    medians = np.median(values[draws], axis=1)
    return [float(value) for value in np.percentile(medians, [2.5, 97.5])]


def build_axial_recto_distance(labels: np.ndarray) -> np.ndarray:
    """Build a stack of independent 2-D distances to recto-band bit 8."""
    labels = np.asarray(labels, dtype=np.uint8)
    if labels.ndim != 3:
        raise ValueError("labels must be a 3-D uint8 array")
    recto = (labels & np.uint8(8)) != 0
    distances = np.full(labels.shape, np.inf, dtype=np.float32)
    for z_index in range(labels.shape[0]):
        if recto[z_index].any():
            distances[z_index] = distance_transform_edt(
                ~recto[z_index]
            ).astype(np.float32)
    return distances


def sample_distances(distance_stack: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """Trilinearly sample the axial distance stack."""
    values = map_coordinates(
        distance_stack,
        np.asarray(coords, dtype=np.float64).T,
        order=1,
        mode="constant",
        cval=np.inf,
        prefilter=False,
    )
    return np.asarray(values, dtype=np.float64)


def read_frozen_label_window(labels_path: Path) -> tuple[np.ndarray, dict]:
    """Read and validate only the preregistered PHerc1203 label window."""
    array = zarr.open_array(str(labels_path), mode="r")
    if tuple(array.shape) != LABEL_SHAPE or np.dtype(array.dtype) != np.uint8:
        raise ValueError(
            f"unexpected label shape/dtype: {array.shape} {array.dtype}"
        )
    attrs = dict(array.attrs)
    if attrs.get("origin_l1") != [3936, 0, 0]:
        raise ValueError("label origin_l1 does not match the frozen reference")
    expected_bits = {
        "valid": 1,
        "material": 2,
        "centerline": 4,
        "recto_band": 8,
        "boundary_poor": 16,
    }
    if attrs.get("bits") != expected_bits:
        raise ValueError("label bit schema does not match the frozen reference")
    labels = np.asarray(array[LABEL_WINDOW], dtype=np.uint8)
    expected_window_shape = tuple(part.stop - part.start for part in LABEL_WINDOW)
    if labels.shape != expected_window_shape:
        raise ValueError(f"short label window: {labels.shape}")
    metadata = {
        "path": labels_path.as_posix(),
        "shape_zyx": list(array.shape),
        "dtype": str(np.dtype(array.dtype)),
        "origin_l1": attrs["origin_l1"],
        "bits": attrs["bits"],
        "registration_heldout_um": attrs.get("registration_heldout_um"),
        "window_local_zyx": [
            [part.start, part.stop] for part in LABEL_WINDOW
        ],
        "window_sha256": hashlib.sha256(
            np.ascontiguousarray(labels).tobytes()
        ).hexdigest(),
        "zarray_sha256": sha256(labels_path / ".zarray"),
        "zattrs_sha256": sha256(labels_path / ".zattrs"),
    }
    return labels, metadata


def validate_meshes(pre: dict, candidate: dict, former: dict) -> None:
    """Enforce the shared-index contract required for a paired comparison."""
    for name, arm in (("candidate", candidate), ("former", former)):
        if arm["verts"].shape != pre["verts"].shape:
            raise ValueError(f"{name}: vertex count differs from pre")
        if not np.array_equal(arm["faces"], pre["faces"]):
            raise ValueError(f"{name}: topology/face indices differ from pre")
        if not np.array_equal(arm["uv"], pre["uv"]):
            raise ValueError(f"{name}: registered UV differs from pre")
        if arm["cubes"] != pre["cubes"]:
            raise ValueError(f"{name}: cube ranges differ from pre")
    for name, mesh in (("pre", pre), ("candidate", candidate), ("former", former)):
        if not np.isfinite(mesh["verts"]).all():
            raise ValueError(f"{name}: non-finite vertex coordinate")


def cube_coordinates(cube_name: str) -> tuple[int, int, int]:
    match = CUBE_NAME_RE.fullmatch(cube_name)
    if not match:
        raise ValueError(f"malformed cube name: {cube_name}")
    return tuple(int(value) for value in match.groups())


def validate_cube_inventory(cubes: list[tuple[str, int, int]]) -> None:
    names = [item[0] for item in cubes]
    if len(names) != len(set(names)):
        raise ValueError("duplicate cube provenance name")
    expected = {
        f"z{EXPECTED_Z:05d}_y{y:05d}_x{x:05d}"
        for y in EXPECTED_YX
        for x in EXPECTED_YX
    }
    extra = sorted(set(names) - expected)
    if extra:
        raise ValueError(f"cube outside frozen inventory: {extra[0]}")


def _mesh_record(mesh: dict) -> dict:
    return {"path": mesh["path"], "sha256": mesh["sha256"]}


def _row_column_summary(per_cube: list[dict], axis: str) -> dict[str, dict]:
    index = 1 if axis == "y" else 2
    out: dict[str, dict] = {}
    for coordinate in EXPECTED_YX:
        rows = [
            row
            for row in per_cube
            if row["qualified"]
            and cube_coordinates(row["cube"])[index] == coordinate
        ]
        effects = np.asarray(
            [row["corrected_effect_l1"] for row in rows], dtype=np.float64
        )
        out[str(coordinate)] = {
            "n_cubes": len(rows),
            "median_corrected_effect_l1": (
                float(np.median(effects)) if len(effects) else None
            ),
        }
    return out


def _all_axis_medians(summary: dict[str, dict], positive: bool) -> bool:
    for row in summary.values():
        value = row["median_corrected_effect_l1"]
        if row["n_cubes"] < MIN_CUBES_PER_ROW_COLUMN or value is None:
            return False
        if positive and not value > 0:
            return False
        if not positive and not value < 0:
            return False
    return True


def _inside_box(points: np.ndarray, margin: float) -> np.ndarray:
    low = BOX_LOW_L0 + margin
    high = BOX_HIGH_L0 - margin
    return np.all((points >= low) & (points <= high), axis=1)


def _finite_crop_safe(values: np.ndarray) -> np.ndarray:
    return np.isfinite(values) & (values < CROP_SAFE_DISTANCE_L1)


def _subset_metrics(
    indices: np.ndarray,
    weights: np.ndarray,
    real: dict[str, np.ndarray],
    null: dict[str, np.ndarray],
) -> dict[str, float | None]:
    keys = (
        "pre_real_weighted_median_l1",
        "iter0_real_weighted_median_l1",
        "iter4_real_weighted_median_l1",
        "iter0_null_weighted_median_l1",
        "iter4_null_weighted_median_l1",
        "real_effect_l1",
        "null_effect_l1",
        "corrected_effect_l1",
        "unweighted_real_effect_l1",
        "unweighted_null_effect_l1",
        "unweighted_corrected_effect_l1",
    )
    if not len(indices):
        return {key: None for key in keys}
    local_weights = weights[indices]
    weighted = {
        "pre_real": weighted_median(real["pre"][indices], local_weights),
        "iter0_real": weighted_median(real["iter0"][indices], local_weights),
        "iter4_real": weighted_median(real["iter4"][indices], local_weights),
        "iter0_null": weighted_median(null["iter0"][indices], local_weights),
        "iter4_null": weighted_median(null["iter4"][indices], local_weights),
    }
    real_effect = weighted["iter4_real"] - weighted["iter0_real"]
    null_effect = weighted["iter4_null"] - weighted["iter0_null"]
    unweighted = {
        f"{arm}_{kind}": float(np.median(source[arm][indices]))
        for kind, source in (("real", real), ("null", null))
        for arm in ("iter0", "iter4")
    }
    unweighted_real = unweighted["iter4_real"] - unweighted["iter0_real"]
    unweighted_null = unweighted["iter4_null"] - unweighted["iter0_null"]
    return {
        "pre_real_weighted_median_l1": weighted["pre_real"],
        "iter0_real_weighted_median_l1": weighted["iter0_real"],
        "iter4_real_weighted_median_l1": weighted["iter4_real"],
        "iter0_null_weighted_median_l1": weighted["iter0_null"],
        "iter4_null_weighted_median_l1": weighted["iter4_null"],
        "real_effect_l1": real_effect,
        "null_effect_l1": null_effect,
        "corrected_effect_l1": real_effect - null_effect,
        "unweighted_real_effect_l1": unweighted_real,
        "unweighted_null_effect_l1": unweighted_null,
        "unweighted_corrected_effect_l1": unweighted_real - unweighted_null,
    }


def _prefix_metrics(prefix: str, values: dict[str, float | None]) -> dict:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def score(args: argparse.Namespace) -> dict:
    labels, label_metadata = read_frozen_label_window(Path(args.labels))
    valid_labels = (labels & np.uint8(1)) != 0
    boundary_poor_labels = (labels & np.uint8(16)) != 0
    distance_stack = build_axial_recto_distance(labels)

    pre = load_obj(Path(args.pre))
    candidate = load_obj(Path(args.candidate))
    former = load_obj(Path(args.former))
    validate_meshes(pre, candidate, former)
    validate_cube_inventory(pre["cubes"])

    meshes = {"pre": pre, "iter0": candidate, "iter4": former}
    coords = {
        name: l0_to_window_coordinates(mesh["verts"])
        for name, mesh in meshes.items()
    }
    weights = fixed_vertex_area_weights(pre["verts"], pre["faces"])
    used = np.zeros(len(pre["verts"]), dtype=bool)
    used[pre["faces"].ravel()] = True
    pre_reference = (
        used
        & np.isfinite(weights)
        & (weights > 0)
        & _inside_box(pre["verts"].astype(np.float64), PRE_BOX_MARGIN_L0)
    )
    pre_target_all_offsets = np.ones(len(pre["verts"]), dtype=bool)
    boundary_poor_any_offset = np.zeros(len(pre["verts"]), dtype=bool)
    offsets = registration_offsets()
    offset_outcomes: list[
        tuple[dict[str, np.ndarray], dict[str, np.ndarray]]
    ] = []
    for offset in offsets:
        real_coords = {name: value + offset for name, value in coords.items()}
        null_coords = {
            name: value + offset + NULL_SHIFT_VECTOR
            for name, value in coords.items()
        }
        real = {
            name: sample_distances(distance_stack, value)
            for name, value in real_coords.items()
        }
        null = {
            name: sample_distances(distance_stack, value)
            for name, value in null_coords.items()
        }
        offset_outcomes.append((real, null))
        pre_reference &= all_trilinear_corners_true(
            valid_labels, real_coords["pre"]
        )
        pre_reference &= all_trilinear_corners_true(
            valid_labels, null_coords["pre"]
        )
        pre_reference &= _finite_crop_safe(real["pre"])
        pre_reference &= _finite_crop_safe(null["pre"])
        pre_target_all_offsets &= real["pre"] <= PRE_DISTANCE_MAX_L1
        boundary_poor_any_offset |= any_trilinear_corner_true(
            boundary_poor_labels, real_coords["pre"]
        )

    physical_target = pre_reference & pre_target_all_offsets
    primary_predeclared = physical_target & ~boundary_poor_any_offset
    boundary_poor_predeclared = physical_target & boundary_poor_any_offset

    arm_support = np.ones(len(pre["verts"]), dtype=bool)
    for name in ("iter0", "iter4"):
        arm_support &= _inside_box(
            meshes[name]["verts"].astype(np.float64), ARM_BOX_MARGIN_L0
        )
        for offset, (real, null) in zip(offsets, offset_outcomes):
            real_coords = coords[name] + offset
            null_coords = real_coords + NULL_SHIFT_VECTOR
            arm_support &= all_trilinear_corners_true(valid_labels, real_coords)
            arm_support &= all_trilinear_corners_true(valid_labels, null_coords)
            arm_support &= _finite_crop_safe(real[name])
            arm_support &= _finite_crop_safe(null[name])

    primary_support_fail = primary_predeclared & ~arm_support
    boundary_poor_support_fail = boundary_poor_predeclared & ~arm_support
    arm_support_gate = not bool(primary_support_fail.any())
    scored_primary = primary_predeclared & arm_support
    scored_boundary_poor = boundary_poor_predeclared & arm_support
    central_real, central_null = offset_outcomes[0]
    sampled_evaluations = 0
    unsafe_distance_evaluations = 0
    max_finite_sampled_distance = None
    for real, null in offset_outcomes:
        for source in (real, null):
            for name in ("pre", "iter0", "iter4"):
                values = source[name][primary_predeclared]
                sampled_evaluations += len(values)
                unsafe_distance_evaluations += int(
                    (~_finite_crop_safe(values)).sum()
                )
                finite = values[np.isfinite(values)]
                if len(finite):
                    local_max = float(finite.max())
                    if (
                        max_finite_sampled_distance is None
                        or local_max > max_finite_sampled_distance
                    ):
                        max_finite_sampled_distance = local_max

    per_cube: list[dict] = []
    for cube_name, lo, hi in pre["cubes"]:
        reference_indices = np.flatnonzero(pre_reference[lo:hi]) + lo
        primary_indices = np.flatnonzero(primary_predeclared[lo:hi]) + lo
        scored_indices = np.flatnonzero(scored_primary[lo:hi]) + lo
        poor_predeclared_indices = (
            np.flatnonzero(boundary_poor_predeclared[lo:hi]) + lo
        )
        poor_scored_indices = np.flatnonzero(scored_boundary_poor[lo:hi]) + lo
        reference_area = float(weights[reference_indices].sum())
        primary_area = float(weights[primary_indices].sum())
        area_fraction = (
            primary_area / reference_area if reference_area > 0 else None
        )
        cube_support_failures = int(primary_support_fail[lo:hi].sum())
        row: dict[str, object] = {
            "cube": cube_name,
            "n_used": int(used[lo:hi].sum()),
            "n_pre_reference": int(len(reference_indices)),
            "n_primary_predeclared": int(len(primary_indices)),
            "n_primary_scored": int(len(scored_indices)),
            "n_arm_support_failures": cube_support_failures,
            "n_boundary_poor_predeclared": int(len(poor_predeclared_indices)),
            "n_boundary_poor_scored": int(len(poor_scored_indices)),
            "n_boundary_poor_arm_support_failures": int(
                boundary_poor_support_fail[lo:hi].sum()
            ),
            "pre_reference_area_l0_sq": reference_area,
            "primary_area_l0_sq": primary_area,
            "primary_area_fraction": area_fraction,
            "boundary_poor_predeclared_area_l0_sq": float(
                weights[poor_predeclared_indices].sum()
            ),
            "boundary_poor_scored_area_l0_sq": float(
                weights[poor_scored_indices].sum()
            ),
        }
        row.update(
            _subset_metrics(
                scored_indices, weights, central_real, central_null
            )
        )
        row.update(
            _prefix_metrics(
                "boundary_poor",
                _subset_metrics(
                    poor_scored_indices, weights, central_real, central_null
                ),
            )
        )
        row["qualified"] = bool(
            len(primary_indices) >= MIN_VERTICES_PER_CUBE
            and len(scored_indices) == len(primary_indices)
            and area_fraction is not None
            and area_fraction >= MIN_CUBE_AREA_FRACTION
            and cube_support_failures == 0
        )
        per_cube.append(row)

    qualified = [row for row in per_cube if row["qualified"]]
    effects = np.asarray(
        [row["corrected_effect_l1"] for row in qualified], dtype=np.float64
    )
    real_effects = np.asarray(
        [row["real_effect_l1"] for row in qualified], dtype=np.float64
    )
    unweighted_effects = np.asarray(
        [row["unweighted_corrected_effect_l1"] for row in qualified],
        dtype=np.float64,
    )
    row_summary = _row_column_summary(per_cube, "y")
    column_summary = _row_column_summary(per_cube, "x")
    row_medians = np.asarray(
        [
            item["median_corrected_effect_l1"]
            for item in row_summary.values()
            if item["median_corrected_effect_l1"] is not None
        ],
        dtype=np.float64,
    )
    column_medians = np.asarray(
        [
            item["median_corrected_effect_l1"]
            for item in column_summary.values()
            if item["median_corrected_effect_l1"] is not None
        ],
        dtype=np.float64,
    )
    row_bootstrap95 = (
        bootstrap_median_ci(row_medians, seed=ROW_BOOTSTRAP_SEED)
        if len(row_medians) == len(EXPECTED_YX)
        else [None, None]
    )
    column_bootstrap95 = (
        bootstrap_median_ci(column_medians, seed=COLUMN_BOOTSTRAP_SEED)
        if len(column_medians) == len(EXPECTED_YX)
        else [None, None]
    )
    reference_area_total = float(weights[pre_reference].sum())
    primary_area_total = float(weights[primary_predeclared].sum())
    global_area_fraction = (
        primary_area_total / reference_area_total
        if reference_area_total > 0 else None
    )
    coverage_gate = (
        MIN_CUBES <= len(pre["cubes"]) <= 25
        and arm_support_gate
        and len(qualified) >= MIN_CUBES
        and global_area_fraction is not None
        and global_area_fraction >= MIN_GLOBAL_AREA_FRACTION
        and all(
            row["n_cubes"] >= MIN_CUBES_PER_ROW_COLUMN
            for row in [*row_summary.values(), *column_summary.values()]
        )
    )

    median_effect = float(np.median(effects)) if len(effects) else None
    median_real_effect = (
        float(np.median(real_effects)) if len(real_effects) else None
    )
    positive_fraction = float((effects > 0).mean()) if len(effects) else None
    negative_fraction = float((effects < 0).mean()) if len(effects) else None
    candidate_large_worse = (
        float((effects < -LARGE_WORSENING_L1).mean()) if len(effects) else None
    )
    former_large_worse = (
        float((effects > LARGE_WORSENING_L1).mean()) if len(effects) else None
    )
    boundary_poor_rows = [
        row
        for row in per_cube
        if row["n_boundary_poor_scored"] >= MIN_VERTICES_PER_CUBE
        and row["boundary_poor_corrected_effect_l1"] is not None
    ]
    boundary_poor_summary = {
        "descriptive_only": True,
        "cubes_with_min_vertices": len(boundary_poor_rows),
        "predeclared_vertices": int(boundary_poor_predeclared.sum()),
        "scored_vertices": int(scored_boundary_poor.sum()),
        "arm_support_failures": int(boundary_poor_support_fail.sum()),
        "predeclared_area_l0_sq": float(
            weights[boundary_poor_predeclared].sum()
        ),
        "scored_area_l0_sq": float(
            weights[scored_boundary_poor].sum()
        ),
        "median_real_effect_l1": (
            float(
                np.median(
                    [
                        row["boundary_poor_real_effect_l1"]
                        for row in boundary_poor_rows
                    ]
                )
            )
            if boundary_poor_rows else None
        ),
        "median_corrected_effect_l1": (
            float(
                np.median(
                    [
                        row["boundary_poor_corrected_effect_l1"]
                        for row in boundary_poor_rows
                    ]
                )
            )
            if boundary_poor_rows else None
        ),
    }

    cube_ranges = {name: (lo, hi) for name, lo, hi in pre["cubes"]}
    sensitivity: list[dict] = []
    for offset, (real, null) in zip(offsets, offset_outcomes):
        offset_real_effects: list[float] = []
        offset_corrected_effects: list[float] = []
        for row in qualified:
            lo, hi = cube_ranges[row["cube"]]
            indices = np.flatnonzero(scored_primary[lo:hi]) + lo
            metrics = _subset_metrics(indices, weights, real, null)
            offset_real_effects.append(float(metrics["real_effect_l1"]))
            offset_corrected_effects.append(
                float(metrics["corrected_effect_l1"])
            )
        sensitivity.append(
            {
                "offset_zyx_l1": [float(value) for value in offset],
                "median_real_effect_l1": (
                    float(np.median(offset_real_effects))
                    if offset_real_effects else None
                ),
                "median_corrected_effect_l1": (
                    float(np.median(offset_corrected_effects))
                    if offset_corrected_effects else None
                ),
            }
        )
    candidate_registration_stable = bool(
        sensitivity
        and all(
            row["median_real_effect_l1"] is not None
            and row["median_real_effect_l1"] > 0
            and row["median_corrected_effect_l1"] is not None
            and row["median_corrected_effect_l1"] > 0
            for row in sensitivity
        )
    )
    former_registration_stable = bool(
        sensitivity
        and all(
            row["median_real_effect_l1"] is not None
            and row["median_real_effect_l1"] < 0
            and row["median_corrected_effect_l1"] is not None
            and row["median_corrected_effect_l1"] < 0
            for row in sensitivity
        )
    )

    candidate_pass = bool(
        coverage_gate
        and median_effect is not None
        and median_effect >= EFFECT_MIN_L1
        and median_real_effect is not None
        and median_real_effect > 0
        and row_bootstrap95[0] is not None
        and row_bootstrap95[0] > 0
        and column_bootstrap95[0] is not None
        and column_bootstrap95[0] > 0
        and positive_fraction is not None
        and positive_fraction >= 0.80
        and _all_axis_medians(row_summary, True)
        and _all_axis_medians(column_summary, True)
        and candidate_large_worse is not None
        and candidate_large_worse <= 0.10
        and candidate_registration_stable
    )
    former_pass = bool(
        coverage_gate
        and median_effect is not None
        and median_effect <= -EFFECT_MIN_L1
        and median_real_effect is not None
        and median_real_effect < 0
        and row_bootstrap95[1] is not None
        and row_bootstrap95[1] < 0
        and column_bootstrap95[1] is not None
        and column_bootstrap95[1] < 0
        and negative_fraction is not None
        and negative_fraction >= 0.80
        and _all_axis_medians(row_summary, False)
        and _all_axis_medians(column_summary, False)
        and former_large_worse is not None
        and former_large_worse <= 0.10
        and former_registration_stable
    )
    if candidate_pass and former_pass:
        raise AssertionError("opposed physical-decision gates both passed")
    if not arm_support_gate:
        metric_decision = "metric_safety_gate_fail_arm_support_attrition"
    elif not coverage_gate:
        metric_decision = "physical_coverage_gate_fail"
    elif candidate_pass:
        metric_decision = "candidate_metric_gate_pass_external_gates_pending"
    elif former_pass:
        metric_decision = "former_metric_gate_pass_external_gates_pending"
    else:
        metric_decision = "tie_or_inconclusive_default_decision"

    result = {
        "tool": "score_snap_physical",
        "protocol": {
            "reference_credit": (
                "7jycwjmbfn-eng/pherc0139-physical-audit and Villa PR #1382"
            ),
            "coordinate_mapping": (
                "p_label_local=(p_L0-0.5)/2-[3936,0,0]-[976,1408,1408]"
            ),
            "distance": "per-label-z-plane 2-D EDT to recto_band bit 8",
            "null": "shift every mesh sample +64 L1 voxels in y",
            "pre_distance_max_l1": PRE_DISTANCE_MAX_L1,
            "boundary_poor": "any pre corner with bit 16; descriptive only",
            "weight": "one third pre-snap triangle area per incident vertex",
            "effect": (
                "(iter4-real minus iter0-real) minus "
                "(iter4-null minus iter0-null); positive favors iter0"
            ),
            "spatial_inference": (
                "five y-row and five x-column medians; one-slab consistency only"
            ),
            "row_bootstrap_seed": ROW_BOOTSTRAP_SEED,
            "column_bootstrap_seed": COLUMN_BOOTSTRAP_SEED,
            "bootstrap_resamples": BOOTSTRAP_N,
            "registration_p95_um": REGISTRATION_P95_UM,
            "registration_radius_l1": REGISTRATION_RADIUS_L1,
            "registration_offsets": len(offsets),
        },
        "label_metadata": label_metadata,
        "inputs": {
            "pre": _mesh_record(pre),
            "candidate_iter0": _mesh_record(candidate),
            "former_iter4": _mesh_record(former),
        },
        "invariants": {
            "vertices": len(pre["verts"]),
            "faces": len(pre["faces"]),
            "cubes_present": len(pre["cubes"]),
            "topology_uv_cube_ranges_identical": True,
            "all_mesh_coordinates_finite": True,
            "pre_reference_vertices": int(pre_reference.sum()),
            "physical_target_vertices": int(physical_target.sum()),
            "primary_predeclared_vertices": int(primary_predeclared.sum()),
            "boundary_poor_vertices": int(boundary_poor_predeclared.sum()),
            "arm_support_failures": int(primary_support_fail.sum()),
            "arm_support_gate": arm_support_gate,
            "sampled_distance_evaluations": sampled_evaluations,
            "unsafe_distance_evaluations": unsafe_distance_evaluations,
            "max_finite_sampled_distance_l1": max_finite_sampled_distance,
            "qualified_cubes": len(qualified),
            "global_primary_area_fraction": global_area_fraction,
            "coverage_gate": coverage_gate,
        },
        "gates": {
            "min_vertices_per_cube": MIN_VERTICES_PER_CUBE,
            "min_qualified_cubes": MIN_CUBES,
            "min_cubes_per_row_column": MIN_CUBES_PER_ROW_COLUMN,
            "min_cube_area_fraction": MIN_CUBE_AREA_FRACTION,
            "min_global_area_fraction": MIN_GLOBAL_AREA_FRACTION,
            "effect_min_l1": EFFECT_MIN_L1,
            "effect_min_um": EFFECT_MIN_L1 * L1_VOXEL_UM,
            "positive_or_negative_fraction_min": 0.80,
            "large_worsening_l1": LARGE_WORSENING_L1,
            "large_worsening_fraction_max": 0.10,
        },
        "summary": {
            "median_corrected_effect_l1": median_effect,
            "median_real_effect_l1": median_real_effect,
            "row_cluster_bootstrap95_l1": row_bootstrap95,
            "column_cluster_bootstrap95_l1": column_bootstrap95,
            "positive_fraction": positive_fraction,
            "negative_fraction": negative_fraction,
            "candidate_large_worsening_fraction": candidate_large_worse,
            "former_large_worsening_fraction": former_large_worse,
            "median_unweighted_corrected_effect_l1": (
                float(np.median(unweighted_effects))
                if len(unweighted_effects) else None
            ),
            "boundary_poor": boundary_poor_summary,
            "row_y": row_summary,
            "column_x": column_summary,
            "registration_sensitivity": sensitivity,
            "candidate_registration_stable": candidate_registration_stable,
            "former_registration_stable": former_registration_stable,
            "candidate_gate_pass": candidate_pass,
            "former_gate_pass": former_pass,
            "metric_decision": metric_decision,
            "external_pipeline_safety_gates_required": True,
        },
        "per_cube": per_cube,
    }
    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--pre", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--former", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    result = score(args)
    print(
        json.dumps(
            {
                "invariants": result["invariants"],
                "summary": result["summary"],
            },
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
