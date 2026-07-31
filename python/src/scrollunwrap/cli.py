"""scrollunwrap CLI: stream a remote OME-Zarr ROI -> mesh -> flatten -> tifxyz."""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from dataclasses import asdict
from pathlib import Path

from . import _paths
from .adaptive_bpa import AdaptiveBpaConfig
from .cube_planner import plan_cubes
from .flatten import run_flatboi
from .geometry_report import analyze_welded, build_report, write_report
from .mesher import run_mesher
from .meshprep import (load_welded_obj, prep_component, reorder_welded_to_xyz,
                       select_components, write_obj_vf)
from .roi_finder import Roi, auto_pick_rois
from .tifxyz import run_obj2tifxyz
from .villa_dataset import (classify_patch, import_scroll_hint, init_dataset,
                            register_patch)
from .zarr_source import open_volume


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="scrollunwrap")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run the full pipeline")
    r.add_argument("--zarr", required=True, help="OME-Zarr URI (s3://... or local)")
    r.add_argument("--level", type=int, default=0, help="meshing resolution level (default 0)")
    g = r.add_mutually_exclusive_group()
    g.add_argument("--bbox", type=int, nargs=6, metavar=("Z0", "Z1", "Y0", "Y1", "X0", "X1"),
                   help="explicit level-0 voxel ROI (overrides --auto-roi)")
    g.add_argument("--auto-roi", type=int, default=3, help="auto-pick N dense ROIs (default 3)")
    r.add_argument("--scan-level", type=int, default=5, help="coarse level for ROI scan (default 5)")
    r.add_argument("--roi-cubes", type=int, nargs=3, default=(2, 2, 2), metavar=("NZ", "NY", "NX"))
    r.add_argument("--halo", type=int, default=13)
    r.add_argument("--threshold", default=None, help="e.g. '>=1', '>128' (default: nonzero=surface)")
    r.add_argument("--iters", type=int, default=20,
                   help="flatboi SLIM iterations (default 20; each iteration is a "
                        "sparse solve, so this dominates per-component runtime)")
    r.add_argument("--energy", default="symmetric_dirichlet", choices=["symmetric_dirichlet", "conformal"])
    r.add_argument("--tol", type=float, default=None)
    r.add_argument("--step-size", type=int, default=20)
    r.add_argument("--components", choices=["all", "largest"], default="all",
                   help="flatten every qualifying sheet (default) or only the largest")
    r.add_argument("--min-faces", type=int, default=200,
                   help="skip components smaller than this many faces (default 200)")
    r.add_argument("--max-components", type=int, default=0,
                   help="cap number of components per ROI (0 = unlimited)")
    r.add_argument("--fill-max-edges", type=int, default=30)
    r.add_argument("--max-concurrent", type=int, default=None)
    r.add_argument("--threads-per-cube", type=int, default=1)
    r.add_argument(
        "--mls-backend",
        choices=["auto", "rust-cpu", "cubecl-cpu", "cubecl-hip", "cubecl-cuda"],
        default=None,
        help="MLS implementation selected by a CubeCL-enabled cube_mesh "
             "(default: binary/environment default)",
    )
    r.add_argument("--cube-timeout", type=float, default=600.0,
                   help="per-cube mesher timeout in seconds (default 600; 0 = no "
                        "timeout). A cube exceeding it is skipped, not allowed to hang")
    r.add_argument("--no-qem", action="store_true")
    r.add_argument("--adaptive-bpa", action="store_true",
                   help="audit BPA winding topology and rerun flagged cubes at "
                        "rho 1.10 then 1.00")
    r.add_argument("--umb-y", type=float, default=None,
                   help="scroll umbilicus Y in level-0 voxels (required with --adaptive-bpa)")
    r.add_argument("--umb-x", type=float, default=None,
                   help="scroll umbilicus X in level-0 voxels (required with --adaptive-bpa)")
    r.add_argument("--wrap-pitch", type=float, default=9.5,
                   help="radial voxels per winding for topology audit (default 9.5)")
    r.add_argument("--bpa-rhos", type=float, nargs="+", default=(1.20, 1.10, 1.00),
                   help="candidate BPA radii, highest/preferred first")
    r.add_argument("--bpa-local-span-tol", type=float, default=1.5,
                   help="max winding span of one outer leaf-cube component")
    r.add_argument("--bpa-min-coverage", type=float, default=0.95)
    r.add_argument("--bpa-max-boundary-ratio", type=float, default=1.25)
    r.add_argument("--keep-bpa-candidates", action="store_true",
                   help="retain rejected candidate dumps (large; diagnostic only)")
    r.add_argument("--mesh-only", action="store_true",
                   help="stop after welding: write welded.obj + report and skip "
                        "the per-component flatten/tifxyz stage")
    r.add_argument("--flatboi-timeout", type=float, default=300.0,
                   help="per-component flatboi timeout (s); a component whose SLIM "
                        "solve exceeds this is skipped (default 300)")
    r.add_argument("--no-render", action="store_true")
    r.add_argument("--villa-dataset", type=Path, default=None,
                   help="also emit a validated dataset consumable by Villa's spiral fitter")
    r.add_argument("--villa-classification", choices=["auto", "verified", "unverified"],
                   default="auto", help="dataset patch trust class (default: UV-gated auto)")
    r.add_argument("--s3-anon", choices=["auto", "yes", "no"], default="auto",
                   help="S3 auth: auto = signed when AWS creds (incl. STS) are in "
                        "the environment else anonymous; yes = force anonymous; "
                        "no = force signed/private")
    r.add_argument("--s3-endpoint", default=None,
                   help="custom S3 endpoint URL for non-AWS/private object stores")
    r.add_argument("--out", required=True, type=Path)
    r.add_argument("--cube-mesh-bin", default=_paths.default_cube_mesh())
    r.add_argument("--grid-weld-bin", default=_paths.default_grid_weld())
    r.add_argument("--wind-audit-bin", default=_paths.default_wind_audit())
    r.add_argument("--flatboi-bin", default=_paths.default_flatboi())
    r.add_argument("--obj2tifxyz-bin", default=_paths.default_obj2tifxyz())

    h = sub.add_parser(
        "villa-hint",
        help="import an existing ScrollFiesta tifxyz output as an unverified hint",
    )
    h.add_argument("source", type=Path, help="directory containing x/y/z/mask.tif")
    h.add_argument("--dataset", required=True, type=Path)
    h.add_argument("--id", default=None)
    return p


def _require_bin(path, name: str):
    if not path or not Path(path).exists():
        sys.exit(f"error: {name} not found at {path!r}; build it or pass the flag.")
    return Path(path)


def _s3_opts(args):
    """Resolve (anon, storage_options) from the S3 auth flags."""
    anon = {"auto": None, "yes": True, "no": False}[args.s3_anon]
    so = {"endpoint_url": args.s3_endpoint} if args.s3_endpoint else None
    return anon, so


def _resolve_rois(args) -> list[Roi]:
    if args.bbox is not None:
        return [Roi(bbox_l0=tuple(args.bbox), surface_voxels=-1, density=-1.0,
                    scan_cell=(-1, -1, -1))]
    anon, so = _s3_opts(args)
    return auto_pick_rois(args.zarr, args.auto_roi, scan_level=args.scan_level,
                          roi_cubes=tuple(args.roi_cubes), target_level=args.level,
                          anon=anon, storage_options=so)


def _process_component(args, comp, cdir: Path, patch_id: str, bins: dict) -> dict:
    """Flatten one connected component -> its own tifxyz + report + renders."""
    cdir.mkdir(parents=True, exist_ok=True)
    entry: dict = {"dir": str(cdir), "faces": int(len(comp.faces))}
    pp = prep_component(comp, args.level, fill_max_edges=args.fill_max_edges)
    prepped_obj = cdir / "prepped.obj"
    write_obj_vf(prepped_obj, pp.mesh.vertices, pp.mesh.faces)
    flat = run_flatboi(prepped_obj, binary=bins["flatboi"], iters=args.iters,
                       energy=args.energy, tol=args.tol, env=_paths.default_env(),
                       timeout=args.flatboi_timeout)
    pre_report = build_report(prepped_info=pp.info, flat_obj=flat)
    classification = classify_patch(pre_report, args.villa_classification)
    if args.villa_dataset:
        tifxyz_dir = (Path(args.villa_dataset)
                      / f"{classification}_patches" / patch_id)
    else:
        tifxyz_dir = cdir / "tifxyz"
    _, meta, _ = run_obj2tifxyz(
        flat, tifxyz_dir, binary=bins["obj2tifxyz"],
        step_size=args.step_size, require_mask=bool(args.villa_dataset))
    crep = build_report(prepped_info=pp.info, flat_obj=flat, tifxyz_meta=meta)
    write_report(cdir / "report.json", crep)
    entry.update({"ok": True, "prepped": pp.info, "uv": crep.get("uv"),
                  "tifxyz_meta": meta, "tifxyz": str(tifxyz_dir)})
    if args.villa_dataset:
        entry["villa_patch"] = register_patch(
            Path(args.villa_dataset), patch_id, classification, tifxyz_dir,
            component_report={"faces": entry["faces"], "prepped": pp.info,
                              "uv": crep.get("uv")},
        )
    if not args.no_render:
        try:
            from . import render
            (cdir / "renders").mkdir(exist_ok=True)
            render.render_mesh_views(comp, cdir / "renders" / "mesh.png")
            render.render_uv_layout(flat, cdir / "renders" / "uv.png")
            render.render_tifxyz_channels(tifxyz_dir, cdir / "renders" / "tifxyz.png")
        except Exception as e:  # rendering must never fail the run
            entry["render_warn"] = str(e)
    return entry


def _process_roi(args, roi: Roi, roi_dir: Path, bins: dict) -> dict:
    roi_dir.mkdir(parents=True, exist_ok=True)
    anon, so = _s3_opts(args)
    vol = open_volume(args.zarr, args.level, anon=anon, storage_options=so)
    cubes = plan_cubes(roi.bbox_l0, args.level)

    adaptive = None
    if args.adaptive_bpa:
        adaptive = AdaptiveBpaConfig(
            wind_audit_bin=bins["wind_audit"],
            umb_y=args.umb_y, umb_x=args.umb_x, pitch=args.wrap_pitch,
            rhos=tuple(args.bpa_rhos),
            local_span_tol=args.bpa_local_span_tol,
            min_point_coverage=args.bpa_min_coverage,
            max_boundary_ratio=args.bpa_max_boundary_ratio,
            keep_candidates=args.keep_bpa_candidates,
        )

    mesh_res = run_mesher(
        vol, cubes, roi_dir / "mesh",
        cube_mesh_bin=bins["cube_mesh"], grid_weld_bin=bins["grid_weld"],
        halo=args.halo, threshold=args.threshold, skip_qem=args.no_qem,
        max_concurrent=args.max_concurrent, threads_per_cube=args.threads_per_cube,
        cube_timeout=(None if args.cube_timeout <= 0 else args.cube_timeout),
        adaptive_bpa=adaptive, mls_backend=args.mls_backend)

    welded = load_welded_obj(mesh_res.welded_obj)
    # Viewer-ready copy: (x,y,z) order + per-vertex color (native welded.obj is z,y,x).
    welded_xyz = Path(mesh_res.welded_obj).with_name("welded_xyz.obj")
    reorder_welded_to_xyz(mesh_res.welded_obj, welded_xyz)
    welded_stats = analyze_welded(welded)
    base = {"roi": asdict(roi), "cubes_planned": len(cubes),
            "mls_backend": args.mls_backend or "binary-default",
            "cubes_meshed": mesh_res.n_obj, "weld_report": mesh_res.weld_report,
            "welded": welded_stats}

    if args.mesh_only:
        rep = {**base, "mesh_only": True, "components": []}
        write_report(roi_dir / "report.json", rep)
        print(f"  mesh-only: {mesh_res.welded_obj} "
              f"({welded_stats['n_vertices']} verts, {welded_stats['n_faces']} faces, "
              f"{welded_stats['n_components']} components)")
        return rep

    comps = select_components(welded, min_faces=args.min_faces,
                              max_components=args.max_components, mode=args.components)
    roi_report: dict = {**base, "n_components_selected": len(comps), "components": []}
    print(f"  welded: {len(welded.faces)} faces, "
          f"{roi_report['welded']['n_components']} components; "
          f"flattening {len(comps)} (>= {args.min_faces} faces)")
    for ci, comp in enumerate(comps):
        cdir = roi_dir / f"comp_{ci:03d}"
        patch_id = f"{roi_dir.name}_comp_{ci:03d}"
        try:
            entry = _process_component(args, comp, cdir, patch_id, bins)
            print(f"    comp {ci}: {entry['faces']} faces -> tifxyz "
                  f"(flipped UV {entry['uv']['flipped_uv_triangles']}/{entry['uv']['n_faces']})")
        except Exception as e:
            entry = {"dir": str(cdir), "faces": int(len(comp.faces)), "ok": False,
                     "error": str(e)}
            print(f"    comp {ci}: {entry['faces']} faces -> skipped: {e}")
        roi_report["components"].append(entry)

    write_report(roi_dir / "report.json", roi_report)
    return roi_report


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    if args.cmd == "villa-hint":
        patch_id = args.id or args.source.name
        entry = import_scroll_hint(
            args.dataset.resolve(), args.source.resolve(), patch_id)
        print(json.dumps(entry, indent=2))
        return 0

    if args.adaptive_bpa and (args.umb_y is None or args.umb_x is None):
        sys.exit("error: --adaptive-bpa requires --umb-y and --umb-x")
    if args.adaptive_bpa and (not args.bpa_rhos or any(r <= 0 for r in args.bpa_rhos)):
        sys.exit("error: --bpa-rhos must contain positive radii")
    bins = {
        "cube_mesh": _require_bin(args.cube_mesh_bin, "cube_mesh"),
        "grid_weld": _require_bin(args.grid_weld_bin, "grid_weld"),
    }
    if args.adaptive_bpa:
        bins["wind_audit"] = _require_bin(args.wind_audit_bin, "wind_audit")
    if not args.mesh_only:
        # flatten/tifxyz binaries are only exercised past the weld stage
        bins["flatboi"] = _require_bin(args.flatboi_bin, "flatboi")
        bins["obj2tifxyz"] = _require_bin(args.obj2tifxyz_bin,
                                          "vc_obj2tifxyz_legacy")
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    if args.villa_dataset:
        args.villa_dataset = init_dataset(Path(args.villa_dataset).resolve())
        print(f"Villa dataset -> {args.villa_dataset}")

    rois = _resolve_rois(args)
    if not rois:
        sys.exit("error: no ROIs found (region may be empty); try --bbox.")
    (out / "rois.json").write_text(json.dumps([asdict(r) for r in rois], indent=2))
    print(f"ROIs: {len(rois)} -> {out/'rois.json'}")

    summary = []
    for i, roi in enumerate(rois):
        roi_dir = out / f"roi_{i:03d}"
        print(f"\n=== ROI {i} bbox(l0)={roi.bbox_l0} density={roi.density:.3f} ===")
        try:
            rep = _process_roi(args, roi, roi_dir, bins)
            if rep.get("mesh_only"):
                w = rep["welded"]
                ok = w["n_faces"] > 0
                summary.append({"roi": i, "ok": ok, "dir": str(roi_dir),
                                "mesh_only": True, "verts": w["n_vertices"],
                                "faces": w["n_faces"], "components": w["n_components"]})
                print(f"  mesh -> {roi_dir} ({w['n_vertices']} verts, {w['n_faces']} faces)")
            else:
                comp_ok = sum(bool(c.get("ok")) for c in rep["components"])
                summary.append({"roi": i, "ok": comp_ok > 0, "dir": str(roi_dir),
                                "components_ok": comp_ok,
                                "components_total": len(rep["components"])})
                print(f"  done -> {roi_dir} ({comp_ok}/{len(rep['components'])} components)")
        except Exception as e:
            traceback.print_exc()
            summary.append({"roi": i, "ok": False, "dir": str(roi_dir), "error": str(e)})
            print(f"  [fail] ROI {i}: {e}")

    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    n_ok = sum(s["ok"] for s in summary)
    msg = f"\nDone: {n_ok}/{len(rois)} ROIs OK"
    if any(not s.get("mesh_only") for s in summary):
        total_comps = sum(s.get("components_ok", 0) for s in summary)
        msg += f" ({total_comps} tifxyz components total)"
    print(msg + f" -> {out/'summary.json'}")
    return 0 if n_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
