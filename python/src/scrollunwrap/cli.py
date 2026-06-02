"""scrollunwrap CLI: stream a remote OME-Zarr ROI -> mesh -> flatten -> tifxyz."""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from dataclasses import asdict
from pathlib import Path

from . import _paths
from .cube_planner import plan_cubes
from .flatten import run_flatboi
from .geometry_report import analyze_welded, build_report, write_report
from .mesher import run_mesher
from .meshprep import (load_welded_obj, prep_component, select_components,
                       write_obj_vf)
from .roi_finder import Roi, auto_pick_rois
from .tifxyz import run_obj2tifxyz
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
    r.add_argument("--cube-timeout", type=float, default=600.0,
                   help="per-cube mesher timeout in seconds (default 600); a "
                        "stalled dense cube is skipped, not allowed to hang")
    r.add_argument("--no-qem", action="store_true")
    r.add_argument("--no-render", action="store_true")
    r.add_argument("--no-anon", action="store_true", help="disable anonymous s3 access")
    r.add_argument("--out", required=True, type=Path)
    r.add_argument("--cube-mesh-bin", default=_paths.default_cube_mesh())
    r.add_argument("--grid-weld-bin", default=_paths.default_grid_weld())
    r.add_argument("--flatboi-bin", default=_paths.default_flatboi())
    r.add_argument("--obj2tifxyz-bin", default=_paths.default_obj2tifxyz())
    return p


def _require_bin(path, name: str):
    if not path or not Path(path).exists():
        sys.exit(f"error: {name} not found at {path!r}; build it or pass the flag.")
    return Path(path)


def _resolve_rois(args) -> list[Roi]:
    if args.bbox is not None:
        return [Roi(bbox_l0=tuple(args.bbox), surface_voxels=-1, density=-1.0,
                    scan_cell=(-1, -1, -1))]
    return auto_pick_rois(args.zarr, args.auto_roi, scan_level=args.scan_level,
                          roi_cubes=tuple(args.roi_cubes), target_level=args.level,
                          anon=not args.no_anon)


def _process_component(args, comp, cdir: Path, bins: dict) -> dict:
    """Flatten one connected component -> its own tifxyz + report + renders."""
    cdir.mkdir(parents=True, exist_ok=True)
    entry: dict = {"dir": str(cdir), "faces": int(len(comp.faces))}
    pp = prep_component(comp, args.level, fill_max_edges=args.fill_max_edges)
    prepped_obj = cdir / "prepped.obj"
    write_obj_vf(prepped_obj, pp.mesh.vertices, pp.mesh.faces)
    flat = run_flatboi(prepped_obj, binary=bins["flatboi"], iters=args.iters,
                       energy=args.energy, tol=args.tol, env=_paths.default_env())
    _, meta, _ = run_obj2tifxyz(flat, cdir / "tifxyz", binary=bins["obj2tifxyz"],
                                step_size=args.step_size)
    crep = build_report(prepped_info=pp.info, flat_obj=flat, tifxyz_meta=meta)
    write_report(cdir / "report.json", crep)
    entry.update({"ok": True, "prepped": pp.info, "uv": crep.get("uv"),
                  "tifxyz_meta": meta})
    if not args.no_render:
        try:
            from . import render
            (cdir / "renders").mkdir(exist_ok=True)
            render.render_mesh_views(comp, cdir / "renders" / "mesh.png")
            render.render_uv_layout(flat, cdir / "renders" / "uv.png")
            render.render_tifxyz_channels(cdir / "tifxyz", cdir / "renders" / "tifxyz.png")
        except Exception as e:  # rendering must never fail the run
            entry["render_warn"] = str(e)
    return entry


def _process_roi(args, roi: Roi, roi_dir: Path, bins: dict) -> dict:
    roi_dir.mkdir(parents=True, exist_ok=True)
    vol = open_volume(args.zarr, args.level, anon=not args.no_anon)
    cubes = plan_cubes(roi.bbox_l0, args.level)

    mesh_res = run_mesher(
        vol, cubes, roi_dir / "mesh",
        cube_mesh_bin=bins["cube_mesh"], grid_weld_bin=bins["grid_weld"],
        halo=args.halo, threshold=args.threshold, skip_qem=args.no_qem,
        max_concurrent=args.max_concurrent, threads_per_cube=args.threads_per_cube,
        cube_timeout=args.cube_timeout)

    welded = load_welded_obj(mesh_res.welded_obj)
    comps = select_components(welded, min_faces=args.min_faces,
                              max_components=args.max_components, mode=args.components)

    roi_report: dict = {
        "roi": asdict(roi),
        "cubes_planned": len(cubes),
        "cubes_meshed": mesh_res.n_obj,
        "weld_report": mesh_res.weld_report,
        "welded": analyze_welded(welded),
        "n_components_selected": len(comps),
        "components": [],
    }
    print(f"  welded: {len(welded.faces)} faces, "
          f"{roi_report['welded']['n_components']} components; "
          f"flattening {len(comps)} (>= {args.min_faces} faces)")
    for ci, comp in enumerate(comps):
        cdir = roi_dir / f"comp_{ci:03d}"
        try:
            entry = _process_component(args, comp, cdir, bins)
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
    bins = {
        "cube_mesh": _require_bin(args.cube_mesh_bin, "cube_mesh"),
        "grid_weld": _require_bin(args.grid_weld_bin, "grid_weld"),
        "flatboi": _require_bin(args.flatboi_bin, "flatboi"),
        "obj2tifxyz": _require_bin(args.obj2tifxyz_bin, "vc_obj2tifxyz_legacy"),
    }
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

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
    total_comps = sum(s.get("components_ok", 0) for s in summary)
    print(f"\nDone: {n_ok}/{len(rois)} ROIs produced surfaces "
          f"({total_comps} tifxyz components total) -> {out/'summary.json'}")
    return 0 if n_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
