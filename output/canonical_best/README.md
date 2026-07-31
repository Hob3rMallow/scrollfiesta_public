# Canonical best pre-ribbon results

This directory freezes the atlas-stage results validated on 2026-07-31. It
contains no ribbon output. The large grids were meshed through the CubeCL CUDA
path and solved with bounded Gauss-Seidel atlas tiles no larger than 5x5x5.

## Start here

| Grid | Textured atlas | Flat atlas | World mesh with UV | Flat weld |
|---|---|---|---|---|
| 4x5x5 | [`texture/atlas_folded_preview.png`](pherc0139_4x5x5/texture/atlas_folded_preview.png) | [`atlas/after_uv.obj`](pherc0139_4x5x5/atlas/after_uv.obj) | [`atlas/atlas_bake.obj`](pherc0139_4x5x5/atlas/atlas_bake.obj) | [`weld/welded.obj`](pherc0139_4x5x5/weld/welded.obj) |
| 10x10x10 | [`texture/atlas_folded_preview.png`](pherc0139_10x10x10/texture/atlas_folded_preview.png) | [`atlas/after_uv.obj`](pherc0139_10x10x10/atlas/after_uv.obj) | [`atlas/atlas_bake.obj`](pherc0139_10x10x10/atlas/atlas_bake.obj) | intentionally omitted |
| 4x21x21 | [`texture/atlas_folded_preview.png`](pherc0139_4x21x21/texture/atlas_folded_preview.png) | [`atlas/after_uv.obj`](pherc0139_4x21x21/atlas/after_uv.obj) | [`atlas/atlas_bake.obj`](pherc0139_4x21x21/atlas/atlas_bake.obj) | intentionally omitted |

The 4x5x5 `after_uv.obj`, `atlas_bake.obj`, and `groups_xyz.obj` are
byte-for-byte identical to private checkpoint
`overlap_fix_ablation/58_ribbon_checkpoint`. Its texture was freshly rebaked
from that exact `atlas_bake.obj`, not copied from an earlier visual probe.

## Artifact contract

- `atlas/after_uv.obj` is the solved flat atlas.
- `atlas/atlas_bake.obj` is the giant world-space triangle mesh with solved
  texture coordinates.
- `atlas/groups_xyz.obj` is the giant world-space mesh colored by atlas group.
- `atlas/atlas_solution.bin` is the checkpoint for a later ribbon pass.
- `texture/atlas_rawtex.{tif,png}` is the full atlas-domain CT readback.
- `texture/atlas_rawtex_strip.{tif,png}` is the same texture folded into rows.
- `texture/atlas_folded_preview.png` is a compact two-column review image.
- `run/` preserves solver manifests, placement audits, and meshing provenance.
- `SHA256SUMS.txt` covers every committed artifact in this directory.

The 10x10x10 and 4x21x21 grids deliberately do not contain a flat whole-grid
weld. That operation is memory-infeasible at these sizes and is not consumed
by the atlas path. Their `atlas_bake.obj` and `groups_xyz.obj` files are the
canonical giant meshes. The 4x5x5 weld is retained because the bounded grid
makes it practical and useful to inspect.

All OBJ, BIN, TIFF, and PNG payloads in this tree are tracked with Git LFS.
