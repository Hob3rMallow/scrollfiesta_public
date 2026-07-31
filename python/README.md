# ScrollFiesta Python utilities

The supported submission workflow is the native pipeline documented in the
repository `README.txt`:

```text
grid_pipeline -> scroll_whole -> scroll_unroll
```

It defaults to CVT/RVD meshing and produces the registered scroll OBJs plus
full-resolution RAW texture TIFs and PNG previews. The Python package is an
optional data-preparation and small-ROI interoperability layer; it is not the
default whole-scroll flattener.

## Prepare a grid from OME-Zarr

`carve_grid_tifs.py` writes the exact `cubes_PRED/` and `cubes_RAW/` layout the
native tools consume. RAW cubes are always emitted because snapping and texture
baking sample outside the predicted surface. The command fails if every carved
RAW cube is zero, catching an incorrect Zarr level or mirror before a long run.

```powershell
py -m uv run --project python python/scripts/carve_grid_tifs.py `
  --pred-zarr s3://.../prediction.zarr `
  --raw-zarr s3://.../volume.zarr `
  --bbox Z0 Z1 Y0 Y1 X0 X1 `
  --umbilicus Y X --out GRID
```

Then run the native flow from the repository root:

```powershell
build\Release\grid_pipeline.exe GRID output\run --halo 13
build\Release\scroll_whole.exe output\run\dump output\run_placed
build\Release\scroll_whole.exe output\run_placed --reregister --audit
build\Release\scroll_unroll.exe output\run_placed output\run_unroll `
  --raw GRID\cubes_RAW --steps 12345 --id run
```

## Optional streamed ROI runner

`scrollunwrap run` streams an OME-Zarr prediction ROI, calls the ScrollFiesta
mesher, and can hand the welded OBJ to external Villa/Volume Cartographer tools.
Use it for focused experiments or interoperability, not as the canonical
whole-scroll entry point.

```text
scrollunwrap run \
  --zarr s3://.../prediction.zarr --level 0 --auto-roi 3 --out RUN_DIR \
  --cube-mesh-bin .../cube_mesh --grid-weld-bin .../grid_weld \
  --flatboi-bin .../flatboi --obj2tifxyz-bin .../vc_obj2tifxyz_legacy
```

For an adaptive BPA audit, provide the scroll umbilicus:

```text
scrollunwrap run ... \
  --adaptive-bpa --umb-y 3405 --umb-x 2878 --wrap-pitch 9.5 \
  --villa-dataset RUN/scrollfiesta_dataset
```

A flagged 1.20-radius result is retried at 1.10 and 1.00. A candidate is selected
only when it is topology-clean, retains at least 95% of baseline BPA points
within 1.5 voxels, and has no more than 1.25 times the baseline boundary edges.
The decision is recorded in each cube's `adaptive_bpa.json` and in
`adaptive_bpa_summary.json`.

An existing tifxyz export can be registered as a low-weight Villa hint:

```text
scrollunwrap villa-hint RUN/tifxyz_patch \
  --dataset RUN/scrollfiesta_dataset --id stage1_patch
```

Hints always go to `unverified_patches`; the command intentionally has no
verified-classification option.

## S3 and development

`--s3-anon auto` signs requests when standard AWS credentials are present and
otherwise uses anonymous access. `yes` forces anonymous access, `no` forces
signed access, and `--s3-endpoint URL` supports compatible object stores.
Credentials are read only through the standard AWS environment/provider chain.

```text
uv sync
uv run pytest -m "not network and not requires_binary"
```