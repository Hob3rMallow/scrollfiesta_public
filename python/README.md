# scrollunwrap

Streams a remote **OME-Zarr** surface-prediction volume over `s3://`, feeds a
region of interest **directly in memory** to the scrollfiesta mesher (no
intermediate TIFF cubes), flattens the welded mesh with villa's C++ `flatboi`
(SLIM), and emits a **tifxyz** surface — with geometry/topology checks and PNG
renders for validation.

```
scrollunwrap run \
  --zarr s3://vesuvius-challenge-open-data/PHerc0139/.../<...>.zarr --level 0 \
  --auto-roi 3 \
  --out RUN_DIR \
  --cube-mesh-bin .../scrollfiesta_public/src/cube_mesh \
  --grid-weld-bin .../scrollfiesta_public/src/grid_weld \
  --flatboi-bin .../volume-cartographer/build-macos/bin/flatboi \
  --obj2tifxyz-bin .../volume-cartographer/build-macos/bin/vc_obj2tifxyz_legacy
```

## Private S3

`--s3-anon auto` (default) signs requests when AWS credentials are present in the
environment, else falls back to anonymous. Credentials are read from the standard
AWS environment (incl. temporary **STS** sessions) by s3fs/botocore — never from
code or flags:

```
export AWS_ACCESS_KEY_ID=... AWS_SECRET_ACCESS_KEY=... AWS_SESSION_TOKEN=...
scrollunwrap run --zarr s3://volpkgs/.../<...>.zarr --auto-roi 1 --out RUN
```

Use `--s3-anon no` to force signed access, `--s3-anon yes` to force anonymous, and
`--s3-endpoint URL` for non-AWS/private object stores (e.g. MinIO).

Dev:

```
uv sync
uv run pytest -m "not network and not requires_binary"   # pure-unit, no deps
```

