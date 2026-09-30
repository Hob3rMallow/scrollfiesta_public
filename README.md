# ScrollFiesta

ScrollFiesta is a C/C++ research pipeline for turning segmented CT surface
predictions into registered scroll meshes, textured unrolls, and review
artifacts.

This repository is the **geometry** half of the work. It takes a completed
surface-probability volume as its input and does not produce or repair one.
Building that volume — cross-resolution distillation of a native-fine predictor
into the coarse model, and additive completion of short dropouts in the
prediction — lives in the companion **Socratic Method** project at
<https://github.com/Hob3rMallow/socratic_method>.

The supported end-to-end path is:

```text
grid_pipeline -> scroll_whole -> scroll_unroll
```

An optional atlas branch continues through `atlas_overlap_fix` and
`atlas_ribbon_fit`. The default path runs on CPU; the per-cube MLS projection
also has an optional CubeCL CUDA/HIP backend.

Since September 2026 the per-cube meshes can also be assembled directly into
flattened, CT-textured sheets by the native sheet assembler:

```text
grid_pipeline -> sheet_assemble -> sheet_strip / sheet_layers / vmesh_tifxyz
```

`sheet_assemble` flattens every chart, joins charts across cube seams, places
them in one winding-aware parameter domain, repairs overlap and distortion,
audits the result against the source meshes, and bakes and scores the sheet.
See [docs/SHEET_ASSEMBLY.md](docs/SHEET_ASSEMBLY.md).

This repository is under active development. The public C API is versioned,
but pipeline outputs and experimental command-line options may still change.
For the updated research write-up, see [submission.pdf](submission.pdf).
The [September release check](submission_update/release/README.md) includes
the PHerc0139 4×5×5 CT sheet, audit, and commands. For
curated example results, start with
[output/canonical_best/README.md](output/canonical_best/README.md). The August 2026 exhibit and its claim ledger are indexed in [submission_update/gallery/README.md](submission_update/gallery/README.md).

## What the pipeline produces

- resumable per-cube meshes and diagnostics;
- global cube registration, placement, and seam audits;
- UV unrolls with overlap ownership and surface snapping;
- full-resolution RAW texture TIFFs and readable PNG previews;
- OBJ geometry plus JSON/CSV metrics for review and regression checks.

Retired experimental `scroll_unwrap`, `scroll_atlas`, and `scroll_ribbon`
front ends are not part of the supported package.

## Build and test

The portable build requires CMake 3.24 or newer and C17/C++17 compilers.
Top-level builds enable the library, command-line tools, tests, TIFF support,
and install rules by default.

On Windows, first stage the vendored libtiff and zlib libraries, then build
with either CMake or Visual Studio 2022:

```powershell
.\build-deps.ps1
cmake -S . -B build -A x64 `
  -DSCROLLFIESTA_BUILD_TOOLS=ON `
  -DSCROLLFIESTA_BUILD_TESTS=ON
cmake --build build --config Release --parallel 4
ctest --test-dir build -C Release --output-on-failure
```

The checked-in solution remains available for native MSBuild workflows:

```powershell
.\build-deps.ps1
msbuild scrollfiesta.sln -p:Configuration=Release -p:Platform=x64 -m
```

On Linux, install a C/C++ toolchain, libtiff and a BLAS/LAPACK
(`liblapack-dev`), then run:

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DSCROLLFIESTA_BUILD_TOOLS=ON \
  -DSCROLLFIESTA_BUILD_TESTS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

For a dependency-light library and API test build, disable the CLI and TIFF:

```bash
cmake -S . -B build \
  -DSCROLLFIESTA_BUILD_TOOLS=OFF \
  -DSCROLLFIESTA_WITH_TIFF=OFF \
  -DSCROLLFIESTA_BUILD_TESTS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

`scroll_whole` and `sheet_assemble` need a BLAS/LAPACK, which backs the TAUCS
Cholesky in their ribbon, metric and repair solves. The build takes the first
it finds: libraries staged by `build-deps`, then a system LAPACK, then the
vendored CLAPACK source under `deps/clapack/` (correct anywhere, but roughly a
minute of extra compilation). If none is available, `scroll_whole` and
`sheet_assemble` are skipped and the rest of the toolchain still builds.
Library-only and embedded builds never look for one.

The low-level [`src/Makefile`](src/Makefile) remains available for the core
GCC per-cube/grid build. CMake is the build of record for the full supported
toolchain and is what CI exercises.

## Supported executables

| Executable | Purpose |
|---|---|
| `cube_mesh` | Mesh one 128³ prediction cube. |
| `grid_pipeline` | Mesh a prediction grid in parallel and write resumable per-cube output. |
| `grid_weld` | Directly weld already-meshed cubes when a flat weld is useful. |
| `scroll_whole` | Register and place usable cubes in one scroll frame. |
| `scroll_unroll` | Unroll, snap, relax, and bake textures and diagnostics. |
| `atlas_overlap_fix` | Repair and audit winding ownership in a registered atlas. |
| `atlas_ribbon_fit` | Fit the terminal collision-free ribbon through a solved atlas. |
| `sheet_assemble` | Assemble a per-cube mesh pile into audited, flattened sheets (September 2026). |
| `obj_bake_raw` | Bake RAW CT texture through a flattened mesh (used by `sheet_assemble`). |
| `ribbon_verdict` | Score a flattened sheet against its acceptance gates (used by `sheet_assemble`). |
| `sheet_strip` | Flow a sheet review's bake into strip pages in spiral order, small pieces in a tray. |
| `sheet_reading` | Flip mirror-imaged pieces of a sheet image or strip page into reading handedness. |
| `sheet_layers` | Render the surface-volume layers along the inward normal on the sheet grid, for ink detection. |
| `vmesh_tifxyz` | Export flattened pieces as VC3D tifxyz segments. |
| `scroll_axis_track` | Derive the scroll axis (umbilicus) table from a mesh pile. |
| `swirl_slim` | Cut radius-gated bridges in a marbled VMESH region and SLIM-refine each freed piece. |

Small diagnostic tools are also built when `SCROLLFIESTA_BUILD_TOOLS=ON`.
Run a tool with no arguments for usage, or with `--selftest` where supported.

## Input grid

A prepared grid contains matching 128 × 128 × 128 multipage `uint8` TIFF
cubes:

```text
<grid>/cubes_PRED/z#####_y#####_x#####.tif
<grid>/cubes_RAW/z#####_y#####_x#####.tif
```

`PRED` is the binary recto-surface mask. `RAW` is the grayscale CT volume used
for snapping and texture baking. Filenames encode the cube's world-space
`(z,y,x)` origin; adjacent cubes differ by 128.

[`python/scripts/carve_grid_tifs.py`](python/scripts/carve_grid_tifs.py) can
carve both trees from OME-Zarr sources and rejects an accidentally all-zero RAW
carve. See [python/README.md](python/README.md) for that optional workflow.

## Minimal end-to-end run

With Windows release binaries and a prepared `GRID` directory:

```powershell
build\Release\grid_pipeline.exe GRID output\run --halo 13
build\Release\scroll_whole.exe output\run\dump output\run_placed
build\Release\scroll_whole.exe output\run_placed --reregister --audit
build\Release\scroll_unroll.exe output\run_placed output\run_unroll `
  --raw GRID\cubes_RAW --steps 12345 --id run
```

Stages `12345` mean base bake, seam join/relax, overlap ownership, CT surface
snap, and final relaxation. Large generated data belongs under `output/` and is
ignored unless deliberately curated into `output/canonical_best/`.

The detailed stages, expected artifacts, overlap repair, ribbon fitting, and
review requirements are documented in
[docs/PIPELINE.md](docs/PIPELINE.md).

## Public configuration

[`configs/default.json`](configs/default.json) is the single configuration for
the native sheet assembler and optional quadribbon tool. It uses PHerc0139
4×5×5 calibration and the supporting [`configs/axis.csv`](configs/axis.csv).
Set `raw_source.path`, `bake.exe`, and `verdict.exe` for your local data and
build; replace the geometry calibration and axis table for other datasets.
Paths are relative to the working directory.

The staged sheet workflow and its reference commands are documented in
[docs/SHEET_ASSEMBLY.md](docs/SHEET_ASSEMBLY.md). On the saved 100-cube fixture,
this configuration completes through audit but is **not geometry-qualified**
(47.4% area coverage, 171 overlapping triangle pairs). It is a starting
configuration, not a validated reproduction of the curated atlas result.

## C library

The supported public API is declared in one header:
[`include/scrollfiesta.h`](include/scrollfiesta.h). It exposes mesh cleanup,
split/repair operations, CVT/RVD remeshing, developability fairing, topology
surgery, BPA reconstruction, MLS/LOP projection, and the in-memory
volume-to-mesh pipeline.

Two consumption models are supported:

1. Runtime loading through `sf_get_api`, which keeps allocator and compiler/CRT
   boundaries explicit. [`tests/api_dlopen.c`](tests/api_dlopen.c) is the
   complete example.
2. Build-time linking with CMake:

   ```cmake
   FetchContent_Declare(
     scrollfiesta
     GIT_REPOSITORY <this-repository>
     GIT_TAG vX.Y.Z
   )
   set(SCROLLFIESTA_BUILD_TOOLS OFF)
   set(SCROLLFIESTA_WITH_TIFF OFF)
   FetchContent_MakeAvailable(scrollfiesta)
   target_link_libraries(your_app PRIVATE scrollfiesta::scrollfiesta)
   ```

The main build switches are `SCROLLFIESTA_BUILD_TOOLS`,
`SCROLLFIESTA_BUILD_TESTS`, `SCROLLFIESTA_WITH_TIFF`,
`SCROLLFIESTA_OPENMP`, and `SCROLLFIESTA_INSTALL`. Top-level builds default
them on; embedded builds default to the library-only shape.

API coordinates use ordinary `(x,y,z)` voxel units. Internally the pipeline
uses `(z,y,x)`, and CLI-written OBJ vertex columns follow that internal order.
Do not interchange those OBJ files with true-XYZ API meshes without converting
the axes and winding as documented in the public header.

A ScrollFiesta operation may be called from any thread, but the library
serializes operations per process. Individual operations can parallelize
internally through `sf_common_opts.n_threads`.

## Repository map

| Path | Contents |
|---|---|
| `include/` | Versioned public C API. |
| `src/` | Core algorithms and command-line front ends. |
| `configs/` | One public configuration and its PHerc0139 reference axis. |
| `tests/` | Direct-link, runtime-load, and embedding smoke tests. |
| `scripts/` | Build helpers, data conversion, diagnostics, and module tests. |
| `python/` | Optional data preparation and ROI interoperability tools. |
| `sample_outputs/` | Small public examples. |
| `output/canonical_best/` | Curated Git LFS result bundle and checksums. |
| `submission_update/` | LaTeX source and figures for the current paper. |
| `deps/` | Vendored dependencies and their upstream notices. |

The documentation map in [docs/README.md](docs/README.md) distinguishes current
interfaces from historical implementation notes.

## Contributing and merge readiness

See [CONTRIBUTING.md](CONTRIBUTING.md) for build profiles, validation
expectations, artifact rules, and the review checklist used while changes are
being integrated for the submission.

## License

ScrollFiesta's own code is MIT-licensed. Vendored dependencies retain their
own terms. In particular, GCoptimization and Triangle carry redistribution or
commercial-use restrictions. Review [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md)
before distributing a binary or derived package.
