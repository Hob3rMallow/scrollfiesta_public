ScrollFiesta -- virtual meshing, registration, and unrolling
=============================================================

ScrollFiesta turns segmented CT surface predictions into registered scroll
meshes and textured unrolls. The supported end-to-end path is:

  grid_pipeline  ->  scroll_whole  ->  scroll_unroll

The default path uses the current CVT/RVD simplifier, keeps cube-face geometry
for whole-scroll registration, performs global placement/re-registration, then
runs weld, overlap ownership, CT surface snap, final relaxation, and RAW texture
baking. The retired experimental scroll_unwrap, scroll_atlas, and scroll_ribbon
front ends are not part of this package.


PRIMARY EXECUTABLES
-------------------
Binaries land in build\Release\ on Windows:

  cube_mesh       Mesh one 128^3 prediction cube.
  grid_pipeline   Normal entry point for a prediction grid. Runs cube_mesh in
                  parallel and writes resumable per-cube meshes and logs.
  grid_weld       Optional direct weld of already-meshed cubes.
  scroll_whole    Register and place every usable cube into one scroll frame;
                  supports global re-registration and seam auditing.
  scroll_unroll   Produce the complete unroll, CT snap, texture TIFs, PNG
                  previews, diagnostic OBJs, and JSON metrics.
  atlas_overlap_fix
                  Optional registered-atlas overlap repair/audit tool.

The Visual Studio solution and the top-level CMake build select these supported
front ends by default. Older prototype drivers and their exclusive flattening
modules have been removed; developer diagnostics remain available as individual
projects where they still exercise live code.


BUILDING
--------
Windows, Visual Studio 2022:

    .\build-deps.ps1
    msbuild scrollfiesta.sln -p:Configuration=Release -p:Platform=x64 -m

CMake (Windows or Linux):

    cmake -S . -B build -DSCROLLFIESTA_BUILD_TOOLS=ON \
          -DSCROLLFIESTA_BUILD_TESTS=ON
    cmake --build build --config Release --parallel
    ctest --test-dir build -C Release --output-on-failure

The low-level src/Makefile remains available for core GCC builds, but CMake is
the portable build of record for the full supported toolchain.


INPUT GRID
----------
A prepared grid has matching 128 x 128 x 128 multipage uint8 TIFF cubes:

    <grid>/cubes_PRED/z#####_y#####_x#####.tif
    <grid>/cubes_RAW/z#####_y#####_x#####.tif

PRED is the binary recto-surface mask. RAW is the grayscale CT volume used by
the snap and texture bake. The filename is the cube's world-space (z,y,x)
origin and adjacent cubes differ by 128. python/scripts/carve_grid_tifs.py can
carve both trees from OME-Zarr sources and rejects an accidentally all-zero RAW
carve.


END-TO-END USAGE
----------------
1. Mesh the grid. CVT/RVD simplification and a zero trim inset are the
   defaults, so steps 1-5 run as written with no extra flags:

    build\Release\grid_pipeline.exe GRID output\run --halo 13

   Important outputs include output\run\dump\, output\run\logs\,
   pipeline_summary.csv, and the optional welded.obj. Use --skip-existing to
   resume. The zero inset matters: step 2 registers cubes by pairing skins
   ACROSS cube seams, so charts have to reach the cube faces. Insetting them
   by a voxel halves the seam evidence on PHerc0139-4x5x5 (2,713 pairs against
   6,441) and step 2's audit then fails at 15.78% whole-turn error instead of
   passing at 2.98%. --trim-inset 1 and --simplify qem are explicit legacy
   choices for a weld-only run, where the seam bridge does want a gap.

   To arm the per-cube multi-winding BPA growth gate as well, add the scroll
   geometry: --umb-y 3405 --umb-x 2878 --wrap-pitch 9.5 for PHerc0139. This
   also arms grid_weld's seam winding/phase gates. It is off by default
   because it is still being characterised: on the 4x5x5 fixture it splits a
   representative cube from 10 pieces into 68.

2. Register and place the per-cube sheets:

    build\Release\scroll_whole.exe output\run\dump output\run_placed
    build\Release\scroll_whole.exe output\run_placed --reregister --audit

   The placed directory contains one *_placed.obj plus UV/skin/group records per
   usable cube, placed_index.json, audit.json, and logs.

3. Unroll, snap to the CT surface, and bake texture data. Stages 12345 are the
   default and mean base bake, seam join/relax, overlap ownership, snap, and
   final relax:

    build\Release\scroll_unroll.exe output\run_placed output\run_unroll \
        --raw GRID\cubes_RAW --steps 12345 --id run

   The unroll directory contains, for every stage:

     *_rawtex.tif             full-resolution grayscale texture
     *_rawtex_strip.tif/png   tiled readable strip
     *_rawtex_preview.png     downsampled preview
     *_diagclass.tif/png      validity/ownership diagnostics
     *_xyzmap*.png            world-coordinate diagnostics
     *_seamzoom.png           seam inspection image
     *_stats.json             stage metrics

   It also contains final scroll/leftover OBJs and pipeline_stats.json. Optional
   --export-tifxyz and --export-atlas outputs are documented by
   scroll_unroll --help.

   `scripts/run_canonical_grid.ps1` runs this five-stage review alongside the
   atlas/ribbon path. A canonical run is considered incomplete unless the
   step4 snapped and step5 final-relaxed full TIFs, readable strips, preview
   PNGs, and metrics JSON files exist under `snap_relax/`. They are linked from
   `CANONICAL_OUTPUTS.md`, and their paths are recorded in `logs/SUMMARY.txt`.

4. Resolve atlas overlaps. Registration fixes where each chart sits WITHIN a
   winding; this decides WHICH winding it sits on. --rounds alternates a tabu
   search over the integer winding depths with a continuous per-chart relayout,
   feeding each round's field to the next and stopping early once a round moves
   no charts:

    build\Release\atlas_overlap_fix.exe output\run_placed output\run_atlas \
        --rounds 7

   Outputs before/after UV OBJs, per-round winding-layer OBJs, charts.csv,
   groups.csv, diag_before/after.png, overlap_fix_stats.json, and
   atlas_solution.bin -- the checkpoint the ribbon stage consumes.

5. Fit the ribbon. This is the terminal stage of the whole-scroll path: it
   slices the registered atlas into constant-v rows and fits one collision-free
   ribbon through them, auditing overlap exactly rather than on a grid:

    build\Release\atlas_ribbon_fit.exe output\run_placed \
        output\run_atlas\atlas_solution.bin output\run_ribbon --mode ribbon

   Outputs ribbon.obj, ribbon_rows_world.obj, observations_world.obj,
   u_overlap_audit.json, the u_* collision/registration ledgers, and
   ribbon_fit_stats.json. --mode also selects the intermediate stages
   (f0/observations .. f6/collision) for ablation, and --raw DIR additionally
   audits adjacent-row RAW texture phase.

Run each primary tool with --selftest for its embedded module suite.


OUTPUT CONVENTIONS
------------------
CLI OBJ vertices are written in internal (z,y,x) order. The C API uses ordinary
(x,y,z) coordinates. Texture rasters use unroll (u,v); their JSON files record
canvas origin, dimensions, fill, multi-coverage, darkness, and seam metrics.
Large validation products belong under output/ and are intentionally ignored by
git.

USING SCROLLFIESTA AS A LIBRARY (C API)
---------------------------------------
Everything the pipeline can do is also exposed as a C library behind ONE
header, include/scrollfiesta.h: mesh cleanup (manifold repair, pinhole/hole
fill, sliver cleanup, component cull), the detangle/split family (depth peel,
developability cut, bridge cut, overlap separation), QEM decimation,
developability fairing, orientation, Seamster cut-to-disk / handle severing,
BPA reconstruction, MLS/LOP projection, and the full in-memory volume->mesh
pipeline. See the header for the full reference; the important conventions:

  - Coordinates at the API boundary are (x,y,z) voxel units. (Internally the
    code is (z,y,x), and OBJ files written by the CLI TOOLS are "v z y x" --
    those are NOT interchangeable with the true-xyz OBJs of the API.)
  - Everything returns sf_status; no exception/longjmp/signal crosses the
    boundary. Long operations take a progress callback; returning nonzero
    from it cancels the call (SF_CANCELLED, outputs untouched).
  - Outputs carry per-vertex provenance (sf_mesh.vmap: output vertex ->
    input vertex, -1 for newly created) wherever the operation permits.
  - Library-filled buffers are released with sf_mesh_free / sf_mesh_list_free
    / sf_free -- never with your own allocator.

Two ways to consume it:

1) Runtime loading (recommended for applications -- how VC3D binds):
       h = LoadLibrary/dlopen("scrollfiesta.dll" / "libscrollfiesta.so");
       sf_get_api_fn get_api = (sf_get_api_fn)GetProcAddress/dlsym(h, "sf_get_api");
       const sf_api *sf = get_api(SCROLLFIESTA_ABI_VERSION);   /* NULL = incompatible */
       ... call everything through the sf-> function-pointer table ...
   The table is append-only within an ABI version, and all memory crosses
   the boundary through the table's malloc/free -- so the library and the
   host can be built with different compilers/CRTs, and upgrading is
   "swap the DLL, restart". tests/api_dlopen.c is a complete example.

2) Build-time linking via CMake:
       FetchContent_Declare(scrollfiesta GIT_REPOSITORY <this repo> GIT_TAG vX.Y.Z)
       set(SCROLLFIESTA_BUILD_TOOLS OFF)
       set(SCROLLFIESTA_WITH_TIFF   OFF)   # library core needs no TIFF
       FetchContent_MakeAvailable(scrollfiesta)
       target_link_libraries(your_app PRIVATE scrollfiesta::scrollfiesta)
   Options (defaults: everything ON when top-level, library-only when
   embedded): SCROLLFIESTA_BUILD_TOOLS, SCROLLFIESTA_BUILD_TESTS,
   SCROLLFIESTA_WITH_TIFF, SCROLLFIESTA_OPENMP, SCROLLFIESTA_INSTALL.
   The vendored deps (Triangle, Clipper2) build from deps/src/ automatically;
   no system packages are required for the library core. NB: the project
   never calls find_package(OpenMP) and no library TU includes <omp.h>, so
   hosts that stub OpenMP (VC3D) embed cleanly.

Threading: any single operation may be called from any thread, but at most
one ScrollFiesta operation runs at a time per process (internal mutex);
operations parallelize internally per sf_common_opts.n_threads.

Licensing: ScrollFiesta's own code is MIT; the vendored dependencies keep
their own licenses -- notably Shewchuk's Triangle (hole filling) restricts
commercial distribution. See THIRD_PARTY_LICENSES.md.


NOTES
-----
- CPU by default: every stage runs on the CPU with no GPU required. The MLS
  projection additionally has an optional GPU backend (the CubeCL CUDA/HIP
  crate in deps/scrollfiesta-mls-cubecl, built with
  scripts/build_mls_cubecl.ps1 and selected at runtime with MLS_BACKEND);
  CPU and GPU backends feed the same geometry and topology gates. The U-Net
  that produces the input predictions is a separate, upstream tool and is not
  part of ScrollFiesta.
- Topology is the priority: the weld is tuned to never merge distinct scroll
  wraps and never split a single sheet, even at the cost of leaving a seam gap.
- Each executable prints its usage when run with no arguments.
- Architecture and design rationale: submission.pdf (LaTeX source and figures
  in submission_update/; superseded submissions in previous_submissions/)
- Example outputs: sample_outputs/, and a full worked run in output/

