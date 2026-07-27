ScrollFiesta -- virtual meshing & unwrapping for the Herculaneum papyri
======================================================================

ScrollFiesta turns segmented (binary) CT voxel data of carbonized scroll
fragments into a clean, single-sided, orientable triangle mesh suitable for
virtual unwrapping. Input is a per-cube binary recto-surface mask (e.g. an
nnU-Net prediction); output is a welded .obj surface mesh stitched across a
grid of cubes. The method is point-cloud -> ball-pivoting, run per cube and
then welded across the grid. It is CPU-only.


WHAT YOU GET (after building)
-----------------------------
Three executables (in build\Release\ on Windows, src/ on Linux):

  cube_mesh      Per-cube mesher. Extracts one cube's prediction TIFF into
                  single-sided sheet meshes: voxel-centres -> MLS smoothing ->
                  Hoppe normal orientation -> ball-pivoting -> sheet split ->
                  hole fill -> QEM decimation -> trim.

  grid_pipeline   Orchestrator -- the normal entry point. Scans a grid of cube
                  TIFFs, runs cube_mesh on each (concurrently), then runs
                  grid_weld to stitch them.

  grid_weld       Stitcher. Concatenates the per-cube meshes at their world
                  origins and welds them across cube seams (ball-pivoting
                  bridge + hole fill + cleanup) into one .obj.


BUILDING (in brief)
-------------------
Windows (Visual Studio 2022): build scrollfiesta.sln as Release|x64, e.g.
    msbuild scrollfiesta.sln -p:Configuration=Release -p:Platform=x64 -m
  Binaries land in build\Release\ . The vendored dependencies (libtiff, zlib,
  Clipper2) build from deps\src\ via CMake into deps\lib\win64\Release\ ;
  Triangle compiles inline. Build those library .libs once before the solution
  -- run build-deps.ps1 (or read it for the exact CMake invocations).

Linux (GCC 13): stage the vendored static libs once, then build the binaries:
    ./build-deps.sh                   # Triangle + Clipper2 -> deps/lib/*.a
    cd src && make release -j$(nproc)  # all six binaries -> src/
  Needs gcc/g++ 13, cmake, system libtiff (libtiff-dev), OpenMP, and pthread.
  Unlike the Windows build, libtiff/zlib are NOT vendored here -- the system
  libtiff is used via -ltiff; build-deps.sh stages only Triangle and Clipper2.
  "make release" builds cube_mesh, grid_pipeline, grid_weld and the three
  diagnostic tools (obj_components, pinhole_verdict, seam_audit); plain "make"
  builds just cube_mesh (debug + AddressSanitizer). The Makefile mirrors the
  .vcxproj source lists.



INPUT FORMAT
------------
The pipeline consumes a GRID of per-cube prediction TIFFs:

    <grid>/cubes_PRED/z#####_y#####_x#####.tif

Each TIFF is a 128 x 128 x 128 multipage uint8 volume (128 pages of 128x128),
a binary mask: 0 = background, nonzero = recto surface. The filename encodes
the cube's voxel-space origin (z, y, x, zero-padded to 5 digits); adjacent
cubes differ by 128 on one axis, and that id doubles as the world offset used
to stitch cubes together. A single cube TIFF on its own also works (see "One
cube" below).


USAGE
-----

1) Whole grid -- the usual case
- - - - - - - - - - - - - - - -
    grid_pipeline <grid_dir> <output_dir> --halo 13 --qem

  <grid_dir> must contain cubes_PRED/*.tif . This meshes every cube and welds
  them, producing:
      <output_dir>/welded.obj            the stitched scroll surface
      <output_dir>/welded.obj.weld_report.json   topology audit
  plus the per-cube mesh dumps.

  Common options:
      --halo N              neighbour voxels loaded for seamless cube borders
                            (default 13; do not go below the MLS kernel radius)
      --qem / --no-qem      enable / skip decimation (keep full resolution)
      --threads-per-cube N  OpenMP threads per cube (default 1)
      --max-concurrent N    how many cubes to mesh at once (default ~#cores)
      --max-cubes N         stop after N cubes (quick smoke test)
      --skip-weld           mesh every cube but don't stitch
      --exe PATH/--weld PATH  override the cube_mesh / grid_weld binary paths
                              (default build/Release/...)

  Example (Windows, a 4x5x5 region under PHerc0139-4x5x5/):
      build\Release\grid_pipeline.exe PHerc0139-4x5x5 output\run1 --halo 13 --qem


2) One cube
- - - - - -
    cube_mesh <cube.tif> <out.obj> --halo 13 --dump-obj <dump_dir>

  Meshes a single cube and writes per-component OBJs under
      <dump_dir>/<cube_id>/<cube_id>_<stage>/
  The final per-cube mesh is the "raw_snap" stage.

  Notes:
    - The positional <out.obj> is a placeholder; meshes are written ONLY when
      you pass --dump-obj.
    - With --halo > 0 the cube's up-to-26 neighbours are read from the SAME
      directory as <cube.tif>, so a single-cube halo run needs the neighbour
      TIFFs present alongside it.
    - --no-qem keeps full resolution; --no-timeout disables the wall-clock guard.
    - Set VESUVIUS_THREADS=N to control OpenMP threads.


3) Weld already-meshed cubes by hand
- - - - - - - - - - - - - - - - - -
    grid_weld <dir_of_per_cube_objs> <welded.obj>

  Reads the per-cube "raw_snap" OBJs from the directory, welds across seams,
  and writes <welded.obj>. Useful flags:
      --no-holefill / --no-cleanup / --no-bridge   turn off individual stages
      --dump-stages <dir>                          write each weld stage as a
                                                   colored OBJ for debugging


OUTPUT
------
A welded, single-sided, orientable .obj triangle mesh of the scroll surface
(vertex order Z Y X, world-offset applied per cube). Hand it to an off-the-
shelf surface parameterizer (e.g. ABF++) to flatten it for ink reading -- that
flattening step is outside this toolkit.


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
- CPU-only: no GPU is used. The U-Net that produces the input predictions is a
  separate, upstream tool and is not part of ScrollFiesta.
- Topology is the priority: the weld is tuned to never merge distinct scroll
  wraps and never split a single sheet, even at the cost of leaving a seam gap.
- Each executable prints its usage when run with no arguments.
- Architecture and design rationale: submission.pdf
- Example outputs: sample_outputs/

