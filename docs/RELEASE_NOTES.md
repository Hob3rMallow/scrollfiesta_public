# Release notes

## October 2026 update: close-ups, fibre tools, a better mesher

### Close-ups and fibre tools

Seven new tools come from the development tree's close-up chain:
`ink_closeup`, `fiber_flatten`, `closeup_strip_turn`, `vmesh_fiber_rigid`,
`closeup_repair`, `tifxyz_snap` and `vcr_gpu`, along with `sheet_strip
--closeup-window`. See [CLOSEUPS.md](CLOSEUPS.md).

- **Close-ups.** Any region of an assembled sheet becomes one continuous,
  isometric surface. Every nearby chart is an observation, and a robust fused
  surface fills holes and borders. SLIM flattens it with the scroll axis down
  the rows. The render writes the 21-layer surface volume at one pixel per
  voxel, with a CT bake and a provenance map.
- **Growth.** A close-up grows past its box by adopting charts from pools,
  including charts no assembly placed, and bridges gaps along the F0
  prediction ridge. On the PHerc0139 4×5×5 reference, a 400 × 350 px seed grew
  from 15,092 to 59,391 data nodes in 38 s.
- **Fibre alignment.** Fibres are measured physically, from the CT's fibre
  cross or from an inferred 3-D fibre field. The close-up is then re-flattened
  so the fibres run along the page axes: angular local/global least squares
  from a checked SLIM start, and XYZ never moves. `fiber_flatten` re-flattens a
  canvas along a fibre field. `closeup_strip_turn` and `vmesh_fiber_rigid`
  apply exact rigid turns.
- **Repair and snap.** `closeup_repair` removes duplicate fill, wrong-side
  fill and flaps from a close-up surface. `tifxyz_snap` snaps any tifxyz
  segment onto the papyrus recto edge read in the RAW CT. On the 4×5×5 demo
  close-up, 82.3% of nodes find the edge, and the median offset from it falls
  from +0.44 to +0.01 voxels.
- **Surface volumes on the GPU.** `vcr_gpu` renders villa's
  `vc_render_tifxyz` surface volume on an NVIDIA GPU and refuses any option it
  does not reproduce. On the demo close-up, its 21 layers are voxel-identical
  to `vc_render_tifxyz`'s, in 0.3 s against 6 s. It loads the CUDA driver and
  NVRTC at run time, so building it needs no CUDA toolkit.
- **Completion in one process.** `ink_closeup <job.json> --complete-chain`
  runs a close-up's completion passes in one process. It reads the pools
  once, keeps each pass's surface in memory and writes only the final pass's
  mesh. Exact speedups come with it: the surface fit's data sums, incremental
  refits, parallel hole solves, KD-tree builds and normals. On the
  development tree's six-pass completion of a PHerc0306B close-up, the stage
  went from 1,794 s at 16 threads to 147 s at 32, and every output stayed
  byte-identical.
- **Back into the strip.** `sheet_strip --closeup-window` places a close-up
  back into its reading strip at native scale, clear of stationary fragments.
  `sheet_strip` also reads its blocks from a list file (`@blocks.txt`), turns
  pieces, and can bound each line's pixels.
- **The development tree's engines.** `ink_closeup` is its engine of record
  of 30 October. It reproduces all eight outputs of the record's hard
  qualification replay byte for byte, and all 36 files of its six-pass
  completion run. `vcr_gpu`'s output is byte-identical to the development
  tree's. All 130 recorded close-up window jobs plan identically to the record
  engine's.
- `vmesh_tifxyz --canvas --native-grid` exports a close-up as one tifxyz at one
  sample per voxel. `meta.json` now records `uv_origin`, which the re-UV reads.

### Mesher: radius 8, clip-at-box and B2G

The per-cube mesher now matches the development tree's mesher of record. Its
meshes change.

- The step-0 MLS kernel radius goes from 12 to 8 voxels. At 12 it reached the
  neighbouring sheet where bands are thin and close together, and projected
  two sheets onto one surface in the gap.
- With a halo, the final trim cuts every chart at the exact cube box, so
  adjacent cubes' charts end on the same plane. `--trim-inset` now only sets
  how far the step-0 cloud reaches into the halo.
- The B2G CVT defaults: three topology retries, early probe aborts and the
  manifold guard. Earlier releases left them out to keep their meshes
  byte-identical.
- `cube_mesh` writes only final meshes unless given `--full-dumps`. An opt-in
  `--surface f0-inner` places the surface on the prediction's inner boundary.
- On the 4×5×5 fixture with production flags (`--halo 13 --trim-inset 0`), all
  1,201 meshes are byte-identical to the record executable's.
- Against the previous mesher on the same fixture: mesh area within 1 voxel
  of a predicted band centre goes from 94.1% to 97.4%, area 2+ voxels off from
  1.08% to 0.13%, and stacked duplicates at cube faces from 0.09% to 0.01%.
- `sheet_assemble` on the new meshes, through audit: coverage goes from 81.3%
  to 88.9%, placed charts from 590 to 645, overlapping triangle pairs from 209
  to 59, strict bad faces from 72,255 to 6,226, and unresolved seams from 1,324
  to 814. One chart still fails the metric audit, and its worst face is now far
  more stretched: the largest singular value goes from 11.4 to 558.

### Shared modules

The close-up engine needs the development tree's current SLIM, sparse solver,
contact walk, RAW, TIFF and Zarr readers and exception handling. Those modules
now match the development tree, with this repository's own fixes kept.

- **On their own, they leave the mesher and weld unchanged.** On the 4×5×5
  fixture, all 1,181 meshes and the weld log were byte-identical to the
  previous build's, both with the default trim and with `--trim-inset 0`.
- **`sheet_assemble` output changes slightly.** The development tree's SLIM
  and sparse-solver changes move chart flattenings at floating-point level.
  They include a stable matrix assembly order: duplicate entries are now summed
  in input order, not in the order an unstable sort left them. On the 4×5×5
  reference, with the previous meshes, the cleaned charts of 13 of the 100
  cubes differ. Through audit, coverage goes from 81.13% to 81.27%,
  overlapping triangle pairs from 213 to 209, strict bad faces from 73,483 to
  72,255, and the singular-value range from [0.023, 13.2] to [0.116, 11.4].
  Unresolved seams go from 1,314 to 1,324. Repeated runs are byte-identical.
- The exception stack is thread-local, so a `RAISE` on an OpenMP worker can no
  longer jump into another thread's frame.
- **Contact tests finish instead of giving up.** A material-boundary test
  whose gather overflowed its 512 segments, or whose band test ran out of
  pieces or witness candidates, used to leave its pair unmeasured. A
  certificate walk then reported its overlaps as a lower bound, and an
  objective walk failed. Now the test is retried with up to 64 times the
  storage and 256 times the work budget, and measured. A pair that still
  cannot finish says which limit stopped it.
- **Off the Windows heap lock.** In MSVC builds, arena chunks, the sparse
  solver's temporaries and AMD's workspace come from each thread's own heap,
  as TAUCS's blocks already did. In the development tree, 57% of a 32-thread
  assembly's thread samples were waiting on the process heap's lock in
  exactly these places. Only allocation changes.
- **`AsmRelation` has no implicit padding.** GCC 16 at `-O3` copied relations
  member by member and left stack bytes in the two padding holes, which
  `AsmStore` and the resume files then wrote. The holes are named, zeroed
  members now, and a compile-time check keeps new padding out. The layout and
  file formats are unchanged.
- The vendored CLAPACK, the CMake build's last-resort LAPACK, compiles with
  GCC 14 and later.
- These last four changes and the completion stage leave the 4×5×5 outputs
  as they were with the new mesher: the 1,201 meshes, the weld log and a full
  `sheet_assemble` run through audit are byte-identical, reports apart from
  timings.

### Tests

- `ctest` runs every new tool's selftest, plus `ink_closeup
  --material-contact-selftest`; 25 tests in all. They pass in a Release build
  and under MSVC AddressSanitizer. `vcr_gpu`'s GPU selftest needs an NVIDIA
  driver, so `ctest` runs its storage selftest unless the build sets
  `-DSCROLLFIESTA_GPU_TESTS=ON`. Its GPU selftest passes on an RTX PRO 6000.

### Known issues

- `grid_weld`'s exact chart disk certificate still fails on the 4×5×5
  fixture with the new mesher (787 of 807 charts are disks; 762 of 781
  before), so the weld writes no `welded.obj`. `sheet_assemble` and
  `scroll_whole` read the per-cube meshes and are unaffected.
- The default `sheet_strip` page mode needs a review layout with
  `ribbon_layout.json`, which this repository writes only from a qualified
  layout (`sheet_assemble --winding-order`). The updated `sheet_strip` was
  therefore checked by its selftest and the 130 close-up window replays, not
  on a real strip page.

## October 2026

### Contributed fixes

Thanks to everyone who built, measured and fixed ScrollFiesta this month.
Merged in this release:

- **#12** (TAUIL-Abd-Elilah): one 2-D Clipper2 ABI in every build. The CMake
  build linked both `Clipper2` and `Clipper2Z` and leaked `USINGZ` into the 2-D
  hole-fill wrapper, so a Linux build meshed 99 of the 100 reference cubes and
  lost 16 of 64 PHerc. 1203 cubes to a crash in hole filling (#26, jo-tud).
- **#14** (tesinclair) and **#21** (evilaliv3): a clean clone builds on Linux
  with GCC 14 and 15, without root; `build-deps.sh` stages the vendored zlib
  and libtiff when the host has no libtiff package.
- **#23** (evilaliv3): the vendored CLAPACK declares `integer` as `int`, which
  is what TAUCS passes on LP64 systems.
- **#24** (rk23): native AppleClang builds with libomp, and Mach memory
  sampling on macOS.
- **#16** (kadenpool): `reorder_welded_to_xyz` keeps texture coordinates, so
  atlas meshes go into villa's `vc_obj2tifxyz`.
- **#22** (anveshrd): an exception from `Clipper2_union` refuses one hole
  instead of ending the process, with a regression test for the loop that
  threw it.
- **#17** (evilaliv3): the bounded nearest-neighbour search in the resplit
  assignment. The exact assignment below replaced that loop in this same
  release.

### Mesher

- **Faster, same meshes.** The development tree's exact speedups are in: the
  KD tree, MLS kernel, resplit assignment, CVT/RVD (an exact double tier for
  the bisector predicates, a whole-triangle filter, a tight candidate radius),
  ball pivoting, the vertex weld, Hoppe orientation and both multicut solvers.
  On the 100-cube PHerc0139 4×5×5 fixture, all 1,181 meshes and the weld log
  are byte-identical to the September release's, in two runs of each build.
  Per-cube work fell 1.9× and the grid's wall time 2.3× (171–178 s → 76 s,
  32 cubes at a time, one thread each, Windows/MSVC).
- `grid_pipeline` starts the longest cubes first.
- **Triangle.** A polygon that passes one vertex twice no longer makes
  Triangle follow an uninitialised pointer, which could end in an internal
  error, an endless loop or heap corruption depending on what the heap held.
  The hole filler has refused such polygons since September; now every caller
  is safe.
- **Watchdog (Windows).** Triangle's hole-fill watchdog bounds the worker's
  CPU time instead of wall time and stops it only inside Triangle. Under load
  the old 10 s wall bound abandoned holes, so meshes depended on machine load,
  and it could hang the process at exit.
- `cube_mesh --selftest` checks every fast path against its reference, plus
  the duplicate-vertex polygon and the watchdog; `ctest` runs it.

### Builds and CI

- CI's ASan job had failed since the September release. The assembler's
  local-failure cache keyed each attempt partly on the addresses of two
  per-call stack objects, which ASan moves on every call, so identical
  attempts never hit the cache. The key ignores those addresses now, and the
  cache's selftest repeats its attempt from a deeper stack to check it.
- The same job found two memory errors, both fixed: the conflict audit read
  its report cells after the arena restore that released them, and Triangle's
  exact arithmetic loaded one element past the end of an expansion (results
  unchanged).
- GCC and Clang no longer fuse `a*b+c` (`-ffp-contract=off`), so every build
  computes what MSVC computes; this is the cause of the macOS contact-gradient
  failure reported in #25. The exact predicates get `-frounding-math
  -fno-fast-math`, as their README requires.

### Known issues

- On the reference 4×5×5 fixture, `grid_weld`'s exact chart disk certificate
  fails (762 of 781 charts are disks), so the weld writes no `welded.obj` and
  its log ends `grid_weld: I/O failure`. `grid_pipeline` reports this as
  non-fatal. The per-cube meshes that `scroll_whole` and `sheet_assemble`
  read are unaffected.
- Not merged yet: #25, whose gradient change alters assembler output on every
  platform (this release removes the FMA contraction it reports instead), and
  #11 and #13, which change unrolling defaults.
