# Release notes

## October 2026

### Close-ups and fibre tools

Four new tools come from the development tree's close-up chain:
`ink_closeup`, `fiber_flatten`, `closeup_strip_turn` and `vmesh_fiber_rigid`.
See [CLOSEUPS.md](CLOSEUPS.md).

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
- **The development tree's engine.** `ink_closeup` is the development tree's
  engine of record. On its hard qualification replay, all eight output files
  are byte-identical to the record's.
- `vmesh_tifxyz --canvas --native-grid` exports a close-up as one tifxyz at one
  sample per voxel. `meta.json` now records `uv_origin`, which the re-UV reads.

### Shared modules

The close-up engine needs the development tree's current SLIM, sparse solver,
contact walk, RAW, TIFF and Zarr readers and exception handling. Those modules
now match the development tree, with this repository's own fixes kept.

- **The mesher and weld are unchanged.** On the 4×5×5 fixture, all 1,181
  meshes and the weld log are byte-identical to the previous build's, both with
  the default trim and with `--trim-inset 0`.
- **`sheet_assemble` output changes slightly.** The development tree's SLIM
  and sparse-solver changes move chart flattenings at floating-point level.
  They include a stable matrix assembly order: duplicate entries are now summed
  in input order, not in the order an unstable sort left them. On the 4×5×5
  reference the cleaned charts of 13 of the 100 cubes differ. Through audit,
  coverage goes from 81.13% to 81.27%, overlapping triangle pairs from 213 to
  209, strict bad faces from 73,483 to 72,255, and the singular-value range
  from [0.023, 13.2] to [0.116, 11.4]. Unresolved seams go from 1,314 to 1,324.
  Repeated runs are byte-identical.
- The exception stack is thread-local, so a `RAISE` on an OpenMP worker can no
  longer jump into another thread's frame.

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
- `ctest` runs the new tools' selftests, plus `ink_closeup
  --material-contact-selftest`. All 22 tests pass in a Release build and under
  MSVC AddressSanitizer.

### Known issues

- On the reference 4×5×5 fixture, `grid_weld`'s exact chart disk certificate
  fails (762 of 781 charts are disks), so the weld writes no `welded.obj` and
  its log ends `grid_weld: I/O failure`. `grid_pipeline` reports this as
  non-fatal. The per-cube meshes that `scroll_whole` and `sheet_assemble`
  read are unaffected.
- Not merged yet: #25, whose gradient change alters assembler output on every
  platform (this release removes the FMA contraction it reports instead), and
  #11 and #13, which change unrolling defaults.
