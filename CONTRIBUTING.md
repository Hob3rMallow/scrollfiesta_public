# Contributing to ScrollFiesta

ScrollFiesta is active research software with a versioned C API and a large
artifact pipeline. The main integration goal is to keep algorithm changes,
review evidence, and generated outputs separable so working branches can be
merged without obscuring regressions.

## Keep changes reviewable

- Prefer one behavioral change per branch or pull request.
- Separate mechanical cleanup from algorithm or parameter changes when
  practical.
- Do not commit local build trees, raw datasets, ordinary `output/` runs, or
  experiment logs.
- Preserve unrelated worktree changes and avoid bulk formatting of untouched
  code.
- Explain which stage is affected and what evidence establishes that the
  result improved or remained equivalent.

## Build profiles

The fastest cross-platform baseline is the library-only profile:

```bash
cmake -S . -B build \
  -DSCROLLFIESTA_BUILD_TOOLS=OFF \
  -DSCROLLFIESTA_WITH_TIFF=OFF \
  -DSCROLLFIESTA_BUILD_TESTS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

On multi-config generators, add `--config Release` to the build command and
`-C Release` to `ctest`.

Changes to CLI tools, TIFF I/O, registration, unrolling, or atlas/ribbon code
also require the full profile documented in [README.md](README.md). On Windows,
run `build-deps.ps1` before configuring or building the full toolchain.

## Validation expectations

At minimum:

1. Build the smallest profile that compiles every changed source file.
2. Run `ctest --output-on-failure` for that profile.
3. Run the relevant executable's `--selftest` when it has one.
4. For algorithm changes, run a representative fixture and retain the audit or
   metrics summary needed to compare before and after.
5. Review the final diff for generated files, accidental rewrites, debug
   output, and unrelated edits.

Changes to installation or CMake exports should also run:

```bash
cmake --install build --prefix build/install-smoke
```

CI covers the full Linux build, sanitizer tests, MSVC and MinGW library builds,
runtime loading, and a hostile embedding-consumer configuration.

## Code and interface rules

- New C code targets C17; C++ wrappers target C++17.
- Keep the public consumer surface in `include/scrollfiesta.h` and use plain C
  types at that boundary.
- API coordinates are `(x,y,z)`; internal pipeline and CLI OBJ coordinates are
  `(z,y,x)`. Axis conversion must preserve face winding.
- Library-owned memory must be released with the matching ScrollFiesta free
  function. Do not let exceptions, signals, or `longjmp` cross the API
  boundary.
- Add new sources to the top-level CMake build. Update the Visual Studio or
  low-level Makefile source lists only when the affected target is supported by
  that build.
- A public API change must update both the CMake project version and the
  `SCROLLFIESTA_VERSION_*` macros. Breaking struct, enum, or signature changes
  must also bump `SCROLLFIESTA_ABI_VERSION` and include a migration note.

## Documentation and artifacts

- Update [README.md](README.md) for entry-point or setup changes and
  [docs/PIPELINE.md](docs/PIPELINE.md) for stage or artifact-contract changes.
- Keep comments focused on non-obvious constraints and current rationale.
- Treat `src/**/*.md` as implementation notebooks; correct or clearly mark
  stale claims when code changes invalidate them.
- Do not promote a run into `output/canonical_best/` as part of an unrelated
  code change. Canonical promotions require the expected review artifacts,
  updated `SHA256SUMS.txt`, and Git LFS tracking.
- Dependency additions or upgrades must update
  [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) and preserve upstream
  notices.

## Pull-request checklist

- [ ] The change has a focused purpose and no unrelated generated files.
- [ ] Every changed source file is compiled by a documented build profile.
- [ ] Relevant CTest and embedded self-tests pass.
- [ ] Algorithm changes include before/after metrics or audit evidence.
- [ ] Public API/version/ABI implications are handled.
- [ ] User-facing commands and artifact contracts are documented.
- [ ] New dependencies and redistributed artifacts have been license-reviewed.
