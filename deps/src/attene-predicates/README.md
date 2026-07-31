# attene-predicates (vendored)

Exact geometric predicates for the CVT / restricted-Voronoi remesher
(`src/remesh/rvd_clip.c`, via the C wrapper `src/common/attene_predicates_wrap.*`).

These files are a clean-room reimplementation of Marco Attene's indirect predicates
("Indirect Predicates for Geometric Constructions", arXiv:2105.09772) plus a
GTE-derived exact-arithmetic backend (BSNumber / UIntegerAP). They were authored for
a private engine and lifted here with the author's permission; **within this
repository they are covered by this repository's license.**

Layout:
- `attene_predicates.{h,cpp}` — interval filter + expansion (error-free) arithmetic
- `attene_implicit_point3d.{h,cpp}`, `attene_indirect_predicates3d.cpp` — implicit
  point constructions (line-plane, three-plane) + indirect orient3d over them
- `bsnumber.h` — exact binary-scientific number (header-only)
- `uinteger.{h,cpp}` — arbitrary-precision unsigned-integer backend
- `ieee_binary.{h,cpp}`, `bit_scan.h` — IEEE-754 decomposition + bit helpers
- `buf.h` (`BUF_STANDALONE`) — stretchy buffer used by `uinteger.cpp`
- `attene_vec.h` — minimal `Vec2d/Vec3d/Vec3f` shim (replaces the engine's vector lib)

Self-contained (no GMP/MPFR), C++17, MSVC + GCC/Clang. Only the RVD-relevant subset
is vendored (the 2D implicit-point files are not needed). The Voronoi "side"
predicates the remesher actually calls are built on this in
`attene_predicates_wrap.cpp`, not in these files.

Build note: on GCC/Clang compile the wrapper TU with `-frounding-math -fno-fast-math`
— the interval filter depends on directed (round-toward-+inf) FP arithmetic, which
MSVC scopes via `#pragma float_control` / `fenv_access` here.
