# Third-party licenses

ScrollFiesta's own code is MIT-licensed (see `LICENSE`). The vendored
dependencies under `deps/src/` are **not** covered by that grant; each carries
its own license, summarized below. The full license texts live inside each
dependency's directory.

## ⚠ Triangle — restricted commercial distribution

**`deps/src/triangle/`** — Jonathan Shewchuk's Triangle (constrained Delaunay
triangulation; used by the hole-filling stage). License terms (from
`deps/src/triangle/README`):

> These programs may be freely redistributed under the condition that the
> copyright notices [...] are not removed, **and no compensation is
> received**. Private, research, and institutional use is free. [...]
> Distribution of this code as part of a commercial system is permissible
> **ONLY BY DIRECT ARRANGEMENT WITH THE AUTHOR.**

This means ScrollFiesta binaries containing Triangle may be redistributed
freely for non-commercial purposes with notices intact, but may **not** be
sold or distributed as part of a commercial system without an arrangement
with the author. If that constraint matters for your use, build without the
CDT hole-fill path or contact the author. (A future build option to swap in a
permissively-licensed CDT is on the roadmap.)

## Compiled into the library / binaries

| Dependency | Path | License | Used for |
|---|---|---|---|
| Triangle | `deps/src/triangle/` | Shewchuk's custom terms (see above) | CDT hole filling |
| Clipper2 | `deps/src/Clipper2/` | Boost Software License 1.0 (`deps/src/Clipper2/LICENSE`) | 2D polygon booleans in hole filling |
| andres/graph | `deps/src/graph/` | No license text in the vendored copy; header-only, from <http://www.andres.sc/graph.html> / <https://github.com/bjoern-andres/graph> — verify upstream terms before commercial redistribution | Lifted-multicut solver (overlap separation) |
| libtiff 4.7.1 | `deps/src/tiff-4.7.1/` | libtiff license (BSD-like; `deps/src/tiff-4.7.1/LICENSE.md`) | TIFF I/O — **Windows builds and CLI tools only**; the library core can be built without TIFF (`SCROLLFIESTA_WITH_TIFF=OFF`), and Linux tool builds link the system libtiff |
| zlib | `deps/src/zlib/` | zlib license (`deps/src/zlib/LICENSE`) | libtiff dependency (Windows builds only) |

## Present in the tree but NOT compiled or linked

| Dependency | Path | License | Status |
|---|---|---|---|
| PoissonRecon | `deps/src/PoissonRecon/` | MIT/BSD-style (`deps/src/PoissonRecon/LICENSE`) | Vestigial — belonged to the retired Marching-Cubes/Poisson pipeline; not built by any current build system |
| probabilistic-quadrics | `deps/src/probabilistic-quadrics/` | MIT (`deps/src/probabilistic-quadrics/LICENSE`) | Reference only — the probabilistic-quadric term in `src/common/qem.c` is an independent C reimplementation of the paper; the header is not included |
