# Third-party licenses

ScrollFiesta's own code is MIT-licensed (see `LICENSE`). Vendored dependencies
retain their own terms. This file is a practical inventory, not legal advice;
use the license/notice files in each dependency directory as the authoritative
text.

## Important redistribution restrictions

### GCoptimization 3.0 — research use only

`docs/research/gco-v3.0/` is compiled into the overlap/ownership solver used by
the supported library and `scroll_unroll`. Its `GCO_README.TXT` permits use and
distribution for **research purposes only**, prohibits commercial use, restricts
redistribution, and requires the listed publications to be cited. Do not ship a
commercial binary containing this code without separate permission or replacing
the solver.

### Triangle — restricted commercial distribution

`deps/src/triangle/` is used for constrained-Delaunay hole filling. Its README
permits free redistribution only when notices remain and no compensation is
received; distribution as part of a commercial system requires a direct
arrangement with Jonathan Shewchuk. Build without that path or obtain permission
if commercial distribution matters.

## Dependencies compiled by supported builds

| Dependency | Path | Terms / notice | Used for |
|---|---|---|---|
| GCoptimization 3.0 | `docs/research/gco-v3.0/` | Research-only custom license and citation requirements in `GCO_README.TXT` | Graph-cut overlap ownership |
| Triangle | `deps/src/triangle/` | Shewchuk custom non-commercial redistribution terms in `README` | CDT hole filling |
| Clipper2 | `deps/src/Clipper2/` | Boost Software License 1.0 in `LICENSE` | 2D polygon booleans |
| andres/graph | `deps/src/graph/` | Vendored copy has no license text; verify upstream terms before redistribution | Lifted-multicut graph structures |
| attene-predicates | `deps/src/attene-predicates/` | Clean-room/project-authored subset; its README states it is covered by this repository's MIT license | Exact CVT/RVD predicates |
| libtiff 4.7.1 | `deps/src/tiff-4.7.1/` | BSD-like libtiff license in `LICENSE.md` | TIFF CLI I/O; optional in library-only builds |
| zlib | `deps/src/zlib/` | zlib license in `LICENSE` | Windows libtiff dependency |
| stb_image_write | `deps/include/stb_image_write.h` | Public-domain/MIT dual offer in the header | PNG diagnostic output |

## Dependencies used by auxiliary atlas tools

| Dependency | Path | Terms / notice | Used for |
|---|---|---|---|
| herculaneum-mls-cubecl | `deps/scrollfiesta-mls-cubecl/` | `MIT OR Apache-2.0` per its `Cargo.toml`; the pinned upstream (`VENDORED_FROM.md`) ships no separate license text file, so that declaration is the authoritative statement. Optional: only built when the GPU MLS backend is requested. | Optional CubeCL CPU/CUDA/HIP MLS projection backend |
| TAUCS | `deps/taucs/` | LGPL-3.0-or-later, GPL-3.0-or-later, or Apache-2.0 choice, plus the research-citation requirement in `LICENSE.txt` | Sparse solves in atlas tooling |
| CLAPACK 3.1.1 | `deps/clapack/` | Netlib CLAPACK/LAPACK distribution; no single top-level license file is present in this vendored release. Preserve its notices and audit the distribution before redistribution. | BLAS/LAPACK support for TAUCS |
| f2c runtime | `deps/clapack/F2CLIBS/` | AT&T/Lucent/Bellcore permission notice in `libf2c/Notice` plus notices in the bundled sources | CLAPACK runtime support |

The unused PoissonRecon and probabilistic-quadrics source drops were removed;
no current build includes them.