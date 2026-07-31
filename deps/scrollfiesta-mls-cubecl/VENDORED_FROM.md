# Vendored provenance

Upstream: https://github.com/altommo/scrollfiesta-mls-cubecl

Pinned commit: `5fa316fdde43fa48158556bfec13e11e14406776`

Vendored on 2026-07-24. The upstream package declares `MIT OR Apache-2.0` in
`Cargo.toml`. Before redistribution, preserve that declaration and obtain the
corresponding license text from upstream if it is not added there.

Local integration changes are deliberately small and opt-in:

- add CubeCL 0.10's CUDA runtime as feature `cubecl-cuda`;
- add backend id 4 / `MLS_BACKEND=cubecl-cuda` to the existing C ABI;
- emit both `cdylib`/import-library and `staticlib` artifacts;
- add a CUDA current-support numerical smoke test.

The tracked `PACKAGE_MANIFEST.sha256` is the upstream package manifest and does
not cover these local changes.
