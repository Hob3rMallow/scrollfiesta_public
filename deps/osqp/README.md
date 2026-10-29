OSQP 1.0.0, QDLDL 0.1.8 and its bundled AMD implementation.

These are the C sources used by the installed reference Python package
`osqp==1.1.3`. All matching core files were checked against the
[upstream v1.0.0 release](https://github.com/osqp/osqp/releases/tag/v1.0.0),
allowing only CRLF/LF conversion. `sources.json` records every copied file
and its hash. QDLDL and generated type headers come from that wheel's codegen
sources. Polishing and AMD come from the upstream release.

The project compiles these sources directly, with double values, 32-bit
indices, the built-in direct solver, and dynamic setup. The executable has
no Python dependency. Licenses are retained in `LICENSE` and `amd/LICENSE`.
`osqp_configure.h` is the local build configuration; upstream sources are
unchanged.
