# Documentation map

ScrollFiesta contains public documentation, implementation notebooks, a paper,
and vendored research code. They do not all have the same stability.

## Current interfaces and workflows

| Document | Use it for |
|---|---|
| [`README.md`](../README.md) | Project overview, build instructions, quick start, and repository map. |
| [`PIPELINE.md`](PIPELINE.md) | Supported whole-scroll stages and required run artifacts. |
| [`CONTRIBUTING.md`](../CONTRIBUTING.md) | Development profiles, validation expectations, and merge checklist. |
| [`include/scrollfiesta.h`](../include/scrollfiesta.h) | Authoritative public C API, ownership, errors, ABI, and coordinate conventions. |
| [`output/canonical_best/README.md`](../output/canonical_best/README.md) | Curated result bundle and artifact contract. |
| [`THIRD_PARTY_LICENSES.md`](../THIRD_PARTY_LICENSES.md) | Dependency inventory and redistribution constraints. |

The generated API behavior, executable `--help` output, tests, and source code
take precedence when an implementation notebook disagrees with a current
interface document.

## Implementation notebooks

The Markdown files under `src/common/`, `src/extract/`, and `src/split/` are
useful design history. They predate parts of the current public API and contain
some references to retired private planning guides. Treat them as explanatory
notes, not as a supported contract, and verify signatures and behavior against
the adjacent headers and source.

## Paper and research dependencies

- [`submission.pdf`](../submission.pdf) is the current rendered paper.
- `submission_update/` contains its LaTeX source and figures.
- `previous_submissions/` contains superseded submissions.
- `docs/research/` contains vendored research code used by the implementation;
  its upstream documentation and licensing terms remain authoritative for that
  code.
