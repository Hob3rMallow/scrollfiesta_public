# Documentation map

ScrollFiesta contains public documentation, implementation notebooks, a paper,
and vendored research code. They do not all have the same stability.

## Current interfaces and workflows

| Document | Use it for |
|---|---|
| [`README.md`](../README.md) | Project overview, build instructions, quick start, and repository map. |
| [`PIPELINE.md`](PIPELINE.md) | Supported whole-scroll stages and required run artifacts. |
| [`SHEET_ASSEMBLY.md`](SHEET_ASSEMBLY.md) | The native sheet assembler (September 2026): stages, configuration, artifacts and finishing tools. |
| [`configs/default.json`](../configs/default.json) | The shared public configuration, calibrated from the PHerc0139 4×5×5 profiles. |
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
- [September release evidence](../submission_update/release/README.md) records
  the public configuration's CT preview, geometry audit, and bake provenance.
- `previous_submissions/` contains superseded submissions.
- `docs/research/` contains vendored research code used by the implementation;
  its upstream documentation and licensing terms remain authoritative for that
  code.

To rebuild the paper with a LaTeX distribution that provides pdfLaTeX and
BibTeX, run these commands from `submission_update/`:

```text
pdflatex -interaction=nonstopmode -halt-on-error -no-shell-escape submission.tex
bibtex submission
pdflatex -interaction=nonstopmode -halt-on-error -no-shell-escape submission.tex
pdflatex -interaction=nonstopmode -halt-on-error -no-shell-escape submission.tex
```

Inspect the rendered pages, then copy `submission_update/submission.pdf` to
the canonical `submission.pdf` at the repository root. Auxiliary LaTeX files
and the build copy of the PDF are ignored. Refresh the release evidence's
`SHA256SUMS.txt` when its PDF, result ledger, or figure changes.
