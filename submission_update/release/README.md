# September 2026 public-release check

This is the evidence for the public-configuration PHerc0139 4x5x5 sheet in
the updated [paper](../../submission.pdf). The CT texture is copied directly
from the native baker, wrapped into four consecutive rows for display.

![Public-configuration sheet](../figures/public_release_4x5x5.png)

The run accounts for all 916,995 source faces and represents 646,199 of them:
**47.385% source-area coverage**, with **171 overlapping triangle pairs**.
Source geometry is preserved and there are no invalid faces or strict metric
failures. The geometry remains unqualified. These numbers are separate from
the [August research exhibits](../gallery/README.md).

## Evidence

[`results.json`](results.json) retains the full assembly audit, CT
coverage report, raster grid, configuration and artifact hashes. The saved
mesh pile is `head_39e72d0_20260826/mesh/dump`; source and field hashes
are in the audit. Large mesh, CT, and per-run outputs remain outside Git.
Configuration hashes normalize CRLF to LF; binary hashes use exact bytes.
`SHA256SUMS.txt` checks the committed figure, result ledger and PDF.

## Reproduce the workflow

Use binaries from the documented CMake Release build. `PILE` below is
the saved 100-cube mesh pile, and `RAW` is its matching directory of
128-cubed multipage CT TIFFs. Paths in the configuration are relative to the
repository root. The reference bounds are half-open: z=[4352,4864),
y=[3072,3712), x=[2560,3200).

```powershell
build\Release\sheet_assemble.exe PILE output\sheet configs\default.json --stop-after audit
build\Release\obj_bake_raw.exe output\sheet\stage6_repaired.vmesh RAW output\sheet\preview `
  --id sheet --chunk 128 --normal-range 4 --normal-samples 9 `
  --raster-du 1 --raster-dv 1 --window 47 193 `
  --skip-faces-reaching-missing-raw --no-3d --no-flat
```

The saved reference was finished by resuming `stage5_state.asr` with
`--resume-placement` under the shipped repair defaults. A fresh run may
differ with platform, parallel scheduling or budget timing; the ledger records
the exact audited result, not a bitwise reproducibility claim.

The fixed [47,193] contrast window is explicit for this preview; the assembler's
normal bake keeps its native window. The original raster is 11,170 by 642
pixels at one pixel per coarse voxel. CT sampling uses nine normal samples
over +/-4 voxels. The baker leaves 18,456 faces unpainted where those samples
would reach outside the available cube grid. Its coverage mask also identifies
38,630 pixels of bounded raster crack or pinhole fill.

The updated release passed the Windows CMake Release build, native quadribbon
build, installation check, and all 16 CTest checks. Geometry qualification is
reported separately from software validation.
