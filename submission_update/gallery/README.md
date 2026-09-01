# August 2026 progress exhibit

This directory is the evidence index for Appendix B, Progress Exhibit Gallery.
The gallery deliberately separates visual quality, geometric scale, and fit
state.

## Start here

- PHerc0139 10x10x10 coherent hybrid strip:
  ../figures/quadribbon_10x10x10_hybrid_strip.jpg
- PHerc0139 21x21x21 certified winding page:
  ../figures/sheetbook_wind50.png
- PHerc1203 complete 5x5x5 geometry views:
  ../figures/prize_pherc1203_5x5x5_views.png
- PHerc1447 complete 5x5x5 geometry views:
  ../figures/prize_pherc1447_5x5x5_views.png
- W043 geometry audit before/after:
  ../figures/w043_original_geometry_audit.png and
  ../figures/w043_native_sheet_fit_audit.png
- Machine-readable claims and boundaries: results.json
- File-integrity ledger: SHA256SUMS.txt

## Result ladder

| Target | Region | Verified result | Boundary |
|---|---:|---|---|
| PHerc0139 | 4x5x5 | Current integrated fitted/baked page | Visual quality-mask faults fall 16%; zero primary-sheet hard conflicts |
| PHerc0139 | 10x10x10 | 16.444M vertices / 32.608M faces | Hybrid champion; not the current general pre-weld fitter |
| PHerc0139 | 21x21x21 | 9,110 meshes; 195.971M vertices / 381.764M faces | Stage-1 concat only; global fit was not run |
| PHerc1203 | 5x5x5 | 125/125 cubes, zero reject/fail; 125/125 manifold | Prediction-derived geometry only |
| PHerc1447 | 5x5x5 | 125/125 cubes, zero reject/fail; 125/125 manifold | Prediction-derived geometry only |

The 21x21x21 stage-1 VMesh is 6,932,816,788 bytes. It is not duplicated into
the paper repository. Its SHA-256 is
398ae1acd722604aa24c0683629a37fe19bbd80a71d8975d09989e30e7102a33,
and its 9,110-row cube ledger independently totals to 195,970,607 vertices and
381,764,120 faces.

## Prize-target render contract

Each prize-target image was generated from exactly one canonical step-12 OBJ
per cube. The 125 objects were concatenated without welding or geometric
changes, then rendered in four views: isometric, top, front, and side. Colour
is view-dependent surface-normal shading; it is not a material, winding, or
confidence label.

Both carved grids contain prediction TIFFs only. No RAW texture, bake, ink
prediction, global weld, or quadribbon fit is claimed for these prize-target
regions.

## W043 downstream check

The native 9.362-micrometer normal-only fit improves median signed geometry
offset from +0.810 to +0.395 voxels, RMSE from 2.027 to 1.806, q95 absolute
offset from 4.082 to 3.590, and supervised AUC from 0.6901 to 0.7098. This is a
real but modest downstream improvement, not oracle geometry.

## Provenance

The new geometry artifacts were generated at experiment commit
d92a623cd0e5b3b3c884eec989c6a297e92885c5. The public report snapshot used to
prepare this exhibit was
cfbc90ea3bba92162f0b532a34fffae1382639ba. Absolute workstation paths and the
multi-gigabyte mesh payload are intentionally excluded.