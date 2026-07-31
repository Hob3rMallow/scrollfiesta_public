# PHerc0139 4x5x5

This is the exact checkpoint-58 atlas export on the pinned 100-cube test bed.

- 171,769 faces, 600 charts, and 16 groups.
- Cross-group overlap: 282,369 to 0.
- Overlapping groups: 12 to 0.
- Atlas exports `after_uv.obj`, `atlas_bake.obj`, and `groups_xyz.obj` match
  checkpoint 58 byte for byte.
- Fresh RAW bake: 41,222 x 511 pixels, 32.0% filled, intensity window
  `[43,191]`, and 96.90% diagnostic-ok pixels.

`weld/welded.obj` is the best retained flat seam weld for the same 100-cube
fixture (716,149 vertices / 475,503 faces). It is an inspection artifact; the
checkpoint-58 atlas consumes the pinned per-cube placement rather than this
flat weld.
