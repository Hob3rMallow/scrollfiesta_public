# PHerc0139 10x10x10

Fresh end-to-end pre-ribbon run using `cubecl-cuda` meshing and 64 bounded
atlas tiles (maximum loaded tile: 5x5x5 / 125 cubes).

- 1,762,542 atlas faces.
- Cross-group overlap: 3,664,494 to 2,371,771 (35.3% reduction).
- Intra-group overlap: 6,249,964 to 4,039,882 (35.4% reduction).
- Fresh RAW bake: 88,560 x 428 pixels, 12.2% filled, intensity window
  `[32,200]`, and 98.23% diagnostic-ok pixels.
- Meshing completed 999 of 1,000 input cubes; 998 placed sidecars entered the
  atlas manifest.

The placement audit has strong join completeness (`|du|<2`: 79.35%) but its
whole-turn error is 5.52%, just above the strict 5% gate. This is retained as
the best terminating 10x10x10 atlas, not presented as a perfect registration.

No flat whole-grid weld is included: it is unnecessary for the atlas and is
not memory-safe at this size. `atlas/atlas_bake.obj` is the giant textured
world mesh.
