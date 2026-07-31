# PHerc0139 4x21x21

Fresh end-to-end pre-ribbon run using the CubeCL CUDA mesh fleet and 49
bounded atlas tiles (maximum loaded tile: 4x5x5 / 100 cubes). The reflected
frame legitimately produces a negative signed ladder
`-0.9882629717598193`.

- 1,753 input cubes: 1,702 accepted and 51 rejected as empty/solid slabs.
- 3,958,802 atlas faces.
- Cross-group overlap: 1,651,025 to 926,023 (43.9% reduction).
- Intra-group overlap: 64,491,422 to 24,141,976 (62.6% reduction).
- Fresh RAW bake: 353,441 x 172 pixels, 11.1% filled, intensity window
  `[40,196]`, and 98.82% diagnostic-ok pixels.

The placement audit reports high join completeness (`|du|<2`: 82.23%) but an
8.10% whole-turn error, above the strict 5% gate. The final atlas records 1,952
torn weld relations, with zero torn lateral or family relations. These
limitations remain visible in the preserved stats rather than being hidden.

No flat whole-grid weld is included: that historical step was the
memory-infeasible bottleneck. `atlas/atlas_bake.obj` and
`atlas/groups_xyz.obj` are the canonical giant meshes.
