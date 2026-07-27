#ifndef RAW_SNAPBACK_INCLUDED
#define RAW_SNAPBACK_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "arena.h"
#include "mesh_types.h"

/* Raw prediction snap-back.
 *
 * After all topology-cleaning mesh operations (Poisson, hole fill, bay fill)
 * but BEFORE revoxelization, pull mesh vertices toward the raw nnUNet
 * prediction surface.  Topology is encoded in connectivity (unchanged by
 * vertex movement), so this is safe.
 *
 * Per-vertex alpha: strong pull for vertices near the raw surface,
 * pure Laplacian smoothing for vertices far away (hole fill regions). */
int RawSnap_process(Arena_T arena,
                    ComponentMesh *meshes, size_t n_meshes,
                    const uint8_t *input_vol,
                    size_t D, size_t H, size_t W,
                    int iterations, float alpha, float dist_max);

#endif
