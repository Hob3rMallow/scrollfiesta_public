#pragma once

#include "mls_cubecl.h"

/* No layout adapter is required. ScrollFiesta's Arena_T is allocator scratch
 * state, and MLS_project_verts intentionally ignores it. Include this header
 * at the existing call site and link the CubeCL library; the exported symbol
 * consumes ScrollFiesta's native z-y-x buffers and additive cell origin.
 *
 * The caller must keep the existing ping-pong loop. Each call performs one
 * tangent-plane projection and rebuilds support from that pass's input. */
