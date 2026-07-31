#ifndef COMPONENT_CULL_INCLUDED
#define COMPONENT_CULL_INCLUDED

#include <stddef.h>
#include "../common/arena.h"
#include "../common/mesh_types.h"

/*
 * "Kibble" removal: a connectivity pass + surface-area filter, run after
 * hole-fill/QEM. Splits every input mesh into its mesh-connectivity-components,
 * computes each component's triangle-area sum, and drops any component whose area
 * is < min_frac of its PARENT input mesh's area. The input array already carries
 * the semantic sheet split; using the whole-cube total catastrophically removes
 * valid wraps when a cube contains more than 1/min_frac sheets. At least the
 * largest connectivity component of every non-empty input is retained.
 *
 * Output is the surviving components, one ComponentMesh per kept connectivity-
 * component (arena-allocated, comp_id reassigned 1..*n_out, pca_normal/centroid
 * recomputed). *n_culled (optional) receives the number of components dropped.
 * Empty / trivial input returns cleanly. Returns 0 on success.
 */
int ComponentCull_by_area(Arena_T arena,
                          const ComponentMesh *in, size_t n_in,
                          float min_frac,
                          ComponentMesh **out, size_t *n_out,
                          size_t *n_culled);

#endif
