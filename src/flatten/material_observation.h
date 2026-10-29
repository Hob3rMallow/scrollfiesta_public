#ifndef MATERIAL_OBSERVATION_INCLUDED
#define MATERIAL_OBSERVATION_INCLUDED

#include "scroll_model.h"

/* Independent crop-local SOURCE/EVIDENCE stage of the revised lane. Streams
 * original cubes, retains all source support in a ledger, and solves only the
 * compact soft phase graph as an initialization diagnostic. Does not yet
 * publish an intrinsic UV map or certify continuation. No full-scroll model,
 * cached chart labels, inherited coordinates or source edits are involved.
 * Bounds are inclusive source cube origins; no hidden halo is loaded. */
int MaterialObservation_build(Arena_T arena,const char *source_dir,
                               const char *output_dir,const long *bounds,
                               const ScrollModelOptions *options);

#endif
