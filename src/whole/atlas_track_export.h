#ifndef ATLAS_TRACK_EXPORT_INCLUDED
#define ATLAS_TRACK_EXPORT_INCLUDED

#include <stddef.h>

#include "../common/arena.h"
#include "../unroll/scaffold.h"
#include "atlas_ribbon_fit.h"
#include "atlas_track_grow.h"

typedef struct {
    size_t pieces;
    size_t written;
    size_t empty;
    size_t over_cap;
    size_t quarantine_faces;
    size_t valid_pixels;
    size_t multi_pixels;
    size_t conflict_pixels;
    double conflict_fraction;
    double seconds;
} AtlasTrackExportStats;

/* Adapt the evidence-grown mesh to the maintained L3 winding organizer.
 * The organizer writes one VC3D tifxyz segment per (winding x axial slab),
 * plus atlas.json.  Local grown UV remains the raster geometry; snapped
 * sample winding supplies the physical piece identity. */
int AtlasTrackExportWinding_run(
    Arena_T arena,
    const AtlasRibbonObservationSet *evidence,
    const ScaffoldCalib *calibration,
    const AtlasTrackGrowResult *result,
    const char *output_root,
    const char *prefix,
    int wraps_per_piece,
    double slab_v,
    double du,
    double dv,
    AtlasTrackExportStats *out);

#endif
