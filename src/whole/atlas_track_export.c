#define _USE_MATH_DEFINES
#include "atlas_track_export.h"

#include "../unroll/export_atlas.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static double ate_mod_tau(double angle)
{
    double value = fmod(angle, SCAFFOLD_2PI);
    if (value < 0.0) value += SCAFFOLD_2PI;
    return value;
}

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
    AtlasTrackExportStats *out)
{
    if (arena == NULL || evidence == NULL || calibration == NULL ||
        result == NULL || output_root == NULL || prefix == NULL ||
        out == NULL || result->nv == 0 || result->nf == 0 ||
        result->xyz == NULL || result->global_uv == NULL ||
        result->faces == NULL || result->source_target == NULL ||
        result->vertex_winding == NULL || wraps_per_piece < 1 ||
        !(slab_v > 0.0) || !(du > 0.0) || !(dv > 0.0))
        return -1;
    memset(out, 0, sizeof *out);

    PieceSet ps;
    memset(&ps, 0, sizeof ps);
    ps.nv = result->nv;
    ps.nf = result->nf;
    ps.n_cubes = 1;
    ps.verts = (float *)ARENA_ALLOC(
        arena, ps.nv * 3 * sizeof(*ps.verts));
    ps.uv = (float *)ARENA_ALLOC(
        arena, ps.nv * 2 * sizeof(*ps.uv));
    ps.phi = (float *)ARENA_ALLOC(
        arena, ps.nv * sizeof(*ps.phi));
    ps.normals = (float *)ARENA_CALLOC(
        arena, ps.nv * 3, sizeof(*ps.normals));
    ps.gid = (int32_t *)ARENA_ALLOC(
        arena, ps.nv * sizeof(*ps.gid));
    ps.faces = (int32_t *)ARENA_ALLOC(
        arena, ps.nf * 3 * sizeof(*ps.faces));
    ps.face_cube = (int32_t *)ARENA_CALLOC(
        arena, ps.nf, sizeof(*ps.face_cube));
    ps.ids = (char (*)[48])ARENA_CALLOC(arena, 1, 48);
    memcpy(ps.ids[0], "atlas_track_grow", 16);
    ps.cube_voff = (size_t *)ARENA_ALLOC(
        arena, 2 * sizeof(*ps.cube_voff));
    ps.cube_voff[0] = 0;
    ps.cube_voff[1] = ps.nv;
    ps.cube_org = (long (*)[3])ARENA_CALLOC(arena, 3, sizeof(long));
    ps.u_min = ps.v_min = DBL_MAX;
    ps.u_max = ps.v_max = -DBL_MAX;

    for (size_t i = 0; i < ps.nv; i++) {
        for (int d = 0; d < 3; d++)
            ps.verts[i * 3 + (size_t)d] =
                (float)result->xyz[i * 3 + (size_t)d];
        double u = result->global_uv[i * 2];
        double v = result->global_uv[i * 2 + 1];
        ps.uv[i * 2] = (float)u;
        ps.uv[i * 2 + 1] = (float)v;
        if (u < ps.u_min) ps.u_min = u;
        if (u > ps.u_max) ps.u_max = u;
        if (v < ps.v_min) ps.v_min = v;
        if (v > ps.v_max) ps.v_max = v;
        ps.gid[i] = result->vertex_component[i];

        int32_t target = result->source_target[i];
        double p[2];
        if (target >= 0 && (size_t)target < evidence->ntarget) {
            const AtlasRibbonTarget *sample = &evidence->target[target];
            p[0] = sample->p[0];
            p[1] = sample->p[1];
        } else if (target == -1) {
            double relative[3] = {
                result->xyz[i * 3] - evidence->axis_point[0],
                result->xyz[i * 3 + 1] - evidence->axis_point[1],
                result->xyz[i * 3 + 2] - evidence->axis_point[2]
            };
            p[0] = relative[0] * evidence->basis0[0] +
                   relative[1] * evidence->basis0[1] +
                   relative[2] * evidence->basis0[2];
            p[1] = relative[0] * evidence->basis1[0] +
                   relative[1] * evidence->basis1[1] +
                   relative[2] * evidence->basis1[2];
        } else {
            return -1;
        }
        double theta = atan2(p[1], p[0]);
        double alpha = ate_mod_tau(
            (double)result->winding_direction * theta) / SCAFFOLD_2PI;
        double turn = (double)result->vertex_winding[i] + alpha;
        /* ExportAtlas bins on sense*phi/(2*pi), so synthesize the registered
         * phi that exactly represents the consumer's snapped turn. */
        ps.phi[i] = (float)((double)calibration->sense *
                            SCAFFOLD_2PI * turn);
    }
    memcpy(ps.faces, result->faces,
           ps.nf * 3 * sizeof(*ps.faces));

    AtlasOpts options;
    AtlasOpts_default(&options);
    options.piece_wraps = wraps_per_piece;
    options.slab_v = slab_v;
    options.du = du;
    options.dv = dv;
    options.write_winding = 1;
    options.verbose = 1;
    AtlasStats stats;
    if (ExportAtlas_run(
            arena, &ps, calibration, output_root, prefix,
            &options, &stats) != 0)
        return -1;
    out->pieces = stats.n_pieces;
    out->written = stats.n_written;
    out->empty = stats.n_empty;
    out->over_cap = stats.n_over_cap;
    out->quarantine_faces = stats.n_quarantine_faces;
    out->valid_pixels = stats.total_valid_px;
    out->multi_pixels = stats.total_multi_px;
    out->conflict_pixels = stats.total_conflict_px;
    out->conflict_fraction = stats.conflict_frac;
    out->seconds = stats.seconds;
    return 0;
}
