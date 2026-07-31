#include "snap_quilt.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../common/gco_wrap.h"

enum {
    SQ_GOOD = 0,
    SQ_FIXABLE = 1,
    SQ_CRACK = 2,
    SQ_TARGET_CLOSEST = 0,
    SQ_TARGET_QUILT_ARGMIN = 1,
    SQ_TARGET_QUILT_MRF = 2
};

#define SQ_INVALID_COST 1000000
#define SQ_MAX_VALID_COST 8000
#define SQ_MAX_HITS 16

static double sq_clamp(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static int sq_make_direction(const float *p, const float *normal,
                             const SnapQuiltOpts *o, double out[3])
{
    double axis[3] = { o->axis_dir[0], o->axis_dir[1], o->axis_dir[2] };
    double al = sqrt(axis[0] * axis[0] + axis[1] * axis[1]
                     + axis[2] * axis[2]);
    double d[3], axial = 0.0, len = 0.0;
    int k = 0;

    if (al < 1e-12) {
        axis[0] = 1.0;
        axis[1] = 0.0;
        axis[2] = 0.0;
    } else {
        for (k = 0; k < 3; k++) axis[k] /= al;
    }
    for (k = 0; k < 3; k++) {
        d[k] = (double)p[k] - (double)o->axis_point[k];
        axial += d[k] * axis[k];
    }
    for (k = 0; k < 3; k++) {
        out[k] = d[k] - axial * axis[k];
        len += out[k] * out[k];
    }
    len = sqrt(len);
    if (len > 1e-9) {
        for (k = 0; k < 3; k++) out[k] /= len;
        return 0;
    }

    /* The radial direction is undefined on the axis.  This is rare, but a
     * usable mesh normal is a safer fallback than manufacturing a direction. */
    if (normal != NULL) {
        len = sqrt((double)normal[0] * normal[0]
                   + (double)normal[1] * normal[1]
                   + (double)normal[2] * normal[2]);
        if (len > 0.5) {
            for (k = 0; k < 3; k++) out[k] = (double)normal[k] / len;
            return 0;
        }
    }
    out[0] = out[1] = out[2] = 0.0;
    return -1;
}

/* Return the free ray length.  The occupancy tree may contain the source's
 * local patch, so hits whose ORIGINAL position is within local_r of p do not
 * block.  KDTree_ball_query returns the total hit count even when the output
 * buffer is shorter; always cap before indexing the fixed buffer. */
static int sq_is_local_one_ring(const int32_t *off, const int32_t *tgt,
                                size_t vi, size_t base, int32_t hit)
{
    size_t gj = hit >= 0 ? (size_t)hit : SIZE_MAX;
    if (gj == base + vi) return 1;
    if (gj < base) return 0;
    gj -= base;
    for (int32_t e = off[vi]; e < off[vi + 1]; e++)
        if ((size_t)tgt[e] == gj) return 1;
    return 0;
}

static double sq_free_reach(KDTree_T occ, const float *occ_verts,
                             const int32_t *off, const int32_t *tgt,
                             size_t vi, size_t base,
                             const double p[3], const double d[3], double sign,
                            double reach, double step, double occ2,
                            double local2, double *out_block)
{
    double t = 0.0;

    *out_block = reach;
    for (t = step; t <= reach + 1e-9; t += step) {
        float q[3];
        int32_t hit[SQ_MAX_HITS];
        size_t nh = 0, h = 0;
        int blocked = 0;

        q[0] = (float)(p[0] + sign * t * d[0]);
        q[1] = (float)(p[1] + sign * t * d[1]);
        q[2] = (float)(p[2] + sign * t * d[2]);
        nh = KDTree_ball_query(occ, q, (float)occ2, hit, SQ_MAX_HITS);
        if (nh > SQ_MAX_HITS) nh = SQ_MAX_HITS;
        for (h = 0; h < nh; h++) {
            size_t bi = (size_t)hit[h];
            double ez = (double)occ_verts[bi * 3 + 0] - p[0];
            double ey = (double)occ_verts[bi * 3 + 1] - p[1];
            double ex = (double)occ_verts[bi * 3 + 2] - p[2];
            double dsq = ez * ez + ey * ey + ex * ex;
            /* Spatial proximity alone cannot identify the source patch: an
             * adjacent recto/verso ply can be closer than local_r.  Ignore only
             * the source and its actual mesh 1-ring; every other hit is another
             * geometric constraint, even when very close. */
            if (dsq > local2
                || !sq_is_local_one_ring(off, tgt, vi, base, hit[h])) {
                blocked = 1;
                break;
            }
        }
        if (blocked) {
            *out_block = t;
            return t > step ? t - step : 0.0;
        }
    }
    return reach;
}

static int sq_is_candidate(CubeTable *ct, const double p[3],
                           const double d[3], double depth, double band,
                           double cv, double min_gain, double *out_sample)
{
    double q[3] = {
        p[0] + depth * d[0],
        p[1] + depth * d[1],
        p[2] + depth * d[2]
    };
    double s = sample_trilinear(ct, q[0], q[1], q[2]);

    *out_sample = s;
    return s >= 0.0 && s > band && s - cv > min_gain;
}

static void sq_propagate_boundary_reference(
    Arena_T arena, size_t nv, size_t nsite, const int32_t *site_of,
    const int32_t *vert_of, const int32_t *off, const int32_t *tgt,
    const double *cv, const uint8_t *has, const uint8_t *dark, double band,
    double *ref, uint8_t *has_ref)
{
    int32_t *queue = (int32_t *)ARENA_ALLOC(
        arena, (long)(nsite * sizeof(int32_t)));
    size_t head = 0, tail = 0, si = 0;

    (void)nv;
    for (si = 0; si < nsite; si++) {
        int32_t vi = vert_of[si];
        double sum = 0.0;
        size_t n = 0;
        int32_t e = 0;
        for (e = off[vi]; e < off[vi + 1]; e++) {
            int32_t j = tgt[e];
            if (!dark[j] && has[j]) {
                sum += cv[j];
                n++;
            }
        }
        if (n > 0) {
            ref[si] = sum / (double)n;
            has_ref[si] = 1;
            queue[tail++] = (int32_t)si;
        }
    }

    /* Multi-source graph propagation carries the nearest clean seam's value
     * through the dark island.  It is deliberately not an intensity maximizer:
     * it asks which candidate continues the clean boundary texture. */
    while (head < tail) {
        int32_t si0 = queue[head++];
        int32_t vi = vert_of[si0];
        int32_t e = 0;
        for (e = off[vi]; e < off[vi + 1]; e++) {
            int32_t sj = site_of[tgt[e]];
            if (sj >= 0 && !has_ref[sj]) {
                ref[sj] = ref[si0];
                has_ref[sj] = 1;
                queue[tail++] = sj;
            }
        }
    }
    for (si = 0; si < nsite; si++)
        if (!has_ref[si]) ref[si] = band;
}

int SnapQuilt_select(Arena_T arena,
                     CubeTable *ct,
                     KDTree_T occ,
                     const float *occ_verts,
                     size_t occ_index_base,
                     const float *verts,
                     const float *normals,
                     size_t nv,
                     CSR_T adj,
                     const double *cv,
                     const uint8_t *has,
                     const uint8_t *dark,
                     const SnapQuiltOpts *opts,
                     uint8_t *vclass,
                     float *voff,
                     float *vdir,
                     float *vgain,
                     float *vblock,
                     SnapQuiltStats *out)
{
    Arena_Mark mark;
    const int32_t *off = NULL, *tgt = NULL;
    int32_t *site_of = NULL, *vert_of = NULL;
    double *ref = NULL;
    float *dirs = NULL, *energy = NULL;
    uint8_t *valid = NULL, *has_ref = NULL;
    int32_t *data = NULL, *smooth = NULL;
    int *labels = NULL, *best_label = NULL, *close_label = NULL;
    size_t nsite = 0, si = 0, ncost = 0;
    int bins = 0, zero = 0, mode = 0;
    double reach = 0.0, step = 0.0, depth_step = 0.0;
    double occ2 = 0.0, local2 = 0.0, tau = 0.0;
    double max_cost = 1.0, scale = 1.0;
    GCO_Handle gc = NULL;

    assert(arena && ct && occ && occ_verts && verts && adj && cv && has
           && dark && opts && vclass && voff && vdir && vgain && vblock && out);
    memset(out, 0, sizeof *out);
    if (nv == 0) return 0;
    if (nv > (size_t)INT_MAX) return -1;

    mark = Arena_save(arena);
    off = CSR_offset(adj);
    tgt = CSR_target(adj);
    reach = opts->reach > 0.0 ? opts->reach : 8.0;
    step = opts->step > 0.0 ? opts->step : 0.25;
    occ2 = opts->occ_thresh > 0.0 ? opts->occ_thresh : 0.9;
    occ2 *= occ2;
    local2 = opts->local_r > 0.0 ? opts->local_r : 2.0;
    local2 *= local2;
    tau = opts->smooth_tau > 0.0 ? opts->smooth_tau : 4.0;
    mode = opts->target_mode;
    if (mode < SQ_TARGET_CLOSEST || mode > SQ_TARGET_QUILT_MRF)
        mode = SQ_TARGET_QUILT_MRF;
    bins = opts->depth_bins >= 3 ? opts->depth_bins : 33;
    if ((bins & 1) == 0) bins++;
    if (bins > 257) bins = 257;
    zero = bins / 2;
    depth_step = 2.0 * reach / (double)(bins - 1);

    site_of = (int32_t *)ARENA_ALLOC(
        arena, (long)(nv * sizeof(int32_t)));
    for (si = 0; si < nv; si++) site_of[si] = -1;
    for (si = 0; si < nv; si++) {
        if (dark[si] && has[si] && off[si] < off[si + 1]) {
            site_of[si] = (int32_t)nsite;
            nsite++;
        }
    }
    if (nsite == 0) {
        Arena_restore(arena, mark);
        return 0;
    }
    if (nsite > (size_t)INT_MAX
        || nsite > SIZE_MAX / (size_t)bins
        || nsite * (size_t)bins > (size_t)LONG_MAX / sizeof(float)) {
        Arena_restore(arena, mark);
        return -1;
    }
    ncost = nsite * (size_t)bins;
    vert_of = (int32_t *)ARENA_ALLOC(
        arena, (long)(nsite * sizeof(int32_t)));
    for (si = 0; si < nv; si++)
        if (site_of[si] >= 0) vert_of[site_of[si]] = (int32_t)si;

    ref = (double *)ARENA_ALLOC(arena, (long)(nsite * sizeof(double)));
    has_ref = (uint8_t *)ARENA_CALLOC(arena, (long)nsite, 1);
    dirs = (float *)ARENA_CALLOC(arena, (long)nsite, 3 * sizeof(float));
    energy = (float *)ARENA_ALLOC(arena, (long)(ncost * sizeof(float)));
    valid = (uint8_t *)ARENA_CALLOC(arena, (long)ncost, 1);
    labels = (int *)ARENA_ALLOC(arena, (long)(nsite * sizeof(int)));
    best_label = (int *)ARENA_ALLOC(arena, (long)(nsite * sizeof(int)));
    close_label = (int *)ARENA_ALLOC(arena, (long)(nsite * sizeof(int)));
    sq_propagate_boundary_reference(arena, nv, nsite, site_of, vert_of,
                                    off, tgt, cv, has, dark, opts->band,
                                    ref, has_ref);

    for (si = 0; si < nsite; si++) {
        int32_t vi = vert_of[si];
        const float *normal = normals != NULL ? &normals[(size_t)vi * 3] : NULL;
        double p[3] = { verts[(size_t)vi * 3 + 0],
                        verts[(size_t)vi * 3 + 1],
                        verts[(size_t)vi * 3 + 2] };
        double dir[3] = { 0.0, 0.0, 0.0 };
        double clear_minus = 0.0, clear_plus = 0.0;
        double block_minus = reach, block_plus = reach;
        double best = 1e30, closest = 1e30, closest_cost = 1e30;
        int best_l = zero, close_l = zero;
        int plus = 0, minus = 0, have_candidate = 0;
        int boundary_edges = 0, l = 0;
        int32_t e = 0;

        out->n_dark++;
        if (sq_make_direction(&verts[(size_t)vi * 3], normal, opts, dir) != 0) {
            valid[si * (size_t)bins + (size_t)zero] = 1;
            energy[si * (size_t)bins + (size_t)zero] = 0.0f;
            labels[si] = best_label[si] = close_label[si] = zero;
            vblock[vi] = (float)reach;
            continue;
        }
        dirs[si * 3 + 0] = (float)dir[0];
        dirs[si * 3 + 1] = (float)dir[1];
        dirs[si * 3 + 2] = (float)dir[2];
        clear_minus = sq_free_reach(occ, occ_verts, off, tgt, (size_t)vi,
                                    occ_index_base, p, dir, -1.0, reach,
                                    step, occ2, local2, &block_minus);
        clear_plus = sq_free_reach(occ, occ_verts, off, tgt, (size_t)vi,
                                   occ_index_base, p, dir, 1.0, reach,
                                   step, occ2, local2, &block_plus);
        vblock[vi] = (float)(block_minus < block_plus ? block_minus : block_plus);
        for (e = off[vi]; e < off[vi + 1]; e++)
            if (!dark[tgt[e]]) boundary_edges++;

        for (l = 0; l < bins; l++) {
            double depth = ((double)l - (double)zero) * depth_step;
            double ad = fabs(depth), s = -1.0, c = 0.0;
            size_t ci = si * (size_t)bins + (size_t)l;
            double free_reach = depth < 0.0 ? clear_minus : clear_plus;

            energy[ci] = 0.0f;
            if (l == zero || ad > free_reach + 1e-9) continue;
            if (!sq_is_candidate(ct, p, dir, depth, opts->band, cv[vi],
                                 opts->min_gain, &s))
                continue;
            /* No clean boundary means no quilting observation.  In that
             * anchorless case distance is the conservative unary: inventing a
             * reference from the threshold can make a farther verso transition
             * look artificially perfect. */
            c = opts->w_close * ad;
            if (has_ref[si]) c += opts->w_match * fabs(s - ref[si]);
            if (boundary_edges > 0 && opts->smooth_mu > 0.0)
                c += (double)boundary_edges * opts->smooth_mu
                     * (ad < tau ? ad : tau);
            if (c < 0.0) c = 0.0;
            energy[ci] = (float)c;
            valid[ci] = 1;
            have_candidate = 1;
            if (depth > 0.0) plus = 1;
            else minus = 1;
            if (c < best) {
                best = c;
                best_l = l;
            }
            if (ad < closest - 1e-9
                || (fabs(ad - closest) <= 1e-9 && c < closest_cost)) {
                closest = ad;
                closest_cost = c;
                close_l = l;
            }
            if (c > max_cost) max_cost = c;
        }

        if (!have_candidate) {
            size_t ci = si * (size_t)bins + (size_t)zero;
            valid[ci] = 1;
            energy[ci] = 0.0f;
            best_l = close_l = zero;
            if (vblock[vi] < (float)(reach - 1e-6)) out->n_blocked++;
        } else if (plus && minus) {
            out->n_bidir++;
        }
        best_label[si] = best_l;
        close_label[si] = close_l;
        labels[si] = mode == SQ_TARGET_CLOSEST ? close_l : best_l;
    }

    if (opts->smooth_mu > 0.0) {
        double sm = opts->smooth_mu * tau;
        if (sm > max_cost) max_cost = sm;
    }
    scale = (double)SQ_MAX_VALID_COST / max_cost;
    scale = sq_clamp(scale, 1e-6, 64.0);

    data = (int32_t *)ARENA_ALLOC(
        arena, (long)(ncost * sizeof(int32_t)));
    smooth = (int32_t *)ARENA_ALLOC(
        arena, (long)((size_t)bins * (size_t)bins * sizeof(int32_t)));
    for (si = 0; si < ncost; si++) {
        if (!valid[si]) {
            data[si] = SQ_INVALID_COST;
        } else {
            double q = floor((double)energy[si] * scale + 0.5);
            if (q > SQ_MAX_VALID_COST) q = SQ_MAX_VALID_COST;
            data[si] = (int32_t)q;
        }
    }
    /* Build an exactly metric integer truncated-linear matrix.  Rounding every
     * floating pair independently can violate the triangle inequality by one,
     * which makes alpha-expansion reject an otherwise metric model. */
    {
        int64_t unit = (int64_t)floor(
            opts->smooth_mu * depth_step * scale + 0.5);
        int64_t cap = (int64_t)floor(
            opts->smooth_mu * tau * scale + 0.5);
        if (unit < 1 && opts->smooth_mu > 0.0) unit = 1;
        if (cap > SQ_MAX_VALID_COST) cap = SQ_MAX_VALID_COST;
        if (cap < 0) cap = 0;
    for (int l0 = 0; l0 < bins; l0++) {
        for (int l1 = 0; l1 < bins; l1++) {
            int64_t q = (int64_t)abs(l0 - l1) * unit;
            if (q > cap) q = cap;
            smooth[(size_t)l0 * (size_t)bins + (size_t)l1] =
                (int32_t)q;
        }
    }
    }

    if (mode == SQ_TARGET_QUILT_MRF && bins >= 2) {
        gc = GCO_create((int)nsite, bins);
        if (gc != NULL) {
            GCO_set_data_cost(gc, data);
            GCO_set_smooth_cost(gc, smooth);
            for (si = 0; si < nsite; si++) {
                int32_t vi = vert_of[si];
                int32_t e = 0;
                GCO_set_label(gc, (int)si, labels[si]);
                for (e = off[vi]; e < off[vi + 1]; e++) {
                    int32_t sj = site_of[tgt[e]];
                    if (sj > (int32_t)si)
                        GCO_set_neighbor(gc, (int)si, sj, 1);
                }
            }
            if (GCO_expansion(gc, -1) >= 0) {
                GCO_get_labels(gc, labels, (int)nsite);
            } else {
                out->n_gco_fallback++;
                for (si = 0; si < nsite; si++) labels[si] = best_label[si];
            }
            GCO_destroy(gc);
            gc = NULL;
        } else {
            out->n_gco_fallback++;
            for (si = 0; si < nsite; si++) labels[si] = best_label[si];
        }
    }

    for (si = 0; si < nsite; si++) {
        int32_t vi = vert_of[si];
        int l = labels[si];
        double depth = 0.0, p[3], dir[3], s = -1.0;
        size_t ci = 0;

        if (mode == SQ_TARGET_CLOSEST) l = close_label[si];
        else if (mode == SQ_TARGET_QUILT_ARGMIN) l = best_label[si];
        if (l < 0 || l >= bins
            || !valid[si * (size_t)bins + (size_t)l])
            l = best_label[si];
        ci = si * (size_t)bins + (size_t)l;
        if (l == zero || !valid[ci]) {
            vclass[vi] = SQ_CRACK;
            out->n_crack++;
            continue;
        }
        depth = ((double)l - (double)zero) * depth_step;
        p[0] = verts[(size_t)vi * 3 + 0];
        p[1] = verts[(size_t)vi * 3 + 1];
        p[2] = verts[(size_t)vi * 3 + 2];
        dir[0] = dirs[si * 3 + 0];
        dir[1] = dirs[si * 3 + 1];
        dir[2] = dirs[si * 3 + 2];
        if (!sq_is_candidate(ct, p, dir, depth, opts->band, cv[vi],
                             opts->min_gain, &s)) {
            vclass[vi] = SQ_CRACK;
            out->n_crack++;
            continue;
        }
        voff[vi] = (float)fabs(depth);
        vdir[(size_t)vi * 3 + 0] =
            (float)((depth < 0.0 ? -1.0 : 1.0) * dir[0]);
        vdir[(size_t)vi * 3 + 1] =
            (float)((depth < 0.0 ? -1.0 : 1.0) * dir[1]);
        vdir[(size_t)vi * 3 + 2] =
            (float)((depth < 0.0 ? -1.0 : 1.0) * dir[2]);
        vgain[vi] = (float)(s - cv[vi]);
        vclass[vi] = SQ_FIXABLE;
        out->n_fixable++;
        out->mean_quilt_cost += energy[ci];
        out->mean_distance += fabs(depth);
    }
    if (out->n_fixable > 0) {
        out->mean_quilt_cost /= (double)out->n_fixable;
        out->mean_distance /= (double)out->n_fixable;
    }

    (void)opts->verbose; /* Logging is owned by the public snap stages. */
    Arena_restore(arena, mark);
    return 0;
}
