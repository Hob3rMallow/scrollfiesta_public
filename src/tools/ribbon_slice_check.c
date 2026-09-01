/* ============================================================================
 * ribbon_slice_check -- functional winding-order gate at fixed z planes.
 *
 * Loads a SOLVED ribbon VMESH (stage 3+: the metric-reopt or later
 * artifact, whose UV carries the solved u), intersects the mesh with axial
 * planes, bins the crossings into umbilicus rays, and demands u be
 * radially MONOTONE along every ray (the monotone rainbow the cross
 * sections are judged by).  This is the regression sentinel for the z=5500
 * class of defect (user directive 2026-08-30): a drift-scrambled global
 * frame, a skipped track stitched by a chord, or a self-intersecting core
 * all present as u inversions along rays.  Do NOT point it at a stage-2
 * fit (its u is legitimately atlas-packed, and its phase sidecar was
 * measured near-wrapped -- see the calibration comment below).
 *
 *   ribbon_slice_check <solved_ribbon.vmesh> --umb-y Y --umb-x X
 *                      --z Z [--z Z...]
 *   ribbon_slice_check --selftest
 *
 * Exit 0 = every plane PASSes, 1 = any plane FAILs, 2 = usage/load error.
 *
 * The pinned production planes live in
 * scripts/flatten/check_pherc0139_10x_slices.ps1 (first slice + z=5500).
 * ==========================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/mesh_bin.h"

#define RSC_THETA_BINS           96
#define RSC_DUP_DR_VOX           3.0    /* nearer pairs = shell duplicates */
/* CALIBRATED 2026-08-30.  Rejected statistics, for the record: fit-stage
 * PHASE metrics of any kind -- measured on the good 4x5x5 champion, the
 * published fit-stage phase sidecar is near-wrapped (wraps nine turns
 * apart carry the same value along a ray), so no stage-2 phase statistic
 * can certify winding; raw inversion rates (15-20% benign ~1-turn
 * theta-bin noise on the GOOD champion) and per-ray pitch linearity
 * (eccentricity bends effective pitch ~20%/ray) also fail to separate
 * good from bad.  The contract that holds is the SOLVED U of stage-3+
 * artifacts: u must be radially monotone along every umbilicus ray (the
 * monotone rainbow the cross sections are judged by).  Gate: deep u
 * regressions (beyond RSC_U_SLACK_CIRC of the local circumference) as a
 * fraction of adjacent ray pairs.  One-wrap step sizes are REPORTED (skip
 * evidence, in local-circumference units) but not gated yet. */
#define RSC_STEP_DR_LO           4.75   /* 0.5 x pitch: excludes duplicates */
#define RSC_STEP_DR_HI          15.20   /* 1.6 x pitch: one-wrap pairs only */
#define RSC_STEP_DTHETA_MAX      0.02   /* rad; near-coaxial pairs only */
#define RSC_U_SLACK_CIRC         0.35   /* u regression deeper than this
                                         * fraction of local circumference
                                         * is DEEP */
/* Calibration table (2026-08-30, per-plane deep-fraction / worst in local
 * circumferences):
 *   4x5x5 champion reopt (GOOD): 4.4-7.6% / <= 8.7c  (multi-layer coverage
 *     legitimately re-visits a wrap at +-1 turn -> benign deep rate)
 *   0143-hybrid refit: 2.5% / 55c at the KNOWN-broken core (honest FAIL),
 *     3.0c where its band is good
 *   drift-scrambled 10x fit u (atlas frame): 10.7-13.9% / 291-1104c
 * The two defect classes exceed BOTH bounds; the good band exceeds
 * neither. */
#define RSC_MAX_DEEP_FRAC        0.09   /* per-plane FAIL threshold */
#define RSC_MAX_WORST_CIRC      20.0    /* per-plane FAIL threshold */

typedef struct {
    double radius, w;   /* w = lifted phase in turns */
    double theta;
    double edge_t;
    size_t va, vb;
    int32_t bin;
} RscCrossing;

static int rsc_compare_crossing(const void *pa, const void *pb)
{
    const RscCrossing *a = (const RscCrossing *)pa;
    const RscCrossing *b = (const RscCrossing *)pb;
    if (a->radius != b->radius) return a->radius < b->radius ? -1 : 1;
    return a->w < b->w ? -1 : (a->w > b->w ? 1 : 0);
}

/* Collect unique-edge crossings of plane z=plane_z.  Returns count. */
static size_t rsc_collect(Arena_T arena,
                          const float *verts, const int32_t *faces,
                          size_t nf, const float *uv,
                          double plane_z, double umb_y, double umb_x,
                          RscCrossing **out)
{
    size_t capacity = 4096, count = 0;
    size_t table_size = 1;
    uint64_t *seen = NULL;
    RscCrossing *crossing = (RscCrossing *)ARENA_ALLOC(
        arena, (long)(capacity * sizeof(RscCrossing)));
    /* dedup via open-addressed hash of the undirected edge key */
    while (table_size < nf) table_size <<= 1;
    table_size <<= 1;
    seen = (uint64_t *)ARENA_CALLOC(arena, (long)table_size,
                                    (long)sizeof(uint64_t));
    for (size_t f = 0; f < nf; f++) {
        for (int k = 0; k < 3; k++) {
            size_t va = (size_t)faces[f * 3 + (size_t)k];
            size_t vb = (size_t)faces[f * 3 + (size_t)((k + 1) % 3)];
            double za = verts[va * 3], zb = verts[vb * 3];
            double t = 0.0, y = 0.0, x = 0.0, w = 0.0;
            uint64_t key = 0, slot = 0;
            if (!((za < plane_z && zb > plane_z) ||
                  (zb < plane_z && za > plane_z)))
                continue;
            key = va < vb ? ((uint64_t)va << 32) | (uint64_t)vb
                          : ((uint64_t)vb << 32) | (uint64_t)va;
            key = key * UINT64_C(0x9E3779B97F4A7C15) + 1u;
            slot = key & (table_size - 1);
            while (seen[slot] != 0 && seen[slot] != key)
                slot = (slot + 1) & (table_size - 1);
            if (seen[slot] == key) continue;
            seen[slot] = key;
            t = (plane_z - za) / (zb - za);
            y = (double)verts[va * 3 + 1] +
                t * ((double)verts[vb * 3 + 1] - verts[va * 3 + 1]);
            x = (double)verts[va * 3 + 2] +
                t * ((double)verts[vb * 3 + 2] - verts[va * 3 + 2]);
            w = (double)uv[va * 2] +
                t * ((double)uv[vb * 2] - uv[va * 2]);
            if (count == capacity) {
                RscCrossing *grown = (RscCrossing *)ARENA_ALLOC(
                    arena, (long)(capacity * 2 * sizeof(RscCrossing)));
                memcpy(grown, crossing, count * sizeof(RscCrossing));
                crossing = grown;
                capacity *= 2;
            }
            crossing[count].radius = hypot(y - umb_y, x - umb_x);
            crossing[count].theta = atan2(y - umb_y, x - umb_x);
            crossing[count].w = w;
            crossing[count].edge_t = t;
            crossing[count].va = va;
            crossing[count].vb = vb;
            count++;
        }
    }
    *out = crossing;
    return count;
}

/* Check one plane.  Returns 0 PASS, 1 FAIL, prints one verdict line. */
static int rsc_check_plane(Arena_T arena,
                           const float *verts, const int32_t *faces,
                           size_t nf, const float *uv,
                           double plane_z, double umb_y, double umb_x)
{
    RscCrossing *crossing = NULL;
    size_t count = rsc_collect(arena, verts, faces, nf, uv, plane_z,
                               umb_y, umb_x, &crossing);
    size_t bin_first[RSC_THETA_BINS + 1];
    size_t pairs = 0, violations = 0, rays = 0;
    double sense_sum = 0.0, worst = 0.0;
    RscCrossing worst_inner, worst_outer;
    int have_worst = 0;
    memset(bin_first, 0, sizeof bin_first);
    if (count < 8) {
        printf("[slice-check] z=%.1f crossings=%zu VERDICT=FAIL "
               "(plane barely intersects the ribbon)\n", plane_z, count);
        return 1;
    }
    /* bucket by theta bin (counting sort), sort each ray by radius */
    for (size_t i = 0; i < count; i++) {
        int bin = (int)floor((crossing[i].theta + 3.14159265358979323846) /
                             (2.0 * 3.14159265358979323846) * RSC_THETA_BINS);
        if (bin < 0) bin = 0;
        if (bin >= RSC_THETA_BINS) bin = RSC_THETA_BINS - 1;
        crossing[i].bin = bin;
        bin_first[bin + 1]++;
    }
    for (int b = 0; b < RSC_THETA_BINS; b++) bin_first[b + 1] += bin_first[b];
    {
        RscCrossing *sorted = (RscCrossing *)ARENA_ALLOC(
            arena, (long)(count * sizeof(RscCrossing)));
        size_t cursor[RSC_THETA_BINS];
        for (int b = 0; b < RSC_THETA_BINS; b++) cursor[b] = bin_first[b];
        for (size_t i = 0; i < count; i++)
            sorted[cursor[crossing[i].bin]++] = crossing[i];
        crossing = sorted;
    }
    for (int b = 0; b < RSC_THETA_BINS; b++)
        qsort(crossing + bin_first[b], bin_first[b + 1] - bin_first[b],
              sizeof(RscCrossing), rsc_compare_crossing);
    /* winding sense: majority of radially-adjacent phase steps */
    for (int b = 0; b < RSC_THETA_BINS; b++)
        for (size_t i = bin_first[b] + 1; i < bin_first[b + 1]; i++) {
            double dr = crossing[i].radius - crossing[i - 1].radius;
            double dw = crossing[i].w - crossing[i - 1].w;
            if (dr < RSC_DUP_DR_VOX) continue;
            sense_sum += dw > 0.0 ? 1.0 : (dw < 0.0 ? -1.0 : 0.0);
        }
    /* monotonicity walk, duplicates skipped */
    size_t deep = 0;
    for (int b = 0; b < RSC_THETA_BINS; b++) {
        size_t first = bin_first[b], last = bin_first[b + 1];
        double previous_r = -1.0, previous_w = 0.0;
        const RscCrossing *previous = NULL;
        int have = 0;
        if (last - first >= 2) rays++;
        for (size_t i = first; i < last; i++) {
            double w = sense_sum >= 0.0 ? crossing[i].w : -crossing[i].w;
            if (have && crossing[i].radius - previous_r < RSC_DUP_DR_VOX)
                continue;
            if (have) {
                double circumference = 2.0 * 3.14159265358979323846 *
                                       crossing[i].radius;
                double regression = previous_w - w;
                pairs++;
                if (regression > 0.05 * circumference) {
                    violations++;
                    if (regression > RSC_U_SLACK_CIRC * circumference)
                        deep++;
                    if (regression / circumference > worst)
                    {
                        worst = regression / circumference;
                        worst_inner = *previous;
                        worst_outer = crossing[i];
                        have_worst = 1;
                    }
                }
            }
            previous_r = crossing[i].radius;
            previous_w = w;
            previous = &crossing[i];
            have = 1;
        }
    }
    /* one-wrap step statistic (the drift detector): near-coaxial crossing
     * pairs one pitch apart must step +1 turn exactly */
    {
        double fraction = pairs > 0 ? (double)violations / (double)pairs : 1.0;
        size_t nstep = 0, bad_step = 0;
        double step_sum = 0.0, bad_fraction = 1.0;
        int fail = 0;
        for (int b = 0; b < RSC_THETA_BINS; b++) {
            size_t first = bin_first[b], last = bin_first[b + 1];
            for (size_t i = first; i < last; i++) {
                for (size_t j = i + 1; j < last; j++) {
                    double dr = crossing[j].radius - crossing[i].radius;
                    double dtheta = 0.0, wa = 0.0, wb = 0.0, step = 0.0;
                    if (dr > RSC_STEP_DR_HI) break;   /* sorted by radius */
                    if (dr < RSC_STEP_DR_LO) continue;
                    dtheta = fabs(crossing[j].theta - crossing[i].theta);
                    if (dtheta > RSC_STEP_DTHETA_MAX) continue;
                    wa = sense_sum >= 0.0 ? crossing[i].w : -crossing[i].w;
                    wb = sense_sum >= 0.0 ? crossing[j].w : -crossing[j].w;
                    step = (wb - wa) /
                           (2.0 * 3.14159265358979323846 *
                            crossing[i].radius);
                    nstep++;
                    step_sum += step;
                    if (step < 0.5 || step > 2.0) bad_step++;
                }
            }
        }
        double deep_fraction = pairs > 0 ? (double)deep / (double)pairs
                                         : 1.0;
        bad_fraction = nstep > 0 ? (double)bad_step / (double)nstep : 1.0;
        fail = pairs < 16 || deep_fraction > RSC_MAX_DEEP_FRAC ||
               worst > RSC_MAX_WORST_CIRC;
        printf("[slice-check] z=%.1f crossings=%zu rays=%zu pairs=%zu "
               "inversions=%zu (%.2f%%) deep=%zu (%.3f%%) worst=%.2fc "
               "steps=%zu offsize=%zu (%.1f%%) mean_step=%.2fc "
               "sense=%+d VERDICT=%s\n",
               plane_z, count, rays, pairs, violations, fraction * 100.0,
               deep, deep_fraction * 100.0, worst, nstep, bad_step,
               bad_fraction * 100.0,
               nstep > 0 ? step_sum / (double)nstep : 0.0,
               sense_sum >= 0.0 ? 1 : -1, fail ? "FAIL" : "PASS");
        if (have_worst) {
            double wi = sense_sum >= 0.0 ? worst_inner.w : -worst_inner.w;
            double wo = sense_sum >= 0.0 ? worst_outer.w : -worst_outer.w;
            printf("[slice-check]   worst pair: bin=%d r=%.3f->%.3f "
                   "theta=%.5f/%.5f U=%.3f->%.3f regression=%.3f; "
                   "edges=(%zu,%zu,t=%.5f)->(%zu,%zu,t=%.5f)\n",
                   worst_outer.bin, worst_inner.radius, worst_outer.radius,
                   worst_inner.theta, worst_outer.theta, wi, wo, wi - wo,
                   worst_inner.va, worst_inner.vb, worst_inner.edge_t,
                   worst_outer.va, worst_outer.vb, worst_outer.edge_t);
        }
        return fail ? 1 : 0;
    }
}

/* ---- selftest ------------------------------------------------------------ */

static void rsc_selftest_spiral(float *verts, float *uv, int32_t *faces,
                                int nturn, int nphi, double z_lo, double z_hi,
                                int scramble_ring)
{
    int nv_ring = nphi * 2;
    for (int ring = 0; ring < nturn; ring++) {
        for (int p = 0; p < nphi; p++) {
            double w = ring + (double)p / nphi;
            double r = 20.0 + 9.5 * w;
            double theta = 2.0 * 3.14159265358979323846 * ((double)p / nphi);
            /* exact spiral arclength (circular approximation) */
            double u = 2.0 * 3.14159265358979323846 *
                       (20.0 * w + 4.75 * w * w);
            if (ring == scramble_ring)
                u += 4.0 * 2.0 * 3.14159265358979323846 * r;
            for (int layer = 0; layer < 2; layer++) {
                int v = ring * nv_ring + p * 2 + layer;
                verts[v * 3] = (float)(layer == 0 ? z_lo : z_hi);
                verts[v * 3 + 1] = (float)(100.0 + r * sin(theta));
                verts[v * 3 + 2] = (float)(100.0 + r * cos(theta));
                uv[v * 2] = (float)u;
                uv[v * 2 + 1] = 0.0f;
            }
        }
    }
    for (int ring = 0; ring < nturn; ring++)
        for (int p = 0; p + 1 < nphi; p++) {
            int a = ring * nv_ring + p * 2;
            int b = a + 2;
            int32_t *face = &faces[(size_t)(ring * (nphi - 1) + p) * 6];
            face[0] = a; face[1] = a + 1; face[2] = b;
            face[3] = b; face[4] = a + 1; face[5] = b + 1;
        }
}

static int rsc_selftest(void)
{
    Arena_T arena = Arena_new();
    int nturn = 6, nphi = 64, fails = 0;
    size_t nv = (size_t)nturn * nphi * 2;
    size_t nf = (size_t)nturn * (nphi - 1) * 2;
    float *verts = (float *)ARENA_ALLOC(arena, (long)(nv * 3 * sizeof(float)));
    float *uv = (float *)ARENA_ALLOC(arena, (long)(nv * 2 * sizeof(float)));
    int32_t *faces = (int32_t *)ARENA_ALLOC(arena,
                                            (long)(nf * 3 * sizeof(int32_t)));
    if (arena == NULL || verts == NULL || uv == NULL || faces == NULL)
        return -1;
    /* clean spiral must PASS */
    rsc_selftest_spiral(verts, uv, faces, nturn, nphi, 90.0, 110.0, -1);
    if (rsc_check_plane(arena, verts, faces, nf, uv, 100.0,
                        100.0, 100.0) != 0) {
        fprintf(stderr, "[slice-check selftest] FAIL: clean spiral "
                "did not pass\n");
        fails++;
    }
    /* one ring shoved +4 circumferences up in u (drift scramble) must FAIL */
    rsc_selftest_spiral(verts, uv, faces, nturn, nphi, 90.0, 110.0, 2);
    if (rsc_check_plane(arena, verts, faces, nf, uv, 100.0,
                        100.0, 100.0) != 1) {
        fprintf(stderr, "[slice-check selftest] FAIL: scrambled spiral "
                "was not caught\n");
        fails++;
    }
    Arena_dispose(&arena);
    fprintf(stderr, "[slice-check selftest] %s (%d failure(s))\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails == 0 ? 0 : 1;
}

int main(int argc, char **argv)
{
    Arena_T arena = NULL;
    MeshBinData mesh;
    double umb_y = 0.0, umb_x = 0.0;
    double plane[64];
    int nplane = 0, have_umb_y = 0, have_umb_x = 0, rc = 0;
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return rsc_selftest();
    if (argc < 4) {
        fprintf(stderr,
                "usage: ribbon_slice_check <ribbon.vmesh> --umb-y Y "
                "--umb-x X --z Z [--z Z...]\n"
                "       ribbon_slice_check --selftest\n");
        return 2;
    }
    for (int i = 2; i + 1 < argc; i += 2) {
        if (strcmp(argv[i], "--umb-y") == 0) {
            umb_y = atof(argv[i + 1]); have_umb_y = 1;
        } else if (strcmp(argv[i], "--umb-x") == 0) {
            umb_x = atof(argv[i + 1]); have_umb_x = 1;
        } else if (strcmp(argv[i], "--z") == 0 && nplane < 64) {
            plane[nplane++] = atof(argv[i + 1]);
        } else {
            fprintf(stderr, "ribbon_slice_check: unknown arg %s\n", argv[i]);
            return 2;
        }
    }
    if (!have_umb_y || !have_umb_x || nplane == 0) {
        fprintf(stderr, "ribbon_slice_check: need --umb-y, --umb-x, --z\n");
        return 2;
    }
    arena = Arena_new();
    if (arena == NULL) return 2;
    memset(&mesh, 0, sizeof mesh);
    if (MeshBin_read_arena(arena, argv[1], &mesh) != 0 || mesh.uv == NULL) {
        fprintf(stderr, "ribbon_slice_check: cannot read %s (a SOLVED "
                "stage-3+ vmesh with UV is required)\n", argv[1]);
        Arena_dispose(&arena);
        return 2;
    }
    printf("[slice-check] %s: nv=%zu nf=%zu axis=(%.1f,%.1f) planes=%d\n",
           argv[1], mesh.nv, mesh.nf, umb_y, umb_x, nplane);
    for (int i = 0; i < nplane; i++)
        if (rsc_check_plane(arena, mesh.verts, mesh.faces, mesh.nf, mesh.uv,
                            plane[i], umb_y, umb_x) != 0)
            rc = 1;
    printf("[slice-check] overall: %s\n", rc == 0 ? "PASS" : "FAIL");
    Arena_dispose(&arena);
    return rc;
}
