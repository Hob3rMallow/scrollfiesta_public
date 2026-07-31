#include "atlas_turn_advance.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef ATA_2PI
#define ATA_2PI 6.283185307179586476925286766559
#endif
#ifndef ATA_PI
#define ATA_PI 3.14159265358979323846
#endif

static int ata_cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

static double ata_median(double *v, size_t n)
{
    if (n == 0) return 0.0;
    qsort(v, n, sizeof(double), ata_cmp_double);
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/* The pinned map, the one thing here that is not a measurement. */
static double ata_spiral_u(double a, double b, double phi)
{
    return a * phi + b * phi * phi / (4.0 * ATA_PI);
}

typedef struct {
    int32_t piece;
    int32_t turn;
    int32_t sample;
} AtaKey;

static int ata_cmp_key(const void *va, const void *vb)
{
    const AtaKey *a = (const AtaKey *)va;
    const AtaKey *b = (const AtaKey *)vb;
    if (a->piece != b->piece) return a->piece < b->piece ? -1 : 1;
    if (a->turn != b->turn) return a->turn < b->turn ? -1 : 1;
    return 0;
}

int AtlasTurnAdvance_bins(Arena_T arena,
                          const double *u, const double *phi,
                          const double *radius, const int32_t *piece,
                          size_t nsamples, size_t npieces, double sense,
                          AtlasTurnAdvanceResult *out)
{
    if (arena == NULL || out == NULL || u == NULL || phi == NULL ||
        radius == NULL || piece == NULL)
        return -1;
    memset(out, 0, sizeof *out);

    AtaKey *key = (AtaKey *)ARENA_ALLOC(
        arena, (nsamples ? nsamples : 1) * sizeof(AtaKey));
    size_t nkey = 0;
    for (size_t i = 0; i < nsamples; i++) {
        if (piece[i] < 0 || (size_t)piece[i] >= npieces) continue;
        key[nkey].piece = piece[i];
        key[nkey].turn = (int32_t)floor(sense * phi[i] / ATA_2PI);
        key[nkey].sample = (int32_t)i;
        nkey++;
    }
    if (nkey == 0) return 0;
    qsort(key, nkey, sizeof(AtaKey), ata_cmp_key);

    AtlasTurnAdvanceBin *bin = (AtlasTurnAdvanceBin *)ARENA_ALLOC(
        arena, nkey * sizeof(AtlasTurnAdvanceBin));
    double *scratch_u = (double *)ARENA_ALLOC(arena, nkey * sizeof(double));
    double *scratch_r = (double *)ARENA_ALLOC(arena, nkey * sizeof(double));
    size_t nbin = 0;
    {
        size_t i = 0;
        while (i < nkey) {
            size_t j = i;
            while (j < nkey && key[j].piece == key[i].piece &&
                   key[j].turn == key[i].turn) j++;
            size_t n = 0;
            for (size_t k = i; k < j; k++) {
                scratch_u[n] = u[key[k].sample];
                scratch_r[n] = radius[key[k].sample];
                n++;
            }
            bin[nbin].piece = key[i].piece;
            bin[nbin].turn = key[i].turn;
            bin[nbin].samples = n;
            bin[nbin].u_median = ata_median(scratch_u, n);
            bin[nbin].radius_median = ata_median(scratch_r, n);
            nbin++;
            i = j;
        }
    }

    /* Steps between ADJACENT turns of the same piece.  A gap in the turn
     * sequence is skipped rather than divided out: a piece that is absent for
     * two turns says nothing about how fast u should have advanced across
     * them. */
    double *step = (double *)ARENA_ALLOC(arena, (nbin ? nbin : 1) *
                                                sizeof(double));
    double *circ = (double *)ARENA_ALLOC(arena, (nbin ? nbin : 1) *
                                                sizeof(double));
    size_t nstep = 0;
    for (size_t i = 0; i + 1 < nbin; i++) {
        if (bin[i].piece != bin[i + 1].piece) continue;
        if (bin[i + 1].turn != bin[i].turn + 1) continue;
        double du = bin[i + 1].u_median - bin[i].u_median;
        double r = 0.5 * (bin[i].radius_median + bin[i + 1].radius_median);
        double owed = ATA_2PI * r;
        /* Sign convention comes from the data, not from us: the piece's own
         * majority direction decides which way is forward. */
        step[nstep] = du;
        circ[nstep] = owed;
        nstep++;
    }
    out->steps = nstep;
    if (nstep > 0) {
        double *mag = (double *)ARENA_ALLOC(arena, nstep * sizeof(double));
        double *cmag = (double *)ARENA_ALLOC(arena, nstep * sizeof(double));
        size_t positive = 0;
        for (size_t i = 0; i < nstep; i++) if (step[i] > 0.0) positive++;
        double dir = positive * 2 >= nstep ? 1.0 : -1.0;
        for (size_t i = 0; i < nstep; i++) {
            if (step[i] * dir < 0.0) out->steps_backward++;
            mag[i] = fabs(step[i]);
            cmag[i] = circ[i];
        }
        out->step_median = ata_median(mag, nstep);
        out->circumference_median = ata_median(cmag, nstep);
        out->advance_fraction = out->circumference_median > 1e-9
            ? out->step_median / out->circumference_median : 0.0;
    }
    out->bins = bin;
    out->nbins = nbin;
    return 0;
}

int AtlasTurnAdvance_edges(Arena_T arena,
                           const double *u, const double *phi,
                           const int32_t *edge_a, const int32_t *edge_b,
                           const double *target, size_t nedges,
                           size_t nsamples,
                           double spiral_a, double spiral_b,
                           int expect_zero,
                           AtlasTurnAdvanceEdgeStats *out)
{
    if (arena == NULL || out == NULL || u == NULL || phi == NULL) return -1;
    memset(out, 0, sizeof *out);
    if (nedges == 0) return 0;
    if (edge_a == NULL || edge_b == NULL) return -1;

    double *dphi = (double *)ARENA_ALLOC(arena, nedges * sizeof(double));
    double *du = (double *)ARENA_ALLOC(arena, nedges * sizeof(double));
    double *dexp = (double *)ARENA_ALLOC(arena, nedges * sizeof(double));
    double *ratio = (double *)ARENA_ALLOC(arena, nedges * sizeof(double));
    size_t nd = 0, nr = 0;

    for (size_t i = 0; i < nedges; i++) {
        int32_t a = edge_a[i], b = edge_b[i];
        if (a < 0 || b < 0 || (size_t)a >= nsamples || (size_t)b >= nsamples)
            continue;
        out->edges++;
        double dp = phi[b] - phi[a];
        double d = u[b] - u[a];
        /*
         * What this relation OWED.  An alignment row owes nothing by design,
         * so it is scored the other way round: against the spiral arc its
         * members span, which is the advance it is destroying if they are not
         * really one ruling.
         */
        double owed = target != NULL && !expect_zero
            ? target[i]
            : ata_spiral_u(spiral_a, spiral_b, phi[b]) -
              ata_spiral_u(spiral_a, spiral_b, phi[a]);
        dphi[nd] = fabs(dp) / ATA_2PI;
        du[nd] = d;
        dexp[nd] = owed;
        nd++;
        out->owed_total += fabs(owed);
        out->delivered_total += fabs(d);
        if (fabs(owed) > 1e-6) {
            out->scored++;
            ratio[nr++] = d / owed;
            if (d * owed < 0.0) out->reversed++;
            /* A quarter turn of real angle with under a quarter of the u it
             * implies: the relation is holding two wraps together. */
            if (fabs(dp) > 0.25 * ATA_2PI && fabs(d) < 0.25 * fabs(owed))
                out->collapsed++;
        }
    }
    out->dphi_turns_median = ata_median(dphi, nd);
    out->du_median = ata_median(du, nd);
    out->du_expected_median = ata_median(dexp, nd);
    out->ratio_median = ata_median(ratio, nr);
    return 0;
}

/* ------------------------------------------------------------------ selftest */

#define ATA_CHECK(cond, msg)                                                  \
    do {                                                                      \
        if (!(cond)) {                                                        \
            fprintf(stderr, "[atlas_turn_advance] FAIL: %s\n", (msg));        \
            failures++;                                                       \
        }                                                                     \
    } while (0)

int AtlasTurnAdvance_selftest(void)
{
    int failures = 0;
    Arena_T arena = Arena_new();

    /* A correct 4-turn spiral: u IS the pinned arclength. */
    const double sa = 200.0, sb = 12.0;   /* r(0) = 200, pitch 12 */
    enum { N = 400 };
    double u[N], phi[N], radius[N];
    int32_t piece[N];
    for (int i = 0; i < N; i++) {
        phi[i] = ATA_2PI * 4.0 * (double)i / (double)(N - 1);
        u[i] = ata_spiral_u(sa, sb, phi[i]);
        radius[i] = sa + sb * phi[i] / ATA_2PI;
        piece[i] = 0;
    }
    AtlasTurnAdvanceResult res;
    ATA_CHECK(AtlasTurnAdvance_bins(arena, u, phi, radius, piece,
                                    (size_t)N, 1, 1.0, &res) == 0,
              "bins on a correct spiral return 0");
    ATA_CHECK(res.nbins == 5, "four turns of phi fall in five turn bins");
    ATA_CHECK(res.steps_backward == 0, "a correct field never goes backward");
    ATA_CHECK(fabs(res.advance_fraction - 1.0) < 0.05,
              "a correct field advances one circumference per turn");

    /* The same geometry with u flattened: every sample gets one turn's worth,
     * which is exactly the observed failure mode. */
    double flat[N];
    for (int i = 0; i < N; i++)
        flat[i] = ata_spiral_u(sa, sb, fmod(phi[i], ATA_2PI));
    AtlasTurnAdvanceResult collapsed;
    ATA_CHECK(AtlasTurnAdvance_bins(arena, flat, phi, radius, piece,
                                    (size_t)N, 1, 1.0, &collapsed) == 0,
              "bins on a collapsed field return 0");
    ATA_CHECK(collapsed.advance_fraction < 0.25,
              "a collapsed field shows a small advance fraction");

    /* Edges: consecutive samples on the correct field deliver what they owe. */
    int32_t ea[N], eb[N];
    size_t ne = 0;
    for (int i = 0; i + 1 < N; i++) { ea[ne] = i; eb[ne] = i + 1; ne++; }
    AtlasTurnAdvanceEdgeStats es;
    ATA_CHECK(AtlasTurnAdvance_edges(arena, u, phi, ea, eb, NULL, ne,
                                     (size_t)N, sa, sb, 0, &es) == 0,
              "edges on a correct field return 0");
    ATA_CHECK(es.collapsed == 0, "no collapsed edges on a correct field");
    ATA_CHECK(es.reversed == 0, "no reversed edges on a correct field");
    ATA_CHECK(fabs(es.ratio_median - 1.0) < 1e-6,
              "each edge delivers exactly the arc it owes");

    /* An alignment row spanning a whole turn while asserting no u change:
     * one edge, and it must be caught. */
    int32_t fa[1] = {0}, fb[1] = {0};
    double same[2];
    /* sample 0 and the sample one full turn later */
    int later = (N - 1) / 4;
    fb[0] = (int32_t)later;
    same[0] = u[0];
    ATA_CHECK(fabs(phi[later] - phi[0] - ATA_2PI) < 0.1,
              "the probe pair really is one turn apart");
    {
        double held[N];
        for (int i = 0; i < N; i++) held[i] = u[i];
        held[later] = held[0];      /* the alignment row's instruction */
        AtlasTurnAdvanceEdgeStats align;
        ATA_CHECK(AtlasTurnAdvance_edges(arena, held, phi, fa, fb, NULL, 1,
                                         (size_t)N, sa, sb, 1, &align) == 0,
                  "alignment scoring returns 0");
        ATA_CHECK(align.collapsed == 1,
                  "an alignment row spanning a turn is flagged collapsed");
        ATA_CHECK(align.owed_total > 1000.0,
                  "and it is charged the full circumference it destroyed");
    }
    (void)same;

    /* Empty input must not crash. */
    AtlasTurnAdvanceEdgeStats empty;
    ATA_CHECK(AtlasTurnAdvance_edges(arena, u, phi, NULL, NULL, NULL, 0,
                                     (size_t)N, sa, sb, 0, &empty) == 0 &&
              empty.edges == 0, "zero edges is not an error");

    Arena_dispose(&arena);
    fprintf(stderr, "[atlas_turn_advance] selftest %s (%d failures)\n",
            failures == 0 ? "PASSED" : "FAILED", failures);
    return failures;
}
