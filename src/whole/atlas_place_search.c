#include "atlas_place_search.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void AtlasPlaceSearchOptions_default(AtlasPlaceSearchOptions *opts)
{
    if (opts == NULL) return;
    opts->cell = 4.0;
    opts->lambda_move = 8.0;
    opts->lambda_pair = 1.0;
    opts->rounds = 6;
}

/* Shift of a piece under a given choice; background (-1) never moves. */
static double aps_shift_of(const double *shift, size_t ncand,
                           const int32_t *choice, int32_t piece,
                           size_t npieces)
{
    if (piece < 0 || (size_t)piece >= npieces) return 0.0;
    return shift[(size_t)piece * ncand + (size_t)choice[piece]];
}

/* A pair is satisfied when the two ends land within half the circumference
 * they owe of that circumference -- the same test the radial-order stage uses
 * to call a pair coincident, so the two numbers are directly comparable. */
static int aps_pair_ok(const AtlasPlacePair *p, double sa, double sb)
{
    double want = fabs(p->target);
    if (!(want > 0.0)) return 1;
    double du = p->du_base + sb - sa;
    return fabs(du - p->target) < 0.5 * want;
}

/*
 * A cell key packs the two raster indices so that a u translation is a plain
 * addition.  v occupies the high bits and never moves, u the low bits, and the
 * bias keeps u non-negative so the packed order matches (v, u) lexicographic
 * order -- which is what lets the intersection below be a linear merge over
 * two sorted arrays instead of a hash lookup per cell.
 */
#define APS_UBITS 26
#define APS_UBIAS (((int64_t)1 << (APS_UBITS - 1)))
#define APS_UMASK ((((int64_t)1) << APS_UBITS) - 1)

static int64_t aps_key(int64_t iv, int64_t iu)
{
    return (iv << APS_UBITS) | ((iu + APS_UBIAS) & APS_UMASK);
}

static int aps_cmp_i64(const void *a, const void *b)
{
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

static size_t aps_dedup(int64_t *v, size_t n)
{
    if (n == 0) return 0;
    qsort(v, n, sizeof(int64_t), aps_cmp_i64);
    size_t m = 1;
    for (size_t i = 1; i < n; i++)
        if (v[i] != v[m - 1]) v[m++] = v[i];
    return m;
}

/* |A shifted by `off` in u| intersect |B|, both sorted and deduped.  Adding
 * `off` preserves order because the u field never overflows into v at the
 * magnitudes involved (26 bits of cells is ~67M vox at cell 1). */
static size_t aps_overlap(const int64_t *a, size_t na, int64_t off,
                          const int64_t *b, size_t nb)
{
    size_t i = 0, j = 0, hits = 0;
    while (i < na && j < nb) {
        int64_t x = a[i] + off;
        if (x < b[j]) i++;
        else if (x > b[j]) j++;
        else { hits++; i++; j++; }
    }
    return hits;
}

int AtlasPlaceSearch_solve(Arena_T arena,
                           const double *u, const double *v, size_t nv,
                           const int32_t *vertex_piece,
                           size_t npieces, size_t ncand,
                           const double *shift, const double *cost,
                           const AtlasPlacePair *pairs, size_t npairs,
                           const AtlasPlaceSearchOptions *opts,
                           int32_t *out_choice,
                           AtlasPlaceSearchStats *stats)
{
    if (arena == NULL || u == NULL || v == NULL || vertex_piece == NULL ||
        shift == NULL || out_choice == NULL || stats == NULL)
        return -1;
    if (npairs > 0 && pairs == NULL) return -1;
    if (npieces == 0 || ncand == 0) return -1;
    AtlasPlaceSearchOptions o;
    if (opts != NULL) o = *opts; else AtlasPlaceSearchOptions_default(&o);
    if (!(o.cell > 0.0)) return -1;

    memset(stats, 0, sizeof *stats);
    stats->pieces = npieces;
    double inv = 1.0 / o.cell;

    /* --- rasterize once: background, and each piece at zero shift --------- */
    size_t *count = (size_t *)ARENA_CALLOC(arena, npieces + 1, sizeof(size_t));
    for (size_t i = 0; i < nv; i++) {
        int32_t p = vertex_piece[i];
        if (p < 0) count[npieces]++;
        else if ((size_t)p < npieces) count[p]++;
    }
    int64_t **cells = (int64_t **)ARENA_ALLOC(
        arena, (npieces + 1) * sizeof(int64_t *));
    size_t *ncells = (size_t *)ARENA_CALLOC(arena, npieces + 1,
                                            sizeof(size_t));
    for (size_t p = 0; p <= npieces; p++)
        cells[p] = (int64_t *)ARENA_ALLOC(
            arena, (count[p] ? count[p] : 1) * sizeof(int64_t));
    for (size_t i = 0; i < nv; i++) {
        int32_t p = vertex_piece[i];
        if (p >= 0 && (size_t)p >= npieces) continue;
        size_t slot = p < 0 ? npieces : (size_t)p;
        int64_t iu = (int64_t)floor(u[i] * inv);
        int64_t iv = (int64_t)floor(v[i] * inv);
        cells[slot][ncells[slot]++] = aps_key(iv, iu);
    }
    for (size_t p = 0; p <= npieces; p++)
        ncells[p] = aps_dedup(cells[p], ncells[p]);
    stats->static_cells = ncells[npieces];
    for (size_t p = 0; p < npieces; p++) stats->moving_cells += ncells[p];

    /* Candidate j of piece p is a whole-cell offset along u. */
    int64_t *off = (int64_t *)ARENA_ALLOC(
        arena, npieces * ncand * sizeof(int64_t));
    for (size_t p = 0; p < npieces; p++)
        for (size_t j = 0; j < ncand; j++)
            off[p * ncand + j] =
                (int64_t)llround(shift[p * ncand + j] * inv);

    /* --- radial pairs, indexed by the piece each one can move -------------- */
    size_t *pdeg = (size_t *)ARENA_CALLOC(arena, npieces + 1, sizeof(size_t));
    for (size_t i = 0; i < npairs; i++) {
        const AtlasPlacePair *pr = &pairs[i];
        int a = pr->piece_a >= 0 && (size_t)pr->piece_a < npieces;
        int b = pr->piece_b >= 0 && (size_t)pr->piece_b < npieces;
        stats->pairs++;
        if (a && b && pr->piece_a == pr->piece_b) {
            /* Both ends ride together: no choice here can change it. */
            stats->pairs_intra++;
            continue;
        }
        if (a) pdeg[(size_t)pr->piece_a + 1]++;
        if (b) pdeg[(size_t)pr->piece_b + 1]++;
    }
    for (size_t p = 0; p < npieces; p++) pdeg[p + 1] += pdeg[p];
    size_t npadj = pdeg[npieces];
    int32_t *padj = (int32_t *)ARENA_ALLOC(
        arena, (npadj ? npadj : 1) * sizeof(int32_t));
    {
        size_t *cur = (size_t *)ARENA_ALLOC(arena, npieces * sizeof(size_t));
        memcpy(cur, pdeg, npieces * sizeof(size_t));
        for (size_t i = 0; i < npairs; i++) {
            const AtlasPlacePair *pr = &pairs[i];
            int a = pr->piece_a >= 0 && (size_t)pr->piece_a < npieces;
            int b = pr->piece_b >= 0 && (size_t)pr->piece_b < npieces;
            if (a && b && pr->piece_a == pr->piece_b) continue;
            if (a) padj[cur[(size_t)pr->piece_a]++] = (int32_t)i;
            if (b) padj[cur[(size_t)pr->piece_b]++] = (int32_t)i;
        }
    }

    /*
     * --- coordinate descent ----------------------------------------------
     * Start every piece at its INCOMING placement -- the candidate with the
     * smallest |shift| (the driver puts the zero-shift wind somewhere in the
     * middle of the table, not at index 0).  Starting at index 0 would begin
     * the descent, and measure every "before" statistic, at a fictitious
     * all-pieces-at-minus-span state.
     */
    for (size_t p = 0; p < npieces; p++) {
        size_t best = 0;
        for (size_t j = 1; j < ncand; j++)
            if (fabs(shift[p * ncand + j]) <
                fabs(shift[p * ncand + best])) best = j;
        out_choice[p] = (int32_t)best;
    }
    for (size_t i = 0; i < npairs; i++) {
        const AtlasPlacePair *pr = &pairs[i];
        double sa = aps_shift_of(shift, ncand, out_choice, pr->piece_a,
                                 npieces);
        double sb = aps_shift_of(shift, ncand, out_choice, pr->piece_b,
                                 npieces);
        if (aps_pair_ok(pr, sa, sb)) stats->pairs_satisfied_before++;
    }

    size_t before = 0;
    for (size_t p = 0; p < npieces; p++) {
        int64_t op = off[p * ncand + (size_t)out_choice[p]];
        before += aps_overlap(cells[p], ncells[p], op,
                              cells[npieces], ncells[npieces]);
        for (size_t q = p + 1; q < npieces; q++) {
            int64_t oq = off[q * ncand + (size_t)out_choice[q]];
            /* shift q's cells into p's frame by the difference */
            before += aps_overlap(cells[p], ncells[p], op - oq,
                                  cells[q], ncells[q]);
        }
    }
    stats->collisions_before = before;

    int converged = 0;
    int r = 0;
    for (r = 0; r < o.rounds; r++) {
        int changed = 0;
        for (size_t p = 0; p < npieces; p++) {
            double best_score = 0.0;
            int32_t best_j = out_choice[p];
            int have = 0;
            for (size_t j = 0; j < ncand; j++) {
                int64_t op = off[p * ncand + j];
                size_t hits = aps_overlap(cells[p], ncells[p], op,
                                          cells[npieces], ncells[npieces]);
                for (size_t q = 0; q < npieces; q++) {
                    if (q == p) continue;
                    int64_t oq = off[q * ncand + (size_t)out_choice[q]];
                    hits += aps_overlap(cells[p], ncells[p], op - oq,
                                        cells[q], ncells[q]);
                }
                double score = (double)hits;
                if (cost != NULL)
                    score += o.lambda_move * cost[p * ncand + j];
                /* Radial disagreement, counted only over the pairs this piece
                 * can actually move; the rest are constant in this decision. */
                if (o.lambda_pair > 0.0) {
                    size_t bad = 0;
                    for (size_t e = pdeg[p]; e < pdeg[p + 1]; e++) {
                        const AtlasPlacePair *pr = &pairs[padj[e]];
                        double sa = pr->piece_a == (int32_t)p
                            ? shift[p * ncand + j]
                            : aps_shift_of(shift, ncand, out_choice,
                                           pr->piece_a, npieces);
                        double sb = pr->piece_b == (int32_t)p
                            ? shift[p * ncand + j]
                            : aps_shift_of(shift, ncand, out_choice,
                                           pr->piece_b, npieces);
                        if (!aps_pair_ok(pr, sa, sb)) bad++;
                    }
                    score += o.lambda_pair * (double)bad;
                }
                stats->candidates_evaluated++;
                if (!have || score < best_score) {
                    best_score = score;
                    best_j = (int32_t)j;
                    have = 1;
                }
            }
            if (best_j != out_choice[p]) { out_choice[p] = best_j; changed = 1; }
        }
        if (!changed) { converged = 1; break; }
    }
    stats->rounds_run = r + (converged ? 1 : 0);
    stats->converged = converged;

    size_t after = 0, after_static = 0;
    for (size_t p = 0; p < npieces; p++) {
        int64_t op = off[p * ncand + (size_t)out_choice[p]];
        size_t s = aps_overlap(cells[p], ncells[p], op,
                               cells[npieces], ncells[npieces]);
        after += s;
        after_static += s;
        for (size_t q = p + 1; q < npieces; q++) {
            int64_t oq = off[q * ncand + (size_t)out_choice[q]];
            after += aps_overlap(cells[p], ncells[p], op - oq,
                                 cells[q], ncells[q]);
        }
        if (out_choice[p] != 0) stats->moved_pieces++;
    }
    stats->collisions_after = after;
    stats->collisions_static_after = after_static;
    for (size_t i = 0; i < npairs; i++) {
        const AtlasPlacePair *pr = &pairs[i];
        double sa = aps_shift_of(shift, ncand, out_choice, pr->piece_a,
                                 npieces);
        double sb = aps_shift_of(shift, ncand, out_choice, pr->piece_b,
                                 npieces);
        if (aps_pair_ok(pr, sa, sb)) stats->pairs_satisfied_after++;
    }
    return 0;
}

/* ------------------------------------------------------------------ selftest */

#define APS_CHECK(cond, msg)                                                  \
    do {                                                                      \
        if (!(cond)) {                                                        \
            fprintf(stderr, "[atlas_place_search] FAIL: %s\n", (msg));        \
            failures++;                                                       \
        }                                                                     \
    } while (0)

int AtlasPlaceSearch_selftest(void)
{
    int failures = 0;
    Arena_T arena = Arena_new();

    /*
     * Three unit patches stacked on the same u, plus a background patch one
     * step to the right.  The only placement with no collision at all puts the
     * three at steps 0, 2 and 3 -- step 1 is taken by the background -- and
     * the travel penalty should keep them from wandering further.
     */
    enum { PER = 25, NP = 3, NV = (NP + 1) * PER, NC = 5 };
    double u[NV], v[NV];
    int32_t piece[NV];
    const double step = 40.0;
    for (int p = 0; p < NP + 1; p++)
        for (int i = 0; i < PER; i++) {
            int idx = p * PER + i;
            u[idx] = (double)(i % 5) * 4.0 + (p == NP ? step : 0.0);
            v[idx] = (double)(i / 5) * 4.0;
            piece[idx] = p == NP ? -1 : p;
        }

    double shift[NP * NC], cost[NP * NC];
    for (int p = 0; p < NP; p++)
        for (int j = 0; j < NC; j++) {
            shift[p * NC + j] = step * (double)j;
            cost[p * NC + j] = (double)j;
        }

    AtlasPlaceSearchOptions opts;
    AtlasPlaceSearchOptions_default(&opts);
    opts.cell = 4.0;
    opts.lambda_move = 1.0;
    opts.lambda_pair = 0.0;
    int32_t choice[NP];
    AtlasPlaceSearchStats st;
    int rc = AtlasPlaceSearch_solve(arena, u, v, (size_t)NV, piece,
                                    (size_t)NP, (size_t)NC, shift, cost,
                                    NULL, 0, &opts, choice, &st);
    APS_CHECK(rc == 0, "solve returns 0");
    APS_CHECK(st.collisions_before > 0, "the stacked input really does collide");
    APS_CHECK(st.collisions_after == 0, "the search finds a collision-free placement");
    APS_CHECK(st.converged, "and it converges rather than running out of rounds");
    int seen[NC];
    memset(seen, 0, sizeof seen);
    for (int p = 0; p < NP; p++) {
        APS_CHECK(choice[p] >= 0 && choice[p] < NC, "choice in range");
        APS_CHECK(!seen[choice[p]], "no two pieces land on the same step");
        seen[choice[p]] = 1;
    }
    APS_CHECK(!seen[1], "the step occupied by background is avoided");

    /* With no travel penalty at all a piece may wander, but it must still
     * never choose a colliding step when a free one exists. */
    opts.lambda_move = 0.0;
    AtlasPlaceSearchStats st2;
    rc = AtlasPlaceSearch_solve(arena, u, v, (size_t)NV, piece, (size_t)NP,
                                (size_t)NC, shift, cost, NULL, 0, &opts,
                                choice, &st2);
    APS_CHECK(rc == 0 && st2.collisions_after == 0,
              "still collision-free with no travel penalty");

    /*
     * Collision-freedom alone is degenerate: several permutations tile just as
     * well.  A radial pair breaks the tie, and it must be able to override a
     * placement the collision term is indifferent to.  Pieces 0 and 2 are
     * asserted to sit exactly two steps apart.
     */
    AtlasPlacePair pr[1];
    pr[0].piece_a = 0;
    pr[0].piece_b = 2;
    pr[0].du_base = 0.0;
    pr[0].target = 2.0 * step;
    opts.lambda_move = 1.0;
    opts.lambda_pair = 100.0;
    AtlasPlaceSearchStats st4;
    rc = AtlasPlaceSearch_solve(arena, u, v, (size_t)NV, piece, (size_t)NP,
                                (size_t)NC, shift, cost, pr, 1, &opts,
                                choice, &st4);
    APS_CHECK(rc == 0, "solve with a radial pair returns 0");
    APS_CHECK(st4.pairs == 1 && st4.pairs_intra == 0, "the pair is counted");
    APS_CHECK(st4.pairs_satisfied_before == 0,
              "the stacked input does not satisfy it");
    APS_CHECK(st4.pairs_satisfied_after == 1,
              "the search moves the two pieces two steps apart");
    APS_CHECK(st4.collisions_after == 0,
              "and it stays collision-free while doing so");

    /* A pair with both ends in one piece is unreachable and must be reported
     * as such rather than counted against any candidate. */
    pr[0].piece_b = 0;
    AtlasPlaceSearchStats st5;
    rc = AtlasPlaceSearch_solve(arena, u, v, (size_t)NV, piece, (size_t)NP,
                                (size_t)NC, shift, cost, pr, 1, &opts,
                                choice, &st5);
    APS_CHECK(rc == 0 && st5.pairs_intra == 1,
              "an intra-piece pair is reported as unreachable");

    /* A single piece with nothing to hit should not move at all. */
    for (int i = 0; i < NV; i++) piece[i] = i < PER ? 0 : -1;
    for (int i = PER; i < NV; i++) u[i] = 10000.0;   /* background far away */
    opts.lambda_pair = 0.0;
    AtlasPlaceSearchStats st3;
    rc = AtlasPlaceSearch_solve(arena, u, v, (size_t)NV, piece, 1,
                                (size_t)NC, shift, cost, NULL, 0, &opts,
                                choice, &st3);
    APS_CHECK(rc == 0 && choice[0] == 0,
              "a piece with nothing to avoid stays where it is");

    Arena_dispose(&arena);
    fprintf(stderr, "[atlas_place_search] selftest %s (%d failures)\n",
            failures == 0 ? "PASSED" : "FAILED", failures);
    return failures;
}
