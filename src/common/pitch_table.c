#include "pitch_table.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PITCH_TABLE_MAX_KNOTS 4096
#define PITCH_TABLE_MIN_PITCH 0.25

struct PitchTable {
    int     nknots;
    double *r;       /* ascending radii */
    double *pitch;   /* pitch at each knot, > 0 */
    double *turns;   /* cumulative integral of dr/pitch from r[0] */
};

/* Analytic winding turns across one linear-pitch segment [r0,r1] evaluated
 * from r0 to r (r within the segment).  pitch(x) = p0 + s*(x - r0). */
static double segment_turns(double r0, double p0, double r1, double p1,
                            double r) {
    double dr = r1 - r0;
    if (!(dr > 0.0)) return 0.0;
    double s = (p1 - p0) / dr;
    double x = r - r0;
    if (fabs(s) < 1e-9) return x / p0;
    return log((p0 + s * x) / p0) / s;
}

static double table_turns_at(PitchTable t, double r) {
    if (r <= t->r[0])
        return (r - t->r[0]) / t->pitch[0];   /* constant extension */
    int hi = t->nknots - 1;
    if (r >= t->r[hi])
        return t->turns[hi] + (r - t->r[hi]) / t->pitch[hi];
    int k = 0;
    while (k + 1 < t->nknots && t->r[k + 1] < r) k++;
    return t->turns[k] + segment_turns(t->r[k], t->pitch[k],
                                       t->r[k + 1], t->pitch[k + 1], r);
}

double PitchTable_pitch_at(PitchTable t, double r) {
    if (t == NULL || t->nknots < 1) return 0.0;
    if (r <= t->r[0]) return t->pitch[0];
    int hi = t->nknots - 1;
    if (r >= t->r[hi]) return t->pitch[hi];
    int k = 0;
    while (k + 1 < t->nknots && t->r[k + 1] < r) k++;
    double dr = t->r[k + 1] - t->r[k];
    double w = dr > 0.0 ? (r - t->r[k]) / dr : 0.0;
    return t->pitch[k] + w * (t->pitch[k + 1] - t->pitch[k]);
}

double PitchTable_dturns(PitchTable t, double r_a, double r_b,
                         double scalar_fallback) {
    if (t == NULL || t->nknots < 1)
        return (r_b - r_a) / scalar_fallback;
    return table_turns_at(t, r_b) - table_turns_at(t, r_a);
}

int PitchTable_knots(PitchTable t) { return t ? t->nknots : 0; }

void PitchTable_free(PitchTable t) {
    if (t == NULL) return;
    free(t->r); free(t->pitch); free(t->turns); free(t);
}

static int knot_cmp(const void *aa, const void *bb) {
    double a = *(const double *)aa, b = *(const double *)bb;
    return a < b ? -1 : a > b ? 1 : 0;
}

int PitchTable_load(const char *path, PitchTable *out) {
    if (out == NULL) return -1;
    *out = NULL;
    if (path == NULL || path[0] == '\0') return -1;
    FILE *f = fopen(path, "rb");
    if (f == NULL) return -1;
    /* rows collected as (r, pitch) pairs for a sort-by-r before use */
    double (*rows)[2] = malloc(PITCH_TABLE_MAX_KNOTS * sizeof *rows);
    if (rows == NULL) { fclose(f); return -1; }
    int n = 0, bad = 0;
    char line[512];
    while (fgets(line, (int)sizeof line, f) != NULL) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\r' || *p == '\n' || *p == '\0') continue;
        double r = 0.0, pitch = 0.0;
        if (sscanf(p, "%lf,%lf", &r, &pitch) != 2) {
            /* tolerate one header line ("r,pitch,...") */
            if (n == 0 && strchr(p, ',') != NULL) continue;
            bad = 1; break;
        }
        if (!isfinite(r) || !isfinite(pitch) ||
            pitch < PITCH_TABLE_MIN_PITCH) { bad = 1; break; }
        if (n >= PITCH_TABLE_MAX_KNOTS) { bad = 1; break; }
        rows[n][0] = r; rows[n][1] = pitch; n++;
    }
    fclose(f);
    if (bad || n < 1) { free(rows); return -1; }
    qsort(rows, (size_t)n, sizeof *rows, knot_cmp);
    for (int i = 1; i < n; i++)
        if (!(rows[i][0] > rows[i - 1][0])) { free(rows); return -1; }
    PitchTable t = calloc(1, sizeof *t);
    if (t == NULL) { free(rows); return -1; }
    t->nknots = n;
    t->r = malloc((size_t)n * sizeof *t->r);
    t->pitch = malloc((size_t)n * sizeof *t->pitch);
    t->turns = malloc((size_t)n * sizeof *t->turns);
    if (t->r == NULL || t->pitch == NULL || t->turns == NULL) {
        PitchTable_free(t); free(rows); return -1;
    }
    for (int i = 0; i < n; i++) { t->r[i] = rows[i][0]; t->pitch[i] = rows[i][1]; }
    free(rows);
    t->turns[0] = 0.0;
    for (int i = 1; i < n; i++)
        t->turns[i] = t->turns[i - 1] +
            segment_turns(t->r[i - 1], t->pitch[i - 1],
                          t->r[i], t->pitch[i], t->r[i]);
    *out = t;
    return 0;
}

PitchTable PitchTable_from_env(void) {
    static PitchTable cached = NULL;
    static int attempted = 0;
    if (!attempted) {
        attempted = 1;
        const char *path = getenv("VES_WRAP_PITCH_TABLE");
        if (path != NULL && path[0] != '\0') {
            if (PitchTable_load(path, &cached) != 0)
                fprintf(stderr,
                        "pitch_table: FAILED to load VES_WRAP_PITCH_TABLE=%s; "
                        "gates fall back to scalar pitch\n", path);
            else
                fprintf(stderr,
                        "pitch_table: loaded %d knots from %s "
                        "(r %.1f..%.1f, pitch %.2f..%.2f)\n",
                        cached->nknots, path,
                        cached->r[0], cached->r[cached->nknots - 1],
                        cached->pitch[0], cached->pitch[cached->nknots - 1]);
        }
    }
    return cached;
}

int PitchTable_selftest(void) {
    int fails = 0;
#define PT_CK(c, m) do { if (!(c)) { fprintf(stderr, "  FAIL: %s\n", (m)); \
                          fails++; } else fprintf(stderr, "  ok: %s\n", (m)); \
                        } while (0)
    /* NULL table == scalar, bit-exact. */
    {
        double a = PitchTable_dturns(NULL, 37.25, 129.5, 9.5);
        double b = (129.5 - 37.25) / 9.5;
        PT_CK(a == b, "NULL table matches scalar expression bit-exactly");
    }
    /* Constant table equals scalar to double rounding. */
    {
        char path[512];
        const char *tmp = getenv("TEMP");
        snprintf(path, sizeof path, "%s\\pitch_table_selftest_const.csv",
                 tmp ? tmp : ".");
        FILE *f = fopen(path, "wb");
        if (f) {
            fputs("# selftest\nr,pitch,count\n50,20,100\n450,20,100\n", f);
            fclose(f);
        }
        PitchTable t = NULL;
        int rc = PitchTable_load(path, &t);
        PT_CK(rc == 0 && PitchTable_knots(t) == 2, "constant table loads");
        if (t) {
            double a = PitchTable_dturns(t, 60.0, 400.0, 9.5);
            PT_CK(fabs(a - (400.0 - 60.0) / 20.0) < 1e-9,
                  "constant table integrates as dr/pitch");
            PT_CK(fabs(PitchTable_pitch_at(t, 10.0) - 20.0) < 1e-12 &&
                  fabs(PitchTable_pitch_at(t, 900.0) - 20.0) < 1e-12,
                  "table clamps to end pitches outside its range");
        }
        PitchTable_free(t);
        remove(path);
    }
    /* Two-segment log integral matches numeric quadrature; antisymmetry;
     * one wrap of radius spans ~one turn. */
    {
        char path[512];
        const char *tmp = getenv("TEMP");
        snprintf(path, sizeof path, "%s\\pitch_table_selftest_ramp.csv",
                 tmp ? tmp : ".");
        FILE *f = fopen(path, "wb");
        if (f) {
            fputs("r,pitch\n100,10\n200,30\n300,20\n", f);
            fclose(f);
        }
        PitchTable t = NULL;
        int rc = PitchTable_load(path, &t);
        PT_CK(rc == 0 && PitchTable_knots(t) == 3, "ramp table loads");
        if (t) {
            double analytic = PitchTable_dturns(t, 120.0, 280.0, 9.5);
            double numeric = 0.0;
            const int steps = 200000;
            for (int i = 0; i < steps; i++) {
                double r = 120.0 + (280.0 - 120.0) * (i + 0.5) / steps;
                numeric += (280.0 - 120.0) / steps /
                           PitchTable_pitch_at(t, r);
            }
            PT_CK(fabs(analytic - numeric) < 1e-6,
                  "analytic segment integral matches quadrature");
            PT_CK(fabs(PitchTable_dturns(t, 280.0, 120.0, 9.5) + analytic)
                      < 1e-12,
                  "dturns is antisymmetric");
            double p = PitchTable_pitch_at(t, 150.0);
            double one_wrap = PitchTable_dturns(t, 150.0, 150.0 + p, 9.5);
            PT_CK(one_wrap > 0.85 && one_wrap < 1.15,
                  "one local wrap of radius spans about one turn");
        }
        PitchTable_free(t);
        remove(path);
    }
    /* Loader rejects garbage. */
    {
        char path[512];
        const char *tmp = getenv("TEMP");
        snprintf(path, sizeof path, "%s\\pitch_table_selftest_bad.csv",
                 tmp ? tmp : ".");
        FILE *f = fopen(path, "wb");
        if (f) { fputs("r,pitch\n100,10\n100,12\n", f); fclose(f); }
        PitchTable t = NULL;
        PT_CK(PitchTable_load(path, &t) != 0 && t == NULL,
              "duplicate radii are rejected");
        remove(path);
        PT_CK(PitchTable_load("definitely_missing_file.csv", &t) != 0,
              "missing file is rejected");
    }
#undef PT_CK
    fprintf(stderr, "[selftest] pitch_table: %s (%d failures)\n",
            fails == 0 ? "ALL PASS" : "FAILURES", fails);
    return fails;
}
