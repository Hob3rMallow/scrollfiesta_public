/* ============================================================================
 * quadribbon_solidify.c -- see quadribbon_solidify.h.
 *
 * One solid quad ribbon per atlas column run: observations exact, holes
 * filled by quad_strip's cylindrical residual harmonic solve (TAUCS), fills
 * capped by lattice distance, every vertex labeled by provenance.  Runs are
 * separate wraps and are never bridged.
 * ==========================================================================*/
#include "quadribbon_solidify.h"
#include "../common/zarr_u8.h"
#include "quad_strip.h"

#include "../common/json_read.h"
#include "../common/mesh_bin.h"
#include "../common/pipeline_constants.h"
#include "ridge_track.h"
#include "../common/ves_platform.h"
#include "../common/ves_png.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define QS_PATH 2048
#define QS_TWO_PI 6.283185307179586476925286766559
#define QS_NONE 255u

/* ---- small file helpers -------------------------------------------------- */

static int qs_replace_file(const char *temporary, const char *destination)
{
#ifdef _WIN32
    return MoveFileExA(temporary, destination,
                       MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)
         ? 0 : -1;
#else
    return rename(temporary, destination);
#endif
}

/* Atomic raw write: <path>.tmp then rename, so a crash never leaves a
 * half-written sidecar beside a complete VMESH. */
static int qs_write_raw(const char *path, const void *data, size_t nbytes)
{
    char tmp[QS_PATH];
    FILE *f = NULL;
    if (snprintf(tmp, sizeof tmp, "%s.tmp", path) >= (int)sizeof tmp) return -1;
    f = fopen(tmp, "wb");
    if (f == NULL) return -1;
    if (nbytes > 0 && fwrite(data, 1, nbytes, f) != nbytes) {
        fclose(f);
        remove(tmp);
        return -1;
    }
    if (fclose(f) != 0) { remove(tmp); return -1; }
    if (qs_replace_file(tmp, path) != 0) { remove(tmp); return -1; }
    return 0;
}

static long long qs_file_size(const char *path)
{
    FILE *f = fopen(path, "rb");
    long long n = -1;
    if (f == NULL) return -1;
#ifdef _WIN32
    if (_fseeki64(f, 0, SEEK_END) == 0) n = _ftelli64(f);
#else
    if (fseek(f, 0, SEEK_END) == 0) n = (long long)ftell(f);
#endif
    fclose(f);
    return n;
}

static int qs_file_has_text(const char *path, const char *needle)
{
    FILE *f = fopen(path, "rb");
    char line[4096];
    if (f == NULL || needle == NULL) {
        if (f != NULL) fclose(f);
        return 0;
    }
    while (fgets(line, sizeof line, f) != NULL) {
        if (strstr(line, needle) != NULL) {
            fclose(f);
            return 1;
        }
    }
    fclose(f);
    return 0;
}

/* Read a raw sidecar that must hold EXACTLY count*esize bytes; NULL otherwise
 * (missing or wrong length -- a stale sidecar fails closed). */
static void *qs_read_exact(Arena_T arena, const char *path, size_t count,
                           size_t esize)
{
    long long sz = qs_file_size(path);
    size_t want = count * esize;
    FILE *f = NULL;
    void *buf = NULL;
    if (sz < 0 || (unsigned long long)sz != (unsigned long long)want) return NULL;
    f = fopen(path, "rb");
    if (f == NULL) return NULL;
    buf = ARENA_ALLOC(arena, want > 0 ? want : 1);
    if (want > 0 && fread(buf, 1, want, f) != want) { fclose(f); return NULL; }
    fclose(f);
    return buf;
}

/* Header peek: vertex count of a VESMESH1 file without reading the payload. */
static uint64_t qs_vmesh_nv(const char *path)
{
    unsigned char hdr[64];
    FILE *f = fopen(path, "rb");
    uint64_t nv = 0;
    if (f == NULL) return 0;
    if (fread(hdr, 1, sizeof hdr, f) != sizeof hdr ||
        memcmp(hdr, "VESMESH1", 8) != 0) {
        fclose(f);
        return 0;
    }
    fclose(f);
    memcpy(&nv, hdr + 24, sizeof nv);
    return nv;
}

static void qs_stem_of(const char *vmesh, char *stem, size_t cap)
{
    size_t n = strlen(vmesh);
    if (n > 6 && strcmp(vmesh + n - 6, ".vmesh") == 0) n -= 6;
    if (n >= cap) n = cap - 1;
    memcpy(stem, vmesh, n);
    stem[n] = '\0';
}

static void qs_dir_of(const char *stem, char *dir, size_t cap)
{
    const char *a = strrchr(stem, '/');
    const char *b = strrchr(stem, '\\');
    const char *slash = a != NULL && (b == NULL || a > b) ? a : b;
    size_t n = slash != NULL ? (size_t)(slash - stem) : 0;
    if (n >= cap) n = cap - 1;
    memcpy(dir, stem, n);
    dir[n] = '\0';
    if (n == 0) snprintf(dir, cap, ".");
}

static int qs_cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

static double qs_quantile(const double *sorted, size_t n, double q)
{
    size_t i = 0;
    if (n == 0) return 0.0;
    i = (size_t)(q * (double)(n - 1) + 0.5);
    if (i >= n) i = n - 1;
    return sorted[i];
}

/* Lattice pitch straight from the deliverable: the median gap between
 * consecutive DISTINCT values of one UV component (ribbon_verdict's rule). */
static double qs_detect_pitch(const float *uv, size_t nv, int comp)
{
    double *vals = NULL, *gaps = NULL, med = 0.0;
    size_t i = 0, nun = 0, ng = 0;
    if (nv < 2) return 1.0;
    vals = (double *)malloc(nv * sizeof *vals);
    if (vals == NULL) return 1.0;
    for (i = 0; i < nv; i++) vals[i] = (double)uv[i * 2 + comp];
    qsort(vals, nv, sizeof *vals, qs_cmp_double);
    for (i = 1, nun = 1; i < nv; i++)
        if (vals[i] > vals[nun - 1] + 1e-6) vals[nun++] = vals[i];
    if (nun < 2) { free(vals); return 1.0; }
    gaps = (double *)malloc((nun - 1) * sizeof *gaps);
    if (gaps == NULL) { free(vals); return 1.0; }
    for (i = 1; i < nun; i++) gaps[ng++] = vals[i] - vals[i - 1];
    qsort(gaps, ng, sizeof *gaps, qs_cmp_double);
    med = gaps[ng / 2];
    free(gaps);
    free(vals);
    return med > 1e-6 ? med : 1.0;
}

/* ---- input ---------------------------------------------------------------- */

typedef struct QsInput {
    MeshBinData mesh;
    const uint8_t *support;    /* [nv] 255 fixed / 0 generated */
    const uint8_t *prov;       /* [nv] QS_PROV_* (derived when absent) */
    const float   *phase;      /* [nv] lifted phase, radians */
    const int32_t *material;   /* [nv] or NULL */
    const int32_t *claimant_chart; /* [nv] stable semantic claimant, or NULL */
    const int32_t *reconstruction; /* [nv] or NULL */
    const int32_t *lane;       /* [nv] or reconstruction fallback */
    int prov_present;          /* provenance came from its exact sidecar */
    int stats_n_slices;        /* from <stem>_stats.json, 0 = unknown */
    double bbox_lo[3], bbox_hi[3];
} QsInput;

static int qs_load_input(Arena_T arena, const char *in_vmesh, QsInput *in,
                         QsReport *rep)
{
    char stem[QS_PATH], path[QS_PATH];
    size_t nv = 0;
    memset(in, 0, sizeof *in);
    if (MeshBin_read_arena(arena, in_vmesh, &in->mesh) != 0) {
        fprintf(stderr, "[solidify] cannot read %s\n", in_vmesh);
        return -1;
    }
    if (in->mesh.uv == NULL) {
        fprintf(stderr, "[solidify] %s carries no UV; the lattice IS the "
                        "parameterization\n", in_vmesh);
        return -1;
    }
    nv = in->mesh.nv;
    qs_stem_of(in_vmesh, stem, sizeof stem);
    snprintf(path, sizeof path, "%s_support.u8", stem);
    in->support = (const uint8_t *)qs_read_exact(arena, path, nv, 1);
    if (in->support == NULL) {
        fprintf(stderr, "[solidify] required sidecar missing or mis-sized: "
                        "%s (need %zu bytes)\n", path, nv);
        return -1;
    }
    snprintf(path, sizeof path, "%s_phase.f32", stem);
    in->phase = (const float *)qs_read_exact(arena, path, nv, sizeof(float));
    if (in->phase == NULL) {
        fprintf(stderr, "[solidify] required sidecar missing or mis-sized: "
                        "%s (the winding-aware fill needs the lifted phase)\n",
                path);
        return -1;
    }
    snprintf(path, sizeof path, "%s_provenance.u8", stem);
    in->prov = (const uint8_t *)qs_read_exact(arena, path, nv, 1);
    in->prov_present = in->prov != NULL;
    if (in->prov == NULL) {
        uint8_t *p = (uint8_t *)ARENA_ALLOC(arena, nv > 0 ? nv : 1);
        for (size_t i = 0; i < nv; i++)
            p[i] = in->support[i] != 0 ? (uint8_t)QS_PROV_OBSERVED_FIXED
                                       : (uint8_t)QS_PROV_OBSERVED_GENERATED;
        in->prov = p;
    }
    snprintf(path, sizeof path, "%s_material_identity.i32", stem);
    in->material = (const int32_t *)qs_read_exact(arena, path, nv,
                                                  sizeof(int32_t));
    snprintf(path, sizeof path, "%s_claimant_chart.i32", stem);
    in->claimant_chart = (const int32_t *)qs_read_exact(
        arena, path, nv, sizeof(int32_t));
    snprintf(path, sizeof path, "%s_reconstruction_component.i32", stem);
    in->reconstruction = (const int32_t *)qs_read_exact(
        arena, path, nv, sizeof(int32_t));
    snprintf(path, sizeof path, "%s_lane.i32", stem);
    in->lane = (const int32_t *)qs_read_exact(arena, path, nv, sizeof(int32_t));
    if (in->lane == NULL) in->lane = in->reconstruction;
    snprintf(path, sizeof path, "%s_stats.json", stem);
    {
        const char *err = NULL;
        Arena_Mark mark = Arena_save(arena);
        const JsonValue *root = Json_parse_file(arena, path, &err);
        if (root != NULL) {
            const JsonValue *sl = Json_object_get(root, "slicing");
            long ns = sl != NULL ? Json_member_long(sl, "n_slices", 0) : 0;
            if (ns > 0 && ns < INT_MAX) in->stats_n_slices = (int)ns;
        }
        Arena_restore(arena, mark);
    }
    /* the box is the OBSERVED extent: fills of an earlier solidify must not
     * widen the box a refit then fills up to, or the fill creeps one margin
     * per round */
    for (int c = 0; c < 3; c++) { in->bbox_lo[c] = 1e300; in->bbox_hi[c] = -1e300; }
    for (size_t i = 0; i < nv; i++) {
        if (in->prov[i] == QS_PROV_FILL) continue;
        for (int c = 0; c < 3; c++) {
            double p = in->mesh.verts[i * 3 + (size_t)c];
            if (p < in->bbox_lo[c]) in->bbox_lo[c] = p;
            if (p > in->bbox_hi[c]) in->bbox_hi[c] = p;
        }
    }
    if (in->bbox_lo[0] > in->bbox_hi[0]) {
        fprintf(stderr, "[solidify] %s carries no observation at all\n", in_vmesh);
        return -1;
    }
    rep->nv_in = nv;
    rep->nf_in = in->mesh.nf;
    for (size_t i = 0; i < nv; i++) {
        if (in->prov[i] == QS_PROV_FILL) rep->input_fills++;
        else if (in->support[i] != 0) rep->observed_fixed++;
        else rep->observed_generated++;
    }
    return 0;
}

/* ---- lattice recovery ---------------------------------------------------- */

typedef struct QsLattice {
    double du, dv, slice_h;
    long col_min, col_max;
    int n_slices, n_bands, k_min;
    int32_t *col;      /* [nv] atlas column */
    int32_t *slice;    /* [nv] slice k in [0, n_slices) */
    int8_t  *band;     /* [nv] peel band */
    float   *row_z;    /* [n_slices] plane z of each slice */
    size_t   row_z_deviations;
    int32_t *col_off;  /* vertices bucketed by column: [ncols+1] */
    int32_t *col_vert; /* [nv] */
} QsLattice;

static int qs_recover_lattice(Arena_T arena, const QsInput *in, QsLattice *L,
                              QsReport *rep)
{
    const MeshBinData *m = &in->mesh;
    size_t nv = m->nv;
    long vrow_max = 0;
    double *z = NULL;
    size_t nz = 0;
    double zmin = 0.0;
    int32_t *kz = NULL;
    memset(L, 0, sizeof *L);
    if (nv == 0) return -1;
    L->du = qs_detect_pitch(m->uv, nv, 0);
    L->dv = qs_detect_pitch(m->uv, nv, 1);
    L->col = (int32_t *)ARENA_ALLOC(arena, nv * sizeof *L->col);
    L->slice = (int32_t *)ARENA_ALLOC(arena, nv * sizeof *L->slice);
    L->band = (int8_t *)ARENA_ALLOC(arena, nv);
    kz = (int32_t *)ARENA_ALLOC(arena, nv * sizeof *kz);
    L->col_min = LONG_MAX;
    L->col_max = LONG_MIN;
    {
        size_t bad = 0;
        for (size_t i = 0; i < nv; i++) {
            double uq = (double)m->uv[i * 2] / L->du;
            double vq = (double)m->uv[i * 2 + 1] / L->dv;
            long c = lround(uq), r = lround(vq);
            if (fabs(uq - (double)c) > QS_LATTICE_EPS ||
                fabs(vq - (double)r) > QS_LATTICE_EPS || r < 0)
                bad++;
            L->col[i] = (int32_t)c;
            L->slice[i] = (int32_t)r;      /* vrow for now */
            if (c < L->col_min) L->col_min = c;
            if (c > L->col_max) L->col_max = c;
            if (r > vrow_max) vrow_max = r;
        }
        if (bad != 0) {
            fprintf(stderr, "[solidify] %zu vertices are off the UV lattice "
                            "(du=%.3f dv=%.3f); refusing a non-lattice ribbon\n",
                    bad, L->du, L->dv);
            return -1;
        }
    }
    /* slice planes from z: distinct values, median gap = slice_h */
    z = (double *)malloc(nv * sizeof *z);
    if (z == NULL) return -1;
    for (size_t i = 0; i < nv; i++) z[i] = (double)m->verts[i * 3];
    qsort(z, nv, sizeof *z, qs_cmp_double);
    nz = 1;
    for (size_t i = 1; i < nv; i++)
        if (z[i] > z[nz - 1] + 1e-3) z[nz++] = z[i];
    zmin = z[0];
    if (nz >= 2) {
        double *gaps = (double *)malloc((nz - 1) * sizeof *gaps);
        if (gaps == NULL) { free(z); return -1; }
        for (size_t i = 1; i < nz; i++) gaps[i - 1] = z[i] - z[i - 1];
        qsort(gaps, nz - 1, sizeof *gaps, qs_cmp_double);
        L->slice_h = gaps[(nz - 1) / 2];
        free(gaps);
    } else {
        L->slice_h = L->dv;
    }
    free(z);
    if (!(L->slice_h > 1e-6)) L->slice_h = L->dv;
    for (size_t i = 0; i < nv; i++)
        kz[i] = (int32_t)lround(((double)m->verts[i * 3] - zmin) / L->slice_h);

    /* offsets vrow - kz: one cluster per peel band, np apart */
    {
        long off_min = LONG_MAX, off_max = LONG_MIN;
        size_t *hist = NULL, range = 0;
        long np = 0;
        for (size_t i = 0; i < nv; i++) {
            long off = (long)L->slice[i] - (long)kz[i];
            if (off < off_min) off_min = off;
            if (off > off_max) off_max = off;
        }
        range = (size_t)(off_max - off_min + 1);
        if (range > (size_t)1 << 22) {
            fprintf(stderr, "[solidify] z/row offsets span %zu values; the "
                            "input is not a slice lattice\n", range);
            return -1;
        }
        hist = (size_t *)calloc(range, sizeof *hist);
        if (hist == NULL) return -1;
        for (size_t i = 0; i < nv; i++)
            hist[(size_t)((long)L->slice[i] - (long)kz[i] - off_min)]++;
        /* cluster nearby offsets (a two-cell tolerance absorbs off-plane
         * generated vertices); each cluster's mode is one band's offset */
        {
            long modes[64];
            int nmodes = 0;
            size_t floor_count = nv / 1000 > 2 ? nv / 1000 : 2;
            size_t i = 0;
            while (i < range && nmodes < 64) {
                size_t j = i, total = 0, best = 0;
                long best_off = off_min + (long)i;
                if (hist[i] == 0) { i++; continue; }
                while (j < range &&
                       (hist[j] != 0 || (j + 1 < range && hist[j + 1] != 0) ||
                        (j + 2 < range && hist[j + 2] != 0))) {
                    total += hist[j];
                    if (hist[j] > best) { best = hist[j]; best_off = off_min + (long)j; }
                    j++;
                }
                if (total >= floor_count) modes[nmodes++] = best_off;
                i = j;
            }
            L->k_min = nmodes > 0 ? (int)modes[0] : 0;
            if (in->stats_n_slices > 0) {
                np = in->stats_n_slices;
                if (nmodes >= 2 && modes[1] - modes[0] != np)
                    fprintf(stderr, "[solidify] WARNING: stats n_slices=%ld but the "
                                    "z lattice suggests %ld rows per band\n",
                            np, modes[1] - modes[0]);
            } else if (nmodes >= 2) {
                np = modes[1] - modes[0];
            } else {
                np = vrow_max + 1;
            }
            if (np <= 0) np = vrow_max + 1;
        }
        free(hist);
        L->n_slices = (int)np;
        L->n_bands = (int)(vrow_max / np) + 1;
        for (size_t i = 0; i < nv; i++) {
            long vrow = L->slice[i];
            long b = vrow / np;
            L->band[i] = (int8_t)(b > 127 ? 127 : b);
            L->slice[i] = (int32_t)(vrow - b * np);
        }
    }
    /* per-slice plane z: the first band-0 vertex of each slice defines it */
    L->row_z = (float *)ARENA_ALLOC(arena, (size_t)L->n_slices * sizeof(float));
    {
        uint8_t *seen = (uint8_t *)ARENA_CALLOC(arena, (size_t)L->n_slices, 1);
        for (size_t i = 0; i < nv; i++) {
            int k = L->slice[i];
            if (!seen[k]) { seen[k] = 1; L->row_z[k] = m->verts[i * 3]; }
            else if (fabs((double)m->verts[i * 3] - (double)L->row_z[k]) > 1e-3)
                L->row_z_deviations++;
        }
        for (int k = 0; k < L->n_slices; k++)
            if (!seen[k])
                L->row_z[k] = (float)(zmin + (double)(k - L->k_min) * L->slice_h);
    }
    /* bucket vertices by column for the per-run scans */
    {
        size_t ncols = (size_t)(L->col_max - L->col_min + 1);
        L->col_off = (int32_t *)ARENA_CALLOC(arena, ncols + 1, sizeof(int32_t));
        L->col_vert = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
        for (size_t i = 0; i < nv; i++)
            L->col_off[(size_t)(L->col[i] - L->col_min) + 1]++;
        for (size_t c = 0; c < ncols; c++) L->col_off[c + 1] += L->col_off[c];
        {
            int32_t *fill = (int32_t *)ARENA_ALLOC(arena, ncols * sizeof(int32_t));
            memcpy(fill, L->col_off, ncols * sizeof(int32_t));
            for (size_t i = 0; i < nv; i++) {
                size_t c = (size_t)(L->col[i] - L->col_min);
                L->col_vert[fill[c]++] = (int32_t)i;
            }
        }
    }
    rep->du = L->du;
    rep->dv = L->dv;
    rep->slice_h = L->slice_h;
    rep->n_slices = L->n_slices;
    rep->n_bands = L->n_bands;
    fprintf(stderr,
            "[solidify] lattice: du=%.3f dv=%.3f slice_h=%.3f cols=[%ld,%ld] "
            "slices=%d bands=%d k_min=%d (%s)%s\n",
            L->du, L->dv, L->slice_h, L->col_min, L->col_max, L->n_slices,
            L->n_bands, L->k_min,
            in->stats_n_slices > 0 ? "n_slices from stats sidecar"
                                   : "n_slices inferred from the z lattice",
            L->row_z_deviations != 0 ? " [row-z deviations present]" : "");
    if (L->row_z_deviations != 0)
        fprintf(stderr, "[solidify] WARNING: %zu vertices sit off their slice "
                        "plane by >1e-3 vox\n", L->row_z_deviations);
    return 0;
}

/* ---- output accumulator -------------------------------------------------- */

typedef struct QsOut {
    size_t nv, nf, cap_v, cap_f;
    float *verts, *uv;
    int32_t *faces;
    uint8_t *support, *prov;
    float *phase;
    int32_t *material, *run, *lane;
    int32_t *atlas_col;   /* the certificate's column before the layout */
} QsOut;

typedef struct QsRunStat {
    long c0, c1;
    long c0_out;              /* first output column after layout (-1 = not emitted) */
    size_t rows_occupied;     /* rows with at least one observation */
    double density_occ;       /* observed / (rows_occupied x width) */
    size_t snap_candidates, snap_moved, snap_no_ridge;
    size_t fold_cells, fold_columns;
    size_t dark_cells, dark_snapped, dark_unrecovered;
    size_t hole_cells, refused_hole_cells;
    size_t midline_probed, midline_outliers;
    size_t midline_isolated, midline_dense;
    size_t holes, holes_refused, holes_filled;   /* hole-level fill votes */
    int width, build;         /* build: 0 observed-only, 1 built, 2 retried */
    size_t observed, adopted, rejected, no_ref, duplicates, filled;
    size_t dropped_bbox, dropped_core, dropped_unref;
    size_t dropped_unbridged, dropped_skip, dropped_chord, dropped_material;
    size_t dropped_order, dropped_confetti;
    size_t dropped_fill_far;          /* FILL CAP: unsupported fill cells not emitted */
    size_t obs_smooth_cells;
    size_t column_probed, column_jumps;
    size_t depth_probed, depth_outliers;
    size_t ridge_cells, ridge_moved, ridge_on_before, ridge_on_after;
    size_t ridge_fill_cells, ridge_fill_on_before, ridge_fill_on_after;
    double ridge_shift_p50, ridge_shift_p90, ridge_shift_max;
    size_t faces, skipped_gate, skipped_material;
    double ramp_pitches_per_turn;      /* RAMP gate: |dr/dturn| / pitch, NAN = not measured */
    double ramp_turns;                 /* turns the run spans */
    size_t scale_faces, scale_over, scale_over_fill;   /* SCALE gate */
    size_t fill_cells, fill_far;                       /* FILL-SUPPORT gate */
    size_t fill_dist_max;
    double t_build;
} QsRunStat;

typedef struct QsDup { double *v; size_t n, cap; } QsDup;

static void qs_dup_push(QsDup *d, double x)
{
    if (d->n == d->cap) {
        size_t nc = d->cap ? d->cap * 2 : 1024;
        double *nb = (double *)realloc(d->v, nc * sizeof *nb);
        if (nb == NULL) return;
        d->v = nb;
        d->cap = nc;
    }
    d->v[d->n++] = x;
}

/* provenance image classes */
enum { QS_IMG_EMPTY = 0, QS_IMG_FIXED, QS_IMG_GENERATED, QS_IMG_BAND1,
       QS_IMG_FILL, QS_IMG_DROPPED, QS_IMG_DUPLICATE, QS_IMG_REJECTED,
       QS_IMG_REFUSED, QS_IMG_CONFETTI, QS_IMG_FOLD, QS_IMG_DARKSNAP,
       QS_IMG_MIDLINE_ISO,      /* MIDLINE isolated suspect, kept (QS_DARK_DEMOTE off) */
       QS_IMG_MIDLINE_DENSE };  /* MIDLINE dense suspect: the fit's depth error, kept */
#define QS_IMG_COUNT 14

typedef struct QsCtx {
    const QsConfig *cfg;
    const QsInput *in;
    const QsLattice *L;
    QsOut out;
    QsDup dup;
    uint8_t *img;         /* [n_slices * ncols] */
    size_t ncols;
    QsRunStat *runs;
    int n_runs;
    QsReport *rep;
    long col_out;             /* LAYOUT: next free output column */
    ZarrU8 *zarr;             /* RAW CT for the fill snap, or NULL */
    QsDup snap;               /* |displacement| samples of snapped fills */
    QsDup midline;            /* |ridge offset| samples of kept observed cells */
    QsDup ramp;               /* |pitches of radius gained per turn| over runs */
    QsDup ridge;              /* per-run median ridge shift */
    QsDup depth;              /* per-run median |radius - local median| */
    QsDup column;             /* per-run median |radius - column median|, pitches */
    QsDup obs_smooth;         /* per-run median |radial move| of the depth smoothing */
    QsDup scale;              /* |log2(area3d / area_uv)| over emitted triangles */
    QsDup scale_obs;          /* ... observed-only triangles */
    QsDup scale_fill;         /* ... triangles touching a fill */
    QsDup fill_dist;          /* BFS distance of emitted fill cells to an observation */
    size_t ct_coordinate_refusals;
} QsCtx;

/* Straightened -> world, for CT probes.  A no-op when the lane did not
 * straighten (axis_curve_n == 0). */
static void qs_ct_world(const QsCtx *ctx, const double *p, double out[3])
{
    const QsConfig *cfg = ctx->cfg;
    size_t lo = 0, hi = 0;
    double t = 0.0, y = 0.0, x = 0.0;
    if (cfg->coordinate_warp) {
        uint32_t flags=0;
        if (AxisWarp_to_world(cfg->coordinate_warp,p,out,&flags)!=0) {
            ((QsCtx*)ctx)->ct_coordinate_refusals++;
            out[0]=out[1]=out[2]=NAN;
        }
        return;
    }
    out[0] = p[0]; out[1] = p[1]; out[2] = p[2];
    if (cfg->axis_curve_n == 0 || cfg->axis_curve_z == NULL) return;
    if (p[0] <= cfg->axis_curve_z[0]) {
        y = cfg->axis_curve_y[0]; x = cfg->axis_curve_x[0];
    } else if (p[0] >= cfg->axis_curve_z[cfg->axis_curve_n - 1]) {
        y = cfg->axis_curve_y[cfg->axis_curve_n - 1];
        x = cfg->axis_curve_x[cfg->axis_curve_n - 1];
    } else {
        hi = 1;
        while (hi < cfg->axis_curve_n && cfg->axis_curve_z[hi] < p[0]) hi++;
        lo = hi - 1;
        t = (cfg->axis_curve_z[hi] - cfg->axis_curve_z[lo]) > 1e-9
          ? (p[0] - cfg->axis_curve_z[lo]) / (cfg->axis_curve_z[hi] - cfg->axis_curve_z[lo])
          : 0.0;
        y = cfg->axis_curve_y[lo] + t * (cfg->axis_curve_y[hi] - cfg->axis_curve_y[lo]);
        x = cfg->axis_curve_x[lo] + t * (cfg->axis_curve_x[hi] - cfg->axis_curve_x[lo]);
    }
    out[1] = p[1] + y - cfg->axis_ref_y;
    out[2] = p[2] + x - cfg->axis_ref_x;
}

/* CT at a straightened position. */
static double qs_ct_at(const QsCtx *ctx, double pz, double py, double px)
{
    double p[3], w[3];
    p[0] = pz; p[1] = py; p[2] = px;
    qs_ct_world(ctx, p, w);
    if (!isfinite(w[0]) || !isfinite(w[1]) || !isfinite(w[2])) return NAN;
    return (double)ZarrU8_sample(ctx->zarr, w[0], w[1], w[2]);
}

/* RidgeTrackSampler over the lattice frame */
int QuadribbonSolidify_probe_ct(const QsConfig *cfg, const double *points,
                                size_t count, double *world, uint8_t *values)
{
    QsCtx context={0}; ZarrU8 volume={0}; int result=-1;
    if (!cfg || !cfg->raw_zarr || !points || !world || !values) return -1;
    if (ZarrU8_open(&volume,cfg->raw_zarr)!=0) return -1;
    context.cfg=cfg; context.zarr=&volume;
    for (size_t i=0;i<count;i++) {
        qs_ct_world(&context,points+3*i,world+3*i);
        if (context.ct_coordinate_refusals) goto done;
        for (int a=0;a<3;a++) {
            double index=floor(world[3*i+a]+.5);
            if (!isfinite(index) || index<0. || index>=volume.shape[a]) goto done;
        }
        double value=qs_ct_at(&context,points[3*i],points[3*i+1],points[3*i+2]);
        if (!isfinite(value) || context.ct_coordinate_refusals || volume.missing) goto done;
        values[i]=(uint8_t)value;
    }
    result=0;
done:
    ZarrU8_close(&volume); return result;
}

/* RidgeTrackSampler over the lattice frame */
static double qs_ridge_sample(void *vctx, double z, double y, double x)
{
    return qs_ct_at((const QsCtx *)vctx, z, y, x);
}

/* RidgeTrackToWorld: the verdict measures the world twin, so edge lengths must
 * be judged there -- a displacement in z changes the straightening offset. */
static void qs_ridge_world(void *vctx, const double in[3], double out[3])
{
    qs_ct_world((const QsCtx *)vctx, in, out);
}

/* Local axis tangent at world z, from the curve's own slope.  Returns 0 when no
 * curve is armed (the caller then uses the z-plane radius). */
static int qs_axis_tangent(const QsCtx *ctx, double z, double t[3], double centre[2])
{
    const QsConfig *cfg = ctx->cfg;
    size_t lo = 0, hi = 0, n = cfg->axis_curve_n;
    double dz = 0.0, dy = 0.0, dx = 0.0, len = 0.0, u = 0.0;
    if (n < 2 || cfg->axis_curve_z == NULL) return 0;
    if (z <= cfg->axis_curve_z[0]) { lo = 0; hi = 1; }
    else if (z >= cfg->axis_curve_z[n - 1]) { lo = n - 2; hi = n - 1; }
    else {
        hi = 1;
        while (hi < n && cfg->axis_curve_z[hi] < z) hi++;
        lo = hi - 1;
    }
    dz = cfg->axis_curve_z[hi] - cfg->axis_curve_z[lo];
    if (!(fabs(dz) > 1e-9)) return 0;
    dy = (cfg->axis_curve_y[hi] - cfg->axis_curve_y[lo]) / dz;
    dx = (cfg->axis_curve_x[hi] - cfg->axis_curve_x[lo]) / dz;
    len = sqrt(1.0 + dy * dy + dx * dx);
    t[0] = 1.0 / len; t[1] = dy / len; t[2] = dx / len;
    u = (z - cfg->axis_curve_z[lo]) / dz;
    centre[0] = cfg->axis_curve_y[lo] + u * (cfg->axis_curve_y[hi] - cfg->axis_curve_y[lo]);
    centre[1] = cfg->axis_curve_x[lo] + u * (cfg->axis_curve_x[hi] - cfg->axis_curve_x[lo]);
    return 1;
}

static double qs_radius(const QsCtx *ctx, const float *p)
{
    if (QS_PERP_RADIUS && ctx->cfg->axis_curve_n >= 2) {
        /* the lattice frame is straightened; go back to world, then measure
         * perpendicular to the local tangent */
        double lp[3], w[3], t[3], ctr[2], d[3], dt = 0.0, perp[3];
        lp[0] = (double)p[0]; lp[1] = (double)p[1]; lp[2] = (double)p[2];
        qs_ct_world(ctx, lp, w);
        if (qs_axis_tangent(ctx, w[0], t, ctr)) {
            d[0] = 0.0; d[1] = w[1] - ctr[0]; d[2] = w[2] - ctr[1];
            dt = d[0] * t[0] + d[1] * t[1] + d[2] * t[2];
            perp[0] = d[0] - dt * t[0];
            perp[1] = d[1] - dt * t[1];
            perp[2] = d[2] - dt * t[2];
            return sqrt(perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2]);
        }
    }
    return hypot((double)p[1] - ctx->cfg->axis_y, (double)p[2] - ctx->cfg->axis_x);
}

static int qs_edge_ok(const float *pa, const float *pb)
{
    double dz = (double)pa[0] - pb[0], dy = (double)pa[1] - pb[1];
    double dx = (double)pa[2] - pb[2];
    return dz * dz + dy * dy + dx * dx <= QS_WRAP_GATE_VOX * QS_WRAP_GATE_VOX;
}

static int qs_mat_ok(int32_t a, int32_t b, int32_t c)
{
    int32_t m = a >= 0 ? a : (b >= 0 ? b : c);
    if (a >= 0 && a != m) return 0;
    if (b >= 0 && b != m) return 0;
    if (c >= 0 && c != m) return 0;
    return 1;
}

/* Legality of one triangle over cells (ca,cb,cc): 1 legal, 0 gate, -1 material */
static int qs_tri_legal(const float *pos, const int32_t *mat, size_t ca,
                        size_t cb, size_t cc)
{
    const float *pa = pos + ca * 3, *pb = pos + cb * 3, *pc = pos + cc * 3;
    if (!qs_edge_ok(pa, pb) || !qs_edge_ok(pb, pc) || !qs_edge_ok(pa, pc))
        return 0;
    if (!qs_mat_ok(mat[ca], mat[cb], mat[cc])) return -1;
    return 1;
}

/* Wrap-skip test for a bridge between observed cells a and b that are
 * `steps` lattice steps apart along an axis (stride 1 = row, W = column):
 * each side's radius is extrapolated across the bridge from the least-squares
 * slope over up to `reach` observed cells behind it; the bridge is consistent
 * when either extrapolation lands within QS_SKIP_DR_PITCH pitches of the
 * other side.  A pure pitch offset (the next wrap) fails both. */
static int qs_bridge_consistent(const QsCtx *ctx, const uint8_t *domain,
                                const float *pos, size_t a, size_t b,
                                size_t stride, int steps, int lo_index,
                                int hi_index, int hi_limit, int reach,
                                double pitch, double tol_pitch)
{
    double ra = qs_radius(ctx, pos + a * 3), rb = qs_radius(ctx, pos + b * 3);
    double sa = 0.0, sb = 0.0;
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int n = 0;
    /* slope behind a (towards lower index) */
    for (int k = 0; k < reach && lo_index - k >= 0; k++) {
        size_t c = a - (size_t)k * stride;
        double x = -(double)k, y = 0.0;
        if (!domain[c]) continue;
        y = qs_radius(ctx, pos + c * 3);
        sx += x; sy += y; sxx += x * x; sxy += x * y; n++;
    }
    if (n >= 2 && n * sxx - sx * sx > 1e-9) sa = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    sx = sy = sxx = sxy = 0; n = 0;
    /* slope beyond b (towards higher index) */
    for (int k = 0; k < reach && hi_index + k < hi_limit; k++) {
        size_t c = b + (size_t)k * stride;
        double x = (double)k, y = 0.0;
        if (!domain[c]) continue;
        y = qs_radius(ctx, pos + c * 3);
        sx += x; sy += y; sxx += x * x; sxy += x * y; n++;
    }
    if (n >= 2 && n * sxx - sx * sx > 1e-9) sb = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    if (fabs(ra + sa * (double)steps - rb) <= tol_pitch * pitch) return 1;
    if (fabs(rb - sb * (double)steps - ra) <= tol_pitch * pitch) return 1;
    return 0;
}

/* Expected lifted-phase change across ncols columns between cells a and b:
 * u is arc length, so the angle per column is du / r at the LOCAL radius.
 * The run-median angle assumed one radius for a run spanning many wraps and
 * refused most gaps for "order" (19,479 cells in one 4x5x5 run, 2026-09-02). */
static double qs_want_dtheta(const QsCtx *ctx, const QsLattice *L, const float *pos,
                             size_t a, size_t b, int ncols, double dtheta_col)
{
    double r = 0.5 * (qs_radius(ctx, pos + a * 3) + qs_radius(ctx, pos + b * 3));
    double sgn = dtheta_col < 0.0 ? -1.0 : 1.0;
    if (!(r > 1.0) || !isfinite(r)) return dtheta_col * (double)ncols;
    return sgn * (double)ncols * L->du / r;
}

/* Unit lattice normal at cell c from its u and v neighbours (0 = degenerate). */
static int qs_lattice_normal(const float *pos, int W, int H, size_t c, double out[3])
{
    int j = (int)(c / (size_t)W), i = (int)(c % (size_t)W);
    const float *pl = pos + (i > 0 ? c - 1 : c) * 3, *pr = pos + (i + 1 < W ? c + 1 : c) * 3;
    const float *pd = pos + (j > 0 ? c - (size_t)W : c) * 3, *pu = pos + (j + 1 < H ? c + (size_t)W : c) * 3;
    double eu[3], ev[3], len = 0.0;
    if (!isfinite((double)pl[1]) || !isfinite((double)pr[1]) ||
        !isfinite((double)pd[1]) || !isfinite((double)pu[1])) return 0;
    for (int a = 0; a < 3; a++) { eu[a] = (double)pr[a] - pl[a]; ev[a] = (double)pu[a] - pd[a]; }
    out[0] = eu[1] * ev[2] - eu[2] * ev[1];
    out[1] = eu[2] * ev[0] - eu[0] * ev[2];
    out[2] = eu[0] * ev[1] - eu[1] * ev[0];
    len = sqrt(out[0] * out[0] + out[1] * out[1] + out[2] * out[2]);
    if (!(len > 1e-9)) return 0;
    for (int a = 0; a < 3; a++) out[a] /= len;
    return 1;
}

static int qs_process_run(Arena_T arena, QsCtx *ctx, int run_index,
                          long c0, long c1)
{
    const QsInput *in = ctx->in;
    const QsLattice *L = ctx->L;
    const MeshBinData *m = &in->mesh;
    QsRunStat *rs = &ctx->runs[run_index];
    int W = (int)(c1 - c0 + 1), H = L->n_slices;
    size_t HW = (size_t)H * (size_t)W;
    Arena_Mark mark = Arena_save(arena);
    int32_t *sel = NULL, *pend = NULL, *dist = NULL, *src = NULL, *queue = NULL;
    uint8_t *is_band0 = NULL, *domain = NULL, *prov = NULL, *keep = NULL, *used = NULL;
    uint8_t *folded = NULL;      /* fold-demoted cells: interpolated without vetoes */
    float *field = NULL, *phase = NULL, *pos = NULL, *bverts = NULL, *bphase = NULL;
    float *pos_obs = NULL;       /* observed positions only (NaN elsewhere): MIDLINE probe */
    uint8_t *bfilled = NULL;
    int32_t *mat = NULL, *lane = NULL, *tri = NULL, *oidx = NULL;
    size_t bnv = 0, ntri = 0, n_obs = 0;
    int built = 0, rc = 0;
    double t_build0 = 0.0;
    const double pitch = ctx->cfg->pitch > 0.0 ? ctx->cfg->pitch : 9.5;

    memset(rs, 0, sizeof *rs);
    rs->c0 = c0; rs->c1 = c1; rs->width = W;
    if (HW == 0 || HW > (size_t)INT32_MAX) { Arena_restore(arena, mark); return -1; }

    sel = (int32_t *)ARENA_ALLOC(arena, HW * sizeof *sel);
    pend = (int32_t *)ARENA_ALLOC(arena, HW * sizeof *pend);
    is_band0 = (uint8_t *)ARENA_CALLOC(arena, HW, 1);
    domain = (uint8_t *)ARENA_CALLOC(arena, HW, 1);
    prov = (uint8_t *)ARENA_ALLOC(arena, HW);
    field = (float *)ARENA_ALLOC(arena, HW * 3 * sizeof *field);
    phase = (float *)ARENA_ALLOC(arena, HW * sizeof *phase);
    mat = (int32_t *)ARENA_ALLOC(arena, HW * sizeof *mat);
    lane = (int32_t *)ARENA_ALLOC(arena, HW * sizeof *lane);
    for (size_t c = 0; c < HW; c++) {
        sel[c] = -1; pend[c] = -1; prov[c] = QS_NONE; mat[c] = -1; lane[c] = -1;
        phase[c] = NAN;
        field[c * 3] = field[c * 3 + 1] = field[c * 3 + 2] = NAN;
    }

    /* pass 1: band-0 observations claim cells; other bands wait */
    for (long col = c0; col <= c1; col++) {
        size_t cc = (size_t)(col - L->col_min);
        for (int32_t e = L->col_off[cc]; e < L->col_off[cc + 1]; e++) {
            size_t v = (size_t)L->col_vert[e];
            size_t cell = (size_t)L->slice[v] * (size_t)W + (size_t)(col - c0);
            if (in->prov[v] == QS_PROV_FILL) continue;   /* refit: unobserved */
            if (L->band[v] == 0) {
                if (sel[cell] < 0) { sel[cell] = (int32_t)v; is_band0[cell] = 1; }
                else ctx->rep->double_cover++;
            } else if (pend[cell] < 0) {
                pend[cell] = (int32_t)v;
            }
        }
    }
    /* pass 2: duplicates (both bands present) and band-1 adoption */
    for (size_t c = 0; c < HW; c++) {
        int j = (int)(c / (size_t)W), i = (int)(c % (size_t)W);
        uint8_t *px = ctx->img + (size_t)j * ctx->ncols + (size_t)(c0 + i - L->col_min);
        if (pend[c] < 0) continue;
        if (sel[c] >= 0) {
            double r0 = qs_radius(ctx, m->verts + (size_t)sel[c] * 3);
            double r1 = qs_radius(ctx, m->verts + (size_t)pend[c] * 3);
            qs_dup_push(&ctx->dup, fabs(r1 - r0) / pitch);
            rs->duplicates++;
            *px = QS_IMG_DUPLICATE;
            continue;
        }
        {
            /* nearest band-0 reference within the window */
            int best_d2 = INT_MAX;
            size_t best = SIZE_MAX;
            for (int dj = -QS_BAND1_NEIGHBOUR_CELLS; dj <= QS_BAND1_NEIGHBOUR_CELLS; dj++) {
                int jj = j + dj;
                if (jj < 0 || jj >= H) continue;
                for (int di = -QS_BAND1_NEIGHBOUR_CELLS; di <= QS_BAND1_NEIGHBOUR_CELLS; di++) {
                    int ii = i + di, d2 = dj * dj + di * di;
                    size_t n = 0;
                    if (ii < 0 || ii >= W) continue;
                    n = (size_t)jj * (size_t)W + (size_t)ii;
                    if (!is_band0[n] || d2 >= best_d2) continue;
                    best_d2 = d2;
                    best = n;
                }
            }
            if (best == SIZE_MAX) {
                rs->no_ref++;
                *px = QS_IMG_REJECTED;
            } else {
                double r0 = qs_radius(ctx, m->verts + (size_t)sel[best] * 3);
                double r1 = qs_radius(ctx, m->verts + (size_t)pend[c] * 3);
                if (fabs(r1 - r0) <= QS_BAND1_ADOPT_DR_PITCH * pitch) {
                    sel[c] = pend[c];
                    prov[c] = QS_PROV_BAND1;
                    rs->adopted++;
                } else {
                    rs->rejected++;
                    *px = QS_IMG_REJECTED;
                }
            }
        }
    }
    /* pass 3: the observation field */
    for (size_t c = 0; c < HW; c++) {
        size_t v = 0;
        int j = (int)(c / (size_t)W);
        if (sel[c] < 0) continue;
        v = (size_t)sel[c];
        domain[c] = 1;
        field[c * 3 + 0] = (float)j;                  /* row-unit z for quad_strip */
        field[c * 3 + 1] = m->verts[v * 3 + 1];
        field[c * 3 + 2] = m->verts[v * 3 + 2];
        phase[c] = in->phase[v];
        if (prov[c] == QS_NONE) prov[c] = in->prov[v];
        mat[c] = in->material != NULL ? in->material[v] : -1;
        lane[c] = in->lane != NULL ? in->lane[v] : -1;
        n_obs++;
        if (!isfinite((double)phase[c]) || !isfinite((double)field[c * 3 + 1]) ||
            !isfinite((double)field[c * 3 + 2])) {
            fprintf(stderr, "[solidify] run %d: non-finite observation at cell "
                            "(row %d, col %ld)\n", run_index, j,
                    c0 + (long)(c % (size_t)W));
            Arena_restore(arena, mark);
            return -1;
        }
    }
    rs->observed = n_obs;
    /* FOLD demotion: a column whose observed cells barely move along u (the
     * lattice folded back on itself; the bake smears them into streaks) is
     * demoted to a hole and re-interpolated from its neighbours. */
    {
        /* per CELL: an observed cell whose u-step to its right neighbour is
         * shorter than QS_FOLD_STEP_FRAC x du is a fold duplicate of that
         * neighbour (the smear spreads across columns, so a per-column rule
         * never fired); the right cell of the pair is demoted */
        const double step_min = QS_FOLD_STEP_FRAC * L->du;
        uint8_t *demote = (uint8_t *)ARENA_CALLOC(arena, HW, 1);
        folded = (uint8_t *)ARENA_CALLOC(arena, HW, 1);
        for (size_t c = 0; c + 1 < HW; c++) {
            const float *p, *q;
            if ((c + 1) % (size_t)W == 0) continue;
            if (!domain[c] || !domain[c + 1]) continue;
            p = m->verts + (size_t)sel[c] * 3;
            q = m->verts + (size_t)sel[c + 1] * 3;
            if (sqrt(pow((double)p[0] - q[0], 2) + pow((double)p[1] - q[1], 2) +
                     pow((double)p[2] - q[2], 2)) < step_min)
                demote[c + 1] = 1;
        }
        for (size_t c = 0; c < HW; c++) {
            int j = (int)(c / (size_t)W), i = (int)(c % (size_t)W);
            uint8_t *px = ctx->img + (size_t)j * ctx->ncols + (size_t)(c0 + i - L->col_min);
            if (!demote[c] || !domain[c]) continue;
            domain[c] = 0; sel[c] = -1; prov[c] = QS_NONE; mat[c] = -1; lane[c] = -1;
            phase[c] = NAN;
            field[c * 3] = field[c * 3 + 1] = field[c * 3 + 2] = NAN;
            *px = QS_IMG_FOLD;
            folded[c] = 1;
            rs->fold_cells++;
            n_obs--;
        }
        rs->fold_columns = 0;
        rs->observed = n_obs;
    }
    rs->c0_out = -1;
    {
        /* density over the rows the run occupies: an arc fragment of a
         * wrap fills its rows well, confetti does not */
        size_t rows_occ = 0;
        for (int j = 0; j < H; j++) {
            size_t row = (size_t)j * (size_t)W;
            for (int i = 0; i < W; i++)
                if (domain[row + (size_t)i]) { rows_occ++; break; }
        }
        rs->rows_occupied = rows_occ;
        rs->density_occ = rows_occ > 0 ? (double)n_obs / ((double)rows_occ * (double)W) : 0.0;
    }
    if (n_obs < (size_t)(ctx->cfg->run_min_observed > 0 ? ctx->cfg->run_min_observed
                                                          : QS_RUN_MIN_OBSERVED) ||
        (n_obs < (size_t)QS_RUN_CONFETTI_MAX_OBS && rs->density_occ < QS_RUN_MIN_DENSITY)) {
        /* confetti: fragments the certificate could not place are not part of
         * the single ribbon (counted, painted, never emitted) */
        for (size_t c = 0; c < HW; c++) {
            int j = (int)(c / (size_t)W), i = (int)(c % (size_t)W);
            if (domain[c])
                ctx->img[(size_t)j * ctx->ncols + (size_t)(c0 + i - L->col_min)] = QS_IMG_CONFETTI;
        }
        rs->dropped_confetti = n_obs;
        rs->build = -1;
        fprintf(stderr, "[solidify] run %d/%d cols [%ld,%ld] w=%d: confetti (%zu observed "
                        "cells < %d, or < %d with density %.3f over %zu occupied rows < %.2f) "
                        "-- not emitted\n", run_index + 1,
                ctx->n_runs, c0, c1, W, n_obs,
                ctx->cfg->run_min_observed > 0 ? ctx->cfg->run_min_observed : QS_RUN_MIN_OBSERVED,
                QS_RUN_CONFETTI_MAX_OBS, rs->density_occ, rs->rows_occupied, QS_RUN_MIN_DENSITY);
        Arena_restore(arena, mark);
        return 0;
    }
    /* phase convention guard: the sidecar must agree with atan2 up to 2*pi */
    for (size_t c = 0; c < HW; c++) {
        double raw = 0.0, branch = 0.0, mm = 0.0;
        if (!domain[c]) continue;
        raw = atan2((double)field[c * 3 + 2] - ctx->cfg->axis_x,
                    (double)field[c * 3 + 1] - ctx->cfg->axis_y);
        branch = nearbyint(((double)phase[c] - raw) / QS_TWO_PI);
        mm = fabs(raw + branch * QS_TWO_PI - (double)phase[c]);
        if (mm > ctx->rep->phase_mismatch_max) ctx->rep->phase_mismatch_max = mm;
        if (mm > QS_PHASE_MISMATCH_MAX) ctx->rep->phase_mismatch_over_limit++;
    }

    /* build the run rectangle: observed exact, holes = cylindrical residual
     * harmonic continuation.  Two facts drive the arguments: quad_strip
     * demands z == z0 + row (so z is the row index and the plane z is
     * restored below) and weights its stencil 1/d^2 in its OWN uv (so
     * grid_du = 1 keeps the Laplacian isotropic; uv is rewritten below). */
    /* RAMP test: u is arc length and phase is the geometric angle about the
     * axis, so a run that follows ONE wrap gains exactly one pitch of radius
     * per turn.  Fit the observed cells' column medians of radius and phase,
     * and report the radius gained per turn in pitches.  A run that gains
     * several is walking across the wrap stack (the prediction fused the wraps
     * there); the reviewer sees a swirl. */
    {
        double sx = 0.0, sxx = 0.0, sr = 0.0, sxr = 0.0, sp = 0.0, sxp = 0.0;
        size_t nc = 0;
        double ramp_pitch = ctx->cfg->pitch > 0.0 ? ctx->cfg->pitch : 9.5;
        rs->ramp_pitches_per_turn = NAN;
        rs->ramp_turns = 0.0;
        for (int i = 0; i < W; i++) {
            double rsum = 0.0, psum = 0.0;
            size_t n = 0;
            for (int j = 0; j < H; j++) {
                size_t c = (size_t)j * (size_t)W + (size_t)i;
                if (sel[c] < 0 || !isfinite(phase[c])) continue;
                rsum += qs_radius(ctx, m->verts + (size_t)sel[c] * 3);
                psum += (double)phase[c];
                n++;
            }
            if (n < 3) continue;
            {
                double xv = (double)i, rv = rsum / (double)n, pv = psum / (double)n;
                sx += xv; sxx += xv * xv; sr += rv; sxr += xv * rv; sp += pv; sxp += xv * pv;
                nc++;
            }
        }
        if (nc >= 20) {
            double den = (double)nc * sxx - sx * sx;
            if (fabs(den) > 1e-9) {
                double dr = ((double)nc * sxr - sx * sr) / den;      /* vox per column */
                double dp = ((double)nc * sxp - sx * sp) / den;      /* rad per column */
                double turns = fabs(dp) * (double)W / QS_TWO_PI;
                rs->ramp_turns = turns;
                if (fabs(dp) > 1e-9 && turns >= QS_RAMP_MIN_TURNS) {
                    double dr_turn = dr / dp * QS_TWO_PI;
                    rs->ramp_pitches_per_turn = fabs(dr_turn) / ramp_pitch;
                    qs_dup_push(&ctx->ramp, rs->ramp_pitches_per_turn);
                }
            }
        }
    }

    /* COLUMN-CONTINUITY repair.  A lattice column is one angular position, so its
     * radius must vary smoothly down the rows: a wrap is a cylinder.  A cell
     * whose radius jumps a pitch or more away from its column's local median is
     * on a DIFFERENT WRAP -- the certificate's lift changed between adjacent z --
     * and every face touching it is then refused by the wrap gate, which is what
     * shatters the sheet.  Demote it; the fill bridges the rows instead. */
    if (QS_COLUMN_JUMP_PITCH > 0.0) {
        double pitch_c = ctx->cfg->pitch > 0.0 ? ctx->cfg->pitch : 9.5;
        double lim = QS_COLUMN_JUMP_PITCH * pitch_c;
        const int hwr = QS_COLUMN_JUMP_WIN / 2;
        double *crad = (double *)ARENA_ALLOC(arena, HW * sizeof *crad);
        double *cmed = (double *)ARENA_ALLOC(arena, HW * sizeof *cmed);
        double *cdev = (double *)ARENA_ALLOC(arena, HW * sizeof *cdev);
        size_t ndev = 0;
        if (crad != NULL && cmed != NULL && cdev != NULL) {
            for (size_t c = 0; c < HW; c++)
                crad[c] = sel[c] >= 0 ? qs_radius(ctx, m->verts + (size_t)sel[c] * 3) : NAN;
            for (int i = 0; i < W; i++)
                for (int j = 0; j < H; j++) {
                    size_t c = (size_t)j * (size_t)W + (size_t)i;
                    double buf[QS_COLUMN_JUMP_WIN];
                    int n = 0, k = 0;
                    if (sel[c] < 0) { cmed[c] = NAN; continue; }
                    for (k = -hwr; k <= hwr; k++) {
                        int jj = j + k;
                        size_t o = 0;
                        if (jj < 0 || jj >= H) continue;
                        o = (size_t)jj * (size_t)W + (size_t)i;
                        if (sel[o] < 0 || !isfinite(crad[o])) continue;
                        buf[n++] = crad[o];
                    }
                    if (n < 3) { cmed[c] = NAN; continue; }
                    qsort(buf, (size_t)n, sizeof *buf, qs_cmp_double);
                    cmed[c] = buf[n / 2];
                }
            if (QS_COLUMN_JUMP_2D) {
                /* second pass ALONG the columns: a wrong-wrap cell is a pitch off
                 * its neighbourhood in both directions, and the column-only test
                 * only sees it where the column carries enough rows */
                double *rowmed = (double *)ARENA_ALLOC(arena, HW * sizeof *rowmed);
                if (rowmed != NULL) {
                    for (int j = 0; j < H; j++)
                        for (int i = 0; i < W; i++) {
                            size_t c = (size_t)j * (size_t)W + (size_t)i;
                            double buf[QS_COLUMN_JUMP_WIN];
                            int n = 0, k = 0;
                            if (sel[c] < 0 || !isfinite(cmed[c])) { rowmed[c] = NAN; continue; }
                            for (k = -hwr; k <= hwr; k++) {
                                int ii = i + k;
                                size_t o = 0;
                                if (ii < 0 || ii >= W) continue;
                                o = (size_t)j * (size_t)W + (size_t)ii;
                                if (sel[o] < 0 || !isfinite(cmed[o])) continue;
                                buf[n++] = cmed[o];
                            }
                            if (n < 3) { rowmed[c] = cmed[c]; continue; }
                            qsort(buf, (size_t)n, sizeof *buf, qs_cmp_double);
                            rowmed[c] = buf[n / 2];
                        }
                    for (size_t c = 0; c < HW; c++)
                        if (isfinite(rowmed[c])) cmed[c] = rowmed[c];
                }
            }
            for (size_t c = 0; c < HW; c++) {
                if (sel[c] < 0 || !isfinite(crad[c]) || !isfinite(cmed[c])) continue;
                rs->column_probed++;
                cdev[ndev++] = fabs(crad[c] - cmed[c]) / pitch_c;
            }
            if (ndev > 0) {
                qsort(cdev, ndev, sizeof *cdev, qs_cmp_double);
                qs_dup_push(&ctx->column, cdev[ndev / 2]);
            }
            for (size_t c = 0; c < HW; c++) {
                int j = 0, i = 0;
                uint8_t *px = NULL;
                if (sel[c] < 0 || !isfinite(crad[c]) || !isfinite(cmed[c])) continue;
                if (fabs(crad[c] - cmed[c]) <= lim) continue;
                j = (int)(c / (size_t)W); i = (int)(c % (size_t)W);
                px = ctx->img + (size_t)j * ctx->ncols + (size_t)(c0 + i - L->col_min);
                if (QS_COLUMN_RELOCATE) {
                    /* move the cell radially onto the radius its column agrees on:
                     * the correct wrap's material is there (the wraps are separated),
                     * so this keeps the cell rather than costing its coverage */
                    const float *sp = m->verts + (size_t)sel[c] * 3;
                    double dy = (double)sp[1] - ctx->cfg->axis_y;
                    double dx = (double)sp[2] - ctx->cfg->axis_x;
                    double len = sqrt(dy * dy + dx * dx);
                    if (len > 1e-9) {
                        double scale = cmed[c] / len;
                        field[c * 3] = sp[0];
                        field[c * 3 + 1] = (float)(ctx->cfg->axis_y + dy * scale);
                        field[c * 3 + 2] = (float)(ctx->cfg->axis_x + dx * scale);
                        prov[c] = QS_PROV_OBSERVED_GENERATED;
                        rs->column_jumps++;
                        *px = QS_IMG_MIDLINE_ISO;
                        continue;
                    }
                }
                domain[c] = 0; sel[c] = -1; prov[c] = QS_NONE; mat[c] = -1; lane[c] = -1;
                phase[c] = NAN;
                field[c * 3] = field[c * 3 + 1] = field[c * 3 + 2] = NAN;
                *px = QS_IMG_MIDLINE_ISO;
                rs->column_jumps++;
                n_obs--;
            }
        }
    }

    /* DEPTH-OUTLIER demotion.  The verdict's TEXTURE/depth-seams gate high-passes
     * the per-cell radius with a 9-cell local median and counts adjacent steps
     * over 0.75 vox: blocks the claim DP placed a voxel or two proud of their
     * surroundings, which bake as rectangles of misregistered texture.  Here the
     * same high-pass is a REPAIR: an observed cell too far from its
     * neighbourhood is demoted, and the harmonic fill puts the surface through
     * its neighbours back. */
    if (QS_DEPTH_OUTLIER_VOX > 0.0) {
        double *rad = (double *)ARENA_ALLOC(arena, HW * sizeof *rad);
        double *med = (double *)ARENA_ALLOC(arena, HW * sizeof *med);
        double *tmp = (double *)ARENA_ALLOC(arena, HW * sizeof *tmp);
        double *hpv = (double *)ARENA_ALLOC(arena, HW * sizeof *hpv);
        const int hw = QS_DEPTH_OUTLIER_WIN / 2;
        size_t nhp = 0;
        if (rad != NULL && med != NULL && tmp != NULL && hpv != NULL) {
            for (size_t c = 0; c < HW; c++)
                rad[c] = sel[c] >= 0 ? qs_radius(ctx, m->verts + (size_t)sel[c] * 3) : NAN;
            /* separable median: along u, then along v */
            for (int j = 0; j < H; j++)
                for (int i = 0; i < W; i++) {
                    size_t c = (size_t)j * (size_t)W + (size_t)i;
                    double buf[QS_DEPTH_OUTLIER_WIN];
                    int n = 0, k = 0;
                    if (sel[c] < 0) { tmp[c] = NAN; continue; }
                    for (k = -hw; k <= hw; k++) {
                        int ii = i + k;
                        size_t o = 0;
                        if (ii < 0 || ii >= W) continue;
                        o = (size_t)j * (size_t)W + (size_t)ii;
                        if (sel[o] < 0) continue;
                        buf[n++] = rad[o];
                    }
                    if (n == 0) { tmp[c] = rad[c]; continue; }
                    qsort(buf, (size_t)n, sizeof *buf, qs_cmp_double);
                    tmp[c] = buf[n / 2];
                }
            for (int j = 0; j < H; j++)
                for (int i = 0; i < W; i++) {
                    size_t c = (size_t)j * (size_t)W + (size_t)i;
                    double buf[QS_DEPTH_OUTLIER_WIN];
                    int n = 0, k = 0;
                    if (sel[c] < 0) { med[c] = NAN; continue; }
                    for (k = -hw; k <= hw; k++) {
                        int jj = j + k;
                        size_t o = 0;
                        if (jj < 0 || jj >= H) continue;
                        o = (size_t)jj * (size_t)W + (size_t)i;
                        if (sel[o] < 0 || !isfinite(tmp[o])) continue;
                        buf[n++] = tmp[o];
                    }
                    if (n == 0) { med[c] = rad[c]; continue; }
                    qsort(buf, (size_t)n, sizeof *buf, qs_cmp_double);
                    med[c] = buf[n / 2];
                }
            for (size_t c = 0; c < HW; c++) {
                if (sel[c] < 0 || !isfinite(rad[c]) || !isfinite(med[c])) continue;
                hpv[nhp++] = fabs(rad[c] - med[c]);
                rs->depth_probed++;
            }
            if (nhp > 0) {
                qsort(hpv, nhp, sizeof *hpv, qs_cmp_double);
                qs_dup_push(&ctx->depth, hpv[nhp / 2]);
            }
            for (size_t c = 0; c < HW; c++) {
                int j = 0, i = 0;
                uint8_t *px = NULL;
                if (sel[c] < 0 || !isfinite(rad[c]) || !isfinite(med[c])) continue;
                if (fabs(rad[c] - med[c]) <= QS_DEPTH_OUTLIER_VOX) continue;
                j = (int)(c / (size_t)W); i = (int)(c % (size_t)W);
                px = ctx->img + (size_t)j * ctx->ncols + (size_t)(c0 + i - L->col_min);
                domain[c] = 0; sel[c] = -1; prov[c] = QS_NONE; mat[c] = -1; lane[c] = -1;
                phase[c] = NAN;
                field[c * 3] = field[c * 3 + 1] = field[c * 3 + 2] = NAN;
                *px = QS_IMG_MIDLINE_DENSE;
                rs->depth_outliers++;
                n_obs--;
            }
        }
    }

    /* MIDLINE test runs BEFORE the harmonic build so that a demoted cell is
     * genuinely refit from its neighbours (after the build its position
     * would already be its own observation). */
    pos_obs = (float *)ARENA_ALLOC(arena, HW * 3 * sizeof *pos_obs);
    for (size_t c = 0; c < HW; c++) {
        if (sel[c] >= 0) {
            const float *p = m->verts + (size_t)sel[c] * 3;
            pos_obs[c * 3] = p[0]; pos_obs[c * 3 + 1] = p[1]; pos_obs[c * 3 + 2] = p[2];
        } else {
            pos_obs[c * 3] = pos_obs[c * 3 + 1] = pos_obs[c * 3 + 2] = NAN;
        }
    }
    /* MIDLINE test (user rule): darkness counts against a cell in proportion
     * to its distance from the papyrus.  Each observed cell probes the CT
     * along its lattice normal: it is DARK when nothing within +-1 vox reaches
     * QS_DARK_U8, and its distance to the papyrus is the first sample at or
     * above QS_FILL_SNAP_MIN_U8 walking outward on both sides within
     * QS_MIDLINE_REACH.  A dark cell further than QS_MIDLINE_OUTLIER_VOX from
     * any papyrus, or any cell with no papyrus in reach, is SUSPECT: demoted
     * to a hole and refit through its neighbours by the harmonic
     * interpolation (the surface through the neighbours).  Dim cells that sit
     * on papyrus are kept.  The kept cells' distance to the papyrus face is
     * reported (p50 / p90). */
    if (ctx->zarr != NULL) {
        uint8_t *suspect_map = (uint8_t *)ARENA_CALLOC(arena, HW, 1);
        for (size_t c = 0; c < HW; c++) {
            double nn[3], near_v = 0.0, face_d = -1.0;
            const float *pc = pos_obs + c * 3;
            int suspect = 0;
            if (sel[c] < 0) continue;
            if (!qs_lattice_normal(pos_obs, W, H, c, nn)) continue;
            rs->midline_probed++;
            for (double t = -1.0; t <= 1.0 + 1e-9; t += 0.5) {
                double v = qs_ct_at(ctx, (double)pc[0] + t * nn[0],
                                    (double)pc[1] + t * nn[1], (double)pc[2] + t * nn[2]);
                if (v > near_v) near_v = v;
            }
            for (double r = 0.0; r <= QS_MIDLINE_REACH + 1e-9 && face_d < 0.0; r += 0.5) {
                for (int sgn = (r > 0.0 ? -1 : 1); sgn <= 1 && face_d < 0.0; sgn += 2) {
                    double t = sgn * r;
                    double v = qs_ct_at(ctx, (double)pc[0] + t * nn[0],
                                        (double)pc[1] + t * nn[1], (double)pc[2] + t * nn[2]);
                    if (v >= QS_FILL_SNAP_MIN_U8) face_d = r;
                }
            }
            if (face_d < 0.0) suspect = 1;                                   /* no papyrus in reach */
            else if (near_v < QS_DARK_U8 && face_d > QS_MIDLINE_OUTLIER_VOX) suspect = 1;   /* dark and far */
            if (suspect) {
                suspect_map[c] = 1;
                rs->dark_cells++;
                rs->midline_outliers++;
            } else {
                qs_dup_push(&ctx->midline, face_d);
            }
        }
        /* pass 2: an ISOLATED suspect (a dash, a stray row) is refit through
         * its neighbours; a suspect inside a dense suspect region is the fit's
         * systematic depth error and stays (reported) */
        for (size_t c = 0; c < HW; c++) {
            int j = (int)(c / (size_t)W), i = (int)(c % (size_t)W), nsus = 0;
            uint8_t *px = NULL;
            if (!suspect_map[c]) continue;
            for (int dj = -1; dj <= 1; dj++)
                for (int di = -1; di <= 1; di++) {
                    int jj = j + dj, ii = i + di;
                    if ((di == 0 && dj == 0) || jj < 0 || jj >= H || ii < 0 || ii >= W) continue;
                    if (suspect_map[(size_t)jj * (size_t)W + (size_t)ii]) nsus++;
                }
            px = ctx->img + (size_t)j * ctx->ncols + (size_t)(c0 + i - L->col_min);
            if (nsus <= QS_MIDLINE_ISOLATED_MAX) {
                rs->midline_isolated++;
                if (QS_DARK_DEMOTE) {
                    domain[c] = 0; sel[c] = -1; prov[c] = QS_NONE; mat[c] = -1; lane[c] = -1;
                    phase[c] = NAN;
                    field[c * 3] = field[c * 3 + 1] = field[c * 3 + 2] = NAN;
                    *px = QS_IMG_DARKSNAP;
                    rs->dark_snapped++;
                    n_obs--;
                } else {
                    rs->dark_unrecovered++;
                    *px = QS_IMG_MIDLINE_ISO;
                }
            } else {
                rs->midline_dense++;
                *px = QS_IMG_MIDLINE_DENSE;
            }
        }
    }
    rs->observed = n_obs;
    t_build0 = ves_clock_sec();
    if (W >= 2 && H >= 2 && n_obs >= 1) {
        for (int attempt = 0; attempt < 2 && !built; attempt++) {
            QuadStripOpts o;
            Arena_Mark bmark = Arena_save(arena);
            QuadStrip_defaults(&o);
            o.mode = QUAD_STRIP_RECT;
            o.cylindrical_fill = 1;
            o.axis_y = ctx->cfg->axis_y;
            o.axis_x = ctx->cfg->axis_x;
            o.exact_fitted = 1;
            o.max_row_gap = 1;
            o.trend_per_row = 1;
            if (attempt == 1) {
                /* retry without the adopted peel-band cells */
                size_t stripped = 0;
                for (size_t c = 0; c < HW; c++)
                    if (domain[c] && prov[c] == QS_PROV_BAND1) {
                        domain[c] = 0; sel[c] = -1; prov[c] = QS_NONE;
                        mat[c] = -1; lane[c] = -1; stripped++;
                    }
                if (stripped == 0) break;
                n_obs -= stripped;
                rs->observed = n_obs;
                rs->rejected += stripped;
                rs->adopted -= stripped;
                fprintf(stderr, "[solidify] run %d: retrying the build without "
                                "%zu adopted peel-band cells\n", run_index, stripped);
                if (n_obs == 0) break;
            }
            rc = QuadStrip_build_topology_with_phase_ex(
                arena, NULL, domain, field, phase, H, W, 0, 0, 1.0, &o,
                &bverts, &bnv, NULL, NULL, NULL, &bfilled, &bphase);
            if (rc == 0 && bnv == HW) {
                built = attempt == 0 ? 1 : 2;
            } else {
                fprintf(stderr, "[solidify] run %d: build %s (rc=%d nv=%zu of %zu)\n",
                        run_index, attempt == 0 ? "failed" : "failed again",
                        rc, bnv, HW);
                Arena_restore(arena, bmark);
                bverts = NULL; bphase = NULL; bfilled = NULL;
            }
        }
    }
    rs->t_build = ves_clock_sec() - t_build0;
    rs->build = built;
    if (built) ctx->rep->runs_built++;
    if (built == 2) ctx->rep->runs_retried++;
    if (!built) ctx->rep->runs_observed_only++;

    /* positions: observed exact (input floats), fills from the build */
    pos = (float *)ARENA_ALLOC(arena, HW * 3 * sizeof *pos);
    for (size_t c = 0; c < HW; c++) {
        int j = (int)(c / (size_t)W);
        if (sel[c] >= 0) {
            const float *p = m->verts + (size_t)sel[c] * 3;
            pos[c * 3] = p[0]; pos[c * 3 + 1] = p[1]; pos[c * 3 + 2] = p[2];
        } else if (built) {
            pos[c * 3] = L->row_z[j];
            pos[c * 3 + 1] = bverts[c * 3 + 1];
            pos[c * 3 + 2] = bverts[c * 3 + 2];
        } else {
            pos[c * 3] = pos[c * 3 + 1] = pos[c * 3 + 2] = NAN;
        }
    }
    /* CT-RIDGE SNAP of the fills: each fill cell probes the RAW CT along its
     * lattice normal and moves to the nearest papyrus ridge within reach;
     * the displacement field is smoothed over fill cells so neighbouring
     * fills do not pick different ridges.  Observed cells never move. */
    if (built && ctx->zarr != NULL && QS_FILL_SNAP_ENABLE) {
        float *tshift = (float *)ARENA_ALLOC(arena, HW * sizeof *tshift);
        float *tnext = (float *)ARENA_ALLOC(arena, HW * sizeof *tnext);
        float *nrm = (float *)ARENA_ALLOC(arena, HW * 3 * sizeof *nrm);
        uint8_t *cand = (uint8_t *)ARENA_CALLOC(arena, HW, 1);
        for (size_t c = 0; c < HW; c++) {
            int j = (int)(c / (size_t)W), i = (int)(c % (size_t)W);
            const float *pl, *pr, *pd, *pu;
            double eu[3], ev[3], nn[3], len = 0.0;
            tshift[c] = 0.0f;
            nrm[c * 3] = nrm[c * 3 + 1] = nrm[c * 3 + 2] = 0.0f;
            if (sel[c] >= 0 || !isfinite((double)pos[c * 3 + 1])) continue;
            pl = pos + (i > 0 ? c - 1 : c) * 3;
            pr = pos + (i + 1 < W ? c + 1 : c) * 3;
            pd = pos + (j > 0 ? c - (size_t)W : c) * 3;
            pu = pos + (j + 1 < H ? c + (size_t)W : c) * 3;
            if (!isfinite((double)pl[1]) || !isfinite((double)pr[1]) ||
                !isfinite((double)pd[1]) || !isfinite((double)pu[1])) continue;
            for (int a = 0; a < 3; a++) { eu[a] = (double)pr[a] - pl[a]; ev[a] = (double)pu[a] - pd[a]; }
            nn[0] = eu[1] * ev[2] - eu[2] * ev[1];
            nn[1] = eu[2] * ev[0] - eu[0] * ev[2];
            nn[2] = eu[0] * ev[1] - eu[1] * ev[0];
            len = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
            if (!(len > 1e-9)) continue;
            for (int a = 0; a < 3; a++) nrm[c * 3 + a] = (float)(nn[a] / len);
            cand[c] = 1;
            rs->snap_candidates++;
            {
                /* nearest bright voxel along the normal, then the ridge maximum */
                const float *pc = pos + c * 3;
                double best_t = 0.0;
                int found = 0;
                for (double r = 0.0; r <= QS_FILL_SNAP_REACH + 1e-9 && !found; r += 0.5) {
                    for (int sgn = (r > 0.0 ? -1 : 1); sgn <= 1 && !found; sgn += 2) {
                        double t = sgn * r;
                        double v = qs_ct_at(ctx, (double)pc[0] + t * nrm[c * 3],
                                            (double)pc[1] + t * nrm[c * 3 + 1],
                                            (double)pc[2] + t * nrm[c * 3 + 2]);
                        if (v >= QS_FILL_SNAP_MIN_U8) { best_t = t; found = 1; }
                    }
                }
                if (!found) { rs->snap_no_ridge++; continue; }
                {
                    double t0 = best_t, bt = best_t, bv = -1.0;
                    for (double t = t0 - QS_FILL_SNAP_RIDGE_HALF; t <= t0 + QS_FILL_SNAP_RIDGE_HALF + 1e-9; t += 0.5) {
                        double v = qs_ct_at(ctx, (double)pc[0] + t * nrm[c * 3],
                                            (double)pc[1] + t * nrm[c * 3 + 1],
                                            (double)pc[2] + t * nrm[c * 3 + 2]);
                        if (v > bv) { bv = v; bt = t; }
                    }
                    tshift[c] = (float)bt;
                }
            }
        }
        for (int it = 0; it < QS_FILL_SNAP_SMOOTH_ITERS; it++) {
            for (size_t c = 0; c < HW; c++) {
                int j = (int)(c / (size_t)W), i = (int)(c % (size_t)W);
                double acc = 0.0;
                int n = 0;
                tnext[c] = tshift[c];
                if (!cand[c]) continue;
                if (i > 0) { acc += cand[c - 1] ? tshift[c - 1] : 0.0; n++; }
                if (i + 1 < W) { acc += cand[c + 1] ? tshift[c + 1] : 0.0; n++; }
                if (j > 0) { acc += cand[c - (size_t)W] ? tshift[c - (size_t)W] : 0.0; n++; }
                if (j + 1 < H) { acc += cand[c + (size_t)W] ? tshift[c + (size_t)W] : 0.0; n++; }
                tnext[c] = (float)(0.5 * tshift[c] + 0.5 * (n > 0 ? acc / n : 0.0));
            }
            memcpy(tshift, tnext, HW * sizeof *tshift);
        }
        for (size_t c = 0; c < HW; c++) {
            if (!cand[c] || fabs((double)tshift[c]) < 1e-6) continue;
            for (int a = 0; a < 3; a++) pos[c * 3 + a] += tshift[c] * nrm[c * 3 + a];
            rs->snap_moved++;
            qs_dup_push(&ctx->snap, fabs((double)tshift[c]));
        }
    }
    /* BFS from observed cells: fill distance (the cap) and label inheritance */
    dist = (int32_t *)ARENA_ALLOC(arena, HW * sizeof *dist);
    src = (int32_t *)ARENA_ALLOC(arena, HW * sizeof *src);
    queue = (int32_t *)ARENA_ALLOC(arena, HW * sizeof *queue);
    {
        size_t head = 0, tail = 0;
        for (size_t c = 0; c < HW; c++) {
            dist[c] = -1; src[c] = -1;
            if (domain[c]) { dist[c] = 0; src[c] = (int32_t)c; queue[tail++] = (int32_t)c; }
        }
        while (head < tail) {
            size_t c = (size_t)queue[head++];
            int d = dist[c], j = (int)(c / (size_t)W), i = (int)(c % (size_t)W);
            size_t nb[4];
            int nn = 0;
            if (i > 0) nb[nn++] = c - 1;
            if (i + 1 < W) nb[nn++] = c + 1;
            if (j > 0) nb[nn++] = c - (size_t)W;
            if (j + 1 < H) nb[nn++] = c + (size_t)W;
            for (int q = 0; q < nn; q++) {
                size_t n = nb[q];
                if (dist[n] >= 0) continue;
                dist[n] = d + 1;
                src[n] = src[c];
                queue[tail++] = (int32_t)n;
            }
        }
    }
    /* Fill domain (user rule 2026-09-01): the V domain is the box's z extent
     * and every wrap spans it, so for any u-range bounded on both sides by
     * GOOD columns (>= QS_GOOD_COLUMN_FRAC of the rows observed) EVERY row is
     * filled -- interpolation along u between the bounds, continuation along v
     * to the box faces (the box and core gates below still apply).  A stretch
     * of non-good columns wider than QS_FILL_MAX_U_GAP stays empty, and so
     * does everything outside the first..last good column: extrapolation in u
     * is never interpolation.  A gap is REFUSED when its two sides sit on
     * different wraps -- the radius extrapolated across the gap from each
     * side's local slope must land within QS_SKIP_DR_PITCH pitches of the
     * other side, row by row -- or carry two material lineages.  Inside the
     * domain a row bridge (rims <= QS_BRIDGE_MAX_CELLS apart) that fails the
     * fit's chord gates vetoes its row's cells: the harmonic solve would
     * otherwise blend a foreign rim in. */
    {
        enum { QS_COL_OUT = 0, QS_COL_FILL = 1, QS_COL_WIDE = 2, QS_COL_SKIP = 3,
               QS_COL_MAT = 4, QS_COL_ORDER = 5 };
        double dtheta_col = 0.0;   /* the run's median angle step per column */
        const double good_frac = ctx->cfg->good_column_frac > 0.0 ? ctx->cfg->good_column_frac
                                                                  : QS_GOOD_COLUMN_FRAC;
        const int max_gap = ctx->cfg->fill_max_u_gap > 0 ? ctx->cfg->fill_max_u_gap
                                                           : QS_FILL_MAX_U_GAP;
        const double skip_tol = ctx->cfg->skip_dr_pitch > 0.0 ? ctx->cfg->skip_dr_pitch
                                                               : QS_SKIP_DR_PITCH;
        int32_t *ncol = (int32_t *)ARENA_CALLOC(arena, (size_t)W, sizeof(int32_t));
        int32_t *colmat = (int32_t *)ARENA_ALLOC(arena, (size_t)W * sizeof(int32_t));
        uint8_t *good = (uint8_t *)ARENA_CALLOC(arena, (size_t)W, 1);
        uint8_t *colstate = (uint8_t *)ARENA_CALLOC(arena, (size_t)W, 1);
        uint8_t *rowveto = (uint8_t *)ARENA_CALLOC(arena, HW, 1);
        int good_min = (int)ceil(good_frac * (double)H);
        int first_good = -1, last_good = -1, prev_good = -1;
        if (good_min < 1) good_min = 1;
        for (int i = 0; i < W; i++) {
            /* observed count and dominant material per column */
            int32_t m0 = -1, m1 = -1;
            size_t n0 = 0, n1 = 0;
            for (int j = 0; j < H; j++) {
                size_t c = (size_t)j * (size_t)W + (size_t)i;
                if (!domain[c]) continue;
                ncol[i]++;
                if (mat[c] < 0) continue;
                if (m0 < 0 || mat[c] == m0) { m0 = mat[c]; n0++; }
                else if (m1 < 0 || mat[c] == m1) { m1 = mat[c]; n1++; }
            }
            colmat[i] = n1 > n0 ? m1 : m0;
            good[i] = ncol[i] >= good_min;
            if (good[i]) { if (first_good < 0) first_good = i; last_good = i; }
        }
        {
            /* median lifted-phase step between adjacent observed cells: u is
             * arc length, so this is the run's angle per column */
            double *dd = (double *)ARENA_ALLOC(arena, HW * sizeof *dd);
            size_t nd = 0;
            for (size_t c = 0; c + 1 < HW; c++) {
                if ((c + 1) % (size_t)W == 0) continue;
                if (!domain[c] || !domain[c + 1]) continue;
                dd[nd++] = (double)phase[c + 1] - (double)phase[c];
            }
            if (nd > 0) {
                qsort(dd, nd, sizeof *dd, qs_cmp_double);
                dtheta_col = dd[nd / 2];
            }
        }
        for (int i = 0; i < W; i++) {
            if (!good[i]) continue;
            colstate[i] = QS_COL_FILL;
            if (prev_good >= 0 && i - prev_good > 1) {
                int gL = prev_good, gR = i, gap = gR - gL - 1;
                uint8_t state = QS_COL_FILL;
                if (gap > max_gap) {
                    state = QS_COL_WIDE;
                } else if (colmat[gL] >= 0 && colmat[gR] >= 0 && colmat[gL] != colmat[gR]) {
                    state = QS_COL_MAT;
                } else {
                    /* wrap-skip test row by row across the gap; a gap the two
                     * sides share too few rows on cannot be verified */
                    size_t comparable = 0, failing = 0, misordered = 0;
                    size_t min_rows = (size_t)ceil(QS_SKIP_MIN_ROWS_FRAC * (double)H);
                    if (min_rows < 4) min_rows = 4;
                    for (int j = 0; j < H; j++) {
                        size_t row = (size_t)j * (size_t)W;
                        double want = 0.0, have = 0.0;
                        if (!domain[row + (size_t)gL] || !domain[row + (size_t)gR]) continue;
                        comparable++;
                        want = qs_want_dtheta(ctx, L, pos, row + (size_t)gL,
                                              row + (size_t)gR, gR - gL, dtheta_col);
                        have = (double)phase[row + (size_t)gR] - (double)phase[row + (size_t)gL];
                        if (fabs(have - want) > QS_ORDER_TOL_FRAC * fabs(want) + QS_ORDER_TOL_RAD)
                            misordered++;
                        else if (!qs_bridge_consistent(ctx, domain, pos, row + (size_t)gL,
                                                       row + (size_t)gR, 1, gR - gL, gL, gR, W,
                                                       QS_BRIDGE_MAX_CELLS, pitch, skip_tol))
                            failing++;
                    }
                    if (misordered * 4 > comparable) state = QS_COL_ORDER;
                    else if (comparable < min_rows || failing * 4 > comparable) state = QS_COL_SKIP;
                }
                for (int k = gL + 1; k < gR; k++) colstate[k] = state;
            }
            prev_good = i;
        }
        /* bridge vetoes inside the domain: every run of non-observed cells
         * between two observed rims at most QS_FILL_MAX_U_GAP apart, along
         * rows AND along columns, must join one material lineage, pass the
         * fit's stretch gate when short (<= QS_BRIDGE_MAX_CELLS), and pass the
         * wrap-skip test; a failing bridge vetoes its cells (1 chord, 2
         * material, 3 skip) */
        for (int axis = 0; axis < 2; axis++) {
            int outer_n = axis == 0 ? H : W, inner_n = axis == 0 ? W : H;
            size_t stride = axis == 0 ? 1u : (size_t)W;
            double cell = axis == 0 ? L->du : L->dv;
            for (int o = 0; o < outer_n; o++) {
                int t = 0;
                while (t < inner_n) {
                    size_t c = axis == 0 ? (size_t)o * (size_t)W + (size_t)t
                                         : (size_t)t * (size_t)W + (size_t)o;
                    int t0 = t, t1 = t;
                    if (domain[c]) { t++; continue; }
                    while (t1 + 1 < inner_n && !domain[c + (size_t)(t1 + 1 - t0) * stride]) t1++;
                    if (t0 > 0 && t1 + 1 < inner_n) {
                        size_t a = c - stride, b = c + (size_t)(t1 + 1 - t0) * stride;
                        uint8_t v = 0;
                        if (!qs_mat_ok(mat[a], mat[b], -1)) {
                            v = 2;
                        } else if (axis == 0) {
                            /* angular order: the bridge must turn by its span
                             * times the run's angle per column */
                            double want = qs_want_dtheta(ctx, L, pos, a, b, t1 - t0 + 2, dtheta_col);
                            double have = (double)phase[b] - (double)phase[a];
                            if (fabs(have - want) > QS_ORDER_TOL_FRAC * fabs(want) + QS_ORDER_TOL_RAD)
                                v = 4;
                        }
                        if (v == 0 && t1 - t0 + 1 <= QS_BRIDGE_MAX_CELLS) {
                            const float *pa = pos + a * 3, *pb = pos + b * 3;
                            double chord = sqrt(pow((double)pa[0] - pb[0], 2) +
                                                pow((double)pa[1] - pb[1], 2) +
                                                pow((double)pa[2] - pb[2], 2));
                            double arc = (double)(t1 - t0 + 2) * cell;
                            /* the stretch reference is what the neighbouring
                             * observed lines actually span over the same cells:
                             * a uniformly stretched sheet is not a foreign rim */
                            double ref = 1.0, rsum = 0.0;
                            int rn = 0;
                            size_t side = axis == 0 ? (size_t)W : 1u;
                            for (int sgn = -1; sgn <= 1; sgn += 2) {
                                size_t a2 = 0, b2 = 0;
                                if (sgn < 0 && o == 0) continue;
                                if (sgn > 0 && o + 1 >= outer_n) continue;
                                a2 = sgn < 0 ? a - side : a + side;
                                b2 = sgn < 0 ? b - side : b + side;
                                if (!domain[a2] || !domain[b2]) continue;
                                {
                                    const float *qa = pos + a2 * 3, *qb = pos + b2 * 3;
                                    double ch2 = sqrt(pow((double)qa[0] - qb[0], 2) +
                                                      pow((double)qa[1] - qb[1], 2) +
                                                      pow((double)qa[2] - qb[2], 2));
                                    rsum += ch2 / arc;
                                    rn++;
                                }
                            }
                            if (rn > 0 && rsum / rn > 1.0) ref = rsum / rn;
                            if (chord > QS_FILL_STRETCH * ref * arc ||
                                chord > QS_FILL_CHORD_ABS * ref) v = 1;
                        }
                        if (v == 0 && t1 - t0 + 1 <= max_gap &&
                            !qs_bridge_consistent(ctx, domain, pos, a, b, stride,
                                                  t1 - t0 + 2, t0 - 1, t1 + 1,
                                                  inner_n, QS_BRIDGE_MAX_CELLS, pitch,
                                                  skip_tol))
                            v = 3;
                        if (v)
                            for (int k = t0; k <= t1; k++) {
                                size_t x = c + (size_t)(k - t0) * stride;
                                if (rowveto[x] == 0 || v < rowveto[x]) rowveto[x] = v;
                            }
                    }
                    t = t1 + 1;
                }
            }
        }
        /* HOLE-LEVEL VOTE (user review 2026-09-02: comb teeth at hole rims).
         * Bridge vetoes are per row / column, so a hole whose rows alternate
         * between passing and failing bridges was baked as a comb.  Every
         * 4-connected hole inside the fill domain is now decided as a whole:
         * refused when >= QS_HOLE_VETO_FRAC of its cells carry any veto or
         * >= QS_HOLE_FUSION_FRAC carry a fusion-class veto (material, wrap
         * skip, angular order); filled entirely otherwise. */
        {
            int32_t *hole = (int32_t *)ARENA_ALLOC(arena, HW * sizeof *hole);
            int32_t *hq = (int32_t *)ARENA_ALLOC(arena, HW * sizeof *hq);
            for (size_t c = 0; c < HW; c++) hole[c] = -1;
            for (size_t c = 0; c < HW; c++) {
                size_t head = 0, tail = 0, n = 0, n_any = 0, n_fusion = 0;
                size_t cnt[5] = {0, 0, 0, 0, 0};
                int refuse = 0;
                uint8_t code = 0;
                if (domain[c] || hole[c] >= 0 || colstate[c % (size_t)W] != QS_COL_FILL) continue;
                hole[c] = (int32_t)c;
                hq[tail++] = (int32_t)c;
                while (head < tail) {
                    size_t x = (size_t)hq[head++];
                    int jj = (int)(x / (size_t)W), ii = (int)(x % (size_t)W);
                    size_t nb[4];
                    int nn = 0;
                    n++;
                    if (rowveto[x]) { n_any++; cnt[rowveto[x]]++; }
                    if (rowveto[x] == 2 || rowveto[x] == 3 || rowveto[x] == 4) n_fusion++;
                    if (ii > 0) nb[nn++] = x - 1;
                    if (ii + 1 < W) nb[nn++] = x + 1;
                    if (jj > 0) nb[nn++] = x - (size_t)W;
                    if (jj + 1 < H) nb[nn++] = x + (size_t)W;
                    for (int q = 0; q < nn; q++) {
                        size_t y = nb[q];
                        if (domain[y] || hole[y] >= 0 || colstate[y % (size_t)W] != QS_COL_FILL) continue;
                        hole[y] = (int32_t)c;
                        hq[tail++] = (int32_t)y;
                    }
                }
                rs->holes++;
                rs->hole_cells += n;
                refuse = (double)n_fusion >= QS_HOLE_FUSION_FRAC * (double)n ||
                         (double)n_any >= QS_HOLE_VETO_FRAC * (double)n;
                if (refuse) {
                    rs->refused_hole_cells += n;
                    /* the dominant veto class labels the whole hole; fusion first */
                    code = cnt[2] >= cnt[3] && cnt[2] >= cnt[4] && cnt[2] > 0 ? 2
                         : (cnt[3] >= cnt[4] && cnt[3] > 0 ? 3 : (cnt[4] > 0 ? 4 : 1));
                    rs->holes_refused++;
                } else if (n_any > 0) {
                    rs->holes_filled++;
                }
                for (size_t k = 0; k < tail; k++) {
                    size_t x = (size_t)hq[k];
                    rowveto[x] = (refuse && !(folded != NULL && folded[x])) ? code : 0;
                }
            }
        }
        keep = (uint8_t *)ARENA_CALLOC(arena, HW, 1);
        for (size_t c = 0; c < HW; c++) {
            int j = (int)(c / (size_t)W), i = (int)(c % (size_t)W);
            uint8_t *px = ctx->img + (size_t)j * ctx->ncols + (size_t)(c0 + i - L->col_min);
            if (domain[c]) { keep[c] = 1; continue; }
            if (!built) continue;
            if (colstate[i] == QS_COL_OUT || colstate[i] == QS_COL_WIDE) {
                rs->dropped_unbridged++;
                if (*px == QS_IMG_EMPTY) *px = QS_IMG_DROPPED;
            } else if (colstate[i] == QS_COL_SKIP) {
                rs->dropped_skip++;
                if (*px == QS_IMG_EMPTY) *px = QS_IMG_REFUSED;
            } else if (colstate[i] == QS_COL_ORDER || rowveto[c] == 4) {
                rs->dropped_order++;
                if (*px == QS_IMG_EMPTY) *px = QS_IMG_REFUSED;
            } else if (colstate[i] == QS_COL_MAT || rowveto[c] == 2) {
                rs->dropped_material++;
                if (*px == QS_IMG_EMPTY) *px = QS_IMG_REFUSED;
            } else if (rowveto[c] == 1) {
                rs->dropped_chord++;
                if (*px == QS_IMG_EMPTY) *px = QS_IMG_REFUSED;
            } else if (rowveto[c] == 3) {
                rs->dropped_skip++;
                if (*px == QS_IMG_EMPTY) *px = QS_IMG_REFUSED;
            } else {
                keep[c] = 2;   /* candidate: box and core gates below */
            }
        }
        (void)first_good; (void)last_good;
    }
    for (size_t c = 0; c < HW; c++) {
        if (keep[c] != 2) continue;
        keep[c] = 0;
        {
            const float *p = pos + c * 3;
            int inside = 1;
            for (int a = 0; a < 3; a++)
                if ((double)p[a] < in->bbox_lo[a] - QS_BBOX_MARGIN_VOX ||
                    (double)p[a] > in->bbox_hi[a] + QS_BBOX_MARGIN_VOX)
                    inside = 0;
            if (!inside) { rs->dropped_bbox++; continue; }
            if (qs_radius(ctx, p) < QS_CORE_RADIUS_VOX) { rs->dropped_core++; continue; }
        }
        keep[c] = 1;
        prov[c] = QS_PROV_FILL;
        mat[c] = mat[(size_t)src[c]];
        lane[c] = lane[(size_t)src[c]];
        phase[c] = bphase != NULL && isfinite((double)bphase[c]) ? bphase[c]
                                                                 : phase[(size_t)src[c]];
    }
    /* triangles: per cell with four kept corners, the split with more legal
     * triangles wins (tie -> main diagonal, the writer's convention) */
    used = (uint8_t *)ARENA_CALLOC(arena, HW, 1);
    tri = (int32_t *)ARENA_ALLOC(arena, (HW * 2 + 1) * 3 * sizeof *tri);
    for (int j = 0; j + 1 < H; j++) {
        for (int i = 0; i + 1 < W; i++) {
            size_t a = (size_t)j * (size_t)W + (size_t)i, b = a + 1;
            size_t cc = a + (size_t)W, d = cc + 1;
            int m1 = 0, m2 = 0, a1 = 0, a2 = 0, nm = 0, na = 0;
            if (!keep[a] || !keep[b] || !keep[cc] || !keep[d]) continue;
            m1 = qs_tri_legal(pos, mat, a, b, cc);
            m2 = qs_tri_legal(pos, mat, b, d, cc);
            a1 = qs_tri_legal(pos, mat, a, b, d);
            a2 = qs_tri_legal(pos, mat, a, d, cc);
            nm = (m1 == 1) + (m2 == 1);
            na = (a1 == 1) + (a2 == 1);
            if (nm >= na) {
                if (m1 == 1) { tri[ntri*3] = (int32_t)a; tri[ntri*3+1] = (int32_t)b; tri[ntri*3+2] = (int32_t)cc; ntri++; }
                else if (m1 == 0) rs->skipped_gate++; else rs->skipped_material++;
                if (m2 == 1) { tri[ntri*3] = (int32_t)b; tri[ntri*3+1] = (int32_t)d; tri[ntri*3+2] = (int32_t)cc; ntri++; }
                else if (m2 == 0) rs->skipped_gate++; else rs->skipped_material++;
            } else {
                if (a1 == 1) { tri[ntri*3] = (int32_t)a; tri[ntri*3+1] = (int32_t)b; tri[ntri*3+2] = (int32_t)d; ntri++; }
                else if (a1 == 0) rs->skipped_gate++; else rs->skipped_material++;
                if (a2 == 1) { tri[ntri*3] = (int32_t)a; tri[ntri*3+1] = (int32_t)d; tri[ntri*3+2] = (int32_t)cc; ntri++; }
                else if (a2 == 0) rs->skipped_gate++; else rs->skipped_material++;
            }
        }
    }
    /* RIDGE TRACK: the lattice is built; now put it ON the papyrus.  Every cell
     * of the run (observed and filled alike) picks an offset along its own
     * normal that lands on a CT ridge, stays near the fit and agrees with its
     * neighbours.  Off by constant until measured. */
    if (QS_RIDGE_SNAP_ENABLE && ctx->zarr != NULL && built) {
        RidgeTrackOpts ro;
        RidgeTrackReport rr;
        uint8_t *cellvalid = (uint8_t *)ARENA_ALLOC(arena, HW);
        uint8_t *cellwide = (uint8_t *)ARENA_ALLOC(arena, HW);
        RidgeTrack_defaults(&ro);
        if (cellvalid != NULL && cellwide != NULL) {
            /* fills have no observation to stay near: they search the wide reach */
            for (size_t c = 0; c < HW; c++) {
                cellvalid[c] = (uint8_t)(keep[c] ? 1 : 0);
                cellwide[c] = (uint8_t)((keep[c] && prov[c] == QS_PROV_FILL) ? 1 : 0);
            }
            /* NULL to_world: the verdict measures the STRAIGHTENED file (measured
             * 2026-09-03 -- the baseline's straightened max edge is 5.999 with 0
             * over the limit, while its world twin already holds 143 over from the
             * straightening shear), so the clamp must judge edges in the lattice
             * frame, which is the frame `pos` is already in. */
            double *celldepth = (double *)ARENA_ALLOC(arena, HW * sizeof *celldepth);
            if (celldepth != NULL)
                for (size_t c = 0; c < HW; c++)
                    celldepth[c] = keep[c] ? qs_radius(ctx, pos + c * 3) : 0.0;
            if (RidgeTrack_solve(arena, &ro, qs_ridge_sample, NULL, (void *)ctx,
                                 pos, cellvalid, cellwide, W, H, tri, ntri, celldepth, &rr) == 0) {
                rs->ridge_cells = rr.cells;
                rs->ridge_moved = rr.moved;
                rs->ridge_on_before = rr.on_ridge_before;
                rs->ridge_on_after = rr.on_ridge;
                rs->ridge_shift_p50 = rr.shift_p50;
                rs->ridge_shift_p90 = rr.shift_p90;
                rs->ridge_shift_max = rr.shift_max;
                rs->ridge_fill_cells = rr.fill_cells;
                rs->ridge_fill_on_before = rr.fill_on_before;
                rs->ridge_fill_on_after = rr.fill_on_after;
                qs_dup_push(&ctx->ridge, rr.shift_p50);
                fprintf(stderr, "      ridge track: %zu cells, %zu moved (|shift| p50 %.2f p90 %.2f max %.2f vox), "
                        "on papyrus %zu -> %zu (%.1f%% -> %.1f%%), CT at the cell p50 %.0f -> %.0f, "
                        "max neighbour disagreement %.2f vox, %zu edge clamp(s), longest edge %.2f vox, "
                        "depth deviation from the local median p90 %.2f -> %.2f vox; fills %zu on papyrus %zu -> %zu\n",
                        rr.cells, rr.moved, rr.shift_p50, rr.shift_p90, rr.shift_max,
                        rr.on_ridge_before, rr.on_ridge,
                        rr.cells ? 100.0 * (double)rr.on_ridge_before / (double)rr.cells : 0.0,
                        rr.cells ? 100.0 * (double)rr.on_ridge / (double)rr.cells : 0.0,
                        rr.ct_before_p50, rr.ct_after_p50, rr.neighbour_max,
                        rr.edge_clamped, rr.edge_max, rr.hp_p90_before, rr.hp_p90_after,
                        rr.fill_cells, rr.fill_on_before, rr.fill_on_after);
            }
        }
    }

    /* DEPTH SMOOTHING.  What is left of the depth seams after everything else
     * joins two real observations: the extracted surface's own roughness at the
     * 2-cell scale, which is what the gate's high-pass measures.  A few bounded
     * Jacobi passes on the radius field remove it.  The bound (well inside the
     * verdict's 3 vox coverage radius) keeps every cell on the source vertex it
     * came from, and cells move along their own lattice normal so the surface
     * stays where it is in u and v. */
    if (QS_OBS_SMOOTH_VOX > 0.0) {
        double *rr0 = (double *)ARENA_ALLOC(arena, HW * sizeof *rr0);
        double *rr1 = (double *)ARENA_ALLOC(arena, HW * sizeof *rr1);
        double *shift = (double *)ARENA_ALLOC(arena, HW * sizeof *shift);
        if (rr0 != NULL && rr1 != NULL && shift != NULL) {
            size_t nsh = 0;
            for (size_t c = 0; c < HW; c++)
                rr0[c] = keep[c] ? qs_radius(ctx, pos + c * 3) : NAN;
            for (int it = 0; it < QS_OBS_SMOOTH_ITERS; it++) {
                for (size_t c = 0; c < HW; c++) rr1[c] = rr0[c];
                for (int j = 0; j < H; j++)
                    for (int i = 0; i < W; i++) {
                        size_t c = (size_t)j * (size_t)W + (size_t)i;
                        double acc = 0.0;
                        int n = 0;
                        if (!keep[c] || !isfinite(rr0[c])) continue;
                        if (i > 0 && keep[c - 1] && isfinite(rr0[c - 1])) { acc += rr0[c - 1]; n++; }
                        if (i + 1 < W && keep[c + 1] && isfinite(rr0[c + 1])) { acc += rr0[c + 1]; n++; }
                        if (j > 0 && keep[c - (size_t)W] && isfinite(rr0[c - (size_t)W])) {
                            acc += rr0[c - (size_t)W]; n++;
                        }
                        if (j + 1 < H && keep[c + (size_t)W] && isfinite(rr0[c + (size_t)W])) {
                            acc += rr0[c + (size_t)W]; n++;
                        }
                        if (n == 0) continue;
                        rr1[c] = rr0[c] + QS_OBS_SMOOTH_LAMBDA * (acc / (double)n - rr0[c]);
                    }
                for (size_t c = 0; c < HW; c++) rr0[c] = rr1[c];
            }
            /* apply as a radial move, bounded */
            for (size_t c = 0; c < HW; c++) {
                double r_old = 0.0, dr = 0.0, dy = 0.0, dx = 0.0, len = 0.0;
                if (!keep[c] || !isfinite(rr0[c])) continue;
                r_old = qs_radius(ctx, pos + c * 3);
                dr = rr0[c] - r_old;
                if (dr > QS_OBS_SMOOTH_VOX) dr = QS_OBS_SMOOTH_VOX;
                if (dr < -QS_OBS_SMOOTH_VOX) dr = -QS_OBS_SMOOTH_VOX;
                if (fabs(dr) < 1e-6) continue;
                dy = (double)pos[c * 3 + 1] - ctx->cfg->axis_y;
                dx = (double)pos[c * 3 + 2] - ctx->cfg->axis_x;
                len = sqrt(dy * dy + dx * dx);
                if (!(len > 1e-9)) continue;
                pos[c * 3 + 1] = (float)((double)pos[c * 3 + 1] + dr * dy / len);
                pos[c * 3 + 2] = (float)((double)pos[c * 3 + 2] + dr * dx / len);
                shift[nsh++] = fabs(dr);
                rs->obs_smooth_cells++;
            }
            if (nsh > 0) {
                qsort(shift, nsh, sizeof *shift, qs_cmp_double);
                qs_dup_push(&ctx->obs_smooth, shift[nsh / 2]);
            }
        }
    }

    /* FILL CAP: a fill cell too far from any observation is invention, not
     * interpolation.  Dropping it here (before the faces are counted) leaves an
     * honest hole instead of a smeared patch.  Off by default. */
    if (QS_FILL_MAX_DIST_CELLS > 0) {
        size_t ntri_keep = 0;
        for (size_t c = 0; c < HW; c++)
            if (!domain[c] && keep[c] && dist[c] > QS_FILL_MAX_DIST_CELLS) {
                keep[c] = 0;
                rs->dropped_fill_far++;
            }
        for (size_t t = 0; t < ntri; t++) {
            size_t a = (size_t)tri[t * 3], b = (size_t)tri[t * 3 + 1], c2 = (size_t)tri[t * 3 + 2];
            if (!keep[a] || !keep[b] || !keep[c2]) continue;
            tri[ntri_keep * 3] = (int32_t)a;
            tri[ntri_keep * 3 + 1] = (int32_t)b;
            tri[ntri_keep * 3 + 2] = (int32_t)c2;
            ntri_keep++;
        }
        ntri = ntri_keep;
    }
    for (size_t t = 0; t < ntri * 3; t++) used[(size_t)tri[t]] = 1;
    /* final keep: observations always, fills only when a face references them */
    oidx = (int32_t *)ARENA_ALLOC(arena, HW * sizeof *oidx);
    {
        size_t nv_before = ctx->out.nv;
        long base = ctx->col_out;
    for (size_t c = 0; c < HW; c++) {
        int j = (int)(c / (size_t)W), i = (int)(c % (size_t)W);
        uint8_t *px = ctx->img + (size_t)j * ctx->ncols + (size_t)(c0 + i - L->col_min);
        oidx[c] = -1;
        if (!keep[c]) {
            if (!domain[c] && built && *px == QS_IMG_EMPTY) *px = QS_IMG_DROPPED;
            continue;
        }
        if (!domain[c] && !used[c]) {
            rs->dropped_unref++;
            if (*px == QS_IMG_EMPTY) *px = QS_IMG_DROPPED;
            continue;
        }
        if (ctx->out.nv >= ctx->out.cap_v) { Arena_restore(arena, mark); return -1; }
        oidx[c] = (int32_t)ctx->out.nv;
        {
            size_t o = ctx->out.nv++;
            QsOut *out = &ctx->out;
            out->verts[o * 3] = pos[c * 3];
            out->verts[o * 3 + 1] = pos[c * 3 + 1];
            out->verts[o * 3 + 2] = pos[c * 3 + 2];
            out->uv[o * 2] = (float)((double)(base + i) * L->du);   /* LAYOUT */
            out->uv[o * 2 + 1] = (float)((double)j * L->dv);
            out->support[o] = domain[c] ? in->support[(size_t)sel[c]] : 0;
            out->prov[o] = prov[c];
            out->phase[o] = phase[c];
            out->material[o] = mat[c];
            out->lane[o] = lane[c];
            out->run[o] = run_index;
            out->atlas_col[o] = (int32_t)(c0 + i);
            if (domain[c]) {
                if (*px != QS_IMG_DUPLICATE && *px != QS_IMG_MIDLINE_ISO &&
                    *px != QS_IMG_MIDLINE_DENSE)
                    *px = prov[c] == QS_PROV_BAND1 ? QS_IMG_BAND1
                        : (prov[c] == QS_PROV_OBSERVED_GENERATED ? QS_IMG_GENERATED
                        : (prov[c] == QS_PROV_OBSERVED_SNAPPED ? QS_IMG_DARKSNAP
                                                                 : QS_IMG_FIXED));
            } else {
                *px = QS_IMG_FILL;
                rs->filled++;
                /* FILL-SUPPORT: how far is this invented cell from real data? */
                rs->fill_cells++;
                if (dist[c] >= 0) {
                    qs_dup_push(&ctx->fill_dist, (double)dist[c]);
                    if ((size_t)dist[c] > rs->fill_dist_max)
                        rs->fill_dist_max = (size_t)dist[c];
                    if (dist[c] > QS_FILL_SUPPORT_CELLS) rs->fill_far++;
                }
            }
        }
    }
        if (ctx->out.nv > nv_before) {
            rs->c0_out = base;
            ctx->col_out = base + (long)W + QS_RUN_GAP_COLS;
        }
    }
    if (ctx->out.nf + ntri > ctx->out.cap_f) { Arena_restore(arena, mark); return -1; }
    for (size_t t = 0; t < ntri; t++) {
        for (int e = 0; e < 3; e++)
            ctx->out.faces[(ctx->out.nf + t) * 3 + (size_t)e] =
                oidx[(size_t)tri[t * 3 + (size_t)e]];
    }
    /* SCALE gate: a lattice triangle should cover as much 3-D area as it covers
     * (u,v) area -- u is arc length and v is z, so the ribbon is metric by
     * construction.  A triangle that covers twice its uv area smears the CT it
     * bakes; one that covers half is crumpled.  Anisotropy (the verdict's gate)
     * cannot see either. */
    for (size_t t = 0; t < ntri; t++) {
        size_t ca = (size_t)tri[t * 3], cb = (size_t)tri[t * 3 + 1], cc2 = (size_t)tri[t * 3 + 2];
        const float *pa = pos + ca * 3, *pb = pos + cb * 3, *pc2 = pos + cc2 * 3;
        double e1[3], e2[3], cx, cy, cz, a3, a2, l2;
        double ua = (double)(ca % (size_t)W) * L->du, va = (double)(ca / (size_t)W) * L->dv;
        double ub = (double)(cb % (size_t)W) * L->du, vb = (double)(cb / (size_t)W) * L->dv;
        double uc = (double)(cc2 % (size_t)W) * L->du, vc = (double)(cc2 / (size_t)W) * L->dv;
        int is_fill = (!domain[ca] || !domain[cb] || !domain[cc2]);
        a2 = 0.5 * fabs((ub - ua) * (vc - va) - (uc - ua) * (vb - va));
        if (!(a2 > 1e-9)) continue;
        e1[0] = (double)pb[0] - pa[0]; e1[1] = (double)pb[1] - pa[1]; e1[2] = (double)pb[2] - pa[2];
        e2[0] = (double)pc2[0] - pa[0]; e2[1] = (double)pc2[1] - pa[1]; e2[2] = (double)pc2[2] - pa[2];
        cx = e1[1] * e2[2] - e1[2] * e2[1];
        cy = e1[2] * e2[0] - e1[0] * e2[2];
        cz = e1[0] * e2[1] - e1[1] * e2[0];
        a3 = 0.5 * sqrt(cx * cx + cy * cy + cz * cz);
        if (!(a3 > 1e-12)) continue;
        l2 = fabs(log2(a3 / a2));
        rs->scale_faces++;
        qs_dup_push(&ctx->scale, l2);
        qs_dup_push(is_fill ? &ctx->scale_fill : &ctx->scale_obs, l2);
        if (l2 > QS_GATE_SCALE_LOG2) {
            rs->scale_over++;
            if (is_fill) rs->scale_over_fill++;
        }
    }
    ctx->out.nf += ntri;
    rs->faces = ntri;
    fprintf(stderr,
            "[solidify] run %d/%d cols [%ld,%ld] w=%d: observed=%zu adopted=%zu "
            "rejected=%zu no-ref=%zu dup=%zu filled=%zu not-filled(unbridged/skip/order/chord/"
            "material/bbox/core/unref)=%zu/%zu/%zu/%zu/%zu/%zu/%zu/%zu holes=%zu (refused %zu, "
            "filled-by-vote %zu) fold=%zu cells/%zu cols dark=%zu (snapped %zu) snap=%zu/%zu (no ridge %zu) faces=%zu skipped(gate/mat)=%zu/%zu scale-over=%zu/%zu fill-far=%zu/%zu (max %zu) ramp=%.2f pitches/turn over %.2f turns build=%s out-cols [%ld,%ld] %.2fs\n",
            run_index + 1, ctx->n_runs, c0, c1, W, rs->observed, rs->adopted,
            rs->rejected, rs->no_ref, rs->duplicates, rs->filled, rs->dropped_unbridged,
            rs->dropped_skip, rs->dropped_order, rs->dropped_chord, rs->dropped_material,
            rs->dropped_bbox, rs->dropped_core, rs->dropped_unref, rs->holes, rs->holes_refused,
            rs->holes_filled, rs->fold_cells, rs->fold_columns, rs->dark_cells, rs->dark_snapped,
            rs->snap_moved, rs->snap_candidates, rs->snap_no_ridge,
            rs->faces, rs->skipped_gate, rs->skipped_material,
            rs->scale_over, rs->scale_faces, rs->fill_far, rs->fill_cells, rs->fill_dist_max,
            rs->ramp_pitches_per_turn, rs->ramp_turns,
            built == 0 ? "observed-only" : (built == 1 ? "ok" : "retry-ok"),
            rs->c0_out, rs->c0_out >= 0 ? rs->c0_out + (long)W - 1 : -1L, rs->t_build);
    Arena_restore(arena, mark);
    return 0;
}

/* ---- outputs -------------------------------------------------------------- */

static int qs_write_provenance_png(const QsCtx *ctx, const char *path)
{
    size_t w = ctx->ncols, h = (size_t)ctx->L->n_slices;
    uint8_t *rgb = (uint8_t *)malloc(w * h * 3);
    int rc = -1;
    static const uint8_t palette[QS_IMG_COUNT][3] = {
        { 0, 0, 0 },        /* empty */
        { 200, 200, 200 },  /* observed fixed */
        { 120, 120, 120 },  /* observed generated */
        { 80, 120, 255 },   /* band-1 adopted */
        { 60, 200, 60 },    /* fill */
        { 110, 20, 20 },    /* dropped (cap / bbox / core / unreferenced) */
        { 230, 150, 40 },   /* duplicate cell (band 0 kept, band 1 conflicts) */
        { 200, 60, 200 },   /* band-1 rejected (no band 0 at the cell) */
        { 120, 40, 160 },   /* bridge refused: chord / material / wrap-skip */
        { 70, 90, 120 },    /* confetti run: observed but not emitted */
        { 255, 140, 0 },    /* FOLD: compressed column demoted to a hole */
        { 40, 200, 220 },   /* dark observed cell snapped onto the CT ridge */
        { 255, 0, 160 },    /* MIDLINE isolated suspect, kept in place */
        { 255, 230, 0 }     /* MIDLINE dense suspect: off the papyrus, kept in place */
    };
    if (rgb == NULL) return -1;
    for (size_t i = 0; i < w * h; i++) {
        const uint8_t *p = palette[ctx->img[i] < QS_IMG_COUNT ? ctx->img[i] : 0];
        rgb[i * 3] = p[0]; rgb[i * 3 + 1] = p[1]; rgb[i * 3 + 2] = p[2];
    }
    rc = VesPng_write_rgb(path, rgb, (int)w, (int)h);
    free(rgb);
    return rc;
}

/* JSON string escape for paths (backslashes, quotes, control characters) */
static void qs_json_escape(const char *in, char *out, size_t cap)
{
    size_t o = 0;
    if (cap == 0) return;
    for (; *in != '\0' && o + 2 < cap; in++) {
        unsigned char ch = (unsigned char)*in;
        if (ch == '"' || ch == '\\') { out[o++] = '\\'; out[o++] = (char)ch; }
        else if (ch < 0x20) out[o++] = ' ';
        else out[o++] = (char)ch;
    }
    out[o] = '\0';
}

static int qs_write_report(const QsCtx *ctx, const char *path,
                           const char *in_vmesh, const char *out_stem)
{
    const QsReport *r = ctx->rep;
    char in_esc[1024] = { 0 }, out_esc[1024] = { 0 };
    FILE *f = fopen(path, "wb");
    if (f == NULL) return -1;
    qs_json_escape(in_vmesh, in_esc, sizeof in_esc);
    qs_json_escape(out_stem, out_esc, sizeof out_esc);
    fprintf(f, "{\n  \"schema\": \"vesuvius-quadribbon-solid-v1\",\n");
    fprintf(f, "  \"input\": \"%s\",\n  \"output\": \"%s.vmesh\",\n", in_esc, out_esc);
    fprintf(f, "  \"lattice\": { \"du\": %.3f, \"dv\": %.3f, \"slice_h\": %.4f, "
               "\"n_slices\": %d, \"n_bands\": %d, \"col_min\": %ld, \"col_max\": %ld, "
               "\"atlas_cols\": %zu, \"runs\": %d, \"row_z_deviations\": %zu },\n",
            r->du, r->dv, r->slice_h, r->n_slices, r->n_bands, ctx->L->col_min,
            ctx->L->col_max, ctx->ncols, ctx->n_runs, ctx->L->row_z_deviations);
    fprintf(f, "  \"input_counts\": { \"nv\": %zu, \"nf\": %zu, \"observed_fixed\": %zu, "
               "\"observed_generated\": %zu, \"fills_treated_unobserved\": %zu },\n",
            r->nv_in, r->nf_in, r->observed_fixed, r->observed_generated, r->input_fills);
    fprintf(f, "  \"selection\": { \"observed\": %zu, \"band1_adopted\": %zu, "
               "\"band1_rejected_other_wrap\": %zu, \"band1_no_reference\": %zu, "
               "\"duplicate_cells\": %zu, \"double_cover\": %zu, "
               "\"duplicate_dr_pitch_p10\": %.3f, \"duplicate_dr_pitch_p50\": %.3f, "
               "\"duplicate_dr_pitch_p90\": %.3f, \"adopt_dr_pitch\": %.3f },\n",
            r->observed, r->band1_adopted, r->band1_rejected, r->band1_no_reference,
            r->duplicates, r->double_cover, r->dup_dr_pitch_p10, r->dup_dr_pitch_p50,
            r->dup_dr_pitch_p90, QS_BAND1_ADOPT_DR_PITCH);
    fprintf(f, "  \"fill\": { \"cells_filled\": %zu, \"dropped_unbridged\": %zu, "
               "\"dropped_skip\": %zu, \"dropped_order\": %zu, \"dropped_chord\": %zu, "
               "\"dropped_material\": %zu, "
               "\"dropped_bbox\": %zu, \"dropped_core\": %zu, \"dropped_unreferenced\": %zu, "
               "\"policy\": \"every row of every column between the first and last good column; "
               "gaps refused by wrap-skip, material or width; row/column bridges vetoed by the "
               "chord, material and wrap-skip gates; confetti runs not emitted\", "
               "\"solid_gates\": { \"fold_cells\": %zu, \"fold_columns\": %zu, \"dark_cells\": %zu, "
               "\"dark_snapped\": %zu, \"dark_unrecovered\": %zu, \"hole_cells\": %zu, "
               "\"refused_hole_cells\": %zu, \"gates_failed\": %d, "
               "\"midline_probed\": %zu, \"midline_outliers\": %zu, \"midline_isolated\": %zu, \"midline_dense\": %zu, "
               "\"midline_p50\": %.3f, \"midline_p90\": %.3f },\n"
               "  \"scale\": { \"faces\": %zu, \"over\": %zu, \"over_fill\": %zu, \"limit_log2\": %.3f, "
               "\"p50\": %.4f, \"p90\": %.4f, \"p99\": %.4f, \"observed_p90\": %.4f, \"fill_p90\": %.4f },\n"
               "  \"fill_support\": { \"cells\": %zu, \"far\": %zu, \"limit_cells\": %d, \"p50\": %.1f, "
               "\"p90\": %.1f, \"max\": %zu },\n"
               "  \"ramp\": { \"runs\": %zu, \"runs_over\": %zu, \"cells\": %zu, \"cells_over\": %zu, "
               "\"limit_pitches_per_turn\": %.2f, \"p50\": %.3f, \"p90\": %.3f, \"max\": %.3f },\n"
               "              \"fill_snap\": { \"candidates\": %zu, \"moved\": %zu, \"no_ridge\": %zu, "
               "\"shift_p50\": %.3f, \"shift_p90\": %.3f, \"reach\": %.1f, \"min_u8\": %d },\n"
               "              \"confetti_observed\": %zu, \"run_min_observed\": %d, "
               "\"good_column_frac\": %.3f, \"max_u_gap\": %d, \"skip_dr_pitch\": %.2f, "
               "\"bridge_max_cells\": %d, \"stretch\": %.2f, \"chord_abs\": %.1f, "
               "\"core_radius\": %.1f, \"bbox_margin\": %.1f, "
               "\"runs_built\": %d, \"runs_retried\": %d, \"runs_observed_only\": %d },\n",
            r->filled, r->dropped_unbridged, r->dropped_skip, r->dropped_order, r->dropped_chord,
            r->dropped_material, r->dropped_bbox, r->dropped_core,
            r->dropped_unreferenced,

            r->fold_cells, r->fold_columns, r->dark_cells, r->dark_snapped, r->dark_unrecovered,

            r->hole_cells, r->refused_hole_cells, r->gates_failed,

            r->midline_probed, r->midline_outliers, r->midline_isolated, r->midline_dense, r->midline_p50, r->midline_p90,
            r->scale_faces, r->scale_over, r->scale_over_fill, (double)QS_GATE_SCALE_LOG2,
            r->scale_p50, r->scale_p90, r->scale_p99, r->scale_obs_p90, r->scale_fill_p90,
            r->fill_cells, r->fill_far, (int)QS_FILL_SUPPORT_CELLS, r->fill_dist_p50, r->fill_dist_p90,
            r->fill_dist_max,
            r->ramp_runs, r->ramp_runs_over, r->ramp_cells, r->ramp_cells_over,
            (double)QS_GATE_RAMP_PITCHES, r->ramp_p50, r->ramp_p90, r->ramp_max,
            r->snap_candidates, r->snap_moved, r->snap_no_ridge, r->snap_shift_p50, r->snap_shift_p90,

            QS_FILL_SNAP_REACH, QS_FILL_SNAP_MIN_U8,

            r->dropped_confetti,
            ctx->cfg->run_min_observed > 0 ? ctx->cfg->run_min_observed : QS_RUN_MIN_OBSERVED,
            ctx->cfg->good_column_frac > 0.0 ? ctx->cfg->good_column_frac : QS_GOOD_COLUMN_FRAC,
            ctx->cfg->fill_max_u_gap > 0 ? ctx->cfg->fill_max_u_gap : QS_FILL_MAX_U_GAP,
            ctx->cfg->skip_dr_pitch > 0.0 ? ctx->cfg->skip_dr_pitch : QS_SKIP_DR_PITCH,
            QS_BRIDGE_MAX_CELLS, QS_FILL_STRETCH, QS_FILL_CHORD_ABS, QS_CORE_RADIUS_VOX,
            QS_BBOX_MARGIN_VOX, r->runs_built, r->runs_retried, r->runs_observed_only);
    fprintf(f, "  \"emission\": { \"nv\": %zu, \"nf\": %zu, \"faces_skipped_gate\": %zu, "
               "\"faces_skipped_material\": %zu, \"wrap_gate_vox\": %.1f },\n",
            r->nv_out, r->nf_out, r->faces_skipped_gate, r->faces_skipped_material,
            QS_WRAP_GATE_VOX);
    fprintf(f, "  \"phase\": { \"mismatch_max_rad\": %.4f, \"mismatch_over_limit_cells\": %zu, "
               "\"limit_rad\": %.4f },\n",
            r->phase_mismatch_max, r->phase_mismatch_over_limit, QS_PHASE_MISMATCH_MAX);
    fprintf(f, "  \"provenance_png\": { \"empty\": [0,0,0], \"observed_fixed\": [200,200,200], "
               "\"observed_generated\": [120,120,120], \"band1_adopted\": [80,120,255], "
               "\"fill\": [60,200,60], \"unbridged\": [110,20,20], \"duplicate_cell\": [230,150,40], "
               "\"band1_rejected\": [200,60,200], \"bridge_refused\": [120,40,160], "
               "\"confetti_run\": [70,90,120] },\n");
    fprintf(f, "  \"runs\": [\n");
    for (int i = 0; i < ctx->n_runs; i++) {
        const QsRunStat *s = &ctx->runs[i];
        fprintf(f, "    { \"run\": %d, \"c0\": %ld, \"c1\": %ld, \"width\": %d, \"c0_out\": %ld, "
                   "\"holes\": %zu, \"holes_refused\": %zu, \"holes_filled_by_vote\": %zu, "
                   "\"observed\": %zu, \"adopted\": %zu, \"rejected\": %zu, \"no_ref\": %zu, "
                   "\"duplicates\": %zu, \"filled\": %zu, \"dropped_unbridged\": %zu, "
                   "\"dropped_skip\": %zu, \"dropped_order\": %zu, \"dropped_chord\": %zu, "
                   "\"dropped_material\": %zu, "
                   "\"dropped_bbox\": %zu, \"dropped_core\": %zu, \"dropped_unref\": %zu, "
                   "\"faces\": %zu, \"skipped_gate\": %zu, \"skipped_material\": %zu, "
                   "\"build\": \"%s\", \"t_build_sec\": %.2f }%s\n",
                i, s->c0, s->c1, s->width, s->c0_out, s->holes, s->holes_refused, s->holes_filled,
                s->observed, s->adopted, s->rejected, s->no_ref,
                s->duplicates, s->filled, s->dropped_unbridged, s->dropped_skip, s->dropped_order,
                s->dropped_chord, s->dropped_material, s->dropped_bbox, s->dropped_core,
                s->dropped_unref, s->faces, s->skipped_gate, s->skipped_material,
                s->build < 0 ? "confetti" : (s->build == 0 ? "observed-only" : (s->build == 1 ? "ok" : "retry-ok")),
                s->t_build, i + 1 < ctx->n_runs ? "," : "");
    }
    fprintf(f, "  ],\n");
    fprintf(f, "  \"timing_sec\": { \"load\": %.2f, \"lattice\": %.2f, \"build\": %.2f, "
               "\"emit\": %.2f, \"write\": %.2f, \"total\": %.2f }\n}\n",
            r->t_load, r->t_lattice, r->t_build, r->t_emit, r->t_write, r->t_total);
    return fclose(f) == 0 ? 0 : -1;
}

static const char *const qs_sidecar_suffix[] = {
    "_support.u8", "_provenance.u8", "_phase.f32", "_material_identity.i32",
    "_reconstruction_component.i32", "_lane.i32", "_atlas_column.i32"
};

static int qs_write_outputs(QsCtx *ctx, const char *out_stem, const char *in_vmesh)
{
    char path[QS_PATH], dir[QS_PATH];
    const QsOut *o = &ctx->out;
    size_t nv = o->nv;
    if (nv == 0 || o->nf == 0) {
        fprintf(stderr, "[solidify] nothing to emit (nv=%zu nf=%zu)\n", nv, o->nf);
        return -1;
    }
    /* companions first, VMESH last: a crash cannot pair a fresh mesh with
     * stale provenance */
    snprintf(path, sizeof path, "%s_support.u8", out_stem);
    if (qs_write_raw(path, o->support, nv) != 0) return -1;
    snprintf(path, sizeof path, "%s_provenance.u8", out_stem);
    if (qs_write_raw(path, o->prov, nv) != 0) return -1;
    snprintf(path, sizeof path, "%s_phase.f32", out_stem);
    if (qs_write_raw(path, o->phase, nv * sizeof(float)) != 0) return -1;
    snprintf(path, sizeof path, "%s_material_identity.i32", out_stem);
    if (qs_write_raw(path, o->material, nv * sizeof(int32_t)) != 0) return -1;
    snprintf(path, sizeof path, "%s_reconstruction_component.i32", out_stem);
    if (qs_write_raw(path, o->run, nv * sizeof(int32_t)) != 0) return -1;
    snprintf(path, sizeof path, "%s_lane.i32", out_stem);
    if (qs_write_raw(path, o->lane, nv * sizeof(int32_t)) != 0) return -1;
    snprintf(path, sizeof path, "%s_atlas_column.i32", out_stem);
    if (qs_write_raw(path, o->atlas_col, nv * sizeof(int32_t)) != 0) return -1;
    /* the metric stages discover sidecars by FIXED names beside their input */
    qs_dir_of(out_stem, dir, sizeof dir);
    snprintf(path, sizeof path, "%s/ribbon_phase.f32", dir);
    if (qs_write_raw(path, o->phase, nv * sizeof(float)) != 0) return -1;
    snprintf(path, sizeof path, "%s/ribbon_material_identity.i32", dir);
    if (qs_write_raw(path, o->material, nv * sizeof(int32_t)) != 0) return -1;
    /* minimal stats sidecar so a refit recovers the band structure exactly */
    snprintf(path, sizeof path, "%s_stats.json", out_stem);
    {
        FILE *f = fopen(path, "wb");
        if (f == NULL) return -1;
        fprintf(f, "{ \"schema\": \"vesuvius-quadribbon-solid-stats-v1\",\n"
                   "  \"ribbon\": { \"nk\": %d, \"nu\": %zu, \"du\": %.3f, \"dv\": %.3f },\n"
                   "  \"slicing\": { \"n_slices\": %d } }\n",
                ctx->L->n_slices, ctx->ncols, ctx->L->du, ctx->L->dv, ctx->L->n_slices);
        if (fclose(f) != 0) return -1;
    }
    snprintf(path, sizeof path, "%s_provenance.png", out_stem);
    if (qs_write_provenance_png(ctx, path) != 0)
        fprintf(stderr, "[solidify] WARNING: cannot write %s\n", path);
    snprintf(path, sizeof path, "%s_report.json", out_stem);
    if (qs_write_report(ctx, path, in_vmesh, out_stem) != 0) return -1;
    snprintf(path, sizeof path, "%s.vmesh", out_stem);
    if (MeshBin_write(path, o->verts, nv, o->faces, o->nf, o->uv) != 0) {
        fprintf(stderr, "[solidify] cannot write %s\n", path);
        return -1;
    }
    return 0;
}

/* ---- projective local authority path ------------------------------------
 *
 * The legacy solidifier above deliberately reasons about whole column runs.
 * That is useful for exploratory repair, but it is not a projective map: a
 * new row can merge runs, change a median/trend, alter a harmonic domain, and
 * move old cells.  The authority path below has a much smaller contract.
 *
 *   - a lattice address is absolute (col,row), from configured du/dv;
 *   - every non-fill input sample is copied bit-for-bit, while preserving
 *     whether Stage 3 obtained it from direct support or bounded continuation;
 *   - a missing address is reconstructed only from the two nearest
 *     directly-supported observations which bracket it on one lattice axis,
 *     within a fixed radius; upstream-generated and reconstructed cells are
 *     never evidence for another cell;
 *   - every upstream sample retains its occupied address; at a genuinely
 *     missing address horizontal brackets win over vertical deterministically;
 *   - faces depend only on the four addresses of one unit lattice square.
 *
 * Thus the dependency certificate is exact, not statistical: a direct cell
 * names itself, an upstream-generated cell preserves that classification and
 * names itself, and a locally reconstructed cell names its two direct endpoint
 * addresses.  Extending an input can affect a cell only through new evidence
 * in its declared Stage-3 or Stage-4 stencil.
 */

enum {
    QS_DEP_OBSERVED = 0,
    QS_DEP_HORIZONTAL = 1,
    QS_DEP_VERTICAL = 2,
    QS_DEP_UPSTREAM_GENERATED = 3
};

/* Stage 3 bridges at most eight missing U columns and then eight missing V
 * rows.  Including the endpoint gives a Cartesian 9x9 evidence radius for a
 * copied Stage-3 generated sample. */
#define QS_UPSTREAM_EVIDENCE_RADIUS_CELLS (QS_BRIDGE_MAX_CELLS + 1)

typedef struct QsProjectiveEntry {
    int32_t col, row;
    int32_t source_a, source_b; /* input vertex indices */
    uint8_t kind;
} QsProjectiveEntry;

typedef struct QsProjectiveOut {
    size_t nv, nf;
    float *verts, *uv, *phase;
    int32_t *faces, *material, *claimant_chart, *run, *lane, *atlas_col;
    int32_t *address;             /* [nv][2] col,row */
    int32_t *dependency;          /* [nv][4] a_col,a_row,b_col,b_row */
    uint8_t *support, *prov, *dependency_kind;
} QsProjectiveOut;

static int qs_projective_kind_rank(uint8_t kind)
{
    /* Preserve the upstream artifact before constructing any new geometry.
     * A local fill is considered only at an address Stage 3 left absent. */
    switch (kind) {
        case QS_DEP_OBSERVED: return 0;
        case QS_DEP_UPSTREAM_GENERATED: return 1;
        case QS_DEP_HORIZONTAL: return 2;
        case QS_DEP_VERTICAL: return 3;
        default: return 4;
    }
}

static int qs_projective_is_copied_input(uint8_t kind)
{
    return kind == QS_DEP_OBSERVED ||
           kind == QS_DEP_UPSTREAM_GENERATED;
}

static int qs_projective_entry_row_compare(const void *pa, const void *pb)
{
    const QsProjectiveEntry *a = (const QsProjectiveEntry *)pa;
    const QsProjectiveEntry *b = (const QsProjectiveEntry *)pb;
    if (a->row != b->row) return a->row < b->row ? -1 : 1;
    if (a->col != b->col) return a->col < b->col ? -1 : 1;
    if (a->kind != b->kind) {
        int ar = qs_projective_kind_rank(a->kind);
        int br = qs_projective_kind_rank(b->kind);
        return ar < br ? -1 : 1;
    }
    if (a->source_a != b->source_a)
        return a->source_a < b->source_a ? -1 : 1;
    if (a->source_b != b->source_b)
        return a->source_b < b->source_b ? -1 : 1;
    return 0;
}

static int qs_projective_entry_col_compare(const void *pa, const void *pb)
{
    const QsProjectiveEntry *a = (const QsProjectiveEntry *)pa;
    const QsProjectiveEntry *b = (const QsProjectiveEntry *)pb;
    if (a->col != b->col) return a->col < b->col ? -1 : 1;
    if (a->row != b->row) return a->row < b->row ? -1 : 1;
    return 0;
}

static int qs_projective_push(QsProjectiveEntry **entry, size_t *n,
                              size_t *cap, QsProjectiveEntry value)
{
    if (*n == *cap) {
        size_t next = *cap == 0 ? 4096 : *cap * 2;
        QsProjectiveEntry *p = NULL;
        if (next < *cap || next > SIZE_MAX / sizeof *p) return -1;
        p = (QsProjectiveEntry *)realloc(*entry, next * sizeof *p);
        if (p == NULL) return -1;
        *entry = p; *cap = next;
    }
    (*entry)[(*n)++] = value;
    return 0;
}

/* Deterministic claimant for the defensive duplicate-address case.  The
 * source ordinal is consulted only after every published field is identical,
 * so input reordering cannot change the output bits. */
static int qs_projective_prefer_source(const QsInput *in, int32_t ia,
                                       int32_t ib)
{
    size_t a = (size_t)ia, b = (size_t)ib;
    int c = 0;
    if (in->support[a] != in->support[b])
        return in->support[a] > in->support[b] ? ia : ib;
    if (in->prov[a] != in->prov[b])
        return in->prov[a] < in->prov[b] ? ia : ib;
    if (in->material != NULL && in->material[a] != in->material[b])
        return in->material[a] < in->material[b] ? ia : ib;
    if (in->claimant_chart != NULL &&
        in->claimant_chart[a] != in->claimant_chart[b])
        return in->claimant_chart[a] < in->claimant_chart[b] ? ia : ib;
    c = memcmp(in->mesh.verts + a * 3, in->mesh.verts + b * 3,
               3 * sizeof(float));
    if (c != 0) return c < 0 ? ia : ib;
    c = memcmp(in->mesh.uv + a * 2, in->mesh.uv + b * 2,
               2 * sizeof(float));
    if (c != 0) return c < 0 ? ia : ib;
    c = memcmp(in->phase + a, in->phase + b, sizeof(float));
    if (c != 0) return c < 0 ? ia : ib;
    return ia < ib ? ia : ib;
}

static int qs_projective_bridge_ok(const QsConfig *cfg, const QsInput *in,
                                   int32_t ia, int32_t ib, int steps,
                                   int horizontal)
{
    const float *a = in->mesh.verts + (size_t)ia * 3;
    const float *b = in->mesh.verts + (size_t)ib * 3;
    double dz = (double)b[0] - a[0], dy = (double)b[1] - a[1];
    double dx = (double)b[2] - a[2];
    double chord = sqrt(dz * dz + dy * dy + dx * dx);
    double pa = (double)in->phase[ia], pb = (double)in->phase[ib];
    double dphase = pb - pa;
    double ra = hypot((double)a[1] - cfg->axis_y,
                      (double)a[2] - cfg->axis_x);
    double rb = hypot((double)b[1] - cfg->axis_y,
                      (double)b[2] - cfg->axis_x);
    if (steps <= 1 || !isfinite(chord) || !isfinite(pa) || !isfinite(pb) ||
        !(ra > 1.0) || !(rb > 1.0))
        return 0;
    if (in->material == NULL || in->material[ia] < 0 ||
        in->material[ia] != in->material[ib])
        return 0;
    if (in->claimant_chart != NULL &&
        in->claimant_chart[ia] != in->claimant_chart[ib])
        return 0;
    if (chord > QS_FILL_CHORD_ABS ||
        chord > (double)steps * QS_WRAP_GATE_VOX)
        return 0;
    if (horizontal) {
        double want = (double)steps * cfg->lattice_du / (0.5 * (ra + rb));
        /* A different turn differs by about 2*pi.  This tolerance admits
         * local geometric noise while retaining a wide deterministic margin
         * to that wrong-wrap case. */
        if (fabs(fabs(dphase) - want) > 0.50 * want + 0.20) return 0;
    } else {
        /* One canonical column follows one angle down z; a turn jump is not a
         * vertical hole. */
        if (fabs(dphase) > 0.75) return 0;
    }
    return 1;
}

static long qs_projective_find(const QsProjectiveEntry *entry, size_t n,
                               int32_t col, int32_t row)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (entry[mid].row < row ||
            (entry[mid].row == row && entry[mid].col < col))
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo < n && entry[lo].row == row && entry[lo].col == col
         ? (long)lo : -1;
}

static void qs_projective_out_dispose(QsProjectiveOut *o)
{
    if (o == NULL) return;
    free(o->verts); free(o->uv); free(o->phase); free(o->faces);
    free(o->material); free(o->claimant_chart); free(o->run); free(o->lane);
    free(o->atlas_col);
    free(o->address); free(o->dependency); free(o->support); free(o->prov);
    free(o->dependency_kind);
    memset(o, 0, sizeof *o);
}

static int qs_write_projective_outputs(const QsProjectiveOut *o,
                                       const QsConfig *cfg,
                                       const char *in_vmesh,
                                       const char *out_stem,
                                       const QsReport *rep,
                                       int fill_radius, size_t nrows)
{
    char path[QS_PATH], dir[QS_PATH], in_esc[1024], out_esc[1024];
    FILE *f = NULL;
    size_t nv = o->nv;
    const char *schema = cfg->projective_copy_only
                       ? "vesuvius-quadribbon-solid-projective-v5"
                       : "vesuvius-quadribbon-solid-projective-v4";
    const char *contract = cfg->projective_copy_only
        ? "lossless semantic atlas-member transport; exact Stage-3 vertices, "
          "faces, claimant identities, coordinates and sample fields; no "
          "storage-color-local reconstruction"
        : "absolute lattice; direct observations immutable; upstream "
          "continuation explicitly classified; direct-only fixed-radius "
          "brackets; non-recursive fill; copied input faces; local "
          "reconstructed faces; stable material-key labels";
    snprintf(path, sizeof path, "%s_support.u8", out_stem);
    if (qs_write_raw(path, o->support, nv) != 0) return -1;
    snprintf(path, sizeof path, "%s_provenance.u8", out_stem);
    if (qs_write_raw(path, o->prov, nv) != 0) return -1;
    snprintf(path, sizeof path, "%s_phase.f32", out_stem);
    if (qs_write_raw(path, o->phase, nv * sizeof(float)) != 0) return -1;
    snprintf(path, sizeof path, "%s_material_identity.i32", out_stem);
    if (qs_write_raw(path, o->material, nv * sizeof(int32_t)) != 0) return -1;
    snprintf(path, sizeof path, "%s_claimant_chart.i32", out_stem);
    if (qs_write_raw(path, o->claimant_chart, nv * sizeof(int32_t)) != 0)
        return -1;
    snprintf(path, sizeof path, "%s_reconstruction_component.i32", out_stem);
    if (qs_write_raw(path, o->run, nv * sizeof(int32_t)) != 0) return -1;
    snprintf(path, sizeof path, "%s_lane.i32", out_stem);
    if (qs_write_raw(path, o->lane, nv * sizeof(int32_t)) != 0) return -1;
    snprintf(path, sizeof path, "%s_atlas_column.i32", out_stem);
    if (qs_write_raw(path, o->atlas_col, nv * sizeof(int32_t)) != 0) return -1;
    snprintf(path, sizeof path, "%s_lattice_address.i32", out_stem);
    if (qs_write_raw(path, o->address, nv * 2 * sizeof(int32_t)) != 0)
        return -1;
    snprintf(path, sizeof path, "%s_dependency_kind.u8", out_stem);
    if (qs_write_raw(path, o->dependency_kind, nv) != 0) return -1;
    snprintf(path, sizeof path, "%s_dependency_endpoints.i32", out_stem);
    if (qs_write_raw(path, o->dependency, nv * 4 * sizeof(int32_t)) != 0)
        return -1;

    qs_dir_of(out_stem, dir, sizeof dir);
    snprintf(path, sizeof path, "%s/ribbon_phase.f32", dir);
    if (qs_write_raw(path, o->phase, nv * sizeof(float)) != 0) return -1;
    snprintf(path, sizeof path, "%s/ribbon_material_identity.i32", dir);
    if (qs_write_raw(path, o->material, nv * sizeof(int32_t)) != 0) return -1;

    snprintf(path, sizeof path, "%s_stats.json", out_stem);
    f = fopen(path, "wb");
    if (f == NULL) return -1;
    fprintf(f, "{ \"schema\": \"vesuvius-quadribbon-projective-stats-v1\",\n"
               "  \"ribbon\": { \"nk\": %zu, \"nu\": %zu, \"du\": %.9g, "
               "\"dv\": %.9g },\n  \"slicing\": { \"n_slices\": %zu } }\n",
            nrows, nv, cfg->lattice_du, cfg->lattice_dv, nrows);
    if (fclose(f) != 0) return -1;

    qs_json_escape(in_vmesh, in_esc, sizeof in_esc);
    qs_json_escape(out_stem, out_esc, sizeof out_esc);
    snprintf(path, sizeof path, "%s_report.json", out_stem);
    f = fopen(path, "wb");
    if (f == NULL) return -1;
    fprintf(f,
            "{\n  \"schema\": \"%s\",\n"
            "  \"input\": \"%s\",\n  \"output\": \"%s.vmesh\",\n"
            "  \"contract\": \"%s\",\n"
            "  \"lattice\": { \"du\": %.9g, \"dv\": %.9g, "
            "\"rows_present\": %zu },\n"
            "  \"dependency\": { \"copy_only\": %s, \"fill_radius\": %d, "
            "\"upstream_generated_radius\": %d, "
            "\"recursive\": false, \"horizontal_precedes_vertical\": true, "
            "\"upstream_generated_precedes_local_fill\": true, "
            "\"upstream_generated_is_evidence\": false, "
            "\"address_sidecar\": \"%s_lattice_address.i32\", "
            "\"endpoint_sidecar\": \"%s_dependency_endpoints.i32\", "
            "\"claimant_sidecar\": \"%s_claimant_chart.i32\" },\n"
            "  \"input_counts\": { \"vertices\": %zu, \"faces\": %zu, "
            "\"prior_fills_excluded\": %zu },\n"
            "  \"output_counts\": { \"vertices\": %zu, \"faces\": %zu, "
            "\"observed\": %zu, \"filled\": %zu, \"duplicate_addresses\": %zu, "
            "\"faces_skipped_gate\": %zu, \"faces_skipped_material\": %zu }\n}\n",
            schema, in_esc, out_esc, contract,
            cfg->lattice_du, cfg->lattice_dv, nrows,
            cfg->projective_copy_only ? "true" : "false", fill_radius,
            QS_UPSTREAM_EVIDENCE_RADIUS_CELLS,
            out_esc, out_esc, out_esc, rep->nv_in, rep->nf_in,
            rep->input_fills, rep->nv_out, rep->nf_out, rep->observed,
            rep->filled, rep->duplicates, rep->faces_skipped_gate,
            rep->faces_skipped_material);
    if (fclose(f) != 0) return -1;

    snprintf(path, sizeof path, "%s.vmesh", out_stem);
    if (MeshBin_write(path, o->verts, o->nv, o->faces, o->nf, o->uv) != 0) {
        fprintf(stderr, "[solidify-projective] cannot write %s\n", path);
        return -1;
    }
    return 0;
}

/* The authoritative projective Stage-4 operation is deliberately an identity
 * on the Stage-3 semantic atlas.  Collision colors are only storage: running a
 * local hole finder independently inside each color would make generated
 * geometry depend on how a larger domain happened to recolor the same
 * claimant.  This path therefore copies one atlas member exactly, adds only
 * absolute-address/dependency certificates, and leaves optional reconstruction
 * to a non-authoritative presentation layer. */
static int qs_projective_copy_run(const QsConfig *cfg, const QsInput *in,
                                  const char *in_vmesh,
                                  const char *out_stem, QsReport *rep,
                                  double t0)
{
    QsProjectiveOut out;
    size_t nv = in->mesh.nv, nf = in->mesh.nf, nrows = 0;
    int32_t row_min = INT32_MAX, row_max = INT32_MIN;
    int rc = -1;
    memset(&out, 0, sizeof out);
    if (!(cfg->lattice_du > 0.0) || !(cfg->lattice_dv > 0.0) ||
        !isfinite(cfg->lattice_du) || !isfinite(cfg->lattice_dv) ||
        in->material == NULL || in->claimant_chart == NULL ||
        nv == 0 || nf == 0 ||
        nv > (size_t)INT32_MAX ||
        nv > SIZE_MAX / (4 * sizeof(int32_t)) ||
        nf > SIZE_MAX / (3 * sizeof(int32_t))) {
        fprintf(stderr,
                "[solidify-projective] copy-only atlas member requires a "
                "finite canonical lattice and exact material, claimant-chart "
                "sidecars (provenance is derived canonically from support "
                "when absent)\n");
        return -1;
    }

    out.nv = nv; out.nf = nf;
    /* Core semantic arrays are immutable arena-backed input data.  Borrow
     * them while writing instead of doubling the peak memory of a 10x10x10
     * chart; only the new certificate arrays need allocation. */
    out.verts = (float *)in->mesh.verts;
    out.uv = (float *)in->mesh.uv;
    out.phase = (float *)in->phase;
    out.faces = (int32_t *)in->mesh.faces;
    out.material = (int32_t *)in->material;
    out.claimant_chart = (int32_t *)in->claimant_chart;
    out.run = (int32_t *)malloc(nv * sizeof(int32_t));
    out.lane = (int32_t *)malloc(nv * sizeof(int32_t));
    out.atlas_col = (int32_t *)malloc(nv * sizeof(int32_t));
    out.address = (int32_t *)malloc(nv * 2 * sizeof(int32_t));
    out.dependency = (int32_t *)malloc(nv * 4 * sizeof(int32_t));
    out.support = (uint8_t *)in->support;
    out.prov = (uint8_t *)in->prov;
    out.dependency_kind = (uint8_t *)malloc(nv);
    if (out.verts == NULL || out.uv == NULL || out.phase == NULL ||
        out.faces == NULL || out.material == NULL ||
        out.claimant_chart == NULL || out.run == NULL || out.lane == NULL ||
        out.atlas_col == NULL || out.address == NULL ||
        out.dependency == NULL || out.support == NULL || out.prov == NULL ||
        out.dependency_kind == NULL)
        goto done;

    for (size_t i = 0; i < nv; i++) {
        double cq = (double)in->mesh.uv[i * 2] / cfg->lattice_du;
        double rq = (double)in->mesh.uv[i * 2 + 1] / cfg->lattice_dv;
        long long cc = llround(cq), rr = llround(rq);
        int32_t semantic_label = in->material[i];
        if (!isfinite(cq) || !isfinite(rq) ||
            !isfinite((double)in->phase[i]) ||
            fabs(cq - (double)cc) > QS_LATTICE_EPS ||
            fabs(rq - (double)rr) > QS_LATTICE_EPS ||
            cc < INT32_MIN || cc > INT32_MAX ||
            rr < INT32_MIN || rr > INT32_MAX) {
            fprintf(stderr, "[solidify-projective] input vertex %zu is not "
                            "finite canonical lattice evidence\n", i);
            goto done;
        }
        out.address[i * 2] = (int32_t)cc;
        out.address[i * 2 + 1] = (int32_t)rr;
        out.atlas_col[i] = (int32_t)cc;
        out.run[i] = in->reconstruction != NULL
                   ? in->reconstruction[i] : semantic_label;
        out.lane[i] = in->lane != NULL ? in->lane[i] : out.run[i];
        out.dependency_kind[i] = in->support[i] != 0
                               ? QS_DEP_OBSERVED
                               : QS_DEP_UPSTREAM_GENERATED;
        out.dependency[i * 4] = (int32_t)cc;
        out.dependency[i * 4 + 1] = (int32_t)rr;
        out.dependency[i * 4 + 2] = (int32_t)cc;
        out.dependency[i * 4 + 3] = (int32_t)rr;
        if ((int32_t)rr < row_min) row_min = (int32_t)rr;
        if ((int32_t)rr > row_max) row_max = (int32_t)rr;
    }
    if (in->stats_n_slices > 0) {
        nrows = (size_t)in->stats_n_slices;
    } else if (row_max >= row_min) {
        uint64_t span = (uint64_t)((int64_t)row_max -
                                   (int64_t)row_min) + 1u;
        nrows = span > (uint64_t)SIZE_MAX ? SIZE_MAX : (size_t)span;
    }
    rep->observed = nv;
    rep->filled = 0;
    rep->nv_out = nv; rep->nf_out = nf;
    rep->du = cfg->lattice_du; rep->dv = cfg->lattice_dv;
    rep->n_slices = (int)(nrows > INT_MAX ? INT_MAX : nrows);
    rep->cols_out = 0;
    rep->t_total = ves_clock_sec() - t0;
    rc = qs_write_projective_outputs(&out, cfg, in_vmesh, out_stem, rep,
                                      0, nrows);
    fprintf(stderr,
            "[solidify-projective] %s copy-only atlas member: %zu samples, "
            "%zu faces, zero storage-color-local fills\n",
            rc == 0 ? "wrote" : "FAILED", nv, nf);
done:
    /* Do not free the borrowed Arena-owned semantic arrays. */
    out.verts = NULL; out.uv = NULL; out.phase = NULL; out.faces = NULL;
    out.material = NULL; out.claimant_chart = NULL;
    out.support = NULL; out.prov = NULL;
    qs_projective_out_dispose(&out);
    return rc;
}

static int qs_projective_run(const QsConfig *cfg, const QsInput *in,
                             const char *in_vmesh, const char *out_stem,
                             QsReport *rep, double t0)
{
    QsProjectiveEntry *entry = NULL, *direct = NULL, *by_col = NULL;
    int32_t *source_out = NULL;
    QsProjectiveOut out;
    size_t n = 0, cap = 0, input_n = 0, direct_n = 0;
    size_t i = 0, w = 0, nrows = 0;
    int fill_radius = cfg->local_fill_radius > 0
                    ? cfg->local_fill_radius : QS_BRIDGE_MAX_CELLS;
    int rc = -1;
    memset(&out, 0, sizeof out);
    if (cfg->projective_copy_only)
        return qs_projective_copy_run(cfg, in, in_vmesh, out_stem, rep, t0);
    if (!(cfg->lattice_du > 0.0) || !(cfg->lattice_dv > 0.0) ||
        !isfinite(cfg->lattice_du) || !isfinite(cfg->lattice_dv) ||
        fill_radius < 1 || in->material == NULL ||
        in->mesh.nv > (size_t)INT32_MAX) {
        fprintf(stderr, "[solidify-projective] invalid canonical lattice, "
                "missing material identity, or too many input vertices\n");
        return -1;
    }

    for (i = 0; i < in->mesh.nv; i++) {
        double cq = (double)in->mesh.uv[i * 2] / cfg->lattice_du;
        double rq = (double)in->mesh.uv[i * 2 + 1] / cfg->lattice_dv;
        long long cc = llround(cq), rr = llround(rq);
        QsProjectiveEntry e;
        if (in->prov[i] == QS_PROV_FILL) continue;
        if (!isfinite(cq) || !isfinite(rq) ||
            fabs(cq - (double)cc) > QS_LATTICE_EPS ||
            fabs(rq - (double)rr) > QS_LATTICE_EPS ||
            cc < INT32_MIN || cc > INT32_MAX ||
            rr < INT32_MIN || rr > INT32_MAX ||
            !isfinite((double)in->phase[i])) {
            fprintf(stderr, "[solidify-projective] input vertex %zu is not "
                    "finite canonical lattice evidence\n", i);
            goto done;
        }
        e.col = (int32_t)cc; e.row = (int32_t)rr;
        e.source_a = e.source_b = (int32_t)i;
        e.kind = in->support[i] != 0
               ? QS_DEP_OBSERVED : QS_DEP_UPSTREAM_GENERATED;
        if (qs_projective_push(&entry, &n, &cap, e) != 0) goto done;
    }
    if (n == 0) goto done;
    qsort(entry, n, sizeof *entry, qs_projective_entry_row_compare);
    /* Exact-address claimant selection, independent of source ordering. */
    for (i = 0; i < n; ) {
        size_t j = i + 1;
        QsProjectiveEntry best = entry[i];
        while (j < n && entry[j].row == entry[i].row &&
               entry[j].col == entry[i].col) {
            if (entry[j].kind == best.kind)
                best.source_a = best.source_b =
                    qs_projective_prefer_source(in, best.source_a,
                                                 entry[j].source_a);
            rep->duplicates++;
            j++;
        }
        entry[w++] = best;
        i = j;
    }
    n = input_n = w;

    /* Only genuinely direct Stage-3 evidence may create Stage-4 fills.  A
     * copied continuation remains output, but cannot recursively widen the
     * evidence closure or displace a direct bracket endpoint. */
    direct = (QsProjectiveEntry *)malloc(input_n * sizeof *direct);
    if (direct == NULL) goto done;
    for (i = 0; i < input_n; i++)
        if (entry[i].kind == QS_DEP_OBSERVED)
            direct[direct_n++] = entry[i];

    /* Horizontal brackets in row-major order. */
    for (i = 1; i < direct_n; i++) {
        const QsProjectiveEntry *a = &direct[i - 1], *b = &direct[i];
        int64_t steps = (int64_t)b->col - (int64_t)a->col;
        if (a->row != b->row || steps <= 1 ||
            steps - 1 > fill_radius || steps > INT_MAX ||
            !qs_projective_bridge_ok(cfg, in, a->source_a, b->source_a,
                                     (int)steps, 1))
            continue;
        for (int64_t d = 1; d < steps; d++) {
            QsProjectiveEntry e;
            e.col = (int32_t)((int64_t)a->col + d); e.row = a->row;
            e.source_a = a->source_a; e.source_b = b->source_a;
            e.kind = QS_DEP_HORIZONTAL;
            if (qs_projective_push(&entry, &n, &cap, e) != 0) goto done;
        }
    }

    /* Vertical brackets use a separate sort of direct observations only. */
    by_col = (QsProjectiveEntry *)malloc(
        (direct_n > 0 ? direct_n : 1) * sizeof *by_col);
    if (by_col == NULL) goto done;
    if (direct_n > 0)
        memcpy(by_col, direct, direct_n * sizeof *by_col);
    qsort(by_col, direct_n, sizeof *by_col, qs_projective_entry_col_compare);
    for (i = 1; i < direct_n; i++) {
        const QsProjectiveEntry *a = &by_col[i - 1], *b = &by_col[i];
        int64_t steps = (int64_t)b->row - (int64_t)a->row;
        if (a->col != b->col || steps <= 1 ||
            steps - 1 > fill_radius || steps > INT_MAX ||
            !qs_projective_bridge_ok(cfg, in, a->source_a, b->source_a,
                                     (int)steps, 0))
            continue;
        for (int64_t d = 1; d < steps; d++) {
            QsProjectiveEntry e;
            e.col = a->col; e.row = (int32_t)((int64_t)a->row + d);
            e.source_a = a->source_a; e.source_b = b->source_a;
            e.kind = QS_DEP_VERTICAL;
            if (qs_projective_push(&entry, &n, &cap, e) != 0) goto done;
        }
    }
    free(by_col); by_col = NULL;
    free(direct); direct = NULL;

    qsort(entry, n, sizeof *entry, qs_projective_entry_row_compare);
    /* Direct > copied upstream continuation > horizontal > vertical. */
    w = 0;
    for (i = 0; i < n; ) {
        size_t j = i + 1;
        entry[w++] = entry[i];
        while (j < n && entry[j].row == entry[i].row &&
               entry[j].col == entry[i].col)
            j++;
        i = j;
    }
    n = w;
    if (n > (size_t)INT32_MAX || n > SIZE_MAX / 6 ||
        in->mesh.nf > SIZE_MAX - 2 * n ||
        in->mesh.nf + 2 * n > SIZE_MAX / (3 * sizeof(int32_t)))
        goto done;
    out.nv = n;
    out.verts = (float *)malloc(n * 3 * sizeof(float));
    out.uv = (float *)malloc(n * 2 * sizeof(float));
    out.phase = (float *)malloc(n * sizeof(float));
    out.material = (int32_t *)malloc(n * sizeof(int32_t));
    out.claimant_chart = (int32_t *)malloc(n * sizeof(int32_t));
    out.run = (int32_t *)malloc(n * sizeof(int32_t));
    out.lane = (int32_t *)malloc(n * sizeof(int32_t));
    out.atlas_col = (int32_t *)malloc(n * sizeof(int32_t));
    out.address = (int32_t *)malloc(n * 2 * sizeof(int32_t));
    out.dependency = (int32_t *)malloc(n * 4 * sizeof(int32_t));
    out.support = (uint8_t *)malloc(n);
    out.prov = (uint8_t *)malloc(n);
    out.dependency_kind = (uint8_t *)malloc(n);
    out.faces = (int32_t *)malloc(
        (in->mesh.nf + 2 * n + 1) * 3 * sizeof(int32_t));
    source_out = (int32_t *)malloc(
        (in->mesh.nv > 0 ? in->mesh.nv : 1) * sizeof(int32_t));
    if (out.verts == NULL || out.uv == NULL || out.phase == NULL ||
        out.material == NULL || out.claimant_chart == NULL ||
        out.run == NULL || out.lane == NULL ||
        out.atlas_col == NULL || out.address == NULL ||
        out.dependency == NULL || out.support == NULL || out.prov == NULL ||
        out.dependency_kind == NULL || out.faces == NULL || source_out == NULL)
        goto done;
    for (i = 0; i < in->mesh.nv; i++) source_out[i] = -1;

    for (i = 0; i < n; i++) {
        const QsProjectiveEntry *e = &entry[i];
        size_t a = (size_t)e->source_a, b = (size_t)e->source_b;
        out.address[i * 2] = e->col; out.address[i * 2 + 1] = e->row;
        out.atlas_col[i] = e->col;
        out.dependency_kind[i] = e->kind;
        if (qs_projective_is_copied_input(e->kind)) {
            memcpy(out.verts + i * 3, in->mesh.verts + a * 3,
                   3 * sizeof(float));
            memcpy(out.uv + i * 2, in->mesh.uv + a * 2,
                   2 * sizeof(float));
            out.phase[i] = in->phase[a];
            out.material[i] = in->material[a];
            out.claimant_chart[i] = in->claimant_chart != NULL
                                  ? in->claimant_chart[a]
                                  : in->material[a];
            /* The input reconstruction/lane label is assigned by a
             * crop-wide Stage-3 traversal and is therefore presentation
             * metadata, not canonical evidence.  Publishing it here made an
             * otherwise identical A subset B differ in a downstream field.
             * Material identity is the stable block-addressed semantic key. */
            out.lane[i] = in->material[a];
            out.support[i] = in->support[a]; out.prov[i] = in->prov[a];
            out.dependency[i * 4] = e->col;
            out.dependency[i * 4 + 1] = e->row;
            out.dependency[i * 4 + 2] = e->col;
            out.dependency[i * 4 + 3] = e->row;
            source_out[a] = (int32_t)i;
            rep->observed++;
        } else {
            int32_t ac = (int32_t)llround(
                (double)in->mesh.uv[a * 2] / cfg->lattice_du);
            int32_t ar = (int32_t)llround(
                (double)in->mesh.uv[a * 2 + 1] / cfg->lattice_dv);
            int32_t bc = (int32_t)llround(
                (double)in->mesh.uv[b * 2] / cfg->lattice_du);
            int32_t br = (int32_t)llround(
                (double)in->mesh.uv[b * 2 + 1] / cfg->lattice_dv);
            double den = e->kind == QS_DEP_HORIZONTAL
                       ? (double)((int64_t)bc - ac)
                       : (double)((int64_t)br - ar);
            double num = e->kind == QS_DEP_HORIZONTAL
                       ? (double)((int64_t)e->col - ac)
                       : (double)((int64_t)e->row - ar);
            double t = num / den;
            double ra = hypot((double)in->mesh.verts[a * 3 + 1] - cfg->axis_y,
                              (double)in->mesh.verts[a * 3 + 2] - cfg->axis_x);
            double rb = hypot((double)in->mesh.verts[b * 3 + 1] - cfg->axis_y,
                              (double)in->mesh.verts[b * 3 + 2] - cfg->axis_x);
            double radius = ra + t * (rb - ra);
            double phase = (double)in->phase[a] +
                           t * ((double)in->phase[b] - in->phase[a]);
            out.verts[i * 3] = (float)((double)in->mesh.verts[a * 3] +
                                t * ((double)in->mesh.verts[b * 3] -
                                     in->mesh.verts[a * 3]));
            out.verts[i * 3 + 1] =
                (float)(cfg->axis_y + radius * cos(phase));
            out.verts[i * 3 + 2] =
                (float)(cfg->axis_x + radius * sin(phase));
            out.uv[i * 2] = (float)((double)e->col * cfg->lattice_du);
            out.uv[i * 2 + 1] = (float)((double)e->row * cfg->lattice_dv);
            out.phase[i] = (float)phase;
            out.material[i] = in->material[a];
            out.claimant_chart[i] = in->claimant_chart != NULL
                                  ? in->claimant_chart[a]
                                  : in->material[a];
            out.lane[i] = in->material[a];
            out.support[i] = 0; out.prov[i] = QS_PROV_FILL;
            out.dependency[i * 4] = ac; out.dependency[i * 4 + 1] = ar;
            out.dependency[i * 4 + 2] = bc; out.dependency[i * 4 + 3] = br;
            rep->filled++;
        }
        /* Stable semantic identity replaces crop-order run numbering. */
        out.run[i] = out.material[i];
        if (i == 0 || entry[i].row != entry[i - 1].row) nrows++;
    }

    /* Every input face is upstream output in its own right.  Preserve it
     * exactly (only vertex ordinals change); its vertices retain the direct or
     * upstream-generated dependency class needed by the nested-domain gate. */
    for (i = 0; i < in->mesh.nf; i++) {
        int32_t sa = in->mesh.faces[i * 3];
        int32_t sb = in->mesh.faces[i * 3 + 1];
        int32_t sc = in->mesh.faces[i * 3 + 2];
        int32_t a = -1, b = -1, c = -1;
        if (sa < 0 || sb < 0 || sc < 0 ||
            (size_t)sa >= in->mesh.nv || (size_t)sb >= in->mesh.nv ||
            (size_t)sc >= in->mesh.nv ||
            in->prov[sa] == QS_PROV_FILL ||
            in->prov[sb] == QS_PROV_FILL ||
            in->prov[sc] == QS_PROV_FILL)
            continue;
        a = source_out[sa]; b = source_out[sb]; c = source_out[sc];
        if (a < 0 || b < 0 || c < 0 || a == b || b == c || a == c)
            continue;
        out.faces[out.nf * 3] = a;
        out.faces[out.nf * 3 + 1] = b;
        out.faces[out.nf * 3 + 2] = c;
        out.nf++;
    }

    /* A reconstructed triangle depends only on the three endpoint-certified
     * vertices of one unit square.  Fully copied-input triangles were copied
     * above and are deliberately not re-gated. */
    for (i = 0; i < n; i++) {
        long ib = -1, ic = -1, id = -1;
        int legal = 0;
        if (entry[i].col == INT32_MAX || entry[i].row == INT32_MAX) continue;
        ib = qs_projective_find(entry, n, entry[i].col + 1, entry[i].row);
        ic = qs_projective_find(entry, n, entry[i].col, entry[i].row + 1);
        id = qs_projective_find(entry, n, entry[i].col + 1, entry[i].row + 1);
        if (ib < 0 || ic < 0 || id < 0) continue;
        if (!qs_projective_is_copied_input(out.dependency_kind[i]) ||
            !qs_projective_is_copied_input(out.dependency_kind[ib]) ||
            !qs_projective_is_copied_input(out.dependency_kind[ic])) {
            legal = qs_tri_legal(out.verts, out.material, i, (size_t)ib,
                                 (size_t)ic);
            if (legal > 0) {
                out.faces[out.nf * 3] = (int32_t)i;
                out.faces[out.nf * 3 + 1] = (int32_t)ib;
                out.faces[out.nf * 3 + 2] = (int32_t)ic;
                out.nf++;
            } else if (legal == 0) rep->faces_skipped_gate++;
            else rep->faces_skipped_material++;
        }
        if (!qs_projective_is_copied_input(out.dependency_kind[ib]) ||
            !qs_projective_is_copied_input(out.dependency_kind[id]) ||
            !qs_projective_is_copied_input(out.dependency_kind[ic])) {
            legal = qs_tri_legal(out.verts, out.material, (size_t)ib,
                                 (size_t)id, (size_t)ic);
            if (legal > 0) {
                out.faces[out.nf * 3] = (int32_t)ib;
                out.faces[out.nf * 3 + 1] = (int32_t)id;
                out.faces[out.nf * 3 + 2] = (int32_t)ic;
                out.nf++;
            } else if (legal == 0) rep->faces_skipped_gate++;
            else rep->faces_skipped_material++;
        }
    }
    if (out.nf == 0) {
        fprintf(stderr, "[solidify-projective] no locally supported faces\n");
        goto done;
    }
    rep->nv_out = out.nv; rep->nf_out = out.nf;
    rep->du = cfg->lattice_du; rep->dv = cfg->lattice_dv;
    rep->n_slices = (int)(nrows > INT_MAX ? INT_MAX : nrows);
    rep->cols_out = 0; /* absolute sparse addresses: no crop-relative width */
    rep->t_total = ves_clock_sec() - t0;
    rc = qs_write_projective_outputs(&out, cfg, in_vmesh, out_stem, rep,
                                      fill_radius, nrows);
    fprintf(stderr,
            "[solidify-projective] %s: %zu copied input samples + %zu "
            "bounded fills -> %zu vertices, %zu local faces; duplicates=%zu, "
            "skipped gate/material=%zu/%zu\n",
            rc == 0 ? "wrote" : "FAILED", rep->observed, rep->filled,
            rep->nv_out, rep->nf_out, rep->duplicates,
            rep->faces_skipped_gate, rep->faces_skipped_material);
done:
    free(entry); free(direct); free(by_col); free(source_out);
    qs_projective_out_dispose(&out);
    return rc;
}

/* ---- public ---------------------------------------------------------------- */

int QuadribbonSolidify_complete(const char *out_stem)
{
    char path[QS_PATH];
    uint64_t nv = 0;
    int projective_v5 = 0;
    static const size_t esize[] = { 1, 1, 4, 4, 4, 4, 4 };   /* + _atlas_column.i32 */
    snprintf(path, sizeof path, "%s.vmesh", out_stem);
    if (!MeshBin_looks_complete(path)) return 0;
    nv = qs_vmesh_nv(path);
    if (nv == 0) return 0;
    for (size_t i = 0; i < sizeof qs_sidecar_suffix / sizeof qs_sidecar_suffix[0]; i++) {
        long long sz = 0;
        snprintf(path, sizeof path, "%s%s", out_stem, qs_sidecar_suffix[i]);
        sz = qs_file_size(path);
        if (sz < 0 || (unsigned long long)sz != nv * esize[i]) return 0;
    }
    snprintf(path, sizeof path, "%s_report.json", out_stem);
    if (qs_file_size(path) <= 0) return 0;
    projective_v5 = qs_file_has_text(
        path, "vesuvius-quadribbon-solid-projective-v5");
    if (qs_file_has_text(path,
                         "vesuvius-quadribbon-solid-projective-v3") ||
        qs_file_has_text(path,
                         "vesuvius-quadribbon-solid-projective-v4") ||
        projective_v5) {
        snprintf(path, sizeof path, "%s_lattice_address.i32", out_stem);
        if (qs_file_size(path) != (long long)(nv * 2 * sizeof(int32_t)))
            return 0;
        snprintf(path, sizeof path, "%s_dependency_kind.u8", out_stem);
        if (qs_file_size(path) != (long long)nv) return 0;
        snprintf(path, sizeof path, "%s_dependency_endpoints.i32", out_stem);
        if (qs_file_size(path) != (long long)(nv * 4 * sizeof(int32_t)))
            return 0;
        if (projective_v5) {
            snprintf(path, sizeof path, "%s_claimant_chart.i32", out_stem);
            if (qs_file_size(path) !=
                (long long)(nv * sizeof(int32_t)))
                return 0;
        }
    }
    return 1;
}

int QuadribbonSolidify_run(Arena_T arena, const QsConfig *cfg,
                           const char *in_vmesh, const char *out_stem,
                           QsReport *rep)
{
    QsInput in;
    QsLattice L;
    QsCtx ctx;
    QsReport local;
    double t0 = ves_clock_sec(), t1 = 0.0;
    int rc = -1;
    if (rep == NULL) rep = &local;
    memset(rep, 0, sizeof *rep);
    memset(&ctx, 0, sizeof ctx);
    if (arena == NULL || cfg == NULL || in_vmesh == NULL || out_stem == NULL) return -1;
    fprintf(stderr, "[solidify] %s -> %s.vmesh (axis %.1f,%.1f pitch %.3f)\n",
            in_vmesh, out_stem, cfg->axis_y, cfg->axis_x, cfg->pitch);
    if (qs_load_input(arena, in_vmesh, &in, rep) != 0) return -1;
    if (cfg->coordinate_warp) for (size_t i=0;i<in.mesh.nv;i++) {
        double point[3]={in.mesh.verts[3*i],in.mesh.verts[3*i+1],in.mesh.verts[3*i+2]},world[3];
        uint32_t flags=0;
        if (AxisWarp_to_world(cfg->coordinate_warp,point,world,&flags)!=0) {
            fprintf(stderr,"[solidify] input coordinate %zu refused (flags %u)\n",i,(unsigned)flags); return -1;
        }
    }
    t1 = ves_clock_sec(); rep->t_load = t1 - t0;
    if (cfg->projective_local)
        return qs_projective_run(cfg, &in, in_vmesh, out_stem, rep, t0);
    if (qs_recover_lattice(arena, &in, &L, rep) != 0) return -1;
    rep->t_lattice = ves_clock_sec() - t1;
    ctx.cfg = cfg; ctx.in = &in; ctx.L = &L; ctx.rep = rep;
    ctx.ncols = (size_t)(L.col_max - L.col_min + 1);
    ctx.img = (uint8_t *)ARENA_CALLOC(arena, ctx.ncols * (size_t)L.n_slices, 1);
    /* runs = maximal consecutive occupied columns */
    {
        uint8_t *occ = (uint8_t *)ARENA_CALLOC(arena, ctx.ncols, 1);
        size_t cells_total = 0;
        int nr = 0;
        for (size_t i = 0; i < in.mesh.nv; i++) occ[(size_t)(L.col[i] - L.col_min)] = 1;
        {   /* LAYOUT: a run spans unoccupied gaps up to run_merge_gap columns
             * when the SAME material lineage sits on both sides (a hole in a
             * continuous sheet); the certificate packs islands with 2-column
             * gaps too, and those keep their own runs (measured 2026-09-02:
             * a lineage-blind merge fused the whole atlas into one run) */
            int32_t *colmat = (int32_t *)ARENA_ALLOC(arena, ctx.ncols * sizeof *colmat);
            int32_t *cm0 = (int32_t *)ARENA_ALLOC(arena, ctx.ncols * sizeof *cm0);
            int32_t *cm1 = (int32_t *)ARENA_ALLOC(arena, ctx.ncols * sizeof *cm1);
            size_t *cn0 = (size_t *)ARENA_CALLOC(arena, ctx.ncols, sizeof *cn0);
            size_t *cn1 = (size_t *)ARENA_CALLOC(arena, ctx.ncols, sizeof *cn1);
            size_t last_occ = SIZE_MAX, merged = 0;
            for (size_t c = 0; c < ctx.ncols; c++) { colmat[c] = -1; cm0[c] = -1; cm1[c] = -1; }
            if (in.material != NULL)
                for (size_t i = 0; i < in.mesh.nv; i++) {
                    size_t c = (size_t)(L.col[i] - L.col_min);
                    int32_t mm = in.material[i];
                    if (mm < 0) continue;
                    if (cm0[c] < 0 || cm0[c] == mm) { cm0[c] = mm; cn0[c]++; }
                    else if (cm1[c] < 0 || cm1[c] == mm) { cm1[c] = mm; cn1[c]++; }
                }
            for (size_t c = 0; c < ctx.ncols; c++) colmat[c] = cn1[c] > cn0[c] ? cm1[c] : cm0[c];
            {
                /* per (column, row): lifted phase and radius of the band-0
                 * observation, for the continuity test across a gap */
                size_t HS = (size_t)L.n_slices;
                float *cph = (float *)ARENA_ALLOC(arena, ctx.ncols * HS * sizeof *cph);
                float *crr = (float *)ARENA_ALLOC(arena, ctx.ncols * HS * sizeof *crr);
                double dtheta_sign = 1.0;
                size_t rejected_phase = 0;
                for (size_t i = 0; i < ctx.ncols * HS; i++) { cph[i] = NAN; crr[i] = 0.0f; }
                for (size_t i = 0; i < in.mesh.nv; i++) {
                    size_t cc = (size_t)(L.col[i] - L.col_min), rr = (size_t)L.slice[i];
                    if (L.band[i] != 0 || rr >= HS) continue;
                    if (!isfinite((double)in.phase[i])) continue;
                    cph[cc * HS + rr] = in.phase[i];
                    crr[cc * HS + rr] = (float)qs_radius(&ctx, in.mesh.verts + i * 3);
                }
                {   /* the winding direction: sign of the median phase step */
                    double pos_n = 0.0, neg_n = 0.0;
                    for (size_t cc = 0; cc + 1 < ctx.ncols; cc++)
                        for (size_t rr = 0; rr < HS; rr++) {
                            double d = (double)cph[(cc + 1) * HS + rr] - (double)cph[cc * HS + rr];
                            if (!isfinite(d)) continue;
                            if (d > 0.0) pos_n += 1.0; else if (d < 0.0) neg_n += 1.0;
                        }
                    dtheta_sign = neg_n > pos_n ? -1.0 : 1.0;
                }
                for (size_t c = 0; c < ctx.ncols; c++) {
                    if (!occ[c]) continue;
                    if (last_occ != SIZE_MAX && cfg->run_merge_gap > 0 &&
                        c - last_occ - 1 <= (size_t)cfg->run_merge_gap &&
                        colmat[c] >= 0 && colmat[c] == colmat[last_occ]) {
                        size_t comparable = 0, agree = 0;
                        for (size_t rr = 0; rr < HS; rr++) {
                            double pa = cph[last_occ * HS + rr], pb = cph[c * HS + rr];
                            double ra = crr[last_occ * HS + rr], rb = crr[c * HS + rr];
                            double want = 0.0, have = 0.0, r = 0.5 * (ra + rb);
                            if (!isfinite(pa) || !isfinite(pb) || !(r > 1.0)) continue;
                            comparable++;
                            want = dtheta_sign * (double)(c - last_occ) * L.du / r;
                            have = pb - pa;
                            if (fabs(have - want) <= QS_ORDER_TOL_FRAC * fabs(want) + QS_ORDER_TOL_RAD) agree++;
                        }
                        if (comparable >= 4 && (double)agree >= QS_MERGE_PHASE_AGREE * (double)comparable) {
                            for (size_t k = last_occ + 1; k < c; k++) occ[k] = 1;
                            merged++;
                        } else {
                            rejected_phase++;
                        }
                    }
                    last_occ = c;
                }
                if (merged > 0 || rejected_phase > 0)
                    fprintf(stderr, "[solidify] layout: %zu same-lineage, phase-continuous column gap(s) bridged inside runs; "
                            "%zu same-lineage gap(s) kept as run breaks (phase jumps)\n", merged, rejected_phase);
            }
        }
        for (size_t c = 0; c < ctx.ncols; c++)
            if (occ[c] && (c == 0 || !occ[c - 1])) nr++;
        ctx.runs = (QsRunStat *)ARENA_CALLOC(arena, nr > 0 ? (size_t)nr : 1, sizeof *ctx.runs);
        ctx.n_runs = nr;
        {
            int k = 0;
            for (size_t c = 0; c < ctx.ncols; c++) {
                if (occ[c] && (c == 0 || !occ[c - 1])) {
                    size_t e = c;
                    while (e + 1 < ctx.ncols && occ[e + 1]) e++;
                    ctx.runs[k].c0 = L.col_min + (long)c;
                    ctx.runs[k].c1 = L.col_min + (long)e;
                    cells_total += (e - c + 1) * (size_t)L.n_slices;
                    k++;
                }
            }
        }
        rep->n_runs = nr;
        ctx.out.cap_v = cells_total;
        ctx.out.cap_f = cells_total * 2 + 1;
        ctx.out.verts = (float *)ARENA_ALLOC(arena, (cells_total + 1) * 3 * sizeof(float));
        ctx.out.uv = (float *)ARENA_ALLOC(arena, (cells_total + 1) * 2 * sizeof(float));
        ctx.out.faces = (int32_t *)ARENA_ALLOC(arena, ctx.out.cap_f * 3 * sizeof(int32_t));
        ctx.out.support = (uint8_t *)ARENA_ALLOC(arena, cells_total + 1);
        ctx.out.prov = (uint8_t *)ARENA_ALLOC(arena, cells_total + 1);
        ctx.out.phase = (float *)ARENA_ALLOC(arena, (cells_total + 1) * sizeof(float));
        ctx.out.material = (int32_t *)ARENA_ALLOC(arena, (cells_total + 1) * sizeof(int32_t));
        ctx.out.run = (int32_t *)ARENA_ALLOC(arena, (cells_total + 1) * sizeof(int32_t));
        ctx.out.lane = (int32_t *)ARENA_ALLOC(arena, (cells_total + 1) * sizeof(int32_t));
        ctx.out.atlas_col = (int32_t *)ARENA_ALLOC(arena, (cells_total + 1) * sizeof(int32_t));
        fprintf(stderr, "[solidify] %d column run(s), %zu rectangle cells, %zu atlas "
                        "columns x %d slices\n", nr, cells_total, ctx.ncols, L.n_slices);
    }
    t1 = ves_clock_sec();
    ctx.col_out = L.col_min;     /* LAYOUT: the sheet starts where the atlas starts */
    {
        static ZarrU8 zarr;      /* one chunk cache per process; the stage is serial */
        if (cfg->raw_zarr != NULL && cfg->raw_zarr[0] != '\0') {   /* the DARK gate always probes */
            if (ZarrU8_open(&zarr, cfg->raw_zarr) == 0) {
                ctx.zarr = &zarr;
                fprintf(stderr, "[solidify] RAW %s armed (chunk %ld): fill snap %s, dark-cell snap %s (reach %.1f vox, ridge >= %d)\n",
                        cfg->raw_zarr, zarr.chunk[0], QS_FILL_SNAP_ENABLE ? "ON" : "off",
                        QS_DARK_CELL_SNAP_ENABLE ? "ON" : "off", QS_DARK_SNAP_REACH, QS_FILL_SNAP_MIN_U8);
            } else {
                fprintf(stderr, "[solidify] WARNING: cannot open RAW zarr %s -- fills are NOT snapped\n",
                        cfg->raw_zarr);
            }
        }
    }
    for (int r = 0; r < ctx.n_runs; r++) {
        if (qs_process_run(arena, &ctx, r, ctx.runs[r].c0, ctx.runs[r].c1) != 0) {
            fprintf(stderr, "[solidify] run %d failed\n", r);
            free(ctx.dup.v);
            free(ctx.snap.v);
            if (ctx.zarr != NULL) ZarrU8_close(ctx.zarr);
            return -1;
        }
    }
    rep->t_build = ves_clock_sec() - t1;
    for (int r = 0; r < ctx.n_runs; r++) {
        const QsRunStat *s = &ctx.runs[r];
        rep->observed += s->observed;
        rep->band1_adopted += s->adopted;
        rep->band1_rejected += s->rejected;
        rep->band1_no_reference += s->no_ref;
        rep->duplicates += s->duplicates;
        rep->filled += s->filled;
        rep->dropped_unbridged += s->dropped_unbridged;
        rep->dropped_skip += s->dropped_skip;
        rep->dropped_order += s->dropped_order;
        rep->dropped_confetti += s->dropped_confetti;
        rep->dropped_chord += s->dropped_chord;
        rep->dropped_material += s->dropped_material;
        rep->dropped_bbox += s->dropped_bbox;
        rep->dropped_core += s->dropped_core;
        rep->dropped_unreferenced += s->dropped_unref;
        rep->faces_skipped_gate += s->skipped_gate;
        rep->faces_skipped_material += s->skipped_material;
        rep->holes += s->holes;
        rep->holes_refused += s->holes_refused;
        rep->holes_filled += s->holes_filled;
        rep->snap_candidates += s->snap_candidates;
        rep->snap_moved += s->snap_moved;
        rep->snap_no_ridge += s->snap_no_ridge;
        rep->fold_cells += s->fold_cells;
        rep->fold_columns += s->fold_columns;
        rep->dark_cells += s->dark_cells;
        rep->dark_snapped += s->dark_snapped;
        rep->dark_unrecovered += s->dark_unrecovered;
        rep->hole_cells += s->hole_cells;
        rep->refused_hole_cells += s->refused_hole_cells;
        rep->midline_probed += s->midline_probed;
        rep->midline_outliers += s->midline_outliers;
        rep->midline_isolated += s->midline_isolated;
        rep->midline_dense += s->midline_dense;
        if (isfinite(s->ramp_pitches_per_turn)) {
            rep->ramp_runs++;
            rep->ramp_cells += s->observed;
            if (s->ramp_pitches_per_turn > QS_GATE_RAMP_PITCHES) {
                rep->ramp_runs_over++;
                rep->ramp_cells_over += s->observed;
            }
            if (s->ramp_pitches_per_turn > rep->ramp_max)
                rep->ramp_max = s->ramp_pitches_per_turn;
        }
        rep->scale_faces += s->scale_faces;
        rep->scale_over += s->scale_over;
        rep->scale_over_fill += s->scale_over_fill;
        rep->dropped_fill_far += s->dropped_fill_far;
        rep->obs_smooth_cells += s->obs_smooth_cells;
        rep->column_probed += s->column_probed;
        rep->column_jumps += s->column_jumps;
        rep->depth_probed += s->depth_probed;
        rep->depth_outliers += s->depth_outliers;
        rep->ridge_cells += s->ridge_cells;
        rep->ridge_moved += s->ridge_moved;
        rep->ridge_on_before += s->ridge_on_before;
        rep->ridge_on_after += s->ridge_on_after;
        rep->ridge_fill_cells += s->ridge_fill_cells;
        rep->ridge_fill_on_before += s->ridge_fill_on_before;
        rep->ridge_fill_on_after += s->ridge_fill_on_after;
        if (s->ridge_shift_max > rep->ridge_shift_max) rep->ridge_shift_max = s->ridge_shift_max;
        rep->fill_cells += s->fill_cells;
        rep->fill_far += s->fill_far;
        if (s->fill_dist_max > rep->fill_dist_max) rep->fill_dist_max = s->fill_dist_max;
    }
    if (ctx.obs_smooth.n > 0) {
        qsort(ctx.obs_smooth.v, ctx.obs_smooth.n, sizeof *ctx.obs_smooth.v, qs_cmp_double);
        rep->obs_smooth_p50 = qs_quantile(ctx.obs_smooth.v, ctx.obs_smooth.n, 0.50);
        rep->obs_smooth_p90 = qs_quantile(ctx.obs_smooth.v, ctx.obs_smooth.n, 0.90);
        rep->obs_smooth_max = ctx.obs_smooth.v[ctx.obs_smooth.n - 1];
    }
    free(ctx.obs_smooth.v);
    if (ctx.column.n > 0) {
        qsort(ctx.column.v, ctx.column.n, sizeof *ctx.column.v, qs_cmp_double);
        rep->column_dev_p50 = qs_quantile(ctx.column.v, ctx.column.n, 0.50);
        rep->column_dev_p90 = qs_quantile(ctx.column.v, ctx.column.n, 0.90);
    }
    free(ctx.column.v);
    if (ctx.depth.n > 0) {
        qsort(ctx.depth.v, ctx.depth.n, sizeof *ctx.depth.v, qs_cmp_double);
        rep->depth_hp_p50 = qs_quantile(ctx.depth.v, ctx.depth.n, 0.50);
        rep->depth_hp_p90 = qs_quantile(ctx.depth.v, ctx.depth.n, 0.90);
    }
    free(ctx.depth.v);
    if (ctx.ridge.n > 0) {
        qsort(ctx.ridge.v, ctx.ridge.n, sizeof *ctx.ridge.v, qs_cmp_double);
        rep->ridge_shift_p50 = qs_quantile(ctx.ridge.v, ctx.ridge.n, 0.50);
        rep->ridge_shift_p90 = qs_quantile(ctx.ridge.v, ctx.ridge.n, 0.90);
    }
    free(ctx.ridge.v);
    if (ctx.ramp.n > 0) {
        qsort(ctx.ramp.v, ctx.ramp.n, sizeof *ctx.ramp.v, qs_cmp_double);
        rep->ramp_p50 = qs_quantile(ctx.ramp.v, ctx.ramp.n, 0.50);
        rep->ramp_p90 = qs_quantile(ctx.ramp.v, ctx.ramp.n, 0.90);
    }
    free(ctx.ramp.v);
    if (ctx.scale.n > 0) {
        qsort(ctx.scale.v, ctx.scale.n, sizeof *ctx.scale.v, qs_cmp_double);
        rep->scale_p50 = qs_quantile(ctx.scale.v, ctx.scale.n, 0.50);
        rep->scale_p90 = qs_quantile(ctx.scale.v, ctx.scale.n, 0.90);
        rep->scale_p99 = qs_quantile(ctx.scale.v, ctx.scale.n, 0.99);
    }
    if (ctx.scale_obs.n > 0) {
        qsort(ctx.scale_obs.v, ctx.scale_obs.n, sizeof *ctx.scale_obs.v, qs_cmp_double);
        rep->scale_obs_p90 = qs_quantile(ctx.scale_obs.v, ctx.scale_obs.n, 0.90);
    }
    if (ctx.scale_fill.n > 0) {
        qsort(ctx.scale_fill.v, ctx.scale_fill.n, sizeof *ctx.scale_fill.v, qs_cmp_double);
        rep->scale_fill_p90 = qs_quantile(ctx.scale_fill.v, ctx.scale_fill.n, 0.90);
    }
    if (ctx.fill_dist.n > 0) {
        qsort(ctx.fill_dist.v, ctx.fill_dist.n, sizeof *ctx.fill_dist.v, qs_cmp_double);
        rep->fill_dist_p50 = qs_quantile(ctx.fill_dist.v, ctx.fill_dist.n, 0.50);
        rep->fill_dist_p90 = qs_quantile(ctx.fill_dist.v, ctx.fill_dist.n, 0.90);
    }
    free(ctx.scale.v); free(ctx.scale_obs.v); free(ctx.scale_fill.v); free(ctx.fill_dist.v);
    if (ctx.midline.n > 0) {
        qsort(ctx.midline.v, ctx.midline.n, sizeof *ctx.midline.v, qs_cmp_double);
        rep->midline_p50 = qs_quantile(ctx.midline.v, ctx.midline.n, 0.50);
        rep->midline_p90 = qs_quantile(ctx.midline.v, ctx.midline.n, 0.90);
    }
    free(ctx.midline.v);
    {   /* SOLID GATES: the artifact classes the reviewer boxed, as fractions */
        double obs_all = (double)(rep->observed + rep->fold_cells);
        double f_fold = obs_all > 0.0 ? (double)rep->fold_cells / obs_all : 0.0;
        double f_dark = obs_all > 0.0 ? (double)rep->midline_isolated / obs_all : 0.0;   /* isolated suspects */
        double f_dense = obs_all > 0.0 ? (double)rep->midline_dense / obs_all : 0.0;
        double f_ref = rep->hole_cells > 0 ? (double)rep->refused_hole_cells / (double)rep->hole_cells : 0.0;
        double f_scale = rep->scale_faces > 0 ? (double)rep->scale_over / (double)rep->scale_faces : 0.0;
        double f_far = rep->fill_cells > 0 ? (double)rep->fill_far / (double)rep->fill_cells : 0.0;
        double f_ramp = rep->ramp_cells > 0 ? (double)rep->ramp_cells_over / (double)rep->ramp_cells : 0.0;
        int fails = 0;
        if (f_fold > QS_GATE_FOLD_MAX) fails++;
        if (f_dark > QS_GATE_DARK_MAX) fails++;
        if (f_ref > QS_GATE_REFUSED_HOLE_MAX) fails++;
        if (f_scale > QS_GATE_SCALE_MAX) fails++;
        if (f_far > QS_GATE_FILL_FAR_MAX) fails++;
        if (f_ramp > QS_GATE_RAMP_MAX) fails++;
        rep->gates_failed = fails;
        fprintf(stderr, "[solidify GATES] FOLD %.3f%% of observed cells demoted (max %.1f%%) %s | MIDLINE isolated "
                "suspects %.3f%% of observed cells (dark and > %.1f vox from any papyrus, or none in reach; max %.1f%%; "
                "%zu refit) %s, dense suspects %.1f%% left in place (%zu; the fit's depth error, reported); kept cells' "
                "distance to the papyrus face p50 %.2f p90 %.2f vox | "
                "REFUSED-HOLES %.1f%% of hole cells refused (max %.0f%%) %s -> %s\n",
                100.0 * f_fold, 100.0 * QS_GATE_FOLD_MAX, f_fold > QS_GATE_FOLD_MAX ? "FAIL" : "ok",
                100.0 * f_dark, QS_MIDLINE_OUTLIER_VOX, 100.0 * QS_GATE_DARK_MAX, rep->midline_isolated,
                f_dark > QS_GATE_DARK_MAX ? "FAIL" : "ok", 100.0 * f_dense, rep->midline_dense,
                rep->midline_p50, rep->midline_p90,
                100.0 * f_ref, 100.0 * QS_GATE_REFUSED_HOLE_MAX, f_ref > QS_GATE_REFUSED_HOLE_MAX ? "FAIL" : "ok",
                fails == 0 ? "PASS" : "NOT_SHIPPABLE");
        fprintf(stderr, "[solidify GATES] SCALE %.2f%% of %zu emitted triangles differ from their (u,v) area by "
                "> %.2f octaves (max %.0f%%; %.0f%% of those touch a fill) %s | |log2 area ratio| p50 %.3f p90 %.3f "
                "p99 %.3f (observed p90 %.3f, fill p90 %.3f) | FILL-SUPPORT %.1f%% of %zu fill cells are more than "
                "%d lattice cells from an observation (max %.0f%%; distance p50 %.0f p90 %.0f max %zu) %s\n",
                100.0 * f_scale, rep->scale_faces, QS_GATE_SCALE_LOG2, 100.0 * QS_GATE_SCALE_MAX,
                rep->scale_over > 0 ? 100.0 * (double)rep->scale_over_fill / (double)rep->scale_over : 0.0,
                f_scale > QS_GATE_SCALE_MAX ? "FAIL" : "ok",
                rep->scale_p50, rep->scale_p90, rep->scale_p99, rep->scale_obs_p90, rep->scale_fill_p90,
                100.0 * f_far, rep->fill_cells, QS_FILL_SUPPORT_CELLS, 100.0 * QS_GATE_FILL_FAR_MAX,
                rep->fill_dist_p50, rep->fill_dist_p90, rep->fill_dist_max,
                f_far > QS_GATE_FILL_FAR_MAX ? "FAIL" : "ok");
        if (rep->obs_smooth_cells > 0)
            fprintf(stderr, "[solidify] DEPTH SMOOTHING: %zu cells moved radially (bound %.2f vox, "
                    "per-run |move| p50 %.3f p90 %.3f max %.3f vox)\n",
                    rep->obs_smooth_cells, (double)QS_OBS_SMOOTH_VOX,
                    rep->obs_smooth_p50, rep->obs_smooth_p90, rep->obs_smooth_max);
        if (rep->column_probed > 0)
            fprintf(stderr, "[solidify] COLUMN CONTINUITY: %zu of %zu observed cells jump more than "
                    "%.2f pitches from their column's local median and were demoted (%.2f%%); per-run "
                    "|radius - column median| p50 %.3f p90 %.3f pitches\n",
                    rep->column_jumps, rep->column_probed, (double)QS_COLUMN_JUMP_PITCH,
                    100.0 * (double)rep->column_jumps / (double)rep->column_probed,
                    rep->column_dev_p50, rep->column_dev_p90);
        if (rep->depth_probed > 0)
            fprintf(stderr, "[solidify] DEPTH OUTLIERS: %zu of %zu observed cells more than %.2f vox "
                    "from their 9-cell local median, demoted to fill (%.2f%%); per-run |radius - median| "
                    "p50 %.2f p90 %.2f vox\n",
                    rep->depth_outliers, rep->depth_probed, (double)QS_DEPTH_OUTLIER_VOX,
                    100.0 * (double)rep->depth_outliers / (double)rep->depth_probed,
                    rep->depth_hp_p50, rep->depth_hp_p90);
        if (rep->ridge_cells > 0)
            fprintf(stderr, "[solidify] RIDGE TRACK: %zu cells, %zu moved (%.1f%%); on papyrus %.1f%% -> %.1f%%; "
                    "per-run |shift| p50 %.2f p90 %.2f, max %.2f vox\n",
                    rep->ridge_cells, rep->ridge_moved,
                    100.0 * (double)rep->ridge_moved / (double)rep->ridge_cells,
                    100.0 * (double)rep->ridge_on_before / (double)rep->ridge_cells,
                    100.0 * (double)rep->ridge_on_after / (double)rep->ridge_cells,
                    rep->ridge_shift_p50, rep->ridge_shift_p90, rep->ridge_shift_max);
        if (rep->ridge_fill_cells > 0)
            fprintf(stderr, "[solidify] RIDGE TRACK fills: %zu fill cells searched to %.1f vox; on papyrus %zu -> %zu (%.1f%% -> %.1f%%)\n",
                    rep->ridge_fill_cells, (double)QS_RIDGE_FILL_REACH, rep->ridge_fill_on_before, rep->ridge_fill_on_after,
                    100.0 * (double)rep->ridge_fill_on_before / (double)rep->ridge_fill_cells,
                    100.0 * (double)rep->ridge_fill_on_after / (double)rep->ridge_fill_cells);
        fprintf(stderr, "[solidify GATES] RAMP %.1f%% of observed cells sit in runs that gain more than %.1f pitches of "
                "radius per turn (max %.0f%%; %zu of %zu runs; pitches/turn p50 %.2f p90 %.2f max %.1f -- one is correct) %s\n",
                100.0 * f_ramp, QS_GATE_RAMP_PITCHES, 100.0 * QS_GATE_RAMP_MAX, rep->ramp_runs_over, rep->ramp_runs,
                rep->ramp_p50, rep->ramp_p90, rep->ramp_max, f_ramp > QS_GATE_RAMP_MAX ? "FAIL" : "ok");
    }
    rep->cols_out = ctx.col_out - L.col_min;
    if (ctx.snap.n > 0) {
        qsort(ctx.snap.v, ctx.snap.n, sizeof *ctx.snap.v, qs_cmp_double);
        rep->snap_shift_p50 = qs_quantile(ctx.snap.v, ctx.snap.n, 0.50);
        rep->snap_shift_p90 = qs_quantile(ctx.snap.v, ctx.snap.n, 0.90);
    }
    free(ctx.snap.v);
    if (ctx.zarr != NULL) {
        fprintf(stderr, "[solidify] fill snap: %zu of %zu fill cells moved (|shift| p50 %.2f p90 %.2f vox), "
                "%zu without a ridge within %.1f vox; zarr chunks read %zu (missing %zu)\n",
                rep->snap_moved, rep->snap_candidates, rep->snap_shift_p50, rep->snap_shift_p90,
                rep->snap_no_ridge, QS_FILL_SNAP_REACH, ctx.zarr->misses, ctx.zarr->missing);
        ZarrU8_close(ctx.zarr);
        ctx.zarr = NULL;
    }
    if (ctx.dup.n > 0) {
        qsort(ctx.dup.v, ctx.dup.n, sizeof *ctx.dup.v, qs_cmp_double);
        rep->dup_dr_pitch_p10 = qs_quantile(ctx.dup.v, ctx.dup.n, 0.10);
        rep->dup_dr_pitch_p50 = qs_quantile(ctx.dup.v, ctx.dup.n, 0.50);
        rep->dup_dr_pitch_p90 = qs_quantile(ctx.dup.v, ctx.dup.n, 0.90);
    }
    free(ctx.dup.v);
    rep->nv_out = ctx.out.nv;
    rep->nf_out = ctx.out.nf;
    t1 = ves_clock_sec();
    if (ctx.ct_coordinate_refusals) {
        fprintf(stderr,"[solidify] refused %zu CT coordinate conversions; no solidified mesh published\n",ctx.ct_coordinate_refusals);
        return -1;
    }
    rc = qs_write_outputs(&ctx, out_stem, in_vmesh);
    rep->t_write = ves_clock_sec() - t1;
    rep->t_total = ves_clock_sec() - t0;
    fprintf(stderr,
            "[solidify] %s: %zu verts (%zu observed, %zu adopted band-1, %zu filled), "
            "%zu faces; band-1 rejected %zu (+%zu without reference), duplicates %zu "
            "(|dr| p10/50/90 = %.2f/%.2f/%.2f pitch); not filled: unbridged %zu, wrap-skip %zu, "
            "misordered %zu, chord-gate %zu, mixed-material %zu, box %zu, core %zu, unreferenced "
            "%zu; confetti "
            "runs not emitted: %zu observed cells; holes %zu (refused %zu, filled by vote %zu); "
            "layout %ld columns (atlas %zu); skipped "
            "faces gate/material = %zu/%zu; phase mismatch max %.3f rad (%zu over limit); "
            "%.1fs\n",
            rc == 0 ? "wrote" : "FAILED", rep->nv_out, rep->observed, rep->band1_adopted,
            rep->filled, rep->nf_out, rep->band1_rejected, rep->band1_no_reference,
            rep->duplicates, rep->dup_dr_pitch_p10, rep->dup_dr_pitch_p50,
            rep->dup_dr_pitch_p90, rep->dropped_unbridged, rep->dropped_skip, rep->dropped_order,
            rep->dropped_chord, rep->dropped_material, rep->dropped_bbox, rep->dropped_core,
            rep->dropped_unreferenced, rep->dropped_confetti, rep->holes, rep->holes_refused,
            rep->holes_filled, rep->cols_out, ctx.ncols, rep->faces_skipped_gate,
            rep->faces_skipped_material,
            rep->phase_mismatch_max, rep->phase_mismatch_over_limit, rep->t_total);
    return rc;
}

/* ============================================================================
 * Self-test: a synthetic multi-run lattice on a cylinder (axis at the origin)
 * with an interior hole, peel-band duplicates, adopted/rejected band-1 cells,
 * a cap-sized notch, a wrap-gate jump, a material border and an umbilicus-core
 * run; the output is read back through the public files.
 * ==========================================================================*/

typedef struct QsTv {
    size_t n, cap;
    int32_t *col, *k, *band;
    float *pos;       /* z,y,x */
    float *phase;
    uint8_t *support, *prov;
    int32_t *mat, *lane;
} QsTv;

static int qs_tv_push(QsTv *t, long col, int k, int band, double z, double y,
                      double x, double phase, uint8_t support, uint8_t prov,
                      int32_t mat, int32_t lane)
{
    if (t->n == t->cap) {
        size_t nc = t->cap ? t->cap * 2 : 1024;
        int32_t *c2 = (int32_t *)realloc(t->col, nc * sizeof *c2);
        int32_t *k2 = (int32_t *)realloc(t->k, nc * sizeof *k2);
        int32_t *b2 = (int32_t *)realloc(t->band, nc * sizeof *b2);
        float *p2 = (float *)realloc(t->pos, nc * 3 * sizeof *p2);
        float *ph2 = (float *)realloc(t->phase, nc * sizeof *ph2);
        uint8_t *s2 = (uint8_t *)realloc(t->support, nc);
        uint8_t *pr2 = (uint8_t *)realloc(t->prov, nc);
        int32_t *m2 = (int32_t *)realloc(t->mat, nc * sizeof *m2);
        int32_t *l2 = (int32_t *)realloc(t->lane, nc * sizeof *l2);
        if (c2) t->col = c2;
        if (k2) t->k = k2;
        if (b2) t->band = b2;
        if (p2) t->pos = p2;
        if (ph2) t->phase = ph2;
        if (s2) t->support = s2;
        if (pr2) t->prov = pr2;
        if (m2) t->mat = m2;
        if (l2) t->lane = l2;
        if (!c2 || !k2 || !b2 || !p2 || !ph2 || !s2 || !pr2 || !m2 || !l2) return -1;
        t->cap = nc;
    }
    t->col[t->n] = (int32_t)col;
    t->k[t->n] = k;
    t->band[t->n] = band;
    t->pos[t->n * 3] = (float)z;
    t->pos[t->n * 3 + 1] = (float)y;
    t->pos[t->n * 3 + 2] = (float)x;
    t->phase[t->n] = (float)phase;
    t->support[t->n] = support;
    t->prov[t->n] = prov;
    t->mat[t->n] = mat;
    t->lane[t->n] = lane;
    t->n++;
    return 0;
}

static void qs_tv_free(QsTv *t)
{
    free(t->col); free(t->k); free(t->band); free(t->pos); free(t->phase);
    free(t->support); free(t->prov); free(t->mat); free(t->lane);
    memset(t, 0, sizeof *t);
}

/* Write the synthetic ribbon as a VMESH + sidecars (lattice quads between
 * neighbouring cells of one band). */
static int qs_tv_write(const QsTv *t, const char *stem, double du, double dv,
                       int np, int write_stats)
{
    char path[QS_PATH];
    size_t n = t->n, nf = 0;
    float *uv = (float *)malloc(n * 2 * sizeof *uv);
    int32_t *faces = (int32_t *)malloc((n * 2 + 1) * 3 * sizeof *faces);
    long cmin = LONG_MAX, cmax = LONG_MIN;
    int bmax = 0;
    int32_t *map = NULL;
    size_t ncols = 0, nrows = 0;
    int rc = -1;
    if (uv == NULL || faces == NULL) { free(uv); free(faces); return -1; }
    for (size_t i = 0; i < n; i++) {
        if (t->col[i] < cmin) cmin = t->col[i];
        if (t->col[i] > cmax) cmax = t->col[i];
        if (t->band[i] > bmax) bmax = t->band[i];
        uv[i * 2] = (float)((double)t->col[i] * du);
        uv[i * 2 + 1] = (float)((double)(t->band[i] * np + t->k[i]) * dv);
    }
    ncols = (size_t)(cmax - cmin + 1);
    nrows = (size_t)np * (size_t)(bmax + 1);
    map = (int32_t *)malloc(ncols * nrows * sizeof *map);
    if (map == NULL) { free(uv); free(faces); return -1; }
    for (size_t i = 0; i < ncols * nrows; i++) map[i] = -1;
    for (size_t i = 0; i < n; i++)
        map[(size_t)(t->band[i] * np + t->k[i]) * ncols + (size_t)(t->col[i] - cmin)] = (int32_t)i;
    for (size_t r = 0; r + 1 < nrows; r++) {
        if ((r + 1) % (size_t)np == 0) continue;   /* never across a band */
        for (size_t c = 0; c + 1 < ncols; c++) {
            int32_t a = map[r * ncols + c], b = map[r * ncols + c + 1];
            int32_t cc = map[(r + 1) * ncols + c], d = map[(r + 1) * ncols + c + 1];
            if (a < 0 || b < 0 || cc < 0 || d < 0) continue;
            faces[nf * 3] = a; faces[nf * 3 + 1] = b; faces[nf * 3 + 2] = cc; nf++;
            faces[nf * 3] = b; faces[nf * 3 + 1] = d; faces[nf * 3 + 2] = cc; nf++;
        }
    }
    free(map);
    snprintf(path, sizeof path, "%s_support.u8", stem);
    if (qs_write_raw(path, t->support, n) != 0) goto done;
    snprintf(path, sizeof path, "%s_provenance.u8", stem);
    if (qs_write_raw(path, t->prov, n) != 0) goto done;
    snprintf(path, sizeof path, "%s_phase.f32", stem);
    if (qs_write_raw(path, t->phase, n * sizeof(float)) != 0) goto done;
    snprintf(path, sizeof path, "%s_material_identity.i32", stem);
    if (qs_write_raw(path, t->mat, n * sizeof(int32_t)) != 0) goto done;
    snprintf(path, sizeof path, "%s_claimant_chart.i32", stem);
    if (qs_write_raw(path, t->lane, n * sizeof(int32_t)) != 0) goto done;
    snprintf(path, sizeof path, "%s_lane.i32", stem);
    if (qs_write_raw(path, t->lane, n * sizeof(int32_t)) != 0) goto done;
    snprintf(path, sizeof path, "%s_stats.json", stem);
    if (write_stats) {
        FILE *f = fopen(path, "wb");
        if (f == NULL) goto done;
        fprintf(f, "{ \"slicing\": { \"n_slices\": %d } }\n", np);
        fclose(f);
    } else {
        remove(path);
    }
    snprintf(path, sizeof path, "%s.vmesh", stem);
    if (nf == 0) goto done;
    rc = MeshBin_write(path, t->pos, n, faces, nf, uv);
done:
    free(uv);
    free(faces);
    return rc;
}

static long qs_test_find_address(const int32_t *address, size_t n,
                                 int32_t col, int32_t row)
{
    /* Copy-only deliberately preserves the source member's vertex order;
     * the larger member may therefore be permuted.  This is self-test code on
     * tiny fixtures, so use an order-independent scan rather than assuming the
     * legacy fill path's row-major output order. */
    for (size_t i = 0; i < n; i++)
        if (address[i * 2] == col && address[i * 2 + 1] == row)
            return (long)i;
    return -1;
}

/* Metamorphic checker used by the native self-test: every A-address must
 * exist in B with byte-identical published values and dependency certificate,
 * and every A triangle must exist in B after remapping indices by address. */
static int qs_test_projective_subset(Arena_T arena, const char *a_stem,
                                     const char *b_stem)
{
    MeshBinData a, b;
    char path[QS_PATH];
    const int32_t *aa = NULL, *ba = NULL, *ad = NULL, *bd = NULL;
    const int32_t *am = NULL, *bm = NULL, *ac = NULL, *bc = NULL;
    const int32_t *al = NULL, *bl = NULL;
    const float *ap = NULL, *bp = NULL;
    const uint8_t *as = NULL, *bs = NULL, *av = NULL, *bv = NULL;
    const uint8_t *ak = NULL, *bk = NULL;
    int fails = 0, saw_fill = 0;
    memset(&a, 0, sizeof a); memset(&b, 0, sizeof b);
    snprintf(path, sizeof path, "%s.vmesh", a_stem);
    if (MeshBin_read_arena(arena, path, &a) != 0) return 1;
    snprintf(path, sizeof path, "%s.vmesh", b_stem);
    if (MeshBin_read_arena(arena, path, &b) != 0) return 1;
#define QS_TEST_READ(dst, stem, suffix, count, type) do {                    \
    snprintf(path, sizeof path, "%s" suffix, stem);                        \
    dst = (const type *)qs_read_exact(arena, path, count, sizeof(type));     \
    if (dst == NULL) return 1;                                                \
} while (0)
    QS_TEST_READ(aa, a_stem, "_lattice_address.i32", a.nv * 2, int32_t);
    QS_TEST_READ(ba, b_stem, "_lattice_address.i32", b.nv * 2, int32_t);
    QS_TEST_READ(ad, a_stem, "_dependency_endpoints.i32", a.nv * 4, int32_t);
    QS_TEST_READ(bd, b_stem, "_dependency_endpoints.i32", b.nv * 4, int32_t);
    QS_TEST_READ(am, a_stem, "_material_identity.i32", a.nv, int32_t);
    QS_TEST_READ(bm, b_stem, "_material_identity.i32", b.nv, int32_t);
    QS_TEST_READ(ac, a_stem, "_claimant_chart.i32", a.nv, int32_t);
    QS_TEST_READ(bc, b_stem, "_claimant_chart.i32", b.nv, int32_t);
    QS_TEST_READ(al, a_stem, "_lane.i32", a.nv, int32_t);
    QS_TEST_READ(bl, b_stem, "_lane.i32", b.nv, int32_t);
    QS_TEST_READ(ap, a_stem, "_phase.f32", a.nv, float);
    QS_TEST_READ(bp, b_stem, "_phase.f32", b.nv, float);
    QS_TEST_READ(as, a_stem, "_support.u8", a.nv, uint8_t);
    QS_TEST_READ(bs, b_stem, "_support.u8", b.nv, uint8_t);
    QS_TEST_READ(av, a_stem, "_provenance.u8", a.nv, uint8_t);
    QS_TEST_READ(bv, b_stem, "_provenance.u8", b.nv, uint8_t);
    QS_TEST_READ(ak, a_stem, "_dependency_kind.u8", a.nv, uint8_t);
    QS_TEST_READ(bk, b_stem, "_dependency_kind.u8", b.nv, uint8_t);
#undef QS_TEST_READ
    for (size_t i = 0; i < a.nv; i++) {
        long j = qs_test_find_address(ba, b.nv, aa[i * 2], aa[i * 2 + 1]);
        if (j < 0) { fails++; continue; }
        if (memcmp(a.verts + i * 3, b.verts + (size_t)j * 3,
                   3 * sizeof(float)) != 0 ||
            memcmp(a.uv + i * 2, b.uv + (size_t)j * 2,
                   2 * sizeof(float)) != 0 ||
            memcmp(ap + i, bp + j, sizeof(float)) != 0 ||
            am[i] != bm[j] || ac[i] != bc[j] ||
            al[i] != bl[j] || as[i] != bs[j] ||
            av[i] != bv[j] || ak[i] != bk[j] ||
            memcmp(ad + i * 4, bd + (size_t)j * 4,
                   4 * sizeof(int32_t)) != 0)
            fails++;
        if (ak[i] != QS_DEP_OBSERVED) saw_fill = 1;
    }
    for (size_t f = 0; f < a.nf; f++) {
        int32_t mapped[3];
        int found = 0;
        for (int q = 0; q < 3; q++) {
            size_t vi = (size_t)a.faces[f * 3 + (size_t)q];
            long j = qs_test_find_address(ba, b.nv, aa[vi * 2],
                                          aa[vi * 2 + 1]);
            mapped[q] = (int32_t)j;
        }
        for (size_t g = 0; g < b.nf && !found; g++)
            found = b.faces[g * 3] == mapped[0] &&
                    b.faces[g * 3 + 1] == mapped[1] &&
                    b.faces[g * 3 + 2] == mapped[2];
        if (!found) fails++;
    }
    if (!saw_fill || !(b.nv > a.nv)) fails++;
    return fails;
}

#define QS_CK(cond, msg) do { if (!(cond)) { fprintf(stderr, "  FAIL: %s\n", msg); fails++; } } while (0)

int QuadribbonSolidify_selftest(void)
{

    int fails = 0;
    char tmp[QS_PATH], dir[QS_PATH], in_stem[QS_PATH], out_stem[QS_PATH];
    char refit_stem[QS_PATH], in_vmesh[QS_PATH], out_vmesh[QS_PATH], path[QS_PATH];
    const double du = 2.0, dv = 2.0, slice_h = 2.0, z0 = 100.0, pitch = 9.5;
    const double R_ref = 100.0;
    const int np = 20;
    QsTv tv;
    QsConfig cfg;
    QsReport rep;
    Arena_T arena = NULL;
    memset(&tv, 0, sizeof tv);
    memset(&cfg, 0, sizeof cfg);
    cfg.axis_y = 0.0; cfg.axis_x = 0.0; cfg.pitch = pitch;
    cfg.run_min_observed = 20;   /* the synthetic runs hold 15..800 cells */
    cfg.good_column_frac = 0.10; /* the scene's sparse columns hold 1 of 20 rows */
    cfg.fill_max_u_gap = 16;
    cfg.skip_dr_pitch = 0.5;
    {   /* CT FRAME: straightened -> world must invert the lane's straightening */
        static const double az[3] = { 0.0, 100.0, 200.0 };
        static const double ay[3] = { 10.0, 30.0, 50.0 };
        static const double ax[3] = { 5.0, 5.0, 25.0 };
        QsConfig cf;
        QsCtx cx;
        double p[3], w[3];
        int bad = 0;
        memset(&cf, 0, sizeof cf);
        memset(&cx, 0, sizeof cx);
        cf.axis_curve_z = az; cf.axis_curve_y = ay; cf.axis_curve_x = ax;
        cf.axis_curve_n = 3; cf.axis_ref_y = 10.0; cf.axis_ref_x = 5.0;
        cx.cfg = &cf;
        /* a point straightened from world (50, 40, 15) at z = 50:
         * axis(50) = (20, 5), so straightened = (50, 40 - 20 + 10, 15 - 5 + 5) */
        p[0] = 50.0; p[1] = 30.0; p[2] = 15.0;
        qs_ct_world(&cx, p, w);
        if (fabs(w[0] - 50.0) > 1e-9 || fabs(w[1] - 40.0) > 1e-9 || fabs(w[2] - 15.0) > 1e-9) {
            fprintf(stderr, "  qs_ct_world -> (%.3f, %.3f, %.3f), want (50, 40, 15)\n", w[0], w[1], w[2]);
            bad++;
        }
        cf.axis_curve_n = 0;                     /* unarmed = identity */
        qs_ct_world(&cx, p, w);
        if (fabs(w[1] - p[1]) > 1e-9 || fabs(w[2] - p[2]) > 1e-9) {
            fprintf(stderr, "  qs_ct_world is not the identity when unarmed\n");
            bad++;
        }
        if (bad) fprintf(stderr, "[selftest] solidify ct-frame FAIL\n");
        else fprintf(stderr, "[selftest] solidify ct-frame PASS\n");
        fails += bad;
    }
    if (ves_temp_dir(tmp, sizeof tmp) != 0) { fprintf(stderr, "  FAIL: no temp dir\n"); return 1; }
    snprintf(dir, sizeof dir, "%sqs_selftest_%d", tmp, ves_getpid());
    ves_mkdir(dir);
    snprintf(in_stem, sizeof in_stem, "%s/in", dir);
    snprintf(out_stem, sizeof out_stem, "%s/solid", dir);
    snprintf(refit_stem, sizeof refit_stem, "%s/refit", dir);
    snprintf(in_vmesh, sizeof in_vmesh, "%s.vmesh", in_stem);
    snprintf(out_vmesh, sizeof out_vmesh, "%s.vmesh", out_stem);

    /* ---- the synthetic scene --------------------------------------------- */
#define QS_TV(col, k, band, r, dphase, sup, prov, mat, lane) do {                \
        double ph_ = (double)(col) * du / R_ref + (dphase);                       \
        double z_ = z0 + (double)(k) * slice_h;                                   \
        if (qs_tv_push(&tv, (col), (k), (band), z_, (r) * cos(ph_), (r) * sin(ph_), \
                       ph_, (uint8_t)(sup), (uint8_t)(prov), (mat), (lane)) != 0) \
            fails++;                                                              \
    } while (0)
    /* run A: cols 0..39, R=100, interior hole rows 7..12 x cols 15..20;
     * band-0 removed at cols 30..31 x rows 2..3 (adopted band-1 there) and at
     * cols 34..35 x rows 2..3 (rejected band-1 there); a +8 vox jump at (10,30);
     * a 4x4 hole at rows 14..17 x cols 24..27 whose right rim (col 28) sits one
     * pitch out -- a chord to the next wrap, which the stretch gate refuses */
    for (long c = 0; c < 40; c++)
        for (int k = 0; k < np; k++) {
            int hole = k >= 7 && k <= 12 && c >= 15 && c <= 20;
            int adopt_cell = (c == 30 || c == 31) && (k == 2 || k == 3);
            int reject_cell = (c == 34 || c == 35) && (k == 2 || k == 3);
            int skip_hole = k >= 14 && k <= 17 && c >= 24 && c <= 27;
            double r = 100.0;
            if (hole || adopt_cell || reject_cell || skip_hole) continue;
            if (k == 10 && c == 30) r = 108.0;               /* the gate jump */
            if (k >= 14 && k <= 17 && c == 28) r = 100.0 + pitch;  /* foreign rim */
            QS_TV(c, k, 0, r, 0.0, (c == 3 && k == 3) ? 0 : 255, 0, 1, 0);
        }
    for (long c = 5; c <= 9; c++)                              /* duplicates */
        for (int k = 3; k <= 5; k++)
            QS_TV(c, k, 1, 100.0 + 1.5 * pitch, 0.0, 255, 0, 1, 0);
    for (long c = 30; c <= 31; c++)                            /* adoptable */
        for (int k = 2; k <= 3; k++)
            QS_TV(c, k, 1, 100.0 + 0.1 * pitch, 0.0, 255, 0, 1, 0);
    for (long c = 34; c <= 35; c++)                            /* other wrap */
        for (int k = 2; k <= 3; k++)
            QS_TV(c, k, 1, 100.0 + 1.5 * pitch, 0.0, 255, 0, 1, 0);
    /* run B: cols 42..61, R=109.5, materials 2 | 3 at col 52; a 2x2 hole
     * straddling the border (rows 10..11 x cols 51..52) has mixed rims */
    for (long c = 42; c < 62; c++)
        for (int k = 0; k < np; k++) {
            if ((k == 10 || k == 11) && (c == 51 || c == 52)) continue;
            QS_TV(c, k, 0, 109.5, 0.0, 255, 0, c < 52 ? 2 : 3, 1);
        }
    /* run C: cols 70..129, R=100: a notch open to the top rows 0..5 x cols
     * 75..80, a closed 8x6 hole rows 3..8 x cols 90..97, six SPARSE columns
     * 100..105 (row 0 only) between good ones, and twenty sparse columns
     * 110..129 (row 0 only) after the last good column */
    for (long c = 70; c < 130; c++)
        for (int k = 0; k < np; k++) {
            if (k <= 5 && c >= 75 && c <= 80) continue;
            if (k >= 3 && k <= 8 && c >= 90 && c <= 97) continue;
            if (k > 0 && ((c >= 100 && c <= 105) || c >= 110)) continue;
            QS_TV(c, k, 0, 100.0, 0.0, 255, 0, 7, 3);
        }
    /* run D: cols 140..159, R=50 (inside the core), hole rows 5..8 x cols 145..148.
     * u is arc length, so at R=50 the angle per column is du/50 (the macro's
     * R_ref phase gets the difference added) */
    for (long c = 140; c < 160; c++)
        for (int k = 0; k < np; k++) {
            if (k >= 5 && k <= 8 && c >= 145 && c <= 148) continue;
            QS_TV(c, k, 0, 50.0, (double)c * du / 50.0 - (double)c * du / R_ref, 255, 0, 9, 4);
        }
    /* run E: cols 170..199: good columns 170..179 at R=100, four sparse
     * columns 180..183 (row 0 only, R=100), good columns 184..199 at R+pitch
     * -- the next wrap fused into the same run: the gap is a wrap skip */
    for (long c = 170; c < 200; c++)
        for (int k = 0; k < np; k++) {
            if (k > 0 && c >= 180 && c <= 183) continue;
            QS_TV(c, k, 0, c >= 184 ? 100.0 + pitch : 100.0, 0.0, 255, 0, 11, 5);
        }
    /* run F: cols 210..214, three rows only -- confetti, never emitted */
    for (long c = 210; c < 215; c++)
        for (int k = 4; k < 7; k++)
            QS_TV(c, k, 0, 100.0, 0.0, 255, 0, 13, 6);
    /* run G: cols 230..259 at R=100: good columns 230..244 in order, five sparse
     * columns 245..249 (row 0), then good columns 250..259 whose angle runs
     * BACKWARDS (a fragment placed out of angular order): the gap fails the
     * order test */
    for (long c = 230; c < 260; c++)
        for (int k = 0; k < np; k++) {
            double dphase = c >= 250 ? (double)(244 - (c - 250) - c) * du / R_ref : 0.0;
            if (k > 0 && c >= 245 && c <= 249) continue;
            QS_TV(c, k, 0, 100.0, dphase, 255, 0, 15, 7);
        }
#undef QS_TV
    QS_CK(qs_tv_write(&tv, in_stem, du, dv, np, 1) == 0, "synthetic ribbon written");
    QS_CK(QuadribbonSolidify_complete(out_stem) == 0, "complete() is false before the run");

    arena = Arena_new();
    QS_CK(arena != NULL, "arena");
    if (arena == NULL) { qs_tv_free(&tv); return fails; }
    QS_CK(QuadribbonSolidify_run(arena, &cfg, in_vmesh, out_stem, &rep) == 0, "solidify run");
    QS_CK(QuadribbonSolidify_complete(out_stem) == 1, "complete() is true after the run");
    QS_CK(rep.n_runs == 7, "seven column runs");
    QS_CK(rep.dropped_order == 5 * 19, "a gap to a fragment placed out of angular order is refused");
    QS_CK(rep.dropped_confetti == 15, "a confetti run is counted and not emitted");
    QS_CK(rep.n_bands == 2 && rep.n_slices == np, "two peel bands of np slices");
    QS_CK(rep.duplicates == 15, "15 duplicate cells measured");
    QS_CK(rep.dup_dr_pitch_p50 > 1.4 && rep.dup_dr_pitch_p50 < 1.6, "duplicate |dr| ~1.5 pitch");
    QS_CK(rep.band1_adopted == 4, "4 same-wrap band-1 cells adopted");
    QS_CK(rep.band1_rejected == 4, "4 other-wrap band-1 cells rejected");
    QS_CK(rep.dropped_unbridged == 20 * 19, "sparse columns after the last good one stay empty");
    /* The foreign rim sits a pitch off its column's local radius, so the
     * COLUMN-CONTINUITY repair demotes it BEFORE the wrap-skip gate sees it
     * (2026-09-03).  What the wrap-skip gate still refuses is the gap between
     * two wraps; the rim's cells become fill on the correct wrap instead of a
     * hole, which is the better outcome and what the emission check below now
     * asserts. */
    QS_CK(rep.dropped_skip == 4 * 19, "the gap between two wraps is refused (wrap-skip)");
    QS_CK(rep.column_jumps > 0, "the column-continuity repair demotes the foreign rim");
    QS_CK(rep.dropped_chord == 0, "the adaptive chord gate defers the foreign rim to the wrap-skip gate");
    QS_CK(rep.dropped_material == 4, "a hole bridging two lineages is refused");
    QS_CK(rep.dropped_core == 16, "core-run hole is not filled");
    QS_CK(rep.filled == 40 + 36 + 48 + 6 * 19 + 20,
          "every row filled between good columns, plus the 20 cells the column repair demoted "
          "(A: 36 + 4, C: notch 36, hole 48, sparse gap 114)");
    QS_CK(rep.faces_skipped_gate > 0, "the +8 vox jump skips triangles");
    QS_CK(rep.faces_skipped_material > 0, "the material border skips triangles");
    QS_CK(rep.phase_mismatch_max < 1e-3, "sidecar phase agrees with atan2");
    QS_CK(rep.holes_refused >= 1, "the two-lineage hole is refused as a whole (the foreign-rim "
                                  "hole is now repaired by column continuity, not refused)");
    QS_CK(rep.cols_out > 0 && rep.cols_out < 260, "layout compacts the atlas gaps and the confetti run");

    /* read back through the public files */
    {
        MeshBinData m;
        uint8_t *prov = NULL;
        int32_t *mat = NULL, *lane = NULL, *run = NULL;
        int32_t *acol = NULL;
        Arena_Mark mark = Arena_save(arena);
        memset(&m, 0, sizeof m);
        QS_CK(MeshBin_read_arena(arena, out_vmesh, &m) == 0 && m.uv != NULL, "output readable");
        if (m.uv != NULL && m.nv > 0) {
            /* LAYOUT: no gap between occupied output columns wider than the separator */
            long cmin = LONG_MAX, cmax = LONG_MIN, worst = 0;
            uint8_t *occ = NULL;
            for (size_t i = 0; i < m.nv; i++) {
                long cc = (long)nearbyint((double)m.uv[i * 2] / du);
                if (cc < cmin) cmin = cc;
                if (cc > cmax) cmax = cc;
            }
            occ = (uint8_t *)ARENA_CALLOC(arena, (size_t)(cmax - cmin + 1), 1);
            for (size_t i = 0; i < m.nv; i++)
                occ[(size_t)((long)nearbyint((double)m.uv[i * 2] / du) - cmin)] = 1;
            for (long cc = cmin, last = cmin; cc <= cmax; cc++)
                if (occ[cc - cmin]) { if (cc - last - 1 > worst) worst = cc - last - 1; last = cc; }
            QS_CK(worst <= QS_RUN_GAP_COLS, "output columns are contiguous up to the run separator");
        }
        snprintf(path, sizeof path, "%s_provenance.u8", out_stem);
        prov = (uint8_t *)qs_read_exact(arena, path, m.nv, 1);
        snprintf(path, sizeof path, "%s_material_identity.i32", out_stem);
        mat = (int32_t *)qs_read_exact(arena, path, m.nv, 4);
        snprintf(path, sizeof path, "%s_lane.i32", out_stem);
        lane = (int32_t *)qs_read_exact(arena, path, m.nv, 4);
        snprintf(path, sizeof path, "%s_reconstruction_component.i32", out_stem);
        run = (int32_t *)qs_read_exact(arena, path, m.nv, 4);
        snprintf(path, sizeof path, "%s_atlas_column.i32", out_stem);
        acol = (int32_t *)qs_read_exact(arena, path, m.nv, 4);
        QS_CK(prov && mat && lane && run && acol, "sidecars sized nv");
        if (prov && mat && lane && run && acol) {
            /* LAYOUT: within a run the output column is the atlas column
             * shifted by one constant */
            size_t layout_bad = 0;
            long shift[8] = {LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN, LONG_MIN};
            for (size_t i = 0; i < m.nv; i++) {
                long lcol = lround((double)m.uv[i * 2] / du);
                long sh = (long)acol[i] - lcol;
                if (run[i] < 0 || run[i] >= 8) { layout_bad++; continue; }
                if (shift[run[i]] == LONG_MIN) shift[run[i]] = sh;
                else if (shift[run[i]] != sh) layout_bad++;
            }
            QS_CK(layout_bad == 0, "the layout shifts every run by one constant");
        }
        if (prov && mat && lane && run && acol) {
            double max_edge = 0.0, fill_r_err = 0.0, fill_phi_err = 0.0;
            size_t cross_run = 0, mixed_mat = 0, uv_mismatch = 0, z_err = 0;
            size_t fills_a = 0, fills_c_bad_label = 0, adopted = 0, obs_missing = 0;
            for (size_t f = 0; f < m.nf; f++) {
                const int32_t *t = m.faces + f * 3;
                for (int e = 0; e < 3; e++) {
                    const float *pa = m.verts + (size_t)t[e] * 3;
                    const float *pb = m.verts + (size_t)t[(e + 1) % 3] * 3;
                    double d = sqrt(pow((double)pa[0] - pb[0], 2) + pow((double)pa[1] - pb[1], 2) +
                                    pow((double)pa[2] - pb[2], 2));
                    if (d > max_edge) max_edge = d;
                }
                if (run[t[0]] != run[t[1]] || run[t[0]] != run[t[2]]) cross_run++;
                if (!qs_mat_ok(mat[t[0]], mat[t[1]], mat[t[2]])) mixed_mat++;
            }
            QS_CK(max_edge <= QS_WRAP_GATE_VOX + 1e-6, "no output edge over the wrap gate");
            QS_CK(cross_run == 0, "no face joins two runs");
            QS_CK(mixed_mat == 0, "no face joins two materials");
            for (size_t i = 0; i < m.nv; i++) {
                long col = (long)acol[i];   /* certificate space */
                int k = (int)lround((double)m.uv[i * 2 + 1] / dv);
                const float *p = m.verts + i * 3;
                double r = hypot((double)p[1], (double)p[2]);
                double phi = atan2((double)p[2], (double)p[1]);
                if (fabs((double)p[0] - (z0 + k * slice_h)) > 1e-4) z_err++;
                if (prov[i] == QS_PROV_FILL && col < 40) {
                    double want = (double)col * du / R_ref;
                    fills_a++;
                    if (fabs(r - 100.0) > fill_r_err) fill_r_err = fabs(r - 100.0);
                    if (fabs(phi - want) > fill_phi_err) fill_phi_err = fabs(phi - want);
                }
                if (prov[i] == QS_PROV_FILL && col >= 70 && col < 130 &&
                    (mat[i] != 7 || lane[i] != 3))
                    fills_c_bad_label++;
                if (prov[i] == QS_PROV_BAND1) adopted++;
                if (prov[i] == QS_PROV_OBSERVED_FIXED || prov[i] == QS_PROV_OBSERVED_GENERATED) {
                    /* the input band-0 vertex at this cell must exist with the
                     * same v bits; u is the layout column (exact lattice) */
                    int found = 0;
                    for (size_t j = 0; j < tv.n; j++) {
                        if (tv.band[j] == 0 && tv.col[j] == col && tv.k[j] == k) {
                            float v_in = (float)((double)tv.k[j] * dv);
                            float u_out = (float)((double)lround((double)m.uv[i * 2] / du) * du);
                            found = 1;
                            if (memcmp(&u_out, &m.uv[i * 2], 4) != 0 ||
                                memcmp(&v_in, &m.uv[i * 2 + 1], 4) != 0)
                                uv_mismatch++;
                            break;
                        }
                    }
                    if (!found) obs_missing++;
                }
            }
            QS_CK(fills_a == 60, "run A emits 60 fills (40 by the domain rule + 20 the column "
                                 "repair handed back)");
            fprintf(stderr, "  run A fills: max |r-100| = %.4f, max |phi err| = %.4f\n",
                    fill_r_err, fill_phi_err);
            QS_CK(fill_r_err < 0.05, "fills stay on the cylinder radius");
            QS_CK(fill_phi_err < 0.02, "fills follow the winding phase");
            QS_CK(z_err == 0, "every output vertex sits on its slice plane");
            QS_CK(fills_c_bad_label == 0, "fills inherit lane and material");
            QS_CK(adopted == 4, "adopted cells carry provenance 2");
            QS_CK(uv_mismatch == 0, "observed v is bit-identical to the input and u sits on the layout lattice");
            QS_CK(obs_missing == 0, "every observed output cell exists in the input");
            /* the notch, the closed hole and the sparse gap present; the
             * trailing sparse columns, the foreign-rim hole and the wrap-skip
             * gap absent */
            {
                int has_notch = 0, has_tail = 0, has_closed = 0, has_skip = 0;
                int has_gap = 0, has_wrapskip = 0;
                for (size_t i = 0; i < m.nv; i++) {
                    long col = (long)acol[i];   /* certificate space */
                    int k = (int)lround((double)m.uv[i * 2 + 1] / dv);
                    if (col == 77 && k == 2) has_notch = 1;
                    if (col == 120 && k == 10) has_tail = 1;
                    if (col == 93 && k == 5) has_closed = 1;
                    if (col == 25 && k == 15) has_skip = 1;
                    if (col == 102 && k == 10) has_gap = 1;
                    if (col == 181 && k == 10) has_wrapskip = 1;
                    if (col >= 210 && col <= 214) has_tail = 1;
                    if (col == 247 && k == 10) has_wrapskip = 1;
                }
                QS_CK(has_notch, "notch open to the border is filled (v-complete)");
                QS_CK(!has_tail, "sparse columns after the last good one not emitted");
                QS_CK(has_closed, "closed hole cell emitted");
                QS_CK(has_gap, "sparse gap between good columns filled on every row");
                /* the foreign-rim cell IS emitted now -- as FILL interpolated from its
                 * correct-wrap neighbours, not as the foreign observation */
                QS_CK(has_skip, "foreign-rim cell emitted as fill on the correct wrap");
                QS_CK(!has_wrapskip, "wrap-skip gap not emitted");
            }
        }
        Arena_restore(arena, mark);
    }
    /* refit idempotence: the output is a valid input; observations unchanged,
     * fills reproduced (band structure inferred, no stats sidecar) */
    {
        QsReport rep2;
        MeshBinData a, b;
        Arena_Mark mark = Arena_save(arena);
        snprintf(path, sizeof path, "%s_stats.json", out_stem);
        remove(path);
        QS_CK(QuadribbonSolidify_run(arena, &cfg, out_vmesh, refit_stem, &rep2) == 0, "refit run");
        QS_CK(rep2.n_bands == 1 && rep2.n_slices == np, "refit sees one band (inferred)");
        QS_CK(rep2.input_fills == rep.filled, "refit treats prior fills as unobserved");
        QS_CK(rep2.nv_out == rep.nv_out && rep2.nf_out == rep.nf_out, "refit reproduces the lattice");
        memset(&a, 0, sizeof a); memset(&b, 0, sizeof b);
        snprintf(path, sizeof path, "%s.vmesh", refit_stem);
        if (MeshBin_read_arena(arena, out_vmesh, &a) == 0 &&
            MeshBin_read_arena(arena, path, &b) == 0 && a.nv == b.nv) {
            double maxd = 0.0;
            for (size_t i = 0; i < a.nv * 3; i++) {
                double d = fabs((double)a.verts[i] - (double)b.verts[i]);
                if (d > maxd) maxd = d;
            }
            QS_CK(maxd < 1e-4, "refit geometry within 1e-4 of the first solid");
        } else {
            QS_CK(0, "refit output readable with the same nv");
        }
        Arena_restore(arena, mark);
    }
    /* Projectivity: A is a three-row absolute-V strip with a one-column hole.
     * B contains A, is deliberately written in reverse order, and adds rows
     * on both sides.  Those rows also offer a vertical bracket for A's holes;
     * the fixed horizontal-precedence rule must keep A byte-identical. */
    {
        QsTv pa, pb;
        QsConfig pcfg;
        QsReport pra, prb;
        Arena_T parena = Arena_new();
        char ai[QS_PATH], bi[QS_PATH], ao[QS_PATH], bo[QS_PATH];
        char aco[QS_PATH], bco[QS_PATH];
        int local_fails = 0;
        memset(&pa, 0, sizeof pa); memset(&pb, 0, sizeof pb);
        memset(&pcfg, 0, sizeof pcfg);
        memset(&pra, 0, sizeof pra); memset(&prb, 0, sizeof prb);
        pcfg.axis_y = 0.0; pcfg.axis_x = 0.0; pcfg.pitch = pitch;
        pcfg.projective_local = 1;
        pcfg.lattice_du = du; pcfg.lattice_dv = dv;
        pcfg.local_fill_radius = 4;
        snprintf(ai, sizeof ai, "%s/proj_a_in", dir);
        snprintf(bi, sizeof bi, "%s/proj_b_in", dir);
        snprintf(ao, sizeof ao, "%s/proj_a_out", dir);
        snprintf(bo, sizeof bo, "%s/proj_b_out", dir);
        snprintf(aco, sizeof aco, "%s/proj_a_copy", dir);
        snprintf(bco, sizeof bco, "%s/proj_b_copy", dir);
        for (int row = 10; row <= 12; row++)
            for (int col = 100; col <= 104; col++) {
                double ph = (double)col * du / R_ref;
                if (col == 102) continue;
                if (qs_tv_push(&pa, col, row, 0, (double)row * dv,
                               R_ref * cos(ph), R_ref * sin(ph), ph,
                               255, QS_PROV_OBSERVED_FIXED, 123456, 123456) != 0)
                    local_fails++;
            }
        for (int row = 13; row >= 8; row--)
            for (int col = 104; col >= 100; col--) {
                double ph = (double)col * du / R_ref;
                if (row >= 10 && row <= 12 && col == 102) continue;
                if (qs_tv_push(&pb, col, row, 0, (double)row * dv,
                               R_ref * cos(ph), R_ref * sin(ph), ph,
                               255, QS_PROV_OBSERVED_FIXED, 123456, 123456) != 0)
                    local_fails++;
            }
        /* A Stage-3 continuation at a locally bracketed address is upstream
         * output and must be preserved.  An isolated continuation is retained
         * too, but must not become an endpoint that fills its neighbour. */
        {
            const int cols[] = { 102, 110, 112 };
            const uint8_t support[] = { 0, 0, 255 };
            const uint8_t prov[] = {
                QS_PROV_OBSERVED_GENERATED,
                QS_PROV_OBSERVED_GENERATED,
                QS_PROV_OBSERVED_FIXED
            };
            for (size_t q = 0; q < sizeof cols / sizeof cols[0]; q++) {
                double ph = (double)cols[q] * du / R_ref;
                if (qs_tv_push(&pa, cols[q], 11, 0, 11.0 * dv,
                               R_ref * cos(ph), R_ref * sin(ph), ph,
                               support[q], prov[q], 123456, 123456) != 0)
                    local_fails++;
            }
            for (size_t q = sizeof cols / sizeof cols[0]; q-- > 0; ) {
                double ph = (double)cols[q] * du / R_ref;
                if (qs_tv_push(&pb, cols[q], 11, 0, 11.0 * dv,
                               R_ref * cos(ph), R_ref * sin(ph), ph,
                               support[q], prov[q], 123456, 123456) != 0)
                    local_fails++;
            }
        }
        if (local_fails == 0 &&
            (qs_tv_write(&pa, ai, du, dv, 32, 1) != 0 ||
             qs_tv_write(&pb, bi, du, dv, 32, 1) != 0))
            local_fails++;
        snprintf(path, sizeof path, "%s.vmesh", ai);
        if (local_fails == 0 &&
            QuadribbonSolidify_run(parena, &pcfg, path, ao, &pra) != 0)
            local_fails++;
        snprintf(path, sizeof path, "%s.vmesh", bi);
        if (local_fails == 0 &&
            QuadribbonSolidify_run(parena, &pcfg, path, bo, &prb) != 0)
            local_fails++;
        if (local_fails == 0 &&
            (!QuadribbonSolidify_complete(ao) ||
             !QuadribbonSolidify_complete(bo)))
            local_fails++;
        if (local_fails == 0)
            local_fails += qs_test_projective_subset(parena, ao, bo);
        if (local_fails == 0) {
            Arena_Mark pmark = Arena_save(parena);
            MeshBinData pm;
            const int32_t *addr = NULL, *dep = NULL;
            const uint8_t *kind = NULL;
            long at_bracket = -1, at_upstream = -1, at_recursive = -1;
            memset(&pm, 0, sizeof pm);
            snprintf(path, sizeof path, "%s.vmesh", ao);
            if (MeshBin_read_arena(parena, path, &pm) != 0) {
                local_fails++;
            } else {
                snprintf(path, sizeof path, "%s_lattice_address.i32", ao);
                addr = (const int32_t *)qs_read_exact(
                    parena, path, pm.nv * 2, sizeof(int32_t));
                snprintf(path, sizeof path, "%s_dependency_kind.u8", ao);
                kind = (const uint8_t *)qs_read_exact(
                    parena, path, pm.nv, sizeof(uint8_t));
                snprintf(path, sizeof path, "%s_dependency_endpoints.i32", ao);
                dep = (const int32_t *)qs_read_exact(
                    parena, path, pm.nv * 4, sizeof(int32_t));
                if (addr == NULL || kind == NULL || dep == NULL) {
                    local_fails++;
                } else {
                    at_bracket = qs_test_find_address(addr, pm.nv, 102, 11);
                    at_upstream = qs_test_find_address(addr, pm.nv, 110, 11);
                    at_recursive = qs_test_find_address(addr, pm.nv, 111, 11);
                    if (at_bracket < 0 ||
                        kind[at_bracket] != QS_DEP_UPSTREAM_GENERATED ||
                        at_upstream < 0 ||
                        kind[at_upstream] != QS_DEP_UPSTREAM_GENERATED ||
                        at_recursive >= 0 ||
                        dep[at_upstream * 4] != 110 ||
                        dep[at_upstream * 4 + 1] != 11 ||
                        dep[at_upstream * 4 + 2] != 110 ||
                        dep[at_upstream * 4 + 3] != 11)
                        local_fails++;
                }
            }
            Arena_restore(parena, pmark);
        }
        /* The production atlas contract is stronger than the legacy local
         * fill mode above: it copies each storage member exactly, so a sample
         * remains identical even if a larger domain assigns its claimant to a
         * different storage color. */
        pcfg.projective_copy_only = 1;
        snprintf(path, sizeof path, "%s.vmesh", ai);
        if (local_fails == 0 &&
            QuadribbonSolidify_run(parena, &pcfg, path, aco, &pra) != 0)
            local_fails++;
        snprintf(path, sizeof path, "%s.vmesh", bi);
        if (local_fails == 0 &&
            QuadribbonSolidify_run(parena, &pcfg, path, bco, &prb) != 0)
            local_fails++;
        if (local_fails == 0 &&
            (!QuadribbonSolidify_complete(aco) ||
             !QuadribbonSolidify_complete(bco) ||
             pra.nv_out != pa.n || prb.nv_out != pb.n ||
             pra.filled != 0 || prb.filled != 0 ||
             qs_test_projective_subset(parena, aco, bco) != 0))
            local_fails++;
        QS_CK(local_fails == 0,
              "projective A subset B is byte-identical; production copy-only atlas transport preserves every sample and legacy fills remain non-recursive");
        qs_tv_free(&pa); qs_tv_free(&pb);
        if (parena != NULL) Arena_dispose(&parena);
        {
            static const char *const ps[] = {
                ".vmesh", "_support.u8", "_provenance.u8", "_phase.f32",
                "_material_identity.i32", "_claimant_chart.i32",
                "_lane.i32", "_stats.json",
                "_reconstruction_component.i32", "_atlas_column.i32",
                "_lattice_address.i32", "_dependency_kind.u8",
                "_dependency_endpoints.i32", "_report.json"
            };
            const char *st[] = { ai, bi, ao, bo, aco, bco };
            for (size_t si = 0; si < sizeof st / sizeof st[0]; si++)
                for (size_t pi = 0; pi < sizeof ps / sizeof ps[0]; pi++) {
                    snprintf(path, sizeof path, "%s%s", st[si], ps[pi]);
                    remove(path);
                }
        }
    }
    /* completeness fails closed on a missing sidecar */
    snprintf(path, sizeof path, "%s_lane.i32", out_stem);
    remove(path);
    QS_CK(QuadribbonSolidify_complete(out_stem) == 0, "complete() false with a sidecar missing");

    /* cleanup (best effort) */
    {
        static const char *const stems[] = { "in", "solid", "refit" };
        static const char *const suf[] = {
            ".vmesh", "_support.u8", "_provenance.u8", "_phase.f32",
            "_material_identity.i32", "_claimant_chart.i32",
            "_reconstruction_component.i32", "_lane.i32",
            "_stats.json", "_provenance.png", "_report.json"
        };
        for (size_t s = 0; s < 3; s++)
            for (size_t i = 0; i < sizeof suf / sizeof suf[0]; i++) {
                snprintf(path, sizeof path, "%s/%s%s", dir, stems[s], suf[i]);
                remove(path);
            }
        snprintf(path, sizeof path, "%s/ribbon_phase.f32", dir); remove(path);
        snprintf(path, sizeof path, "%s/ribbon_material_identity.i32", dir); remove(path);
#ifdef _WIN32
        RemoveDirectoryA(dir);
#else
        rmdir(dir);
#endif
    }
    Arena_dispose(&arena);
    qs_tv_free(&tv);
    fprintf(stderr, "QuadribbonSolidify_selftest: %s (%d failure(s))\n",
            fails ? "FAIL" : "PASS", fails);
    return fails;
}
