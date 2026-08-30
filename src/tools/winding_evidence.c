/*
 * winding_evidence.c -- mine Villa winding constraints from tifxyz patches.
 *
 * Native geometry policy lives here: robust Archimedean phase estimates,
 * nearest-patch correspondence, confidence scoring, and accepted/candidate
 * decisions.  The Python package is only a launcher/artifact registrar.
 *
 * usage: winding_evidence <villa_dataset> <out_dir>
 *          --axis-point Z Y X [--axis-dir DZ DY DX] [--pitch P]
 *          [--winding-offset W] [--neighbor-distance D]
 *          [--max-phase-residual R] [--min-confidence C]
 *          [--max-relative-winding K] [--max-samples N] [--dry-run]
 *        winding_evidence --selftest
 *
 * out_dir receives winding_evidence.json plus abs_winding.json,
 * relative_windings.json, and same_windings.json in VC PCL v1 format.
 */
#include "../common/ves_platform.h"

#include <float.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
  #include <windows.h>
#else
  #include <dirent.h>
  #include <sys/stat.h>
#endif

#include "../common/arena.h"
#include "../common/kdtree.h"
#include "../common/tiff_io.h"

#define WE_PI 3.14159265358979323846264338327950288
#define WE_FORMAT "scrollfiesta-winding-evidence-v1"

static int regular_file_exists(const char *path)
{
#ifdef _WIN32
    DWORD attributes = GetFileAttributesA(path);
    return attributes != INVALID_FILE_ATTRIBUTES &&
           !(attributes & FILE_ATTRIBUTE_DIRECTORY);
#else
    struct stat status;
    return stat(path, &status) == 0 && S_ISREG(status.st_mode);
#endif
}

typedef struct {
    const char *dataset, *out_dir;
    double point_zyx[3], dir_zyx[3];
    double origin[3], axis[3], radial_x[3], radial_y[3]; /* xyz */
    double pitch, winding_offset, neighbor_distance;
    double max_phase_residual, min_confidence;
    int axis_point_set, max_relative_winding, dry_run;
    size_t max_samples;
} Params;

typedef struct {
    char id[256], rel[512], dir[1200];
    int verified, W, H;
    size_t valid_cells, raster_cells, n;
    float *xyz;                         /* n x 3, xyz */
    double *phase;                      /* n */
    double bbox_lo[3], bbox_hi[3], representative[3];
    double median, mad, p05, p95;
    int winding, absolute_accepted;
    double absolute_residual, absolute_confidence;
    Arena_T tree_arena;
    KDTree_T tree;
} Patch;

typedef struct {
    int a, b, delta, accepted;
    char id[560];
    double measured_delta, residual, delta_mad;
    double nearest_distance, bbox_gap, confidence;
    double point_a[3], point_b[3];
    size_t correspondences;
} Relation;

typedef struct {
    double distance, delta;
    size_t ia, ib;
} Correspondence;

static double clamp01(double value)
{
    if (value < 0.0) return 0.0;
    if (value > 1.0) return 1.0;
    return value;
}

static double dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void cross3(const double a[3], const double b[3], double out[3])
{
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

static int normalize3(double value[3])
{
    double length = sqrt(dot3(value, value));
    if (!isfinite(length) || length <= 1e-12) return -1;
    value[0] /= length; value[1] /= length; value[2] /= length;
    return 0;
}

static int prepare_axis(Params *P)
{
    double reference[3] = {1.0, 0.0, 0.0};
    /* Public z,y,x -> internal x,y,z. */
    for (int k = 0; k < 3; k++) {
        P->origin[k] = P->point_zyx[2-k];
        P->axis[k] = P->dir_zyx[2-k];
    }
    if (normalize3(P->axis) != 0) return -1;
    if (fabs(dot3(reference, P->axis)) > 0.95) {
        reference[0] = 0.0; reference[1] = 1.0;
    }
    double projection = dot3(reference, P->axis);
    for (int k = 0; k < 3; k++)
        P->radial_x[k] = reference[k] - projection * P->axis[k];
    if (normalize3(P->radial_x) != 0) return -1;
    cross3(P->axis, P->radial_x, P->radial_y);
    return normalize3(P->radial_y);
}

static double point_phase(const float point[3], const Params *P)
{
    double relative[3], radial[3];
    for (int k = 0; k < 3; k++) relative[k] = point[k] - P->origin[k];
    double axial = dot3(relative, P->axis);
    for (int k = 0; k < 3; k++) radial[k] = relative[k] - axial * P->axis[k];
    double rx = dot3(radial, P->radial_x);
    double ry = dot3(radial, P->radial_y);
    return hypot(rx, ry) / P->pitch - atan2(ry, rx) / (2.0 * WE_PI) +
           P->winding_offset;
}

static int compare_double(const void *aa, const void *bb)
{
    double a = *(const double *)aa, b = *(const double *)bb;
    return (a > b) - (a < b);
}

static int compare_patch(const void *aa, const void *bb)
{
    const Patch *a = (const Patch *)aa, *b = (const Patch *)bb;
    return strcmp(a->id, b->id);
}

static int compare_correspondence(const void *aa, const void *bb)
{
    const Correspondence *a = (const Correspondence *)aa;
    const Correspondence *b = (const Correspondence *)bb;
    if (a->distance < b->distance) return -1;
    if (a->distance > b->distance) return 1;
    return 0;
}

static double median_sorted(const double *values, size_t n)
{
    if (n & 1) return values[n/2];
    return 0.5 * (values[n/2-1] + values[n/2]);
}

static uint64_t fnv_update(uint64_t hash, const void *data, size_t n)
{
    const unsigned char *bytes = (const unsigned char *)data;
    for (size_t i = 0; i < n; i++) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

static uint64_t fnv_string(uint64_t hash, const char *text)
{
    return fnv_update(hash, text, strlen(text));
}

static void json_string(FILE *stream, const char *text)
{
    fputc('"', stream);
    for (const unsigned char *p = (const unsigned char *)text; *p; p++) {
        switch (*p) {
        case '"': fputs("\\\"", stream); break;
        case '\\': fputs("\\\\", stream); break;
        case '\b': fputs("\\b", stream); break;
        case '\f': fputs("\\f", stream); break;
        case '\n': fputs("\\n", stream); break;
        case '\r': fputs("\\r", stream); break;
        case '\t': fputs("\\t", stream); break;
        default:
            if (*p < 0x20) fprintf(stream, "\\u%04x", (unsigned)*p);
            else fputc(*p, stream);
        }
    }
    fputc('"', stream);
}

static int append_patch(Patch **patches, size_t *n, size_t *cap,
                        const char *root, const char *subdir,
                        const char *name, int verified)
{
    if (strlen(name) >= sizeof((*patches)[0].id)) {
        fprintf(stderr, "winding_evidence: patch id too long: %s\n", name);
        return -1;
    }
    if (*n == *cap) {
        size_t next = *cap ? *cap * 2 : 64;
        if (next < *cap || next > SIZE_MAX / sizeof(Patch)) return -1;
        Patch *grown = (Patch *)realloc(*patches, next * sizeof(*grown));
        if (grown == NULL) return -1;
        *patches = grown; *cap = next;
    }
    Patch *patch = &(*patches)[(*n)++];
    memset(patch, 0, sizeof(*patch));
    if (snprintf(patch->id, sizeof(patch->id), "%s", name) < 0 ||
        snprintf(patch->rel, sizeof(patch->rel), "%s/%s", subdir, name) < 0 ||
        snprintf(patch->dir, sizeof(patch->dir), "%s/%s/%s", root, subdir, name) < 0 ||
        strlen(patch->rel) >= sizeof(patch->rel)-1 ||
        strlen(patch->dir) >= sizeof(patch->dir)-1) {
        (*n)--;
        fprintf(stderr, "winding_evidence: patch path too long: %s\n", name);
        return -1;
    }
    patch->verified = verified;
    return 0;
}

static int scan_class(const char *root, const char *subdir, int verified,
                      Patch **patches, size_t *n, size_t *cap)
{
    char base[1200];
    snprintf(base, sizeof(base), "%s/%s", root, subdir);
#ifdef _WIN32
    char pattern[1300]; WIN32_FIND_DATAA fd;
    snprintf(pattern, sizeof(pattern), "%s/*", base);
    HANDLE handle = FindFirstFileA(pattern, &fd);
    if (handle == INVALID_HANDLE_VALUE) return 0;
    do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) ||
            fd.cFileName[0] == '.') continue;
        if (append_patch(patches, n, cap, root, subdir, fd.cFileName,
                         verified) != 0) {
            FindClose(handle); return -1;
        }
    } while (FindNextFileA(handle, &fd));
    FindClose(handle);
#else
    DIR *directory = opendir(base);
    if (directory == NULL) return 0;
    struct dirent *entry;
    while ((entry = readdir(directory)) != NULL) {
        if (entry->d_name[0] == '.') continue;
        char entry_path[1400]; struct stat status;
        int path_length = snprintf(entry_path, sizeof(entry_path), "%s/%s",
                                   base, entry->d_name);
        if (path_length < 0 || (size_t)path_length >= sizeof(entry_path) ||
            stat(entry_path, &status) != 0) {
            closedir(directory); return -1;
        }
        if (!S_ISDIR(status.st_mode)) continue;
        if (append_patch(patches, n, cap, root, subdir, entry->d_name,
                         verified) != 0) {
            closedir(directory); return -1;
        }
    }
    closedir(directory);
#endif
    return 0;
}

static int valid_coordinate(float x, float y, float z, uint8_t mask)
{
    return mask != 0 && isfinite(x) && isfinite(y) && isfinite(z) &&
           x >= 0.0f && y >= 0.0f && z >= 0.0f;
}

static int load_patch(Patch *patch, const Params *P)
{
    Arena_T arena = Arena_new();
    float *x = NULL, *y = NULL, *z = NULL;
    uint8_t *mask = NULL, *provenance = NULL;
    int wx = 0, hx = 0, wy = 0, hy = 0, wz = 0, hz = 0;
    int wm = 0, hm = 0, depth = 0;
    char path[1400];
    snprintf(path, sizeof(path), "%s/x.tif", patch->dir);
    if (TiffIO_load_float2d(arena, path, &x, &wx, &hx) != 0) goto invalid;
    snprintf(path, sizeof(path), "%s/y.tif", patch->dir);
    if (TiffIO_load_float2d(arena, path, &y, &wy, &hy) != 0) goto invalid;
    snprintf(path, sizeof(path), "%s/z.tif", patch->dir);
    if (TiffIO_load_float2d(arena, path, &z, &wz, &hz) != 0) goto invalid;
    snprintf(path, sizeof(path), "%s/mask.tif", patch->dir);
    if (TiffIO_load(arena, path, &mask, &depth, &hm, &wm) != 0) goto invalid;
    if (wx < 2 || hx < 2 || wy != wx || hy != hx || wz != wx || hz != hx ||
        wm != wx || hm != hx || depth != 1) goto invalid;
    snprintf(path, sizeof(path), "%s/provenance.tif", patch->dir);
    if (regular_file_exists(path)) {
        int provenance_depth = 0, provenance_h = 0, provenance_w = 0;
        if (TiffIO_load(arena, path, &provenance, &provenance_depth,
                        &provenance_h, &provenance_w) != 0 ||
            provenance_depth != 1 || provenance_w != wx || provenance_h != hx)
            goto invalid;
        for (size_t i = 0; i < (size_t)wx * (size_t)hx; i++) {
            if (provenance[i] > 3 ||
                ((provenance[i] == 0 || provenance[i] == 3) && mask[i] != 0) ||
                ((provenance[i] == 1 || provenance[i] == 2) && mask[i] == 0))
                goto invalid;
        }
    }
    patch->W = wx; patch->H = hx;
    patch->raster_cells = (size_t)wx * (size_t)hx;
    for (int k = 0; k < 3; k++) {
        patch->bbox_lo[k] = DBL_MAX;
        patch->bbox_hi[k] = -DBL_MAX;
    }
    size_t valid = 0, valid_quads = 0;
    for (size_t i = 0; i < patch->raster_cells; i++)
        if (valid_coordinate(x[i], y[i], z[i], mask[i])) valid++;
    for (int row = 0; row + 1 < hx; row++) {
        for (int col = 0; col + 1 < wx; col++) {
            size_t i0 = (size_t)row*(size_t)wx+(size_t)col;
            size_t i1 = i0+1, i2 = i0+(size_t)wx, i3 = i2+1;
            if (valid_coordinate(x[i0],y[i0],z[i0],mask[i0]) &&
                valid_coordinate(x[i1],y[i1],z[i1],mask[i1]) &&
                valid_coordinate(x[i2],y[i2],z[i2],mask[i2]) &&
                valid_coordinate(x[i3],y[i3],z[i3],mask[i3])) valid_quads++;
        }
    }
    if (valid == 0 || valid_quads == 0) goto invalid;
    patch->valid_cells = valid;
    size_t stride = (valid + P->max_samples - 1) / P->max_samples;
    size_t allocation = (valid + stride - 1) / stride;
    patch->xyz = (float *)malloc(allocation * 3 * sizeof(float));
    patch->phase = (double *)malloc(allocation * sizeof(double));
    if (patch->xyz == NULL || patch->phase == NULL) goto oom;
    size_t seen = 0;
    for (size_t i = 0; i < patch->raster_cells; i++) {
        if (!valid_coordinate(x[i], y[i], z[i], mask[i])) continue;
        const double coordinate[3] = {x[i], y[i], z[i]};
        for (int k = 0; k < 3; k++) {
            if (coordinate[k] < patch->bbox_lo[k]) patch->bbox_lo[k] = coordinate[k];
            if (coordinate[k] > patch->bbox_hi[k]) patch->bbox_hi[k] = coordinate[k];
        }
        if (seen % stride == 0 && patch->n < allocation) {
            patch->xyz[patch->n*3+0] = x[i];
            patch->xyz[patch->n*3+1] = y[i];
            patch->xyz[patch->n*3+2] = z[i];
            patch->phase[patch->n] = point_phase(&patch->xyz[patch->n*3], P);
            if (!isfinite(patch->phase[patch->n])) goto invalid;
            patch->n++;
        }
        seen++;
    }
    if (patch->n == 0) goto invalid;

    double *sorted = (double *)malloc(patch->n * sizeof(double));
    double *deviation = (double *)malloc(patch->n * sizeof(double));
    if (sorted == NULL || deviation == NULL) {
        free(sorted); free(deviation); goto oom;
    }
    memcpy(sorted, patch->phase, patch->n * sizeof(double));
    qsort(sorted, patch->n, sizeof(double), compare_double);
    patch->median = median_sorted(sorted, patch->n);
    size_t i05 = (size_t)floor(0.05 * (double)(patch->n - 1));
    size_t i95 = (size_t)floor(0.95 * (double)(patch->n - 1));
    patch->p05 = sorted[i05]; patch->p95 = sorted[i95];
    size_t representative = 0;
    double best = DBL_MAX;
    for (size_t i = 0; i < patch->n; i++) {
        deviation[i] = fabs(patch->phase[i] - patch->median);
        if (deviation[i] < best) { best = deviation[i]; representative = i; }
    }
    qsort(deviation, patch->n, sizeof(double), compare_double);
    patch->mad = median_sorted(deviation, patch->n);
    free(sorted); free(deviation);
    for (int k = 0; k < 3; k++)
        patch->representative[k] = patch->xyz[representative*3+(size_t)k];
    if (patch->median < (double)INT_MIN + 1.0 ||
        patch->median > (double)INT_MAX - 1.0) goto invalid;
    patch->winding = (int)floor(patch->median + 0.5);
    patch->absolute_residual = fabs(patch->median - (double)patch->winding);
    double trust = patch->verified ? 1.0 : 0.72;
    double coherence = exp(-patch->mad / 0.10);
    double phase_score = fmax(0.0, 1.0 - patch->absolute_residual / 0.25);
    patch->absolute_confidence = clamp01(trust * coherence * phase_score);
    /* Unverified scroll hints remain useful evidence, but cannot silently turn
     * themselves into solver constraints. Classification is a native gate. */
    patch->absolute_accepted = patch->verified && patch->winding > 0 &&
                               patch->absolute_confidence >= P->min_confidence;
    patch->tree_arena = Arena_new();
    patch->tree = KDTree_new(patch->tree_arena, patch->xyz, patch->n);
    Arena_dispose(&arena);
    return 0;

invalid:
    fprintf(stderr, "winding_evidence: invalid tifxyz patch %s\n", patch->dir);
    Arena_dispose(&arena);
    return -1;
oom:
    fprintf(stderr, "winding_evidence: out of memory loading %s\n", patch->id);
    Arena_dispose(&arena);
    return -1;
}

static double bbox_gap(const Patch *a, const Patch *b)
{
    double sum = 0.0;
    for (int k = 0; k < 3; k++) {
        double gap = 0.0;
        if (a->bbox_lo[k] > b->bbox_hi[k]) gap = a->bbox_lo[k] - b->bbox_hi[k];
        else if (b->bbox_lo[k] > a->bbox_hi[k]) gap = b->bbox_lo[k] - a->bbox_hi[k];
        sum += gap * gap;
    }
    return sqrt(sum);
}

static void keep_correspondence(Correspondence best[128], size_t *n,
                                double distance, double delta,
                                size_t ia, size_t ib)
{
    if (*n < 128) {
        best[*n].distance = distance; best[*n].delta = delta;
        best[*n].ia = ia; best[*n].ib = ib; (*n)++;
        return;
    }
    size_t worst = 0;
    for (size_t i = 1; i < *n; i++)
        if (best[i].distance > best[worst].distance) worst = i;
    if (distance < best[worst].distance) {
        best[worst].distance = distance; best[worst].delta = delta;
        best[worst].ia = ia; best[worst].ib = ib;
    }
}

static int measure_pair(const Patch *a, const Patch *b, const Params *P,
                        Relation *relation)
{
    Correspondence best[128]; size_t nbest = 0;
    for (size_t i = 0; i < a->n; i++) {
        float distance2 = 0.0f;
        size_t j = KDTree_nearest(b->tree, &a->xyz[i*3], &distance2);
        keep_correspondence(best, &nbest, sqrt((double)distance2),
                            b->phase[j] - a->phase[i], i, j);
    }
    for (size_t j = 0; j < b->n; j++) {
        float distance2 = 0.0f;
        size_t i = KDTree_nearest(a->tree, &b->xyz[j*3], &distance2);
        keep_correspondence(best, &nbest, sqrt((double)distance2),
                            b->phase[j] - a->phase[i], i, j);
    }
    if (nbest == 0) return 0;
    qsort(best, nbest, sizeof(best[0]), compare_correspondence);
    if (best[0].distance > P->neighbor_distance) return 0;
    size_t keep = nbest < 64 ? nbest : 64;
    while (keep > 1 && best[keep-1].distance > P->neighbor_distance) keep--;
    double delta[64], deviation[64];
    for (size_t i = 0; i < keep; i++) delta[i] = best[i].delta;
    qsort(delta, keep, sizeof(double), compare_double);
    relation->measured_delta = median_sorted(delta, keep);
    for (size_t i = 0; i < keep; i++)
        deviation[i] = fabs(delta[i] - relation->measured_delta);
    qsort(deviation, keep, sizeof(double), compare_double);
    relation->delta_mad = median_sorted(deviation, keep);
    relation->nearest_distance = best[0].distance;
    relation->correspondences = keep;
    for (int k = 0; k < 3; k++) {
        relation->point_a[k] = a->xyz[best[0].ia*3+(size_t)k];
        relation->point_b[k] = b->xyz[best[0].ib*3+(size_t)k];
    }
    return 1;
}

static void relation_id(char *out, size_t capacity, const Patch *a, const Patch *b)
{
    uint64_t hash = UINT64_C(1469598103934665603);
    hash = fnv_string(hash, a->id); hash = fnv_update(hash, "\0", 1);
    hash = fnv_string(hash, b->id);
    snprintf(out, capacity, "%s__%s__%016llx", a->id, b->id,
             (unsigned long long)hash);
}

static int append_relation(Relation **relations, size_t *n, size_t *cap,
                           const Patch *patches, int ia, int ib,
                           const Params *P)
{
    if (*n == *cap) {
        size_t next = *cap ? *cap * 2 : 128;
        if (next < *cap || next > SIZE_MAX / sizeof(Relation)) return -1;
        Relation *grown = (Relation *)realloc(*relations, next * sizeof(*grown));
        if (grown == NULL) return -1;
        *relations = grown; *cap = next;
    }
    Relation relation; memset(&relation, 0, sizeof(relation));
    relation.a = ia; relation.b = ib;
    relation.bbox_gap = bbox_gap(&patches[ia], &patches[ib]);
    if (relation.bbox_gap > P->neighbor_distance ||
        !measure_pair(&patches[ia], &patches[ib], P, &relation)) return 0;
    if (!isfinite(relation.measured_delta) ||
        relation.measured_delta < (double)INT_MIN + 1.0 ||
        relation.measured_delta > (double)INT_MAX - 1.0) return 0;
    relation.delta = (int)floor(relation.measured_delta + 0.5);
    if (abs(relation.delta) > P->max_relative_winding) return 0;
    relation.residual = fabs(relation.measured_delta - (double)relation.delta);
    double coherence = exp(-fmax(fmax(patches[ia].mad, patches[ib].mad),
                                 relation.delta_mad) / 0.12);
    double distance_score = exp(-relation.nearest_distance / (1.5 * P->pitch));
    double phase_score = fmax(0.0, 1.0 - relation.residual / P->max_phase_residual);
    double trust = patches[ia].verified && patches[ib].verified ? 1.0 : 0.78;
    relation.confidence = clamp01(trust * phase_score *
                                  sqrt(coherence * distance_score));
    relation.accepted = patches[ia].verified && patches[ib].verified &&
                         relation.residual <= P->max_phase_residual &&
                         relation.confidence >= P->min_confidence;
    relation_id(relation.id, sizeof(relation.id), &patches[ia], &patches[ib]);
    (*relations)[(*n)++] = relation;
    return 0;
}

static int promote_file(const char *temporary, const char *destination)
{
#ifdef _WIN32
    return MoveFileExA(temporary, destination,
                       MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) ? 0 : -1;
#else
    return rename(temporary, destination);
#endif
}

static FILE *open_output(const Params *P, const char *name,
                         char temporary[1400], char destination[1400])
{
    int destination_length = snprintf(destination, 1400, "%s/%s", P->out_dir, name);
    int temporary_length = snprintf(temporary, 1400, "%s.tmp", destination);
    if (destination_length < 0 || destination_length >= 1400 ||
        temporary_length < 0 || temporary_length >= 1400) return NULL;
    if (ves_ensure_parent_dir(destination) != 0) return NULL;
    return fopen(temporary, "wb");
}

static uint64_t evidence_hash(const Patch *patches, size_t npatch,
                              const Relation *relations, size_t nrelation,
                              const Params *P)
{
    uint64_t hash = UINT64_C(1469598103934665603);
    hash = fnv_update(hash, P->point_zyx, sizeof(P->point_zyx));
    hash = fnv_update(hash, P->dir_zyx, sizeof(P->dir_zyx));
    hash = fnv_update(hash, &P->pitch, sizeof(P->pitch));
    hash = fnv_update(hash, &P->winding_offset, sizeof(P->winding_offset));
    hash = fnv_update(hash, &P->neighbor_distance, sizeof(P->neighbor_distance));
    hash = fnv_update(hash, &P->max_phase_residual, sizeof(P->max_phase_residual));
    hash = fnv_update(hash, &P->min_confidence, sizeof(P->min_confidence));
    hash = fnv_update(hash, &P->max_relative_winding,
                      sizeof(P->max_relative_winding));
    hash = fnv_update(hash, &P->max_samples, sizeof(P->max_samples));
    for (size_t i = 0; i < npatch; i++) {
        hash = fnv_string(hash, patches[i].id);
        hash = fnv_update(hash, &patches[i].verified, sizeof(int));
        hash = fnv_update(hash, &patches[i].valid_cells, sizeof(size_t));
        hash = fnv_update(hash, &patches[i].median, sizeof(double));
        hash = fnv_update(hash, &patches[i].mad, sizeof(double));
        hash = fnv_update(hash, &patches[i].winding, sizeof(int));
        hash = fnv_update(hash, &patches[i].absolute_confidence, sizeof(double));
        hash = fnv_update(hash, &patches[i].absolute_accepted, sizeof(int));
    }
    for (size_t i = 0; i < nrelation; i++) {
        hash = fnv_string(hash, relations[i].id);
        hash = fnv_update(hash, &relations[i].delta, sizeof(int));
        hash = fnv_update(hash, &relations[i].confidence, sizeof(double));
        hash = fnv_update(hash, &relations[i].accepted, sizeof(int));
    }
    return hash;
}

static void write_point(FILE *stream, const double point[3], double winding)
{
    fprintf(stream,
            "{\"p\":[%.9g,%.9g,%.9g],\"creation_time\":0,\"wind_a\":%.9g}",
            point[0], point[1], point[2], winding);
}

static int write_evidence(const Patch *patches, size_t npatch,
                          const Relation *relations, size_t nrelation,
                          const Params *P, uint64_t generation,
                          size_t abs_ok, size_t same_ok, size_t relative_ok)
{
    char temporary[1400], destination[1400];
    FILE *stream = open_output(P, "winding_evidence.json", temporary, destination);
    if (stream == NULL) return -1;
    fprintf(stream, "{\n  \"format\": \"%s\", \"version\": 1,\n", WE_FORMAT);
    fprintf(stream, "  \"coordinate_order\": \"xyz_level0_voxels\",\n");
    fprintf(stream, "  \"generation_id\": \"%016llx\",\n",
            (unsigned long long)generation);
    fprintf(stream, "  \"source\": {\"dataset\": ");
    json_string(stream, P->dataset);
    fprintf(stream, "},\n  \"axis\": {\"point_zyx\":[%.9g,%.9g,%.9g],"
                    "\"direction_zyx\":[%.9g,%.9g,%.9g],\"pitch\":%.9g,"
                    "\"winding_offset\":%.9g,\"angular_zero_xyz\":[%.9g,%.9g,%.9g],"
                    "\"angular_positive_xyz\":[%.9g,%.9g,%.9g],"
                    "\"phase_formula\":\"radius/pitch - theta/(2*pi) + winding_offset\"},\n",
            P->point_zyx[0],P->point_zyx[1],P->point_zyx[2],
            P->dir_zyx[0],P->dir_zyx[1],P->dir_zyx[2],P->pitch,P->winding_offset,
            P->radial_x[0],P->radial_x[1],P->radial_x[2],
            P->radial_y[0],P->radial_y[1],P->radial_y[2]);
    fprintf(stream, "  \"parameters\": {\"max_samples_per_patch\":%zu,"
                    "\"neighbor_distance\":%.9g,\"max_phase_residual\":%.9g,"
                    "\"min_confidence\":%.9g,\"max_relative_winding\":%d},\n",
            P->max_samples, P->neighbor_distance, P->max_phase_residual,
            P->min_confidence, P->max_relative_winding);
    fprintf(stream, "  \"nodes\": {\n");
    for (size_t i = 0; i < npatch; i++) {
        const Patch *patch = &patches[i];
        fprintf(stream, "    "); json_string(stream, patch->id); fprintf(stream, ": {");
        fprintf(stream, "\"id\":"); json_string(stream, patch->id);
        fprintf(stream, ",\"classification\":\"%s\",\"path\":",
                patch->verified ? "verified" : "unverified");
        json_string(stream, patch->rel);
        fprintf(stream,
                ",\"bbox_xyz\":[[%.9g,%.9g,%.9g],[%.9g,%.9g,%.9g]],"
                "\"representative_xyz\":[%.9g,%.9g,%.9g],"
                "\"phase\":{\"median\":%.9g,\"mad\":%.9g,\"p05\":%.9g,\"p95\":%.9g},"
                "\"sampling\":{\"samples\":%zu,\"valid_cells\":%zu,\"raster_cells\":%zu},"
                "\"absolute\":{\"winding\":%d,\"residual\":%.9g,\"confidence\":%.9g,"
                "\"status\":\"%s\",\"provenance\":{\"method\":"
                "\"archimedean_phase_robust_median\",\"source_patch\":",
                patch->bbox_lo[0],patch->bbox_lo[1],patch->bbox_lo[2],
                patch->bbox_hi[0],patch->bbox_hi[1],patch->bbox_hi[2],
                patch->representative[0],patch->representative[1],patch->representative[2],
                patch->median,patch->mad,patch->p05,patch->p95,
                patch->n,patch->valid_cells,patch->raster_cells,
                patch->winding,patch->absolute_residual,patch->absolute_confidence,
                patch->absolute_accepted ? "accepted" : "candidate");
        json_string(stream, patch->id); fprintf(stream, "}}}%s\n",
                i + 1 < npatch ? "," : "");
    }
    fprintf(stream, "  },\n  \"relations\": {\n");
    for (size_t i = 0; i < nrelation; i++) {
        const Relation *relation = &relations[i];
        const Patch *a = &patches[relation->a], *b = &patches[relation->b];
        fprintf(stream, "    "); json_string(stream, relation->id); fprintf(stream, ": {");
        fprintf(stream, "\"id\":"); json_string(stream, relation->id);
        fprintf(stream, ",\"a\":"); json_string(stream, a->id);
        fprintf(stream, ",\"b\":"); json_string(stream, b->id);
        fprintf(stream,
                ",\"type\":\"%s\",\"delta_winding\":%d,\"measured_delta\":%.9g,"
                "\"residual\":%.9g,\"delta_mad\":%.9g,\"nearest_distance\":%.9g,"
                "\"bbox_gap\":%.9g,\"confidence\":%.9g,\"status\":\"%s\","
                "\"point_a_xyz\":[%.9g,%.9g,%.9g],\"point_b_xyz\":[%.9g,%.9g,%.9g],"
                "\"provenance\":{\"method\":\"nearest_tifxyz_phase_correspondence\","
                "\"correspondences\":%zu,\"source_patches\":[",
                relation->delta == 0 ? "same" : "relative", relation->delta,
                relation->measured_delta, relation->residual, relation->delta_mad,
                relation->nearest_distance, relation->bbox_gap, relation->confidence,
                relation->accepted ? "accepted" : "candidate",
                relation->point_a[0],relation->point_a[1],relation->point_a[2],
                relation->point_b[0],relation->point_b[1],relation->point_b[2],
                relation->correspondences);
        json_string(stream, a->id); fputc(',', stream); json_string(stream, b->id);
        fprintf(stream, "]}}%s\n", i + 1 < nrelation ? "," : "");
    }
    fprintf(stream,
            "  },\n  \"summary\":{\"patches\":%zu,\"absolute_accepted\":%zu,"
            "\"same_accepted\":%zu,\"relative_accepted\":%zu,"
            "\"relation_candidates\":%zu}\n}\n",
            npatch, abs_ok, same_ok, relative_ok,
            nrelation - same_ok - relative_ok);
    if (fclose(stream) != 0 || promote_file(temporary, destination) != 0) {
        ves_unlink(temporary); return -1;
    }
    return 0;
}

enum PclKind { PCL_ABSOLUTE, PCL_SAME, PCL_RELATIVE };

static int write_pcl(const Patch *patches, size_t npatch,
                     const Relation *relations, size_t nrelation,
                     const Params *P, uint64_t generation,
                     enum PclKind kind, const char *filename)
{
    char temporary[1400], destination[1400];
    FILE *stream = open_output(P, filename, temporary, destination);
    if (stream == NULL) return -1;
    fprintf(stream, "{\n  \"vc_pointcollections_json_version\":\"1\",\n  \"collections\":{\n");
    size_t total = 0;
    if (kind == PCL_ABSOLUTE) {
        for (size_t i = 0; i < npatch; i++) if (patches[i].absolute_accepted) total++;
    } else {
        for (size_t i = 0; i < nrelation; i++)
            if (relations[i].accepted &&
                ((kind == PCL_SAME) == (relations[i].delta == 0))) total++;
    }
    size_t emitted = 0;
    if (kind == PCL_ABSOLUTE) {
        for (size_t i = 0; i < npatch; i++) {
            const Patch *patch = &patches[i];
            if (!patch->absolute_accepted) continue;
            fprintf(stream, "    \"%zu\":{\"name\":", emitted);
            char name[600];
            snprintf(name, sizeof(name), "scrollfiesta__absolute__%s", patch->id);
            json_string(stream, name);
            fprintf(stream, ",\"points\":{\"0\":");
            write_point(stream, patch->representative, (double)patch->winding);
            fprintf(stream,
                    "},\"metadata\":{\"winding_is_absolute\":true},"
                    "\"color\":[0.95,0.42,0.12],\"tags\":{"
                    "\"scrollfiesta_managed\":\"1\",\"scrollfiesta_kind\":\"absolute\","
                    "\"scrollfiesta_patch\":");
            json_string(stream, patch->id);
            fprintf(stream, ",\"scrollfiesta_confidence\":\"%.6f\","
                    "\"scrollfiesta_generation\":\"%016llx\","
                    "\"scrollfiesta_provenance\":\"archimedean_phase_robust_median\"}}%s\n",
                    patch->absolute_confidence, (unsigned long long)generation,
                    ++emitted < total ? "," : "");
        }
    } else {
        for (size_t i = 0; i < nrelation; i++) {
            const Relation *relation = &relations[i];
            if (!relation->accepted ||
                ((kind == PCL_SAME) != (relation->delta == 0))) continue;
            const Patch *a = &patches[relation->a], *b = &patches[relation->b];
            fprintf(stream, "    \"%zu\":{\"name\":", emitted);
            char name[900];
            snprintf(name, sizeof(name), "scrollfiesta__%s__%s__%s",
                     relation->delta == 0 ? "same" : "relative", a->id, b->id);
            json_string(stream, name);
            fprintf(stream, ",\"points\":{\"0\":");
            write_point(stream, relation->point_a, 0.0);
            fprintf(stream, ",\"1\":");
            write_point(stream, relation->point_b, (double)relation->delta);
            fprintf(stream,
                    "},\"metadata\":{\"winding_is_absolute\":false},"
                    "\"color\":[%.2f,%.2f,%.2f],\"tags\":{"
                    "\"scrollfiesta_managed\":\"1\",\"scrollfiesta_kind\":\"%s\","
                    "\"scrollfiesta_relation\":",
                    relation->delta == 0 ? 0.15 : 0.18,
                    relation->delta == 0 ? 0.78 : 0.46,
                    relation->delta == 0 ? 0.34 : 0.95,
                    relation->delta == 0 ? "same" : "relative");
            json_string(stream, relation->id);
            fprintf(stream, ",\"scrollfiesta_patches\":");
            char pair[560]; snprintf(pair, sizeof(pair), "%s|%s", a->id, b->id);
            json_string(stream, pair);
            fprintf(stream, ",\"scrollfiesta_delta_winding\":\"%d\","
                    "\"scrollfiesta_confidence\":\"%.6f\","
                    "\"scrollfiesta_generation\":\"%016llx\","
                    "\"scrollfiesta_provenance\":\"nearest_tifxyz_phase_correspondence\"}}%s\n",
                    relation->delta, relation->confidence,
                    (unsigned long long)generation, ++emitted < total ? "," : "");
        }
    }
    fprintf(stream, "  }\n}\n");
    if (fclose(stream) != 0 || promote_file(temporary, destination) != 0) {
        ves_unlink(temporary); return -1;
    }
    return 0;
}

static void free_patches(Patch *patches, size_t n)
{
    if (patches == NULL) return;
    for (size_t i = 0; i < n; i++) {
        free(patches[i].xyz); free(patches[i].phase);
        if (patches[i].tree_arena != NULL)
            Arena_dispose(&patches[i].tree_arena);
    }
    free(patches);
}

static void synth_patch(Patch *patch, const Params *P, const char *id,
                        double winding, double z0)
{
    memset(patch, 0, sizeof(*patch));
    snprintf(patch->id, sizeof(patch->id), "%s", id);
    patch->verified = 1; patch->n = 16;
    patch->xyz = (float *)malloc(patch->n * 3 * sizeof(float));
    patch->phase = (double *)malloc(patch->n * sizeof(double));
    for (size_t i = 0; i < patch->n; i++) {
        double theta = 0.20 + 0.015*(double)i;
        double radius = P->pitch * (winding + theta/(2.0*WE_PI));
        patch->xyz[i*3+0] = (float)(100.0 + radius*cos(theta));
        patch->xyz[i*3+1] = (float)(100.0 + radius*sin(theta));
        patch->xyz[i*3+2] = (float)(z0 + 0.25*(double)i);
        patch->phase[i] = point_phase(&patch->xyz[i*3], P);
    }
    patch->tree_arena = Arena_new();
    patch->tree = KDTree_new(patch->tree_arena, patch->xyz, patch->n);
}

static int selftest(void)
{
    Params P; memset(&P, 0, sizeof(P));
    P.point_zyx[0] = 0.0; P.point_zyx[1] = 100.0; P.point_zyx[2] = 100.0;
    P.dir_zyx[0] = 1.0; P.pitch = 10.0; P.neighbor_distance = 16.0;
    P.max_phase_residual = 0.22; P.min_confidence = 0.30;
    P.max_relative_winding = 12;
    int fails = 0;
    if (prepare_axis(&P) != 0) return 1;
    float point[3];
    double theta = 0.73, winding = 5.0;
    double radius = P.pitch * (winding + theta/(2.0*WE_PI));
    point[0] = (float)(100.0 + radius*cos(theta));
    point[1] = (float)(100.0 + radius*sin(theta)); point[2] = 7.0f;
    if (fabs(point_phase(point, &P) - winding) > 2e-6) {
        fprintf(stderr, "selftest: axis phase mismatch\n"); fails++;
    }
    Patch patches[2];
    synth_patch(&patches[0], &P, "a", 3.0, 10.0);
    synth_patch(&patches[1], &P, "b", 4.0, 10.0);
    patches[0].mad = patches[1].mad = 0.0;
    Relation *relations = NULL; size_t nr = 0, cap = 0;
    if (append_relation(&relations, &nr, &cap, patches, 0, 1, &P) != 0 ||
        nr != 1 || relations[0].delta != 1 ||
        fabs(relations[0].measured_delta - 1.0) > 2e-5 ||
        !relations[0].accepted) {
        fprintf(stderr, "selftest: relative relation mismatch\n"); fails++;
    }
    for (int i = 0; i < 2; i++) {
        free(patches[i].xyz); free(patches[i].phase);
        Arena_dispose(&patches[i].tree_arena);
    }
    free(relations);
    fprintf(stderr, "[winding_evidence selftest] %s (%d failures)\n",
            fails ? "FAILED" : "PASSED", fails);
    return fails ? 1 : 0;
}

static void usage(const char *program)
{
    fprintf(stderr,
        "usage: %s <villa_dataset> <out_dir> --axis-point Z Y X\n"
        "         [--axis-dir DZ DY DX] [--pitch P] [--winding-offset W]\n"
        "         [--neighbor-distance D] [--max-phase-residual R]\n"
        "         [--min-confidence C] [--max-relative-winding K]\n"
        "         [--max-samples N] [--dry-run]\n"
        "       %s --selftest\n", program, program);
}

int main(int argc, char **argv)
{
    setlocale(LC_NUMERIC, "C");
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc < 4) { usage(argv[0]); return 2; }
    Params P; memset(&P, 0, sizeof(P));
    P.dataset = argv[1]; P.out_dir = argv[2];
    P.dir_zyx[0] = 1.0; P.pitch = 9.5;
    P.max_phase_residual = 0.22; P.min_confidence = 0.55;
    P.max_relative_winding = 12; P.max_samples = 4096;
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--axis-point") == 0 && i+3 < argc) {
            for (int k = 0; k < 3; k++) P.point_zyx[k] = atof(argv[++i]);
            P.axis_point_set = 1;
        } else if (strcmp(argv[i], "--axis-dir") == 0 && i+3 < argc) {
            for (int k = 0; k < 3; k++) P.dir_zyx[k] = atof(argv[++i]);
        } else if (strcmp(argv[i], "--pitch") == 0 && i+1 < argc) {
            P.pitch = atof(argv[++i]);
        } else if (strcmp(argv[i], "--winding-offset") == 0 && i+1 < argc) {
            P.winding_offset = atof(argv[++i]);
        } else if (strcmp(argv[i], "--neighbor-distance") == 0 && i+1 < argc) {
            P.neighbor_distance = atof(argv[++i]);
        } else if (strcmp(argv[i], "--max-phase-residual") == 0 && i+1 < argc) {
            P.max_phase_residual = atof(argv[++i]);
        } else if (strcmp(argv[i], "--min-confidence") == 0 && i+1 < argc) {
            P.min_confidence = atof(argv[++i]);
        } else if (strcmp(argv[i], "--max-relative-winding") == 0 && i+1 < argc) {
            P.max_relative_winding = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--max-samples") == 0 && i+1 < argc) {
            P.max_samples = (size_t)strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--dry-run") == 0) {
            P.dry_run = 1;
        } else {
            fprintf(stderr, "winding_evidence: unknown/incomplete arg %s\n", argv[i]);
            return 2;
        }
    }
    if (P.neighbor_distance == 0.0) P.neighbor_distance = 2.5 * P.pitch;
    int finite_axis = 1;
    for (int k = 0; k < 3; k++)
        finite_axis = finite_axis && isfinite(P.point_zyx[k]) && isfinite(P.dir_zyx[k]);
    if (!P.axis_point_set || !finite_axis || prepare_axis(&P) != 0 ||
        !isfinite(P.pitch) || P.pitch <= 0.0 ||
        !isfinite(P.winding_offset) ||
        !isfinite(P.neighbor_distance) || P.neighbor_distance <= 0.0 ||
        !isfinite(P.max_phase_residual) || P.max_phase_residual <= 0.0 ||
        P.max_phase_residual >= 0.5 || !isfinite(P.min_confidence) ||
        P.min_confidence < 0.0 || P.min_confidence > 1.0 ||
        P.max_relative_winding < 0 || P.max_samples < 8) {
        fprintf(stderr, "winding_evidence: invalid axis or scoring parameter\n");
        return 2;
    }

    Patch *patches = NULL; size_t npatch = 0, patch_cap = 0;
    if (scan_class(P.dataset, "verified_patches", 1,
                   &patches, &npatch, &patch_cap) != 0 ||
        scan_class(P.dataset, "unverified_patches", 0,
                   &patches, &npatch, &patch_cap) != 0) {
        free_patches(patches, npatch); return 1;
    }
    if (npatch == 0) {
        fprintf(stderr, "winding_evidence: no tifxyz patches under %s\n", P.dataset);
        free(patches); return 1;
    }
    qsort(patches, npatch, sizeof(patches[0]), compare_patch);
    for (size_t i = 1; i < npatch; i++) {
        if (strcmp(patches[i-1].id, patches[i].id) == 0) {
            fprintf(stderr, "winding_evidence: duplicate patch id %s\n", patches[i].id);
            free_patches(patches, npatch); return 1;
        }
    }
    for (size_t i = 0; i < npatch; i++) {
        if (load_patch(&patches[i], &P) != 0) {
            free_patches(patches, npatch); return 1;
        }
    }

    Relation *relations = NULL; size_t nrelation = 0, relation_cap = 0;
    for (size_t i = 0; i < npatch; i++) {
        for (size_t j = i + 1; j < npatch; j++) {
            if (bbox_gap(&patches[i], &patches[j]) > P.neighbor_distance) continue;
            if (append_relation(&relations, &nrelation, &relation_cap,
                                patches, (int)i, (int)j, &P) != 0) {
                free(relations); free_patches(patches, npatch); return 1;
            }
        }
    }
    size_t abs_ok = 0, same_ok = 0, relative_ok = 0;
    for (size_t i = 0; i < npatch; i++) abs_ok += patches[i].absolute_accepted != 0;
    for (size_t i = 0; i < nrelation; i++) if (relations[i].accepted) {
        if (relations[i].delta == 0) same_ok++; else relative_ok++;
    }
    uint64_t generation = evidence_hash(patches, npatch, relations, nrelation, &P);
    printf("winding_evidence: patches=%zu absolute=%zu same=%zu relative=%zu candidates=%zu generation=%016llx\n",
           npatch, abs_ok, same_ok, relative_ok,
           nrelation - same_ok - relative_ok, (unsigned long long)generation);
    int rc = 0;
    if (!P.dry_run) {
        if (write_evidence(patches, npatch, relations, nrelation, &P,
                           generation, abs_ok, same_ok, relative_ok) != 0 ||
            write_pcl(patches, npatch, relations, nrelation, &P, generation,
                      PCL_ABSOLUTE, "abs_winding.json") != 0 ||
            write_pcl(patches, npatch, relations, nrelation, &P, generation,
                      PCL_SAME, "same_windings.json") != 0 ||
            write_pcl(patches, npatch, relations, nrelation, &P, generation,
                      PCL_RELATIVE, "relative_windings.json") != 0) {
            fprintf(stderr, "winding_evidence: failed writing %s\n", P.out_dir);
            rc = 1;
        } else {
            printf("winding_evidence: wrote native evidence and Villa PCLs -> %s\n",
                   P.out_dir);
        }
    }
    free(relations); free_patches(patches, npatch);
    return rc;
}
