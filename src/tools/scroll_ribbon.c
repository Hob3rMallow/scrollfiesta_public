/*
 * scroll_ribbon.c -- slice / arc-length / StrokeStrip unwrap of a scroll sheet.
 *
 * Reads a welded (ideally single-component; see obj_biggest) .obj, computes the
 * winding field phi (src/flatten/unwrap.c, keep_phi), then runs the ribbon
 * parameterization (src/flatten/ribbon.c): constant-t slicing -> per-slice
 * polylines -> joint arc-length solve -> fitted quad ribbon.  u is true
 * per-slice arc length (vox), v is axial (vox). Writes:
 *
 *   <out_dir>/<id>_ribbon.obj/.vmesh  -- fitted ribbon grid
 *   <out_dir>/<id>_ribbon_stats.json  -- diagnostics
 *   <out_dir>/<id>_winding_index.f32  -- registered turns, input-vertex order
 *   <out_dir>/<id>_winding_confidence.f32 -- MRF/jump confidence
 *
 * Usage:
 *   scroll_ribbon <input.obj> <out_dir> [--id <id>]
 *                 [--axis-point z y x] [--axis-dir z y x] [--wrap-spacing B]
 *                 [--slice-h F] [--sample-h F] [--match-r F] [--match-ang F]
 *                 [--iters N] [--final-iters N] [--grid-u F]
 *   scroll_ribbon --selftest
 *
 * Default axis: Z line through the PHerc0139 manifest umbilicus
 * (y,x) = (3405, 2878) -- override with --axis-point/--axis-dir for other data.
 */
#include "../common/ves_platform.h"

#ifdef _OPENMP
#include "../common/ves_omp.h"
#endif

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/except.h"
#include "../common/mesh_bin.h"
#include "../common/obj_io.h"
#include "../flatten/lane_turn_mrf.h"
#include "../flatten/seam_cut.h"
#include "../flatten/unwrap.h"
#include "../flatten/ribbon.h"
#include "../flatten/sparse_solve.h"

/* Derive an id from the input path: basename minus extension. */
static void derive_id(const char *input_path, char *out, size_t cap)
{
    const char *base = ves_path_basename(input_path);
    size_t n = strlen(base);
    const char *dot = strrchr(base, '.');
    if (dot != NULL && dot != base) n = (size_t)(dot - base);
    if (n >= cap) n = cap - 1;
    memcpy(out, base, n);
    out[n] = '\0';
}

/* Escape a UTF-8 path/id for a JSON string value. In particular, absolute
 * Windows paths contain backslashes which otherwise turn sequences such as
 * `\v` into invalid JSON escapes. */
static void json_escape_copy(char *dst, size_t cap, const char *src)
{
    size_t n = 0;
    if (cap == 0) return;
    for (const unsigned char *p = (const unsigned char *)src; *p; p++) {
        const char *esc = NULL;
        char unicode[7];
        switch (*p) {
        case '\\': esc = "\\\\"; break;
        case '"':  esc = "\\\""; break;
        case '\b': esc = "\\b"; break;
        case '\f': esc = "\\f"; break;
        case '\n': esc = "\\n"; break;
        case '\r': esc = "\\r"; break;
        case '\t': esc = "\\t"; break;
        default:
            if (*p < 0x20) {
                snprintf(unicode, sizeof unicode, "\\u%04x", (unsigned)*p);
                esc = unicode;
            }
            break;
        }
        if (esc != NULL) {
            size_t len = strlen(esc);
            if (n + len >= cap) break;
            memcpy(dst + n, esc, len);
            n += len;
        } else {
            if (n + 1 >= cap) break;
            dst[n++] = (char)*p;
        }
    }
    dst[n] = '\0';
}

static int sibling_path(const char *input, const char *name,
                        char *out, size_t cap)
{
    const char *a = strrchr(input, '/');
    const char *b = strrchr(input, '\\');
    const char *slash = a != NULL && (b == NULL || a > b) ? a : b;
    size_t prefix = slash != NULL ? (size_t)(slash - input + 1) : 0;
    size_t n = strlen(name);
    if (prefix > cap || n >= cap - prefix) return -1;
    if (prefix != 0) memcpy(out, input, prefix);
    memcpy(out + prefix, name, n + 1);
    return 0;
}

/* Replace the final extension while retaining the input stem and directory. */
static int stem_suffix_path(const char *input, const char *suffix,
                            char *out, size_t cap)
{
    const char *a = strrchr(input, '/');
    const char *b = strrchr(input, '\\');
    const char *slash = a != NULL && (b == NULL || a > b) ? a : b;
    const char *dot = strrchr(input, '.');
    if (dot == NULL || (slash != NULL && dot < slash)) dot = input + strlen(input);
    size_t prefix = (size_t)(dot - input);
    size_t n = strlen(suffix);
    if (prefix >= cap || n >= cap - prefix) return -1;
    memcpy(out, input, prefix);
    memcpy(out + prefix, suffix, n + 1);
    return 0;
}

static void *read_raw_exact(Arena_T arena, const char *path,
                            size_t element_size, size_t count)
{
    FILE *file = fopen(path, "rb");
    void *data;
    int extra;
    if (file == NULL || element_size == 0 ||
        count > SIZE_MAX / element_size) {
        if (file != NULL) fclose(file);
        return NULL;
    }
    data = ARENA_ALLOC(arena, count * element_size);
    if (fread(data, element_size, count, file) != count) {
        fclose(file);
        return NULL;
    }
    extra = fgetc(file);
    if (extra != EOF || ferror(file) || fclose(file) != 0) return NULL;
    return data;
}

static int write_raw_exact(const char *path, const void *data,
                           size_t element_size, size_t count)
{
    if (path == NULL || data == NULL || element_size == 0 ||
        count > SIZE_MAX / element_size)
        return -1;
    ves_ensure_parent_dir(path);
    FILE *file = fopen(path, "wb");
    if (file == NULL) return -1;
    int ok = fwrite(data, element_size, count, file) == count && !ferror(file);
    if (fclose(file) != 0) ok = 0;
    return ok ? 0 : -1;
}

typedef struct {
    double q;
    double u;
    double local_u;
} ProjectiveUKnot;

static int projective_u_knot_compare(const void *pa, const void *pb)
{
    const ProjectiveUKnot *a = (const ProjectiveUKnot *)pa;
    const ProjectiveUKnot *b = (const ProjectiveUKnot *)pb;
    if (a->q != b->q) return a->q < b->q ? -1 : 1;
    if (a->u != b->u) return a->u < b->u ? -1 : 1;
    return 0;
}

/* Extend an immutable parent's metric chart through an exact overlap.  The
 * parent U(q) curve is copied on its supported q interval; beyond either end
 * the child's local arc-length map continues with a constant translation.
 * Boundary vertices themselves are copied bit-for-bit.  V has an exact
 * translation because both domains use the same axial geometry.
 *
 * This operation has no fitted scalar and no crop-wide vote: every datum is an
 * exact shared input vertex.  It is therefore compositional along a canonical
 * block chain. */
static int apply_projective_boundary_uv(
    Arena_T arena, float *uv, const float *verts, size_t nv,
    const UnwrapResult *ures, const UnwrapOpts *uopts)
{
    (void)verts;
    if (uopts == NULL || uopts->boundary_winding == NULL) return 0;
    if (uopts->boundary_u == NULL || uopts->boundary_v == NULL ||
        uopts->boundary_material == NULL || ures == NULL ||
        ures->winding_index == NULL)
        return -1;
    size_t nsample = 0;
    for (size_t i = 0; i < nv; i++)
        if (isfinite((double)uopts->boundary_winding[i])) nsample++;
    if (nsample == 0) return -1;
    ProjectiveUKnot *knot = (ProjectiveUKnot *)ARENA_ALLOC(
        arena, nsample * sizeof *knot);
    size_t at = 0;
    double vshift = 0.0;
    int have_vshift = 0;
    for (size_t i = 0; i < nv; i++) {
        double q = (double)uopts->boundary_winding[i];
        if (!isfinite(q)) continue;
        double pu = (double)uopts->boundary_u[i];
        double pv = (double)uopts->boundary_v[i];
        double vs = pv - (double)uv[i * 2 + 1];
        if (!isfinite(pu) || !isfinite(pv) || uopts->boundary_material[i] < 0)
            return -1;
        if (!have_vshift) { vshift = vs; have_vshift = 1; }
        else if (fabs(vs - vshift) > 0.01) {
            fprintf(stderr, "  projective UV: non-translational V boundary "
                    "at vertex %zu (%+.9g vs %+.9g)\n", i, vs, vshift);
            return -1;
        }
        knot[at].q = q;
        knot[at].u = pu;
        knot[at].local_u = (double)uv[i * 2];
        at++;
    }
    qsort(knot, at, sizeof *knot, projective_u_knot_compare);
    /* Collapse identical float-q samples.  Parent U is mathematically a
     * function of q; tolerate only float roundoff, never a material-scale
     * disagreement hidden behind an average. */
    size_t nknot = 0;
    for (size_t i = 0; i < at; ) {
        size_t j = i + 1;
        double su = knot[i].u, sl = knot[i].local_u;
        while (j < at && knot[j].q == knot[i].q) {
            /* q is persisted as float while U can change by several thousand
             * vox/turn.  One q ulp at the outer scroll is therefore a few
             * hundredths of a voxel in U; 0.25 remains far below one raster
             * sample and catches only a genuinely multivalued chart. */
            if (fabs(knot[j].u - knot[i].u) > 0.25) {
                fprintf(stderr, "  projective UV: parent U is multivalued at "
                        "q=%.9g (%.9g vs %.9g)\n",
                        knot[i].q, knot[i].u, knot[j].u);
                return -1;
            }
            su += knot[j].u;
            sl += knot[j].local_u;
            j++;
        }
        knot[nknot].q = knot[i].q;
        knot[nknot].u = su / (double)(j - i);
        knot[nknot].local_u = sl / (double)(j - i);
        if (nknot > 0 && knot[nknot].u + 0.25 < knot[nknot - 1].u) {
            fprintf(stderr, "  projective UV: parent U(q) is not monotone at "
                    "q=%.9g (%.9g after %.9g)\n", knot[nknot].q,
                    knot[nknot].u, knot[nknot - 1].u);
            return -1;
        }
        nknot++;
        i = j;
    }
    for (size_t i = 0; i < nv; i++) {
        double q = (double)ures->winding_index[i];
        double local_u = (double)uv[i * 2];
        double u = 0.0;
        if (!isfinite(q) || !isfinite(local_u)) return -1;
        if (nknot == 1 || q <= knot[0].q) {
            u = local_u + knot[0].u - knot[0].local_u;
        } else if (q >= knot[nknot - 1].q) {
            u = local_u + knot[nknot - 1].u -
                            knot[nknot - 1].local_u;
        } else {
            size_t lo = 0, hi = nknot - 1;
            while (hi - lo > 1) {
                size_t mid = lo + (hi - lo) / 2;
                if (knot[mid].q <= q) lo = mid;
                else hi = mid;
            }
            double d = knot[hi].q - knot[lo].q;
            double f = d > 0.0 ? (q - knot[lo].q) / d : 0.0;
            u = knot[lo].u + f * (knot[hi].u - knot[lo].u);
        }
        uv[i * 2] = (float)u;
        uv[i * 2 + 1] = (float)((double)uv[i * 2 + 1] + vshift);
    }
    /* Shared vertices are not merely close: the recursive certificate is the
     * parent's exact IEEE-754 sample at the overlap. */
    for (size_t i = 0; i < nv; i++) {
        if (!isfinite((double)uopts->boundary_winding[i])) continue;
        uv[i * 2] = uopts->boundary_u[i];
        uv[i * 2 + 1] = uopts->boundary_v[i];
    }
    fprintf(stderr, "  projective UV boundary: %zu exact samples, %zu U(q) "
            "knots, V shift=%+.6f\n", nsample, nknot, vshift);
    return 0;
}

/* Materialize the probabilistic winding certificate as a topology-preserving
 * VMESH that the fixed-topology metric solve can consume directly.  U is a
 * generously scaled winding coordinate (one turn = 1200 units); V is axial Z.
 * Geometry and triangle indices are unchanged. */
static int write_winding_ready_vmesh(
        Arena_T arena, const char *input_path, const char *out_dir,
        const char *id, const float *verts, size_t nv,
         const int32_t *faces, size_t nf, const UnwrapResult *ures,
         const UnwrapOpts *uopts,
         double winding_seconds, int conflict_correction_enabled)
{
    const double pixels_per_turn = 1200.0;   /* legacy, reported only */
    double effective_ppt = 1200.0;
    char vmesh_path[4096], report_path[4096];
    char json_input[4096], json_output[4096];
    float *uv;
    double winding_lo = 1e300, winding_hi = -1e300;
    double z_lo = 1e300, z_hi = -1e300;

    if (ures == NULL || ures->winding_index == NULL || nv == 0) return -1;
    uv = (float *)ARENA_ALLOC(arena, nv * 2 * sizeof(*uv));
    for (size_t i = 0; i < nv; i++) {
        double winding = (double)ures->winding_index[i];
        double z = (double)verts[i*3];
        if (!isfinite(winding) || !isfinite(z)) return -1;
        if (winding < winding_lo) winding_lo = winding;
        if (winding > winding_hi) winding_hi = winding;
        if (z < z_lo) z_lo = z;
        if (z > z_hi) z_hi = z;
    }
    /* U is ARC LENGTH along the sheet, not turns x a constant.
     *
     * This used to be (winding - winding_lo) * pixels_per_turn with
     * pixels_per_turn = 1200 -- the same number of voxels per turn at every
     * radius.  On the 4x5x5 the material runs from r = 13 to r = 438, so a
     * turn at the umbilicus really is ~82 vox of papyrus and a turn at the rim
     * ~2,750, and pricing both at 1200 oversampled the core by ~15x.  Measured
     * on the emitted ribbon: 2.20% of triangles carried under a quarter of the
     * sheet's median 3-D area per unit of UV area, ALL of them squashed, none
     * stretched, at a median radius of 45 vox against 242 for the sheet.
     *
     * Unwrap already integrates the empirical radius profile against lifted
     * phase, so carry that.  It is a monotone reparameterization of the same
     * certificate: equal winding still gives equal U and the turn ordering is
     * untouched.  Orientation is taken from the winding index so a
     * negative-sense scroll still runs U in the certificate's direction. */
    {
        int flip = 0;
        double a_lo = 0.0, a_hi = 0.0, u_hi = 0.0;
        size_t i_lo = 0, i_hi = 0;
        for (size_t i = 0; i < nv; i++) {
            double w = (double)ures->winding_index[i];
            if (w <= winding_lo) i_lo = i;
            if (w >= winding_hi) i_hi = i;
            if ((double)ures->uv[i*2] > u_hi) u_hi = (double)ures->uv[i*2];
        }
        a_lo = (double)ures->uv[i_lo*2];
        a_hi = (double)ures->uv[i_hi*2];
        flip = a_lo > a_hi;
        for (size_t i = 0; i < nv; i++) {
            double a = (double)ures->uv[i*2];
            if (!isfinite(a)) return -1;
            uv[i*2] = (float)(flip ? u_hi - a : a);
            uv[i*2+1] = (float)((double)verts[i*3] - z_lo);
        }
        effective_ppt = winding_hi > winding_lo
                      ? u_hi / (winding_hi - winding_lo) : pixels_per_turn;
        fprintf(stderr,
                "  winding certificate U: arc length, span %.1f vox over "
                "%.3f turns (mean %.1f vox/turn; the fixed %.0f vox/turn this "
                "replaces oversampled the umbilicus)%s\n",
                u_hi, winding_hi - winding_lo, effective_ppt,
                 pixels_per_turn, flip ? " [flipped to winding sense]" : "");
    }
    if (apply_projective_boundary_uv(arena, uv, verts, nv, ures, uopts) != 0)
        return -1;

    snprintf(vmesh_path, sizeof vmesh_path, "%s/%s_winding.vmesh",
             out_dir, id);
    ves_ensure_parent_dir(vmesh_path);
    if (MeshBin_write(vmesh_path, verts, nv, faces, nf, uv) != 0) return -1;

    snprintf(report_path, sizeof report_path,
             "%s/%s_winding_visualization.json", out_dir, id);
    json_escape_copy(json_input, sizeof json_input, input_path);
    json_escape_copy(json_output, sizeof json_output, vmesh_path);
    FILE *report = fopen(report_path, "w");
    if (report == NULL) return -1;
    fprintf(report,
        "{\n"
        "  \"schema\": \"vesuvius-winding-certificate-v2\",\n"
        "  \"source\": \"%s\",\n"
        "  \"output\": \"%s\",\n"
        "  \"vertices\": %zu,\n"
        "  \"faces\": %zu,\n"
        "  \"pixels_per_turn\": %.1f,\n"
        "  \"winding_range\": [%.9g, %.9g],\n"
        "  \"winding_span\": %.9g,\n"
        "  \"z_range\": [%.9g, %.9g],\n"
        "  \"registration\": {\n"
        "    \"mesh_components\": %d,\n"
        "    \"continuation_components\": %zu,\n"
        "    \"relation_components\": %zu,\n"
        "    \"relations\": %zu,\n"
        "    \"eligible_relations\": %zu,\n"
        "    \"residual_conflicts\": %zu,\n"
        "    \"continuation_satisfaction\": %.9g,\n"
        "    \"order_satisfaction\": %.9g,\n"
        "    \"mrf_rounds\": %zu,\n"
        "    \"mrf_label_changes\": %zu,\n"
        "    \"mrf_abstained_sites\": %zu,\n"
        "    \"mrf_energy_before\": %.9g,\n"
        "    \"mrf_energy_after\": %.9g,\n"
         "    \"mrf_mean_confidence\": %.9g,\n"
         "    \"conflict_resolution\": { \"enabled\": %s, "
         "\"rounds\": %zu, "
        "\"bins_before\": %zu, \"bins_after\": %zu, "
        "\"losing_claims\": %zu, \"forbidden_assignments\": %zu, "
        "\"winner_locks\": %zu, "
        "\"label_changes\": %zu, \"converged\": %s }\n"
        "  },\n"
        "  \"timing_sec\": { \"winding\": %.3f }\n"
        "}\n",
        json_input, json_output, nv, nf, effective_ppt,
        winding_lo, winding_hi, winding_hi - winding_lo, z_lo, z_hi,
        ures->n_components, ures->winding_continuation_components,
        ures->winding_relation_components, ures->winding_relations,
        ures->winding_eligible_relations, ures->winding_relation_conflicts,
        ures->continuation_satisfaction, ures->order_satisfaction,
        ures->winding_mrf_rounds, ures->winding_mrf_label_changes,
        ures->winding_mrf_abstained_sites, ures->winding_mrf_energy_before,
         ures->winding_mrf_energy_after, ures->winding_mrf_mean_confidence,
         conflict_correction_enabled ? "true" : "false",
         ures->winding_mrf_conflict_rounds,
        ures->winding_mrf_conflict_bins_before,
        ures->winding_mrf_conflict_bins_after,
        ures->winding_mrf_conflict_losing_claims,
        ures->winding_mrf_conflict_exclusions,
        ures->winding_mrf_conflict_winner_locks,
        ures->winding_mrf_conflict_label_changes,
        ures->winding_mrf_conflict_converged ? "true" : "false",
        winding_seconds);
    if (fclose(report) != 0) return -1;
    fprintf(stderr,
            "  winding-only: wrote %s + %s (topology preserved)\n",
            vmesh_path, report_path);
    return 0;
}
/* Load the concat's per-cube table (mesh_cubes.tsv beside the input) and
 * expand it to a per-vertex source-cube id.  Row order IS concat order, so
 * the cumulative nv column maps vertex ranges to cubes; a table whose total
 * does not match the mesh is stale and refused.  Arms the register's exact
 * OVERLAP relation family on halo-overlap piles; harmless (NULL) elsewhere. */
static int32_t *load_vertex_cube_table(Arena_T arena, const char *mesh_path,
                                       size_t nv)
{
    char path[2048];
    size_t len = strlen(mesh_path);
    const char *slash = strrchr(mesh_path, '/');
    const char *bslash = strrchr(mesh_path, 92);
    const char *cut = slash > bslash ? slash : bslash;
    size_t dirlen = cut != NULL ? (size_t)(cut - mesh_path) + 1 : 0;
    FILE *f = NULL;
    int32_t *out = NULL;
    size_t at = 0, row = 0;
    char line[4096];
    if (dirlen + 16 >= sizeof path) return NULL;
    memcpy(path, mesh_path, dirlen);
    memcpy(path + dirlen, "mesh_cubes.tsv", 15);
    f = fopen(path, "rb");
    if (f == NULL) return NULL;
    out = (int32_t *)ARENA_ALLOC(arena, (size_t)(nv * sizeof(int32_t)));
    if (fgets(line, sizeof line, f) == NULL) { fclose(f); return NULL; }
    size_t empty_rows = 0;
    while (fgets(line, sizeof line, f) != NULL) {
        unsigned long long cnv = 0;
        const char *tab = strchr(line, 9);
        if (tab == NULL) continue;
        cnv = strtoull(tab + 1, NULL, 10);
        /* A cube whose vertices were all routed away (core curl) is a
         * legitimate EMPTY row, not a broken table: skip it but keep its
         * row id so provenance stays aligned.  Refusing here silently
         * disarmed the overlap family on every fixed-axis pile lane
         * (measured 2026-09-03: 4 of 100 empty rows -> ovl=0). */
        if (cnv == 0) { empty_rows++; row++; continue; }
        if (at + cnv > nv) {
            fprintf(stderr, "  vertex cube table: row %zu overruns the mesh "
                    "(%llu verts at %zu of %zu); overlap family NOT armed\n",
                    row, cnv, at, nv);
            fclose(f);
            return NULL;
        }
        for (size_t v = 0; v < cnv; v++) out[at + v] = (int32_t)row;
        at += cnv;
        row++;
    }
    fclose(f);
    if (at != nv) {
        fprintf(stderr, "  vertex cube table: %zu verts listed vs %zu in the "
                "mesh; overlap family NOT armed\n", at, nv);
        return NULL;
    }
    fprintf(stderr, "  vertex cube table: %zu cubes (%zu empty) over %zu verts "
            "-> overlap family ARMED (%s)",
            row, empty_rows, nv, path);
    fputc(10, stderr);
    return out;
}



int main(int argc, char *argv[])
{
    /* Unbuffered stderr: this is a long-running diagnostic tool whose per-phase
     * progress lines must appear live (block-buffering hid them until exit,
     * blinding memory/hang monitoring of the multi-minute phases). */
    setvbuf(stderr, NULL, _IONBF, 0);

    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) {
        int f0 = Sparse_selftest() != 0;
        int f1 = Ribbon_selftest();
        int f2 = Unwrap_selftest();
        int f3 = SeamCut_selftest();
        int f4 = LaneTurn_selftest();
        int f5 = RibbonGridIO_selftest();
        int fails = f0 + f1 + f2 + f3 + f4 + f5;
        fprintf(stderr, "\nscroll_ribbon --selftest: %s (sparse=%d ribbon=%d unwrap=%d seamcut=%d laneturn=%d gridio=%d failures)\n",
                fails == 0 ? "PASS" : "FAIL", f0, f1, f2, f3, f4, f5);
        return fails == 0 ? 0 : 1;
    }

    if (argc < 3) {
        fprintf(stderr,
            "Usage: %s <input.obj> <out_dir> [--id <id>]\n"
            "          [--axis-point z y x] [--axis-dir z y x] [--wrap-spacing B]\n"
            "          [--no-winding-field] [--winding-field-eps F] [--winding-field-beta F]\n"
            "          [--no-winding-conflict-correction]\n"
            "                            (diagnostic: retain initial MRF certificate)\n"
            "          [--boundary-winding FILE --boundary-u FILE\n"
            "           --boundary-v FILE --boundary-material FILE]\n"
            "                            (exact canonical-parent overlap)\n"
            "          [--slice-h F] [--sample-h F] [--match-r F] [--match-ang F]\n"
            "          [--iters N] [--final-iters N] [--metric-iters N]\n"
            "          [--metric-weight F] [--threads N] [--no-amg]\n"
            "          [--winding-only] (emit certificate VMESH + sidecars; skip ribbon solve)\n"
            "          [--dump-gmg-levels] [--dump-gmg-level-objs]\n"
            "                            (fitted VMESH immediately after every\n"
            "                             robust solve round, then per final\n"
            "                             coarse-to-fine level; optional text OBJ)\n"
            "          [--grid-u F]\n"
            "          [--trust-gauge | --solve-carried-u]\n"
            "                            (reuse source winding/material sidecars;\n"
            "                             trust freezes U, solve refines U with\n"
            "                             the geometry-driven Stage-C solve)\n"
            "          [--coarse-chain-seed] (with --solve-carried-u: solve\n"
            "                             chain offsets before fine samples)\n"
            "          [--preserve-input-topology]\n"
            "                            (parameterize the input ribbon in place:\n"
            "                             identical vertices/faces, new UV only)\n"
            "          [--metric-project-only]\n"
            "                            (diagnostic legacy exact-row projection;\n"
            "                             requires --preserve-input-topology)\n"
            "          [--vmesh-only]    (skip the redundant text OBJ output)\n"
            "          [--peel-layers N] (claims-mode collision depth, 1..64; default 2; checked memory budget)\n"
            "          [--peel-layer N] (store only zero-based projective-owner layer N)\n"
            "          [--promote-peel-min-share F]\n"
            "                            (pack collision-peel U runs with at least\n"
            "                             F times the primary vertex count beside it)\n"
            "          [--promote-peel-min-vertices N]\n"
            "                            (also promote runs with at least N vertices;\n"
            "                             the smaller positive threshold wins)\n"
             "          [--fit-cover-width] (diagnostic: retain cover grid extent)\n"
            "          [--verify-fit-width] (diagnostic: compare compact/padded fits)\n"
            "          [--no-sever]     (skip fusion-handle severing)\n"
            "          [--component-global] (whole welded-component fitted-grid\n"
            "                            consensus; resolve duplicate u claims by\n"
            "                            authoritative winding instead of chain order)\n"
            "          [--coherent-claims] (default with --trust-gauge + claims:\n"
            "                            retain UV frame, use measured cross-row\n"
            "                            continuity and balanced U/V break costs)\n"
            "          [--legacy-owner-claims] (diagnostic: reproduce the old\n"
            "                            slice-independent owner selection)\n"
            "          [--ownership claims|construction] (component-global mode:\n"
            "                            'claims' = frozen whole-row claimant path,\n"
            "                            DEFAULT;\n"
            "                            'construction' = certified merges + frozen\n"
            "                            row-invariant fusion, 2026-08-13)\n"
            "          [--discard-conflicting-claims]\n"
            "                            (claims mode: every multiply claimed UV cell\n"
            "                             becomes unsupported minimal-surface input)\n"
            "          [--radial-bridge-cut | --no-radial-bridge-cut]\n"
            "                            (force legacy radial slice splitting on/off;\n"
            "                             default auto: off with a winding reference)\n"
            "       %s --selftest\n", argv[0], argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    const char *out_dir    = argv[2];
    char id[512];
    derive_id(input_path, id, sizeof id);

    /* default: Z axis through the PHerc0139 manifest umbilicus */
    UnwrapOpts uopts;
    memset(&uopts, 0, sizeof uopts);
    uopts.axis_mode = UNWRAP_AXIS_EXTERNAL;
    uopts.axis_point[0] = 0.0f; uopts.axis_point[1] = 3405.0f; uopts.axis_point[2] = 2878.0f;
    uopts.axis_dir[0] = 1.0f;   uopts.axis_dir[1] = 0.0f;      uopts.axis_dir[2] = 0.0f;
    uopts.keep_phi = 1;          /* Ribbon consumes this graph-traced scaffold. */

    RibbonOpts ropts;
    RibbonOpts_default(&ropts);
    ropts.direct_ribbon = 1;
    int sever_handles = 1;        /* genus > 0 collapses the winding chart */
    int dump_gmg_levels = 0;
    int dump_gmg_level_objs = 0;
    int preserve_input_topology = 0;
    int metric_iters_explicit = 0;
    int metric_project_only_requested = 0;
    int trust_gauge = 0;   /* --trust-gauge: parameterize the carried certificate
                            * U directly (solve_reference_u==0 fast path) --
                            * skip the pair graph / winding / metric ADMM solve.
                            * The streaming block fit and M1 use this; the slow
                            * metric solve is being moved to geometric multigrid
                            * on coarse proxies (see project_fullscroll_block_
                            * architecture). */
    int winding_only = 0;
    int solve_carried_u = 0; /* same measured source/correspondence as trust,
                             * but U is an initializer, not a frozen result */
    int legacy_owner_claims = 0;
    int write_text_obj = 1;
    double promote_peel_min_share = 0.0;
    size_t promote_peel_min_vertices = 0;
    char gmg_level_prefix[4096];
    const char *boundary_winding_path = NULL;
    const char *boundary_u_path = NULL;
    const char *boundary_v_path = NULL;
    const char *boundary_material_path = NULL;

    int no_cube_table = 0;
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--id") == 0 && i + 1 < argc) {
            snprintf(id, sizeof id, "%s", argv[++i]);
        } else if (strcmp(argv[i], "--axis-point") == 0 && i + 3 < argc) {
            uopts.axis_point[0]=(float)atof(argv[++i]);
            uopts.axis_point[1]=(float)atof(argv[++i]);
            uopts.axis_point[2]=(float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--axis-dir") == 0 && i + 3 < argc) {
            uopts.axis_dir[0]=(float)atof(argv[++i]);
            uopts.axis_dir[1]=(float)atof(argv[++i]);
            uopts.axis_dir[2]=(float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--wrap-spacing") == 0 && i + 1 < argc) {
            uopts.wrap_spacing = atof(argv[++i]);
        } else if (strcmp(argv[i], "--no-cube-table") == 0) {
            no_cube_table = 1;    /* A/B switch: leave the OVERLAP family unarmed */
        } else if (strcmp(argv[i], "--boundary-winding") == 0 &&
                   i + 1 < argc) {
            boundary_winding_path = argv[++i];
        } else if (strcmp(argv[i], "--boundary-u") == 0 && i + 1 < argc) {
            boundary_u_path = argv[++i];
        } else if (strcmp(argv[i], "--boundary-v") == 0 && i + 1 < argc) {
            boundary_v_path = argv[++i];
        } else if (strcmp(argv[i], "--boundary-material") == 0 &&
                   i + 1 < argc) {
            boundary_material_path = argv[++i];
        } else if (strcmp(argv[i], "--winding-sense") == 0 && i + 1 < argc) {
            int ws = atoi(argv[++i]);
            uopts.winding_sense = ws > 0 ? 1 : (ws < 0 ? -1 : 0);
            ropts.winding_sense = uopts.winding_sense;
        } else if (strcmp(argv[i], "--no-winding-field") == 0) {
            uopts.winding_field_mode = -1;
        } else if (strcmp(argv[i], "--winding-field-eps") == 0 &&
                   i + 1 < argc) {
            uopts.winding_field_epsilon = atof(argv[++i]);
        } else if (strcmp(argv[i], "--winding-field-beta") == 0 &&
                   i + 1 < argc) {
            uopts.winding_field_beta = atof(argv[++i]);
        } else if (strcmp(argv[i], "--no-winding-conflict-correction") == 0) {
            uopts.winding_conflict_mode = -1;
        } else if (strcmp(argv[i], "--slice-h") == 0 && i + 1 < argc) {
            ropts.slice_h = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--sample-h") == 0 && i + 1 < argc) {
            ropts.sample_h = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--match-r") == 0 && i + 1 < argc) {
            ropts.match_r = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--match-ang") == 0 && i + 1 < argc) {
            ropts.match_ang_deg = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--iters") == 0 && i + 1 < argc) {
            ropts.relax_iters = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--final-iters") == 0 && i + 1 < argc) {
            ropts.final_iters = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--metric-iters") == 0 && i + 1 < argc) {
            ropts.metric_iters = atoi(argv[++i]);
            metric_iters_explicit = 1;
        } else if (strcmp(argv[i], "--metric-weight") == 0 && i + 1 < argc) {
            ropts.metric_weight = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--scaffold-solve") == 0) {
            /* fit-for-reparameterization: Stage-C u is a scaffold the
             * pipeline re-solves, so one round each and a capped ADMM */
            ropts.scaffold_solve = 1;
            ropts.relax_iters = 1;
            ropts.final_iters = 1;
            ropts.metric_iters = 0;
            metric_iters_explicit = 1;
        } else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            ropts.solve_threads = atoi(argv[++i]);
            /* --threads used to bound only the Stage-C solve workers, so
             * every other parallel region in the ribbon lane still opened
             * an OpenMP team of one thread per logical core.  A `--threads 8`
             * 10x fit measured 52 cores busy.  Cap the whole process. */
#ifdef _OPENMP
            if (ropts.solve_threads > 0)
                omp_set_num_threads(ropts.solve_threads);
#endif
        } else if (strcmp(argv[i], "--no-amg") == 0) {
            ropts.solve_amg = 0;
        } else if (strcmp(argv[i], "--dump-gmg-levels") == 0) {
            dump_gmg_levels = 1;
        } else if (strcmp(argv[i], "--dump-gmg-level-objs") == 0) {
            dump_gmg_levels = 1;
            dump_gmg_level_objs = 1;
        } else if (strcmp(argv[i], "--grid-u") == 0 && i + 1 < argc) {
            ropts.grid_u = (float)atof(argv[++i]);
        } else if (strcmp(argv[i], "--preserve-input-topology") == 0) {
            preserve_input_topology = 1;
        } else if (strcmp(argv[i], "--metric-project-only") == 0) {
            metric_project_only_requested = 1;
        } else if (strcmp(argv[i], "--trust-gauge") == 0) {
            trust_gauge = 1;
        } else if (strcmp(argv[i], "--solve-carried-u") == 0) {
            solve_carried_u = 1;
        } else if (strcmp(argv[i], "--coarse-chain-seed") == 0) {
            ropts.coarse_chain_seed = 1;
        } else if (strcmp(argv[i], "--coherent-claims") == 0) {
            ropts.coherent_claims = 1;
        } else if (strcmp(argv[i], "--legacy-owner-claims") == 0) {
            legacy_owner_claims = 1;
        } else if (strcmp(argv[i], "--winding-only") == 0) {
            winding_only = 1;
        } else if (strcmp(argv[i], "--radial-site-split") == 0) {
            uopts.radial_site_split = 1;   /* R1 experiment; see unwrap.h */
        } else if (strcmp(argv[i], "--vmesh-only") == 0) {
            write_text_obj = 0;
        } else if (strcmp(argv[i], "--peel-layers") == 0 &&
                   i + 1 < argc) {
            char *end = NULL;
            const char *value = argv[++i];
            long layers = strtol(value, &end, 10);
            if (end == value || *end != '\0' || layers < 1 || layers > RIBBON_MAX_PEEL_LAYERS) {
                fprintf(stderr,
                        "ERROR: --peel-layers must be an integer in [1,%d]\n", RIBBON_MAX_PEEL_LAYERS);
                return 1;
            }
            ropts.peel_layers = (int)layers;
        } else if (strcmp(argv[i], "--peel-layer") == 0 && i + 1 < argc) {
            char *end = NULL;
            const char *value = argv[++i];
            long layer = strtol(value, &end, 10);
            if (end == value || *end != '\0' || layer < 0 || layer >= RIBBON_MAX_PEEL_LAYERS) {
                fprintf(stderr, "ERROR: --peel-layer must be an integer in [0,%d]\n", RIBBON_MAX_PEEL_LAYERS-1);
                return 1;
            }
            ropts.peel_layer = (int)layer;
        } else if (strcmp(argv[i], "--promote-peel-min-share") == 0 &&
                   i + 1 < argc) {
            promote_peel_min_share = strtod(argv[++i], NULL);
            if (!isfinite(promote_peel_min_share) ||
                promote_peel_min_share < 0.0) {
                fprintf(stderr,
                        "ERROR: --promote-peel-min-share must be finite and "
                        "non-negative\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--promote-peel-min-vertices") == 0 &&
                   i + 1 < argc) {
            char *end = NULL;
            unsigned long long value = strtoull(argv[++i], &end, 10);
            if (end == argv[i] || *end != '\0' || value == 0 ||
                value > (unsigned long long)SIZE_MAX) {
                fprintf(stderr,
                        "ERROR: --promote-peel-min-vertices must be a "
                        "positive integer\n");
                return 1;
            }
            promote_peel_min_vertices = (size_t)value;
        } else if (strcmp(argv[i], "--fit-cover-width") == 0) {
            ropts.fit_cover_width = 1;
        } else if (strcmp(argv[i], "--verify-fit-width") == 0) {
            ropts.verify_fit_width = 1;
        } else if (strcmp(argv[i], "--no-sever") == 0) {
            sever_handles = 0;
        } else if (strcmp(argv[i], "--component-global") == 0) {
            ropts.component_global = 1;
        } else if (strcmp(argv[i], "--discard-conflicting-claims") == 0) {
            ropts.discard_conflicting_claims = 1;
        } else if (strcmp(argv[i], "--ownership") == 0 && i + 1 < argc) {
            ++i;
            if (strcmp(argv[i], "claims") == 0) {
                ropts.ownership_construction = 0;
            } else if (strcmp(argv[i], "construction") == 0) {
                ropts.ownership_construction = 1;
            } else {
                fprintf(stderr,
                        "ERROR: --ownership expects claims|construction, "
                        "got '%s'\n", argv[i]);
                return 1;
            }
        } else if (strcmp(argv[i], "--radial-bridge-cut") == 0) {
            ropts.radial_bridge_cut = 1;
        } else if (strcmp(argv[i], "--no-radial-bridge-cut") == 0) {
            ropts.radial_bridge_cut = -1;
        } else {
            fprintf(stderr, "WARN: ignoring unknown arg '%s'\n", argv[i]);
        }
    }
    memcpy(ropts.axis_point, uopts.axis_point, sizeof ropts.axis_point);
    memcpy(ropts.axis_dir,   uopts.axis_dir,   sizeof ropts.axis_dir);
    if (ropts.coarse_chain_seed && !solve_carried_u) {
        fprintf(stderr, "ERROR: --coarse-chain-seed requires --solve-carried-u\n");
        return 1;
    }
    if (solve_carried_u) {
        if (trust_gauge || winding_only || preserve_input_topology ||
            metric_project_only_requested || legacy_owner_claims) {
            fprintf(stderr, "ERROR: --solve-carried-u is a fitted-grid metric "
                    "solve, incompatible with frozen-U, winding-only, "
                    "fixed-topology, or legacy-owner modes\n");
            return 1;
        }
        /* Reuse the complete sidecar loader and preserve source topology.
         * The only coordinate authority changed below is solve_reference_u. */
        trust_gauge = 1;
    }
    {
        int boundary_args = (boundary_winding_path != NULL) +
                            (boundary_u_path != NULL) +
                            (boundary_v_path != NULL) +
                            (boundary_material_path != NULL);
        if (boundary_args != 0 && boundary_args != 4) {
            fprintf(stderr, "ERROR: projective boundary requires all four of "
                    "--boundary-winding/--boundary-u/--boundary-v/"
                    "--boundary-material\n");
            return 1;
        }
        if (boundary_args == 4) {
            /* Boundary arrays index the input topology.  A topology-changing
             * sever pass would invalidate their exact vertex correspondence. */
            sever_handles = 0;
        }
    }
    if (preserve_input_topology) {
        /* This is the controlled quad-ribbon parameterization experiment.
         * Geometry and connectivity are the independent variable; only
         * per-vertex UV may change.  Do not sever the input, rebuild a fitted
         * grid, select claimants, or emit reconstructed faces.
         *
         * The input constant-V rows become the StrokeStrip observations.  The
         * carried coordinate is the global initializer/correspondence evidence;
         * the staged solve first establishes cross-sections, then freezes them
         * while ramping physical stroke-edge unit speed. */
        sever_handles = 0;
        ropts.direct_ribbon = 0;
        ropts.fit_ribbon = 0;
        ropts.preserve_input_rows = 1;
        ropts.pin_orient = 1;
        ropts.stitch_solve_gauges = 0;
        ropts.metric_project_only = metric_project_only_requested;
        if (!metric_iters_explicit && !metric_project_only_requested)
            ropts.metric_iters = 3;
        if (metric_project_only_requested)
            ropts.metric_iters = 0;
        if (dump_gmg_levels) {
            fprintf(stderr,
                    "ERROR: solve checkpoints currently rebuild fitted grids; "
                    "they are deliberately disabled by "
                    "--preserve-input-topology\n");
            return 1;
        }
    }
    /* A trusted winding certificate already fixed the source topology and
     * integer gauge.  Cutting handles here would change its vertex order and
     * invalidate every carried sidecar before the fit even starts. */
    if (trust_gauge) {
        sever_handles = 0;
        ropts.projective_grid = 1;
        if (promote_peel_min_share > 0.0 ||
            promote_peel_min_vertices > 0) {
            fprintf(stderr,
                    "  trusted projective atlas: peel promotion disabled; "
                    "secondary claims remain in extras\n");
            promote_peel_min_share = 0.0;
            promote_peel_min_vertices = 0;
        }
    }
    if (legacy_owner_claims && ropts.coherent_claims) {
        fprintf(stderr, "ERROR: coherent and legacy-owner claims are mutually exclusive\n");
        return 1;
    }
    /* Trusting the coordinate frame does not imply independent row ownership.
     * Keep that obsolete policy available only as an explicit control. */
    if (trust_gauge && ropts.component_global && !ropts.ownership_construction &&
        !ropts.discard_conflicting_claims && !preserve_input_topology &&
        !winding_only && !legacy_owner_claims)
        ropts.coherent_claims = 1;
    if ((ropts.coherent_claims || legacy_owner_claims) &&
        (!trust_gauge || !ropts.component_global || ropts.ownership_construction ||
         ropts.discard_conflicting_claims || preserve_input_topology || winding_only)) {
        fprintf(stderr, "ERROR: coherent/legacy-owner selection requires --trust-gauge "
                        "--component-global with retained claims and fitted-grid output\n");
        return 1;
    }
    if (metric_project_only_requested && !preserve_input_topology) {
        fprintf(stderr,
                "ERROR: --metric-project-only requires "
                "--preserve-input-topology\n");
        return 1;
    }
    if (ropts.discard_conflicting_claims && !ropts.component_global) {
        fprintf(stderr,
                "ERROR: --discard-conflicting-claims requires --component-global\n");
        return 1;
    }
    if (ropts.discard_conflicting_claims && ropts.ownership_construction) {
        fprintf(stderr,
                "ERROR: --discard-conflicting-claims applies to --ownership claims\n");
        return 1;
    }
    ropts.wrap_spacing = (float)uopts.wrap_spacing;   /* --wrap-spacing override */
    if (dump_gmg_levels) {
        snprintf(gmg_level_prefix, sizeof gmg_level_prefix,
                 "%s/%s", out_dir, id);
        ropts.solve_level_prefix = gmg_level_prefix;
        ropts.solve_level_write_obj = dump_gmg_level_objs;
    }

    fprintf(stderr, "scroll_ribbon: input=%s out_dir=%s id=%s\n"
            "  axis point(zyx)=(%.1f,%.1f,%.1f) dir(zyx)=(%.2f,%.2f,%.2f)\n"
            "  slice_h=%.2f sample_h=%.2f match_r=%.2f match_ang=%.1f "
            "iters=%d+%d+%d metric_weight=%.3g threads=%d amg=%d "
            "grid_u=%.2f component_global=%d "
            "ownership=%s conflict_policy=%s winding_conflict=%s "
            "radial_bridge_cut=%s "
            "wrap_spacing=%.3f%s gmg_levels=%d%s output=%s\n",
            input_path, out_dir, id,
            (double)uopts.axis_point[0], (double)uopts.axis_point[1], (double)uopts.axis_point[2],
            (double)uopts.axis_dir[0], (double)uopts.axis_dir[1], (double)uopts.axis_dir[2],
            (double)ropts.slice_h, (double)ropts.sample_h, (double)ropts.match_r,
            (double)ropts.match_ang_deg, ropts.relax_iters, ropts.final_iters,
            ropts.metric_iters, (double)ropts.metric_weight,
            ropts.solve_threads, ropts.solve_amg,
            (double)ropts.grid_u, ropts.component_global,
            ropts.ownership_construction ? "construction" : "claims",
            ropts.discard_conflicting_claims ? "discard-to-minimal" : "select",
            uopts.winding_conflict_mode < 0 ? "initial-only" : "corrected",
            ropts.radial_bridge_cut > 0 ? "forced" :
            ropts.radial_bridge_cut < 0 ? "off" : "auto",
            (double)ropts.wrap_spacing,
            ropts.wrap_spacing > 0.0f ? " (pinned)" : " (auto)",
            dump_gmg_levels,
            dump_gmg_level_objs ? " (text OBJ too)" : "",
            winding_only ? "winding certificate only" :
             preserve_input_topology ? "input topology + solved UV" :
                                       "fitted ribbon grid");
    if (promote_peel_min_share > 0.0)
        fprintf(stderr,
                "  collision-peel promotion: minimum %.6g of primary "
                "fitted vertices\n", promote_peel_min_share);
    if (promote_peel_min_vertices > 0)
        fprintf(stderr,
                "  collision-peel absolute promotion criterion: %zu vertices\n",
                promote_peel_min_vertices);
    if (ropts.peel_layers > 0)
        fprintf(stderr, "  collision-peel layers: %d (explicit)\n",
                ropts.peel_layers);
    if (ropts.coherent_claims)
        fprintf(stderr, "  claim policy: carried UV frame; measured rolling-row "
                        "continuity; balanced U/V breaks; crossing-aware U fill\n");

    double t0 = ves_clock_sec();
    Arena_T arena = Arena_new();
    int ok = 0;
    double t_load = 0, t_sever = 0, t_unwrap = 0,
           t_ribbon = 0, t_write = 0;
    UnwrapResult ures;
    RibbonResult rres;
    memset(&ures, 0, sizeof ures);
    memset(&rres, 0, sizeof rres);

    TRY
        /* --- Load. --- */
        float *verts = NULL; int32_t *faces = NULL;
        double *input_reference_u = NULL;
        double *input_reference_v = NULL;
        float *input_reference_u_confidence = NULL;
        float *certificate_winding = NULL;
        float *certificate_phi = NULL;
        int32_t *certificate_material = NULL;
        int32_t *certificate_relation = NULL;
        int32_t *certificate_mesh_component = NULL;
        int32_t *certificate_owner_block = NULL;
        size_t certificate_material_count = 0;
        size_t certificate_relation_count = 0;
        size_t certificate_mesh_component_count = 0;
        float *sidecar_phase = NULL;
        int32_t *sidecar_material = NULL;
        char sidecar_phase_path[4096], sidecar_material_path[4096];
        size_t nv = 0, nf = 0;
        double ta = ves_clock_sec();
        {
            size_t input_len = strlen(input_path);
            int is_vmesh = input_len >= 6 &&
                _stricmp(input_path + input_len - 6, ".vmesh") == 0;
            if (is_vmesh) {
                MeshBinData mesh;
                memset(&mesh, 0, sizeof mesh);
                if (MeshBin_read_arena(arena, input_path, &mesh) != 0) {
                    fprintf(stderr, "ERROR: cannot read %s\n", input_path);
                    RAISE(IO_Failed);
                }
                verts = mesh.verts; nv = mesh.nv;
                faces = mesh.faces; nf = mesh.nf;
                /* Both fixed-topology metric refinement and a quadribbon
                 * scaffold fit consume carried winding U.  The latter used
                 * to advertise Stage-C U as a scaffold while silently
                 * discarding it here, because only preserve-input-topology
                 * armed the reference.  That made a direct fit from a
                 * winding-only VMESH recompute its own gauge instead of
                 * fitting the supplied certificate. */
                if ((preserve_input_topology || ropts.scaffold_solve ||
                     trust_gauge) &&
                    mesh.uv != NULL) {
                    double ulo = 1e300, uhi = -1e300;
                    input_reference_u = (double *)ARENA_ALLOC(
                        arena, nv * sizeof(*input_reference_u));
                    input_reference_v = (double *)ARENA_ALLOC(
                        arena, nv * sizeof(*input_reference_v));
                    for (size_t vi = 0; vi < nv; vi++) {
                        double u = (double)mesh.uv[vi * 2];
                        double v = (double)mesh.uv[vi * 2 + 1];
                        if (!isfinite(u) || !isfinite(v)) {
                            fprintf(stderr,
                                    "ERROR: non-finite input UV at vertex %zu\n",
                                    vi);
                            RAISE(IO_Failed);
                        }
                        input_reference_u[vi] = u;
                        input_reference_v[vi] = v;
                        if (u < ulo) ulo = u;
                        if (u > uhi) uhi = u;
                    }
                    ropts.reference_u = input_reference_u;
                    ropts.reference_v = input_reference_v;
                    ropts.solve_reference_u = !trust_gauge || solve_carried_u;
                    fprintf(stderr,
                            "  quadribbon U scaffold: input VMESH U=[%.3f,%.3f] "
                            "orients/gates Stage-C correspondences\n",
                            ulo, uhi);
                    {
                        char confidence_path[4096];
                        if (stem_suffix_path(input_path, "_confidence.f32",
                                             confidence_path,
                                             sizeof confidence_path) == 0) {
                            input_reference_u_confidence = (float *)read_raw_exact(
                                arena, confidence_path, sizeof(float), nv);
                        }
                        if (input_reference_u_confidence != NULL) {
                            size_t armed = 0;
                            for (size_t vi = 0; vi < nv; vi++) {
                                double q = input_reference_u_confidence[vi];
                                if (isfinite(q) && q >= 0.75 && q <= 1.0)
                                    armed++;
                            }
                            ropts.reference_u_confidence =
                                input_reference_u_confidence;
                            fprintf(stderr,
                                    "  quadribbon winding posterior: %s "
                                    "(%zu/%zu samples arm hard order)\n",
                                    confidence_path, armed, nv);
                        }
                    }
                }
                if (trust_gauge) {
                    char cert_path[4096];
                    int32_t material_hi = -1, relation_hi = -1, mesh_hi = -1;
                    if (stem_suffix_path(input_path, "_index.f32", cert_path,
                                         sizeof cert_path) == 0)
                        certificate_winding = (float *)read_raw_exact(
                            arena, cert_path, sizeof(float), nv);
                    if (stem_suffix_path(input_path,
                                         "_material_identity.i32", cert_path,
                                         sizeof cert_path) == 0)
                        certificate_material = (int32_t *)read_raw_exact(
                            arena, cert_path, sizeof(int32_t), nv);
                    if (stem_suffix_path(input_path, "_relation_island.i32",
                                         cert_path, sizeof cert_path) == 0)
                        certificate_relation = (int32_t *)read_raw_exact(
                            arena, cert_path, sizeof(int32_t), nv);
                    if (stem_suffix_path(input_path, "_mesh_component.i32",
                                         cert_path, sizeof cert_path) == 0)
                        certificate_mesh_component = (int32_t *)read_raw_exact(
                            arena, cert_path, sizeof(int32_t), nv);
                    if (nv <= SIZE_MAX / 3 &&
                        stem_suffix_path(input_path, "_owner_block.i32",
                                         cert_path, sizeof cert_path) == 0)
                        certificate_owner_block = (int32_t *)read_raw_exact(
                            arena, cert_path, sizeof(int32_t), nv * 3);
                    if (certificate_winding == NULL ||
                        input_reference_u_confidence == NULL ||
                        certificate_material == NULL ||
                        certificate_relation == NULL ||
                        certificate_mesh_component == NULL) {
                        fprintf(stderr,
                                "ERROR: --trust-gauge requires the complete "
                                "winding certificate sidecars beside %s\n",
                                input_path);
                        RAISE(IO_Failed);
                    }
                    certificate_phi = (float *)ARENA_ALLOC(
                        arena, nv * sizeof(*certificate_phi));
                    for (size_t vi = 0; vi < nv; vi++) {
                        double q = (double)certificate_winding[vi];
                        double phi = (uopts.winding_sense < 0 ? -1.0 : 1.0) *
                                     q * 6.283185307179586476925286766559;
                        if (!isfinite(q) || !isfinite(phi) ||
                            certificate_material[vi] < 0 ||
                            certificate_relation[vi] < 0 ||
                            certificate_mesh_component[vi] < 0) {
                            fprintf(stderr,
                                    "ERROR: invalid trusted certificate at "
                                    "vertex %zu\n", vi);
                            RAISE(IO_Failed);
                        }
                        certificate_phi[vi] = (float)phi;
                        if (certificate_material[vi] > material_hi)
                            material_hi = certificate_material[vi];
                        if (certificate_relation[vi] > relation_hi)
                            relation_hi = certificate_relation[vi];
                        if (certificate_mesh_component[vi] > mesh_hi)
                            mesh_hi = certificate_mesh_component[vi];
                        if (certificate_owner_block != NULL) {
                            for (int d = 0; d < 3; d++) {
                                if (certificate_owner_block[vi * 3 +
                                                            (size_t)d] ==
                                    INT32_MIN) {
                                    fprintf(stderr,
                                            "ERROR: invalid canonical owner "
                                            "coordinate at vertex %zu\n", vi);
                                    RAISE(IO_Failed);
                                }
                            }
                        }
                    }
                    certificate_material_count = (size_t)material_hi + 1;
                    certificate_relation_count = (size_t)relation_hi + 1;
                    certificate_mesh_component_count = (size_t)mesh_hi + 1;
                    /* The canonical block merger gives every physical source
                     * component a stable global id.  Do not throw that away
                     * and renumber charts from crop-local face order: claimant
                     * ties would then select a different sheet merely because
                     * lower-coordinate cubes were prepended. */
                    ropts.reference_chart = certificate_mesh_component;
                    ropts.reference_chart_count =
                        certificate_mesh_component_count;
                    ropts.reference_owner_block = certificate_owner_block;
                    fprintf(stderr,
                            "  trusted certificate: winding/material/relation/"
                            "component sidecars authoritative (%zu/%zu/%zu "
                            "labels); canonical ownership %s\n",
                            certificate_material_count,
                            certificate_relation_count,
                            certificate_mesh_component_count,
                            certificate_owner_block != NULL
                                ? "armed" : "unavailable (legacy certificate)");
                }
            } else if (ObjIO_read(
                           arena, input_path, &verts, &nv, &faces, &nf) != 0) {
                fprintf(stderr, "ERROR: cannot read %s\n", input_path);
                RAISE(IO_Failed);
            }
        }
        t_load = ves_clock_sec() - ta;
        fprintf(stderr, "  loaded %zu verts, %zu faces (%.2fs)\n", nv, nf, t_load);
        if (boundary_winding_path != NULL) {
            float *bw = (float *)read_raw_exact(
                arena, boundary_winding_path, sizeof(float), nv);
            float *bu = (float *)read_raw_exact(
                arena, boundary_u_path, sizeof(float), nv);
            float *bv = (float *)read_raw_exact(
                arena, boundary_v_path, sizeof(float), nv);
            int32_t *bm = (int32_t *)read_raw_exact(
                arena, boundary_material_path, sizeof(int32_t), nv);
            size_t supported = 0;
            if (bw == NULL || bu == NULL || bv == NULL || bm == NULL) {
                fprintf(stderr, "ERROR: cannot read exact projective boundary "
                        "arrays for %zu vertices\n", nv);
                RAISE(IO_Failed);
            }
            for (size_t vi = 0; vi < nv; vi++) {
                int valid = isfinite((double)bw[vi]);
                if ((valid && (!isfinite((double)bu[vi]) ||
                               !isfinite((double)bv[vi]) || bm[vi] < 0)) ||
                    (!valid && (isfinite((double)bu[vi]) ||
                                isfinite((double)bv[vi]) || bm[vi] >= 0))) {
                    fprintf(stderr, "ERROR: inconsistent projective boundary "
                            "validity at vertex %zu\n", vi);
                    RAISE(IO_Failed);
                }
                supported += valid;
            }
            if (supported == 0) {
                fprintf(stderr, "ERROR: projective boundary has no supported "
                        "vertices\n");
                RAISE(IO_Failed);
            }
            uopts.boundary_winding = bw;
            uopts.boundary_u = bu;
            uopts.boundary_v = bv;
            uopts.boundary_material = bm;
            fprintf(stderr, "  projective boundary: %zu/%zu exact parent "
                    "vertices loaded\n", supported, nv);
        }
        if (preserve_input_topology &&
            sibling_path(input_path, "ribbon_phase.f32",
                         sidecar_phase_path, sizeof sidecar_phase_path) == 0 &&
            sibling_path(input_path, "ribbon_material_identity.i32",
                         sidecar_material_path,
                         sizeof sidecar_material_path) == 0) {
            sidecar_phase = (float *)read_raw_exact(
                arena, sidecar_phase_path, sizeof(*sidecar_phase), nv);
            sidecar_material = (int32_t *)read_raw_exact(
                arena, sidecar_material_path, sizeof(*sidecar_material), nv);
            if (sidecar_phase == NULL || sidecar_material == NULL) {
                sidecar_phase = NULL;
                sidecar_material = NULL;
            }
        }
        /* --- Sever fusion handles (genus reduction). A residual inter-wrap
         * fusion bridge is a genus handle; the winding BFS shortcuts through
         * it and collapses whole turns (distinct wraps land on the same u).
         * Cut every non-separating loop so the sheet is genus 0 first. --- */
        long handles_cut = 0;
        if (sever_handles) {
            ta = ves_clock_sec();
            float *sv = NULL; int32_t *sf = NULL;
            size_t snv = 0, snf = 0;
            if (SeamCut_sever_short_handles(arena, verts, nv, faces, nf,
                                            0.0 /* all lengths */,
                                            &sv, &snv, &sf, &snf,
                                            &handles_cut) == 0) {
                fprintf(stderr, "  sever: %ld handle loop(s) cut, verts %zu -> %zu (%.2fs)\n",
                        handles_cut, nv, snv, ves_clock_sec() - ta);
                verts = sv; nv = snv; faces = sf; nf = snf;
            } else {
                fprintf(stderr, "  WARN: handle severing failed; continuing on input mesh\n");
            }
        }
        t_sever = ves_clock_sec() - ta;

        /* --- Reference winding scaffold. The graph trace supplies absolute
         * turn placement; Ribbon retains its own slice-domain arc-length solve.
         * Runs on the FULL mesh -- there is no internal coarse proxy: remeshing
         * the crumpled core pinched thin necks apart (fragmenting a 1-component
         * sheet into hundreds), so parameterization uses the input resolution. */
        ta = ves_clock_sec();
        int urc = 0;
        if (trust_gauge && certificate_winding != NULL) {
            double qlo = INFINITY, qhi = -INFINITY;
            double ulo = INFINITY, uhi = -INFINITY;
            for (size_t vi = 0; vi < nv; vi++) {
                double q = certificate_winding[vi];
                double u = input_reference_u[vi];
                if (q < qlo) qlo = q;
                if (q > qhi) qhi = q;
                if (u < ulo) ulo = u;
                if (u > uhi) uhi = u;
            }
            ures.phi = certificate_phi;
            ures.winding_index = certificate_winding;
            ures.winding_confidence = input_reference_u_confidence;
            ures.continuation_island = certificate_material;
            ures.island = certificate_relation;
            ures.mesh_component = certificate_mesh_component;
            ures.n_components = certificate_mesh_component_count <= INT_MAX
                              ? (int)certificate_mesh_component_count : INT_MAX;
            ures.winding_continuation_components = certificate_material_count;
            ures.winding_relation_components = certificate_relation_count;
            ures.winding_sense = uopts.winding_sense;
            ures.turns = qhi - qlo;
            ures.u_span = uhi - ulo;
            ures.continuation_satisfaction = 1.0;
            ures.order_satisfaction = 1.0;
            fprintf(stderr,
                    "  unwrap ref: BYPASSED crop-wide solve; consuming "
                    "immutable trusted certificate (%.3f turns)\n",
                    ures.turns);
        } else {
            uopts.vertex_cube = no_cube_table ? NULL :
                load_vertex_cube_table(arena, input_path, nv);
            urc = Unwrap_run(arena, verts, nv, faces, nf, &uopts, &ures);
        }
        t_unwrap = ves_clock_sec() - ta;
        if (urc == 0) {
            ropts.reference_phi = ures.phi;
            {
                char winding_path[4096];
                snprintf(winding_path, sizeof winding_path,
                         "%s/%s_winding_index.f32", out_dir, id);
                int sidecar_rc = write_raw_exact(
                    winding_path, ures.winding_index, sizeof(float), nv);
                snprintf(winding_path, sizeof winding_path,
                         "%s/%s_winding_confidence.f32", out_dir, id);
                sidecar_rc |= write_raw_exact(
                    winding_path, ures.winding_confidence, sizeof(float), nv);
                /* Streaming-block registration needs the certificate gauge
                 * per physical sheet, not one vote over every sheet in a
                 * shared cube.  continuation_island is exactly the material
                 * identity consumed by Ribbon below, so persist it beside the
                 * topology-preserving winding VMESH.  The other two labels
                 * are diagnostics: relation-island is an observability class,
                 * while mesh-component records the input topology. */
                snprintf(winding_path, sizeof winding_path,
                         "%s/%s_winding_material_identity.i32", out_dir, id);
                sidecar_rc |= write_raw_exact(
                    winding_path, ures.continuation_island,
                    sizeof(int32_t), nv);
                snprintf(winding_path, sizeof winding_path,
                         "%s/%s_winding_relation_island.i32", out_dir, id);
                sidecar_rc |= write_raw_exact(
                    winding_path, ures.island, sizeof(int32_t), nv);
                snprintf(winding_path, sizeof winding_path,
                         "%s/%s_winding_mesh_component.i32", out_dir, id);
                sidecar_rc |= write_raw_exact(
                    winding_path, ures.mesh_component,
                    sizeof(int32_t), nv);
                if (certificate_owner_block != NULL) {
                    snprintf(winding_path, sizeof winding_path,
                             "%s/%s_winding_owner_block.i32", out_dir, id);
                    sidecar_rc |= write_raw_exact(
                        winding_path, certificate_owner_block,
                        sizeof(int32_t), nv * 3);
                }
                if (ures.winding_field_used && ures.field_winding != NULL &&
                    ures.field_jump != NULL) {
                    snprintf(winding_path, sizeof winding_path,
                             "%s/%s_winding_field.f32", out_dir, id);
                    sidecar_rc |= write_raw_exact(
                        winding_path, ures.field_winding, sizeof(float), nv);
                    snprintf(winding_path, sizeof winding_path,
                             "%s/%s_winding_jump.f32", out_dir, id);
                    sidecar_rc |= write_raw_exact(
                        winding_path, ures.field_jump, sizeof(float), nv);
                }
                if (sidecar_rc != 0) {
                    fprintf(stderr,
                            "  ERROR: winding-index sidecar emission failed\n");
                    RAISE(IO_Failed);
                }
                fprintf(stderr,
                    "  winding index: wrote %s/%s_winding_{index,confidence}.f32"
                    " + material/relation/mesh identity i32%s%s\n",
                    out_dir, id,
                    certificate_owner_block != NULL ? " + owner-block i32x3" : "",
                    ures.winding_field_used ? " (+field,jump)" : "");
            }
            /* Positive same-sheet continuation is material identity, not an
             * atlas layout.  Use it only to prevent adjacent physical branches
             * from entering one StrokeStrip run; every input face still emits. */
            ropts.reference_material_island = ures.continuation_island;
            ropts.reference_material_island_count =
                ures.winding_continuation_components;
            fprintf(stderr,
                "  unwrap ref: comps=%d material_candidates=%zu gauge_islands=%zu "
                "turns=%.3f r_ref=%.1f spiral b=%.3f r2=%.3f (%.2fs)\n",
                ures.n_components, ures.winding_continuation_components,
                ures.winding_relation_components, ures.turns, ures.r_ref,
                ures.spiral_b, ures.spiral_r2, t_unwrap);
            fprintf(stderr,
                "              field=%s clean=%.1f%% supported=%zu/%d; "
                "MRF rounds=%zu changes=%zu abstain=%zu conf=%.3f; "
                "exclude=%zu rounds %zu->%zu bins (%s)\n",
                ures.winding_field_used ? "on" : "off",
                100.0 * ures.winding_field_clean_fraction,
                ures.winding_field_supported_components, ures.n_components,
                ures.winding_mrf_rounds, ures.winding_mrf_label_changes,
                ures.winding_mrf_abstained_sites,
                ures.winding_mrf_mean_confidence,
                ures.winding_mrf_conflict_rounds,
                ures.winding_mrf_conflict_bins_before,
                ures.winding_mrf_conflict_bins_after,
                ures.winding_mrf_conflict_converged
                    ? "converged" : "fixed-point residual");
            if (ures.n_components != 1)
                fprintf(stderr, "  note: %d components -- unwrapping ALL sheets in one "
                        "global (u,v) frame\n",
                        ures.n_components);
            /* The winding field supplies ordering; the continuation label is a
             * fail-closed candidate-identity gate.  Neither supplies metric U. */
            if (sidecar_phase != NULL && sidecar_material != NULL) {
                double *offset_sum = (double *)ARENA_CALLOC(
                    arena, (size_t)ures.n_components, sizeof(*offset_sum));
                size_t *offset_count = (size_t *)ARENA_CALLOC(
                    arena, (size_t)ures.n_components, sizeof(*offset_count));
                float *repaired = (float *)ARENA_ALLOC(
                    arena, nv * sizeof(*repaired));
                size_t repaired_count = 0;
                size_t ungauged_repair_count = 0;
                int32_t max_material = -1;
                int sidecar_ok = ures.mesh_component != NULL &&
                                 ures.phi != NULL && ures.n_components > 0;
                double phase_lo = 1e300, phase_hi = -1e300;

                for (size_t i = 0; sidecar_ok && i < nv; i++) {
                    int32_t component = ures.mesh_component[i];
                    int32_t material = sidecar_material[i];
                    if (component < 0 || component >= ures.n_components ||
                        material < 0) {
                        sidecar_ok = 0;
                        break;
                    }
                    if (material > max_material) max_material = material;
                    if (isfinite((double)sidecar_phase[i]) &&
                        isfinite((double)ures.phi[i])) {
                        offset_sum[component] +=
                            (double)sidecar_phase[i] - (double)ures.phi[i];
                        offset_count[component]++;
                    }
                }
                for (size_t i = 0; sidecar_ok && i < nv; i++) {
                    double phase;
                    int32_t component = ures.mesh_component[i];
                    if (isfinite((double)sidecar_phase[i])) {
                        phase = (double)sidecar_phase[i];
                    } else if (offset_count[component] != 0 &&
                               isfinite((double)ures.phi[i])) {
                        phase = (double)ures.phi[i] +
                            offset_sum[component] /
                            (double)offset_count[component];
                        repaired_count++;
                    } else if (isfinite((double)ures.phi[i])) {
                        /* A tiny generated-only connected component can have
                         * no finite value in an old sidecar.  Keep the finite
                         * graph lift for that component; never propagate the
                         * sidecar's non-finite validity sentinel. */
                        phase = (double)ures.phi[i];
                        repaired_count++;
                        ungauged_repair_count++;
                    } else {
                        sidecar_ok = 0;
                        break;
                    }
                    if (!isfinite(phase)) {
                        sidecar_ok = 0;
                        break;
                    }
                    repaired[i] = (float)phase;
                    if (phase < phase_lo) phase_lo = phase;
                    if (phase > phase_hi) phase_hi = phase;
                }
                if (sidecar_ok && max_material < INT32_MAX) {
                    ropts.reference_phi = repaired;
                    ropts.reference_phi_authoritative = 1;
                    ropts.reference_material_island = sidecar_material;
                    ropts.reference_material_island_count =
                        (size_t)max_material + 1;
                    fprintf(stderr,
                        "  quadribbon scaffold: %s + %s; turns=%.3f, "
                        "material=%d, repaired_nonfinite=%zu "
                        "(ungauged=%zu)\n",
                        sidecar_phase_path, sidecar_material_path,
                        (phase_hi - phase_lo) /
                            6.283185307179586476925286766559,
                        max_material + 1, repaired_count,
                        ungauged_repair_count);
                } else {
                    fprintf(stderr,
                        "  WARN: quadribbon sidecars failed finite/label "
                        "validation; using recomputed unwrap scaffold\n");
                }
            }
        } else {
            ropts.reference_phi = NULL;
            ropts.reference_material_island = NULL;
            ropts.reference_material_island_count = 0;
            ropts.reference_island = NULL;
            ropts.reference_island_count = 0;
        }

        if (winding_only) {
            if (urc != 0 || write_winding_ready_vmesh(
                      arena, input_path, out_dir, id, verts, nv, faces, nf,
                      &ures, &uopts, t_unwrap,
                     uopts.winding_conflict_mode >= 0) != 0) {
                fprintf(stderr,
                        "ERROR: winding-only certificate emission failed\n");
                RAISE(IO_Failed);
            }
            ok = 1;
        } else {
        /* --- Ribbon parameterization (on the full mesh). --- */
        ta = ves_clock_sec();
        int rrc = Ribbon_run(arena, verts, nv, faces, nf, &ropts, &rres);
        t_ribbon = ves_clock_sec() - ta;
        if (rrc != 0) {
            fprintf(stderr, "ERROR: ribbon parameterization failed\n");
            RAISE(IO_Failed);
        }
        if (preserve_input_topology)
            fprintf(stderr,
                    "  fixed-topology contract: V is axial; U is the staged "
                    "StrokeStrip solution on input rows%s\n",
                    ropts.metric_project_only
                        ? " (legacy exact-row projection diagnostic)"
                        : " (final C(u) frozen during metric continuation)");
        fprintf(stderr,
            "  ribbon: slices=%d chains=%d (closed=%d, fragmented slices=%d)\n"
            "          samples=%zu pairs=%zu(+%zu cont) cover=%.1f%% bridge_cuts=%zu\n"
            "          winding: %d groups (radial-placed=%d unreached=%d) pitch=%.2f "
            "conflicts=%zu; spiral b=%.2f r2=%.3f\n"
            "          solve: %d components, registration max shift %.1f vox\n"
            "          StrokeStrip: runs=%zu members=%zu links=%zu pruned=%zu "
            "length_rms=%.5f align_rms=%.5f\n"
            "          mono_repairs=%zu\n"
            "          u_span=%.1f v_span=%.1f turns=%.2f\n"
            "          |du/ds-1|: mean=%.4f max=%.3f hist[<1 <2 <5 <10 <20 >=20 %%]="
            "%ld %ld %ld %ld %ld %ld  (%.2fs)\n",
            rres.n_slices, rres.n_chains, rres.n_closed, rres.n_multi_slices,
            rres.n_samples, rres.n_pairs, rres.n_cont_pairs,
            100.0 * rres.match_cover, rres.bridge_cuts,
            rres.w_groups, rres.w_prior_groups, rres.w_unreached, rres.pitch_used,
            rres.w_conflicts, rres.spiral_b, rres.spiral_r2,
            rres.n_qp_comps, rres.reg_max_shift,
            rres.n_strip_runs, rres.n_strip_members, rres.n_strip_links,
            rres.n_strip_links_pruned, rres.strip_length_rms,
            rres.strip_align_rms,
            rres.mono_repairs,
            rres.u_span, rres.v_span, rres.phi_span_turns,
            rres.duds_err_mean, rres.duds_err_max,
            rres.duds_hist[0], rres.duds_hist[1], rres.duds_hist[2],
            rres.duds_hist[3], rres.duds_hist[4], rres.duds_hist[5], t_ribbon);

        /* --- Emit the parameterized ribbon. --- */
        ta = ves_clock_sec();
        char path[4096];

        size_t rib_nv = 0, rib_nf = 0;
        size_t rib_atlas_cols = 0, rib_atlas_runs = 0;
        size_t rib_empty_cols_removed = 0;
        RibbonWriteStats wstats;
        memset(&wstats, 0, sizeof wstats);
        if (preserve_input_topology) {
            size_t finite_uv = 0;
            if (rres.uv == NULL) {
                fprintf(stderr,
                        "  ERROR: topology-preserving parameterization "
                        "produced no vertex UV\n");
                RAISE(IO_Failed);
            }
            for (size_t i = 0; i < nv; i++)
                if (isfinite((double)rres.uv[i*2]) &&
                    isfinite((double)rres.uv[i*2+1]))
                    finite_uv++;
            if (finite_uv != nv) {
                fprintf(stderr,
                        "  ERROR: topology-preserving parameterization has "
                        "%zu/%zu finite UV vertices\n", finite_uv, nv);
                RAISE(IO_Failed);
            }

            snprintf(path, sizeof path, "%s/%s_ribbon.obj", out_dir, id);
            ves_ensure_parent_dir(path);
            char vmesh_path[4096];
            if (MeshBin_companion_path(
                    path, vmesh_path, sizeof vmesh_path) != 0) {
                fprintf(stderr,
                        "  ERROR: cannot derive authoritative VMESH path from %s\n",
                        path);
                RAISE(IO_Failed);
            }
            if (write_text_obj &&
                ObjIO_write_uv(path, verts, nv, faces, nf, rres.uv) != 0) {
                fprintf(stderr, "  ERROR: UV OBJ emission failed: %s\n", path);
                RAISE(IO_Failed);
            }
            if (MeshBin_write(vmesh_path, verts, nv, faces, nf, rres.uv) != 0) {
                fprintf(stderr,
                        "  ERROR: authoritative VMESH emission failed: %s\n",
                        vmesh_path);
                RAISE(IO_Failed);
            }
            rib_nv = nv;
            rib_nf = nf;
            if (write_text_obj)
                fprintf(stderr,
                        "  wrote %s and %s (%zu verts, %zu faces; topology "
                        "identical to input)\n",
                        path, vmesh_path, rib_nv, rib_nf);
            else
                fprintf(stderr,
                        "  wrote %s (%zu verts, %zu faces; topology "
                        "identical to input; text OBJ skipped)\n",
                        vmesh_path, rib_nv, rib_nf);
            fprintf(stderr,
                    "          UV transfer: %zu direct, %zu neighbor-filled, "
                    "%zu unmapped\n",
                    nv - rres.uv_filled - rres.uv_fallback,
                    rres.uv_filled, rres.uv_fallback);
        } else {
            snprintf(path, sizeof path, "%s/%s_ribbon.obj", out_dir, id);
            ves_ensure_parent_dir(path);
            if (Ribbon_write_obj_promote_peels(
                    path, &rres, write_text_obj, 0,
                    promote_peel_min_share, promote_peel_min_vertices,
                    &wstats) != 0) {
                fprintf(stderr,"  ERROR: fitted ribbon emission failed: %s\n",path);
                RAISE(IO_Failed);
            }
            rib_nv = wstats.vertices;
            rib_nf = wstats.faces;
            rib_atlas_cols = wstats.atlas_columns;
            rib_atlas_runs = wstats.atlas_runs;
            rib_empty_cols_removed = wstats.empty_columns_removed;
            if (rib_nf == 0 && rres.grid_projective) {
                fprintf(stderr, "  wrote empty/point-only layer inventory; no surface file\n");
            } else if (write_text_obj) {
                fprintf(stderr,
                        "  wrote %s (%zu verts, %zu faces, direct grid %zux%zu)\n",
                        path, rib_nv, rib_nf, rres.nu, rres.nk);
            } else {
                char vmesh_path[4096];
                if (MeshBin_companion_path(
                        path, vmesh_path, sizeof vmesh_path) != 0) {
                    fprintf(stderr,
                            "  ERROR: cannot derive fitted VMESH path from %s\n",
                            path);
                    RAISE(IO_Failed);
                }
                fprintf(stderr,
                        "  wrote %s (%zu verts, %zu faces, direct grid %zux%zu; "
                        "text OBJ skipped)\n",
                        vmesh_path, rib_nv, rib_nf, rres.nu, rres.nk);
            }
            fprintf(stderr,
                    "          confidence: %zu fixed support, %zu generated continuation\n",
                    wstats.supported_vertices, wstats.generated_vertices);
            fprintf(stderr,
                    "          emitted atlas: %zu -> %zu columns in %zu "
                    "occupied run(s), %zu empty column(s) removed\n",
                    rres.nu, rib_atlas_cols, rib_atlas_runs,
                    rib_empty_cols_removed);
            if (wstats.promoted_peel_runs > 0 ||
                promote_peel_min_share > 0.0)
                fprintf(stderr,
                        "          peel partition: %zu run(s)/%zu vertices "
                        "promoted, %zu run(s)/%zu vertices retained as extras "
                        "(threshold %zu of %zu primary-layer vertices)\n",
                        wstats.promoted_peel_runs,
                        wstats.promoted_peel_vertices,
                        wstats.extra_peel_runs,
                        wstats.extra_peel_vertices,
                        wstats.peel_promotion_threshold,
                        wstats.primary_peel_vertices);
            if (ropts.component_global) {
                fprintf(stderr, "          reconstruction components: %zu "
                        "(%zu far path conflict(s), %zu branch split(s), "
                        "%zu relation cut(s), support %.1f)\n",
                        rres.grid_reconstruction_components,
                        rres.grid_branch_conflicts,
                        rres.grid_branch_splits,
                        rres.grid_branch_relation_cuts,
                        rres.grid_branch_relation_cut_support);
                fprintf(stderr, "          component-global grid claims: "
                        "%zu extra claim(s) in %zu cell(s), %zu discarded, "
                        "%zu later candidate(s) selected\n",
                        rres.grid_claim_conflicts,
                        rres.grid_claim_conflict_cells,
                        rres.grid_claim_discarded_cells,
                        rres.grid_claim_replaced);
            }
        }
        t_write = ves_clock_sec() - ta;

        /* --- Stats JSON. --- */
        double t_total = ves_clock_sec() - t0;
        snprintf(path, sizeof path, "%s/%s_ribbon_stats.json", out_dir, id);
        FILE *js = fopen(path, "w");
        if (js != NULL) {
            char json_id[512], json_input[4096];
            json_escape_copy(json_id, sizeof json_id, id);
            json_escape_copy(json_input, sizeof json_input, input_path);
            fprintf(js,
                "{\n"
                "  \"id\": \"%s\",\n"
                "  \"input\": \"%s\",\n"
                "  \"n_verts\": %zu,\n"
                "  \"n_faces\": %zu,\n"
                "  \"axis_point_zyx\": [%.4f, %.4f, %.4f],\n"
                "  \"axis_dir_zyx\": [%.6f, %.6f, %.6f],\n"
                "  \"opts\": { \"slice_h\": %.3f, \"sample_h\": %.3f, \"match_r\": %.3f,\n"
                "            \"match_ang_deg\": %.1f, \"relax_iters\": %d, \"final_iters\": %d,\n"
                "            \"metric_iters\": %d, \"metric_weight\": %.6g,\n"
                "            \"solve_carried_u\": %s,\n"
                "            \"coarse_chain_seed\": %s,\n"
                "            \"grid_u\": %.3f, \"peel_layers_requested\": %d, \"peel_layer_selected\": %d, "
                "\"component_global\": %s,\n"
                "            \"ownership\": \"%s\", \"discard_conflicting_claims\": %s,\n"
                "            \"coherent_claims\": %s, \"ufill_rejected_crossings\": %zu,\n"
                "            \"wrap_spacing_requested\": %.4f,\n"
                "            \"promote_peel_min_share\": %.9g,\n"
                "            \"promote_peel_min_vertices\": %zu },\n"
                "  \"winding\": { \"n_components\": %d, \"turns\": %.4f, \"r_ref\": %.3f,\n"
                "               \"spiral_b\": %.4f, \"spiral_r2\": %.4f, \"u_span_cyl\": %.2f },\n"
                "  \"winding_registration\": { \"sense\": %d, \"sense_pinned\": %d,\n"
                "              \"sense_spiral\": %d, \"sense_ray_vote\": %d,\n"
                "              \"sense_ray_vote_agree\": %.3f, \"sense_ray_vote_pairs\": %zu,\n"
                "              \"bins\": %zu, \"strands\": %zu,\n"
                "    \"continuation_observations\": %zu, \"order_observations\": %zu, \"order_observations_suppressed\": %zu,\n"
                "    \"relations\": %zu, \"eligible_relations\": %zu, \"forest_relations\": %zu, \"order_relations_suppressed\": %zu,\n"
                "    \"continuation_components\": %zu, \"relation_components\": %zu, \"packed_relation_components\": %zu, \"packed_mesh_components\": %zu,\n"
                "    \"relation_conflicts\": %zu, \"observations_dropped\": %zu,\n"
                "    \"continuation_satisfaction\": %.9g, \"order_satisfaction\": %.9g,\n"
                "    \"turn_correction_min\": %d, \"turn_correction_max\": %d,\n"
                "    \"repair_closers\": %zu, \"repair_conflicts_pre\": %zu, \"repair_shifts\": %zu, \"repair_capped_roots\": %zu,\n"
                "    \"anchor_span_pre_turns\": %.6g, \"anchor_span_turns\": %.6g,\n"
                "    \"mrf_rounds\": %zu, \"mrf_label_changes\": %zu, \"mrf_abstained_sites\": %zu,\n"
                "    \"mrf_energy_before\": %.9g, \"mrf_energy_after\": %.9g, \"mrf_mean_confidence\": %.9g,\n"
                "    \"conflict_resolution\": { \"rounds\": %zu, \"bins_before\": %zu, \"bins_after\": %zu,\n"
                "      \"losing_claims\": %zu, \"forbidden_assignments\": %zu, \"winner_locks\": %zu, \"label_changes\": %zu, \"converged\": %s } },\n"
                "  \"winding_field\": { \"used\": %s, \"backend\": %d, \"samples\": %zu,\n"
                "    \"supported_components\": %zu, \"clean_vertices\": %zu, \"invalid_vertices\": %zu,\n"
                "    \"clean_fraction\": %.9g },\n"
                "  \"handles_cut\": %ld,\n"
                "  \"slicing\": { \"n_slices\": %d, \"n_chains\": %d, \"n_closed\": %d,\n"
                "               \"fragmented_slices\": %d, \"n_samples\": %zu, \"bridge_cuts\": %zu },\n"
                "  \"winding_index\": { \"groups\": %d, \"radial_placed\": %d, \"unreached\": %d,\n"
                "                     \"pitch_used\": %.4f, \"pitch_source\": \"%s\",\n"
                "                     \"conflicts\": %zu,\n"
                "                     \"solve_components\": %d, \"reg_max_shift\": %.3f,\n"
                "                     \"spiral_a\": %.4f, \"spiral_b\": %.4f, \"spiral_r2\": %.4f },\n"
                "  \"pairs\": { \"cross_slice\": %zu, \"continuation\": %zu, \"match_cover\": %.4f },\n"
                "  \"strokestrip\": { \"runs\": %zu, \"members\": %zu, \"links\": %zu,\n"
                "                   \"links_pruned\": %zu, \"length_rms\": %.9g, \"align_rms\": %.9g,\n"
                "                   \"gauge_components\": %zu, \"gauge_span\": %.9g, \"gauge_max_shift\": %.9g },\n"
                 "  \"solve\": { \"mono_repairs\": %zu, \"duds_err_mean\": %.6f, \"duds_err_max\": %.6f,\n"
                 "             \"duds_hist_pct\": [%ld, %ld, %ld, %ld, %ld, %ld] },\n"
                "  \"fixed_topology\": { \"enabled\": %s, "
                "\"post_layout\": false, \"v_coordinate\": \"axial\" },\n"
                "  \"uv\": { \"u_span\": %.3f, \"v_span\": %.3f, \"turns\": %.4f,\n"
                "          \"filled_verts\": %zu, \"fallback_verts\": %zu,\n"
                "          \"phase_rejected_candidates\": %zu, \"phase_blocked_queries\": %zu,\n"
                "          \"component_rejected_candidates\": %zu, \"component_blocked_queries\": %zu,\n"
                "          \"gauge_observations\": %zu, \"gauge_relations\": %zu, \"gauge_max_shift\": %.6f,\n"
                "          \"atlas_islands\": %zu, \"atlas_packed_islands\": %zu,\n"
                "          \"atlas_gutter\": %.6f, \"atlas_pack_saved\": %.6f },\n"
                "  \"ribbon\": { \"nu\": %zu, \"nk\": %zu, \"grid_layers\": %zu, "
                "\"du\": %.3f, \"dv\": %.3f,\n"
                "              \"obj_verts\": %zu, \"obj_faces\": %zu,\n"
                "              \"atlas_cols\": %zu, \"atlas_column_runs\": %zu, "
                "\"empty_columns_removed\": %zu,\n"
                "              \"primary_peel_vertices\": %zu, "
                "\"peel_promotion_threshold\": %zu,\n"
                "              \"promoted_peel_runs\": %zu, "
                "\"promoted_peel_vertices\": %zu,\n"
                "              \"extra_peel_runs\": %zu, "
                "\"extra_peel_vertices\": %zu,\n"
                "              \"long_edges_rowfit\": %zu, \"claim_repairs\": %zu, "
                "\"long_edges_repaired\": %zu, \"long_edges_vfill\": %zu,\n"
                "              \"local_outlier_repairs\": %zu, "
                "\"local_pair_repairs\": %zu, \"local_pair_slots\": %zu, "
                "\"local_run_repairs\": %zu, \"local_run_slots\": %zu, "
                "\"local_patch_repairs\": %zu, \"local_patch_slots\": %zu, "
                "\"local_patch_small_components\": %zu, "
                "\"local_patch_enclosed_components\": %zu, "
                "\"local_patch_claim_slots\": %zu, "
                "\"local_supported_replacements\": %zu, "
                "\"long_edges_local\": %zu,\n"
                "              \"long_edges_smooth\": %zu, \"both_diagonals_long\": %zu,\n"
                "              \"claim_conflicts\": %zu, \"claim_conflict_cells\": %zu,\n"
                "              \"claim_discarded_cells\": %zu, \"claim_replaced\": %zu,\n"
                "              \"claim_candidates\": %zu, \"claim_required_peels\": %zu, "
                "\"claim_overflow\": %zu, \"claim_stored\": %zu,\n"
                "              \"subcell_claim_chains\": %zu,\n"
                "              \"reconstruction_components\": %zu, "
                "\"branch_conflicts\": %zu, \"branch_splits\": %zu,\n"
                "              \"branch_relation_cuts\": %zu, "
                "\"branch_relation_cut_support\": %.6f,\n"
                "              \"row_wrap_splits\": %zu,\n"
                "              \"row_wrap_same_mesh\": %zu,\n"
                "              \"row_wrap_same_solve\": %zu, \"row_wrap_same_island\": %zu },\n"
                "  \"ribbon_topology\": {\n"
                "              \"long_edges_diag_main\": %zu, "
                "\"long_edges_diag_anti\": %zu,\n"
                "              \"cell_lane_hist\": "
                "[%zu, %zu, %zu, %zu, %zu, %zu, %zu, %zu],\n"
                 "              \"vfill_runs\": %zu, \"vfill_rows\": %zu,\n"
                 "              \"vfill_reject_crossing\": %zu,\n"
                 "              \"material_faces_cut\": %zu,\n"
                 "              \"vfill_run_hist\": "
                "[%zu, %zu, %zu, %zu, %zu, %zu] },\n"
                "  \"timing_sec\": { \"load\": %.3f, \"sever\": %.3f,\n"
                "                  \"winding\": %.3f, \"ribbon\": %.3f, \"write\": %.3f,\n"
                "                  \"total\": %.3f }\n"
                "}\n",
                json_id, json_input, nv, nf,
                (double)uopts.axis_point[0], (double)uopts.axis_point[1], (double)uopts.axis_point[2],
                (double)uopts.axis_dir[0], (double)uopts.axis_dir[1], (double)uopts.axis_dir[2],
                (double)ropts.slice_h, (double)ropts.sample_h, (double)ropts.match_r,
                (double)ropts.match_ang_deg, ropts.relax_iters, ropts.final_iters,
                 ropts.metric_iters, (double)ropts.metric_weight,
                 solve_carried_u ? "true" : "false",
                 ropts.coarse_chain_seed ? "true" : "false",
                 (double)ropts.grid_u, ropts.peel_layers, ropts.peel_layer,
                 ropts.component_global ? "true" : "false",
                 ropts.ownership_construction ? "construction" : "claims",
                 ropts.discard_conflicting_claims ? "true" : "false",
                 ropts.coherent_claims ? "true" : "false",
                 rres.grid_ufill_reject_crossing,
                 (double)ropts.wrap_spacing,
                 promote_peel_min_share,
                 promote_peel_min_vertices,
                ures.n_components, ures.turns, ures.r_ref,
                ures.spiral_b, ures.spiral_r2, ures.u_span,
                ures.winding_sense, uopts.winding_sense,
                ures.winding_sense_spiral, ures.winding_sense_vote,
                ures.winding_sense_vote_agree, ures.winding_sense_vote_n,
                ures.winding_bins,
                ures.winding_strands,
                ures.continuation_observations, ures.order_observations,
                ures.order_observations_suppressed,
                ures.winding_relations, ures.winding_eligible_relations,
                ures.winding_forest_relations,
                ures.winding_order_relations_suppressed,
                ures.winding_continuation_components,
                ures.winding_relation_components,
                ures.winding_packed_relation_components,
                ures.winding_packed_mesh_components,
                ures.winding_relation_conflicts,
                ures.winding_observations_dropped,
                ures.continuation_satisfaction, ures.order_satisfaction,
                ures.turn_correction_min, ures.turn_correction_max,
                ures.winding_repair_closers, ures.winding_repair_conflicts_pre,
                ures.winding_repair_shifts, ures.winding_repair_capped_roots,
                ures.winding_anchor_span_pre_turns,
                ures.winding_anchor_span_turns,
                ures.winding_mrf_rounds,
                ures.winding_mrf_label_changes,
                ures.winding_mrf_abstained_sites,
                ures.winding_mrf_energy_before,
                ures.winding_mrf_energy_after,
                ures.winding_mrf_mean_confidence,
                ures.winding_mrf_conflict_rounds,
                ures.winding_mrf_conflict_bins_before,
                ures.winding_mrf_conflict_bins_after,
                ures.winding_mrf_conflict_losing_claims,
                ures.winding_mrf_conflict_exclusions,
                ures.winding_mrf_conflict_winner_locks,
                ures.winding_mrf_conflict_label_changes,
                ures.winding_mrf_conflict_converged ? "true" : "false",
                ures.winding_field_used ? "true" : "false",
                ures.winding_field_backend,
                ures.winding_field_samples,
                ures.winding_field_supported_components,
                ures.winding_field_clean_vertices,
                ures.winding_field_invalid_vertices,
                ures.winding_field_clean_fraction,
                handles_cut,
                rres.n_slices, rres.n_chains, rres.n_closed,
                rres.n_multi_slices, rres.n_samples, rres.bridge_cuts,
                rres.w_groups, rres.w_prior_groups, rres.w_unreached,
                rres.pitch_used,
                rres.pitch_source == RIB_PITCH_PINNED ? "pinned" :
                rres.pitch_source == RIB_PITCH_ESTIMATED ? "estimated" : "fallback",
                rres.w_conflicts,
                rres.n_qp_comps, rres.reg_max_shift,
                rres.spiral_a, rres.spiral_b, rres.spiral_r2,
                rres.n_pairs, rres.n_cont_pairs, rres.match_cover,
                rres.n_strip_runs, rres.n_strip_members, rres.n_strip_links,
                rres.n_strip_links_pruned, rres.strip_length_rms,
                rres.strip_align_rms,
                rres.strip_gauge_components, rres.strip_gauge_span,
                rres.strip_gauge_max_shift,
                rres.mono_repairs, rres.duds_err_mean, rres.duds_err_max,
                rres.duds_hist[0], rres.duds_hist[1], rres.duds_hist[2],
                rres.duds_hist[3], rres.duds_hist[4], rres.duds_hist[5],
                preserve_input_topology ? "true" : "false",
                rres.u_span, rres.v_span, rres.phi_span_turns,
                rres.uv_filled, rres.uv_fallback,
                rres.uv_phase_rejects, rres.uv_phase_blocked,
                rres.uv_component_rejects, rres.uv_component_blocked,
                rres.uv_gauge_observations, rres.uv_gauge_relations,
                rres.uv_gauge_max_shift,
                rres.uv_atlas_islands, rres.uv_atlas_packed_islands,
                 rres.uv_atlas_gutter, rres.uv_atlas_pack_saved,
                 rres.nu, rres.nk, rres.grid_layers,
                 (double)rres.grid_du, (double)rres.grid_dv,
                 rib_nv, rib_nf, rib_atlas_cols, rib_atlas_runs,
                 rib_empty_cols_removed,
                 wstats.primary_peel_vertices,
                 wstats.peel_promotion_threshold,
                 wstats.promoted_peel_runs, wstats.promoted_peel_vertices,
                 wstats.extra_peel_runs, wstats.extra_peel_vertices,
                  rres.grid_long_edges_rowfit, rres.grid_claim_repairs,
                  rres.grid_long_edges_repaired,rres.grid_long_edges_vfill,
                  rres.grid_local_outlier_repairs,
                  rres.grid_local_pair_repairs,
                  rres.grid_local_pair_slots,
                  rres.grid_local_run_repairs,
                  rres.grid_local_run_slots,
                  rres.grid_local_patch_repairs,
                  rres.grid_local_patch_slots,
                  rres.grid_local_patch_small_components,
                  rres.grid_local_patch_enclosed_components,
                  rres.grid_local_patch_claim_slots,
                  rres.grid_local_supported_replacements,
                  rres.grid_long_edges_local,
                  rres.grid_long_edges_smooth,
                  rres.grid_both_diagonals_long,
                  rres.grid_claim_conflicts,
                  rres.grid_claim_conflict_cells,
                  rres.grid_claim_discarded_cells,
                  rres.grid_claim_replaced,
                  rres.grid_claim_candidates,
                  rres.grid_claim_required_peels,
                  rres.grid_claim_overflow,
                  rres.grid_claim_stored,
                  rres.grid_subcell_claim_chains,
                  rres.grid_reconstruction_components,
                  rres.grid_branch_conflicts,
                  rres.grid_branch_splits,
                  rres.grid_branch_relation_cuts,
                  rres.grid_branch_relation_cut_support,
                 rres.grid_row_wrap_splits,
                 rres.grid_row_wrap_same_mesh,
                 rres.grid_row_wrap_same_solve,rres.grid_row_wrap_same_island,
                 rres.grid_long_edges_diag_main,
                 rres.grid_long_edges_diag_anti,
                 rres.grid_cell_lane_hist[0], rres.grid_cell_lane_hist[1],
                 rres.grid_cell_lane_hist[2], rres.grid_cell_lane_hist[3],
                 rres.grid_cell_lane_hist[4], rres.grid_cell_lane_hist[5],
                 rres.grid_cell_lane_hist[6], rres.grid_cell_lane_hist[7],
                  rres.grid_vfill_runs, rres.grid_vfill_rows,
                  rres.grid_vfill_reject_crossing,
                  rres.grid_material_faces_cut,
                 rres.grid_vfill_run_hist[0], rres.grid_vfill_run_hist[1],
                 rres.grid_vfill_run_hist[2], rres.grid_vfill_run_hist[3],
                 rres.grid_vfill_run_hist[4], rres.grid_vfill_run_hist[5],
                t_load, t_sever, t_unwrap, t_ribbon, t_write, t_total);
            fclose(js);
            fprintf(stderr, "  wrote %s\n", path);
        }
        ok = 1;
        }
    EXCEPT(IO_Failed)
        fprintf(stderr, "scroll_ribbon: I/O failure\n");
    EXCEPT(Arena_Failed)
        fprintf(stderr, "scroll_ribbon: out of memory\n");
    END_TRY;

    double total = ves_clock_sec() - t0;
    fprintf(stderr, "scroll_ribbon: %s (total %.2fs)\n", ok ? "OK" : "FAILED", total);

    Arena_dispose(&arena);
    return ok ? 0 : 1;
}
