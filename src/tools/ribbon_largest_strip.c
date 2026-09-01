/* ribbon_largest_strip.c -- select the largest connected strip of a fitted
 * quadribbon.
 *
 * POLICY (2026-09-01).  A fitted ribbon is one lattice in (u, v), but its
 * quads exist only where the fit had support, so the lattice carries several
 * disconnected strips.  Only the largest is shippable as a big sheet: the
 * rest are slivers a reader cannot follow and that no downstream stage can
 * join without inventing a relation.  Rather than arbitrate the small pieces,
 * drop them and publish the strip that is connected by construction.
 *
 * Connectivity is the MESH's own face adjacency, not a lattice-neighbour
 * rule: two cells that merely sit next to each other in (u, v) are not the
 * same sheet unless the fit actually emitted a quad joining them.  Using the
 * faces is what keeps this from silently bridging two wraps that happen to
 * land on adjacent columns.
 *
 * Components are ranked by OCCUPIED LATTICE CELLS (readable area), not by
 * vertex or face count, so a finely tessellated patch cannot outrank a coarse
 * long strip.  Each candidate's u span is reported alongside, because
 * "longest strip" is the property a reader cares about and the two can
 * disagree.
 *
 *   ribbon_largest_strip <in.vmesh> <out.vmesh> [--du F=2] [--dv F=1]
 *                        [--top N=10] [--report out.json]
 *   ribbon_largest_strip --selftest
 * Exit: 0 ok / selftest pass, 1 IO, 2 usage, 3 selftest fail.
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/union_find.h"

typedef struct {
    int32_t root;      /* an arbitrary vertex of the component */
    size_t  nv;        /* vertices */
    size_t  nf;        /* faces */
    size_t  ncell;     /* distinct occupied (col,row) lattice cells */
    double  u_lo, u_hi;
    double  v_lo, v_hi;
} StripStat;

static int rls_cmp_desc(const void *a, const void *b)
{
    const StripStat *x = (const StripStat *)a;
    const StripStat *y = (const StripStat *)b;
    double sx = 0.0, sy = 0.0;
    if (x->ncell != y->ncell) return x->ncell < y->ncell ? 1 : -1;
    /* Ties break on u span so the longer strip wins, then on root so the
     * order does not depend on the qsort implementation. */
    sx = x->u_hi - x->u_lo;
    sy = y->u_hi - y->u_lo;
    if (sx < sy) return 1;
    if (sx > sy) return -1;
    return x->root < y->root ? -1 : (x->root > y->root ? 1 : 0);
}

typedef struct { int32_t comp; int32_t col; int32_t row; } CellKey;

static int rls_cmp_cell(const void *a, const void *b)
{
    const CellKey *x = (const CellKey *)a;
    const CellKey *y = (const CellKey *)b;
    if (x->comp != y->comp) return x->comp < y->comp ? -1 : 1;
    if (x->col  != y->col)  return x->col  < y->col  ? -1 : 1;
    if (x->row  != y->row)  return x->row  < y->row  ? -1 : 1;
    return 0;
}

/* Component statistics for a parameterized mesh.  out_label receives a dense
 * 0-based component id per vertex; out_stat a table sorted largest-first.
 * Returns 0 on success, nonzero on empty or oversized input. */
static int rls_components(Arena_T arena,
                          const float *uv, size_t nv,
                          const int32_t *faces, size_t nf,
                          double du, double dv,
                          int32_t **out_label, size_t *out_ncomp,
                          StripStat **out_stat)
{
    UnionFind uf;
    int32_t *label = NULL, *dense = NULL;
    StripStat *st = NULL;
    CellKey *key = NULL;
    size_t ncomp = 0, i = 0, f = 0, g = 0;
    double dus = 0.0, dvs = 0.0;

    if (nv == 0 || nv > (size_t)INT32_MAX) return -1;

    uf = UF_new(arena, (int32_t)nv);
    for (f = 0; f < nf; f++) {
        uf_union(&uf, faces[f*3+0], faces[f*3+1]);
        uf_union(&uf, faces[f*3+1], faces[f*3+2]);
    }

    label = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
    dense = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
    for (i = 0; i < nv; i++) dense[i] = -1;
    for (i = 0; i < nv; i++) {
        int32_t r = uf_find(&uf, (int32_t)i);
        if (dense[r] < 0) dense[r] = (int32_t)ncomp++;
        label[i] = dense[r];
    }

    st = (StripStat *)ARENA_ALLOC(arena, (long)(ncomp * sizeof(StripStat)));
    for (g = 0; g < ncomp; g++) {
        st[g].root = -1; st[g].nv = 0; st[g].nf = 0; st[g].ncell = 0;
        st[g].u_lo = 1e300; st[g].u_hi = -1e300;
        st[g].v_lo = 1e300; st[g].v_hi = -1e300;
    }
    for (i = 0; i < nv; i++) {
        StripStat *s = &st[label[i]];
        double u = (double)uv[i*2+0];
        double v = (double)uv[i*2+1];
        if (s->root < 0) s->root = (int32_t)i;
        s->nv++;
        if (u < s->u_lo) s->u_lo = u;
        if (u > s->u_hi) s->u_hi = u;
        if (v < s->v_lo) s->v_lo = v;
        if (v > s->v_hi) s->v_hi = v;
    }
    for (f = 0; f < nf; f++) st[label[faces[f*3+0]]].nf++;

    /* One sort plus a linear scan dedups every component's cells at once: the
     * cells of a component are contiguous once the key is (comp, col, row). */
    key = (CellKey *)ARENA_ALLOC(arena, (long)(nv * sizeof(CellKey)));
    dus = du > 1e-9 ? du : 1.0;
    dvs = dv > 1e-9 ? dv : 1.0;
    for (i = 0; i < nv; i++) {
        key[i].comp = label[i];
        key[i].col  = (int32_t)lround((double)uv[i*2+0] / dus);
        key[i].row  = (int32_t)lround((double)uv[i*2+1] / dvs);
    }
    qsort(key, nv, sizeof(CellKey), rls_cmp_cell);
    for (i = 0; i < nv; i++)
        if (i == 0 || rls_cmp_cell(&key[i], &key[i-1]) != 0)
            st[key[i].comp].ncell++;

    qsort(st, ncomp, sizeof(StripStat), rls_cmp_desc);
    *out_label = label;
    *out_ncomp = ncomp;
    *out_stat = st;
    return 0;
}

static int rls_run(const char *in_path, const char *out_path,
                   double du, double dv, int top, const char *report)
{
    Arena_T arena = Arena_new();
    MeshBinData mesh;
    int32_t *label = NULL, *remap = NULL, *of = NULL;
    StripStat *st = NULL;
    float *ov = NULL, *ouv = NULL;
    size_t ncomp = 0, total_cell = 0, onv = 0, onf = 0, g = 0, i = 0, f = 0;
    int32_t keep = 0;
    int rc = 1;

    memset(&mesh, 0, sizeof mesh);
    if (MeshBin_read_arena(arena, in_path, &mesh) != 0) {
        fprintf(stderr, "ERROR: cannot read %s\n", in_path);
        goto done;
    }
    if (mesh.uv == NULL) {
        fprintf(stderr, "ERROR: %s carries no UV; a strip needs (u,v)\n",
                in_path);
        goto done;
    }
    fprintf(stderr, "  loaded %zu verts, %zu faces\n", mesh.nv, mesh.nf);

    if (rls_components(arena, mesh.uv, mesh.nv, mesh.faces, mesh.nf,
                       du, dv, &label, &ncomp, &st) != 0) {
        fprintf(stderr, "ERROR: component analysis failed\n");
        goto done;
    }

    for (g = 0; g < ncomp; g++) total_cell += st[g].ncell;
    fprintf(stderr, "  strips: %zu component(s), %zu occupied cell(s)\n",
            ncomp, total_cell);
    for (g = 0; g < ncomp && (int)g < top; g++)
        fprintf(stderr,
                "    #%zu cells=%zu (%.2f%%) verts=%zu faces=%zu "
                "u=[%.1f,%.1f] span=%.1f v=[%.1f,%.1f]\n",
                g, st[g].ncell,
                total_cell ? 100.0 * (double)st[g].ncell / (double)total_cell
                           : 0.0,
                st[g].nv, st[g].nf,
                st[g].u_lo, st[g].u_hi, st[g].u_hi - st[g].u_lo,
                st[g].v_lo, st[g].v_hi);

    keep = label[st[0].root];   /* st is sorted; recover the winner's label */

    remap = (int32_t *)ARENA_ALLOC(arena, (long)(mesh.nv * sizeof(int32_t)));
    for (i = 0; i < mesh.nv; i++)
        remap[i] = label[i] == keep ? (int32_t)onv++ : -1;

    ov  = (float *)ARENA_ALLOC(arena, (long)(onv * 3 * sizeof(float)));
    ouv = (float *)ARENA_ALLOC(arena, (long)(onv * 2 * sizeof(float)));
    for (i = 0; i < mesh.nv; i++) {
        size_t o = 0;
        int d = 0;
        if (remap[i] < 0) continue;
        o = (size_t)remap[i];
        for (d = 0; d < 3; d++) ov[o*3+(size_t)d] = mesh.verts[i*3+(size_t)d];
        ouv[o*2+0] = mesh.uv[i*2+0];
        ouv[o*2+1] = mesh.uv[i*2+1];
    }
    of = (int32_t *)ARENA_ALLOC(arena, (long)(mesh.nf * 3 * sizeof(int32_t)));
    for (f = 0; f < mesh.nf; f++) {
        int d = 0;
        if (label[mesh.faces[f*3+0]] != keep) continue;
        for (d = 0; d < 3; d++)
            of[onf*3+(size_t)d] = remap[mesh.faces[f*3+(size_t)d]];
        onf++;
    }

    fprintf(stderr,
            "  kept strip #0: %zu/%zu verts (%.2f%%), %zu/%zu faces, "
            "%zu cells, u span %.1f\n",
            onv, mesh.nv, mesh.nv ? 100.0*(double)onv/(double)mesh.nv : 0.0,
            onf, mesh.nf, st[0].ncell, st[0].u_hi - st[0].u_lo);

    if (MeshBin_write(out_path, ov, onv, of, onf, ouv) != 0) {
        fprintf(stderr, "ERROR: cannot write %s\n", out_path);
        goto done;
    }
    fprintf(stderr, "  wrote %s\n", out_path);

    if (report != NULL) {
        FILE *fp = fopen(report, "wb");
        if (fp == NULL) {
            fprintf(stderr, "ERROR: cannot write %s\n", report);
            goto done;
        }
        fprintf(fp, "{\n  \"input\": \"%s\",\n", in_path);
        fprintf(fp, "  \"n_components\": %zu,\n", ncomp);
        fprintf(fp, "  \"total_cells\": %zu,\n", total_cell);
        fprintf(fp, "  \"kept\": { \"cells\": %zu, \"verts\": %zu, "
                    "\"faces\": %zu, \"u_span\": %.4f, \"v_span\": %.4f },\n",
                st[0].ncell, onv, onf,
                st[0].u_hi - st[0].u_lo, st[0].v_hi - st[0].v_lo);
        fprintf(fp, "  \"components\": [");
        for (g = 0; g < ncomp && (int)g < top; g++)
            fprintf(fp, "%s\n    { \"cells\": %zu, \"verts\": %zu, "
                        "\"faces\": %zu, \"u_lo\": %.4f, \"u_hi\": %.4f }",
                    g ? "," : "", st[g].ncell, st[g].nv, st[g].nf,
                    st[g].u_lo, st[g].u_hi);
        fprintf(fp, "\n  ]\n}\n");
        fclose(fp);
        fprintf(stderr, "  wrote %s\n", report);
    }
    rc = 0;
done:
    Arena_dispose(&arena);
    return rc;
}

/* ---- selftest ---------------------------------------------------------- */

static int rls_check(int cond, const char *what, int *fails)
{
    fprintf(stderr, "[ribbon_largest_strip selftest] %s -> %s\n",
            what, cond ? "ok" : "FAIL");
    if (!cond) (*fails)++;
    return cond;
}

static int rls_selftest(void)
{
    enum { NV = 8*2 + 2*3 + 4 };
    int fails = 0;
    Arena_T arena = Arena_new();
    float uv[NV*2];
    float xyz[NV*3];
    int32_t faces[64*3];
    size_t nv = 0, nf = 0, base_a = 0, base_b = 0, base_c = 0;
    int32_t *label = NULL;
    StripStat *st = NULL;
    size_t ncomp = 0;
    int c = 0, r = 0, i = 0, rc = 0;

    /* Strip A: 8 columns x 2 rows, at du = 2.  Strip B: 2 columns x 3 rows --
     * more rows, fewer cells.  Strip C: one quad.  A must win on cells. */
    base_a = nv;
    for (c = 0; c < 8; c++)
        for (r = 0; r < 2; r++) {
            uv[nv*2+0] = (float)(c * 2);
            uv[nv*2+1] = (float)r;
            xyz[nv*3+0] = (float)c; xyz[nv*3+1] = (float)r; xyz[nv*3+2] = 0.0f;
            nv++;
        }
    for (c = 0; c + 1 < 8; c++) {
        int32_t a = (int32_t)(base_a + (size_t)c * 2);
        int32_t b = a + 2;
        faces[nf*3+0] = a; faces[nf*3+1] = b;   faces[nf*3+2] = a+1; nf++;
        faces[nf*3+0] = b; faces[nf*3+1] = b+1; faces[nf*3+2] = a+1; nf++;
    }

    base_b = nv;
    for (c = 0; c < 2; c++)
        for (r = 0; r < 3; r++) {
            uv[nv*2+0] = (float)(1000 + c * 2);
            uv[nv*2+1] = (float)r;
            xyz[nv*3+0] = (float)c; xyz[nv*3+1] = (float)r; xyz[nv*3+2] = 9.0f;
            nv++;
        }
    for (r = 0; r + 1 < 3; r++) {
        int32_t a = (int32_t)(base_b + (size_t)r);
        int32_t b = a + 3;
        faces[nf*3+0] = a; faces[nf*3+1] = b;   faces[nf*3+2] = a+1; nf++;
        faces[nf*3+0] = b; faces[nf*3+1] = b+1; faces[nf*3+2] = a+1; nf++;
    }

    base_c = nv;
    for (i = 0; i < 4; i++) {
        uv[nv*2+0] = (float)(2000 + (i/2) * 2);
        uv[nv*2+1] = (float)(i % 2);
        xyz[nv*3+0] = 0.0f; xyz[nv*3+1] = 0.0f; xyz[nv*3+2] = 20.0f;
        nv++;
    }
    faces[nf*3+0] = (int32_t)base_c;
    faces[nf*3+1] = (int32_t)base_c + 2;
    faces[nf*3+2] = (int32_t)base_c + 1;
    nf++;
    faces[nf*3+0] = (int32_t)base_c + 2;
    faces[nf*3+1] = (int32_t)base_c + 3;
    faces[nf*3+2] = (int32_t)base_c + 1;
    nf++;

    rc = rls_components(arena, uv, nv, faces, nf, 2.0, 1.0,
                        &label, &ncomp, &st);
    rls_check(rc == 0, "component analysis returns 0", &fails);
    rls_check(ncomp == 3, "three disjoint strips found", &fails);
    if (rc == 0 && ncomp == 3) {
        rls_check(st[0].ncell == 16, "widest strip ranks first on cells",
                  &fails);
        rls_check(st[1].ncell == 6, "taller-but-narrower strip is second",
                  &fails);
        rls_check(st[2].ncell == 4, "single quad is last", &fails);
        rls_check(fabs((st[0].u_hi - st[0].u_lo) - 14.0) < 1e-6,
                  "winner u span measured in UV units", &fails);
        rls_check(label[st[0].root] == label[base_a],
                  "winner root maps back to the wide strip", &fails);
    }

    /* Coincident UVs must not inflate a cell count. */
    {
        float uv2[8];
        int32_t f2[6];
        int32_t *l2 = NULL;
        StripStat *s2 = NULL;
        size_t nc2 = 0;
        int j = 0, rc2 = 0;
        for (j = 0; j < 4; j++) {
            uv2[j*2+0] = (float)((j/2) * 2);
            uv2[j*2+1] = (float)(j % 2);
        }
        uv2[3*2+0] = uv2[2*2+0];
        uv2[3*2+1] = uv2[2*2+1];
        f2[0]=0; f2[1]=1; f2[2]=2; f2[3]=1; f2[4]=3; f2[5]=2;
        rc2 = rls_components(arena, uv2, 4, f2, 2, 2.0, 1.0, &l2, &nc2, &s2);
        rls_check(rc2 == 0 && nc2 == 1 && s2[0].ncell == 3,
                  "coincident UV vertices count as one cell", &fails);
    }

    /* Empty input is a clean failure, never a crash. */
    {
        int32_t *l3 = NULL;
        StripStat *s3 = NULL;
        size_t nc3 = 0;
        rls_check(rls_components(arena, uv, 0, faces, 0, 2.0, 1.0,
                                 &l3, &nc3, &s3) != 0,
                  "zero-vertex input rejected cleanly", &fails);
    }

    (void)xyz;
    Arena_dispose(&arena);
    fprintf(stderr, "\nribbon_largest_strip --selftest: %s (%d failure%s)\n",
            fails ? "FAIL" : "PASS", fails, fails == 1 ? "" : "s");
    return fails ? 3 : 0;
}

static void usage(const char *argv0)
{
    fprintf(stderr,
        "Usage: %s <in.vmesh> <out.vmesh> [--du F=2] [--dv F=1]\n"
        "          [--top N=10] [--report out.json]\n"
        "       %s --selftest\n"
        "Keeps the connected strip with the most occupied (u,v) cells.\n",
        argv0, argv0);
}

int main(int argc, char **argv)
{
    double du = 2.0, dv = 1.0;
    int top = 10, i = 0;
    const char *report = NULL;

    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return rls_selftest();
    if (argc < 3) { usage(argv[0]); return 2; }

    for (i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--du") == 0 && i + 1 < argc)
            du = atof(argv[++i]);
        else if (strcmp(argv[i], "--dv") == 0 && i + 1 < argc)
            dv = atof(argv[++i]);
        else if (strcmp(argv[i], "--top") == 0 && i + 1 < argc)
            top = atoi(argv[++i]);
        else if (strcmp(argv[i], "--report") == 0 && i + 1 < argc)
            report = argv[++i];
        else { usage(argv[0]); return 2; }
    }
    return rls_run(argv[1], argv[2], du, dv, top, report);
}
