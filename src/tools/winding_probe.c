/* ============================================================================
 * winding_probe -- exact generalized-winding-number staircase probe.
 *
 * THE experiment (2026-08-31): the winding authority our pipeline lacks may
 * simply BE the generalized winding number field of the oriented prediction
 * soup (Jacobson 2013; Barill 2018 for soups; Xie/Hafner/Wojtan 2026 fig. 5
 * "swiss roll").  Between wraps, w(q) is the single-valued universal-cover
 * coordinate of the inter-wrap spiral channel; across the sheet it jumps by
 * exactly +-1.  Per vertex, w+- = GWN(v +- eps*n) gives a turn coordinate
 * ((w+ + w-)/2) and a built-in confidence (jump = w- - w+, expected 1;
 * 2 marks a same-orientation duplicate shell, 0 an orientation cancel).
 *
 * BRUTE FORCE ON PURPOSE: exact Van Oosterom solid-angle sum over every
 * face for every query, OpenMP over vertices.  No octrees until the
 * staircase is proven.
 *
 *   winding_probe <vmesh_list.txt> <out_dir> --umb-y Y --umb-x X
 *                 [--u-per-turn F=1200] [--eps F=1.5]
 *                 [--backend auto|direct|exact|fast] [--beta F=2]
 *   winding_probe --selftest
 *
 * Writes <out_dir>/probe_ribbon.vmesh with UV = (w * u_per_turn, z - zmin),
 * sidecars probe_winding.f32 / probe_jump.f32, and prints the jump
 * histogram (the staircase verdict).  Components are first orientation-
 * voted radially (a wrongly flipped component SUBTRACTS from the field).
 * ==========================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/union_find.h"
#include "../flatten/winding_field.h"
#include "../flatten/winding_mrf.h"

#define WP_MAX_MESHES 512
#define WP_PI 3.14159265358979323846

typedef struct {
    float *verts;      /* [nv*3] (z,y,x) */
    int32_t *faces;    /* [nf*3] */
    size_t nv, nf;
} WpSoup;

static int wp_load_list(Arena_T arena, const char *list_path, WpSoup *soup)
{
    char line[1024];
    FILE *handle = fopen(list_path, "r");
    size_t nv = 0, nf = 0, count = 0;
    MeshBinData mesh[WP_MAX_MESHES];
    if (handle == NULL) {
        fprintf(stderr, "winding_probe: cannot open %s\n", list_path);
        return -1;
    }
    while (fgets(line, sizeof line, handle) != NULL &&
           count < WP_MAX_MESHES) {
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r' ||
                           line[len-1] == ' '))
            line[--len] = '\0';
        if (len == 0) continue;
        memset(&mesh[count], 0, sizeof mesh[count]);
        if (MeshBin_read_arena(arena, line, &mesh[count]) != 0) {
            fprintf(stderr, "winding_probe: cannot read %s\n", line);
            fclose(handle);
            return -1;
        }
        nv += mesh[count].nv;
        nf += mesh[count].nf;
        count++;
    }
    fclose(handle);
    if (count == 0 || nv == 0 || nf == 0) {
        fprintf(stderr, "winding_probe: empty input list\n");
        return -1;
    }
    soup->verts = (float *)ARENA_ALLOC(arena,
                                       (long)(nv * 3 * sizeof(float)));
    soup->faces = (int32_t *)ARENA_ALLOC(arena,
                                         (long)(nf * 3 * sizeof(int32_t)));
    soup->nv = nv;
    soup->nf = nf;
    {
        size_t vat = 0, fat = 0;
        for (size_t m = 0; m < count; m++) {
            memcpy(soup->verts + vat * 3, mesh[m].verts,
                   mesh[m].nv * 3 * sizeof(float));
            for (size_t f = 0; f < mesh[m].nf * 3; f++)
                soup->faces[fat * 3 + f] =
                    mesh[m].faces[f] + (int32_t)vat;
            vat += mesh[m].nv;
            fat += mesh[m].nf;
        }
    }
    fprintf(stderr, "[winding-probe] loaded %zu mesh(es): nv=%zu nf=%zu\n",
            count, soup->nv, soup->nf);
    return 0;
}

/* radial orientation vote per connected component; flips losing faces */
static void wp_orient_radially(Arena_T arena, WpSoup *soup,
                               double umb_y, double umb_x)
{
    UnionFind graph = UF_new(arena, (int32_t)soup->nv);
    double *vote = NULL;
    int32_t flipped_components = 0;
    size_t flipped_faces = 0;
    for (size_t f = 0; f < soup->nf; f++) {
        uf_union(&graph, soup->faces[f*3], soup->faces[f*3+1]);
        uf_union(&graph, soup->faces[f*3], soup->faces[f*3+2]);
    }
    vote = (double *)ARENA_CALLOC(arena, (long)soup->nv,
                                  (long)sizeof(double));
    for (size_t f = 0; f < soup->nf; f++) {
        const float *a = &soup->verts[(size_t)soup->faces[f*3] * 3];
        const float *b = &soup->verts[(size_t)soup->faces[f*3+1] * 3];
        const float *c = &soup->verts[(size_t)soup->faces[f*3+2] * 3];
        double e1[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
        double e2[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
        /* normal in (z,y,x) */
        double n[3] = { e1[1]*e2[2] - e1[2]*e2[1],
                        e1[2]*e2[0] - e1[0]*e2[2],
                        e1[0]*e2[1] - e1[1]*e2[0] };
        double cy = (a[1]+b[1]+c[1])/3.0 - umb_y;
        double cx = (a[2]+b[2]+c[2])/3.0 - umb_x;
        double rr = hypot(cy, cx);
        int32_t root = uf_find(&graph, soup->faces[f*3]);
        if (rr > 1e-9)
            vote[root] += (n[1]*cy + n[2]*cx) / rr;
    }
    for (size_t f = 0; f < soup->nf; f++) {
        int32_t root = uf_find(&graph, soup->faces[f*3]);
        if (vote[root] < 0.0) {
            int32_t swap = soup->faces[f*3+1];
            soup->faces[f*3+1] = soup->faces[f*3+2];
            soup->faces[f*3+2] = swap;
            flipped_faces++;
        }
    }
    for (size_t v = 0; v < soup->nv; v++)
        if (uf_find(&graph, (int32_t)v) == (int32_t)v && vote[v] < 0.0)
            flipped_components++;
    fprintf(stderr,
            "[winding-probe] radial orientation: %d component(s) flipped "
            "(%zu faces)\n", flipped_components, flipped_faces);
}

/* exact GWN at q over the whole soup (Van Oosterom & Strackee solid angle) */
static double wp_gwn(const WpSoup *soup, const double q[3])
{
    double total = 0.0;
    for (size_t f = 0; f < soup->nf; f++) {
        const float *pa = &soup->verts[(size_t)soup->faces[f*3] * 3];
        const float *pb = &soup->verts[(size_t)soup->faces[f*3+1] * 3];
        const float *pc = &soup->verts[(size_t)soup->faces[f*3+2] * 3];
        double a[3] = { pa[0]-q[0], pa[1]-q[1], pa[2]-q[2] };
        double b[3] = { pb[0]-q[0], pb[1]-q[1], pb[2]-q[2] };
        double c[3] = { pc[0]-q[0], pc[1]-q[1], pc[2]-q[2] };
        double la = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
        double lb = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
        double lc = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
        double det = a[0]*(b[1]*c[2] - b[2]*c[1]) -
                     a[1]*(b[0]*c[2] - b[2]*c[0]) +
                     a[2]*(b[0]*c[1] - b[1]*c[0]);
        double denom = la*lb*lc + (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])*lc +
                       (b[0]*c[0]+b[1]*c[1]+b[2]*c[2])*la +
                       (c[0]*a[0]+c[1]*a[1]+c[2]*a[2])*lb;
        total += 2.0 * atan2(det, denom);
    }
    return total / (4.0 * WP_PI);
}

static int wp_run(const char *list_path, const char *field_list_path,
                  const char *out_dir,
                  double umb_y, double umb_x, double u_per_turn, double eps,
                  WindingFieldBackend backend, double beta)
{
    Arena_T arena = Arena_new();
    WpSoup soup;
    WpSoup field;
    const WpSoup *field_soup = NULL;
    float *normal = NULL, *uv = NULL, *w_out = NULL, *jump_out = NULL;
    double zmin = 1e300, wmin = 1e300, wmax = -1e300;
    size_t clean = 0, doubled = 0, cancelled = 0, other = 0;
    WindingField_T winding_field = NULL;
    WindingFieldOptions field_options;
    WindingFieldStats field_stats;
    WindingFieldEvalStats eval_stats;
    char path[1024];
    if (arena == NULL) return -1;
    memset(&soup, 0, sizeof soup);
    memset(&field, 0, sizeof field);
    if (wp_load_list(arena, list_path, &soup) != 0) {
        Arena_dispose(&arena);
        return -1;
    }
    wp_orient_radially(arena, &soup, umb_y, umb_x);
    /* The FIELD may integrate over a larger soup than the queried block:
     * GWN's global validity needs enclosing coverage (measured 2026-08-31:
     * a 3x3x3 block's own arcs compress ~25 wraps into ~6 apparent turns
     * of w while the local jump stays 1 at 98.2% -- field global, queries
     * local is the working configuration). */
    if (field_list_path != NULL) {
        if (wp_load_list(arena, field_list_path, &field) != 0) {
            Arena_dispose(&arena);
            return -1;
        }
        wp_orient_radially(arena, &field, umb_y, umb_x);
        field_soup = &field;
        fprintf(stderr, "[winding-probe] field soup: nv=%zu nf=%zu\n",
                field.nv, field.nf);
    } else {
        field_soup = &soup;
    }
    normal = (float *)ARENA_CALLOC(arena, (long)(soup.nv * 3),
                                   (long)sizeof(float));
    for (size_t f = 0; f < soup.nf; f++) {
        const float *a = &soup.verts[(size_t)soup.faces[f*3] * 3];
        const float *b = &soup.verts[(size_t)soup.faces[f*3+1] * 3];
        const float *c = &soup.verts[(size_t)soup.faces[f*3+2] * 3];
        double e1[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
        double e2[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
        double n[3] = { e1[1]*e2[2] - e1[2]*e2[1],
                        e1[2]*e2[0] - e1[0]*e2[2],
                        e1[0]*e2[1] - e1[1]*e2[0] };
        for (int k = 0; k < 3; k++) {
            normal[(size_t)soup.faces[f*3] * 3 + (size_t)k] += (float)n[k];
            normal[(size_t)soup.faces[f*3+1] * 3 + (size_t)k] += (float)n[k];
            normal[(size_t)soup.faces[f*3+2] * 3 + (size_t)k] += (float)n[k];
        }
    }
    w_out = (float *)ARENA_ALLOC(arena, (long)(soup.nv * sizeof(float)));
    jump_out = (float *)ARENA_ALLOC(arena, (long)(soup.nv * sizeof(float)));
    uv = (float *)ARENA_ALLOC(arena, (long)(soup.nv * 2 * sizeof(float)));
    for (size_t v = 0; v < soup.nv; v++)
        if (soup.verts[v*3] < zmin) zmin = soup.verts[v*3];
    WindingField_default_options(&field_options);
    field_options.backend = backend;
    field_options.beta = beta;
    if (WindingField_build(
            arena, field_soup->verts, field_soup->nv,
            field_soup->faces, field_soup->nf, &field_options,
            &winding_field, &field_stats) != 0) {
        fprintf(stderr, "winding_probe: field build failed\n");
        Arena_dispose(&arena);
        return -1;
    }
    fprintf(stderr,
            "[winding-probe] field engine: nodes=%zu leaves=%zu "
            "exterior=%zu (multiplicity=%zu, %.3f/face) beta=%.2f; "
            "evaluating %zu two-sided queries on <= half CPU...\n",
            field_stats.bvh_nodes, field_stats.bvh_leaves,
            field_stats.exterior_edges, field_stats.exterior_multiplicity,
            (double)field_stats.exterior_multiplicity /
                (double)field_stats.faces,
            field_options.beta, soup.nv);
    if (WindingField_evaluate_sides(
            winding_field, soup.verts, normal, soup.nv, eps, backend,
            w_out, jump_out, &eval_stats) != 0) {
        fprintf(stderr, "winding_probe: field evaluation failed\n");
        Arena_dispose(&arena);
        return -1;
    }
    fprintf(stderr,
            "[winding-probe] field evaluation: backend=%s queries=%zu "
            "exact_tri=%zu approximate_nodes=%zu ray_fallback=%zu "
            "invalid_normals=%zu\n",
            eval_stats.backend_used == WINDING_FIELD_DIRECT ? "direct" :
            eval_stats.backend_used == WINDING_FIELD_BOUNDARY_EXACT
                ? "boundary-exact" : "fast",
            eval_stats.queries, eval_stats.exact_triangle_evaluations,
            eval_stats.approximate_node_evaluations,
            eval_stats.degenerate_ray_fallbacks,
            eval_stats.invalid_normals);
    for (size_t v = 0; v < soup.nv; v++) {
        double j = jump_out[v];
        if (!isfinite(j)) { other++; continue; }
        if (fabs(j - 1.0) <= 0.5) clean++;
        else if (fabs(j - 2.0) <= 0.5) doubled++;
        else if (fabs(j) <= 0.5) cancelled++;
        else other++;
        if (isfinite(w_out[v])) {
            if (w_out[v] < wmin) wmin = w_out[v];
            if (w_out[v] > wmax) wmax = w_out[v];
        }
        uv[v*2] = (float)((double)w_out[v] * u_per_turn);
        uv[v*2+1] = (float)((double)soup.verts[v*3] - zmin);
    }
    /* re-zero u */
    {
        double umin = 1e300;
        for (size_t v = 0; v < soup.nv; v++)
            if (isfinite(uv[v*2]) && uv[v*2] < umin) umin = uv[v*2];
        for (size_t v = 0; v < soup.nv; v++)
            uv[v*2] = isfinite(uv[v*2]) ? (float)(uv[v*2] - umin) : 0.0f;
    }
    fprintf(stderr,
            "[winding-probe] STAIRCASE VERDICT: jump~1 (clean sheet) "
            "%zu/%zu = %.1f%%; jump~2 (duplicate shell) %.1f%%; jump~0 "
            "(orientation cancel) %.1f%%; other %.1f%%\n",
            clean, soup.nv, 100.0*(double)clean/(double)soup.nv,
            100.0*(double)doubled/(double)soup.nv,
            100.0*(double)cancelled/(double)soup.nv,
            100.0*(double)other/(double)soup.nv);
    fprintf(stderr, "[winding-probe] w range: %.2f .. %.2f turns\n",
            wmin, wmax);
    snprintf(path, sizeof path, "%s/probe_ribbon.vmesh", out_dir);
    if (MeshBin_write(path, soup.verts, soup.nv, soup.faces, soup.nf,
                      uv) != 0) {
        fprintf(stderr, "winding_probe: cannot write %s\n", path);
        Arena_dispose(&arena);
        return -1;
    }
    {
        FILE *handle = NULL;
        snprintf(path, sizeof path, "%s/probe_winding.f32", out_dir);
        handle = fopen(path, "wb");
        if (handle) {
            fwrite(w_out, sizeof(float), soup.nv, handle);
            fclose(handle);
        }
        snprintf(path, sizeof path, "%s/probe_jump.f32", out_dir);
        handle = fopen(path, "wb");
        if (handle) {
            fwrite(jump_out, sizeof(float), soup.nv, handle);
            fclose(handle);
        }
    }
    fprintf(stderr, "[winding-probe] wrote %s/probe_ribbon.vmesh "
            "(+winding/jump sidecars)\n", out_dir);
    Arena_dispose(&arena);
    return 0;
}

/* ---- selftest ------------------------------------------------------------ */

static int wp_selftest(void)
{
    /* closed unit cube (12 tris, outward-oriented): w(center)=1, w(far)=0 */
    static const float cube_verts[8*3] = {
        0,0,0, 0,0,1, 0,1,0, 0,1,1, 1,0,0, 1,0,1, 1,1,0, 1,1,1 };
    static const int32_t cube_faces[12*3] = {
        0,1,3, 0,3,2,  4,6,7, 4,7,5,  0,4,5, 0,5,1,
        2,3,7, 2,7,6,  0,2,6, 0,6,4,  1,5,7, 1,7,3 };
    WpSoup soup;
    double center[3] = { 0.5, 0.5, 0.5 };
    double outside[3] = { 5.0, 5.0, 5.0 };
    double w_in = 0.0, w_out2 = 0.0;
    int fails = 0;
    soup.verts = (float *)cube_verts;
    soup.faces = (int32_t *)cube_faces;
    soup.nv = 8;
    soup.nf = 12;
    w_in = wp_gwn(&soup, center);
    w_out2 = wp_gwn(&soup, outside);
    if (fabs(fabs(w_in) - 1.0) > 1e-9 || fabs(w_out2) > 1e-9) {
        fprintf(stderr,
                "[winding-probe selftest] FAIL: cube w_in=%.6f w_out=%.6f\n",
                w_in, w_out2);
        fails++;
    }
    /* jump across one face: +-eps around the face center along +z-ish */
    {
        double qa[3] = { 1.2, 0.5, 0.5 };  /* just outside the z=1 face */
        double qb[3] = { 0.8, 0.5, 0.5 };  /* just inside */
        double jump = fabs(wp_gwn(&soup, qb) - wp_gwn(&soup, qa));
        if (fabs(jump - 1.0) > 1e-9) {
            fprintf(stderr,
                    "[winding-probe selftest] FAIL: face jump=%.6f\n", jump);
            fails++;
        }
    }
    fprintf(stderr, "[winding-probe selftest] %s (%d failure(s))\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    fails += WindingField_selftest();
    fails += WindingMRF_selftest();
    return fails == 0 ? 0 : 1;
}

int main(int argc, char **argv)
{
    double umb_y = 0.0, umb_x = 0.0, u_per_turn = 1200.0, eps = 1.5;
    double beta = 2.0;
    const char *field_list = NULL;
    WindingFieldBackend backend = WINDING_FIELD_AUTO;
    int have_y = 0, have_x = 0;
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return wp_selftest();
    if (argc < 4) {
        fprintf(stderr,
                "usage: winding_probe <vmesh_list.txt> <out_dir> "
                "--umb-y Y --umb-x X [--field-list F.txt] "
                "[--u-per-turn F] [--eps F] "
                "[--backend auto|direct|exact|fast] [--beta F]\n"
                "       winding_probe --selftest\n");
        return 2;
    }
    for (int i = 3; i + 1 < argc; i += 2) {
        if (strcmp(argv[i], "--umb-y") == 0) {
            umb_y = atof(argv[i+1]); have_y = 1;
        } else if (strcmp(argv[i], "--umb-x") == 0) {
            umb_x = atof(argv[i+1]); have_x = 1;
        } else if (strcmp(argv[i], "--field-list") == 0) {
            field_list = argv[i+1];
        } else if (strcmp(argv[i], "--u-per-turn") == 0) {
            u_per_turn = atof(argv[i+1]);
        } else if (strcmp(argv[i], "--eps") == 0) {
            eps = atof(argv[i+1]);
        } else if (strcmp(argv[i], "--beta") == 0) {
            beta = atof(argv[i+1]);
        } else if (strcmp(argv[i], "--backend") == 0) {
            if (strcmp(argv[i+1], "auto") == 0)
                backend = WINDING_FIELD_AUTO;
            else if (strcmp(argv[i+1], "direct") == 0)
                backend = WINDING_FIELD_DIRECT;
            else if (strcmp(argv[i+1], "exact") == 0)
                backend = WINDING_FIELD_BOUNDARY_EXACT;
            else if (strcmp(argv[i+1], "fast") == 0)
                backend = WINDING_FIELD_FAST;
            else {
                fprintf(stderr, "winding_probe: bad backend %s\n", argv[i+1]);
                return 2;
            }
        } else {
            fprintf(stderr, "winding_probe: unknown arg %s\n", argv[i]);
            return 2;
        }
    }
    if (!have_y || !have_x) {
        fprintf(stderr, "winding_probe: need --umb-y and --umb-x\n");
        return 2;
    }
    return wp_run(argv[1], field_list, argv[2], umb_y, umb_x, u_per_turn,
                  eps, backend, beta) == 0 ? 0 : 1;
}
