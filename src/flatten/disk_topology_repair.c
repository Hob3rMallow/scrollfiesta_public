#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "disk_topology_repair.h"

#include "mesh_topo.h"
#include "seam_cut.h"
#include "topology_invariants.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int32_t dtr_find(int32_t *parent, int32_t v)
{
    while (parent[v] != v) {
        parent[v] = parent[parent[v]];
        v = parent[v];
    }
    return v;
}

static void dtr_union(int32_t *parent, int32_t *size, int32_t a, int32_t b)
{
    a = dtr_find(parent, a);
    b = dtr_find(parent, b);
    if (a == b) return;
    if (size[a] < size[b]) {
        int32_t t = a; a = b; b = t;
    }
    parent[b] = a;
    size[a] += size[b];
}

static int dtr_append_component(float *dst_v, int32_t *dst_f,
                                int32_t *dst_source,
                                size_t vertex_capacity, size_t face_capacity,
                                size_t *dst_nv, size_t *dst_nf,
                                const float *src_v, size_t src_nv,
                                const int32_t *src_f, size_t src_nf,
                                const int32_t *src_to_input)
{
    size_t vb = *dst_nv, fb = *dst_nf;
    if (src_nv > vertex_capacity - vb || src_nf > face_capacity - fb)
        return -1;
    if (vb > (size_t)INT32_MAX) return -1;
    memcpy(&dst_v[vb*3], src_v, src_nv*3*sizeof(float));
    for (size_t v = 0; v < src_nv; v++)
        dst_source[vb+v] = src_to_input[v];
    for (size_t f = 0; f < src_nf; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t local = src_f[f*3+(size_t)k];
            if (local < 0 || (size_t)local >= src_nv ||
                vb + (size_t)local > (size_t)INT32_MAX)
                return -1;
            dst_f[(fb+f)*3+(size_t)k] = (int32_t)(vb+(size_t)local);
        }
    }
    *dst_nv += src_nv;
    *dst_nf += src_nf;
    return 0;
}

int DiskTopologyRepair_process(Arena_T arena,
                               const float *verts, size_t nv,
                               const int32_t *faces, size_t nf,
                               float **out_verts, size_t *out_nv,
                               int32_t **out_faces, size_t *out_nf,
                               int32_t **out_vertex_source,
                               DiskTopologyRepairStats *stats)
{
    DiskTopologyRepairStats st;
    int32_t *parent = NULL, *usize = NULL, *root_to_component = NULL;
    int32_t *vertex_component = NULL, *global_to_local = NULL;
    uint8_t *used = NULL;
    size_t *component_faces = NULL, *offset = NULL, *cursor = NULL;
    int32_t *face_index = NULL;
    float *result_v = NULL;
    int32_t *result_f = NULL, *result_source = NULL;
    size_t result_nv = 0, result_nf = 0, ncomponents = 0;
    size_t vertex_capacity;
    int rc = -1;

    memset(&st, 0, sizeof st);
    st.vertices_in = nv; st.faces_in = nf;
    if (stats) *stats = st;
    if (!arena || !verts || !faces || !out_verts || !out_nv ||
        !out_faces || !out_nf || !out_vertex_source ||
        nv == 0 || nf == 0 || nv > (size_t)INT32_MAX ||
        nf > SIZE_MAX/3)
        return -1;
    *out_verts = NULL; *out_nv = 0;
    *out_faces = NULL; *out_nf = 0; *out_vertex_source = NULL;

    parent = (int32_t *)malloc(nv*sizeof(int32_t));
    usize = (int32_t *)malloc(nv*sizeof(int32_t));
    used = (uint8_t *)calloc(nv, 1);
    root_to_component = (int32_t *)malloc(nv*sizeof(int32_t));
    vertex_component = (int32_t *)malloc(nv*sizeof(int32_t));
    global_to_local = (int32_t *)malloc(nv*sizeof(int32_t));
    if (!parent || !usize || !used || !root_to_component ||
        !vertex_component || !global_to_local)
        goto cleanup;
    for (size_t v = 0; v < nv; v++) {
        parent[v] = (int32_t)v; usize[v] = 1;
        root_to_component[v] = -1;
        vertex_component[v] = -1;
        global_to_local[v] = -1;
    }
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f*3], b = faces[f*3+1], c = faces[f*3+2];
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= nv || (size_t)b >= nv || (size_t)c >= nv ||
            a == b || b == c || c == a)
            goto cleanup;
        used[a] = used[b] = used[c] = 1;
        dtr_union(parent, usize, a, b);
        dtr_union(parent, usize, b, c);
    }
    for (size_t v = 0; v < nv; v++) {
        if (!used[v]) continue;
        int32_t root = dtr_find(parent, (int32_t)v);
        if (root_to_component[root] < 0)
            root_to_component[root] = (int32_t)ncomponents++;
        vertex_component[v] = root_to_component[root];
    }
    if (ncomponents == 0) goto cleanup;
    st.components = ncomponents;

    component_faces = (size_t *)calloc(ncomponents, sizeof(size_t));
    offset = (size_t *)malloc((ncomponents+1)*sizeof(size_t));
    cursor = (size_t *)malloc(ncomponents*sizeof(size_t));
    face_index = (int32_t *)malloc(nf*sizeof(int32_t));
    if (!component_faces || !offset || !cursor || !face_index) goto cleanup;
    for (size_t f = 0; f < nf; f++) {
        int32_t component = vertex_component[faces[f*3]];
        if (component < 0 ||
            component != vertex_component[faces[f*3+1]] ||
            component != vertex_component[faces[f*3+2]])
            goto cleanup;
        component_faces[component]++;
    }
    offset[0] = 0;
    for (size_t c = 0; c < ncomponents; c++) {
        offset[c+1] = offset[c] + component_faces[c];
        cursor[c] = offset[c];
    }
    for (size_t f = 0; f < nf; f++) {
        int32_t component = vertex_component[faces[f*3]];
        face_index[cursor[component]++] = (int32_t)f;
    }

    vertex_capacity = nf*3; /* one independent vertex per face corner */
    result_v = (float *)ARENA_ALLOC(arena,
        vertex_capacity*3*sizeof(float));
    result_f = (int32_t *)ARENA_ALLOC(arena, nf*3*sizeof(int32_t));
    result_source = (int32_t *)ARENA_ALLOC(arena,
        vertex_capacity*sizeof(int32_t));

    for (size_t component = 0; component < ncomponents; component++) {
        size_t cnf = component_faces[component], cnv = 0;
        size_t corner_capacity = cnf*3;
        float *cv = (float *)malloc(corner_capacity*3*sizeof(float));
        int32_t *cf = (int32_t *)malloc(cnf*3*sizeof(int32_t));
        int32_t *local_to_global =
            (int32_t *)malloc(corner_capacity*sizeof(int32_t));
        Arena_T component_arena = NULL;
        MeshTopoInfo before;
        TopologyAuditReport exact_before, exact_after;
        const float *append_v = NULL;
        const int32_t *append_f = NULL, *append_source = NULL;
        size_t append_nv = 0, append_nf = 0;
        SeamCutResult cut;
        int component_ok = 0;

        memset(&exact_before, 0, sizeof exact_before);
        memset(&exact_after, 0, sizeof exact_after);

        if (!cv || !cf || !local_to_global) {
            free(cv); free(cf); free(local_to_global);
            goto cleanup;
        }
        for (size_t q = offset[component]; q < offset[component+1]; q++) {
            size_t f = (size_t)face_index[q];
            for (int k = 0; k < 3; k++) {
                int32_t global = faces[f*3+(size_t)k];
                int32_t local = global_to_local[global];
                if (local < 0) {
                    if (cnv >= corner_capacity || cnv > (size_t)INT32_MAX)
                        goto component_cleanup;
                    local = (int32_t)cnv;
                    global_to_local[global] = local;
                    local_to_global[cnv] = global;
                    memcpy(&cv[cnv*3], &verts[(size_t)global*3],
                           3*sizeof(float));
                    cnv++;
                }
                cf[(q-offset[component])*3+(size_t)k] = local;
            }
        }
        component_arena = Arena_new();
        if (TopologyAudit_analyze(
                cv, cnv, cf, cnf, NULL, &exact_before) != 0) {
            fprintf(stderr,
                    "  [disk-repair] component %zu exact topology audit "
                    "failed: %s\n", component, exact_before.error);
            goto component_cleanup;
        }
        if (exact_before.face_components != 1 ||
            exact_before.isolated_vertices != 0 ||
            exact_before.invalid_surface_components != 0 ||
            exact_before.nonorientable_components != 0 ||
            exact_before.inconsistent_winding_components != 0 ||
            !exact_before.betti_complete || exact_before.component == NULL ||
            !exact_before.component[0].surface_valid ||
            !exact_before.component[0].orientable ||
            !exact_before.component[0].input_winding_consistent) {
            const TopologyComponentInvariant *q =
                exact_before.component ? &exact_before.component[0] : NULL;
            fprintf(stderr,
                    "  [disk-repair] component %zu rejected before SeamCut: "
                    "face-components=%zu isolated=%zu invalid=%zu "
                    "nonorientable=%zu winding=%zu nonmanifold=%zuE/%zuV "
                    "irregular-boundary=%zu\n",
                    component, exact_before.face_components,
                    exact_before.isolated_vertices,
                    exact_before.invalid_surface_components,
                    exact_before.nonorientable_components,
                    exact_before.inconsistent_winding_components,
                    q ? q->nonmanifold_edges : 0,
                    q ? q->nonmanifold_vertices : 0,
                    q ? q->boundary_irregular_vertices : 0);
            goto component_cleanup;
        }
        const TopologyComponentInvariant *exact = &exact_before.component[0];
        if (MeshTopo_analyze(component_arena, cv, cnv, cf, cnf, &before) != 0 ||
            before.n_nonmanifold_edges != 0 || before.n_unref_verts != 0 ||
            before.genus < -1e-6 || exact->orientable_genus < 0 ||
            before.n_boundary_loops != (long)exact->boundary_loops ||
            fabs(before.genus - (double)exact->orientable_genus) > 1e-6 ||
            before.is_disk != exact->homeomorphic_to_disk) {
            fprintf(stderr,
                    "  [disk-repair] component %zu lightweight/exact "
                    "topology disagreement (loops=%ld/%zu genus=%.3f/%lld "
                    "disk=%d/%d)\n",
                    component, before.n_boundary_loops,
                    exact->boundary_loops, before.genus,
                    (long long)exact->orientable_genus,
                    before.is_disk, exact->homeomorphic_to_disk);
            goto component_cleanup;
        }
        if (exact->boundary_loops > 0)
            st.boundary_loops_before += exact->boundary_loops;

        if (exact->homeomorphic_to_disk) {
            st.already_disks++;
            st.boundary_loops_after++;
            append_v = cv; append_nv = cnv;
            append_f = cf; append_nf = cnf;
            append_source = local_to_global;
        } else {
            SeamCutOpts options;
            int cut_rc;
            options.accept_frac = 1.0; /* topology only; no distortion terminals */
            options.r_max = 5;
            memset(&cut, 0, sizeof cut);
            cut_rc = SeamCut_run(component_arena, cv, cnv, cf, cnf,
                                 &options, &cut);
            if (cut_rc != 0 || !cut.is_disk_after ||
                cut.loops_after != 1 || cut.genus_after != 0 ||
                TopologyAudit_analyze(cut.verts, cut.nv,
                                      cut.faces, cut.nf, NULL,
                                      &exact_after) != 0 ||
                exact_after.face_components != 1 ||
                exact_after.isolated_vertices != 0 ||
                !exact_after.all_components_are_disks ||
                exact_after.invalid_surface_components != 0 ||
                exact_after.nonorientable_components != 0 ||
                exact_after.inconsistent_winding_components != 0) {
                fprintf(stderr,
                        "  [disk-repair] component %zu failed: "
                        "V=%zu F=%zu loops=%ld genus=%.3f -> "
                        "rc=%d loops=%ld genus=%ld disk=%d seam=%ld "
                        "exact-disk=%d invalid=%zu winding=%zu\n",
                        component, cnv, cnf, before.n_boundary_loops,
                        before.genus, cut_rc, cut.loops_after,
                        cut.genus_after, cut.is_disk_after, cut.seam_edges,
                        exact_after.all_components_are_disks,
                        exact_after.invalid_surface_components,
                        exact_after.inconsistent_winding_components);
                goto component_cleanup;
            }
            st.repaired_components++;
            if (cut.genus_before > 0)
                st.handles_opened += (size_t)cut.genus_before;
            st.seam_edges += cut.seam_edges > 0 ? (size_t)cut.seam_edges : 0;
            st.boundary_loops_after++;
            append_v = cut.verts; append_nv = cut.nv;
            append_f = cut.faces; append_nf = cut.nf;
            /* SeamCut maps new vertices into the local component. */
            int32_t *mapped = (int32_t *)ARENA_ALLOC(
                component_arena, append_nv*sizeof(int32_t));
            for (size_t v = 0; v < append_nv; v++) {
                int32_t local = cut.vmap[v];
                if (local < 0 || (size_t)local >= cnv)
                    goto component_cleanup;
                mapped[v] = local_to_global[local];
            }
            append_source = mapped;
        }
        if (dtr_append_component(result_v, result_f, result_source,
                                 vertex_capacity, nf,
                                 &result_nv, &result_nf,
                                 append_v, append_nv, append_f, append_nf,
                                 append_source) != 0)
            goto component_cleanup;
        component_ok = 1;

component_cleanup:
        for (size_t v = 0; v < cnv; v++)
            global_to_local[local_to_global[v]] = -1;
        TopologyAudit_dispose(&exact_after);
        TopologyAudit_dispose(&exact_before);
        if (component_arena) Arena_dispose(&component_arena);
        free(local_to_global); free(cf); free(cv);
        if (!component_ok) {
            st.failed_components++;
            goto cleanup;
        }
    }

    st.vertices_out = result_nv; st.faces_out = result_nf;
    *out_verts = result_v; *out_nv = result_nv;
    *out_faces = result_f; *out_nf = result_nf;
    *out_vertex_source = result_source;
    if (stats) *stats = st;
    rc = 0;

cleanup:
    if (rc != 0 && stats) *stats = st;
    free(face_index); free(cursor); free(offset); free(component_faces);
    free(global_to_local); free(vertex_component); free(root_to_component);
    free(used); free(usize); free(parent);
    return rc;
}

static void dtr_build_fixture(float *verts, int32_t *faces,
                              size_t *out_nv, size_t *out_nf)
{
    const int n = 8;
    size_t nv = 0, nf = 0;
    /* One disk. */
    const float tri[9] = { -4,0,0, -3,1,0, -3,-1,0 };
    memcpy(verts, tri, sizeof tri);
    faces[0] = 0; faces[1] = 1; faces[2] = 2;
    nv = 3; nf = 1;
    /* One annulus: a two-ring open cylinder. */
    size_t base = nv;
    for (int row = 0; row < 2; row++)
        for (int i = 0; i < n; i++) {
            double a = 6.2831853071795864769*(double)i/(double)n;
            verts[nv*3] = (float)row;
            verts[nv*3+1] = (float)cos(a);
            verts[nv*3+2] = (float)sin(a);
            nv++;
        }
    for (int i = 0; i < n; i++) {
        int j = (i+1)%n;
        int32_t a = (int32_t)(base+(size_t)i);
        int32_t b = (int32_t)(base+(size_t)j);
        int32_t c = (int32_t)(base+(size_t)n+(size_t)i);
        int32_t d = (int32_t)(base+(size_t)n+(size_t)j);
        faces[nf*3] = a; faces[nf*3+1] = b; faces[nf*3+2] = d; nf++;
        faces[nf*3] = a; faces[nf*3+1] = d; faces[nf*3+2] = c; nf++;
    }
    *out_nv = nv; *out_nf = nf;
}

int DiskTopologyRepair_selftest(void)
{
    float verts[19*3];
    int32_t faces[17*3];
    size_t nv, nf, onv, onf, onv2, onf2;
    float *ov = NULL, *ov2 = NULL;
    int32_t *of = NULL, *map = NULL, *of2 = NULL, *map2 = NULL;
    DiskTopologyRepairStats first = {0}, second = {0};
    Arena_T a = Arena_new();
    int ok, invalid_ok;

    dtr_build_fixture(verts, faces, &nv, &nf);
    ok = DiskTopologyRepair_process(a, verts, nv, faces, nf,
                                    &ov, &onv, &of, &onf, &map, &first) == 0 &&
         first.components == 2 && first.already_disks == 1 &&
         first.repaired_components == 1 && first.failed_components == 0 &&
         onf == nf;
    if (ok)
        ok = DiskTopologyRepair_process(a, ov, onv, of, onf,
                                        &ov2, &onv2, &of2, &onf2, &map2,
                                        &second) == 0 &&
             second.components == 2 && second.already_disks == 2 &&
             second.repaired_components == 0 && onf2 == onf;
    fprintf(stderr,
            "[selftest] disk topology repair: components=2 repaired=%zu "
            "idempotent_disks=%zu -> %s\n",
            first.repaired_components, ok ? second.already_disks : 0,
            ok ? "ok" : "FAIL");
    {
        const float bowtie_v[5*3] = {
            0,0,0, 0,1,0, 0,0,1, 0,-1,0, 0,0,-1
        };
        const int32_t bowtie_f[2*3] = {0,1,2, 0,3,4};
        DiskTopologyRepairStats invalid = {0};
        float *bad_v = NULL;
        int32_t *bad_f = NULL, *bad_map = NULL;
        size_t bad_nv = 0, bad_nf = 0;
        invalid_ok = DiskTopologyRepair_process(
                         a, bowtie_v, 5, bowtie_f, 2,
                         &bad_v, &bad_nv, &bad_f, &bad_nf,
                         &bad_map, &invalid) != 0 &&
                     invalid.failed_components == 1;
        fprintf(stderr,
                "[selftest] disk topology repair rejects vertex bowtie -> %s\n",
                invalid_ok ? "ok" : "FAIL");
    }
    Arena_dispose(&a);
    return ok && invalid_ok ? 0 : 1;
}
