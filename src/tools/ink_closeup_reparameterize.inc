/* Metric recovery after completion. SLIM changes only UV: the fitted material,
 * face order, provenance and oriented normals are kept exactly. Completion's
 * fixed seed coordinates are an initialization, not a metric certificate. */
#include "ink_closeup_cones.inc"
#include "ink_closeup_material.inc"

static int ic_slim_linear(Arena_T arena, void *context, size_t nv, size_t ne,
                          const double *diagonal, const int32_t *edge_a,
                          const int32_t *edge_b, const double *blocks,
                          const double *rhs, double *result)
{
    (void)context;
    if (nv > (size_t)INT_MAX / 3 || ne > ((size_t)INT_MAX - 3 * nv) / 4) return -1;
    Arena_Mark mark = Arena_save(arena);
    size_t count = 3 * nv + 4 * ne, at = 0;
    int *row = ARENA_ALLOC(arena, count * sizeof *row);
    int *col = ARENA_ALLOC(arena, count * sizeof *col);
    double *value = ARENA_ALLOC(arena, count * sizeof *value);
    for (size_t v = 0; v < nv; v++) {
        int x = (int)(2 * v), y = x + 1;
        row[at] = col[at] = x; value[at++] = diagonal[3 * v];
        row[at] = y; col[at] = x; value[at++] = diagonal[3 * v + 1];
        row[at] = col[at] = y; value[at++] = diagonal[3 * v + 2];
    }
    for (size_t e = 0; e < ne; e++) {
        if (edge_a[e] < 0 || edge_b[e] < 0 || (size_t)edge_a[e] >= nv ||
            (size_t)edge_b[e] >= nv || edge_a[e] == edge_b[e]) {
            Arena_restore(arena, mark); return -1;
        }
        for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) {
            int i = 2 * edge_a[e] + a, j = 2 * edge_b[e] + b;
            row[at] = i > j ? i : j; col[at] = i > j ? j : i;
            value[at++] = blocks[3 * e + (size_t)(a + b)];
        }
    }
    int rc = Sparse_solve_sym((int)(2 * nv), (int)count, row, col, value, rhs, result, SPARSE_SPD);
    Arena_restore(arena, mark);
    return rc;
}

#include "../flatten/closeup_align.h"
#include "ink_closeup_projection.inc"
#include "ink_closeup_projection_test.inc"
#include "ink_closeup_physical_align.inc"
#include "ink_closeup_global_contact_test.inc"
#include "ink_closeup_align_plan.inc"

static int ic_reparameterize_dense(const MeshBinData *m, const SlimRefineDirection *directions,
                                    double *uv, SlimRefineStats *stats)
{
    Arena_T arena = Arena_new();
    double *seed = malloc(m->nv * 2 * sizeof *seed);
    double *conformal = NULL;
    float *guide = NULL;
    SlimRefineOpts options;
    int rc = -1;
    if (!arena || !seed || !m->nv || !m->nf || (!m->uv && !m->uv64)) goto out;
    for (size_t i = 0; i < 2 * m->nv; i++) seed[i] = MeshBin_uv(m, i);
    if (ic_reparameterize_mode) {
        conformal=malloc(2*m->nv*sizeof *conformal);guide=malloc(2*m->nv*sizeof *guide);
        if (!conformal || !guide) goto out;
        memcpy(conformal,seed,2*m->nv*sizeof *conformal);
        if (ic_conformal_seed(m,conformal)) { fprintf(stderr,"ERROR: conformal initialization failed\n");goto out; }
        for (size_t i=0;i<2*m->nv;i++) guide[i]=(float)conformal[i];
    }
    SlimRefine_defaults(&options);
    options.iterations = CLOSEUP_FLAT_ITERS;
    options.energy_tolerance = CLOSEUP_FLAT_TOL;
    options.strip_pack = 0;
    options.guard_boundary = 1;
    options.local_steps = 1;
    options.local_step_min_faces = 1;
    options.local_line_search = 24;
    options.local_relaxation = 4;
    options.local_relax_boundary = 1;
    options.linear_solve = ic_slim_linear;
    options.verbose = 1;
    options.directions = directions;
    options.direction_weight = directions ? CLOSEUP_FIBER_DIRECTION_WEIGHT : 0;
    options.direction_rigid_gauge = directions != NULL;
    options.guide_uv=guide;
    options.guide_warm_start=guide!=NULL;
    rc = SlimRefine_run_double(arena, m->verts, m->nv, m->faces, m->nf, seed, &options, uv, stats);
    if (!rc && (stats->reverted || stats->flips_after || stats->boundary_intersections_after ||
                !isfinite(stats->energy_after) || stats->energy_after > stats->energy_before + 1e-8)) rc = -1;
    if (!rc) {
        /* SLIM normalizes each component's orientation. Restore the input's
         * reading handedness component-wise before a subsequent rigid turn. */
        int32_t *parent = ARENA_ALLOC(arena, m->nv * sizeof *parent);
        int8_t *sign = ARENA_CALLOC(arena, m->nv, sizeof *sign);
        for (size_t v = 0; v < m->nv; v++) parent[v] = (int32_t)v;
        for (size_t f = 0; f < m->nf; f++) for (int k = 1; k < 3; k++) {
            int32_t a = m->faces[3*f], b = m->faces[3*f+(size_t)k];
            while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
            while (parent[b] != b) { parent[b] = parent[parent[b]]; b = parent[b]; }
            if (a != b) parent[b] = a;
        }
        for (size_t f = 0; f < m->nf; f++) {
            size_t a = (size_t)m->faces[3*f], b = (size_t)m->faces[3*f+1], c = (size_t)m->faces[3*f+2];
            int32_t p = (int32_t)a;
            while (parent[p] != p) p = parent[p];
            double det = (MeshBin_uv(m,2*b)-MeshBin_uv(m,2*a))*(MeshBin_uv(m,2*c+1)-MeshBin_uv(m,2*a+1)) -
                         (MeshBin_uv(m,2*b+1)-MeshBin_uv(m,2*a+1))*(MeshBin_uv(m,2*c)-MeshBin_uv(m,2*a));
            sign[p] = det < 0 ? -1 : 1;
        }
        for (size_t v = 0; v < m->nv; v++) {
            int32_t p = (int32_t)v;
            while (parent[p] != p) p = parent[p];
            if (sign[p] < 0) uv[2*v+1] = -uv[2*v+1];
        }
        double lo[2] = {INFINITY, INFINITY};
        for (size_t i = 0; i < 2*m->nv; i++) {
            if (!isfinite(uv[i])) { rc = -1; break; }
            lo[i%2] = fmin(lo[i%2], uv[i]);
        }
        for (size_t i = 0; !rc && i < 2*m->nv; i++) uv[i] += CLOSEUP_RASTER_MARGIN_PX - lo[i%2];
    }
out:
    free(seed);free(conformal);free(guide);Arena_dispose(&arena);
    return rc;
}

static int ic_reparameterize_uv_directed(const MeshBinData *m, const SlimRefineDirection *directions,
                                          double *uv, SlimRefineStats *stats)
{
    /* Completion retains grid samples with no incident face. They remain in
     * the handoff, but have no metric and must not form empty SLIM components. */
    MeshBinData dense = {0};
    int32_t *map = malloc(m->nv * sizeof *map);
    double *result = malloc(2*m->nv*sizeof *result);
    int rc = -1;
    dense.verts = malloc(3*m->nv*sizeof *dense.verts);
    dense.uv64 = malloc(2*m->nv*sizeof *dense.uv64);
    dense.faces = malloc(3*m->nf*sizeof *dense.faces);
    if (!map || !result || !dense.verts || !dense.uv64 || !dense.faces) goto out;
    for (size_t v = 0; v < m->nv; v++) map[v] = -1;
    for (size_t i = 0; i < 3*m->nf; i++) map[m->faces[i]] = 0;
    for (size_t v = 0; v < m->nv; v++) if (map[v] >= 0) {
        size_t d = dense.nv++;
        map[v] = (int32_t)d;
        memcpy(dense.verts+3*d,m->verts+3*v,3*sizeof(float));
        dense.uv64[2*d] = MeshBin_uv(m,2*v); dense.uv64[2*d+1] = MeshBin_uv(m,2*v+1);
    }
    dense.nf = m->nf;
    for (size_t i = 0; i < 3*m->nf; i++) dense.faces[i] = map[m->faces[i]];
    fprintf(stderr,"[ink_closeup] SLIM: %zu incident vertices; %zu isolated samples retained\n",dense.nv,m->nv-dense.nv);
    rc = ic_reparameterize_dense(&dense,directions,result,stats);
    if (!rc) for (size_t v = 0; v < m->nv; v++) for (size_t k = 0; k < 2; k++)
        uv[2*v+k] = map[v] >= 0 ? result[2*(size_t)map[v]+k] : MeshBin_uv(m,2*v+k);
out:
    free(map); free(result); MeshBin_dispose(&dense);
    return rc;
}

static int ic_reparameterize_uv(const MeshBinData *m, double *uv, SlimRefineStats *stats)
{ return ic_reparameterize_uv_directed(m,NULL,uv,stats); }

static int ic_fiber_qualified(const SlimRefineStats *s)
{
    double rms = sqrt(s->direction_angle_after)*180/3.14159265358979323846;
    return s->direction_faces > 0 && s->direction_area_fraction >= CLOSEUP_FIBER_MIN_AREA &&
           isfinite(rms) && rms <= CLOSEUP_FIBER_RMS_DEG &&
           isfinite(s->objective_after) && s->objective_after <= s->objective_before+1e-8;
}

static void ic_growth_report(FILE *report, const SlimRefineStats *stats)
{
    fprintf(report," \"feasible_growth\": {\"component_steps\": %zu, \"vertex_steps\": %zu, "
            "\"alpha_mean\": %.12g, \"alpha_max\": %.12g},\n",
            stats->growth_component_steps_accepted,stats->growth_vertex_steps,
            stats->growth_vertex_alpha_mean,stats->growth_vertex_alpha_max);
}

/* Emits one JSON member plus a trailing comma for the surrounding report. */
static void ic_fiber_report(FILE *fp, const SlimRefineStats *s, const CloseupFiberStats *m, const IcConfig *c)
{
    char field_json[4800] = "null", identity_json[5000] = "null", identity_path[2200];
    char angle_before[64], angle_after[64], rms[64];
    ic_align_json_number(angle_before,s->direction_angle_before);
    ic_align_json_number(angle_after,s->direction_angle_after);
    ic_align_json_number(rms,sqrt(s->direction_angle_after)*57.2957795130823208768);
    if (c->fiber_volume[0]) {
        ic_json_str(field_json,sizeof field_json,c->fiber_volume);
        snprintf(identity_path,sizeof identity_path,"%s/inference_identity.json",c->fiber_volume);
        ic_json_str(identity_json,sizeof identity_json,identity_path);
    }
    fprintf(fp," \"fiber_alignment\": {\"reference_kind\": \"%s\", \"field_root\": %s, \"field_identity\": %s,"
            " \"qualified\": %s, \"direction_faces\": %zu, \"direction_area_fraction\": %.12g,"
            " \"direction_weight\": %.12g, \"rigid_gauge_accepted\": %s, \"rigid_gauge_radians\": %.12g,"
            " \"direction_energy_before\": %.12g, \"direction_energy_after\": %.12g,"
            " \"direction_angle_before\": %s, \"direction_angle_after\": %s,"
            " \"objective_before\": %.12g, \"objective_after\": %.12g,"
            " \"rms_degrees\": %s, \"rms_limit_degrees\": %.12g, \"minimum_area_fraction\": %.12g,"
            " \"measurement\": {\"patches\": %zu, \"measured_patches\": %zu, \"faces\": %zu,"
            " \"samples\": %zu, \"area\": %.12g, \"measured_area\": %.12g,"
            " \"chunk_misses\": %zu, \"chunk_missing\": %zu, \"seconds\": %.3f}},\n",
            c->fiber_volume[0]?"inferred_3d_fiber_cross":"ct_physical_tangent_patches",field_json,identity_json,
            ic_fiber_qualified(s)?"true":"false",s->direction_faces,s->direction_area_fraction,
            CLOSEUP_FIBER_DIRECTION_WEIGHT,s->direction_gauge_accepted?"true":"false",s->direction_gauge_radians,
            s->direction_energy_before,s->direction_energy_after,angle_before,angle_after,
            s->objective_before,s->objective_after,rms,
            CLOSEUP_FIBER_RMS_DEG,CLOSEUP_FIBER_MIN_AREA,m->patches,m->measured_patches,m->faces,m->samples,
            m->area,m->measured_area,m->chunk_misses,m->chunk_missing,m->seconds);
}

static int ic_measure_directions(const IcConfig *c, const MeshBinData *m, const float *normal,
                                  const double *frame, SlimRefineDirection *directions, CloseupFiberStats *stats)
{
    if (c->fiber_volume[0]) {
        char path[2200]; snprintf(path,sizeof path,"%s/inference_identity.json",c->fiber_volume);
        FILE *identity = fopen(path,"rb");
        if (!identity) { fprintf(stderr,"ERROR: physical fiber volume has no inference identity\n"); return -1; }
        fclose(identity);
        return CloseupFiberVolume_measure(m->verts,normal,m->nv,m->faces,m->nf,c->fiber_volume,directions,stats);
    }
    return CloseupFiber_measure(m->verts,normal,frame,m->nv,m->faces,m->nf,c->raw,directions,stats);
}

static int ic_freeze_fiber_identity(const IcConfig *c, const char *target)
{
    char path[2800]; size_t size = 0;
    if (!c->fiber_volume[0]) return 0;
    snprintf(path,sizeof path,"%s/inference_identity.json",c->fiber_volume);
    void *identity = ic_read_file(path,0,&size);
    if (!identity || !size) { free(identity); return -1; }
    snprintf(path,sizeof path,"%s_surface_inference_identity.json",target);
    int rc = ic_write_file(path,identity,size);
    free(identity); return rc;
}

/* Read-only measurement probe; it produces no repaired surface or eligibility
 * claim. This makes CT support inspectable before the expensive global solve. */
static int ic_stage_measure_fibers(const IcConfig *c)
{
    MeshBinData m = {0}; CloseupFiberStats stats = {0};
    char path[2600], target[2600], json_source[4800], json_raw[4800], json_field[4800] = "null";
    float *normal = NULL; double *frame = NULL;
    SlimRefineDirection *directions = NULL;
    size_t got = 0; int rc = 1; FILE *fp = NULL;
    snprintf(path,sizeof path,"%s_surface.vmesh",c->surface_source);
    if (MeshBin_read_precise_malloc(path,&m)) goto done;
    snprintf(path,sizeof path,"%s_surface_nrm.f32",c->surface_source);
    normal = ic_read_file(path,3*m.nv*sizeof(float),&got);
    frame = malloc(2*m.nv*sizeof *frame); directions = calloc(m.nf,sizeof *directions);
    if (!normal || got != 3*m.nv*sizeof(float) || !frame || !directions) goto done;
    for (size_t k = 0; k < 2*m.nv; k++) frame[k] = MeshBin_uv(&m,k);
    if (ic_measure_directions(c,&m,normal,frame,directions,&stats)) goto done;
    snprintf(target,sizeof target,"%s/%s",c->out_dir,c->id);
    snprintf(path,sizeof path,"%s_surface_fiber_arms.f64",target); ves_ensure_parent_dir(path);
    if (ic_write_file(path,directions,m.nf*sizeof *directions)) goto done;
    if (ic_freeze_fiber_identity(c,target)) goto done;
    snprintf(path,sizeof path,"%s_fiber_measure.json",target); fp = fopen(path,"wb");
    if (!fp) goto done;
    ic_json_str(json_source,sizeof json_source,c->surface_source); ic_json_str(json_raw,sizeof json_raw,c->raw);
    if (c->fiber_volume[0]) ic_json_str(json_field,sizeof json_field,c->fiber_volume);
    fprintf(fp,"{\"schema\":\"closeup-physical-fiber-measure-v1\",\"surface_source\":%s,\"raw\":%s,"
            "\"reference_kind\":\"%s\",\"field_root\":%s,"
            "\"surface_faces\":%zu,\"direction_faces\":%zu,\"patches\":%zu,\"measured_patches\":%zu,"
            "\"samples\":%zu,\"area\":%.12g,\"measured_area\":%.12g,\"seconds\":%.3f}\n",
            json_source,json_raw,c->fiber_volume[0]?"inferred_3d_fiber_cross":"ct_physical_tangent_patches",json_field,
            m.nf,stats.faces,stats.patches,stats.measured_patches,stats.samples,
            stats.area,stats.measured_area,stats.seconds);
    rc = fclose(fp) != 0; fp = NULL;
done:
    if (fp) fclose(fp);
    free(normal); free(frame); free(directions); MeshBin_dispose(&m);
    return rc;
}

static int ic_stage_reparameterize(const IcConfig *c)
{
    char path[2600], target[2600], source_json[4800];
    MeshBinData m = {0};
    SlimRefineStats stats = {0};
    IcPhysicalAlignStats physical = {0};
    IcAlignPlanStats plan = {0};
    AsmContactStats contact = {0};
    CloseupFiberStats measured = {0};
    SlimRefineDirection *directions = NULL;
    double *uv = NULL;
    void *side[3] = {NULL,NULL,NULL};
    int32_t *source_map = NULL,*original_aliases=NULL,*output_aliases=NULL;
    size_t original_nv = 0, cones = 0, cut_edges = 0, fiber_cones = 0;
    const char *suffix[3] = {"_surface_nrm.f32", "_surface_attr.u8", "_surface_chart.i32"};
    size_t length[3], got = 0;
    FILE *report = NULL;
    int rc = 1;
    snprintf(path, sizeof path, "%s_surface.vmesh", c->surface_source);
    if (MeshBin_read_precise_malloc(path, &m)) goto out;
    original_nv = m.nv;
    if(ic_material_read_aliases(c->surface_source,original_nv,&original_aliases))goto out;
    length[0] = m.nv * 3 * sizeof(float); length[1] = m.nv; length[2] = m.nv * sizeof(int32_t);
    for (int k = 0; k < 3; k++) {
        snprintf(path, sizeof path, "%s%s", c->surface_source, suffix[k]);
        side[k] = ic_read_file(path, length[k], &got);
        if (!side[k] || got != length[k]) goto out;
    }
    if (c->fiber_constraints) {
        double *frame = malloc(2*m.nv*sizeof *frame);
        directions = calloc(m.nf,sizeof *directions);
        if (!frame || !directions) { free(frame); goto out; }
        for (size_t k = 0; k < 2*m.nv; k++) frame[k] = MeshBin_uv(&m,k);
        int measured_rc = ic_measure_directions(c,&m,side[0],frame,directions,&measured);
        free(frame); if (measured_rc) goto out;
    }
    if (directions) {
        if (ic_reparameterize_physical_plan(c,&m,side[1],directions,original_aliases,ic_reparameterize_mode==2,
                                            &uv,&source_map,&stats,&physical,&contact,
                                            &cones,&cut_edges,&fiber_cones,&plan)) goto out;
    } else {
        if (ic_reparameterize_mode==2 && ic_cone_cut_directed(&m,NULL,&source_map,&cones,&cut_edges,&fiber_cones)) goto out;
        uv = malloc(2*m.nv*sizeof *uv);
        if (!uv || ic_reparameterize_uv_directed(&m,NULL,uv,&stats)) goto out;
    }
    if(ic_material_compose_aliases(original_nv,original_aliases,m.nv,source_map,&output_aliases))goto out;
    if(!directions)(void)ic_material_measure(&m,uv,output_aliases,&contact);
    if (source_map) for (int k=0;k<3;k++) {
        size_t stride=length[k]/original_nv;
        uint8_t *mapped=malloc(m.nv*stride);
        if (!mapped) goto out;
        for (size_t v=0;v<m.nv;v++) memcpy(mapped+v*stride,(uint8_t *)side[k]+(size_t)source_map[v]*stride,stride);
        free(side[k]);side[k]=mapped;length[k]=m.nv*stride;
    }
    snprintf(target, sizeof target, "%s/%s", c->out_dir, c->id);
    snprintf(path, sizeof path, "%s_surface.vmesh", target);
    ves_ensure_parent_dir(path);
    if (MeshBin_write_uv64(path, m.verts, m.nv, m.faces, m.nf, uv)) goto out;
    if (ic_material_write_aliases(target,m.nv,output_aliases)) goto out;
    for (int k = 0; k < 3; k++) {
        snprintf(path, sizeof path, "%s%s", target, suffix[k]);
        if (ic_write_file(path, side[k], length[k])) goto out;
    }
    if (source_map) {
        snprintf(path,sizeof path,"%s_surface_source_vertex.i32",target);
        if (ic_write_file(path,source_map,m.nv*sizeof *source_map)) goto out;
    }
    if (directions) {
        snprintf(path,sizeof path,"%s_surface_fiber_arms.f64",target);
        if (ic_write_file(path,directions,m.nf*sizeof *directions)) goto out;
        if (ic_freeze_fiber_identity(c,target)) goto out;
    }
    snprintf(path, sizeof path, "%s_surface_report.json", target);
    report = fopen(path, "wb");
    if (!report) goto out;
    ic_json_str(source_json, sizeof source_json, c->surface_source);
    fprintf(report,"{\n");
    ic_material_report(report,&contact);
    ic_growth_report(report,&stats);
    if (directions) { ic_fiber_report(report,&stats,&measured,c); ic_physical_alignment_report(report,&physical); ic_alignment_plan_report(report,&plan); }
    fprintf(report, " \"schema\": \"ink-closeup-reparameterize-v1\",\n \"surface_source\": %s,\n"
            " \"method\": \"slim_full_mesh\",\n \"world_geometry_preserved\": true,\n"
            " \"source_vertices\": %zu, \"cone_terminals\": %zu, \"fiber_cone_terminals\": %zu, \"cut_edges\": %zu, \"initialization\": \"%s\",\n"
            " \"surface\": {\"verts\": %zu, \"faces\": %zu},\n"
            " \"iterations\": %zu, \"accepted_steps\": %zu,\n"
            " \"energy_before\": %.12g, \"energy_after\": %.12g,\n"
            " \"flips_after\": %zu, \"boundary_intersections_after\": %zu,\n \"seconds\": %.3f\n}\n",
            source_json, original_nv, cones, fiber_cones, cut_edges, directions ? "physical_field_and_original_seed_trials" : (ic_reparameterize_mode ? "conformal" : "source"),
            m.nv, m.nf, stats.iterations_run, stats.accepted_steps,
            stats.energy_before, stats.energy_after, stats.flips_after, stats.boundary_intersections_after, stats.total_seconds);
    rc = fclose(report) != 0 || !ic_material_qualified(&contact) ||
         (directions && !ic_fiber_qualified(&stats)); report = NULL;
out:
    if (report) fclose(report);
    if (rc) fprintf(stderr, "ERROR: close-up native proposal did not qualify; source retained\n");
    for (int k = 0; k < 3; k++) free(side[k]);
    free(source_map);free(original_aliases);free(output_aliases);
    free(directions);
    free(uv); MeshBin_dispose(&m);
    return rc;
}

static int ic_reparameterize_selftest(void)
{
    if (SeamCut_selftest()) return 1;
    if (ic_cone_feasibility_selftest()) return 1;
    float xyz[30];
    double seed[20], result[20];
    int32_t faces[24] = {0,1,4,0,4,3,1,2,5,1,5,4,3,4,7,3,7,6,4,5,8,4,8,7};
    for (int k = 0; k < 9; k++) {
        xyz[3*k] = (float)(k%3); xyz[3*k+1] = (float)(k/3); xyz[3*k+2] = 0;
        seed[2*k] = -1.3*(k%3); seed[2*k+1] = 0.8*(k/3);
    }
    xyz[27] = 99; xyz[28] = 88; xyz[29] = 77; seed[18] = 9; seed[19] = 8;
    MeshBinData m = {0}; SlimRefineStats stats = {0};
    m.nv = 10; m.nf = 8; m.verts = xyz; m.faces = faces; m.uv64 = seed;
    if (ic_reparameterize_uv(&m,result,&stats) || fabs(stats.energy_after-16.0) > 1e-6 || stats.flips_after) return 1;
    for (int f = 0; f < 8; f++) {
        int a=faces[3*f], b=faces[3*f+1], c=faces[3*f+2];
        double det=(result[2*b]-result[2*a])*(result[2*c+1]-result[2*a+1])-
                   (result[2*b+1]-result[2*a+1])*(result[2*c]-result[2*a]);
        if (!(det < 0)) return 1;
    }
    for (int k = 0; k < 9; k++) if (xyz[3*k] != k%3 || xyz[3*k+1] != k/3 || xyz[3*k+2] != 0) return 1;
    if (result[18] != 9 || result[19] != 8 || xyz[27] != 99 || xyz[28] != 88 || xyz[29] != 77) return 1;
    fprintf(stderr,"[closeup reparameterize selftest] PASS\n");
    return 0;
}

#include "ink_closeup_physical_align_test.inc"
