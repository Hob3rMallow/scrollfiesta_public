/* asm_emit.c -- sheet.vmesh (+uv, sidecars), extras, and the RAW bake. */
#include "asm_emit.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/mesh_bin.h"
#include "../common/ves_platform.h"
#include "asm_metric.h"
#include "asm_audit.h"

void AsmEmit_default_opts(AsmEmitOpts *o)
{
    o->raw_source = NULL;
    o->bake_exe = "build/Release/obj_bake_raw.exe";
    o->normal_range = 4.0;
    o->raster_du = 1.0;
    o->raster_dv = 1.0;
    o->level_v_to_z = 1;
}

static void ae_pose(const AsmChart *c, double u, double v, double *gx, double *gy)
{
    if (c->flags & ASM_CHART_MIRROR) u = -u;
    double ct = cos(c->pose_theta), st = sin(c->pose_theta);
    *gx = ct * u - st * v + c->pose_x;
    *gy = st * u + ct * v + c->pose_y;
}

/* Encoding is part of the metric contract. Check the exact array which the
 * writer will publish, after both rigid transforms and the selected encoding. */
static int ae_encode_chart(Arena_T arena, const AsmChart *c,
                            double gauge_c, double gauge_s, double off_u, double off_v,
                            void *uv, int precise, AsmMetricStats *before, AsmMetricStats *after)
{
    Arena_Mark mark = Arena_save(arena);
    AsmMetricFace *ref = ARENA_ALLOC(arena,c->nf*sizeof *ref);
    memset(before,0,sizeof *before); memset(after,0,sizeof *after);
    for (size_t k = 0; k < c->nv; k++) {
        double p[2]; AsmChart_point(c,k,p); double gx = p[0], gy = p[1];
        double u = gauge_c*gx-gauge_s*gy+off_u, v = gauge_s*gx+gauge_c*gy+off_v;
        if (precise) { ((double *)uv)[k*2] = u; ((double *)uv)[k*2+1] = v; }
        else { ((float *)uv)[k*2] = (float)u; ((float *)uv)[k*2+1] = (float)v; }
    }
    int ok = AsmMetric_prepare(c,ref) == 0;
    if (ok) {
        if (c->placed_uv) AsmMetric_measure_encoded_uv64(c,ref,c->placed_uv,before);
        else AsmMetric_measure(c,ref,c->uv,before);
        if (precise) AsmMetric_measure_encoded_uv64(c,ref,(const double *)uv,after);
        else AsmMetric_measure_encoded(c,ref,(const float *)uv,after);
        ok = AsmMetric_preserved(before,after);
        if (ok && precise) ok = AsmMetric_encoded_faces_preserved(c,ref,(const double *)uv);
    }
    Arena_restore(arena,mark);
    return ok;
}

static int ae_write_i32(const char *path, const int32_t *v, size_t n)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;
    size_t w = fwrite(v, sizeof(int32_t), n, fp);
    fclose(fp);
    return w == n ? 0 : -1;
}

static int ae_write_obj(const char *path, const float *verts, const double *uv, size_t nv,
                        const int32_t *faces, size_t nf)
{
    FILE *fp = fopen(path, "wb");
    if (!fp) return -1;
    fprintf(fp, "# sheet_assemble: v z y x (pipeline order); scientific uv in the companion .vmesh\n");
    for (size_t i = 0; i < nv; i++) fprintf(fp, "v %.3f %.3f %.3f\n", verts[i*3], verts[i*3+1], verts[i*3+2]);
    for (size_t i = 0; i < nv; i++) fprintf(fp, "vt %.17g %.17g\n", uv[i*2], uv[i*2+1]);
    for (size_t f = 0; f < nf; f++)
        fprintf(fp, "f %d/%d %d/%d %d/%d\n", faces[f*3]+1, faces[f*3]+1, faces[f*3+1]+1, faces[f*3+1]+1, faces[f*3+2]+1, faces[f*3+2]+1);
    fclose(fp);
    return 0;
}

/* concatenate the charts selected by `select` (1 = primary sheet, 0 = extras) */
static int ae_write_set(AsmRun *run, const char *out_dir, const char *stem, int32_t primary, int want_primary,
                        double gauge_c, double gauge_s, double off_u, double off_v,
                        size_t *out_charts, size_t *out_nv, size_t *out_nf, double *out_area,
                        size_t *metric_checked, size_t *metric_regressions,
                        double *umin, double *umax, double *vmin, double *vmax)
{
    (void)primary;
    Arena_Mark mark = Arena_save(run->arena);
    size_t nv = 0, nf = 0, ncharts = 0;
    double area = 0.0;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || !c->placed) continue;
        if (AsmChart_registered(c) != want_primary) continue;
        nv += c->nv; nf += c->nf; ncharts++; area += c->area3d;
    }
    *out_charts = ncharts; *out_nv = nv; *out_nf = nf; *out_area = area;
    if (nv == 0) { Arena_restore(run->arena, mark); return 0; }
    float *verts = ARENA_ALLOC(run->arena, nv * 3 * sizeof(float));
    double *uv = ARENA_ALLOC(run->arena, nv * 2 * sizeof(double));
    int32_t *faces = ARENA_ALLOC(run->arena, nf * 3 * sizeof(int32_t));
    int32_t *side_comp = ARENA_ALLOC(run->arena, nv * sizeof(int32_t));
    int32_t *side_chart = ARENA_ALLOC(run->arena, nv * sizeof(int32_t));
    int32_t *side_cube = ARENA_ALLOC(run->arena, nv * sizeof(int32_t));
    size_t pv = 0, pf = 0;
    char path[2048];
    snprintf(path,sizeof path,"%s/%s_encoded_metric.csv",out_dir,stem);
    FILE *metric_log = fopen(path,"wb");
    if (!metric_log) { Arena_restore(run->arena,mark); return -1; }
    int metric_io = fputs("chart,faces,area,local10,local25,local_invalid,encoded10,encoded25,encoded_invalid,preserved\n",metric_log) >= 0;
    size_t regressions = 0;
    *umin = *vmin = 1e300; *umax = *vmax = -1e300;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || !c->placed) continue;
        if (AsmChart_registered(c) != want_primary) continue;
        AsmMetricStats before, after;
        int preserved = ae_encode_chart(run->arena,c,gauge_c,gauge_s,off_u,off_v,uv+pv*2,1,&before,&after);
        (*metric_checked)++;
        if (!preserved) { regressions++; (*metric_regressions)++; }
        double area_scale = before.area > 0 ? 1/before.area : 0;
        if (fprintf(metric_log,"%d,%zu,%.17g,%.17g,%.17g,%zu,%.17g,%.17g,%zu,%d\n",c->id,c->nf,
            before.area,before.within10*area_scale,before.within25*area_scale,before.invalid,
            after.within10*area_scale,after.within25*area_scale,after.invalid,preserved) < 0) metric_io = 0;
        for (size_t k = 0; k < c->nv; k++) {
            memcpy(&verts[(pv+k)*3], &c->xyz[k*3], 3 * sizeof(float));
            double u = uv[(pv+k)*2], v = uv[(pv+k)*2+1];
            if (u < *umin) *umin = u; if (u > *umax) *umax = u;
            if (v < *vmin) *vmin = v; if (v > *vmax) *vmax = v;
            side_comp[pv+k] = c->component; side_chart[pv+k] = c->id; side_cube[pv+k] = c->cube;
        }
        for (size_t f = 0; f < c->nf; f++)
            for (int e = 0; e < 3; e++) faces[(pf+f)*3 + (size_t)e] = c->faces[f*3 + (size_t)e] + (int32_t)pv;
        pv += c->nv; pf += c->nf;
    }
    if (fclose(metric_log) != 0) metric_io = 0;
    if (regressions || !metric_io) {
        fprintf(stderr,"[assemble emit] %s refused: %zu encoded chart metric regressions; diagnostic %s\n",stem,regressions,path);
        Arena_restore(run->arena,mark);
        return -1;
    }
    snprintf(path, sizeof path, "%s/%s.vmesh", out_dir, stem);
    int rc = MeshBin_write_uv64(path, verts, nv, faces, nf, uv);
    snprintf(path, sizeof path, "%s/%s.obj", out_dir, stem);
    rc |= ae_write_obj(path, verts, uv, nv, faces, nf);
    snprintf(path, sizeof path, "%s/%s_component.i32", out_dir, stem);
    rc |= ae_write_i32(path, side_comp, nv);
    snprintf(path, sizeof path, "%s/%s_material_identity.i32", out_dir, stem);
    rc |= ae_write_i32(path, side_comp, nv);
    snprintf(path, sizeof path, "%s/%s_chart.i32", out_dir, stem);
    rc |= ae_write_i32(path, side_chart, nv);
    snprintf(path, sizeof path, "%s/%s_cube.i32", out_dir, stem);
    rc |= ae_write_i32(path, side_cube, nv);
    if (!rc) rc = AsmAudit_encoded(run->arena,out_dir,stem,verts,uv,nv,faces,nf,side_chart,side_comp,side_cube);
    Arena_restore(run->arena, mark);
    return rc;
}

int AsmEmit_run(AsmRun *run, const AsmEmitOpts *o, const char *out_dir, AsmEmitStats *st)
{
    memset(st, 0, sizeof *st);
    st->bake_rc = -1;
    double t0 = ves_clock_sec();
    /* primary = the component with the largest area */
    int32_t ncomp = 0;
    for (size_t i = 0; i < run->n_charts; i++) if (AsmChart_registered(run->charts+i) && run->charts[i].component + 1 > ncomp) ncomp = run->charts[i].component + 1;
    if (ncomp == 0) return -1;
    double *area = ARENA_CALLOC(run->arena, (size_t)ncomp, sizeof(double));
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (AsmChart_registered(c) && c->component >= 0) area[c->component] += c->area3d;
    }
    int32_t primary = 0;
    for (int32_t c = 1; c < ncomp; c++) if (area[c] > area[primary]) primary = c;

    /* gauge: rotate so +v follows the world-z gradient of the primary sheet */
    double gc = 1.0, gs = 0.0, theta = 0.0;
    if (o->level_v_to_z) {
        /* least squares z ~ a*gx + b*gy + c over primary vertices */
        double sxx = 0, sxy = 0, syy = 0, sx = 0, sy = 0, sxz = 0, syz = 0, sz = 0, n = 0;
        for (size_t i = 0; i < run->n_charts; i++) {
            const AsmChart *c = &run->charts[i];
            if (!AsmChart_in_layout(c) || !c->placed || c->component != primary) continue;
            for (size_t k = 0; k < c->nv; k += 4) {
                double gx, gy;
                double point[2]; AsmChart_point(c,k,point); gx = point[0]; gy = point[1];
                double z = c->xyz[k*3];
                sxx += gx*gx; sxy += gx*gy; syy += gy*gy; sx += gx; sy += gy; sxz += gx*z; syz += gy*z; sz += z; n += 1;
            }
        }
        if (n > 10) {
            /* solve the 3x3 normal equations [sxx sxy sx; sxy syy sy; sx sy n] [a b c] = [sxz syz sz] */
            double A[9] = { sxx, sxy, sx, sxy, syy, sy, sx, sy, n }, B[3] = { sxz, syz, sz }, X[3] = { 0, 0, 0 };
            /* Gaussian elimination */
            for (int col = 0; col < 3; col++) {
                int piv = col;
                for (int r = col + 1; r < 3; r++) if (fabs(A[r*3+col]) > fabs(A[piv*3+col])) piv = r;
                for (int k = 0; k < 3; k++) { double t = A[col*3+k]; A[col*3+k] = A[piv*3+k]; A[piv*3+k] = t; }
                { double t = B[col]; B[col] = B[piv]; B[piv] = t; }
                if (fabs(A[col*3+col]) < 1e-12) continue;
                for (int r = 0; r < 3; r++) {
                    if (r == col) continue;
                    double f = A[r*3+col] / A[col*3+col];
                    for (int k = 0; k < 3; k++) A[r*3+k] -= f * A[col*3+k];
                    B[r] -= f * B[col];
                }
            }
            for (int k = 0; k < 3; k++) X[k] = fabs(A[k*3+k]) > 1e-12 ? B[k] / A[k*3+k] : 0.0;
            double ga = X[0], gb = X[1];
            if (hypot(ga, gb) > 1e-9) {
                /* rotate the gradient direction onto +v (angle pi/2) */
                theta = 0.5 * 3.14159265358979323846 - atan2(gb, ga);
                gc = cos(theta); gs = sin(theta);
            }
        }
    }
    st->gauge_theta = theta;
    /* offsets so the primary sheet starts at (1, 1) */
    double umin = 1e300, umax = -1e300, vmin = 1e300, vmax = -1e300;
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = &run->charts[i];
        if (!AsmChart_in_layout(c) || !c->placed || c->component != primary) continue;
        for (size_t k = 0; k < c->nv; k++) {
            double gx, gy;
            double point[2]; AsmChart_point(c,k,point); gx = point[0]; gy = point[1];
            double u = gc * gx - gs * gy, v = gs * gx + gc * gy;
            if (u < umin) umin = u; if (u > umax) umax = u;
            if (v < vmin) vmin = v; if (v > vmax) vmax = v;
        }
    }
    double off_u = 1.0 - umin, off_v = 1.0 - vmin;
    double a, b, c2, d;
    if (ae_write_set(run, out_dir, "sheet", primary, 1, gc, gs, off_u, off_v,
                     &st->sheet_charts, &st->sheet_verts, &st->sheet_faces, &st->sheet_area,
                     &st->metric_checked_charts,&st->metric_regressions,&a,&b,&c2,&d) != 0) return -1;
    st->u_span = b - a; st->v_span = d - c2;
    if (ae_write_set(run, out_dir, "sheet_extras", primary, 0, 1.0, 0.0, 0.0, 0.0,
                     &st->extras_charts, &st->extras_verts, &st->extras_faces, &st->extras_area,
                     &st->metric_checked_charts,&st->metric_regressions,&a,&b,&c2,&d) != 0) return -1;

    /* bake */
    if (o->raw_source && o->raw_source[0] && st->sheet_verts > 0) {
        char mesh[2048], bake_dir[2048], log[2048], du[64], dv[64], nr[64];
        snprintf(mesh, sizeof mesh, "%s/sheet.vmesh", out_dir);
        snprintf(bake_dir, sizeof bake_dir, "%s/bake", out_dir);
        snprintf(log, sizeof log, "%s/bake.log", out_dir);
        snprintf(du, sizeof du, "%.17g", o->raster_du);
        snprintf(dv, sizeof dv, "%.17g", o->raster_dv);
        snprintf(nr, sizeof nr, "%.17g", o->normal_range);
        ves_mkdir(bake_dir);
        const char *argv[] = { o->bake_exe, mesh, o->raw_source, bake_dir, "--id", "sheet",
                               "--raster-du", du, "--raster-dv", dv, "--raster-auto",
                               "--normal-range", nr, "--no-diag", "--require-complete-raw", NULL };
        st->bake_rc = ves_run_subprocess_logged(o->bake_exe, argv, 0.0, log);
        if (st->bake_rc != 0) {
            st->sec = ves_clock_sec() - t0;
            fprintf(stderr,"[assemble emit] requested RAW bake failed (rc %d); see %s\n",st->bake_rc,log);
            return -1;
        }
    }
    st->sec = ves_clock_sec() - t0;
    return 0;
}

static int ae_bake_failure_selftest(void)
{
    Arena_T arena = Arena_new();
    float xyz[] = {0,0,0, 1,0,0, 0,1,0}, uv[] = {0,0, 1,0, 0,1};
    int32_t faces[] = {0,1,2};
    AsmChart chart = {0}; AsmRun run = {0}; AsmEmitOpts opts; AsmEmitStats st;
    chart.xyz = xyz; chart.uv = uv; chart.faces = faces;
    chart.nv = 3; chart.nf = 1; chart.area3d = .5; chart.placed = 1; chart.placement_state = ASM_PLACE_ROOT;
    run.arena = arena; run.charts = &chart; run.n_charts = 1;
    AsmEmit_default_opts(&opts); opts.level_v_to_z = 0;
    char dir[256], path[320], missing[320];
    snprintf(dir,sizeof dir,"output/asm_emit_bake_selftest_%d",ves_getpid());
    snprintf(path,sizeof path,"%s/sheet.vmesh",dir);
    snprintf(missing,sizeof missing,"%s/baker_does_not_exist",dir);
    int fails = ves_ensure_parent_dir(path) != 0;
    if (AsmEmit_run(&run,&opts,dir,&st) != 0 || st.bake_rc != -1 || st.sheet_faces != 1) {
        fprintf(stderr,"  assembly bake selftest FAIL: explicitly skipped bake\n"); fails++;
    }
    {
        double actual_uv[]={1,1,2,1,1,2};int32_t labels[]={0,0,0};
        if(AsmAudit_encoded(arena,dir,"sheet",xyz,actual_uv,3,faces,1,labels,labels,labels))fails++;
        char sidecar[320];snprintf(sidecar,sizeof sidecar,"%s/sheet_chart.i32",dir);FILE *bad=fopen(sidecar,"r+b");
        if(!bad)fails++;else{int32_t wrong=13;fwrite(&wrong,sizeof wrong,1,bad);fclose(bad);
            if(!AsmAudit_encoded(arena,dir,"sheet",xyz,actual_uv,3,faces,1,labels,labels,labels))fails++;}
        /* A registered chart with a distinct material label belongs in the
         * sheet; an unregistered local pose remains explicitly in extras. */
        AsmChart charts[2]={chart,chart}; charts[1].id=1;charts[1].component=1;charts[1].pose_x=3;
        run.charts=charts;run.n_charts=2;
        if(AsmEmit_run(&run,&opts,dir,&st) || st.sheet_charts!=2 || st.extras_charts)fails++;
        charts[1].placement_state=ASM_PLACE_NONE;
        if(AsmEmit_run(&run,&opts,dir,&st) || st.sheet_charts!=1 || st.extras_charts!=1)fails++;
        run.charts=&chart;run.n_charts=1;
    }
    opts.raw_source = "missing_raw_for_bake_failure_control";
    opts.bake_exe = missing;
    if (AsmEmit_run(&run,&opts,dir,&st) == 0 || st.bake_rc == 0 || st.sheet_faces != 1) {
        fprintf(stderr,"  assembly bake selftest FAIL: failed child accepted as successful emission\n"); fails++;
    }
    Arena_dispose(&arena);
    fprintf(stderr,"  assembly bake failure propagation: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails;
}

int AsmEmit_selftest(void)
{
    Arena_T arena = Arena_new();
    float xyz[] = {0,0,0,1,0,0,0,1,0}, uv[] = {0,0,1,0,0,1}, encoded[6];
    int32_t face[] = {0,1,2}; AsmChart c = {0}; AsmMetricStats before, after;
    c.nv = 3; c.nf = 1; c.xyz = xyz; c.uv = uv; c.faces = face;
    int fails = ae_bake_failure_selftest();
    for (int mirror = 0; mirror < 2; mirror++) {
        c.flags = mirror ? ASM_CHART_MIRROR : 0; c.pose_theta = .713;
        c.pose_x = 341.2; c.pose_y = -42.3;
        if (!ae_encode_chart(arena,&c,cos(.217),sin(.217),5245.0003,1,encoded,0,&before,&after)) fails++;
        if (before.invalid || after.invalid) fails++;
    }
    c.flags = 0; c.pose_x = c.pose_y = c.pose_theta = 0;
    uv[2] = .90000005f;
    if (ae_encode_chart(arena,&c,1,0,5245.0003,1,encoded,0,&before,&after) ||
        before.within10 != before.area || after.within10 != 0 || after.invalid) fails++;
    double precise[6];
    if (!ae_encode_chart(arena,&c,1,0,5245.0003,1,precise,1,&before,&after) || after.invalid) fails++;
    /* A perfectly shaped thin original face collapses when absolute-frame
     * spacing exceeds its height. Rigid export may not silently accept it. */
    uv[2] = 1; xyz[7] = uv[5] = 1e-6f;
    if (ae_encode_chart(arena,&c,1,0,1,100,encoded,0,&before,&after) ||
        before.invalid || !after.invalid) fails++;
    if (!ae_encode_chart(arena,&c,1,0,1,100,precise,1,&before,&after) || after.invalid) fails++;
    if (uv[2] != 1 || uv[5] != 1e-6f || xyz[7] != 1e-6f || face[2] != 2) fails++;
    Arena_dispose(&arena);
    fprintf(stderr,"  assembly encoded export metric: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails;
}
