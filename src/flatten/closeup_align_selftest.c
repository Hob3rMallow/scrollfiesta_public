#include "closeup_align.h"
#include "sparse_solve.h"
#include "../common/mesh_manifold.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct CaTestSolve { int calls, mode; } CaTestSolve;

static int cat_solve(Arena_T arena, void *context, size_t nv, size_t ne,
                      const double *diagonal, const int32_t *edge_a,
                      const int32_t *edge_b, const double *blocks,
                      const double *rhs, double *result)
{
    CaTestSolve *test = context;
    if (test) {
        test->calls++;
        if (test->mode == 1) { (void)ARENA_ALLOC(arena,100000); RAISE(Arena_Failed); }
        if (test->mode == 2) return -1;
        if (test->mode == 3 && test->calls <= 8) {
            /* Poison only the eight optional projection solves: collapsed
             * triangles must never become a field seed passed to SLIM. */
            memset(result,0,2*nv*sizeof *result); return 0;
        }
    }
    if (nv > (size_t)INT_MAX/3 || ne > ((size_t)INT_MAX-3*nv)/4) return -1;
    Arena_Mark mark = Arena_save(arena);
    volatile int rc = -1;
    TRY {
        size_t count = 3*nv+4*ne, at = 0;
        int *row = ARENA_ALLOC(arena,count*sizeof *row), *col = ARENA_ALLOC(arena,count*sizeof *col);
        double *value = ARENA_ALLOC(arena,count*sizeof *value);
        for (size_t v = 0; v < nv; v++) {
            int x = (int)(2*v), y = x+1;
            row[at] = col[at] = x; value[at++] = diagonal[3*v];
            row[at] = y; col[at] = x; value[at++] = diagonal[3*v+1];
            row[at] = col[at] = y; value[at++] = diagonal[3*v+2];
        }
        for (size_t e = 0; e < ne; e++) for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) {
            int x = 2*edge_a[e]+a, y = 2*edge_b[e]+b;
            row[at] = x > y ? x : y; col[at] = x < y ? x : y;
            value[at++] = blocks[3*e+(size_t)(a+b)];
        }
        rc = Sparse_solve_sym((int)(2*nv),(int)count,row,col,value,rhs,result,SPARSE_SPD);
    } FINALLY { Arena_restore(arena,mark); } END_TRY;
    return rc;
}

static int cat_check(int ok, const char *name)
{
    fprintf(stderr,"[closeup align selftest] %s %s\n",ok ? "PASS" : "FAIL",name);
    return !ok;
}

static int cat_same_mark(Arena_Mark a, Arena_Mark b)
{ return a.avail == b.avail && a.chunk == b.chunk; }

int CloseupAlign_selftest(void)
{
    enum { NX = 9, NY = 7, NV = NX*NY, NF = 2*(NX-1)*(NY-1) };
    float xyz[NV*3], saved_xyz[NV*3];
    int32_t faces[NF*3], saved_faces[NF*3];
    double uv[NV*2], original[NV*2], out[NV*2], reflected[NV*2], sentinel[NV*2];
    SlimRefineDirection arms[NF], saved_arms[NF];
    CloseupAlignOpts opts; CloseupAlign_defaults(&opts);
    CloseupAlignStats stats = {0}, saved_stats;
    CaTestSolve solver = {0};
    Arena_T arena = Arena_new_sized(4096);
    int failures = 0;
    opts.slim.linear_solve = cat_solve; opts.slim.linear_context = &solver;
    opts.slim.threads = 2; opts.slim.iterations = 20; opts.slim.verbose = 0;
    for (int y = 0; y < NY; y++) for (int x = 0; x < NX; x++) {
        int v = y*NX+x;
        xyz[3*v] = (float)x; xyz[3*v+1] = (float)y; xyz[3*v+2] = 0;
        double u = 1.08*x+.12*sin(.5*y), w = .95*y+.11*sin(.4*x), a = .36;
        uv[2*v] = 137+cos(a)*u-sin(a)*w; uv[2*v+1] = -62+sin(a)*u+cos(a)*w;
    }
    for (int y = 0, f = 0; y < NY-1; y++) for (int x = 0; x < NX-1; x++) {
        int v = y*NX+x;
        faces[3*f] = v; faces[3*f+1] = v+1; faces[3*f+2] = v+NX+1; f++;
        faces[3*f] = v; faces[3*f+1] = v+NX+1; faces[3*f+2] = v+NX; f++;
    }
    memset(arms,0,sizeof arms);
    for (int f = 0; f < NF; f++) { arms[f].h[0] = 1; arms[f].v[1] = 1; arms[f].weight = 2.5; }
    memcpy(original,uv,sizeof uv); memcpy(saved_xyz,xyz,sizeof xyz);
    memcpy(saved_faces,faces,sizeof faces); memcpy(saved_arms,arms,sizeof arms);
    for (int k = 0; k < NV*2; k++) out[k] = sentinel[k] = -12345;
    Arena_Mark mark = Arena_save(arena);
    int rc = CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&opts,out,&stats);
    failures += cat_check(rc == CLOSEUP_ALIGN_PROPOSED && stats.target_met && stats.scale_nonregression &&
        stats.after.fiber_rms_degrees < 5 && !stats.after.inverted_faces,"warped grid improves under physical metric guards");
    failures += cat_check(!memcmp(xyz,saved_xyz,sizeof xyz) && !memcmp(faces,saved_faces,sizeof faces) &&
        !memcmp(uv,original,sizeof uv) && !memcmp(arms,saved_arms,sizeof arms) &&
        cat_same_mark(mark,Arena_save(arena)),"immutable inputs and borrowed arena restored");
    double before[2] = {0}, after[2] = {0};
    for (int f = 0; f < NF; f++) for (int k = 0; k < 3; k++) for (int a = 0; a < 2; a++) {
        before[a] += uv[2*faces[3*f+k]+a]/(NF*3.0);
        after[a] += out[2*faces[3*f+k]+a]/(NF*3.0);
    }
    int corners[4] = {0,NX-1,NX*(NY-1),NV-1}, moved = 0;
    for (int k = 0; k < 4; k++) {
        int v = corners[k]; moved += hypot(out[2*v]-uv[2*v],out[2*v+1]-uv[2*v+1]) > .01;
    }
    failures += cat_check(moved == 4 && hypot(after[0]-before[0],after[1]-before[1]) < 1e-10,
        "all four actual boundary corners free; physical-area centroid preserved");
    for (int v = 0; v < NV; v++) uv[2*v+1] = -uv[2*v+1];
    solver.calls = 0;
    rc = CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&opts,reflected,&stats);
    double mismatch = 0;
    for (int v = 0; v < NV; v++) mismatch = fmax(mismatch,
        hypot(out[2*v]-reflected[2*v],out[2*v+1]+reflected[2*v+1]));
    failures += cat_check(rc == 0 && mismatch < 1e-7,"mirrored input preserves handedness and equivalent solution");
    memcpy(uv,original,sizeof uv);
    for (int f = 0; f < NF; f++) arms[f].weight = 0;
    solver.calls = 0;
    rc = CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&opts,reflected,&stats);
    failures += cat_check(rc == CLOSEUP_ALIGN_NO_EVIDENCE && !solver.calls &&
        !memcmp(uv,reflected,sizeof uv) && isnan(stats.after.fiber_rms_degrees) && !stats.after.measured_faces,
        "zero evidence is exact unchanged abstention and unknown RMS");
    memcpy(arms,saved_arms,sizeof arms);
    solver.mode = 3; solver.calls = 0;
    rc = CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&opts,reflected,&stats);
    int only_old = 1;
    for (size_t k = 0; k < stats.trials; k++) if (stats.trial[k].seed != 0) only_old = 0;
    failures += cat_check(rc >= 0 && stats.field_attempted && !stats.field_feasible && stats.field_status < 0 && only_old,
        "noninjective optional seed refused; old-seed control retained");
    solver.mode = 2; solver.calls = 0; memcpy(reflected,sentinel,sizeof reflected);
    rc = CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&opts,reflected,&stats);
    failures += cat_check(rc == CLOSEUP_ALIGN_SOLVER_FAILED && !memcmp(reflected,sentinel,sizeof reflected) &&
        cat_same_mark(mark,Arena_save(arena)),"linear backend failure leaves UV unchanged and arena restored");
    solver.mode = 1; solver.calls = 0; memset(&stats,0x5a,sizeof stats); memcpy(&saved_stats,&stats,sizeof stats);
    volatile int caught = 0;
    TRY { (void)CloseupAlign_propose(arena,xyz,NV,faces,NF,uv,arms,NULL,&opts,reflected,&stats); }
    EXCEPT(Arena_Failed) { caught = 1; } END_TRY;
    failures += cat_check(caught && !memcmp(reflected,sentinel,sizeof reflected) &&
        !memcmp(&stats,&saved_stats,sizeof stats) && cat_same_mark(mark,Arena_save(arena)),
        "allocation exception restores borrowed arena and both outputs atomically");
    solver.mode = 0;
    {
        const float p[12] = {0,0,0, 1,0,0, 1,1,0, 0,1,0};
        const int32_t f[6] = {0,1,2, 0,2,3};
        const double shear[8] = {0,0, 1,0, 2,.001, 1,.001};
        CloseupAlignMetric m;
        int measured = CloseupAlign_measure(arena,p,4,f,2,shear,arms,NULL,1,&m);
        failures += cat_check(!measured && m.fiber_rms_degrees > 63 && m.fiber_rms_degrees < 64 &&
            fabs(m.confidence_mass-2.5) < 1e-12,"joint complementary axes reject collapsed-axis false alignment; confidence unclipped");
        double invalid[8] = {0,0, 1,0, 0,1, 1,1}, output[8] = {99,99,99,99,99,99,99,99};
        rc = CloseupAlign_propose(arena,p,4,f,2,invalid,arms,NULL,&opts,output,&stats);
        failures += cat_check(rc == CLOSEUP_ALIGN_INFEASIBLE && output[0] == 99 && output[7] == 99,
            "mixed orientation input is refused atomically");
    }
    {
        const float p[15] = {0,0,0, 1,0,0, 0,1,0, -1,0,0, 0,-1,0};
        const int32_t f[6] = {0,1,2, 0,3,4};
        const double seed[10] = {0,0, 1,0, 0,1, -1,0, 0,-1};
        double output[10]; for (int k = 0; k < 10; k++) output[k] = 99;
        rc = CloseupAlign_propose(arena,p,5,f,2,seed,arms,NULL,&opts,output,&stats);
        failures += cat_check(rc == CLOSEUP_ALIGN_TOPOLOGY && !stats.input_manifold_disk &&
            output[0] == 99 && output[9] == 99,"pinched vertex refuses despite positive faces, chi=1 and no edge multiplicity");
        const float ring[24] = {-2,-2,0, 2,-2,0, 2,2,0, -2,2,0, -1,-1,0, 1,-1,0, 1,1,0, -1,1,0};
        const int32_t ring_faces[24] = {0,1,5,0,5,4,1,2,6,1,6,5,2,3,7,2,7,6,3,0,4,3,4,7};
        double ring_uv[16], ring_original[16], ring_output[16];
        for (int v = 0; v < 8; v++) {
            double x = 1.08*ring[3*v], y = .96*ring[3*v+1];
            ring_uv[2*v] = 100+cos(.32)*x-sin(.32)*y;
            ring_uv[2*v+1] = 200+sin(.32)*x+cos(.32)*y;
            ring_output[2*v] = ring_output[2*v+1] = 99;
        }
        memcpy(ring_original,ring_uv,sizeof ring_uv);
        rc = CloseupAlign_propose(arena,ring,8,ring_faces,8,ring_uv,arms,NULL,&opts,ring_output,&stats);
        int boundary_moved = 0;
        for (int v = 0; v < 8; v++) boundary_moved +=
            hypot(ring_uv[2*v]-ring_output[2*v],ring_uv[2*v+1]-ring_output[2*v+1]) > .01;
        failures += cat_check(rc == CLOSEUP_ALIGN_PROPOSED && stats.input_topology_supported &&
            !stats.input_manifold_disk && stats.boundary_loops == 2 && stats.target_met &&
            !stats.after.inverted_faces && boundary_moved == 8 && !memcmp(ring_uv,ring_original,sizeof ring_uv),
            "valid annulus aligns with both boundary loops free, no cuts and unchanged source");
        for (int v = 0; v < 8; v++) ring_uv[2*v+1] = -ring_uv[2*v+1];
        rc = CloseupAlign_propose(arena,ring,8,ring_faces,8,ring_uv,arms,NULL,&opts,ring_output,&stats);
        failures += cat_check(rc == CLOSEUP_ALIGN_PROPOSED && stats.boundary_loops == 2 && stats.target_met,
            "mirrored annulus preserves supported boundary orientation");
    }
    {
        /* Conflicting measured arms require strain on an already isometric
         * square. Every native trial rolls back safely; that is not a solver
         * failure and must not abort a whole multi-component close-up. */
        const float x[12] = {0,0,0, 1,0,0, 1,1,0, 0,1,0};
        const int32_t f[6] = {0,1,2, 0,2,3};
        const double seed[8] = {0,0, 1,0, 1,1, 0,1};
        const uint8_t protected_faces[2] = {1,1};
        SlimRefineDirection cross[2] = {{0}};
        double result[8]; CloseupAlignOpts control = opts;
        control.field_initialize = 0; control.slim.iterations = 60;
        for (int k = 0; k < 2; k++) {
            double angle = k ? -.45 : .45;
            cross[k].h[0] = cos(angle); cross[k].h[1] = sin(angle);
            cross[k].v[0] = -sin(angle); cross[k].v[1] = cos(angle);
            cross[k].weight = 1;
        }
        solver.mode = 0; solver.calls = 0;
        rc = CloseupAlign_propose(arena,x,4,f,2,seed,cross,protected_faces,&control,result,&stats);
        int rollbacks = stats.trials == 4;
        for (size_t k = 0; k < stats.trials; k++)
            rollbacks &= stats.trial[k].solver_status == 1 && stats.trial[k].feasible;
        failures += cat_check(rc == CLOSEUP_ALIGN_UNCHANGED && rollbacks && stats.input_feasible &&
            !memcmp(seed,result,sizeof seed) && stats.scale_nonregression && !stats.target_met &&
            stats.before.sd_mean == stats.after.sd_mean && !stats.after.inverted_faces &&
            cat_same_mark(mark,Arena_save(arena)),
            "certified metric rollbacks retain exact input and try every lower weight");
    }
    failures += MeshManifold_selftest(arena);
    failures += CloseupAlign_rigid_selftest();
    Arena_dispose(&arena);
    fprintf(stderr,"[closeup align selftest] %s: %d failure(s)\n",failures ? "FAIL" : "PASS",failures);
    return failures;
}
