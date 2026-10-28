#include "closeup_fiber.h"
#include "../assemble/asm_fibre_frame.h"
#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct CfiberFace { int32_t component, row, column, face; } CfiberFace;
typedef struct CfiberPatch { size_t first, end; } CfiberPatch;

static int cfiber_compare(const void *a, const void *b)
{
    const CfiberFace *x = a, *y = b;
    if (x->component != y->component) return x->component < y->component ? -1 : 1;
    if (x->row != y->row) return x->row < y->row ? -1 : 1;
    if (x->column != y->column) return x->column < y->column ? -1 : 1;
    return 0;
}

static double cfiber_dot(const double a[3], const double b[3])
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

static int cfiber_unit(double a[3])
{
    double n = sqrt(cfiber_dot(a,a));
    if (!(n > 1e-12) || !isfinite(n)) return -1;
    for (int k = 0; k < 3; k++) a[k] /= n;
    return 0;
}

static void cfiber_cross(const double a[3], const double b[3], double c[3])
{
    for (int k = 0; k < 3; k++) c[k] = a[(k+1)%3]*b[(k+2)%3] - a[(k+2)%3]*b[(k+1)%3];
}

/* Invert the projection from this face to the sampling plane. Orthogonally
 * projecting a plane arm back onto a tilted face is its adjoint, not its
 * inverse, and introduces a second foreshortening of the same direction. */
static int cfiber_lift(const double arm[3], const double patch_normal[3],
                        const double face_normal[3], double out[3])
{
    double agreement=cfiber_dot(patch_normal,face_normal);
    if (!(agreement>0)) return -1;
    double offset=cfiber_dot(face_normal,arm)/agreement;
    for (int k=0;k<3;k++) out[k]=arm[k]-offset*patch_normal[k];
    return cfiber_unit(out);
}

/* A physical basis from a face's indexing frame. Used to choose a sampling
 * plane only; the CT measurement below determines the fiber cross. */
static double cfiber_face(const float *xyz, const double *uv, const int32_t *f,
                           double center[3], double h[3], double n[3])
{
    size_t a = (size_t)f[0], b = (size_t)f[1], c = (size_t)f[2];
    double e[3], g[3];
    double du1 = uv[2*b]-uv[2*a], dv1 = uv[2*b+1]-uv[2*a+1];
    double du2 = uv[2*c]-uv[2*a], dv2 = uv[2*c+1]-uv[2*a+1];
    double det = du1*dv2-dv1*du2;
    if (!(fabs(det) > 1e-12) || !isfinite(det)) return 0;
    for (int k = 0; k < 3; k++) {
        e[k] = (double)xyz[3*b+(size_t)k]-xyz[3*a+(size_t)k];
        g[k] = (double)xyz[3*c+(size_t)k]-xyz[3*a+(size_t)k];
        center[k] = ((double)xyz[3*a+(size_t)k]+xyz[3*b+(size_t)k]+xyz[3*c+(size_t)k])/3;
        h[k] = (e[k]*dv2-g[k]*dv1)/det;
    }
    cfiber_cross(e,g,n);
    double area = .5*sqrt(cfiber_dot(n,n));
    if (!(area > 0) || !isfinite(area) || cfiber_unit(h) || cfiber_unit(n)) return 0;
    if (det < 0) for (int k = 0; k < 3; k++) n[k] = -n[k];
    return area;
}

static size_t cfiber_find(const CfiberFace *order, const CfiberPatch *patches,
                           size_t np, const CfiberFace *key)
{
    size_t lo = 0, hi = np;
    while (lo < hi) {
        size_t mid = lo+(hi-lo)/2;
        if (cfiber_compare(order+patches[mid].first,key) < 0) lo = mid+1; else hi = mid;
    }
    return lo < np && !cfiber_compare(order+patches[lo].first,key) ? lo : SIZE_MAX;
}

static int cfiber_measure_volume(Arena_T arena, ZarrU8 *volume,
                                  const float *xyz, const float *normal,
                                  const double *frame_uv, size_t nv,
                                  const int32_t *faces, size_t nf,
                                  SlimRefineDirection *out, CloseupFiberStats *stats)
{
    double started = ves_clock_sec();
    CloseupFiberStats s = {0};
    int rc = -1;
    if (stats) memset(stats,0,sizeof *stats);
    if (!volume || !xyz || !frame_uv || !faces || !out || nv > INT_MAX || nf > INT_MAX) return -1;
    memset(out,0,nf*sizeof *out);
    if (!nv || !nf) return 0;
    size_t misses_before = volume->misses, missing_before = volume->missing;
    UnionFind uf = UF_new(arena,(int32_t)nv);
    for (size_t f = 0; f < nf; f++) {
        const int32_t *tri = faces+3*f;
        for (int k = 0; k < 3; k++) if (tri[k] < 0 || (size_t)tri[k] >= nv) goto done;
        uf_union(&uf,tri[0],tri[1]); uf_union(&uf,tri[0],tri[2]);
    }
    CfiberFace *order = ARENA_ALLOC(arena,nf*sizeof *order);
    CfiberPatch *patches = ARENA_ALLOC(arena,nf*sizeof *patches);
    for (size_t f = 0; f < nf; f++) {
        const int32_t *tri = faces+3*f;
        double p[3], h[3], n[3];
        double area = cfiber_face(xyz,frame_uv,tri,p,h,n);
        if (!(area > 0)) goto done;
        s.area += area;
        double u = 0, v = 0;
        for (int k = 0; k < 3; k++) { u += frame_uv[2*(size_t)tri[k]]/3; v += frame_uv[2*(size_t)tri[k]+1]/3; }
        u = floor(u/CLOSEUP_FIBER_PATCH_PX); v = floor(v/CLOSEUP_FIBER_PATCH_PX);
        if (!isfinite(u) || !isfinite(v) || u < INT_MIN+1.0 || u > INT_MAX-1.0 || v < INT_MIN+1.0 || v > INT_MAX-1.0) goto done;
        order[f].component = uf_find(&uf,tri[0]); order[f].column = (int32_t)u;
        order[f].row = (int32_t)v; order[f].face = (int32_t)f;
    }
    qsort(order,nf,sizeof *order,cfiber_compare);
    size_t np = 0;
    for (size_t first = 0; first < nf;) {
        size_t end = first+1;
        while (end < nf && !cfiber_compare(order+first,order+end)) end++;
        patches[np].first = first; patches[np++].end = end; first = end;
    }
    s.patches = np;
    for (size_t p = 0; p < np; p++) {
        Arena_Mark mark = Arena_save(arena);
        const CfiberFace *key = order+patches[p].first;
        double center[3] = {0}, h[3] = {0}, n[3] = {0}, v[3] = {0}, mass = 0;
        for (size_t i = patches[p].first; i < patches[p].end; i++) {
            double fc[3], fh[3], fn[3];
            double area = cfiber_face(xyz,frame_uv,faces+3*(size_t)order[i].face,fc,fh,fn);
            for (int k = 0; k < 3; k++) { center[k] += area*fc[k]; h[k] += area*fh[k]; n[k] += area*fn[k]; }
            mass += area;
        }
        if (!(mass > 0) || cfiber_unit(n)) { Arena_restore(arena,mark); continue; }
        double hn = cfiber_dot(h,n);
        for (int k = 0; k < 3; k++) { center[k] /= mass; h[k] -= hn*n[k]; }
        if (cfiber_unit(h)) { Arena_restore(arena,mark); continue; }
        cfiber_cross(n,h,v);
        size_t neighbors[9], count = 0, cap = 0;
        for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
            CfiberFace neighbor_key = *key; neighbor_key.row += dy; neighbor_key.column += dx;
            size_t at = cfiber_find(order,patches,np,&neighbor_key);
            if (at != SIZE_MAX) { neighbors[count++] = at; cap += patches[at].end-patches[at].first; }
        }
        if (cap > INT_MAX/3) goto done;
        AsmChart chart = {0};
        chart.xyz = ARENA_ALLOC(arena,9*cap*sizeof(float));
        chart.nrm = ARENA_ALLOC(arena,9*cap*sizeof(float));
        chart.uv = ARENA_ALLOC(arena,6*cap*sizeof(float));
        chart.faces = ARENA_ALLOC(arena,3*cap*sizeof(int32_t));
        size_t *face_order = ARENA_ALLOC(arena,cap*sizeof *face_order);
        size_t *face_samples = ARENA_CALLOC(arena,cap,sizeof *face_samples);
        uint8_t *supported = ARENA_CALLOC(arena,patches[p].end-patches[p].first,sizeof *supported);
        double u0 = ((double)key->column+.5)*CLOSEUP_FIBER_PATCH_PX;
        double v0 = ((double)key->row+.5)*CLOSEUP_FIBER_PATCH_PX;
        for (size_t q = 0; q < count; q++) {
            const CfiberPatch *neighbor_patch = patches+neighbors[q];
            for (size_t i = neighbor_patch->first; i < neighbor_patch->end; i++) {
                const int32_t *tri = faces+3*(size_t)order[i].face;
                double uc = 0, vc = 0, fc[3], fh[3], fn[3];
                for (int k = 0; k < 3; k++) { uc += frame_uv[2*(size_t)tri[k]]/3; vc += frame_uv[2*(size_t)tri[k]+1]/3; }
                if (fabs(uc-u0) > CLOSEUP_FIBER_PATCH_PX || fabs(vc-v0) > CLOSEUP_FIBER_PATCH_PX) continue;
                if (!(cfiber_face(xyz,frame_uv,tri,fc,fh,fn) > 0) || cfiber_dot(fn,n) < CLOSEUP_FIBER_PATCH_NORMAL) continue;
                for (int k = 0; k < 3; k++) {
                    size_t source = (size_t)tri[k], dst = chart.nv++;
                    double d[3];
                    for (int axis = 0; axis < 3; axis++) {
                        chart.xyz[3*dst+(size_t)axis] = xyz[3*source+(size_t)axis];
                        chart.nrm[3*dst+(size_t)axis] = normal ? normal[3*source+(size_t)axis] : (float)fn[axis];
                        d[axis] = (double)xyz[3*source+(size_t)axis]-center[axis];
                    }
                    chart.uv[2*dst] = (float)cfiber_dot(d,h); chart.uv[2*dst+1] = (float)cfiber_dot(d,v);
                    chart.faces[3*chart.nf+(size_t)k] = (int32_t)dst;
                }
                face_order[chart.nf++]=i;
            }
        }
        AsmFibreMoments moments = {0}; AsmFibreTurn turn = {0};
        AsmFibre_measure_patch(volume,arena,&chart,&moments,face_samples); AsmFibre_turn(&moments,&turn);
        for (size_t f=0;f<chart.nf;f++) if (face_samples[f] &&
            face_order[f]>=patches[p].first && face_order[f]<patches[p].end)
            supported[face_order[f]-patches[p].first]=1;
        s.samples += moments.px;
        if (turn.decided_by) {
            double angle = turn.cross_deg*3.14159265358979323846/180;
            double ca = cos(angle), sa = sin(angle);
            s.measured_patches++;
            for (size_t i = patches[p].first; i < patches[p].end; i++) {
                if (!supported[i-patches[p].first]) continue;
                size_t f = (size_t)order[i].face;
                double fc[3], fh[3], fn[3];
                double area = cfiber_face(xyz,frame_uv,faces+3*f,fc,fh,fn);
                if (cfiber_dot(fn,n) < CLOSEUP_FIBER_PATCH_NORMAL) continue;
                double ph[3],pv[3];
                for (int k=0;k<3;k++) { ph[k]=ca*h[k]+sa*v[k];pv[k]=-sa*h[k]+ca*v[k]; }
                if (cfiber_lift(ph,n,fn,out[f].h) || cfiber_lift(pv,n,fn,out[f].v)) continue;
                out[f].weight = fmin(turn.coherence,1.0);
                s.faces++; s.measured_area += area;
            }
        }
        Arena_restore(arena,mark);
        if ((p+1)%64 == 0 || p+1 == np) {
            fprintf(stderr,"[closeup fiber] patches %zu/%zu measured=%zu area=%.2f%% samples=%zu\n",
                    p+1,np,s.measured_patches,s.area > 0 ? 100*s.measured_area/s.area : 0,s.samples);
            fflush(stderr);
        }
    }
    s.chunk_misses = volume->misses-misses_before;
    s.chunk_missing = volume->missing-missing_before;
    rc = 0;
done:
    s.seconds = ves_clock_sec()-started;
    if (stats) *stats = s;
    if (rc) fprintf(stderr,"ERROR: close-up physical fiber measurement failed\n");
    return rc;
}

int CloseupFiber_measure_volume(ZarrU8 *volume, const float *xyz, const float *normal,
                                 const double *frame_uv, size_t nv,
                                 const int32_t *faces, size_t nf,
                                 SlimRefineDirection *out, CloseupFiberStats *stats)
{
    if (stats) memset(stats,0,sizeof *stats);
    if (!volume || !xyz || !frame_uv || !faces || !out || nv > INT_MAX || nf > INT_MAX) return -1;
    if (!nv || !nf) { memset(out,0,nf*sizeof *out); return 0; }
    /* Allocate before TRY so this handle is unchanged across longjmp. A failed
     * Arena_new has no owned scratch; every later exception disposes it here. */
    Arena_T arena = Arena_new();
    volatile int rc = -1;
    TRY {
        rc = cfiber_measure_volume(arena,volume,xyz,normal,frame_uv,nv,faces,nf,out,stats);
    } FINALLY {
        Arena_dispose(&arena);
    } END_TRY;
    return rc;
}

int CloseupFiber_measure(const float *xyz, const float *normal,
                          const double *frame_uv, size_t nv,
                          const int32_t *faces, size_t nf, const char *raw,
                          SlimRefineDirection *out, CloseupFiberStats *stats)
{
    double started = ves_clock_sec();
    if (stats) memset(stats,0,sizeof *stats);
    if (!xyz || !frame_uv || !faces || !raw || !out || nv > INT_MAX || nf > INT_MAX) return -1;
    memset(out,0,nf*sizeof *out);
    if (!nv || !nf) return 0;
    ZarrU8 *volume = calloc(1,sizeof *volume);
    if (!volume || ZarrU8_open(volume,raw)) {
        free(volume);
        fprintf(stderr,"ERROR: close-up physical fiber measurement failed\n");
        return -1;
    }
    volatile int rc = -1;
    TRY {
        rc = CloseupFiber_measure_volume(volume,xyz,normal,frame_uv,nv,faces,nf,out,stats);
    } FINALLY {
        if (stats) stats->seconds = ves_clock_sec()-started;
        ZarrU8_close(volume); free(volume);
    } END_TRY;
    return rc;
}

#include "closeup_fiber_test.inc"
