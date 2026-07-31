/* cvt_remesh_obj — CVT/RVD remesh an OBJ to a target vertex count.
 *   cvt_remesh_obj <in.obj> <out.obj> <target_verts> [--iters N] [--seed S]
 *   cvt_remesh_obj --selftest
 * M1: pure Lloyd. Mirrors qslim_obj's shape. */
#include "../remesh/cvt_remesh.h"
#include "../common/obj_io.h"
#include "../common/arena.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Graded-density selftest: a flat sheet inside a synthetic owned box; assert the
 * output edge length near the box planes ~ h_seam and deep interior ~ h_interior
 * (validates the sizing field + energy_exp=4 gradient preservation + derived N). */
static int graded_selftest(void) {
    Arena_T a = Arena_new();
    int G = 40, side = G+1;
    size_t nv = (size_t)side*side, nf = (size_t)G*G*2;
    float   *V = (float *)ARENA_ALLOC(a, (long)(nv*3*sizeof(float)));
    int32_t *F = (int32_t *)ARENA_ALLOC(a, (long)(nf*3*sizeof(int32_t)));
    for (int i=0;i<side;i++) for (int j=0;j<side;j++) {
        size_t idx=(size_t)i*side+j;
        V[idx*3+0]=20.0f;            /* z = 20, well inside the owned box */
        V[idx*3+1]=(float)i;         /* y in [0,40] */
        V[idx*3+2]=(float)j;         /* x in [0,40] */
    }
    size_t fi=0;
    for (int i=0;i<G;i++) for (int j=0;j<G;j++) {
        int32_t v00=(int32_t)((size_t)i*side+j), v01=(int32_t)((size_t)i*side+j+1);
        int32_t v10=(int32_t)((size_t)(i+1)*side+j), v11=(int32_t)((size_t)(i+1)*side+j+1);
        F[fi*3+0]=v00;F[fi*3+1]=v10;F[fi*3+2]=v11;fi++;
        F[fi*3+0]=v00;F[fi*3+1]=v11;F[fi*3+2]=v01;fi++;
    }
    CvtField fld;
    fld.lo=2.0; fld.hi=38.0; fld.h_seam=3.0; fld.h_interior=10.0;
    fld.band_lo=3.0; fld.band_hi=12.0; fld.energy_exp=4;
    CvtOpts o; CVT_default_opts(&o); o.n_iters=20; o.seed=7;
    float *OV; size_t onv; int32_t *OF; size_t onf;
    int rc = CVT_remesh(a, V, nv, F, nf, 100, &o, &fld, &OV, &onv, &OF, &onf);

    double near_sum=0.0, far_sum=0.0; int near_n=0, far_n=0;
    for (size_t t=0; t<onf; t++) for (int e=0;e<3;e++) {
        int32_t va=OF[t*3+e], vb=OF[t*3+(e+1)%3];
        double mz=(OV[va*3]+OV[vb*3])*0.5, my=(OV[va*3+1]+OV[vb*3+1])*0.5, mx=(OV[va*3+2]+OV[vb*3+2])*0.5;
        double dz=OV[va*3]-OV[vb*3], dy=OV[va*3+1]-OV[vb*3+1], dx=OV[va*3+2]-OV[vb*3+2];
        double len=sqrt(dz*dz+dy*dy+dx*dx);
        double d0=mz-fld.lo,d1=fld.hi-mz,d2=my-fld.lo,d3=fld.hi-my,d4=mx-fld.lo,d5=fld.hi-mx;
        double d=d0; if(d1<d)d=d1; if(d2<d)d=d2; if(d3<d)d=d3; if(d4<d)d=d4; if(d5<d)d=d5;
        if (d < fld.band_lo) { near_sum+=len; near_n++; }
        else if (d > fld.band_hi) { far_sum+=len; far_n++; }
    }
    double nm = near_n ? near_sum/near_n : 0.0, fm = far_n ? far_sum/far_n : 0.0;
    int ok = (rc==0 && onf>0 && near_n>0 && far_n>0
              && nm>2.0 && nm<5.0            /* ~h_seam=3 */
              && fm>6.5 && fm<13.5           /* ~h_interior=10 */
              && fm > nm*1.6);               /* a real gradient, not washed out */
    fprintf(stderr, "cvt graded selftest: rc=%d onv=%zu onf=%zu  near_edge=%.2f (h_seam 3)  "
            "far_edge=%.2f (h_int 10) -> %s\n", rc, onv, onf, nm, fm, ok?"OK":"FAIL");
    Arena_dispose(&a);
    return ok ? 0 : 1;
}

static int selftest(void) {
    if (CVT_projection_selftest() != 0) return 1;   /* KD-tree projection vs brute */
    if (graded_selftest() != 0) return 1;           /* graded sizing field */
    /* 16x16 grid on the unit square (z=0); remesh to 64 verts, expect a dual. */
    Arena_T a = Arena_new();
    int M = 16, side = M+1;
    size_t nv = (size_t)side*side, nf = (size_t)M*M*2;
    float   *V = (float *)ARENA_ALLOC(a, (long)(nv*3*sizeof(float)));
    int32_t *F = (int32_t *)ARENA_ALLOC(a, (long)(nf*3*sizeof(int32_t)));
    for (int i=0;i<side;i++) for (int j=0;j<side;j++) {
        size_t idx=(size_t)i*side+j;
        V[idx*3+0]=0.0f; V[idx*3+1]=(float)i/M; V[idx*3+2]=(float)j/M;
    }
    size_t fi=0;
    for (int i=0;i<M;i++) for (int j=0;j<M;j++) {
        int32_t v00=(int32_t)((size_t)i*side+j), v01=(int32_t)((size_t)i*side+j+1);
        int32_t v10=(int32_t)((size_t)(i+1)*side+j), v11=(int32_t)((size_t)(i+1)*side+j+1);
        F[fi*3+0]=v00;F[fi*3+1]=v10;F[fi*3+2]=v11;fi++;
        F[fi*3+0]=v00;F[fi*3+1]=v11;F[fi*3+2]=v01;fi++;
    }
    CvtOpts o; CVT_default_opts(&o); o.n_iters=20; o.seed=7;
    float *OV; size_t onv; int32_t *OF; size_t onf;
    int rc = CVT_remesh(a, V, nv, F, nf, 64, &o, NULL, &OV, &onv, &OF, &onf);
    int ok = (rc==0 && onv==64 && onf>0);
    /* index range check */
    for (size_t t=0; t<onf*3 && ok; t++) if (OF[t]<0 || (size_t)OF[t]>=onv) ok=0;
    fprintf(stderr, "cvt selftest: rc=%d onv=%zu onf=%zu -> %s\n",
            rc, onv, onf, ok?"OK":"FAIL");
    Arena_dispose(&a);
    return ok ? 0 : 1;
}

int main(int argc, char **argv) {
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc < 4) {
        fprintf(stderr, "usage: %s <in.obj> <out.obj> <target_verts> [--iters N] [--seed S]\n", argv[0]);
        fprintf(stderr, "       %s --selftest\n", argv[0]);
        return 2;
    }
    const char *in = argv[1], *out = argv[2];
    size_t target = (size_t)strtoull(argv[3], NULL, 10);
    CvtOpts o; CVT_default_opts(&o); o.verbose = 1;
    for (int i = 4; i < argc; i++) {
        if (!strcmp(argv[i], "--iters") && i+1 < argc) o.n_iters = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i+1 < argc) o.seed = (uint32_t)strtoul(argv[++i], NULL, 10);
    }

    Arena_T a = Arena_new();
    float *V; size_t nv; int32_t *F; size_t nf;
    if (ObjIO_read(a, in, &V, &nv, &F, &nf) != 0) {
        fprintf(stderr, "error: cannot read %s\n", in); return 1;
    }
    fprintf(stderr, "in:  %zu verts, %zu faces; target %zu sites, %d iters\n",
            nv, nf, target, o.n_iters);

    float *OV; size_t onv; int32_t *OF; size_t onf;
    int rc = CVT_remesh(a, V, nv, F, nf, target, &o, NULL, &OV, &onv, &OF, &onf);
    if (rc != 0) { fprintf(stderr, "error: CVT_remesh rc=%d\n", rc); return 1; }
    fprintf(stderr, "out: %zu verts, %zu faces\n", onv, onf);

    if (ObjIO_write(out, OV, onv, OF, onf) != 0) {
        fprintf(stderr, "error: cannot write %s\n", out); return 1;
    }
    fprintf(stderr, "wrote %s\n", out);
    Arena_dispose(&a);
    return 0;
}
