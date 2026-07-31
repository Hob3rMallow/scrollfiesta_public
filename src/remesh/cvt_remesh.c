#include "cvt_remesh.h"
#include "rvd_clip.h"
#include "../common/closest_tri.h"
#include "../common/kdtree.h"

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void CVT_default_opts(CvtOpts *o) {
    if (!o) return;
    o->n_iters        = 50;
    o->lambda_na      = 1.0f;
    o->lambda_cvt0    = 1.0f;
    o->cvt_decay      = 0.95f;
    o->lambda_cvt_min = 1e-3f;
    o->seed           = 0;
    o->verbose        = 0;
}

/* Solve a symmetric 3x3 system A x = b, A stored as {a00,a01,a02,a11,a12,a22}.
 * A is positive-definite here (lambda_CVT*m*I regularizes N), so the cofactor
 * inverse is stable; falls back to a diagonal solve if the determinant underflows. */
static void sym3_solve(const double A[6], const double b[3], double x[3]) {
    double a0=A[0],a1=A[1],a2=A[2],a3=A[3],a4=A[4],a5=A[5];
    double c00=a3*a5-a4*a4, c01=a2*a4-a1*a5, c02=a1*a4-a2*a3;
    double det=a0*c00+a1*c01+a2*c02;
    if (fabs(det) < 1e-20) {
        x[0]= a0>0?b[0]/a0:b[0]; x[1]= a3>0?b[1]/a3:b[1]; x[2]= a5>0?b[2]/a5:b[2];
        return;
    }
    double c11=a0*a5-a2*a2, c12=a1*a2-a0*a4, c22=a0*a3-a1*a1, inv=1.0/det;
    x[0]=inv*(c00*b[0]+c01*b[1]+c02*b[2]);
    x[1]=inv*(c01*b[0]+c11*b[1]+c12*b[2]);
    x[2]=inv*(c02*b[0]+c12*b[1]+c22*b[2]);
}

/* xorshift32 -> uniform [0,1) */
static uint32_t xs32(uint32_t *s) { uint32_t x=*s; x^=x<<13; x^=x>>17; x^=x<<5; *s=x; return x; }
static double   rnd01(uint32_t *s) { return (double)(xs32(s) >> 8) * (1.0/16777216.0); }

/* Closest point on the mesh to p (brute force over triangles — M1/M2). */
static void project_to_mesh(const double *V, const int32_t *F, size_t nf,
                            const double p[3], double out[3]) {
    double best = 1e300;
    out[0]=p[0]; out[1]=p[1]; out[2]=p[2];
    for (size_t f = 0; f < nf; f++) {
        const double *A=&V[(size_t)F[f*3+0]*3];
        const double *B=&V[(size_t)F[f*3+1]*3];
        const double *C=&V[(size_t)F[f*3+2]*3];
        double q[3], u, v;
        ClosestTri_point(p, A, B, C, q, &u, &v);
        double dz=p[0]-q[0], dy=p[1]-q[1], dx=p[2]-q[2];
        double d2 = dz*dz+dy*dy+dx*dx;
        if (d2 < best) { best=d2; out[0]=q[0]; out[1]=q[1]; out[2]=q[2]; }
    }
}

/* Closest point on the mesh via a triangle-centroid KD-tree. Any triangle that
 * could be the closest has its centroid within (nearest-centroid distance +
 * 2*max_tri_r) of p, so one ball query returns all candidates. Replaces the O(nf)
 * brute force -- the dominant cost on the dense pre-QEM meshes. */
static void project_to_mesh_kd(KDTree_T kd, const double *V, const int32_t *F,
                               double max_tri_r, const double p[3], double out[3],
                               int32_t *cand, size_t cand_cap) {
    float pf[3] = {(float)p[0],(float)p[1],(float)p[2]};
    float d0sq = 0.0f;
    size_t i0 = KDTree_nearest(kd, pf, &d0sq);
    double R = sqrt((double)d0sq) + 2.0*max_tri_r + 1e-6;
    size_t nc = KDTree_ball_query(kd, pf, (float)(R*R), cand, cand_cap);
    if (nc == 0) { cand[0] = (int32_t)i0; nc = 1; }
    double best = 1e300;
    out[0]=p[0]; out[1]=p[1]; out[2]=p[2];
    for (size_t k = 0; k < nc; k++) {
        size_t f = (size_t)cand[k];
        const double *A=&V[(size_t)F[f*3+0]*3], *B=&V[(size_t)F[f*3+1]*3], *C=&V[(size_t)F[f*3+2]*3];
        double q[3], u, v;
        ClosestTri_point(p, A, B, C, q, &u, &v);
        double dz=p[0]-q[0], dy=p[1]-q[1], dx=p[2]-q[2], d2=dz*dz+dy*dy+dx*dx;
        if (d2 < best) { best=d2; out[0]=q[0]; out[1]=q[1]; out[2]=q[2]; }
    }
}

static double closest_on_seg(const double p[3], const double a[3], const double b[3], double out[3]) {
    double ab[3]={b[0]-a[0],b[1]-a[1],b[2]-a[2]};
    double t = (p[0]-a[0])*ab[0]+(p[1]-a[1])*ab[1]+(p[2]-a[2])*ab[2];
    double d = ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2];
    t = d>0.0 ? t/d : 0.0; if (t<0) t=0; if (t>1) t=1;
    out[0]=a[0]+t*ab[0]; out[1]=a[1]+t*ab[1]; out[2]=a[2]+t*ab[2];
    double dz=p[0]-out[0], dy=p[1]-out[1], dx=p[2]-out[2];
    return dz*dz+dy*dy+dx*dx;
}

/* Closest point on the union of boundary segments (seg[nseg*6] = a.xyz,b.xyz). */
static void project_to_boundary(const double *seg, size_t nseg, const double p[3], double out[3]) {
    double best=1e300;
    out[0]=p[0]; out[1]=p[1]; out[2]=p[2];
    for (size_t s=0; s<nseg; s++) {
        double q[3];
        double d2 = closest_on_seg(p, &seg[s*6], &seg[s*6+3], q);
        if (d2<best) { best=d2; out[0]=q[0]; out[1]=q[1]; out[2]=q[2]; }
    }
}

typedef struct {
    int32_t *loopv;    /* boundary vertex indices, concatenated per loop */
    int32_t *loopoff;  /* [nloops+1] offsets into loopv */
    int      nloops;
    double  *seg;      /* [nseg*6] boundary segment endpoint positions */
    size_t   nseg;
} Boundary;

static int cmp_edge(const void *pa, const void *pb) {
    const int32_t *a=(const int32_t*)pa, *b=(const int32_t*)pb;
    if (a[0]!=b[0]) return a[0]<b[0]?-1:1;
    if (a[1]!=b[1]) return a[1]<b[1]?-1:1;
    return 0;
}

/* Extract boundary loops (edges incident to exactly one face) and chain them.
 * Assumes a manifold boundary (each boundary vertex has degree 2). */
static void extract_boundary(Arena_T a, const double *V, const int32_t *F,
                             size_t nf, size_t nv, Boundary *bd) {
    memset(bd, 0, sizeof(*bd));
    size_t ne = nf*3;
    int32_t *E = (int32_t*)ARENA_ALLOC(a, (size_t)(ne*2*sizeof(int32_t)));
    for (size_t f=0; f<nf; f++)
        for (int e=0; e<3; e++) {
            int32_t u=F[f*3+e], w=F[f*3+(e+1)%3];
            int32_t lo=u<w?u:w, hi=u<w?w:u;
            E[(f*3+e)*2]=lo; E[(f*3+e)*2+1]=hi;
        }
    qsort(E, ne, 2*sizeof(int32_t), cmp_edge);

    /* boundary edges = runs of length 1 */
    int32_t *nbr0=(int32_t*)ARENA_ALLOC(a,(size_t)(nv*sizeof(int32_t)));
    int32_t *nbr1=(int32_t*)ARENA_ALLOC(a,(size_t)(nv*sizeof(int32_t)));
    for (size_t i=0;i<nv;i++){ nbr0[i]=-1; nbr1[i]=-1; }
    size_t nbnd=0;
    /* first pass: count boundary edges */
    for (size_t i=0;i<ne;) {
        size_t j=i+1;
        while (j<ne && E[j*2]==E[i*2] && E[j*2+1]==E[i*2+1]) j++;
        if (j-i==1) nbnd++;
        i=j;
    }
    bd->nseg=nbnd;
    bd->seg=(double*)ARENA_ALLOC(a,(size_t)(nbnd*6*sizeof(double)));
    size_t si=0;
    for (size_t i=0;i<ne;) {
        size_t j=i+1;
        while (j<ne && E[j*2]==E[i*2] && E[j*2+1]==E[i*2+1]) j++;
        if (j-i==1) {
            int32_t u=E[i*2], w=E[i*2+1];
            const double *Pu=&V[(size_t)u*3], *Pw=&V[(size_t)w*3];
            bd->seg[si*6+0]=Pu[0]; bd->seg[si*6+1]=Pu[1]; bd->seg[si*6+2]=Pu[2];
            bd->seg[si*6+3]=Pw[0]; bd->seg[si*6+4]=Pw[1]; bd->seg[si*6+5]=Pw[2];
            si++;
            /* adjacency for chaining */
            if (nbr0[u]<0) nbr0[u]=w; else nbr1[u]=w;
            if (nbr0[w]<0) nbr0[w]=u; else nbr1[w]=u;
        }
        i=j;
    }

    /* chain into loops */
    uint8_t *vis=(uint8_t*)ARENA_CALLOC(a,(size_t)nv,1);
    bd->loopv=(int32_t*)ARENA_ALLOC(a,(size_t)((nbnd+1)*sizeof(int32_t)));
    bd->loopoff=(int32_t*)ARENA_ALLOC(a,(size_t)((nbnd+2)*sizeof(int32_t)));
    int nloops=0; int32_t w=0; bd->loopoff[0]=0;
    for (size_t s=0; s<nv; s++) {
        if (nbr0[s]<0 || vis[s]) continue;
        int32_t prev=-1, cur=(int32_t)s;
        while (cur>=0 && !vis[cur]) {
            vis[cur]=1; bd->loopv[w++]=cur;
            int32_t nx = (nbr0[cur]!=prev) ? nbr0[cur] : nbr1[cur];
            prev=cur; cur=nx;
        }
        bd->loopoff[++nloops]=w;
    }
    bd->nloops=nloops;
}

/* Place m points evenly by arc length along a loop of vertices; writes to out[m*3]. */
static void place_on_loop(const double *V, const int32_t *idx, int n, int m, double *out) {
    if (n < 2 || m < 1) return;
    double L=0.0;
    for (int i=0;i<n;i++) {
        const double *a=&V[(size_t)idx[i]*3], *b=&V[(size_t)idx[(i+1)%n]*3];
        double dz=b[0]-a[0],dy=b[1]-a[1],dx=b[2]-a[2];
        L += sqrt(dz*dz+dy*dy+dx*dx);
    }
    if (L<=0.0) { for (int k=0;k<m;k++){ const double*p=&V[(size_t)idx[0]*3]; out[k*3]=p[0];out[k*3+1]=p[1];out[k*3+2]=p[2]; } return; }
    double step=L/m, target=0.0, acc=0.0; int seg=0;
    const double *a=&V[(size_t)idx[0]*3], *b=&V[(size_t)idx[1%n]*3];
    double dz=b[0]-a[0],dy=b[1]-a[1],dx=b[2]-a[2];
    double seglen=sqrt(dz*dz+dy*dy+dx*dx);
    for (int k=0;k<m;k++) {
        while (acc+seglen < target && seg < n-1) {
            acc+=seglen; seg++;
            a=&V[(size_t)idx[seg]*3]; b=&V[(size_t)idx[(seg+1)%n]*3];
            dz=b[0]-a[0];dy=b[1]-a[1];dx=b[2]-a[2]; seglen=sqrt(dz*dz+dy*dy+dx*dx);
        }
        double t = seglen>0 ? (target-acc)/seglen : 0.0; if (t<0)t=0; if (t>1)t=1;
        out[k*3+0]=a[0]+t*dz; out[k*3+1]=a[1]+t*dy; out[k*3+2]=a[2]+t*dx;
        target+=step;
    }
}

/* No CVT-level Arena_save/restore (QEM model; avoids nested-mark aliasing with the
 * dual `tris` Rvd_accumulate returns). */
int CVT_remesh(Arena_T arena,
               const float *in_verts, size_t in_nv,
               const int32_t *in_faces, size_t in_nf,
               size_t target_nsites,
               const CvtOpts *opts,
               const CvtField *field,
               float **out_verts, size_t *out_nv,
               int32_t **out_faces, size_t *out_nf) {
    if (!arena || !in_verts || !in_faces || in_nv < 3 || in_nf < 1 || target_nsites < 3)
        return -1;
    if (out_verts) *out_verts = NULL;
    if (out_faces) *out_faces = NULL;
    if (out_nv) *out_nv = 0;
    if (out_nf) *out_nf = 0;

    CvtOpts o; if (opts) o = *opts; else CVT_default_opts(&o);
    uint32_t rng = o.seed ? o.seed : 12345u;
    size_t N = target_nsites;

    double *V = (double *)ARENA_ALLOC(arena, (size_t)(in_nv*3*sizeof(double)));
    for (size_t i = 0; i < in_nv*3; i++) V[i] = (double)in_verts[i];

    Rvd_T r = Rvd_new(arena, V, in_nv, in_faces, in_nf);
    if (!r) return -2;
    double surf_area = Rvd_mesh_area(r);
    if (surf_area <= 0.0) return -3;

    /* Projection acceleration: KD-tree over triangle centroids (+ the max
     * centroid-to-vertex radius, so the ball query is a provable superset). */
    float *tc = (float *)ARENA_ALLOC(arena, (size_t)(in_nf*3*sizeof(float)));
    double max_tri_r = 0.0;
    for (size_t f = 0; f < in_nf; f++) {
        const double *A=&V[(size_t)in_faces[f*3+0]*3];
        const double *B=&V[(size_t)in_faces[f*3+1]*3];
        const double *C=&V[(size_t)in_faces[f*3+2]*3];
        double c0=(A[0]+B[0]+C[0])/3.0, c1=(A[1]+B[1]+C[1])/3.0, c2=(A[2]+B[2]+C[2])/3.0;
        tc[f*3+0]=(float)c0; tc[f*3+1]=(float)c1; tc[f*3+2]=(float)c2;
        double ra=sqrt((A[0]-c0)*(A[0]-c0)+(A[1]-c1)*(A[1]-c1)+(A[2]-c2)*(A[2]-c2));
        double rb=sqrt((B[0]-c0)*(B[0]-c0)+(B[1]-c1)*(B[1]-c1)+(B[2]-c2)*(B[2]-c2));
        double rc2=sqrt((C[0]-c0)*(C[0]-c0)+(C[1]-c1)*(C[1]-c1)+(C[2]-c2)*(C[2]-c2));
        double rr=ra>rb?(ra>rc2?ra:rc2):(rb>rc2?rb:rc2);
        if (rr>max_tri_r) max_tri_r=rr;
    }
    KDTree_T proj_kd = KDTree_new(arena, tc, in_nf);
    int32_t *proj_cand = (int32_t *)ARENA_ALLOC(arena, (size_t)(1024*sizeof(int32_t)));

    /* Seeding CDF over triangles. Uniform (field==NULL): weight = area, so sampling
     * is area-weighted. Graded (field!=NULL): weight = area / h(centroid)^2, so seed
     * density ~ 1/h^2 (spacing ~ h); the weighted total then DERIVES the site count
     * N = round(2 * sum(area/h^2)) -- the generalization of the uniform N = 2A/h^2. */
    double *cum = (double *)ARENA_ALLOC(arena, (size_t)(in_nf*sizeof(double)));
    double tot = 0.0;
    for (size_t f = 0; f < in_nf; f++) {
        const double *A=&V[(size_t)in_faces[f*3+0]*3];
        const double *B=&V[(size_t)in_faces[f*3+1]*3];
        const double *C=&V[(size_t)in_faces[f*3+2]*3];
        double ab[3]={B[0]-A[0],B[1]-A[1],B[2]-A[2]};
        double ac[3]={C[0]-A[0],C[1]-A[1],C[2]-A[2]};
        double cr[3]={ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0]};
        double af = 0.5*sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
        if (field) {
            double c0=(A[0]+B[0]+C[0])/3.0, c1=(A[1]+B[1]+C[1])/3.0, c2=(A[2]+B[2]+C[2])/3.0;
            af *= sizing_seed_rho(field, c0, c1, c2);
        }
        tot += af;
        cum[f] = tot;
    }
    if (field) {
        long long nn = (long long)(2.0*tot + 0.5);   /* derived site count */
        if (nn < 4) nn = 4;
        if ((size_t)nn > in_nf) nn = (long long)in_nf;
        N = (size_t)nn;
    }

    /* boundary loops + (uniform-path) target spacing */
    Boundary bd; extract_boundary(arena, V, in_faces, in_nf, in_nv, &bd);
    double spacing = sqrt(2.0*surf_area/(double)N);   /* approx target edge length */

    double *sites = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    uint8_t *bflag = (uint8_t *)ARENA_CALLOC(arena, (size_t)N, 1);
    size_t ns = 0;

    /* boundary seeds first (constrained to their loop) */
    for (int k = 0; k < bd.nloops && ns < N; k++) {
        int lo = bd.loopoff[k], hi = bd.loopoff[k+1], n = hi - lo;
        if (n < 2) continue;
        double L=0.0, hc=0.0;
        for (int i=0;i<n;i++){
            const double *aa=&V[(size_t)bd.loopv[lo+i]*3], *bb=&V[(size_t)bd.loopv[lo+(i+1)%n]*3];
            double dz=bb[0]-aa[0],dy=bb[1]-aa[1],dx=bb[2]-aa[2];
            double sl=sqrt(dz*dz+dy*dy+dx*dx); L+=sl;
            if (field) {
                double mz=(aa[0]+bb[0])*0.5, my=(aa[1]+bb[1])*0.5, mx=(aa[2]+bb[2])*0.5;
                hc += sl / sizing_h(field, mz, my, mx);   /* local-h arc-length count */
            }
        }
        /* Count by LOCAL edge length: seam rims (h==h_seam) get dense verts the weld
         * bridge can span; interior-hole rims (h==h_interior) stay coarse. */
        int m = field ? (int)(hc + 0.5) : (int)(L/spacing + 0.5);
        if (m < 3) m = 3;
        if ((size_t)m > N - ns) m = (int)(N - ns);
        if (m < 1) break;
        place_on_loop(V, &bd.loopv[lo], n, m, &sites[ns*3]);
        for (int i=0;i<m;i++) bflag[ns+i]=1;
        ns += (size_t)m;
    }

    size_t n_bnd = ns;   /* boundary-constrained seeds precede the interior ones */

    /* interior seeds (area-weighted barycentric), kept clear of the boundary strip
     * so the boundary seeds own the rim cleanly. Without this an interior cell
     * wedges to the edge, breaks two boundary seeds' adjacency, and the dual
     * boundary saw-tooths inward. */
    double clear = 0.85 * spacing;
    for (; ns < N; ns++) {
        double pos[3] = {0,0,0};
        for (int tries = 0; ; tries++) {
            double rr = rnd01(&rng) * tot;
            size_t loo=0, hii=in_nf-1;
            while (loo < hii) { size_t mid=(loo+hii)/2; if (cum[mid] < rr) loo=mid+1; else hii=mid; }
            size_t f = loo;
            double u=rnd01(&rng), v=rnd01(&rng);
            if (u+v > 1.0) { u=1.0-u; v=1.0-v; }
            const double *A=&V[(size_t)in_faces[f*3+0]*3];
            const double *B=&V[(size_t)in_faces[f*3+1]*3];
            const double *C=&V[(size_t)in_faces[f*3+2]*3];
            pos[0]=A[0]+u*(B[0]-A[0])+v*(C[0]-A[0]);
            pos[1]=A[1]+u*(B[1]-A[1])+v*(C[1]-A[1]);
            pos[2]=A[2]+u*(B[2]-A[2])+v*(C[2]-A[2]);
            if (bd.nseg == 0 || tries >= 8) break;
            double q[3]; project_to_boundary(bd.seg, bd.nseg, pos, q);
            double dz=pos[0]-q[0],dy=pos[1]-q[1],dx=pos[2]-q[2];
            double cl = field ? 0.85*sizing_h(field, pos[0],pos[1],pos[2]) : clear;
            if (dz*dz+dy*dy+dx*dx >= cl*cl) break;   /* far enough from the rim */
        }
        sites[ns*3+0]=pos[0]; sites[ns*3+1]=pos[1]; sites[ns*3+2]=pos[2];
    }

    double *area = (double *)ARENA_ALLOC(arena, (size_t)(N*sizeof(double)));
    double *cent = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    double *NA   = (double *)ARENA_ALLOC(arena, (size_t)(N*6*sizeof(double)));
    double *bNA  = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    double lna = (double)o.lambda_na;

    for (int it = 0; it < o.n_iters; it++) {
        double lcvt = (double)o.lambda_cvt0 * pow((double)o.cvt_decay, (double)it);
        if (lcvt < (double)o.lambda_cvt_min) lcvt = (double)o.lambda_cvt_min;
        if (Rvd_accumulate(r, sites, N, area, cent, NA, bNA, NULL, NULL, field) != 0) return -4;
        double moved = 0.0;
        for (size_t i = 0; i < N; i++) {
            double np[3];
            if (bflag[i] && bd.nseg > 0) {
                /* boundary site: CVT-relax along its loop */
                project_to_boundary(bd.seg, bd.nseg, &cent[i*3], np);
            } else if (area[i] <= 0.0) {
                np[0]=sites[i*3+0]; np[1]=sites[i*3+1]; np[2]=sites[i*3+2];  /* empty cell */
            } else {
                /* interior CWF solve: (lna*N + lcvt*m*I) x = lna*b + lcvt*m*c */
                double m = area[i];
                const double *Ni=&NA[i*6], *bi=&bNA[i*3], *ci=&cent[i*3];
                double A6[6] = { lna*Ni[0]+lcvt*m, lna*Ni[1],        lna*Ni[2],
                                 lna*Ni[3]+lcvt*m, lna*Ni[4],        lna*Ni[5]+lcvt*m };
                double rhs[3] = { lna*bi[0]+lcvt*m*ci[0],
                                  lna*bi[1]+lcvt*m*ci[1],
                                  lna*bi[2]+lcvt*m*ci[2] };
                double x[3]; sym3_solve(A6, rhs, x);
                project_to_mesh_kd(proj_kd, V, in_faces, max_tri_r, x, np, proj_cand, 1024);
                /* keep interior sites off the boundary strip every iteration (not
                 * just at seeding), so boundary seeds keep sole ownership of the
                 * rim and the dual boundary stays the clean boundary-seed loop. */
                if (bd.nseg > 0) {
                    double q[3]; project_to_boundary(bd.seg, bd.nseg, np, q);
                    double dz=np[0]-q[0], dy=np[1]-q[1], dx=np[2]-q[2];
                    double dd=sqrt(dz*dz+dy*dy+dx*dx);
                    double cl = field ? 0.85*sizing_h(field, np[0],np[1],np[2]) : clear;
                    if (dd < cl && dd > 1e-9) {
                        double sc=cl/dd;
                        double push[3]={ q[0]+dz*sc, q[1]+dy*sc, q[2]+dx*sc };
                        project_to_mesh_kd(proj_kd, V, in_faces, max_tri_r, push, np, proj_cand, 1024);
                    }
                }
            }
            double dz=np[0]-sites[i*3+0], dy=np[1]-sites[i*3+1], dx=np[2]-sites[i*3+2];
            moved += sqrt(dz*dz+dy*dy+dx*dx);
            sites[i*3+0]=np[0]; sites[i*3+1]=np[1]; sites[i*3+2]=np[2];
        }
        if (o.verbose) fprintf(stderr, "  cvt iter %2d: move %.5g  lcvt %.4f  (%zu bnd / %zu)\n",
                               it, moved, lcvt, n_bnd, N);
    }

    int32_t *tris = NULL; size_t nt = 0;
    /* Dual is the geometric Restricted-Delaunay of the converged sites -- density
     * doesn't affect it, so pass NULL (the area/centroid outputs here are scratch). */
    if (Rvd_accumulate(r, sites, N, area, cent, NULL, NULL, &tris, &nt, NULL) != 0) return -5;

    float *ov = (float *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(float)));
    for (size_t i = 0; i < N*3; i++) ov[i] = (float)sites[i];

    *out_verts = ov; *out_nv = N;
    *out_faces = tris; *out_nf = nt;
    return 0;
}

/* project_to_boundary + the index of the nearest segment (for LOCAL spacing). */
static void project_to_boundary_seg(const double *seg, size_t nseg,
                                    const double p[3], double out[3],
                                    size_t *out_seg) {
    double best = 1e300;
    out[0]=p[0]; out[1]=p[1]; out[2]=p[2];
    *out_seg = 0;
    for (size_t s=0; s<nseg; s++) {
        double q[3];
        double d2 = closest_on_seg(p, &seg[s*6], &seg[s*6+3], q);
        if (d2<best) { best=d2; out[0]=q[0]; out[1]=q[1]; out[2]=q[2]; *out_seg=s; }
    }
}

/* Pinned-boundary variant: see cvt_remesh.h. Mirrors CVT_remesh with three
 * deltas -- boundary sites are the EXACT input boundary verts (no arc-length
 * resampling, no along-loop relax, bit-exact positions), the interior count is
 * derived from target_h, and the interior clearance uses the LOCAL boundary
 * segment length so coarse junction rings and fine hole rims each get a
 * proportionate berth. */
int CVT_remesh_pinned(Arena_T arena,
                      const float *in_verts, size_t in_nv,
                      const int32_t *in_faces, size_t in_nf,
                      double target_h,
                      const CvtOpts *opts,
                      float **out_verts, size_t *out_nv,
                      int32_t **out_faces, size_t *out_nf,
                      int32_t **out_site_src, size_t *out_n_pin) {
    if (!arena || !in_verts || !in_faces || in_nv < 3 || in_nf < 1 || target_h <= 0.0)
        return -1;
    if (out_verts) *out_verts = NULL;
    if (out_faces) *out_faces = NULL;
    if (out_nv) *out_nv = 0;
    if (out_nf) *out_nf = 0;
    if (out_site_src) *out_site_src = NULL;
    if (out_n_pin) *out_n_pin = 0;

    CvtOpts o; if (opts) o = *opts; else CVT_default_opts(&o);
    uint32_t rng = o.seed ? o.seed : 12345u;

    double *V = (double *)ARENA_ALLOC(arena, (size_t)(in_nv*3*sizeof(double)));
    for (size_t i = 0; i < in_nv*3; i++) V[i] = (double)in_verts[i];

    Rvd_T r = Rvd_new(arena, V, in_nv, in_faces, in_nf);
    if (!r) return -2;
    double surf_area = Rvd_mesh_area(r);
    if (surf_area <= 0.0) return -3;

    /* Boundary loops: every loop vert is a pinned site, in loop order. */
    Boundary bd; extract_boundary(arena, V, in_faces, in_nf, in_nv, &bd);
    if (bd.nloops == 0 || bd.nseg == 0) return -6;   /* closed patch: not our case */
    size_t n_pin = (size_t)bd.loopoff[bd.nloops];
    if (n_pin < 3) return -6;

    /* Per-segment length, for the scale-aware clearance. */
    double *seglen = (double *)ARENA_ALLOC(arena, (size_t)(bd.nseg*sizeof(double)));
    for (size_t s = 0; s < bd.nseg; s++) {
        double dz=bd.seg[s*6+3]-bd.seg[s*6+0];
        double dy=bd.seg[s*6+4]-bd.seg[s*6+1];
        double dx=bd.seg[s*6+5]-bd.seg[s*6+2];
        seglen[s] = sqrt(dz*dz+dy*dy+dx*dx);
    }

    /* Interior count from the target spacing (uniform N = 2A/h^2). */
    size_t n_int = (size_t)(2.0*surf_area/(target_h*target_h) + 0.5);
    size_t N = n_pin + n_int;

    /* Projection acceleration (as CVT_remesh). */
    float *tc = (float *)ARENA_ALLOC(arena, (size_t)(in_nf*3*sizeof(float)));
    double max_tri_r = 0.0;
    for (size_t f = 0; f < in_nf; f++) {
        const double *A=&V[(size_t)in_faces[f*3+0]*3];
        const double *B=&V[(size_t)in_faces[f*3+1]*3];
        const double *C=&V[(size_t)in_faces[f*3+2]*3];
        double c0=(A[0]+B[0]+C[0])/3.0, c1=(A[1]+B[1]+C[1])/3.0, c2=(A[2]+B[2]+C[2])/3.0;
        tc[f*3+0]=(float)c0; tc[f*3+1]=(float)c1; tc[f*3+2]=(float)c2;
        double ra=sqrt((A[0]-c0)*(A[0]-c0)+(A[1]-c1)*(A[1]-c1)+(A[2]-c2)*(A[2]-c2));
        double rb=sqrt((B[0]-c0)*(B[0]-c0)+(B[1]-c1)*(B[1]-c1)+(B[2]-c2)*(B[2]-c2));
        double rc2=sqrt((C[0]-c0)*(C[0]-c0)+(C[1]-c1)*(C[1]-c1)+(C[2]-c2)*(C[2]-c2));
        double rr=ra>rb?(ra>rc2?ra:rc2):(rb>rc2?rb:rc2);
        if (rr>max_tri_r) max_tri_r=rr;
    }
    KDTree_T proj_kd = KDTree_new(arena, tc, in_nf);
    int32_t *proj_cand = (int32_t *)ARENA_ALLOC(arena, (size_t)(1024*sizeof(int32_t)));

    /* Area CDF for interior seeding (uniform). */
    double *cum = (double *)ARENA_ALLOC(arena, (size_t)(in_nf*sizeof(double)));
    double tot = 0.0;
    for (size_t f = 0; f < in_nf; f++) {
        const double *A=&V[(size_t)in_faces[f*3+0]*3];
        const double *B=&V[(size_t)in_faces[f*3+1]*3];
        const double *C=&V[(size_t)in_faces[f*3+2]*3];
        double ab[3]={B[0]-A[0],B[1]-A[1],B[2]-A[2]};
        double ac[3]={C[0]-A[0],C[1]-A[1],C[2]-A[2]};
        double cr[3]={ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0]};
        tot += 0.5*sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
        cum[f] = tot;
    }

    double *sites = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    int32_t *src  = (int32_t *)ARENA_ALLOC(arena, (size_t)(N*sizeof(int32_t)));

    /* Pinned sites: exact boundary verts, loop order, bit-exact positions. */
    for (size_t i = 0; i < n_pin; i++) {
        int32_t vid = bd.loopv[i];
        sites[i*3+0]=V[(size_t)vid*3+0];
        sites[i*3+1]=V[(size_t)vid*3+1];
        sites[i*3+2]=V[(size_t)vid*3+2];
        src[i] = vid;
    }

    /* Interior seeds, scale-aware clearance from the boundary polyline. */
    for (size_t ns = n_pin; ns < N; ns++) {
        double pos[3] = {0,0,0};
        for (int tries = 0; ; tries++) {
            double rr = rnd01(&rng) * tot;
            size_t loo=0, hii=in_nf-1;
            while (loo < hii) { size_t mid=(loo+hii)/2; if (cum[mid] < rr) loo=mid+1; else hii=mid; }
            size_t f = loo;
            double u=rnd01(&rng), v=rnd01(&rng);
            if (u+v > 1.0) { u=1.0-u; v=1.0-v; }
            const double *A=&V[(size_t)in_faces[f*3+0]*3];
            const double *B=&V[(size_t)in_faces[f*3+1]*3];
            const double *C=&V[(size_t)in_faces[f*3+2]*3];
            pos[0]=A[0]+u*(B[0]-A[0])+v*(C[0]-A[0]);
            pos[1]=A[1]+u*(B[1]-A[1])+v*(C[1]-A[1]);
            pos[2]=A[2]+u*(B[2]-A[2])+v*(C[2]-A[2]);
            if (tries >= 8) break;
            double q[3]; size_t sidx;
            project_to_boundary_seg(bd.seg, bd.nseg, pos, q, &sidx);
            double dz=pos[0]-q[0],dy=pos[1]-q[1],dx=pos[2]-q[2];
            /* Wedge rule: an interior site closer to a boundary segment than
             * HALF its length can capture the segment's midpoint and sever
             * the two pinned endpoints' dual adjacency (= a dropped boundary
             * edge = a crack). 0.85 * seglen clears the half-length threshold
             * for ANY segment; deliberately NOT capped by target_h -- long
             * segments (despur stair-step diagonals) need proportionally wide
             * berths, which is geometrically required, not a density knob. */
            double cl = 0.85 * seglen[sidx];
            if (dz*dz+dy*dy+dx*dx >= cl*cl) break;
        }
        sites[ns*3+0]=pos[0]; sites[ns*3+1]=pos[1]; sites[ns*3+2]=pos[2];
        src[ns] = -1;
    }

    double *area = (double *)ARENA_ALLOC(arena, (size_t)(N*sizeof(double)));
    double *cent = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    double *NA   = (double *)ARENA_ALLOC(arena, (size_t)(N*6*sizeof(double)));
    double *bNA  = (double *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(double)));
    double lna = (double)o.lambda_na;

    for (int it = 0; it < o.n_iters; it++) {
        double lcvt = (double)o.lambda_cvt0 * pow((double)o.cvt_decay, (double)it);
        if (lcvt < (double)o.lambda_cvt_min) lcvt = (double)o.lambda_cvt_min;
        if (Rvd_accumulate(r, sites, N, area, cent, NA, bNA, NULL, NULL, NULL) != 0) return -4;
        for (size_t i = n_pin; i < N; i++) {   /* pinned sites NEVER move */
            double np[3];
            if (area[i] <= 0.0) {
                np[0]=sites[i*3+0]; np[1]=sites[i*3+1]; np[2]=sites[i*3+2];
            } else {
                double m = area[i];
                const double *Ni=&NA[i*6], *bi=&bNA[i*3], *ci=&cent[i*3];
                double A6[6] = { lna*Ni[0]+lcvt*m, lna*Ni[1],        lna*Ni[2],
                                 lna*Ni[3]+lcvt*m, lna*Ni[4],        lna*Ni[5]+lcvt*m };
                double rhs[3] = { lna*bi[0]+lcvt*m*ci[0],
                                  lna*bi[1]+lcvt*m*ci[1],
                                  lna*bi[2]+lcvt*m*ci[2] };
                double x[3]; sym3_solve(A6, rhs, x);
                project_to_mesh_kd(proj_kd, V, in_faces, max_tri_r, x, np, proj_cand, 1024);
                /* scale-aware push off the boundary strip, every iteration
                 * (same uncapped wedge rule as at seeding) */
                {
                    double q[3]; size_t sidx;
                    project_to_boundary_seg(bd.seg, bd.nseg, np, q, &sidx);
                    double dz=np[0]-q[0], dy=np[1]-q[1], dx=np[2]-q[2];
                    double dd=sqrt(dz*dz+dy*dy+dx*dx);
                    double cl = 0.85 * seglen[sidx];
                    if (dd < cl && dd > 1e-9) {
                        double sc=cl/dd;
                        double push[3]={ q[0]+dz*sc, q[1]+dy*sc, q[2]+dx*sc };
                        project_to_mesh_kd(proj_kd, V, in_faces, max_tri_r, push, np, proj_cand, 1024);
                    }
                }
            }
            sites[i*3+0]=np[0]; sites[i*3+1]=np[1]; sites[i*3+2]=np[2];
        }
    }

    int32_t *tris = NULL; size_t nt = 0;
    if (Rvd_accumulate(r, sites, N, area, cent, NULL, NULL, &tris, &nt, NULL) != 0) return -5;
    if (nt == 0) return -7;

    float *ov = (float *)ARENA_ALLOC(arena, (size_t)(N*3*sizeof(float)));
    for (size_t i = 0; i < N*3; i++) ov[i] = (float)sites[i];

    *out_verts = ov; *out_nv = N;
    *out_faces = tris; *out_nf = nt;
    if (out_site_src) *out_site_src = src;
    if (out_n_pin) *out_n_pin = n_pin;
    return 0;
}

/* ---- projection self-test: project_to_mesh_kd MUST equal the O(nf) brute force on
 * every query, including on meshes with wildly-varying triangle sizes (a giant
 * triangle whose closest point is far from any centroid) -- that is the case the
 * "nearest-centroid + 2*max_tri_r" ball radius must still capture. ---- */
static uint32_t g_ps = 0xB0Bu;
static double psr(double lo, double hi) { uint32_t x=g_ps; x^=x<<13; x^=x>>17; x^=x<<5; g_ps=x;
                                          return lo + (hi-lo)*((double)(x>>8)*(1.0/16777216.0)); }

int CVT_projection_selftest(void) {
    Arena_T a = Arena_new();
    int ok = 1;
    g_ps = 0xB0Bu;

    for (int kind = 0; kind < 3; kind++) {
        int M = 16, side = M+1;
        size_t base_nv = (size_t)side*side;
        size_t base_nf = (size_t)M*M*2;
        int giant = (kind == 2);
        size_t nv = base_nv + (giant ? 3 : 0);
        size_t nf = base_nf + (giant ? 1 : 0);
        double  *V = (double *)ARENA_ALLOC(a, (size_t)(nv*3*sizeof(double)));
        int32_t *F = (int32_t *)ARENA_ALLOC(a, (size_t)(nf*3*sizeof(int32_t)));

        for (int i=0;i<side;i++) for (int j=0;j<side;j++) {
            size_t idx=(size_t)i*side+j;
            double u=(double)i/M, v=(double)j/M;
            if (kind==1) { u=pow(u,1.7); v=pow(v,1.7); }          /* graded: mixed triangle sizes */
            V[idx*3+0]= (kind==0||kind==2)?0.0 : 0.4*sin(6.2831*u)*sin(6.2831*v); /* wavy for kind 1 */
            V[idx*3+1]=10.0*u; V[idx*3+2]=10.0*v;
        }
        size_t fi=0;
        for (int i=0;i<M;i++) for (int j=0;j<M;j++) {
            int32_t v00=(int32_t)((size_t)i*side+j), v01=(int32_t)((size_t)i*side+j+1);
            int32_t v10=(int32_t)((size_t)(i+1)*side+j), v11=(int32_t)((size_t)(i+1)*side+j+1);
            F[fi*3+0]=v00;F[fi*3+1]=v10;F[fi*3+2]=v11;fi++;
            F[fi*3+0]=v00;F[fi*3+1]=v11;F[fi*3+2]=v01;fi++;
        }
        if (giant) { /* one huge triangle far to the side */
            int32_t g0=(int32_t)base_nv, g1=g0+1, g2=g0+2;
            V[(size_t)g0*3+0]=0; V[(size_t)g0*3+1]=-40; V[(size_t)g0*3+2]=-40;
            V[(size_t)g1*3+0]=0; V[(size_t)g1*3+1]= 60; V[(size_t)g1*3+2]=-40;
            V[(size_t)g2*3+0]=0; V[(size_t)g2*3+1]=-40; V[(size_t)g2*3+2]= 60;
            F[base_nf*3+0]=g0; F[base_nf*3+1]=g1; F[base_nf*3+2]=g2;
        }

        /* build the projection KD-tree exactly as CVT_remesh does */
        float *tc = (float *)ARENA_ALLOC(a, (size_t)(nf*3*sizeof(float)));
        double max_tri_r = 0.0;
        for (size_t f=0; f<nf; f++) {
            const double *A=&V[(size_t)F[f*3+0]*3], *B=&V[(size_t)F[f*3+1]*3], *C=&V[(size_t)F[f*3+2]*3];
            double c0=(A[0]+B[0]+C[0])/3.0, c1=(A[1]+B[1]+C[1])/3.0, c2=(A[2]+B[2]+C[2])/3.0;
            tc[f*3+0]=(float)c0; tc[f*3+1]=(float)c1; tc[f*3+2]=(float)c2;
            double ra=sqrt((A[0]-c0)*(A[0]-c0)+(A[1]-c1)*(A[1]-c1)+(A[2]-c2)*(A[2]-c2));
            double rb=sqrt((B[0]-c0)*(B[0]-c0)+(B[1]-c1)*(B[1]-c1)+(B[2]-c2)*(B[2]-c2));
            double rc=sqrt((C[0]-c0)*(C[0]-c0)+(C[1]-c1)*(C[1]-c1)+(C[2]-c2)*(C[2]-c2));
            double rr=ra>rb?(ra>rc?ra:rc):(rb>rc?rb:rc);
            if (rr>max_tri_r) max_tri_r=rr;
        }
        KDTree_T kd = KDTree_new(a, tc, nf);
        int32_t *cand = (int32_t *)ARENA_ALLOC(a, (size_t)((nf+1)*sizeof(int32_t)));

        for (int q=0; q<3000; q++) {
            double p[3] = { psr(-6.0,6.0), psr(-45.0,65.0), psr(-45.0,65.0) };
            double okd[3], obf[3];
            project_to_mesh_kd(kd, V, F, max_tri_r, p, okd, cand, nf+1);
            project_to_mesh(V, F, nf, p, obf);
            double dk=(p[0]-okd[0])*(p[0]-okd[0])+(p[1]-okd[1])*(p[1]-okd[1])+(p[2]-okd[2])*(p[2]-okd[2]);
            double db=(p[0]-obf[0])*(p[0]-obf[0])+(p[1]-obf[1])*(p[1]-obf[1])+(p[2]-obf[2])*(p[2]-obf[2]);
            if (fabs(sqrt(dk)-sqrt(db)) > 1e-3*(sqrt(db)+1.0)) {
                fprintf(stderr, "CVT_projection_selftest FAIL: kind=%d q=(%.2f,%.2f,%.2f) kd_d=%.5f brute_d=%.5f\n",
                        kind, p[0],p[1],p[2], sqrt(dk), sqrt(db));
                ok = 0; break;
            }
        }
        if (!ok) break;
    }

    Arena_dispose(&a);
    fprintf(stderr, "CVT_projection_selftest: %s\n", ok ? "OK" : "FAIL");
    return ok ? 0 : -1;
}
