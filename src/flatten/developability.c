#define _USE_MATH_DEFINES
#include "developability.h"
#include "../common/eig3.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- vec3 (double) -------------------------------------------------------- */
static void  sub3(const double a[3], const double b[3], double o[3])
{ o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2]; }
static double dot3(const double a[3], const double b[3])
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static void  cross3(const double a[3], const double b[3], double o[3])
{ o[0]=a[1]*b[2]-a[2]*b[1]; o[1]=a[2]*b[0]-a[0]*b[2]; o[2]=a[0]*b[1]-a[1]*b[0]; }

static double corner_angle3(const double eA[3], const double eB[3])
{
    double la = sqrt(dot3(eA,eA)), lb = sqrt(dot3(eB,eB));
    if (la < 1e-20 || lb < 1e-20) return 0.0;
    double q = dot3(eA,eB)/(la*lb);
    if (q < -1.0) q = -1.0; if (q > 1.0) q = 1.0;
    return acos(q);
}

/* ---- boundary detection (undirected edge multiplicity) -------------------- */
static int cmp_u64d(const void *a, const void *b)
{ uint64_t x=*(const uint64_t*)a, y=*(const uint64_t*)b; return (x<y)?-1:(x>y)?1:0; }

static void mark_boundary(Arena_T arena, const int32_t *faces, size_t nf, size_t nv,
                          char *is_boundary)
{
    Arena_Mark mark = Arena_save(arena);
    size_t ne = nf*3;
    uint64_t *keys = (uint64_t *)ARENA_ALLOC(arena, (long)(ne*sizeof(uint64_t)));
    size_t m = 0;
    for (size_t t = 0; t < nf; t++) {
        int32_t v[3] = { faces[t*3+0], faces[t*3+1], faces[t*3+2] };
        for (int k = 0; k < 3; k++) {
            int32_t a = v[k], b = v[(k+1)%3];
            int32_t lo = a<b?a:b, hi = a<b?b:a;
            keys[m++] = (uint64_t)lo * (uint64_t)nv + (uint64_t)hi;
        }
    }
    qsort(keys, m, sizeof(uint64_t), cmp_u64d);
    for (size_t i = 0; i < nv; i++) is_boundary[i] = 0;
    size_t i = 0;
    while (i < m) {
        size_t j = i+1;
        while (j < m && keys[j] == keys[i]) j++;
        if (j - i == 1) {                          /* boundary edge */
            uint64_t key = keys[i];
            int32_t lo = (int32_t)(key / nv), hi = (int32_t)(key % nv);
            is_boundary[lo] = 1; is_boundary[hi] = 1;
        }
        i = j;
    }
    Arena_restore(arena, mark);
}

/* ---- energy + gradient ---------------------------------------------------- */

/* Assemble A_v = sum_corner theta * N N^T per interior vertex, return
 * E = sum lambda_min(A_v). Optionally store per-vertex eigenpair. */
static double assemble_energy(const double *X, size_t nv, const int32_t *faces,
                              size_t nf, const char *interior,
                              double *Asym /*[nv*6]*/, double *lam, double *uvec)
{
    for (size_t v = 0; v < nv*6; v++) Asym[v] = 0.0;
    for (size_t t = 0; t < nf; t++) {
        int32_t f[3] = { faces[t*3+0], faces[t*3+1], faces[t*3+2] };
        if (f[0]==f[1] || f[1]==f[2] || f[2]==f[0]) continue;
        const double *Xa=&X[f[0]*3], *Xb=&X[f[1]*3], *Xc=&X[f[2]*3];
        double e1[3], e2[3], n[3];
        sub3(Xb,Xa,e1); sub3(Xc,Xa,e2); cross3(e1,e2,n);
        double area2 = sqrt(dot3(n,n));
        if (area2 < 1e-20) continue;
        double N[3] = { n[0]/area2, n[1]/area2, n[2]/area2 };
        double NN[6] = { N[0]*N[0], N[1]*N[1], N[2]*N[2],
                         N[0]*N[1], N[0]*N[2], N[1]*N[2] };
        /* corner angles */
        double eab[3], eac[3], ebc[3];
        sub3(Xb,Xa,eab); sub3(Xc,Xa,eac); sub3(Xc,Xb,ebc);
        double th[3];
        { double t0[3],t1[3]; sub3(Xb,Xa,t0); sub3(Xc,Xa,t1); th[0]=corner_angle3(t0,t1); }
        { double t0[3],t1[3]; sub3(Xc,Xb,t0); sub3(Xa,Xb,t1); th[1]=corner_angle3(t0,t1); }
        { double t0[3],t1[3]; sub3(Xa,Xc,t0); sub3(Xb,Xc,t1); th[2]=corner_angle3(t0,t1); }
        for (int s = 0; s < 3; s++) {
            if (!interior[f[s]]) continue;
            double *A = &Asym[(size_t)f[s]*6];
            for (int k = 0; k < 6; k++) A[k] += th[s]*NN[k];
        }
        (void)eab; (void)eac; (void)ebc;
    }
    double E = 0.0;
    for (size_t v = 0; v < nv; v++) {
        if (!interior[v]) continue;
        const double *A = &Asym[v*6];
        double M[3][3] = { { A[0], A[3], A[4] },
                           { A[3], A[1], A[5] },
                           { A[4], A[5], A[2] } };
        double l, u[3];
        Eig3_smallest(M, &l, u);
        E += l;
        if (lam) lam[v] = l;
        if (uvec) { uvec[v*3+0]=u[0]; uvec[v*3+1]=u[1]; uvec[v*3+2]=u[2]; }
    }
    return E;
}

/* grad[nv*3] += d E / d x  (caller zeroes grad first). uvec holds the
 * lambda_min eigenvector per interior vertex (from assemble_energy). */
static void accumulate_grad(const double *X, size_t nv, const int32_t *faces,
                            size_t nf, const char *interior, const double *uvec,
                            double *grad)
{
    (void)nv;
    for (size_t t = 0; t < nf; t++) {
        int32_t f[3] = { faces[t*3+0], faces[t*3+1], faces[t*3+2] };
        if (f[0]==f[1] || f[1]==f[2] || f[2]==f[0]) continue;
        const double *X0=&X[f[0]*3], *X1=&X[f[1]*3], *X2=&X[f[2]*3];
        double e1[3], e2[3], n[3];
        sub3(X1,X0,e1); sub3(X2,X0,e2); cross3(e1,e2,n);
        double area2 = sqrt(dot3(n,n));
        if (area2 < 1e-20) continue;
        double N[3] = { n[0]/area2, n[1]/area2, n[2]/area2 };
        const double *Xp[3] = { X0, X1, X2 };

        for (int pa = 0; pa < 3; pa++) {
            if (!interior[f[pa]]) continue;
            const double *u = &uvec[(size_t)f[pa]*3];
            double s = dot3(u, N);
            int q1 = (pa+1)%3, q2 = (pa+2)%3;
            double eA[3], eB[3];
            sub3(Xp[q1], Xp[pa], eA);     /* apex -> q1 */
            sub3(Xp[q2], Xp[pa], eB);     /* apex -> q2 */
            double th = corner_angle3(eA, eB);

            /* corner-angle gradient per slot */
            double cAN[3], cBN[3];
            cross3(eA, N, cAN); cross3(eB, N, cBN);
            double la2 = dot3(eA,eA), lb2 = dot3(eB,eB);
            if (la2 < 1e-20 || lb2 < 1e-20) continue;
            double dth[3][3];
            for (int k = 0; k < 3; k++) {
                dth[q1][k] =  cAN[k]/la2;     /* d theta / d x_q1 = (eA x N)/|eA|^2 */
                dth[q2][k] = -cBN[k]/lb2;     /* d theta / d x_q2 = -(eB x N)/|eB|^2 */
                dth[pa][k] = -(dth[q1][k] + dth[q2][k]);
            }

            /* (u . N) gradient per FACE slot (fixed face e1,e2) */
            double w[3] = { (u[0]-s*N[0])/area2, (u[1]-s*N[1])/area2, (u[2]-s*N[2])/area2 };
            double e1xw[3], e2xw[3];
            cross3(e1, w, e1xw); cross3(e2, w, e2xw);
            double gUN[3][3];
            for (int k = 0; k < 3; k++) {
                gUN[0][k] = e1xw[k] - e2xw[k];
                gUN[1][k] = e2xw[k];
                gUN[2][k] = -e1xw[k];
            }

            double s2 = s*s, c2 = 2.0*th*s;
            for (int sl = 0; sl < 3; sl++) {
                double *g = &grad[(size_t)f[sl]*3];
                for (int k = 0; k < 3; k++) {
                    g[k] += s2*dth[sl][k] + c2*gUN[sl][k];
                }
            }
        }
    }
}

int Develop_optimize(Arena_T arena, float *verts, size_t nv,
                     const int32_t *faces, size_t nf,
                     const uint8_t *pin_mask, float *vert_normals,
                     const DevelopOpts *opts, DevelopStats *out)
{
    assert(arena);
    if (out) { out->e_initial=0; out->e_final=0; out->iters=0; out->max_disp=0; out->converged=1; }
    if (!verts || !faces) return -1;
    if (nv < 3 || nf < 1) return 0;                 /* trivial: no-op */

    int    max_iters = (opts && opts->max_iters>0)    ? opts->max_iters    : 100;
    double tol_grad  = (opts && opts->tol_grad>0)     ? opts->tol_grad     : 1e-6;
    double tol_e     = (opts && opts->tol_e>0)        ? opts->tol_e        : 1e-6;
    double max_disp  = (opts && opts->max_disp_vox>0) ? opts->max_disp_vox : 1.0;
    int    pin_bd    = (opts) ? opts->pin_boundary : 1;
    double c1        = (opts && opts->c1>0)           ? opts->c1           : 1e-4;
    double t0        = (opts && opts->t0>0)           ? opts->t0           : 1.0;
    int    max_ls    = (opts && opts->max_ls>0)       ? opts->max_ls       : 30;

    Arena_Mark mark = Arena_save(arena);
    double *X     = (double *)ARENA_ALLOC(arena, (long)(nv*3*sizeof(double)));
    double *X0    = (double *)ARENA_ALLOC(arena, (long)(nv*3*sizeof(double)));
    double *Xtry  = (double *)ARENA_ALLOC(arena, (long)(nv*3*sizeof(double)));
    for (size_t i = 0; i < nv*3; i++) { X[i] = (double)verts[i]; X0[i] = X[i]; }

    char *is_b = (char *)ARENA_ALLOC(arena, (long)nv);
    mark_boundary(arena, faces, nf, nv, is_b);
    char *interior = (char *)ARENA_ALLOC(arena, (long)nv);
    char *movable  = (char *)ARENA_ALLOC(arena, (long)nv);
    for (size_t v = 0; v < nv; v++) {
        interior[v] = (char)!is_b[v];
        int mv = interior[v] || !pin_bd;            /* boundary movable only if not pinned */
        if (pin_mask && pin_mask[v]) mv = 0;
        if (is_b[v] && pin_bd) mv = 0;
        movable[v] = (char)mv;
    }

    double *Asym = (double *)ARENA_ALLOC(arena, (long)(nv*6*sizeof(double)));
    double *uvec = (double *)ARENA_ALLOC(arena, (long)(nv*3*sizeof(double)));
    double *grad = (double *)ARENA_ALLOC(arena, (long)(nv*3*sizeof(double)));

    double E = assemble_energy(X, nv, faces, nf, interior, Asym, NULL, uvec);
    double e_init = E;
    int iter = 0, converged = 0;
    for (iter = 0; iter < max_iters; iter++) {
        for (size_t i = 0; i < nv*3; i++) grad[i] = 0.0;
        accumulate_grad(X, nv, faces, nf, interior, uvec, grad);
        /* zero non-movable, measure */
        double gmax = 0.0, gd = 0.0;
        for (size_t v = 0; v < nv; v++) {
            if (!movable[v]) { grad[v*3]=grad[v*3+1]=grad[v*3+2]=0.0; continue; }
            for (int k = 0; k < 3; k++) {
                double g = grad[v*3+k];
                if (fabs(g) > gmax) gmax = fabs(g);
                gd += g*g;
            }
        }
        if (gmax < tol_grad) { converged = 1; break; }

        /* Armijo backtracking on direction d = -grad, with displacement cap. */
        double t = t0;
        double Et = E;
        int ls_ok = 0;
        for (int ls = 0; ls < max_ls; ls++) {
            for (size_t v = 0; v < nv; v++) {
                for (int k = 0; k < 3; k++) {
                    double xv = X[v*3+k] - t*grad[v*3+k];
                    Xtry[v*3+k] = xv;
                }
                /* radial displacement clamp to X0 */
                double dx=Xtry[v*3]-X0[v*3], dy=Xtry[v*3+1]-X0[v*3+1], dz=Xtry[v*3+2]-X0[v*3+2];
                double d = sqrt(dx*dx+dy*dy+dz*dz);
                if (d > max_disp && d > 1e-20) {
                    double f = max_disp/d;
                    Xtry[v*3]   = X0[v*3]   + dx*f;
                    Xtry[v*3+1] = X0[v*3+1] + dy*f;
                    Xtry[v*3+2] = X0[v*3+2] + dz*f;
                }
            }
            Et = assemble_energy(Xtry, nv, faces, nf, interior, Asym, NULL, NULL);
            if (Et <= E - c1*t*gd) { ls_ok = 1; break; }
            t *= 0.5;
        }
        if (!ls_ok) break;                           /* stuck */
        memcpy(X, Xtry, nv*3*sizeof(double));
        double improve = E - Et;
        E = assemble_energy(X, nv, faces, nf, interior, Asym, NULL, uvec); /* refresh uvec */
        if (improve < tol_e*(fabs(E)+1.0)) { converged = 1; iter++; break; }
    }

    /* write back + measure displacement */
    double maxd = 0.0;
    for (size_t v = 0; v < nv; v++) {
        double dx=X[v*3]-X0[v*3], dy=X[v*3+1]-X0[v*3+1], dz=X[v*3+2]-X0[v*3+2];
        double d = sqrt(dx*dx+dy*dy+dz*dz);
        if (d > maxd) maxd = d;
        verts[v*3+0]=(float)X[v*3+0]; verts[v*3+1]=(float)X[v*3+1]; verts[v*3+2]=(float)X[v*3+2];
    }

    /* refresh vertex normals (angle-weighted face normals) if requested */
    if (vert_normals) {
        for (size_t i = 0; i < nv*3; i++) vert_normals[i] = 0.0f;
        for (size_t t = 0; t < nf; t++) {
            int32_t f[3] = { faces[t*3+0], faces[t*3+1], faces[t*3+2] };
            const double *Xa=&X[f[0]*3], *Xb=&X[f[1]*3], *Xc=&X[f[2]*3];
            double e1[3],e2[3],n[3]; sub3(Xb,Xa,e1); sub3(Xc,Xa,e2); cross3(e1,e2,n);
            double a2=sqrt(dot3(n,n)); if (a2<1e-20) continue;
            double N[3]={n[0]/a2,n[1]/a2,n[2]/a2};
            double th[3];
            { double p[3],q[3]; sub3(Xb,Xa,p); sub3(Xc,Xa,q); th[0]=corner_angle3(p,q); }
            { double p[3],q[3]; sub3(Xc,Xb,p); sub3(Xa,Xb,q); th[1]=corner_angle3(p,q); }
            { double p[3],q[3]; sub3(Xa,Xc,p); sub3(Xb,Xc,q); th[2]=corner_angle3(p,q); }
            for (int s=0;s<3;s++) for(int k=0;k<3;k++)
                vert_normals[f[s]*3+k] += (float)(th[s]*N[k]);
        }
        for (size_t v=0; v<nv; v++) {
            double nx=vert_normals[v*3],ny=vert_normals[v*3+1],nz=vert_normals[v*3+2];
            double ln=sqrt(nx*nx+ny*ny+nz*nz);
            if (ln>1e-20){vert_normals[v*3]=(float)(nx/ln);vert_normals[v*3+1]=(float)(ny/ln);vert_normals[v*3+2]=(float)(nz/ln);}
        }
    }

    if (out) { out->e_initial=e_init; out->e_final=E; out->iters=iter; out->max_disp=maxd; out->converged=converged; }
    if (opts && opts->verbose)
        fprintf(stderr, "  [develop] E %.4e -> %.4e in %d iters, max_disp=%.3f vox\n",
                e_init, E, iter, maxd);
    Arena_restore(arena, mark);
    return 0;
}

/* Per-vertex energy only (no descent). Reuses the same assemble_energy that
 * Develop_optimize runs, so the splitter's seam detector and the optimizer share
 * one definition of the Crane covariance energy. */
int Develop_vertex_energy(Arena_T arena, const float *verts, size_t nv,
                          const int32_t *faces, size_t nf, double *lam_out)
{
    assert(arena);
    if (!verts || !faces || !lam_out) return -1;
    for (size_t i = 0; i < nv; i++) lam_out[i] = 0.0;
    if (nv < 3 || nf < 1) return 0;                 /* trivial: all developable */

    Arena_Mark mark = Arena_save(arena);
    double *X = (double *)ARENA_ALLOC(arena, (long)(nv*3*sizeof(double)));
    for (size_t i = 0; i < nv*3; i++) X[i] = (double)verts[i];

    char *is_b = (char *)ARENA_ALLOC(arena, (long)nv);
    mark_boundary(arena, faces, nf, nv, is_b);
    char *interior = (char *)ARENA_ALLOC(arena, (long)nv);
    for (size_t v = 0; v < nv; v++) interior[v] = (char)!is_b[v];

    /* assemble_energy writes lam[v] only for interior verts; boundary stays 0
     * from the zero-fill above. uvec=NULL -> eigenvectors are not needed here. */
    double *Asym = (double *)ARENA_ALLOC(arena, (long)(nv*6*sizeof(double)));
    (void)assemble_energy(X, nv, faces, nf, interior, Asym, lam_out, NULL);

    Arena_restore(arena, mark);
    return 0;
}

/* ============================================================================
 * Self-test: FD gradient check + monotone energy decrease on a bumped grid.
 * ==========================================================================*/

static void build_bump(Arena_T arena, int N, double amp,
                       float **ov, size_t *onv, int32_t **of, size_t *onf)
{
    size_t nv=(size_t)N*N, nf=(size_t)(N-1)*(N-1)*2;
    float   *v=(float*)ARENA_ALLOC(arena,(long)(nv*3*sizeof(float)));
    int32_t *f=(int32_t*)ARENA_ALLOC(arena,(long)(nf*3*sizeof(int32_t)));
    for (int j=0;j<N;j++) for(int i=0;i<N;i++){
        double x=(double)i/(N-1), y=(double)j/(N-1);
        size_t idx=(size_t)j*N+i;
        v[idx*3+0]=(float)(10.0*x);
        v[idx*3+1]=(float)(10.0*y);
        /* smooth interior bump, zero on the boundary */
        v[idx*3+2]=(float)(amp*sin(M_PI*x)*sin(M_PI*y));
    }
    size_t fi=0;
    for(int j=0;j<N-1;j++)for(int i=0;i<N-1;i++){
        int32_t a=(int32_t)((size_t)j*N+i),b=(int32_t)((size_t)j*N+i+1);
        int32_t c=(int32_t)((size_t)(j+1)*N+i),d=(int32_t)((size_t)(j+1)*N+i+1);
        f[fi*3+0]=a;f[fi*3+1]=b;f[fi*3+2]=c;fi++;
        f[fi*3+0]=b;f[fi*3+1]=d;f[fi*3+2]=c;fi++;
    }
    *ov=v;*onv=nv;*of=f;*onf=fi;
}

int Develop_selftest(void)
{
    int fails = 0;
    Arena_T arena = Arena_new();

    /* (1) FD gradient check on a bumped grid. */
    {
        float *v; int32_t *f; size_t nv, nf;
        build_bump(arena, 9, 1.5, &v, &nv, &f, &nf);
        double *X = (double*)ARENA_ALLOC(arena,(long)(nv*3*sizeof(double)));
        for (size_t i=0;i<nv*3;i++) X[i]=(double)v[i];
        char *is_b=(char*)ARENA_ALLOC(arena,(long)nv);
        mark_boundary(arena,f,nf,nv,is_b);
        char *interior=(char*)ARENA_ALLOC(arena,(long)nv);
        for (size_t i=0;i<nv;i++) interior[i]=(char)!is_b[i];
        double *Asym=(double*)ARENA_ALLOC(arena,(long)(nv*6*sizeof(double)));
        double *uvec=(double*)ARENA_ALLOC(arena,(long)(nv*3*sizeof(double)));
        double *grad=(double*)ARENA_ALLOC(arena,(long)(nv*3*sizeof(double)));
        assemble_energy(X,nv,f,nf,interior,Asym,NULL,uvec);
        for (size_t i=0;i<nv*3;i++) grad[i]=0.0;
        accumulate_grad(X,nv,f,nf,interior,uvec,grad);

        double h=1e-6, worst=0.0, wfd=0, wan=0;
        int checked=0;
        for (size_t vtx=0; vtx<nv; vtx++) {
            if (!interior[vtx]) continue;
            for (int k=0;k<3;k++) {
                double save=X[vtx*3+k];
                X[vtx*3+k]=save+h; double Ep=assemble_energy(X,nv,f,nf,interior,Asym,NULL,NULL);
                X[vtx*3+k]=save-h; double Em=assemble_energy(X,nv,f,nf,interior,Asym,NULL,NULL);
                X[vtx*3+k]=save;
                double fd=(Ep-Em)/(2*h), an=grad[vtx*3+k];
                /* only judge components with a meaningful gradient (skip near-flat
                 * vertices where lambda_min is degenerate and the eigvec arbitrary) */
                if (fabs(fd) < 1e-2 && fabs(an) < 1e-2) continue;
                double err=fabs(fd-an)/(fabs(an)+fabs(fd)+1e-9);
                if (err>worst) { worst=err; wfd=fd; wan=an; }
                checked++;
            }
        }
        int ok = (worst < 1e-3);
        fprintf(stderr, "  [develop FD] checked=%d worst_rel_err=%.3e (fd=%.5f an=%.5f) -> %s\n",
                checked, worst, wfd, wan, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    /* (2) Optimizer reduces energy on the bump, respects displacement cap. */
    {
        float *v; int32_t *f; size_t nv, nf;
        build_bump(arena, 11, 2.0, &v, &nv, &f, &nf);
        DevelopOpts o; memset(&o,0,sizeof(o));
        o.max_iters=60; o.max_disp_vox=3.0; o.pin_boundary=1;
        DevelopStats st;
        int rc=Develop_optimize(arena, v, nv, f, nf, NULL, NULL, &o, &st);
        int ok = (rc==0) && (st.e_final < st.e_initial) && (st.max_disp <= 3.0+1e-6);
        fprintf(stderr, "  [develop opt] E %.4e -> %.4e disp=%.3f -> %s\n",
                st.e_initial, st.e_final, st.max_disp, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    /* (3) Flat plane: already developable -> E ~ 0, ~no movement. */
    {
        float *v; int32_t *f; size_t nv, nf;
        build_bump(arena, 8, 0.0, &v, &nv, &f, &nf);   /* amp 0 = flat */
        DevelopOpts o; memset(&o,0,sizeof(o)); o.max_iters=20;
        DevelopStats st;
        Develop_optimize(arena, v, nv, f, nf, NULL, NULL, &o, &st);
        int ok = (st.e_initial < 1e-9) && (st.max_disp < 1e-6);
        fprintf(stderr, "  [develop flat] E0=%.3e disp=%.3e -> %s\n",
                st.e_initial, st.max_disp, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    /* (4) pin_mask holds vertices fixed. */
    {
        float *v; int32_t *f; size_t nv, nf;
        build_bump(arena, 7, 1.0, &v, &nv, &f, &nf);
        uint8_t *pin=(uint8_t*)ARENA_ALLOC(arena,(long)nv);
        for (size_t i=0;i<nv;i++) pin[i]=1;            /* pin everything -> no-op */
        float *v0=(float*)ARENA_ALLOC(arena,(long)(nv*3*sizeof(float)));
        memcpy(v0,v,nv*3*sizeof(float));
        DevelopOpts o; memset(&o,0,sizeof(o)); o.max_iters=20;
        DevelopStats st;
        Develop_optimize(arena, v, nv, f, nf, pin, NULL, &o, &st);
        int ok = (st.max_disp < 1e-9);
        for (size_t i=0;i<nv*3;i++) if (v[i]!=v0[i]) ok=0;
        fprintf(stderr, "  [develop pin] disp=%.3e -> %s\n", st.max_disp, ok?"ok":"FAIL");
        if (!ok) fails++;
    }

    if (fails==0) fprintf(stderr, "[develop selftest] ok\n");
    Arena_dispose(&arena);
    return fails;
}
