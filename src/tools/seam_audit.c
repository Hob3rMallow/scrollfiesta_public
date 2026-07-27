/*
 * seam_audit.c -- detect invalid/overlapping geometry near a cube-boundary seam.
 *
 * The manifold audit in grid_weld only counts edge multiplicity (>2 incident
 * faces). It is BLIND to triangles that overlap or interpenetrate without
 * sharing an edge -- exactly the "invalid/overlapping geometry on one side of
 * the bridge" a thrashing BPA seam-weld leaves behind.
 *
 * The dominant failure mode is NEAR-COINCIDENT geometry: a bridge triangle and a
 * sheet triangle lying in ~the same plane, a fraction of a voxel apart -- a
 * doubled surface that Z-FIGHTS while rendering. That is NOT clean "stabbing"
 * interpenetration, and a naive triangle-triangle sign test misses it (parallel
 * triangles with a small gap test as "all on one side -> disjoint"). So the
 * primary detector here is:
 *
 *   NEAR-COINCIDENT OVERLAP -- two non-adjacent triangles whose normals are
 *     near-parallel (within --angle), whose perpendicular separation is small
 *     (within --gap voxels), and whose in-plane projections overlap. The gap
 *     distribution is reported so you can tell true Z-fight (sub-voxel) from a
 *     merely-close pair.
 *
 * plus, for completeness, the non-parallel case:
 *
 *   INTERPENETRATION -- two non-parallel non-adjacent triangles that actually
 *     cross (Moller 1997 triangle-triangle interval test, double precision).
 *
 * Pairs that share any vertex (legitimate fan/edge adjacency) are excluded, so a
 * hit is a real defect. Parallel sheets farther apart than --gap (e.g. adjacent
 * scroll wraps at >=7 vox) are deliberately NOT flagged.
 *
 * Coordinate convention matches the pipeline OBJs: "v c0 c1 c2" with c0=z, c1=y,
 * c2=x. Plane axis 0=z, 1=y, 2=x (so the x=2944 seam is axis 2, coord 2944).
 *
 * Usage:
 *   seam_audit <mesh.obj> [--cube <size=128>] [--plane <axis> <coord>]
 *              [--band <vox=4>] [--gap <vox=1.0>] [--angle <deg=20>]
 *              [--out <bad.obj>]
 *   seam_audit --selftest
 *
 * Exit: 0 = clean (or selftest pass), 1 = offenders found (or selftest fail),
 *       2 = usage/IO error.
 *
 * Build: standalone (seam_audit.vcxproj compiles only this file).
 */
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ===================================================================
 * Vector helpers (double precision; coords can be ~thousands).
 * =================================================================== */
static void vsub(double r[3], const double a[3], const double b[3])
{ r[0]=a[0]-b[0]; r[1]=a[1]-b[1]; r[2]=a[2]-b[2]; }
static void vcross(double r[3], const double a[3], const double b[3])
{ r[0]=a[1]*b[2]-a[2]*b[1]; r[1]=a[2]*b[0]-a[0]*b[2]; r[2]=a[0]*b[1]-a[1]*b[0]; }
static double vdot(const double a[3], const double b[3])
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

/* Unit normal of triangle (p0,p1,p2); 0 if degenerate (area ~ 0). */
static int tri_unit_normal(const double p0[3], const double p1[3],
                           const double p2[3], double n[3])
{
    double e1[3], e2[3];
    vsub(e1,p1,p0); vsub(e2,p2,p0);
    vcross(n,e1,e2);
    double len = sqrt(vdot(n,n));
    if (len < 1e-12) return 0;
    n[0]/=len; n[1]/=len; n[2]/=len;
    return 1;
}

/* ===================================================================
 * 2D primitives for the in-plane overlap test (project out the dominant
 * normal axis -- a parallel projection that preserves overlap topology for
 * near-parallel triangles).
 * =================================================================== */
static double cross2(double ax,double ay,double bx,double by){ return ax*by - ay*bx; }

/* Strict-interior point in triangle: all three edge cross-products the same
 * sign with magnitude above eps (so vertices lying exactly on an edge -- a bare
 * touch -- do NOT count as overlap). */
static int pt_in_tri2d(double px,double py,
                       double ax,double ay,double bx,double by,double cx,double cy)
{
    double d1 = cross2(bx-ax, by-ay, px-ax, py-ay);
    double d2 = cross2(cx-bx, cy-by, px-bx, py-by);
    double d3 = cross2(ax-cx, ay-cy, px-cx, py-cy);
    const double e = 1e-9;
    if (fabs(d1) < e || fabs(d2) < e || fabs(d3) < e) return 0;  /* on an edge */
    int neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    int pos = (d1 > 0) || (d2 > 0) || (d3 > 0);
    return !(neg && pos);   /* all strictly same side -> interior */
}

/* Proper segment-segment crossing (interiors only). */
static int seg_cross2d(double p1x,double p1y,double p2x,double p2y,
                       double q1x,double q1y,double q2x,double q2y)
{
    double d1x=p2x-p1x, d1y=p2y-p1y, d2x=q2x-q1x, d2y=q2y-q1y;
    double denom = cross2(d1x,d1y,d2x,d2y);
    if (fabs(denom) < 1e-12) return 0;          /* parallel/collinear */
    double sx=q1x-p1x, sy=q1y-p1y;
    double t = cross2(sx,sy,d2x,d2y) / denom;
    double s = cross2(sx,sy,d1x,d1y) / denom;
    return (t>1e-9 && t<1.0-1e-9 && s>1e-9 && s<1.0-1e-9);
}

/* In-plane area overlap of two near-parallel triangles, projected with n. */
static int inplane_overlap(const double n[3],
                           const double v0[3],const double v1[3],const double v2[3],
                           const double u0[3],const double u1[3],const double u2[3])
{
    double ax=fabs(n[0]),ay=fabs(n[1]),az=fabs(n[2]);
    int i0,i1;
    if (ax>ay && ax>az) { i0=1; i1=2; }
    else if (ay>az)     { i0=0; i1=2; }
    else                { i0=0; i1=1; }
    const double *T1[3]={v0,v1,v2}, *T2[3]={u0,u1,u2};
    double a[3][2], b[3][2];
    for (int k=0;k<3;k++){ a[k][0]=T1[k][i0]; a[k][1]=T1[k][i1];
                           b[k][0]=T2[k][i0]; b[k][1]=T2[k][i1]; }
    /* edge crossings */
    for (int e=0;e<3;e++)
        for (int f=0;f<3;f++)
            if (seg_cross2d(a[e][0],a[e][1],a[(e+1)%3][0],a[(e+1)%3][1],
                            b[f][0],b[f][1],b[(f+1)%3][0],b[(f+1)%3][1]))
                return 1;
    /* containment (one triangle's vertex strictly inside the other) */
    for (int k=0;k<3;k++) {
        if (pt_in_tri2d(a[k][0],a[k][1], b[0][0],b[0][1],b[1][0],b[1][1],b[2][0],b[2][1])) return 1;
        if (pt_in_tri2d(b[k][0],b[k][1], a[0][0],a[0][1],a[1][0],a[1][1],a[2][0],a[2][1])) return 1;
    }
    return 0;
}

/* ===================================================================
 * Non-parallel interpenetration: Moller 1997 triangle-triangle interval test.
 * (Used only when the triangles are NOT near-parallel; the parallel/near-
 * coincident case is handled separately and more precisely below.)
 * =================================================================== */
static void isect_interval(const double pv[3], const double d[3], double *t0, double *t1)
{
    int i0,i1,i2;
    double d0d1=d[0]*d[1], d0d2=d[0]*d[2];
    if (d0d1>0.0)      { i0=2; i1=0; i2=1; }
    else if (d0d2>0.0) { i0=1; i1=0; i2=2; }
    else               { i0=0; i1=1; i2=2; }
    double a = pv[i1] + (pv[i0]-pv[i1]) * (d[i1]/(d[i1]-d[i0]));
    double b = pv[i2] + (pv[i0]-pv[i2]) * (d[i2]/(d[i2]-d[i0]));
    if (a<=b){ *t0=a; *t1=b; } else { *t0=b; *t1=a; }
}

static int moller_interpenetrate(const double v0[3],const double v1[3],const double v2[3],
                                 const double u0[3],const double u1[3],const double u2[3])
{
    double n1[3], n2[3];
    if (!tri_unit_normal(v0,v1,v2,n1)) return 0;
    if (!tri_unit_normal(u0,u1,u2,n2)) return 0;
    double d2 = -vdot(n2,u0);
    double dv[3] = { vdot(n2,v0)+d2, vdot(n2,v1)+d2, vdot(n2,v2)+d2 };
    for (int k=0;k<3;k++) if (fabs(dv[k])<1e-7) dv[k]=0.0;
    if ((dv[0]>0&&dv[1]>0&&dv[2]>0)||(dv[0]<0&&dv[1]<0&&dv[2]<0)) return 0;
    if (dv[0]==0&&dv[1]==0&&dv[2]==0) return 0;        /* coplanar -> parallel path */
    double d1 = -vdot(n1,v0);
    double du[3] = { vdot(n1,u0)+d1, vdot(n1,u1)+d1, vdot(n1,u2)+d1 };
    for (int k=0;k<3;k++) if (fabs(du[k])<1e-7) du[k]=0.0;
    if ((du[0]>0&&du[1]>0&&du[2]>0)||(du[0]<0&&du[1]<0&&du[2]<0)) return 0;
    double D[3]; vcross(D,n1,n2);
    int idx=0; double m=fabs(D[0]);
    if (fabs(D[1])>m){ m=fabs(D[1]); idx=1; }
    if (fabs(D[2])>m){ idx=2; }
    double pv[3]={v0[idx],v1[idx],v2[idx]}, pu[3]={u0[idx],u1[idx],u2[idx]};
    double a0,a1,b0,b1;
    isect_interval(pv,dv,&a0,&a1);
    isect_interval(pu,du,&b0,&b1);
    if (a1 < b0+1e-9 || b1 < a0+1e-9) return 0;
    return 1;
}

/* ===================================================================
 * The pair test. Near-parallel -> slab + in-plane overlap (the Z-fight case);
 * else -> Moller interpenetration. cos_par = cos(angle threshold); gap_max is
 * the max perpendicular separation (vox) to consider near-coincident.
 * =================================================================== */
#define HIT_NONE    0
#define HIT_OVERLAP 1   /* near-coincident, parallel, in-plane overlap */
#define HIT_STAB    2   /* non-parallel interpenetration */

static int tri_pair_test(const double v0[3],const double v1[3],const double v2[3],
                         const double u0[3],const double u1[3],const double u2[3],
                         double cos_par, double gap_max, double *out_gap)
{
    double n1[3], n2[3];
    if (!tri_unit_normal(v0,v1,v2,n1)) return HIT_NONE;
    if (!tri_unit_normal(u0,u1,u2,n2)) return HIT_NONE;

    if (fabs(vdot(n1,n2)) >= cos_par) {
        /* near-parallel: perpendicular separation of the FAR triangle from this
         * triangle's plane (whole far-tri must sit inside the slab -> a genuine
         * doubled surface, not a tilted partial clip). */
        double d1 = -vdot(n1,v0);
        double s0=fabs(vdot(n1,u0)+d1), s1=fabs(vdot(n1,u1)+d1), s2=fabs(vdot(n1,u2)+d1);
        double gmax = s0>s1?(s0>s2?s0:s2):(s1>s2?s1:s2);
        if (gmax > gap_max) return HIT_NONE;        /* too far apart (e.g. wraps) */
        if (inplane_overlap(n1, v0,v1,v2, u0,u1,u2)) {
            double gmin = s0<s1?(s0<s2?s0:s2):(s1<s2?s1:s2);
            if (out_gap) *out_gap = 0.5*(gmin+gmax);
            return HIT_OVERLAP;
        }
        return HIT_NONE;
    }
    return moller_interpenetrate(v0,v1,v2,u0,u1,u2) ? HIT_STAB : HIT_NONE;
}

/* ===================================================================
 * Minimal OBJ reader (self-contained). First three floats of each "v" line;
 * 0-based indices of each triangular "f" line (tolerates v/vt/vn slashes).
 * =================================================================== */
typedef struct { float *v; size_t nv; int32_t *f; size_t nf; } Mesh;

static int obj_read(const char *path, Mesh *m)
{
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;
    size_t vcap=1024, fcap=1024;
    m->v=(float*)malloc(vcap*3*sizeof(float));
    m->f=(int32_t*)malloc(fcap*3*sizeof(int32_t));
    m->nv=0; m->nf=0;
    char line[512];
    while (fgets(line,sizeof line,fp)) {
        if (line[0]=='v' && line[1]==' ') {
            if (m->nv>=vcap){ vcap*=2; m->v=(float*)realloc(m->v,vcap*3*sizeof(float)); }
            double a=0,b=0,c=0;
            sscanf(line+2,"%lf %lf %lf",&a,&b,&c);
            m->v[m->nv*3+0]=(float)a; m->v[m->nv*3+1]=(float)b; m->v[m->nv*3+2]=(float)c;
            m->nv++;
        } else if (line[0]=='f' && line[1]==' ') {
            int ia=0,ib=0,ic=0;
            if (sscanf(line+2,"%d/%*d/%*d %d/%*d/%*d %d/%*d/%*d",&ia,&ib,&ic)==3 ||
                sscanf(line+2,"%d//%*d %d//%*d %d//%*d",&ia,&ib,&ic)==3 ||
                sscanf(line+2,"%d/%*d %d/%*d %d/%*d",&ia,&ib,&ic)==3 ||
                sscanf(line+2,"%d %d %d",&ia,&ib,&ic)==3) {
                if (m->nf>=fcap){ fcap*=2; m->f=(int32_t*)realloc(m->f,fcap*3*sizeof(int32_t)); }
                m->f[m->nf*3+0]=ia-1; m->f[m->nf*3+1]=ib-1; m->f[m->nf*3+2]=ic-1;
                m->nf++;
            }
        }
    }
    fclose(fp);
    return 0;
}
static void mesh_free(Mesh *m){ free(m->v); free(m->f); m->v=NULL; m->f=NULL; }

/* ===================================================================
 * Seam-plane detection (mirrors seam_weld.c::detect_planes).
 * =================================================================== */
#define SEAM_MIN_SIDE 2.0
typedef struct { int axis; double coord; } Plane;

static size_t detect_planes(const Mesh *m, double cube, double band, Plane *planes, size_t cap)
{
    double mn[3]={1e30,1e30,1e30}, mx[3]={-1e30,-1e30,-1e30};
    for (size_t v=0; v<m->nv; v++)
        for (int a=0;a<3;a++){ double c=m->v[v*3+(size_t)a]; if(c<mn[a])mn[a]=c; if(c>mx[a])mx[a]=c; }
    size_t np=0;
    for (int a=0;a<3 && np<cap;a++){
        long k0=(long)ceil(mn[a]/cube), k1=(long)floor(mx[a]/cube);
        for (long k=k0;k<=k1 && np<cap;k++){
            double c=(double)k*cube; int lo=0,hi=0;
            for (size_t v=0; v<m->nv; v++){
                double cc=m->v[v*3+(size_t)a];
                if (cc>c-band && cc<c-SEAM_MIN_SIDE) lo=1;
                if (cc>c+SEAM_MIN_SIDE && cc<c+band) hi=1;
                if (lo&&hi) break;
            }
            if (lo&&hi){ planes[np].axis=a; planes[np].coord=c; np++; }
        }
    }
    return np;
}

/* ===================================================================
 * Audit one plane: gather band triangles, sweep-prune on an in-plane axis,
 * AABB-reject, then run the pair test on non-adjacent pairs.
 * =================================================================== */
typedef struct { size_t fi; double lo, hi; } SeamTri;
static int cmp_seamtri(const void *pa,const void *pb)
{ double a=((const SeamTri*)pa)->lo, b=((const SeamTri*)pb)->lo; return a<b?-1:(a>b?1:0); }

static void tri_pts(const Mesh *m, size_t fi, double p0[3], double p1[3], double p2[3])
{
    int32_t a=m->f[fi*3+0], b=m->f[fi*3+1], c=m->f[fi*3+2];
    for (int k=0;k<3;k++){ p0[k]=m->v[(size_t)a*3+(size_t)k];
                           p1[k]=m->v[(size_t)b*3+(size_t)k];
                           p2[k]=m->v[(size_t)c*3+(size_t)k]; }
}
static int shared_count(const Mesh *m, size_t fa, size_t fb)
{
    int n=0;
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
        if (m->f[fa*3+i]==m->f[fb*3+j]) n++;
    return n;   /* 0=disjoint, 1=hinge(shared vertex), 2=edge-adjacent */
}

static double tri_area(const double p0[3], const double p1[3], const double p2[3])
{
    double e1[3], e2[3], c[3];
    vsub(e1,p1,p0); vsub(e2,p2,p0); vcross(c,e1,e2);
    return 0.5*sqrt(vdot(c,c));
}

/* Minimum altitude of a triangle = 2*area / longest edge. A sliver (near-
 * degenerate triangle) has a tiny min altitude even when its area is not quite
 * zero -- this is the quantity an edge-collapse / QEM pass keys on, and the
 * better "is this a sliver?" signal than area alone. Returns 0 for a point. */
static double tri_min_alt(const double p0[3], const double p1[3], const double p2[3])
{
    double e0[3], e1[3], e2[3];
    vsub(e0,p1,p0); vsub(e1,p2,p1); vsub(e2,p0,p2);
    double l0=sqrt(vdot(e0,e0)), l1=sqrt(vdot(e1,e1)), l2=sqrt(vdot(e2,e2));
    double lmax = l0>l1?(l0>l2?l0:l2):(l1>l2?l1:l2);
    if (lmax < 1e-12) return 0.0;
    return 2.0*tri_area(p0,p1,p2)/lmax;
}

/* Squared distance from point p to segment [a,b]; foot parameter t (0=a, 1=b)
 * returned via t_out (clamped to [0,1]). */
static double pt_seg_d2(const double p[3], const double a[3], const double b[3], double *t_out)
{
    double ab[3], ap[3];
    vsub(ab,b,a); vsub(ap,p,a);
    double denom = vdot(ab,ab);
    double t = denom > 1e-18 ? vdot(ap,ab)/denom : 0.0;
    if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
    double q[3] = { a[0]+t*ab[0], a[1]+t*ab[1], a[2]+t*ab[2] };
    double d[3]; vsub(d,p,q);
    if (t_out) *t_out = t;
    return vdot(d,d);
}

/* T-junction test: a NON-shared vertex of one triangle lies within `eps` of the
 * INTERIOR of an edge of the other (foot parameter strictly between the edge's
 * endpoints, so a shared corner does not count). This is the signature the bridge
 * leaves when it welds a vertex onto the middle of a coarse triangle's edge
 * rather than at its corner -- the two faces then read as "overlapping" to the
 * coplanar test even though they only touch along the split edge. */
static int is_tjunction(const Mesh *m, size_t fa, size_t fb, double eps)
{
    const double tedge = 0.02;             /* keep clear of the endpoints */
    const int32_t *F[2] = { &m->f[fa*3], &m->f[fb*3] };
    for (int pass=0; pass<2; pass++) {
        const int32_t *FP = F[pass], *FE = F[1-pass];   /* point-tri vs edge-tri */
        for (int i=0;i<3;i++) {
            int32_t vp = FP[i];
            if (vp==FE[0] || vp==FE[1] || vp==FE[2]) continue;   /* shared corner */
            double p[3] = { m->v[(size_t)vp*3+0], m->v[(size_t)vp*3+1], m->v[(size_t)vp*3+2] };
            for (int e=0;e<3;e++) {
                int32_t ea=FE[e], eb=FE[(e+1)%3];
                double a[3] = { m->v[(size_t)ea*3+0], m->v[(size_t)ea*3+1], m->v[(size_t)ea*3+2] };
                double b[3] = { m->v[(size_t)eb*3+0], m->v[(size_t)eb*3+1], m->v[(size_t)eb*3+2] };
                double t; double d2 = pt_seg_d2(p,a,b,&t);
                if (d2 <= eps*eps && t > tedge && t < 1.0-tedge) return 1;
            }
        }
    }
    return 0;
}

/* A face is "cross-region" if it uses verts from both sides of the --split
 * boundary -- i.e. it spans the two welded cubes, so it is a bridge face the
 * weld created (not original cube geometry). split<0 -> unknown. */
static int face_cross_region(const Mesh *m, size_t f, long split)
{
    if (split < 0) return 0;
    int below=0, above=0;
    for (int k=0;k<3;k++){ long v=m->f[f*3+k]; if (v<split) below=1; else above=1; }
    return below && above;
}

/* True if face f comes within `band` of any detected seam plane (the weld zone). */
static int face_near_seam(const Mesh *m, size_t f, const Plane *planes,
                          size_t np, double band)
{
    double p0[3],p1[3],p2[3]; tri_pts(m,f,p0,p1,p2);
    for (size_t p=0;p<np;p++){
        int ax=planes[p].axis; double co=planes[p].coord;
        double mn=p0[ax], mx=p0[ax];
        if(p1[ax]<mn)mn=p1[ax]; if(p1[ax]>mx)mx=p1[ax];
        if(p2[ax]<mn)mn=p2[ax]; if(p2[ax]>mx)mx=p2[ax];
        if (mx >= co-band && mn <= co+band) return 1;
    }
    return 0;
}

/* Label a vertex by the --split boundary: 'A' below split, 'B' at/above.
 * (split is the first cube's vertex count for a 2-component concat; lets the
 * --explain dump show whether a coincident pair is cross-component.) */
static char region_of(int32_t vi, long split)
{ return split < 0 ? '?' : ((long)vi < split ? 'A' : 'B'); }

/* Nearest-vertex distance between the two faces (min over the 9 vertex pairs):
 * if a doubled face is just a near-coincident copy, this is ~the vertex jitter. */
static double min_vertex_gap(const Mesh *m, size_t fa, size_t fb)
{
    double best = 1e30;
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
        int32_t a = m->f[fa*3+i], b = m->f[fb*3+j];
        double dz=(double)m->v[(size_t)a*3+0]-m->v[(size_t)b*3+0];
        double dy=(double)m->v[(size_t)a*3+1]-m->v[(size_t)b*3+1];
        double dx=(double)m->v[(size_t)a*3+2]-m->v[(size_t)b*3+2];
        double d=sqrt(dz*dz+dy*dy+dx*dx);
        if (d<best) best=d;
    }
    return best;
}

typedef struct { size_t overlap, stab, fold; double gmin, gsum, gmax; } Stats;

static size_t audit_plane(const Mesh *m, const Plane *pl, double band,
                          double cos_par, double gap_max, int min_excl,
                          int explain, long split, uint8_t *bad, Stats *st_acc)
{
    int axis=pl->axis, sweep=(axis+1)%3;
    SeamTri *st=(SeamTri*)malloc((m->nf?m->nf:1)*sizeof(SeamTri));
    size_t k=0;
    for (size_t f=0; f<m->nf; f++){
        double p0[3],p1[3],p2[3]; tri_pts(m,f,p0,p1,p2);
        double amin=p0[axis],amax=p0[axis];
        if(p1[axis]<amin)amin=p1[axis]; if(p1[axis]>amax)amax=p1[axis];
        if(p2[axis]<amin)amin=p2[axis]; if(p2[axis]>amax)amax=p2[axis];
        if (amax < pl->coord-band || amin > pl->coord+band) continue;
        double lo=p0[sweep],hi=p0[sweep];
        if(p1[sweep]<lo)lo=p1[sweep]; if(p1[sweep]>hi)hi=p1[sweep];
        if(p2[sweep]<lo)lo=p2[sweep]; if(p2[sweep]>hi)hi=p2[sweep];
        st[k].fi=f; st[k].lo=lo; st[k].hi=hi; k++;
    }
    qsort(st,k,sizeof(SeamTri),cmp_seamtri);

    size_t pairs=0;
    for (size_t i=0;i<k;i++){
        double va0[3],va1[3],va2[3]; tri_pts(m,st[i].fi,va0,va1,va2);
        double amn[3],amx[3];
        for (int c=0;c<3;c++){ amn[c]=va0[c]; amx[c]=va0[c];
            if(va1[c]<amn[c])amn[c]=va1[c]; if(va1[c]>amx[c])amx[c]=va1[c];
            if(va2[c]<amn[c])amn[c]=va2[c]; if(va2[c]>amx[c])amx[c]=va2[c]; }
        for (size_t j=i+1; j<k && st[j].lo <= st[i].hi + gap_max; j++){
            int sc = shared_count(m, st[i].fi, st[j].fi);
            /* sc==1 hinge (fan corner): skip unless --hinge. sc==2 edge-adjacent:
             * KEEP -- a fold-back double (two triangles sharing an edge, interiors
             * overlapping) is a self-intersection the blanket shared-vertex skip
             * used to miss. Only its near-parallel HIT_OVERLAP is real; the Moller
             * path is rejected below (it trivially fires on the shared edge). */
            if (sc >= min_excl && sc != 2) continue;
            double vb0[3],vb1[3],vb2[3]; tri_pts(m,st[j].fi,vb0,vb1,vb2);
            double bmn[3],bmx[3]; int sep=0;
            for (int c=0;c<3;c++){ bmn[c]=vb0[c]; bmx[c]=vb0[c];
                if(vb1[c]<bmn[c])bmn[c]=vb1[c]; if(vb1[c]>bmx[c])bmx[c]=vb1[c];
                if(vb2[c]<bmn[c])bmn[c]=vb2[c]; if(vb2[c]>bmx[c])bmx[c]=vb2[c]; }
            for (int c=0;c<3;c++) if (amx[c] < bmn[c]-gap_max || bmx[c] < amn[c]-gap_max){ sep=1; break; }
            if (sep) continue;
            double gap=0.0;
            int r = tri_pair_test(va0,va1,va2, vb0,vb1,vb2, cos_par, gap_max, &gap);
            if (sc >= 1 && r == HIT_STAB) continue;   /* shared geom: Moller bogus */
            if (r){
                pairs++;
                bad[st[i].fi]=1; bad[st[j].fi]=1;
                if (sc >= 1){
                    st_acc->fold++;                    /* edge/hinge fold-back double */
                } else if (r==HIT_OVERLAP){
                    st_acc->overlap++;
                    if (gap<st_acc->gmin) st_acc->gmin=gap;
                    if (gap>st_acc->gmax) st_acc->gmax=gap;
                    st_acc->gsum+=gap;
                } else st_acc->stab++;
                if (explain){
                    size_t fa=st[i].fi, fb=st[j].fi;
                    int32_t *A=&m->f[fa*3], *B=&m->f[fb*3];
                    double na[3], nb[3];
                    tri_unit_normal(va0,va1,va2,na); tri_unit_normal(vb0,vb1,vb2,nb);
                    printf("  PAIR %s  plane-gap=%.4f  min-vtx-gap=%.4f  shared=%d  n.dot=%.4f\n",
                           r==HIT_OVERLAP?"OVERLAP":"STAB", gap,
                           min_vertex_gap(m,fa,fb), shared_count(m,fa,fb), vdot(na,nb));
                    printf("       minalt A=%.4f B=%.4f   T-junction=%s\n",
                           tri_min_alt(va0,va1,va2), tri_min_alt(vb0,vb1,vb2),
                           is_tjunction(m,fa,fb,0.25) ? "YES (vtx on opposite edge)" : "no");
                    printf("    A f%zu v[%d%c %d%c %d%c] area=%.4f\n", fa,
                           A[0],region_of(A[0],split), A[1],region_of(A[1],split),
                           A[2],region_of(A[2],split), tri_area(va0,va1,va2));
                    printf("       (%.3f %.3f %.3f)(%.3f %.3f %.3f)(%.3f %.3f %.3f)\n",
                           va0[0],va0[1],va0[2], va1[0],va1[1],va1[2], va2[0],va2[1],va2[2]);
                    printf("    B f%zu v[%d%c %d%c %d%c] area=%.4f\n", fb,
                           B[0],region_of(B[0],split), B[1],region_of(B[1],split),
                           B[2],region_of(B[2],split), tri_area(vb0,vb1,vb2));
                    printf("       (%.3f %.3f %.3f)(%.3f %.3f %.3f)(%.3f %.3f %.3f)\n",
                           vb0[0],vb0[1],vb0[2], vb1[0],vb1[1],vb1[2], vb2[0],vb2[1],vb2[2]);
                }
            }
        }
    }
    free(st);
    return pairs;
}

/* ===================================================================
 * Self-test: synthetic pairs with known verdicts, focused on the near-coplanar
 * regime (parallel + small gap) the real defect lives in.
 * =================================================================== */
static int selftest(void)
{
    int fail=0;
    double cp = cos(20.0*M_PI/180.0), g;
    struct { const char *name; double v[3][3], u[3][3]; int expect; double gap; double cos_par; } T[] = {
      /* 1. exact coplanar overlap */
      {"coplanar overlap", {{0,0,0},{4,0,0},{0,4,0}}, {{1,1,0},{5,1,0},{1,5,0}}, HIT_OVERLAP, 1.0, 0},
      /* 2. near-coincident, 0.1 vox gap, overlapping in-plane (THE z-fight case) */
      {"near-coincident 0.1v", {{0,0,0},{4,0,0},{0,4,0}}, {{1,1,0.1},{5,1,0.1},{1,5,0.1}}, HIT_OVERLAP, 1.0, 0},
      /* 3. parallel but 5 vox apart (e.g. adjacent wrap) -> NOT flagged */
      {"parallel 5v apart", {{0,0,0},{4,0,0},{0,4,0}}, {{1,1,5},{5,1,5},{1,5,5}}, HIT_NONE, 1.0, 0},
      /* 4. coplanar but disjoint in-plane */
      {"coplanar disjoint", {{0,0,0},{1,0,0},{0,1,0}}, {{10,10,0},{11,10,0},{10,11,0}}, HIT_NONE, 1.0, 0},
      /* 5. coplanar edge-touch only */
      {"coplanar edge-touch", {{0,0,0},{4,0,0},{0,4,0}}, {{4,0,0},{4,4,0},{0,4,0}}, HIT_NONE, 1.0, 0},
      /* 6. non-parallel true stab */
      {"interpenetration", {{0,0,0},{4,0,0},{0,4,0}}, {{1,1,-2},{1,1,2},{3,1,0}}, HIT_STAB, 1.0, 0},
      /* 7. near-coincident but 0.6v gap with gap_max 0.5 -> NOT flagged */
      {"gap 0.6 > max 0.5", {{0,0,0},{4,0,0},{0,4,0}}, {{1,1,0.6},{5,1,0.6},{1,5,0.6}}, HIT_NONE, 0.5, 0},
      /* 8. EDGE-ADJACENT FOLD-BACK: T2 shares edge (0,0,0)-(4,0,0); its apex
       *    (0.5,0.5) lands strictly inside T1 (same side, ~coplanar) -> a fold-back
       *    double whose interiors overlap. This is the seam-bridge self-intersect. */
      {"edge-adj fold-back", {{0,0,0},{4,0,0},{0,4,0}}, {{0,0,0},{4,0,0},{0.5,0.5,0.05}}, HIT_OVERLAP, 1.0, 0},
      /* 9. EDGE-ADJACENT FAN: same shared edge but apex on the OPPOSITE side -> no
       *    AREA overlap (touch only along the edge) -> a legit crease, NOT flagged. */
      {"edge-adj flat fan", {{0,0,0},{4,0,0},{0,4,0}}, {{0,0,0},{4,0,0},{1,-4,0}}, HIT_NONE, 1.0, 0},
    };
    for (size_t i=0;i<sizeof T/sizeof T[0];i++){
        double gmax = T[i].gap;
        g=0;
        int r = tri_pair_test(T[i].v[0],T[i].v[1],T[i].v[2],
                              T[i].u[0],T[i].u[1],T[i].u[2], cp, gmax, &g);
        const char *kind = r==HIT_OVERLAP?"OVERLAP":(r==HIT_STAB?"STAB":"none");
        printf("  [selftest] %-22s -> %s (gap=%.2f) expect %s\n",
               T[i].name, kind, g,
               T[i].expect==HIT_OVERLAP?"OVERLAP":(T[i].expect==HIT_STAB?"STAB":"none"));
        if (r != T[i].expect) fail++;
    }

    /* tri_min_alt: a fat triangle has a large min altitude; a sliver tiny. */
    {
        double eq0[3]={0,0,0}, eq1[3]={1,0,0}, eq2[3]={0,1,0};
        double sl0[3]={0,0,0}, sl1[3]={10,0,0}, sl2[3]={5,0.02,0};
        double ma_eq=tri_min_alt(eq0,eq1,eq2), ma_sl=tri_min_alt(sl0,sl1,sl2);
        printf("  [selftest] tri_min_alt  fat=%.4f sliver=%.4f (expect fat>0.5, sliver<0.10)\n",
               ma_eq, ma_sl);
        if (!(ma_eq > 0.5)) fail++;
        if (!(ma_sl < 0.10)) fail++;
    }

    /* is_tjunction: B's apex sits on the interior of A's edge -> YES; a disjoint
     * far triangle -> no. */
    {
        float V[9*3] = {
            0.0f,0.0f,0.0f,   10.0f,0.0f,0.0f,   0.0f,10.0f,0.0f,   /* 0,1,2  tri A   */
            5.0f,0.0f,0.0f,    5.0f,-3.0f,0.0f,   8.0f,-3.0f,0.0f,   /* 3,4,5  tri B   */
            100.0f,100.0f,0.0f,101.0f,100.0f,0.0f,100.0f,101.0f,0.0f /* 6,7,8  tri C   */
        };
        int32_t F[3*3] = { 0,1,2,  3,4,5,  6,7,8 };
        Mesh tm; tm.v=V; tm.nv=9; tm.f=F; tm.nf=3;
        int tj_AB=is_tjunction(&tm,0,1,0.25), tj_AC=is_tjunction(&tm,0,2,0.25);
        printf("  [selftest] is_tjunction  A-B=%d (expect 1)  A-C=%d (expect 0)\n", tj_AB, tj_AC);
        if (tj_AB!=1) fail++;
        if (tj_AC!=0) fail++;
    }

    printf("SEAM AUDIT SELFTEST %s\n", fail?"FAILED":"PASSED");
    return fail?1:0;
}

/* ===================================================================
 * Whole-mesh degeneracy census (--degen). Counts degenerate (near-zero-area)
 * and sliver (tiny min-altitude) triangles, split by whether they sit in the
 * weld zone (within --band of a seam plane) and whether they are bridge faces
 * (cross-region; needs --split). Answers "is the bridge welding against /
 * creating zero-area triangles?" and sizes the edge-collapse cleanup.
 * Returns the total degenerate+sliver count.
 * =================================================================== */
static int degen_scan(const Mesh *m, const Plane *planes, size_t np, double band,
                      long split, const char *out_path)
{
    const double AREA_EPS = 1e-3;   /* vox^2: effectively zero area  */
    const double ALT_EPS  = 0.10;   /* vox:   thinner than 1/10 voxel -> sliver */
    size_t n_degen=0, n_sliver=0, n_degen_seam=0, n_sliver_seam=0;
    size_t n_degen_bridge=0, n_sliver_bridge=0;
    uint8_t *flag = (uint8_t *)calloc(m->nf?m->nf:1,1);

    enum { WORST=32 };
    size_t widx[WORST]; double walt[WORST]; size_t nworst=0;

    for (size_t f=0; f<m->nf; f++){
        double p0[3],p1[3],p2[3]; tri_pts(m,f,p0,p1,p2);
        double ar = tri_area(p0,p1,p2);
        double ma = tri_min_alt(p0,p1,p2);
        int degen  = (ar < AREA_EPS);
        int sliver = (!degen && ma < ALT_EPS);
        if (!degen && !sliver) continue;
        flag[f]=1;
        int seam   = face_near_seam(m,f,planes,np,band);
        int bridge = face_cross_region(m,f,split);
        if (degen){ n_degen++;  if(seam)n_degen_seam++;  if(bridge)n_degen_bridge++; }
        else      { n_sliver++; if(seam)n_sliver_seam++; if(bridge)n_sliver_bridge++; }
        if (nworst<WORST || ma<walt[nworst-1]){
            size_t pos = nworst<WORST ? nworst : WORST-1;
            if (nworst<WORST) nworst++;
            while (pos>0 && walt[pos-1]>ma){ walt[pos]=walt[pos-1]; widx[pos]=widx[pos-1]; pos--; }
            walt[pos]=ma; widx[pos]=f;
        }
    }

    printf("------------------------------------------------------------\n");
    printf("DEGENERACY CENSUS (area<%.0e -> degenerate; min-alt<%.2f vox -> sliver):\n", AREA_EPS, ALT_EPS);
    printf("  degenerate (~zero-area): %zu   (%zu near seam", n_degen, n_degen_seam);
    if (split>=0) printf(", %zu bridge", n_degen_bridge);
    printf(")\n");
    printf("  sliver (thin):           %zu   (%zu near seam", n_sliver, n_sliver_seam);
    if (split>=0) printf(", %zu bridge", n_sliver_bridge);
    printf(")\n");
    printf("  worst %zu by min-altitude:\n", nworst);
    for (size_t i=0;i<nworst;i++){
        size_t f=widx[i];
        double p0[3],p1[3],p2[3]; tri_pts(m,f,p0,p1,p2);
        int32_t *F=&m->f[f*3];
        printf("    f%zu minalt=%.5f area=%.5f  v[%d%c %d%c %d%c]%s%s\n",
               f, walt[i], tri_area(p0,p1,p2),
               F[0],region_of(F[0],split), F[1],region_of(F[1],split), F[2],region_of(F[2],split),
               face_near_seam(m,f,planes,np,band)?"  [seam]":"",
               face_cross_region(m,f,split)?"  [bridge]":"");
    }

    if (out_path){
        FILE *fp=fopen(out_path,"w");
        if (fp){
            uint8_t *bv=(uint8_t*)calloc(m->nv?m->nv:1,1);
            for (size_t f=0;f<m->nf;f++) if(flag[f]){ bv[m->f[f*3+0]]=1; bv[m->f[f*3+1]]=1; bv[m->f[f*3+2]]=1; }
            fprintf(fp,"# seam_audit --degen: %zu degenerate + %zu sliver triangles in RED, rest grey\n",
                    n_degen, n_sliver);
            for (size_t v=0;v<m->nv;v++)
                fprintf(fp, bv[v] ? "v %.6f %.6f %.6f 1.0 0.0 0.0\n" : "v %.6f %.6f %.6f 0.70 0.70 0.70\n",
                        (double)m->v[v*3+0],(double)m->v[v*3+1],(double)m->v[v*3+2]);
            for (size_t f=0;f<m->nf;f++)
                fprintf(fp,"f %d %d %d\n", m->f[f*3+0]+1,m->f[f*3+1]+1,m->f[f*3+2]+1);
            free(bv); fclose(fp);
            printf("wrote %s (full mesh; degenerate+sliver triangles red)\n", out_path);
        }
    }
    free(flag);
    return (int)(n_degen + n_sliver);
}

int main(int argc, char **argv)
{
    if (argc>=2 && !strcmp(argv[1],"--selftest")) return selftest();
    if (argc<2){
        fprintf(stderr,
          "Usage: %s <mesh.obj> [--cube <size=128>] [--plane <axis> <coord>]\n"
          "       [--band <vox=4>] [--gap <vox=1.0>] [--angle <deg=20>] [--hinge]\n"
          "       [--explain] [--degen] [--split <nvA>] [--out <bad.obj>]\n"
          "       %s --selftest\n"
          "  --degen   whole-mesh sliver/degenerate census (T-junction source)\n"
          "  --explain per-pair detail incl. min-altitude + T-junction verdict\n", argv[0], argv[0]);
        return 2;
    }
    const char *path=argv[1], *out_path=NULL;
    double cube=128.0, band=4.0, gap_max=1.0, angle_deg=20.0;
    int have_plane=0, plane_axis=0; double plane_coord=0.0;
    int min_excl=1;   /* exclude pairs sharing >= this many verts (1 = any shared) */
    int explain=0; long split=-1; long clones=-1; int degen=0;
    int near_set=0; double ncz=0, ncy=0, ncx=0, nr=0;  /* --near crop */
    for (int i=2;i<argc;i++){
        if (!strcmp(argv[i],"--cube") && i+1<argc) cube=atof(argv[++i]);
        else if (!strcmp(argv[i],"--band") && i+1<argc) band=atof(argv[++i]);
        else if (!strcmp(argv[i],"--gap") && i+1<argc) gap_max=atof(argv[++i]);
        else if (!strcmp(argv[i],"--angle") && i+1<argc) angle_deg=atof(argv[++i]);
        else if (!strcmp(argv[i],"--hinge")) min_excl=2;  /* also test 1-vertex-sharing folds */
        else if (!strcmp(argv[i],"--degen")) degen=1;     /* whole-mesh sliver/degenerate census */
        else if (!strcmp(argv[i],"--explain")) explain=1;
        else if (!strcmp(argv[i],"--split") && i+1<argc) split=atol(argv[++i]);
        else if (!strcmp(argv[i],"--clones") && i+1<argc) clones=atol(argv[++i]);
        else if (!strcmp(argv[i],"--near") && i+4<argc) {
            near_set=1; ncz=atof(argv[++i]); ncy=atof(argv[++i]);
            ncx=atof(argv[++i]); nr=atof(argv[++i]);
        }
        else if (!strcmp(argv[i],"--out") && i+1<argc) out_path=argv[++i];
        else if (!strcmp(argv[i],"--plane") && i+2<argc){ have_plane=1; plane_axis=atoi(argv[++i]); plane_coord=atof(argv[++i]); }
        else { fprintf(stderr,"seam_audit: unknown arg %s\n", argv[i]); return 2; }
    }
    double cos_par = cos(angle_deg * M_PI / 180.0);

    Mesh m;
    if (obj_read(path,&m)!=0){ fprintf(stderr,"seam_audit: cannot read %s\n", path); return 2; }
    printf("seam_audit: %s -- %zu verts, %zu faces (gap<=%.2f vox, angle<=%.0f deg)\n",
           path, m.nv, m.nf, gap_max, angle_deg);

    /* --near: crop to faces with any vertex within radius of a point, list them,
     * and write a small colored OBJ (green=comp A <split, blue=comp B, RED=clone
     * vert >= --clones). Lets you eyeball the exact before/after local geometry. */
    if (near_set) {
        uint8_t *keep = (uint8_t *)calloc(m.nf ? m.nf : 1, 1);
        size_t nk = 0, nk_clone = 0;
        for (size_t f = 0; f < m.nf; f++) {
            int hit = 0, uses_clone = 0;
            for (int k = 0; k < 3; k++) {
                int32_t vi = m.f[f*3+k];
                double dz=m.v[(size_t)vi*3+0]-ncz, dy=m.v[(size_t)vi*3+1]-ncy, dx=m.v[(size_t)vi*3+2]-ncx;
                if (sqrt(dz*dz+dy*dy+dx*dx) <= nr) hit = 1;
                if (clones >= 0 && (long)vi >= clones) uses_clone = 1;
            }
            if (hit) { keep[f]=1; nk++; if (uses_clone) nk_clone++; }
        }
        printf("--near (%.2f %.2f %.2f) r=%.2f: %zu faces (%zu use a clone vert >= %ld)\n",
               ncz, ncy, ncx, nr, nk, nk_clone, clones);
        size_t shown = 0;
        for (size_t f = 0; f < m.nf && shown < 120; f++) {
            if (!keep[f]) continue;
            int32_t a=m.f[f*3+0], b=m.f[f*3+1], c=m.f[f*3+2];
            printf("  f%zu  v[%d%c%s %d%c%s %d%c%s]\n", f,
                a, region_of(a,split), (clones>=0&&(long)a>=clones)?"*CLONE":"",
                b, region_of(b,split), (clones>=0&&(long)b>=clones)?"*CLONE":"",
                c, region_of(c,split), (clones>=0&&(long)c>=clones)?"*CLONE":"");
            shown++;
        }
        if (out_path) {
            FILE *fp = fopen(out_path, "w");
            if (fp) {
                for (size_t v = 0; v < m.nv; v++) {
                    float cr,cg,cb;
                    if (clones>=0 && (long)v>=clones)      { cr=1.0f; cg=0.0f; cb=0.0f; }
                    else if (split>=0 && (long)v<split)    { cr=0.0f; cg=0.8f; cb=0.0f; }
                    else                                   { cr=0.2f; cg=0.4f; cb=1.0f; }
                    fprintf(fp, "v %.6f %.6f %.6f %.3f %.3f %.3f\n",
                            (double)m.v[v*3+0],(double)m.v[v*3+1],(double)m.v[v*3+2],
                            (double)cr,(double)cg,(double)cb);
                }
                for (size_t f = 0; f < m.nf; f++)
                    if (keep[f]) fprintf(fp, "f %d %d %d\n",
                                         m.f[f*3+0]+1, m.f[f*3+1]+1, m.f[f*3+2]+1);
                fclose(fp);
                printf("wrote %s (%zu cropped faces; green=A blue=B RED=clone)\n", out_path, nk);
            }
        }
        free(keep); mesh_free(&m);
        return 0;
    }

    Plane planes[64]; size_t np=0;
    if (have_plane){ planes[0].axis=plane_axis; planes[0].coord=plane_coord; np=1; }
    else np = detect_planes(&m, cube, band, planes, 64);

    /* --degen: whole-mesh degeneracy census (independent of seam pairs; planes
     * only used to label which slivers sit in the weld zone). Diagnostic only. */
    if (degen){
        int n = degen_scan(&m, planes, np, band, split, out_path);
        printf("VERDICT: %s\n", n ? "degenerate/sliver triangles present (cleanable by edge-collapse)"
                                  : "no degenerate/sliver triangles");
        mesh_free(&m);
        return 0;
    }

    if (np==0){ printf("seam_audit: no seam plane found (cube=%.0f) -- nothing to audit\n", cube); mesh_free(&m); return 0; }

    static const char *AX="zyx";
    uint8_t *bad=(uint8_t*)calloc(m.nf?m.nf:1,1);
    Stats st = {0,0,0,1e30,0.0,0.0};
    size_t total_pairs=0;
    for (size_t p=0;p<np;p++){
        size_t pr = audit_plane(&m,&planes[p],band,cos_par,gap_max,min_excl,explain,split,bad,&st);
        total_pairs += pr;
        printf("  plane axis=%d (%c) coord=%.1f: %zu offending pair(s)\n",
               planes[p].axis, AX[planes[p].axis], planes[p].coord, pr);
    }

    /* Offender bbox (which side / where) + face count. */
    size_t nbad=0; double mn[3]={1e30,1e30,1e30}, mx[3]={-1e30,-1e30,-1e30};
    for (size_t f=0; f<m.nf; f++){
        if (!bad[f]) continue;
        nbad++;
        for (int k=0;k<3;k++){
            double a=m.v[(size_t)m.f[f*3+0]*3+(size_t)k];
            double b=m.v[(size_t)m.f[f*3+1]*3+(size_t)k];
            double c=m.v[(size_t)m.f[f*3+2]*3+(size_t)k];
            double lo=a<b?(a<c?a:c):(b<c?b:c), hi=a>b?(a>c?a:c):(b>c?b:c);
            if(lo<mn[k])mn[k]=lo; if(hi>mx[k])mx[k]=hi;
        }
    }
    printf("------------------------------------------------------------\n");
    printf("NEAR-COINCIDENT OVERLAP pairs (doubled surface / Z-fight): %zu\n", st.overlap);
    if (st.overlap) printf("  gap (vox): min %.3f  mean %.3f  max %.3f\n",
                           st.gmin, st.gsum/(double)st.overlap, st.gmax);
    printf("INTERPENETRATION pairs (non-parallel stab): %zu\n", st.stab);
    printf("FOLD-BACK self-intersections (edge-adjacent, interiors overlap): %zu\n", st.fold);
    printf("offending triangles: %zu", nbad);
    if (nbad) printf("  bbox z[%.1f,%.1f] y[%.1f,%.1f] x[%.1f,%.1f]",
                     mn[0],mx[0], mn[1],mx[1], mn[2],mx[2]);
    printf("\n");

    if (out_path && nbad){
        FILE *fp=fopen(out_path,"w");
        if (fp){
            /* Whole mesh IN CONTEXT: verts of a flagged triangle are RED, the
             * rest light grey; ALL faces are emitted so you can judge whether the
             * flagged triangles are genuine doubles or normal neighbours. */
            uint8_t *bad_vert = (uint8_t *)calloc(m.nv ? m.nv : 1, 1);
            for (size_t f=0; f<m.nf; f++)
                if (bad[f]) { bad_vert[m.f[f*3+0]]=1; bad_vert[m.f[f*3+1]]=1; bad_vert[m.f[f*3+2]]=1; }
            fprintf(fp,"# seam_audit: %zu flagged triangles in RED, rest grey\n", nbad);
            for (size_t v=0;v<m.nv;v++)
                fprintf(fp, bad_vert[v] ? "v %.6f %.6f %.6f 1.0 0.0 0.0\n"
                                        : "v %.6f %.6f %.6f 0.70 0.70 0.70\n",
                        (double)m.v[v*3+0],(double)m.v[v*3+1],(double)m.v[v*3+2]);
            for (size_t f=0;f<m.nf;f++)
                fprintf(fp,"f %d %d %d\n", m.f[f*3+0]+1,m.f[f*3+1]+1,m.f[f*3+2]+1);
            free(bad_vert);
            fclose(fp);
            printf("wrote %s (full mesh; %zu flagged triangles red)\n", out_path, nbad);
        }
    }
    printf("VERDICT: %s\n", total_pairs ? "OVERLAP/INVALID GEOMETRY near seam" : "clean near seam");

    free(bad); mesh_free(&m);
    return total_pairs ? 1 : 0;
}
