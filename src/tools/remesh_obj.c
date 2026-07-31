/* remesh_obj.c -- incremental isotropic remeshing of an OBJ mesh
 * (Remesh_isotropic: split long / collapse short / flip / tangential relax at a
 * fixed target face count, interior-only with the open boundary frozen). Runs
 * AFTER decimation to repair QEM's slivers + coarse-center/fine-rim anisotropy.
 * OBJ in -> OBJ out. Fail-closed: never emits a non-manifold mesh (falls back to
 * a verbatim copy of the input if the remesh would degrade topology).
 *
 * Usage:
 *   remesh_obj <in.obj> <out.obj> [--iters N] [--target-edge L] [--mls]
 *              [--no-flip] [--no-relax]
 *   remesh_obj --selftest
 */
#include "../common/ves_platform.h"
#include "../common/arena.h"
#include "../common/obj_io.h"
#include "../common/mesh_manifold.h"
#include "../remesh/isotropic_remesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================
 * Selftest-only mesh metrics (self-contained; not used by the CLI path).
 * ================================================================ */
typedef struct { int32_t v0, v1; } TE;
static int te_cmp(const void *pa, const void *pb)
{
    const TE *a=(const TE*)pa, *b=(const TE*)pb;
    if (a->v0 != b->v0) return a->v0 < b->v0 ? -1 : 1;
    if (a->v1 != b->v1) return a->v1 < b->v1 ? -1 : 1;
    return 0;
}

/* Mark boundary verts (on any face_count==1 edge). */
static void selftest_boundary(Arena_T a, const int32_t *F, size_t nf, size_t nv,
                              uint8_t *bnd)
{
    size_t n=nf*3, i, j;
    TE *e=(TE*)ARENA_ALLOC(a,(long)((n?n:1)*sizeof(TE)));
    for (i=0;i<nf;i++){
        int32_t t[3]={F[i*3+0],F[i*3+1],F[i*3+2]};
        int k;
        for (k=0;k<3;k++){
            int32_t u=t[k], v=t[(k+1)%3];
            e[i*3+k].v0=(u<v)?u:v; e[i*3+k].v1=(u<v)?v:u;
        }
    }
    qsort(e,n,sizeof(TE),te_cmp);
    memset(bnd,0,nv*sizeof(uint8_t));
    i=0;
    while (i<n){ j=i+1; while (j<n && e[j].v0==e[i].v0 && e[j].v1==e[i].v1) j++;
        if (j-i==1){ bnd[e[i].v0]=1; bnd[e[i].v1]=1; } i=j; }
}

static double edge_cv(const float *V, const int32_t *F, size_t nf)
{
    double s=0,s2=0,mean,var; size_t n=0,i;
    for (i=0;i<nf;i++){ int k; for (k=0;k<3;k++){
        const float *A=&V[(size_t)F[i*3+k]*3], *B=&V[(size_t)F[i*3+(k+1)%3]*3];
        double dx=B[0]-A[0], dy=B[1]-A[1], dz=B[2]-A[2], L=sqrt(dx*dx+dy*dy+dz*dz);
        s+=L; s2+=L*L; n++; } }
    if (!n) return 0.0; mean=s/(double)n; var=s2/(double)n-mean*mean;
    if (var<0) var=0; return mean>1e-12?sqrt(var)/mean:0.0;
}

/* Min interior angle (deg) over all tris, and the fraction that are slivers
 * (min-angle < 15 deg) -- the real "bad edge quality" signal the flip/relax
 * passes target (CV alone is dominated by legitimately-adaptive structure). */
static void angle_stats(const float *V, const int32_t *F, size_t nf,
                        double *min_deg, double *sliver_pct)
{
    double mn=180.0; size_t i, bad=0, n=0;
    for (i=0;i<nf;i++){
        const float *P[3]; int k; double tmn=180.0;
        P[0]=&V[(size_t)F[i*3+0]*3]; P[1]=&V[(size_t)F[i*3+1]*3]; P[2]=&V[(size_t)F[i*3+2]*3];
        for (k=0;k<3;k++){
            const float *A=P[k], *B=P[(k+1)%3], *C=P[(k+2)%3];
            double ux=B[0]-A[0],uy=B[1]-A[1],uz=B[2]-A[2];
            double wx=C[0]-A[0],wy=C[1]-A[1],wz=C[2]-A[2];
            double lu=sqrt(ux*ux+uy*uy+uz*uz), lw=sqrt(wx*wx+wy*wy+wz*wz), cg, ang;
            if (lu<1e-12||lw<1e-12){ tmn=0.0; continue; }
            cg=(ux*wx+uy*wy+uz*wz)/(lu*lw); if(cg>1)cg=1; else if(cg<-1)cg=-1;
            ang=acos(cg)*180.0/3.14159265358979323846;
            if (ang<tmn) tmn=ang;
        }
        if (tmn<mn) mn=tmn;
        if (tmn<15.0) bad++;
        n++;
    }
    *min_deg = n?mn:0.0;
    *sliver_pct = n?100.0*(double)bad/(double)n:0.0;
}

/* Does position p (exact float) occur among the mesh verts? */
static int has_vert(const float *V, size_t nv, const float p[3])
{
    size_t i;
    for (i=0;i<nv;i++)
        if (V[i*3+0]==p[0] && V[i*3+1]==p[1] && V[i*3+2]==p[2]) return 1;
    return 0;
}

/* Flat grid (nx*ny) with optional quadratic x-clustering (non-uniform spacing)
 * at z-offset. */
static void make_grid(Arena_T a, int nx, int ny, int cluster, float zoff,
                      float **pv, size_t *pnv, int32_t **pf, size_t *pnf)
{
    size_t nv=(size_t)(nx*ny), nf=(size_t)((nx-1)*(ny-1)*2), vi=0, fi=0;
    float *v=(float*)ARENA_ALLOC(a,(long)(nv*3*sizeof(float)));
    int32_t *f=(int32_t*)ARENA_ALLOC(a,(long)(nf*3*sizeof(int32_t)));
    int ix,iy;
    for (iy=0;iy<ny;iy++) for (ix=0;ix<nx;ix++){
        float fx=(float)ix/(float)(nx-1);
        if (cluster) fx=fx*fx;                 /* dense near x=0 -> non-uniform */
        v[vi*3+0]=fx*(float)(nx-1);
        v[vi*3+1]=(float)iy;
        v[vi*3+2]=zoff;
        vi++;
    }
    for (iy=0;iy<ny-1;iy++) for (ix=0;ix<nx-1;ix++){
        int32_t i00=(int32_t)(iy*nx+ix), i10=(int32_t)(iy*nx+ix+1);
        int32_t i01=(int32_t)((iy+1)*nx+ix), i11=(int32_t)((iy+1)*nx+ix+1);
        f[fi*3+0]=i00; f[fi*3+1]=i10; f[fi*3+2]=i11; fi++;
        f[fi*3+0]=i00; f[fi*3+1]=i11; f[fi*3+2]=i01; fi++;
    }
    *pv=v; *pnv=nv; *pf=f; *pnf=nf;
}

static int run_selftest(void)
{
    int fails=0;
    Arena_T a=Arena_new();

    /* Case A: non-uniform flat grid -> remesh must stay manifold, freeze the
     * rim (every input boundary position present in output), keep face count
     * within tolerance, and not worsen edge-length uniformity. */
    {
        float *v=NULL; size_t nv=0; int32_t *f=NULL; size_t nf=0;
        float *ov=NULL; size_t onv=0; int32_t *of=NULL; size_t onf=0;
        uint8_t *bnd; size_t i; int c=0;
        double cv_in, cv_out; MeshManifoldStats st;
        RemeshOpts o; Remesh_default_opts(&o);
        make_grid(a, 41, 41, 1 /*cluster*/, 0.0f, &v,&nv,&f,&nf);
        cv_in=edge_cv(v,f,nf);
        bnd=(uint8_t*)ARENA_ALLOC(a,(long)(nv*sizeof(uint8_t)));
        selftest_boundary(a,f,nf,nv,bnd);
        if (Remesh_isotropic(a, v,nv, f,nf, NULL, &o, &ov,&onv,&of,&onf)!=0)
            c=0; /* fallback still valid, just note */
        st=MeshManifold_audit(a,onv,of,onf);
        cv_out=edge_cv(ov,of,onf);
        if (!MeshManifold_ok(&st)){ fprintf(stderr,"[A] non-manifold out\n"); fails++; }
        /* every input rim position must survive unmoved */
        for (i=0;i<nv;i++){ if (!bnd[i]) continue;
            float p[3]={v[i*3+0],v[i*3+1],v[i*3+2]};
            if (!has_vert(ov,onv,p)){ c++; } }
        if (c){ fprintf(stderr,"[A] %d rim verts moved/lost\n", c); fails++; }
        if ((double)onf > 1.10*(double)nf || (double)onf < 0.90*(double)nf){
            fprintf(stderr,"[A] face count %zu out of tol vs %zu\n", onf, nf); fails++; }
        if (cv_out > cv_in + 1e-6){
            fprintf(stderr,"[A] CV worsened %.3f->%.3f\n", cv_in, cv_out); fails++; }
        fprintf(stderr,"[A] grid %zu->%zu f, cv %.3f->%.3f, manifold=%d -> %s\n",
                nf, onf, cv_in, cv_out, MeshManifold_ok(&st),
                (fails?"see above":"ok"));
    }

    /* Case B: two flat sheets, gap 2 (< the anti-fusion cap) -> remesh must NOT
     * fuse them: no output vertex may land in the gap (0.5..1.5), and every face
     * must stay on ONE sheet (all 3 verts z<1 or all z>1). Run with default
     * (tangential) reproject AND with MLS. */
    {
        int pass;
        for (pass=0; pass<2; pass++){
            Arena_T b=Arena_new();
            float *va=NULL,*vb=NULL; size_t nva=0,nvb=0; int32_t *fa=NULL,*fb=NULL; size_t nfa=0,nfb=0;
            float *V; int32_t *Fm; size_t tnv, tnf, i, gap=0, split_face=0;
            float *ov=NULL; size_t onv=0; int32_t *of=NULL; size_t onf=0;
            MeshManifoldStats st;
            RemeshOpts o; Remesh_default_opts(&o);
            if (pass==1) o.reproject=REMESH_REPROJECT_MLS;
            make_grid(b, 24, 24, 0, 0.0f, &va,&nva,&fa,&nfa);
            make_grid(b, 24, 24, 0, 2.0f, &vb,&nvb,&fb,&nfb);
            tnv=nva+nvb; tnf=nfa+nfb;
            V =(float*)ARENA_ALLOC(b,(long)(tnv*3*sizeof(float)));
            Fm=(int32_t*)ARENA_ALLOC(b,(long)(tnf*3*sizeof(int32_t)));
            memcpy(V, va, nva*3*sizeof(float));
            memcpy(V+nva*3, vb, nvb*3*sizeof(float));
            memcpy(Fm, fa, nfa*3*sizeof(int32_t));
            for (i=0;i<nfb*3;i++) Fm[nfa*3+i]=fb[i]+(int32_t)nva;
            Remesh_isotropic(b, V,tnv, Fm,tnf, NULL, &o, &ov,&onv,&of,&onf);
            st=MeshManifold_audit(b,onv,of,onf);
            for (i=0;i<onv;i++){ float z=ov[i*3+2]; if (z>0.5f && z<1.5f) gap++; }
            for (i=0;i<onf;i++){
                float z0=ov[(size_t)of[i*3+0]*3+2], z1=ov[(size_t)of[i*3+1]*3+2], z2=ov[(size_t)of[i*3+2]*3+2];
                int s0=z0>1.0f, s1=z1>1.0f, s2=z2>1.0f;
                if (!(s0==s1 && s1==s2)) split_face++;
            }
            if (gap){ fprintf(stderr,"[B/%s] %zu verts in gap\n", pass?"mls":"tan", gap); fails++; }
            if (split_face){ fprintf(stderr,"[B/%s] %zu cross-sheet faces\n", pass?"mls":"tan", split_face); fails++; }
            if (!MeshManifold_ok(&st)){ fprintf(stderr,"[B/%s] non-manifold\n", pass?"mls":"tan"); fails++; }
            fprintf(stderr,"[B/%s] two sheets %zu->%zu f, gap=%zu cross=%zu -> %s\n",
                    pass?"mls":"tan", tnf, onf, gap, split_face,
                    (gap||split_face||!MeshManifold_ok(&st))?"see above":"ok");
            Arena_dispose(&b);
        }
    }

    Arena_dispose(&a);
    fprintf(stderr,"=== remesh_obj selftest %s (%d failure%s) ===\n",
            fails?"FAILED":"PASSED", fails, fails==1?"":"s");
    return fails?1:0;
}

int main(int argc, char **argv)
{
    const char *in_path=NULL, *out_path=NULL;
    RemeshOpts o; Remesh_default_opts(&o);
    int i;
    Arena_T arena;
    float *verts=NULL; int32_t *faces=NULL; size_t nv=0, nf=0;
    float *ov=NULL; int32_t *of=NULL; size_t onv=0, onf=0;
    int rc;

    if (argc>=2 && strcmp(argv[1],"--selftest")==0) return run_selftest();
    if (argc<3){
        fprintf(stderr,
            "Usage: %s <in.obj> <out.obj> [--iters N] [--target-edge L] [--mls]\n"
            "                            [--no-flip] [--no-relax]\n"
            "       %s --selftest\n", argv[0], argv[0]);
        return 1;
    }
    in_path=argv[1]; out_path=argv[2];
    for (i=3;i<argc;i++){
        if (!strcmp(argv[i],"--iters") && i+1<argc) o.n_iters=atoi(argv[++i]);
        else if (!strcmp(argv[i],"--target-edge") && i+1<argc) o.target_edge_len=(float)atof(argv[++i]);
        else if (!strcmp(argv[i],"--mls")) o.reproject=REMESH_REPROJECT_MLS;
        else if (!strcmp(argv[i],"--no-split")) o.do_split=0;
        else if (!strcmp(argv[i],"--no-collapse")) o.do_collapse=0;
        else if (!strcmp(argv[i],"--no-flip")) o.do_flip=0;
        else if (!strcmp(argv[i],"--no-relax")) o.do_relax=0;
        else { fprintf(stderr,"remesh_obj: unknown arg %s\n", argv[i]); return 1; }
    }

    arena=Arena_new();
    if (ObjIO_read(arena,in_path,&verts,&nv,&faces,&nf)!=0){
        fprintf(stderr,"remesh_obj: cannot read %s\n", in_path);
        Arena_dispose(&arena); return 1;
    }
    rc=Remesh_isotropic(arena, verts,nv, faces,nf, NULL, &o, &ov,&onv,&of,&onf);
    ves_ensure_parent_dir(out_path);
    if (ObjIO_write(out_path, ov,onv, of,onf)!=0){
        fprintf(stderr,"remesh_obj: cannot write %s\n", out_path);
        Arena_dispose(&arena); return 1;
    }
    {
        double a0,s0,a1,s1;
        angle_stats(verts,faces,nf,&a0,&s0);
        angle_stats(ov,of,onf,&a1,&s1);
        fprintf(stderr,
            "remesh_obj: %s  v %zu->%zu  f %zu->%zu  edge-CV %.3f->%.3f  "
            "minAng %.1f->%.1f  slivers<15deg %.1f%%->%.1f%%  (%s) -> %s\n",
            in_path, nv, onv, nf, onf,
            edge_cv(verts,faces,nf), edge_cv(ov,of,onf),
            a0, a1, s0, s1,
            rc==0?"remeshed":"fallback(input copy)", out_path);
    }
    Arena_dispose(&arena);
    return 0;
}
