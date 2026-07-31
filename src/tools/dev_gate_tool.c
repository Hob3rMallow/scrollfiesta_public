/* dev_gate_tool.c -- run the end-of-pipeline developability "bad sheet" gate on
 * one component OBJ (the same DevGate_classify the pipeline calls). Loads the
 * mesh, derives its PCA normal (the oracle's ray direction), and prints a
 * one-line machine-parseable verdict so a script can sweep the whole grid's
 * step12_final components without a multi-hour re-mesh.
 *
 * Usage:
 *   dev_gate_tool <component.obj> [--kthresh R] [--frac F] [--minbad N]
 *                                 [--grid G] [--no-oracle]
 *   dev_gate_tool --selftest
 */
#include "../common/ves_platform.h"
#include "../common/arena.h"
#include "../common/obj_io.h"
#include "../common/mesh_types.h"
#include "../common/pca.h"
#include "../remesh/dev_gate.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *reason_str(int r){
    if((r&DEVGATE_NONDEV) && (r&DEVGATE_ORACLE)) return "nondev+oracle";
    if(r&DEVGATE_NONDEV) return "nondev";
    if(r&DEVGATE_ORACLE) return "oracle";
    return "none";
}

static int run_selftest(void){
    int fail=0;
    Arena_T arena=Arena_new();
    /* 7x7 flat grid: developable, gate must pass */
    int N=7; size_t nv=(size_t)N*N;
    float *V=(float*)malloc(nv*3*sizeof(float));
    for(int j=0;j<N;j++) for(int i=0;i<N;i++){ size_t k=(size_t)j*N+i; V[k*3]=0.0f; V[k*3+1]=(float)j; V[k*3+2]=(float)i; }
    size_t nf=(size_t)(N-1)*(N-1)*2; int32_t *F=(int32_t*)malloc(nf*3*sizeof(int32_t));
    size_t k=0;
    for(int j=0;j<N-1;j++) for(int i=0;i<N-1;i++){
        int32_t a=(int32_t)(j*N+i),b=(int32_t)(j*N+i+1),c=(int32_t)((j+1)*N+i),d=(int32_t)((j+1)*N+i+1);
        F[k*3]=a;F[k*3+1]=b;F[k*3+2]=d;k++; F[k*3]=a;F[k*3+1]=d;F[k*3+2]=c;k++;
    }
    ComponentMesh cm; memset(&cm,0,sizeof cm);
    cm.verts=V; cm.faces=F; cm.nv=nv; cm.nf=nf; cm.self=&cm;
    PCA_normal(V,nv,cm.pca_normal,cm.centroid);
    DevGateParams p; DevGate_default_params(&p);
    p.min_bad=1; p.frac_thresh=0.0; p.oracle_grid=0;  /* criterion-1 only, sensitive */
    DevGateVerdict v;
    DevGate_classify(arena,&cm,&p,&v);
    int f1=(v.bad!=0); fail|=f1;
    fprintf(stderr,"[selftest] flat grid bad=%d (int_bad=%zu) -> %s\n", v.bad, v.n_interior_bad, f1?"FAIL":"ok");
    /* raise a cluster of interior verts -> non-developable, gate must flag */
    for(int j=2;j<=4;j++) for(int i=2;i<=4;i++){ size_t idx=(size_t)j*N+i; V[idx*3]= (i+j)%2 ? 0.9f : -0.9f; }
    DevGate_classify(arena,&cm,&p,&v);
    int f2=!(v.bad && (v.reason&DEVGATE_NONDEV)); fail|=f2;
    fprintf(stderr,"[selftest] crumpled grid bad=%d reason=%s int_bad=%zu maxK=%.2f -> %s\n",
            v.bad, reason_str(v.reason), v.n_interior_bad, v.max_abs_k, f2?"FAIL":"ok");
    free(V); free(F); Arena_dispose(&arena);
    fprintf(stderr,"=== dev_gate_tool selftest %s ===\n", fail?"FAILED":"PASSED");
    return fail?1:0;
}

int main(int argc, char **argv){
    if(argc>=2 && !strcmp(argv[1],"--selftest")) return run_selftest();
    if(argc<2){
        fprintf(stderr,"Usage: %s <component.obj> [--kthresh R] [--frac F] [--minbad N] [--grid G] [--no-oracle]\n"
                       "       %s --selftest\n", argv[0], argv[0]);
        return 1;
    }
    const char *in=argv[1];
    DevGateParams p; DevGate_default_params(&p);
    for(int i=2;i<argc;i++){
        if(!strcmp(argv[i],"--kthresh")&&i+1<argc) p.k_thresh=atof(argv[++i]);
        else if(!strcmp(argv[i],"--frac")&&i+1<argc) p.frac_thresh=atof(argv[++i]);
        else if(!strcmp(argv[i],"--minbad")&&i+1<argc) p.min_bad=(size_t)atol(argv[++i]);
        else if(!strcmp(argv[i],"--grid")&&i+1<argc) p.oracle_grid=atoi(argv[++i]);
        else if(!strcmp(argv[i],"--no-oracle")) p.oracle_grid=0;
        else { fprintf(stderr,"unknown arg %s\n",argv[i]); return 1; }
    }
    Arena_T arena=Arena_new();
    float *V=NULL; int32_t *F=NULL; size_t nv=0,nf=0;
    if(ObjIO_read(arena,in,&V,&nv,&F,&nf)!=0){
        fprintf(stderr,"dev_gate_tool: cannot read %s\n",in); Arena_dispose(&arena); return 1;
    }
    ComponentMesh cm; memset(&cm,0,sizeof cm);
    cm.verts=V; cm.faces=F; cm.nv=nv; cm.nf=nf; cm.self=&cm;
    if(nv) PCA_normal(V,nv,cm.pca_normal,cm.centroid);
    DevGateVerdict v;
    DevGate_classify(arena,&cm,&p,&v);
    /* one machine-parseable line */
    printf("GATE bad=%d reason=%-13s nv=%-7zu nf=%-7zu int=%-7zu int_bad=%-5zu frac=%.4f maxK=%.3f oracle=%d | %s\n",
           v.bad, reason_str(v.reason), nv, nf, v.n_interior, v.n_interior_bad,
           v.frac_interior_bad, v.max_abs_k, v.oracle_sheets, in);
    Arena_dispose(&arena);
    return 0;
}
