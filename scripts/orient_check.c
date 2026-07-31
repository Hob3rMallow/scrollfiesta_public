/* Decide whether a mesh is orientable: run a CORRECT per-face BFS
 * winding-propagation (recomputing each face's live half-edges at
 * visit time, so flips are never stale) and report same-direction
 * interior edges before and after. If after==0 the mesh is orientable
 * and an in-pipeline orientation pass will fix it; residual>0 means
 * genuine non-orientable triangulation (folds/overlaps) BPA must avoid.
 *
 * For each residual same-dir edge also classify the two incident
 * triangles by normal dot: dot<0 => fold/overlap (back-to-back layers),
 * dot>0 => pure winding contradiction (odd cycle).
 *
 * Plain C99, single file. Reads "v z y x" and "f a b c" (1-based). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

typedef struct { float z, y, x; } V3;

static V3 *V; static int *F;

static void tri_normal(int f, float n[3]) {
    int a=F[f*3+0], b=F[f*3+1], c=F[f*3+2];
    float e1[3]={V[b].z-V[a].z, V[b].y-V[a].y, V[b].x-V[a].x};
    float e2[3]={V[c].z-V[a].z, V[c].y-V[a].y, V[c].x-V[a].x};
    n[0]=e1[1]*e2[2]-e1[2]*e2[1];
    n[1]=e1[2]*e2[0]-e1[0]*e2[2];
    n[2]=e1[0]*e2[1]-e1[1]*e2[0];
    float L=sqrtf(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]); if(L>0){n[0]/=L;n[1]/=L;n[2]/=L;}
}

/* count interior edges whose two faces traverse the shared edge the SAME
 * way (same_dir). Returns same_dir; *out_interior set. */
static int count_same_dir(int nf, int *out_interior, int *out_bdry) {
    #define HSZ (1u<<23)
    static int *bk=NULL; static int cap=0;
    if(!bk){ bk=malloc(HSZ*sizeof(int)); }
    for(size_t i=0;i<HSZ;i++) bk[i]=-1;
    /* store directed half-edges: key by undirected (lo,hi); remember dir */
    typedef struct { int lo, hi, dir0, cnt, next; } HE; /* dir0: 1 if first face went lo->hi */
    if(cap<nf*3){ cap=nf*3; }
    HE *he=malloc((size_t)nf*3*sizeof(HE)); int en=0;
    int same=0, interior=0, bdry;
    for(int f=0; f<nf; f++){
        int v[3]={F[f*3],F[f*3+1],F[f*3+2]};
        for(int e=0;e<3;e++){
            int a=v[e], b=v[(e+1)%3]; int lo=a,hi=b,dir= (a<b)?1:0; if(a>b){lo=b;hi=a;}
            uint64_t h=(uint64_t)(uint32_t)lo*0x9E3779B97F4A7C15ULL+(uint64_t)(uint32_t)hi; h^=h>>33; h*=0xff51afd7ed558ccdULL; h^=h>>33;
            size_t bb=h&(HSZ-1); int idx=bk[bb];
            while(idx>=0){ if(he[idx].lo==lo&&he[idx].hi==hi) break; idx=he[idx].next; }
            if(idx<0){ idx=en++; he[idx].lo=lo; he[idx].hi=hi; he[idx].dir0=dir; he[idx].cnt=1; he[idx].next=bk[bb]; bk[bb]=idx; }
            else { he[idx].cnt++; if(he[idx].cnt==2){ interior++; if(he[idx].dir0==dir) same++; } }
        }
    }
    bdry=0; for(int i=0;i<en;i++) if(he[i].cnt==1) bdry++;
    free(he);
    if(out_interior)*out_interior=interior; if(out_bdry)*out_bdry=bdry;
    return same;
}

int main(int argc, char **argv){
    if(argc<2){ fprintf(stderr,"usage: %s mesh.obj\n",argv[0]); return 1; }
    FILE *f=fopen(argv[1],"r"); if(!f){perror(argv[1]);return 1;}
    int vcap=1<<16,nv=0,fcap=1<<16,nf=0;
    V=malloc((size_t)vcap*sizeof(V3)); F=malloc((size_t)fcap*3*sizeof(int));
    char line[256];
    while(fgets(line,sizeof line,f)){
        if(line[0]=='v'&&line[1]==' '){ V3 v; if(sscanf(line,"v %f %f %f",&v.z,&v.y,&v.x)==3){ if(nv>=vcap){vcap*=2;V=realloc(V,(size_t)vcap*sizeof(V3));} V[nv++]=v; } }
        else if(line[0]=='f'&&line[1]==' '){ int a,b,c; if(sscanf(line,"f %d %d %d",&a,&b,&c)==3){ if(nf>=fcap){fcap*=2;F=realloc(F,(size_t)fcap*3*sizeof(int));} F[nf*3+0]=a-1;F[nf*3+1]=b-1;F[nf*3+2]=c-1;nf++; } }
    }
    fclose(f);
    fprintf(stderr,"loaded %d v %d f\n",nv,nf);

    int interior=0,bdry=0;
    int same_before=count_same_dir(nf,&interior,&bdry);
    printf("interior edges: %d   boundary edges: %d\n",interior,bdry);
    printf("same_dir BEFORE: %d  (%.2f%% of interior)\n",same_before,100.0*same_before/(interior?interior:1));

    /* Build manifold-edge face adjacency (run==2 edges). */
    #define HSZ2 (1u<<23)
    int *bk=malloc(HSZ2*sizeof(int)); for(size_t i=0;i<HSZ2;i++)bk[i]=-1;
    typedef struct { int lo,hi,f0,f1,next; } EF;
    EF *ef=malloc((size_t)nf*3*sizeof(EF)); int en=0;
    for(int fc=0;fc<nf;fc++){ int v[3]={F[fc*3],F[fc*3+1],F[fc*3+2]};
        for(int e=0;e<3;e++){ int a=v[e],b=v[(e+1)%3],lo=a<b?a:b,hi=a<b?b:a;
            uint64_t h=(uint64_t)(uint32_t)lo*0x9E3779B97F4A7C15ULL+(uint64_t)(uint32_t)hi; h^=h>>33;h*=0xff51afd7ed558ccdULL;h^=h>>33;
            size_t bb=h&(HSZ2-1); int idx=bk[bb];
            while(idx>=0){ if(ef[idx].lo==lo&&ef[idx].hi==hi)break; idx=ef[idx].next; }
            if(idx<0){ idx=en++; ef[idx].lo=lo;ef[idx].hi=hi;ef[idx].f0=fc;ef[idx].f1=-1;ef[idx].next=bk[bb];bk[bb]=idx; }
            else { if(ef[idx].f1<0) ef[idx].f1=fc; else ef[idx].f1=-2; /* >2 faces: non-manifold edge, mark */ }
        }
    }
    /* face adjacency CSR over manifold (exactly-2) edges */
    int *deg=calloc((size_t)nf,sizeof(int));
    for(int i=0;i<en;i++) if(ef[i].f1>=0){ deg[ef[i].f0]++; deg[ef[i].f1]++; }
    int *off=malloc((size_t)(nf+1)*sizeof(int)); off[0]=0; for(int i=0;i<nf;i++)off[i+1]=off[i]+deg[i];
    int *adjf=malloc((size_t)off[nf]*sizeof(int)); int *cur=malloc((size_t)nf*sizeof(int)); memcpy(cur,off,(size_t)nf*sizeof(int));
    for(int i=0;i<en;i++) if(ef[i].f1>=0){ adjf[cur[ef[i].f0]++]=ef[i].f1; adjf[cur[ef[i].f1]++]=ef[i].f0; }

    /* CORRECT BFS orientation: at visit time read LIVE windings. */
    uint8_t *vis=calloc((size_t)nf,1); int *q=malloc((size_t)nf*sizeof(int));
    long flips=0; int comps=0;
    for(int s=0;s<nf;s++){ if(vis[s])continue; comps++; int qh=0,qt=0; q[qt++]=s; vis[s]=1;
        while(qh<qt){ int fc=q[qh++];
            int a=F[fc*3+0],b=F[fc*3+1],c=F[fc*3+2];
            for(int k=off[fc];k<off[fc+1];k++){ int fn=adjf[k]; if(vis[fn])continue;
                /* find shared edge between fc and fn; check if fn traverses it same dir as fc */
                int na=F[fn*3+0],nb=F[fn*3+1],nc=F[fn*3+2];
                int pe[3][2]={{a,b},{b,c},{c,a}}; int ne[3][2]={{na,nb},{nb,nc},{nc,na}};
                int needflip=0;
                for(int pi=0;pi<3;pi++)for(int ni=0;ni<3;ni++){
                    if((pe[pi][0]==ne[ni][0]&&pe[pi][1]==ne[ni][1])){ needflip=1; }
                }
                if(needflip){ int t=F[fn*3+1]; F[fn*3+1]=F[fn*3+2]; F[fn*3+2]=t; flips++; }
                vis[fn]=1; q[qt++]=fn;
            }
        }
    }
    int interior2=0,same_after=count_same_dir(nf,&interior2,NULL);
    printf("BFS orient: %d components, %ld flips\n",comps,flips);
    printf("same_dir AFTER:  %d  (%.2f%% of interior)\n",same_after,100.0*same_after/(interior2?interior2:1));

    /* Classify residual same_dir edges by normal dot. */
    if(same_after>0){
        for(size_t i=0;i<HSZ2;i++)bk[i]=-1; en=0;
        int fold=0,wind=0;
        for(int fc=0;fc<nf;fc++){ int v[3]={F[fc*3],F[fc*3+1],F[fc*3+2]};
            for(int e=0;e<3;e++){ int a=v[e],b=v[(e+1)%3],lo=a<b?a:b,hi=a<b?b:a,dir=a<b?1:0;
                uint64_t h=(uint64_t)(uint32_t)lo*0x9E3779B97F4A7C15ULL+(uint64_t)(uint32_t)hi;h^=h>>33;h*=0xff51afd7ed558ccdULL;h^=h>>33;
                size_t bb=h&(HSZ2-1); int idx=bk[bb];
                while(idx>=0){ if(ef[idx].lo==lo&&ef[idx].hi==hi)break; idx=ef[idx].next; }
                if(idx<0){ idx=en++; ef[idx].lo=lo;ef[idx].hi=hi;ef[idx].f0=fc;ef[idx].f1=dir;/*reuse f1=dir0*/ ef[idx].next=bk[bb];bk[bb]=idx; }
                else { if(ef[idx].f1==dir){ float n0[3],n1[3]; tri_normal(ef[idx].f0,n0); tri_normal(fc,n1); float d=n0[0]*n1[0]+n0[1]*n1[1]+n0[2]*n1[2]; if(d<0)fold++; else wind++; } }
            }
        }
        printf("residual classification: fold(dot<0)=%d  winding(dot>=0)=%d\n",fold,wind);
    }
    return 0;
}
