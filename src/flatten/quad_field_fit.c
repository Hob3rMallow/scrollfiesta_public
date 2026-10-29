#include "quad_field.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* One compact root directory, not a UV rectangle. Root identities depend
 * only on the fixed global chart and absolute dyadic lattice. */
typedef struct {
    QuadFieldCell *cell;
    int32_t *hash;
    size_t count,capacity,slots;
} QffRoots;

static uint64_t qff_hash(int32_t chart,int32_t x,int32_t y)
{
    uint64_t h=(uint32_t)chart;
    h=(h^(uint32_t)x)*UINT64_C(1099511628211);
    h=(h^(uint32_t)y)*UINT64_C(1099511628211);
    h^=h>>30; h*=UINT64_C(0xbf58476d1ce4e5b9); h^=h>>27;
    return h;
}
static int qff_insert(QffRoots *r,QuadFieldCell cell)
{
    size_t slot=0;
    if (r->slots==0 || 2*r->count>=r->slots) {
        size_t wanted=r->slots ? 2*r->slots : 16384;
        int32_t *hash=(int32_t *)calloc(wanted,sizeof(int32_t));
        if (hash==NULL) return -1;
        for (size_t i=0;i<r->count;i++) {
            const QuadFieldCell *c=r->cell+i;
            size_t at=(size_t)qff_hash(c->chart,c->x,c->y)&(wanted-1);
            while (hash[at]) at=(at+1)&(wanted-1);
            hash[at]=(int32_t)i+1;
        }
        free(r->hash); r->hash=hash; r->slots=wanted;
    }
    slot=(size_t)qff_hash(cell.chart,cell.x,cell.y)&(r->slots-1);
    while (r->hash[slot]) {
        const QuadFieldCell *c=r->cell+r->hash[slot]-1;
        if (c->chart==cell.chart && c->x==cell.x && c->y==cell.y) return 0;
        slot=(slot+1)&(r->slots-1);
    }
    if (r->count>=QUAD_FIELD_MAX_CELLS) return -1;
    if (r->count==r->capacity) {
        size_t wanted=r->capacity ? 2*r->capacity : 8192;
        QuadFieldCell *grown=(QuadFieldCell *)realloc(r->cell,wanted*sizeof *grown);
        if (grown==NULL) return -1;
        r->cell=grown; r->capacity=wanted;
    }
    r->cell[r->count]=cell; r->hash[slot]=(int32_t)++r->count;
    return 0;
}
static int qff_compare(const void *pa,const void *pb)
{
    const QuadFieldCell *a=(const QuadFieldCell *)pa,*b=(const QuadFieldCell *)pb;
    if (a->chart!=b->chart) return a->chart<b->chart ? -1 : 1;
    if (a->y!=b->y) return a->y<b->y ? -1 : 1;
    return a->x<b->x ? -1 : a->x>b->x;
}

/* Surface area quadrature gives dense and sparse input triangulations the
 * same units of evidence. Isolated points are retained as unit observations;
 * their lack of supported faces remains a topology qualification issue. */
static int qff_weights(Arena_T arena,const MeshBinData *m,double **out)
{
    if ((m->nv && m->verts==NULL) || (m->nf && m->faces==NULL)) return -1;
    double *w=(double *)ARENA_CALLOC(arena,m->nv,sizeof(double));
    for (size_t i=0;i<m->nf;i++) {
        double u[3]={0},v[3]={0},n[3]={0},area=0;
        int32_t a=m->faces[3*i],b=m->faces[3*i+1],c=m->faces[3*i+2];
        if (a<0 || b<0 || c<0 || (size_t)a>=m->nv || (size_t)b>=m->nv || (size_t)c>=m->nv) return -1;
        for (int k=0;k<3;k++) { u[k]=(double)m->verts[3*(size_t)b+k]-m->verts[3*(size_t)a+k]; v[k]=(double)m->verts[3*(size_t)c+k]-m->verts[3*(size_t)a+k]; }
        n[0]=u[1]*v[2]-u[2]*v[1]; n[1]=u[2]*v[0]-u[0]*v[2]; n[2]=u[0]*v[1]-u[1]*v[0];
        area=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])/6;
        if (!isfinite(area)) return -1;
        w[a]+=area; w[b]+=area; w[c]+=area;
    }
    for (size_t i=0;i<m->nv;i++) if (w[i]==0) w[i]=1;
    *out=w; return 0;
}

int QuadField_fit_stream(Arena_T arena,size_t chunks,QuadFieldReader reader,
                          void *context,int root_cells,double step,double tolerance,
                          QuadFieldCheckpoint checkpoint,void *checkpoint_context,
                          QuadField_T *out)
{
    QffRoots roots={0};
    QuadField_T field=NULL,coarse=NULL;
    double started=ves_clock_sec();
    int rc=-1;
    if (!arena || !reader || !out || root_cells<1 || root_cells>(1<<28) ||
        (root_cells&(root_cells-1)) || !isfinite(step) || step<=0 ||
        !isfinite(tolerance) || tolerance<=0) return -1;
    *out=NULL;
    for (size_t chunk=0;chunk<chunks;chunk++) {
        Arena_T scratch=Arena_new();
        MeshBinData mesh={0};
        int32_t *chart=NULL;
        int ok=reader(context,scratch,chunk,&mesh,&chart)==0;
        if (ok && mesh.nv && (!mesh.verts || !mesh.uv || !chart)) ok=0;
        for (size_t i=0;ok && i<mesh.nv;i++) if (chart[i]>=0) {
            double x=floor(mesh.uv[2*i]/(step*root_cells))*root_cells;
            double y=floor(mesh.uv[2*i+1]/(step*root_cells))*root_cells;
            if (!isfinite(x) || !isfinite(y) || x<INT32_MIN || y<INT32_MIN ||
                x>INT32_MAX-root_cells || y>INT32_MAX-root_cells) { ok=0; break; }
            ok=qff_insert(&roots,(QuadFieldCell){chart[i],(int32_t)x,(int32_t)y,root_cells})==0;
        }
        Arena_dispose(&scratch);
        if (!ok) goto done;
        if (chunk%100==0 || chunk+1==chunks) fprintf(stderr,"[quad field] root inventory %zu/%zu: %zu sparse roots, %.1fs\n",chunk+1,chunks,roots.count,ves_clock_sec()-started);
    }
    qsort(roots.cell,roots.count,sizeof *roots.cell,qff_compare);
    if (QuadField_new(arena,roots.cell,roots.count,step,&field)!=0) goto done;
    free(roots.cell); roots.cell=NULL; free(roots.hash); roots.hash=NULL;
    for (int level=0;level<=28;level++) {
        QuadFieldReport report={0};
        size_t refine=0,unresolved=0;
        uint8_t *marked=NULL;
        double maximum=0,square=0,mass=0;
        for (size_t chunk=0;chunk<chunks;chunk++) {
            Arena_T scratch=Arena_new();
            MeshBinData mesh={0};
            int32_t *chart=NULL;
            double *weight=NULL;
            int ok=reader(context,scratch,chunk,&mesh,&chart)==0;
            if (ok && mesh.nv && (!mesh.verts || !mesh.uv || !chart)) ok=0;
            if (ok) ok=qff_weights(scratch,&mesh,&weight)==0;
            for (size_t i=0;ok && i<mesh.nv;i++) if (chart[i]>=0) {
                size_t c=0;
                double xyz[3]={mesh.verts[3*i],mesh.verts[3*i+1],mesh.verts[3*i+2]};
                ok=QuadField_locate(field,chart[i],mesh.uv[2*i],mesh.uv[2*i+1],&c)==0 &&
                   QuadField_observe(field,c,mesh.uv[2*i],mesh.uv[2*i+1],xyz,weight[i])==0;
            }
            Arena_dispose(&scratch);
            if (!ok) goto done;
        }
        if (QuadField_solve(field,coarse,&report)!=0) goto done;
        marked=(uint8_t *)ARENA_CALLOC(arena,QuadField_cell_count(field),1);
        for (size_t chunk=0;chunk<chunks;chunk++) {
            Arena_T scratch=Arena_new();
            MeshBinData mesh={0};
            int32_t *chart=NULL;
            double *weight=NULL;
            int ok=reader(context,scratch,chunk,&mesh,&chart)==0;
            if (ok && mesh.nv && (!mesh.verts || !mesh.uv || !chart)) ok=0;
            if (ok) ok=qff_weights(scratch,&mesh,&weight)==0;
            for (size_t i=0;ok && i<mesh.nv;i++) if (chart[i]>=0) {
                size_t c=0;
                double xyz[3]={0},err=0;
                ok=QuadField_locate(field,chart[i],mesh.uv[2*i],mesh.uv[2*i+1],&c)==0 &&
                   QuadField_evaluate(field,c,mesh.uv[2*i],mesh.uv[2*i+1],xyz)==0;
                if (!ok) break;
                for (int k=0;k<3;k++) err+=(xyz[k]-mesh.verts[3*i+k])*(xyz[k]-mesh.verts[3*i+k]);
                square+=weight[i]*err; mass+=weight[i]; err=sqrt(err);
                if (err>maximum) maximum=err;
                if (err>tolerance) {
                    if (QuadField_cell(field,c)->size>1) { if (!marked[c]) refine++; marked[c]=1; }
                    else unresolved++;
                }
            }
            Arena_dispose(&scratch);
            if (!ok) goto done;
        }
        report.maximum_error=maximum; report.rms_error=mass>0 ? sqrt(square/mass) : 0;
        fprintf(stderr,"[quad field] level %d: %zu quads, %zu true DOFs, %zu hanging nodes, rms %.6g max %.6g, refine %zu, unresolved %zu, %.2fs\n",
                level,report.cells,report.unknowns,report.hanging_nodes,report.rms_error,maximum,refine,unresolved,ves_clock_sec()-started);
        if (checkpoint && checkpoint(checkpoint_context,level,field,&report)!=0) goto done;
        if (!refine) {
            *out=field;
            /* A geometric residual failure is not turned into a success by
             * reaching the finest allowed level. Checkpoints remain useful. */
            rc=unresolved ? 1 : 0; goto done;
        }
        coarse=field;
        if (QuadField_refine(arena,coarse,marked,&field)!=0) goto done;
    }
done:
    free(roots.cell); free(roots.hash);
    return rc;
}
