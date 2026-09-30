#include "quad_ribbon_fit.h"
#include "ribbon_lod.h"
#include "ribbon_domains.h"
#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"
#include "../common/json_read.h"
#include "../remesh/intersection_cleanup.h"

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/ves_omp.h"

typedef struct { int32_t chart,x,y,face; } QrfFace;
typedef struct { int32_t a,b; } QrfEdge;
typedef struct { double u,v; } QrfPoint;
typedef struct {
    const MeshBinData *mesh;
    const double *budget; /* optional original-face clearance envelope */
    double step;
    int32_t *chart;
    uint8_t *retained,*anchor,*used;
    QrfFace *faces;
    size_t nfaces,nretained;
    int32_t *contact_faces;
    size_t ncontact_faces;
} QrfSource;
typedef struct {
    QrfSource *source;
    QuadField_T field;
    int32_t face,chart;
    size_t cell;
    QrfPoint uv[3];
    double corner[4][3]; /* current original unit cell, shared by its two faces */
    double origin[3],du[3],dv[3],normal[3],area_ratio;
    uint8_t *marked;
    double *mixed,*rounding;
    double tolerance,maximum,mesh_maximum;
    size_t bad_orientation,unresolved,clearance_failures;
    int observe;
} QrfIntegral;

static int qrf_face_compare(const void *pa,const void *pb)
{
    const QrfFace *a=(const QrfFace *)pa,*b=(const QrfFace *)pb;
    if (a->chart!=b->chart) return a->chart<b->chart ? -1 : 1;
    if (a->y!=b->y) return a->y<b->y ? -1 : 1;
    if (a->x!=b->x) return a->x<b->x ? -1 : 1;
    return a->face<b->face ? -1 : a->face>b->face;
}
static int qrf_edge_compare(const void *pa,const void *pb)
{
    const QrfEdge *a=(const QrfEdge *)pa,*b=(const QrfEdge *)pb;
    if (a->a!=b->a) return a->a<b->a ? -1 : 1;
    return a->b<b->b ? -1 : a->b>b->b;
}
static double qrf_norm(const double a[3])
{ return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]); }
static void qrf_cross(const double a[3],const double b[3],double c[3])
{
    c[0]=a[1]*b[2]-a[2]*b[1]; c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}

static int qrf_affine(QrfIntegral *p,int32_t face)
{
    const MeshBinData *m=p->source->mesh;
    const int32_t *tri=m->faces+3*(size_t)face;
    double e[3]={0},g[3]={0},a=0,b=0,c=0,d=0,det=0;
    p->face=face; p->chart=p->source->chart[tri[0]];
    for (int j=0;j<3;j++) {
        p->uv[j].u=m->uv[2*(size_t)tri[j]];
        p->uv[j].v=m->uv[2*(size_t)tri[j]+1];
    }
    a=p->uv[1].u-p->uv[0].u; b=p->uv[1].v-p->uv[0].v;
    c=p->uv[2].u-p->uv[0].u; d=p->uv[2].v-p->uv[0].v;
    det=a*d-b*c;
    if (!(det>0)) return -1;
    for (int k=0;k<3;k++) {
        p->origin[k]=m->verts[3*(size_t)tri[0]+k];
        e[k]=(double)m->verts[3*(size_t)tri[1]+k]-p->origin[k];
        g[k]=(double)m->verts[3*(size_t)tri[2]+k]-p->origin[k];
        p->du[k]=(d*e[k]-b*g[k])/det;
        p->dv[k]=(a*g[k]-c*e[k])/det;
    }
    qrf_cross(e,g,p->normal); p->area_ratio=qrf_norm(p->normal)/det;
    return isfinite(p->area_ratio) && p->area_ratio>0 ? 0 : -1;
}
static void qrf_target(const QrfIntegral *p,double u,double v,double xyz[3])
{
    for (int k=0;k<3;k++) xyz[k]=p->origin[k]+(u-p->uv[0].u)*p->du[k]+(v-p->uv[0].v)*p->dv[k];
}
static int qrf_barycentric(const MeshBinData *m,int32_t face,double u,double v,double w[3])
{
    const int32_t *t=m->faces+3*(size_t)face;
    double a=(double)m->uv[2*(size_t)t[1]]-m->uv[2*(size_t)t[0]];
    double b=(double)m->uv[2*(size_t)t[1]+1]-m->uv[2*(size_t)t[0]+1];
    double c=(double)m->uv[2*(size_t)t[2]]-m->uv[2*(size_t)t[0]];
    double d=(double)m->uv[2*(size_t)t[2]+1]-m->uv[2*(size_t)t[0]+1];
    double x=u-m->uv[2*(size_t)t[0]],y=v-m->uv[2*(size_t)t[0]+1],det=a*d-b*c;
    if (!(det>0)) return -1;
    w[1]=(x*d-y*c)/det; w[2]=(a*y-b*x)/det; w[0]=1-w[1]-w[2];
    return w[0]>=-1e-10 && w[1]>=-1e-10 && w[2]>=-1e-10 ? 0 : 1;
}

/* Source queries use the original face and its barycentric coordinates,
 * never a nearest-strand choice. Only already-supported regular cells enter
 * this index; duplicate/cut occupancy was excluded by RibbonLod. */
static int qrf_reference(const QrfSource *s,int32_t chart,double u,double v,
                         int32_t *face,double w[3],int32_t *vertex)
{
    double x=u/s->step,y=v/s->step;
    int32_t ix=(int32_t)floor(x),iy=(int32_t)floor(y);
    int nx=x==(double)ix ? 2 : 1,ny=y==(double)iy ? 2 : 1;
    for (int j=0;j<ny;j++) for (int i=0;i<nx;i++) {
        QrfFace key={chart,ix-i,iy-j,INT32_MIN};
        size_t lo=0,hi=s->nfaces;
        while (lo<hi) { size_t mid=lo+(hi-lo)/2; if (qrf_face_compare(s->faces+mid,&key)<0) lo=mid+1; else hi=mid; }
        for (size_t k=lo;k<s->nfaces;k++) {
            const QrfFace *r=s->faces+k;
            if (r->chart!=chart || r->x!=key.x || r->y!=key.y) break;
            if (qrf_barycentric(s->mesh,r->face,u,v,w)==0) {
                *face=r->face; *vertex=-1;
                for (int q=0;q<3;q++) if (fabs(w[q]-1)<1e-12)
                    *vertex=s->mesh->faces[3*(size_t)*face+q];
                return 0;
            }
        }
    }
    return -1;
}

static int qrf_source(Arena_T arena,const MeshBinData *m,double grid,
                       QrfSource *s,QuadFieldCell **cells,size_t *nc)
{
    Arena_T scratch=Arena_new();
    UnionFind uf={0};
    MeshBinData lod={0};
    RibbonLodReport report={0};
    int32_t *original=NULL;
    QrfEdge *edges=NULL;
    int rc=-1;
    if (!m || m->nv>INT32_MAX || m->nf>INT32_MAX ||
        (m->nv && (!m->verts || !m->uv)) || (m->nf && !m->faces)) goto done;
    memset(s,0,sizeof *s); s->mesh=m; s->step=grid;
    s->chart=(int32_t *)ARENA_ALLOC(arena,m->nv*sizeof(int32_t));
    s->retained=(uint8_t *)ARENA_CALLOC(arena,m->nf,1);
    s->anchor=(uint8_t *)ARENA_CALLOC(arena,m->nv,1);
    s->used=(uint8_t *)ARENA_CALLOC(arena,m->nv,1);
    s->faces=(QrfFace *)ARENA_ALLOC(arena,m->nf*sizeof(QrfFace));
    uf=UF_new(scratch,(int32_t)m->nv);
    edges=(QrfEdge *)ARENA_ALLOC(scratch,3*m->nf*sizeof(QrfEdge));
    for (size_t f=0;f<m->nf;f++) for (int k=0;k<3;k++) {
        int32_t a=m->faces[3*f+k],b=m->faces[3*f+(k+1)%3];
        if (a<0 || b<0 || (size_t)a>=m->nv || (size_t)b>=m->nv) goto done;
        uf_union(&uf,a,b); s->used[a]=1;
        edges[3*f+k]=(QrfEdge){a<b ? a : b,a<b ? b : a};
    }
    for (size_t i=0;i<m->nv;i++) s->chart[i]=uf_find(&uf,(int32_t)i);
    qsort(edges,3*m->nf,sizeof *edges,qrf_edge_compare);
    for (size_t i=0;i<3*m->nf;) {
        size_t end=i+1;
        while (end<3*m->nf && qrf_edge_compare(edges+i,edges+end)==0) end++;
        if (end-i!=2)
            s->anchor[edges[i].a]=s->anchor[edges[i].b]=1;
        i=end;
    }
    if (RibbonLod_build(scratch,m,grid,grid,QUAD_RIBBON_ROOT_CELLS,
            QUAD_RIBBON_INITIAL_ERROR,&lod,&original,&report)!=0) goto done;
    *nc=report.nquads;
    *cells=(QuadFieldCell *)ARENA_ALLOC(arena,*nc*sizeof(QuadFieldCell));
    for (size_t i=0;i<*nc;i++) {
        const RibbonLodQuad *q=report.quads+i;
        int32_t chart=s->chart[original[q->corner[0]]];
        for (int k=1;k<4;k++) if (s->chart[original[q->corner[k]]]!=chart) goto done;
        (*cells)[i]=(QuadFieldCell){chart,q->x,q->y,q->cells};
    }
    rc=0;
done:
    Arena_dispose(&scratch); return rc;
}

static int qrf_support(QrfSource *s,QuadField_T field)
{
    const MeshBinData *m=s->mesh;
    for (size_t f=0;f<m->nf;f++) {
        const int32_t *t=m->faces+3*f;
        double u=0,v=0,xmin=DBL_MAX,ymin=DBL_MAX;
        size_t c=0;
        int supported=1;
        for (int j=0;j<3;j++) {
            u+=m->uv[2*(size_t)t[j]]/3.0; v+=m->uv[2*(size_t)t[j]+1]/3.0;
            xmin=fmin(xmin,m->uv[2*(size_t)t[j]]); ymin=fmin(ymin,m->uv[2*(size_t)t[j]+1]);
        }
        if (QuadField_locate(field,s->chart[t[0]],u,v,&c)!=0) supported=0;
        if (supported) {
            const QuadFieldCell *cell=QuadField_cell(field,c);
            for (int j=0;j<3;j++) {
                double x=m->uv[2*(size_t)t[j]]/s->step,y=m->uv[2*(size_t)t[j]+1]/s->step;
                if (x<cell->x || x>cell->x+cell->size || y<cell->y || y>cell->y+cell->size) supported=0;
            }
        }
        if (!supported) {
            s->retained[f]=1; s->nretained++;
            for (int j=0;j<3;j++) s->anchor[t[j]]=1;
        } else {
            QrfIntegral p={0}; p.source=s;
            /* The supported domain is a union of exact original lattice
             * cells. Verify that contract before using the unclipped overlay;
             * irregular/cut faces must have been excluded by RibbonLod. */
            if (xmin/s->step!=floor(xmin/s->step) || ymin/s->step!=floor(ymin/s->step)) return -1;
            for (int j=0;j<3;j++) {
                double x=(m->uv[2*(size_t)t[j]]-xmin)/s->step;
                double y=(m->uv[2*(size_t)t[j]+1]-ymin)/s->step;
                if ((x!=0 && x!=1) || (y!=0 && y!=1)) return -1;
            }
            if (qrf_affine(&p,(int32_t)f)!=0) return -1;
            s->faces[s->nfaces++]=(QrfFace){s->chart[t[0]],
                (int32_t)floor(xmin/s->step),(int32_t)floor(ymin/s->step),(int32_t)f};
        }
    }
    qsort(s->faces,s->nfaces,sizeof(QrfFace),qrf_face_compare);
    return 0;
}

static int qrf_refine_at(QrfSource *s,QuadField_T f,int32_t chart,
                         double u,double v,uint8_t *marked,size_t *count)
{
    int refined=0;
    for (int y=-1;y<=1;y+=2) for (int x=-1;x<=1;x+=2) {
        size_t cell=0;
        if (QuadField_locate(f,chart,u+x*.25*s->step,v+y*.25*s->step,&cell)!=0) continue;
        if (QuadField_cell(f,cell)->size<=1) continue;
        if (!marked[cell]) { marked[cell]=1; (*count)++; }
        refined=1;
    }
    return refined ? 0 : -1;
}

/* All original physical-boundary and exceptional-face vertices remain
 * explicit masters. Refine the adjacent cell when an anchor would otherwise
 * be a hanging constraint; never force mutually inconsistent constraints. */
static int qrf_pins(QrfSource *s,QuadField_T f,uint8_t *marked,size_t *count)
{
    *count=0;
    /* Preserve existing contact triangles, not just their three positions.
     * A coarser triangle spanning them would change the contact signature. */
    for (size_t i=0;i<s->ncontact_faces;i++) {
        int32_t face=s->contact_faces[i];
        const int32_t *tri=s->mesh->faces+3*(size_t)face;
        double u=0,v=0;
        size_t cell=0;
        if (s->retained[face]) continue;
        for (int k=0;k<3;k++) {
            u+=s->mesh->uv[2*(size_t)tri[k]]/3.0;
            v+=s->mesh->uv[2*(size_t)tri[k]+1]/3.0;
        }
        if (QuadField_locate(f,s->chart[tri[0]],u,v,&cell)!=0) return -1;
        if (QuadField_cell(f,cell)->size>1 && !marked[cell]) {
            marked[cell]=1; (*count)++;
        }
    }
    for (size_t i=0;i<s->mesh->nv;i++) if (s->anchor[i]) {
        double u=s->mesh->uv[2*i],v=s->mesh->uv[2*i+1];
        size_t cell=0,node=0;
        QuadFieldNodeKey key={0};
        int found=0;
        if (QuadField_locate(f,s->chart[i],u,v,&cell)!=0) continue;
        if (u/s->step!=round(u/s->step) || v/s->step!=round(v/s->step)) return -1;
        found=QuadField_find_node(f,s->chart[i],(int32_t)round(u/s->step),(int32_t)round(v/s->step),&node)==0;
        if (found && QuadField_node_key(f,node,&key)!=0) return -1;
        if (found && !key.constrained) {
            double xyz[3]={s->mesh->verts[3*i],s->mesh->verts[3*i+1],s->mesh->verts[3*i+2]};
            if (QuadField_fix_node(f,node,xyz)!=0) return -1;
        } else {
            if (qrf_refine_at(s,f,s->chart[i],u,v,marked,count)!=0) return -1;
        }
    }
    return 0;
}

/* At the source sampling limit use its actual two affine triangles. A Q1
 * patch is not the same surface at a diagonal crease, no matter how accurately
 * its four corner samples were fitted. This is a basis choice, not an omitted
 * observation, changed material identity, or post-solve vertex snap. */
static int qrf_elements(const QrfSource *s,QuadField_T f)
{
    for (size_t c=0;c<QuadField_cell_count(f);c++) {
        const QuadFieldCell *key=QuadField_cell(f,c);
        int32_t face=-1,vertex=-1;
        double w[3]={0};
        int mask=0;
        if (key->size!=1) continue;
        if (QuadField_element(f,c)>0) continue;
        if (qrf_reference(s,key->chart,(key->x+.5)*s->step,
                (key->y+.5)*s->step,&face,w,&vertex)!=0) return -1;
        for (int k=0;k<3;k++) {
            size_t id=(size_t)s->mesh->faces[3*(size_t)face+k];
            double x=s->mesh->uv[2*id]/s->step-key->x;
            double y=s->mesh->uv[2*id+1]/s->step-key->y;
            int corner=x==0 ? (y==0 ? 0 : 3) : (y==0 ? 1 : 2);
            if ((x!=0 && x!=1) || (y!=0 && y!=1)) return -1;
            mask|=1<<corner;
        }
        if (mask!=7 && mask!=13 && mask!=11 && mask!=14) return -1;
        if (QuadField_set_element(f,c,mask==7 || mask==13 ? 1 : 2)!=0) return -1;
    }
    return 0;
}

/* At the terminal source resolution the data itself is the detail, not an
 * extrapolated least-squares estimate. Keep independent detail nodes at their
 * original observations. A hanging node stays algebraically constrained (it
 * must not force refinement across the entire connected sheet); the error and
 * orientation checks refine its coarse neighbour only when actually needed.
 * This fixes source data, never a previous coarse iterate. */
static int qrf_details(const QrfSource *s,QuadField_T f)
{
    for (size_t c=0;c<QuadField_cell_count(f);c++) {
        const QuadFieldCell *cell=QuadField_cell(f,c);
        if (!QuadField_element(f,c)) continue;
        for (int y=0;y<2;y++) for (int x=0;x<2;x++) {
            size_t node=0;
            QuadFieldNodeKey key={0};
            int32_t face=-1,vertex=-1;
            double w[3]={0},xyz[3]={0};
            if (QuadField_find_node(f,cell->chart,cell->x+x,cell->y+y,&node)!=0 ||
                QuadField_node_key(f,node,&key)!=0) return -1;
            if (key.constrained) continue;
            if (qrf_reference(s,key.chart,key.x*s->step,key.y*s->step,&face,w,&vertex)!=0 || vertex<0) return -1;
            for (int k=0;k<3;k++) xyz[k]=s->mesh->verts[3*(size_t)vertex+k];
            if (QuadField_fix_node(f,node,xyz)!=0) return -1;
        }
    }
    return 0;
}

/* Six-point degree-four triangle quadrature integrates Q1 normal equations
 * exactly on the original PL/leaf overlay (the area-density is constant on
 * an original source triangle). No refinement-level vertex-vote dependence. */
static int qrf_quadrature(QrfIntegral *p,size_t cell,const QrfPoint *polygon,int n)
{
    static const double a[2]={.445948490915965,.091576213509771};
    static const double weight[2]={.223381589678011,.109951743655322};
    for (int t=1;t+1<n;t++) {
        QrfPoint tri[3]={polygon[0],polygon[t],polygon[t+1]};
        double area=.5*fabs((tri[1].u-tri[0].u)*(tri[2].v-tri[0].v)-(tri[1].v-tri[0].v)*(tri[2].u-tri[0].u));
        if (area<=1e-18) continue;
        for (int family=0;family<2;family++) for (int j=0;j<3;j++) {
            double b[3]={a[family],a[family],a[family]},u=0,v=0,xyz[3]={0};
            double w=weight[family]*area*p->area_ratio;
            b[j]=1-2*a[family];
            for (int k=0;k<3;k++) { u+=b[k]*tri[k].u; v+=b[k]*tri[k].v; }
            qrf_target(p,u,v,xyz);
            if (QuadField_observe(p->field,cell,u,v,xyz,w)!=0) return -1;
        }
    }
    return 0;
}

/* Source cells and fitted leaves share one integer lattice. Every source
 * triangle lies wholly in one leaf; the overlay is the source triangle itself.
 * No sub-grid clipping, recursive descent or extra boundary samples are needed. */
static int qrf_triangle(QrfIntegral *p,int32_t x,int32_t y)
{
    double step=p->source->step,x0=x*step,y0=y*step,x1=(x+1)*step,y1=(y+1)*step;
    /* Retain the old clipper's cyclic ordering so accumulation stays bitwise
     * reproducible, including the six-point quadrature's rounding order. */
    QrfPoint polygon[3]={p->uv[1],p->uv[2],p->uv[0]};
    size_t cell=p->cell;
    const QuadFieldCell *key=QuadField_cell(p->field,cell);
    if (p->observe && qrf_quadrature(p,cell,polygon,3)!=0) return -1;
    if (!p->observe) {
        double maximum=0,du[3]={0},dv[3]={0},cross[3]={0};
        const double (*corner)[3]=p->corner;
        int bad=0;
        if (QuadField_element(p->field,cell)>0) {
            double fitted[3][3]={{0}},dot=0;
            /* The element diagonal matches this source triangle. Their
             * difference is affine: its exact norm bound is at a vertex. */
            for (int i=0;i<3;i++) {
                double target[3]={0},error[3]={0};
                int at=(p->uv[i].u==x1 ? 1 : 0)+(p->uv[i].v==y1 ? 2 : 0);
                memcpy(fitted[i],corner[at],sizeof fitted[i]);
                qrf_target(p,p->uv[i].u,p->uv[i].v,target);
                for (int k=0;k<3;k++) error[k]=fitted[i][k]-target[k];
                maximum=fmax(maximum,qrf_norm(error));
            }
            for (int k=0;k<3;k++) {
                du[k]=fitted[1][k]-fitted[0][k]; dv[k]=fitted[2][k]-fitted[0][k];
            }
            qrf_cross(du,dv,cross);
            for (int k=0;k<3;k++) dot+=cross[k]*p->normal[k];
            bad=!(dot>0);
        } else {
        /* Q1 minus the source triangle's affine extension is bilinear.
         * Its norm on this enclosing rectangle is bounded by its corners. */
        for (int j=0;j<2;j++) for (int i=0;i<2;i++) {
            double u=i ? x1 : x0,v=j ? y1 : y0,target[3]={0},error[3]={0};
            int k=i+2*j;
            qrf_target(p,u,v,target);
            for (int d=0;d<3;d++) error[d]=corner[k][d]-target[d];
            maximum=fmax(maximum,qrf_norm(error));
        }
        /* The Jacobian normal of Q1 is affine in UV. Its dot with the
         * original triangle normal attains its minimum at an original triangle
         * vertex, not at a box corner outside that triangle's support. */
        for (int q=0;q<3;q++) {
            double dot=0;
            double u=(polygon[q].u-x0)/(x1-x0),v=(polygon[q].v-y0)/(y1-y0);
            for (int d=0;d<3;d++) {
                du[d]=(1-v)*(corner[1][d]-corner[0][d])+v*(corner[3][d]-corner[2][d]);
                dv[d]=(1-u)*(corner[2][d]-corner[0][d])+u*(corner[3][d]-corner[1][d]);
            }
            qrf_cross(du,dv,cross);
            for (int d=0;d<3;d++) dot+=cross[d]*p->normal[d];
            if (!(dot>0)) bad=1;
        }
        }
        p->maximum=fmax(p->maximum,maximum);
        maximum+=p->mixed[cell]+p->rounding[cell];
        p->mesh_maximum=fmax(p->mesh_maximum,maximum);
        if (bad) p->bad_orientation++;
        if (p->source->budget && maximum>p->source->budget[p->face]) p->clearance_failures++;
        if (maximum>p->tolerance || bad ||
            (p->source->budget && maximum>p->source->budget[p->face])) {
            if (key->size>1) p->marked[cell]=1;
            else {
                size_t marked=0;
                for (int j=0;j<2;j++) for (int i=0;i<2;i++)
                    (void)qrf_refine_at(p->source,p->field,key->chart,
                        (key->x+i)*step,(key->y+j)*step,p->marked,&marked);
                if (bad && p->unresolved<8) {
                    double target[3]={0};
                    qrf_target(p,(p->uv[0].u+p->uv[1].u+p->uv[2].u)/3,
                        (p->uv[0].v+p->uv[1].v+p->uv[2].v)/3,target);
                    fprintf(stderr,"[global ribbon fit] orientation witness face=%d leaf=(%d,%d,%d,%d) original_center=(%.9g,%.9g,%.9g) mesh_bound=%.9g\n",
                        p->face,key->chart,key->x,key->y,key->size,target[0],target[1],target[2],maximum);
                }
                p->unresolved++;
            }
        }
    }
    return 0;
}

static int qrf_integrate(QrfSource *s,QuadField_T field,int observe,double tolerance,
                          Arena_T scratch,uint8_t *marked,QuadRibbonFitReport *report)
{
    QrfIntegral p={0};
    uint8_t *cached=NULL;
    p.source=s; p.field=field; p.observe=observe; p.marked=marked; p.tolerance=tolerance;
    if (observe) {
        cached=(uint8_t *)ARENA_CALLOC(scratch,QuadField_cell_count(field),1);
        for (size_t c=0;c<QuadField_cell_count(field);c++)
            cached[c]=(uint8_t)(QuadField_observation_count(field,c)>0);
    }
    if (!observe) {
        size_t nc=QuadField_cell_count(field);
        p.mixed=(double *)ARENA_ALLOC(scratch,nc*sizeof(double));
        p.rounding=(double *)ARENA_ALLOC(scratch,nc*sizeof(double));
        for (size_t c=0;c<nc;c++) {
            const QuadFieldCell *key=QuadField_cell(field,c);
            double value[4][3]={{0}},mixed[3]={0},rounding[3]={0};
            for (int j=0;j<2;j++) for (int i=0;i<2;i++) {
                if (QuadField_evaluate(field,c,(key->x+i*key->size)*s->step,
                        (key->y+j*key->size)*s->step,value[i+2*j])!=0) return -1;
            }
            for (int k=0;k<3;k++) {
                mixed[k]=value[0][k]-value[1][k]-value[2][k]+value[3][k];
                if (QuadField_element(field,c)>0) {
                    /* P1 has exactly these vertices and no finer hanging
                     * edge. Exact original float samples have ZERO export
                     * error, essential when preserving an original contact. */
                    for (int i=0;i<4;i++)
                        rounding[k]=fmax(rounding[k],fabs((double)(float)value[i][k]-value[i][k]));
                } else {
                    for (int i=0;i<4;i++) rounding[k]=fmax(rounding[k],fabs(value[i][k]));
                    rounding[k]*=.5*FLT_EPSILON;
                }
            }
            /* A PL interpolation of st on any template triangle in the unit
             * square differs from st by at most 1/4. Include float export. */
            p.mixed[c]=QuadField_element(field,c)>0 ? 0 : .25*qrf_norm(mixed);
            p.rounding[c]=qrf_norm(rounding);
        }
    }
    for (size_t i=0;i<s->nfaces;i++) {
        const QrfFace *face=s->faces+i;
        int new_cell=!i || face->chart!=s->faces[i-1].chart ||
            face->x!=s->faces[i-1].x || face->y!=s->faces[i-1].y;
        if (new_cell) {
            const QuadFieldCell *key=NULL;
            if (QuadField_locate(field,face->chart,(face->x+.5)*s->step,
                    (face->y+.5)*s->step,&p.cell)!=0) return -1;
            key=QuadField_cell(field,p.cell);
            if (key->x>face->x || key->y>face->y ||
                key->x+key->size<face->x+1 || key->y+key->size<face->y+1) return -1;
            if (!observe) for (int y=0;y<2;y++) for (int x=0;x<2;x++)
                if (QuadField_evaluate(field,p.cell,(face->x+x)*s->step,
                        (face->y+y)*s->step,p.corner[x+2*y])!=0) return -1;
        }
        if (observe && cached[p.cell]) continue;
        if (qrf_affine(&p,face->face)!=0 ||
            qrf_triangle(&p,face->x,face->y)!=0) return -1;
    }
    if (!observe) {
        report->continuous_bound=p.maximum; report->mesh_bound=p.mesh_maximum;
        report->rms=report->solve.rms_error;
        report->orientation_failures=p.bad_orientation; report->unresolved_cells=p.unresolved;
        report->clearance_failures=p.clearance_failures;
    }
    return 0;
}

static int qrf_export(Arena_T arena,QrfSource *s,QuadField_T field,QuadRibbonFitOutput *out)
{
    MeshBinData quad={0};
    int32_t *nodes=NULL,*map=NULL;
    size_t extra=0,nv=0,nf=0;
    memset(out,0,sizeof *out);
    if (QuadField_mesh_selection(arena,field,NULL,&quad,&nodes)!=0) return -1;
    map=(int32_t *)ARENA_ALLOC(arena,s->mesh->nv*sizeof(int32_t));
    for (size_t i=0;i<s->mesh->nv;i++) map[i]=-1;
    /* Reserve enough for original fallback/isolated vertices. Unused capacity
     * is temporary arena storage, never serialized as fabricated samples. */
    if (quad.nv+s->mesh->nv>INT32_MAX || quad.nf+s->nretained>INT32_MAX) return -1;
    out->mesh.verts=(float *)ARENA_ALLOC(arena,3*(quad.nv+s->mesh->nv)*sizeof(float));
    out->mesh.uv=(float *)ARENA_ALLOC(arena,2*(quad.nv+s->mesh->nv)*sizeof(float));
    out->mesh.faces=(int32_t *)ARENA_ALLOC(arena,3*(quad.nf+s->nretained)*sizeof(int32_t));
    out->reference_vertex=(int32_t *)ARENA_ALLOC(arena,(quad.nv+s->mesh->nv)*sizeof(int32_t));
    out->reference_face=(int32_t *)ARENA_ALLOC(arena,(quad.nv+s->mesh->nv)*sizeof(int32_t));
    out->reference_barycentric=(double *)ARENA_CALLOC(arena,3*(quad.nv+s->mesh->nv),sizeof(double));
    out->face_cell=(int32_t *)ARENA_ALLOC(arena,(quad.nf+s->nretained)*sizeof(int32_t));
    out->reference_triangle=(int32_t *)ARENA_ALLOC(arena,(quad.nf+s->nretained)*sizeof(int32_t));
    memcpy(out->mesh.verts,quad.verts,3*quad.nv*sizeof(float));
    memcpy(out->mesh.uv,quad.uv,2*quad.nv*sizeof(float));
    memcpy(out->mesh.faces,quad.faces,3*quad.nf*sizeof(int32_t));
    for (size_t i=0;i<quad.nv;i++) {
        QuadFieldNodeKey key={0};
        int32_t face=-1,vertex=-1;
        double w[3]={0};
        if (QuadField_node_key(field,(size_t)nodes[i],&key)!=0 ||
            qrf_reference(s,key.chart,key.x*s->step,key.y*s->step,&face,w,&vertex)!=0) return -1;
        out->reference_face[i]=face; out->reference_vertex[i]=vertex;
        memcpy(out->reference_barycentric+3*i,w,3*sizeof(double));
        if (vertex>=0) {
            if (map[vertex]>=0) return -1;
            map[vertex]=(int32_t)i;
        }
    }
    for (size_t f=0;f<quad.nf;f++) {
        double u=0,v=0;
        size_t c=0;
        int32_t source_face=-1,source_vertex=-1;
        double bary[3]={0};
        int exact=1;
        QuadFieldNodeKey key={0};
        if (QuadField_node_key(field,(size_t)nodes[quad.faces[3*f]],&key)!=0) return -1;
        for (int k=0;k<3;k++) { u+=quad.uv[2*(size_t)quad.faces[3*f+k]]/3.0; v+=quad.uv[2*(size_t)quad.faces[3*f+k]+1]/3.0; }
        if (QuadField_locate(field,key.chart,u,v,&c)!=0) return -1;
        out->face_cell[f]=(int32_t)c;
        if (qrf_reference(s,key.chart,u,v,&source_face,bary,&source_vertex)!=0) return -1;
        for (int k=0;k<3;k++) {
            size_t at=(size_t)quad.faces[3*f+k];
            int32_t original=out->reference_vertex[at];
            int found=0;
            if (original<0) { exact=0; break; }
            for (int j=0;j<3;j++) found|=original==s->mesh->faces[3*(size_t)source_face+j];
            if (!found || memcmp(out->mesh.verts+3*at,s->mesh->verts+3*(size_t)original,3*sizeof(float))) exact=0;
        }
        out->reference_triangle[f]=exact ? source_face : -1;
    }
    for (size_t f=0;f<s->mesh->nf;f++) if (s->retained[f])
        for (int k=0;k<3;k++) if (map[s->mesh->faces[3*f+k]]==-1) map[s->mesh->faces[3*f+k]]=-2;
    for (size_t i=0;i<s->mesh->nv;i++) if (!s->used[i] && map[i]==-1) map[i]=-2;
    for (size_t i=0;i<s->mesh->nv;i++) if (map[i]==-2) {
        size_t at=quad.nv+extra++;
        map[i]=(int32_t)at;
        memcpy(out->mesh.verts+3*at,s->mesh->verts+3*i,3*sizeof(float));
        memcpy(out->mesh.uv+2*at,s->mesh->uv+2*i,2*sizeof(float));
        out->reference_vertex[at]=(int32_t)i; out->reference_face[at]=-1;
    }
    nv=quad.nv+extra; nf=quad.nf;
    for (size_t f=0;f<s->mesh->nf;f++) if (s->retained[f]) {
        for (int k=0;k<3;k++) {
            size_t at=(size_t)map[s->mesh->faces[3*f+k]];
            out->mesh.faces[3*nf+k]=(int32_t)at;
            if (out->reference_face[at]<0) {
                out->reference_face[at]=(int32_t)f;
                out->reference_barycentric[3*at+k]=1;
            }
        }
        out->reference_triangle[nf]=(int32_t)f;
        out->face_cell[nf++]=-1;
    }
    out->mesh.nv=nv; out->mesh.nf=nf;
    return 0;
}

static void qrf_intersection_params(IntersectionCleanupParams *params)
{
    IntersectionCleanup_default_params(params);
    params->gap_max=0; params->parallel_angle_deg=0; params->include_hinges=0;
    params->hit_kind_mask=INTERSECTION_HIT_MASK_STAB;
    /* The parallel visitor streams bounded candidate chunks. No conflict
     * graph or accepted prefix is stored, even on a contact-rich source. */
    params->max_conflicts=SIZE_MAX;
}

static int qrf_mark_contact(size_t a,size_t b,int kind,void *context)
{
    uint8_t *marked=(uint8_t *)context;
    if (kind!=INTERSECTION_HIT_STAB) return -1;
    marked[a]=marked[b]=1;
    return 0;
}

static int qrf_original_contacts(Arena_T arena,QrfSource *s,const double *exact_budget,
                                   QuadRibbonFitReport *report)
{
    Arena_T scratch=Arena_new();
    IntersectionCleanupParams params={0};
    IntersectionCleanupStats stats={0};
    const MeshBinData *m=s->mesh;
    uint8_t *marked=NULL;
    size_t count=0;
    double started=ves_clock_sec();
    int rc=-1;
    if (!m->nf) { rc=0; goto done; }
    marked=(uint8_t *)ARENA_CALLOC(scratch,m->nf,1);
    qrf_intersection_params(&params);
    if (m->nf>=2 && (IntersectionCleanup_audit_visit_parallel(m->verts,m->nv,m->faces,m->nf,
            NULL,&params,NULL,qrf_mark_contact,marked,&stats)!=0 || stats.budget_rejected)) goto done;
    /* Other work batches are absent from the local audit. A zero clearance
     * envelope therefore retains the entire original triangle, not merely
     * corner positions on a potentially different diagonal. */
    if (exact_budget) for (size_t f=0;f<m->nf;f++) if (exact_budget[f]==0) marked[f]=1;
    for (size_t f=0;f<m->nf;f++) count+=marked[f]!=0;
    s->contact_faces=(int32_t *)ARENA_ALLOC(arena,count*sizeof(int32_t));
    for (size_t f=0;f<m->nf;f++) if (marked[f]) {
        s->contact_faces[s->ncontact_faces++]=(int32_t)f;
        for (int k=0;k<3;k++) s->anchor[m->faces[3*f+k]]=1;
    }
    report->source_intersections=stats.stab_pairs;
    fprintf(stderr,"[global ribbon fit] original proper contacts: %zu pairs, %zu exact triangles, %.2fs\n",
        report->source_intersections,count,ves_clock_sec()-started);
    rc=0;
done:
    report->intersection_seconds+=ves_clock_sec()-started;
    Arena_dispose(&scratch); return rc;
}

/* Identity is original vertex provenance AND triangle membership AND XYZ
 * bits, never material/component-label equality or a nearest-UV match. */
static int qrf_original_triangle(const QrfSource *s,QuadField_T field,
                                   const QuadRibbonFitOutput *out,size_t face)
{
    const int32_t *tri=out->mesh.faces+3*face;
    int32_t original[3]={0},cell=out->face_cell[face];
    double u=0,v=0,w[3]={0};
    int32_t reference=-1,vertex=-1;
    for (int k=0;k<3;k++) {
        size_t id=(size_t)tri[k];
        original[k]=out->reference_vertex[id];
        if (original[k]<0 || (size_t)original[k]>=s->mesh->nv ||
            memcmp(out->mesh.verts+3*id,s->mesh->verts+3*(size_t)original[k],3*sizeof(float))) return 0;
        u+=out->mesh.uv[2*id]/3.0; v+=out->mesh.uv[2*id+1]/3.0;
    }
    /* Negative cells are verbatim exceptional source triangles from export. */
    if (cell<0) return 1;
    if (qrf_reference(s,QuadField_cell(field,(size_t)cell)->chart,u,v,
            &reference,w,&vertex)!=0) return 0;
    for (int k=0;k<3;k++) {
        int found=0;
        for (int j=0;j<3;j++) found|=original[k]==s->mesh->faces[3*(size_t)reference+j];
        if (!found) return 0;
    }
    return 1;
}

typedef struct {
    QrfSource *source;
    QuadField_T field;
    const QuadRibbonFitOutput *output;
    uint8_t *marked;
    size_t failures;
} QrfCollision;

static int qrf_refine_collision(size_t a,size_t b,int kind,void *context)
{
    QrfCollision *p=(QrfCollision *)context;
    const QuadRibbonFitOutput *out=p->output;
    size_t faces[2]={a,b};
    int changed[2]={0};
    if (kind!=INTERSECTION_HIT_STAB) return -1;
    for (int j=0;j<2;j++) changed[j]=!qrf_original_triangle(p->source,p->field,out,faces[j]);
    if (!changed[0] && !changed[1]) return 0;
    p->failures++;
    if (p->failures<=8)
        fprintf(stderr,"[global ribbon fit] new-contact witness faces=%zu/%zu cells=%d/%d changed=%d/%d\n",
            a,b,out->face_cell[a],out->face_cell[b],changed[0],changed[1]);
    for (int j=0;j<2;j++) if (changed[j]) {
        int32_t cell=out->face_cell[faces[j]];
        const QuadFieldCell *key=NULL;
        if (cell<0) return -1; /* retained source data must never have moved */
        key=QuadField_cell(p->field,(size_t)cell);
        if (key->size>1) p->marked[cell]=1;
        else for (int k=0;k<3;k++) {
            size_t id=(size_t)out->mesh.faces[3*faces[j]+k],count=0;
            (void)qrf_refine_at(p->source,p->field,key->chart,
                out->mesh.uv[2*id],out->mesh.uv[2*id+1],p->marked,&count);
        }
    }
    return 0;
}

static int qrf_collisions(QrfSource *source,QuadField_T field,
                            const QuadRibbonFitOutput *output,uint8_t *marked,
                            QuadRibbonFitReport *report)
{
    IntersectionCleanupParams params={0};
    IntersectionCleanupStats stats={0};
    QrfCollision context={source,field,output,marked,0};
    double started=ves_clock_sec();
    int rc=0;
    qrf_intersection_params(&params);
    if (output->mesh.nf>=2) rc=IntersectionCleanup_audit_visit_parallel(output->mesh.verts,output->mesh.nv,
        output->mesh.faces,output->mesh.nf,NULL,&params,NULL,
        qrf_refine_collision,&context,&stats);
    report->intersection_seconds+=ves_clock_sec()-started;
    report->intersection_passes++;
    report->output_intersections=stats.stab_pairs;
    report->intersection_failures=context.failures;
    report->intersection_refinements+=context.failures;
    /* Original contacts were pinned before solving. Their disappearance is
     * a failure, not an improvement in the fitted pair count. */
    if (!context.failures && stats.stab_pairs!=report->source_intersections) rc=-1;
    fprintf(stderr,"[global ribbon fit] emitted proper contacts: %zu original, %zu fitted, %zu changed pairs, %.2fs\n",
        report->source_intersections,stats.stab_pairs,context.failures,ves_clock_sec()-started);
    return rc || stats.budget_rejected ? -1 : 0;
}

static int qrf_run(const MeshBinData *input,double grid,const double *budget,int partitioned,
                     QuadRibbonFitCheckpoint checkpoint,void *context)
{
    Arena_T arena=Arena_new(),current_arena=Arena_new(),previous_arena=NULL;
    QrfSource source={0};
    QuadFieldCell *cells=NULL;
    QuadField_T field=NULL,previous=NULL;
    QuadRibbonFitReport report={0};
    size_t nc=0,iteration=0;
    double started=ves_clock_sec();
    int rc=-1,solved=0;
    if (!isfinite(grid) || grid<=0) goto done;
    if (qrf_source(arena,input,grid,&source,&cells,&nc)!=0) {
        fprintf(stderr,"[global ribbon fit] source hierarchy construction failed\n"); goto done;
    }
    fprintf(stderr,"[global ribbon fit] initial supported domain: %zu cells, %.2fs\n",nc,ves_clock_sec()-started);
    if (QuadField_new(current_arena,cells,nc,source.step,&field)!=0) {
        if (nc>QUAD_FIELD_MAX_CELLS) {
            QrfEdge *domains=(QrfEdge *)ARENA_ALLOC(arena,nc*sizeof *domains);
            size_t ndomains=0,largest=0,oversized=0;
            int32_t largest_chart=-1;
            for (size_t i=0;i<nc;i++) domains[i]=(QrfEdge){cells[i].chart,0};
            qsort(domains,nc,sizeof *domains,qrf_edge_compare);
            for (size_t i=0;i<nc;) {
                size_t end=i+1;
                while (end<nc && domains[end].a==domains[i].a) end++;
                ndomains++; oversized+=end-i>QUAD_FIELD_MAX_CELLS;
                if (end-i>largest) { largest=end-i; largest_chart=domains[i].a; }
                i=end;
            }
            fprintf(stderr,"[global ribbon fit] capacity diagnostic: %zu supported source-connected domains; largest has %zu cells (source root %d); %zu domains exceed the cell cap\n",
                ndomains,largest,largest_chart,oversized);
        }
        fprintf(stderr,"[global ribbon fit] initial finite-element construction failed\n"); goto done;
    }
    if (qrf_support(&source,field)!=0) {
        fprintf(stderr,"[global ribbon fit] original support indexing failed\n"); goto done;
    }
    if (qrf_original_contacts(arena,&source,partitioned ? budget : NULL,&report)!=0) {
        fprintf(stderr,"[global ribbon fit] original contact inventory failed\n"); goto done;
    }
    source.budget=budget;
    if (budget) for (size_t f=0;f<input->nf;f++) {
        if (!isfinite(budget[f]) || budget[f]<0) goto done;
        report.clearance_faces+=budget[f]<QUAD_RIBBON_INITIAL_ERROR;
        report.exact_contact_faces+=budget[f]==0;
    }
    fprintf(stderr,"[global ribbon fit] support: %zu original triangles, %zu quad leaves, %zu retained exceptional triangles, %.2fs\n",
        input->nf,nc,source.nretained,ves_clock_sec()-started);
    report.input_faces=input->nf; report.retained_faces=source.nretained;
    for (int level=0;level<QUAD_RIBBON_LEVELS;level++) {
        double tolerance=ldexp(QUAD_RIBBON_INITIAL_ERROR,-level);
        for (;;) {
            Arena_T scratch=Arena_new();
            uint8_t *marked=(uint8_t *)ARENA_CALLOC(scratch,QuadField_cell_count(field),1);
            size_t count=0;
            int ok=1;
            if (!solved) {
                ok=qrf_elements(&source,field)==0 && qrf_pins(&source,field,marked,&count)==0;
                if (ok && count) {
                    Arena_T next=Arena_new(); QuadField_T refined=NULL;
                    ok=QuadField_refine(next,field,marked,&refined)==0;
                    if (ok) { Arena_dispose(&current_arena); current_arena=next; field=refined; }
                    else Arena_dispose(&next);
                    Arena_dispose(&scratch);
                    if (!ok) goto done;
                    continue;
                }
                if (ok) ok=qrf_details(&source,field)==0;
                {
                    double started_assembly=ves_clock_sec();
                    if (ok && previous) ok=QuadField_reuse_observations(field,previous)==0;
                    if (ok) ok=qrf_integrate(&source,field,1,tolerance,scratch,marked,&report)==0;
                    report.assembly_seconds+=ves_clock_sec()-started_assembly;
                }
                {
                    double started_solve=ves_clock_sec();
                    if (ok) {
                        /* GENMMD and other legacy TAUCS helpers keep mutable
                         * static state. Keep this complete factor/solve/free
                         * transaction serial while independent source jobs
                         * assemble, refine, validate and publish concurrently. */
#ifdef _OPENMP
#pragma omp critical(qrf_taucs_transaction)
#endif
                        { ok=QuadField_solve(field,previous,&report.solve)==0; }
                    }
                    report.solve_seconds+=ves_clock_sec()-started_solve;
                }
                if (ok) { iteration++; solved=1; }
            }
            {
                double started_validation=ves_clock_sec();
                if (ok) ok=qrf_integrate(&source,field,0,tolerance,scratch,marked,&report)==0;
                report.validation_seconds+=ves_clock_sec()-started_validation;
            }
            if (!ok) { Arena_dispose(&scratch); goto done; }
            for (size_t i=0;i<QuadField_cell_count(field);i++) if (marked[i]) count++;
            report.iteration=iteration; report.tolerance=tolerance; report.elapsed_seconds=ves_clock_sec()-started;
            fprintf(stderr,"[global ribbon fit] level %d solve %zu: %zu quads, %zu free + %zu fixed DOFs, Q1 bound %.6g, mesh bound %.6g, rms %.6g, refine %zu, orientation %zu, clearance %zu, unresolved %zu, %.2fs\n",
                level,iteration,report.solve.cells,report.solve.unknowns,report.solve.fixed_nodes,
                report.continuous_bound,report.mesh_bound,report.rms,count,report.orientation_failures,
                report.clearance_failures,report.unresolved_cells,report.elapsed_seconds);
            if (!count) {
                QuadRibbonFitOutput output={0};
                if (report.unresolved_cells || report.orientation_failures || report.clearance_failures) {
                    ok=qrf_export(scratch,&source,field,&output)==0;
                    report.output_faces=output.mesh.nf; report.output_vertices=output.mesh.nv;
                    if (ok && checkpoint) ok=checkpoint(context,-level-1,field,&output,&report)==0;
                    Arena_dispose(&scratch); rc=1; goto done;
                }
                ok=qrf_export(scratch,&source,field,&output)==0;
                report.output_faces=output.mesh.nf; report.output_vertices=output.mesh.nv;
                if (ok) ok=qrf_collisions(&source,field,&output,marked,&report)==0;
                for (size_t i=0;i<QuadField_cell_count(field);i++) if (marked[i]) count++;
                report.elapsed_seconds=ves_clock_sec()-started;
                if (ok && !count && report.intersection_failures) ok=0;
                if (ok && !count && checkpoint) ok=checkpoint(context,level,field,&output,&report)==0;
                if (!ok || !count) {
                    Arena_dispose(&scratch);
                    if (!ok) goto done;
                    break;
                }
            }
            {
                Arena_T next=Arena_new(); QuadField_T refined=NULL;
                ok=QuadField_refine(next,field,marked,&refined)==0;
                if (ok) {
                    if (previous_arena) Arena_dispose(&previous_arena);
                    previous_arena=current_arena; previous=field;
                    current_arena=next; field=refined; solved=0;
                } else Arena_dispose(&next);
            }
            Arena_dispose(&scratch);
            if (!ok) goto done;
        }
    }
    rc=0;
done:
    fprintf(stderr,"[global ribbon fit] timings: assembly %.3fs, solve %.3fs, validation %.3fs, intersections %.3fs\n",
        report.assembly_seconds,report.solve_seconds,report.validation_seconds,report.intersection_seconds);
    if (rc) fprintf(stderr,"[global ribbon fit] incomplete rc=%d after %.2fs; no source or accepted checkpoint changed\n",rc,ves_clock_sec()-started);
    if (previous_arena) Arena_dispose(&previous_arena);
    Arena_dispose(&current_arena); Arena_dispose(&arena);
    return rc;
}

int QuadRibbonFit_run(const MeshBinData *input,double grid,
    QuadRibbonFitCheckpoint checkpoint,void *context)
{
    return qrf_run(input,grid,NULL,0,checkpoint,context);
}

typedef struct {
    const char *directory;
    uint64_t fingerprint;
    const MeshBinData *source;
    const char *input;
    uint64_t budget_fingerprint;
    const int32_t *source_vertices,*source_faces; /* part-local -> original */
    size_t original_vertices;
    QuadRibbonFitReport *reports; /* optional four-checkpoint inventory */
    size_t *p1_cells;
    unsigned *completed;
    const uint8_t *const *sidecars; /* optional original-global read cache */
} QrfWriter;
static const char *qrf_sidecar_suffix[]={"_material_identity.i32","_phase.f32","_support.u8",
    "_reconstruction_component.i32","_claimant_chart.i32","_claimant_island.i32"};
static int qrf_write_array(const char *path,const void *data,size_t width,size_t count)
{
    FILE *f=fopen(path,"wb");
    int ok=0;
    if (!f) return -1;
    ok=fwrite(data,width,count,f)==count;
    if (fclose(f)!=0) ok=0;
    return ok ? 0 : -1;
}

static int qrf_write_indices(const char *path,const int32_t *data,size_t count,
                               const int32_t *map,size_t map_count)
{
    int32_t buffer[1024]={0};
    FILE *f=fopen(path,"wb");
    int ok=f!=NULL;
    if (!f) return -1;
    for (size_t start=0;ok && start<count;) {
        size_t n=count-start<1024 ? count-start : 1024;
        for (size_t i=0;i<n;i++) {
            int32_t id=data[start+i];
            if (id < -1 || (id>=0 && (size_t)id>=map_count)) { ok=0; break; }
            buffer[i]=id>=0 && map ? map[id] : id;
        }
        if (ok) ok=fwrite(buffer,sizeof *buffer,n,f)==n;
        start+=n;
    }
    if (fclose(f)!=0) ok=0;
    return ok ? 0 : -1;
}

/* Preserve carried observation identity and support at the exact source
 * vertex. No label is inferred from spatial proximity or field chart number. */
static int qrf_input_stem(const char *input,char *stem,size_t capacity)
{
    size_t len=strlen(input);
    if (len<6 || len>=capacity || strcmp(input+len-6,".vmesh")) return -1;
    memcpy(stem,input,len-6); stem[len-6]=0; len-=6;
    if (len>=6 && strcmp(stem+len-6,"_world")==0) stem[len-6]=0;
    return 0;
}

static int qrf_sidecars(const QrfWriter *writer,const char *dir,const QuadRibbonFitOutput *out)
{
    Arena_T arena=Arena_new();
    char stem[1200],path[1400];
    int rc=-1;
    size_t source_vertices=writer->source_vertices ? writer->original_vertices : writer->source->nv;
    if (qrf_input_stem(writer->input,stem,sizeof stem)!=0) goto done;
    for (size_t j=0;j<sizeof qrf_sidecar_suffix/sizeof qrf_sidecar_suffix[0];j++) {
        size_t width=j==2 ? 1 : 4;
        const uint8_t *source=NULL;
        uint8_t *selected=NULL;
        Arena_Mark mark=Arena_save(arena);
        if (writer->sidecars) {
            source=writer->sidecars[j];
            if (!source) continue;
        } else {
            uint8_t *loaded=NULL;
            FILE *f=NULL;
            int ok=1;
            snprintf(path,sizeof path,"%s%s",stem,qrf_sidecar_suffix[j]);
            f=fopen(path,"rb");
            if (!f) { if (errno==ENOENT) continue; goto done; }
            loaded=(uint8_t *)ARENA_ALLOC(arena,width*source_vertices);
            ok=fread(loaded,width,source_vertices,f)==source_vertices && fgetc(f)==EOF;
            if (fclose(f)!=0) ok=0;
            if (!ok) goto done;
            source=loaded;
        }
        selected=(uint8_t *)ARENA_ALLOC(arena,width*out->mesh.nv);
        for (size_t i=0;i<out->mesh.nv;i++) {
            int32_t original=out->reference_vertex[i];
            if (original<0 || (size_t)original>=writer->source->nv) goto done;
            if (writer->source_vertices) original=writer->source_vertices[original];
            if (original<0 || (size_t)original>=source_vertices) goto done;
            memcpy(selected+width*i,source+width*(size_t)original,width);
        }
        snprintf(path,sizeof path,"%s/fit_ribbon%s",dir,qrf_sidecar_suffix[j]);
        if (qrf_write_array(path,selected,width,out->mesh.nv)!=0) goto done;
        Arena_restore(arena,mark);
    }
    rc=0;
done:
    Arena_dispose(&arena); return rc;
}
static int qrf_checkpoint(void *context,int level,QuadField_T field,
                           const QuadRibbonFitOutput *output,const QuadRibbonFitReport *r)
{
    QrfWriter *writer=(QrfWriter *)context;
    char dir[1200],path[1400];
    FILE *f=NULL;
    size_t p1=0;
    for (size_t c=0;c<QuadField_cell_count(field);c++) if (QuadField_element(field,c)>0) p1++;
    if (level<0) snprintf(dir,sizeof dir,"%s/failed_level%d",writer->directory,-level-1);
    else snprintf(dir,sizeof dir,"%s/level%d",writer->directory,level);
    if (ves_mkdir(dir)!=0 && errno!=EEXIST) return -1;
    snprintf(path,sizeof path,"%s/fit_ribbon_world.vmesh",dir);
    if (MeshBin_write(path,output->mesh.verts,output->mesh.nv,output->mesh.faces,output->mesh.nf,output->mesh.uv)!=0) return -1;
    snprintf(path,sizeof path,"%s/field.qfield",dir);
    if (QuadField_save(field,path,writer->fingerprint)!=0) return -1;
    snprintf(path,sizeof path,"%s/reference_vertex.i32",dir);
    if (qrf_write_indices(path,output->reference_vertex,output->mesh.nv,writer->source_vertices,writer->source->nv)!=0) return -1;
    snprintf(path,sizeof path,"%s/reference_face.i32",dir);
    if (qrf_write_indices(path,output->reference_face,output->mesh.nv,writer->source_faces,writer->source->nf)!=0) return -1;
    snprintf(path,sizeof path,"%s/reference_triangle.i32",dir);
    if (qrf_write_indices(path,output->reference_triangle,output->mesh.nf,writer->source_faces,writer->source->nf)!=0) return -1;
    snprintf(path,sizeof path,"%s/reference_barycentric.f64",dir);
    if (qrf_write_array(path,output->reference_barycentric,3*sizeof(double),output->mesh.nv)!=0) return -1;
    snprintf(path,sizeof path,"%s/face_cell.i32",dir);
    if (qrf_write_array(path,output->face_cell,sizeof(int32_t),output->mesh.nf)!=0) return -1;
    if (qrf_sidecars(writer,dir,output)!=0) return -1;
    snprintf(path,sizeof path,"%s/report.json",dir); f=fopen(path,"wb");
    if (!f) return -1;
    fprintf(f,"{\"schema\":\"supported-global-quad-fit-v2\",\"level\":%d,\"input_fingerprint\":\"%016llx\",\"input_faces\":%zu,\"retained_faces\":%zu,\"vertices\":%zu,\"faces\":%zu,\"quads\":%zu,\"free_dofs\":%zu,\"fixed_dofs\":%zu,\"solves\":%zu,\"tolerance\":%.17g,\"continuous_bound\":%.17g,\"mesh_bound\":%.17g,\"rms\":%.17g,\"orientation_failures\":%zu,\"unresolved\":%zu,\"elapsed_seconds\":%.6f,",
        level,(unsigned long long)writer->fingerprint,r->input_faces,r->retained_faces,r->output_vertices,r->output_faces,
        r->solve.cells,r->solve.unknowns,r->solve.fixed_nodes,r->iteration,r->tolerance,r->continuous_bound,r->mesh_bound,r->rms,
        r->orientation_failures,r->unresolved_cells,r->elapsed_seconds);
    fprintf(f,"\"clearance_budget_fingerprint\":\"%016llx\",\"clearance_faces\":%zu,\"exact_contact_faces\":%zu,\"clearance_failures\":%zu,",
        (unsigned long long)writer->budget_fingerprint,r->clearance_faces,r->exact_contact_faces,r->clearance_failures);
    fprintf(f,"\"source_intersections\":%zu,\"output_intersections\":%zu,\"intersection_failures\":%zu,\"intersection_passes\":%zu,\"intersection_refinements\":%zu,\"intersection_seconds\":%.6f,",
        r->source_intersections,r->output_intersections,r->intersection_failures,
        r->intersection_passes,r->intersection_refinements,r->intersection_seconds);
    fprintf(f,"\"field_nodes\":%zu,\"node_storage_bytes\":%zu,",
        r->solve.nodes,r->solve.node_storage_bytes);
    fprintf(f,"\"source_triangle_cells\":%zu,\"bilinear_cells\":%zu,\"timings\":{\"assembly\":%.6f,\"solve\":%.6f,\"validation\":%.6f},\"finest_policy\":\"original source diagonal and independent source samples; hanging nodes remain constrained\",\"qualification\":\"geometric approximation only; physical intersections and material correspondence require independent audit\"}\n",
        p1,QuadField_cell_count(field)-p1,r->assembly_seconds,r->solve_seconds,r->validation_seconds);
    if (fclose(f)!=0) return -1;
    if (level>=0 && level<QUAD_RIBBON_LEVELS && writer->reports) {
        writer->reports[level]=*r;
        writer->p1_cells[level]=p1;
        *writer->completed|=1u<<level;
    }
    return 0;
}
static uint64_t qrf_hash(uint64_t h,const void *data,size_t bytes)
{
    const uint8_t *p=(const uint8_t *)data;
    for (size_t i=0;i<bytes;i++) { h^=p[i]; h*=UINT64_C(1099511628211); }
    return h;
}

typedef struct {
    QuadRibbonFitReport report[QUAD_RIBBON_LEVELS];
    size_t p1[QUAD_RIBBON_LEVELS];
    size_t input_bytes;
    uint64_t fingerprint;
    unsigned completed;
} QrfPartRecord;

static int qrf_read_array(const char *path,void *data,size_t width,size_t count)
{
    FILE *f=fopen(path,"rb");
    int ok=0;
    if (!f) return -1;
    ok=fread(data,width,count,f)==count && fgetc(f)==EOF;
    if (fclose(f)!=0) ok=0;
    return ok ? 0 : -1;
}

static int qrf_cache_sidecars(Arena_T arena,const QrfWriter *writer,const uint8_t **cache)
{
    char stem[1200],path[1400];
    size_t bytes=0;
    if (qrf_input_stem(writer->input,stem,sizeof stem)!=0) return -1;
    for (size_t j=0;j<sizeof qrf_sidecar_suffix/sizeof qrf_sidecar_suffix[0];j++) {
        size_t width=j==2 ? 1 : 4;
        uint8_t *loaded=NULL;
        FILE *f=NULL;
        int ok=1;
        cache[j]=NULL;
        snprintf(path,sizeof path,"%s%s",stem,qrf_sidecar_suffix[j]);
        f=fopen(path,"rb");
        if (!f) { if (errno==ENOENT) continue; return -1; }
        loaded=(uint8_t *)ARENA_ALLOC(arena,width*writer->source->nv);
        ok=fread(loaded,width,writer->source->nv,f)==writer->source->nv && fgetc(f)==EOF;
        if (fclose(f)!=0) ok=0;
        if (!ok) return -1;
        cache[j]=loaded; bytes+=width*writer->source->nv;
    }
    fprintf(stderr,"[global ribbon bank] source sidecar cache: %zu bytes, one read per member\n",bytes);
    return 0;
}

static int qrf_domain_capacity(RibbonDomains_T domains,size_t cap)
{
    RibbonDomainsReport report={0};
    size_t refused=0,maximum=0;
    if (RibbonDomains_report(domains,&report)!=0) return -1;
    /* A one-part fit borrows the original arrays. Multi-part fits copy one
     * intact part; reject its measured storage BEFORE expensive clearance. */
    if (report.parts>1) for (size_t p=0;p<report.parts;p++) {
        RibbonDomainPart part={0};
        if (RibbonDomains_part(domains,p,&part)!=0) return -1;
        if (part.input_bytes>maximum) maximum=part.input_bytes;
        if (part.input_bytes>cap) {
            fprintf(stderr,"[global ribbon fit] input-copy preflight refused part %zu/%zu: %zu faces / %zu vertices require %zu bytes, cap %zu; no clearance or solve attempted\n",
                p+1,report.parts,part.faces,part.vertices,part.input_bytes,
                cap);
            refused++;
        }
    }
    fprintf(stderr,"[global ribbon fit] input-copy preflight: %zu parts, %zu refused, %zu maximum active bytes, %zu cap\n",
        report.parts,refused,maximum,cap);
    /* Positive means a measured capacity refusal, not invalid source data.
     * The atlas can finish inventorying its other members without solving. */
    return refused ? 1 : 0;
}

static int qrf_domain_plan(Arena_T arena,const MeshBinData *mesh,RibbonDomains_T *out)
{
    RibbonDomainsReport report={0};
    if (RibbonDomains_new(arena,mesh,QUAD_RIBBON_BATCH_FACES,
            QUAD_RIBBON_BATCH_VERTICES,out)!=0 || RibbonDomains_report(*out,&report)!=0) return -1;
    fprintf(stderr,"[global ribbon fit] whole-component work plan: %zu components, %zu parts, %zu oversized; largest component %zu faces / %zu vertices\n",
        report.components,report.parts,report.oversized_parts,
        report.largest_component_faces,report.largest_component_vertices);
    return qrf_domain_capacity(*out,QUAD_RIBBON_BATCH_INPUT_BYTES);
}

static int qrf_partition_clearance(const MeshBinData *mesh,RibbonDomains_T domains,double *budget)
{
    RibbonDomainsReport report={0};
    size_t tested=0;
    double started=ves_clock_sec();
    if (RibbonDomains_report(domains,&report)!=0) return -1;
    if (report.parts<=1) return 0;
    if (IntersectionCleanup_partition_budget(mesh->verts,mesh->nv,mesh->faces,mesh->nf,
            RibbonDomains_face_owners(domains),QUAD_RIBBON_INITIAL_ERROR,budget,&tested)!=0) return -1;
    fprintf(stderr,"[global ribbon fit] cross-part clearance: %zu tested pairs, %.2fs\n",tested,ves_clock_sec()-started);
    return 0;
}

/* Join one fixed-width channel without holding all output samples in RAM.
 * Optional carried channels must be present in EVERY part or absent in all.
 * Only field-cell IDs are rebased; source IDs are already original-global. */
static int qrf_merge_channel(const char *root,int level,const QrfPartRecord *parts,
                               size_t count,const char *name,size_t width,
                               int per_face,int optional,int rebase)
{
    uint8_t buffer[4096]={0};
    char path[1400];
    FILE *output=NULL,*input=NULL;
    size_t offset=0;
    int present=-1,rc=-1;
    if (!width || width>sizeof buffer || (rebase && width!=sizeof(int32_t))) return -1;
    for (size_t p=0;p<count;p++) {
        size_t n=per_face ? parts[p].report[level].output_faces : parts[p].report[level].output_vertices;
        snprintf(path,sizeof path,"%s/parts/part%05zu/level%d/%s",root,p,level,name);
        input=fopen(path,"rb");
        if (!input) {
            if (optional && errno==ENOENT && present!=1) { present=0; continue; }
            goto done;
        }
        if (present==0) goto done;
        if (present<0) {
            present=1;
            snprintf(path,sizeof path,"%s/level%d/%s",root,level,name);
            output=fopen(path,"wb"); if (!output) goto done;
        }
        while (n) {
            size_t take=n<sizeof buffer/width ? n : sizeof buffer/width;
            if (fread(buffer,width,take,input)!=take) goto done;
            if (rebase) for (size_t i=0;i<take;i++) {
                int32_t id=0;
                memcpy(&id,buffer+i*width,sizeof id);
                if (id < -1 || (id>=0 && ((size_t)id>=parts[p].report[level].solve.cells ||
                        offset>INT32_MAX-(size_t)id))) goto done;
                if (id>=0) id+=(int32_t)offset;
                memcpy(buffer+i*width,&id,sizeof id);
            }
            if (fwrite(buffer,width,take,output)!=take) goto done;
            n-=take;
        }
        if (fgetc(input)!=EOF || ferror(input)) goto done;
        if (fclose(input)!=0) { input=NULL; goto done; }
        input=NULL;
        if (rebase) offset+=parts[p].report[level].solve.cells;
    }
    rc=0;
done:
    if (input && fclose(input)!=0) rc=-1;
    if (output && fclose(output)!=0) rc=-1;
    return rc;
}

typedef struct {
    const MeshBinData *source,*mesh;
    const int32_t *vertex,*triangle;
    size_t failures;
} QrfBankAudit;

static int qrf_bank_contact(size_t a,size_t b,int kind,void *context)
{
    QrfBankAudit *audit=(QrfBankAudit *)context;
    size_t pair[2]={a,b};
    int exact=kind==INTERSECTION_HIT_STAB;
    for (int j=0;j<2 && exact;j++) {
        int32_t source_face=audit->triangle[pair[j]];
        if (source_face<0 || (size_t)source_face>=audit->source->nf) { exact=0; break; }
        for (int k=0;k<3;k++) {
            size_t vertex=(size_t)audit->mesh->faces[3*pair[j]+k];
            int32_t original=audit->vertex[vertex];
            int found=0;
            if (original<0 || (size_t)original>=audit->source->nv) { exact=0; break; }
            for (int t=0;t<3;t++) found|=original==audit->source->faces[3*(size_t)source_face+t];
            if (!found || memcmp(audit->mesh->verts+3*vertex,
                    audit->source->verts+3*(size_t)original,3*sizeof(float))) exact=0;
        }
    }
    if (!exact) {
        audit->failures++;
        if (audit->failures<=8) fprintf(stderr,"[global ribbon bank] changed contact: output triangles %zu/%zu, source triangles %d/%d\n",
            a,b,audit->triangle[a],audit->triangle[b]);
    }
    return 0;
}

static int qrf_bank_audit(const QrfWriter *writer,int level,size_t original_pairs,size_t *pairs)
{
    Arena_T arena=Arena_new();
    MeshBinData mesh={0};
    IntersectionCleanupParams params={0};
    IntersectionCleanupStats stats={0};
    QrfBankAudit audit={writer->source,&mesh,NULL,NULL,0};
    int32_t *vertices=NULL,*triangles=NULL;
    uint8_t *seen=NULL;
    char path[1400];
    double started=ves_clock_sec();
    int rc=-1;
    snprintf(path,sizeof path,"%s/level%d/fit_ribbon_world.vmesh",writer->directory,level);
    if (MeshBin_read_arena(arena,path,&mesh)!=0) goto done;
    vertices=(int32_t *)ARENA_ALLOC(arena,mesh.nv*sizeof *vertices);
    triangles=(int32_t *)ARENA_ALLOC(arena,mesh.nf*sizeof *triangles);
    seen=(uint8_t *)ARENA_CALLOC(arena,writer->source->nf,1);
    snprintf(path,sizeof path,"%s/level%d/reference_vertex.i32",writer->directory,level);
    if (qrf_read_array(path,vertices,sizeof *vertices,mesh.nv)!=0) goto done;
    snprintf(path,sizeof path,"%s/level%d/reference_triangle.i32",writer->directory,level);
    if (qrf_read_array(path,triangles,sizeof *triangles,mesh.nf)!=0) goto done;
    for (size_t f=0;f<mesh.nf;f++) {
        int32_t id=triangles[f];
        if (id < -1 || (id>=0 && ((size_t)id>=writer->source->nf || seen[id]))) goto done;
        if (id>=0) seen[id]=1;
    }
    audit.vertex=vertices; audit.triangle=triangles;
    qrf_intersection_params(&params);
    if (mesh.nf>=2 && (IntersectionCleanup_audit_visit_parallel(mesh.verts,mesh.nv,mesh.faces,mesh.nf,
            NULL,&params,NULL,qrf_bank_contact,&audit,&stats)!=0 || stats.budget_rejected)) goto done;
    *pairs=stats.stab_pairs;
    /* Every fitted pair is an exact original pair, and original triangles
     * occur at most once. Equal cardinality therefore proves set equality
     * without retaining a potentially unbounded pair graph. */
    if (audit.failures || stats.stab_pairs!=original_pairs) goto done;
    rc=0;
done:
    fprintf(stderr,"[global ribbon bank] merged level %d contacts: %zu original, %zu fitted, %zu changed; %s, %.2fs\n",
        level,original_pairs,stats.stab_pairs,audit.failures,rc ? "FAILED" : "preserved",ves_clock_sec()-started);
    Arena_dispose(&arena); return rc;
}

static int qrf_bank_manifest(const QrfWriter *writer,const QrfPartRecord *part,size_t count,
                              int level,size_t pairs,double elapsed)
{
    QuadRibbonFitReport total={0};
    size_t p1=0,max_cells=0,max_free=0,max_nodes=0,max_input=0,offset=0;
    double weight=0,error=0;
    char path[1400];
    FILE *f=NULL;
    for (size_t p=0;p<count;p++) {
        const QuadRibbonFitReport *r=part[p].report+level;
        total.input_faces+=r->input_faces; total.retained_faces+=r->retained_faces;
        total.output_vertices+=r->output_vertices; total.output_faces+=r->output_faces;
        total.solve.cells+=r->solve.cells; total.solve.nodes+=r->solve.nodes;
        total.solve.unknowns+=r->solve.unknowns; total.solve.fixed_nodes+=r->solve.fixed_nodes;
        total.iteration+=r->iteration; total.clearance_faces+=r->clearance_faces;
        total.exact_contact_faces+=r->exact_contact_faces;
        total.intersection_refinements+=r->intersection_refinements;
        total.continuous_bound=fmax(total.continuous_bound,r->continuous_bound);
        total.mesh_bound=fmax(total.mesh_bound,r->mesh_bound);
        total.assembly_seconds+=r->assembly_seconds; total.solve_seconds+=r->solve_seconds;
        total.validation_seconds+=r->validation_seconds; p1+=part[p].p1[level];
        if (r->solve.cells>max_cells) max_cells=r->solve.cells;
        if (r->solve.unknowns>max_free) max_free=r->solve.unknowns;
        if (r->solve.node_storage_bytes>max_nodes) max_nodes=r->solve.node_storage_bytes;
        if (part[p].input_bytes>max_input) max_input=part[p].input_bytes;
        weight+=r->solve.observation_weight;
        error+=r->solve.observation_weight*r->rms*r->rms;
    }
    snprintf(path,sizeof path,"%s/level%d/fields.json",writer->directory,level);
    f=fopen(path,"wb"); if (!f) return -1;
    fprintf(f,"{\"schema\":\"supported-quad-field-bank-v1\",\"input_fingerprint\":\"%016llx\",\"cells\":%zu,\"chart_namespace\":\"part-local source vertex; use source_vertices.i32 to resolve\",\"parts\":[",
        (unsigned long long)writer->fingerprint,total.solve.cells);
    for (size_t p=0;p<count;p++) {
        fprintf(f,"%s{\"part\":%zu,\"field\":\"../parts/part%05zu/level%d/field.qfield\",\"source_vertices\":\"../parts/part%05zu/source_vertices.i32\",\"source_faces\":\"../parts/part%05zu/source_faces.i32\",\"input_fingerprint\":\"%016llx\",\"cell_offset\":%zu,\"cells\":%zu,\"nodes\":%zu}",
            p ? "," : "",p,p,level,p,p,(unsigned long long)part[p].fingerprint,offset,
            part[p].report[level].solve.cells,part[p].report[level].solve.nodes);
        offset+=part[p].report[level].solve.cells;
    }
    fprintf(f,"]}\n");
    if (fclose(f)!=0) return -1;
    snprintf(path,sizeof path,"%s/level%d/report.json",writer->directory,level);
    f=fopen(path,"wb"); if (!f) return -1;
    fprintf(f,"{\"schema\":\"supported-global-quad-fit-v3\",\"level\":%d,\"input_fingerprint\":\"%016llx\",\"parts\":%zu,\"input_faces\":%zu,\"retained_faces\":%zu,\"vertices\":%zu,\"faces\":%zu,\"quads\":%zu,\"field_nodes\":%zu,\"free_dofs\":%zu,\"fixed_dofs\":%zu,\"solves\":%zu,\"tolerance\":%.17g,\"continuous_bound\":%.17g,\"mesh_bound\":%.17g,\"rms\":%.17g,\"observation_weight\":%.17g,",
        level,(unsigned long long)writer->fingerprint,count,total.input_faces,total.retained_faces,
        total.output_vertices,total.output_faces,total.solve.cells,total.solve.nodes,total.solve.unknowns,
        total.solve.fixed_nodes,total.iteration,ldexp(QUAD_RIBBON_INITIAL_ERROR,-level),
        total.continuous_bound,total.mesh_bound,weight>0 ? sqrt(error/weight) : 0,weight);
    fprintf(f,"\"orientation_failures\":0,\"unresolved\":0,\"clearance_failures\":0,\"intersection_failures\":0,\"source_intersections\":%zu,\"output_intersections\":%zu,\"intersection_refinements\":%zu,\"clearance_faces\":%zu,\"exact_contact_faces\":%zu,\"clearance_budget_fingerprint\":\"%016llx\",",
        pairs,pairs,total.intersection_refinements,total.clearance_faces,total.exact_contact_faces,
        (unsigned long long)writer->budget_fingerprint);
    fprintf(f,"\"max_active_part_input_bytes\":%zu,\"part_input_byte_cap\":%zu,",
        max_input,(size_t)QUAD_RIBBON_BATCH_INPUT_BYTES);
    fprintf(f,"\"source_triangle_cells\":%zu,\"bilinear_cells\":%zu,\"max_active_part_cells\":%zu,\"max_active_part_free_dofs\":%zu,\"max_active_part_node_storage_bytes\":%zu,\"elapsed_seconds\":%.6f,\"timings\":{\"assembly\":%.6f,\"solve\":%.6f,\"validation\":%.6f},\"qualification\":\"whole-source proper triangle contacts preserved; not a coplanar, continuous embedding or material-correspondence certificate; input inventory and merged audit are not yet out-of-core\"}\n",
        p1,total.solve.cells-p1,max_cells,max_free,max_nodes,elapsed,total.assembly_seconds,
        total.solve_seconds,total.validation_seconds);
    return fclose(f)==0 ? 0 : -1;
}

static int qrf_bank_merge(const QrfWriter *writer,const QrfPartRecord *parts,size_t count,
                           int level,size_t original_pairs,double started)
{
    static const char *names[]={"reference_vertex.i32","reference_face.i32",
        "reference_barycentric.f64","face_cell.i32","reference_triangle.i32",
        "fit_ribbon_material_identity.i32","fit_ribbon_phase.f32","fit_ribbon_support.u8",
        "fit_ribbon_reconstruction_component.i32","fit_ribbon_claimant_chart.i32","fit_ribbon_claimant_island.i32"};
    Arena_T arena=Arena_new(),scratch=Arena_new();
    MeshBinStream_T stream=NULL;
    size_t nv=0,nf=0,nc=0,pairs=0;
    char path[1400];
    int rc=-1;
    for (size_t p=0;p<count;p++) {
        const QuadRibbonFitReport *r=parts[p].report+level;
        if (!(parts[p].completed & (1u<<level)) || r->orientation_failures || r->unresolved_cells ||
            r->clearance_failures || r->intersection_failures) goto done;
        nv+=r->output_vertices; nf+=r->output_faces; nc+=r->solve.cells;
    }
    if (nv>INT32_MAX || nf>INT32_MAX || nc>INT32_MAX) goto done;
    snprintf(path,sizeof path,"%s/level%d",writer->directory,level);
    if (ves_mkdir(path)!=0 && errno!=EEXIST) goto done;
    snprintf(path,sizeof path,"%s/level%d/fit_ribbon_world.vmesh",writer->directory,level);
    if (MeshBin_stream_open(arena,path,nv,nf,1,&stream)!=0) goto done;
    for (size_t p=0;p<count;p++) {
        MeshBinData mesh={0};
        snprintf(path,sizeof path,"%s/parts/part%05zu/level%d/fit_ribbon_world.vmesh",writer->directory,p,level);
        if (MeshBin_read_arena(scratch,path,&mesh)!=0 ||
            mesh.nv!=parts[p].report[level].output_vertices || mesh.nf!=parts[p].report[level].output_faces ||
            MeshBin_stream_append(stream,&mesh)!=0) goto done;
        Arena_free(scratch);
    }
    if (MeshBin_stream_close(&stream,1)!=0) goto done;
    for (size_t j=0;j<sizeof names/sizeof names[0];j++) {
        size_t width=j==2 ? 3*sizeof(double) : j==7 ? 1 : sizeof(int32_t);
        if (qrf_merge_channel(writer->directory,level,parts,count,names[j],width,
                j==3 || j==4,j>=5,j==3)!=0) goto done;
    }
    if (qrf_bank_audit(writer,level,original_pairs,&pairs)!=0 ||
        qrf_bank_manifest(writer,parts,count,level,pairs,ves_clock_sec()-started)!=0) goto done;
    rc=0;
done:
    if (stream) (void)MeshBin_stream_close(&stream,0);
    Arena_dispose(&scratch); Arena_dispose(&arena); return rc;
}

static int qrf_write_part(const QrfWriter *writer,RibbonDomains_T domains,
    const double *budget,const uint8_t *const *sidecars,size_t p,size_t count,QrfPartRecord *record)
{
    Arena_T scratch=Arena_new();
    MeshBinData mesh={0};
    RibbonDomainPart part={0};
    QrfWriter child=*writer;
    double *local_budget=NULL;
    char dir[1200],path[1400];
    int rc=-1;
    if (RibbonDomains_part(domains,p,&part)!=0 ||
        RibbonDomains_extract(scratch,domains,p,QUAD_RIBBON_BATCH_INPUT_BYTES,&mesh)!=0) goto done;
    record->input_bytes=part.input_bytes;
    snprintf(dir,sizeof dir,"%s/parts/part%05zu",writer->directory,p);
    if (ves_mkdir(dir)!=0 && errno!=EEXIST) goto done;
    child.directory=dir; child.source=&mesh; child.original_vertices=writer->source->nv;
    child.source_vertices=part.source_vertices; child.source_faces=part.source_faces;
    child.sidecars=sidecars;
    child.reports=record->report; child.p1_cells=record->p1; child.completed=&record->completed;
    child.fingerprint=qrf_hash(UINT64_C(14695981039346656037),mesh.verts,3*mesh.nv*sizeof(float));
    child.fingerprint=qrf_hash(child.fingerprint,mesh.uv,2*mesh.nv*sizeof(float));
    child.fingerprint=qrf_hash(child.fingerprint,mesh.faces,3*mesh.nf*sizeof(int32_t));
    record->fingerprint=child.fingerprint;
    snprintf(path,sizeof path,"%s/source_vertices.i32",dir);
    if (qrf_write_array(path,part.source_vertices,sizeof(int32_t),mesh.nv)!=0) goto done;
    snprintf(path,sizeof path,"%s/source_faces.i32",dir);
    if (qrf_write_array(path,part.source_faces,sizeof(int32_t),mesh.nf)!=0) goto done;
    local_budget=(double *)ARENA_ALLOC(scratch,mesh.nf*sizeof *local_budget);
    for (size_t f=0;f<mesh.nf;f++) local_budget[f]=budget ? budget[part.source_faces[f]] : QUAD_RIBBON_INITIAL_ERROR;
    child.budget_fingerprint=qrf_hash(UINT64_C(14695981039346656037),local_budget,mesh.nf*sizeof(double));
    fprintf(stderr,"[global ribbon bank] part %zu/%zu: %zu whole components, %zu faces / %zu vertices, %zu input-copy bytes%s\n",
        p+1,count,part.components,mesh.nf,mesh.nv,part.input_bytes,
        part.oversized ? " (oversized, existing guards still apply)" : "");
    rc=qrf_run(&mesh,QUAD_RIBBON_GRID_STEP,local_budget,1,qrf_checkpoint,&child);
done:
    fprintf(stderr,"[global ribbon bank] part %zu/%zu %s; %u completed levels\n",
        p+1,count,rc ? "INCOMPLETE" : "complete",record->completed);
    Arena_dispose(&scratch); return rc;
}

static int qrf_domain_workers(size_t parts,int requested)
{
    int workers=1;
#ifdef _OPENMP
    /* Never move a caller's Hanson exception frame to another thread, or
     * introduce nested work teams. The legacy exception stack is global. */
    if (!Except_stack && !omp_in_parallel()) {
        workers=requested<omp_get_max_threads() ? requested : omp_get_max_threads();
        if (workers>(int)QUAD_RIBBON_BATCH_WORKERS) workers=(int)QUAD_RIBBON_BATCH_WORKERS;
        if (parts<(size_t)workers) workers=(int)parts;
    }
#else
    (void)parts; (void)requested;
#endif
    return workers>0 ? workers : 1;
}

static int qrf_write_domains_workers(const QrfWriter *writer,RibbonDomains_T domains,
    const double *budget,int requested)
{
    Arena_T arena=Arena_new();
    RibbonDomainsReport plan={0};
    QrfPartRecord *records=NULL;
    const uint8_t *sidecars[sizeof qrf_sidecar_suffix/sizeof qrf_sidecar_suffix[0]]={0};
    IntersectionCleanupParams params={0};
    IntersectionCleanupStats stats={0};
    char dir[1200];
    double started=ves_clock_sec();
    int rc=-1,failures=0,workers=1;
    if (requested<1 || RibbonDomains_report(domains,&plan)!=0 || plan.parts>INT_MAX) goto done;
    if (plan.parts<=1) {
        rc=qrf_run(writer->source,QUAD_RIBBON_GRID_STEP,budget,0,qrf_checkpoint,(void *)writer);
        goto done;
    }
    qrf_intersection_params(&params);
    if (IntersectionCleanup_audit_visit_parallel(writer->source->verts,writer->source->nv,
            writer->source->faces,writer->source->nf,NULL,&params,NULL,NULL,NULL,&stats)!=0 || stats.budget_rejected) goto done;
    fprintf(stderr,"[global ribbon bank] whole-source baseline: %zu proper pairs, %.2fs\n",stats.stab_pairs,ves_clock_sec()-started);
    if (qrf_cache_sidecars(arena,writer,sidecars)!=0) goto done;
    records=(QrfPartRecord *)ARENA_CALLOC(arena,plan.parts,sizeof *records);
    snprintf(dir,sizeof dir,"%s/parts",writer->directory);
    if (ves_mkdir(dir)!=0 && errno!=EEXIST) goto done;
    workers=qrf_domain_workers(plan.parts,requested);
    fprintf(stderr,"[global ribbon bank] scheduler: %d workers, private part arenas, serial TAUCS transactions and ordered merge\n",workers);
    if (workers==1) {
        for (size_t p=0;p<plan.parts;p++)
            if (qrf_write_part(writer,domains,budget,sidecars,p,plan.parts,records+p)!=0) goto done;
    } else {
        int p=0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(workers) schedule(dynamic,1) reduction(+:failures)
#endif
        for (p=0;p<(int)plan.parts;p++)
            failures+=qrf_write_part(writer,domains,budget,sidecars,(size_t)p,plan.parts,records+p)!=0;
        if (failures) goto done;
    }
    for (int level=0;level<QUAD_RIBBON_LEVELS;level++)
        if (qrf_bank_merge(writer,records,plan.parts,level,stats.stab_pairs,started)!=0) goto done;
    rc=0;
done:
    Arena_dispose(&arena); return rc;
}

static int qrf_write_domains(const QrfWriter *writer,RibbonDomains_T domains,const double *budget)
{
    return qrf_write_domains_workers(writer,domains,budget,(int)QUAD_RIBBON_BATCH_WORKERS);
}

int QuadRibbonFit_write(const char *input,const char *directory)
{
    Arena_T arena=Arena_new();
    MeshBinData mesh={0};
    RibbonDomains_T domains=NULL;
    double *budget=NULL;
    QrfWriter writer={directory,UINT64_C(14695981039346656037),&mesh,input,0};
    int rc=-1;
    if (!input || !directory || (ves_mkdir(directory)!=0 && errno!=EEXIST) ||
        MeshBin_read_arena(arena,input,&mesh)!=0 || !mesh.uv) goto done;
    writer.fingerprint=qrf_hash(writer.fingerprint,mesh.verts,3*mesh.nv*sizeof(float));
    writer.fingerprint=qrf_hash(writer.fingerprint,mesh.uv,2*mesh.nv*sizeof(float));
    writer.fingerprint=qrf_hash(writer.fingerprint,mesh.faces,3*mesh.nf*sizeof(int32_t));
    if (qrf_domain_plan(arena,&mesh,&domains)!=0) goto done;
    budget=(double *)ARENA_ALLOC(arena,mesh.nf*sizeof *budget);
    for (size_t f=0;f<mesh.nf;f++) budget[f]=QUAD_RIBBON_INITIAL_ERROR;
    if (qrf_partition_clearance(&mesh,domains,budget)!=0) goto done;
    writer.budget_fingerprint=qrf_hash(UINT64_C(14695981039346656037),budget,mesh.nf*sizeof(double));
    rc=qrf_write_domains(&writer,domains,budget);
done:
    Arena_dispose(&arena); return rc;
}

typedef struct {
    MeshBinData mesh;
    RibbonDomains_T domains;
    double *budget;
    char input[1200],directory[1200];
    size_t limited,exact,point_only;
    int present;
} QrfAtlasLayer;

static void qrf_json_string(FILE *f,const char *s)
{
    fputc('"',f);
    for (const unsigned char *p=(const unsigned char *)s;*p;p++) {
        if (*p=='"' || *p=='\\') fputc('\\',f);
        if (*p<32) fprintf(f,"\\u%04x",(unsigned)*p);
        else fputc(*p,f);
    }
    fputc('"',f);
}

int QuadRibbonFit_write_atlas(const char *stem,const char *directory,const AxisWarp *axis)
{
    Arena_T arena=Arena_new();
    const JsonValue *root=NULL,*ribbon=NULL,*inventory=NULL,*inventory_layers=NULL;
    const char *error=NULL;
    QrfAtlasLayer *layers=NULL;
    char path[1400];
    long declared=0;
    size_t surfaces=0,tests=0,refused_members=0;
    double started=ves_clock_sec(),clearance_seconds=0;
    FILE *manifest=NULL;
    int rc=-1,written=0;
    if (!stem || !directory || (axis && !AxisWarp_valid(axis)) ||
        (ves_mkdir(directory)!=0 && errno!=EEXIST)) goto done;
    written=snprintf(path,sizeof path,"%s/atlas.json",directory);
    if (written<0 || (size_t)written>=sizeof path) goto done;
    manifest=fopen(path,"rb");
    if (manifest) { fclose(manifest); manifest=NULL; fprintf(stderr,"[global ribbon atlas] use a fresh output directory; atlas.json already exists\n"); goto done; }
    if (errno!=ENOENT) goto done;
    written=snprintf(path,sizeof path,"%s_stats.json",stem);
    if (written<0 || (size_t)written>=sizeof path) goto done;
    root=Json_parse_file(arena,path,&error);
    ribbon=Json_object_get(root,"ribbon");
    declared=Json_member_long(ribbon,"grid_layers",-1);
    if (declared<1 || declared>100000 || (size_t)declared>SIZE_MAX/sizeof *layers) goto done;
    written=snprintf(path,sizeof path,"%s_layers.json",stem);
    if (written<0 || (size_t)written>=sizeof path) goto done;
    manifest=fopen(path,"rb");
    if (manifest) {
        const char *schema=NULL;
        if (fclose(manifest)!=0) { manifest=NULL; goto done; }
        manifest=NULL;
        inventory=Json_parse_file(arena,path,&error);
        if (Json_as_bool(Json_object_get(inventory,"partial_atlas"),0)) {
            fprintf(stderr,"[global ribbon atlas] selected-peel inventory is a partial atlas; collect every declared peel before joint clearance\n");
            goto done;
        }
        schema=Json_as_string(Json_object_get(inventory,"schema"));
        inventory_layers=Json_object_get(inventory,"layers");
        if (!schema || strcmp(schema,"ribbon-emitted-layers-v1") ||
            !Json_as_bool(Json_object_get(inventory,"complete"),0) ||
            Json_member_long(inventory,"grid_layers",-1)!=declared ||
            Json_array_len(inventory_layers)!=(size_t)declared) goto done;
    } else if (errno!=ENOENT) goto done;
    layers=(QrfAtlasLayer *)ARENA_CALLOC(arena,(size_t)declared,sizeof *layers);
    for (long i=0;i<declared;i++) {
        QrfAtlasLayer *l=layers+i;
        long expected_vertices=-1,expected_faces=-1;
        FILE *check=NULL;
        int planned=0;
        if (inventory) {
            const JsonValue *record=Json_array_get(inventory_layers,(size_t)i);
            long points=Json_member_long(record,"point_only_vertices",-1);
            expected_vertices=Json_member_long(record,"vertices",-1);
            expected_faces=Json_member_long(record,"faces",-1);
            if (Json_member_long(record,"layer",-1)!=i || expected_vertices<0 ||
                expected_faces<0 || points<0 || (!expected_faces && expected_vertices)) goto done;
            l->point_only=(size_t)points;
        }
        if (i==0) written=snprintf(l->input,sizeof l->input,"%s.vmesh",stem);
        else if (i==1) written=snprintf(l->input,sizeof l->input,"%s_extras.vmesh",stem);
        else written=snprintf(l->input,sizeof l->input,"%s_extras_layer%02ld.vmesh",stem,i);
        if (written<0 || (size_t)written>=sizeof l->input) goto done;
        check=fopen(l->input,"rb");
        if (!check) {
            if (i>0 && errno==ENOENT && inventory && expected_faces==0) continue;
            fprintf(stderr,"[global ribbon atlas] missing declared surface %s (requires explicit zero-face emission inventory)\n",l->input);
            goto done;
        }
        if (fclose(check)!=0 || MeshBin_read_arena(arena,l->input,&l->mesh)!=0 || !l->mesh.uv) goto done;
        if (inventory && ((size_t)expected_vertices!=l->mesh.nv ||
                          (size_t)expected_faces!=l->mesh.nf || !expected_faces)) goto done;
        if (axis) for (size_t v=0;v<l->mesh.nv;v++) {
            double metric[3]={l->mesh.verts[3*v],l->mesh.verts[3*v+1],l->mesh.verts[3*v+2]},world[3];
            uint32_t flags=0;
            if (AxisWarp_to_world(axis,metric,world,&flags)!=0) {
                fprintf(stderr,"[global ribbon atlas] inverse coordinate %zu refused (flags %u)\n",v,(unsigned)flags); goto done;
            }
            for (int a=0;a<3;a++) l->mesh.verts[3*v+a]=(float)world[a];
        }
        fprintf(stderr,"[global ribbon atlas] preflight layer %ld: %zu faces / %zu vertices\n",i,l->mesh.nf,l->mesh.nv);
        planned=qrf_domain_plan(arena,&l->mesh,&l->domains);
        if (planned<0) goto done;
        refused_members+=planned>0;
        written=snprintf(l->directory,sizeof l->directory,"%s/layer%02ld",directory,i);
        if (written<0 || (size_t)written>=sizeof l->directory) goto done;
        l->present=1; surfaces++;
    }
    if (refused_members) {
        fprintf(stderr,"[global ribbon atlas] complete input preflight refused %zu of %zu surface members; no clearance or solve attempted\n",
            refused_members,surfaces);
        goto done;
    }
    for (long i=0;i<declared;i++) if (layers[i].present) {
        QrfAtlasLayer *l=layers+i;
        l->budget=(double *)ARENA_ALLOC(arena,l->mesh.nf*sizeof(double));
        for (size_t f=0;f<l->mesh.nf;f++) l->budget[f]=QUAD_RIBBON_INITIAL_ERROR;
    }
    clearance_seconds=ves_clock_sec();
    for (long i=0;i<declared;i++) if (layers[i].present)
        for (long j=i+1;j<declared;j++) if (layers[j].present) {
            QrfAtlasLayer *a=layers+i,*b=layers+j;
            size_t tested=0;
            if (IntersectionCleanup_cross_budget(a->mesh.verts,a->mesh.nv,a->mesh.faces,a->mesh.nf,
                    b->mesh.verts,b->mesh.nv,b->mesh.faces,b->mesh.nf,QUAD_RIBBON_INITIAL_ERROR,
                    a->budget,b->budget,&tested)!=0 || SIZE_MAX-tests<tested) goto done;
            tests+=tested;
            fprintf(stderr,"[global ribbon atlas] clearance layers %ld/%ld: %zu tested pairs, %.2fs\n",i,j,tested,ves_clock_sec()-clearance_seconds);
        }
    for (long i=0;i<declared;i++) if (layers[i].present)
        if (qrf_partition_clearance(&layers[i].mesh,layers[i].domains,layers[i].budget)!=0) goto done;
    clearance_seconds=ves_clock_sec()-clearance_seconds;
    for (long i=0;i<declared;i++) if (layers[i].present) {
        QrfAtlasLayer *l=layers+i;
        QrfWriter writer={l->directory,UINT64_C(14695981039346656037),&l->mesh,l->input,0};
        if (ves_mkdir(l->directory)!=0 && errno!=EEXIST) goto done;
        writer.fingerprint=qrf_hash(writer.fingerprint,l->mesh.verts,3*l->mesh.nv*sizeof(float));
        writer.fingerprint=qrf_hash(writer.fingerprint,l->mesh.uv,2*l->mesh.nv*sizeof(float));
        writer.fingerprint=qrf_hash(writer.fingerprint,l->mesh.faces,3*l->mesh.nf*sizeof(int32_t));
        writer.budget_fingerprint=qrf_hash(UINT64_C(14695981039346656037),l->budget,l->mesh.nf*sizeof(double));
        snprintf(path,sizeof path,"%s/clearance_budget.f64",l->directory);
        if (qrf_write_array(path,l->budget,sizeof(double),l->mesh.nf)!=0) goto done;
        for (size_t f=0;f<l->mesh.nf;f++) { l->limited+=l->budget[f]<QUAD_RIBBON_INITIAL_ERROR; l->exact+=l->budget[f]==0; }
        fprintf(stderr,"[global ribbon atlas] layer %ld: %zu clearance-limited faces, %zu require exact source detail\n",i,l->limited,l->exact);
        if (qrf_write_domains(&writer,l->domains,l->budget)!=0) goto done;
    }
    snprintf(path,sizeof path,"%s/atlas.json",directory);
    manifest=fopen(path,"wb"); if (!manifest) goto done;
    fprintf(manifest,"{\"schema\":\"supported-global-ribbon-atlas-v1\",\"source_stem\":"); qrf_json_string(manifest,stem);
    fprintf(manifest,",\"world_frame\":true,\"axis_unwarped\":%s,\"declared_layers\":%ld,\"surface_layers\":%zu,\"clearance_pairs_tested\":%zu,\"clearance_seconds\":%.6f,\"elapsed_seconds\":%.6f,\"layers\":[",
        axis ? "true" : "false",declared,surfaces,tests,clearance_seconds,ves_clock_sec()-started);
    for (long i=0;i<declared;i++) {
        QrfAtlasLayer *l=layers+i;
        fprintf(manifest,"%s{\"layer\":%ld,\"present\":%s,\"input\":",i ? "," : "",i,l->present ? "true" : "false"); qrf_json_string(manifest,l->input);
        fprintf(manifest,",\"faces\":%zu,\"clearance_limited_faces\":%zu,\"exact_contact_faces\":%zu,\"point_only_vertices\":%zu}",l->mesh.nf,l->limited,l->exact,l->point_only);
    }
    fprintf(manifest,"],\"qualification\":\"preserves original inter-layer contact geometry; existing intersections, duplicate interpretations and within-layer global geometry still require independent audit\"}\n");
    rc=ferror(manifest) ? -1 : 0;
    if (fclose(manifest)!=0) rc=-1;
    manifest=NULL;
done:
    if (manifest) fclose(manifest);
    fprintf(stderr,"[global ribbon atlas] %s after %.2fs (%zu surface members)\n",rc ? "INCOMPLETE" : "completed geometric envelopes",ves_clock_sec()-started,surfaces);
    Arena_dispose(&arena); return rc;
}

typedef struct { int calls; size_t faces; double bound; } QrfTest;
static int qrf_test_checkpoint(void *context,int level,QuadField_T f,
    const QuadRibbonFitOutput *out,const QuadRibbonFitReport *r)
{
    QrfTest *t=(QrfTest *)context;
    if (level!=t->calls || !f || r->mesh_bound>r->tolerance || r->orientation_failures || r->unresolved_cells) return -1;
    for (size_t i=0;i<out->mesh.nv;i++) {
        double sum=out->reference_barycentric[3*i]+out->reference_barycentric[3*i+1]+out->reference_barycentric[3*i+2];
        if (out->reference_face[i]>=0 && fabs(sum-1)>1e-9) return -1;
    }
    t->calls++; t->faces=out->mesh.nf; t->bound=r->mesh_bound;
    return 0;
}

static int qrf_test_equal_file(const char *a,const char *b)
{
    unsigned char left[4096],right[4096];
    FILE *fa=fopen(a,"rb"),*fb=fopen(b,"rb");
    int rc=-1;
    if (!fa || !fb) goto done;
    for (;;) {
        size_t na=fread(left,1,sizeof left,fa),nb=fread(right,1,sizeof right,fb);
        if (na!=nb || memcmp(left,right,na)!=0 || ferror(fa) || ferror(fb)) goto done;
        if (na<sizeof left) { rc=0; break; }
    }
done:
    if (fa && fclose(fa)!=0) rc=-1;
    if (fb && fclose(fb)!=0) rc=-1;
    return rc;
}

static int qrf_test_worker_stacks(void)
{
    int failures=0;
#if defined(_WIN32) && defined(_OPENMP)
    int workers=qrf_domain_workers(2,2);
#pragma omp parallel num_threads(workers) reduction(+:failures)
    {
        int marker=0;
        MEMORY_BASIC_INFORMATION memory={0};
        size_t reserved=0;
        if (VirtualQuery(&marker,&memory,sizeof memory)!=sizeof memory) failures++;
        else {
            reserved=(size_t)((uintptr_t)memory.BaseAddress+memory.RegionSize-
                              (uintptr_t)memory.AllocationBase);
            failures+=reserved<((size_t)512*1024*1024);
        }
        fprintf(stderr,"[selftest] qrf worker %d stack reservation: %zu bytes\n",
            omp_get_thread_num(),reserved);
    }
#endif
    return failures ? -1 : 0;
}

static int qrf_test_bank(const MeshBinData *base)
{
    Arena_T arena=Arena_new();
    MeshBinData mesh={0};
    RibbonDomains_T domains=NULL;
    RibbonDomainsReport plan={0};
    char dir[1100],parallel[1200],input[1200],path[1400],other_path[1400];
    QrfWriter writer={0};
    double *budget=NULL;
    int rc=-1;
    {
        Except_Frame frame={0},*saved=Except_stack;
        int workers=0;
        Except_stack=&frame;
        workers=qrf_domain_workers(8,4);
        Except_stack=saved;
        if (workers!=1 || qrf_domain_workers(1,4)!=1 || qrf_test_worker_stacks()!=0) goto done;
    }
    mesh.nv=2*base->nv; mesh.nf=2*base->nf;
    mesh.verts=(float *)ARENA_ALLOC(arena,3*mesh.nv*sizeof(float));
    mesh.uv=(float *)ARENA_ALLOC(arena,2*mesh.nv*sizeof(float));
    mesh.faces=(int32_t *)ARENA_ALLOC(arena,3*mesh.nf*sizeof(int32_t));
    for (size_t i=0;i<mesh.nv;i++) {
        size_t original=i%base->nv;
        int other=i>=base->nv;
        for (int k=0;k<3;k++) mesh.verts[3*i+k]=base->verts[3*original+k];
        for (int k=0;k<2;k++) mesh.uv[2*i+k]=2*base->uv[2*original+k];
        /* Distinct indexed components with identical UVs and a genuine
         * crossing. The work plan must not join them or average their field. */
        if (other) {
            mesh.verts[3*i]+=.375f; mesh.verts[3*i+1]+=.25f;
            mesh.verts[3*i+2]+=.15f*(base->uv[2*original]-8.5f);
        }
    }
    for (size_t f=0;f<mesh.nf;f++) for (int k=0;k<3;k++)
        mesh.faces[3*f+k]=base->faces[3*(f%base->nf)+k]+(f>=base->nf ? (int32_t)base->nv : 0);
    if (RibbonDomains_new(arena,&mesh,base->nf,base->nv,&domains)!=0 ||
        RibbonDomains_report(domains,&plan)!=0 || plan.parts!=2 || plan.components!=2) goto done;
    {
        RibbonDomainPart part={0};
        RibbonDomains_T borrowed=NULL;
        if (RibbonDomains_part(domains,0,&part)!=0 || !part.input_bytes ||
            qrf_domain_capacity(domains,part.input_bytes-1)!=1 ||
            qrf_domain_capacity(domains,part.input_bytes)!=0 ||
            qrf_domain_capacity(NULL,part.input_bytes)!=-1 ||
            RibbonDomains_new(arena,&mesh,mesh.nf,mesh.nv,&borrowed)!=0 ||
            qrf_domain_capacity(borrowed,0)!=0) goto done;
    }
    budget=(double *)ARENA_ALLOC(arena,mesh.nf*sizeof(double));
    for (size_t f=0;f<mesh.nf;f++) budget[f]=QUAD_RIBBON_INITIAL_ERROR;
    if (qrf_partition_clearance(&mesh,domains,budget)!=0) goto done;
    /* Retain this small, uniquely named native regression for independent
     * readers. Never replace an older experiment or remove user artifacts. */
    snprintf(dir,sizeof dir,"output/qrf_bank_selftest_%d_%llu",ves_getpid(),
        (unsigned long long)(ves_clock_sec()*1000000));
    if (ves_mkdir(dir)!=0) goto done;
    snprintf(input,sizeof input,"%s/source.vmesh",dir);
    if (MeshBin_write(input,mesh.verts,mesh.nv,mesh.faces,mesh.nf,mesh.uv)!=0) goto done;
    /* Exercise the shared read-only provenance cache as well as geometry.
     * Missing optional sidecars are not files to compare; this fixture now
     * supplies all six with distinct nonzero source-indexed payloads. */
    for (size_t j=0;j<sizeof qrf_sidecar_suffix/sizeof qrf_sidecar_suffix[0];j++) {
        size_t width=j==2 ? 1 : 4;
        void *values=ARENA_ALLOC(arena,mesh.nv*width);
        for (size_t v=0;v<mesh.nv;v++) {
            if (j==2) ((uint8_t *)values)[v]=(uint8_t)(v%251+1);
            else if (j==1) ((float *)values)[v]=(float)(v%97)*.25f;
            else ((int32_t *)values)[v]=(int32_t)(1000*j+v+1);
        }
        snprintf(path,sizeof path,"%s/source%s",dir,qrf_sidecar_suffix[j]);
        if (qrf_write_array(path,values,width,mesh.nv)!=0) goto done;
    }
    writer.directory=dir; writer.source=&mesh; writer.input=input;
    writer.fingerprint=qrf_hash(UINT64_C(14695981039346656037),mesh.verts,3*mesh.nv*sizeof(float));
    writer.fingerprint=qrf_hash(writer.fingerprint,mesh.uv,2*mesh.nv*sizeof(float));
    writer.fingerprint=qrf_hash(writer.fingerprint,mesh.faces,3*mesh.nf*sizeof(int32_t));
    if (qrf_write_domains_workers(&writer,domains,budget,1)!=0) goto done;
    snprintf(parallel,sizeof parallel,"%s/parallel",dir);
    if (ves_mkdir(parallel)!=0) goto done;
    writer.directory=parallel;
    if (qrf_write_domains_workers(&writer,domains,budget,2)!=0) goto done;
    for (int level=0;level<QUAD_RIBBON_LEVELS;level++) {
        static const char *payload[]={"fit_ribbon_world.vmesh","face_cell.i32",
            "reference_vertex.i32","reference_face.i32","reference_barycentric.f64",
            "reference_triangle.i32","fit_ribbon_material_identity.i32","fit_ribbon_phase.f32",
            "fit_ribbon_support.u8","fit_ribbon_reconstruction_component.i32",
            "fit_ribbon_claimant_chart.i32","fit_ribbon_claimant_island.i32"};
        const JsonValue *report=NULL;
        const char *error=NULL;
        snprintf(path,sizeof path,"%s/level%d/report.json",dir,level);
        report=Json_parse_file(arena,path,&error);
        if (!report || Json_member_long(report,"parts",-1)!=2 ||
            Json_member_long(report,"source_intersections",-1)<1 ||
            Json_member_long(report,"source_intersections",-1)!=Json_member_long(report,"output_intersections",-2)) goto done;
        for (int part=-1;part<2;part++) {
            char relative[128];
            if (part<0) snprintf(relative,sizeof relative,"level%d",level);
            else snprintf(relative,sizeof relative,"parts/part%05d/level%d",part,level);
            for (size_t j=0;j<=sizeof payload/sizeof payload[0];j++) {
                const char *name=j<sizeof payload/sizeof payload[0] ? payload[j] :
                    part<0 ? "fields.json" : "field.qfield";
                snprintf(path,sizeof path,"%s/%s/%s",dir,relative,name);
                snprintf(other_path,sizeof other_path,"%s/%s/%s",parallel,relative,name);
                if (qrf_test_equal_file(path,other_path)!=0) {
                    fprintf(stderr,"[selftest] serial/parallel mismatch: %s/%s\n",relative,name);
                    goto done;
                }
            }
        }
    }
    rc=0;
done:
    fprintf(stderr,"[selftest] bounded field bank with cross-part contacts and exact serial/parallel checkpoints %s\n",rc ? "FAIL" : "PASS");
    Arena_dispose(&arena); return rc;
}
int QuadRibbonFit_selftest(void)
{
    enum { N=16,NV=(N+1)*(N+1),NF=2*N*N };
    Arena_T arena=Arena_new();
    MeshBinData m={0}; QrfTest test={0};
    size_t nf=0;
    int fail=IntersectionCleanup_selftest()!=0;
    fail+=RibbonDomains_selftest()!=0;
    m.nv=NV; m.nf=NF;
    m.verts=(float *)ARENA_ALLOC(arena,3*NV*sizeof(float));
    m.uv=(float *)ARENA_ALLOC(arena,2*NV*sizeof(float));
    m.faces=(int32_t *)ARENA_ALLOC(arena,3*NF*sizeof(int32_t));
    for (int y=0;y<=N;y++) for (int x=0;x<=N;x++) {
        size_t i=(size_t)y*(N+1)+x;
        m.uv[2*i]=(float)x; m.uv[2*i+1]=(float)y;
        m.verts[3*i]=(float)y; m.verts[3*i+1]=(float)x;
        m.verts[3*i+2]=(float)(.04*x*x+.03*y*y);
    }
    for (int y=0;y<N;y++) for (int x=0;x<N;x++) {
        int a=y*(N+1)+x,b=a+1,c=a+N+2,d=a+N+1;
        int32_t tri[6]={a,b,c,a,c,d};
        memcpy(m.faces+3*nf,tri,sizeof tri); nf+=2;
    }
    fail+=QuadRibbonFit_run(&m,1,qrf_test_checkpoint,&test)!=0;
    fail+=test.calls!=QUAD_RIBBON_LEVELS || test.bound>.5 || test.faces>=NF;
    fail+=qrf_test_bank(&m)!=0;
    fprintf(stderr,"[selftest] supported global ribbon fit %s (%d failures)\n",fail ? "FAIL" : "PASS",fail);
    Arena_dispose(&arena); return fail ? -1 : 0;
}
