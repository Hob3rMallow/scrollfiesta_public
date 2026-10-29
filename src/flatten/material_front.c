#include "material_front.h"
#include "../common/pipeline_constants.h"

#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct { int32_t a,b,third; size_t face; int edge; } MfEdge;
typedef struct { float lo[3],hi[3],center[3]; size_t front; } MfBox;
typedef struct { float lo[3],hi[3]; size_t first,last,left,right; } MfNode;
typedef struct { int32_t chart,vertex; size_t front; int end; } MfEndpoint;
typedef struct { uint64_t edge; size_t front; } MfVisit;

static int mf_edge_compare(const void *lhs,const void *rhs)
{
    const MfEdge *a=lhs,*b=rhs;
    int32_t al=a->a<a->b ? a->a : a->b, ah=a->a<a->b ? a->b : a->a;
    int32_t bl=b->a<b->b ? b->a : b->b, bh=b->a<b->b ? b->b : b->a;
    if (al!=bl) return al<bl ? -1 : 1;
    if (ah!=bh) return ah<bh ? -1 : 1;
    return a->face<b->face ? -1 : a->face>b->face;
}

static double mf_dot(const double a[3],const double b[3])
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

static void mf_cross(const double a[3],const double b[3],double c[3])
{ c[0]=a[1]*b[2]-a[2]*b[1]; c[1]=a[2]*b[0]-a[0]*b[2]; c[2]=a[0]*b[1]-a[1]*b[0]; }

static int mf_endpoint_compare(const void *lhs,const void *rhs)
{
    const MfEndpoint *a=lhs,*b=rhs;
    if (a->chart!=b->chart) return a->chart<b->chart ? -1 : 1;
    if (a->vertex!=b->vertex) return a->vertex<b->vertex ? -1 : 1;
    if (a->end!=b->end) return a->end<b->end ? -1 : 1;
    return a->front<b->front ? -1 : a->front>b->front;
}

static int mf_visit_compare(const void *lhs,const void *rhs)
{
    const MfVisit *a=lhs,*b=rhs;
    return a->edge<b->edge ? -1 : a->edge>b->edge;
}

int MaterialFront_trace(Arena_T arena,MaterialFront *fronts,size_t n)
{
    if (!arena || (n && !fronts) || n>SIZE_MAX/(2*sizeof(MfEndpoint))) return -1;
    if (!n) return 0;
    Arena_Mark mark=Arena_save(arena);
    MfEndpoint *ends=ARENA_ALLOC(arena,2*n*sizeof *ends);
    MfVisit *order=ARENA_ALLOC(arena,n*sizeof *order);
    size_t *next=ARENA_ALLOC(arena,n*sizeof *next);
    size_t *previous=ARENA_ALLOC(arena,n*sizeof *previous);
    uint8_t *seen=ARENA_CALLOC(arena,n,sizeof *seen);
    uint64_t prefix=fronts[0].source_edge>>32, regions=0;
    int rc=-1;
    for (size_t i=0;i<n;i++) {
        const MaterialFront *f=fronts+i;
        if (f->chart<0 || f->source_vertices[0]<0 || f->source_vertices[1]<0 ||
            f->source_vertices[0]==f->source_vertices[1] || f->source_edge>>32!=prefix) goto done;
        ends[2*i]=(MfEndpoint){f->chart,f->source_vertices[0],i,0};
        ends[2*i+1]=(MfEndpoint){f->chart,f->source_vertices[1],i,1};
        order[i]=(MfVisit){f->source_edge,i}; next[i]=previous[i]=SIZE_MAX;
    }
    qsort(ends,2*n,sizeof *ends,mf_endpoint_compare);
    qsort(order,n,sizeof *order,mf_visit_compare);
    for (size_t i=1;i<n;i++) if (order[i].edge==order[i-1].edge) goto done;
    for (size_t i=0;i<2*n;) {
        size_t j=i+1;
        while (j<2*n && ends[j].chart==ends[i].chart && ends[j].vertex==ends[i].vertex) j++;
        if (j-i>2 || (j-i==2 && ends[i].end==ends[i+1].end)) goto done;
        if (j-i==2) {
            next[ends[i+1].front]=ends[i].front;
            previous[ends[i].front]=ends[i+1].front;
        }
        i=j;
    }
    /* Trace open chains first, then closed loops. Source identities, not input
     * enumeration, choose the origin; no spatial neighbor can change a curve. */
    for (int pass=0;pass<2;pass++) for (size_t i=0;i<n;i++) {
        size_t seed=order[i].front, at=seed;
        double arc=0;
        if (seen[seed] || (!pass && previous[seed]!=SIZE_MAX)) continue;
        do {
            MaterialFront *f=fronts+at;
            double length=0;
            if (seen[at]) goto done;
            for (int k=0;k<3;k++) {
                double delta=(double)f->xyz[k+3]-f->xyz[k];
                length+=delta*delta;
            }
            length=sqrt(length);
            if (!isfinite(length) || !(length>0)) goto done;
            f->curve=fronts[seed].source_edge; f->arc=arc;
            f->region_base=(prefix<<32)|regions;
            arc+=length; seen[at]=1; at=next[at];
        } while (at!=SIZE_MAX && at!=seed);
        double blocks=ceil(arc/MATERIAL_FRONT_CORRELATION_LENGTH);
        if (!isfinite(blocks) || blocks<1 || blocks>(double)(UINT32_MAX-regions)) goto done;
        regions+=(uint64_t)blocks;
    }
    rc=0;
done:
    Arena_restore(arena,mark);
    return rc;
}

int MaterialFront_extract(Arena_T arena,const MeshBinData *mesh,
                          const ScrollSource *source,uint32_t source_prefix,
                          int32_t chart_base,MaterialFront **out,size_t *count)
{
    float *direction=NULL;
    size_t n=0;
    if (!out || !count) return -1;
    *out=NULL; *count=0;
    if (!arena || !mesh || !source || chart_base<0 || mesh->nf>UINT32_MAX/3 ||
        mesh->nf>SIZE_MAX/(3*sizeof(MfEdge)) ||
        source->ncomponents>(size_t)(INT32_MAX-chart_base) ||
        (mesh->nv && (!source->component || !source->q))) return -1;
    if (ScrollSource_boundary_directions(arena,mesh,&direction,NULL)!=0) return -1;
    for (size_t i=0;i<mesh->nv;i++)
        if (source->component[i]<0 || (size_t)source->component[i]>=source->ncomponents ||
            !isfinite(source->q[i])) return -1;
    MfEdge *edges=ARENA_ALLOC(arena,3*mesh->nf*sizeof *edges);
    for (size_t f=0;f<mesh->nf;f++) for (int k=0;k<3;k++) {
        edges[3*f+(size_t)k]=(MfEdge){mesh->faces[3*f+k],mesh->faces[3*f+(k+1)%3],
                                      mesh->faces[3*f+(k+2)%3],f,k};
    }
    qsort(edges,3*mesh->nf,sizeof *edges,mf_edge_compare);
    for (size_t i=0;i<3*mesh->nf;) {
        size_t j=i+1;
        while (j<3*mesh->nf && ((edges[j].a==edges[i].a && edges[j].b==edges[i].b) ||
               (edges[j].b==edges[i].a && edges[j].a==edges[i].b))) j++;
        if (j==i+1) n++;
        i=j;
    }
    MaterialFront *fronts=ARENA_CALLOC(arena,n,sizeof *fronts);
    n=0;
    for (size_t i=0;i<3*mesh->nf;) {
        size_t j=i+1;
        const MfEdge *e=edges+i;
        while (j<3*mesh->nf && ((edges[j].a==e->a && edges[j].b==e->b) ||
               (edges[j].b==e->a && edges[j].a==e->b))) j++;
        if (j==i+1) {
            double u[3]={0},v[3]={0},normal[3]={0},outward[3]={0};
            const float *a=mesh->verts+3*(size_t)e->a,*b=mesh->verts+3*(size_t)e->b;
            const float *c=mesh->verts+3*(size_t)e->third;
            const float *da=direction+3*(size_t)e->a,*db=direction+3*(size_t)e->b;
            if ((da[0]==0 && da[1]==0 && da[2]==0) ||
                (db[0]==0 && db[1]==0 && db[2]==0)) { i=j; continue; }
            for (int k=0;k<3;k++) { u[k]=(double)b[k]-a[k]; v[k]=(double)c[k]-a[k]; }
            mf_cross(u,v,normal);
            double nl=sqrt(mf_dot(normal,normal)), el=sqrt(mf_dot(u,u));
            if (!(nl>0) || !(el>0)) { i=j; continue; }
            for (int k=0;k<3;k++) normal[k]/=nl;
            mf_cross(u,normal,outward);
            MaterialFront *f=fronts+n++;
            memcpy(f->xyz,a,3*sizeof(float)); memcpy(f->xyz+3,b,3*sizeof(float));
            for (int k=0;k<3;k++) { f->normal[k]=(float)normal[k]; f->outward[k]=(float)(outward[k]/el); }
            f->phase[0]=source->q[e->a]; f->phase[1]=source->q[e->b];
            f->chart=chart_base+source->component[e->a];
            f->face=((uint64_t)source_prefix<<32)|(uint64_t)e->face;
            f->source_edge=((uint64_t)source_prefix<<32)|(uint64_t)(3*e->face+(size_t)e->edge);
            f->source_vertices[0]=e->a; f->source_vertices[1]=e->b;
        }
        i=j;
    }
    if (MaterialFront_trace(arena,fronts,n)!=0) return -1;
    *out=fronts; *count=n;
    return 0;
}

static int mf_box_compare(const MfBox *a,const MfBox *b,int axis)
{
    if (a->center[axis]!=b->center[axis]) return a->center[axis]<b->center[axis] ? -1 : 1;
    return a->front<b->front ? -1 : a->front>b->front;
}
static int mf_box0(const void *a,const void *b) { return mf_box_compare(a,b,0); }
static int mf_box1(const void *a,const void *b) { return mf_box_compare(a,b,1); }
static int mf_box2(const void *a,const void *b) { return mf_box_compare(a,b,2); }

static size_t mf_tree(MfBox *box,size_t first,size_t last,MfNode *nodes,size_t *nn)
{
    size_t id=(*nn)++;
    MfNode *node=nodes+id;
    node->first=first; node->last=last; node->left=node->right=SIZE_MAX;
    memcpy(node->lo,box[first].lo,sizeof node->lo); memcpy(node->hi,box[first].hi,sizeof node->hi);
    for (size_t i=first+1;i<last;i++) for (int k=0;k<3;k++) {
        node->lo[k]=fminf(node->lo[k],box[i].lo[k]); node->hi[k]=fmaxf(node->hi[k],box[i].hi[k]);
    }
    if (last-first>8) {
        int axis=0;
        for (int k=1;k<3;k++) if (node->hi[k]-node->lo[k]>node->hi[axis]-node->lo[axis]) axis=k;
        qsort(box+first,last-first,sizeof *box,axis==0 ? mf_box0 : axis==1 ? mf_box1 : mf_box2);
        size_t mid=first+(last-first)/2;
        node->left=mf_tree(box,first,mid,nodes,nn); node->right=mf_tree(box,mid,last,nodes,nn);
    }
    return id;
}

static int mf_overlap_box(const MfBox *a,const MfNode *b)
{
    for (int k=0;k<3;k++)
        if ((double)a->hi[k]+MATERIAL_FRONT_GAP<b->lo[k] ||
            (double)b->hi[k]+MATERIAL_FRONT_GAP<a->lo[k]) return 0;
    return 1;
}

/* Integrate overlap along A, while testing the corresponding B interval.
 * Opposite oriented source halfedges represent a consistently oriented seam.
 * The normal sign itself may differ between independently oriented charts. */
static int mf_interval(const MaterialFront *a,const MaterialFront *b,
                        MaterialEvidenceInterval *answer,MaterialFrontReport *report)
{
    double u[3]={0},v[3]={0},delta[3]={0},al=0,bl=0,along=0;
    for (int k=0;k<3;k++) { u[k]=(double)a->xyz[3+k]-a->xyz[k]; v[k]=(double)b->xyz[3+k]-b->xyz[k]; }
    al=sqrt(mf_dot(u,u)); bl=sqrt(mf_dot(v,v));
    if (!(al>0) || !(bl>0)) return 0;
    for (int k=0;k<3;k++) { u[k]/=al; v[k]/=bl; delta[k]=(double)b->xyz[k]-a->xyz[k]; }
    along=mf_dot(u,v);
    if (fabs(along)<MATERIAL_FRONT_TANGENT_COS) return 0;
    double normal_agreement=0;
    for (int k=0;k<3;k++) normal_agreement+=(double)a->normal[k]*b->normal[k];
    if (fabs(normal_agreement)<MATERIAL_FRONT_TANGENT_COS) { report->contacts++; return 0; }
    double p=mf_dot(delta,u), q=p+bl*along;
    double lo=fmax(0,fmin(p,q)), hi=fmin(al,fmax(p,q));
    if (!(hi>lo)) return 0;
    /* Local differences of float source coordinates are exact in double.
     * Norms, normalization and projected endpoints are not. An apparent
     * overlap below their roundoff envelope is point contact, not physical
     * support; it can even round to ZERO length on the paired source edge.
     * This numerical guard is not a voxel-scale material threshold. */
    double length_error=64*DBL_EPSILON*(al+bl+fabs(delta[0])+fabs(delta[1])+fabs(delta[2]));
    if (hi-lo<=length_error || (hi-p)/along==(lo-p)/along) {
        report->degenerate_intervals++;
        return 0;
    }
    double residual_max=0, gap_max=0;
    int32_t target=0;
    for (int sample=0;sample<3;sample++) {
        double s=lo+(hi-lo)*.5*sample;
        double t=(s-p)/along, gap[3]={0},toward_a=0,toward_b=0;
        for (int k=0;k<3;k++) {
            gap[k]=(double)b->xyz[k]+t*v[k]-a->xyz[k]-s*u[k];
            toward_a+=gap[k]*a->outward[k]; toward_b-=gap[k]*b->outward[k];
        }
        double dist=sqrt(mf_dot(gap,gap));
        if (dist>MATERIAL_FRONT_GAP) return 0;
        double opposed=0;
        for (int k=0;k<3;k++) opposed+=(double)a->outward[k]*b->outward[k];
        if (opposed>=0 || (dist>1e-8 && (toward_a<=0 || toward_b<=0))) {
            report->contacts++; return 0;
        }
        double qa=a->phase[0]+(s/al)*(a->phase[1]-a->phase[0]);
        double qb=b->phase[0]+(t/bl)*(b->phase[1]-b->phase[0]);
        double d=qa-qb, rounded=round(d);
        if (rounded<=INT32_MIN || rounded>INT32_MAX ||
            (sample && (int32_t)rounded!=target)) { report->ambiguous_phase++; return 0; }
        target=(int32_t)rounded;
        residual_max=fmax(residual_max,fabs(d-rounded)); gap_max=fmax(gap_max,dist);
    }
    if (residual_max>MATERIAL_FRONT_PHASE_TOL) { report->ambiguous_phase++; return 0; }
    memset(answer,0,sizeof *answer);
    answer->a=a->chart; answer->b=b->chart; answer->target=target;
    answer->kind=WINDING_MRF_EQUAL; answer->begin=lo; answer->end=hi;
    answer->begin_b=(lo-p)/along; answer->end_b=(hi-p)/along;
    answer->source_face_a=a->face; answer->source_face_b=b->face;
    answer->source_edge_a=a->source_edge; answer->source_edge_b=b->source_edge;
    answer->quality=exp(-.5*(gap_max*gap_max)/(MATERIAL_FRONT_GAP*MATERIAL_FRONT_GAP)) *
                    exp(-.5*(residual_max*residual_max)/(MATERIAL_FRONT_PHASE_TOL*MATERIAL_FRONT_PHASE_TOL));
    return 1;
}

/* Merge both sources' block-boundary events in overlap parameter t. The paired
 * B interval can run backwards. Midpoint classification avoids attributing an
 * interval endpoint to a neighboring block. Splitting carries the original
 * source face/halfedge and adds no observations of its own. */
static int mf_split_interval(const MaterialFront *a,const MaterialFront *b,
                              const MaterialEvidenceInterval *v,
                              MaterialEvidenceInterval *result,size_t capacity,
                              size_t *count)
{
    double a0=a->arc+v->begin, da=v->end-v->begin;
    double b0=b->arc+v->begin_b, db=v->end_b-v->begin_b;
    double width=MATERIAL_FRONT_CORRELATION_LENGTH, previous=0;
    double na=floor(a0/width)+1;
    double nb=db>0 ? floor(b0/width)+1 : ceil(b0/width)-1;
    if (!isfinite(a0) || !isfinite(b0) || da<=0 || db==0) {
        fprintf(stderr,"[material fronts] invalid mapped interval a0=%.17g b0=%.17g da=%.17g db=%.17g\n",a0,b0,da,db);
        return -1;
    }
    for (size_t step=0;previous<1;step++) {
        if (step>MATERIAL_FRONT_MAX_BLOCK_SPLITS) {
            fprintf(stderr,"[material fronts] block subdivision limit at t=%.17g a0=%.17g b0=%.17g da=%.17g db=%.17g\n",previous,a0,b0,da,db);
            return -1;
        }
        double ta=(na*width-a0)/da, tb=(nb*width-b0)/db;
        double end=fmin(1,fmin(ta,tb));
        if (end>previous) {
            double mid=.5*(previous+end);
            double ba=floor((a0+mid*da)/width), bb=floor((b0+mid*db)/width);
            if (ba<0 || bb<0 || ba>(double)(UINT32_MAX-(uint32_t)a->region_base) ||
                bb>(double)(UINT32_MAX-(uint32_t)b->region_base) || *count>=capacity) {
                fprintf(stderr,"[material fronts] block refusal count=%zu capacity=%zu a_block=%.17g b_block=%.17g t=[%.17g,%.17g] a0=%.17g b0=%.17g da=%.17g db=%.17g\n",
                    *count,capacity,ba,bb,previous,end,a0,b0,da,db);
                return -1;
            }
            if (result) {
                MaterialEvidenceInterval *part=result+*count;
                *part=*v; part->unit=(uint64_t)*count;
                part->begin=v->begin+previous*da; part->end=v->begin+end*da;
                part->begin_b=v->begin_b+previous*db; part->end_b=v->begin_b+end*db;
                part->region_a=a->region_base+(uint64_t)ba;
                part->region_b=b->region_base+(uint64_t)bb;
            }
            (*count)++;
            previous=end;
        }
        if (ta<=tb) na++;
        if (tb<=ta) nb+=db>0 ? 1 : -1;
    }
    return 0;
}

static int mf_interval_bytes(size_t n,size_t cap,size_t *bytes)
{
    *bytes=0;
    if (n>SIZE_MAX/sizeof(MaterialEvidenceInterval)) return -1;
    *bytes=n*sizeof(MaterialEvidenceInterval);
    return *bytes>cap ? -1 : 0;
}

int MaterialFront_collect(Arena_T arena,const MaterialFront *fronts,size_t n,
                          MaterialEvidenceInterval **out,size_t *count,
                          MaterialFrontReport *report)
{
    MaterialFrontReport local={0};
    if (!out || !count) return -1;
    *out=NULL; *count=0;
    if (!report) report=&local;
    memset(report,0,sizeof *report); report->fronts=n;
    if (!arena || (n && !fronts) || n>SIZE_MAX/(2*sizeof(MfNode))) return -1;
    if (!n) return 0;
    Arena_T scratch=Arena_new();
    int rc=-1;
    MfBox *boxes=ARENA_ALLOC(scratch,n*sizeof *boxes);
    MfNode *nodes=ARENA_ALLOC(scratch,2*n*sizeof *nodes);
    size_t nn=0,total=0;
    for (size_t i=0;i<n;i++) {
        if (fronts[i].chart<0 || !isfinite(fronts[i].phase[0]) || !isfinite(fronts[i].phase[1])) goto done;
        boxes[i].front=i;
        for (int k=0;k<3;k++) {
            float a=fronts[i].xyz[k], b=fronts[i].xyz[3+k];
            if (!isfinite(a) || !isfinite(b) || !isfinite(fronts[i].outward[k]) ||
                !isfinite(fronts[i].normal[k])) goto done;
            boxes[i].lo[k]=fminf(a,b); boxes[i].hi[k]=fmaxf(a,b);
            boxes[i].center[k]=(float)(.5*((double)a+b));
        }
    }
    mf_tree(boxes,0,n,nodes,&nn);
    MaterialEvidenceInterval *result=NULL;
    /* Count then emit: bounded allocation, no silent prefix on overflow. */
    for (int pass=0;pass<2;pass++) {
        size_t at=0;
        MaterialFrontReport stats={0}; stats.fronts=n;
        for (size_t i=0;i<n;i++) {
            size_t stack[64]={0},top=1;
            MfBox query={0};
            for (int k=0;k<3;k++) {
                query.lo[k]=fminf(fronts[i].xyz[k],fronts[i].xyz[k+3]);
                query.hi[k]=fmaxf(fronts[i].xyz[k],fronts[i].xyz[k+3]);
            }
            while (top) {
                const MfNode *node=nodes+stack[--top];
                if (!mf_overlap_box(&query,node)) continue;
                if (node->left!=SIZE_MAX) {
                    if (top>62) goto done;
                    stack[top++]=node->right; stack[top++]=node->left;
                    continue;
                }
                for (size_t j=node->first;j<node->last;j++) {
                    size_t other=boxes[j].front;
                    MaterialEvidenceInterval interval={0};
                    if (other<=i || fronts[i].chart==fronts[other].chart) continue;
                    stats.candidates++;
                    if (!mf_interval(fronts+i,fronts+other,&interval,&stats)) continue;
                    if (mf_split_interval(fronts+i,fronts+other,&interval,pass ? result : NULL,
                                          pass ? total : SIZE_MAX/sizeof *result,&at)!=0) {
                        fprintf(stderr,"[material fronts] interval failed pass=%d fronts=%zu/%zu edges=%llu/%llu arc=%.17g/%.17g begin=%.17g/%.17g end=%.17g/%.17g; no evidence prefix published\n",
                            pass,i,other,(unsigned long long)fronts[i].source_edge,
                            (unsigned long long)fronts[other].source_edge,fronts[i].arc,fronts[other].arc,
                            interval.begin,interval.begin_b,interval.end,interval.end_b);
                        goto done;
                    }
                    stats.supported_length+=interval.end-interval.begin;
                }
            }
        }
        stats.supported_intervals=at;
        if (!pass) {
            size_t bytes=0;
            total=at;
            int fits=mf_interval_bytes(total,(size_t)MATERIAL_FRONT_INTERVAL_BYTES,&bytes)==0;
            fprintf(stderr,"[material fronts] complete count: %zu fronts, %zu intervals, %zu interval bytes; cap=%zu bytes\n",
                    n,total,bytes,(size_t)MATERIAL_FRONT_INTERVAL_BYTES);
            if (!fits) {
                *report=stats;
                fprintf(stderr,"[material fronts] complete count exceeds resident interval cap; no evidence prefix published\n");
                goto done;
            }
            result=ARENA_ALLOC(arena,bytes);
        }
        else { if (at!=total) goto done; *report=stats; }
    }
    *out=result; *count=total;
    rc=0;
done:
    Arena_dispose(&scratch);
    return rc;
}

int MaterialFront_selftest(void)
{
    Arena_T arena=Arena_new();
    MaterialFront f[2]={0};
    MaterialEvidenceInterval *v=NULL;
    MaterialFrontReport report={0};
    size_t n=0;
    int fail=0;
    {
        size_t bytes=0, cap=(size_t)MATERIAL_FRONT_INTERVAL_BYTES;
        size_t boundary=cap/sizeof(MaterialEvidenceInterval);
        fail+=mf_interval_bytes(boundary,cap,&bytes)!=0 || bytes>cap;
        fail+=mf_interval_bytes(boundary+1,cap,&bytes)==0 || bytes<=cap;
        fail+=mf_interval_bytes(SIZE_MAX,cap,&bytes)==0 || bytes!=0;
        fail+=mf_interval_bytes(8076891,cap,&bytes)!=0 || bytes!=904611792;
    }
    f[0].xyz[3]=8; f[1].xyz[0]=8; f[1].xyz[1]=1; f[1].xyz[4]=1;
    f[0].outward[1]=1; f[1].outward[1]=-1;
    f[0].normal[2]=f[1].normal[2]=1;
    f[1].chart=1; f[1].face=1; f[1].phase[0]=f[1].phase[1]=-1;
    f[0].source_vertices[1]=1;
    f[1].source_vertices[0]=2; f[1].source_vertices[1]=3; f[1].source_edge=3;
    fail+=MaterialFront_trace(arena,f,2)!=0;
    fail+=MaterialFront_collect(arena,f,2,&v,&n,&report)!=0;
    fail+=n!=1 || v[0].target!=1 || fabs(v[0].end-v[0].begin-8)>1e-12;
    /* Same spatial extent, but face-to-face instead of an open continuation. */
    f[1].outward[1]=1;
    fail+=MaterialFront_collect(arena,f,2,&v,&n,&report)!=0 || n!=0 || report.contacts==0;
    f[1].outward[1]=-1; f[1].xyz[1]=f[1].xyz[4]=20;
    fail+=MaterialFront_collect(arena,f,2,&v,&n,&report)!=0 || n!=0;
    fail+=MaterialFront_collect(arena,NULL,0,&v,&n,&report)!=0 || n!=0 || v!=NULL;
    {
        MaterialFront chain[8]={0}, original[8]={0};
        MaterialEvidenceFactor *factors=NULL;
        MaterialEvidenceReport evidence={0};
        size_t nf=0;
        for (int side=0;side<2;side++) for (int i=0;i<4;i++) {
            MaterialFront *edge=chain+4*side+i;
            edge->xyz[0]=(float)(side ? 32-8*i : 8*i);
            edge->xyz[3]=(float)(side ? 24-8*i : 8*i+8);
            edge->xyz[1]=edge->xyz[4]=(float)side;
            edge->normal[2]=1; edge->outward[1]=side ? -1.0f : 1.0f;
            edge->phase[0]=edge->phase[1]=(double)-side;
            edge->chart=side; edge->face=(uint64_t)(4*side+i);
            edge->source_edge=3*edge->face;
            edge->source_vertices[0]=5*side+i; edge->source_vertices[1]=5*side+i+1;
        }
        fail+=MaterialFront_trace(arena,chain,8)!=0;
        memcpy(original,chain,sizeof chain);
        for (int i=0;i<8;i++) chain[i]=original[7-i];
        fail+=MaterialFront_trace(arena,chain,8)!=0;
        for (int i=0;i<8;i++) {
            fail+=chain[i].curve!=original[7-i].curve || chain[i].region_base!=original[7-i].region_base ||
                  fabs(chain[i].arc-original[7-i].arc)>1e-12;
        }
        fail+=MaterialFront_collect(arena,chain,8,&v,&n,&report)!=0 || n!=4;
        fail+=MaterialEvidence_reduce(arena,v,n,MATERIAL_FRONT_LENGTH_SCALE,&factors,&nf,&evidence)!=0;
        fail+=nf!=1 || factors[0].units!=2 || fabs(factors[0].effective_support-2)>1e-12 || factors[0].significant;
        /* Paired intervals preserve orientation, including opposite halfedges. */
        for (size_t i=0;i<n;i++) fail+=v[i].begin_b<=v[i].end_b;
        /* Nonbranching original topology only; never choose a geometric next
         * edge when an original vertex has two successors. */
        chain[1].source_vertices[0]=chain[0].source_vertices[0];
        fail+=MaterialFront_trace(arena,chain,8)==0;
    }
    {
        /* One long original edge crosses several physical blocks. Both
         * source block boundaries split the interval, but its face cap is 1. */
        f[0].xyz[3]=40; f[1].xyz[0]=40;
        f[1].xyz[1]=f[1].xyz[4]=1;
        fail+=MaterialFront_trace(arena,f,2)!=0;
        fail+=MaterialFront_collect(arena,f,2,&v,&n,&report)!=0 || n!=5;
        MaterialEvidenceFactor *factors=NULL;
        size_t nf=0;
        fail+=MaterialEvidence_reduce(arena,v,n,MATERIAL_FRONT_LENGTH_SCALE,&factors,&nf,NULL)!=0;
        fail+=nf!=1 || fabs(factors[0].effective_support-1)>1e-12 || factors[0].significant;
    }
    {
        /* Frozen PHerc0139 z05760_y03328_x02560 edges 32298 / 142999.
         * Distinct source vertices share exactly one geometric endpoint.
         * The old overlap was 2.22e-16 on A, zero on B: it aborted 21x5x5. */
        MaterialFront point[2]={0};
        const float xyz[2][6]={
            {5801.3505859375f,3395.193603515625f,2658.986328125f,
             5801.12646484375f,3395.456787109375f,2658.85498046875f},
            {5797.21875f,3398.60009765625f,2658.55419921875f,
             5801.12646484375f,3395.456787109375f,2658.85498046875f}};
        const float normals[2][3]={{.5925931930541992f,.7018807530403137f,.39521729946136475f},
                                   {.5196273922920227f,.6937170624732971f,.49874240159988403f}};
        const float outward[2][3]={{.5305785536766052f,.029045334085822105f,-.847137987613678f},
                                   {-.35357099771499634f,-.35681331157684326f,.8646801710128784f}};
        for (int side=0;side<2;side++) {
            memcpy(point[side].xyz,xyz[side],sizeof xyz[side]);
            memcpy(point[side].normal,normals[side],sizeof normals[side]);
            memcpy(point[side].outward,outward[side],sizeof outward[side]);
            point[side].chart=side; point[side].face=(uint64_t)side;
            point[side].source_edge=(uint64_t)(3*side);
            point[side].source_vertices[0]=2*side; point[side].source_vertices[1]=2*side+1;
        }
        point[0].arc=167.75320979782359; point[1].arc=3.6435917893572989;
        fail+=MaterialFront_collect(arena,point,2,&v,&n,&report)!=0 || n!=0 ||
              report.degenerate_intervals!=1;
        /* Small but real float-representable overlap is retained, not
         * confused with the roundoff-scale zero-support witness. */
        point[0]=f[0]; point[1]=f[1];
        point[0].xyz[3]=8; point[1].xyz[0]=16;
        point[1].xyz[3]=nextafterf(8,0);
        point[0].arc=point[1].arc=0;
        fail+=MaterialFront_collect(arena,point,2,&v,&n,&report)!=0 || n!=1 ||
              report.degenerate_intervals!=0;
    }
    fprintf(stderr,"[selftest] material fronts %s (%d failures)\n",fail ? "FAIL" : "PASS",fail);
    Arena_dispose(&arena);
    return fail ? -1 : 0;
}
