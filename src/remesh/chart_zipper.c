/*
 * chart_zipper.c -- chart-level trim-boundary zipper for grid_weld.
 *
 * This is the deliberately small seam operation described by Turk/Levoy and
 * the contour-triangulation literature: recover ordered boundary paths, match
 * paths across one seam, and triangulate one complete strip.  It replaces the
 * former BPA -> patch repair -> fragment forest sequence.
 */
#include "chart_zipper.h"

#include "seam_planes.h"
#include "intersection_cleanup.h"
#include "../common/kdtree.h"
#include "../common/pipeline_constants.h"
#include "../common/union_find.h"
#include "../flatten/topology_invariants.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/pitch_table.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define CZ_PARALLEL_EPS 0.75
#define CZ_INF (DBL_MAX / 16.0)

typedef struct {
    int32_t lo, hi;
    int32_t a, b, opposite;
    int32_t face;
} CzHalfEdge;

typedef struct {
    int32_t a, b;                 /* directed source boundary half-edge */
    int32_t opposite;
    int32_t chart;
    int32_t plane;
    uint8_t side;
} CzBoundaryEdge;

typedef struct {
    int32_t *verts;               /* ordered path vertices */
    int32_t *edges;               /* boundary edge ids, n-1 for open paths */
    size_t n;
    int32_t chart;
    int32_t plane;
    uint8_t side;
    uint8_t closed;
    double length;
    double centroid[3];
    double normal[3];
} CzChain;

typedef struct {
    int32_t chain_a, chain_b;     /* a is low side, b is high side */
    size_t a0, a1, b0, b1;       /* inclusive vertex intervals */
    size_t hits;
    uint8_t driver;               /* 0: a nominated b; 1: b nominated a */
    double full_coverage;
    double span_coverage;
    double support;
    double mean_gap;
    double score;
} CzCandidate;

typedef struct {
    size_t candidate;
    size_t face0;
    size_t nf;
} CzAcceptedTx;

typedef struct {
    CzChain *v;
    size_t n, cap;
} CzChainVec;

typedef struct {
    CzCandidate *v;
    size_t n, cap;
} CzCandidateVec;

typedef struct {
    float *points;
    int32_t *chain;
    int32_t *pos;
    size_t n;
    KDTree_T tree;
} CzPointIndex;

static int cz_halfedge_cmp(const void *aa, const void *bb)
{
    const CzHalfEdge *a = (const CzHalfEdge *)aa;
    const CzHalfEdge *b = (const CzHalfEdge *)bb;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    return 0;
}

static int cz_candidate_cmp(const void *aa, const void *bb)
{
    const CzCandidate *a = (const CzCandidate *)aa;
    const CzCandidate *b = (const CzCandidate *)bb;
    if (a->full_coverage != b->full_coverage)
        return a->full_coverage > b->full_coverage ? -1 : 1;
    if (a->span_coverage != b->span_coverage)
        return a->span_coverage > b->span_coverage ? -1 : 1;
    if (a->mean_gap != b->mean_gap)
        return a->mean_gap < b->mean_gap ? -1 : 1;
    if (a->support != b->support)
        return a->support > b->support ? -1 : 1;
    if (a->score != b->score) return a->score > b->score ? -1 : 1;
    if (a->chain_a != b->chain_a) return a->chain_a < b->chain_a ? -1 : 1;
    return a->chain_b < b->chain_b ? -1 : (a->chain_b > b->chain_b ? 1 : 0);
}

static int cz_chain_push(CzChainVec *v, const CzChain *c)
{
    if (v->n == v->cap) {
        size_t cap = v->cap ? v->cap * 2 : 64;
        CzChain *p = (CzChain *)realloc(v->v, cap * sizeof(*p));
        if (!p) return -1;
        v->v = p; v->cap = cap;
    }
    v->v[v->n++] = *c;
    return 0;
}

static int cz_candidate_push(CzCandidateVec *v, const CzCandidate *c)
{
    if (v->n == v->cap) {
        size_t cap = v->cap ? v->cap * 2 : 128;
        CzCandidate *p = (CzCandidate *)realloc(v->v, cap * sizeof(*p));
        if (!p) return -1;
        v->v = p; v->cap = cap;
    }
    v->v[v->n++] = *c;
    return 0;
}

static int cz_face_reserve(int32_t **faces, size_t *capacity, size_t need)
{
    size_t cap;
    int32_t *grown;
    if (need <= *capacity) return 0;
    if (need > SIZE_MAX / (3 * sizeof(**faces))) return -1;
    cap = *capacity ? *capacity : 1024;
    while (cap < need) {
        size_t next = cap <= SIZE_MAX / 2 ? cap * 2 : need;
        if (next < cap || next > SIZE_MAX / (3 * sizeof(**faces))) {
            cap = need;
            break;
        }
        cap = next;
    }
    grown = (int32_t *)realloc(*faces, cap * 3 * sizeof(**faces));
    if (!grown) return -1;
    *faces = grown;
    *capacity = cap;
    return 0;
}

static double cz_dist2(const float *verts, int32_t a, int32_t b)
{
    double d0 = (double)verts[(size_t)a*3+0] - verts[(size_t)b*3+0];
    double d1 = (double)verts[(size_t)a*3+1] - verts[(size_t)b*3+1];
    double d2 = (double)verts[(size_t)a*3+2] - verts[(size_t)b*3+2];
    return d0*d0 + d1*d1 + d2*d2;
}

static double cz_dot(const double *a, const double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void cz_compute_normals(const float *verts, size_t nv,
                               const int32_t *faces, size_t nf,
                               double *normal)
{
    memset(normal, 0, nv * 3 * sizeof(*normal));
    for (size_t f = 0; f < nf; f++) {
        int32_t a=faces[f*3], b=faces[f*3+1], c=faces[f*3+2];
        const float *pa=&verts[(size_t)a*3];
        const float *pb=&verts[(size_t)b*3];
        const float *pc=&verts[(size_t)c*3];
        double u0=pb[0]-pa[0],u1=pb[1]-pa[1],u2=pb[2]-pa[2];
        double v0=pc[0]-pa[0],v1=pc[1]-pa[1],v2=pc[2]-pa[2];
        double n[3]={u1*v2-u2*v1,u2*v0-u0*v2,u0*v1-u1*v0};
        int32_t q[3]={a,b,c};
        for(int k=0;k<3;k++) for(int d=0;d<3;d++)
            normal[(size_t)q[k]*3+(size_t)d]+=n[d];
    }
    for(size_t v=0;v<nv;v++){
        double *n=&normal[v*3];
        double m=sqrt(cz_dot(n,n));
        if(m>1e-20){n[0]/=m;n[1]/=m;n[2]/=m;}
    }
}

static int cz_gate_armed(const BpaBridgeGate *g)
{
    return g && g->pitch > 0.0 &&
           ((g->axis_z && g->axis_y && g->axis_x && g->axis_n >= 2) ||
            g->umb_y != 0.0 || g->umb_x != 0.0);
}

static double cz_wrap_pmpi(double angle)
{
    while(angle>M_PI) angle-=2.0*M_PI;
    while(angle<=-M_PI) angle+=2.0*M_PI;
    return angle;
}

/* Branch-cut-free winding difference along one proposed zipper edge.  Taking
 * two absolute atan2-based winding coordinates and subtracting them creates a
 * false one-turn jump whenever a perfectly valid edge crosses -pi/+pi.  The
 * local edge invariant is instead
 *
 *   dw = (r_b-r_a)/pitch - wrap(theta_b-theta_a)/(2*pi).
 *
 * Evaluate the curved axis independently at each endpoint's z. */
static double cz_winding_delta(const BpaBridgeGate *g,
                               const float *pa,const float *pb)
{
    double uya,uxa,uyb,uxb;
    BpaBridgeGate_axis_at(g,(double)pa[0],&uya,&uxa);
    BpaBridgeGate_axis_at(g,(double)pb[0],&uyb,&uxb);
    double ya=(double)pa[1]-uya, xa=(double)pa[2]-uxa;
    double yb=(double)pb[1]-uyb, xb=(double)pb[2]-uxb;
    double dtheta=cz_wrap_pmpi(atan2(yb,xb)-atan2(ya,xa));
    return PitchTable_dturns(PitchTable_from_env(),hypot(ya,xa),
                             hypot(yb,xb),g->pitch)-dtheta/(2.0*M_PI);
}

static int cz_gate_pair(const BpaBridgeGate *g,
                        const float *verts, int32_t a, int32_t b)
{
    if(!cz_gate_armed(g)) return 1;
    double dw=cz_winding_delta(g,&verts[(size_t)a*3],
                               &verts[(size_t)b*3]);
    double limit=g->hard>0.0?g->hard:(g->tol>0.0?g->tol:0.35);
    return fabs(dw)<=limit;
}

void ChartZipper_default_params(ChartZipperParams *p)
{
    if(!p) return;
    p->cube_size=128.0f;
    p->band=6.0f;
    p->max_cross_edge=2.0f*(float)BRIDGE_RHO_MAX;
    p->normal_dot_min=0.35f;
    p->min_coverage=0.45f;
    p->min_support_length=2.0f;
    p->ambiguity_ratio=0.80f;
    p->min_arc_ratio=0.75f;
    p->skinny_weight=0.20f;
    p->min_triangle_altitude=0.10f;
    /* Exact embeddedness is the hard transaction gate. A positive slab is a
     * useful diagnostic, but on noisy papyrus it also labels close, distinct
     * layers and would reject otherwise safe welds. */
    p->conflict_gap=0.0f;
    p->conflict_parallel_angle_deg=20.0f;
    p->min_hits=2;
    p->ball_results=96;
    p->trace=0;
}

static int cz_best_plane(const float *verts, const CzHalfEdge *h,
                         const SeamPlane *planes, size_t np, double band,
                         int *out_plane, int *out_side)
{
    int best=-1;
    double best_d=DBL_MAX;
    for(size_t p=0;p<np;p++){
        int ax=planes[p].axis;
        double ca=verts[(size_t)h->a*3+(size_t)ax];
        double cb=verts[(size_t)h->b*3+(size_t)ax];
        double da=fabs(ca-planes[p].coord), db=fabs(cb-planes[p].coord);
        if(da>=band||db>=band||fabs(ca-cb)>=CZ_PARALLEL_EPS) continue;
        double d=da>db?da:db;
        if(d<best_d){best=(int)p;best_d=d;}
    }
    if(best<0) return 0;
    {
        int ax=planes[best].axis;
        double co=planes[best].coord;
        double oc=verts[(size_t)h->opposite*3+(size_t)ax];
        double mc=0.5*((double)verts[(size_t)h->a*3+(size_t)ax]+
                       (double)verts[(size_t)h->b*3+(size_t)ax]);
        double side_coord=fabs(oc-co)>1e-9?oc:mc;
        *out_plane=best;
        *out_side=side_coord>co?1:0;
    }
    return 1;
}

static void cz_reverse_chain(CzChain *c)
{
    for(size_t i=0;i<c->n/2;i++){
        int32_t t=c->verts[i];c->verts[i]=c->verts[c->n-1-i];c->verts[c->n-1-i]=t;
    }
    if(c->n>1) for(size_t i=0;i<(c->n-1)/2;i++){
        int32_t t=c->edges[i];
        c->edges[i]=c->edges[c->n-2-i];
        c->edges[c->n-2-i]=t;
    }
}

static void cz_finish_chain(const float *verts, const double *normal,
                            const CzBoundaryEdge *be, CzChain *c)
{
    size_t forward=0,backward=0;
    for(size_t i=0;i+1<c->n;i++){
        const CzBoundaryEdge *e=&be[c->edges[i]];
        if(e->a==c->verts[i]&&e->b==c->verts[i+1]) forward++;
        else backward++;
    }
    if(backward>forward) cz_reverse_chain(c);
    for(size_t i=0;i<c->n;i++){
        int32_t v=c->verts[i];
        for(int d=0;d<3;d++){
            c->centroid[d]+=verts[(size_t)v*3+(size_t)d];
            c->normal[d]+=normal[(size_t)v*3+(size_t)d];
        }
        if(i+1<c->n) c->length+=sqrt(cz_dist2(verts,v,c->verts[i+1]));
    }
    if(c->n) for(int d=0;d<3;d++) c->centroid[d]/=(double)c->n;
    {
        double m=sqrt(cz_dot(c->normal,c->normal));
        if(m>1e-20) for(int d=0;d<3;d++) c->normal[d]/=m;
    }
}

static int cz_extract_chains(const float *verts,size_t nv,
                             const double *normal,
                             const CzBoundaryEdge *be,size_t nbe,
                             size_t plane,int side,
                             CzChainVec *chains,
                             ChartZipperStats *stats)
{
    int32_t *adj0=NULL,*adj1=NULL;
    uint8_t *visited=NULL,*branched=NULL;
    int rc=-1;
    adj0=(int32_t*)malloc(nv*sizeof(*adj0));
    adj1=(int32_t*)malloc(nv*sizeof(*adj1));
    visited=(uint8_t*)calloc(nbe?nbe:1,1);
    branched=(uint8_t*)calloc(nv?nv:1,1);
    if(!adj0||!adj1||!visited||!branched) goto done;
    for(size_t v=0;v<nv;v++) adj0[v]=adj1[v]=-1;
    for(size_t e=0;e<nbe;e++) if(be[e].plane==(int32_t)plane&&be[e].side==(uint8_t)side){
        int32_t q[2]={be[e].a,be[e].b};
        for(int k=0;k<2;k++){
            int32_t v=q[k];
            if(adj0[v]<0) adj0[v]=(int32_t)e;
            else if(adj1[v]<0) adj1[v]=(int32_t)e;
            else branched[v]=1;
        }
    }
    for(int pass=0;pass<2;pass++) for(size_t seed=0;seed<nbe;seed++){
        if(visited[seed]||be[seed].plane!=(int32_t)plane||be[seed].side!=(uint8_t)side) continue;
        int32_t ends[2]={be[seed].a,be[seed].b};
        int32_t start=-1;
        for(int k=0;k<2;k++){
            int32_t v=ends[k];
            int deg=(adj0[v]>=0)+(adj1[v]>=0);
            if(pass==0&&deg==1) start=v;
            if(pass==1&&deg==2) start=v;
        }
        if(start<0) continue;
        size_t cap=32,nvpath=0,nepath=0;
        int32_t *pv=(int32_t*)malloc(cap*sizeof(*pv));
        int32_t *pe=(int32_t*)malloc(cap*sizeof(*pe));
        int32_t cur=start,prev_edge=-1;
        int closed=0,bad=0;
        if(!pv||!pe){free(pv);free(pe);goto done;}
        pv[nvpath++]=cur;
        for(;;){
            if(branched[cur]){bad=1;break;}
            int32_t e0=adj0[cur],e1=adj1[cur],e=-1;
            if(e0>=0&&e0!=prev_edge&&!visited[e0]) e=e0;
            else if(e1>=0&&e1!=prev_edge&&!visited[e1]) e=e1;
            if(e<0) break;
            visited[e]=1;
            int32_t nxt=be[e].a==cur?be[e].b:be[e].a;
            if(nvpath+1>=cap){
                cap*=2;
                int32_t *qv=(int32_t*)realloc(pv,cap*sizeof(*pv));
                int32_t *qe=(int32_t*)realloc(pe,cap*sizeof(*pe));
                if(!qv||!qe){free(qv?qv:pv);free(qe?qe:pe);goto done;}
                pv=qv;pe=qe;
            }
            pe[nepath++]=e;
            if(nxt==start){closed=1;break;}
            pv[nvpath++]=nxt;
            prev_edge=e;cur=nxt;
            if(nepath>nbe){bad=1;break;}
        }
        if(bad||nvpath<2){
            if(bad&&stats) stats->branched_chains_rejected++;
            free(pv);free(pe);continue;
        }
        {
            CzChain c;
            memset(&c,0,sizeof c);
            c.verts=pv;c.edges=pe;c.n=nvpath;c.closed=(uint8_t)closed;
            c.chart=be[seed].chart;c.plane=(int32_t)plane;c.side=(uint8_t)side;
            for(size_t i=0;i<nepath;i++) if(be[pe[i]].chart!=c.chart) bad=1;
            if(bad){free(pv);free(pe);continue;}
            cz_finish_chain(verts,normal,be,&c);
            if(cz_chain_push(chains,&c)!=0){free(pv);free(pe);goto done;}
            if(stats){stats->chains++;if(closed)stats->closed_chains_rejected++;else stats->open_chains++;}
        }
    }
    rc=0;
done:
    free(branched);free(visited);free(adj1);free(adj0);
    return rc;
}

static double cz_arc_length(const float *verts,const CzChain *c,size_t i0,size_t i1)
{
    double s=0.0;
    if(i0>i1){size_t t=i0;i0=i1;i1=t;}
    for(size_t i=i0;i<i1;i++) s+=sqrt(cz_dist2(verts,c->verts[i],c->verts[i+1]));
    return s;
}

/* Endpoint cost for a candidate pair of boundary intervals.  The winding and
 * distance checks are the same hard gates used by the final strip. */
static double cz_interval_endpoint_cost(const float *verts,
                                        const CzChain *a,size_t a0,size_t a1,
                                        const CzChain *b,size_t b0,size_t b1,
                                        const BpaBridgeGate *gate,
                                        const ChartZipperParams *p)
{
    double max2=(double)p->max_cross_edge*p->max_cross_edge;
    int32_t av[2]={a->verts[a0],a->verts[a1]};
    int32_t bv[2]={b->verts[b0],b->verts[b1]};
    double best=DBL_MAX;
    for(int reverse=0;reverse<2;reverse++){
        int32_t x=bv[reverse?1:0],y=bv[reverse?0:1];
        double d0=cz_dist2(verts,av[0],x),d1=cz_dist2(verts,av[1],y);
        if(d0>max2||d1>max2||!cz_gate_pair(gate,verts,av[0],x)||
           !cz_gate_pair(gate,verts,av[1],y)) continue;
        if(d0+d1<best)best=d0+d1;
    }
    return best;
}

/* Nearest-point voting can collapse a long, curled child contour onto one edge
 * of a smoother parent contour.  Zipping those raw bounds would require a fan
 * of long diagonals.  Grow only the shorter interval, one existing boundary
 * edge at a time, until the arc lengths are comparable.  Every proposed growth
 * must keep both paired endpoints inside the ordinary cross-edge and winding
 * gates; source geometry is never split or moved. */
static int cz_balance_intervals(const float *verts,
                                const CzChain *a,size_t *a0,size_t *a1,
                                const CzChain *b,size_t *b0,size_t *b1,
                                const BpaBridgeGate *gate,
                                const ChartZipperParams *p)
{
    int changed=0;
    for(size_t round=0;round<a->n+b->n;round++){
        double la=cz_arc_length(verts,a,*a0,*a1);
        double lb=cz_arc_length(verts,b,*b0,*b1);
        double small=la<lb?la:lb,large=la<lb?lb:la;
        if(large<=1e-12||small>=large*(double)p->min_arc_ratio)break;
        const CzChain *s=la<lb?a:b;
        size_t *s0=la<lb?a0:b0,*s1=la<lb?a1:b1;
        double best=DBL_MAX;int which=0;
        if(*s0>0){
            size_t q0=*s0-1,q1=*s1;
            double ec=la<lb?
                cz_interval_endpoint_cost(verts,a,q0,q1,b,*b0,*b1,gate,p):
                cz_interval_endpoint_cost(verts,a,*a0,*a1,b,q0,q1,gate,p);
            if(ec<DBL_MAX){
                double nl=cz_arc_length(verts,s,q0,q1);
                double cost=ec+0.10*fabs(nl-large);
                if(cost<best){best=cost;which=-1;}
            }
        }
        if(*s1+1<s->n){
            size_t q0=*s0,q1=*s1+1;
            double ec=la<lb?
                cz_interval_endpoint_cost(verts,a,q0,q1,b,*b0,*b1,gate,p):
                cz_interval_endpoint_cost(verts,a,*a0,*a1,b,q0,q1,gate,p);
            if(ec<DBL_MAX){
                double nl=cz_arc_length(verts,s,q0,q1);
                double cost=ec+0.10*fabs(nl-large);
                if(cost<best){best=cost;which=1;}
            }
        }
        if(which==0)break;
        if(which<0)(*s0)--;else (*s1)++;
        changed=1;
    }
    return changed;
}

static int cz_point_index(Arena_T arena,const float *verts,
                          const CzChain *chains,size_t nc,
                          int32_t plane,int side,CzPointIndex *out)
{
    size_t n=0,k=0;
    memset(out,0,sizeof(*out));
    for(size_t c=0;c<nc;c++) if(!chains[c].closed&&chains[c].plane==plane&&chains[c].side==(uint8_t)side) n+=chains[c].n;
    if(n==0) return 0;
    out->points=(float*)ARENA_ALLOC(arena,(n*3*sizeof(float)));
    out->chain=(int32_t*)ARENA_ALLOC(arena,(n*sizeof(int32_t)));
    out->pos=(int32_t*)ARENA_ALLOC(arena,(n*sizeof(int32_t)));
    for(size_t c=0;c<nc;c++) if(!chains[c].closed&&chains[c].plane==plane&&chains[c].side==(uint8_t)side){
        for(size_t i=0;i<chains[c].n;i++){
            int32_t v=chains[c].verts[i];
            memcpy(&out->points[k*3],&verts[(size_t)v*3],3*sizeof(float));
            out->chain[k]=(int32_t)c;out->pos[k]=(int32_t)i;k++;
        }
    }
    out->n=n;out->tree=KDTree_new(arena,out->points,n);
    return 0;
}

static int cz_match_compatible(const float *verts,const double *normal,
                               int32_t a,int32_t b,const BpaBridgeGate *gate,
                               const ChartZipperParams *p)
{
    const double *na=&normal[(size_t)a*3],*nb=&normal[(size_t)b*3];
    double ma=cz_dot(na,na),mb=cz_dot(nb,nb);
    if(ma>0.25&&mb>0.25&&fabs(cz_dot(na,nb))<(double)p->normal_dot_min) return 0;
    return cz_gate_pair(gate,verts,a,b);
}

/* A small chart clipped wholly inside the seam band presents its complete
 * boundary loop, rather than an open path.  Rejecting every such loop loses a
 * safe continuation at block corners.  Turn it into one ordinary boundary
 * interval by omitting the loop edge farthest from all geometry on the
 * opposite side.  This changes only contour bookkeeping: source vertices,
 * source faces, and even the omitted boundary edge remain untouched.  The
 * ordinary forest, embeddedness, and disk certificates still gate the entire
 * zipper transaction.
 *
 * At least min_hits vertices must have a compatible opposite-side neighbour
 * inside max_cross_edge.  Thus a merely nearby closed island remains rejected.
 */
static int cz_open_supported_closed_chains(const float *verts,
                                           const double *normal,
                                           CzChain *chains,size_t nc,
                                           int32_t plane,int side,
                                           const CzPointIndex *target,
                                           const BpaBridgeGate *gate,
                                           const ChartZipperParams *p,
                                           const CzBoundaryEdge *be,
                                           ChartZipperStats *stats)
{
    int32_t *results=NULL;
    if(target->n==0) return 0;
    results=(int32_t*)malloc(p->ball_results*sizeof(*results));
    if(!results) return -1;
    for(size_t ci=0;ci<nc;ci++){
        CzChain *c=&chains[ci];
        double *near2=NULL;
        size_t compatible=0,break_edge=0;
        double best_min=-1.0,best_sum=-1.0;
        int32_t *new_verts=NULL,*new_edges=NULL;
        if(!c->closed||c->plane!=plane||c->side!=(uint8_t)side||c->n<3)
            continue;
        near2=(double*)malloc(c->n*sizeof(*near2));
        if(!near2) { free(results); return -1; }
        for(size_t i=0;i<c->n;i++){
            int32_t dv=c->verts[i];
            float q[3]={verts[(size_t)dv*3+0],verts[(size_t)dv*3+1],
                        verts[(size_t)dv*3+2]};
            float nearest_sq=0.0f;
            (void)KDTree_nearest(target->tree,q,&nearest_sq);
            near2[i]=(double)nearest_sq;
            {
                size_t nr=KDTree_ball_query(
                    target->tree,q,p->max_cross_edge*p->max_cross_edge,
                    results,p->ball_results);
                int found=0;
                for(size_t r=0;r<nr&&r<p->ball_results;r++){
                    int32_t pi=results[r];
                    int32_t tc,tv;
                    if(pi<0||(size_t)pi>=target->n) continue;
                    tc=target->chain[pi];
                    if(tc<0||(size_t)tc>=nc||chains[tc].chart==c->chart)
                        continue;
                    tv=chains[tc].verts[target->pos[pi]];
                    if(cz_match_compatible(verts,normal,dv,tv,gate,p)){
                        found=1;break;
                    }
                }
                if(found) compatible++;
            }
        }
        if(compatible<p->min_hits){
            free(near2);continue;
        }
        for(size_t e=0;e<c->n;e++){
            double da=near2[e],db=near2[(e+1)%c->n];
            double mn=da<db?da:db,sum=da+db;
            if(mn>best_min||(mn==best_min&&sum>best_sum)){
                best_min=mn;best_sum=sum;break_edge=e;
            }
        }
        new_verts=(int32_t*)malloc(c->n*sizeof(*new_verts));
        new_edges=(int32_t*)malloc(c->n*sizeof(*new_edges));
        if(!new_verts||!new_edges){
            free(new_edges);free(new_verts);free(near2);free(results);
            return -1;
        }
        for(size_t i=0;i<c->n;i++)
            new_verts[i]=c->verts[(break_edge+1+i)%c->n];
        for(size_t i=0;i+1<c->n;i++)
            new_edges[i]=c->edges[(break_edge+1+i)%c->n];
        free(c->verts);free(c->edges);
        c->verts=new_verts;c->edges=new_edges;c->closed=0;
        c->length=0.0;
        memset(c->centroid,0,sizeof c->centroid);
        memset(c->normal,0,sizeof c->normal);
        cz_finish_chain(verts,normal,be,c);
        if(stats){
            stats->closed_chains_opened++;
            stats->open_chains++;
            if(stats->closed_chains_rejected>0)
                stats->closed_chains_rejected--;
        }
        if(p->trace)fprintf(stderr,
            "    [zip trace] open closed path=%zu chart=%d plane=%d "
            "break_edge=%zu far_gap=%.3f compatible=%zu/%zu\n",
            ci,c->chart,c->plane,break_edge,
            best_min>=0.0?sqrt(best_min):0.0,compatible,c->n);
        free(near2);
    }
    free(results);
    return 0;
}

static int cz_generate_direction(Arena_T arena,const float *verts,
                                 const double *normal,
                                 const CzChain *chains,size_t nc,
                                 int32_t plane,int driver_side,
                                 const CzPointIndex *target,
                                 const BpaBridgeGate *gate,
                                 const ChartZipperParams *p,
                                 CzCandidateVec *candidate,
                                 ChartZipperStats *stats)
{
    (void)arena;
    size_t *hits=NULL,*dmin=NULL,*dmax=NULL,*tmin=NULL,*tmax=NULL;
    double *sumgap=NULL;
    int32_t *results=NULL,*touched=NULL;
    int rc=-1;
    if(target->n==0) return 0;
    hits=(size_t*)calloc(nc?nc:1,sizeof(*hits));
    dmin=(size_t*)malloc((nc?nc:1)*sizeof(*dmin));
    dmax=(size_t*)calloc(nc?nc:1,sizeof(*dmax));
    tmin=(size_t*)malloc((nc?nc:1)*sizeof(*tmin));
    tmax=(size_t*)calloc(nc?nc:1,sizeof(*tmax));
    sumgap=(double*)calloc(nc?nc:1,sizeof(*sumgap));
    results=(int32_t*)malloc(p->ball_results*sizeof(*results));
    touched=(int32_t*)malloc((nc?nc:1)*sizeof(*touched));
    if(!hits||!dmin||!dmax||!tmin||!tmax||!sumgap||!results||!touched) goto done;
    for(size_t c=0;c<nc;c++){
        const CzChain *d=&chains[c];
        if(d->closed||d->plane!=plane||d->side!=(uint8_t)driver_side||d->n<2) continue;
        size_t ntouched=0;
        for(size_t i=0;i<d->n;i++){
            int32_t dv=d->verts[i];
            float q[3]={verts[(size_t)dv*3],verts[(size_t)dv*3+1],verts[(size_t)dv*3+2]};
            size_t nr=KDTree_ball_query(target->tree,q,p->max_cross_edge*p->max_cross_edge,
                                        results,p->ball_results);
            double best=DBL_MAX;int32_t best_chain=-1,best_pos=-1;
            for(size_t r=0;r<nr&&r<p->ball_results;r++){
                int32_t pi=results[r];
                if(pi<0||(size_t)pi>=target->n) continue;
                int32_t tc=target->chain[pi];
                if(tc<0||(size_t)tc>=nc||chains[tc].chart==d->chart) continue;
                int32_t tv=chains[tc].verts[target->pos[pi]];
                double dd=cz_dist2(verts,dv,tv);
                if(dd<best&&cz_match_compatible(verts,normal,dv,tv,gate,p)){
                    best=dd;best_chain=tc;best_pos=target->pos[pi];
                }
            }
            if(best_chain<0) continue;
            size_t tc=(size_t)best_chain;
            if(hits[tc]==0){
                touched[ntouched++]=best_chain;
                dmin[tc]=dmax[tc]=i;tmin[tc]=tmax[tc]=(size_t)best_pos;
            }else{
                if(i<dmin[tc])dmin[tc]=i;if(i>dmax[tc])dmax[tc]=i;
                if((size_t)best_pos<tmin[tc])tmin[tc]=(size_t)best_pos;
                if((size_t)best_pos>tmax[tc])tmax[tc]=(size_t)best_pos;
            }
            hits[tc]++;sumgap[tc]+=sqrt(best);
        }
        int32_t best_tc=-1,second_tc=-1;
        double best_score=-DBL_MAX,second_score=-DBL_MAX;
        for(size_t z=0;z<ntouched;z++){
            size_t tc=(size_t)touched[z];
            size_t dspan=dmax[tc]-dmin[tc]+1;
            size_t tspan=tmax[tc]-tmin[tc]+1;
            double full=(double)hits[tc]/(double)d->n;
            double span=(double)hits[tc]/(double)(dspan>tspan?dspan:tspan);
            double support=cz_arc_length(verts,d,dmin[tc],dmax[tc]);
            double ts=cz_arc_length(verts,&chains[tc],tmin[tc],tmax[tc]);
            if(ts<support)support=ts;
            double mean=sumgap[tc]/(double)hits[tc];
            double score=100.0*full+25.0*span+support-2.0*mean;
            if(score>best_score){second_score=best_score;second_tc=best_tc;best_score=score;best_tc=(int32_t)tc;}
            else if(score>second_score){second_score=score;second_tc=(int32_t)tc;}
        }
        if(best_tc>=0){
            size_t tc=(size_t)best_tc;
            size_t dspan=dmax[tc]-dmin[tc]+1;
            size_t tspan=tmax[tc]-tmin[tc]+1;
            size_t di0=dmin[tc],di1=dmax[tc];
            size_t ti0=tmin[tc],ti1=tmax[tc];
            double full=(double)hits[tc]/(double)d->n;
            double span=(double)hits[tc]/(double)(dspan>tspan?dspan:tspan);
            int balanced=cz_balance_intervals(verts,d,&di0,&di1,
                                               &chains[tc],&ti0,&ti1,
                                               gate,p);
            double support=cz_arc_length(verts,d,di0,di1);
            double ts=cz_arc_length(verts,&chains[tc],ti0,ti1);
            if(ts<support)support=ts;
            int ambiguous=second_tc>=0&&second_score>=best_score*(double)p->ambiguity_ratio;
            if(ambiguous){
                if(stats)stats->ambiguous_candidates++;
                if(p->trace) fprintf(stderr,
                    "    [zip trace] ambiguous driver chart=%d path=%zu side=%d "
                    "best=%d score=%.3f second=%d score=%.3f\n",
                    d->chart,c,driver_side,best_tc,best_score,
                    second_tc,second_score);
            }
            else if(hits[tc]>=p->min_hits&&span>=p->min_coverage&&
                    support>=p->min_support_length&&dspan>=2&&tspan>=2){
                CzCandidate cc;
                memset(&cc,0,sizeof cc);
                if(driver_side==0){
                    cc.chain_a=(int32_t)c;cc.a0=di0;cc.a1=di1;
                    cc.chain_b=best_tc;cc.b0=ti0;cc.b1=ti1;
                }else{
                    cc.chain_a=best_tc;cc.a0=ti0;cc.a1=ti1;
                    cc.chain_b=(int32_t)c;cc.b0=di0;cc.b1=di1;
                }
                cc.hits=hits[tc];cc.full_coverage=full;cc.span_coverage=span;
                cc.driver=(uint8_t)driver_side;
                cc.support=support;cc.mean_gap=sumgap[tc]/(double)hits[tc];cc.score=best_score;
                if(cz_candidate_push(candidate,&cc)!=0) goto done;
                if(balanced&&stats)stats->balanced_candidates++;
                if(p->trace) fprintf(stderr,
                    "    [zip trace] candidate charts=%d:%d paths=%d:%d "
                    "span=%zu/%zu:%zu/%zu hits=%zu full=%.3f local=%.3f "
                    "support=%.3f gap=%.3f score=%.3f\n",
                    chains[cc.chain_a].chart,chains[cc.chain_b].chart,
                    cc.chain_a,cc.chain_b,cc.a0,cc.a1,cc.b0,cc.b1,
                    cc.hits,cc.full_coverage,cc.span_coverage,
                    cc.support,cc.mean_gap,cc.score);
            }else if(p->trace){
                fprintf(stderr,
                    "    [zip trace] weak driver chart=%d path=%zu -> chart=%d "
                    "hits=%zu spans=%zu:%zu full=%.3f local=%.3f "
                    "support=%.3f gap=%.3f\n",
                    d->chart,c,chains[tc].chart,hits[tc],dspan,tspan,
                    full,span,support,sumgap[tc]/(double)hits[tc]);
            }
        }
        for(size_t z=0;z<ntouched;z++){
            size_t tc=(size_t)touched[z];hits[tc]=0;sumgap[tc]=0.0;
        }
    }
    rc=0;
done:
    free(touched);free(results);free(sumgap);free(tmax);free(tmin);free(dmax);free(dmin);free(hits);
    return rc;
}

typedef struct { size_t a,b; } CzMatchPosition;

static int cz_match_position_cmp(const void *pa,const void *pb)
{
    const CzMatchPosition *a=(const CzMatchPosition*)pa;
    const CzMatchPosition *b=(const CzMatchPosition*)pb;
    if(a->a!=b->a)return a->a<b->a?-1:1;
    return a->b<b->b?-1:a->b>b->b?1:0;
}

/* Fraction of a driver's nearest-neighbour votes that agree with one of the
 * two possible boundary orientations.  A zipper can only follow a monotone
 * walk through the two ordered contours; this is the lineage certificate we
 * need from welding, without constructing its triangle strip. */
static double cz_order_coverage(const size_t *value,size_t n,size_t *tails)
{
    size_t best=0;
    if(n==0)return 0.0;
    for(int reverse=0;reverse<2;reverse++){
        size_t length=0;
        for(size_t i=0;i<n;i++){
            size_t x=reverse?SIZE_MAX-value[i]:value[i];
            size_t lo=0,hi=length;
            /* upper_bound admits repeated target vertices: advancing several
             * driver vertices against one target vertex is a valid zipper
             * walk. */
            while(lo<hi){
                size_t mid=lo+(hi-lo)/2;
                if(tails[mid]<=x)lo=mid+1;else hi=mid;
            }
            tails[lo]=x;
            if(lo==length)length++;
        }
        if(length>best)best=length;
    }
    return (double)best/(double)n;
}

static int cz_zip_dp(const float *verts,const double *normal,
                     const CzChain *a,size_t a0,size_t a1,
                     const CzChain *b,size_t b0,size_t b1,
                     const BpaBridgeGate *gate,const ChartZipperParams *p,
                     int32_t *out,size_t *out_nf);

static int cz_zip_trimmed(const float *verts,const double *normal,
                          const CzChain *a,size_t *a0,size_t *a1,
                          const CzChain *b,size_t *b0,size_t *b1,
                          const BpaBridgeGate *gate,
                          const ChartZipperParams *p,
                          int32_t *out,size_t *out_nf);

/* Turn a scored interval pair into compact reciprocal nearest-neighbour
 * samples.  This deliberately stops before triangulation.  Using the actual
 * compatible hits, rather than the candidate's conservative min/max bounds,
 * keeps disjoint child intervals on one large chart port disjoint. */
static int cz_emit_matches(Arena_T arena,const float *verts,
                           const double *normal,const CzChain *chains,
                           const CzCandidate *candidate,size_t ncandidate,
                           const BpaBridgeGate *gate,
                           const ChartZipperParams *p,
                           ChartZipperMatch **out_match,size_t *out_nmatch,
                           ChartZipperMatchSample **out_sample,
                           size_t *out_nsample)
{
    ChartZipperMatch *match=NULL;
    ChartZipperMatchSample *sample=NULL;
    size_t maximum_samples=0,nmatch=0,nsample=0;
    if(!out_match||!out_nmatch||!out_sample||!out_nsample)return -1;
    *out_match=NULL;*out_nmatch=0;*out_sample=NULL;*out_nsample=0;
    for(size_t ci=0;ci<ncandidate;ci++){
        const CzCandidate *cc=&candidate[ci];
        if(cc->a1<cc->a0||cc->b1<cc->b0)return -1;
        size_t na=cc->a1-cc->a0+1,nb=cc->b1-cc->b0+1;
        if(na>SIZE_MAX-nb||na+nb>SIZE_MAX-maximum_samples)return -1;
        maximum_samples+=na+nb;
    }
    match=(ChartZipperMatch*)ARENA_ALLOC(
        arena,(ncandidate?ncandidate:1)*sizeof(*match));
    sample=(ChartZipperMatchSample*)ARENA_ALLOC(
        arena,(maximum_samples?maximum_samples:1)*sizeof(*sample));
    for(size_t ci=0;ci<ncandidate;ci++){
        const CzCandidate *cc=&candidate[ci];
        const CzChain *a=&chains[cc->chain_a],*b=&chains[cc->chain_b];
        size_t a0=cc->a0,a1=cc->a1,b0=cc->b0,b1=cc->b1,walk_faces=0;
        if(cz_zip_dp(verts,normal,a,a0,a1,b,b0,b1,gate,p,NULL,
                     &walk_faces)!=0||walk_faces==0){
            walk_faces=0;
            if(cz_zip_trimmed(verts,normal,a,&a0,&a1,b,&b0,&b1,gate,p,
                              NULL,&walk_faces)!=0||walk_faces==0)
                continue;
        }
        size_t na=a1-a0+1,nb=b1-b0+1;
        CzMatchPosition *position=(CzMatchPosition*)malloc(
            (na+nb)*sizeof(*position));
        size_t driver_n=cc->driver?nb:na;
        size_t *order=(size_t*)malloc((driver_n?driver_n:1)*sizeof(*order));
        size_t *tails=(size_t*)malloc((driver_n?driver_n:1)*sizeof(*tails));
        size_t norder=0;
        size_t nposition=0,amin=na,amax=0,bmin=nb,bmax=0;
        double max2=(double)p->max_cross_edge*p->max_cross_edge;
        if(!position||!order||!tails){free(tails);free(order);free(position);return -1;}
        for(size_t i=0;i<na;i++){
            int32_t av=a->verts[a0+i];
            size_t best=nb;double best_d=DBL_MAX;
            for(size_t j=0;j<nb;j++){
                int32_t bv=b->verts[b0+j];
                double d=cz_dist2(verts,av,bv);
                if(d<=max2&&d<best_d&&cz_match_compatible(
                        verts,normal,av,bv,gate,p)){
                    best=j;best_d=d;
                }
            }
            if(best<nb){
                position[nposition++]=(CzMatchPosition){i,best};
                if(cc->driver==0)order[norder++]=best;
            }
        }
        for(size_t j=0;j<nb;j++){
            int32_t bv=b->verts[b0+j];
            size_t best=na;double best_d=DBL_MAX;
            for(size_t i=0;i<na;i++){
                int32_t av=a->verts[a0+i];
                double d=cz_dist2(verts,av,bv);
                if(d<=max2&&d<best_d&&cz_match_compatible(
                        verts,normal,av,bv,gate,p)){
                    best=i;best_d=d;
                }
            }
            if(best<na){
                if(cc->driver==1)order[norder++]=best;
                int duplicate=0;
                for(size_t k=0;k<nposition;k++)
                    if(position[k].a==best&&position[k].b==j){
                        duplicate=1;break;
                    }
                if(!duplicate)
                    position[nposition++]=(CzMatchPosition){best,j};
            }
        }
        double order_coverage=cz_order_coverage(order,norder,tails);
        free(tails);free(order);
        if(nposition<2){free(position);continue;}
        qsort(position,nposition,sizeof(*position),cz_match_position_cmp);
        for(size_t k=0;k<nposition;k++){
            if(position[k].a<amin)amin=position[k].a;
            if(position[k].a>amax)amax=position[k].a;
            if(position[k].b<bmin)bmin=position[k].b;
            if(position[k].b>bmax)bmax=position[k].b;
        }
        if(amin>=amax||bmin>=bmax){free(position);continue;}
        size_t first=nsample;
        for(size_t k=0;k<nposition;k++){
            sample[nsample].vertex_a=a->verts[a0+position[k].a];
            sample[nsample].vertex_b=b->verts[b0+position[k].b];
            nsample++;
        }
        free(position);
        if(nsample-first<2){
            nsample=first;
            continue;
        }
        match[nmatch].chart_a=a->chart;
        match[nmatch].chart_b=b->chart;
        match[nmatch].port_a=cc->chain_a;
        match[nmatch].port_b=cc->chain_b;
        match[nmatch].interval_a_first=a0+amin;
        match[nmatch].interval_a_last=a0+amax;
        match[nmatch].interval_b_first=b0+bmin;
        match[nmatch].interval_b_last=b0+bmax;
        match[nmatch].sample_first=first;
        match[nmatch].sample_count=nsample-first;
        match[nmatch].full_coverage=cc->full_coverage;
        match[nmatch].span_coverage=cc->span_coverage;
        match[nmatch].mean_gap=cc->mean_gap;
        match[nmatch].support=cc->support;
        match[nmatch].score=cc->score;
        match[nmatch].order_coverage=order_coverage;
        nmatch++;
    }
    *out_match=match;*out_nmatch=nmatch;
    *out_sample=sample;*out_nsample=nsample;
    return 0;
}

static double cz_triangle_cost(const float *verts,int32_t a,int32_t b,int32_t c,
                               const double refn[3],
                               const ChartZipperParams *p)
{
    const float *pa=&verts[(size_t)a*3],*pb=&verts[(size_t)b*3],*pc=&verts[(size_t)c*3];
    double u[3]={pb[0]-pa[0],pb[1]-pa[1],pb[2]-pa[2]};
    double v[3]={pc[0]-pa[0],pc[1]-pa[1],pc[2]-pa[2]};
    double n[3]={u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
    double area2=sqrt(cz_dot(n,n));
    if(area2<1e-10) return CZ_INF;
    double l0=cz_dist2(verts,a,b),l1=cz_dist2(verts,b,c),l2=cz_dist2(verts,c,a);
    double lm=l0>l1?(l0>l2?l0:l2):(l1>l2?l1:l2);
    if (lm<=1e-20 ||
        area2/sqrt(lm)<(double)p->min_triangle_altitude)
        return CZ_INF;
    double align=0.0,rm=sqrt(cz_dot(refn,refn));
    if(rm>1e-12) align=1.0-fabs(cz_dot(n,refn))/(area2*rm);
    return sqrt(lm)+(double)p->skinny_weight*lm/(area2+1e-12)+0.5*align;
}

static int cz_segment_edges_free(const CzChain *c,size_t i0,size_t i1,
                                 const uint8_t *edge_used)
{
    if(i0>i1){size_t t=i0;i0=i1;i1=t;}
    for(size_t i=i0;i<i1;i++) if(edge_used[c->edges[i]]) return 0;
    return 1;
}

static void cz_mark_segment(const CzChain *c,size_t i0,size_t i1,uint8_t *edge_used)
{
    if(i0>i1){size_t t=i0;i0=i1;i1=t;}
    for(size_t i=i0;i<i1;i++) edge_used[c->edges[i]]=1;
}

static int cz_pair_already_used(const CzChain *a,const CzChain *b,
                                const int64_t *pairs,size_t npairs)
{
    int32_t lo=a->chart<b->chart?a->chart:b->chart;
    int32_t hi=a->chart<b->chart?b->chart:a->chart;
    int64_t key=((int64_t)lo<<32)|(uint32_t)hi;
    for(size_t i=0;i<npairs;i++) if(pairs[i]==key) return 1;
    return 0;
}

static int cz_zip_dp(const float *verts,const double *normal,
                     const CzChain *a,size_t a0,size_t a1,
                     const CzChain *b,size_t b0,size_t b1,
                     const BpaBridgeGate *gate,const ChartZipperParams *p,
                     int32_t *out,size_t *out_nf)
{
    (void)normal;
    size_t na=a1-a0+1,nb=b1-b0+1;
    int reverse_b=0;
    double same=cz_dist2(verts,a->verts[a0],b->verts[b0])+cz_dist2(verts,a->verts[a1],b->verts[b1]);
    double cross=cz_dist2(verts,a->verts[a0],b->verts[b1])+cz_dist2(verts,a->verts[a1],b->verts[b0]);
    if(cross<same) reverse_b=1;
    double *prevrow=NULL,*row=NULL;
    uint8_t *step=NULL;
    int32_t *temp=NULL;
    int rc=-1;
    prevrow=(double*)malloc(nb*sizeof(*prevrow));
    row=(double*)malloc(nb*sizeof(*row));
    if(out){
        step=(uint8_t*)calloc(na*nb,1);
        temp=(int32_t*)malloc((na+nb-2)*3*sizeof(*temp));
    }
    if(!out_nf||!prevrow||!row||(out&&(!step||!temp))) goto done;
    for(size_t j=0;j<nb;j++) prevrow[j]=CZ_INF;
    int32_t av0=a->verts[a0];
    int32_t bv0=b->verts[reverse_b?b1:b0];
    if(cz_dist2(verts,av0,bv0)>(double)p->max_cross_edge*p->max_cross_edge||
       !cz_gate_pair(gate,verts,av0,bv0)) goto done;
    prevrow[0]=sqrt(cz_dist2(verts,av0,bv0));
    double refn[3]={a->normal[0]+b->normal[0],a->normal[1]+b->normal[1],a->normal[2]+b->normal[2]};
    for(size_t i=0;i<na;i++){
        if(i>0) for(size_t j=0;j<nb;j++) row[j]=CZ_INF;
        double *cur=i==0?prevrow:row;
        for(size_t j=0;j<nb;j++){
            if(i==0&&j==0) continue;
            int32_t av=a->verts[a0+i];
            size_t bj=reverse_b?b1-j:b0+j;
            int32_t bv=b->verts[bj];
            if(cz_dist2(verts,av,bv)>(double)p->max_cross_edge*p->max_cross_edge||
               !cz_gate_pair(gate,verts,av,bv)) continue;
            double best=CZ_INF;uint8_t dir=0;
            if(i>0&&prevrow[j]<CZ_INF){
                int32_t ap=a->verts[a0+i-1];
                double q=cz_triangle_cost(verts,av,ap,bv,refn,p);
                if(q<CZ_INF&&prevrow[j]+q<best){best=prevrow[j]+q;dir=1;}
            }
            if(j>0&&cur[j-1]<CZ_INF){
                size_t bpj=reverse_b?b1-(j-1):b0+(j-1);
                int32_t bp=b->verts[bpj];
                double q=cz_triangle_cost(verts,av,bp,bv,refn,p);
                if(q<CZ_INF&&cur[j-1]+q<best){best=cur[j-1]+q;dir=2;}
            }
            cur[j]=best;if(step)step[i*nb+j]=dir;
        }
        if(i>0){double *t=prevrow;prevrow=row;row=t;}
    }
    if(prevrow[nb-1]>=CZ_INF) goto done;
    if(!out){
        *out_nf=na+nb-2;
    }else{
        size_t i=na-1,j=nb-1,nf=0;
        while(i>0||j>0){
            uint8_t dir=step[i*nb+j];
            int32_t av=a->verts[a0+i];
            size_t bj=reverse_b?b1-j:b0+j;
            int32_t bv=b->verts[bj];
            if(dir==1){
                int32_t ap=a->verts[a0+i-1];
                temp[nf*3]=av;temp[nf*3+1]=ap;temp[nf*3+2]=bv;i--;
            }else if(dir==2){
                size_t bpj=reverse_b?b1-(j-1):b0+(j-1);
                int32_t bp=b->verts[bpj];
                temp[nf*3]=av;temp[nf*3+1]=bp;temp[nf*3+2]=bv;j--;
            }else goto done;
            nf++;
        }
        for(size_t f=0;f<nf;f++) memcpy(&out[f*3],&temp[(nf-1-f)*3],3*sizeof(int32_t));
        *out_nf=nf;
    }
    rc=0;
done:
    free(temp);free(step);free(row);free(prevrow);
    return rc;
}

/* A valid chart join needs one supported boundary interval, not necessarily
 * every seam-facing edge of a curled micro-chart.  If the full voted interval
 * cannot be triangulated inside the hard edge gate, trim only the longer arc
 * and take the longest geometrically valid sub-interval.  This is contour
 * selection, not mesh repair: no source face or vertex changes. */
static int cz_zip_trimmed(const float *verts,const double *normal,
                          const CzChain *a,size_t *a0,size_t *a1,
                          const CzChain *b,size_t *b0,size_t *b1,
                          const BpaBridgeGate *gate,
                          const ChartZipperParams *p,
                          int32_t *out,size_t *out_nf)
{
    size_t na=*a1-*a0+1,nb=*b1-*b0+1;
    int trim_a=na>nb;
    size_t long_n=trim_a?na:nb;
    int32_t *trial=NULL,*best_faces=NULL;
    int rc=-1;
    size_t cap=na+nb-2;
    if(out){
        trial=(int32_t*)malloc((cap?cap:1)*3*sizeof(*trial));
        best_faces=(int32_t*)malloc((cap?cap:1)*3*sizeof(*best_faces));
    }
    if(!out_nf||(out&&(!trial||!best_faces)))goto done;
    if(long_n<=2)goto done;
    for(size_t win=long_n-1;win>=2;win--){
        double best_cost=DBL_MAX;size_t best_start=0,best_nf=0;
        for(size_t start=0;start+win<=long_n;start++){
            size_t qa0=*a0,qa1=*a1,qb0=*b0,qb1=*b1;
            if(trim_a){qa0=*a0+start;qa1=qa0+win-1;}
            else{qb0=*b0+start;qb1=qb0+win-1;}
            double la=cz_arc_length(verts,a,qa0,qa1);
            double lb=cz_arc_length(verts,b,qb0,qb1);
            double support=la<lb?la:lb;
            if(support<(double)p->min_support_length)continue;
            double ec=cz_interval_endpoint_cost(verts,a,qa0,qa1,
                                                b,qb0,qb1,gate,p);
            if(ec>=DBL_MAX)continue;
            size_t tnf=0;
            if(cz_zip_dp(verts,normal,a,qa0,qa1,b,qb0,qb1,
                         gate,p,out?trial:NULL,&tnf)!=0||tnf==0)continue;
            if(ec<best_cost){
                best_cost=ec;best_start=start;best_nf=tnf;
                if(out)memcpy(best_faces,trial,tnf*3*sizeof(*trial));
            }
        }
        if(best_nf>0){
            if(trim_a){*a0+=best_start;*a1=*a0+win-1;}
            else{*b0+=best_start;*b1=*b0+win-1;}
            if(out)memcpy(out,best_faces,best_nf*3*sizeof(*out));
            *out_nf=best_nf;rc=0;goto done;
        }
        if(win==2)break;
    }
done:
    free(best_faces);free(trial);
    return rc;
}

static int cz_process_transactions(Arena_T arena,
                         const float *verts,size_t nv,
                        const int32_t *faces,size_t nf,
                        const int32_t *chart_component,
                        size_t n_chart_components,
                        const BpaBridgeGate *gate,
                        const ChartZipperParams *params_in,
                         int32_t **out_faces,size_t *out_nf,
                         size_t *out_n_bridge,
                         ChartZipperStats *stats,
                         ChartZipperTransaction **out_transactions,
                         size_t *out_n_transactions,
                         int mode,
                         ChartZipperMatch **out_matches,
                         size_t *out_n_matches,
                         ChartZipperMatchSample **out_samples,
                         size_t *out_n_samples)
{
    ChartZipperParams defaults;
    const ChartZipperParams *p=params_in;
    CzHalfEdge *he=NULL;
    CzBoundaryEdge *be=NULL;
    double *normal=NULL;
    uint8_t *used=NULL,*edge_used=NULL;
    SeamPlane planes[64];size_t np=0,nbe=0;
    CzChainVec chains={0};CzCandidateVec candidates={0};
    int32_t *result=NULL,*ziptemp=NULL;
    size_t result_capacity=0;
    int64_t *accepted_pairs=NULL;size_t naccepted_pairs=0;
    CzAcceptedTx *accepted_tx=NULL;size_t naccepted_tx=0;
    size_t *conflict_degree=NULL;
    UnionFind chart_forest;
    int enumerate_all=mode==1;
    int match_only=mode==2;
    int rc=-1;
    if(stats)memset(stats,0,sizeof(*stats));
    if(out_transactions)*out_transactions=NULL;
    if(out_n_transactions)*out_n_transactions=0;
    if(out_matches)*out_matches=NULL;
    if(out_n_matches)*out_n_matches=0;
    if(out_samples)*out_samples=NULL;
    if(out_n_samples)*out_n_samples=0;
    if(out_faces)*out_faces=NULL;if(out_nf)*out_nf=0;if(out_n_bridge)*out_n_bridge=0;
    if(!arena||!verts||!faces||!chart_component||!out_faces||!out_nf||
       !out_n_bridge||((out_transactions==NULL)!=(out_n_transactions==NULL))||
       ((out_matches==NULL)!=(out_n_matches==NULL))||
       ((out_samples==NULL)!=(out_n_samples==NULL))||
       (match_only&&(!out_matches||!out_samples))||
       n_chart_components==0||nv>(size_t)INT32_MAX||
       n_chart_components>(size_t)INT32_MAX) return -1;
    if(!p){ChartZipper_default_params(&defaults);p=&defaults;}
    used=(uint8_t*)calloc(nv?nv:1,1);
    normal=(double*)calloc((nv?nv:1)*3,sizeof(*normal));
    he=(CzHalfEdge*)malloc((nf?nf*3:1)*sizeof(*he));
    if(!used||!normal||!he) goto done;
    cz_compute_normals(verts,nv,faces,nf,normal);
    for(size_t f=0;f<nf;f++){
        int32_t q[3]={faces[f*3],faces[f*3+1],faces[f*3+2]};
        for(int k=0;k<3;k++){
            int32_t a=q[k],b=q[(k+1)%3],o=q[(k+2)%3];
            if(a<0||b<0||o<0||(size_t)a>=nv||(size_t)b>=nv||(size_t)o>=nv) goto done;
            used[a]=used[b]=used[o]=1;
            CzHalfEdge *h=&he[f*3+(size_t)k];
            h->lo=a<b?a:b;h->hi=a<b?b:a;h->a=a;h->b=b;h->opposite=o;h->face=(int32_t)f;
        }
    }
    np=SeamPlanes_detect(verts,nv,used,p->cube_size,p->band,planes,64);
    if(stats)stats->planes=np;
    qsort(he,nf*3,sizeof(*he),cz_halfedge_cmp);
    for(size_t i=0;i<nf*3;){size_t j=i+1;while(j<nf*3&&he[j].lo==he[i].lo&&he[j].hi==he[i].hi)j++;if(j-i==1)nbe++;i=j;}
    be=(CzBoundaryEdge*)malloc((nbe?nbe:1)*sizeof(*be));
    if(!be)goto done;
    nbe=0;
    for(size_t i=0;i<nf*3;){
        size_t j=i+1;while(j<nf*3&&he[j].lo==he[i].lo&&he[j].hi==he[i].hi)j++;
        if(j-i==1){
            int32_t ca=chart_component[he[i].a],cb=chart_component[he[i].b];
            if(ca>=0&&ca==cb&&(size_t)ca<n_chart_components){
                int plane=-1,side=0;
                if(cz_best_plane(verts,&he[i],planes,np,p->band,&plane,&side)){
                    be[nbe].a=he[i].a;be[nbe].b=he[i].b;be[nbe].opposite=he[i].opposite;
                    be[nbe].chart=ca;be[nbe].plane=plane;be[nbe].side=(uint8_t)side;nbe++;
                }
            }
        }
        i=j;
    }
    if(stats){stats->boundary_edges=0;for(size_t i=0;i<nf*3;){size_t j=i+1;while(j<nf*3&&he[j].lo==he[i].lo&&he[j].hi==he[i].hi)j++;if(j-i==1)stats->boundary_edges++;i=j;}stats->seam_boundary_edges=nbe;}
    for(size_t plane=0;plane<np;plane++)for(int side=0;side<2;side++)
        if(cz_extract_chains(verts,nv,normal,be,nbe,plane,side,&chains,stats)!=0)goto done;
    /* Resolve only closed contours that already have sustained compatible
     * support on an ordinary opposite-side path.  Build the seed indexes
     * before opening either side so two closed islands cannot nominate one
     * another into existence. */
    for(size_t plane=0;plane<np;plane++){
        CzPointIndex low_seed,high_seed;
        if(cz_point_index(arena,verts,chains.v,chains.n,(int32_t)plane,0,
                          &low_seed)!=0||
           cz_point_index(arena,verts,chains.v,chains.n,(int32_t)plane,1,
                          &high_seed)!=0)goto done;
        if(cz_open_supported_closed_chains(
                verts,normal,chains.v,chains.n,(int32_t)plane,0,&high_seed,
                gate,p,be,stats)!=0||
           cz_open_supported_closed_chains(
                verts,normal,chains.v,chains.n,(int32_t)plane,1,&low_seed,
                gate,p,be,stats)!=0)goto done;
    }
    if(p->trace){
        for(size_t ci=0;ci<chains.n;ci++){
            const CzChain *c=&chains.v[ci];
            double mn[3]={DBL_MAX,DBL_MAX,DBL_MAX};
            double mx[3]={-DBL_MAX,-DBL_MAX,-DBL_MAX};
            for(size_t j=0;j<c->n;j++) for(int d=0;d<3;d++){
                double q=verts[(size_t)c->verts[j]*3+(size_t)d];
                if(q<mn[d])mn[d]=q;if(q>mx[d])mx[d]=q;
            }
            fprintf(stderr,
                "    [zip trace] path=%zu chart=%d plane=%d axis=%d "
                "coord=%.3f side=%u verts=%zu length=%.3f closed=%u "
                "bbox=[%.3f %.3f %.3f]-[%.3f %.3f %.3f]\n",
                ci,c->chart,c->plane,planes[c->plane].axis,
                planes[c->plane].coord,(unsigned)c->side,c->n,c->length,
                (unsigned)c->closed,mn[0],mn[1],mn[2],mx[0],mx[1],mx[2]);
        }
    }
    for(size_t plane=0;plane<np;plane++){
        CzPointIndex low,high;
        if(cz_point_index(arena,verts,chains.v,chains.n,(int32_t)plane,0,&low)!=0||
           cz_point_index(arena,verts,chains.v,chains.n,(int32_t)plane,1,&high)!=0)goto done;
        if(cz_generate_direction(arena,verts,normal,chains.v,chains.n,(int32_t)plane,0,&high,gate,p,&candidates,stats)!=0||
           cz_generate_direction(arena,verts,normal,chains.v,chains.n,(int32_t)plane,1,&low,gate,p,&candidates,stats)!=0)goto done;
    }
    if(stats)stats->candidates=candidates.n;
    qsort(candidates.v,candidates.n,sizeof(*candidates.v),cz_candidate_cmp);
    if(match_only){
        if(cz_emit_matches(arena,verts,normal,chains.v,candidates.v,
                           candidates.n,gate,p,out_matches,out_n_matches,
                           out_samples,out_n_samples)!=0)
            goto done;
        rc=0;
        goto done;
    }
    if(!enumerate_all){
        edge_used=(uint8_t*)calloc(nbe?nbe:1,1);
        accepted_pairs=(int64_t*)malloc(
            (candidates.n?candidates.n:1)*sizeof(*accepted_pairs));
    }
    accepted_tx=(CzAcceptedTx*)malloc((candidates.n?candidates.n:1)*sizeof(*accepted_tx));
    ziptemp=(int32_t*)malloc((nbe?nbe:1)*3*sizeof(*ziptemp));
    if(cz_face_reserve(&result,&result_capacity,nf)!=0||
       (!enumerate_all&&(!edge_used||!accepted_pairs))||
       !accepted_tx||!ziptemp)goto done;
    if(!enumerate_all)chart_forest=UF_new(arena,(int32_t)n_chart_components);
    memcpy(result,faces,nf*3*sizeof(int32_t));
    size_t write=nf;
    for(size_t ci=0;ci<candidates.n;ci++){
        CzCandidate *cc=&candidates.v[ci];
        CzChain *a=&chains.v[cc->chain_a],*b=&chains.v[cc->chain_b];
        if(!enumerate_all&&
           cz_pair_already_used(a,b,accepted_pairs,naccepted_pairs)){
            if(stats)stats->duplicate_pair_candidates++;
            if(p->trace)fprintf(stderr,
                "    [zip trace] reject duplicate charts=%d:%d paths=%d:%d\n",
                a->chart,b->chart,cc->chain_a,cc->chain_b);
            continue;
        }
        /* Gluing two disk components along one boundary interval preserves a
         * disk.  Joining charts already connected by earlier zippers would
         * close a graph cycle and can create an annulus or handle, so reject
         * that transaction before touching the face array. */
        if(!enumerate_all&&
           uf_find(&chart_forest,a->chart)==uf_find(&chart_forest,b->chart)){
            if(stats)stats->cycle_candidates++;
            if(p->trace)fprintf(stderr,
                "    [zip trace] reject cycle charts=%d:%d paths=%d:%d\n",
                a->chart,b->chart,cc->chain_a,cc->chain_b);
            continue;
        }
        if(!enumerate_all&&
           (!cz_segment_edges_free(a,cc->a0,cc->a1,edge_used)||
            !cz_segment_edges_free(b,cc->b0,cc->b1,edge_used))){
            if(stats)stats->overlap_candidates++;
            if(p->trace)fprintf(stderr,
                "    [zip trace] reject overlap charts=%d:%d paths=%d:%d\n",
                a->chart,b->chart,cc->chain_a,cc->chain_b);
            continue;
        }
        size_t znf=0;
        size_t za0=cc->a0,za1=cc->a1,zb0=cc->b0,zb1=cc->b1;
        if(cz_zip_dp(verts,normal,a,cc->a0,cc->a1,b,cc->b0,cc->b1,
                     gate,p,ziptemp,&znf)!=0||znf==0){
            znf=0;
            if(cz_zip_trimmed(verts,normal,a,&za0,&za1,b,&zb0,&zb1,
                              gate,p,ziptemp,&znf)!=0||znf==0){
                if(stats)stats->geometry_candidates++;
                if(p->trace)fprintf(stderr,
                    "    [zip trace] reject geometry charts=%d:%d paths=%d:%d "
                    "span=%zu:%zu/%zu:%zu\n",
                    a->chart,b->chart,cc->chain_a,cc->chain_b,
                    cc->a0,cc->a1,cc->b0,cc->b1);
                continue;
            }
            if(stats)stats->trimmed_candidates++;
            if(p->trace)fprintf(stderr,
                "    [zip trace] trim charts=%d:%d paths=%d:%d "
                "span=%zu:%zu/%zu:%zu -> %zu:%zu/%zu:%zu\n",
                a->chart,b->chart,cc->chain_a,cc->chain_b,
                cc->a0,cc->a1,cc->b0,cc->b1,za0,za1,zb0,zb1);
        }
        accepted_tx[naccepted_tx].candidate=ci;
        accepted_tx[naccepted_tx].face0=write;
        accepted_tx[naccepted_tx].nf=znf;
        naccepted_tx++;
        if(znf>SIZE_MAX-write||
           cz_face_reserve(&result,&result_capacity,write+znf)!=0)goto done;
        memcpy(&result[write*3],ziptemp,znf*3*sizeof(int32_t));write+=znf;
        if(!enumerate_all){
            cz_mark_segment(a,za0,za1,edge_used);
            cz_mark_segment(b,zb0,zb1,edge_used);
            int32_t lo=a->chart<b->chart?a->chart:b->chart;
            int32_t hi=a->chart<b->chart?b->chart:a->chart;
            accepted_pairs[naccepted_pairs++]=((int64_t)lo<<32)|(uint32_t)hi;
            uf_union(&chart_forest,a->chart,b->chart);
        }
        if(stats){stats->accepted_transactions++;stats->bridge_faces+=znf;}
        if(p->trace)fprintf(stderr,
            "    [zip trace] ACCEPT charts=%d:%d paths=%d:%d faces=%zu\n",
            a->chart,b->chart,cc->chain_a,cc->chain_b,znf);
    }
    /*
     * Geometry is an atomic transaction too. Audit the provisional source +
     * zipper complex with the same robust triangle predicate used by the
     * standalone seam diagnostic. If any face of a zipper strip overlaps,
     * stabs, or folds into existing geometry, reject the WHOLE strip. Source
     * faces remain immutable; no cleanup face deletion or vertex cut occurs.
     */
    if(!enumerate_all&&naccepted_tx>0){
        IntersectionCleanupParams ip;
        IntersectionCleanupStats is;
        size_t out_write=nf,kept_tx=0,kept_faces=0;
        conflict_degree=(size_t*)calloc(write?write:1,sizeof(*conflict_degree));
        if(!conflict_degree)goto done;
        IntersectionCleanup_default_params(&ip);
        ip.gap_max=(double)p->conflict_gap;
        ip.parallel_angle_deg=(double)p->conflict_parallel_angle_deg;
        ip.include_hinges=1;
        if(IntersectionCleanup_audit(verts,nv,result,write,NULL,&ip,
                                     conflict_degree,&is)!=0)
            goto done;
        if(stats)stats->conflict_pairs=is.conflicts;
        for(size_t t=0;t<naccepted_tx;t++){
            const CzAcceptedTx *tx=&accepted_tx[t];
            int bad=0;
            for(size_t f=tx->face0;f<tx->face0+tx->nf;f++)
                if(conflict_degree[f]){bad=1;break;}
            if(bad){
                if(stats)stats->conflict_transactions++;
                if(p->trace)fprintf(stderr,
                    "    [zip trace] reject conflict candidate=%zu "
                    "faces=%zu range=%zu:%zu\n",
                    tx->candidate,tx->nf,tx->face0,tx->face0+tx->nf);
                continue;
            }
            if(out_write!=tx->face0)
                memmove(&result[out_write*3],&result[tx->face0*3],
                        tx->nf*3*sizeof(int32_t));
            accepted_tx[kept_tx]=*tx;
            accepted_tx[kept_tx].face0=out_write;
            out_write+=tx->nf;kept_faces+=tx->nf;kept_tx++;
        }
        naccepted_tx=kept_tx;
        if(stats){
            stats->accepted_transactions=kept_tx;
            stats->bridge_faces=kept_faces;
        }
        write=out_write;
        /*
         * Certify the COMPACTED transaction too.  The first audit addresses all
         * provisional strips at their original ranges; this second pass is the
         * authoritative output-layout check and guards range-compaction and
         * subsequent face-order regressions.  Any residual still removes its
         * entire strip atomically.  A final zero-conflict audit then commits.
         */
        if(naccepted_tx>0){
            if(IntersectionCleanup_audit(verts,nv,result,write,NULL,&ip,
                                         conflict_degree,&is)!=0)
                goto done;
            if(stats){
                stats->postcompact_conflict_pairs=is.conflicts;
                stats->conflict_pairs+=is.conflicts;
            }
            if(is.conflicts>0){
                size_t compact_write=nf,compact_kept=0,compact_faces=0;
                for(size_t t=0;t<naccepted_tx;t++){
                    CzAcceptedTx tx=accepted_tx[t];
                    int bad=0;
                    for(size_t f=tx.face0;f<tx.face0+tx.nf;f++)
                        if(conflict_degree[f]){bad=1;break;}
                    if(bad){
                        if(stats){
                            stats->conflict_transactions++;
                            stats->postcompact_conflict_transactions++;
                        }
                        if(p->trace)fprintf(stderr,
                            "    [zip trace] reject postcompact conflict "
                            "candidate=%zu faces=%zu range=%zu:%zu\n",
                            tx.candidate,tx.nf,tx.face0,tx.face0+tx.nf);
                        continue;
                    }
                    if(compact_write!=tx.face0)
                        memmove(&result[compact_write*3],
                                &result[tx.face0*3],
                                tx.nf*3*sizeof(int32_t));
                    accepted_tx[compact_kept]=tx;
                    accepted_tx[compact_kept].face0=compact_write;
                    compact_write+=tx.nf;
                    compact_faces+=tx.nf;
                    compact_kept++;
                }
                if(compact_kept==naccepted_tx)goto done;
                naccepted_tx=compact_kept;
                write=compact_write;
                if(stats){
                    stats->accepted_transactions=compact_kept;
                    stats->bridge_faces=compact_faces;
                }
                if(IntersectionCleanup_audit(
                        verts,nv,result,write,NULL,&ip,
                        conflict_degree,&is)!=0||is.conflicts!=0)
                    goto done;
            }
        }
        if(stats)stats->embedded_certificate=1;
    }else if(!enumerate_all&&stats){
        stats->embedded_certificate=1;
    }
    if(out_transactions){
        ChartZipperTransaction *tx=(ChartZipperTransaction*)ARENA_ALLOC(
            arena,(naccepted_tx?naccepted_tx:1)*sizeof(*tx));
        for(size_t t=0;t<naccepted_tx;t++){
            const CzAcceptedTx *at=&accepted_tx[t];
            const CzCandidate *cc=&candidates.v[at->candidate];
            const CzChain *a=&chains.v[cc->chain_a];
            const CzChain *b=&chains.v[cc->chain_b];
            tx[t].chart_a=a->chart;tx[t].chart_b=b->chart;
            tx[t].face_first=at->face0;tx[t].face_count=at->nf;
            tx[t].orient_a=0;tx[t].orient_b=0;
            tx[t].flags=enumerate_all?0:CHART_ZIPPER_TRANSACTION_PHYSICAL;
            tx[t].full_coverage=cc->full_coverage;
            tx[t].span_coverage=cc->span_coverage;
            tx[t].mean_gap=cc->mean_gap;
            tx[t].support=cc->support;tx[t].score=cc->score;
        }
        *out_transactions=tx;*out_n_transactions=naccepted_tx;
    }
    {
        int32_t *arena_result=(int32_t*)ARENA_ALLOC(
            arena,(write?write:1)*3*sizeof(*arena_result));
        memcpy(arena_result,result,write*3*sizeof(*arena_result));
        *out_faces=arena_result;
    }
    *out_nf=write;*out_n_bridge=write-nf;rc=0;
done:
    if(rc!=0&&out_faces)*out_faces=NULL;
    free(conflict_degree);free(accepted_tx);free(accepted_pairs);
    free(result);free(ziptemp);free(edge_used);
    for(size_t i=0;i<chains.n;i++){free(chains.v[i].verts);free(chains.v[i].edges);}
    free(chains.v);free(candidates.v);free(be);free(he);free(normal);free(used);
    return rc;
}

int ChartZipper_process_with_transactions(Arena_T arena,
                         const float *verts,size_t nv,
                         const int32_t *faces,size_t nf,
                         const int32_t *chart_component,
                         size_t n_chart_components,
                         const BpaBridgeGate *gate,
                         const ChartZipperParams *params,
                         int32_t **out_faces,size_t *out_nf,
                         size_t *out_n_bridge,
                         ChartZipperStats *stats,
                         ChartZipperTransaction **out_transactions,
                         size_t *out_n_transactions)
{
    return cz_process_transactions(
        arena,verts,nv,faces,nf,chart_component,n_chart_components,
        gate,params,out_faces,out_nf,out_n_bridge,stats,
        out_transactions,out_n_transactions,0,NULL,NULL,NULL,NULL);
}

int ChartZipper_enumerate_transactions(Arena_T arena,
                         const float *verts,size_t nv,
                         const int32_t *faces,size_t nf,
                         const int32_t *chart_component,
                         size_t n_chart_components,
                         const BpaBridgeGate *gate,
                         const ChartZipperParams *params,
                         int32_t **out_faces,size_t *out_nf,
                         size_t *out_n_bridge,
                         ChartZipperStats *stats,
                         ChartZipperTransaction **out_transactions,
                         size_t *out_n_transactions)
{
    if(!out_transactions||!out_n_transactions)return -1;
    return cz_process_transactions(
        arena,verts,nv,faces,nf,chart_component,n_chart_components,
        gate,params,out_faces,out_nf,out_n_bridge,stats,
        out_transactions,out_n_transactions,1,NULL,NULL,NULL,NULL);
}

int ChartZipper_match(Arena_T arena,
                         const float *verts,size_t nv,
                         const int32_t *faces,size_t nf,
                         const int32_t *chart_component,
                         size_t n_chart_components,
                         const BpaBridgeGate *gate,
                         const ChartZipperParams *params,
                         ChartZipperMatch **out_matches,size_t *out_n_matches,
                         ChartZipperMatchSample **out_samples,
                         size_t *out_n_samples,
                         ChartZipperStats *stats)
{
    int32_t *unused_faces=NULL;
    size_t unused_nf=0,unused_bridge=0;
    return cz_process_transactions(
        arena,verts,nv,faces,nf,chart_component,n_chart_components,
        gate,params,&unused_faces,&unused_nf,&unused_bridge,stats,
        NULL,NULL,2,out_matches,out_n_matches,out_samples,out_n_samples);
}

int ChartZipper_process(Arena_T arena,
                        const float *verts,size_t nv,
                        const int32_t *faces,size_t nf,
                        const int32_t *chart_component,
                        size_t n_chart_components,
                        const BpaBridgeGate *gate,
                        const ChartZipperParams *params,
                        int32_t **out_faces,size_t *out_nf,
                        size_t *out_n_bridge,
                        ChartZipperStats *stats)
{
    return ChartZipper_process_with_transactions(
        arena,verts,nv,faces,nf,chart_component,n_chart_components,
        gate,params,out_faces,out_nf,out_n_bridge,stats,NULL,NULL);
}

/* --- native fixtures ----------------------------------------------------- */

static int cz_fixture(int two_children,int far_apart)
{
    /* Low chart is a 3x2 grid strip with seam path x=-1, y=0..10.  The high
     * side is either one matching strip or two disjoint small strips sharing
     * the low path.  All charts are disks and all coordinates are (z,y,x). */
    float v[30*3];int32_t f[32*3],chart[30];size_t nv=0,nf=0;
    #define ADDV(z,y,x,ch) (v[nv*3]=(z),v[nv*3+1]=(y),v[nv*3+2]=(x),chart[nv]=(ch),(int32_t)nv++)
    #define ADDF(a,b,c) do{f[nf*3]=(a);f[nf*3+1]=(b);f[nf*3+2]=(c);nf++;}while(0)
    int32_t l0=ADDV(0,0,-4,0),l1=ADDV(0,5,-4,0),l2=ADDV(0,10,-4,0);
    int32_t s0=ADDV(0,0,-1,0),s1=ADDV(0,5,-1,0),s2=ADDV(0,10,-1,0);
    ADDF(l0,s0,l1);ADDF(s0,s1,l1);ADDF(l1,s1,l2);ADDF(s1,s2,l2);
    if(!two_children){
        float hx=far_apart?20.0f:1.0f,ox=far_apart?23.0f:4.0f;
        int32_t h0=ADDV(0,0,hx,1),h1=ADDV(0,5,hx,1),h2=ADDV(0,10,hx,1);
        int32_t r0=ADDV(0,0,ox,1),r1=ADDV(0,5,ox,1),r2=ADDV(0,10,ox,1);
        ADDF(h0,r0,h1);ADDF(r0,r1,h1);ADDF(h1,r1,h2);ADDF(r1,r2,h2);
    }else{
        int32_t h0=ADDV(0,0,1,1),h1=ADDV(0,4,1,1),r0=ADDV(0,0,4,1),r1=ADDV(0,4,4,1);
        ADDF(h0,r0,h1);ADDF(r0,r1,h1);
        int32_t h2=ADDV(0,6,1,2),h3=ADDV(0,10,1,2),r2=ADDV(0,6,4,2),r3=ADDV(0,10,4,2);
        ADDF(h2,r2,h3);ADDF(r2,r3,h3);
    }
    Arena_T a=Arena_new(),ga=NULL;ChartZipperParams p;ChartZipperStats s,gs;
    int32_t *of=NULL;size_t onf=0,nb=0;TopologyAuditReport tr;memset(&tr,0,sizeof tr);
    ChartZipperTransaction *gtx=NULL;size_t gntx=0,gnf=0,gnb=0;
    ChartZipperMatch *gm=NULL;ChartZipperMatchSample *gms=NULL;
    size_t ngm=0,ngms=0;
    int32_t *gf=NULL;
    ChartZipper_default_params(&p);p.cube_size=10.0f;p.band=6.0f;
    /* Refined production paths have short longitudinal edges.  This fixture's
     * five-voxel sampling needs the corresponding sqrt(5^2 + 2^2) zipper
     * diagonal, while the opposing contour itself remains only two voxels
     * away. */
    p.max_cross_edge=6.0f;p.min_support_length=1.0f;
    int rc=ChartZipper_process(a,v,nv,f,nf,chart,two_children?3:2,NULL,&p,&of,&onf,&nb,&s);
    int fail=rc!=0;
    if(!fail&&!far_apart){
        fail|=s.accepted_transactions!=(size_t)(two_children?2:1);
        fail|=TopologyAudit_analyze(v,nv,of,onf,NULL,&tr)!=0||!tr.all_components_are_disks;
    }else if(!fail&&far_apart) fail|=nb!=0;
    ga=Arena_new();
    if(ChartZipper_enumerate_transactions(
            ga,v,nv,f,nf,chart,two_children?3:2,NULL,&p,
            &gf,&gnf,&gnb,&gs,&gtx,&gntx)!=0)
        fail=1;
    else if(far_apart)
        fail|=gntx!=0||gnb!=0;
    else
        fail|=gntx<=s.accepted_transactions||gntx!=gs.accepted_transactions;
    if(ChartZipper_match(ga,v,nv,f,nf,chart,two_children?3:2,NULL,&p,
                         &gm,&ngm,&gms,&ngms,NULL)!=0)
        fail=1;
    else if(far_apart)
        fail|=ngm!=0||ngms!=0;
    else{
        /* Match-only is the transaction graph with its contour-walk
         * certificate retained and every triangle omitted. */
        fail|=ngm!=gntx||ngms<ngm*2;
        for(size_t m=0;m<ngm;m++){
            fail|=gm[m].sample_count<2||
                   gm[m].sample_first+gm[m].sample_count>ngms||
                   gm[m].chart_a==gm[m].chart_b||
                   gm[m].order_coverage<=0.0||gm[m].order_coverage>1.0;
            for(size_t k=gm[m].sample_first;
                k<gm[m].sample_first+gm[m].sample_count;k++)
                fail|=gms[k].vertex_a<0||gms[k].vertex_b<0||
                      (size_t)gms[k].vertex_a>=nv||
                      (size_t)gms[k].vertex_b>=nv;
        }
    }
    fprintf(stderr,"[selftest] chart zipper %s -> %s "
            "(planes=%zu seam_edges=%zu chains=%zu open=%zu closed=%zu "
            "candidates=%zu ambiguous=%zu geometry=%zu tx=%zu graph=%zu "
            "matches=%zu samples=%zu faces=%zu)\n",
            two_children?"micro absorption":(far_apart?"distance gate":"two disks"),
            fail?"FAIL":"ok",s.planes,s.seam_boundary_edges,s.chains,
            s.open_chains,s.closed_chains_rejected,s.candidates,
            s.ambiguous_candidates,s.geometry_candidates,
            s.accepted_transactions,gntx,ngm,ngms,nb);
    TopologyAudit_dispose(&tr);Arena_dispose(&ga);Arena_dispose(&a);
    #undef ADDF
    #undef ADDV
    return fail?1:0;
}

static int cz_trim_fixture(void)
{
    /* A curled tail can make the complete voted contour fail the hard cross-
     * edge gate even though its longest supported sub-interval is a clean
     * zipper.  The transaction must trim the contour interval, never split or
     * alter either source chart. */
    const float verts[] = {
        0, 0,-1,  0, 6,-1,
        0,-10, 1, 0, 0, 1, 0, 1, 1, 0, 2, 1,
        0, 3, 1,  0, 4, 1, 0, 5, 1, 0, 6, 1
    };
    int32_t av[] = {0,1};
    int32_t bv[] = {2,3,4,5,6,7,8,9};
    CzChain a,b;
    ChartZipperParams p;
    int32_t out[8*3];
    size_t nf=0,a0=0,a1=1,b0=0,b1=7;
    memset(&a,0,sizeof a);memset(&b,0,sizeof b);
    a.verts=av;a.n=2;a.normal[0]=1.0;
    b.verts=bv;b.n=8;b.normal[0]=1.0;
    ChartZipper_default_params(&p);
    p.max_cross_edge=4.0f;p.min_support_length=1.0f;
    int full_rc=cz_zip_dp(verts,NULL,&a,a0,a1,&b,b0,b1,NULL,&p,out,&nf);
    nf=0;
    int trim_rc=cz_zip_trimmed(verts,NULL,&a,&a0,&a1,&b,&b0,&b1,
                               NULL,&p,out,&nf);
    int fail=full_rc==0||trim_rc!=0||nf==0||a0!=0||a1!=1||b0!=1||b1!=7;
    fprintf(stderr,"[selftest] chart zipper curled interval trim -> %s "
            "(full=%d trim=%d span=%zu:%zu/%zu:%zu faces=%zu)\n",
            fail?"FAIL":"ok",full_rc,trim_rc,a0,a1,b0,b1,nf);
    return fail?1:0;
}

static int cz_winding_gate_fixture(void)
{
    BpaBridgeGate gate;
    float v[4*3];
    const double eps=0.01;
    memset(&gate,0,sizeof gate);
    gate.umb_y=1.0; /* arm the constant-axis form without special sentinels */
    gate.pitch=10.0;
    gate.hard=0.40;

    /* Same Archimedean wrap on opposite sides of atan2's branch cut. */
    {
        double t0=M_PI-eps,t1=-M_PI+eps;
        double r0=100.0+gate.pitch*t0/(2.0*M_PI);
        double r1=100.0+gate.pitch*(t0+2.0*eps)/(2.0*M_PI);
        v[0]=0;v[1]=(float)(gate.umb_y+r0*sin(t0));v[2]=(float)(r0*cos(t0));
        v[3]=0;v[4]=(float)(gate.umb_y+r1*sin(t1));v[5]=(float)(r1*cos(t1));
    }
    /* One radial pitch at the same angle is a real adjacent-wrap jump. */
    v[6]=0;v[7]=(float)(gate.umb_y+100.0);v[8]=0;
    v[9]=0;v[10]=(float)(gate.umb_y+110.0);v[11]=0;
    {
        int branch_ok=cz_gate_pair(&gate,v,0,1);
        int jump_ok=cz_gate_pair(&gate,v,2,3);
        int fail=!branch_ok||jump_ok;
        fprintf(stderr,
                "[selftest] chart zipper branch-free winding gate -> %s "
                "(branch=%s radial-jump=%s)\n",
                fail?"FAIL":"ok",branch_ok?"accept":"reject",
                jump_ok?"accept":"reject");
        return fail?1:0;
    }
}

static int cz_closed_contour_fixture(void)
{
    /* The low disk lies wholly inside the x=0 seam band, so all four of its
     * boundary edges form one closed contour.  Only its near side is supported
     * by the ordinary open path of the high disk.  The far edge must be omitted
     * from contour bookkeeping and the supported side zipped without changing
     * either source disk. */
    const float v[]={
         0,0,-3,   0,4,-3,  -1,0,-3.5f, -1,4,-3.5f,
         0,0, 3,   0,4, 3,  -5,0, 5.5f, -5,4, 5.5f
    };
    const int32_t f[]={
        0,1,2, 1,3,2,
        4,6,5, 6,7,5
    };
    const int32_t chart[]={0,0,0,0,1,1,1,1};
    Arena_T a=Arena_new();
    ChartZipperParams p;ChartZipperStats s;
    TopologyAuditReport tr;int32_t *of=NULL;
    size_t onf=0,nb=0;int fail=0;
    memset(&tr,0,sizeof tr);
    ChartZipper_default_params(&p);
    p.cube_size=10.0f;p.band=6.0f;p.max_cross_edge=8.0f;
    p.min_support_length=1.0f;
    if(ChartZipper_process(a,v,8,f,4,chart,2,NULL,&p,
                           &of,&onf,&nb,&s)!=0)
        fail=1;
    if(!fail){
        fail|=s.closed_chains_opened!=1;
        fail|=s.accepted_transactions!=1||nb==0;
        fail|=TopologyAudit_analyze(v,8,of,onf,NULL,&tr)!=0;
        fail|=!tr.all_components_are_disks||tr.face_components!=1;
    }
    fprintf(stderr,
        "[selftest] chart zipper closed-contour continuation -> %s "
        "(planes=%zu chains=%zu open=%zu opened=%zu rejected=%zu "
        "tx=%zu faces=%zu disks=%zu/%zu)\n",
        fail?"FAIL":"ok",s.planes,s.chains,s.open_chains,
        s.closed_chains_opened,
        s.closed_chains_rejected,s.accepted_transactions,nb,
        tr.disk_components,tr.face_components);
    TopologyAudit_dispose(&tr);Arena_dispose(&a);
    return fail?1:0;
}

int ChartZipper_selftest(void)
{
    int fail=0;
    fail+=cz_fixture(0,0);
    fail+=cz_fixture(0,1);
    fail+=cz_fixture(1,0);
    fail+=cz_trim_fixture();
    fail+=cz_winding_gate_fixture();
    fail+=cz_closed_contour_fixture();
    return fail;
}
