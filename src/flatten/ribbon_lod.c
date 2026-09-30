#include "ribbon_lod.h"
#include "../common/ves_platform.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Work is bounded per dyadic patch. The only global arrays are one compact
 * face reference per input triangle and the eventual emitted mesh. */
typedef struct {
    int32_t bx, by, x, y, face;
} RlFace;

typedef struct {
    int bx, by;
    size_t begin, end;
    uint8_t *level; /* leaf edge length, repeated over its finest cells */
} RlTile;

typedef struct {
    int x[8], y[8], tri[6][3], nt;
} RlTemplate;

typedef struct {
    const MeshBinData *mesh;
    const RlFace *face;
    size_t nf;
    int32_t *grid;
    int32_t *cell;
    uint8_t *level;
    int n, bx, by;
    int32_t *out_faces;
    size_t nout;
    double error;
    RibbonLodReport *report;
} RlPatch;

static int rl_div(int x, int n)
{ int q = x/n; return x%n < 0 ? q-1 : q; }

static int rl_compare(const void *pa, const void *pb)
{
    const RlFace *a = (const RlFace *)pa, *b = (const RlFace *)pb;
    if (a->by != b->by) return a->by < b->by ? -1 : 1;
    if (a->bx != b->bx) return a->bx < b->bx ? -1 : 1;
    if (a->y != b->y) return a->y < b->y ? -1 : 1;
    if (a->x != b->x) return a->x < b->x ? -1 : 1;
    return a->face < b->face ? -1 : a->face > b->face;
}

static int rl_vertex(const RlPatch *p, int x, int y)
{ return p->grid[y*(p->n+1)+x]; }

static void rl_xyz(const RlPatch *p, int x, int y, double v[3])
{
    int id = rl_vertex(p, x, y);
    for (int k = 0; k < 3; k++) v[k] = p->mesh->verts[3*(size_t)id+k];
}

/* Evaluate the original triangulated cell. cell records which diagonal is
 * present: 0 for (00,11), 1 for (10,01); -1 means incomplete/ambiguous. */
static void rl_original(const RlPatch *p, double x, double y, double out[3])
{
    int ix = (int)floor(x), iy = (int)floor(y);
    double u = 0.0, v = 0.0, w[3] = {0}, a[3] = {0}, b[3] = {0}, c[3] = {0};
    if (ix == p->n) ix--;
    if (iy == p->n) iy--;
    if (fabs(x-round(x))<1e-12 && fabs(y-round(y))<1e-12) {
        rl_xyz(p,(int)round(x),(int)round(y),out);
        return;
    }
    if (fabs(x-round(x))<1e-12) {
        double edge_a[3]={0}, edge_b[3]={0}, t=y-iy;
        rl_xyz(p,(int)round(x),iy,edge_a); rl_xyz(p,(int)round(x),iy+1,edge_b);
        for (int k=0; k<3; k++) out[k]=(1-t)*edge_a[k]+t*edge_b[k];
        return;
    }
    if (fabs(y-round(y))<1e-12) {
        double edge_a[3]={0}, edge_b[3]={0}, t=x-ix;
        rl_xyz(p,ix,(int)round(y),edge_a); rl_xyz(p,ix+1,(int)round(y),edge_b);
        for (int k=0; k<3; k++) out[k]=(1-t)*edge_a[k]+t*edge_b[k];
        return;
    }
    u = x-ix; v = y-iy;
    if (p->cell[iy*p->n+ix] == 0) {
        rl_xyz(p, ix, iy, a); rl_xyz(p, ix+1, iy+1, c);
        if (u >= v) { rl_xyz(p, ix+1, iy, b); w[0]=1-u; w[1]=u-v; w[2]=v; }
        else { rl_xyz(p, ix, iy+1, b); w[0]=1-v; w[1]=v-u; w[2]=u; }
    } else {
        rl_xyz(p, ix+1, iy, b); rl_xyz(p, ix, iy+1, c);
        if (u+v <= 1) { rl_xyz(p, ix, iy, a); w[0]=1-u-v; w[1]=u; w[2]=v; }
        else { rl_xyz(p, ix+1, iy+1, a); w[0]=u+v-1; w[1]=1-v; w[2]=1-u; }
    }
    for (int k = 0; k < 3; k++) out[k] = w[0]*a[k] + w[1]*b[k] + w[2]*c[k];
}

static double rl_distance(const double a[3], const double b[3])
{ double d = 0.0; for (int k = 0; k < 3; k++) d += (a[k]-b[k])*(a[k]-b[k]); return sqrt(d); }

/* Corner numbering is CCW from (0,0); mids are B,R,T,L. Rotate six
 * published canonical cases to cover all 16 masks. No center vertex and no
 * fan whose degree/aspect ratio grows with the refinement gap. */
static int rl_template(int x, int y, int n, int mask, RlTemplate *t)
{
    static const int masks[6] = {0,8,12,10,14,15};
    static const int counts[6] = {2,3,4,4,5,6};
    static const int tris[6][6][3] = {
        {{0,1,2},{0,2,3}},
        {{0,1,7},{7,1,2},{7,2,3}},
        {{0,1,7},{7,1,6},{1,2,6},{7,6,3}},
        {{0,1,5},{0,5,7},{7,5,3},{3,5,2}},
        {{0,1,6},{1,5,6},{5,2,6},{0,6,7},{7,6,3}},
        {{0,4,7},{4,1,5},{5,2,6},{6,3,7},{4,5,7},{5,6,7}}
    };
    static const int px[8]={0,2,2,0,1,2,1,0};
    static const int py[8]={0,0,2,2,0,1,2,1};
    memset(t,0,sizeof *t);
    for (int i=0; i<8; i++) { t->x[i]=x+px[i]*n/2; t->y[i]=y+py[i]*n/2; }
    for (int c=0; c<6; c++) for (int r=0; r<4; r++) {
        int rotated=((masks[c]<<r)|(masks[c]>>(4-r)))&15;
        if (rotated!=mask) continue;
        t->nt=counts[c];
        for (int f=0; f<t->nt; f++) for (int k=0; k<3; k++) {
            int v=tris[c][f][k];
            t->tri[f][k]=v<4 ? (v+r)%4 : 4+(v-4+r)%4;
        }
        return 0;
    }
    return -1;
}

int RibbonLod_transition(int mask, int triangles[6][3], int *count)
{
    RlTemplate t={0};
    if (mask<0 || mask>15 || triangles==NULL || count==NULL ||
        rl_template(0,0,2,mask,&t)!=0) return -1;
    memcpy(triangles,t.tri,sizeof t.tri); *count=t.nt;
    return 0;
}

static int rl_value(const RlPatch *p, const RlTemplate *t,
                     double x, double y, double out[3])
{
    for (int f=0; f<t->nt; f++) {
        int a=t->tri[f][0], b=t->tri[f][1], c=t->tri[f][2];
        double bx=t->x[b]-t->x[a], by=t->y[b]-t->y[a];
        double cx=t->x[c]-t->x[a], cy=t->y[c]-t->y[a];
        double dx=x-t->x[a], dy=y-t->y[a], den=bx*cy-by*cx;
        double u=(dx*cy-dy*cx)/den, v=(bx*dy-by*dx)/den;
        double av[3]={0}, bv[3]={0}, cv[3]={0};
        if (u < -1e-10 || v < -1e-10 || u+v > 1+1e-10) continue;
        rl_xyz(p,t->x[a],t->y[a],av); rl_xyz(p,t->x[b],t->y[b],bv);
        rl_xyz(p,t->x[c],t->y[c],cv);
        for (int k=0; k<3; k++) out[k]=(1-u-v)*av[k]+u*bv[k]+v*cv[k];
        return 0;
    }
    return -1;
}

static int rl_test_point(const RlPatch *p, const RlTemplate *t,
                         double x, double y, double *maximum)
{
    double original[3]={0}, coarse[3]={0}, error=0.0;
    rl_original(p,x,y,original);
    if (rl_value(p,t,x,y,coarse)!=0) return 0;
    error=rl_distance(original,coarse);
    if (error>*maximum) *maximum=error;
    return error<=p->error-1e-9;
}

static int rl_cross(const RlPatch *p, const RlTemplate *t,
                    double ax, double ay, double bx, double by,
                    double cx, double cy, double dx, double dy, double *maximum)
{
    double ex=bx-ax, ey=by-ay, fx=dx-cx, fy=dy-cy, den=ex*fy-ey*fx;
    double u=0.0, v=0.0;
    if (fabs(den)<1e-15) return 1;
    u=((cx-ax)*fy-(cy-ay)*fx)/den;
    v=((cx-ax)*ey-(cy-ay)*ex)/den;
    if (u<=1e-12 || u>=1-1e-12 || v<=1e-12 || v>=1-1e-12) return 1;
    return rl_test_point(p,t,ax+u*ex,ay+u*ey,maximum);
}

static int rl_complete(const RlPatch *p, int x, int y, int n)
{
    for (int j=y; j<y+n; j++) for (int i=x; i<x+n; i++)
        if (p->cell[j*p->n+i]<0) return 0;
    for (int j=y; j<=y+n; j++) for (int i=x; i<=x+n; i++)
        if (rl_vertex(p,i,j)<0) return 0;
    return 1;
}

/* Overlay vertices = fine vertices + intersections of coarse edges with
 * fine horizontal, vertical and diagonal edges. Difference is affine on
 * every overlay polygon; convexity of its norm bounds the whole surface.
 * Always audit the FINAL transition template, not only its unsplit quad. */
static int rl_accept(const RlPatch *p, int x, int y, int n, int mask, double *maximum)
{
    RlTemplate t={0};
    *maximum=0;
    if (n<=1 || !rl_complete(p,x,y,n) || rl_template(x,y,n,mask,&t)!=0) return 0;
    for (int j=y; j<=y+n; j++) for (int i=x; i<=x+n; i++)
        if (!rl_test_point(p,&t,i,j,maximum)) return 0;
    for (int f=0; f<t.nt; f++) for (int e=0; e<3; e++) {
        int a=t.tri[f][e], b=t.tri[f][(e+1)%3], seen=0;
        double ax=t.x[a], ay=t.y[a], bx=t.x[b], by=t.y[b];
        for (int g=0; g<f; g++) for (int h=0; h<3; h++) {
            int c=t.tri[g][h], d=t.tri[g][(h+1)%3];
            if ((c==a && d==b) || (c==b && d==a)) seen=1;
        }
        if (seen) continue;
        for (int k=x+1; k<x+n; k++) {
            double u=bx!=ax ? (k-ax)/(bx-ax) : -1;
            if (u>0 && u<1 && !rl_test_point(p,&t,k,ay+u*(by-ay),maximum)) return 0;
        }
        for (int k=y+1; k<y+n; k++) {
            double u=by!=ay ? (k-ay)/(by-ay) : -1;
            if (u>0 && u<1 && !rl_test_point(p,&t,ax+u*(bx-ax),k,maximum)) return 0;
        }
        for (int j=y; j<y+n; j++) for (int i=x; i<x+n; i++) {
            if (p->cell[j*p->n+i]==0) {
                if (!rl_cross(p,&t,ax,ay,bx,by,i,j,i+1,j+1,maximum)) return 0;
            } else if (!rl_cross(p,&t,ax,ay,bx,by,i+1,j,i,j+1,maximum)) return 0;
        }
    }
    return 1;
}

static void rl_level(RlPatch *p, int x, int y, int n, int level)
{
    for (int j=y; j<y+n; j++) for (int i=x; i<x+n; i++)
        p->level[j*p->n+i]=(uint8_t)level;
}

static void rl_select(RlPatch *p, int x, int y, int n)
{
    double error=0;
    if (n==1 || rl_accept(p,x,y,n,0,&error)) rl_level(p,x,y,n,n);
    else {
        rl_select(p,x,y,n/2); rl_select(p,x+n/2,y,n/2);
        rl_select(p,x,y+n/2,n/2); rl_select(p,x+n/2,y+n/2,n/2);
    }
}

static int rl_tile_find(const RlTile *tiles, size_t nt, int bx, int by)
{
    size_t lo=0, hi=nt;
    while (lo<hi) {
        size_t mid=lo+(hi-lo)/2;
        if (tiles[mid].by<by || (tiles[mid].by==by && tiles[mid].bx<bx)) lo=mid+1;
        else hi=mid;
    }
    return lo<nt && tiles[lo].bx==bx && tiles[lo].by==by ? (int)lo : -1;
}

static int rl_neighbor(const RlTile *tiles, size_t nt, int tile, int n, int x, int y)
{
    int bx=tiles[tile].bx, by=tiles[tile].by, found=tile;
    if (x<0) { x+=n; bx--; } else if (x>=n) { x-=n; bx++; }
    if (y<0) { y+=n; by--; } else if (y>=n) { y-=n; by++; }
    if (bx!=tiles[tile].bx || by!=tiles[tile].by) found=rl_tile_find(tiles,nt,bx,by);
    /* Missing support is finest level. It preserves every physical boundary
     * knot and creates a graded collar, not a skirt across absent geometry. */
    return found<0 ? 1 : tiles[found].level[y*n+x];
}

static int rl_mask(const RlTile *tiles, size_t nt, int ti, int root,
                    int x, int y, int n, int *must_split)
{
    int mask=0;
    *must_split=0;
    for (int k=0; k<n; k++) {
        int neighbors[4]={rl_neighbor(tiles,nt,ti,root,x+k,y-1),
                          rl_neighbor(tiles,nt,ti,root,x+n,y+k),
                          rl_neighbor(tiles,nt,ti,root,x+k,y+n),
                          rl_neighbor(tiles,nt,ti,root,x-1,y+k)};
        for (int e=0; e<4; e++) {
            if (neighbors[e]<n) mask|=1<<e;
            if (2*neighbors[e]<n) *must_split=1;
        }
    }
    return mask;
}

static void rl_grid(Arena_T arena, RlPatch *p, const RlTile *tile,
                     const int32_t *vx, const int32_t *vy, const uint8_t *pin)
{
    int n=p->n;
    int *masks=(int *)ARENA_CALLOC(arena,n*n,sizeof(int));
    int *counts=(int *)ARENA_CALLOC(arena,n*n,sizeof(int));
    p->bx=tile->bx; p->by=tile->by; p->level=tile->level;
    p->grid=(int32_t *)ARENA_ALLOC(arena,(n+1)*(n+1)*sizeof(int32_t));
    p->cell=(int32_t *)ARENA_ALLOC(arena,n*n*sizeof(int32_t));
    for (int j=0; j<(n+1)*(n+1); j++) p->grid[j]=-1;
    for (int j=0; j<n*n; j++) p->cell[j]=-1;
    for (size_t j=tile->begin; j<tile->end; j++) {
        const RlFace *r=p->face+j;
        const int32_t *tri=p->mesh->faces+3*(size_t)r->face;
        int mask=0, ci=r->y*n+r->x;
        counts[ci]++;
        for (int k=0; k<3; k++) {
            int id=tri[k], px=vx[id]-p->bx*n, py=vy[id]-p->by*n;
            int gi=py*(n+1)+px, corner=(py-r->y)*2+px-r->x;
            mask|=1<<corner;
            if (p->grid[gi]==-1) p->grid[gi]=id;
            else if (p->grid[gi]!=id) p->grid[gi]=-2;
            if (pin[id]) p->grid[gi]=-2;
        }
        if (mask==11) masks[ci]+=1;
        else if (mask==13) masks[ci]+=4;
        else if (mask==7) masks[ci]+=16;
        else if (mask==14) masks[ci]+=64;
        else masks[ci]+=256;
    }
    for (int j=0; j<n*n; j++) {
        if (counts[j]==2 && masks[j]==5) p->cell[j]=0;
        else if (counts[j]==2 && masks[j]==80) p->cell[j]=1;
    }
}

static int rl_emit(RlPatch *p, const RlTile *tiles, size_t nt, int ti)
{
    const RlTile *tile=tiles+ti;
    for (int y=0; y<p->n; y++) for (int x=0; x<p->n; x++) {
        int n=p->level[y*p->n+x], split=0, mask=0;
        double error=0;
        RlTemplate t={0};
        if (x%n || y%n) continue;
        if (n==1) {
            for (size_t i=tile->begin; i<tile->end; i++) if (p->face[i].x==x && p->face[i].y==y) {
                memcpy(p->out_faces+3*p->nout,p->mesh->faces+3*(size_t)p->face[i].face,3*sizeof(int32_t));
                p->nout++; p->report->retained_faces++;
            }
        } else {
            mask=rl_mask(tiles,nt,ti,p->n,x,y,n,&split);
            if (split || !rl_accept(p,x,y,n,mask,&error) || rl_template(x,y,n,mask,&t)!=0) return -1;
            for (int f=0; f<t.nt; f++) {
                for (int k=0; k<3; k++) {
                    int v=t.tri[f][k];
                    p->out_faces[3*p->nout+k]=rl_vertex(p,t.x[v],t.y[v]);
                }
                p->nout++;
            }
            p->report->patches++; p->report->transition_cases[mask]++;
            if (error>p->report->max_error) p->report->max_error=error;
        }
        if (n>1 || rl_complete(p,x,y,n)) {
            RibbonLodQuad *q=p->report->quads+p->report->nquads++;
            q->x=p->bx*p->n+x; q->y=p->by*p->n+y; q->cells=n; q->transition_mask=mask;
            q->corner[0]=rl_vertex(p,x,y); q->corner[1]=rl_vertex(p,x+n,y);
            q->corner[2]=rl_vertex(p,x+n,y+n); q->corner[3]=rl_vertex(p,x,y+n);
        }
    }
    return 0;
}

int RibbonLod_build(Arena_T arena, const MeshBinData *in,
                     double du, double dv, int max_cells, double max_error,
                     MeshBinData *out, int32_t **source_vertex,
                     RibbonLodReport *report)
{
    RlFace *ref = NULL;
    RlTile *tiles = NULL;
    int32_t *edges = NULL;
    int32_t *vx = NULL, *vy = NULL, *faces = NULL, *map = NULL;
    uint8_t *pin = NULL;
    size_t nref = 0, nout = 0, nv = 0, nt = 0;
    if (arena==NULL || in==NULL || out==NULL || source_vertex==NULL || report==NULL ||
        !(du>0) || !(dv>0) || !isfinite(du) || !isfinite(dv) ||
        !(max_error>0) || !isfinite(max_error) || max_cells<1 || max_cells>64 ||
        (max_cells&(max_cells-1)) || in->nv>INT32_MAX || in->nf>INT32_MAX ||
        (in->nv && (in->verts==NULL || in->uv==NULL)) || (in->nf && in->faces==NULL)) return -1;
    memset(out,0,sizeof *out); memset(report,0,sizeof *report);
    report->input_faces=in->nf;
    ref=(RlFace *)ARENA_ALLOC(arena,in->nf*sizeof *ref);
    faces=(int32_t *)ARENA_ALLOC(arena,3*in->nf*sizeof(int32_t));
    vx=(int32_t *)ARENA_ALLOC(arena,in->nv*sizeof(int32_t));
    vy=(int32_t *)ARENA_ALLOC(arena,in->nv*sizeof(int32_t));
    map=(int32_t *)ARENA_ALLOC(arena,in->nv*sizeof(int32_t));
    pin=(uint8_t *)ARENA_CALLOC(arena,in->nv,1);
    for (size_t i=0; i<in->nv; i++) {
        double x=in->uv[2*i]/du, y=in->uv[2*i+1]/dv;
        double ix=round(x), iy=round(y);
        vx[i]=vy[i]=INT32_MIN; map[i]=-2;
        if (!isfinite(x) || !isfinite(y) || !isfinite(in->verts[3*i]) ||
            !isfinite(in->verts[3*i+1]) || !isfinite(in->verts[3*i+2])) return -1;
        if (ix<=INT32_MIN+64.0 || ix>=INT32_MAX-64.0 ||
            iy<=INT32_MIN+64.0 || iy>=INT32_MAX-64.0 || x!=ix || y!=iy) continue;
        /* Exactness here is an ADDRESS check, not a geometric tolerance:
         * snapping a near-lattice UV would invalidate the overlay bound. */
        vx[i]=(int32_t)ix; vy[i]=(int32_t)iy;
    }
    for (size_t f=0; f<in->nf; f++) {
        int x=INT_MAX,y=INT_MAX,xhi=INT_MIN,yhi=INT_MIN,ok=1;
        double signed_area=0.0;
        const int32_t *tri=in->faces+3*f;
        for (int k=0; k<3; k++) {
            int id=tri[k];
            if (id<0 || (size_t)id>=in->nv) return -1;
            map[id]=-1;
            if (vx[id]==INT32_MIN || vy[id]==INT32_MIN) ok=0;
            if (vx[id]<x) x=vx[id]; if (vx[id]>xhi) xhi=vx[id];
            if (vy[id]<y) y=vy[id]; if (vy[id]>yhi) yhi=vy[id];
        }
        if (ok) signed_area=((double)vx[tri[1]]-vx[tri[0]])*((double)vy[tri[2]]-vy[tri[0]]) -
                           ((double)vy[tri[1]]-vy[tri[0]])*((double)vx[tri[2]]-vx[tri[0]]);
        if (!ok || (int64_t)xhi-x!=1 || (int64_t)yhi-y!=1 || signed_area!=1) {
            for (int k=0; k<3; k++) pin[tri[k]]=1;
            memcpy(faces+3*nout,tri,3*sizeof(int32_t)); nout++; report->retained_faces++; continue;
        }
        ref[nref].bx=rl_div(x,max_cells); ref[nref].by=rl_div(y,max_cells);
        ref[nref].x=x-ref[nref].bx*max_cells; ref[nref].y=y-ref[nref].by*max_cells;
        ref[nref].face=(int32_t)f; nref++;
    }
    qsort(ref,nref,sizeof *ref,rl_compare);
    /* Sparse root directory. Fine occupancy is one byte/cell in occupied
     * roots only; geometry/error scratch remains bounded to one root. */
    for (size_t i=0; i<nref; i++)
        if (i==0 || ref[i].bx!=ref[i-1].bx || ref[i].by!=ref[i-1].by) nt++;
    tiles=(RlTile *)ARENA_CALLOC(arena,nt,sizeof *tiles);
    edges=(int32_t *)ARENA_ALLOC(arena,nt*4*(max_cells+1)*sizeof(int32_t));
    report->quads=(RibbonLodQuad *)ARENA_ALLOC(arena,(in->nf/2)*sizeof(RibbonLodQuad));
    for (size_t i=0, ti=0; i<nref; ti++) {
        size_t end=i+1;
        Arena_Mark mark;
        RlPatch p={0};
        while (end<nref && ref[end].bx==ref[i].bx && ref[end].by==ref[i].by) end++;
        tiles[ti].bx=ref[i].bx; tiles[ti].by=ref[i].by;
        tiles[ti].begin=i; tiles[ti].end=end;
        tiles[ti].level=(uint8_t *)ARENA_CALLOC(arena,max_cells*max_cells,1);
        mark=Arena_save(arena);
        p.mesh=in; p.face=ref; p.n=max_cells;
        rl_grid(arena,&p,tiles+ti,vx,vy,pin);
        for (int k=0; k<=max_cells; k++) {
            size_t at=ti*4*(max_cells+1)+(size_t)k;
            edges[at]=rl_vertex(&p,k,0);
            edges[at+max_cells+1]=rl_vertex(&p,max_cells,k);
            edges[at+2*(max_cells+1)]=rl_vertex(&p,k,max_cells);
            edges[at+3*(max_cells+1)]=rl_vertex(&p,0,k);
        }
        Arena_restore(arena,mark);
        i=end;
    }
    /* Distinct source vertices at coincident UV across root edges are cuts,
     * not weld permission. Pin both sides before making ANY selection. */
    for (size_t ti=0; ti<nt; ti++) for (int e=1; e<=2; e++) {
        int other=rl_tile_find(tiles,nt,tiles[ti].bx+(e==1),tiles[ti].by+(e==2));
        if (other<0) continue;
        for (int k=0; k<=max_cells; k++) {
            int a=edges[(ti*4+e)*(max_cells+1)+k];
            int b=edges[((size_t)other*4+(e+2)%4)*(max_cells+1)+k];
            if (a!=b) { if (a>=0) pin[a]=1; if (b>=0) pin[b]=1; }
        }
    }
    for (size_t ti=0; ti<nt; ti++) {
        Arena_Mark mark=Arena_save(arena);
        RlPatch p={0};
        p.mesh=in; p.face=ref; p.n=max_cells; p.error=max_error;
        rl_grid(arena,&p,tiles+ti,vx,vy,pin);
        rl_select(&p,0,0,max_cells);
        Arena_restore(arena,mark);
    }
    /* Refinement-only closure terminates at the unchanged source grid.
     * Balance first, then audit actual transition templates, then rebalance.
     * No view/camera-dependent choice or coordinate smoothing occurs here. */
    for (;;) {
        size_t changed=0, balancing=0;
        do {
            balancing=0;
            for (size_t ti=0; ti<nt; ti++) {
                RlPatch p={0};
                p.n=max_cells; p.level=tiles[ti].level;
                for (int y=0; y<max_cells; y++) for (int x=0; x<max_cells; x++) {
                    int n=p.level[y*max_cells+x], split=0;
                    if (n<=1 || x%n || y%n) continue;
                    (void)rl_mask(tiles,nt,(int)ti,max_cells,x,y,n,&split);
                    if (split) {
                        rl_level(&p,x,y,n,n/2);
                        balancing++; report->balance_splits++;
                    }
                }
            }
        } while (balancing);
        for (size_t ti=0; ti<nt; ti++) {
            Arena_Mark mark=Arena_save(arena);
            RlPatch p={0};
            p.mesh=in; p.face=ref; p.n=max_cells; p.error=max_error;
            rl_grid(arena,&p,tiles+ti,vx,vy,pin);
            for (int y=0; y<max_cells; y++) for (int x=0; x<max_cells; x++) {
                int n=p.level[y*max_cells+x], split=0, mask=0;
                double error=0;
                if (n<=1 || x%n || y%n) continue;
                mask=rl_mask(tiles,nt,(int)ti,max_cells,x,y,n,&split);
                if (!rl_accept(&p,x,y,n,mask,&error)) {
                    rl_level(&p,x,y,n,n/2);
                    changed++; report->error_splits++;
                }
            }
            Arena_restore(arena,mark);
        }
        if (!changed) break;
    }
    for (size_t ti=0; ti<nt; ti++) {
        Arena_Mark mark=Arena_save(arena);
        RlPatch p={0};
        p.mesh=in; p.face=ref; p.n=max_cells; p.error=max_error;
        p.out_faces=faces; p.nout=nout; p.report=report;
        rl_grid(arena,&p,tiles+ti,vx,vy,pin);
        if (rl_emit(&p,tiles,nt,(int)ti)!=0) return -1;
        nout=p.nout;
        Arena_restore(arena,mark);
    }
    if (nout>in->nf) return -1;
    for (size_t i=0; i<3*nout; i++) map[faces[i]]=0;
    /* Existing isolated source points are not part of any replaced disk. */
    for (size_t i=0; i<in->nv; i++) if (map[i]==-2) { map[i]=0; report->isolated_vertices++; }
    for (size_t i=0; i<in->nv; i++) if (map[i]==0) map[i]=(int32_t)nv++;
    out->verts=(float *)ARENA_ALLOC(arena,3*nv*sizeof(float));
    out->uv=(float *)ARENA_ALLOC(arena,2*nv*sizeof(float));
    out->faces=faces; out->nv=nv; out->nf=nout;
    *source_vertex=(int32_t *)ARENA_ALLOC(arena,nv*sizeof(int32_t));
    for (size_t i=0; i<in->nv; i++) if (map[i]>=0) {
        memcpy(out->verts+3*(size_t)map[i],in->verts+3*i,3*sizeof(float));
        memcpy(out->uv+2*(size_t)map[i],in->uv+2*i,2*sizeof(float));
        (*source_vertex)[map[i]]=(int32_t)i;
    }
    for (size_t i=0; i<3*nout; i++) faces[i]=map[faces[i]];
    for (size_t i=0; i<report->nquads; i++) for (int k=0; k<4; k++)
        report->quads[i].corner[k]=map[report->quads[i].corner[k]];
    report->output_faces=nout;
    return 0;
}

int RibbonLod_write(const char *input_path, const char *output_path,
                     double du, double dv, int max_cells, double max_error)
{
    Arena_T arena=Arena_new();
    MeshBinData input={0}, output={0};
    RibbonLodReport report={0};
    int32_t *source=NULL;
    char path[2048];
    FILE *f=NULL;
    int rc=-1;
    double started=ves_clock_sec();
    if (MeshBin_read_arena(arena,input_path,&input)!=0 ||
        RibbonLod_build(arena,&input,du,dv,max_cells,max_error,&output,&source,&report)!=0 ||
        MeshBin_write(output_path,output.verts,output.nv,output.faces,output.nf,output.uv)!=0) goto done;
    snprintf(path,sizeof path,"%s.source_vertex.i32",output_path);
    f=fopen(path,"wb");
    if (f==NULL) goto done;
    if (fwrite(source,sizeof(int32_t),output.nv,f)!=output.nv) goto done;
    if (fclose(f)!=0) { f=NULL; goto done; } f=NULL;
    snprintf(path,sizeof path,"%s.quads.bin",output_path);
    f=fopen(path,"wb");
    if (f==NULL) goto done;
    {
        const uint64_t header[3]={UINT64_C(0x3153444155514c52),1,(uint64_t)report.nquads};
        if (fwrite(header,sizeof header,1,f)!=1 ||
            fwrite(report.quads,sizeof(RibbonLodQuad),report.nquads,f)!=report.nquads) goto done;
    }
    if (fclose(f)!=0) { f=NULL; goto done; } f=NULL;
    snprintf(path,sizeof path,"%s.lod.json",output_path);
    f=fopen(path,"wb");
    if (f==NULL) goto done;
    fprintf(f,"{\"schema\":\"ribbon-lod-v2\",\"scheme\":\"restricted-quadtree-16-primal-transition-templates\","
            "\"input_faces\":%zu,\"output_faces\":%zu,\"patches\":%zu,\"quad_leaves\":%zu,\"retained_faces\":%zu,"
            "\"max_continuous_error_voxels\":%.12g,\"error_limit_voxels\":%.12g,\"max_patch_cells\":%d,\"step_uv\":[%.17g,%.17g],"
            "\"max_edge_neighbor_ratio\":2,\"balance_splits\":%zu,\"transition_error_splits\":%zu,"
            "\"uv_and_boundary_vertices\":\"bitwise retained\",\"retained_isolated_vertices\":%zu,\"elapsed_seconds\":%.6f,\"transition_cases\":[",
            report.input_faces,report.output_faces,report.patches,report.nquads,report.retained_faces,
            report.max_error,max_error,max_cells,du,dv,report.balance_splits,report.error_splits,
            report.isolated_vertices,ves_clock_sec()-started);
    for (int mask=0; mask<16; mask++) fprintf(f,"%s%zu",mask ? "," : "",report.transition_cases[mask]);
    fprintf(f,"]}\n");
    if (fclose(f)!=0) { f=NULL; goto done; } f=NULL;
    fprintf(stderr,"[ribbon LOD] %zu -> %zu faces, %zu patches, max continuous error %.6g / %.6g vox\n",
            input.nf,output.nf,report.patches,report.max_error,max_error);
    rc=0;
done:
    if (f!=NULL) fclose(f);
    Arena_dispose(&arena);
    return rc;
}

typedef struct { int32_t a,b; } RlEdge;

static int rl_edge_compare(const void *pa, const void *pb)
{
    const RlEdge *a=(const RlEdge *)pa, *b=(const RlEdge *)pb;
    if (a->a!=b->a) return a->a<b->a ? -1 : 1;
    return a->b<b->b ? -1 : a->b>b->b;
}

static size_t rl_boundary(Arena_T arena, const MeshBinData *m,
                          const int32_t *source, RlEdge **out)
{
    RlEdge *edges=(RlEdge *)ARENA_ALLOC(arena,m->nf*3*sizeof(RlEdge));
    size_t count=0;
    for (size_t f=0; f<m->nf; f++) for (int k=0; k<3; k++) {
        int a=m->faces[3*f+k], b=m->faces[3*f+(k+1)%3];
        if (source!=NULL) { a=source[a]; b=source[b]; }
        edges[3*f+k].a=a<b ? a : b; edges[3*f+k].b=a<b ? b : a;
    }
    qsort(edges,3*m->nf,sizeof(RlEdge),rl_edge_compare);
    for (size_t i=0; i<3*m->nf; ) {
        size_t end=i+1;
        while (end<3*m->nf && rl_edge_compare(edges+i,edges+end)==0) end++;
        if (end-i!=2) edges[count++]=edges[i];
        i=end;
    }
    *out=edges;
    return count;
}

static int rl_audit(Arena_T arena, const MeshBinData *in, const MeshBinData *out,
                     const int32_t *source)
{
    RlEdge *a=NULL, *b=NULL;
    size_t na=rl_boundary(arena,in,NULL,&a), nb=rl_boundary(arena,out,source,&b);
    int fails=na!=nb;
    if (na==nb) fails+=memcmp(a,b,na*sizeof(RlEdge))!=0;
    for (size_t i=0; i<out->nv; i++) {
        fails+=memcmp(out->uv+2*i,in->uv+2*(size_t)source[i],2*sizeof(float))!=0;
        fails+=memcmp(out->verts+3*i,in->verts+3*(size_t)source[i],3*sizeof(float))!=0;
    }
    return fails;
}

int RibbonLod_selftest(void)
{
    enum { N=32, NV=(N+1)*(N+1), NF=N*N*2 };
    Arena_T arena=Arena_new();
    float v[3*NV]={0}, uv[2*NV]={0};
    int32_t f[3*NF]={0}, *source=NULL;
    MeshBinData in={v,uv,f,NV,NF}, out={0}, repeat={0}, empty={0};
    RibbonLodReport report={0}, again={0};
    int fails=0;
    size_t planar_faces=0;
    for (int mask=0; mask<16; mask++) {
        RlTemplate t={0};
        int area=0, boundary[8]={0};
        fails+=rl_template(0,0,2,mask,&t)!=0;
        for (int i=0; i<t.nt; i++) {
            int a=t.tri[i][0], b=t.tri[i][1], c=t.tri[i][2];
            int cross=(t.x[b]-t.x[a])*(t.y[c]-t.y[a])-(t.y[b]-t.y[a])*(t.x[c]-t.x[a]);
            fails+=cross<=0; area+=cross;
            for (int k=0; k<3; k++) {
                int j=t.tri[i][k], l=t.tri[i][(k+1)%3], h=t.tri[i][(k+2)%3];
                double ux=t.x[l]-t.x[j], uy=t.y[l]-t.y[j];
                double vx=t.x[h]-t.x[j], vy=t.y[h]-t.y[j];
                double cosine=(ux*vx+uy*vy)/sqrt((ux*ux+uy*uy)*(vx*vx+vy*vy));
                fails+=cosine>3/sqrt(10.0)+1e-12; /* min angle >= atan(1/3) */
                boundary[j]=1;
            }
        }
        fails+=area!=8;
        for (int e=0; e<4; e++) fails+=boundary[4+e]!=((mask>>e)&1);
    }
    for (int y=0; y<=N; y++) for (int x=0; x<=N; x++) {
        int i=y*(N+1)+x;
        v[3*i]=uv[2*i+1]=(float)(y-16);
        v[3*i+1]=uv[2*i]=(float)(x-16); v[3*i+2]=2;
    }
    for (int y=0; y<N; y++) for (int x=0; x<N; x++) {
        int a=y*(N+1)+x, at=6*(y*N+x);
        f[at]=a; f[at+1]=a+1; f[at+2]=a+N+2;
        f[at+3]=a; f[at+4]=a+N+2; f[at+5]=a+N+1;
    }
    fails+=RibbonLod_build(arena,&in,1,1,16,.1,&out,&source,&report)!=0;
    planar_faces=out.nf;
    fails+=!(out.nf<NF/2 && report.patches>0 && report.max_error<1e-12);
    fails+=rl_audit(arena,&in,&out,source);
    fails+=RibbonLod_build(arena,&in,1,1,16,.1,&repeat,&source,&again)!=0;
    fails+=repeat.nf!=out.nf || repeat.nv!=out.nv;
    if (repeat.nf==out.nf) fails+=memcmp(repeat.faces,out.faces,3*out.nf*sizeof(int32_t))!=0;
    /* Off-diagonal bump exercises the overlay, not only coarse corners.
     * A boundary bump exercises conformity across absolute root edges. */
    v[3*(13*(N+1)+16)+2]=10;
    v[3*(19*(N+1)+21)+2]=7;
    fails+=RibbonLod_build(arena,&in,1,1,16,.1,&out,&source,&report)!=0;
    fails+=!(out.nf>planar_faces && report.max_error<=.1);
    fails+=rl_audit(arena,&in,&out,source);
    in.nf--; /* one open half-cell, never patched across */
    fails+=RibbonLod_build(arena,&in,1,1,16,.1,&out,&source,&report)!=0;
    fails+=!(report.retained_faces>0 && out.nf<=in.nf);
    fails+=rl_audit(arena,&in,&out,source);
    fails+=RibbonLod_build(arena,&empty,1,1,16,.1,&out,&source,&report)!=0;
    fails+=out.nf!=0 || out.nv!=0;
    empty.verts=v; empty.uv=uv; empty.nv=1;
    fails+=RibbonLod_build(arena,&empty,1,1,16,.1,&out,&source,&report)!=0;
    fails+=out.nv!=1 || out.nf!=0 || report.isolated_vertices!=1;
    fprintf(stderr,"[selftest] restricted quadtree LOD %s (%d failures): 16 templates, overlay error, boundaries, holes, negative origins, repeatability\n",
            fails ? "FAIL" : "PASS",fails);
    Arena_dispose(&arena);
    return fails ? -1 : 0;
}
