#define _USE_MATH_DEFINES
#include "seam_cut.h"
#include "mesh_topo.h"
#include "../common/ves_platform.h"   /* ves_clock_sec for phase timing */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===================================================================
 * Mesh connectivity (half-edge "corner" model).
 *   corner c lives in face c/3, local slot c%3; vertex AT c = faces[c];
 *   directed edge of c = (faces[c] -> faces[NEXT(c)]); twin = opposite.
 * =================================================================== */
#define NEXT(c) (((c)/3)*3 + (((c)%3)+1)%3)
#define PREV(c) (((c)/3)*3 + (((c)%3)+2)%3)

typedef struct {
    size_t         nv, nf, nc;
    const float   *verts;     /* [nv*3] */
    const int32_t *faces;     /* [nf*3] */
    int32_t       *twin;      /* [nc], -1 = boundary */
    int32_t       *voff;      /* [nv+1] CSR offsets into vcor */
    int32_t       *vcor;      /* [nc] incident corners grouped by vertex */
    int            nonmanifold;
} Conn;

typedef struct { uint64_t key; int32_t corner; } EKey;
static int cmp_ekey(const void *a, const void *b)
{
    uint64_t x = ((const EKey *)a)->key, y = ((const EKey *)b)->key;
    return (x < y) ? -1 : (x > y ? 1 : 0);
}

static void conn_build(Arena_T arena, const float *verts, size_t nv,
                       const int32_t *faces, size_t nf, Conn *c)
{
    c->nv = nv; c->nf = nf; c->nc = 3*nf; c->verts = verts; c->faces = faces;
    c->nonmanifold = 0;
    c->twin = (int32_t *)ARENA_ALLOC(arena, (long)(c->nc * sizeof(int32_t)));
    for (size_t i = 0; i < c->nc; i++) c->twin[i] = -1;

    /* twins via sorted undirected-edge keys */
    EKey *ek = (EKey *)ARENA_ALLOC(arena, (long)(c->nc * sizeof(EKey)));
    for (size_t cc = 0; cc < c->nc; cc++) {
        int32_t a = faces[cc], b = faces[NEXT(cc)];
        uint64_t lo = (uint64_t)(a < b ? a : b), hi = (uint64_t)(a < b ? b : a);
        ek[cc].key = lo * (uint64_t)nv + hi;
        ek[cc].corner = (int32_t)cc;
    }
    qsort(ek, c->nc, sizeof(EKey), cmp_ekey);
    size_t i = 0;
    while (i < c->nc) {
        size_t j = i; while (j < c->nc && ek[j].key == ek[i].key) j++;
        if (j - i == 2) {
            c->twin[ek[i].corner] = ek[i+1].corner;
            c->twin[ek[i+1].corner] = ek[i].corner;
        } else if (j - i > 2) {
            c->nonmanifold = 1;   /* >2 faces on an edge: leave as boundary */
        }
        i = j;
    }

    /* vertex -> incident corners (CSR) */
    c->voff = (int32_t *)ARENA_CALLOC(arena, (long)(nv + 1), (long)sizeof(int32_t));
    for (size_t cc = 0; cc < c->nc; cc++) c->voff[faces[cc] + 1]++;
    for (size_t v = 0; v < nv; v++) c->voff[v+1] += c->voff[v];
    c->vcor = (int32_t *)ARENA_ALLOC(arena, (long)(c->nc * sizeof(int32_t)));
    int32_t *cur = (int32_t *)ARENA_ALLOC(arena, (long)(nv * sizeof(int32_t)));
    for (size_t v = 0; v < nv; v++) cur[v] = c->voff[v];
    for (size_t cc = 0; cc < c->nc; cc++) { int32_t v = faces[cc]; c->vcor[cur[v]++] = (int32_t)cc; }
}

/* ---- geometry helpers (double) ---- */
static void vsub(const float *p, const float *q, double o[3])
{ o[0]=(double)p[0]-q[0]; o[1]=(double)p[1]-q[1]; o[2]=(double)p[2]-q[2]; }
static double vdot(const double a[3], const double b[3])
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static double vlen(const double a[3]) { return sqrt(vdot(a,a)); }
/* angle at vertex p0 between p1 and p2 */
static double corner_angle(const float *p0, const float *p1, const float *p2)
{
    double a[3], b[3]; vsub(p1,p0,a); vsub(p2,p0,b);
    double la=vlen(a), lb=vlen(b);
    if (la<1e-20 || lb<1e-20) return 0.0;
    double c = vdot(a,b)/(la*lb);
    if (c>1.0)c=1.0; if(c<-1.0)c=-1.0;
    return acos(c);
}

/* ===================================================================
 * Seam edge set: hash of undirected vertex pairs (open addressing).
 * =================================================================== */
typedef struct { uint64_t *slot; size_t cap, n; size_t nv; } SeamSet;
static void seam_init(Arena_T arena, SeamSet *s, size_t nv, size_t hint)
{
    size_t cap = 16; while (cap < hint*2) cap <<= 1;
    s->slot = (uint64_t *)ARENA_ALLOC(arena, (long)(cap*sizeof(uint64_t)));
    for (size_t i=0;i<cap;i++) s->slot[i] = UINT64_MAX;
    s->cap = cap; s->n = 0; s->nv = nv;
}
static uint64_t seam_key(const SeamSet *s, int32_t a, int32_t b)
{ uint64_t lo=(uint64_t)(a<b?a:b), hi=(uint64_t)(a<b?b:a); return lo*(uint64_t)s->nv+hi; }
static int seam_has(const SeamSet *s, int32_t a, int32_t b)
{
    uint64_t k=seam_key(s,a,b), m=s->cap-1, i=k&m;
    while (s->slot[i]!=UINT64_MAX){ if(s->slot[i]==k) return 1; i=(i+1)&m; }
    return 0;
}
static void seam_add(SeamSet *s, int32_t a, int32_t b)
{
    if (a==b) return;
    uint64_t k=seam_key(s,a,b), m=s->cap-1, i=k&m;
    while (s->slot[i]!=UINT64_MAX){ if(s->slot[i]==k) return; i=(i+1)&m; }
    s->slot[i]=k; s->n++;
}

/* ===================================================================
 * Section 5: region distortion D_r(n).
 * =================================================================== */

/* D'_0(n): standard angle deficit over the one-ring (interior). */
static double distortion_r0(const Conn *c, int32_t n)
{
    double sum = 0.0;
    for (int32_t e = c->voff[n]; e < c->voff[n+1]; e++) {
        int32_t cc = c->vcor[e];               /* corner at n in some face */
        const float *p0 = c->verts + 3*(size_t)n;
        const float *p1 = c->verts + 3*(size_t)c->faces[NEXT(cc)];
        const float *p2 = c->verts + 3*(size_t)c->faces[PREV(cc)];
        sum += corner_angle(p0, p1, p2);
    }
    return (2.0*M_PI - sum) / (2.0*M_PI);
}

/* longest edge emanating from n */
static double longest_edge(const Conn *c, int32_t n)
{
    double best = 0.0;
    const float *pn = c->verts + 3*(size_t)n;
    for (int32_t e = c->voff[n]; e < c->voff[n+1]; e++) {
        int32_t cc = c->vcor[e];
        double d[3];
        vsub(c->verts + 3*(size_t)c->faces[NEXT(cc)], pn, d); double l=vlen(d);
        if (l>best) best=l;
        vsub(c->verts + 3*(size_t)c->faces[PREV(cc)], pn, d); l=vlen(d);
        if (l>best) best=l;
    }
    return best;
}

/* D'_i(n), i>=1: grow region (all-verts within R), sum subtended angles at n
 * over region boundary edges. Returns the deficit, or NAN if the region is not
 * a topological disk (multi-loop -> Section 5 concern (a): skip). Scratch
 * arrays (size nf / nv) are caller-provided and reset on use. */
static double distortion_ri(Arena_T arena, const Conn *c, int32_t n, double R,
                            int32_t *facemark, int32_t mark, int32_t *queue,
                            int32_t *vstamp, int32_t *vid)
{
    const float *pn = c->verts + 3*(size_t)n;
    /* BFS faces within R (all 3 verts within R of n), seeded from n's faces. */
    size_t head=0, tail=0;
    for (int32_t e = c->voff[n]; e < c->voff[n+1]; e++) {
        int32_t f = c->vcor[e] / 3;
        if (facemark[f] != mark) { facemark[f]=mark; queue[tail++]=f; }
    }
    size_t region_start = 0; (void)region_start;
    size_t qn = tail;
    while (head < qn) {
        int32_t f = queue[head++];
        for (int k=0;k<3;k++){
            int32_t cc = f*3+k; int32_t tw = c->twin[cc];
            if (tw < 0) continue;
            int32_t nf2 = tw/3;
            if (facemark[nf2]==mark) continue;
            /* include nf2 if all its verts within R of n */
            int ok=1;
            for (int t=0;t<3;t++){ double d[3]; vsub(c->verts+3*(size_t)c->faces[nf2*3+t], pn, d);
                if (vlen(d) > R) { ok=0; break; } }
            if (ok){ facemark[nf2]=mark; queue[tail++]=nf2; qn=tail; }
        }
    }
    /* region faces = queue[0..tail). Boundary edges + subtended angle sum.
     * Also collect boundary endpoints to verify a single loop (disk). */
    double sum=0.0; long nbe=0;
    /* small union-find over boundary endpoints to count loops */
    Arena_Mark m2 = Arena_save(arena);
    int32_t *bp = (int32_t *)ARENA_ALLOC(arena, (long)(tail*3*2*sizeof(int32_t)));
    long bpn=0;
    for (size_t qi=0; qi<tail; qi++){
        int32_t f=queue[qi];
        for (int k=0;k<3;k++){
            int32_t cc=f*3+k; int32_t tw=c->twin[cc];
            int is_bnd = (tw<0) || (facemark[tw/3]!=mark);
            if (!is_bnd) continue;
            int32_t a=c->faces[cc], b=c->faces[NEXT(cc)];
            if (a==n || b==n) continue;     /* edge touching n: skip subtend */
            sum += corner_angle(pn, c->verts+3*(size_t)a, c->verts+3*(size_t)b);
            bp[bpn++]=a; bp[bpn++]=b; nbe++;
        }
    }
    /* count boundary loops: compact the endpoints (stamp = mark), union the two
     * endpoints of each boundary edge, count roots. O(bpn). A disk has 1 loop. */
    int loops_ok = 1;
    if (nbe >= 1) {
        int ncomp2 = 0;
        for (long i=0;i<bpn;i++){ int32_t v=bp[i]; if(vstamp[v]!=mark){vstamp[v]=mark; vid[v]=ncomp2++;} }
        int32_t *par = (int32_t *)ARENA_ALLOC(arena, (long)((size_t)ncomp2*sizeof(int32_t)));
        for (int k=0;k<ncomp2;k++) par[k]=k;
        for (long i=0;i<bpn;i+=2){
            int32_t x=vid[bp[i]], y=vid[bp[i+1]];
            while(par[x]!=x)x=par[x]; while(par[y]!=y)y=par[y];
            if(x!=y)par[x]=y;
        }
        long nloops=0; for(int k=0;k<ncomp2;k++){ int r=k; while(par[r]!=r)r=par[r]; if(r==k)nloops++; }
        if (nloops != 1) loops_ok = 0;
    } else loops_ok = 0;
    Arena_restore(arena, m2);
    if (!loops_ok) return NAN;
    return (2.0*M_PI - sum) / (2.0*M_PI);
}

/* Compute D_r(n) for all vertices (Section 5), with overlap filtering (b). */
static void compute_distortion(Arena_T arena, const Conn *c, int r_max,
                               double *Dr, int32_t *best_r)
{
    size_t nv=c->nv, nf=c->nf;
    int32_t *facemark = (int32_t *)ARENA_CALLOC(arena, (long)nf, (long)sizeof(int32_t));
    int32_t *queue = (int32_t *)ARENA_ALLOC(arena, (long)(nf*sizeof(int32_t)));
    int32_t *vstamp = (int32_t *)ARENA_CALLOC(arena, (long)nv, (long)sizeof(int32_t));
    int32_t *vid = (int32_t *)ARENA_ALLOC(arena, (long)(nv*sizeof(int32_t)));
    int32_t markctr = 0;

    for (size_t v=0; v<nv; v++) {
        /* boundary vertex with angle sum < 2pi -> developable -> 0 */
        double d0 = distortion_r0(c, (int32_t)v);
        double best = d0; int32_t br = 0;
        double le = longest_edge(c, (int32_t)v);
        for (int i=1; i<=r_max; i++) {
            markctr++;
            double di = distortion_ri(arena, c, (int32_t)v, (double)i*le,
                                      facemark, markctr, queue, vstamp, vid);
            if (di != di) break;          /* NAN -> region not a disk, stop */
            if (di > best) { best = di; br = i; }
        }
        Dr[v] = best; best_r[v] = br;
    }

    /* Section 5 (b): overlap filtering. For nearby high-distortion vertices,
     * demote the lesser to its r=0 value so a single region isn't counted
     * twice. We approximate "n inside M_i(m)" by Euclidean proximity within
     * best_r(m)*longest_edge(m). One pass, larger D_r wins. */
    /* (Cheap version: only adjust along mesh edges, which dominate proximity.) */
    for (size_t m=0; m<nv; m++) {
        if (Dr[m] <= 0.0) continue;
        for (int32_t e=c->voff[m]; e<c->voff[m+1]; e++) {
            int32_t cc=c->vcor[e];
            int32_t nb[2] = { c->faces[NEXT(cc)], c->faces[PREV(cc)] };
            for (int t=0;t<2;t++){
                int32_t n=nb[t];
                if (Dr[n] <= 0.0) continue;
                if (Dr[n] > Dr[m]) { Dr[m] = distortion_r0(c,(int32_t)m); best_r[m]=0; }
                else if (Dr[n] < Dr[m]) { Dr[n] = distortion_r0(c,n); best_r[n]=0; }
            }
        }
    }
}

/* ===================================================================
 * Section 6.2: approximate minimal Steiner tree via front propagation.
 * Lazy binary min-heap Dijkstra from terminal seeds; when fronts from two
 * different terminals meet, mark the connecting path as seam and merge fronts.
 * =================================================================== */
typedef struct { double d; int32_t v; } HItem;
typedef struct { HItem *a; size_t n, cap; } Heap;
static void heap_push(Heap *h, double d, int32_t v){
    if (h->n >= h->cap) return;          /* safety: never overflow */
    size_t i=h->n++; h->a[i].d=d; h->a[i].v=v;
    while(i>0){ size_t p=(i-1)/2; if(h->a[p].d<=h->a[i].d)break; HItem t=h->a[p];h->a[p]=h->a[i];h->a[i]=t; i=p; }
}
static HItem heap_pop(Heap *h){
    HItem top=h->a[0]; h->a[0]=h->a[--h->n]; size_t i=0;
    for(;;){ size_t l=2*i+1,r=2*i+2,s=i;
        if(l<h->n&&h->a[l].d<h->a[s].d)s=l; if(r<h->n&&h->a[r].d<h->a[s].d)s=r;
        if(s==i)break; HItem t=h->a[s];h->a[s]=h->a[i];h->a[i]=t; i=s; }
    return top;
}
static int32_t uf_find2(int32_t *p, int32_t x){ while(p[x]!=x){p[x]=p[p[x]];x=p[x];} return x; }

/* candidate terminal-connection (Voronoi region-crossing edge) */
typedef struct { double w; int32_t tu, tw, u, vv; } SCand;
static int cmp_scand(const void *a, const void *b)
{ double x=((const SCand*)a)->w, y=((const SCand*)b)->w; return (x<y)?-1:(x>y?1:0); }

/* Mehlhorn's 2-approx Steiner tree (paper 6.2 steps 1-3): one multi-source
 * Dijkstra from all terminals -> per-vertex nearest terminal (Voronoi) + dist +
 * parent; then for each edge crossing two Voronoi regions a candidate connecting
 * the two terminals (weight = dist[u]+w(e)+dist[v]); Kruskal MST over terminals;
 * backtrack each chosen candidate to its two terminals to recover the seam path.
 * Boundary edges weigh EPS so the tree runs along boundaries for free (those
 * seam edges are already open and are ignored by the cut). */
static void steiner_tree(Arena_T arena, const Conn *c,
                         const char *is_terminal, SeamSet *seam)
{
    size_t nv=c->nv, nc=c->nc;
    const double EPS_BND = 1e-4;
    double  *dist = (double *)ARENA_ALLOC(arena,(long)(nv*sizeof(double)));
    int32_t *nrst = (int32_t *)ARENA_ALLOC(arena,(long)(nv*sizeof(int32_t)));
    int32_t *parent = (int32_t *)ARENA_ALLOC(arena,(long)(nv*sizeof(int32_t)));
    for (size_t v=0;v<nv;v++){ dist[v]=1e300; nrst[v]=-1; parent[v]=-1; }

    Heap h; h.cap=nc*4+1024; h.n=0;
    h.a=(HItem*)ARENA_ALLOC(arena,(long)(h.cap*sizeof(HItem)));
    for (size_t v=0;v<nv;v++) if(is_terminal[v]){ dist[v]=0; nrst[v]=(int32_t)v; heap_push(&h,0.0,(int32_t)v); }

    while (h.n>0){
        HItem it=heap_pop(&h); int32_t u=it.v;
        if (it.d>dist[u]+1e-12) continue;
        for (int32_t e=c->voff[u]; e<c->voff[u+1]; e++){
            int32_t cc=c->vcor[e];
            int32_t ws[2]={c->faces[NEXT(cc)], c->faces[PREV(cc)]};
            int32_t es[2]={cc, PREV(cc)};
            for (int t=0;t<2;t++){
                int32_t w=ws[t]; int is_bnd=(c->twin[es[t]]<0);
                double dd[3]; vsub(c->verts+3*(size_t)u, c->verts+3*(size_t)w, dd);
                double wt = is_bnd ? EPS_BND : vlen(dd);
                double nd = dist[u]+wt;
                if (nd < dist[w]-1e-12){ dist[w]=nd; nrst[w]=nrst[u]; parent[w]=u; heap_push(&h,nd,w); }
            }
        }
    }

    /* collect candidates (each undirected edge once) */
    Arena_Mark m = Arena_save(arena);
    SCand *cand = (SCand *)ARENA_ALLOC(arena,(long)(nc*sizeof(SCand)));
    long ncand=0;
    for (size_t cc=0; cc<nc; cc++){
        int32_t tw=c->twin[cc];
        if (tw>=0 && (int32_t)cc>tw) continue;        /* interior edge once */
        int32_t u=c->faces[cc], w=c->faces[NEXT(cc)];
        if (nrst[u]<0 || nrst[w]<0 || nrst[u]==nrst[w]) continue;
        double wt;
        if (tw<0) wt=EPS_BND;
        else { double dd[3]; vsub(c->verts+3*(size_t)u, c->verts+3*(size_t)w, dd); wt=vlen(dd); }
        cand[ncand].w=dist[u]+wt+dist[w];
        cand[ncand].tu=nrst[u]; cand[ncand].tw=nrst[w];
        cand[ncand].u=u; cand[ncand].vv=w; ncand++;
    }
    qsort(cand,(size_t)ncand,sizeof(SCand),cmp_scand);

    int32_t *uf = (int32_t *)ARENA_ALLOC(arena,(long)(nv*sizeof(int32_t)));
    for (size_t v=0;v<nv;v++) uf[v]=(int32_t)v;
    for (long i=0;i<ncand;i++){
        int32_t ru=uf_find2(uf,cand[i].tu), rw=uf_find2(uf,cand[i].tw);
        if (ru==rw) continue;
        uf[ru]=rw;
        int32_t x=cand[i].u;  while(parent[x]>=0){ seam_add(seam,x,parent[x]); x=parent[x]; }
        x=cand[i].vv;         while(parent[x]>=0){ seam_add(seam,x,parent[x]); x=parent[x]; }
        seam_add(seam, cand[i].u, cand[i].vv);
    }
    Arena_restore(arena, m);
}

/* ===================================================================
 * Cut surgery: open the mesh along seam edges + existing boundary.
 * Per vertex, incident faces split into sectors by cut edges; one new vertex
 * per sector. Drops unreferenced vertices automatically.
 * =================================================================== */
static int cut_surgery(Arena_T arena, const Conn *c, const SeamSet *seam,
                       SeamCutResult *out)
{
    size_t nc=c->nc, nv=c->nv, nf=c->nf;
    int32_t *sp = (int32_t *)ARENA_ALLOC(arena,(long)(nc*sizeof(int32_t)));
    for (size_t i=0;i<nc;i++) sp[i]=(int32_t)i;
    /* union corners sharing a non-cut edge at their common vertex */
    for (size_t cc=0; cc<nc; cc++){
        int32_t tw=c->twin[cc];
        if (tw<0) continue;
        int32_t v=c->faces[cc], w=c->faces[NEXT(cc)];
        if (seam_has(seam,v,w)) continue;        /* cut edge -> sector break */
        /* the v-corner of the twin face. Look it up by vertex rather than
         * assuming NEXT(tw): QEM can leave adjacent faces inconsistently wound
         * (both traversing the edge the same way), and NEXT(tw) would then land
         * on the wrong vertex and mis-merge the sectors. */
        int32_t tf = (tw/3)*3;
        int32_t d = (c->faces[tf]==v) ? tf : (c->faces[tf+1]==v) ? tf+1 : tf+2;
        /* union cc and d (both are corners at v) */
        int32_t ra=cc; while(sp[ra]!=ra)ra=sp[ra];
        int32_t rb=d;  while(sp[rb]!=rb)rb=sp[rb];
        if(ra!=rb) sp[ra]=rb;
    }
    /* assign new vertex ids per sector root */
    int32_t *newid = (int32_t *)ARENA_ALLOC(arena,(long)(nc*sizeof(int32_t)));
    for (size_t i=0;i<nc;i++) newid[i]=-1;
    int32_t nvnew=0;
    for (size_t cc=0; cc<nc; cc++){ int32_t r=cc; while(sp[r]!=r)r=sp[r];
        if(newid[r]<0) newid[r]=nvnew++; }
    /* build outputs */
    float *V = (float *)ARENA_ALLOC(arena,(long)((size_t)nvnew*3*sizeof(float)));
    int32_t *vmap = (int32_t *)ARENA_ALLOC(arena,(long)((size_t)nvnew*sizeof(int32_t)));
    int32_t *F = (int32_t *)ARENA_ALLOC(arena,(long)(nf*3*sizeof(int32_t)));
    for (size_t cc=0; cc<nc; cc++){
        int32_t r=cc; while(sp[r]!=r)r=sp[r];
        int32_t id=newid[r];
        int32_t v=c->faces[cc];
        V[id*3+0]=c->verts[3*(size_t)v+0]; V[id*3+1]=c->verts[3*(size_t)v+1]; V[id*3+2]=c->verts[3*(size_t)v+2];
        vmap[id]=v;
        F[cc]=id;
    }
    (void)nv;
    out->verts=V; out->faces=F; out->nv=(size_t)nvnew; out->nf=nf; out->vmap=vmap;
    return 0;
}

/* (Dr, idx) pair for terminal-selection sort (Algorithm 1, descending Dr). */
typedef struct { double d; int32_t i; } SeamDI;
static int cmp_di_desc(const void *a, const void *b)
{
    double x=((const SeamDI*)a)->d, y=((const SeamDI*)b)->d;
    return (x<y)?1:(x>y?-1:0);
}

/* ===================================================================
 * Section 6.3: genus reduction (Erickson & Har-Peled 2002, as adopted by
 * Seamster). Repeatedly cut the shortest *non-separating* loop -- one that
 * drops the genus by 1 instead of merely splitting the surface -- until the
 * genus is 0. Visibility is uniform for a buried scroll (paper's V(e)=1), so
 * the loop weight is pure edge length.
 *
 *   per loop:  shortest-path tree (Dijkstra) from a base vertex; each interior
 *              non-tree edge (a,b) closes a fundamental cycle
 *              path(a..lca) + (a,b) + path(b..lca); a cut is genus-reducing iff
 *              it leaves the surface connected (dual face graph stays one
 *              component). Cut the shortest such loop, re-open the mesh, repeat.
 *
 * Each cut re-opens the mesh by vertex duplication (cut_surgery), so the two
 * new banks become ordinary boundary loops that the later Steiner step joins.
 * =================================================================== */

/* BFS over a primal spanning tree (edges marked in tree[]) from src, filling
 * parent[], depth[] (hop count, for LCA) and dist[] (summed edge length, for
 * loop ranking). The tree spans every vertex, so all are reached. */
static void tree_bfs(const Conn *c, const char *tree, int32_t src,
                     double *dist, int32_t *parent, int32_t *depth, int32_t *queue)
{
    for (size_t v=0; v<c->nv; v++){ dist[v]=1e300; parent[v]=-1; depth[v]=-1; }
    size_t h=0, t=0; dist[src]=0.0; depth[src]=0; queue[t++]=src;
    while (h<t){
        int32_t u=queue[h++];
        for (int32_t e=c->voff[u]; e<c->voff[u+1]; e++){
            int32_t cc=c->vcor[e];
            int32_t ws[2]={c->faces[NEXT(cc)], c->faces[PREV(cc)]};
            int32_t es[2]={cc, PREV(cc)};        /* corner-id of each emanating edge */
            for (int k=0;k<2;k++){
                if (!tree[es[k]]) continue;
                int32_t w=ws[k];
                if (depth[w]!=-1) continue;      /* already in BFS */
                double dd[3]; vsub(c->verts+3*(size_t)u, c->verts+3*(size_t)w, dd);
                parent[w]=u; depth[w]=depth[u]+1; dist[w]=dist[u]+vlen(dd); queue[t++]=w;
            }
        }
    }
}

/* Add the fundamental cycle of non-tree edge (a,b) to seam set eg:
 * the two tree paths a..lca and b..lca plus edge (a,b). */
static void extract_cycle(const int32_t *parent, const int32_t *depth,
                          int32_t a, int32_t b, SeamSet *eg)
{
    seam_add(eg, a, b);
    int32_t x=a, y=b;
    while (depth[x]>depth[y]){ seam_add(eg,x,parent[x]); x=parent[x]; }
    while (depth[y]>depth[x]){ seam_add(eg,y,parent[y]); y=parent[y]; }
    while (x!=y){
        seam_add(eg,x,parent[x]); x=parent[x];
        seam_add(eg,y,parent[y]); y=parent[y];
    }
}

/* Non-separating test: with eg's edges cut, the surface stays connected iff the
 * dual face graph (linked across interior, non-cut edges) is one component. */
static int is_nonseparating(Arena_T arena, const Conn *c, const SeamSet *eg)
{
    if (c->nf == 0) return 0;
    Arena_Mark m=Arena_save(arena);
    int32_t *uf=(int32_t*)ARENA_ALLOC(arena,(long)(c->nf*sizeof(int32_t)));
    for (size_t f=0; f<c->nf; f++) uf[f]=(int32_t)f;
    for (size_t cc=0; cc<c->nc; cc++){
        int32_t tw=c->twin[cc];
        if (tw<0 || (int32_t)cc>tw) continue;        /* interior edge, once */
        int32_t a=c->faces[cc], b=c->faces[NEXT(cc)];
        if (seam_has(eg,a,b)) continue;              /* cut -> no dual link */
        int32_t ra=uf_find2(uf,(int32_t)(cc/3)), rb=uf_find2(uf,tw/3);
        if (ra!=rb) uf[ra]=rb;
    }
    int32_t r0=uf_find2(uf,0); int conn=1;
    for (size_t f=1; f<c->nf; f++) if (uf_find2(uf,(int32_t)f)!=r0){ conn=0; break; }
    Arena_restore(arena,m);
    return conn;
}

/* loop candidate: interior non-tree edge ranked by closed-loop length. */
typedef struct { double w; int32_t a, b; } LCand;
static int cmp_lcand(const void *p, const void *q)
{ double x=((const LCand*)p)->w, y=((const LCand*)q)->w; return (x<y)?-1:(x>y?1:0); }

/* Summed edge length of a loop edge-set (cv = current vertex coords). */
static double seam_loop_length(const SeamSet *eg, size_t cnv, const float *cv)
{
    double len=0.0;
    for (size_t si=0; si<eg->cap; si++){
        uint64_t k=eg->slot[si]; if (k==UINT64_MAX) continue;
        int32_t lo=(int32_t)(k/(uint64_t)cnv), hi=(int32_t)(k%(uint64_t)cnv);
        double dd[3]; vsub(cv+3*(size_t)lo, cv+3*(size_t)hi, dd); len+=vlen(dd);
    }
    return len;
}

/* Keep the largest EDGE-connected component (faces sharing an undirected edge).
 * A genus cut can leave a face attached only by a vertex; that island walls the
 * dual BFS and corrupts the genus read, so we strip it between cuts. Returns the
 * new mesh and a new->old vertex map. */
static void keep_largest_edge(Arena_T arena, const float *v, size_t nv,
                              const int32_t *f, size_t nf,
                              float **ov, size_t *onv, int32_t **of, size_t *onf,
                              int32_t **ovmap)
{
    size_t nc = 3*nf;
    EKey *ek = (EKey *)ARENA_ALLOC(arena, (long)(nc*sizeof(EKey)));
    for (size_t cc=0; cc<nc; cc++){
        int32_t a=f[cc], b=f[NEXT(cc)];
        uint64_t lo=(uint64_t)(a<b?a:b), hi=(uint64_t)(a<b?b:a);
        ek[cc].key=lo*(uint64_t)nv+hi; ek[cc].corner=(int32_t)cc;
    }
    qsort(ek, nc, sizeof(EKey), cmp_ekey);
    int32_t *uf=(int32_t*)ARENA_ALLOC(arena,(long)(nf*sizeof(int32_t)));
    for (size_t i=0;i<nf;i++) uf[i]=(int32_t)i;
    size_t i=0;
    while (i<nc){ size_t j=i; while(j<nc && ek[j].key==ek[i].key) j++;
        for (size_t p=i+1;p<j;p++){ int32_t ra=uf_find2(uf,ek[i].corner/3), rb=uf_find2(uf,ek[p].corner/3); if(ra!=rb) uf[ra]=rb; }
        i=j; }
    long *cnt=(long*)ARENA_CALLOC(arena,(long)nf,(long)sizeof(long));
    for (size_t t=0;t<nf;t++) cnt[uf_find2(uf,(int32_t)t)]++;
    int32_t best=-1; long bestn=-1;
    for (size_t r=0;r<nf;r++) if(cnt[r]>bestn){bestn=cnt[r];best=(int32_t)r;}

    int32_t *vm=(int32_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(int32_t)));
    int32_t *rm=(int32_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(int32_t)));
    for (size_t k=0;k<nv;k++) vm[k]=-1;
    int32_t *F=(int32_t*)ARENA_ALLOC(arena,bestn*3*(long)sizeof(int32_t));
    long onv2=0, fi=0;
    for (size_t t=0;t<nf;t++){
        if (uf_find2(uf,(int32_t)t)!=best) continue;
        for (int k=0;k<3;k++){ int32_t vg=f[t*3+k];
            if(vm[vg]==-1){vm[vg]=(int32_t)onv2; rm[onv2]=vg; onv2++;}
            F[fi*3+k]=vm[vg]; }
        fi++;
    }
    float *V=(float*)ARENA_ALLOC(arena,onv2*3*(long)sizeof(float));
    int32_t *vmap=(int32_t*)ARENA_ALLOC(arena,onv2*(long)sizeof(int32_t));
    for (long k=0;k<onv2;k++){ V[k*3+0]=v[rm[k]*3+0]; V[k*3+1]=v[rm[k]*3+1]; V[k*3+2]=v[rm[k]*3+2]; vmap[k]=rm[k]; }
    *ov=V; *onv=(size_t)onv2; *of=F; *onf=(size_t)fi; *ovmap=vmap;
}

/* Reduce the mesh genus to 0 by cutting non-separating loops. On return
 * *pv/*pnv/*pf/*pnf describe the genus-0 (multi-boundary) mesh, and *paccum
 * maps each new vertex back to an ORIGINAL vertex (NULL if no cut was made). */
static void genus_reduce(Arena_T arena, const float **pv, size_t *pnv,
                         const int32_t **pf, size_t *pnf,
                         int32_t **paccum, long genus0, long *ploops_cut,
                         HandleLoop *collect, size_t *ncol, double max_loop_len)
{
    const float   *cv=*pv; size_t cnv=*pnv;
    const int32_t *cf=*pf; size_t cnf=*pnf;
    int32_t *accum=NULL;                 /* current-vertex -> original-vertex */
    long cut=0, cap=genus0*4+16;         /* safety bound on the loop count */

    int dbg = (getenv("SEAMCUT_DEBUG") != NULL);
    for (;;){
        double t_it = ves_clock_sec();
        MeshTopoInfo ti; MeshTopo_analyze(arena, NULL, cnv, cf, cnf, &ti);
        long g=lround(ti.genus);
        if (g<=0) break;
        if (cut>=cap){ fprintf(stderr,"  [seamcut] genus-reduce cap %ld hit (g=%ld left)\n",cap,g); break; }

        Conn c; conn_build(arena, cv, cnv, cf, cnf, &c);

        /* Tree-cotree decomposition, COTREE FIRST. Build a dual spanning tree
         * over the face graph using ALL interior edges -- the dual of a
         * connected mesh is connected, so this reaches every face and marks
         * F-1 cotree edges. Then build the primal spanning tree over the
         * NON-cotree edges (cotree is acyclic in the dual, so primal-minus-
         * cotree stays connected and spans every vertex). The interior edges in
         * NEITHER tree are exactly the 2-chi = 2G + (B-1) homology generators --
         * the only non-contractible fundamental cycles. (Building the primal
         * tree first and excluding it from the dual disconnects the dual on a
         * bounded surface -- that was the dualF=292/3851, gens=3236 bug.) */
        char *cotree=(char*)ARENA_CALLOC(arena,(long)c.nc,1L);
        size_t dual_faces=0;
        {
            char    *fseen=(char*)ARENA_CALLOC(arena,(long)c.nf,1L);
            int32_t *fq=(int32_t*)ARENA_ALLOC(arena,(long)(c.nf*sizeof(int32_t)));
            size_t qh=0,qt=0; fseen[0]=1; fq[qt++]=0;
            while (qh<qt){
                int32_t f=fq[qh++];
                for (int k=0;k<3;k++){
                    int32_t cc=f*3+k, tw=c.twin[cc];
                    if (tw<0) continue;                          /* boundary: no dual edge */
                    int32_t nf2=tw/3;
                    if (fseen[nf2]) continue;
                    fseen[nf2]=1; cotree[cc]=1; cotree[tw]=1; fq[qt++]=nf2;
                }
            }
            dual_faces=qt;
        }
        /* Primal spanning tree, completed so it spans EVERY vertex (on a bounded
         * surface the cotree can otherwise wall off a vertex): union-find takes
         * non-cotree edges first, then promotes cotree edges to reconnect any
         * remaining components (those move out of the cotree). Tree and cotree
         * stay disjoint; generators are the edges in neither. */
        char *tree=(char*)ARENA_CALLOC(arena,(long)c.nc,1L);
        {
            int32_t *uf=(int32_t*)ARENA_ALLOC(arena,(long)(cnv*sizeof(int32_t)));
            for (size_t v=0;v<cnv;v++) uf[v]=(int32_t)v;
            /* pass A1: BOUNDARY edges first -- absorb the boundary into the
             * spanning tree (Seamster's eps-weight boundary trick) so the
             * surviving interior generators are exactly the 2G handle loops, not
             * boundary-to-boundary arcs the interior-only enumeration would miss
             * (that left the last handle stranded). */
            for (size_t cc=0; cc<c.nc; cc++){
                if (c.twin[cc]>=0) continue;                 /* boundary corner */
                int32_t ra=uf_find2(uf,c.faces[cc]), rb=uf_find2(uf,c.faces[NEXT(cc)]);
                if (ra!=rb){ uf[ra]=rb; tree[cc]=1; }
            }
            /* pass A2: interior non-cotree edges. */
            for (size_t cc=0; cc<c.nc; cc++){
                int32_t tw=c.twin[cc];
                if (tw<0 || (int32_t)cc>tw) continue;        /* interior once */
                if (cotree[cc]) continue;
                int32_t ra=uf_find2(uf,c.faces[cc]), rb=uf_find2(uf,c.faces[NEXT(cc)]);
                if (ra!=rb){ uf[ra]=rb; tree[cc]=1; tree[tw]=1; }
            }
            /* pass B: promote cotree edges to finish spanning every vertex. */
            for (size_t cc=0; cc<c.nc; cc++){
                int32_t tw=c.twin[cc];
                if (tw<0 || (int32_t)cc>tw) continue;        /* interior cotree edges */
                if (!cotree[cc]) continue;
                int32_t ra=uf_find2(uf,c.faces[cc]), rb=uf_find2(uf,c.faces[NEXT(cc)]);
                if (ra!=rb){ uf[ra]=rb; tree[cc]=1; tree[tw]=1; cotree[cc]=0; cotree[tw]=0; }
            }
        }
        double  *dist=(double *)ARENA_ALLOC(arena,(long)(cnv*sizeof(double)));
        int32_t *par =(int32_t*)ARENA_ALLOC(arena,(long)(cnv*sizeof(int32_t)));
        int32_t *dep =(int32_t*)ARENA_ALLOC(arena,(long)(cnv*sizeof(int32_t)));
        int32_t *bq =(int32_t*)ARENA_ALLOC(arena,(long)(cnv*sizeof(int32_t)));
        tree_bfs(&c, tree, 0, dist, par, dep, bq);
        /* rank the generators (non-tree, non-cotree) by fundamental-loop length */
        LCand *cand=(LCand*)ARENA_ALLOC(arena,(long)(c.nc*sizeof(LCand)));
        long ncd=0;
        for (size_t cc=0; cc<c.nc; cc++){
            int32_t tw=c.twin[cc];
            if (tw<0 || (int32_t)cc>tw) continue;            /* interior edge once */
            if (tree[cc] || cotree[cc]) continue;            /* in a tree -> not a generator */
            int32_t a=c.faces[cc], b=c.faces[NEXT(cc)];
            if (dist[a]>=1e299 || dist[b]>=1e299) continue;  /* unreached (shouldn't happen) */
            double dd[3]; vsub(cv+3*(size_t)a, cv+3*(size_t)b, dd);
            cand[ncd].w=dist[a]+vlen(dd)+dist[b]; cand[ncd].a=a; cand[ncd].b=b; ncd++;
        }
        qsort(cand,(size_t)ncd,sizeof(LCand),cmp_lcand);

        int found=0; long tested=0, skipped_long=0;
        for (long i=0;i<ncd;i++){
            Arena_Mark mk=Arena_save(arena);
            SeamSet eg; seam_init(arena,&eg,cnv,2*cnv);
            extract_cycle(par, dep, cand[i].a, cand[i].b, &eg);
            tested++;
            if (is_nonseparating(arena,&c,&eg)){
                /* Length bound: skip non-separating loops longer than the cap so
                 * only SHORT handles (thin BPA bridges) are cut, never scroll-scale
                 * loops. Ranking is by fundamental-cycle length, so we test the
                 * ACTUAL loop length here and keep scanning for a short one. */
                double loop_len = seam_loop_length(&eg, cnv, cv);
                if (max_loop_len > 0.0 && loop_len > max_loop_len){
                    skipped_long++;
                    Arena_restore(arena,mk);
                    continue;
                }
                /* Record the loop (in ORIGINAL/input vertex space via accum)
                 * before the cut renumbers vertices. */
                if (collect){
                    size_t ne=0;
                    for (size_t si=0; si<eg.cap; si++) if (eg.slot[si]!=UINT64_MAX) ne++;
                    int32_t *lv=(int32_t*)ARENA_ALLOC(arena,(long)(2*ne*sizeof(int32_t)));
                    size_t kk=0; double len=0.0;
                    for (size_t si=0; si<eg.cap; si++){
                        uint64_t key=eg.slot[si]; if (key==UINT64_MAX) continue;
                        int32_t lo=(int32_t)(key/(uint64_t)cnv), hi=(int32_t)(key%(uint64_t)cnv);
                        double dd[3]; vsub(cv+3*(size_t)lo, cv+3*(size_t)hi, dd); len+=vlen(dd);
                        lv[kk++]= accum?accum[lo]:lo;
                        lv[kk++]= accum?accum[hi]:hi;
                    }
                    collect[*ncol].verts=lv; collect[*ncol].n=kk; collect[*ncol].length=len;
                    (*ncol)++;
                }
                SeamCutResult tmp; memset(&tmp,0,sizeof tmp);
                cut_surgery(arena,&c,&eg,&tmp);
                int32_t *na=(int32_t*)ARENA_ALLOC(arena,(long)(tmp.nv*sizeof(int32_t)));
                for (size_t k=0;k<tmp.nv;k++) na[k]=accum?accum[tmp.vmap[k]]:tmp.vmap[k];
                cv=tmp.verts; cnv=tmp.nv; cf=tmp.faces; cnf=tmp.nf; accum=na;
                /* strip any vertex-only island the cut created, so the next
                 * iteration's dual BFS spans and the genus read stays honest. */
                {
                    float *lv=NULL; int32_t *lf=NULL; size_t lnv=0,lnf=0; int32_t *lvm=NULL;
                    keep_largest_edge(arena, cv,cnv,cf,cnf, &lv,&lnv,&lf,&lnf,&lvm);
                    if (lnv>0 && lnv<cnv){
                        int32_t *na2=(int32_t*)ARENA_ALLOC(arena,(long)(lnv*sizeof(int32_t)));
                        for (size_t k=0;k<lnv;k++) na2[k]=accum[lvm[k]];
                        cv=lv; cnv=lnv; cf=lf; cnf=lnf; accum=na2;
                    }
                }
                cut++; found=1;
                break;                            /* keep arena: tmp/na must live */
            }
            Arena_restore(arena,mk);              /* discard eg, try next candidate */
        }
        if (dbg) fprintf(stderr,"  [seamcut] g=%ld B=%ld nm=%ld V=%zu gens=%ld dualF=%zu/%zu tested=%ld -> %s (%.2fs)\n",
                         g, ti.n_boundary_loops, ti.n_nonmanifold_edges, cnv, ncd,
                         dual_faces, c.nf, tested, found?"cut":"none", ves_clock_sec()-t_it);
        if (!found){
            /* Clean stop: every remaining non-separating loop is longer than the
             * cap -- only scroll-scale handles are left, which we intentionally
             * keep. Not a stall. */
            if (max_loop_len > 0.0 && skipped_long > 0) break;
            if (dbg){
                long ne_int=0, nt=0, nco=0, nboth=0;
                for (size_t cc=0; cc<c.nc; cc++){
                    int32_t tw=c.twin[cc];
                    if (tw<0 || (int32_t)cc>tw) continue;
                    ne_int++; if(tree[cc])nt++; if(cotree[cc])nco++; if(tree[cc]&&cotree[cc])nboth++;
                }
                fprintf(stderr,"  [seamcut] STALL g=%ld E_int=%ld tree=%ld cotree=%ld both=%ld unref=%ld V=%zu F=%zu\n",
                        g, ne_int, nt, nco, nboth, ti.n_unref_verts, cnv, c.nf);
            }
            fprintf(stderr,"  [seamcut] no genus-reducing loop found (g=%ld)\n",g); break;
        }
    }
    *pv=cv; *pnv=cnv; *pf=cf; *pnf=cnf; *paccum=accum; *ploops_cut=cut;
}

/* Enumerate the handle loops (one per genus) without modifying the input. */
int SeamCut_handle_loops(Arena_T arena, const float *verts, size_t nv,
                         const int32_t *faces, size_t nf,
                         HandleLoop **out_loops, size_t *out_n)
{
    assert(arena && out_loops && out_n);
    *out_loops = NULL; *out_n = 0;
    if (nv < 3 || nf < 1 || !verts || !faces) return -1;

    MeshTopoInfo ti; MeshTopo_analyze(arena, NULL, nv, faces, nf, &ti);
    long g0 = lround(ti.genus);
    if (g0 <= 0) return 0;                       /* no handles -> no loops */

    HandleLoop *loops = (HandleLoop *)ARENA_CALLOC(arena, (long)g0,
                                                   (long)sizeof(HandleLoop));
    size_t ncol = 0;
    const float   *cv = verts; size_t cnv = nv;
    const int32_t *cf = faces; size_t cnf = nf;
    int32_t *accum = NULL; long loops_cut = 0;
    genus_reduce(arena, &cv, &cnv, &cf, &cnf, &accum, g0, &loops_cut, loops, &ncol, 0.0);

    *out_loops = loops;
    *out_n = ncol;
    return 0;
}

/* Sever only the SHORT handles: open every non-separating loop shorter than
 * max_loop_len (the thin BPA bridges), leaving longer scroll-scale handles
 * intact. Reuses the tested genus reducer (cut_surgery is manifold-preserving),
 * just length-bounded. *out_* describe the opened mesh; *out_loops_cut counts
 * the handles severed. Outputs arena-allocated. Returns 0; -1 on bad input. */
int SeamCut_sever_short_handles(Arena_T arena, const float *verts, size_t nv,
                                const int32_t *faces, size_t nf, double max_loop_len,
                                float **out_verts, size_t *out_nv,
                                int32_t **out_faces, size_t *out_nf,
                                long *out_loops_cut)
{
    assert(arena && out_verts && out_nv && out_faces && out_nf);
    if (nv < 3 || nf < 1 || !verts || !faces) return -1;

    MeshTopoInfo ti; MeshTopo_analyze(arena, NULL, nv, faces, nf, &ti);
    long g0 = lround(ti.genus);

    const float   *cv = verts; size_t cnv = nv;
    const int32_t *cf = faces; size_t cnf = nf;
    int32_t *accum = NULL; long loops_cut = 0;
    if (g0 > 0)
        genus_reduce(arena, &cv, &cnv, &cf, &cnf, &accum, g0, &loops_cut, NULL, NULL,
                     max_loop_len);

    /* Copy out (cv/cf may alias the inputs if nothing was cut). */
    float   *ov = (float *)  ARENA_ALLOC(arena, (long)(cnv * 3 * sizeof(float)));
    int32_t *of = (int32_t *)ARENA_ALLOC(arena, (long)(cnf * 3 * sizeof(int32_t)));
    memcpy(ov, cv, cnv * 3 * sizeof(float));
    memcpy(of, cf, cnf * 3 * sizeof(int32_t));
    *out_verts = ov; *out_nv = cnv;
    *out_faces = of; *out_nf = cnf;
    if (out_loops_cut) *out_loops_cut = loops_cut;
    return 0;
}

/* ===================================================================
 * Manifold repair (Seamster precondition): split non-manifold vertices/edges.
 * cut_surgery with an EMPTY seam groups each vertex's corners into fans by
 * shared (twin-paired) edges; a bowtie vertex has >1 fan -> split, and a
 * non-manifold edge (twin=-1, left unpaired by conn_build) is never glued ->
 * its faces separate. Clean manifolds pass through unchanged.
 * =================================================================== */
int SeamCut_repair_manifold(Arena_T arena, const float *verts, size_t nv,
                            const int32_t *faces, size_t nf,
                            float **out_verts, size_t *out_nv,
                            int32_t **out_faces, size_t *out_nf,
                            int32_t **out_vmap)
{
    assert(arena && out_verts && out_nv && out_faces && out_nf && out_vmap);
    if (nv < 3 || nf < 1 || !verts || !faces) return -1;

    /* Step 1: drop DEGENERATE faces (a repeated vertex, e.g. QEM collapse
     * artifacts -- their self-edges corrupt the topology) and faces on
     * NON-MANIFOLD EDGES (>2 incident faces -- the vertex-sector split below
     * can't separate them, and they wall off the dual BFS, breaking genus
     * reduction). The small holes left behind are rejoined by the Steiner step.*/
    size_t nc = 3*nf;
    char *drop = (char *)ARENA_CALLOC(arena, (long)nf, 1L);
    long ndrop = 0, ndegen = 0;
    for (size_t t = 0; t < nf; t++) {
        int32_t a = faces[t*3+0], b = faces[t*3+1], c = faces[t*3+2];
        if (a == b || b == c || c == a) { drop[t] = 1; ndrop++; ndegen++; }
    }
    {
        Arena_Mark m = Arena_save(arena);
        EKey *ek = (EKey *)ARENA_ALLOC(arena, (long)(nc * sizeof(EKey)));
        size_t nek = 0;
        for (size_t cc = 0; cc < nc; cc++) {
            if (drop[cc/3]) continue;                    /* skip degenerate faces */
            int32_t a = faces[cc], b = faces[NEXT(cc)];
            uint64_t lo = (uint64_t)(a < b ? a : b), hi = (uint64_t)(a < b ? b : a);
            ek[nek].key = lo*(uint64_t)nv + hi; ek[nek].corner = (int32_t)cc; nek++;
        }
        qsort(ek, nek, sizeof(EKey), cmp_ekey);
        size_t i = 0;
        while (i < nek) {
            size_t j = i; while (j < nek && ek[j].key == ek[i].key) j++;
            if (j - i > 2)                                /* non-manifold edge */
                for (size_t p = i; p < j; p++) { int32_t f = ek[p].corner/3;
                    if (!drop[f]) { drop[f] = 1; ndrop++; } }
            i = j;
        }
        Arena_restore(arena, m);
    }
    const int32_t *use_faces = faces; size_t use_nf = nf;
    if (ndrop > 0) {
        fprintf(stderr, "  [seamcut] manifold repair: dropped %ld face(s) (%ld degenerate, %ld on non-manifold edges)\n",
                ndrop, ndegen, ndrop - ndegen);
        int32_t *ff = (int32_t *)ARENA_ALLOC(arena, (long)((nf - (size_t)ndrop) * 3 * sizeof(int32_t)));
        size_t fnf = 0;
        for (size_t t = 0; t < nf; t++) if (!drop[t]) {
            ff[fnf*3+0]=faces[t*3+0]; ff[fnf*3+1]=faces[t*3+1]; ff[fnf*3+2]=faces[t*3+2]; fnf++;
        }
        use_faces = ff; use_nf = fnf;
        if (use_nf < 1) return -1;
    }

    /* Step 2: split non-manifold VERTICES (bowties) via cut_surgery with an
     * empty seam -- each vertex's corners group into fans by shared edges. */
    Conn c; conn_build(arena, verts, nv, use_faces, use_nf, &c);
    SeamSet seam; seam_init(arena, &seam, nv, 16);   /* empty: nothing extra cut */
    SeamCutResult tmp; memset(&tmp, 0, sizeof tmp);
    if (cut_surgery(arena, &c, &seam, &tmp) != 0) return -1;
    *out_verts = tmp.verts; *out_nv = tmp.nv;
    *out_faces = tmp.faces; *out_nf = tmp.nf; *out_vmap = tmp.vmap;
    return 0;
}

/* ===================================================================
 * Orchestrator.
 * =================================================================== */
int SeamCut_run(Arena_T arena, const float *verts, size_t nv,
                const int32_t *faces, size_t nf,
                const SeamCutOpts *opts, SeamCutResult *out)
{
    assert(arena && out);
    memset(out, 0, sizeof(*out));
    if (nv < 3 || nf < 1 || !verts || !faces) return -1;
    double accept = opts ? opts->accept_frac : 0.0;
    int r_max = (opts && opts->r_max>0) ? opts->r_max : 5;

    Arena_Mark mark = Arena_save(arena);

    MeshTopoInfo ti0; MeshTopo_analyze(arena, NULL, nv, faces, nf, &ti0);
    out->loops_before = ti0.n_boundary_loops;
    out->genus_before = lround(ti0.genus);

    int dbg = (getenv("SEAMCUT_DEBUG") != NULL);

    /* Section 6.3: genus reduction -> genus-0 (multi-boundary) mesh. The new
     * banks become ordinary boundary that the Steiner step below joins. */
    const float   *cv = verts; size_t cnv = nv;
    const int32_t *cf = faces; size_t cnf = nf;
    int32_t *accum = NULL; long loops_cut = 0;
    double tph = ves_clock_sec();
    if (out->genus_before > 0) {
        genus_reduce(arena, &cv, &cnv, &cf, &cnf, &accum, out->genus_before, &loops_cut,
                     NULL, NULL, 0.0);
        if (loops_cut > 0)
            fprintf(stderr, "  [seamcut] genus %ld -> 0 via %ld non-separating loop cut(s)\n",
                    out->genus_before, loops_cut);
    }
    if (dbg){ fprintf(stderr,"  [seamcut] phase genus_reduce %.2fs (nv %zu->%zu)\n", ves_clock_sec()-tph, nv, cnv); tph=ves_clock_sec(); }

    Conn c; conn_build(arena, cv, cnv, cf, cnf, &c);

    /* Section 5: distortion */
    double  *Dr = (double *)ARENA_ALLOC(arena,(long)(cnv*sizeof(double)));
    int32_t *br = (int32_t *)ARENA_ALLOC(arena,(long)(cnv*sizeof(int32_t)));
    compute_distortion(arena, &c, r_max, Dr, br);
    if (dbg){ fprintf(stderr,"  [seamcut] phase distortion %.2fs\n", ves_clock_sec()-tph); tph=ves_clock_sec(); }

    /* boundary vertices (incident to a boundary edge) */
    char *is_bnd = (char *)ARENA_CALLOC(arena,(long)cnv,1L);
    for (size_t cc=0; cc<c.nc; cc++) if (c.twin[cc]<0){ is_bnd[c.faces[cc]]=1; is_bnd[c.faces[NEXT(cc)]]=1; }

    /* Section 6.1: terminal selection (Algorithm 1) + all boundary verts */
    double Dtot=0.0; for (size_t v=0;v<cnv;v++) if(Dr[v]>0.0) Dtot+=Dr[v];
    double Dacc = accept*Dtot;          /* leave this much distortion uncut */
    /* sort vertices by Dr desc */
    int32_t *order=(int32_t*)ARENA_ALLOC(arena,(long)(cnv*sizeof(int32_t)));
    for (size_t v=0;v<cnv;v++) order[v]=(int32_t)v;
    {
        SeamDI *di=(SeamDI*)ARENA_ALLOC(arena,(long)(cnv*sizeof(SeamDI)));
        for(size_t v=0;v<cnv;v++){di[v].d=Dr[v];di[v].i=(int32_t)v;}
        qsort(di, cnv, sizeof(SeamDI), cmp_di_desc);
        for(size_t v=0;v<cnv;v++) order[v]=di[v].i;
    }
    char *is_term=(char*)ARENA_CALLOC(arena,(long)cnv,1L);
    long nterm=0;
    double s=0.0;
    for (size_t k=0; k<cnv; k++){
        int32_t v=order[k];
        if (Dr[v]<=0.0) break;
        if (s >= Dtot - Dacc) break;
        is_term[v]=1; nterm++; s+=Dr[v];
    }
    for (size_t v=0;v<cnv;v++) if(is_bnd[v] && !is_term[v]){ is_term[v]=1; nterm++; }
    out->n_terminals=nterm;

    /* Section 6.2: Steiner tree connecting terminals (boundary loops + high-D) */
    SeamSet seam; seam_init(arena, &seam, cnv, c.nc);
    if (nterm >= 2) steiner_tree(arena, &c, is_term, &seam);
    out->seam_edges=(long)seam.n;
    if (dbg){
        long si=0, sb=0;
        for (size_t cc=0; cc<c.nc; cc++){
            int32_t tw=c.twin[cc];
            if (tw>=0 && (int32_t)cc>tw) continue;
            if (!seam_has(&seam, c.faces[cc], c.faces[NEXT(cc)])) continue;
            if (tw<0) sb++; else si++;
        }
        fprintf(stderr,"  [seamcut] phase steiner %.2fs (nterm=%ld seam=%ld: interior=%ld boundary=%ld)\n",
                ves_clock_sec()-tph, nterm, (long)seam.n, si, sb); tph=ves_clock_sec();
    }

    /* cut surgery (final open into a single-boundary disk) */
    if (cut_surgery(arena, &c, &seam, out)!=0){ Arena_restore(arena,mark); return -1; }
    if (dbg){ fprintf(stderr,"  [seamcut] phase surgery %.2fs\n", ves_clock_sec()-tph); tph=ves_clock_sec(); }

    /* Compose the vmap back to ORIGINAL vertices through the genus-reduction
     * map: out->vmap currently indexes the (possibly genus-reduced) mesh. */
    if (accum) for (size_t i=0;i<out->nv;i++) out->vmap[i]=accum[out->vmap[i]];

    /* re-analyze. out arrays are arena-allocated and must survive, so DO NOT
     * restore the arena -- the caller disposes it later. */
    MeshTopoInfo ti1; MeshTopo_analyze(arena, NULL, out->nv, out->faces, out->nf, &ti1);
    out->loops_after = ti1.n_boundary_loops;
    out->genus_after = lround(ti1.genus);
    out->is_disk_after = ti1.is_disk;
    (void)mark;  /* out arrays live in the arena; caller disposes it later */
    return 0;
}

/* ===================================================================
 * Self-test.
 * =================================================================== */
static void build_grid(Arena_T arena, int nu, int nh, int wrap, double R,
                       float **out_v, size_t *out_nv, int32_t **out_f, size_t *out_nf)
{
    size_t nvv=(size_t)nu*(size_t)nh; int uc=wrap?nu:nu-1;
    size_t nff=(size_t)uc*(size_t)(nh-1)*2;
    float *v=(float*)ARENA_ALLOC(arena,(long)(nvv*3*sizeof(float)));
    int32_t *f=(int32_t*)ARENA_ALLOC(arena,(long)(nff*3*sizeof(int32_t)));
    for(int j=0;j<nh;j++)for(int i=0;i<nu;i++){
        double th=2.0*M_PI*(double)i/(double)nu, z=2.0*(double)j;
        size_t idx=(size_t)j*(size_t)nu+(size_t)i;
        v[idx*3+0]=(float)z; v[idx*3+1]=(float)(R*sin(th)); v[idx*3+2]=(float)(R*cos(th));
    }
    size_t fi=0;
    for(int j=0;j<nh-1;j++)for(int i=0;i<uc;i++){
        int i2=wrap?(i+1)%nu:i+1;
        int32_t a=(int32_t)((size_t)j*(size_t)nu+(size_t)i), b=(int32_t)((size_t)j*(size_t)nu+(size_t)i2);
        int32_t cc=(int32_t)((size_t)(j+1)*(size_t)nu+(size_t)i), d=(int32_t)((size_t)(j+1)*(size_t)nu+(size_t)i2);
        f[fi*3+0]=a;f[fi*3+1]=b;f[fi*3+2]=cc;fi++; f[fi*3+0]=b;f[fi*3+1]=d;f[fi*3+2]=cc;fi++;
    }
    *out_v=v;*out_nv=nvv;*out_f=f;*out_nf=fi;
}

/* Closed torus: nu around the tube, nr around the ring (both wrapped) -> genus
 * 1, no boundary. skip_quad=1 removes one quad to leave a single boundary loop
 * (genus 1 + 1 boundary). */
static void build_torus(Arena_T arena, int nu, int nr, double Rb, double rt,
                        int skip_quad, float **out_v, size_t *out_nv,
                        int32_t **out_f, size_t *out_nf)
{
    size_t NV=(size_t)nu*(size_t)nr;
    float *v=(float*)ARENA_ALLOC(arena,(long)(NV*3*sizeof(float)));
    int32_t *f=(int32_t*)ARENA_ALLOC(arena,(long)((size_t)nu*(size_t)nr*2*3*sizeof(int32_t)));
    for(int j=0;j<nr;j++)for(int i=0;i<nu;i++){
        double phi=2.0*M_PI*(double)i/(double)nu, th=2.0*M_PI*(double)j/(double)nr;
        size_t idx=(size_t)j*(size_t)nu+(size_t)i;
        v[idx*3+0]=(float)((Rb+rt*cos(phi))*cos(th));
        v[idx*3+1]=(float)((Rb+rt*cos(phi))*sin(th));
        v[idx*3+2]=(float)(rt*sin(phi));
    }
    size_t fi=0;
    for(int j=0;j<nr;j++)for(int i=0;i<nu;i++){
        if(skip_quad && i==0 && j==0) continue;       /* leave a hole */
        int i2=(i+1)%nu, j2=(j+1)%nr;
        int32_t a=(int32_t)((size_t)j *(size_t)nu+(size_t)i);
        int32_t b=(int32_t)((size_t)j *(size_t)nu+(size_t)i2);
        int32_t cc=(int32_t)((size_t)j2*(size_t)nu+(size_t)i);
        int32_t d=(int32_t)((size_t)j2*(size_t)nu+(size_t)i2);
        f[fi*3+0]=a;f[fi*3+1]=b;f[fi*3+2]=cc;fi++;
        f[fi*3+0]=b;f[fi*3+1]=d;f[fi*3+2]=cc;fi++;
    }
    *out_v=v;*out_nv=NV;*out_f=f;*out_nf=fi;
}

int SeamCut_selftest(void)
{
    int fails=0; Arena_T arena=Arena_new();
    /* Annulus (closed cylinder strip): 2 boundary loops -> cut -> 1 loop disk. */
    {
        float *v; int32_t *f; size_t nvv,nff;
        build_grid(arena, 24, 10, 1, 10.0, &v,&nvv,&f,&nff);
        SeamCutOpts o={0.0,5}; SeamCutResult res;
        int rc=SeamCut_run(arena, v,nvv,f,nff,&o,&res);
        int ok=(rc==0) && res.loops_before==2 && res.loops_after==1 && res.is_disk_after;
        fprintf(stderr,"  [seamcut annulus] rc=%d loops %ld->%ld seam=%ld disk=%d nv %zu->%zu -> %s\n",
                rc,res.loops_before,res.loops_after,res.seam_edges,res.is_disk_after,nvv,res.nv, ok?"ok":"FAIL");
        fails+=ok?0:1;
    }
    /* Disk (open strip): already 1 loop -> stays a disk. */
    {
        float *v; int32_t *f; size_t nvv,nff;
        build_grid(arena, 16, 10, 0, 10.0, &v,&nvv,&f,&nff);
        SeamCutOpts o={0.0,5}; SeamCutResult res;
        int rc=SeamCut_run(arena, v,nvv,f,nff,&o,&res);
        int ok=(rc==0) && res.loops_after==1 && res.is_disk_after;
        fprintf(stderr,"  [seamcut disk] rc=%d loops %ld->%ld disk=%d -> %s\n",
                rc,res.loops_before,res.loops_after,res.is_disk_after, ok?"ok":"FAIL");
        fails+=ok?0:1;
    }
    /* Torus (genus 1, closed): genus-reduce + Steiner -> single-loop disk. */
    {
        float *v; int32_t *f; size_t nvv,nff;
        build_torus(arena, 20, 14, 10.0, 3.0, 0, &v,&nvv,&f,&nff);
        SeamCutOpts o={0.0,5}; SeamCutResult res;
        int rc=SeamCut_run(arena, v,nvv,f,nff,&o,&res);
        int vmap_ok=(rc==0);
        if (rc==0) for (size_t k=0;k<res.nv;k++)
            if (res.vmap[k]<0 || (size_t)res.vmap[k]>=nvv){ vmap_ok=0; break; }
        int ok=(rc==0) && res.genus_before==1 && res.genus_after==0
                && res.loops_after==1 && res.is_disk_after && vmap_ok;
        fprintf(stderr,"  [seamcut torus] rc=%d genus %ld->%ld loops ->%ld disk=%d vmap=%d nv %zu->%zu -> %s\n",
                rc,res.genus_before,res.genus_after,res.loops_after,res.is_disk_after,vmap_ok,nvv,res.nv, ok?"ok":"FAIL");
        fails+=ok?0:1;
    }
    /* Torus with a hole (genus 1 + 1 boundary loop): -> single-loop disk. */
    {
        float *v; int32_t *f; size_t nvv,nff;
        build_torus(arena, 20, 14, 10.0, 3.0, 1, &v,&nvv,&f,&nff);
        SeamCutOpts o={0.0,5}; SeamCutResult res;
        int rc=SeamCut_run(arena, v,nvv,f,nff,&o,&res);
        int ok=(rc==0) && res.genus_before==1 && res.genus_after==0
                && res.loops_after==1 && res.is_disk_after;
        fprintf(stderr,"  [seamcut torus+hole] rc=%d genus %ld->%ld loops %ld->%ld disk=%d -> %s\n",
                rc,res.genus_before,res.genus_after,res.loops_before,res.loops_after,res.is_disk_after, ok?"ok":"FAIL");
        fails+=ok?0:1;
    }
    /* Handle-loop enumeration: a genus-1 torus has exactly ONE handle loop, with
     * vertices in the input index space. A genus-0 disk has none. */
    {
        float *v; int32_t *f; size_t nvv,nff;
        build_torus(arena, 20, 14, 10.0, 3.0, 0, &v,&nvv,&f,&nff);
        HandleLoop *loops=NULL; size_t nl=0;
        int rc=SeamCut_handle_loops(arena, v,nvv,f,nff, &loops,&nl);
        int idx_ok = (rc==0 && nl>=1 && loops[0].n>0);
        if (idx_ok) for (size_t k=0;k<loops[0].n;k++)
            if (loops[0].verts[k]<0 || (size_t)loops[0].verts[k]>=nvv){ idx_ok=0; break; }
        int ok=(rc==0 && nl==1 && idx_ok && loops[0].length>0.0);
        fprintf(stderr,"  [seamcut handle-loops torus] rc=%d nloops=%zu len=%.2f verts=%zu -> %s\n",
                rc, nl, nl?loops[0].length:0.0, nl?loops[0].n:0, ok?"ok":"FAIL");
        fails+=ok?0:1;

        float *vd; int32_t *fd; size_t nvd,nfd;
        build_grid(arena, 16, 10, 0, 10.0, &vd,&nvd,&fd,&nfd);
        HandleLoop *ld=NULL; size_t nld=0;
        int rc2=SeamCut_handle_loops(arena, vd,nvd,fd,nfd, &ld,&nld);
        int ok2=(rc2==0 && nld==0);
        fprintf(stderr,"  [seamcut handle-loops disk] rc=%d nloops=%zu -> %s\n",
                rc2, nld, ok2?"ok":"FAIL");
        fails+=ok2?0:1;
    }
    /* Length-bounded severing: a genus-1 torus (handle loop ~54 vox) is cut to
     * genus 0 with a generous cap, but LEFT INTACT with a cap below the loop
     * length (the false-positive guard against slicing long/legit handles). */
    {
        float *v; int32_t *f; size_t nvv,nff;
        build_torus(arena, 20, 14, 10.0, 3.0, 0, &v,&nvv,&f,&nff);
        float *ov; int32_t *of; size_t onv,onf; long cut=0;
        int rc=SeamCut_sever_short_handles(arena, v,nvv,f,nff, 1000.0,
                                           &ov,&onv,&of,&onf,&cut);
        MeshTopoInfo ti; ti.genus=-1.0;
        if (rc==0) MeshTopo_analyze(arena, NULL, onv, of, onf, &ti);
        int ok=(rc==0 && cut==1 && lround(ti.genus)==0);
        fprintf(stderr,"  [seamcut sever short cap=1000] rc=%d cut=%ld genus->%ld -> %s\n",
                rc, cut, rc==0?lround(ti.genus):-1, ok?"ok":"FAIL");
        fails+=ok?0:1;

        float *ov2; int32_t *of2; size_t onv2,onf2; long cut2=0;
        int rc3=SeamCut_sever_short_handles(arena, v,nvv,f,nff, 10.0,
                                            &ov2,&onv2,&of2,&onf2,&cut2);
        MeshTopoInfo ti2; ti2.genus=-1.0;
        if (rc3==0) MeshTopo_analyze(arena, NULL, onv2, of2, onf2, &ti2);
        int ok3=(rc3==0 && cut2==0 && lround(ti2.genus)==1);
        fprintf(stderr,"  [seamcut sever short cap=10] rc=%d cut=%ld genus->%ld (expect kept) -> %s\n",
                rc3, cut2, rc3==0?lround(ti2.genus):-1, ok3?"ok":"FAIL");
        fails+=ok3?0:1;
    }
    Arena_dispose(&arena);
    return fails;
}
