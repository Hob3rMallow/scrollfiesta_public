/* ============================================================================
 * overlap_quilt.c -- see overlap_quilt.h.
 *
 * Stage 1 (this file's core): detect UV overlaps, group them into regions, and
 * split each region's co-located material into connected LAYERS with a lifted
 * multicut (reusing LiftedMulticut_kernighan_lin -- the same solver the cube
 * mesher's OverlapSep uses). The one net-new signal vs OverlapSep: overlap
 * edges are gated by a RADIAL gap (two faces at the same UV but > radius_gate
 * apart from the axis are distinct layers; a single folded sheet at one radius
 * is not split).
 *
 * Stage 2-3 (quilting energy + layer pick) and Stage 4 (relocation) are added
 * incrementally; until then the pick is a radial-prior stub.
 * ==========================================================================*/

#include "overlap_quilt.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../common/pca.h"
#include "../common/raw_sample.h"
#include "../split/multicut_wrap.h"

/* ------------------------------------------------------------------ util */

static void *xmalloc(size_t n) {
    void *p = malloc(n > 0 ? n : 1);
    if (!p) { fprintf(stderr, "overlap_quilt: OOM %zu\n", n); exit(1); }
    return p;
}
static void *xcalloc(size_t c, size_t s) {
    void *p = calloc(c > 0 ? c : 1, s);
    if (!p) { fprintf(stderr, "overlap_quilt: OOM %zux%zu\n", c, s); exit(1); }
    return p;
}

/* simple union-find over ints (malloc, path-halving) */
typedef struct { int32_t *p; } UF;
static void uf_init(UF *uf, size_t n) {
    uf->p = (int32_t *)xmalloc(n * sizeof(int32_t));
    for (size_t i = 0; i < n; i++) uf->p[i] = (int32_t)i;
}
static int32_t uf_find(UF *uf, int32_t x) {
    while (uf->p[x] != x) { uf->p[x] = uf->p[uf->p[x]]; x = uf->p[x]; }
    return x;
}
static void uf_union(UF *uf, int32_t a, int32_t b) {
    int32_t ra = uf_find(uf, a), rb = uf_find(uf, b);
    if (ra != rb) uf->p[rb] = ra;
}
static void uf_free(UF *uf) { free(uf->p); uf->p = NULL; }

void QuiltOpts_default(QuiltOpts *o) {
    memset(o, 0, sizeof *o);
    o->axis_point[0] = 0.0f; o->axis_point[1] = 3405.0f; o->axis_point[2] = 2878.0f;
    o->axis_dir[0] = 1.0f;   o->axis_dir[1] = 0.0f;      o->axis_dir[2] = 0.0f;
    o->raw_dir = NULL;
    o->chunk = 128;
    o->radius_gate = 0.0;    /* => pitch/2 */
    o->cell = 2.0;
    o->tex_du = 1.0; o->tex_dv = 1.0;
    o->normal_range = 2.0; o->normal_samples = 5;
    o->spiral_a = 0.0; o->spiral_b = 0.0; o->pitch = 9.5;
    o->relocate = 1;
    o->cross_component_only = 0;
    o->verbose = 0;
}

/* ------------------------------------------------------------- geometry */

/* axis-perpendicular radius of a point p (z,y,x). e1,e2 span the plane. */
double Quilt_point_radius(const float *p, const float ap[3],
                          const double e1[3], const double e2[3]) {
    double d[3] = { (double)p[0]-ap[0], (double)p[1]-ap[1], (double)p[2]-ap[2] };
    double c1 = d[0]*e1[0] + d[1]*e1[1] + d[2]*e1[2];
    double c2 = d[0]*e2[0] + d[1]*e2[1] + d[2]*e2[2];
    return sqrt(c1*c1 + c2*c2);
}

/* ------------------------------------------------- adjacency (edge->faces) */

/* Build face adjacency as an edge-keyed pass: sort (v_lo,v_hi,face) and pair
 * faces sharing an edge. Returns malloc'd pair arrays; caller frees.
 * Also returns, per pair, the shared 2D(UV) edge length for the weight. */
typedef QuiltAdjEdge AdjEdge;   /* exported in the header for SeamOwn */

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x < y) ? -1 : (x > y);
}

/* key stream: (edgekey<<? ) -- we sort a struct instead */
typedef struct { uint64_t ekey; int32_t face; int32_t va, vb; } EdgeRec;
static int cmp_edgerec(const void *a, const void *b) {
    const EdgeRec *x = (const EdgeRec *)a, *y = (const EdgeRec *)b;
    if (x->ekey != y->ekey) return (x->ekey < y->ekey) ? -1 : 1;
    return (x->face < y->face) ? -1 : (x->face > y->face);
}

QuiltAdjEdge *Quilt_build_adjacency(const int32_t *faces, size_t nf,
                                    const float *uv, size_t *out_n) {
    size_t ne = nf * 3, i = 0;
    EdgeRec *er = (EdgeRec *)xmalloc(ne * sizeof(EdgeRec));
    size_t k = 0;
    for (i = 0; i < nf; i++) {
        int32_t v[3] = { faces[i*3+0], faces[i*3+1], faces[i*3+2] };
        for (int e = 0; e < 3; e++) {
            int32_t a = v[e], b = v[(e+1)%3];
            int32_t lo = a < b ? a : b, hi = a < b ? b : a;
            er[k].ekey = ((uint64_t)(uint32_t)lo << 32) | (uint32_t)hi;
            er[k].face = (int32_t)i; er[k].va = lo; er[k].vb = hi;
            k++;
        }
    }
    qsort(er, ne, sizeof(EdgeRec), cmp_edgerec);
    /* runs of equal ekey -> pair consecutive faces (manifold: run len 2) */
    size_t cap = nf + 16, cnt = 0;
    AdjEdge *adj = (AdjEdge *)xmalloc(cap * sizeof(AdjEdge));
    i = 0;
    while (i < ne) {
        size_t j = i + 1;
        while (j < ne && er[j].ekey == er[i].ekey) j++;
        /* pair every two faces in the run (usually exactly 2) */
        for (size_t a = i; a + 1 < j; a += 2) {
            double du = (double)uv[(size_t)er[a].va*2+0] - uv[(size_t)er[a].vb*2+0];
            double dv = (double)uv[(size_t)er[a].va*2+1] - uv[(size_t)er[a].vb*2+1];
            if (cnt >= cap) { cap *= 2; adj = (AdjEdge *)realloc(adj, cap*sizeof(AdjEdge)); if(!adj){fprintf(stderr,"OOM\n");exit(1);} }
            adj[cnt].fa = er[a].face; adj[cnt].fb = er[a+1].face;
            adj[cnt].uvlen = (float)sqrt(du*du + dv*dv);
            cnt++;
        }
        i = j;
    }
    free(er);
    *out_n = cnt;
    return adj;
}

/* ------------------------------------------------- overlap detection */

/* Bin faces into UV cells by their UV bbox; a cell holds candidate faces.
 * Overlap pair = two faces sharing a cell with |r_fa - r_fb| > radius_gate and
 * not mesh-adjacent. We union both into regions and record the pair. */
typedef struct { int32_t fa, fb; } OvPair;

/* ---------------------------------------------- per-region multicut split */

/* Split one region's faces into connected layers via the lifted multicut.
 * rfaces[nrf] = global face ids; g2l maps global->local (caller sets/clears).
 * adj/nadj = GLOBAL adjacency; ov/nov = GLOBAL overlap pairs. Writes a 0-based
 * layer id per region face into out_layer (indexed by local id). Returns the
 * number of layers, or -1 on failure (caller falls back to 1 layer). */
static int split_region_multicut(const int32_t *rfaces, size_t nrf,
                                  const int32_t *g2l,
                                  const AdjEdge *adj, size_t nadj,
                                  const OvPair *ov, size_t nov,
                                  const int32_t *region_of, int32_t rid,
                                  int32_t *out_layer) {
    if (nrf == 0) return -1;
    if (nrf == 1) { out_layer[0] = 0; return 1; }

    /* local adjacency edges (both endpoints in this region) */
    size_t cap_a = nrf + 8, na = 0;
    int32_t *af = (int32_t *)xmalloc(cap_a * sizeof(int32_t));
    int32_t *at = (int32_t *)xmalloc(cap_a * sizeof(int32_t));
    double  *aw = (double  *)xmalloc(cap_a * sizeof(double));
    double avg = 0.0; size_t navg = 0;
    for (size_t i = 0; i < nadj; i++) {
        int32_t fa = adj[i].fa, fb = adj[i].fb;
        if (region_of[fa] != rid || region_of[fb] != rid) continue;
        if (na >= cap_a) { cap_a *= 2;
            af = (int32_t*)realloc(af, cap_a*sizeof(int32_t));
            at = (int32_t*)realloc(at, cap_a*sizeof(int32_t));
            aw = (double *)realloc(aw, cap_a*sizeof(double));
            if(!af||!at||!aw){fprintf(stderr,"OOM\n");exit(1);} }
        af[na] = g2l[fa]; at[na] = g2l[fb]; aw[na] = adj[i].uvlen;
        avg += adj[i].uvlen; navg++; na++;
    }
    if (navg > 0) avg /= (double)navg; else avg = 1.0;
    if (avg < 1e-9) avg = 1.0;

    /* overlap edges (both endpoints in this region) */
    size_t cap_o = 8, no = 0;
    int32_t *of = (int32_t *)xmalloc(cap_o * sizeof(int32_t));
    int32_t *ot = (int32_t *)xmalloc(cap_o * sizeof(int32_t));
    for (size_t i = 0; i < nov; i++) {
        int32_t fa = ov[i].fa, fb = ov[i].fb;
        if (region_of[fa] != rid || region_of[fb] != rid) continue;
        if (no >= cap_o) { cap_o *= 2;
            of = (int32_t*)realloc(of, cap_o*sizeof(int32_t));
            ot = (int32_t*)realloc(ot, cap_o*sizeof(int32_t));
            if(!of||!ot){fprintf(stderr,"OOM\n");exit(1);} }
        of[no] = g2l[fa]; ot[no] = g2l[fb]; no++;
    }

    /* lifted arrays = adjacency (positive, weight uvlen/avg) THEN overlap
     * (negative, -BIG). adjacency passed twice: as the original graph and as
     * the first na lifted edges (OverlapSep convention). */
    size_t nl = na + no;
    int32_t *lf = (int32_t *)xmalloc(nl * sizeof(int32_t));
    int32_t *lt = (int32_t *)xmalloc(nl * sizeof(int32_t));
    double  *lw = (double  *)xmalloc(nl * sizeof(double));
    for (size_t i = 0; i < na; i++) { lf[i]=af[i]; lt[i]=at[i]; lw[i]=aw[i]/avg; }
    for (size_t i = 0; i < no; i++) { lf[na+i]=of[i]; lt[na+i]=ot[i]; lw[na+i]=-1.0e6; }

    int32_t nlayers = 0;
    int rc = LiftedMulticut_kernighan_lin(
        (int32_t)nrf, (int32_t)na, af, at,
        (int32_t)nl, lf, lt, lw,
        /*min_clusters*/ 2, out_layer, &nlayers);

    free(af); free(at); free(aw); free(of); free(ot); free(lf); free(lt); free(lw);
    if (rc != 0) return -1;
    return (int)nlayers;
}

/* ------------------------------------------------ texture / quilting energy */

/* Rasterize faces flist[nfl] (global ids) into a local grid, sampling raw CT
 * (normal-max) at the barycentric surface point per covered pixel. img is set
 * to the sample (>=0) where covered, cov[]=1 there. NAN/0 elsewhere. */
void Quilt_raster_faces(const int32_t *flist, size_t nfl,
                        const float *verts, const float *uv, const float *nrm,
                        const int32_t *faces,
                        struct CubeTable *ct, double range, int nsteps,
                        double u0, double v0, double du, double dv,
                        int gw, int gh, float *img, uint8_t *cov) {
    for (int i = 0; i < gw*gh; i++) { img[i] = 0.0f; cov[i] = 0; }
    for (size_t fi = 0; fi < nfl; fi++) {
        int32_t f = flist[fi];
        size_t a=(size_t)faces[f*3+0], b=(size_t)faces[f*3+1], c=(size_t)faces[f*3+2];
        double ua=((double)uv[a*2+0]-u0)/du, va=((double)uv[a*2+1]-v0)/dv;
        double ub=((double)uv[b*2+0]-u0)/du, vb=((double)uv[b*2+1]-v0)/dv;
        double uc=((double)uv[c*2+0]-u0)/du, vc=((double)uv[c*2+1]-v0)/dv;
        double A2=(ub-ua)*(vc-va)-(vb-va)*(uc-ua);
        if (fabs(A2)<1e-12) continue;
        double lox=ua<ub?(ua<uc?ua:uc):(ub<uc?ub:uc), hix=ua>ub?(ua>uc?ua:uc):(ub>uc?ub:uc);
        double loy=va<vb?(va<vc?va:vc):(vb<vc?vb:vc), hiy=va>vb?(va>vc?va:vc):(vb>vc?vb:vc);
        int x0=(int)ceil(lox-0.001),x1=(int)floor(hix+0.001),y0=(int)ceil(loy-0.001),y1=(int)floor(hiy+0.001);
        if(x0<0)x0=0; if(y0<0)y0=0; if(x1>=gw)x1=gw-1; if(y1>=gh)y1=gh-1;
        for (int yy=y0;yy<=y1;yy++) for (int xx=x0;xx<=x1;xx++){
            double pu=xx,pv=yy;
            double l1=((pu-ua)*(vc-va)-(pv-va)*(uc-ua))/A2;
            double l2=((ub-ua)*(pv-va)-(vb-va)*(pu-ua))/A2;
            double l0=1.0-l1-l2;
            if(l0<-1e-6||l1<-1e-6||l2<-1e-6) continue;
            float p3[3],n3[3]; double nn=0,s;
            for(int k=0;k<3;k++){
                p3[k]=(float)(l0*verts[a*3+k]+l1*verts[b*3+k]+l2*verts[c*3+k]);
                n3[k]=nrm?(float)(l0*nrm[a*3+k]+l1*nrm[b*3+k]+l2*nrm[c*3+k]):0.0f;
            }
            nn=sqrt((double)n3[0]*n3[0]+(double)n3[1]*n3[1]+(double)n3[2]*n3[2]);
            if(nn>1e-6){n3[0]/=(float)nn;n3[1]/=(float)nn;n3[2]/=(float)nn;}
            s=sample_vertex(ct,p3,nn>0.5?n3:NULL,range,nsteps);
            if(s<0) continue;
            size_t pi=(size_t)yy*gw+xx;
            img[pi]=(float)s; cov[pi]=1;   /* last writer wins; layers are single-cover here */
        }
    }
}

/* Quilting boundary energy of a layer against the anchor (Efros min-error-
 * boundary / Zomet false-edge): mean |vL - vA| in INTENSITY across the seam --
 * layer pixels 4-adjacent to an anchor pixel. The layer and anchor are disjoint
 * in UV (adjacent, not overlapping), so this is the value discontinuity the
 * seam would introduce; papyrus is near-uniform in mean brightness, so this is
 * dominated by whether the fiber pattern (dark fibers / bright gaps) LINES UP
 * across the boundary. To stay robust to a per-layer brightness offset, each
 * side is compared after subtracting its own seam-band mean. Returns a large
 * value when there is no real shared boundary. */
double Quilt_boundary_energy(const float *vL, const uint8_t *covL,
                             const float *vA, const uint8_t *covA,
                             int gw, int gh) {
    int dxs[4]={1,-1,0,0}, dys[4]={0,0,1,-1};
    /* pass 1: seam-band means for offset removal */
    double mL=0, mA=0; size_t nb=0;
    for (int y=1;y<gh-1;y++) for (int x=1;x<gw-1;x++){
        size_t pi=(size_t)y*gw+x; if(!covL[pi]) continue;
        for(int d=0;d<4;d++){ size_t qi=(size_t)(y+dys[d])*gw+(x+dxs[d]);
            if(covA[qi]){ mL+=vL[pi]; mA+=vA[qi]; nb++; break; } }
    }
    if (nb<4) return 1e18;
    mL/=(double)nb; mA/=(double)nb;
    double off = mL - mA;   /* layer-vs-anchor brightness offset */
    /* pass 2: mean |(vL - off) - vA| over seam pairs */
    double sum=0.0; size_t n=0;
    for (int y=1;y<gh-1;y++) for (int x=1;x<gw-1;x++){
        size_t pi=(size_t)y*gw+x; if(!covL[pi]) continue;
        for(int d=0;d<4;d++){ int nx=x+dxs[d],ny=y+dys[d]; size_t qi=(size_t)ny*gw+nx;
            if(covA[qi]){ sum+=fabs(((double)vL[pi]-off)-(double)vA[qi]); n++; } }
    }
    return sum/(double)n;
}

/* ---------------------------------------------------------------- run */

int Quilt_run(Arena_T arena,
              const float *verts, size_t nv,
              const int32_t *faces, size_t nf,
              const float *uv,
              const QuiltOpts *opts, QuiltResult *out) {
    memset(out, 0, sizeof *out);
    if (nv < 3 || nf < 1 || uv == NULL) return -1;

    QuiltOpts o = *opts;
    double radius_gate = o.radius_gate > 0.0 ? o.radius_gate
                       : (o.pitch > 0.0 ? o.pitch * 0.5 : 4.75);

    /* axis frame + per-face mean radius */
    float e1f[3], e2f[3];
    PCA_orthonormal_basis(o.axis_dir, e1f, e2f);
    double e1[3] = { e1f[0], e1f[1], e1f[2] }, e2[3] = { e2f[0], e2f[1], e2f[2] };
    double *rf = (double *)xmalloc(nf * sizeof(double));
    for (size_t f = 0; f < nf; f++) {
        double r = 0.0;
        for (int k = 0; k < 3; k++)
            r += Quilt_point_radius(&verts[(size_t)faces[f*3+k]*3], o.axis_point, e1, e2);
        rf[f] = r / 3.0;
    }

    /* per-face UV centroid + bbox for cell binning */
    double *fu = (double *)xmalloc(nf * sizeof(double));
    double *fv = (double *)xmalloc(nf * sizeof(double));
    double umin = 1e30, vmin = 1e30, umax = -1e30, vmax = -1e30;
    for (size_t f = 0; f < nf; f++) {
        double su = 0, sv = 0;
        for (int k = 0; k < 3; k++) {
            su += uv[(size_t)faces[f*3+k]*2+0];
            sv += uv[(size_t)faces[f*3+k]*2+1];
        }
        fu[f] = su/3.0; fv[f] = sv/3.0;
        if (fu[f]<umin)umin=fu[f]; if(fu[f]>umax)umax=fu[f];
        if (fv[f]<vmin)vmin=fv[f]; if(fv[f]>vmax)vmax=fv[f];
    }

    /* UV cell grid (centroid binning) */
    double cell = o.cell > 0 ? o.cell : 2.0;
    int gu = (int)((umax-umin)/cell) + 1, gv = (int)((vmax-vmin)/cell) + 1;
    if (gu < 1) gu = 1; if (gv < 1) gv = 1;
    size_t ncell = (size_t)gu * (size_t)gv;
    /* CSR cell->faces */
    int32_t *coff = (int32_t *)xcalloc(ncell + 1, sizeof(int32_t));
    int32_t *fcell = (int32_t *)xmalloc(nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) {
        int cu = (int)((fu[f]-umin)/cell), cv = (int)((fv[f]-vmin)/cell);
        if (cu<0)cu=0; if(cu>=gu)cu=gu-1; if(cv<0)cv=0; if(cv>=gv)cv=gv-1;
        fcell[f] = cv*gu + cu; coff[fcell[f]+1]++;
    }
    for (size_t i = 0; i < ncell; i++) coff[i+1] += coff[i];
    int32_t *ccur = (int32_t *)xmalloc(ncell * sizeof(int32_t));
    memcpy(ccur, coff, ncell * sizeof(int32_t));
    int32_t *clist = (int32_t *)xmalloc(nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) clist[ccur[fcell[f]]++] = (int32_t)f;
    free(ccur);

    /* adjacency (for region-building + multicut) */
    size_t nadj = 0;
    AdjEdge *adj = Quilt_build_adjacency(faces, nf, uv, &nadj);
    /* adjacency lookup as a set for the "not adjacent" test in overlap detect:
     * hash the (min,max) face pairs */
    uint64_t *akeys = (uint64_t *)xmalloc((nadj + 1) * sizeof(uint64_t));
    for (size_t i = 0; i < nadj; i++) {
        int32_t a = adj[i].fa, b = adj[i].fb, lo=a<b?a:b, hi=a<b?b:a;
        akeys[i] = ((uint64_t)(uint32_t)lo << 32) | (uint32_t)hi;
    }
    qsort(akeys, nadj, sizeof(uint64_t), cmp_u64);

    /* A local delamination diagnostic starts with already-disconnected input
     * sheets. In that mode, compare only cross-component faces so a low radial
     * gate cannot mistake ordinary curvature within either sheet for layers. */
    int32_t *face_comp = NULL;
    if (o.cross_component_only) {
        UF comp_uf;
        uf_init(&comp_uf, nf);
        for (size_t i = 0; i < nadj; i++) {
            uf_union(&comp_uf, adj[i].fa, adj[i].fb);
        }
        face_comp = (int32_t *)xmalloc(nf * sizeof *face_comp);
        for (size_t f = 0; f < nf; f++) {
            face_comp[f] = uf_find(&comp_uf, (int32_t)f);
        }
        uf_free(&comp_uf);
    }

    /* overlap pairs: same cell, |dr|>gate, not adjacent */
    size_t cap_ov = 64, nov = 0;
    OvPair *ov = (OvPair *)xmalloc(cap_ov * sizeof(OvPair));
    uint8_t *is_ov = (uint8_t *)xcalloc(nf, 1);
    size_t multi_before = 0;
    for (size_t c = 0; c < ncell; c++) {
        int32_t a0 = coff[c], a1 = coff[c+1];
        int has_multi = 0;
        for (int32_t i = a0; i < a1; i++) {
            for (int32_t j = i+1; j < a1; j++) {
                int32_t fa = clist[i], fb = clist[j];
                if (face_comp != NULL && face_comp[fa] == face_comp[fb]) continue;
                if (fabs(rf[fa]-rf[fb]) <= radius_gate) continue;
                /* adjacency test */
                int32_t lo=fa<fb?fa:fb, hi=fa<fb?fb:fa;
                uint64_t key = ((uint64_t)(uint32_t)lo<<32)|(uint32_t)hi;
                /* binary search akeys */
                size_t los=0, his=nadj, adjq=0;
                while (los<his){ size_t m=(los+his)/2; if(akeys[m]<key)los=m+1; else his=m; }
                if (los<nadj && akeys[los]==key) adjq=1;
                if (adjq) continue;
                if (nov>=cap_ov){cap_ov*=2; ov=(OvPair*)realloc(ov,cap_ov*sizeof(OvPair)); if(!ov){fprintf(stderr,"OOM\n");exit(1);}}
                ov[nov].fa=fa; ov[nov].fb=fb; nov++;
                is_ov[fa]=1; is_ov[fb]=1; has_multi=1;
            }
        }
        if (has_multi) multi_before++;
    }

    /* regions: union overlap pairs + adjacency edges between two overlap faces */
    UF uf; uf_init(&uf, nf);
    for (size_t i = 0; i < nov; i++) uf_union(&uf, ov[i].fa, ov[i].fb);
    for (size_t i = 0; i < nadj; i++)
        if (is_ov[adj[i].fa] && is_ov[adj[i].fb]) uf_union(&uf, adj[i].fa, adj[i].fb);

    /* dense region ids over faces that carry an overlap */
    int32_t *region_of = (int32_t *)xmalloc(nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) region_of[f] = -1;
    int32_t nreg = 0;
    int32_t *root2rid = (int32_t *)xmalloc(nf * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) root2rid[f] = -1;
    for (size_t f = 0; f < nf; f++) {
        if (!is_ov[f]) continue;
        int32_t r = uf_find(&uf, (int32_t)f);
        if (root2rid[r] < 0) root2rid[r] = nreg++;
        region_of[f] = root2rid[r];
    }
    /* region = overlap faces only; the surrounding clean sheet (the anchor) is
     * detected separately per region from adjacency below. */

    /* per-region face lists (CSR) */
    int32_t *roff = (int32_t *)xcalloc((size_t)nreg + 1, sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) if (region_of[f]>=0) roff[region_of[f]+1]++;
    for (int32_t i = 0; i < nreg; i++) roff[i+1] += roff[i];
    int32_t *rcur = (int32_t *)xmalloc(((size_t)nreg+1) * sizeof(int32_t));
    memcpy(rcur, roff, ((size_t)nreg+1)*sizeof(int32_t));
    int32_t *rfaces = (int32_t *)xmalloc((roff[nreg] > 0 ? (size_t)roff[nreg] : 1) * sizeof(int32_t));
    for (size_t f = 0; f < nf; f++) if (region_of[f]>=0) rfaces[rcur[region_of[f]]++] = (int32_t)f;
    free(rcur);

    /* outputs */
    out->uv = (float *)ARENA_ALLOC(arena, (long)(nv*2*sizeof(float)));
    memcpy(out->uv, uv, nv*2*sizeof(float));
    out->face_keep = (uint8_t *)ARENA_ALLOC(arena, (long)nf);
    memset(out->face_keep, 1, nf);
    out->face_dec = (uint8_t *)ARENA_CALLOC(arena, (long)nf, 1);
    out->face_region = (int32_t *)ARENA_ALLOC(arena, (long)(nf*sizeof(int32_t)));
    out->face_layer = (int32_t *)ARENA_ALLOC(arena, (long)(nf*sizeof(int32_t)));
    out->face_energy = (double *)ARENA_ALLOC(arena, (long)(nf*sizeof(double)));
    for (size_t f = 0; f < nf; f++) {
        out->face_region[f] = region_of[f];
        out->face_layer[f] = -1;
        out->face_energy[f] = -1.0;
    }

    /* per-region ANCHOR faces (non-overlap faces adjacent to the region) --
     * the surrounding clean sheet the quilting energy matches against. */
    int32_t *anoff = (int32_t *)xcalloc((size_t)nreg + 1, sizeof(int32_t));
    for (size_t i = 0; i < nadj; i++) {
        int32_t fa = adj[i].fa, fb = adj[i].fb;
        if (region_of[fa]>=0 && region_of[fb]<0) anoff[region_of[fa]+1]++;
        else if (region_of[fb]>=0 && region_of[fa]<0) anoff[region_of[fb]+1]++;
    }
    for (int32_t i = 0; i < nreg; i++) anoff[i+1] += anoff[i];
    int32_t *ancur = (int32_t *)xmalloc(((size_t)nreg+1)*sizeof(int32_t));
    memcpy(ancur, anoff, ((size_t)nreg+1)*sizeof(int32_t));
    int32_t *afaces = (int32_t *)xmalloc((anoff[nreg]>0?(size_t)anoff[nreg]:1)*sizeof(int32_t));
    for (size_t i = 0; i < nadj; i++) {
        int32_t fa = adj[i].fa, fb = adj[i].fb;
        if (region_of[fa]>=0 && region_of[fb]<0) afaces[ancur[region_of[fa]]++] = fb;
        else if (region_of[fb]>=0 && region_of[fa]<0) afaces[ancur[region_of[fb]]++] = fa;
    }
    free(ancur);

    /* texture: cube table + per-vertex normals (only if RAW given) */
    CubeTable ct; int have_ct = 0; float *nrm = NULL;
    if (o.raw_dir != NULL &&
        cubetable_init(&ct, arena, o.raw_dir, o.chunk, verts, nv, o.normal_range+2.0)==0) {
        nrm = vertex_normals(verts, nv, faces, nf); have_ct = 1;
    }
    out->have_texture = have_ct;

    /* per-region multicut split + quilting-energy layer pick */
    int32_t *g2l = (int32_t *)xmalloc(nf * sizeof(int32_t));
    int32_t *comp2layer = NULL;
    for (size_t f = 0; f < nf; f++) g2l[f] = -1;
    if (face_comp != NULL) {
        comp2layer = (int32_t *)xmalloc(nf * sizeof *comp2layer);
        for (size_t f = 0; f < nf; f++) comp2layer[f] = -1;
    }
    int32_t *loclayer = (int32_t *)xmalloc((roff[nreg]>0?(size_t)roff[nreg]:1) * sizeof(int32_t));
    size_t n_layers_total = 0, max_region_faces = 0;
    const size_t REGION_CAP = 40000;   /* skip the fused core (no anchor anyway) */
    const int    GRID_CAP   = 3000;
    double du = o.tex_du>0?o.tex_du:1.0, dv = o.tex_dv>0?o.tex_dv:1.0;
    double e_sum = 0.0; size_t e_n = 0;
    size_t n_resolved = 0, n_unresolved = 0;
    size_t why_onelayer=0, why_manylayer=0, why_noanchor=0, why_toobig=0, why_grid=0, why_noseam=0;

    for (int32_t r = 0; r < nreg; r++) {
        int32_t b0 = roff[r], b1 = roff[r+1];
        size_t nrf = (size_t)(b1 - b0);
        if (nrf > max_region_faces) max_region_faces = nrf;
        for (int32_t i = b0; i < b1; i++) g2l[rfaces[i]] = i - b0;

        int nlayers = 0;
        if (face_comp != NULL) {
            for (size_t i = 0; i < nrf; i++) {
                int32_t root = face_comp[rfaces[b0+(int32_t)i]];
                if (comp2layer[root] < 0) comp2layer[root] = nlayers++;
                loclayer[i] = comp2layer[root];
            }
            for (size_t i = 0; i < nrf; i++) {
                int32_t root = face_comp[rfaces[b0+(int32_t)i]];
                comp2layer[root] = -1;
            }
        } else {
            nlayers = split_region_multicut(&rfaces[b0], nrf, g2l, adj, nadj,
                                             ov, nov, region_of, r, loclayer);
        }
        if (nlayers < 1) { nlayers = 1; for (size_t i=0;i<nrf;i++) loclayer[i]=0; }
        n_layers_total += (size_t)nlayers;
        for (size_t i = 0; i < nrf; i++)
            out->face_layer[rfaces[b0+(int32_t)i]] = loclayer[i];

        int32_t n_anchor = anoff[r+1]-anoff[r];
        int resolvable = have_ct && nlayers >= 2 && nrf <= REGION_CAP && n_anchor > 0;
        if (!resolvable) {
            if (nlayers < 2) why_onelayer++;
            else if (nrf > REGION_CAP) why_toobig++;
            else if (n_anchor == 0) why_noanchor++;
            n_unresolved++;
            for (int32_t i = b0; i < b1; i++) g2l[rfaces[i]] = -1;
            continue;
        }

        /* local grid over region + anchor UV bbox (+2px margin) */
        double gu0=1e30,gv0=1e30,gu1=-1e30,gv1=-1e30;
        for (int32_t i = b0; i < b1; i++) { int32_t f=rfaces[i];
            for(int k=0;k<3;k++){ double u=uv[(size_t)faces[f*3+k]*2+0],v=uv[(size_t)faces[f*3+k]*2+1];
                if(u<gu0)gu0=u; if(u>gu1)gu1=u; if(v<gv0)gv0=v; if(v>gv1)gv1=v; } }
        for (int32_t i = anoff[r]; i < anoff[r+1]; i++) { int32_t f=afaces[i];
            for(int k=0;k<3;k++){ double u=uv[(size_t)faces[f*3+k]*2+0],v=uv[(size_t)faces[f*3+k]*2+1];
                if(u<gu0)gu0=u; if(u>gu1)gu1=u; if(v<gv0)gv0=v; if(v>gv1)gv1=v; } }
        gu0-=2*du; gv0-=2*dv; gu1+=2*du; gv1+=2*dv;
        int gw=(int)((gu1-gu0)/du)+1, gh=(int)((gv1-gv0)/dv)+1;
        if (gw<3||gh<3||gw>GRID_CAP||gh>GRID_CAP) { n_unresolved++; why_grid++;
            for (int32_t i=b0;i<b1;i++) g2l[rfaces[i]]=-1; continue; }

        size_t np=(size_t)gw*gh;
        float *imgA=(float*)xmalloc(np*sizeof(float));
        uint8_t *covA=(uint8_t*)xmalloc(np),*covL=(uint8_t*)xmalloc(np);
        float *imgL=(float*)xmalloc(np*sizeof(float));
        Quilt_raster_faces(&afaces[anoff[r]], (size_t)n_anchor, verts, uv, nrm, faces,
                     &ct, o.normal_range, o.normal_samples, gu0,gv0,du,dv,gw,gh,imgA,covA);

        /* keep ONE layer, drop the rest (de-blur). Prefer the SUBSTANTIAL
         * layer (>=15% of region) with the lowest quilting seam energy against
         * the anchor; if none of the substantial layers has a usable seam, keep
         * the LARGEST layer so we never leave a hole. Slivers can win an
         * energy contest by luck, so they are excluded from the pick. */
        int32_t *lsz=(int32_t*)xcalloc((size_t)nlayers,sizeof(int32_t));
        double *layer_energy=(double*)xmalloc((size_t)nlayers*sizeof(double));
        for (size_t i=0;i<nrf;i++) lsz[loclayer[i]]++;
        for (int L=0;L<nlayers;L++) layer_energy[L]=1e18;
        int largest=0; for(int L=1;L<nlayers;L++) if(lsz[L]>lsz[largest]) largest=L;
        size_t min_layer=(size_t)(0.15*(double)nrf); if(min_layer<4)min_layer=4;
        int keep=-1; double bestE=1e19;
        int32_t *lf=(int32_t*)xmalloc(nrf*sizeof(int32_t));
        for (int L=0; L<nlayers; L++) {
            size_t nl=0;
            for (size_t i=0;i<nrf;i++) if(loclayer[i]==L) lf[nl++]=rfaces[b0+(int32_t)i];
            Quilt_raster_faces(lf, nl, verts, uv, nrm, faces, &ct, o.normal_range,
                         o.normal_samples, gu0,gv0,du,dv,gw,gh,imgL,covL);
            double E = Quilt_boundary_energy(imgL, covL, imgA, covA, gw, gh);
            layer_energy[L] = E;
            if ((size_t)lsz[L] >= min_layer && E < bestE) {
                bestE=E;
                keep=L;
            }
        }
        free(lf); free(imgA);free(covA);free(covL);free(imgL);
        if (keep < 0 || bestE >= 1e17) { keep = largest; why_noseam++; }  /* fallback: largest */
        else { e_sum += bestE; e_n++; }

        n_resolved++;
        for (size_t i = 0; i < nrf; i++) {
            int32_t gf = rfaces[b0 + (int32_t)i];
            out->face_energy[gf] = layer_energy[loclayer[i]];
            if (loclayer[i] == keep) { out->face_dec[gf]=1; out->face_keep[gf]=1; out->n_kept++; }
            else { out->face_dec[gf]=3; out->face_keep[gf]=0; out->n_dropped++; }  /* Stage 4 relocates */
        }
        free(lsz);
        free(layer_energy);
        for (int32_t i = b0; i < b1; i++) g2l[rfaces[i]] = -1;
    }

    /* ---- Stage 4: relocate fully-dropped SEPARATE components ----
     * A dropped layer that is a separate connected component (never adjacent to
     * a kept face) is a distinct piece we can move rather than lose -- shift its
     * u into an empty band past the sheet (v=z stays), so it survives as a
     * compact chart. Dropped plies still attached to kept geometry stay dropped
     * (they are the occluded second layer of a delamination). */
    if (o.relocate) {
        UF df; uf_init(&df, nf);
        for (size_t i=0;i<nadj;i++){ int32_t a=adj[i].fa,b=adj[i].fb;
            if (out->face_dec[a]==3 && out->face_dec[b]==3) uf_union(&df,a,b); }
        /* which dropped-CC roots touch a kept face -> not relocatable */
        uint8_t *touch_keep=(uint8_t*)xcalloc(nf,1);
        for (size_t i=0;i<nadj;i++){ int32_t a=adj[i].fa,b=adj[i].fb;
            if (out->face_dec[a]==3 && out->face_keep[b]) touch_keep[uf_find(&df,a)]=1;
            if (out->face_dec[b]==3 && out->face_keep[a]) touch_keep[uf_find(&df,b)]=1; }
        /* umax of the whole sheet -> pack relocated pieces to its right */
        double umax_all=-1e30;
        for (size_t i=0;i<nv;i++) if(out->uv[i*2+0]>umax_all) umax_all=out->uv[i*2+0];
        double place = umax_all + 20.0;
        /* per dropped-CC root: gather its faces, if not touch_keep and >= min,
         * shift its verts' u to `place` (packed left-aligned). */
        int32_t *rroot=(int32_t*)xmalloc(nf*sizeof(int32_t));
        for (size_t f=0;f<nf;f++) rroot[f]=-1;
        for (size_t f=0;f<nf;f++) if(out->face_dec[f]==3){ int32_t r=uf_find(&df,(int32_t)f); if(!touch_keep[r]) rroot[f]=r; }
        uint8_t *done=(uint8_t*)xcalloc(nf,1);   /* process each root once */
        for (size_t f=0;f<nf;f++){
            int32_t r=rroot[f]; if(r<0||done[r]) continue; done[r]=1;
            /* first pass: count faces + u,v bbox of this component */
            double cu0=1e30,cu1=-1e30; size_t cn=0;
            for (size_t g=0; g<nf; g++) if(rroot[g]==r){ cn++;
                for(int k=0;k<3;k++){ double u=out->uv[(size_t)faces[g*3+k]*2+0]; if(u<cu0)cu0=u; if(u>cu1)cu1=u; } }
            if (cn < 20) continue;   /* skip debris */
            double shift = place - cu0;
            /* shift verts (component-exclusive since not touching kept) */
            for (size_t g=0; g<nf; g++) if(rroot[g]==r)
                for(int k=0;k<3;k++){ size_t vi=(size_t)faces[g*3+k];
                    out->uv[vi*2+0] += (float)shift; }
            for (size_t g=0; g<nf; g++) if(rroot[g]==r){ out->face_dec[g]=2; out->face_keep[g]=1; out->n_relocated++; out->n_dropped--; }
            place += (cu1-cu0) + 15.0;   /* pack next piece to the right */
        }
        uf_free(&df); free(touch_keep); free(rroot); free(done);
    }

    /* multi-cover cells AFTER repair (kept faces only) */
    size_t multi_after = 0;
    for (size_t c = 0; c < ncell; c++) {
        int32_t a0=coff[c],a1=coff[c+1]; int hm=0;
        for (int32_t i=a0;i<a1&&!hm;i++){ if(!out->face_keep[clist[i]])continue;
            for (int32_t j=i+1;j<a1;j++){ if(!out->face_keep[clist[j]])continue;
                if (face_comp != NULL &&
                    face_comp[clist[i]] == face_comp[clist[j]]) continue;
                if (fabs(rf[clist[i]]-rf[clist[j]])>radius_gate){hm=1;break;} } }
        if (hm) multi_after++;
    }

    out->n_overlap_faces = 0;
    for (size_t f=0; f<nf; f++) if (is_ov[f]) out->n_overlap_faces++;
    out->n_regions = (size_t)nreg;
    out->n_layers_total = n_layers_total;
    out->max_region_faces = max_region_faces;
    out->multi_cells_before = multi_before;
    out->multi_cells_after = multi_after;
    out->energy_mean = e_n ? e_sum/(double)e_n : 0.0;

    if (o.verbose) {
        fprintf(stderr, "[quilt] regions=%d (resolved=%zu unresolved=%zu) overlap_faces=%zu "
                "layers=%zu max_region=%zu kept=%zu dropped=%zu multi_cells %zu->%zu "
                "E=%.2f tex=%d gate=%.2f\n",
                nreg, n_resolved, n_unresolved, out->n_overlap_faces, n_layers_total,
                max_region_faces, out->n_kept, out->n_dropped, multi_before, multi_after,
                out->energy_mean, have_ct, radius_gate);
        fprintf(stderr, "[quilt] unresolved reasons: one_layer=%zu many_layer=%zu "
                "no_anchor=%zu too_big=%zu grid=%zu no_seam=%zu\n",
                why_onelayer, why_manylayer, why_noanchor, why_toobig, why_grid, why_noseam);
    }

    free(nrm);
    free(rf); free(fu); free(fv); free(coff); free(fcell); free(clist);
    free(adj); free(akeys); free(face_comp); free(ov); free(is_ov);
    uf_free(&uf); free(region_of); free(root2rid);
    free(roff); free(rfaces); free(anoff); free(afaces);
    free(g2l); free(comp2layer); free(loclayer);
    return 0;
}

/* ============================================================================
 * Self-test
 * ==========================================================================*/

#define CK(c,m) do{ if(!(c)){ fprintf(stderr,"  FAIL: %s\n",(m)); nf++; } \
                    else fprintf(stderr,"  ok: %s\n",(m)); }while(0)

/* Build a flat rectangular sheet patch at a given radius, laid at UV
 * (u0..u0+W, v0..v0+H), 3D placed on a cylinder of radius `rad` around Z through
 * (0,0): vertex (z=v, y=rad*cos(u/rad), x=rad*sin(u/rad)). Appends. */
static void st_patch(double rad, double u0, double v0, int nu, int nv_,
                     float **V, int32_t **F, size_t *nv, size_t *nf_,
                     size_t *cv, size_t *cf) {
    size_t base = *nv;
    for (int j = 0; j <= nv_; j++)
        for (int i = 0; i <= nu; i++) {
            double u = u0 + i, v = v0 + j, th = u / rad;
            if (*nv+1 > *cv){ *cv = (*cv?*cv*2:1024); *V=(float*)realloc(*V,*cv*3*sizeof(float)); }
            (*V)[*nv*3+0]=(float)v; (*V)[*nv*3+1]=(float)(rad*cos(th)); (*V)[*nv*3+2]=(float)(rad*sin(th));
            (*nv)++;
        }
    /* uv stored separately by caller via a parallel array is awkward; instead
     * we return geometry and let the caller build uv = (u,v). We reconstruct u
     * from angle is lossy, so caller passes uv explicitly. Here we ALSO fill a
     * global uv via a static — simpler: caller builds uv. */
    for (int j = 0; j < nv_; j++)
        for (int i = 0; i < nu; i++) {
            int32_t a=(int32_t)(base + (size_t)j*(nu+1)+i), b=a+1, c=a+nu+1, d=c+1;
            if (*nf_+2 > *cf){ *cf=(*cf?*cf*2:1024); *F=(int32_t*)realloc(*F,*cf*3*sizeof(int32_t)); }
            (*F)[*nf_*3+0]=a;(*F)[*nf_*3+1]=b;(*F)[*nf_*3+2]=c;(*nf_)++;
            (*F)[*nf_*3+0]=b;(*F)[*nf_*3+1]=d;(*F)[*nf_*3+2]=c;(*nf_)++;
        }
}

/* Build the matching uv (u,v) for a patch grid appended by st_patch. */
static void st_uv(double u0, double v0, int nu, int nv_, float **T, size_t *nvt, size_t *ct) {
    for (int j = 0; j <= nv_; j++)
        for (int i = 0; i <= nu; i++) {
            if (*nvt+1 > *ct){ *ct=(*ct?*ct*2:1024); *T=(float*)realloc(*T,*ct*2*sizeof(float)); }
            (*T)[*nvt*2+0]=(float)(u0+i); (*T)[*nvt*2+1]=(float)(v0+j); (*nvt)++;
        }
}

int Quilt_selftest(void) {
    int nf = 0;
    Arena_T arena = Arena_new();
    fprintf(stderr, "[selftest] overlap_quilt\n");

    QuiltOpts o; QuiltOpts_default(&o);
    o.axis_point[0]=0; o.axis_point[1]=0; o.axis_point[2]=0;
    o.axis_dir[0]=1; o.axis_dir[1]=0; o.axis_dir[2]=0;   /* axis Z through origin */
    o.raw_dir = NULL; o.pitch = 10.0; o.radius_gate = 5.0; o.cell = 2.0; o.verbose = 0;

    /* t1: single-layer sheet -> no overlaps, no split */
    {
        float *V=NULL,*T=NULL; int32_t *F=NULL; size_t nv=0,nvt=0,nfc=0,cv=0,ct=0,cf=0;
        st_patch(50.0, 40.0, 0.0, 20, 8, &V,&F,&nv,&nfc,&cv,&cf);
        st_uv(40.0, 0.0, 20, 8, &T,&nvt,&ct);
        QuiltResult R;
        int rc = Quilt_run(arena, V, nv, F, nfc, T, &o, &R);
        CK(rc==0, "t1 run ok");
        CK(R.n_regions==0 && R.n_overlap_faces==0, "t1 single sheet: 0 overlap regions");
        CK(R.n_dropped==0, "t1 nothing dropped");
        free(V);free(T);free(F);
    }

    /* t2: two layers at SAME uv, radial gap 10 (> gate) -> overlap, 2 layers,
     * one dropped. Layer A rad 50, layer B rad 60, both at uv u=[40,60]. */
    {
        float *V=NULL,*T=NULL; int32_t *F=NULL; size_t nv=0,nvt=0,nfc=0,cv=0,ct=0,cf=0;
        st_patch(50.0, 40.0, 0.0, 20, 8, &V,&F,&nv,&nfc,&cv,&cf);
        st_uv(40.0, 0.0, 20, 8, &T,&nvt,&ct);
        st_patch(60.0, 40.0, 0.0, 20, 8, &V,&F,&nv,&nfc,&cv,&cf);
        st_uv(40.0, 0.0, 20, 8, &T,&nvt,&ct);
        QuiltResult R;
        int rc = Quilt_run(arena, V, nv, F, nfc, T, &o, &R);
        CK(rc==0, "t2 run ok");
        CK(R.n_regions>=1, "t2 overlap region formed");
        CK(R.n_layers_total>=2, "t2 multicut split into >=2 layers");
        /* no RAW here -> region is unresolvable -> nothing dropped (blend kept) */
        CK(R.n_dropped==0 && R.have_texture==0, "t2 no texture -> left intact");
        o.cross_component_only = 1;
        rc = Quilt_run(arena, V, nv, F, nfc, T, &o, &R);
        CK(rc==0 && R.n_regions==1 && R.n_layers_total==2,
           "t2 cross-component diagnostic preserves the two input sheets");
        o.cross_component_only = 0;
        free(V);free(T);free(F);
    }

    /* t5: quilting energy unit -- anchor holds a horizontal-fiber intensity
     * pattern (value depends on y); the layer that CONTINUES that pattern
     * across the seam scores lower than a phase-shifted (mis-aligned) one, and
     * a pure brightness offset must NOT fool it (offset removed). */
    {
        int gw=16, gh=16; size_t np=(size_t)gw*gh;
        float *vA=(float*)calloc(np,sizeof(float)), *vGood=(float*)calloc(np,sizeof(float)),
              *vBad=(float*)calloc(np,sizeof(float));
        uint8_t *covA=(uint8_t*)calloc(np,1), *covL=(uint8_t*)calloc(np,1);
        for(int y=0;y<gh;y++)for(int x=0;x<gw;x++){ size_t p=(size_t)y*gw+x;
            double fib = 100.0 + 40.0*sin((double)y);        /* horizontal fibers */
            if(x<8){ covA[p]=1; vA[p]=(float)fib; }
            else { covL[p]=1;
                vGood[p]=(float)(fib + 25.0);                 /* same pattern, +offset */
                vBad[p] =(float)(100.0 + 40.0*sin((double)y+2.2)); } }  /* phase-shifted */
        double Egood=Quilt_boundary_energy(vGood,covL,vA,covA,gw,gh);
        double Ebad =Quilt_boundary_energy(vBad, covL,vA,covA,gw,gh);
        CK(Egood<1e17 && Ebad<1e17, "t5 both layers share a seam");
        CK(Egood<Ebad, "t5 pattern-continuing layer wins (offset-robust)");
        free(vA);free(vGood);free(vBad);free(covA);free(covL);
    }

    /* t3: single folded sheet at SAME radius (gap 0) must NOT split (radius gate
     * self-gates) -- two patches at rad 50 sharing uv but same radius. */
    {
        float *V=NULL,*T=NULL; int32_t *F=NULL; size_t nv=0,nvt=0,nfc=0,cv=0,ct=0,cf=0;
        st_patch(50.0, 40.0, 0.0, 20, 8, &V,&F,&nv,&nfc,&cv,&cf);
        st_uv(40.0, 0.0, 20, 8, &T,&nvt,&ct);
        st_patch(50.0, 40.0, 0.0, 20, 8, &V,&F,&nv,&nfc,&cv,&cf);   /* same rad */
        st_uv(40.0, 0.0, 20, 8, &T,&nvt,&ct);
        QuiltResult R;
        int rc = Quilt_run(arena, V, nv, F, nfc, T, &o, &R);
        CK(rc==0, "t3 run ok");
        CK(R.n_overlap_faces==0 && R.n_dropped==0, "t3 same-radius: radius gate self-gates (no split)");
        free(V);free(T);free(F);
    }

    /* t4: degenerate */
    {
        QuiltResult R;
        CK(Quilt_run(arena, NULL, 0, NULL, 0, NULL, &o, &R) != 0, "t4 empty -> error");
    }

    Arena_dispose(&arena);
    fprintf(stderr, "[selftest] %s (%d failures)\n", nf==0?"ALL PASS":"FAILURES", nf);
    return nf;
}
