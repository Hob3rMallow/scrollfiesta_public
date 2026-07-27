/* seam_hole_fill.c -- see seam_hole_fill.h.
 *
 * Traces simple boundary loops, keeps the small ones that straddle a cube seam
 * plane and are winding-phase-coherent (same wrap), and closes them with a
 * manifold-guarded fan triangulation. Vertices are never added or moved.
 */
#include "seam_hole_fill.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "../common/arena.h"

void SeamHoleFill_default_params(SeamHoleFillParams *p)
{
    p->cube = 128.0;
    p->band = 6.0;
    p->max_extent = 6.5;   /* < 7-vox inter-wrap clearance (phase gate is the real guard) */
    p->max_loop = 8;
    p->umb_y = 0.0; p->umb_x = 0.0; p->pitch = 0.0;  /* phase gate off until set */
    p->wind_tol_turns = 0.25;
}

/* ---- undirected-edge -> incident-face-count hash (open addressing) ---- */
typedef struct { uint64_t *key; int32_t *cnt; size_t cap, mask; } EHash;

static uint64_t ekey(int32_t a, int32_t b)
{
    uint32_t lo = (uint32_t)(a < b ? a : b), hi = (uint32_t)(a < b ? b : a);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}
static void eh_init(Arena_T ar, EHash *h, size_t n)
{
    size_t cap = 16; while (cap < n * 2) cap <<= 1;
    h->cap = cap; h->mask = cap - 1;
    h->key = (uint64_t *)ARENA_ALLOC(ar, (long)(cap * sizeof(uint64_t)));
    h->cnt = (int32_t *)ARENA_ALLOC(ar, (long)(cap * sizeof(int32_t)));
    for (size_t i = 0; i < cap; i++) { h->key[i] = UINT64_MAX; h->cnt[i] = 0; }
}
static int32_t *eh_slot(EHash *h, uint64_t k)
{
    size_t i = (size_t)((k * 1099511628211ull) >> 20) & h->mask;
    for (;;) {
        if (h->key[i] == k) return &h->cnt[i];
        if (h->key[i] == UINT64_MAX) { h->key[i] = k; h->cnt[i] = 0; return &h->cnt[i]; }
        i = (i + 1) & h->mask;
    }
}
static int32_t eh_get(EHash *h, int32_t a, int32_t b)
{
    uint64_t k = ekey(a, b);
    size_t i = (size_t)((k * 1099511628211ull) >> 20) & h->mask;
    for (;;) {
        if (h->key[i] == k) return h->cnt[i];
        if (h->key[i] == UINT64_MAX) return 0;
        i = (i + 1) & h->mask;
    }
}

/* triangle unit normal (double); returns 0 if degenerate */
static int tri_normal(const float *v, int32_t a, int32_t b, int32_t c, double n[3])
{
    double e1[3] = { v[(size_t)b*3+0]-v[(size_t)a*3+0], v[(size_t)b*3+1]-v[(size_t)a*3+1], v[(size_t)b*3+2]-v[(size_t)a*3+2] };
    double e2[3] = { v[(size_t)c*3+0]-v[(size_t)a*3+0], v[(size_t)c*3+1]-v[(size_t)a*3+1], v[(size_t)c*3+2]-v[(size_t)a*3+2] };
    n[0] = e1[1]*e2[2]-e1[2]*e2[1];
    n[1] = e1[2]*e2[0]-e1[0]*e2[2];
    n[2] = e1[0]*e2[1]-e1[1]*e2[0];
    double L = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    if (L < 1e-9) return 0;
    n[0]/=L; n[1]/=L; n[2]/=L;
    return 1;
}
static double tri_min_alt(const float *v, int32_t a, int32_t b, int32_t c)
{
    double e0[3] = { v[(size_t)b*3+0]-v[(size_t)a*3+0], v[(size_t)b*3+1]-v[(size_t)a*3+1], v[(size_t)b*3+2]-v[(size_t)a*3+2] };
    double e1[3] = { v[(size_t)c*3+0]-v[(size_t)b*3+0], v[(size_t)c*3+1]-v[(size_t)b*3+1], v[(size_t)c*3+2]-v[(size_t)b*3+2] };
    double e2[3] = { v[(size_t)a*3+0]-v[(size_t)c*3+0], v[(size_t)a*3+1]-v[(size_t)c*3+1], v[(size_t)a*3+2]-v[(size_t)c*3+2] };
    double l0 = sqrt(e0[0]*e0[0]+e0[1]*e0[1]+e0[2]*e0[2]);
    double l1 = sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
    double l2 = sqrt(e2[0]*e2[0]+e2[1]*e2[1]+e2[2]*e2[2]);
    double lmax = l0>l1?(l0>l2?l0:l2):(l1>l2?l1:l2);
    if (lmax < 1e-12) return 0.0;
    double cr[3]; cr[0]=e0[1]*e2[2]-e0[2]*e2[1]; cr[1]=e0[2]*e2[0]-e0[0]*e2[2]; cr[2]=e0[0]*e2[1]-e0[1]*e2[0];
    double area = 0.5*sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
    return 2.0*area/lmax;
}

/* directed boundary edge (src->dst) + owning face, sorted by src for the walk */
typedef struct { int32_t s, d, face; } DEdge;

/* union-find over faces (path-halving) */
static int32_t uf_find(int32_t *par, int32_t x)
{
    while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; }
    return x;
}
static void uf_union(int32_t *par, int32_t a, int32_t b)
{
    a = uf_find(par, a); b = uf_find(par, b);
    if (a != b) par[a] = b;
}
static int cmp_dedge_und(const void *pa, const void *pb)
{
    const DEdge *a = pa, *b = pb;
    int32_t a0 = a->s<a->d?a->s:a->d, a1 = a->s<a->d?a->d:a->s;
    int32_t b0 = b->s<b->d?b->s:b->d, b1 = b->s<b->d?b->d:b->s;
    if (a0!=b0) return a0<b0?-1:1;
    if (a1!=b1) return a1<b1?-1:1;
    return 0;
}
static int cmp_dedge_src(const void *pa, const void *pb)
{
    const DEdge *a = pa, *b = pb;
    return a->s<b->s?-1:(a->s>b->s?1:0);
}

/* 2D orientation of (a,b,c): >0 CCW, <0 CW. */
static double cross2(const double a[2], const double b[2], const double c[2])
{
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
}
/* Is p inside triangle (a,b,c)? All same-sign barycentric cross products. */
static int pt_in_tri2(const double p[2], const double a[2],
                      const double b[2], const double c[2])
{
    double d1 = cross2(p,a,b), d2 = cross2(p,b,c), d3 = cross2(p,c,a);
    int neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    int pos = (d1 > 0) || (d2 > 0) || (d3 > 0);
    return !(neg && pos);
}

/* Ear-clip a small (3..8) boundary loop into a GENERAL (non-fan) triangulation.
 * Projects to the loop's Newell-normal plane, clips ears (convex, empty, non-
 * sliver), and -- crucially at a seam between sheets -- refuses any ear whose new
 * diagonal already exists as a mesh edge (`eh`), so the result is manifold-safe
 * by construction. This replaces the removed multi-apex fan: it is deterministic
 * and picks a triangulation that FITS the surrounding mesh rather than forcing a
 * star fan from an arbitrary apex. Triangles are reverse-wound to -sn (the
 * boundary-edge reversal). Returns the triangle count (n-2), or 0 if no manifold-
 * safe, non-folding triangulation exists (caller leaves the loop open). */
static int earclip_loop(const float *V, const int32_t *loop, int n,
                        EHash *eh, const double sn[3], int32_t *out)
{
    if (n < 3 || n > 8) return 0;
    int ax = (fabs(sn[0]) <= fabs(sn[1]) && fabs(sn[0]) <= fabs(sn[2])) ? 0
           : (fabs(sn[1]) <= fabs(sn[2]) ? 1 : 2);
    double t[3] = {0,0,0}; t[ax] = 1.0;
    double u[3] = { t[1]*sn[2]-t[2]*sn[1], t[2]*sn[0]-t[0]*sn[2], t[0]*sn[1]-t[1]*sn[0] };
    double ul = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
    if (ul < 1e-12) return 0;
    u[0]/=ul; u[1]/=ul; u[2]/=ul;
    double w[3] = { sn[1]*u[2]-sn[2]*u[1], sn[2]*u[0]-sn[0]*u[2], sn[0]*u[1]-sn[1]*u[0] };
    double P2[8][2];
    for (int k = 0; k < n; k++) {
        const float *pk = &V[(size_t)loop[k]*3];
        P2[k][0] = pk[0]*u[0]+pk[1]*u[1]+pk[2]*u[2];
        P2[k][1] = pk[0]*w[0]+pk[1]*w[1]+pk[2]*w[2];
    }
    double area2 = 0.0;
    for (int k = 0; k < n; k++) { int j=(k+1)%n; area2 += P2[k][0]*P2[j][1]-P2[j][0]*P2[k][1]; }
    double ori = area2 >= 0 ? 1.0 : -1.0;
    int rem[8], m = n;
    for (int k = 0; k < n; k++) rem[k] = k;
    int nt = 0, guard = 0;
    while (m > 3 && guard++ < 64) {
        int clipped = 0;
        for (int pos = 0; pos < m; pos++) {
            int ip = rem[(pos+m-1)%m], ii = rem[pos], iq = rem[(pos+1)%m];
            if (cross2(P2[ip], P2[ii], P2[iq]) * ori <= 0.0) continue;   /* reflex */
            if (eh_get(eh, loop[ip], loop[iq]) != 0) continue;          /* diagonal exists */
            if (tri_min_alt(V, loop[ip], loop[ii], loop[iq]) < 0.05) continue; /* sliver */
            int inside = 0;
            for (int s = 0; s < m; s++) {
                int iv = rem[s];
                if (iv == ip || iv == ii || iv == iq) continue;
                if (pt_in_tri2(P2[iv], P2[ip], P2[ii], P2[iq])) { inside = 1; break; }
            }
            if (inside) continue;
            out[nt*3+0] = loop[ip]; out[nt*3+1] = loop[ii]; out[nt*3+2] = loop[iq]; nt++;
            for (int s = pos; s < m-1; s++) rem[s] = rem[s+1];
            m--; clipped = 1; break;
        }
        if (!clipped) return 0;   /* no manifold-safe ear -> leave the loop open */
    }
    if (m == 3) { out[nt*3+0]=loop[rem[0]]; out[nt*3+1]=loop[rem[1]]; out[nt*3+2]=loop[rem[2]]; nt++; }
    double nref[3] = { -sn[0], -sn[1], -sn[2] };
    for (int t2 = 0; t2 < nt; t2++) {
        double nn[3];
        if (!tri_normal(V, out[t2*3+0], out[t2*3+1], out[t2*3+2], nn)) return 0;
        double d = nn[0]*nref[0]+nn[1]*nref[1]+nn[2]*nref[2];
        if (d < 0.0) { int32_t tmp = out[t2*3+1]; out[t2*3+1] = out[t2*3+2]; out[t2*3+2] = tmp; d = -d; }
        if (d < 0.2) return 0;    /* triangle folds relative to the surface */
    }
    return nt;
}

int SeamHoleFill_process(Arena_T arena, ComponentMesh *cm,
                         const SeamHoleFillParams *p, SeamHoleFillStats *st)
{
    SeamHoleFillStats z; memset(&z, 0, sizeof z);
    if (st) *st = z;
    if (!cm || cm->nf < 1 || cm->nv < 3) return 0;

    const float *V = cm->verts;
    const int32_t *F = cm->faces;
    size_t nf = cm->nf, nv = cm->nv;
    size_t hn = nf * 3;

    Arena_Mark mark = Arena_save(arena);

    /* edge -> incident-face count (live: updated as we commit fills) */
    EHash eh; eh_init(arena, &eh, hn);
    for (size_t f = 0; f < nf; f++) {
        int32_t a = F[f*3+0], b = F[f*3+1], c = F[f*3+2];
        (*eh_slot(&eh, ekey(a,b)))++;
        (*eh_slot(&eh, ekey(b,c)))++;
        (*eh_slot(&eh, ekey(c,a)))++;
    }

    /* directed boundary edges (undirected mult == 1) + per-vertex boundary degree.
     * Also union faces that share an interior (mult==2) edge, so we can later
     * refuse to cap the SOLE boundary loop of a component (a would-be bubble). */
    DEdge *he = (DEdge *)ARENA_ALLOC(arena, (long)(hn * sizeof(DEdge)));
    for (size_t f = 0; f < nf; f++) {
        int32_t a = F[f*3+0], b = F[f*3+1], c = F[f*3+2];
        he[f*3+0].s=a; he[f*3+0].d=b; he[f*3+0].face=(int32_t)f;
        he[f*3+1].s=b; he[f*3+1].d=c; he[f*3+1].face=(int32_t)f;
        he[f*3+2].s=c; he[f*3+2].d=a; he[f*3+2].face=(int32_t)f;
    }
    qsort(he, hn, sizeof(DEdge), cmp_dedge_und);
    int32_t *fpar = (int32_t *)ARENA_ALLOC(arena, (long)(nf * sizeof(int32_t)));
    for (size_t f = 0; f < nf; f++) fpar[f] = (int32_t)f;
    DEdge *bnd = (DEdge *)ARENA_ALLOC(arena, (long)(hn * sizeof(DEdge)));
    size_t nb = 0;
    int32_t *bdeg = (int32_t *)ARENA_CALLOC(arena, (long)nv, (long)sizeof(int32_t));
    for (size_t i = 0; i < hn;) {
        size_t j = i + 1;
        int32_t a0 = he[i].s<he[i].d?he[i].s:he[i].d, a1 = he[i].s<he[i].d?he[i].d:he[i].s;
        while (j < hn) {
            int32_t b0 = he[j].s<he[j].d?he[j].s:he[j].d, b1 = he[j].s<he[j].d?he[j].d:he[j].s;
            if (b0!=a0 || b1!=a1) break; j++;
        }
        if (j - i == 1) { bnd[nb++] = he[i]; bdeg[he[i].s]++; bdeg[he[i].d]++; }
        else if (j - i == 2) { uf_union(fpar, he[i].face, he[i+1].face); }
        i = j;
    }
    if (nb == 0) { Arena_restore(arena, mark); if (st) *st = z; return 0; }

    /* sort boundary edges by src for the walk */
    qsort(bnd, nb, sizeof(DEdge), cmp_dedge_src);
    int32_t *src = (int32_t *)ARENA_ALLOC(arena, (long)(nb * sizeof(int32_t)));
    for (size_t i = 0; i < nb; i++) src[i] = bnd[i].s;
    uint8_t *used = (uint8_t *)ARENA_CALLOC(arena, (long)nb, 1L);

    /* scratch for one loop's ordered vertices + the appended fill triangles */
    int32_t *loop = (int32_t *)ARENA_ALLOC(arena, (long)((p->max_loop + 2) * (long)sizeof(int32_t)));
    /* fill accumulator: at most (max_loop-2) tris per loop; grow generously */
    size_t fill_cap = nf + 4096, fill_n = 0;
    int32_t *fill = (int32_t *)ARENA_ALLOC(arena, (long)(fill_cap * 3 * (long)sizeof(int32_t)));

    const double tol = p->wind_tol_turns;
    const int phase_on = (p->pitch > 0.0 && (p->umb_y != 0.0 || p->umb_x != 0.0));
    const double ext_cap = phase_on ? p->max_extent : (p->max_extent < 4.0 ? p->max_extent : 4.0);

    #define FIND_SRC(key, LO, HI) do {                                       \
        size_t _l=0,_r=nb; while(_l<_r){size_t _m=(_l+_r)/2; if(src[_m]<(key))_l=_m+1; else _r=_m;} (LO)=_l; \
        size_t _l2=0,_r2=nb; while(_l2<_r2){size_t _m=(_l2+_r2)/2; if(src[_m]<=(key))_l2=_m+1; else _r2=_m;} (HI)=_l2; \
    } while (0)

    /* Count boundary loops per component (sole-loop bubble guard). A genuine
     * puncture lives on a component that also owns a larger boundary (>=2 loops);
     * a tiny isolated patch whose ONLY boundary is this loop must NOT be capped
     * into a closed bubble. Cheap tracer pass over a scratch used-map. */
    int32_t *loops_per_comp = (int32_t *)ARENA_CALLOC(arena, (long)nf, (long)sizeof(int32_t));
    uint8_t *used2 = (uint8_t *)ARENA_CALLOC(arena, (long)nb, 1L);
    for (size_t s0 = 0; s0 < nb; s0++) {
        if (used2[s0]) continue;
        size_t cur = s0;
        int32_t root = uf_find(fpar, bnd[s0].face);
        int32_t start_v = bnd[s0].s;
        for (;;) {
            if (used2[cur]) break;
            used2[cur] = 1;
            int32_t dv = bnd[cur].d;
            if (dv == start_v) break;
            size_t lo, hi; FIND_SRC(dv, lo, hi);
            size_t nxt = (size_t)-1;
            for (size_t t = lo; t < hi; t++) if (!used2[t]) { nxt = t; break; }
            if (nxt == (size_t)-1) break;
            cur = nxt;
        }
        loops_per_comp[root]++;
    }

    for (size_t s0 = 0; s0 < nb; s0++) {
        if (used[s0]) continue;
        int32_t comp_root = uf_find(fpar, bnd[s0].face);
        /* walk the loop, collecting ordered verts; bail if non-simple or too long */
        size_t cur = s0, n = 0; int simple = 1;
        int32_t start_v = bnd[s0].s;
        for (;;) {
            if (used[cur]) { break; }
            used[cur] = 1;
            int32_t sv = bnd[cur].s, dv = bnd[cur].d;
            if (bdeg[sv] != 2) simple = 0;
            if (n < (size_t)p->max_loop + 1) loop[n] = sv;
            n++;
            if (dv == start_v) break;
            size_t lo, hi; FIND_SRC(dv, lo, hi);
            size_t nxt = (size_t)-1;
            for (size_t t = lo; t < hi; t++) if (!used[t]) { nxt = t; break; }
            if (nxt == (size_t)-1) { simple = 0; break; }
            cur = nxt;
            if (n > (size_t)p->max_loop + 1) break;   /* too long */
        }
        if (!simple) continue;
        z.loops++;
        if (n < 3 || n > (size_t)p->max_loop) continue;

        /* bbox + extent */
        double mn[3] = { 1e30,1e30,1e30 }, mx[3] = { -1e30,-1e30,-1e30 };
        for (size_t k = 0; k < n; k++) for (int a = 0; a < 3; a++) {
            double c = V[(size_t)loop[k]*3+(size_t)a];
            if (c < mn[a]) mn[a] = c; if (c > mx[a]) mx[a] = c;
        }
        double ext = sqrt((mx[0]-mn[0])*(mx[0]-mn[0]) + (mx[1]-mn[1])*(mx[1]-mn[1]) + (mx[2]-mn[2])*(mx[2]-mn[2]));

        /* near-seam + straddle on some axis */
        int on_seam = 0;
        for (int a = 0; a < 3; a++) {
            double mid = 0.5*(mn[a]+mx[a]);
            double coord = p->cube * floor(mid / p->cube + 0.5);
            if (mn[a] >= coord - p->band && mx[a] <= coord + p->band &&
                mn[a] < coord && mx[a] > coord) { on_seam = 1; break; }
        }
        if (!on_seam) continue;
        z.seam_candidates++;

        /* sole-loop bubble guard: never cap a component's only boundary loop. */
        if (loops_per_comp[comp_root] < 2) { z.skip_bubble++; continue; }

        if (ext > ext_cap) { z.skip_extent++; continue; }

        /* winding-phase gate: consecutive-edge wrapped dw, same as the seam gate */
        if (phase_on) {
            double span = 0.0;
            for (size_t k = 0; k < n; k++) {
                int32_t a = loop[k], b = loop[(k+1)%n];
                double dya = V[(size_t)a*3+1]-p->umb_y, dxa = V[(size_t)a*3+2]-p->umb_x;
                double dyb = V[(size_t)b*3+1]-p->umb_y, dxb = V[(size_t)b*3+2]-p->umb_x;
                double ra = hypot(dya,dxa), rb = hypot(dyb,dxb);
                double dth = atan2(dya,dxa) - atan2(dyb,dxb);
                while (dth >  M_PI) dth -= 2.0*M_PI;
                while (dth < -M_PI) dth += 2.0*M_PI;
                double dw = fabs((ra-rb)/p->pitch - dth/(2.0*M_PI));
                if (dw > span) span = dw;
            }
            if (span > tol) { z.skip_phase++; continue; }
        }

        /* surface normal estimate: the loop polygon's area vector about its
         * centroid (Newell). For a nearly-planar small hole this is the surface
         * normal up to sign; the fan winding is validated against it below. */
        double sn[3] = {0,0,0};
        {
            double C[3] = {0,0,0};
            for (size_t k = 0; k < n; k++) for (int a=0;a<3;a++) C[a]+=V[(size_t)loop[k]*3+(size_t)a];
            for (int a=0;a<3;a++) C[a]/=(double)n;
            for (size_t k = 0; k < n; k++) {
                const float *pa = &V[(size_t)loop[k]*3], *pb = &V[(size_t)loop[(k+1)%n]*3];
                double qa[3]={pa[0]-C[0],pa[1]-C[1],pa[2]-C[2]}, qb[3]={pb[0]-C[0],pb[1]-C[1],pb[2]-C[2]};
                sn[0]+=qa[1]*qb[2]-qa[2]*qb[1];
                sn[1]+=qa[2]*qb[0]-qa[0]*qb[2];
                sn[2]+=qa[0]*qb[1]-qa[1]*qb[0];
            }
            double L = sqrt(sn[0]*sn[0]+sn[1]*sn[1]+sn[2]*sn[2]);
            if (L < 1e-9) { z.skip_geom++; continue; }
            sn[0]/=L; sn[1]/=L; sn[2]/=L;
        }

        /* Ear-clip triangulation (reverse-wound to -sn, the boundary-edge
         * reversal, so the fill is manifold-consistent with the neighbour faces).
         * The earlier MULTI-APEX FAN -- try EVERY vertex as apex, keep the first
         * fan that clears the gates -- was removed: it forced a star fill on any
         * loop shape, the wrong instinct at a seam between sheets. earclip_loop
         * instead picks a general triangulation that FITS the surrounding mesh
         * (it refuses any diagonal that already exists), and leaves the loop OPEN
         * if none is manifold-safe. */
        size_t ntri = n - 2;
        if (ntri > 6) { continue; }  /* guard the fixed buffers (max_loop<=8) */
        int32_t tri[3 * 6];
        int ntri_ec = earclip_loop(V, loop, (int)n, &eh, sn, tri);
        if (ntri_ec <= 0 || (size_t)ntri_ec != ntri) { z.skip_geom++; continue; }
        /* manifold check: for every edge the fill touches, require
         * existing_incident_faces + fill_uses <= 2 (boundary edge 1+1=2 ok;
         * internal diagonal 0+2=2 ok, but a diagonal already in the mesh ->
         * 1+2=3 rejected). Tally over the fill (<=18 edges) vs the live hash. */
        {
            uint64_t ek[3*6]; int euse[3*6]; size_t ne = 0;
            for (size_t t = 0; t < ntri; t++) {
                int32_t vtx[3] = { tri[t*3+0], tri[t*3+1], tri[t*3+2] };
                for (int e = 0; e < 3; e++) {
                    uint64_t kk = ekey(vtx[e], vtx[(e+1)%3]);
                    size_t g = 0; for (; g < ne; g++) if (ek[g]==kk) break;
                    if (g==ne) { ek[ne]=kk; euse[ne]=1; ne++; } else euse[g]++;
                }
            }
            int bad = 0;
            for (size_t g = 0; g < ne && !bad; g++) {
                int32_t a = (int32_t)(ek[g] >> 32), b = (int32_t)(ek[g] & 0xffffffffu);
                if (eh_get(&eh, a, b) + euse[g] > 2) bad = 1;
            }
            if (bad) { z.skip_manifold++; continue; }   /* non-manifold; leave open */
        }

        /* commit: append tris + update live edge counts */
        if (fill_n + ntri > fill_cap) {
            size_t ncap = (fill_n + ntri) * 2;
            int32_t *nf2 = (int32_t *)ARENA_ALLOC(arena, (long)(ncap * 3 * (long)sizeof(int32_t)));
            memcpy(nf2, fill, fill_n * 3 * sizeof(int32_t));
            fill = nf2; fill_cap = ncap;
        }
        for (size_t t = 0; t < ntri; t++) {
            int32_t a=tri[t*3+0], b=tri[t*3+1], c=tri[t*3+2];
            fill[(fill_n)*3+0]=a; fill[(fill_n)*3+1]=b; fill[(fill_n)*3+2]=c; fill_n++;
            (*eh_slot(&eh, ekey(a,b)))++;
            (*eh_slot(&eh, ekey(b,c)))++;
            (*eh_slot(&eh, ekey(c,a)))++;
        }
        z.filled++; z.tris_added += ntri;
    }
    #undef FIND_SRC

    /* append fills to the mesh (arena-realloc; must outlive `mark`, so copy the
     * combined face array into a fresh alloc that we do NOT restore). */
    if (fill_n > 0) {
        size_t nf_new = nf + fill_n;
        int32_t *out = (int32_t *)ARENA_ALLOC(arena, (long)(nf_new * 3 * (long)sizeof(int32_t)));
        memcpy(out, F, nf * 3 * sizeof(int32_t));
        memcpy(out + nf * 3, fill, fill_n * 3 * sizeof(int32_t));
        cm->faces = out;
        cm->nf = nf_new;
        /* scratch (eh, he, bnd, ...) is below `out` in the arena? No -- `out` was
         * allocated last, after all scratch, so we cannot restore past it without
         * freeing `out`. Leave the arena as-is (scratch leaks until arena dispose);
         * grid_weld disposes the arena per weld, so this is bounded. */
    } else {
        Arena_restore(arena, mark);
    }

    if (st) *st = z;
    return 0;
}
