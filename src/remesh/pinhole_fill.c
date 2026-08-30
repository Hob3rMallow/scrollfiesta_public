#include "pinhole_fill.h"

#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/* See pinhole_fill.h for the algorithm. Two phases per component:
 *   (1) split_pinch_verts  — duplicate bowtie verts, one copy per fan;
 *   (2) fill_small_loops   — close EXACT single-triangle (3-loop) holes. */

static size_t next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) { p <<= 1; }
    return p;
}

static uint64_t edge_hash(int32_t lo, int32_t hi) {
    uint64_t h = (uint64_t)(uint32_t)lo * 0x9E3779B97F4A7C15ULL
               + (uint64_t)(uint32_t)hi;
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33;
    return h;
}

/* Unordered-triple key for a triangle-existence set. The stored value IS the
 * hash (a hash-set, not a key/value map); a 0 hash is bumped to 1 so 0 can be
 * the empty sentinel. Collisions only ever cause a fill to be conservatively
 * skipped, never a bad fill, so a 64-bit hash is safe. */
static uint64_t tri_key(int32_t a, int32_t b, int32_t c) {
    int32_t t;
    if (a > b) { t = a; a = b; b = t; }
    if (b > c) { t = b; b = c; c = t; }
    if (a > b) { t = a; a = b; b = t; }
    uint64_t h = 1469598103934665603ULL;
    h = (h ^ (uint64_t)(uint32_t)a) * 1099511628211ULL;
    h = (h ^ (uint64_t)(uint32_t)b) * 1099511628211ULL;
    h = (h ^ (uint64_t)(uint32_t)c) * 1099511628211ULL;
    return h ? h : 1;
}

static int tri_has(const uint64_t *tk, uint64_t mask,
                   int32_t a, int32_t b, int32_t c) {
    uint64_t key = tri_key(a, b, c), s = key & mask;
    while (tk[s] != 0) {
        if (tk[s] == key) { return 1; }
        s = (s + 1) & mask;
    }
    return 0;
}

static void tri_add(uint64_t *tk, uint64_t mask,
                    int32_t a, int32_t b, int32_t c) {
    uint64_t key = tri_key(a, b, c), s = key & mask;
    while (tk[s] != 0) {
        if (tk[s] == key) { return; }
        s = (s + 1) & mask;
    }
    tk[s] = key;
}

static int uf_find(int *uf, int a) {
    while (uf[a] != a) { uf[a] = uf[uf[a]]; a = uf[a]; }
    return a;
}

/* Undirected edge -> face count, plus first incident face (for the face
 * connected-component union-find). Open-chained in a hash table. */
typedef struct { int32_t lo, hi, cnt, next, f0; } E;

/* Combined source + raw-bridge edge state used by the chart-level bowtie
 * closer.  This is deliberately separate from E above: the regular pinhole
 * pass owns a complete mesh, whereas a bridge repair must preserve the source
 * chart boundary that the bridge is going to attach to. */
typedef struct {
    int32_t lo, hi;
    uint32_t count;
    int32_t orient_sum;
    uint32_t source_count;
    uint8_t used;
} BowtieGuardEdge;

typedef struct {
    BowtieGuardEdge *slot;
    size_t cap;
} BowtieGuard;

static int bowtie_guard_init(BowtieGuard *g, size_t edge_hint)
{
    size_t cap = next_pow2(edge_hint > 4 ? edge_hint * 2 : 8);
    g->slot = (BowtieGuardEdge *)calloc(cap, sizeof(*g->slot));
    if (g->slot == NULL) { g->cap = 0; return -1; }
    g->cap = cap;
    return 0;
}

static BowtieGuardEdge *bowtie_guard_slot(BowtieGuard *g,
                                           int32_t a, int32_t b)
{
    int32_t lo = a < b ? a : b, hi = a < b ? b : a;
    size_t s = (size_t)(edge_hash(lo, hi) & (uint64_t)(g->cap - 1));
    while (g->slot[s].used &&
           (g->slot[s].lo != lo || g->slot[s].hi != hi))
        s = (s + 1) & (g->cap - 1);
    return &g->slot[s];
}

static void bowtie_guard_add_edge(BowtieGuard *g, int32_t a, int32_t b,
                                  int is_source)
{
    BowtieGuardEdge *e = bowtie_guard_slot(g, a, b);
    if (!e->used) {
        e->used = 1;
        e->lo = a < b ? a : b;
        e->hi = a < b ? b : a;
    }
    e->count++;
    e->orient_sum += a < b ? 1 : -1;
    if (is_source) e->source_count++;
}

static void bowtie_guard_add_faces(BowtieGuard *g,
                                    const int32_t *faces, size_t nf,
                                    int is_source)
{
    for (size_t f = 0; f < nf; f++)
        for (int k = 0; k < 3; k++)
            bowtie_guard_add_edge(g, faces[f*3+(size_t)k],
                                  faces[f*3+(size_t)((k+1)%3)],is_source);
}

/* Test one proposed triangle against the immutable source/raw-bridge edge map
 * plus triangles tentatively selected at this vertex.  An occupied edge must
 * be used exactly once in the opposite direction; a two-face edge is sealed. */
static int bowtie_guard_triangle_ok(BowtieGuard *g,
                                    const int32_t *pending,
                                    size_t pending_begin, size_t pending_end,
                                    int32_t a, int32_t b, int32_t c,
                                    int defer_source_orientation)
{
    int32_t tri[3] = {a,b,c};
    for (int k = 0; k < 3; k++) {
        int32_t x = tri[k], y = tri[(k+1)%3];
        int32_t lo = x < y ? x : y, hi = x < y ? y : x;
        int32_t dir = x < y ? 1 : -1;
        BowtieGuardEdge *e = bowtie_guard_slot(g, lo, hi);
        uint32_t count = e->used ? e->count : 0;
        int32_t orient_sum = e->used ? e->orient_sum : 0;
        uint32_t pending_count = 0;
        for (size_t f = pending_begin; f < pending_end; f++) {
            for (int q = 0; q < 3; q++) {
                int32_t px = pending[f*3+(size_t)q];
                int32_t py = pending[f*3+(size_t)((q+1)%3)];
                int32_t plo = px < py ? px : py;
                int32_t phi = px < py ? py : px;
                if (plo == lo && phi == hi) {
                    count++;
                    pending_count++;
                    orient_sum += px < py ? 1 : -1;
                }
            }
        }
        if (count >= 2) return 0;
        if (count == 1 && orient_sum != -dir) {
            /* Raw BPA patches are not oriented against the source until the
             * chart transaction is admitted.  A candidate triangle must keep
             * the raw bridge fan internally orientable, but it may traverse a
             * source-only boundary edge in the same direction: the forest can
             * then flip the complete patch as one unit.  It will reject the
             * patch if different source attachments demand inconsistent flips.
             * Never waive a bridge/pending conflict or edge multiplicity. */
            if (!(defer_source_orientation && pending_count == 0 &&
                  e->used && e->source_count == 1 && e->count == 1))
                return 0;
        }
    }
    return 1;
}

/* Two vertices at the same position (< 1e-4 vox apart). The bowtie split in
 * phase 1 duplicates a vertex's position, so a pinched loop carries such a
 * pair; a fill triangle using BOTH would be zero-area. */
static int coincident_v(const float *vp, int a, int b) {
    float dz = vp[a * 3 + 0] - vp[b * 3 + 0];
    float dy = vp[a * 3 + 1] - vp[b * 3 + 1];
    float dx = vp[a * 3 + 2] - vp[b * 3 + 2];
    return (dz * dz + dy * dy + dx * dx) < 1e-8f;
}

/* Min altitude (2*area / longest edge) of triangle (a,b,c), in voxels: how thin
 * the triangle is. ~0 for a (near-)collinear triple; a healthy ~0.7-vox-spaced
 * fill triangle is ~0.4. The exact coincident_v test (1e-4 vox) misses a
 * collinear-but-not-coincident triple — the "chewed" boundary along a cube face
 * is a row of ~1-vox-spaced collinear verts whose fill is a zero-area sliver.
 * Reject below HOLEFILL_MIN_ALT_VOX. */
static float tri_min_altitude(const float *vp, int a, int b, int c) {
    float az = vp[a*3+0], ay = vp[a*3+1], ax = vp[a*3+2];
    float bz = vp[b*3+0], by = vp[b*3+1], bx = vp[b*3+2];
    float cz = vp[c*3+0], cy = vp[c*3+1], cx = vp[c*3+2];
    float e1z = bz-az, e1y = by-ay, e1x = bx-ax;
    float e2z = cz-az, e2y = cy-ay, e2x = cx-ax;
    float crz = e1y*e2x - e1x*e2y;
    float cry = e1x*e2z - e1z*e2x;
    float crx = e1z*e2y - e1y*e2z;
    float area2 = sqrtf(crz*crz + cry*cry + crx*crx);     /* = 2 * area */
    float l1 = e1z*e1z + e1y*e1y + e1x*e1x;
    float l2 = e2z*e2z + e2y*e2y + e2x*e2x;
    float e3z = cz-bz, e3y = cy-by, e3x = cx-bx;
    float l3 = e3z*e3z + e3y*e3y + e3x*e3x;
    float lmax2 = l1 > l2 ? (l1 > l3 ? l1 : l3) : (l2 > l3 ? l2 : l3);
    if (lmax2 < 1e-12f) { return 0.0f; }
    return area2 / sqrtf(lmax2);
}

/* Unit normal of triangle (a,b,c) into n[3]; returns 2*area (0 if degenerate). */
static double face_normal_d(const float *vp, int32_t a, int32_t b, int32_t c,
                            double n[3]) {
    double az=vp[(size_t)a*3+0], ay=vp[(size_t)a*3+1], ax=vp[(size_t)a*3+2];
    double bz=vp[(size_t)b*3+0], by=vp[(size_t)b*3+1], bx=vp[(size_t)b*3+2];
    double cz=vp[(size_t)c*3+0], cy=vp[(size_t)c*3+1], cx=vp[(size_t)c*3+2];
    double e1[3]={bz-az,by-ay,bx-ax}, e2[3]={cz-az,cy-ay,cx-ax};
    n[0]=e1[1]*e2[2]-e1[2]*e2[1];
    n[1]=e1[2]*e2[0]-e1[0]*e2[2];
    n[2]=e1[0]*e2[1]-e1[1]*e2[0];
    double L=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    if (L>1e-20){ n[0]/=L; n[1]/=L; n[2]/=L; }
    return L;
}

/* Build vertex -> incident-face CSR (off[nv+1], inc[3*nf]). */
static void build_vert_faces(Arena_T arena, const int32_t *faces, size_t nf,
                             size_t nv, size_t **out_off, int **out_inc) {
    size_t *off = (size_t *)ARENA_CALLOC(arena, nv + 1, sizeof(size_t));
    for (size_t f = 0; f < nf; f++) {
        for (int k = 0; k < 3; k++) { off[(size_t)faces[f * 3 + k] + 1]++; }
    }
    for (size_t v = 0; v < nv; v++) { off[v + 1] += off[v]; }
    int *inc = (int *)ARENA_ALLOC(arena, off[nv] * sizeof(int));
    size_t *cur = (size_t *)ARENA_ALLOC(arena, nv * sizeof(size_t));
    memcpy(cur, off, nv * sizeof(size_t));
    for (size_t f = 0; f < nf; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t v = faces[f * 3 + k];
            inc[cur[v]++] = (int)f;
        }
    }
    *out_off = off; *out_inc = inc;
}

/* Direction face fi (its 3 verts) traverses undirected edge (v,w) in its winding:
 * +1 if it goes v->w, -1 if w->v, 0 if the edge is not in the face. Used by the
 * bowtie orientation gate to decide whether a gap-filling fan can reverse both
 * fans' boundary half-edges (benign pinch) or cannot (a fold). */
static int face_edge_dir(const int32_t *faces, int fi, int32_t v, int32_t w)
{
    for (int k = 0; k < 3; k++) {
        int32_t x = faces[fi*3+k], y = faces[fi*3+(k+1)%3];
        if (x == v && y == w) return +1;
        if (x == w && y == v) return -1;
    }
    return 0;
}

/* PHASE 0 — close bowtie pinches IN PLACE (keep the vertex).
 *
 * A non-manifold "bowtie" vertex has incident faces forming >=2 fans that share
 * only v. split_pinch_verts (phase 1) resolves it by DUPLICATING v, one copy per
 * fan -- but at a near-coincident cross-cube seam that spawns a duplicate vertex
 * at an existing position, and phase 2 then lays a fill triangle over the
 * neighbour's coincident surface => a doubled, z-fighting patch (see seam_audit).
 *
 * Instead, when the fans form ONE well-defined surface (every incident face's
 * normal agrees with the area-weighted vertex normal within BOWTIE_COHERENCE_COS)
 * and meet across only NARROW gaps (< BOWTIE_GAP_MAX_VOX), join the fans with a
 * shortest safe spanning forest of triangles fanned from v.  The forest merges
 * the fans into one manifold disk while retaining an exterior boundary; filling
 * every cyclic gap would instead create an annulus.  No vertex is duplicated
 * and no coincident triangle is created. A bowtie that is incoherent (divergent
 * sheets) or has no safe spanning set (for example two wraps touching at a
 * point) fails a gate and is left untouched for split_pinch_verts -- so wraps
 * are never welded. Adds fill triangles to cm->faces; nv is unchanged. */
static int close_bowtie_gaps(Arena_T arena, ComponentMesh *cm,
                             const int32_t *source_faces, size_t source_nf,
                             const int32_t *vertex_chart, size_t n_charts,
                             const size_t *chart_face_count,
                             int bridge_boundary_mode,
                             int defer_source_orientation,
                             size_t *out_closed) {
    size_t nv = cm->nv, nf = cm->nf;
    int32_t *faces = cm->faces;
    const float *vp = cm->verts;

    size_t *off = NULL; int *inc = NULL;
    build_vert_faces(arena, faces, nf, nv, &off, &inc);
    int *uf = (int *)ARENA_ALLOC(arena, nf * sizeof(int));

    /* fill triangles (original vertex indices); <= total incidence slots. */
    int32_t *fill = (int32_t *)ARENA_ALLOC(
        arena, (off[nv] ? off[nv] : 1) * 3 * sizeof(int32_t));
    size_t nfill = 0, closed = 0;
    int dbg = (getenv("PINHOLE_DEBUG") != NULL);
    BowtieGuard guard;
    memset(&guard, 0, sizeof guard);
    if (bridge_boundary_mode) {
        size_t edge_hint =
            (bridge_boundary_mode==1 ? source_nf+nf : nf)*3+1;
        if (bowtie_guard_init(&guard, edge_hint) != 0) return -1;
        if(bridge_boundary_mode==1){
            bowtie_guard_add_faces(&guard,source_faces,source_nf,1);
            bowtie_guard_add_faces(&guard,faces,nf,0);
        }else if(bridge_boundary_mode==2){
            bowtie_guard_add_faces(&guard,faces,source_nf,1);
            bowtie_guard_add_faces(&guard,faces+source_nf*3,nf-source_nf,0);
        }
    }

    /* The weld arms a wider gap cap (via SEAM_WRAP_PITCH); wide gaps are then
     * gated by the ORIENTATION test below, not size. Per-cube (unarmed) keeps the
     * conservative 2.0 cap. */
    int bt_armed = 0;
    { const char *e = getenv("SEAM_WRAP_PITCH"); if (e && atof(e) > 0.0) bt_armed = 1; }
    double bt_gapmax = bt_armed ? (double)BOWTIE_GAP_MAX_ARMED_VOX
                                : (double)BOWTIE_GAP_MAX_VOX;

    enum { NBR_CAP = 128 };
    for (size_t v = 0; v < nv; v++) {
        size_t s = off[v], e = off[v + 1];
        double coherence_gate=(double)BOWTIE_COHERENCE_COS;
        if (e - s < 2) { continue; }
        if(bridge_boundary_mode==2&&vertex_chart!=NULL&&
           chart_face_count!=NULL){
            size_t smallest=SIZE_MAX;
            for(size_t j=s;j<e;j++)for(int k=0;k<3;k++){
                int32_t q=faces[(size_t)inc[j]*3+(size_t)k];
                int32_t c=vertex_chart[q];
                if(c>=0&&(size_t)c<n_charts&&
                   chart_face_count[c]<smallest)
                    smallest=chart_face_count[c];
            }
            if(smallest<=BRIDGE_MICRO_CHART_MAX_FACES)
                coherence_gate=
                    (double)BRIDGE_MICRO_ATTACHMENT_COHERENCE_COS;
        }

        /* fan union-find: union incident faces that share a non-v vertex (i.e.
         * an edge (v,w)). Same fan detection as split_pinch_verts. */
        for (size_t j = s; j < e; j++) { uf[inc[j]] = inc[j]; }
        for (size_t j = s; j < e; j++) {
            int fi = inc[j], o1 = -1, o2 = -1;
            for (int k = 0; k < 3; k++) {
                int32_t w = faces[fi * 3 + k];
                if ((size_t)w != v) { if (o1 < 0) { o1 = w; } else { o2 = w; } }
            }
            for (size_t j2 = j + 1; j2 < e; j2++) {
                int gj = inc[j2], p1 = -1, p2 = -1;
                for (int k = 0; k < 3; k++) {
                    int32_t w = faces[gj * 3 + k];
                    if ((size_t)w != v) { if (p1 < 0) { p1 = w; } else { p2 = w; } }
                }
                if (o1 == p1 || o1 == p2 || o2 == p1 || o2 == p2) {
                    int ra = uf_find(uf, fi), rb = uf_find(uf, gj);
                    if (ra != rb) { uf[ra] = rb; }
                }
            }
        }
        size_t fans = 0;
        for (size_t j = s; j < e; j++) { if (uf_find(uf, inc[j]) == inc[j]) { fans++; } }
        if (fans < 2) { continue; }                 /* manifold vertex */

        /* area-weighted vertex normal + coherence (min face-normal . n_hat). */
        double n[3] = {0,0,0};
        for (size_t j = s; j < e; j++) {
            int fi = inc[j]; double fn[3];
            double a2 = face_normal_d(vp, faces[fi*3+0], faces[fi*3+1], faces[fi*3+2], fn);
            n[0]+=fn[0]*a2; n[1]+=fn[1]*a2; n[2]+=fn[2]*a2;
        }
        double nl = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if (nl < 1e-12) {                            /* normal cancels -> split */
            if (dbg) fprintf(stderr, "  [bowtie] v%zu (%.2f %.2f %.2f) fans=%zu SKIP normal-cancel\n",
                             v, (double)vp[v*3+0], (double)vp[v*3+1], (double)vp[v*3+2], fans);
            continue;
        }
        n[0]/=nl; n[1]/=nl; n[2]/=nl;
        double cohmin = 1.0;
        for (size_t j = s; j < e; j++) {
            int fi = inc[j]; double fn[3];
            face_normal_d(vp, faces[fi*3+0], faces[fi*3+1], faces[fi*3+2], fn);
            double d = fn[0]*n[0]+fn[1]*n[1]+fn[2]*n[2];
            if (d < cohmin) { cohmin = d; }
        }
        if (cohmin < coherence_gate) {              /* divergent -> split */
            if (dbg) fprintf(stderr, "  [bowtie] v%zu (%.2f %.2f %.2f) fans=%zu SKIP incoherent cohmin=%.3f\n",
                             v, (double)vp[v*3+0], (double)vp[v*3+1], (double)vp[v*3+2], fans, cohmin);
            continue;
        }

        /* tangent basis (u,w) perpendicular to n, for the angular edge order. */
        int axmin = (fabs(n[0])<=fabs(n[1]) && fabs(n[0])<=fabs(n[2])) ? 0
                  : (fabs(n[1])<=fabs(n[2]) ? 1 : 2);
        double t[3] = {0,0,0}; t[axmin] = 1.0;
        double u[3] = { t[1]*n[2]-t[2]*n[1], t[2]*n[0]-t[0]*n[2], t[0]*n[1]-t[1]*n[0] };
        double ul = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
        if (ul < 1e-12) { continue; }
        u[0]/=ul; u[1]/=ul; u[2]/=ul;
        double wv3[3] = { n[1]*u[2]-n[2]*u[1], n[2]*u[0]-n[0]*u[2], n[0]*u[1]-n[1]*u[0] };

        /* tally non-v neighbours among v's faces: count + the fan they belong to. */
        int nbr[NBR_CAP], ncnt[NBR_CAP], nfan[NBR_CAP], nn = 0, overflow = 0;
        for (size_t j = s; j < e && !overflow; j++) {
            int fi = inc[j], froot = uf_find(uf, fi);
            for (int k = 0; k < 3; k++) {
                int32_t w = faces[fi*3+k];
                if ((size_t)w == v) { continue; }
                int idx = -1;
                for (int q = 0; q < nn; q++) { if (nbr[q] == w) { idx = q; break; } }
                if (idx >= 0) { ncnt[idx]++; }
                else if (nn < NBR_CAP) { nbr[nn]=w; ncnt[nn]=1; nfan[nn]=froot; nn++; }
                else { overflow = 1; break; }
            }
        }
        if (overflow) {
            if (dbg) fprintf(stderr, "  [bowtie] v%zu fans=%zu SKIP neighbor-overflow\n", v, fans);
            continue;
        }

        /* boundary edges at v = neighbours used by exactly ONE incident face.
         * bdir = the direction that owning face traverses (v,w): +1 v->w, -1 w->v.
         * The orientation gate uses it to tell a benign pinch from a fold. */
        int bn[NBR_CAP], bfan[NBR_CAP], bdir[NBR_CAP]; double bang[NBR_CAP]; int nb = 0;
        for (int q = 0; q < nn; q++) {
            if (ncnt[q] != 1) { continue; }
            int32_t w = nbr[q];
            int owner = -1;
            for (size_t j = s; j < e; j++) {
                int fi = inc[j];
                if (faces[fi*3+0]==w || faces[fi*3+1]==w || faces[fi*3+2]==w) { owner = fi; break; }
            }
            double dz=vp[(size_t)w*3+0]-vp[v*3+0];
            double dy=vp[(size_t)w*3+1]-vp[v*3+1];
            double dx=vp[(size_t)w*3+2]-vp[v*3+2];
            double pu=dz*u[0]+dy*u[1]+dx*u[2], pw=dz*wv3[0]+dy*wv3[1]+dx*wv3[2];
            bn[nb]=w; bfan[nb]=nfan[q]; bang[nb]=atan2(pw,pu);
            bdir[nb]= owner>=0 ? face_edge_dir(faces, owner, (int32_t)v, w) : 0;
            nb++;
        }
        if (nb < 2) {
            if (dbg) fprintf(stderr, "  [bowtie] v%zu (%.2f %.2f %.2f) fans=%zu SKIP nb=%d<2\n",
                             v, (double)vp[v*3+0], (double)vp[v*3+1], (double)vp[v*3+2], fans, nb);
            continue;
        }
        for (int i = 1; i < nb; i++) {        /* insertion sort by angle */
            int kn=bn[i], kf=bfan[i], kd=bdir[i]; double ka=bang[i]; int jj=i-1;
            while (jj>=0 && bang[jj]>ka) { bn[jj+1]=bn[jj]; bfan[jj+1]=bfan[jj]; bdir[jj+1]=bdir[jj]; bang[jj+1]=bang[jj]; jj--; }
            bn[jj+1]=kn; bfan[jj+1]=kf; bdir[jj+1]=kd; bang[jj+1]=ka;
        }

        if (bridge_boundary_mode) {
            /* The bridge suffix is an open patch, not a complete surface.  A
             * two-fan bowtie has two angular gaps, but closing both consumes
             * the source-facing boundary and frequently makes a three-face
             * edge.  Build a shortest winding-compatible spanning forest of
             * fan gaps instead.  It joins all bridge fans while deliberately
             * leaving one boundary arc for the source chart attachment. */
            typedef struct {
                int a, b, fa, fb;
                double gap;
            } BridgeGap;
            BridgeGap gap[NBR_CAP];
            int fan_root[NBR_CAP], fan_parent[NBR_CAP];
            uint8_t fan_source[NBR_CAP],fan_bridge[NBR_CAP];
            int ng = 0, nfr = 0, fan_overflow = 0;
            int have_source=0,have_bridge=0;
            size_t fstart = nfill;

            for (size_t j = s; j < e; j++) {
                int fr = uf_find(uf, inc[j]), q = 0;
                while (q < nfr && fan_root[q] != fr) q++;
                if (q == nfr) {
                    if (nfr == NBR_CAP) { fan_overflow = 1; break; }
                    fan_root[nfr] = fr;
                    fan_parent[nfr] = nfr;
                    fan_source[nfr]=fan_bridge[nfr]=0;
                    nfr++;
                }
                if(bridge_boundary_mode==2){
                    int fi=inc[j];
                    if((size_t)fi<source_nf){
                        fan_source[q]=1;have_source=1;
                    }else{
                        fan_bridge[q]=1;have_bridge=1;
                    }
                }
            }
            if (fan_overflow) {
                if (dbg) fprintf(stderr,
                    "  [bridge-bowtie] v%zu SKIP fan-overflow\n", v);
                continue;
            }
            if(bridge_boundary_mode==2&&(!have_source||!have_bridge))
                continue;
            for (int i = 0; i < nb; i++) {
                int a = bn[i], b = bn[(i+1)%nb];
                double dz,dy,dx,glen;
                if (bfan[i] == bfan[(i+1)%nb]) continue;
                if (vertex_chart != NULL) {
                    int32_t cv[3]={
                        vertex_chart[v],vertex_chart[a],vertex_chart[b]};
                    int nc=0;
                    int allow_source_shoulder=0;
                    for(int k=0;k<3;k++){
                        int seen=0;
                        if(cv[k]<0||(size_t)cv[k]>=n_charts){nc=4;break;}
                        for(int q=0;q<k;q++)if(cv[q]==cv[k])seen=1;
                        if(!seen)nc++;
                    }
                    /* A same-strip bowtie closure has one original chart pair.
                     * Three distinct labels means the two bridge fans belong
                     * to different source-chart transactions; joining them
                     * here creates exactly the premature corner complex this
                     * chart-level pass is meant to avoid. */
                    /* In combined source+bridge mode, a missing triangle may
                     * be the final shoulder between a source-chart fan and a
                     * raw two-chart strip even though all three of the
                     * triangle's vertices belong to the source chart.  That
                     * is not an anonymous self-weld: require exactly one fan
                     * to contain source faces, require the other fan to carry
                     * one unambiguous two-chart pair containing this (small)
                     * chart, and tag the appended face by its creation range
                     * for ChartBridgeForest's initial support promotion.
                     * The forest will still admit or reject the complete
                     * patch transaction against live edge/link topology. */
                    if(nc==1&&bridge_boundary_mode==2&&
                       chart_face_count!=NULL&&
                       chart_face_count[cv[0]]<=
                           BRIDGE_MICRO_CHART_MAX_FACES){
                        int ia=0,ib=0,bridge_fr=-1,other=-1,ambiguous=0;
                        while(ia<nfr&&fan_root[ia]!=bfan[i])ia++;
                        while(ib<nfr&&fan_root[ib]!=bfan[(i+1)%nb])ib++;
                        if(ia<nfr&&ib<nfr&&
                           fan_source[ia]!=fan_source[ib]){
                            bridge_fr=fan_source[ia]?
                                fan_root[ib]:fan_root[ia];
                            for(size_t jj=s;jj<e;jj++){
                                int fi=inc[jj],fr=uf_find(uf,fi);
                                int32_t labels[3];int nlbl=0;
                                if(fr!=bridge_fr||(size_t)fi<source_nf)
                                    continue;
                                for(int k=0;k<3;k++){
                                    int32_t c=vertex_chart[
                                        faces[(size_t)fi*3+(size_t)k]];
                                    int seen=0;
                                    for(int q=0;q<nlbl;q++)
                                        if(labels[q]==c)seen=1;
                                    if(!seen&&nlbl<3)labels[nlbl++]=c;
                                }
                                if(nlbl==2){
                                    int32_t candidate=-1;
                                    if(labels[0]==cv[0])candidate=labels[1];
                                    else if(labels[1]==cv[0])candidate=labels[0];
                                    if(candidate>=0){
                                        if(other<0)other=candidate;
                                        else if(other!=candidate)ambiguous=1;
                                    }
                                }else if(nlbl>2){
                                    ambiguous=1;
                                }
                            }
                            allow_source_shoulder=other>=0&&!ambiguous;
                        }
                    }
                    if(nc!=2&&!allow_source_shoulder){
                        if(dbg)fprintf(stderr,
                            "  [bridge-bowtie] v%zu gap SKIP chart-corner "
                            "(%d distinct chart(s))\n",v,nc);
                        continue;
                    }else if(allow_source_shoulder&&dbg){
                        fprintf(stderr,
                            "  [bridge-bowtie] v%zu gap source shoulder "
                            "chart=%d\n",v,cv[0]);
                    }
                }
                dz=vp[(size_t)a*3+0]-vp[(size_t)b*3+0];
                dy=vp[(size_t)a*3+1]-vp[(size_t)b*3+1];
                dx=vp[(size_t)a*3+2]-vp[(size_t)b*3+2];
                glen=sqrt(dz*dz+dy*dy+dx*dx);
                if (glen > bt_gapmax || coincident_v(vp,a,b) ||
                    tri_min_altitude(vp,(int)v,a,b) < HOLEFILL_MIN_ALT_VOX) {
                    if (dbg) fprintf(stderr,
                        "  [bridge-bowtie] v%zu gap %.3f SKIP geometry\n",
                        v,glen);
                    continue;
                }
                gap[ng].a=a; gap[ng].b=b;
                gap[ng].fa=bfan[i]; gap[ng].fb=bfan[(i+1)%nb];
                gap[ng].gap=glen; ng++;
            }
            for (int i = 1; i < ng; i++) {
                BridgeGap key=gap[i]; int j=i-1;
                while (j>=0 && gap[j].gap>key.gap) {
                    gap[j+1]=gap[j]; j--;
                }
                gap[j+1]=key;
            }
            for (int i = 0; i < ng; i++) {
                int ia=0,ib=0,ra,rb,va,vb;
                double fn[3], dot;
                while (ia<nfr && fan_root[ia]!=gap[i].fa) ia++;
                while (ib<nfr && fan_root[ib]!=gap[i].fb) ib++;
                if (ia==nfr || ib==nfr) continue;
                ra=uf_find(fan_parent,ia); rb=uf_find(fan_parent,ib);
                if (ra==rb) continue;
                if (face_normal_d(vp,(int32_t)v,gap[i].a,gap[i].b,fn)
                    <= 1e-20) continue;
                dot=fn[0]*n[0]+fn[1]*n[1]+fn[2]*n[2];
                va=gap[i].a; vb=gap[i].b;
                if (dot<0.0) { int tmp=va; va=vb; vb=tmp; dot=-dot; }
                if (dot < coherence_gate ||
                    !bowtie_guard_triangle_ok(&guard,fill,fstart,nfill,
                                               (int32_t)v,va,vb,
                                               defer_source_orientation)) {
                    if (dbg) fprintf(stderr,
                        "  [bridge-bowtie] v%zu gap %.3f SKIP source-edge/winding\n",
                        v,gap[i].gap);
                    continue;
                }
                fill[nfill*3+0]=(int32_t)v;
                fill[nfill*3+1]=va;
                fill[nfill*3+2]=vb;
                nfill++;
                fan_parent[ra]=rb;
            }
            if (nfr > 0) {
                int root=uf_find(fan_parent,0), connected=1;
                for (int i=1;i<nfr;i++)
                    if (uf_find(fan_parent,i)!=root) { connected=0; break; }
                if (!connected) {
                    if (dbg) fprintf(stderr,
                        "  [bridge-bowtie] v%zu (%.2f %.2f %.2f) fans=%d "
                        "ABORT incomplete spanning closure\n",v,
                        (double)vp[v*3+0],(double)vp[v*3+1],
                        (double)vp[v*3+2],nfr);
                    nfill=fstart;
                    continue;
                }
            }
            if (nfill > fstart) {
                for (size_t f=fstart;f<nfill;f++)
                    for (int k=0;k<3;k++)
                        bowtie_guard_add_edge(
                            &guard,fill[f*3+(size_t)k],
                            fill[f*3+(size_t)((k+1)%3)],0);
                closed++;
                if (dbg) fprintf(stderr,
                    "  [bridge-bowtie] v%zu (%.2f %.2f %.2f) fans=%d "
                    "CLOSED %zu spanning gap(s)\n",v,
                    (double)vp[v*3+0],(double)vp[v*3+1],
                    (double)vp[v*3+2],nfr,nfill-fstart);
            }
            continue;
        }

        /* Join the fans with a shortest safe spanning tree of angular gaps.
         * Filling every cyclic gap is topologically wrong: two disk fans have
         * two gaps around v, and capping both creates an annulus.  One accepted
         * gap joins them into a disk while retaining the other as outer
         * boundary.  Kruskal also handles 3+ fans without creating a cycle. */
        typedef struct {
            int a, b, fa, fb, va, vb;
            double gap;
        } SourceGap;
        SourceGap source_gap[NBR_CAP];
        int fan_root[NBR_CAP], fan_parent[NBR_CAP];
        int nsource_gap = 0, nfr = 0, fan_overflow = 0;
        size_t fstart = nfill;

        for (size_t j = s; j < e; j++) {
            int fr = uf_find(uf, inc[j]), q = 0;
            while (q < nfr && fan_root[q] != fr) q++;
            if (q == nfr) {
                if (nfr == NBR_CAP) { fan_overflow = 1; break; }
                fan_root[nfr] = fr;
                fan_parent[nfr] = nfr;
                nfr++;
            }
        }
        if (fan_overflow) {
            if (dbg) fprintf(stderr,
                "  [bowtie] v%zu fans=%zu SKIP fan-overflow\n", v, fans);
            continue;
        }
        for (int i = 0; i < nb; i++) {
            int a = bn[i], b = bn[(i+1)%nb];
            int da = bdir[i], db = bdir[(i+1)%nb];
            if (bfan[i] == bfan[(i+1)%nb]) { continue; }
            double dz=vp[(size_t)a*3+0]-vp[(size_t)b*3+0];
            double dy=vp[(size_t)a*3+1]-vp[(size_t)b*3+1];
            double dx=vp[(size_t)a*3+2]-vp[(size_t)b*3+2];
            double glen = sqrt(dz*dz+dy*dy+dx*dx);
            const char *reject = NULL;
            if (glen > bt_gapmax) reject = "wide-gap";
            else if (coincident_v(vp, a, b)) reject = "coincident";
            else if (tri_min_altitude(vp, (int)v, a, b) <
                     HOLEFILL_MIN_ALT_VOX) reject = "sliver";
            int va, vb;
            if (reject == NULL && glen > (double)BOWTIE_GAP_MAX_VOX) {
                /* WIDE (armed) gap: orientation gate. The fan must reverse both
                 * fans' boundary half-edges (v,a) and (v,b). That is possible iff
                 * the two owning faces traverse them in OPPOSITE senses at v
                 * (da != db) -- a benign same-sheet pinch. Equal senses (or an
                 * unknown dir) is a FOLD (fans wound oppositely); abort to split.
                 * The winding is then fixed by the reversal, NOT by the normal. */
                if (da == 0 || db == 0 || da == db) {
                    reject = "fold";
                    va = a; vb = b;
                } else if (da == -1) {
                    va = a; vb = b;       /* (v,a,b): v->a, b->v */
                } else {
                    va = b; vb = a;       /* (v,b,a): a->v, v->b */
                }
            } else if (reject == NULL) {
                /* NARROW gap: unchanged -- wind (v,a,b) so its normal agrees with n. */
                double e1[3]={vp[(size_t)a*3+0]-vp[v*3+0], vp[(size_t)a*3+1]-vp[v*3+1], vp[(size_t)a*3+2]-vp[v*3+2]};
                double e2[3]={vp[(size_t)b*3+0]-vp[v*3+0], vp[(size_t)b*3+1]-vp[v*3+1], vp[(size_t)b*3+2]-vp[v*3+2]};
                double cr[3]={e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]};
                va = a; vb = b;
                if (cr[0]*n[0]+cr[1]*n[1]+cr[2]*n[2] < 0) { va = b; vb = a; }
            } else {
                va = a; vb = b;
            }
            if (reject != NULL) {
                if (dbg) fprintf(stderr,
                    "  [bowtie] v%zu candidate gap %.3f SKIP %s\n",
                    v, glen, reject);
                continue;
            }
            if (nsource_gap == NBR_CAP) { fan_overflow = 1; break; }
            source_gap[nsource_gap].a = a;
            source_gap[nsource_gap].b = b;
            source_gap[nsource_gap].fa = bfan[i];
            source_gap[nsource_gap].fb = bfan[(i+1)%nb];
            source_gap[nsource_gap].va = va;
            source_gap[nsource_gap].vb = vb;
            source_gap[nsource_gap].gap = glen;
            nsource_gap++;
        }
        if (fan_overflow) {
            if (dbg) fprintf(stderr,
                "  [bowtie] v%zu fans=%zu SKIP candidate-overflow\n", v, fans);
            continue;
        }
        for (int i = 1; i < nsource_gap; i++) {
            SourceGap key = source_gap[i];
            int j = i - 1;
            while (j >= 0 && source_gap[j].gap > key.gap) {
                source_gap[j + 1] = source_gap[j];
                j--;
            }
            source_gap[j + 1] = key;
        }
        for (int i = 0; i < nsource_gap; i++) {
            int ia = 0, ib = 0;
            while (ia < nfr && fan_root[ia] != source_gap[i].fa) ia++;
            while (ib < nfr && fan_root[ib] != source_gap[i].fb) ib++;
            if (ia == nfr || ib == nfr) continue;
            int ra = uf_find(fan_parent, ia);
            int rb = uf_find(fan_parent, ib);
            if (ra == rb) continue;
            fill[nfill*3+0] = (int32_t)v;
            fill[nfill*3+1] = source_gap[i].va;
            fill[nfill*3+2] = source_gap[i].vb;
            nfill++;
            fan_parent[ra] = rb;
        }
        int connected = nfr > 0;
        if (connected) {
            int root = uf_find(fan_parent, 0);
            for (int i = 1; i < nfr; i++)
                if (uf_find(fan_parent, i) != root) { connected = 0; break; }
        }
        if (!connected) {
            nfill = fstart;
            if (dbg) fprintf(stderr,
                "  [bowtie] v%zu (%.2f %.2f %.2f) fans=%zu ABORT "
                "no safe spanning closure\n",
                v, (double)vp[v*3+0], (double)vp[v*3+1],
                (double)vp[v*3+2], fans);
            continue;
        }
        if (nfill > fstart) {
            closed++;
            if (dbg) fprintf(stderr, "  [bowtie] v%zu (%.2f %.2f %.2f) fans=%zu CLOSED %zu spanning gap(s)\n",
                             v, (double)vp[v*3+0], (double)vp[v*3+1], (double)vp[v*3+2], fans, nfill-fstart);
        } else if (dbg) {
            fprintf(stderr, "  [bowtie] v%zu (%.2f %.2f %.2f) fans=%zu nb=%d candidates=%d NO-FILL\n",
                    v, (double)vp[v*3+0], (double)vp[v*3+1], (double)vp[v*3+2], fans, nb, nsource_gap);
        }
    }

    if (nfill > 0) {
        size_t nf_new = nf + nfill;
        int32_t *nff = (int32_t *)ARENA_ALLOC(
            arena, nf_new * 3 * sizeof(int32_t));
        memcpy(nff, faces, nf * 3 * sizeof(int32_t));
        memcpy(nff + nf * 3, fill, nfill * 3 * sizeof(int32_t));
        cm->faces = nff; cm->nf = nf_new;
    }
    if (closed && getenv("PINHOLE_DEBUG")) {
        fprintf(stderr, "  [pinhole] close_bowtie_gaps: closed %zu bowtie(s) "
                "in place, +%zu fill tris (no vertex split)\n", closed, nfill);
    }
    free(guard.slot);
    if (out_closed) { *out_closed = closed; }
    return 0;
}

int PinholeFill_close_bowties(Arena_T arena,
                              ComponentMesh *meshes, size_t n_meshes,
                              size_t *out_closed)
{
    size_t total = 0;
    for (size_t i = 0; i < n_meshes; i++) {
        size_t closed = 0;
        if (meshes[i].nf == 0) continue;
        if (close_bowtie_gaps(arena, &meshes[i], NULL, 0, NULL, 0, NULL, 0, 0,
                              &closed) != 0) return -1;
        total += closed;
    }
    if (out_closed) *out_closed = total;
    return 0;
}

int PinholeFill_close_bridge_bowties(Arena_T arena,
                                     ComponentMesh *bridge,
                                     const int32_t *source_faces,
                                     size_t source_nf,
                                     const int32_t *vertex_chart,
                                     size_t n_charts,
                                     int defer_source_orientation,
                                     size_t *out_closed)
{
    if (bridge == NULL || (source_nf > 0 && source_faces == NULL) ||
        vertex_chart == NULL || n_charts == 0) return -1;
    if (bridge->nf == 0) {
        if (out_closed) *out_closed = 0;
        return 0;
    }
    return close_bowtie_gaps(arena, bridge, source_faces, source_nf,
                             vertex_chart, n_charts, NULL, 1,
                             defer_source_orientation,
                             out_closed);
}

int PinholeFill_close_bridge_attachments(Arena_T arena,
                                         ComponentMesh *combined,
                                         size_t source_nf,
                                         const int32_t *vertex_chart,
                                         size_t n_charts,
                                         const size_t *chart_face_count,
                                         int defer_source_orientation,
                                         size_t *out_closed)
{
    if(combined==NULL||source_nf>combined->nf||vertex_chart==NULL||
       n_charts==0)return -1;
    if(combined->nf==source_nf){
        if(out_closed)*out_closed=0;
        return 0;
    }
    /* Mode 2 inspects vertex fans in the combined source+bridge mesh.  The
     * source_nf prefix marks fan provenance; the edge guard is built from the
     * combined mesh exactly once. */
    return close_bowtie_gaps(arena,combined,NULL,source_nf,
                             vertex_chart,n_charts,chart_face_count,2,
                             defer_source_orientation,
                             out_closed);
}

/* PHASE 1 — split bowtie / pinch vertices. */
static size_t split_pinch_verts(Arena_T arena, ComponentMesh *cm,
                                size_t *out_splits,
                                int32_t **out_vertex_source) {
    size_t nv = cm->nv, nf = cm->nf;
    int32_t *faces = cm->faces;

    if (out_vertex_source) *out_vertex_source = NULL;

    size_t *off = NULL; int *inc = NULL;
    build_vert_faces(arena, faces, nf, nv, &off, &inc);

    int *uf = (int *)ARENA_ALLOC(arena, nf * sizeof(int));
    /* root_of_slot[j] = fan root (a face id) for incident slot j */
    int *root_of_slot = (int *)ARENA_ALLOC(arena, off[nv] * sizeof(int));

    /* sub-pass 1: fan roots per slot + count extra verts */
    size_t extra = 0;
    for (size_t v = 0; v < nv; v++) {
        size_t s = off[v], e = off[v + 1];
        for (size_t j = s; j < e; j++) { uf[inc[j]] = inc[j]; }
        for (size_t j = s; j < e; j++) {
            int fi = inc[j];
            int o1 = -1, o2 = -1;
            for (int k = 0; k < 3; k++) {
                int32_t w = faces[fi * 3 + k];
                if ((size_t)w != v) { if (o1 < 0) { o1 = w; } else { o2 = w; } }
            }
            for (size_t j2 = j + 1; j2 < e; j2++) {
                int gj = inc[j2];
                int p1 = -1, p2 = -1;
                for (int k = 0; k < 3; k++) {
                    int32_t w = faces[gj * 3 + k];
                    if ((size_t)w != v) {
                        if (p1 < 0) { p1 = w; } else { p2 = w; }
                    }
                }
                if (o1 == p1 || o1 == p2 || o2 == p1 || o2 == p2) {
                    int ra = uf_find(uf, fi), rb = uf_find(uf, gj);
                    if (ra != rb) { uf[ra] = rb; }
                }
            }
        }
        size_t fans = 0;
        for (size_t j = s; j < e; j++) {
            int r = uf_find(uf, inc[j]);
            root_of_slot[j] = r;
            if (r == inc[j]) { fans++; }
        }
        if (fans >= 2) { extra += fans - 1; }
    }

    if (extra == 0) {
        if (out_vertex_source) {
            int32_t *source = (int32_t *)ARENA_ALLOC(
                arena, (nv ? nv : 1) * sizeof(*source));
            for (size_t v = 0; v < nv; v++) source[v] = (int32_t)v;
            *out_vertex_source = source;
        }
        if (out_splits) { *out_splits = 0; }
        return nv;
    }

    /* sub-pass 2: assign output ids per fan root; fill newid[slot] */
    int *newid = (int *)ARENA_ALLOC(arena, off[nv] * sizeof(int));
    int *src_of_new = (int *)ARENA_ALLOC(arena, extra * sizeof(int));
    size_t next_new = nv, nsplit = 0;
    for (size_t v = 0; v < nv; v++) {
        size_t s = off[v], e = off[v + 1];
        /* map of (root -> output id) for this vertex; small, linear scan */
        int seen_root[64]; int seen_id[64]; int nseen = 0;
        for (size_t j = s; j < e; j++) {
            int r = root_of_slot[j];
            int oid = -1;
            for (int t = 0; t < nseen; t++) {
                if (seen_root[t] == r) { oid = seen_id[t]; break; }
            }
            if (oid < 0) {
                /* first fan, or >64-fan overflow (impossible in practice):
                 * keep v so next_new can never exceed the `extra` reserve. */
                if (nseen == 0 || nseen >= 64) {
                    oid = (int)v;
                } else {
                    oid = (int)next_new++;         /* extra fan -> new id */
                    src_of_new[oid - (int)nv] = (int)v;
                    nsplit++;
                }
                if (nseen < 64) {
                    seen_root[nseen] = r; seen_id[nseen] = oid; nseen++;
                }
            }
            newid[j] = oid;
        }
    }

    /* grow vertex-keyed arrays */
    size_t nv_new = nv + extra;
    float *nv_verts = (float *)ARENA_ALLOC(
        arena, nv_new * 3 * sizeof(float));
    memcpy(nv_verts, cm->verts, nv * 3 * sizeof(float));
    float *nv_norm = NULL;
    if (cm->vert_normals) {
        nv_norm = (float *)ARENA_ALLOC(
            arena, nv_new * 3 * sizeof(float));
        memcpy(nv_norm, cm->vert_normals, nv * 3 * sizeof(float));
    }
    uint8_t *nv_pin = NULL;
    if (cm->pin_mask) {
        nv_pin = (uint8_t *)ARENA_ALLOC(arena, nv_new * sizeof(uint8_t));
        memcpy(nv_pin, cm->pin_mask, nv * sizeof(uint8_t));
    }
    for (size_t i = 0; i < extra; i++) {
        int src = src_of_new[i];
        size_t dst = nv + i;
        nv_verts[dst * 3 + 0] = cm->verts[(size_t)src * 3 + 0];
        nv_verts[dst * 3 + 1] = cm->verts[(size_t)src * 3 + 1];
        nv_verts[dst * 3 + 2] = cm->verts[(size_t)src * 3 + 2];
        if (nv_norm) {
            nv_norm[dst * 3 + 0] = cm->vert_normals[(size_t)src * 3 + 0];
            nv_norm[dst * 3 + 1] = cm->vert_normals[(size_t)src * 3 + 1];
            nv_norm[dst * 3 + 2] = cm->vert_normals[(size_t)src * 3 + 2];
        }
        if (nv_pin) { nv_pin[dst] = cm->pin_mask[(size_t)src]; }
    }

    /* rewrite faces: repoint each split corner to its fan's output id */
    for (size_t v = 0; v < nv; v++) {
        for (size_t j = off[v]; j < off[v + 1]; j++) {
            if ((size_t)newid[j] == v) { continue; }
            int fi = inc[j];
            for (int k = 0; k < 3; k++) {
                if ((size_t)faces[fi * 3 + k] == v) {
                    faces[fi * 3 + k] = newid[j];
                    break;
                }
            }
        }
    }

    cm->verts = nv_verts;
    cm->vert_normals = nv_norm ? nv_norm : cm->vert_normals;
    cm->pin_mask = nv_pin ? nv_pin : cm->pin_mask;
    cm->nv = nv_new;
    if (out_vertex_source) {
        int32_t *source = (int32_t *)ARENA_ALLOC(
            arena, nv_new * sizeof(*source));
        for (size_t v = 0; v < nv; v++) source[v] = (int32_t)v;
        for (size_t i = 0; i < extra; i++) source[nv+i] = src_of_new[i];
        *out_vertex_source = source;
    }
    if (out_splits) { *out_splits = nsplit; }
    return nv_new;
}

/* Public phase-1-only entry: split pinch vertices across an array of meshes.
 * Used by the post-trim manifold guard. (split_pinch_verts grows cm->verts from
 * the arena and repoints faces in place; counts accumulate.) */
int PinholeFill_split_pinches(Arena_T arena,
                              ComponentMesh *meshes, size_t n_meshes,
                              size_t *out_splits) {
    size_t total = 0;
    for (size_t i = 0; i < n_meshes; i++) {
        ComponentMesh *cm = &meshes[i];
        if (cm->nf == 0) { continue; }
        size_t sp = 0;
        split_pinch_verts(arena, cm, &sp, NULL);
        total += sp;
    }
    if (out_splits) { *out_splits = total; }
    return 0;
}

int PinholeFill_split_pinches_mapped(Arena_T arena,
                                     ComponentMesh *mesh,
                                     size_t *out_splits,
                                     int32_t **out_vertex_source)
{
    if (out_splits) *out_splits = 0;
    if (out_vertex_source) *out_vertex_source = NULL;
    if (!arena || !mesh || !mesh->verts || !mesh->faces ||
        mesh->nv == 0 || mesh->nf == 0 || !out_vertex_source)
        return -1;
    split_pinch_verts(arena, mesh, out_splits, out_vertex_source);
    return *out_vertex_source ? 0 : -1;
}

/* PHASE 2 — fill small boundary loops. */
static void fill_small_loops(Arena_T arena, ComponentMesh *cm,
                             int respect_pins,
                             size_t *out_filled, size_t *out_added,
                             size_t *out_skipped) {
    size_t nv = cm->nv, nf = cm->nf;
    int32_t *faces = cm->faces;
    const uint8_t *pin = (respect_pins && cm->pin_mask) ? cm->pin_mask : NULL;

    /* undirected edge -> face count, plus first incident face (for component
     * union-find below) */
    size_t hsz = next_pow2(nf * 3 * 2);
    if (hsz < 1024) { hsz = 1024; }
    uint64_t hmask = (uint64_t)(hsz - 1);
    /* Keep allocation arithmetic in size_t.  On Win64, long is still 32-bit;
     * the 4x21x21 source has hsz=2^29, so the old cast wrapped the exact
     * 2-GiB bucket table to LONG_MIN and then sign-extended to SIZE_MAX-2GiB. */
    int *bk = (int *)ARENA_ALLOC(arena, hsz * sizeof(int));
    memset(bk, 0xFF, hsz * sizeof(int));
    E *ed = (E *)ARENA_ALLOC(arena, nf * 3 * sizeof(E));
    size_t en = 0;

    /* face -> connected component (union across every shared edge). Lets us
     * tell a genuine interior pinhole (whose component ALSO owns the big outer
     * rim, so the component has >=2 boundary loops) from the sole rim of a tiny
     * island (its component's ONLY loop) — the latter must never be capped into
     * a closed bubble. */
    int *cuf = (int *)ARENA_ALLOC(arena, nf * sizeof(int));
    for (size_t f = 0; f < nf; f++) { cuf[f] = (int)f; }

    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f * 3 + e], b = faces[f * 3 + (e + 1) % 3];
            int32_t lo = (a < b) ? a : b, hi = (a < b) ? b : a;
            size_t bb = (size_t)(edge_hash(lo, hi) & hmask);
            int idx = bk[bb];
            while (idx >= 0) {
                if (ed[idx].lo == lo && ed[idx].hi == hi) { break; }
                idx = ed[idx].next;
            }
            if (idx < 0) {
                idx = (int)en++;
                ed[idx].lo = lo; ed[idx].hi = hi; ed[idx].cnt = 1;
                ed[idx].next = bk[bb]; bk[bb] = idx; ed[idx].f0 = (int)f;
            } else {
                ed[idx].cnt++;
                int ra = uf_find(cuf, ed[idx].f0), rb = uf_find(cuf, (int)f);
                if (ra != rb) { cuf[ra] = rb; }
            }
        }
    }

    /* directed boundary half-edges (a->b in their single face) */
    size_t nb = 0;
    for (size_t i = 0; i < en; i++) { if (ed[i].cnt == 1) { nb++; } }
    if (nb == 0) {
        if (out_filled) { *out_filled = 0; }
        if (out_added) { *out_added = 0; }
        if (out_skipped) { *out_skipped = 0; }
        return;
    }
    int *he_src = (int *)ARENA_ALLOC(arena, nb * sizeof(int));
    int *he_dst = (int *)ARENA_ALLOC(arena, nb * sizeof(int));
    int *he_face = (int *)ARENA_ALLOC(arena, nb * sizeof(int));
    size_t hn = 0;
    for (size_t f = 0; f < nf; f++) {
        for (int e = 0; e < 3; e++) {
            int32_t a = faces[f * 3 + e], b = faces[f * 3 + (e + 1) % 3];
            int32_t lo = (a < b) ? a : b, hi = (a < b) ? b : a;
            size_t bb = (size_t)(edge_hash(lo, hi) & hmask);
            int idx = bk[bb];
            while (idx >= 0) {
                if (ed[idx].lo == lo && ed[idx].hi == hi) { break; }
                idx = ed[idx].next;
            }
            if (idx >= 0 && ed[idx].cnt == 1) {
                he_src[hn] = a; he_dst[hn] = b; he_face[hn] = (int)f; hn++;
            }
        }
    }

    /* per-vertex outgoing boundary half-edges (CSR) */
    size_t *voff = (size_t *)ARENA_CALLOC(arena, nv + 1, sizeof(size_t));
    for (size_t i = 0; i < nb; i++) { voff[(size_t)he_src[i] + 1]++; }
    for (size_t v = 0; v < nv; v++) { voff[v + 1] += voff[v]; }
    int *vhe = (int *)ARENA_ALLOC(arena, nb * sizeof(int));
    size_t *vcur = (size_t *)ARENA_ALLOC(arena, nv * sizeof(size_t));
    memcpy(vcur, voff, nv * sizeof(size_t));
    for (size_t i = 0; i < nb; i++) { vhe[vcur[(size_t)he_src[i]]++] = (int)i; }

    uint8_t *used = (uint8_t *)ARENA_CALLOC(arena, nb, 1);
    int loop[PINHOLE_MAX_LOOP + 2];

    /* fill triangles buffer: <= nb triangles added (each boundary edge
     * becomes interior, contributing at most one new fan triangle). */
    int32_t *fill = (int32_t *)ARENA_ALLOC(
        arena, nb * 3 * sizeof(int32_t));
    size_t nfill = 0, filled = 0, skipped = 0;

    /* triangle-existence set over the existing faces, kept live as fills are
     * committed. A fan triangle that duplicates an existing/just-filled tri
     * (re-covering a flap's back side) would create a doubled-triangle pocket
     * -> a 4-fan non-manifold edge; reject the whole loop instead. */
    size_t tsz = next_pow2((nf + nb) * 2);
    if (tsz < 1024) { tsz = 1024; }
    uint64_t tmask = (uint64_t)(tsz - 1);
    uint64_t *tset = (uint64_t *)ARENA_CALLOC(arena, tsz, sizeof(uint64_t));
    for (size_t f = 0; f < nf; f++) {
        tri_add(tset, tmask, faces[f * 3 + 0], faces[f * 3 + 1],
                faces[f * 3 + 2]);
    }

    /* First pass: tally boundary loops per face-component. A real interior
     * pinhole's component also carries the large outer rim (>=2 loops); a tiny
     * island's component has only its own rim (1 loop) and must not be capped
     * into a closed bubble. Walk consumes the same half-edges the fill pass
     * will, so each distinct boundary cycle/chain is counted once. */
    size_t *loops_per_comp = (size_t *)ARENA_CALLOC(
        arena, nf, sizeof(size_t));
    for (size_t i0 = 0; i0 < nb; i0++) {
        if (used[i0]) { continue; }
        loops_per_comp[(size_t)uf_find(cuf, he_face[i0])]++;
        int start = he_src[i0], cur = he_dst[i0];
        used[i0] = 1;
        size_t steps = 0;
        while (cur != start && steps <= nb) {
            int nxt = -1;
            for (size_t k = voff[(size_t)cur]; k < voff[(size_t)cur + 1]; k++) {
                if (!used[vhe[k]]) { nxt = vhe[k]; break; }
            }
            if (nxt < 0) { break; }
            used[nxt] = 1; cur = he_dst[nxt]; steps++;
        }
    }
    memset(used, 0, nb);

    for (size_t i0 = 0; i0 < nb; i0++) {
        if (used[i0]) { continue; }
        int start = he_src[i0];
        int len = 0;
        loop[len++] = start;
        int cur = he_dst[i0];
        used[i0] = 1;
        int closed = 0, bad = 0;
        while (1) {
            if (cur == start) { closed = 1; break; }
            if (len > PINHOLE_MAX_LOOP) { bad = 1; break; } /* too big */
            /* duplicate-vertex guard (figure-8 / unsplit pinch) */
            for (int t = 0; t < len; t++) {
                if (loop[t] == cur) { bad = 1; break; }
            }
            if (bad) { break; }
            loop[len++] = cur;
            int nxt = -1;
            for (size_t k = voff[(size_t)cur]; k < voff[(size_t)cur + 1]; k++) {
                if (!used[vhe[k]]) { nxt = vhe[k]; break; }
            }
            if (nxt < 0) { bad = 1; break; }   /* open chain */
            used[nxt] = 1;
            cur = he_dst[nxt];
        }
        if (!closed || bad) { continue; }      /* leave for HoleFill */
        if (len < 3 || len > PINHOLE_MAX_LOOP) { continue; }

        if (pin) {
            int touches = 0;
            for (int t = 0; t < len; t++) {
                if (pin[loop[t]]) { touches = 1; break; }
            }
            if (touches) { skipped++; continue; }
        }

        /* No-merger diameter gate: a genuine pinhole is a few missing triangles
         * spanning a couple of voxels; a wider loop is a real opening (possibly
         * the gap between two wraps) and MUST stay open — filling it could weld
         * two sheets. Refuse anything that does not fit in the diameter ball.
         * A 3-loop is exempt to the wider PINHOLE_TRI_DIAM_VOX cap: its three
         * edges already exist, so capping it adds no new edge and cannot merge
         * two wraps at any diameter (see PINHOLE_TRI_DIAM_VOX). */
        {
            const float *vp = cm->verts;
            float diam2 = 0.0f;
            for (int a = 0; a < len; a++) {
                for (int b = a + 1; b < len; b++) {
                    float dz = vp[loop[a] * 3 + 0] - vp[loop[b] * 3 + 0];
                    float dy = vp[loop[a] * 3 + 1] - vp[loop[b] * 3 + 1];
                    float dx = vp[loop[a] * 3 + 2] - vp[loop[b] * 3 + 2];
                    float d2 = dz * dz + dy * dy + dx * dx;
                    if (d2 > diam2) { diam2 = d2; }
                }
            }
            float cap = (len == 3) ? PINHOLE_TRI_DIAM_VOX : PINHOLE_MAX_DIAM_VOX;
            if (diam2 > cap * cap) {
                skipped++;
#ifdef PINHOLE_DEBUG
                fprintf(stderr, "  [pinhole] skip DIAMETER len=%d diam=%.2f "
                        "cap=%.2f v0=%d\n", len, (double)sqrtf(diam2),
                        (double)cap, loop[0]);
#endif
                continue;
            }
        }

        /* Component gate: never cap the SOLE boundary loop of a component — an
         * isolated island would become a closed bubble. Genuine pinholes live
         * on components that also own the big outer rim (>=2 boundary loops). */
        if (loops_per_comp[(size_t)uf_find(cuf, he_face[i0])] < 2) {
            skipped++;
#ifdef PINHOLE_DEBUG
            fprintf(stderr, "  [pinhole] skip COMPONENT len=%d v0=%d "
                    "loops_in_comp=%zu\n", len, loop[0],
                    loops_per_comp[(size_t)uf_find(cuf, he_face[i0])]);
#endif
            continue;
        }

        /* PinholeFill closes ONLY exact single-triangle (3-loop) holes — a
         * larger loop is a real hole left to HoleFill's CDT path. The 3-loop's
         * unique fill is the triangle on its three verts, wound (loop0,loop2,
         * loop1) = the boundary-edge reversal (consistent with the neighbours,
         * the same winding the old fan used). Skip it if it duplicates a face,
         * is zero-area from a coincident (phase-1 split) pair, or is a
         * degenerate collinear sliver — the "chewed" boundary along a cube face
         * is a row of near-collinear verts whose fill would overlap its own edge
         * and corrupt the weld; that notch is left open for the seam weld. */
        if (len != 3) { continue; }                 /* 4+ -> HoleFill */
        {
            int A = loop[0], B = loop[1], C = loop[2];
            const float *vp = cm->verts;
            if (tri_has(tset, tmask, A, B, C)) { skipped++; continue; }
            if (coincident_v(vp, A, B) || coincident_v(vp, B, C) ||
                coincident_v(vp, A, C)) { skipped++; continue; }
            if (tri_min_altitude(vp, A, B, C) < HOLEFILL_MIN_ALT_VOX) {
                skipped++;
#ifdef PINHOLE_DEBUG
                fprintf(stderr, "  [pinhole] skip DEGENERATE 3-loop "
                        "alt=%.4f v0=%d\n",
                        (double)tri_min_altitude(vp, A, B, C), A);
#endif
                continue;
            }
            fill[nfill * 3 + 0] = A;
            fill[nfill * 3 + 1] = C;
            fill[nfill * 3 + 2] = B;
            nfill++;
            tri_add(tset, tmask, A, B, C);
            filled++;
        }
    }

    if (nfill > 0) {
        size_t nf_new = nf + nfill;
        int32_t *nf_faces = (int32_t *)ARENA_ALLOC(
            arena, nf_new * 3 * sizeof(int32_t));
        memcpy(nf_faces, faces, nf * 3 * sizeof(int32_t));
        memcpy(nf_faces + nf * 3, fill, nfill * 3 * sizeof(int32_t));
        cm->faces = nf_faces;
        cm->nf = nf_new;
    }
    if (out_filled) { *out_filled = filled; }
    if (out_added) { *out_added = nfill; }
    if (out_skipped) { *out_skipped = skipped; }
}

int PinholeFill_fill_small_loops(Arena_T arena,
                                 ComponentMesh *meshes, size_t n_meshes,
                                 int respect_pins,
                                 size_t *out_loops_filled,
                                 size_t *out_tris_added,
                                 size_t *out_loops_skipped)
{
    size_t total_filled = 0, total_added = 0, total_skipped = 0;
    for (size_t i = 0; i < n_meshes; i++) {
        size_t filled = 0, added = 0, skipped = 0;
        if (meshes[i].nf == 0 || meshes[i].nv < 3) continue;
        fill_small_loops(arena, &meshes[i], respect_pins,
                         &filled, &added, &skipped);
        total_filled += filled;
        total_added += added;
        total_skipped += skipped;
    }
    if (out_loops_filled) *out_loops_filled = total_filled;
    if (out_tris_added) *out_tris_added = total_added;
    if (out_loops_skipped) *out_loops_skipped = total_skipped;
    return 0;
}

int PinholeFill_process(Arena_T arena,
                        ComponentMesh *meshes, size_t n_meshes,
                        int respect_pins,
                        size_t *out_pinch_splits,
                        size_t *out_loops_filled,
                        size_t *out_tris_added,
                        size_t *out_loops_skipped) {
    size_t tot_splits = 0, tot_filled = 0, tot_added = 0, tot_skipped = 0;
    for (size_t i = 0; i < n_meshes; i++) {
        ComponentMesh *cm = &meshes[i];
        if (cm->nf < 1 || cm->nv < 3) { continue; }
        Arena_Mark mark = Arena_save(arena);
        size_t splits = 0;
        /* Phase 0: close coherent narrow bowties in place (keep the vertex), so
         * phase 1 only splits the genuinely divergent/wide ones. */
        close_bowtie_gaps(arena, cm, NULL, 0, NULL, 0, NULL, 0, 0, NULL);
        split_pinch_verts(arena, cm, &splits, NULL);
        size_t filled = 0, added = 0, skipped = 0;
        fill_small_loops(arena, cm, respect_pins, &filled, &added, &skipped);
        /* cm->verts/faces/etc now point into arena past `mark`; do NOT
         * restore — the new buffers must outlive this call. */
        (void)mark;
        tot_splits += splits; tot_filled += filled;
        tot_added += added; tot_skipped += skipped;
    }
    if (out_pinch_splits) { *out_pinch_splits = tot_splits; }
    if (out_loops_filled) { *out_loops_filled = tot_filled; }
    if (out_tris_added) { *out_tris_added = tot_added; }
    if (out_loops_skipped) { *out_loops_skipped = tot_skipped; }
    return 0;
}

int HoleFill_meshes(Arena_T arena,
                    ComponentMesh *meshes, size_t n_meshes,
                    int respect_pins,
                    HoleFill_cdt_fn cdt_fill, void *cdt_user,
                    size_t *out_pinch_splits,
                    size_t *out_loops_filled,
                    size_t *out_tris_added,
                    size_t *out_loops_skipped,
                    size_t *out_cdt_filled) {
    /* Phase 1+2: pinch-split + EXACT single-triangle (3-loop) fills. This is the
     * CDT-free core, also used standalone by grid_weld. */
    PinholeFill_process(arena, meshes, n_meshes, respect_pins,
                        out_pinch_splits, out_loops_filled,
                        out_tris_added, out_loops_skipped);
    /* Phase 3: hand each mesh's remaining 4+ closed loops to the injected
     * CDT/Liepa filler. NULL (grid_weld) => 4+ loops stay open for the seam
     * weld. The filler owns its own loop tracing and gates (size, no-merger);
     * the 3-loops are already resolved above, so it must skip len==3. */
    size_t cdt_filled = 0;
    if (cdt_fill) {
        for (size_t i = 0; i < n_meshes; i++) {
            ComponentMesh *cm = &meshes[i];
            if (cm->nf < 1 || cm->nv < 3) { continue; }
            if (cdt_fill(arena, cm, cdt_user) == 0) { cdt_filled++; }
        }
    }
    if (out_cdt_filled) { *out_cdt_filled = cdt_filled; }
    return 0;
}
