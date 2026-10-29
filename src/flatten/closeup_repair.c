#include "closeup_repair.h"
#include "../common/arena.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int cr_observed(uint8_t s) { return s == 1 || s == 3; }
static int cr_fill(uint8_t s) { return s == 2; }

typedef struct CrCell { uint64_t key; int32_t v; } CrCell;

static int cr_cell_order(const void *x, const void *y)
{
    const CrCell *a = x, *b = y;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return a->v < b->v ? -1 : a->v > b->v;
}

static uint64_t cr_key(int64_t cz, int64_t cy, int64_t cx)
{
    const int64_t off = (int64_t)1 << 20;
    return (((uint64_t)(cz + off) & 0x1FFFFFu) << 42) | (((uint64_t)(cy + off) & 0x1FFFFFu) << 21) |
           ((uint64_t)(cx + off) & 0x1FFFFFu);
}

static size_t cr_lower(const CrCell *cells, size_t n, uint64_t key)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (cells[mid].key < key) lo = mid + 1; else hi = mid;
    }
    return lo;
}

static int32_t cr_find(int32_t *par, int32_t a)
{
    while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; }
    return a;
}

static void cr_union(int32_t *par, int32_t a, int32_t b)
{
    a = cr_find(par, a); b = cr_find(par, b);
    if (a != b) { if (a < b) par[b] = a; else par[a] = b; }
}

/* an observed vertex facing away from the axis */
static int cr_away(const uint8_t *state, const uint8_t *defined, const float *facing, size_t v)
{
    return cr_observed(state[v]) && defined[v] && facing[v] < 0.0f;
}

static int64_t cr_cell(float p) { return (int64_t)floor((double)p / CLOSEUP_REPAIR_DUP_VOX); }

/* The first kept vertex w (not v) within DUP_VOX of v in 3-D and more than DUP_PX away on the surface that `accept`
 * admits; -1 when none. */
typedef int (*CrAccept)(size_t v, size_t w, const void *ctx);

static long cr_duplicate_of(const CrCell *cells, size_t nc, const float *verts, const double *uv, const uint8_t *reason,
                            size_t v, CrAccept accept, const void *ctx)
{
    const double r2 = CLOSEUP_REPAIR_DUP_VOX * CLOSEUP_REPAIR_DUP_VOX, f2 = CLOSEUP_REPAIR_DUP_PX * CLOSEUP_REPAIR_DUP_PX;
    const float *p = verts + 3*v;
    int64_t c0 = cr_cell(p[0]), c1 = cr_cell(p[1]), c2 = cr_cell(p[2]);
    for (int dz = -1; dz <= 1; dz++) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
        uint64_t key = cr_key(c0 + dz, c1 + dy, c2 + dx);
        for (size_t i = cr_lower(cells, nc, key); i < nc && cells[i].key == key; i++) {
            size_t w = (size_t)cells[i].v;
            if (w == v || reason[w] != CLOSEUP_REPAIR_KEEP) continue;
            const float *q = verts + 3*w;
            double a = (double)p[0] - q[0], b = (double)p[1] - q[1], c = (double)p[2] - q[2];
            if (a*a + b*b + c*c > r2) continue;
            double du = uv[2*v] - uv[2*w], dv = uv[2*v+1] - uv[2*w+1];
            if (du*du + dv*dv <= f2) continue;
            if (accept(v, w, ctx)) return (long)w;
        }
    }
    return -1;
}

typedef struct CrDupCtx { const uint8_t *state; const uint32_t *depth; } CrDupCtx;

/* a fill vertex duplicates observed papyrus, or fill at least as close to the observed papyrus as itself */
static int cr_accept_fill(size_t v, size_t w, const void *ctx)
{
    const CrDupCtx *c = ctx;
    return cr_observed(c->state[w]) || (cr_fill(c->state[w]) && c->depth[v] >= c->depth[w]);
}

/* doubled observed papyrus (each pair once) */
static int cr_accept_observed(size_t v, size_t w, const void *ctx)
{
    const CrDupCtx *c = ctx;
    return w > v && cr_observed(c->state[w]);
}

int CloseupRepair_mark(const float *verts, const double *uv, size_t nv, const int32_t *faces, size_t nf,
                       const uint8_t *state, const float *nrm, const float *inward, uint8_t *reason,
                       CloseupRepairStats *st)
{
    if (!st) return -1;
    memset(st, 0, sizeof *st);
    double started = ves_clock_sec();
    if ((nv && (!verts || !uv || !state || !nrm || !inward || !reason)) || (nf && !faces) || nv > (size_t)INT32_MAX)
        return -1;
    for (size_t f = 0; f < 3*nf; f++) if (faces[f] < 0 || (size_t)faces[f] >= nv) return -1;
    st->verts = nv;
    if (!nv) return 0;
    memset(reason, CLOSEUP_REPAIR_KEEP, nv);
    Arena_T a = Arena_new();
    /* adjacency: each vertex lists the two other corners of each of its faces */
    size_t *first = ARENA_CALLOC(a, nv + 1, sizeof *first);
    for (size_t f = 0; f < nf; f++) for (int k = 0; k < 3; k++) first[(size_t)faces[3*f + (size_t)k] + 1] += 2;
    for (size_t v = 0; v < nv; v++) first[v + 1] += first[v];
    int32_t *adj = ARENA_ALLOC(a, (first[nv] ? first[nv] : 1) * sizeof *adj);
    size_t *at = ARENA_CALLOC(a, nv, sizeof *at);
    for (size_t f = 0; f < nf; f++) for (int k = 0; k < 3; k++) {
        size_t v = (size_t)faces[3*f + (size_t)k];
        adj[first[v] + at[v]++] = faces[3*f + (size_t)((k + 1) % 3)];
        adj[first[v] + at[v]++] = faces[3*f + (size_t)((k + 2) % 3)];
    }
    /* facing = the normal against the direction to the axis, read with the bulk's sign */
    float *facing = ARENA_ALLOC(a, nv * sizeof *facing);
    uint8_t *defined = ARENA_ALLOC(a, nv);
    size_t obs_def = 0, obs_pos = 0;
    for (size_t v = 0; v < nv; v++) {
        const float *n = nrm + 3*v, *d = inward + 3*v;
        double dd = (double)d[0]*d[0] + (double)d[1]*d[1] + (double)d[2]*d[2];
        double nn = (double)n[0]*n[0] + (double)n[1]*n[1] + (double)n[2]*n[2];
        defined[v] = (uint8_t)(isfinite(dd) && isfinite(nn) && dd > 0.25 && nn > 1e-12);
        facing[v] = defined[v] ? (float)(((double)n[0]*d[0] + (double)n[1]*d[1] + (double)n[2]*d[2]) / sqrt(dd*nn)) : 0.0f;
        if (cr_observed(state[v])) { st->observed++; if (defined[v]) obs_def++; }
        else if (cr_fill(state[v])) st->fill++;
        if (!defined[v]) st->undefined++;
        if (cr_observed(state[v]) && defined[v] && facing[v] > 0.0f) obs_pos++;
    }
    if (obs_def && 2*obs_pos < obs_def) {
        st->inward_flipped = 1;
        obs_pos = 0;
        for (size_t v = 0; v < nv; v++) {
            facing[v] = -facing[v];
            if (cr_observed(state[v]) && defined[v] && facing[v] > 0.0f) obs_pos++;
        }
    }
    st->observed_facing_axis = obs_def ? (double)obs_pos / (double)obs_def : 0.0;
    /* FLAP: connected observed vertices facing away from the axis, large and clearly away */
    int32_t *par = ARENA_ALLOC(a, nv * sizeof *par);
    for (size_t v = 0; v < nv; v++) par[v] = (int32_t)v;
    for (size_t v = 0; v < nv; v++) if (cr_away(state, defined, facing, v))
        for (size_t t = first[v]; t < first[v + 1]; t++)
            if (cr_away(state, defined, facing, (size_t)adj[t])) cr_union(par, (int32_t)v, adj[t]);
    size_t *count = ARENA_CALLOC(a, nv, sizeof *count);
    double *sum = ARENA_CALLOC(a, nv, sizeof *sum);
    for (size_t v = 0; v < nv; v++) if (cr_away(state, defined, facing, v)) {
        size_t r = (size_t)cr_find(par, (int32_t)v);
        count[r]++; sum[r] += facing[v];
    }
    for (size_t v = 0; v < nv; v++) if (cr_away(state, defined, facing, v)) {
        size_t r = (size_t)cr_find(par, (int32_t)v);
        int flap = count[r] >= (size_t)CLOSEUP_REPAIR_FLAP_MIN_VERTS && sum[r] / (double)count[r] <= -CLOSEUP_REPAIR_FLAP_FACING;
        if (r == v) { if (flap) st->flap_components++; else st->flap_specks++; }
        if (flap) { reason[v] = CLOSEUP_REPAIR_FLAP; st->flap++; }
    }
    /* WRONG FILL */
    for (size_t v = 0; v < nv; v++)
        if (cr_fill(state[v]) && defined[v] && facing[v] < -(float)CLOSEUP_REPAIR_FILL_FACING) {
            reason[v] = CLOSEUP_REPAIR_WRONG_FILL; st->wrong_fill++;
        }
    /* fill depth: hops from the nearest observed vertex */
    uint32_t *depth = ARENA_ALLOC(a, nv * sizeof *depth);
    int32_t *queue = ARENA_ALLOC(a, nv * sizeof *queue);
    size_t head = 0, tail = 0;
    for (size_t v = 0; v < nv; v++) {
        depth[v] = UINT32_MAX;
        if (cr_observed(state[v])) { depth[v] = 0; queue[tail++] = (int32_t)v; }
    }
    while (head < tail) {
        size_t v = (size_t)queue[head++];
        for (size_t t = first[v]; t < first[v + 1]; t++) {
            size_t w = (size_t)adj[t];
            if (depth[w] == UINT32_MAX) { depth[w] = depth[v] + 1; queue[tail++] = (int32_t)w; }
        }
    }
    /* DUPLICATE: kept vertices in DUP_VOX cells; fill against everything kept, in vertex order */
    size_t nc = 0;
    CrCell *cells = ARENA_ALLOC(a, nv * sizeof *cells);
    for (size_t v = 0; v < nv; v++) if (reason[v] == CLOSEUP_REPAIR_KEEP && state[v]) {
        const float *p = verts + 3*v;
        cells[nc].key = cr_key(cr_cell(p[0]), cr_cell(p[1]), cr_cell(p[2])); cells[nc].v = (int32_t)v; nc++;
    }
    qsort(cells, nc, sizeof *cells, cr_cell_order);
    CrDupCtx ctx = { state, depth };
    for (size_t v = 0; v < nv; v++) {
        if (!cr_fill(state[v]) || reason[v] != CLOSEUP_REPAIR_KEEP) continue;
        long w = cr_duplicate_of(cells, nc, verts, uv, reason, v, cr_accept_fill, &ctx);
        if (w < 0) continue;
        reason[v] = CLOSEUP_REPAIR_DUPLICATE; st->duplicate++;
        if (cr_observed(state[w])) st->pairs_fill_observed++; else st->pairs_fill_fill++;
    }
    for (size_t v = 0; v < nv; v++)
        if (cr_observed(state[v]) && reason[v] == CLOSEUP_REPAIR_KEEP &&
            cr_duplicate_of(cells, nc, verts, uv, reason, v, cr_accept_observed, &ctx) >= 0)
            st->pairs_observed_observed++;
    /* ORPHAN: kept fill whose component over the kept faces has no observed vertex */
    uint8_t *in_face = ARENA_CALLOC(a, nv, 1), *has_obs = ARENA_CALLOC(a, nv, 1);
    for (size_t v = 0; v < nv; v++) par[v] = (int32_t)v;
    for (size_t f = 0; f < nf; f++) {
        const int32_t *t = faces + 3*f;
        if (reason[t[0]] || reason[t[1]] || reason[t[2]]) continue;
        cr_union(par, t[0], t[1]); cr_union(par, t[1], t[2]);
        in_face[t[0]] = in_face[t[1]] = in_face[t[2]] = 1;
    }
    for (size_t v = 0; v < nv; v++) if (in_face[v] && cr_observed(state[v])) has_obs[cr_find(par, (int32_t)v)] = 1;
    for (size_t v = 0; v < nv; v++) {
        if (!in_face[v] || reason[v] != CLOSEUP_REPAIR_KEEP) continue;
        size_t r = (size_t)cr_find(par, (int32_t)v);
        if (has_obs[r]) continue;
        if (r == v) st->orphan_components++;
        reason[v] = CLOSEUP_REPAIR_ORPHAN; st->orphan++;
    }
    Arena_dispose(&a);
    st->seconds = ves_clock_sec() - started;
    return 0;
}

/* ---- selftest -------------------------------------------------------------------------------------------- */

typedef struct CrMesh {
    size_t nv, nf, cap_v, cap_f;
    float *p, *n, *in;
    double *uv;
    int32_t *f;
    uint8_t *s;
} CrMesh;

typedef void (*CrPlace)(size_t i, size_t j, const void *ctx, float p[3], float n[3], uint8_t *state);

/* an ni x nj grid patch at surface position (u0 + i, v0 + j), two faces per cell; the direction to the axis is
 * the horizontal one to the line y = x = 0 */
static size_t cr_patch(CrMesh *m, size_t ni, size_t nj, double u0, double v0, CrPlace place, const void *ctx)
{
    size_t base = m->nv;
    for (size_t j = 0; j < nj; j++) for (size_t i = 0; i < ni; i++) {
        size_t v = m->nv++;
        place(i, j, ctx, m->p + 3*v, m->n + 3*v, m->s + v);
        m->uv[2*v] = u0 + (double)i; m->uv[2*v+1] = v0 + (double)j;
        double y = m->p[3*v+1], x = m->p[3*v+2], r = sqrt(x*x + y*y);
        m->in[3*v] = 0.0f; m->in[3*v+1] = r > 0 ? (float)(-y / r) : 0.0f; m->in[3*v+2] = r > 0 ? (float)(-x / r) : 0.0f;
    }
    for (size_t j = 0; j + 1 < nj; j++) for (size_t i = 0; i + 1 < ni; i++) {
        int32_t a = (int32_t)(base + j*ni + i), b = a + 1, c = (int32_t)(a + (int32_t)ni), d = c + 1;
        int32_t *f = m->f + 3*m->nf;
        f[0] = a; f[1] = b; f[2] = c; f[3] = b; f[4] = d; f[5] = c;
        m->nf += 2;
    }
    return base;
}

/* a patch in the plane y = y0 (facing -y, the axis side, or +y), x = x0 + i, z = z0 + j */
typedef struct CrFlat { float x0, y0, z0, ny; uint8_t state; } CrFlat;

static void cr_flat(size_t i, size_t j, const void *ctx, float p[3], float n[3], uint8_t *state)
{
    const CrFlat *c = ctx;
    p[0] = c->z0 + (float)j; p[1] = c->y0; p[2] = c->x0 + (float)i;
    n[0] = 0.0f; n[1] = c->ny; n[2] = 0.0f; *state = c->state;
}

/* a sheet (observed, i < 100) that runs on as a fill tongue folded back 1 vox above it */
static void cr_tongue(size_t i, size_t j, const void *ctx, float p[3], float n[3], uint8_t *state)
{
    (void)ctx;
    p[0] = (float)j;
    if (i < 100) { p[1] = 500.0f; p[2] = (float)i; *state = 1; }
    else { p[1] = 501.0f; p[2] = (float)(199 - (long)i); *state = 2; }
    n[0] = 0.0f; n[1] = -1.0f; n[2] = 0.0f;
}

/* an observed strip (i < 10) that runs on as fill, x = x0 + sign * i */
typedef struct CrStrip { float x0, y0, sign; } CrStrip;

static void cr_strip(size_t i, size_t j, const void *ctx, float p[3], float n[3], uint8_t *state)
{
    const CrStrip *c = ctx;
    p[0] = 300.0f + (float)j; p[1] = c->y0; p[2] = c->x0 + c->sign * (float)i;
    n[0] = 0.0f; n[1] = -1.0f; n[2] = 0.0f; *state = (uint8_t)(i < 10 ? 1 : 2);
}

static CrMesh cr_mesh(Arena_T a, size_t cap_v)
{
    CrMesh m; memset(&m, 0, sizeof m);
    m.cap_v = cap_v; m.cap_f = 2*cap_v;
    m.p = ARENA_ALLOC(a, 3*cap_v * sizeof *m.p); m.n = ARENA_ALLOC(a, 3*cap_v * sizeof *m.n);
    m.in = ARENA_ALLOC(a, 3*cap_v * sizeof *m.in); m.uv = ARENA_ALLOC(a, 2*cap_v * sizeof *m.uv);
    m.f = ARENA_ALLOC(a, 3*m.cap_f * sizeof *m.f); m.s = ARENA_ALLOC(a, cap_v);
    return m;
}

static size_t cr_count(const uint8_t *reason, size_t from, size_t to, uint8_t want)
{
    size_t c = 0;
    for (size_t v = from; v < to; v++) c += reason[v] == want;
    return c;
}

int CloseupRepair_selftest(void)
{
    int fails = 0;
    Arena_T a = Arena_new();
    CloseupRepairStats st;
    /* 1: a fill tongue folded back over its own sheet: the part more than DUP_PX from its twin goes */
    {
        CrMesh m = cr_mesh(a, 160*30);
        cr_patch(&m, 160, 30, 0.0, 0.0, cr_tongue, NULL);
        uint8_t *reason = ARENA_ALLOC(a, m.nv);
        int rc = CloseupRepair_mark(m.p, m.uv, m.nv, m.f, m.nf, m.s, m.n, m.in, reason, &st);
        /* |uv - twin uv| = |2i - 199|: far past 64 from i 140 on, well inside it up to i 120 (the band between
         * depends on the diagonal partners); the sheet keeps every vertex */
        size_t wrong = 0, dup = 0;
        for (size_t v = 0; v < m.nv; v++) {
            size_t i = v % 160;
            int d = reason[v] == CLOSEUP_REPAIR_DUPLICATE;
            dup += (size_t)d;
            if ((i < 100 && reason[v] != CLOSEUP_REPAIR_KEEP) || (i >= 140 && !d) || (i >= 100 && i <= 120 && d)) wrong++;
        }
        if (rc || wrong || st.duplicate != dup || st.pairs_fill_observed != dup || st.inward_flipped || st.orphan) {
            fprintf(stderr, "  closeup_repair tongue FAIL: rc %d, %zu wrong labels, duplicate %zu (%zu labelled), orphan %zu\n",
                    rc, wrong, st.duplicate, dup, st.orphan);
            fails++;
        }
    }
    /* 2 + 3: a facing sheet, an observed flap and speck facing away, wrong-side fill, an orphan fill island;
     * then the same with every normal negated (the bulk's sign) */
    for (int negate = 0; negate < 2; negate++) {
        CrMesh m = cr_mesh(a, 6400 + 3000 + 100 + 400 + 400);
        CrFlat sheet = { 0.0f, 500.0f, 0.0f, -1.0f, 1 }, flap = { 0.0f, 540.0f, 200.0f, 1.0f, 1 },
               speck = { 0.0f, 580.0f, 400.0f, 1.0f, 1 }, wfill = { 0.0f, 620.0f, 500.0f, 1.0f, 2 },
               orphan = { 0.0f, 660.0f, 600.0f, -1.0f, 2 };
        size_t s0 = cr_patch(&m, 80, 80, 0.0, 0.0, cr_flat, &sheet);
        size_t f0 = cr_patch(&m, 60, 50, 1000.0, 0.0, cr_flat, &flap);
        size_t k0 = cr_patch(&m, 10, 10, 2000.0, 0.0, cr_flat, &speck);
        size_t w0 = cr_patch(&m, 20, 20, 3000.0, 0.0, cr_flat, &wfill);
        size_t o0 = cr_patch(&m, 20, 20, 4000.0, 0.0, cr_flat, &orphan);
        if (negate) for (size_t v = 0; v < 3*m.nv; v++) m.n[v] = -m.n[v];
        uint8_t *reason = ARENA_ALLOC(a, m.nv);
        int rc = CloseupRepair_mark(m.p, m.uv, m.nv, m.f, m.nf, m.s, m.n, m.in, reason, &st);
        int ok = rc == 0 && st.inward_flipped == negate &&
                 cr_count(reason, s0, f0, CLOSEUP_REPAIR_KEEP) == f0 - s0 &&
                 cr_count(reason, f0, k0, CLOSEUP_REPAIR_FLAP) == k0 - f0 &&
                 cr_count(reason, k0, w0, CLOSEUP_REPAIR_KEEP) == w0 - k0 &&
                 cr_count(reason, w0, o0, CLOSEUP_REPAIR_WRONG_FILL) == o0 - w0 &&
                 cr_count(reason, o0, m.nv, CLOSEUP_REPAIR_ORPHAN) == m.nv - o0 &&
                 st.flap_components == 1 && st.flap_specks == 1 && st.orphan_components == 1 && st.duplicate == 0;
        if (!ok) {
            fprintf(stderr, "  closeup_repair flap/fill%s FAIL: rc %d flipped %d, flap %zu (%zu comps, %zu specks), "
                    "wrong fill %zu, orphan %zu, duplicate %zu\n", negate ? " (negated normals)" : "", rc,
                    st.inward_flipped, st.flap, st.flap_components, st.flap_specks, st.wrong_fill, st.orphan, st.duplicate);
            fails++;
        }
    }
    /* 4: two fill strips, each run out from its own observed end, overlapping: the deeper copy goes */
    {
        CrMesh m = cr_mesh(a, 70*10 + 160*10);
        CrStrip one = { 0.0f, 500.0f, 1.0f }, two = { 209.0f, 500.5f, -1.0f };
        size_t a0 = cr_patch(&m, 70, 10, 0.0, 0.0, cr_strip, &one);     /* x 0..69, fill depth x - 9 */
        size_t b0 = cr_patch(&m, 160, 10, 100.0, 0.0, cr_strip, &two);  /* x 209..50, uv = x + 100 - 2 (x - 50) */
        uint8_t *reason = ARENA_ALLOC(a, m.nv);
        int rc = CloseupRepair_mark(m.p, m.uv, m.nv, m.f, m.nf, m.s, m.n, m.in, reason, &st);
        /* the deeper strip's copy goes where the strips overlap (x <= 69, plus the 1-vox reach past its end), and
         * nothing of it where they are apart (x >= 75); the shallower strip keeps every vertex */
        size_t wrong = 0, dup = 0;
        for (size_t v = b0; v < m.nv; v++) {
            double x = 209.0 - (double)((v - b0) % 160);
            int d = reason[v] == CLOSEUP_REPAIR_DUPLICATE;
            dup += (size_t)d;
            if ((x <= 69.0 && !d) || (x >= 75.0 && d)) wrong++;
        }
        if (rc || wrong || cr_count(reason, a0, b0, CLOSEUP_REPAIR_KEEP) != b0 - a0 || st.pairs_fill_fill != dup ||
            st.duplicate != dup || dup < 200) {
            fprintf(stderr, "  closeup_repair fill over fill FAIL: rc %d, %zu wrong labels on the deeper strip, fill-fill pairs "
                    "%zu, duplicates %zu (%zu labelled)\n", rc, wrong, st.pairs_fill_fill, st.duplicate, dup);
            fails++;
        }
    }
    /* 5: a clean sheet is left alone */
    {
        CrMesh m = cr_mesh(a, 50*50);
        CrFlat sheet = { 0.0f, 500.0f, 0.0f, -1.0f, 1 };
        cr_patch(&m, 50, 50, 0.0, 0.0, cr_flat, &sheet);
        uint8_t *reason = ARENA_ALLOC(a, m.nv);
        int rc = CloseupRepair_mark(m.p, m.uv, m.nv, m.f, m.nf, m.s, m.n, m.in, reason, &st);
        if (rc || cr_count(reason, 0, m.nv, CLOSEUP_REPAIR_KEEP) != m.nv || st.observed_facing_axis < 0.999) {
            fprintf(stderr, "  closeup_repair clean sheet FAIL: rc %d\n", rc);
            fails++;
        }
    }
    Arena_dispose(&a);
    fprintf(stderr, "[selftest] closeup_repair (duplicate fill, wrong-side fill, flaps and specks, the bulk's sign, "
            "fill over fill, orphans, a clean sheet) %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails;
}
