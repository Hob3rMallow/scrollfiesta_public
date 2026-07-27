/*
 * weld_cleanup.c -- post-weld sliver / T-junction cleanup (see weld_cleanup.h).
 *
 * Self-contained on purpose: grid_weld does NOT link the QEM solver stack
 * (quadric / cg / cotan / kdtree), so the Surazhsky-Gotsman flip pass here is a
 * faithful copy of the canonical one in src/common/qem.c (qem_edge_flip_pass),
 * kept separate rather than dragging that whole dependency tree into the
 * terminal weld step. The collapse pass is new: a guarded short-edge collapse
 * with the Dey-Edelsbrunner-Guha link condition + a normal-flip guard, run on a
 * freshly-rebuilt adjacency each round (so there is no stale-adjacency
 * bookkeeping to get wrong).
 */
#include "weld_cleanup.h"
#include "../common/pipeline_constants.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Candidate edge of a target triangle, sorted shortest-first for collapse. */
typedef struct { int32_t u, v; double len; } EdgeCand;

/* ===================================================================
 * Geometry helpers (double accumulation off float positions).
 * =================================================================== */
static double tri_area3(const float *V, int32_t a, int32_t b, int32_t c)
{
    double e1x=(double)V[(size_t)b*3+0]-V[(size_t)a*3+0];
    double e1y=(double)V[(size_t)b*3+1]-V[(size_t)a*3+1];
    double e1z=(double)V[(size_t)b*3+2]-V[(size_t)a*3+2];
    double e2x=(double)V[(size_t)c*3+0]-V[(size_t)a*3+0];
    double e2y=(double)V[(size_t)c*3+1]-V[(size_t)a*3+1];
    double e2z=(double)V[(size_t)c*3+2]-V[(size_t)a*3+2];
    double cx=e1y*e2z-e1z*e2y, cy=e1z*e2x-e1x*e2z, cz=e1x*e2y-e1y*e2x;
    return 0.5*sqrt(cx*cx+cy*cy+cz*cz);
}

static double edge_len3(const float *V, int32_t a, int32_t b)
{
    double dx=(double)V[(size_t)a*3+0]-V[(size_t)b*3+0];
    double dy=(double)V[(size_t)a*3+1]-V[(size_t)b*3+1];
    double dz=(double)V[(size_t)a*3+2]-V[(size_t)b*3+2];
    return sqrt(dx*dx+dy*dy+dz*dz);
}

/* Minimum altitude = 2*area / longest edge. Tiny for any sliver (needle OR
 * cap), 0 for a degenerate triangle. The quantity collapse keys on. */
static double tri_min_alt3(const float *V, int32_t a, int32_t b, int32_t c)
{
    double l0=edge_len3(V,a,b), l1=edge_len3(V,b,c), l2=edge_len3(V,c,a);
    double lmax = l0>l1?(l0>l2?l0:l2):(l1>l2?l1:l2);
    if (lmax < 1e-12) return 0.0;
    return 2.0*tri_area3(V,a,b,c)/lmax;
}

static int is_target(const float *V, int32_t a, int32_t b, int32_t c,
                     const WeldCleanupParams *p)
{
    if (tri_area3(V,a,b,c) < p->degen_area) return 1;
    if (tri_min_alt3(V,a,b,c) < p->sliver_min_alt) return 1;
    return 0;
}

/* ===================================================================
 * Boundary detection (verts on any face-count==1 edge). Mirrors
 * qem.c::qem_detect_boundary.
 * =================================================================== */
typedef struct { int32_t v0, v1, face; } HE;
static int he_cmp(const void *pa, const void *pb)
{
    const HE *a=(const HE*)pa, *b=(const HE*)pb;
    if (a->v0 != b->v0) return a->v0 < b->v0 ? -1 : 1;
    if (a->v1 != b->v1) return a->v1 < b->v1 ? -1 : 1;
    return 0;
}

static void detect_boundary(Arena_T arena, const int32_t *faces, size_t nf,
                            size_t nv, uint8_t *is_boundary)
{
    Arena_Mark mark = Arena_save(arena);
    size_t n_he = nf*3, i = 0;
    HE *he = (HE *)ARENA_ALLOC(arena, (long)(n_he*sizeof(HE)));
    for (i=0;i<nf;i++){
        int32_t f0=faces[i*3+0], f1=faces[i*3+1], f2=faces[i*3+2];
        int32_t tri[3]={f0,f1,f2};
        size_t k=0;
        for (k=0;k<3;k++){
            int32_t a=tri[k], b=tri[(k+1)%3];
            if (a>b){ int32_t t=a; a=b; b=t; }
            he[i*3+k].v0=a; he[i*3+k].v1=b; he[i*3+k].face=(int32_t)i;
        }
    }
    qsort(he, n_he, sizeof(HE), he_cmp);
    memset(is_boundary, 0, nv*sizeof(uint8_t));
    i=0;
    while (i < n_he){
        size_t j=i+1;
        while (j<n_he && he[j].v0==he[i].v0 && he[j].v1==he[i].v1) j++;
        if (j-i == 1){ is_boundary[he[i].v0]=1; is_boundary[he[i].v1]=1; }
        i=j;
    }
    Arena_restore(arena, mark);
}

/* ===================================================================
 * Surazhsky-Gotsman min-angle edge flip (faithful copy of
 * qem.c::qem_edge_flip_pass). Vertices never move; boundary edges never flip;
 * a flip is taken only if it strictly improves the min angle AND keeps both new
 * faces' normals on the same side as the original (no fold). One locked pass
 * flips only non-adjacent edges, so iterate to convergence.
 * =================================================================== */
typedef struct { int32_t v0,v1,face,opposite; int8_t fwd; } MHE;
static int mhe_cmp(const void *pa, const void *pb)
{
    const MHE *a=(const MHE*)pa, *b=(const MHE*)pb;
    if (a->v0 != b->v0) return a->v0 < b->v0 ? -1 : 1;
    if (a->v1 != b->v1) return a->v1 < b->v1 ? -1 : 1;
    return 0;
}

static void face_normal(const float *V, int32_t a, int32_t b, int32_t c, float o[3])
{
    float e1[3], e2[3]; int k;
    for (k=0;k<3;k++){ e1[k]=V[(size_t)b*3+(size_t)k]-V[(size_t)a*3+(size_t)k];
                       e2[k]=V[(size_t)c*3+(size_t)k]-V[(size_t)a*3+(size_t)k]; }
    o[0]=e1[1]*e2[2]-e1[2]*e2[1];
    o[1]=e1[2]*e2[0]-e1[0]*e2[2];
    o[2]=e1[0]*e2[1]-e1[1]*e2[0];
}

static float angle_at(const float *V, int32_t a, int32_t b, int32_t c)
{
    float ab[3], ac[3]; int k;
    float dot=0,la=0,lb=0,denom,cs;
    for (k=0;k<3;k++){ ab[k]=V[(size_t)b*3+(size_t)k]-V[(size_t)a*3+(size_t)k];
                       ac[k]=V[(size_t)c*3+(size_t)k]-V[(size_t)a*3+(size_t)k]; }
    dot=ab[0]*ac[0]+ab[1]*ac[1]+ab[2]*ac[2];
    la=sqrtf(ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2]);
    lb=sqrtf(ac[0]*ac[0]+ac[1]*ac[1]+ac[2]*ac[2]);
    denom=la*lb;
    if (denom < 1e-12f) return 0.0f;
    cs=dot/denom;
    if (cs>1.0f) cs=1.0f; if (cs<-1.0f) cs=-1.0f;
    return acosf(cs);
}

static float min_angle(const float *V, int32_t a, int32_t b, int32_t c)
{
    float a1=angle_at(V,a,b,c), a2=angle_at(V,b,c,a), a3=angle_at(V,c,a,b);
    float m=a1; if (a2<m)m=a2; if (a3<m)m=a3; return m;
}

/* True if edge (u,w) already exists in the sorted-by-(v0,v1) MHE list. A flip
 * onto an existing edge gives that edge a 3rd/4th face -> non-manifold. Same
 * guard as qem.c::maint_edge_exists. */
static int mhe_edge_exists(const MHE *mhe, size_t n_he, int32_t u, int32_t w)
{
    int32_t v0 = (u < w) ? u : w;
    int32_t v1 = (u < w) ? w : u;
    size_t lo = 0, hi = n_he;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (mhe[mid].v0 < v0 || (mhe[mid].v0 == v0 && mhe[mid].v1 < v1))
            lo = mid + 1;
        else
            hi = mid;
    }
    return (lo < n_he && mhe[lo].v0 == v0 && mhe[lo].v1 == v1) ? 1 : 0;
}

static size_t flip_pass(Arena_T arena, const float *V, int32_t *faces,
                        size_t nf, size_t nv, const uint8_t *bnd)
{
    Arena_Mark mark = Arena_save(arena);
    size_t n_he=nf*3, f=0, i=0, n_flipped=0;
    MHE *mhe = (MHE *)ARENA_ALLOC(arena, (long)(n_he*sizeof(MHE)));
    for (f=0;f<nf;f++){
        int32_t v[3]={faces[f*3+0],faces[f*3+1],faces[f*3+2]};
        int e;
        for (e=0;e<3;e++){
            int32_t a=v[e], b=v[(e+1)%3], opp=v[(e+2)%3];
            size_t idx=f*3+(size_t)e;
            mhe[idx].v0=(a<b)?a:b; mhe[idx].v1=(a<b)?b:a;
            mhe[idx].face=(int32_t)f; mhe[idx].opposite=opp; mhe[idx].fwd=(int8_t)((a<b)?1:0);
        }
    }
    qsort(mhe, n_he, sizeof(MHE), mhe_cmp);
    uint8_t *vused=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    uint8_t *fused=(uint8_t*)ARENA_ALLOC(arena,(long)(nf*sizeof(uint8_t)));
    memset(vused,0,nv*sizeof(uint8_t));
    memset(fused,0,nf*sizeof(uint8_t));
    i=0;
    while (i+1 < n_he){
        if (mhe[i].v0==mhe[i+1].v0 && mhe[i].v1==mhe[i+1].v1){
            int32_t a=mhe[i].v0, b=mhe[i].v1, c=0, d=0, fc=0, fd=0;
            float n_orig[3], n1[3], n2[3], dot1, dot2, cur, cur2, flp, flp2;
            if (mhe[i].fwd && !mhe[i+1].fwd){
                c=mhe[i].opposite; fc=mhe[i].face; d=mhe[i+1].opposite; fd=mhe[i+1].face;
            } else if (!mhe[i].fwd && mhe[i+1].fwd){
                c=mhe[i+1].opposite; fc=mhe[i+1].face; d=mhe[i].opposite; fd=mhe[i].face;
            } else { i+=2; continue; }   /* inconsistent winding */
            if (bnd[a] && bnd[b]){ i+=2; continue; }
            if (vused[a]||vused[b]||vused[c]||vused[d]||fused[fc]||fused[fd]){ i+=2; continue; }
            if (mhe_edge_exists(mhe, n_he, c, d)){ i+=2; continue; }  /* flip target already an edge -> non-manifold */
            cur=min_angle(V,a,b,c); cur2=min_angle(V,b,a,d); if (cur2<cur)cur=cur2;
            flp=min_angle(V,c,a,d); flp2=min_angle(V,d,b,c); if (flp2<flp)flp=flp2;
            if (flp <= cur + 0.001f){ i+=2; continue; }
            /* Orientation-preserving flip: (c,a,d)+(d,b,c) keeps the quad
             * a->d->b->c boundary directions. Guard and write must match -- a
             * reversed (c,d,a)+(d,c,b) write injects same_dir on folds (see
             * qem.c::qem_edge_flip_pass). */
            face_normal(V,a,b,c,n_orig); face_normal(V,c,a,d,n1); face_normal(V,d,b,c,n2);
            dot1=n_orig[0]*n1[0]+n_orig[1]*n1[1]+n_orig[2]*n1[2];
            dot2=n_orig[0]*n2[0]+n_orig[1]*n2[1]+n_orig[2]*n2[2];
            if (dot1<=0.0f || dot2<=0.0f){ i+=2; continue; }
            faces[fc*3+0]=c; faces[fc*3+1]=a; faces[fc*3+2]=d;
            faces[fd*3+0]=d; faces[fd*3+1]=b; faces[fd*3+2]=c;
            vused[a]=vused[b]=vused[c]=vused[d]=1; fused[fc]=fused[fd]=1;
            n_flipped++; i+=2;
        } else i++;
    }
    Arena_restore(arena, mark);
    return n_flipped;
}

static size_t flip_rounds(Arena_T arena, const float *V, int32_t *faces,
                          size_t nf, size_t nv, int max_rounds)
{
    Arena_Mark mark = Arena_save(arena);
    uint8_t *bnd=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    size_t total=0; int r;
    detect_boundary(arena, faces, nf, nv, bnd);  /* boundary stable across flips */
    for (r=0;r<max_rounds;r++){
        size_t k=flip_pass(arena, V, faces, nf, nv, bnd);
        total+=k;
        if (k==0) break;
    }
    Arena_restore(arena, mark);
    return total;
}

/* ===================================================================
 * Vertex -> incident-face adjacency (CSR), rebuilt each collapse round on the
 * current compacted face array.
 * =================================================================== */
static void build_vf(Arena_T arena, const int32_t *faces, size_t nf, size_t nv,
                     int32_t **out_off, int32_t **out_idx)
{
    int32_t *off=(int32_t*)ARENA_ALLOC(arena,(long)((nv+1)*sizeof(int32_t)));
    int32_t *idx=(int32_t*)ARENA_ALLOC(arena,(long)((nf?nf*3:1)*sizeof(int32_t)));
    size_t i;
    memset(off, 0, (nv+1)*sizeof(int32_t));
    for (i=0;i<nf*3;i++) off[(size_t)faces[i]+1]++;
    for (i=0;i<nv;i++) off[i+1]+=off[i];
    {
        int32_t *cur=(int32_t*)ARENA_ALLOC(arena,(long)((nv)*sizeof(int32_t)));
        size_t f;
        memcpy(cur, off, nv*sizeof(int32_t));
        for (f=0;f<nf;f++)
            for (i=0;i<3;i++){
                int32_t v=faces[f*3+i];
                idx[cur[v]++]=(int32_t)f;
            }
    }
    *out_off=off; *out_idx=idx;
}

#define LINK_MAXDEG 64

/* Dey-Edelsbrunner-Guha link condition: collapsing edge (surv,mov) is manifold-
 * safe iff every vertex adjacent to BOTH endpoints is an apex of a face sharing
 * the edge. Returns 1 if the collapse would create a non-manifold edge (reject),
 * conservatively 1 if a neighbourhood is too large to buffer. */
static int collapse_violates_link(const int32_t *faces,
                                  const int32_t *off, const int32_t *idx,
                                  int32_t surv, int32_t mov)
{
    int32_t nbr_mov[LINK_MAXDEG]; int n_mov=0;
    int32_t apex[2]; int n_apex=0;
    int32_t j;
    /* neighbours of mov + apexes of shared faces */
    for (j=off[mov]; j<off[mov+1]; j++){
        int32_t fi=idx[j];
        int32_t r0=faces[fi*3+0], r1=faces[fi*3+1], r2=faces[fi*3+2];
        int has_surv, k;
        int32_t fr[3];
        if (r0==r1||r1==r2||r0==r2) continue;
        has_surv = (r0==surv)||(r1==surv)||(r2==surv);
        if (has_surv){
            int32_t ap = (r0!=surv && r0!=mov) ? r0 : (r1!=surv && r1!=mov) ? r1 : r2;
            if (n_apex<2) apex[n_apex++]=ap;
        }
        fr[0]=r0; fr[1]=r1; fr[2]=r2;
        for (k=0;k<3;k++){
            int32_t w=fr[k], t, dup=0;
            if (w==surv||w==mov) continue;
            for (t=0;t<n_mov;t++) if (nbr_mov[t]==w){ dup=1; break; }
            if (!dup){ if (n_mov>=LINK_MAXDEG) return 1; nbr_mov[n_mov++]=w; }
        }
    }
    /* any neighbour of surv that is also a neighbour of mov must be an apex */
    for (j=off[surv]; j<off[surv+1]; j++){
        int32_t fi=idx[j];
        int32_t r0=faces[fi*3+0], r1=faces[fi*3+1], r2=faces[fi*3+2];
        int32_t fr[3]; int k;
        if (r0==r1||r1==r2||r0==r2) continue;
        fr[0]=r0; fr[1]=r1; fr[2]=r2;
        for (k=0;k<3;k++){
            int32_t w=fr[k], t, in_mov=0, is_apex=0;
            if (w==surv||w==mov) continue;
            for (t=0;t<n_mov;t++) if (nbr_mov[t]==w){ in_mov=1; break; }
            if (!in_mov) continue;
            for (t=0;t<n_apex;t++) if (apex[t]==w){ is_apex=1; break; }
            if (!is_apex) return 1;   /* common neighbour that is not an apex */
        }
    }
    return 0;
}

/* Normal-flip guard: moving `mov` onto `surv`'s position must not reverse or
 * degenerate any face incident to mov (other than the dying shared faces).
 * Returns 1 if any flip/degeneracy would occur (reject). */
static int collapse_would_flip(const float *V, const int32_t *faces,
                               const int32_t *off, const int32_t *idx,
                               int32_t surv, int32_t mov)
{
    int32_t j;
    for (j=off[mov]; j<off[mov+1]; j++){
        int32_t fi=idx[j];
        int32_t r0=faces[fi*3+0], r1=faces[fi*3+1], r2=faces[fi*3+2];
        float p[3][3], np[3][3], no[3], nn[3], dot, oldlen, newlen;
        int k, m;
        if (r0==r1||r1==r2||r0==r2) continue;
        if (r0==surv||r1==surv||r2==surv) continue;   /* shared face -> dies */
        {
            int32_t fr[3]={r0,r1,r2};
            for (m=0;m<3;m++) for (k=0;k<3;k++){
                p[m][k]=V[(size_t)fr[m]*3+(size_t)k];
                np[m][k]=(fr[m]==mov)?V[(size_t)surv*3+(size_t)k]:V[(size_t)fr[m]*3+(size_t)k];
            }
        }
        {
            float e1[3],e2[3];
            for (k=0;k<3;k++){ e1[k]=p[1][k]-p[0][k]; e2[k]=p[2][k]-p[0][k]; }
            no[0]=e1[1]*e2[2]-e1[2]*e2[1]; no[1]=e1[2]*e2[0]-e1[0]*e2[2]; no[2]=e1[0]*e2[1]-e1[1]*e2[0];
            for (k=0;k<3;k++){ e1[k]=np[1][k]-np[0][k]; e2[k]=np[2][k]-np[0][k]; }
            nn[0]=e1[1]*e2[2]-e1[2]*e2[1]; nn[1]=e1[2]*e2[0]-e1[0]*e2[2]; nn[2]=e1[0]*e2[1]-e1[1]*e2[0];
        }
        oldlen=no[0]*no[0]+no[1]*no[1]+no[2]*no[2];
        newlen=nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2];
        if (oldlen < 1e-20f) continue;          /* was already degenerate */
        if (newlen < 1e-20f) return 1;          /* would become degenerate */
        dot=no[0]*nn[0]+no[1]*nn[1]+no[2]*nn[2];
        if (dot < 0.0f) return 1;               /* normal reversed */
    }
    return 0;
}

static void lock_ring(uint8_t *locked, const int32_t *off, const int32_t *idx,
                      const int32_t *faces, int32_t x)
{
    int32_t j;
    locked[x]=1;
    for (j=off[x]; j<off[x+1]; j++){
        int32_t fi=idx[j];
        locked[faces[fi*3+0]]=1; locked[faces[fi*3+1]]=1; locked[faces[fi*3+2]]=1;
    }
}

/* One collapse round. Scans target faces, accepts independent guarded collapses
 * (1-ring locking keeps them non-interacting), then applies the remap and
 * compacts `faces` in place. Returns collapses accepted; *nf updated. */
static size_t collapse_round(Arena_T arena, const float *V,
                             int32_t *faces, size_t *nf, size_t nv,
                             const WeldCleanupParams *p)
{
    Arena_Mark mark = Arena_save(arena);
    int32_t *off=NULL, *idx=NULL;
    uint8_t *bnd=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    uint8_t *locked=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    int32_t *remap=(int32_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(int32_t)));
    size_t f, ncoll=0, w=0, i;
    build_vf(arena, faces, *nf, nv, &off, &idx);
    detect_boundary(arena, faces, *nf, nv, bnd);
    memset(locked, 0, nv*sizeof(uint8_t));
    for (i=0;i<nv;i++) remap[i]=(int32_t)i;

    for (f=0; f<*nf; f++){
        int32_t a=faces[f*3+0], b=faces[f*3+1], c=faces[f*3+2];
        EdgeCand E[3];
        int t, s;
        if (a==b||b==c||a==c) continue;
        if (!is_target(V,a,b,c,p)) continue;
        E[0].u=a; E[0].v=b; E[1].u=b; E[1].v=c; E[2].u=c; E[2].v=a;
        for (t=0;t<3;t++) E[t].len=edge_len3(V,E[t].u,E[t].v);
        /* sort 3 edges ascending by length */
        for (t=0;t<2;t++) for (s=t+1;s<3;s++)
            if (E[s].len < E[t].len){ EdgeCand tmp=E[t]; E[t]=E[s]; E[s]=tmp; }
        /* ONLY the shortest edge is a collapse candidate -- it is the needle's
         * short edge. Never fall through to a sliver's LONG edges: those are
         * legitimate geometry (e.g. a bipyramid apex), and collapsing one when
         * the short edge is blocked (boundary / link / fold) would delete real
         * surface. A cap with no genuinely-short edge is left for the flip pass. */
        {
            int32_t u=E[0].u, v=E[0].v, surv, mov;
            if (locked[u] || locked[v]) continue;
            if (E[0].len > p->max_collapse_len) continue;
            if (bnd[u] && bnd[v]) continue;
            surv = bnd[u] ? u : (bnd[v] ? v : u);
            mov  = (surv==u) ? v : u;
            if (collapse_violates_link(faces, off, idx, surv, mov)) continue;
            if (collapse_would_flip(V, faces, off, idx, surv, mov)) continue;
            remap[mov]=surv;
            lock_ring(locked, off, idx, faces, mov);
            lock_ring(locked, off, idx, faces, surv);
            ncoll++;
        }
    }

    /* apply remap + compact in place (write index never overtakes read index) */
    for (f=0; f<*nf; f++){
        int32_t a=remap[faces[f*3+0]], b=remap[faces[f*3+1]], c=remap[faces[f*3+2]];
        if (a==b||b==c||a==c) continue;
        faces[w*3+0]=a; faces[w*3+1]=b; faces[w*3+2]=c; w++;
    }
    *nf=w;
    Arena_restore(arena, mark);
    return ncoll;
}

/* ===================================================================
 * Public API.
 * =================================================================== */
void WeldCleanup_default_params(WeldCleanupParams *p)
{
    assert(p);
    p->sliver_min_alt     = WELD_CLEANUP_SLIVER_MIN_ALT;
    p->degen_area         = WELD_CLEANUP_DEGEN_AREA;
    p->max_collapse_len   = WELD_CLEANUP_MAX_COLLAPSE_LEN;
    p->flip_max_rounds    = WELD_CLEANUP_FLIP_ROUNDS;
    p->collapse_max_rounds= WELD_CLEANUP_COLLAPSE_ROUNDS;
}

static size_t count_targets(const float *V, const int32_t *faces, size_t nf,
                            const WeldCleanupParams *p)
{
    size_t f, n=0;
    for (f=0;f<nf;f++){
        int32_t a=faces[f*3+0], b=faces[f*3+1], c=faces[f*3+2];
        if (a==b||b==c||a==c) continue;
        if (is_target(V,a,b,c,p)) n++;
    }
    return n;
}

int WeldCleanup_process(Arena_T arena, ComponentMesh *cm,
                        const WeldCleanupParams *params,
                        WeldCleanupStats *stats)
{
    WeldCleanupParams p;
    int32_t *wf = NULL;
    size_t nf, nv, total_flips=0, total_collapses=0;
    int r;

    assert(arena);
    assert(cm);
    assert(cm->self == cm);

    if (params) p = *params; else WeldCleanup_default_params(&p);

    nf = cm->nf; nv = cm->nv;
    if (stats){ memset(stats, 0, sizeof *stats); stats->faces_in=nf; stats->faces_out=nf; }
    if (nf == 0 || nv == 0 || cm->faces == NULL || cm->verts == NULL) return 0;

    if (stats) stats->targets_in = count_targets(cm->verts, cm->faces, nf, &p);

    /* Working face array (flips mutate in place; collapse compacts in place). */
    wf = (int32_t *)ARENA_ALLOC(arena, (long)(nf*3*sizeof(int32_t)));
    memcpy(wf, cm->faces, nf*3*sizeof(int32_t));

    /* Pass 1: flip-first (clear cap slivers, no vertex removed). */
    total_flips += flip_rounds(arena, cm->verts, wf, nf, nv, p.flip_max_rounds);

    /* Pass 2: guarded collapse for the residue (needles, zero-area, T-caps). */
    for (r=0; r<p.collapse_max_rounds; r++){
        size_t k = collapse_round(arena, cm->verts, wf, &nf, nv, &p);
        total_collapses += k;
        if (k == 0) break;
    }

    /* Pass 3: one more flip-to-convergence to tidy caps the collapse exposed. */
    if (total_collapses > 0)
        total_flips += flip_rounds(arena, cm->verts, wf, nf, nv, p.flip_max_rounds);

    cm->faces = wf;
    cm->nf = nf;

    if (stats){
        stats->n_flips = total_flips;
        stats->n_collapses = total_collapses;
        stats->faces_out = nf;
        stats->targets_out = count_targets(cm->verts, wf, nf, &p);
    }
    return 0;
}
