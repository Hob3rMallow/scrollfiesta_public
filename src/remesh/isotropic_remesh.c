/*
 * isotropic_remesh.c -- incremental isotropic remeshing (see isotropic_remesh.h).
 *
 * Self-contained on purpose (same rationale as weld_cleanup.c): the flat-array
 * flip / collapse machinery + Dey-Edelsbrunner-Guha link condition + normal-fold
 * guard are faithful copies of the canonical ones in src/common/qem.c and
 * src/remesh/weld_cleanup.c, rebuilt from the flat arrays each sub-pass so there
 * is no stale-adjacency bookkeeping. The genuinely new mechanic here is the
 * winding-safe edge SPLIT (2 tris -> 4) with array growth.
 *
 * SAFETY POSTURE (the removed BoundaryRemesh fragmented meshes via collapse):
 *   - boundary is FROZEN (pin_mask | topological open edge); never split /
 *     collapsed / flipped / moved.
 *   - collapse only touches edges shorter than 4/5 L (never hole rims), gated by
 *     link + fold + a max-length cap < the ~7 vox inter-wrap clearance.
 *   - FAIL-CLOSED: a combinatorial manifold audit runs last; on any regression
 *     the input is returned verbatim. Worst case the pass is a no-op.
 */
#include "isotropic_remesh.h"
#include "../common/pipeline_constants.h"
#include "../common/mesh_manifold.h"
#include "../common/mls_project.h"
#include "../common/csr.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ===================================================================
 * Geometry helpers (double accumulation off float positions).
 * =================================================================== */
static double edge_len3(const float *V, int32_t a, int32_t b)
{
    double dx=(double)V[(size_t)a*3+0]-V[(size_t)b*3+0];
    double dy=(double)V[(size_t)a*3+1]-V[(size_t)b*3+1];
    double dz=(double)V[(size_t)a*3+2]-V[(size_t)b*3+2];
    return sqrt(dx*dx+dy*dy+dz*dz);
}

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

static double mesh_area(const float *V, const int32_t *F, size_t nf)
{
    double a=0.0; size_t f;
    for (f=0;f<nf;f++){
        int32_t i0=F[f*3+0], i1=F[f*3+1], i2=F[f*3+2];
        if (i0==i1||i1==i2||i0==i2) continue;
        a += tri_area3(V,i0,i1,i2);
    }
    return a;
}

/* ===================================================================
 * Boundary detection (verts on any face-count==1 edge). Copy of
 * weld_cleanup.c::detect_boundary.
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
    HE *he = (HE *)ARENA_ALLOC(arena, (long)((n_he?n_he:1)*sizeof(HE)));
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

/* frozen[v] = pin[v] || boundary[v]. Recomputed from the current faces. */
static void compute_frozen(Arena_T arena, const int32_t *faces, size_t nf,
                           size_t nv, const uint8_t *pin, uint8_t *frozen)
{
    size_t i;
    detect_boundary(arena, faces, nf, nv, frozen);
    if (pin) for (i=0;i<nv;i++) if (pin[i]) frozen[i]=1;
}

/* ===================================================================
 * Surazhsky-Gotsman min-angle edge flip. Copy of weld_cleanup.c::flip_pass,
 * extended to skip any edge whose BOTH endpoints are frozen (pin|boundary).
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
                        size_t nf, size_t nv, const uint8_t *frozen)
{
    Arena_Mark mark = Arena_save(arena);
    size_t n_he=nf*3, f=0, i=0, n_flipped=0;
    MHE *mhe = (MHE *)ARENA_ALLOC(arena, (long)((n_he?n_he:1)*sizeof(MHE)));
    uint8_t *vused, *fused;
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
    vused=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    fused=(uint8_t*)ARENA_ALLOC(arena,(long)(nf*sizeof(uint8_t)));
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
            if (frozen[a] && frozen[b]){ i+=2; continue; }
            if (vused[a]||vused[b]||vused[c]||vused[d]||fused[fc]||fused[fd]){ i+=2; continue; }
            if (mhe_edge_exists(mhe, n_he, c, d)){ i+=2; continue; }
            cur=min_angle(V,a,b,c); cur2=min_angle(V,b,a,d); if (cur2<cur)cur=cur2;
            flp=min_angle(V,c,a,d); flp2=min_angle(V,d,b,c); if (flp2<flp)flp=flp2;
            if (flp <= cur + 0.001f){ i+=2; continue; }
            /* Orientation-preserving flip of quad a->d->b->c->a across diagonal
             * c-d: faces (c,a,d)+(d,b,c) keep every boundary-edge direction. Guard
             * normals AND the write must both use these -- a reversed write
             * (c,d,a)+(d,c,b) negates them, injecting same_dir edges on folded
             * geometry (see qem.c::qem_edge_flip_pass, same bug). */
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
                          size_t nf, size_t nv, const uint8_t *pin, int max_rounds)
{
    Arena_Mark mark = Arena_save(arena);
    uint8_t *frozen=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    size_t total=0; int r;
    compute_frozen(arena, faces, nf, nv, pin, frozen);  /* stable across flips */
    for (r=0;r<max_rounds;r++){
        size_t k=flip_pass(arena, V, faces, nf, nv, frozen);
        total+=k;
        if (k==0) break;
    }
    Arena_restore(arena, mark);
    return total;
}

/* ===================================================================
 * Vertex -> incident-face adjacency (CSR). Copy of weld_cleanup.c::build_vf.
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

/* Dey-Edelsbrunner-Guha link condition. Copy of weld_cleanup.c. */
static int collapse_violates_link(const int32_t *faces,
                                  const int32_t *off, const int32_t *idx,
                                  int32_t surv, int32_t mov)
{
    int32_t nbr_mov[LINK_MAXDEG]; int n_mov=0;
    int32_t apex[2]; int n_apex=0;
    int32_t j;
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
            if (!is_apex) return 1;
        }
    }
    return 0;
}

/* Normal-flip guard. Copy of weld_cleanup.c. */
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
        if (r0==surv||r1==surv||r2==surv) continue;
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
        if (oldlen < 1e-20f) continue;
        if (newlen < 1e-20f) return 1;
        dot=no[0]*nn[0]+no[1]*nn[1]+no[2]*nn[2];
        if (dot < 0.0f) return 1;
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

/* One length-based collapse round: collapse every independent edge shorter than
 * `thresh` (and <= max_len), gated by frozen / link / fold. A frozen endpoint is
 * always the survivor (never removed / moved); both-frozen edges are skipped.
 * Compacts `faces` in place; *nf updated. Orphaned verts are left in place and
 * reclaimed by compact_verts at end of iteration. Returns collapses accepted. */
static size_t collapse_round(Arena_T arena, const float *V,
                             int32_t *faces, size_t *nf, size_t nv,
                             const uint8_t *pin, double thresh, double max_len,
                             size_t min_faces)
{
    Arena_Mark mark = Arena_save(arena);
    int32_t *off=NULL, *idx=NULL;
    uint8_t *frozen=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    uint8_t *locked=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    int32_t *remap=(int32_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(int32_t)));
    size_t f, ncoll=0, w=0, i;
    build_vf(arena, faces, *nf, nv, &off, &idx);
    compute_frozen(arena, faces, *nf, nv, pin, frozen);
    memset(locked, 0, nv*sizeof(uint8_t));
    for (i=0;i<nv;i++) remap[i]=(int32_t)i;

    for (f=0; f<*nf; f++){
        int32_t tri[3]={faces[f*3+0],faces[f*3+1],faces[f*3+2]};
        int e;
        if (*nf - 2*ncoll <= min_faces) break;   /* face floor reached */
        if (tri[0]==tri[1]||tri[1]==tri[2]||tri[0]==tri[2]) continue;
        for (e=0;e<3;e++){
            int32_t u=tri[e], v=tri[(e+1)%3], surv, mov;
            double L;
            if (locked[u] || locked[v]) continue;
            L = edge_len3(V,u,v);
            if (L >= thresh) continue;
            if (L > max_len) continue;               /* anti-fusion length cap */
            if (frozen[u] && frozen[v]) continue;
            surv = frozen[u] ? u : (frozen[v] ? v : (u<v?u:v));
            mov  = (surv==u) ? v : u;
            if (collapse_violates_link(faces, off, idx, surv, mov)) continue;
            if (collapse_would_flip(V, faces, off, idx, surv, mov)) continue;
            remap[mov]=surv;
            lock_ring(locked, off, idx, faces, mov);
            lock_ring(locked, off, idx, faces, surv);
            ncoll++;
            break;                                   /* one collapse per face per round */
        }
    }

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
 * Edge split (the new mechanic). Splits every independent interior edge longer
 * than `thresh` at its midpoint: 2 triangles -> 4, +1 vert / +2 faces each.
 *
 * Winding-safe by construction. For a manifold interior edge (a<b) the face
 * that traverses a->b (fwd) is fc, the one traversing b->a is fd. Splitting at
 * m:  fc  in place b->m   (-> the a-side sub-tri),  copy a->m (-> the b-side);
 *     fd  in place a->m   (-> the b-side sub-tri),  copy b->m (-> the a-side).
 * This makes every one of the four new shared edges (m-a,m-b,m-c,m-d) carry two
 * oppositely-wound faces, so m is a clean degree-4 interior manifold vertex.
 *
 * Two-phase and face-locked so all splits in a round are independent (a face is
 * split at most once/round; a triangle needing >1 split resolves over rounds).
 * Grows the working arrays; the previous arrays become arena garbage. Boundary
 * edges appear as MHE singletons and are structurally never split. Both-frozen
 * edges are skipped so the frozen skeleton is never subdivided. */
typedef struct { int32_t a, b, fc, fd; } SplitReq;

static size_t split_round(Arena_T arena,
                          float **pV, size_t *pnv,
                          int32_t **pF, size_t *pnf,
                          uint8_t **ppin,
                          double thresh, size_t max_faces)
{
    Arena_Mark mark = Arena_save(arena);
    float   *V  = *pV;
    int32_t *F  = *pF;
    uint8_t *pin= *ppin;
    size_t nv=*pnv, nf=*pnf, n_he=nf*3, f, i, nsplit=0;
    MHE *mhe;
    uint8_t *floc;                 /* frozen (pin|boundary) */
    uint8_t *fdone;                /* per-face split lock    */
    SplitReq *req;
    size_t nreq=0;
    float   *nV; int32_t *nF; uint8_t *npin;

    if (nf==0){ Arena_restore(arena, mark); return 0; }

    /* frozen skeleton (do not subdivide an edge between two frozen verts) */
    floc=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    compute_frozen(arena, F, nf, nv, pin, floc);

    mhe=(MHE*)ARENA_ALLOC(arena,(long)(n_he*sizeof(MHE)));
    for (f=0;f<nf;f++){
        int32_t v[3]={F[f*3+0],F[f*3+1],F[f*3+2]};
        int e;
        for (e=0;e<3;e++){
            int32_t a=v[e], b=v[(e+1)%3], opp=v[(e+2)%3];
            size_t idx=f*3+(size_t)e;
            mhe[idx].v0=(a<b)?a:b; mhe[idx].v1=(a<b)?b:a;
            mhe[idx].face=(int32_t)f; mhe[idx].opposite=opp; mhe[idx].fwd=(int8_t)((a<b)?1:0);
        }
    }
    qsort(mhe, n_he, sizeof(MHE), mhe_cmp);

    fdone=(uint8_t*)ARENA_ALLOC(arena,(long)(nf*sizeof(uint8_t)));
    memset(fdone,0,nf*sizeof(uint8_t));
    req=(SplitReq*)ARENA_ALLOC(arena,(long)((n_he/2+1)*sizeof(SplitReq)));

    i=0;
    while (i+1 < n_he){
        if (nf + 2*nreq >= max_faces) break;   /* face budget reached */
        if (mhe[i].v0==mhe[i+1].v0 && mhe[i].v1==mhe[i+1].v1){
            int32_t a=mhe[i].v0, b=mhe[i].v1, fc, fd;
            if (mhe[i].fwd==mhe[i+1].fwd){ i+=2; continue; }   /* same_dir: skip */
            if (floc[a] && floc[b]){ i+=2; continue; }
            if (edge_len3(V,a,b) <= thresh){ i+=2; continue; }
            if (mhe[i].fwd){ fc=mhe[i].face; fd=mhe[i+1].face; }
            else           { fc=mhe[i+1].face; fd=mhe[i].face; }
            if (fdone[fc] || fdone[fd]){ i+=2; continue; }
            req[nreq].a=a; req[nreq].b=b; req[nreq].fc=fc; req[nreq].fd=fd; nreq++;
            fdone[fc]=1; fdone[fd]=1;
            i+=2;
        } else i++;
    }

    if (nreq==0){ Arena_restore(arena, mark); return 0; }

    /* grow */
    {
        size_t new_nv=nv+nreq, new_nf=nf+2*nreq, r;
        nV  =(float*)  ARENA_ALLOC(arena,(long)(new_nv*3*sizeof(float)));
        nF  =(int32_t*)ARENA_ALLOC(arena,(long)(new_nf*3*sizeof(int32_t)));
        npin=(uint8_t*)ARENA_ALLOC(arena,(long)(new_nv*sizeof(uint8_t)));
        memcpy(nV, V, nv*3*sizeof(float));
        memcpy(nF, F, nf*3*sizeof(int32_t));
        if (pin) memcpy(npin, pin, nv*sizeof(uint8_t)); else memset(npin,0,nv*sizeof(uint8_t));

        for (r=0;r<nreq;r++){
            int32_t a=req[r].a, b=req[r].b, fc=req[r].fc, fd=req[r].fd;
            int32_t m=(int32_t)(nv+r);
            int32_t Fc2[3], Fd2[3]; int k;
            nV[(size_t)m*3+0]=0.5f*(nV[(size_t)a*3+0]+nV[(size_t)b*3+0]);
            nV[(size_t)m*3+1]=0.5f*(nV[(size_t)a*3+1]+nV[(size_t)b*3+1]);
            nV[(size_t)m*3+2]=0.5f*(nV[(size_t)a*3+2]+nV[(size_t)b*3+2]);
            npin[m]=0;                                   /* new interior vertex */
            /* fc traverses a->b: copy a->m (b-side sub-tri), in place b->m */
            for (k=0;k<3;k++) Fc2[k]=(nF[(size_t)fc*3+(size_t)k]==a)?m:nF[(size_t)fc*3+(size_t)k];
            for (k=0;k<3;k++) if (nF[(size_t)fc*3+(size_t)k]==b) nF[(size_t)fc*3+(size_t)k]=m;
            /* fd traverses b->a: copy b->m (a-side sub-tri), in place a->m */
            for (k=0;k<3;k++) Fd2[k]=(nF[(size_t)fd*3+(size_t)k]==b)?m:nF[(size_t)fd*3+(size_t)k];
            for (k=0;k<3;k++) if (nF[(size_t)fd*3+(size_t)k]==a) nF[(size_t)fd*3+(size_t)k]=m;
            for (k=0;k<3;k++) nF[(nf+2*r+0)*3+(size_t)k]=Fc2[k];
            for (k=0;k<3;k++) nF[(nf+2*r+1)*3+(size_t)k]=Fd2[k];
        }
        *pV=nV; *pF=nF; *ppin=npin; *pnv=new_nv; *pnf=new_nf;
        nsplit=nreq;
    }
    /* Do NOT Arena_restore: nV/nF/npin must survive (bump-allocated on arena). */
    (void)mark;
    return nsplit;
}

/* ===================================================================
 * Tangential Laplacian relax (adapt of qem.c::qem_smooth_pass). Interior verts
 * move to the tangent-plane projection of their 1-ring mean; frozen verts are
 * held fixed. Displacement clamped to keep it stable.
 * =================================================================== */
static void compute_vert_normals(const float *V, size_t nv,
                                 const int32_t *F, size_t nf, float *N)
{
    size_t f, i;
    memset(N, 0, nv*3*sizeof(float));
    for (f=0;f<nf;f++){
        int32_t i0=F[f*3+0], i1=F[f*3+1], i2=F[f*3+2];
        float e1[3], e2[3], c[3]; int k;
        if (i0==i1||i1==i2||i0==i2) continue;
        for (k=0;k<3;k++){ e1[k]=V[(size_t)i1*3+(size_t)k]-V[(size_t)i0*3+(size_t)k];
                           e2[k]=V[(size_t)i2*3+(size_t)k]-V[(size_t)i0*3+(size_t)k]; }
        c[0]=e1[1]*e2[2]-e1[2]*e2[1];
        c[1]=e1[2]*e2[0]-e1[0]*e2[2];
        c[2]=e1[0]*e2[1]-e1[1]*e2[0];   /* area-weighted (unnormalized) */
        for (k=0;k<3;k++){ N[(size_t)i0*3+(size_t)k]+=c[k];
                           N[(size_t)i1*3+(size_t)k]+=c[k];
                           N[(size_t)i2*3+(size_t)k]+=c[k]; }
    }
    for (i=0;i<nv;i++){
        float *n=&N[i*3];
        float l=sqrtf(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        if (l>1e-20f){ n[0]/=l; n[1]/=l; n[2]/=l; }
        else { n[0]=n[1]=n[2]=0.0f; }
    }
}

static void relax_pass(Arena_T arena, float *V, size_t nv,
                       const int32_t *F, size_t nf,
                       const uint8_t *pin, float lambda, float max_disp)
{
    Arena_Mark mark = Arena_save(arena);
    uint8_t *frozen=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    float   *N =(float*)ARENA_ALLOC(arena,(long)(nv*3*sizeof(float)));
    float   *disp=(float*)ARENA_ALLOC(arena,(long)(nv*3*sizeof(float)));
    CSR_T adj; const int32_t *off; const int32_t *tgt;
    size_t i;
    compute_frozen(arena, F, nf, nv, pin, frozen);
    compute_vert_normals(V, nv, F, nf, N);
    adj = CSR_from_faces(arena, F, nf, nv);
    off = CSR_offset(adj); tgt = CSR_target(adj);
    memset(disp, 0, nv*3*sizeof(float));

    for (i=0;i<nv;i++){
        int32_t degree, j; float mean[3]={0,0,0}, d[3], dot, nx,ny,nz, inv_d, dl;
        if (frozen[i]) continue;
        degree = off[(int32_t)i+1]-off[(int32_t)i];
        if (degree==0) continue;
        for (j=off[(int32_t)i]; j<off[(int32_t)i+1]; j++){
            int32_t nb=tgt[j];
            mean[0]+=V[nb*3+0]; mean[1]+=V[nb*3+1]; mean[2]+=V[nb*3+2];
        }
        inv_d=1.0f/(float)degree;
        d[0]=mean[0]*inv_d - V[i*3+0];
        d[1]=mean[1]*inv_d - V[i*3+1];
        d[2]=mean[2]*inv_d - V[i*3+2];
        nx=N[i*3+0]; ny=N[i*3+1]; nz=N[i*3+2];
        dot=d[0]*nx+d[1]*ny+d[2]*nz;       /* project out normal component */
        d[0]-=dot*nx; d[1]-=dot*ny; d[2]-=dot*nz;
        d[0]*=lambda; d[1]*=lambda; d[2]*=lambda;
        dl=sqrtf(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
        if (max_disp>0.0f && dl>max_disp){ float s=max_disp/dl; d[0]*=s; d[1]*=s; d[2]*=s; }
        disp[i*3+0]=d[0]; disp[i*3+1]=d[1]; disp[i*3+2]=d[2];
    }
    for (i=0;i<nv;i++){
        if (frozen[i]) continue;
        V[i*3+0]+=disp[i*3+0]; V[i*3+1]+=disp[i*3+1]; V[i*3+2]+=disp[i*3+2];
    }
    Arena_restore(arena, mark);
}

/* Opt-in MLS reprojection onto the current point set (LOP midpoint), frozen
 * verts restored to their pre-projection positions. Keeps verts on the sheet;
 * radius kept local so it cannot reach across the ~7 vox inter-wrap gap. */
static void mls_reproject(Arena_T arena, float *V, size_t nv,
                          const uint8_t *pin, const int32_t *F, size_t nf,
                          float radius)
{
    Arena_Mark mark = Arena_save(arena);
    uint8_t *frozen=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    float *save=(float*)ARENA_ALLOC(arena,(long)(nv*3*sizeof(float)));
    float *outn=(float*)ARENA_ALLOC(arena,(long)(nv*3*sizeof(float)));
    float cell[3]={0.0f,0.0f,0.0f};
    size_t i;
    compute_frozen(arena, F, nf, nv, pin, frozen);
    memcpy(save, V, nv*3*sizeof(float));
    MLS_project_verts(arena, V, nv, radius, cell, V, outn);  /* in-place OK */
    for (i=0;i<nv;i++) if (frozen[i]){
        V[i*3+0]=save[i*3+0]; V[i*3+1]=save[i*3+1]; V[i*3+2]=save[i*3+2];
    }
    Arena_restore(arena, mark);
}

/* Drop orphaned (unreferenced) verts; remap faces + pin to a dense range. */
static void compact_verts(Arena_T arena, float **pV, size_t *pnv,
                          int32_t *F, size_t nf, uint8_t **ppin)
{
    float *V=*pV; uint8_t *pin=*ppin; size_t nv=*pnv, i;
    int32_t *remap=(int32_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(int32_t)));
    uint8_t *used=(uint8_t*)ARENA_ALLOC(arena,(long)(nv*sizeof(uint8_t)));
    float *nV; uint8_t *npin; size_t nn=0;
    memset(used,0,nv*sizeof(uint8_t));
    for (i=0;i<nf*3;i++) used[(size_t)F[i]]=1;
    for (i=0;i<nv;i++) remap[i]=-1;
    nV  =(float*)  ARENA_ALLOC(arena,(long)((nv?nv:1)*3*sizeof(float)));
    npin=(uint8_t*)ARENA_ALLOC(arena,(long)((nv?nv:1)*sizeof(uint8_t)));
    for (i=0;i<nv;i++){
        if (!used[i]) continue;
        remap[i]=(int32_t)nn;
        nV[nn*3+0]=V[i*3+0]; nV[nn*3+1]=V[i*3+1]; nV[nn*3+2]=V[i*3+2];
        npin[nn]=pin?pin[i]:0;
        nn++;
    }
    for (i=0;i<nf*3;i++) F[i]=remap[(size_t)F[i]];
    *pV=nV; *ppin=npin; *pnv=nn;
}

/* ===================================================================
 * Public API.
 * =================================================================== */
void Remesh_default_opts(RemeshOpts *o)
{
    assert(o);
    o->target_edge_len    = 0.0f;                       /* auto */
    o->n_iters            = REMESH_ITERS;
    o->split_ratio        = REMESH_SPLIT_RATIO;
    o->collapse_ratio     = REMESH_COLLAPSE_RATIO;
    o->do_split           = 1;
    o->do_collapse        = 1;
    o->do_flip            = 1;
    o->do_relax           = 1;
    o->reproject          = REMESH_REPROJECT_TANGENTIAL;
    o->relax_lambda       = REMESH_RELAX_LAMBDA;
    o->max_collapse_len   = REMESH_MAX_COLLAPSE_LEN;
    o->split_max_rounds   = REMESH_SPLIT_MAX_ROUNDS;
    o->collapse_max_rounds= REMESH_COLLAPSE_MAX_ROUNDS;
    o->flip_max_rounds    = REMESH_FLIP_MAX_ROUNDS;
    o->face_count_tol     = REMESH_FACE_COUNT_TOL;
}

/* Verbatim manifold copy of the input into fresh arena arrays. */
static int emit_copy(Arena_T arena,
                     const float *in_verts, size_t in_nv,
                     const int32_t *in_faces, size_t in_nf,
                     float **out_verts, size_t *out_nv,
                     int32_t **out_faces, size_t *out_nf)
{
    float *ov=(float*)ARENA_ALLOC(arena,(long)((in_nv?in_nv:1)*3*sizeof(float)));
    int32_t *of=(int32_t*)ARENA_ALLOC(arena,(long)((in_nf?in_nf:1)*3*sizeof(int32_t)));
    if (in_nv) memcpy(ov, in_verts, in_nv*3*sizeof(float));
    if (in_nf) memcpy(of, in_faces, in_nf*3*sizeof(int32_t));
    *out_verts=ov; *out_nv=in_nv; *out_faces=of; *out_nf=in_nf;
    return 1;
}

int Remesh_isotropic(Arena_T arena,
                     const float *in_verts, size_t in_nv,
                     const int32_t *in_faces, size_t in_nf,
                     const uint8_t *pin_mask,
                     const RemeshOpts *opts,
                     float **out_verts, size_t *out_nv,
                     int32_t **out_faces, size_t *out_nf)
{
    RemeshOpts o;
    double A, L, Lhi, Llo, maxlen;
    float *wv; int32_t *wf; uint8_t *wpin;
    size_t wnv, wnf, F, it, hi, lo;
    MeshManifoldStats st;
    int verbose = (getenv("REMESH_VERBOSE") != NULL);

    if (!arena || !out_verts || !out_nv || !out_faces || !out_nf) return -1;
    if (opts) o=*opts; else Remesh_default_opts(&o);

    /* Trivial / too-small: nothing to gain, return input verbatim. */
    if (in_nv==0 || in_nf<4 || in_verts==NULL || in_faces==NULL)
        return emit_copy(arena,in_verts,in_nv,in_faces,in_nf,out_verts,out_nv,out_faces,out_nf);

    F = in_nf;
    A = mesh_area(in_verts, in_faces, in_nf);
    if (A <= 0.0)
        return emit_copy(arena,in_verts,in_nv,in_faces,in_nf,out_verts,out_nv,out_faces,out_nf);

    L = (o.target_edge_len > 0.0f) ? (double)o.target_edge_len
                                   : sqrt(4.0*A / (1.7320508075688772*(double)F));
    if (!(L > 0.0))
        return emit_copy(arena,in_verts,in_nv,in_faces,in_nf,out_verts,out_nv,out_faces,out_nf);
    Lhi = (double)o.split_ratio * L;
    Llo = (double)o.collapse_ratio * L;
    maxlen = (double)o.max_collapse_len;
    if (Llo > maxlen) Llo = maxlen;   /* never collapse past the anti-fusion cap */
    if (verbose) fprintf(stderr,
        "[remesh] nv=%zu nf=%zu A=%.2f L=%.4f (split>%.4f collapse<%.4f)\n",
        in_nv, in_nf, A, L, Lhi, Llo);

    /* Face budget: hold the count in a soft band around F (= in_nf) so an
     * adaptively non-uniform input (coarse interior + frozen fine rim, whose
     * short edges can't collapse to offset the interior splits) can't inflate
     * past the decimation target. Split stops at hi, collapse stops at lo; both
     * sit inside the hard face_count_tol so the pass never falls back on drift. */
    hi = (size_t)((double)F * (1.0 + 0.8*(double)o.face_count_tol));
    lo = (size_t)((double)F * (1.0 - 0.8*(double)o.face_count_tol));
    if (lo < 4) lo = 4;
    if (hi < F) hi = F;

    /* Working copies (mutable, growable). */
    wnv=in_nv; wnf=in_nf;
    wv =(float*)  ARENA_ALLOC(arena,(long)(wnv*3*sizeof(float)));
    wf =(int32_t*)ARENA_ALLOC(arena,(long)(wnf*3*sizeof(int32_t)));
    wpin=(uint8_t*)ARENA_ALLOC(arena,(long)(wnv*sizeof(uint8_t)));
    memcpy(wv, in_verts, wnv*3*sizeof(float));
    memcpy(wf, in_faces, wnf*3*sizeof(int32_t));
    if (pin_mask) memcpy(wpin, pin_mask, wnv*sizeof(uint8_t));
    else memset(wpin, 0, wnv*sizeof(uint8_t));

    for (it=0; it<(size_t)o.n_iters; it++){
        size_t n_split=0, n_coll=0; int r;
        /* split long interior edges (capped at the budget ceiling) */
        if (o.do_split) for (r=0; r<o.split_max_rounds; r++){
            size_t k;
            if (wnf >= hi) break;
            k = split_round(arena, &wv,&wnv, &wf,&wnf, &wpin, Lhi, hi);
            n_split += k;
            if (k==0) break;
        }
        /* collapse short interior edges (capped at the budget floor) */
        if (o.do_collapse) for (r=0; r<o.collapse_max_rounds; r++){
            size_t k;
            if (wnf <= lo) break;
            k = collapse_round(arena, wv, wf,&wnf, wnv, wpin, Llo, maxlen, lo);
            n_coll += k;
            if (k==0) break;
        }
        /* flip for max-min-angle */
        if (o.do_flip) flip_rounds(arena, wv, wf, wnf, wnv, wpin, o.flip_max_rounds);
        /* tangential relax (+ optional MLS reproject) */
        if (o.do_relax && o.reproject!=REMESH_REPROJECT_NONE)
            relax_pass(arena, wv, wnv, wf, wnf, wpin, o.relax_lambda, (float)(0.5*L));
        if (o.reproject==REMESH_REPROJECT_MLS)
            mls_reproject(arena, wv, wnv, wpin, wf, wnf, (float)(2.0*L));
        /* reclaim orphaned verts from collapse */
        compact_verts(arena, &wv,&wnv, wf, wnf, &wpin);

        if (verbose) fprintf(stderr, "[remesh]  it=%zu +%zu split -%zu coll -> nv=%zu nf=%zu\n",
                             it, n_split, n_coll, wnv, wnf);
        if (n_split==0 && n_coll==0) break;   /* converged */
    }

    /* Fail-closed: manifold audit + face-count tolerance, else return input. */
    st = MeshManifold_audit(arena, wnv, wf, wnf);
    if (!MeshManifold_ok(&st)) {
        if (verbose) fprintf(stderr, "[remesh] FALLBACK non-manifold: nm_edges=%zu nm_verts=%zu (nf=%zu)\n",
                             st.nm_edges, st.nm_verts, wnf);
        return emit_copy(arena,in_verts,in_nv,in_faces,in_nf,out_verts,out_nv,out_faces,out_nf);
    }
    {
        double drift = fabs((double)wnf - (double)F);
        if (drift > (double)o.face_count_tol * (double)F) {
            if (verbose) fprintf(stderr, "[remesh] FALLBACK face drift: nf=%zu vs F=%zu (%.1f%%)\n",
                                 wnf, F, 100.0*drift/(double)F);
            return emit_copy(arena,in_verts,in_nv,in_faces,in_nf,out_verts,out_nv,out_faces,out_nf);
        }
    }
    if (verbose) fprintf(stderr, "[remesh] OK nf %zu -> %zu\n", F, wnf);

    *out_verts=wv; *out_nv=wnv; *out_faces=wf; *out_nf=wnf;
    return 0;
}
