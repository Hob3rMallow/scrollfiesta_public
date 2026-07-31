#include "rvd_clip.h"
#include "../common/kdtree.h"
#include "../common/attene_predicates_wrap.h"

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* Clip-polygon vertex: an AttenePt id + a tag for the edge that FOLLOWS it
 * (vertex i -> i+1). etag 0 = edge lies on a mesh edge with explicit endpoints
 * (me0,me1); etag 1 = edge lies on bisector(owner, bsj). The tag selects the exact
 * construction when the edge is cut. (fi,fj,fk) name the three sites when the
 * vertex itself is a facet-bisector2 Voronoi vertex (fi >= 0), else fi = -1. */
typedef struct {
    int32_t v; int etag; int32_t me0, me1, bsj;
    int32_t fi, fj, fk;
} PV;

#define RVD_MAX_POLY 256
#define RVD_MAX_CAND 1024

struct Rvd_T {
    Arena_T  arena;
    double  *verts;   /* [nv*3] */
    int32_t *faces;   /* [nf*3] */
    size_t   nv, nf;
    double   mesh_area;
};

static inline void vsub3(const double a[3], const double b[3], double o[3]) {
    o[0]=a[0]-b[0]; o[1]=a[1]-b[1]; o[2]=a[2]-b[2];
}
static inline void vcross3(const double a[3], const double b[3], double o[3]) {
    o[0]=a[1]*b[2]-a[2]*b[1]; o[1]=a[2]*b[0]-a[0]*b[2]; o[2]=a[0]*b[1]-a[1]*b[0];
}
static inline double vdot3(const double a[3], const double b[3]) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static inline double vlen3(const double a[3]) { return sqrt(vdot3(a,a)); }
static inline double tri_area(const double A[3], const double B[3], const double C[3]) {
    double ab[3], ac[3], cr[3];
    vsub3(B,A,ab); vsub3(C,A,ac); vcross3(ab,ac,cr);
    return 0.5*vlen3(cr);
}

Rvd_T Rvd_new(Arena_T arena, const double *verts, size_t nv,
              const int32_t *faces, size_t nf) {
    if (!arena || !verts || !faces || nv < 3 || nf < 1) return NULL;
    Rvd_T r = (Rvd_T)ARENA_ALLOC(arena, (size_t)sizeof(*r));
    r->arena = arena; r->nv = nv; r->nf = nf;
    r->verts = (double *)ARENA_ALLOC(arena, (size_t)(nv*3*sizeof(double)));
    memcpy(r->verts, verts, nv*3*sizeof(double));
    r->faces = (int32_t *)ARENA_ALLOC(arena, (size_t)(nf*3*sizeof(int32_t)));
    memcpy(r->faces, faces, nf*3*sizeof(int32_t));
    r->mesh_area = 0.0;
    for (size_t f=0; f<nf; f++) {
        const double *A=&verts[(size_t)faces[f*3+0]*3];
        const double *B=&verts[(size_t)faces[f*3+1]*3];
        const double *C=&verts[(size_t)faces[f*3+2]*3];
        r->mesh_area += tri_area(A,B,C);
    }
    return r;
}

double Rvd_mesh_area(Rvd_T r) { return r ? r->mesh_area : 0.0; }

/* Intersection of the edge following `cur` with bisector(si,sk); sets (fi,fj,fk)
 * to the three sites when the result is a facet-bisector2 Voronoi vertex, else -1. */
static int32_t make_cross(AttenePredCtx *ctx, const PV *cur, int32_t si, int32_t sk,
                          int32_t va, int32_t vb, int32_t vc,
                          int32_t *fi, int32_t *fj, int32_t *fk) {
    if (cur->etag == 0) {
        *fi = *fj = *fk = -1;
        return AttenePred_edge_bisector(ctx, cur->me0, cur->me1, si, sk);
    }
    *fi = si; *fj = cur->bsj; *fk = sk;
    return AttenePred_facet_bisector2(ctx, va, vb, vc, si, cur->bsj, sk);
}

/* Clip mesh triangle (A,B,C) to the Voronoi cell of owner `si`, against the
 * candidate sites. Returns the vertex count of the clipped convex polygon in
 * `poly` (>= 3), or a value < 3 if the cell is empty on this triangle. */
static int clip_cell(AttenePredCtx *ctx, const double *sites,
                     const int32_t *cand, size_t ncand, int32_t si,
                     const double A[3], const double B[3], const double C[3],
                     int32_t va, int32_t vb, int32_t vc, PV *poly) {
    const double *SI = &sites[(size_t)si*3];
    double da[3], db[3], dc[3];
    vsub3(A,SI,da); vsub3(B,SI,db); vsub3(C,SI,dc);
    double la=vlen3(da), lb=vlen3(db), lc=vlen3(dc);
    double Ri = la>lb ? (la>lc?la:lc) : (lb>lc?lb:lc);
    double cut2 = (2.0*Ri)*(2.0*Ri);

    PV tmp[RVD_MAX_POLY];
    int np = 3;
    poly[0] = (PV){ va,0, va,vb,-1, -1,-1,-1 };
    poly[1] = (PV){ vb,0, vb,vc,-1, -1,-1,-1 };
    poly[2] = (PV){ vc,0, vc,va,-1, -1,-1,-1 };

    for (size_t di = 0; di < ncand; di++) {
        int32_t sk = cand[di];
        if (sk == si) continue;
        const double *SK = &sites[(size_t)sk*3];
        double d[3]; vsub3(SK,SI,d);
        if (vdot3(d,d) > cut2) continue;

        int no = 0;
        for (int i = 0; i < np; i++) {
            PV cur = poly[i], nxt = poly[(i+1)%np];
            int cin = (AttenePred_side_bisector(ctx, cur.v, si, sk) <= 0);
            int nin = (AttenePred_side_bisector(ctx, nxt.v, si, sk) <= 0);
            if (cin) {
                if (no >= RVD_MAX_POLY) return -1;
                tmp[no++] = cur;
                if (!nin) {
                    int32_t fi, fj, fk;
                    int32_t X = make_cross(ctx, &cur, si, sk, va, vb, vc, &fi, &fj, &fk);
                    if (no >= RVD_MAX_POLY) return -1;
                    tmp[no++] = (PV){ X,1, -1,-1,sk, fi,fj,fk };   /* follow: clip line */
                }
            } else if (nin) {
                int32_t fi, fj, fk;
                int32_t X = make_cross(ctx, &cur, si, sk, va, vb, vc, &fi, &fj, &fk);
                if (no >= RVD_MAX_POLY) return -1;
                tmp[no++] = (PV){ X, cur.etag, cur.me0, cur.me1, cur.bsj, fi,fj,fk }; /* resume */
            }
        }
        memcpy(poly, tmp, (size_t)no*sizeof(PV));
        np = no;
        if (np < 3) return np;
    }
    return np;
}

/* Per-triangle candidate set: a provable superset of owners + cutting neighbours.
 * Returns the count written into cand. */
static size_t triangle_candidates(KDTree_T kt, const double A[3], const double B[3],
                                  const double C[3], int32_t *cand) {
    double cz[3] = { (A[0]+B[0]+C[0])/3.0, (A[1]+B[1]+C[1])/3.0, (A[2]+B[2]+C[2])/3.0 };
    double e0[3], e1[3], e2[3];
    vsub3(B,A,e0); vsub3(C,B,e1); vsub3(A,C,e2);
    double le0=vlen3(e0), le1=vlen3(e1), le2=vlen3(e2);
    double maxe = le0>le1 ? (le0>le2?le0:le2) : (le1>le2?le1:le2);
    float czf[3] = { (float)cz[0], (float)cz[1], (float)cz[2] };
    float dnn_sq = 0.0f;
    (void)KDTree_nearest(kt, czf, &dnn_sq);
    double d_near = sqrt((double)dnn_sq);
    double R = 3.0*d_near + 8.0*maxe;
    R = R*1.01 + 1.0;                        /* float-kdtree safety */
    return KDTree_ball_query(kt, czf, (float)(R*R), cand, RVD_MAX_CAND);
}

/* qsort helper: compare a "sorted triple" key stored as int32[6] = {key3, tri3}. */
static int cmp_triple(const void *pa, const void *pb) {
    const int32_t *a = (const int32_t *)pa, *b = (const int32_t *)pb;
    for (int i = 0; i < 3; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

/* BFS owner-queue push: enqueue site s for the current triangle if it is unvisited
 * (per-triangle generation stamp) and the queue has room. */
static inline void owner_push(int32_t s, int32_t cur_gen, int32_t *gen,
                              int32_t *q, int *qt) {
    if (s >= 0 && *qt < RVD_MAX_CAND && gen[s] != cur_gen) {
        gen[s] = cur_gen;
        q[(*qt)++] = s;
    }
}

int Rvd_accumulate(Rvd_T r, const double *sites, size_t nsites,
                   double *out_area, double *out_centroid,
                   double *out_NA, double *out_bNA,
                   int32_t **out_dual_tris, size_t *out_ndual,
                   const CvtField *dens) {
    if (!r || !sites || nsites == 0 || !out_area || !out_centroid) return -1;
    if (out_dual_tris) *out_dual_tris = NULL;
    if (out_ndual) *out_ndual = 0;
    int want_moments = (out_NA != NULL && out_bNA != NULL);

    for (size_t i = 0; i < nsites; i++) {
        out_area[i] = 0.0;
        out_centroid[i*3+0] = out_centroid[i*3+1] = out_centroid[i*3+2] = 0.0;
    }
    if (want_moments) {
        memset(out_NA,  0, nsites*6*sizeof(double));
        memset(out_bNA, 0, nsites*3*sizeof(double));
    }

    /* Raw dual triples ({sortedkey[3], tri[3]}) persist past the scratch restore. */
    size_t dual_cap = out_dual_tris ? (nsites * 8 + 16) : 0;
    int32_t *raw = NULL;
    size_t   nraw = 0;
    if (out_dual_tris)
        raw = (int32_t *)ARENA_ALLOC(r->arena, (size_t)(dual_cap * 6 * sizeof(int32_t)));

    Arena_Mark mark = Arena_save(r->arena);

    float *sf = (float *)ARENA_ALLOC(r->arena, (size_t)(nsites*3*sizeof(float)));
    for (size_t i = 0; i < nsites; i++) {
        sf[i*3+0]=(float)sites[i*3+0]; sf[i*3+1]=(float)sites[i*3+1]; sf[i*3+2]=(float)sites[i*3+2];
    }
    KDTree_T kt = KDTree_new(r->arena, sf, nsites);

    AttenePredCtx *ctx = AttenePred_new();
    for (size_t i = 0; i < nsites; i++)
        (void)AttenePred_explicit(ctx, sites[i*3+0], sites[i*3+1], sites[i*3+2]);

    int32_t *cand = (int32_t *)ARENA_ALLOC(r->arena, (size_t)(RVD_MAX_CAND*sizeof(int32_t)));
    PV poly[RVD_MAX_POLY];
    /* Propagation scratch: per-triangle visited marker (generation stamp, never reset)
     * + owner queue. Owners are a connected subset of `cand`, bounded by RVD_MAX_CAND. */
    int32_t *owner_gen = (int32_t *)ARENA_CALLOC(r->arena, (size_t)nsites, (size_t)sizeof(int32_t));
    int32_t *owner_q   = (int32_t *)ARENA_ALLOC(r->arena, (size_t)(RVD_MAX_CAND*sizeof(int32_t)));
    int32_t cur_gen = 0;

    for (size_t f = 0; f < r->nf; f++) {
        const double *A=&r->verts[(size_t)r->faces[f*3+0]*3];
        const double *B=&r->verts[(size_t)r->faces[f*3+1]*3];
        const double *C=&r->verts[(size_t)r->faces[f*3+2]*3];

        double nrm[3]; { double ab[3],ac[3]; vsub3(B,A,ab); vsub3(C,A,ac); vcross3(ab,ac,nrm); }
        double un[3]={nrm[0],nrm[1],nrm[2]}, nnT[6]={0,0,0,0,0,0};
        if (want_moments) {
            double nl=sqrt(un[0]*un[0]+un[1]*un[1]+un[2]*un[2]);
            if (nl>0.0){ un[0]/=nl; un[1]/=nl; un[2]/=nl; }
            nnT[0]=un[0]*un[0]; nnT[1]=un[0]*un[1]; nnT[2]=un[0]*un[2];
            nnT[3]=un[1]*un[1]; nnT[4]=un[1]*un[2]; nnT[5]=un[2]*un[2];
        }

        /* Nearest site of each vertex. If all three agree, the whole convex planar
         * triangle lies in that Voronoi cell -- accumulate it directly, no exact clip
         * (fast path). Otherwise these three sites seed the propagation BFS below. */
        float af[3]={(float)A[0],(float)A[1],(float)A[2]};
        float bf[3]={(float)B[0],(float)B[1],(float)B[2]};
        float cf[3]={(float)C[0],(float)C[1],(float)C[2]};
        float d0,d1,d2;
        size_t na=KDTree_nearest(kt,af,&d0), nb=KDTree_nearest(kt,bf,&d1), nc2=KDTree_nearest(kt,cf,&d2);
        if (na==nb && nb==nc2) {
            int32_t si=(int32_t)na;
            double at=tri_area(A,B,C);
            double ct0=(A[0]+B[0]+C[0])/3.0, ct1=(A[1]+B[1]+C[1])/3.0, ct2=(A[2]+B[2]+C[2])/3.0;
            double w = dens ? at*Rvd_density_rho(dens,ct0,ct1,ct2) : at;  /* NULL => w==at (bit-identical) */
            out_area[si]+=w;
            out_centroid[si*3+0]+=w*ct0; out_centroid[si*3+1]+=w*ct1; out_centroid[si*3+2]+=w*ct2;
            if (want_moments) {
                double *Ni=&out_NA[(size_t)si*6];
                Ni[0]+=w*nnT[0]; Ni[1]+=w*nnT[1]; Ni[2]+=w*nnT[2];
                Ni[3]+=w*nnT[3]; Ni[4]+=w*nnT[4]; Ni[5]+=w*nnT[5];
                double ndotc=un[0]*ct0+un[1]*ct1+un[2]*ct2;
                double *bi=&out_bNA[(size_t)si*3];
                bi[0]+=w*ndotc*un[0]; bi[1]+=w*ndotc*un[1]; bi[2]+=w*ndotc*un[2];
            }
            continue;
        }

        int32_t va=AttenePred_explicit(ctx, A[0],A[1],A[2]);
        int32_t vb=AttenePred_explicit(ctx, B[0],B[1],B[2]);
        int32_t vc=AttenePred_explicit(ctx, C[0],C[1],C[2]);

        size_t ncand = triangle_candidates(kt, A, B, C, cand);
        if (ncand >= RVD_MAX_CAND)
            fprintf(stderr, "  rvd: face %zu saturated candidates (%zu)\n", f, ncand);

        /* Propagation: clip only ACTUAL owners. The owners of a convex triangle are
         * connected in the restricted-Voronoi adjacency graph, so a BFS seeded at the
         * three vertices' nearest sites reaches them all; each clipped cell's bisector
         * edges (etag==1) name the neighbouring site (bsj) to enqueue. `cand` stays the
         * membership superset clip_cell tests against. The queue is sized RVD_MAX_CAND
         * (>= ncand), so it never overflows -- worst case degrades to the old bound. */
        cur_gen++;
        int qh = 0, qt = 0;
        owner_push((int32_t)na,  cur_gen, owner_gen, owner_q, &qt);
        owner_push((int32_t)nb,  cur_gen, owner_gen, owner_q, &qt);
        owner_push((int32_t)nc2, cur_gen, owner_gen, owner_q, &qt);
        while (qh < qt) {
            int32_t si = owner_q[qh++];
            int np = clip_cell(ctx, sites, cand, ncand, si, A, B, C, va, vb, vc, poly);
            if (np < 3) continue;
            for (int e = 0; e < np; e++)
                if (poly[e].etag == 1)
                    owner_push(poly[e].bsj, cur_gen, owner_gen, owner_q, &qt);

            double P0[3];
            if (AttenePred_xyz(ctx, poly[0].v, &P0[0],&P0[1],&P0[2]) != 0) continue;
            for (int i = 1; i+1 < np; i++) {
                double Pi[3], Pj[3];
                if (AttenePred_xyz(ctx, poly[i].v,   &Pi[0],&Pi[1],&Pi[2]) != 0) continue;
                if (AttenePred_xyz(ctx, poly[i+1].v, &Pj[0],&Pj[1],&Pj[2]) != 0) continue;
                double at = tri_area(P0,Pi,Pj);
                double ct0=(P0[0]+Pi[0]+Pj[0])/3.0, ct1=(P0[1]+Pi[1]+Pj[1])/3.0, ct2=(P0[2]+Pi[2]+Pj[2])/3.0;
                double w = dens ? at*Rvd_density_rho(dens,ct0,ct1,ct2) : at;  /* NULL => w==at (bit-identical) */
                out_area[si] += w;
                out_centroid[si*3+0] += w*ct0;
                out_centroid[si*3+1] += w*ct1;
                out_centroid[si*3+2] += w*ct2;
                if (want_moments) {
                    double *Ni=&out_NA[(size_t)si*6];
                    Ni[0]+=w*nnT[0]; Ni[1]+=w*nnT[1]; Ni[2]+=w*nnT[2];
                    Ni[3]+=w*nnT[3]; Ni[4]+=w*nnT[4]; Ni[5]+=w*nnT[5];
                    double ndotc=un[0]*ct0+un[1]*ct1+un[2]*ct2;
                    double *bi=&out_bNA[(size_t)si*3];
                    bi[0]+=w*ndotc*un[0]; bi[1]+=w*ndotc*un[1]; bi[2]+=w*ndotc*un[2];
                }
            }

            /* Dual: each facet-bisector2 corner is a Voronoi vertex of 3 sites. */
            if (out_dual_tris) {
                for (int i = 0; i < np; i++) {
                    if (poly[i].fi < 0) continue;
                    int32_t t[3] = { poly[i].fi, poly[i].fj, poly[i].fk };
                    /* orient by local surface normal */
                    const double *Pa=&sites[(size_t)t[0]*3], *Pb=&sites[(size_t)t[1]*3], *Pc=&sites[(size_t)t[2]*3];
                    double u[3], v[3], cr[3];
                    vsub3(Pb,Pa,u); vsub3(Pc,Pa,v); vcross3(u,v,cr);
                    if (vdot3(cr,nrm) < 0.0) { int32_t tt=t[1]; t[1]=t[2]; t[2]=tt; }
                    int32_t key[3] = { t[0],t[1],t[2] };
                    for (int p=0;p<2;p++) for (int q=0;q<2-p;q++)
                        if (key[q]>key[q+1]) { int32_t s=key[q]; key[q]=key[q+1]; key[q+1]=s; }
                    if (nraw < dual_cap) {
                        int32_t *rec = &raw[nraw*6];
                        rec[0]=key[0]; rec[1]=key[1]; rec[2]=key[2];
                        rec[3]=t[0];   rec[4]=t[1];   rec[5]=t[2];
                        nraw++;
                    }
                }
            }
        }
        AttenePred_truncate(ctx, (int32_t)nsites);
    }

    AttenePred_free(ctx);
    Arena_restore(r->arena, mark);

    for (size_t i = 0; i < nsites; i++) {
        if (out_area[i] > 0.0) {
            out_centroid[i*3+0] /= out_area[i];
            out_centroid[i*3+1] /= out_area[i];
            out_centroid[i*3+2] /= out_area[i];
        } else {
            out_centroid[i*3+0]=sites[i*3+0]; out_centroid[i*3+1]=sites[i*3+1]; out_centroid[i*3+2]=sites[i*3+2];
        }
    }

    /* Dedup the dual triples (each Voronoi vertex is seen once per incident cell). */
    if (out_dual_tris && nraw > 0) {
        qsort(raw, nraw, 6*sizeof(int32_t), cmp_triple);
        int32_t *tris = (int32_t *)ARENA_ALLOC(r->arena, (size_t)(nraw*3*sizeof(int32_t)));
        size_t nt = 0;
        for (size_t i = 0; i < nraw; i++) {
            if (i > 0 && cmp_triple(&raw[i*6], &raw[(i-1)*6]) == 0) continue;
            /* skip degenerate (repeated site index) */
            int32_t a=raw[i*6+3], b=raw[i*6+4], c=raw[i*6+5];
            if (a==b || b==c || a==c) continue;
            tris[nt*3+0]=a; tris[nt*3+1]=b; tris[nt*3+2]=c; nt++;
        }
        *out_dual_tris = tris;
        *out_ndual = nt;
    }

    return 0;
}
