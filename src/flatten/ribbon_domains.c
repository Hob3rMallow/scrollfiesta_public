#include "ribbon_domains.h"
#include "../common/union_find.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct { size_t vertices,faces,components; } RdSize;
struct RibbonDomains_T {
    const MeshBinData *source;
    RibbonDomainsReport report;
    RibbonDomainPart *part;
    int32_t *vertex_owner,*face_owner,*local_vertex;
    int32_t *vertices,*faces;
};

static int rd_input_bytes(size_t vertices,size_t faces,size_t *out)
{
    size_t xyz=0,uv=0,triangles=0;
    if (!out) return -1;
    *out=0;
    if (vertices>(SIZE_MAX-15)/(3*sizeof(float)) ||
        faces>(SIZE_MAX-15)/(3*sizeof(int32_t))) return -1;
    xyz=(vertices*3*sizeof(float)+15)&~(size_t)15;
    uv=(vertices*2*sizeof(float)+15)&~(size_t)15;
    triangles=(faces*3*sizeof(int32_t)+15)&~(size_t)15;
    if (uv>SIZE_MAX-xyz || triangles>SIZE_MAX-xyz-uv) return -1;
    *out=xyz+uv+triangles;
    return 0;
}

int RibbonDomains_new(Arena_T arena,const MeshBinData *source,
    size_t target_faces,size_t target_vertices,RibbonDomains_T *out)
{
    Arena_T scratch=NULL;
    UnionFind uf={0};
    RibbonDomains_T d=NULL;
    size_t *nv=NULL,*nf=NULL,*voffset=NULL,*foffset=NULL;
    int32_t *seed=NULL,*batch=NULL;
    RdSize *sizes=NULL;
    size_t ncomponents=0,nparts=0;
    int rc=-1;
    if (!out) return -1;
    *out=NULL;
    if (!arena || !source || !target_faces || !target_vertices ||
        source->nv>INT32_MAX || source->nf>INT32_MAX ||
        (source->nv && (!source->verts || !source->uv)) ||
        (source->nf && !source->faces)) return -1;
    scratch=Arena_new();
    uf=UF_new(scratch,(int32_t)source->nv);
    seed=(int32_t *)ARENA_ALLOC(scratch,source->nv*sizeof *seed);
    batch=(int32_t *)ARENA_ALLOC(scratch,source->nv*sizeof *batch);
    nv=(size_t *)ARENA_CALLOC(scratch,source->nv,sizeof *nv);
    nf=(size_t *)ARENA_CALLOC(scratch,source->nv,sizeof *nf);
    for (size_t i=0;i<source->nv;i++) {
        seed[i]=-1; batch[i]=-1;
        for (int k=0;k<3;k++) if (!isfinite(source->verts[3*i+k])) goto done;
        for (int k=0;k<2;k++) if (!isfinite(source->uv[2*i+k])) goto done;
    }
    for (size_t f=0;f<source->nf;f++) {
        const int32_t *t=source->faces+3*f;
        for (int k=0;k<3;k++) if (t[k]<0 || (size_t)t[k]>=source->nv) goto done;
        uf_union(&uf,t[0],t[1]); uf_union(&uf,t[0],t[2]);
    }
    for (size_t i=0;i<source->nv;i++) {
        int32_t root=uf_find(&uf,(int32_t)i);
        if (seed[root]<0) { seed[root]=(int32_t)i; ncomponents++; }
        nv[seed[root]]++;
    }
    for (size_t f=0;f<source->nf;f++) nf[seed[uf_find(&uf,source->faces[3*f])]]++;
    sizes=(RdSize *)ARENA_CALLOC(scratch,ncomponents,sizeof *sizes);
    d=(RibbonDomains_T)ARENA_CALLOC(arena,1,sizeof *d);
    d->source=source; d->report.components=ncomponents;
    for (size_t i=0;i<source->nv;i++) if (nv[i]) {
        RdSize *s=nparts ? sizes+nparts-1 : NULL;
        if (!s || s->faces>target_faces || s->vertices>target_vertices ||
            nf[i]>target_faces-s->faces || nv[i]>target_vertices-s->vertices) {
            s=sizes+nparts++;
        }
        s->vertices+=nv[i]; s->faces+=nf[i]; s->components++;
        batch[i]=(int32_t)(nparts-1);
        if (nv[i]>d->report.largest_component_vertices) d->report.largest_component_vertices=nv[i];
        if (nf[i]>d->report.largest_component_faces) d->report.largest_component_faces=nf[i];
    }
    d->report.parts=nparts;
    d->part=(RibbonDomainPart *)ARENA_CALLOC(arena,nparts,sizeof *d->part);
    d->vertex_owner=(int32_t *)ARENA_ALLOC(arena,source->nv*sizeof(int32_t));
    d->face_owner=(int32_t *)ARENA_ALLOC(arena,source->nf*sizeof(int32_t));
    d->local_vertex=(int32_t *)ARENA_ALLOC(arena,source->nv*sizeof(int32_t));
    d->vertices=(int32_t *)ARENA_ALLOC(arena,source->nv*sizeof(int32_t));
    d->faces=(int32_t *)ARENA_ALLOC(arena,source->nf*sizeof(int32_t));
    voffset=(size_t *)ARENA_CALLOC(scratch,nparts+1,sizeof *voffset);
    foffset=(size_t *)ARENA_CALLOC(scratch,nparts+1,sizeof *foffset);
    for (size_t p=0;p<nparts;p++) {
        RibbonDomainPart *part=d->part+p;
        part->vertices=sizes[p].vertices; part->faces=sizes[p].faces;
        if (rd_input_bytes(part->vertices,part->faces,&part->input_bytes)!=0) goto done;
        part->components=sizes[p].components;
        part->oversized=part->vertices>target_vertices || part->faces>target_faces;
        d->report.oversized_parts+=(size_t)part->oversized;
        if (part->vertices>d->report.largest_part_vertices) d->report.largest_part_vertices=part->vertices;
        if (part->faces>d->report.largest_part_faces) d->report.largest_part_faces=part->faces;
        voffset[p+1]=voffset[p]+part->vertices; foffset[p+1]=foffset[p]+part->faces;
        part->source_vertices=d->vertices+voffset[p]; part->source_faces=d->faces+foffset[p];
    }
    if (voffset[nparts]!=source->nv || foffset[nparts]!=source->nf) goto done;
    for (size_t i=0;i<source->nv;i++) {
        int32_t owner=batch[seed[uf_find(&uf,(int32_t)i)]];
        size_t at=voffset[owner]++;
        d->vertex_owner[i]=owner;
        d->local_vertex[i]=(int32_t)(at-(size_t)(d->part[owner].source_vertices-d->vertices));
        d->vertices[at]=(int32_t)i;
    }
    for (size_t f=0;f<source->nf;f++) {
        int32_t owner=d->vertex_owner[source->faces[3*f]];
        for (int k=1;k<3;k++) if (d->vertex_owner[source->faces[3*f+k]]!=owner) goto done;
        d->face_owner[f]=owner; d->faces[foffset[owner]++]=(int32_t)f;
    }
    *out=d; rc=0;
done:
    Arena_dispose(&scratch); return rc;
}

int RibbonDomains_report(RibbonDomains_T d,RibbonDomainsReport *out)
{
    if (!d || !out) return -1;
    *out=d->report; return 0;
}
int RibbonDomains_part(RibbonDomains_T d,size_t p,RibbonDomainPart *out)
{
    if (!d || !out || p>=d->report.parts) return -1;
    *out=d->part[p]; return 0;
}
const int32_t *RibbonDomains_vertex_owners(RibbonDomains_T d)
{ return d ? d->vertex_owner : NULL; }
const int32_t *RibbonDomains_face_owners(RibbonDomains_T d)
{ return d ? d->face_owner : NULL; }

int RibbonDomains_extract(Arena_T scratch,RibbonDomains_T d,size_t p,
    size_t maximum_bytes,MeshBinData *out)
{
    const RibbonDomainPart *part=NULL;
    const MeshBinData *source=NULL;
    if (!out) return -1;
    memset(out,0,sizeof *out);
    if (!scratch || !d || p>=d->report.parts) return -1;
    part=d->part+p; source=d->source;
    if (part->input_bytes>maximum_bytes) return -1;
    out->verts=(float *)ARENA_ALLOC(scratch,3*part->vertices*sizeof(float));
    out->uv=(float *)ARENA_ALLOC(scratch,2*part->vertices*sizeof(float));
    out->faces=(int32_t *)ARENA_ALLOC(scratch,3*part->faces*sizeof(int32_t));
    out->nv=part->vertices; out->nf=part->faces;
    for (size_t i=0;i<part->vertices;i++) {
        size_t v=(size_t)part->source_vertices[i];
        memcpy(out->verts+3*i,source->verts+3*v,3*sizeof(float));
        memcpy(out->uv+2*i,source->uv+2*v,2*sizeof(float));
    }
    for (size_t f=0;f<part->faces;f++) for (int k=0;k<3;k++)
        out->faces[3*f+k]=d->local_vertex[source->faces[3*(size_t)part->source_faces[f]+k]];
    return 0;
}

int RibbonDomains_selftest(void)
{
    Arena_T arena=Arena_new();
    float xyz[36]={0},uv[24]={0};
    int32_t faces[15]={0,1,2,0,2,3,4,5,6,4,6,7,8,9,10};
    MeshBinData source={xyz,uv,faces,12,5};
    RibbonDomains_T d=NULL,over=NULL,empty=NULL;
    RibbonDomainsReport r={0};
    uint8_t seen_v[12]={0},seen_f[5]={0};
    int fail=0;
    {
        size_t bytes=0;
        /* Real 10x10x10 component: size arithmetic, without a giant fixture. */
        fail+=rd_input_bytes(4720858,9364994,&bytes)!=0 || bytes!=206797104;
        fail+=rd_input_bytes(9716848,19255933,&bytes)!=0 || bytes!=425408160;
        fail+=rd_input_bytes(0,0,&bytes)!=0 || bytes!=0;
        fail+=rd_input_bytes(SIZE_MAX,1,&bytes)==0;
        fail+=rd_input_bytes(1,SIZE_MAX,&bytes)==0;
    }
    for (size_t i=0;i<12;i++) {
        /* The first two disconnected squares deliberately coincide in XYZ/UV. */
        xyz[3*i]=(float)(i%4); xyz[3*i+1]=(float)(i%2);
        uv[2*i]=(float)(i%4); uv[2*i+1]=(float)(i%2);
    }
    if (RibbonDomains_new(arena,&source,2,5,&d)!=0 || RibbonDomains_report(d,&r)!=0) { fail++; goto done; }
    fail+=r.components!=4 || r.parts!=3 || r.oversized_parts!=0;
    for (size_t p=0;p<r.parts;p++) {
        Arena_T scratch=Arena_new();
        RibbonDomainPart part={0};
        MeshBinData m={0};
        if (RibbonDomains_part(d,p,&part)!=0 || RibbonDomains_extract(scratch,d,p,1024,&m)!=0) fail++;
        else {
            for (size_t i=0;i<m.nv;i++) {
                size_t v=(size_t)part.source_vertices[i];
                fail+=seen_v[v]++!=0 || RibbonDomains_vertex_owners(d)[v]!=(int32_t)p;
                fail+=memcmp(m.verts+3*i,xyz+3*v,3*sizeof(float))!=0 || memcmp(m.uv+2*i,uv+2*v,2*sizeof(float))!=0;
            }
            for (size_t f=0;f<m.nf;f++) {
                size_t original=(size_t)part.source_faces[f];
                fail+=seen_f[original]++!=0 || RibbonDomains_face_owners(d)[original]!=(int32_t)p;
                for (int k=0;k<3;k++) fail+=part.source_vertices[m.faces[3*f+k]]!=faces[3*original+k];
            }
            fail+=RibbonDomains_extract(scratch,d,p,part.input_bytes-1,&m)==0;
            fail+=m.verts!=NULL || m.uv!=NULL || m.faces!=NULL || m.nv!=0 || m.nf!=0;
            fail+=RibbonDomains_extract(scratch,d,p,part.input_bytes,&m)!=0;
        }
        Arena_dispose(&scratch);
    }
    for (size_t i=0;i<12;i++) fail+=seen_v[i]!=1;
    for (size_t i=0;i<5;i++) fail+=seen_f[i]!=1;
    fail+=RibbonDomains_new(arena,&source,1,3,&over)!=0;
    fail+=RibbonDomains_report(over,&r)!=0 || r.parts!=4 || r.oversized_parts!=2;
    {
        MeshBinData m={0};
        fail+=RibbonDomains_new(arena,&m,1,1,&empty)!=0;
        fail+=RibbonDomains_report(empty,&r)!=0 || r.parts || r.components;
    }
    faces[0]=-1;
    fail+=RibbonDomains_new(arena,&source,2,5,&empty)==0;
    faces[0]=0; uv[0]=NAN;
    fail+=RibbonDomains_new(arena,&source,2,5,&empty)==0;
done:
    fprintf(stderr,"[selftest] ribbon whole-domain batching %s (%d failures)\n",fail ? "FAIL" : "PASS",fail);
    Arena_dispose(&arena); return fail ? -1 : 0;
}
