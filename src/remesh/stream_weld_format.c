#include "stream_weld_format.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SW_SHARD_VERSION 1u
#define SW_SHARD_LINEAGE 1u
#define SW_PATCH_VERSION 3u
#define SW_PATCH_HAS_CONFLICTS 1u

static int sw_mul_ok(size_t a, size_t b, size_t *out)
{
    if (a != 0 && b > SIZE_MAX / a) return 0;
    *out = a * b;
    return 1;
}

static int sw_write(FILE *f, const void *p, size_t size, size_t n)
{
    return n == 0 || fwrite(p, size, n, f) == n;
}

static int sw_read(FILE *f, void *p, size_t size, size_t n)
{
    return n == 0 || fread(p, size, n, f) == n;
}

int StreamWeld_write_shard(const char *path, const char *cube_id,
                           const float *verts, size_t nv,
                           const int32_t *faces, size_t nf,
                           const int32_t *vertex_chart, size_t n_charts,
                           const int32_t *parent0, const int32_t *parent1)
{
    static const char magic[8] = {'G','W','S','H','R','D','1','\n'};
    FILE *f = NULL;
    uint32_t version = SW_SHARD_VERSION, flags = SW_SHARD_LINEAGE;
    uint64_t nv64 = (uint64_t)nv, nf64 = (uint64_t)nf;
    uint64_t nc64 = (uint64_t)n_charts;
    char id[48] = {0};
    int ok = 0;
    if (!path || !cube_id || !verts || !faces || !vertex_chart ||
        !parent0 || !parent1 || strlen(cube_id) >= sizeof(id) ||
        nv > (size_t)INT32_MAX || n_charts > (size_t)INT32_MAX)
        return -1;
    memcpy(id, cube_id, strlen(cube_id));
    f = fopen(path, "wb");
    if (!f) return -1;
    ok = sw_write(f, magic, 1, sizeof(magic)) &&
         sw_write(f, &version, sizeof(version), 1) &&
         sw_write(f, &flags, sizeof(flags), 1) &&
         sw_write(f, id, 1, sizeof(id)) &&
         sw_write(f, &nv64, sizeof(nv64), 1) &&
         sw_write(f, &nf64, sizeof(nf64), 1) &&
         sw_write(f, &nc64, sizeof(nc64), 1) &&
         sw_write(f, verts, sizeof(*verts), nv * 3) &&
         sw_write(f, faces, sizeof(*faces), nf * 3) &&
         sw_write(f, vertex_chart, sizeof(*vertex_chart), nv) &&
         sw_write(f, parent0, sizeof(*parent0), nv) &&
         sw_write(f, parent1, sizeof(*parent1), nv) &&
         fflush(f) == 0;
    if (fclose(f) != 0) ok = 0;
    return ok ? 0 : -1;
}

int StreamWeld_read_shard(Arena_T arena, const char *path,
                          StreamWeldShard *out)
{
    static const char expected[8] = {'G','W','S','H','R','D','1','\n'};
    char magic[8], id[48];
    uint32_t version = 0, flags = 0;
    uint64_t nv64 = 0, nf64 = 0, nc64 = 0;
    size_t nv = 0, nf = 0, n3 = 0;
    FILE *f = NULL;
    int ok = 0, tail = 0;
    if (!arena || !path || !out) return -1;
    memset(out, 0, sizeof(*out));
    f = fopen(path, "rb");
    if (!f) return -1;
    if (!sw_read(f, magic, 1, sizeof(magic)) ||
        memcmp(magic, expected, sizeof(magic)) != 0 ||
        !sw_read(f, &version, sizeof(version), 1) ||
        !sw_read(f, &flags, sizeof(flags), 1) ||
        version != SW_SHARD_VERSION || !(flags & SW_SHARD_LINEAGE) ||
        !sw_read(f, id, 1, sizeof(id)) || memchr(id, '\0', sizeof(id)) == NULL ||
        !sw_read(f, &nv64, sizeof(nv64), 1) ||
        !sw_read(f, &nf64, sizeof(nf64), 1) ||
        !sw_read(f, &nc64, sizeof(nc64), 1) ||
        nv64 == 0 || nf64 == 0 || nv64 > INT32_MAX || nc64 > INT32_MAX ||
        nv64 > SIZE_MAX || nf64 > SIZE_MAX)
        goto done;
    nv = (size_t)nv64; nf = (size_t)nf64;
    if (!sw_mul_ok(nv, 3, &n3)) goto done;
    memcpy(out->cube_id, id, sizeof(id));
    out->nv = nv; out->nf = nf; out->n_charts = (size_t)nc64;
    out->verts = (float *)ARENA_ALLOC(arena, n3 * sizeof(*out->verts));
    if (!sw_mul_ok(nf, 3, &n3)) goto done;
    out->faces = (int32_t *)ARENA_ALLOC(arena, n3 * sizeof(*out->faces));
    out->vertex_chart = (int32_t *)ARENA_ALLOC(
        arena, nv * sizeof(*out->vertex_chart));
    out->parent0 = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(*out->parent0));
    out->parent1 = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(*out->parent1));
    ok = sw_read(f, out->verts, sizeof(*out->verts), nv * 3) &&
         sw_read(f, out->faces, sizeof(*out->faces), nf * 3) &&
         sw_read(f, out->vertex_chart, sizeof(*out->vertex_chart), nv) &&
         sw_read(f, out->parent0, sizeof(*out->parent0), nv) &&
         sw_read(f, out->parent1, sizeof(*out->parent1), nv);
    if (!ok) goto done;
    tail = fgetc(f);
    ok = tail == EOF;
done:
    fclose(f);
    if (!ok) memset(out, 0, sizeof(*out));
    return ok ? 0 : -1;
}

int StreamWeld_write_patch(const char *path,
                           const StreamWeldShard *a,
                           const StreamWeldShard *b,
                           const int32_t *faces, size_t nf,
                           const ChartZipperTransaction *tx, size_t ntx,
                           const StreamWeldConflict *conflicts,
                           size_t nconflicts)
{
    static const char magic[8] = {'G','W','P','A','T','C','3','\n'};
    FILE *f = NULL;
    uint32_t version = SW_PATCH_VERSION, flags = SW_PATCH_HAS_CONFLICTS;
    uint64_t nva, nvb, nca, ncb, ntx64, total_faces = 0, nconflicts64;
    char ida[48] = {0}, idb[48] = {0};
    int ok = 0;
    if (!path || !a || !b || !faces || (ntx && !tx) ||
        (nconflicts && !conflicts) || ntx > UINT32_MAX ||
        strlen(a->cube_id) >= sizeof(ida) ||
        strlen(b->cube_id) >= sizeof(idb)) return -1;
    memcpy(ida, a->cube_id, strlen(a->cube_id));
    memcpy(idb, b->cube_id, strlen(b->cube_id));
    nva=(uint64_t)a->nv;nvb=(uint64_t)b->nv;
    nca=(uint64_t)a->n_charts;ncb=(uint64_t)b->n_charts;
    ntx64=(uint64_t)ntx;nconflicts64=(uint64_t)nconflicts;
    for(size_t t=0;t<ntx;t++){
        if(tx[t].face_first > nf || tx[t].face_count > nf-tx[t].face_first)
            return -1;
        if(tx[t].face_count > UINT64_MAX-total_faces) return -1;
        total_faces += (uint64_t)tx[t].face_count;
    }
    for(size_t i=0;i<nconflicts;i++){
        uint32_t ta=conflicts[i].transaction_a;
        uint32_t tb=conflicts[i].transaction_b;
        if(ta>=tb||tb>=ntx||conflicts[i].reason_mask==0)return -1;
    }
    f=fopen(path,"wb");if(!f)return -1;
    ok=sw_write(f,magic,1,sizeof(magic))&&
       sw_write(f,&version,sizeof(version),1)&&
       sw_write(f,&flags,sizeof(flags),1)&&
       sw_write(f,ida,1,sizeof(ida))&&sw_write(f,idb,1,sizeof(idb))&&
       sw_write(f,&nva,sizeof(nva),1)&&sw_write(f,&nvb,sizeof(nvb),1)&&
       sw_write(f,&nca,sizeof(nca),1)&&sw_write(f,&ncb,sizeof(ncb),1)&&
       sw_write(f,&ntx64,sizeof(ntx64),1)&&
       sw_write(f,&total_faces,sizeof(total_faces),1);
    for(size_t t=0;ok&&t<ntx;t++){
        int32_t ca=tx[t].chart_a,cb=tx[t].chart_b;
        int32_t la=-1,lb=-1;
        int8_t oa=tx[t].orient_a,ob=tx[t].orient_b;
        uint16_t transaction_flags=tx[t].flags;
        uint64_t count=(uint64_t)tx[t].face_count;
        double rank[5]={tx[t].full_coverage,tx[t].span_coverage,
                        tx[t].mean_gap,tx[t].support,tx[t].score};
        if(ca>=0&&(size_t)ca<a->n_charts&&cb>=(int32_t)a->n_charts&&
           (size_t)(cb-(int32_t)a->n_charts)<b->n_charts){
            la=ca;lb=cb-(int32_t)a->n_charts;
        }else if(cb>=0&&(size_t)cb<a->n_charts&&
                 ca>=(int32_t)a->n_charts&&
                 (size_t)(ca-(int32_t)a->n_charts)<b->n_charts){
            la=cb;lb=ca-(int32_t)a->n_charts;
        }else{ok=0;break;}
        if(!((oa==1||oa==-1)&&(ob==1||ob==-1))||
           (transaction_flags&~CHART_ZIPPER_TRANSACTION_PHYSICAL)){
            ok=0;break;
        }
        ok=sw_write(f,&la,sizeof(la),1)&&sw_write(f,&lb,sizeof(lb),1)&&
           sw_write(f,&oa,sizeof(oa),1)&&sw_write(f,&ob,sizeof(ob),1)&&
            sw_write(f,&transaction_flags,sizeof(transaction_flags),1)&&
           sw_write(f,&count,sizeof(count),1)&&
           sw_write(f,rank,sizeof(*rank),5)&&
           sw_write(f,&faces[tx[t].face_first*3],sizeof(*faces),
                     tx[t].face_count*3);
    }
    if(ok)ok=sw_write(f,&nconflicts64,sizeof(nconflicts64),1);
    for(size_t i=0;ok&&i<nconflicts;i++){
        ok=sw_write(f,&conflicts[i].transaction_a,
                    sizeof(conflicts[i].transaction_a),1)&&
           sw_write(f,&conflicts[i].transaction_b,
                    sizeof(conflicts[i].transaction_b),1)&&
           sw_write(f,&conflicts[i].reason_mask,
                    sizeof(conflicts[i].reason_mask),1);
    }
    if(ok)ok=fflush(f)==0;
    if(fclose(f)!=0)ok=0;
    return ok?0:-1;
}

int StreamWeldFormat_selftest(void)
{
    const char *path="stream_weld_format_selftest.gwshard";
    const char *patch_path="stream_weld_format_selftest.gwpatch";
    Arena_T a=Arena_new();StreamWeldShard r;
    float v[9]={0,0,0,1,0,0,0,1,0};
    int32_t f[3]={0,1,2},ch[3]={0,0,0};
    int32_t p0[3]={-2,-2,-2},p1[3]={0,1,2};
    int32_t patch_faces[6]={0,3,4,1,3,5};
    ChartZipperTransaction tx[2];
    StreamWeldConflict conflict={0,1,STREAM_WELD_CONFLICT_GEOMETRY};
    int fail=0;
    memset(&r,0,sizeof(r));
    memset(tx,0,sizeof(tx));
    if(StreamWeld_write_shard(path,"z00000_y00000_x00000",v,3,f,1,ch,1,
                              p0,p1)!=0||
       StreamWeld_read_shard(a,path,&r)!=0||r.nv!=3||r.nf!=1||
       r.n_charts!=1||strcmp(r.cube_id,"z00000_y00000_x00000")!=0||
       memcmp(r.verts,v,sizeof(v))!=0||memcmp(r.faces,f,sizeof(f))!=0||
       memcmp(r.parent1,p1,sizeof(p1))!=0)fail=1;
    if(!fail){
        StreamWeldShard b=r;FILE *pf=NULL;char magic[8];uint32_t version,flags;
        strcpy(b.cube_id,"z00000_y00000_x00128");
        for(size_t i=0;i<2;i++){
            tx[i].chart_a=0;tx[i].chart_b=1;
            tx[i].face_first=i;tx[i].face_count=1;
            tx[i].orient_a=1;tx[i].orient_b=1;
        }
        if(StreamWeld_write_patch(patch_path,&r,&b,patch_faces,2,tx,2,
                                  &conflict,1)!=0)
            fail=1;
        else{
            pf=fopen(patch_path,"rb");
            if(!pf||!sw_read(pf,magic,1,sizeof(magic))||
               !sw_read(pf,&version,sizeof(version),1)||
               !sw_read(pf,&flags,sizeof(flags),1)||
               memcmp(magic,"GWPATC3\n",8)!=0||version!=3||!(flags&1))
                fail=1;
            if(pf)fclose(pf);
        }
    }
    remove(patch_path);remove(path);Arena_dispose(&a);
    fprintf(stderr,"[selftest] stream weld shard/graph-patch round-trip -> %s\n",
            fail?"FAIL":"PASS");
    return fail;
}
