#include "closeup_import.h"
#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/mesh_normals.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int ci_file(const char *path, const void *data, size_t bytes)
{
    FILE *fp=fopen(path,"wb");
    if(!fp)return -1;
    int ok=fwrite(data,1,bytes,fp)==bytes;
    if(fclose(fp))ok=0;
    return ok?0:-1;
}

static void ci_json(FILE *fp,const char *s)
{
    fputc('"',fp);
    for(;*s;s++){
        unsigned char c=(unsigned char)*s;
        if(c=='"' || c=='\\')fputc('\\',fp);
        if(c<32)fprintf(fp,"\\u%04x",(unsigned)c);else fputc(c,fp);
    }
    fputc('"',fp);
}

typedef struct CiChart { int32_t local,source; } CiChart;

static int ci_order(const void *a,const void *b)
{
    const CiChart *x=a,*y=b;
    return (x->local>y->local)-(x->local<y->local);
}

/* Exact subset exporters retain local/source correspondence in charts.csv. */
static int ci_remap(Arena_T arena,const char *path,int32_t *chart,size_t nv)
{
    FILE *fp=fopen(path,"rb");char line[4096];size_t n=0,k=0;
    if(!fp)return -1;
    if(!fgets(line,sizeof line,fp) || strncmp(line,"local,source,",13)){fclose(fp);return -1;}
    while(fgets(line,sizeof line,fp))n++;
    if(!n || n>SIZE_MAX/sizeof(CiChart) || fseek(fp,0,SEEK_SET)){fclose(fp);return -1;}
    CiChart *map=ARENA_ALLOC(arena,n*sizeof *map);
    if(!fgets(line,sizeof line,fp)){fclose(fp);return -1;}
    while(k<n && fgets(line,sizeof line,fp)){
        long long a=0,b=0;
        if(sscanf(line,"%lld,%lld,",&a,&b)!=2 || a<0 || a>INT32_MAX || b<0 || b>INT32_MAX){fclose(fp);return -1;}
        map[k++]=(CiChart){(int32_t)a,(int32_t)b};
    }
    fclose(fp);if(k!=n)return -1;
    qsort(map,n,sizeof *map,ci_order);
    for(k=1;k<n;k++)if(map[k-1].local==map[k].local)return -1;
    for(size_t v=0;v<nv;v++){
        CiChart key={chart[v],0};CiChart *found=bsearch(&key,map,n,sizeof *map,ci_order);
        if(!found)return -1;
        chart[v]=found->source;
    }
    return 0;
}

typedef struct CiNormalTip {
    size_t vertex, faces[2], count;
    float normal[3];
    double min_neighbor_cos;
} CiNormalTip;

/* Some original meshes contain a doubled-back tip: two nondegenerate faces
 * with identical world coordinates and opposite reading-oriented normals.
 * Their shared tip has a zero area-weighted normal. Retain both original
 * faces, but extend the unanimous surrounding sampling direction to this
 * isolated vertex. This does not repair or certify the source geometry.
 * Other undefined normals (unused vertices, broad folds, ambiguous neighbors)
 * remain an import failure. Never propagate through another undefined vertex. */
static int ci_tip_normals(Arena_T arena,const MeshBinData *m,const int32_t *oriented,
                         const int32_t *chart,float *normals,CiNormalTip **tips,size_t *ntips)
{
    *tips=NULL;*ntips=0;
    for(size_t v=0;v<m->nv;v++){
        double n=0;
        for(size_t d=0;d<3;d++)n+=(double)normals[3*v+d]*normals[3*v+d];
        if(!isfinite(n))return -1;
        if(n<0.5)(*ntips)++;
    }
    if(!*ntips)return 0;
    if(*ntips>SIZE_MAX/sizeof(CiNormalTip) || m->nv>SIZE_MAX/sizeof(size_t))return -1;
    CiNormalTip *tip=ARENA_ALLOC(arena,*ntips*sizeof *tip);
    size_t *index=ARENA_ALLOC(arena,m->nv*sizeof *index),next=0;
    memset(tip,0,*ntips*sizeof *tip);
    for(size_t v=0;v<m->nv;v++){
        double n=0;
        for(size_t d=0;d<3;d++)n+=(double)normals[3*v+d]*normals[3*v+d];
        index[v]=SIZE_MAX;
        if(n<0.5){index[v]=next;tip[next++].vertex=v;}
    }
    for(size_t f=0;f<m->nf;f++)for(size_t k=0;k<3;k++){
        size_t v=(size_t)oriented[3*f+k],i=index[v];
        if(i==SIZE_MAX)continue;
        if(tip[i].count<2)tip[i].faces[tip[i].count]=f;
        tip[i].count++;
    }
    for(size_t i=0;i<*ntips;i++){
        CiNormalTip *t=&tip[i];
        if(t->count!=2){
            fprintf(stderr,"[closeup import] undefined normal at vertex %zu: %zu incident faces\n",t->vertex,t->count);
            return -1;
        }
        const int32_t *a=oriented+3*t->faces[0],*b=oriented+3*t->faces[1];
        int matched[3]={0,0,0};
        for(size_t k=0;k<3;k++){
            int found=0;
            for(size_t j=0;j<3;j++)if(!matched[j] &&
                !memcmp(m->verts+3*(size_t)a[k],m->verts+3*(size_t)b[j],3*sizeof(float))){
                matched[j]=1;found=1;break;
            }
            if(!found){
                fprintf(stderr,"[closeup import] undefined normal at vertex %zu: not an exact source tip\n",t->vertex);
                return -1;
            }
        }
        /* The already computed zero sum plus two identical nondegenerate
         * world triangles establishes opposite directions. Reject a collapsed
         * 3-D face even if the UVs happen to have a nonzero area. */
        double e1[3]={0},e2[3]={0},cross[3]={0};
        for(size_t d=0;d<3;d++){
            e1[d]=(double)m->verts[3*(size_t)a[1]+d]-m->verts[3*(size_t)a[0]+d];
            e2[d]=(double)m->verts[3*(size_t)a[2]+d]-m->verts[3*(size_t)a[0]+d];
        }
        cross[0]=e1[1]*e2[2]-e1[2]*e2[1];cross[1]=e1[2]*e2[0]-e1[0]*e2[2];cross[2]=e1[0]*e2[1]-e1[1]*e2[0];
        double area2=cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2];
        if(!isfinite(area2) || area2<=1e-24)return -1;
        size_t neighbors[4]={0},nn=0;
        for(size_t face=0;face<2;face++)for(size_t k=0;k<3;k++){
            size_t v=(size_t)oriented[3*t->faces[face]+k];
            if(v==t->vertex)continue;
            if(chart[v]!=chart[t->vertex] || index[v]!=SIZE_MAX)return -1;
            size_t j=0;for(;j<nn;j++)if(neighbors[j]==v)break;
            if(j==nn)neighbors[nn++]=v;
        }
        if(nn!=3)return -1;
        double sum[3]={0},min_cos=1;
        for(size_t k=0;k<nn;k++){
            const float *n=normals+3*neighbors[k];
            for(size_t d=0;d<3;d++)sum[d]+=n[d];
            for(size_t j=0;j<k;j++){
                const float *p=normals+3*neighbors[j];
                double dot=(double)n[0]*p[0]+(double)n[1]*p[1]+(double)n[2]*p[2];
                min_cos=fmin(min_cos,dot);
            }
        }
        if(min_cos<CLOSEUP_IMPORT_NEIGHBOR_NORMAL_COS){
            fprintf(stderr,"[closeup import] source tip %zu has ambiguous neighboring normals (cos %.9g)\n",t->vertex,min_cos);
            return -1;
        }
        double length=sqrt(sum[0]*sum[0]+sum[1]*sum[1]+sum[2]*sum[2]);
        if(!isfinite(length) || length<1e-12)return -1;
        for(size_t d=0;d<3;d++)t->normal[d]=(float)(sum[d]/length);
        t->min_neighbor_cos=min_cos;
    }
    for(size_t i=0;i<*ntips;i++)memcpy(normals+3*tip[i].vertex,tip[i].normal,3*sizeof(float));
    *tips=tip;
    fprintf(stderr,"[closeup import] retained and flagged %zu exact source tips; sampling normals from agreeing neighbors\n",*ntips);
    return 0;
}

int CloseupImport_write(const char *mesh,const char *charts,const char *stem,const char *chart_map)
{
    const char *suffix[]={"_surface.vmesh","_surface_chart.i32","_surface_attr.u8",
                          "_surface_nrm.f32","_surface_report.json"};
    char paths[5][4096];
    Arena_T arena=Arena_new();
    MeshBinData m={0};
    float *normals=NULL;
    CiNormalTip *tips=NULL;size_t ntips=0;
    FILE *fp=NULL;
    int rc=-1;
    const char *stage="output preflight";
    for(size_t k=0;k<5;k++){
        if(snprintf(paths[k],sizeof paths[k],"%s%s",stem,suffix[k])>=(int)sizeof paths[k])goto done;
        fp=fopen(paths[k],"rb");
        if(fp){fclose(fp);fp=NULL;fprintf(stderr,"[closeup import] existing output refused: %s\n",paths[k]);goto done;}
    }
    stage="mesh read";
    if(MeshBin_read_precise_arena(arena,mesh,&m) || !m.nv || !m.nf || (!m.uv && !m.uv64) ||
       m.nv>SIZE_MAX/(3*sizeof(double)) || m.nf>SIZE_MAX/(3*sizeof(int32_t)))goto done;
    double *uv=ARENA_ALLOC(arena,2*m.nv*sizeof *uv);
    int32_t *chart=ARENA_ALLOC(arena,m.nv*sizeof *chart);
    int32_t *oriented=ARENA_ALLOC(arena,3*m.nf*sizeof *oriented);
    uint8_t *attr=ARENA_ALLOC(arena,m.nv);
    stage="chart labels";fp=fopen(charts,"rb");
    if(!fp)goto done;
    if(fread(chart,sizeof *chart,m.nv,fp)!=m.nv || fgetc(fp)!=EOF)goto done;
    fclose(fp);fp=NULL;
    stage="chart remap";
    if(chart_map && ci_remap(arena,chart_map,chart,m.nv))goto done;
    stage="finite coordinates";
    double low[2]={INFINITY,INFINITY},shift[2]={0,0};
    for(size_t v=0;v<m.nv;v++){
        for(size_t d=0;d<3;d++)if(!isfinite(m.verts[3*v+d]))goto done;
        for(size_t d=0;d<2;d++){
            double x=MeshBin_uv(&m,2*v+d);
            if(!isfinite(x))goto done;
            low[d]=fmin(low[d],x);
        }
    }
    for(size_t d=0;d<2;d++)shift[d]=(double)CLOSEUP_RASTER_MARGIN_PX-floor(low[d]);
    for(size_t v=0;v<m.nv;v++)for(size_t d=0;d<2;d++)uv[2*v+d]=MeshBin_uv(&m,2*v+d)+shift[d];
    /* Normal direction follows the reviewed reading canvas, independently of
     * each chart's stored triangle winding. Original face order is retained. */
    stage="UV face orientation";
    for(size_t f=0;f<m.nf;f++){
        const int32_t *t=m.faces+3*f;
        size_t a=(size_t)t[0],b=(size_t)t[1],c=(size_t)t[2];
        if(a>=m.nv || b>=m.nv || c>=m.nv)goto done;
        double det=(uv[2*b]-uv[2*a])*(uv[2*c+1]-uv[2*a+1])-
                   (uv[2*b+1]-uv[2*a+1])*(uv[2*c]-uv[2*a]);
        if(!isfinite(det) || fabs(det)<1e-20)goto done;
        oriented[3*f]=t[0];oriented[3*f+1]=t[det>0?1:2];oriented[3*f+2]=t[det>0?2:1];
    }
    stage="sampling normals";
    normals=MeshNormals_compute(m.verts,m.nv,oriented,m.nf);
    if(!normals)goto done;
    if(ci_tip_normals(arena,&m,oriented,chart,normals,&tips,&ntips))goto done;
    memset(attr,1,m.nv);
    stage="output write";
    ves_ensure_parent_dir(paths[0]);
    if(ci_file(paths[1],chart,m.nv*sizeof *chart) || ci_file(paths[2],attr,m.nv) ||
       ci_file(paths[3],normals,3*m.nv*sizeof *normals) ||
       MeshBin_write_uv64(paths[0],m.verts,m.nv,m.faces,m.nf,uv))goto done;
    fp=fopen(paths[4],"wb");if(!fp)goto done;
    fputs("{\"schema\":\"ink-closeup-import-v1\",\"source_mesh\":",fp);ci_json(fp,mesh);
    fputs(",\"source_charts\":",fp);ci_json(fp,charts);
    fputs(",\"chart_map\":",fp);if(chart_map)ci_json(fp,chart_map);else fputs("null",fp);
    fprintf(fp,",\"vertices\":%zu,\"faces\":%zu,\"uv_translation\":[%.17g,%.17g],"
               "\"world_geometry_preserved\":true,\"face_order_preserved\":true,"
               "\"normal_orientation\":\"reviewed positive UV reading canvas\",\"fill_vertices\":0,"
               "\"retained_source_tip_count\":%zu,\"retained_source_tips\":[",
               m.nv,m.nf,shift[0],shift[1],ntips);
    for(size_t i=0;i<ntips;i++){
        const CiNormalTip *t=&tips[i];
        fprintf(fp,"%s{\"vertex\":%zu,\"source_chart\":%d,\"faces\":[%zu,%zu],"
                   "\"neighbor_min_cos\":%.17g,\"sampling_normal\":[%.9g,%.9g,%.9g],"
                   "\"geometry_repaired\":false,\"flag\":\"exact doubled-back source tip retained\"}",
                   i?",":"",t->vertex,chart[t->vertex],t->faces[0],t->faces[1],t->min_neighbor_cos,
                   (double)t->normal[0],(double)t->normal[1],(double)t->normal[2]);
    }
    fputs("]}\n",fp);
    int ok=!ferror(fp);if(fclose(fp))ok=0;fp=NULL;
    if(!ok)goto done;
    fprintf(stderr,"[closeup import] retained %zu vertices, %zu faces; UV translation [%.17g, %.17g]\n",
            m.nv,m.nf,shift[0],shift[1]);rc=0;
done:
    if(fp)fclose(fp);
    free(normals);Arena_dispose(&arena);
    if(rc)fprintf(stderr,"[closeup import] failed at %s; source and any partial evidence retained\n",stage);
    return rc;
}

static int ci_tip_selftest(const char *dir)
{
    /* Two coincident world triangles meet at vertex 0 but unfold to distinct
     * UV triangles. Three neighboring vertices have independent +Z support. */
    float xyz[24]={0,0,0, 1,0,0, 0,1,0, 1,0,0, 1,1,0, 2,0,0, 1,2,0, 5,5,0};
    double uv[16]={0,0, 1,0, 0,1, -1,0, 1,1, 0,0, -1,2, 5,5};
    int32_t faces[12]={0,1,2, 0,2,3, 1,4,2, 3,5,6},ids[8]={17,17,17,17,17,17,17,17};
    char mesh[2400],charts[2400],stem[2400],path[2500],report[2048]={0};
    MeshBinData got={0};int failures=0;float normal[21]={0};FILE *fp=NULL;
    snprintf(mesh,sizeof mesh,"%s/tip.vmesh",dir);snprintf(charts,sizeof charts,"%s/tip.i32",dir);
    snprintf(stem,sizeof stem,"%s/tip_imported",dir);
    if(MeshBin_write_uv64(mesh,xyz,7,faces,4,uv) || ci_file(charts,ids,7*sizeof *ids) ||
       CloseupImport_write(mesh,charts,stem,NULL))return 1;
    snprintf(path,sizeof path,"%s_surface.vmesh",stem);
    if(MeshBin_read_precise_malloc(path,&got) || got.nv!=7 || got.nf!=4 ||
       memcmp(got.verts,xyz,21*sizeof(float)) || memcmp(got.faces,faces,sizeof faces))failures++;
    MeshBin_dispose(&got);
    snprintf(path,sizeof path,"%s_surface_nrm.f32",stem);fp=fopen(path,"rb");
    if(!fp || fread(normal,sizeof(float),21,fp)!=21)failures++;
    if(fp)fclose(fp);
    for(size_t v=0;v<7;v++)if(fabs(normal[3*v])+fabs(normal[3*v+1])+fabs(normal[3*v+2]-1)>1e-6)failures++;
    snprintf(path,sizeof path,"%s_surface_report.json",stem);fp=fopen(path,"rb");
    if(!fp)failures++;
    else{fread(report,1,sizeof report-1,fp);fclose(fp);}
    if(!strstr(report,"\"retained_source_tip_count\":1") || !strstr(report,"\"geometry_repaired\":false"))failures++;
    /* A reversed supporting patch makes the neighbors disagree: refuse. */
    uv[13]=-2;
    snprintf(stem,sizeof stem,"%s/tip_ambiguous",dir);
    if(MeshBin_write_uv64(mesh,xyz,7,faces,4,uv) || !CloseupImport_write(mesh,charts,stem,NULL))failures++;
    snprintf(path,sizeof path,"%s_surface.vmesh",stem);if(MeshBin_looks_complete(path))failures++;
    /* Unsupported vertices are not silently given an invented direction. */
    uv[13]=2;snprintf(stem,sizeof stem,"%s/tip_unused",dir);
    if(MeshBin_write_uv64(mesh,xyz,8,faces,4,uv) || ci_file(charts,ids,sizeof ids) ||
       !CloseupImport_write(mesh,charts,stem,NULL))failures++;
    snprintf(path,sizeof path,"%s_surface.vmesh",stem);if(MeshBin_looks_complete(path))failures++;
    fprintf(stderr,"[closeup import] exact source tip retained, flagged; ambiguous/unused normals refused: %s\n",failures?"FAIL":"passed");
    return failures;
}

int CloseupImport_selftest(void)
{
    char dir[2048],mesh[2300],charts[2300],stem[2300],path[2400];
    float xyz[12]={0,0,0, 1,0,0, 1,1,0, 0,1,0};
    double uv[8]={-100,300,-99,300,-99,301,-100,301};
    int32_t faces[6]={0,1,2, 0,3,2},ids[4]={17,17,17,17};
    MeshBinData got={0};int failures=0;
    if(ves_selftest_dir(dir,sizeof dir,"closeup_import"))return 1;
    snprintf(mesh,sizeof mesh,"%s/input.vmesh",dir);snprintf(charts,sizeof charts,"%s/charts.i32",dir);
    snprintf(stem,sizeof stem,"%s/imported",dir);
    if(MeshBin_write_uv64(mesh,xyz,4,faces,2,uv) || ci_file(charts,ids,sizeof ids) ||
       CloseupImport_write(mesh,charts,stem,NULL)){failures++;goto done;}
    snprintf(path,sizeof path,"%s_surface.vmesh",stem);
    if(MeshBin_read_precise_malloc(path,&got) || got.nv!=4 || got.nf!=2 ||
       memcmp(got.verts,xyz,sizeof xyz) || memcmp(got.faces,faces,sizeof faces)){failures++;goto done;}
    for(size_t v=0;v<4;v++)for(size_t d=0;d<2;d++)
        if(fabs((got.uv64[2*v+d]-got.uv64[d])-(uv[2*v+d]-uv[d]))>1e-12)failures++;
    snprintf(path,sizeof path,"%s_surface_nrm.f32",stem);
    FILE *fp=fopen(path,"rb");float n[12]={0};
    if(!fp || fread(n,sizeof(float),12,fp)!=12)failures++;
    if(fp)fclose(fp);
    for(size_t v=0;v<4;v++)if(fabs(n[3*v+2]-1)>1e-6)failures++;
    if(!CloseupImport_write(mesh,charts,stem,NULL))failures++;
    snprintf(path,sizeof path,"%s/map.csv",dir);
    const char mapping[]="local,source,vertices\n17,9107,4\n";
    snprintf(stem,sizeof stem,"%s/mapped",dir);
    if(ci_file(path,mapping,sizeof mapping-1) || CloseupImport_write(mesh,charts,stem,path))failures++;
    snprintf(path,sizeof path,"%s_surface_chart.i32",stem);
    fp=fopen(path,"rb");int32_t mapped[4]={0};
    if(!fp || fread(mapped,sizeof(int32_t),4,fp)!=4)failures++;
    if(fp)fclose(fp);
    for(size_t v=0;v<4;v++)if(mapped[v]!=9107)failures++;
    snprintf(stem,sizeof stem,"%s/truncated",dir);
    if(ci_file(charts,ids,sizeof ids-sizeof ids[0]) || !CloseupImport_write(mesh,charts,stem,NULL))failures++;
    snprintf(path,sizeof path,"%s_surface.vmesh",stem);
    if(MeshBin_looks_complete(path))failures++;
    uv[2]=uv[0];uv[3]=uv[1];
    snprintf(stem,sizeof stem,"%s/degenerate",dir);
    if(MeshBin_write_uv64(mesh,xyz,4,faces,2,uv) || ci_file(charts,ids,sizeof ids) ||
       !CloseupImport_write(mesh,charts,stem,NULL))failures++;
done:
    MeshBin_dispose(&got);
    failures+=ci_tip_selftest(dir);
    fprintf(stderr,"[closeup import] preservation, mixed winding, overwrite and invalid-source controls: %s\n",
            failures?"FAIL":"passed");
    if(!failures)ves_selftest_clean(dir);
    return failures;
}
