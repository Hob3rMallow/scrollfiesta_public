#include "textured_obj.h"

#include "../common/pipeline_constants.h"
#include "../common/tiff_io.h"
#include "../common/ves_png.h"
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int x0,y0,x1,y1;
    size_t count,begin;
} ToTile;

static void to_pixel(const RawtexPlan *p, double du, double dv,
                      double u, double v, double *x, double *y)
{ *x=(u-p->umin)/du; *y=(v-p->vmin)/dv; }

static void to_uv(const ToTile *t, double x, double y, double *s, double *v)
{
    /* Raster bounds are PIXEL EDGES. At a pixel center, s=(i+.5)/W.
     * One V flip converts top-left image rows to OBJ's bottom-left UVs. */
    *s=(x-t->x0)/(t->x1-t->x0);
    *v=1-(y-t->y0)/(t->y1-t->y0);
}

static int to_preamble(FILE *f, const char *name, const float *verts,
                        RawtexUv uv, const float *normals, size_t nv, int flat)
{
    fprintf(f,"# Presentation texture UVs; authoritative voxel UVs remain in VMESH.\n"
              "# Face vertex and texture-coordinate indices are independent.\n"
              "mtllib %s_rawtex.mtl\ns 1\n",name);
    for (size_t i=0; i<nv; i++) {
        if (flat) fprintf(f,"v %.9g %.9g 0\n",(double)Rawtex_uv(uv,2*i),(double)Rawtex_uv(uv,2*i+1));
        else fprintf(f,"v %.9g %.9g %.9g\n",(double)verts[3*i],
                     (double)verts[3*i+1],(double)verts[3*i+2]);
    }
    if (normals!=NULL) for (size_t i=0; i<nv; i++) {
        if (flat) fprintf(f,"vn 0 0 1\n");
        else fprintf(f,"vn %.9g %.9g %.9g\n",(double)normals[3*i],
                     (double)normals[3*i+1],(double)normals[3*i+2]);
    }
    return ferror(f) ? -1 : 0;
}

int TexturedObj_write_field(const char *prefix, const float *verts, RawtexUv uv,
                       const float *normals, size_t nv,
                       const int32_t *faces, size_t nf, const RawtexPlan *plan,
                       double du, double dv, int do_3d, int do_flat)
{
    Arena_T arena=NULL;
    TiffRowReader reader={0};
    ToTile *tiles=NULL;
    int32_t *face_tile=NULL, *order=NULL;
    size_t *vt_index=NULL, *cursor=NULL;
    FILE *obj[2]={NULL,NULL}, *mtl=NULL, *meta=NULL;
    char path[4096], filename[1024];
    const char *name=prefix;
    size_t nt=0, no_data=0, no_data_begin=0, next_vt=0, occupied=0;
    int nx=0, ny=0, rc=-1, max_width=0, max_height=0;
    if (prefix==NULL || plan==NULL || !plan->ok || !(du>0) || !(dv>0) ||
        !isfinite(du) || !isfinite(dv) || plan->W==0 || plan->H==0 ||
        plan->W>INT_MAX || plan->H>INT_MAX || nv>INT32_MAX || nf>INT32_MAX ||
        (nv && (verts==NULL || !Rawtex_has_uv(uv))) || (nf && faces==NULL) ||
        strlen(prefix)>sizeof path-128) return -1;
    if (!do_3d && !do_flat) return 0;
    for (const char *c=prefix; *c; c++) if (*c=='/' || *c=='\\') name=c+1;
    if (*name=='\0' || strlen(name)>sizeof filename-128) return -1;
    for (const char *c=name; *c; c++)
        if (!isalnum((unsigned char)*c) && *c!='_' && *c!='-' && *c!='.') return -1;
    nx=(int)((plan->W-1)/TEXOBJ_TILE_PIXELS+1);
    ny=(int)((plan->H-1)/TEXOBJ_TILE_PIXELS+1);
    nt=(size_t)nx*ny;
    if (nt>INT32_MAX) return -1;
    arena=Arena_new();
    tiles=(ToTile *)ARENA_CALLOC(arena,nt,sizeof *tiles);
    face_tile=(int32_t *)ARENA_ALLOC(arena,nf*sizeof(int32_t));
    order=(int32_t *)ARENA_ALLOC(arena,nf*sizeof(int32_t));
    vt_index=(size_t *)ARENA_CALLOC(arena,nv,sizeof(size_t));
    cursor=(size_t *)ARENA_CALLOC(arena,nt,sizeof(size_t));
    for (int y=0; y<ny; y++) for (int x=0; x<nx; x++) {
        ToTile *t=tiles+(size_t)y*nx+x;
        t->x0=x*TEXOBJ_TILE_PIXELS; t->y0=y*TEXOBJ_TILE_PIXELS;
        t->x1=t->x0+TEXOBJ_TILE_PIXELS; t->y1=t->y0+TEXOBJ_TILE_PIXELS;
        if (t->x1>(int)plan->W) t->x1=(int)plan->W;
        if (t->y1>(int)plan->H) t->y1=(int)plan->H;
    }
    for (size_t f=0; f<nf; f++) {
        double xmin=HUGE_VAL,ymin=HUGE_VAL,xmax=-HUGE_VAL,ymax=-HUGE_VAL,cx=0,cy=0;
        int tx=0,ty=0,tid=0;
        ToTile *t=NULL;
        for (int k=0; k<3; k++) {
            int id=faces[3*f+k];
            double x=0,y=0;
            if (id<0 || (size_t)id>=nv) goto done;
            to_pixel(plan,du,dv,Rawtex_uv(uv,2*(size_t)id),Rawtex_uv(uv,2*(size_t)id+1),&x,&y);
            if (!isfinite(x) || !isfinite(y)) goto done;
            xmin=fmin(xmin,x); xmax=fmax(xmax,x); ymin=fmin(ymin,y); ymax=fmax(ymax,y);
            cx+=x/3; cy+=y/3;
        }
        /* A partially uncovered face is visibly neutral, never clamped into
         * unrelated pixels. No geometry is removed by a display window. */
        if (xmin<0 || ymin<0 || xmax>plan->W || ymax>plan->H) {
            face_tile[f]=-1; no_data++; continue;
        }
        tx=(int)(cx/TEXOBJ_TILE_PIXELS); ty=(int)(cy/TEXOBJ_TILE_PIXELS);
        if (tx>=nx) tx=nx-1; if (ty>=ny) ty=ny-1;
        tid=ty*nx+tx; t=tiles+tid; face_tile[f]=tid; t->count++;
        /* Assign entire triangles to a material; extend its pixel rectangle
         * past every face corner plus a filtering gutter. Overlap pixels are
         * copied from ONE bake, so neither geometry splits nor new CT samples
         * are introduced at texture tile boundaries. */
        t->x0=(int)fmin(t->x0,floor(xmin)-TEXOBJ_GUTTER_PIXELS);
        t->y0=(int)fmin(t->y0,floor(ymin)-TEXOBJ_GUTTER_PIXELS);
        t->x1=(int)fmax(t->x1,ceil(xmax)+TEXOBJ_GUTTER_PIXELS);
        t->y1=(int)fmax(t->y1,ceil(ymax)+TEXOBJ_GUTTER_PIXELS);
        if (t->x1-t->x0>TEXOBJ_MAX_TILE_PIXELS || t->y1-t->y0>TEXOBJ_MAX_TILE_PIXELS) {
            fprintf(stderr,"ERROR: texture tile requires %dx%d pixels: subdivide long UV edges; refusing to downsample.\n",
                    t->x1-t->x0,t->y1-t->y0);
            goto done;
        }
    }
    for (size_t t=0; t<nt; t++) {
        tiles[t].begin=no_data_begin; cursor[t]=no_data_begin; no_data_begin+=tiles[t].count;
    }
    for (size_t f=0, at=no_data_begin; f<nf; f++) {
        if (face_tile[f]>=0) order[cursor[face_tile[f]]++]=(int32_t)f;
        else order[at++]=(int32_t)f;
    }
    snprintf(path,sizeof path,"%s_rawtex.tif",prefix);
    if (TiffIO_rows_open(arena,path,&reader)!=0 ||
        reader.W!=(int)plan->W || reader.H!=(int)plan->H) goto done;
    snprintf(path,sizeof path,"%s_rawtex.mtl",prefix); mtl=fopen(path,"wb");
    if (mtl==NULL) goto done;
    snprintf(path,sizeof path,"%s_rawtex_tiles.json",prefix); meta=fopen(path,"wb");
    if (meta==NULL) goto done;
    for (int flat=0; flat<2; flat++) if (flat ? do_flat : do_3d) {
        snprintf(path,sizeof path,"%s_%s.obj",prefix,flat ? "rawflat" : "raw3d");
        obj[flat]=fopen(path,"wb");
        if (obj[flat]==NULL || to_preamble(obj[flat],name,verts,uv,normals,nv,flat)!=0) goto done;
    }
    fprintf(meta,"{\"schema\":\"textured-obj-v1\",\"origin_uv\":[%.17g,%.17g],\"step_uv\":[%.17g,%.17g],"
                 "\"shape_hw\":[%zu,%zu],\"uv_rule\":\"s=(pixel_x-x0)/width; t=1-(pixel_y-y0)/height\","
                 "\"outer_gutter\":\"replicate-edge\",\"source_uv_modified\":false,\"vertex_colors\":false,\"faces\":%zu,\"no_data_faces\":%zu,\"tiles\":[\n",
            plan->umin,plan->vmin,du,dv,plan->H,plan->W,nf,no_data);
    for (size_t tid=0; tid<nt; tid++) {
        ToTile *t=tiles+tid;
        Arena_Mark mark;
        uint8_t *pixels=NULL;
        int w=t->x1-t->x0,h=t->y1-t->y0;
        if (t->count==0) continue;
        mark=Arena_save(arena);
        pixels=(uint8_t *)ARENA_ALLOC(arena,(size_t)w*h);
        for (int y=0; y<h; y++) {
            int sy=t->y0+y, x0=t->x0<0 ? 0 : t->x0;
            int x1=t->x1>reader.W ? reader.W : t->x1;
            int left=x0-t->x0, middle=x1-x0, right=t->x1-x1;
            uint8_t *row=pixels+(size_t)y*w;
            if (sy<0) sy=0; if (sy>=reader.H) sy=reader.H-1;
            if (middle<=0 || TiffIO_rows_read(&reader,sy,x0,middle,row+left)!=0) goto done;
            /* Give outermost canvas edges the same guard as internal tiles.
             * Then all face UVs lie INSIDE the tile and plain map_Kd works
             * even in readers that mistake '-clamp on' for a filename. */
            if (left) memset(row,row[left],(size_t)left);
            if (right) memset(row+left+middle,row[left+middle-1],(size_t)right);
        }
        snprintf(filename,sizeof filename,"%s_tile_%zu_%zu.png",name,tid/nx,tid%nx);
        snprintf(path,sizeof path,"%s_tile_%zu_%zu.png",prefix,tid/nx,tid%nx);
        if (VesPng_write_gray(path,pixels,w,h)!=0) goto done;
        Arena_restore(arena,mark);
        fprintf(mtl,"newmtl tile_%zu\nKa 1 1 1\nKd 1 1 1\nKs 0 0 0\nd 1\nillum 0\nmap_Kd %s\n\n",tid,filename);
        fprintf(meta,"%s{\"material\":\"tile_%zu\",\"file\":\"%s\",\"pixel_rect\":[%d,%d,%d,%d],\"faces\":%zu}",
                occupied ? ",\n" : "",tid,filename,t->x0,t->y0,t->x1,t->y1,t->count);
        occupied++; if (w>max_width) max_width=w; if (h>max_height) max_height=h;
        for (int flat=0; flat<2; flat++) if (obj[flat]) fprintf(obj[flat],"usemtl tile_%zu\n",tid);
        for (size_t j=t->begin; j<t->begin+t->count; j++) for (int k=0; k<3; k++) {
            int id=faces[3*(size_t)order[j]+k];
            double x=0,y=0,s=0,v=0;
            if (vt_index[id]) continue;
            vt_index[id]=++next_vt;
            to_pixel(plan,du,dv,Rawtex_uv(uv,2*(size_t)id),Rawtex_uv(uv,2*(size_t)id+1),&x,&y);
            to_uv(t,x,y,&s,&v);
            if (s<0 || s>1 || v<0 || v>1) goto done;
            for (int flat=0; flat<2; flat++) if (obj[flat]) fprintf(obj[flat],"vt %.17g %.17g\n",s,v);
        }
        for (size_t j=t->begin; j<t->begin+t->count; j++) {
            const int32_t *tri=faces+3*(size_t)order[j];
            for (int flat=0; flat<2; flat++) if (obj[flat]) {
                fprintf(obj[flat],"f");
                for (int k=0; k<3; k++) {
                    int id=tri[k];
                    if (normals) fprintf(obj[flat]," %d/%zu/%d",id+1,vt_index[id],id+1);
                    else fprintf(obj[flat]," %d/%zu",id+1,vt_index[id]);
                }
                fprintf(obj[flat],"\n");
            }
        }
        for (size_t j=t->begin; j<t->begin+t->count; j++) for (int k=0; k<3; k++)
            vt_index[faces[3*(size_t)order[j]+k]]=0;
    }
    if (no_data) {
        fprintf(mtl,"newmtl no_data\nKa 0.3 0.3 0.3\nKd 0.3 0.3 0.3\nKs 0 0 0\nillum 0\n");
        for (int flat=0; flat<2; flat++) if (obj[flat]) {
            fprintf(obj[flat],"usemtl no_data\nvt 0.5 0.5\n");
            for (size_t j=no_data_begin; j<nf; j++) {
                const int32_t *tri=faces+3*(size_t)order[j];
                fprintf(obj[flat],"f");
                for (int k=0; k<3; k++) {
                    if (normals) fprintf(obj[flat]," %d/%zu/%d",tri[k]+1,next_vt+1,tri[k]+1);
                    else fprintf(obj[flat]," %d/%zu",tri[k]+1,next_vt+1);
                }
                fprintf(obj[flat],"\n");
            }
        }
    }
    fprintf(meta,"\n]}\n");
    if (ferror(meta) || ferror(mtl)) goto done;
    for (int flat=0; flat<2; flat++) if (obj[flat] && ferror(obj[flat])) goto done;
    rc=0;
done:
    for (int flat=0; flat<2; flat++) if (obj[flat] && fclose(obj[flat])!=0) rc=-1;
    if (mtl && fclose(mtl)!=0) rc=-1;
    if (meta && fclose(meta)!=0) rc=-1;
    TiffIO_rows_close(&reader);
    if (arena) Arena_dispose(&arena);
    if (rc==0) fprintf(stderr,"[textured OBJ] %zu faces, %zu full-density tiles (max %dx%d), %zu no-data faces; no per-vertex colors\n",
                      nf,occupied,max_width,max_height,no_data);
    else fprintf(stderr,"ERROR: textured OBJ export failed; output is incomplete.\n");
    return rc;
}

#include "texture_pyramid.inc"

int TexturedObj_selftest(void)
{
    RawtexPlan p={0};
    ToTile t={4,2,10,7,0,0};
    double x=0,y=0,s=0,v=0;
    int fails=tp_selftest();
    p.umin=-7; p.vmin=20; p.W=13; p.H=11;
    /* Nonzero/negative origin, unequal pixel steps and a padded canvas. */
    to_pixel(&p,2,.5,2,21.25,&x,&y);
    to_uv(&t,x,y,&s,&v);
    fails+=fabs(s-.5/6)>1e-15 || fabs(v-(1-.5/5))>1e-15;
    to_uv(&t,10,7,&s,&v); fails+=s!=1 || v!=0;
    to_uv(&t,4,2,&s,&v); fails+=s!=0 || v!=1;
    fprintf(stderr,"[selftest] textured OBJ coordinate transform %s (%d failures)\n",fails ? "FAIL" : "PASS",fails);
    return fails ? -1 : 0;
}

int TexturedObj_write(const char *prefix, const float *verts, const float *uv,
                       const float *normals, size_t nv,
                       const int32_t *faces, size_t nf, const RawtexPlan *plan,
                       double du, double dv, int do_3d, int do_flat)
{
    return TexturedObj_write_field(prefix,verts,(RawtexUv){uv,NULL},normals,nv,faces,nf,plan,du,dv,do_3d,do_flat);
}
