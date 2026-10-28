#include "fiber_map.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "sparse_solve.h"
#include "../common/pipeline_constants.h"
#include "../common/uv_guard.h"
#include "../common/ves_platform.h"

typedef struct FmSample {
    size_t node;
    double x, y;
} FmSample;

typedef struct FmMetric {
    double length, along, height;
    double lower, upper;
} FmMetric;

int FiberMap_read_directions(Arena_T arena, const char *path, FiberMapDirections *out)
{
    if (!arena || !path || !out) return -1;
    *out=(FiberMapDirections){0};
    FILE *stream=fopen(path,"r");
    if (!stream) return -1;
    const char *header="fiber-directions-v1,column,row,h_col,h_row,v_col,v_row,weight";
    char line[512];size_t count=0;
    if (!fgets(line,sizeof line,stream) || strncmp(line,header,strlen(header)) ||
        (line[strlen(header)]!='\n' && line[strlen(header)]!='\r')) { fclose(stream); return -1; }
    for (const char *p=line+strlen(header);*p;p++) if (*p!='\r' && *p!='\n') { fclose(stream); return -1; }
    while (fgets(line,sizeof line,stream)) count++;
    if (ferror(stream) || !count || count>SIZE_MAX/sizeof(FiberMapDirection)) { fclose(stream); return -1; }
    rewind(stream);
    if (!fgets(line,sizeof line,stream)) { fclose(stream); return -1; }
    FiberMapDirection *data=ARENA_ALLOC(arena,count*sizeof *data);
    for (size_t i=0;i<count;i++) {
        FiberMapDirection *s=data+i;int used=0;
        if (!fgets(line,sizeof line,stream) || sscanf(line,"%lf,%lf,%lf,%lf,%lf,%lf,%lf%n",
            &s->column,&s->row,s->h,s->h+1,s->v,s->v+1,&s->weight,&used)!=7) { fclose(stream); return -1; }
        for (const char *p=line+used;*p;p++) if (*p!='\r' && *p!='\n' && *p!=' ' && *p!='\t') { fclose(stream); return -1; }
        double determinant=s->h[0]*s->v[1]-s->h[1]*s->v[0];
        if (!isfinite(s->column) || !isfinite(s->row) || !isfinite(s->h[0]) || !isfinite(s->h[1]) ||
            !isfinite(s->v[0]) || !isfinite(s->v[1]) || !(s->weight>0) || !isfinite(s->weight) ||
            !(determinant>0) || !isfinite(determinant)) { fclose(stream); return -1; }
    }
    if (fclose(stream)) return -1;
    *out=(FiberMapDirections){count,data};return 0;
}

/* Intrinsic coordinates of the ORIGINAL 3-D face. Cross-product height avoids
 * losing thin source triangles to subtraction of nearly equal Gram terms. */
static int fm_reference(const float *xyz, const int32_t face[3], FmMetric *ref)
{
    double e[3], d[3], cross[3], length2=0, dot=0;
    for(size_t k=0;k<3;k++){
        e[k]=(double)xyz[3*(size_t)face[1]+k]-xyz[3*(size_t)face[0]+k];
        d[k]=(double)xyz[3*(size_t)face[2]+k]-xyz[3*(size_t)face[0]+k];
        length2+=e[k]*e[k];dot+=e[k]*d[k];
    }
    for(size_t k=0;k<3;k++)cross[k]=e[(k+1)%3]*d[(k+2)%3]-e[(k+2)%3]*d[(k+1)%3];
    double length=sqrt(length2), area2=hypot(hypot(cross[0],cross[1]),cross[2]);
    if(!(length>0) || !(area2>0) || !isfinite(length) || !isfinite(area2) || !isfinite(dot))return -1;
    ref->length=length;ref->along=dot/length;ref->height=area2/length;return 0;
}

static void fm_metric(const FmMetric *ref, const int32_t face[3], const double *uv,
                      double sx, double sy, double *lo, double *hi)
{
    size_t p=2*(size_t)face[0], q=2*(size_t)face[1], r=2*(size_t)face[2];
    double a=(uv[q]-uv[p])*sx/ref->length, c=(uv[q+1]-uv[p+1])*sy/ref->length;
    double b=((uv[r]-uv[p])*sx-ref->along*a)/ref->height;
    double d=((uv[r+1]-uv[p+1])*sy-ref->along*c)/ref->height;
    *hi=.5*(hypot(a+d,c-b)+hypot(a-d,c+b));
    *lo=*hi>0?fabs(a*d-b*c)/ *hi:0;
}

/* Same per-face contract as AsmMetric_face_preserved: a passing face keeps the
 * physical band; a pre-existing failure may not worsen either singular value.
 * Improvement of stretch cannot pay for new compression elsewhere on the face. */
static int fm_metric_preserved(const FmMetric *ref, double lo, double hi)
{
    return lo>0 && isfinite(lo) && isfinite(hi) && lo<=hi+ASM_METRIC_ROUNDOFF &&
           lo>=ref->lower-ASM_METRIC_ROUNDOFF &&
           hi<=ref->upper+ASM_METRIC_ROUNDOFF*fmax(1,ref->upper);
}

static void fm_affine(const double a[6], double x, double y, double out[2])
{
    out[0]=a[0]*x+a[1]*y+a[2];
    out[1]=a[3]*x+a[4]*y+a[5];
}

static size_t fm_neighbors(size_t q, size_t rows, size_t cols, size_t out[4])
{
    size_t n=0, x=q%cols, y=q/cols;
    if(x)out[n++]=q-1;
    if(x+1<cols)out[n++]=q+1;
    if(y)out[n++]=q-cols;
    if(y+1<rows)out[n++]=q+cols;
    return n;
}

/* Dirichlet constraints keep every unreleased measured node exactly. Solve
 * the displacement from the affine frame, so an affine field extends exactly
 * at the exterior Neumann boundary too. The rectangular graph is connected;
 * at least one retained constraint makes this Laplacian positive definite. */
static int fm_extend(Arena_T scratch, FiberMapResult *map, const double *base, const uint8_t *fixed)
{
    Arena_Mark mark=Arena_save(scratch);
    size_t nn=(size_t)map->rows*(size_t)map->cols, n=0, nt=0;
    int rc=-1;
    SparseFactor_T factor=NULL;
    int *id=ARENA_ALLOC(scratch,nn*sizeof *id);
    for(size_t q=0;q<nn;q++)id[q]=fixed[q]?-1:(int)n++;
    if(!n){rc=0;goto done;}
    if(n==nn || n>(size_t)INT_MAX/3)goto done;
    int *row=ARENA_ALLOC(scratch,3*n*sizeof *row), *col=ARENA_ALLOC(scratch,3*n*sizeof *col);
    double *val=ARENA_ALLOC(scratch,3*n*sizeof *val);
    double *rhs=ARENA_CALLOC(scratch,2*n,sizeof *rhs), *solution=ARENA_ALLOC(scratch,2*n*sizeof *solution);
    for(size_t q=0;q<nn;q++)if(id[q]>=0){
        size_t nb[4], count=fm_neighbors(q,(size_t)map->rows,(size_t)map->cols,nb);
        int a=id[q];
        row[nt]=col[nt]=a;val[nt++]=(double)count;
        for(size_t k=0;k<count;k++){
            size_t r=nb[k];int b=id[r];
            if(b<0){
                rhs[(size_t)a]+=map->uv[2*r]-base[2*r];
                rhs[n+(size_t)a]+=map->uv[2*r+1]-base[2*r+1];
            }else if(b<a){row[nt]=a;col[nt]=b;val[nt++]=-1;}
        }
    }
    if(nt>3*n || Sparse_factor_spd((int)n,(int)nt,row,col,val,&factor) ||
       Sparse_factor_solve_multi(factor,rhs,solution,2))goto done;
    for(size_t q=0;q<nn;q++)if(id[q]>=0){
        map->uv[2*q]=base[2*q]+solution[(size_t)id[q]];
        map->uv[2*q+1]=base[2*q+1]+solution[n+(size_t)id[q]];
        if(!isfinite(map->uv[2*q]) || !isfinite(map->uv[2*q+1]))goto done;
    }
    rc=0;
done:
    Sparse_factor_free(&factor);Arena_restore(scratch,mark);return rc;
}

static void fm_value(const FiberMapResult *map, const FmSample *s, double out[2])
{
    size_t a=s->node, b=a+1, c=a+(size_t)map->cols, d=c+1;
    for(size_t k=0;k<2;k++){
        double top=(1-s->x)*map->uv[2*a+k]+s->x*map->uv[2*b+k];
        double bottom=(1-s->x)*map->uv[2*c+k]+s->x*map->uv[2*d+k];
        out[k]=(1-s->y)*top+s->y*bottom;
    }
}

static void fm_seed_cell(uint8_t *seed, size_t a, size_t cols)
{seed[a]=seed[a+1]=seed[a+cols]=seed[a+cols+1]=1;}

/* A bilinear cell's determinant is affine in x,y. Its four corner
 * determinants are these oriented triples. Certifying their entire paths
 * also protects the whole cell throughout the affine-to-map homotopy. */
static size_t fm_cells(const FiberMapResult *map, const double *base, uint8_t *seed)
{
    size_t bad=0;
    for(int y=0;y+1<map->rows;y++)for(int x=0;x+1<map->cols;x++){
        int32_t a=y*map->cols+x, b=a+1, c=a+map->cols, d=c+1;
        int32_t corners[12]={a,b,c, a,b,d, a,d,c, b,d,c};
        int okay=1;
        for(size_t t=0;t<4;t++)if(!UvGuard_interval(base,map->uv,corners+3*t,1))okay=0;
        if(!okay){bad++;fm_seed_cell(seed,(size_t)a,(size_t)map->cols);}
    }
    return bad;
}

/* Min-plus distance transform on the rectangular four-neighbour grid. Two
 * raster sweeps give the exact Manhattan envelope for arbitrary finite seed
 * costs, including the negative reach costs used below. */
static void fm_distance(int64_t *cost, size_t rows, size_t cols)
{
    size_t nn=rows*cols;
    for(size_t q=0;q<nn;q++){
        if(q>=cols && cost[q-cols]+1<cost[q])cost[q]=cost[q-cols]+1;
        if(q%cols && cost[q-1]+1<cost[q])cost[q]=cost[q-1]+1;
    }
    for(size_t end=nn;end>0;end--){
        size_t q=end-1;
        if(q+cols<nn && cost[q+cols]+1<cost[q])cost[q]=cost[q+cols]+1;
        if(q%cols+1<cols && cost[q+1]+1<cost[q])cost[q]=cost[q+1]+1;
    }
}

/* Every failure must reach ITS first fixed boundary, plus a small halo. A
 * shared BFS stopping depth lets a shallow failure starve a separate failure
 * deep inside an already free region. Compute each seed's nearest-boundary
 * distance, then the union of their individual reach balls, in linear time.
 * Constraints outside that union stay exact. The whole field is recertified. */
static size_t fm_release(Arena_T scratch, const FiberMapResult *map, uint8_t *fixed, uint8_t *seed)
{
    Arena_Mark mark=Arena_save(scratch);
    size_t rows=(size_t)map->rows, cols=(size_t)map->cols, nn=rows*cols, count=0;
    int have_seed=0,have_fixed=0;
    for(size_t q=0;q<nn;q++){have_seed|=seed[q]!=0;have_fixed|=fixed[q]==1;}
    if(!have_seed || !have_fixed)return 0;
    const int64_t inf=INT64_MAX/4;
    int64_t *cost=ARENA_ALLOC(scratch,nn*sizeof *cost);
    for(size_t q=0;q<nn;q++)cost[q]=fixed[q]==1?0:inf;
    fm_distance(cost,rows,cols);
    for(size_t q=0;q<nn;q++)cost[q]=seed[q]?-cost[q]-FIBER_MAP_REPAIR_HALO_CELLS:inf;
    fm_distance(cost,rows,cols);
    for(size_t q=0;q<nn;q++)if(fixed[q]==1 && cost[q]<=0){fixed[q]=0;count++;}
    Arena_restore(scratch,mark);return count;
}

#include "fiber_map_metric.inc"
#include "fiber_map_pinned.inc"

static int fm_fit_grid(Arena_T arena, const FiberMapGrid *grid, const float *xyz, const double *canvas, size_t nv,
                 const int32_t *faces, size_t nf, double *output, FiberMapResult *result, FiberMapStats *st)
{
    double started=ves_clock_sec();int rc=-1;
    if(!arena || !grid || !result || !st || !grid->u || !grid->v ||
       (nv && (!xyz || !canvas || !output)) || (nf && !faces) || grid->rows<2 || grid->cols<2 ||
       !(grid->cell>0) || !isfinite(grid->cell) || !isfinite(grid->row0) || !isfinite(grid->col0))return -1;
    memset(result,0,sizeof *result);memset(st,0,sizeof *st);
    for(size_t k=0;k<6;k++)if(!isfinite(grid->affine[k]))return -1;
    for(size_t k=0;k<2;k++)if(!(grid->source_scale[k]>0) || !isfinite(grid->source_scale[k]))return -1;
    if(!(grid->affine[0]*grid->affine[4]-grid->affine[1]*grid->affine[3]>0))return -1;
    int pad=FIBER_MAP_PAD_CELLS;
    if(grid->rows>INT_MAX-2*pad || grid->cols>INT_MAX-2*pad)return -1;
    FiberMapResult map={0};map.rows=grid->rows+2*pad;map.cols=grid->cols+2*pad;
    size_t nn=(size_t)map.rows*(size_t)map.cols;
    if(nn>(size_t)INT_MAX || nv>(size_t)INT64_MAX || nf>(size_t)INT64_MAX ||
       nv>SIZE_MAX/2/sizeof(double) || nv>SIZE_MAX/sizeof(FmSample) || nf>SIZE_MAX/sizeof(FmMetric))return -1;
    map.row0=grid->row0-pad*grid->cell;map.col0=grid->col0-pad*grid->cell;map.cell=grid->cell;
    map.uv=ARENA_ALLOC(arena,2*nn*sizeof *map.uv);
    Arena_Mark keep=Arena_save(arena);
    double *base=ARENA_ALLOC(arena,2*nn*sizeof *base), *safe=ARENA_ALLOC(arena,2*(nv?nv:1)*sizeof *safe);
    uint8_t *fixed=ARENA_CALLOC(arena,nn,1), *seed=ARENA_CALLOC(arena,nn,1);
    uint8_t *reversed=ARENA_CALLOC(arena,nf?nf:1,1), *bad_face=ARENA_CALLOC(arena,nf?nf:1,1);
    FmSample *sample=ARENA_ALLOC(arena,(nv?nv:1)*sizeof *sample);
    FmMetric *metric=ARENA_ALLOC(arena,(nf?nf:1)*sizeof *metric);
    double baseline[6];memcpy(baseline,grid->affine,sizeof baseline);
    for(size_t v=0;v<nv;v++){
        double x=(canvas[2*v]-map.col0)/map.cell, y=(canvas[2*v+1]-map.row0)/map.cell;
        if(!isfinite(x) || !isfinite(y) || x<0 || y<0 || x>=(double)(map.cols-1) || y>=(double)(map.rows-1))goto done;
        int ix=(int)floor(x), iy=(int)floor(y);
        sample[v]=(FmSample){(size_t)iy*(size_t)map.cols+(size_t)ix,x-ix,y-iy};
        fm_affine(grid->affine,canvas[2*v],canvas[2*v+1],safe+2*v);
    }
    for(size_t f=0;f<nf;f++){
        for(size_t k=0;k<3;k++)if(faces[3*f+k]<0 || (size_t)faces[3*f+k]>=nv)goto done;
        double lo,hi;
        if(fm_reference(xyz,faces+3*f,metric+f)){
            fprintf(stderr,"[fiber_map] invalid original 3-D reference at face %zu\n",f);goto done;
        }
        fm_metric(metric+f,faces+3*f,canvas,1/grid->source_scale[0],1/grid->source_scale[1],&lo,&hi);
        if(!(lo>0) || !isfinite(lo) || !isfinite(hi) || lo>hi+ASM_METRIC_ROUNDOFF)goto done;
        metric[f].lower=fmin(FIBER_MAP_METRIC_MIN,lo);metric[f].upper=fmax(FIBER_MAP_METRIC_MAX,hi);
        st->source_metric_bad+=(lo<FIBER_MAP_METRIC_MIN-ASM_METRIC_ROUNDOFF || hi>FIBER_MAP_METRIC_MAX+ASM_METRIC_ROUNDOFF);
        fm_metric(metric+f,faces+3*f,safe,1,1,&lo,&hi);
        st->affine_metric_bad+=!fm_metric_preserved(metric+f,lo,hi);
    }
    if(st->affine_metric_bad){
        /* Polar rotation of the affine map in PHYSICAL source coordinates.
         * Keep its mapped grid centre; remove only the baseline's strain.
         * This is a safe origin for the solver, not an accepted fallback output. */
        double sx=grid->source_scale[0], sy=grid->source_scale[1];
        double cs=baseline[0]*sx+baseline[4]*sy, sn=baseline[3]*sx-baseline[1]*sy, norm=hypot(cs,sn);
        double cx=grid->col0+.5*(grid->cols-1)*grid->cell, cy=grid->row0+.5*(grid->rows-1)*grid->cell, centre[2];
        if(!(norm>0) || !isfinite(norm))goto done;
        fm_affine(baseline,cx,cy,centre);cs/=norm;sn/=norm;
        baseline[0]=cs/sx;baseline[1]=-sn/sy;baseline[3]=sn/sx;baseline[4]=cs/sy;
        baseline[2]=centre[0]-baseline[0]*cx-baseline[1]*cy;
        baseline[5]=centre[1]-baseline[3]*cx-baseline[4]*cy;
        st->rigid_baseline=1;
        for(size_t v=0;v<nv;v++)fm_affine(baseline,canvas[2*v],canvas[2*v+1],safe+2*v);
    }
    for(int y=0;y<map.rows;y++)for(int x=0;x<map.cols;x++){
        size_t q=(size_t)y*(size_t)map.cols+(size_t)x;
        fm_affine(baseline,map.col0+x*map.cell,map.row0+y*map.cell,base+2*q);
        map.uv[2*q]=base[2*q];map.uv[2*q+1]=base[2*q+1];
        if(y>=pad && y<map.rows-pad && x>=pad && x<map.cols-pad){
            size_t p=(size_t)(y-pad)*(size_t)grid->cols+(size_t)(x-pad);
            if(isfinite(grid->u[p]) && isfinite(grid->v[p])){
                fixed[q]=1;st->source_nodes++;map.uv[2*q]=grid->u[p];map.uv[2*q+1]=grid->v[p];
            }
        }
    }
    st->extended_nodes=nn-st->source_nodes;
    if(!st->source_nodes)goto done;
    for(size_t f=0;f<nf;f++){
        double lo,hi;fm_metric(metric+f,faces+3*f,safe,1,1,&lo,&hi);
        if(!fm_metric_preserved(metric+f,lo,hi)){
            fprintf(stderr,"[fiber_map] baseline metric failed at face %zu: %.17g %.17g, bounds %.17g %.17g\n",
                    f,lo,hi,metric[f].lower,metric[f].upper);goto done;
        }
        int32_t rev[3]={faces[3*f],faces[3*f+2],faces[3*f+1]};
        reversed[f]=(uint8_t)!UvGuard_positive(safe,faces+3*f);
        if(reversed[f] && !UvGuard_positive(safe,rev))goto done;
    }
    fprintf(stderr,"[orch] status Fibre map physical reference: %zu source failures, %zu affine failures; %s baseline\n",
            st->source_metric_bad,st->affine_metric_bad,st->rigid_baseline?"rigid":"affine");fflush(stderr);
    for(int pass=0;!grid->directions && pass<FIBER_MAP_REPAIR_PASSES;pass++){
        if(fm_extend(arena,&map,base,fixed))goto done;
        memset(seed,0,nn);
        size_t cells=fm_cells(&map,base,seed), bad=0, orientation_bad=0, metric_bad=0;
        int64_t vi=0, fi=0;
        #pragma omp parallel for schedule(static) if(nv > 4096)
        for(vi=0;vi<(int64_t)nv;vi++)fm_value(&map,sample+(size_t)vi,output+2*(size_t)vi);
        #pragma omp parallel for schedule(static) reduction(+:bad,orientation_bad,metric_bad) if(nf > 4096)
        for(fi=0;fi<(int64_t)nf;fi++){
            size_t f=(size_t)fi;
            int32_t idx[3]={faces[3*f],faces[3*f+1],faces[3*f+2]};
            if(reversed[f]){int32_t t=idx[1];idx[1]=idx[2];idx[2]=t;}
            int orientation=!UvGuard_interval(safe,output,idx,1);
            double lo,hi;fm_metric(metric+f,faces+3*f,output,1,1,&lo,&hi);
            int physical=!fm_metric_preserved(metric+f,lo,hi);
            bad_face[f]=(uint8_t)(orientation || physical);bad+=bad_face[f];
            orientation_bad+=orientation;metric_bad+=physical;
        }
        st->passes=pass+1;
        st->metric_faces_after=metric_bad;
        if(!pass){st->unsafe_cells_before=cells;st->unsafe_faces_before=orientation_bad;st->metric_faces_before=metric_bad;}
        fprintf(stderr,"[orch] status Continuous fibre map pass %d: %zu unsafe cells, %zu unsafe triangle paths, %zu metric failures; %zu released nodes\n",
                pass+1,cells,orientation_bad,metric_bad,st->relaxed_nodes);fflush(stderr);
        if(!cells && !bad){
            if(fm_harmonic_useful(grid,&map,base,st)){rc=0;break;}
            fprintf(stderr,"[fiber_map] feasible released map has no useful target gain (%.12g -> %.12g); fitting constrained map\n",
                    st->metric_fit_before,st->metric_fit_after);fflush(stderr);
            break;
        }
        for(size_t f=0;f<nf;f++)if(bad_face[f])for(size_t k=0;k<3;k++)
            fm_seed_cell(seed,sample[(size_t)faces[3*f+k]].node,(size_t)map.cols);
        size_t released=fm_release(arena,&map,fixed,seed);
        if(!released || released>=st->source_nodes-st->relaxed_nodes)break;
        st->relaxed_nodes+=released;
    }
    if(rc)rc=fm_metric_fit(arena,grid,&map,base,safe,sample,metric,reversed,faces,nv,nf,output,st);
    if(rc && !grid->directions && st->source_metric_bad)
        rc=fm_pinned_fit(arena,grid,&map,base,safe,sample,metric,reversed,faces,nv,nf,output,st);
    if(!rc){
        double sum=0;
        for(int y=0;y<grid->rows;y++)for(int x=0;x<grid->cols;x++){
            size_t p=(size_t)y*(size_t)grid->cols+(size_t)x;
            if(!isfinite(grid->u[p]) || !isfinite(grid->v[p]))continue;
            size_t q=(size_t)(y+pad)*(size_t)map.cols+(size_t)(x+pad);
            double d=hypot(map.uv[2*q]-grid->u[p],map.uv[2*q+1]-grid->v[p]);
            sum+=d*d;if(d>st->node_correction_max)st->node_correction_max=d;
        }
        st->node_correction_rms=sqrt(sum/(double)st->source_nodes);*result=map;
    }
done:
    st->seconds=ves_clock_sec()-started;Arena_restore(arena,keep);return rc;
}

/* Refine the field basis without inventing additional fibre observations.
 * The extra unobserved border keeps the FINAL padded domain identical: an
 * original map node (x,y) becomes (factor*x,factor*y) in the refined map.
 * Thus every coarse bilinear map remains representable by the finer field. */
static int fm_refine_grid(Arena_T arena, const FiberMapGrid *grid, int factor, FiberMapGrid *fine)
{
    int pad=FIBER_MAP_PAD_CELLS;
    if(!arena || !grid || !fine || !grid->u || !grid->v || factor<1 ||
       grid->rows<2 || grid->cols<2 || !(grid->cell>0) || !isfinite(grid->cell) ||
       !isfinite(grid->row0) || !isfinite(grid->col0))return -1;
    int64_t rows=(int64_t)grid->rows+2*pad,cols=(int64_t)grid->cols+2*pad;
    if(rows>INT_MAX || cols>INT_MAX)return -1;
    rows=(rows-1)*factor+1;cols=(cols-1)*factor+1;
    if(rows>INT_MAX || cols>INT_MAX || rows*cols>INT_MAX)return -1;
    *fine=*grid;fine->rows=(int)rows-2*pad;fine->cols=(int)cols-2*pad;
    fine->cell=grid->cell/factor;
    fine->row0=grid->row0-pad*(factor-1)*fine->cell;
    fine->col0=grid->col0-pad*(factor-1)*fine->cell;
    if(!(fine->cell>0) || !isfinite(fine->row0) || !isfinite(fine->col0))return -1;
    if(factor==1)return 0;
    size_t nn=(size_t)fine->rows*(size_t)fine->cols;
    if(nn>SIZE_MAX/sizeof(float))return -1;
    float *u=ARENA_ALLOC(arena,nn*sizeof *u),*v=ARENA_ALLOC(arena,nn*sizeof *v);
    for(size_t q=0;q<nn;q++)u[q]=v[q]=NAN;
    for(int y=0;y<grid->rows;y++)for(int x=0;x<grid->cols;x++){
        size_t p=(size_t)y*(size_t)grid->cols+(size_t)x;
        size_t q=(size_t)(factor*y+pad*(factor-1))*(size_t)fine->cols+(size_t)(factor*x+pad*(factor-1));
        u[q]=grid->u[p];v[q]=grid->v[p];
    }
    fine->u=u;fine->v=v;return 0;
}

int FiberMap_fit(Arena_T arena, const FiberMapGrid *grid, const float *xyz, const double *canvas, size_t nv,
                 const int32_t *faces, size_t nf, double *output, FiberMapResult *result, FiberMapStats *st)
{
    if(FIBER_MAP_PROJECT_REFINEMENT==1)return fm_fit_grid(arena,grid,xyz,canvas,nv,faces,nf,output,result,st);
    FiberMapGrid fine;
    if(!result || !st || fm_refine_grid(arena,grid,FIBER_MAP_PROJECT_REFINEMENT,&fine))return -1;
    fprintf(stderr,"[fiber_map] basis refinement %d: %d x %d -> %d x %d, cell %.9g -> %.9g; same observations and padded domain\n",
            FIBER_MAP_PROJECT_REFINEMENT,grid->rows,grid->cols,fine.rows,fine.cols,grid->cell,fine.cell);
    return fm_fit_grid(arena,&fine,xyz,canvas,nv,faces,nf,output,result,st);
}

#include "fiber_map_release_test.inc"
#include "fiber_map_metric_test.inc"
#include "fiber_map_direction_test.inc"

int FiberMap_selftest(void)
{
    Arena_T arena=Arena_new();int failures=fm_refine_selftest(arena)+fm_release_selftest(arena)+fm_metric_selftest(arena)+fm_direction_selftest(arena)+fm_harmonic_quality_selftest()+fm_polish_center_selftest()+fm_pinned_selftest(arena);
    enum { N=17, NV=N*N, NF=2*(N-1)*(N-1) };
    float *u=ARENA_ALLOC(arena,NV*sizeof *u), *v=ARENA_ALLOC(arena,NV*sizeof *v), *xyz=ARENA_ALLOC(arena,3*NV*sizeof *xyz);
    double *canvas=ARENA_ALLOC(arena,2*NV*sizeof *canvas), *out=ARENA_ALLOC(arena,2*NV*sizeof *out);
    int32_t *faces=ARENA_ALLOC(arena,3*NF*sizeof *faces);size_t nf=0;
    FiberMapGrid grid={N,N,0,0,4,u,v,{1,.125,7,0,1,-3},{1,1}};
    for(int y=0;y<N;y++)for(int x=0;x<N;x++){
        size_t q=(size_t)(y*N+x);double expected[2];
        canvas[2*q]=4*x;canvas[2*q+1]=4*y;fm_affine(grid.affine,4*x,4*y,expected);
        u[q]=(float)expected[0];v[q]=(float)expected[1];
        if(x+1<N && y+1<N){
            int32_t a=(int32_t)q,b=a+1,c=a+N,d=c+1;
            faces[3*nf]=a;faces[3*nf+1]=b;faces[3*nf+2]=d;nf++;
            faces[3*nf]=a;faces[3*nf+1]=d;faces[3*nf+2]=c;nf++;
        }
    }
    for(int trial=0;trial<8;trial++){
        Arena_Mark mark=Arena_save(arena);FiberMapResult map={0};FiberMapStats st={0};
        grid=(FiberMapGrid){N,N,0,0,4,u,v,{1,.125,7,0,1,-3},{1,1}};
        if(trial==5){grid.source_scale[0]=2;grid.source_scale[1]=.5;grid.affine[0]=.5;grid.affine[1]=0;grid.affine[4]=2;}
        if(trial>=6){grid.affine[0]=trial==6?.6:.9;grid.affine[1]=grid.affine[2]=grid.affine[5]=0;}
        for(size_t q=0;q<NV;q++){
            double expected[2];fm_affine(grid.affine,canvas[2*q],canvas[2*q+1],expected);
            u[q]=(float)(trial>=6?canvas[2*q]:expected[0]);v[q]=(float)(trial>=6?canvas[2*q+1]:expected[1]);
            xyz[3*q]=(float)(canvas[2*q]/grid.source_scale[0]*(trial==7?2:1));
            xyz[3*q+1]=(float)(canvas[2*q+1]/grid.source_scale[1]);xyz[3*q+2]=0;
        }
        if(trial==1)for(int y=0;y<N;y++)for(int x=6;x<=9;x++)u[y*N+x]=v[y*N+x]=NAN;
        if(trial==2){
            for(int y=0;y<N;y++)for(int x=0;x<N;x++){
                size_t q=(size_t)(y*N+x);double p[2];fm_affine(grid.affine,4*x,4*y,p);
                u[q]=(float)(p[0]+2*sin((double)y/3));v[q]=(float)(p[1]+sin((double)x/4));
            }
            v[8*N+8]-=30;
        }
        if(trial==3)for(size_t q=0;q<NV;q++){u[q]=(float)-canvas[2*q];v[q]=(float)-canvas[2*q+1];}
        if(trial==4)v[8*N+8]-=2.8f; /* positive everywhere, but much too compressed */
        int rc=FiberMap_fit(arena,&grid,xyz,canvas,NV,faces,nf,out,&map,&st);
        if(trial==3){
            /* The requested half-turn crosses a determinant zero. A constrained
             * fit may return a different, certified map, never that unsafe path. */
            double from[2*NV],target[2*NV];
            for(size_t q=0;q<NV;q++){
                fm_affine(grid.affine,canvas[2*q],canvas[2*q+1],from+2*q);
                target[2*q]=u[q];target[2*q+1]=v[q];
            }
            if(UvGuard_interval(from,target,faces,1))failures++;
            if(!rc){
                if(st.metric_faces_after)failures++;
                for(size_t f=0;f<nf;f++)if(!UvGuard_interval(from,out,faces+3*f,1)){failures++;break;}
            }
        }
        else{
            if(rc || !UvGuard_mesh(out,faces,nf) || st.metric_faces_after)failures++;
            if(trial<2 || trial>=5){
                if(st.relaxed_nodes)failures++;
                for(size_t q=0;q<NV;q++){
                    double expected[2];fm_affine(grid.affine,canvas[2*q],canvas[2*q+1],expected);
                    if(trial>=6){expected[0]=canvas[2*q];expected[1]=canvas[2*q+1];}
                    if(hypot(out[2*q]-expected[0],out[2*q+1]-expected[1])>1e-7){failures++;break;}
                }
                if(trial>=6 && (!st.rigid_baseline || st.affine_metric_bad!=nf))failures++;
                if(trial==7 && st.source_metric_bad!=nf)failures++;
            }else if((trial==2 && !st.unsafe_cells_before) || !st.relaxed_nodes || st.relaxed_nodes>NV/4 ||
                     fabs(out[0]-u[0])>1e-8 || fabs(out[1]-v[0])>1e-8)failures++;
            if(trial==4 && (!st.metric_faces_before || st.unsafe_faces_before))failures++;
        }
        Arena_restore(arena,mark);
    }
    FmMetric old_bad={1,0,1,.6,2};
    if(fm_metric_preserved(&old_bad,.5,1.5) || fm_metric_preserved(&old_bad,.6,20) ||
       !fm_metric_preserved(&old_bad,.65,1.8) || fm_metric_preserved(&old_bad,NAN,1))failures++;
    Arena_dispose(&arena);
    fprintf(stderr,"[fiber_map] continuous-map selftest: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    return failures;
}
