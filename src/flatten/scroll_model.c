#include "scroll_model.h"
#include "scroll_source.h"
#include "scroll_coordinate.h"
#include "winding_register.h"
#include "winding_mrf.h"
#include "quad_field.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"
#include "../common/pipeline_constants.h"

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define SM_TAU 6.283185307179586476925286766559
#define SM_MAGIC UINT64_C(0x3152435353455656)

typedef struct {
    uint64_t magic, version, source_hash, options_hash, nv, nf, nc;
} SmCubeHeader;

typedef struct {
    SmCubeHeader header;
    int32_t base;
} SmCube;

struct ScrollModel_T {
    ScrollModelOptions options;
    AxisWarp *owned_axis;
    char directory[1024];
    MeshPileEntry *pile;
    SmCube *cube;
    size_t npile, ncharts;
    uint64_t options_hash, model_hash;
    int64_t *winding;
    int32_t *material;
    int32_t profile_begin;
    size_t profile_count;
    double *profile_u, *profile_radius;
};

void ScrollModel_release(ScrollModel_T model)
{
    if (model && model->owned_axis) {
        AxisWarp_dispose(model->owned_axis); model->owned_axis=NULL; model->options.axis=NULL;
    }
}

uint64_t ScrollModel_coordinate_fingerprint(ScrollModel_T model)
{ return AxisWarp_fingerprint(model ? model->options.axis : NULL); }

static int sm_coordinate_point(ScrollModel_T m,const float source[3],double metric[3])
{
    const AxisWarp *axis=m->options.axis;
    double world[3]={source[0],source[1],source[2]}; uint32_t flags=0;
    if (AxisWarp_to_metric(axis,world,metric,&flags)!=0) return -1;
    if (AxisWarp_valid(axis)) {
        metric[1]+=m->options.axis_y-axis->reference_y;
        metric[2]+=m->options.axis_x-axis->reference_x;
    }
    return 0;
}

static uint64_t sm_hash(uint64_t h, const void *data, size_t n)
{
    const uint8_t *p = (const uint8_t *)data;
    for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= UINT64_C(1099511628211); }
    return h;
}

static int sm_directory(const char *path)
{
    char copy[1024];
    size_t n = strlen(path);
    if (n == 0 || n >= sizeof copy) return -1;
    memcpy(copy, path, n+1);
    for (size_t i = 1; i <= n; i++) {
        char c = copy[i];
        if (c != '/' && c != '\\' && c != 0) continue;
        if (i > 0 && copy[i-1] == ':') continue;
        copy[i] = 0;
        if (ves_mkdir(copy) != 0 && errno != EEXIST) return -1;
        copy[i] = c;
    }
    return 0;
}

static void sm_cube_path(ScrollModel_T m, size_t cube, char path[1200])
{
    snprintf(path, 1200, "%s/source/%s.chart", m->directory, m->pile[cube].cube_id);
}

static int sm_read_cube(Arena_T arena, ScrollModel_T m, size_t cube,
                         int32_t **component, int32_t **branch,
                         int32_t **sizes, uint8_t **core)
{
    char path[1200];
    FILE *f = NULL;
    SmCubeHeader h = {0};
    size_t nv = (size_t)m->cube[cube].header.nv, nc = (size_t)m->cube[cube].header.nc;
    int ok = 0;
    sm_cube_path(m, cube, path);
    f = fopen(path, "rb");
    if (f == NULL) return -1;
    *component = (int32_t *)ARENA_ALLOC(arena, nv*sizeof(int32_t));
    *branch = (int32_t *)ARENA_ALLOC(arena, nv*sizeof(int32_t));
    *sizes = (int32_t *)ARENA_ALLOC(arena, nc*sizeof(int32_t));
    *core = (uint8_t *)ARENA_ALLOC(arena, nc);
    ok = fread(&h, sizeof h, 1, f) == 1 &&
         memcmp(&h, &m->cube[cube].header, sizeof h) == 0 &&
         fread(*component, sizeof(int32_t), nv, f) == nv &&
         fread(*branch, sizeof(int32_t), nv, f) == nv &&
         fread(*sizes, sizeof(int32_t), nc, f) == nc &&
         fread(*core, 1, nc, f) == nc && fgetc(f) == EOF;
    if (fclose(f) != 0) ok = 0;
    if (!ok) return -1;
    for (size_t i = 0; i < nv; i++) if ((*component)[i] < 0 || (size_t)(*component)[i] >= nc) return -1;
    return 0;
}

static int sm_prepare_source(ScrollModel_T m)
{
    double start = ves_clock_sec();
    uint64_t total_nv = 0, total_nf = 0;
    char table_path[1200];
    FILE *table = NULL;
    snprintf(table_path, sizeof table_path, "%s/source_catalog.tsv", m->directory);
    table = fopen(table_path, "wb");
    if (table == NULL) return -1;
    fprintf(table, "cube\tvertices\tfaces\tchart_base\tcharts\tsource_hash\n");
    for (size_t i = 0; i < m->npile; i++) {
        Arena_T scratch = Arena_new();
        MeshBinData mesh = {0};
        ScrollSource source = {0};
        SmCubeHeader h = {0}, old = {0};
        char path[1200];
        FILE *f = NULL;
        int ok = 0, cached = 0;
        if (MeshBin_read_arena(scratch, m->pile[i].path, &mesh) != 0) goto cube_done;
        h.magic = SM_MAGIC;
        h.version = SCROLL_SOURCE_VERSION;
        h.options_hash = m->options_hash;
        h.nv = mesh.nv; h.nf = mesh.nf;
        h.source_hash = sm_hash(UINT64_C(14695981039346656037), mesh.verts, mesh.nv*3*sizeof(float));
        h.source_hash = sm_hash(h.source_hash, mesh.faces, mesh.nf*3*sizeof(int32_t));
        sm_cube_path(m, i, path);
        f = fopen(path, "rb");
        if (f != NULL) {
            cached = fread(&old, sizeof old, 1, f) == 1 && old.magic == h.magic &&
                     old.version == h.version && old.source_hash == h.source_hash &&
                     old.options_hash == h.options_hash && old.nv == h.nv && old.nf == h.nf &&
                     old.nc <= old.nv;
            if (cached) {
                long bytes = 0;
                if (fseek(f, 0, SEEK_END) != 0) cached = 0;
                bytes = ftell(f);
                if (bytes < 0 || (uint64_t)bytes != sizeof old + 8*old.nv + 5*old.nc) cached = 0;
            }
            fclose(f); f = NULL;
        }
        if (!cached) {
            if (ScrollSource_build(scratch, &mesh, m->options.axis,
                                   m->options.axis_y, m->options.axis_x,
                                   m->options.sense, m->options.core_radius, &source) != 0)
                goto cube_done;
            if (source.contradictory_edges) {
                fprintf(stderr, "[scroll-model] UNQUALIFIED %s: %zu inconsistent source-edge lifts; preserved as unresolved charts, not registered or filled\n",
                        m->pile[i].cube_id, source.contradictory_edges);
            }
            h.nc = source.ncomponents;
            f = fopen(path, "wb");
            if (f == NULL) goto cube_done;
            ok = fwrite(&h, sizeof h, 1, f) == 1 &&
                 fwrite(source.component, sizeof(int32_t), mesh.nv, f) == mesh.nv &&
                 fwrite(source.branch, sizeof(int32_t), mesh.nv, f) == mesh.nv &&
                 fwrite(source.component_size, sizeof(int32_t), source.ncomponents, f) == source.ncomponents &&
                 fwrite(source.core, 1, source.ncomponents, f) == source.ncomponents;
            if (fclose(f) != 0) ok = 0;
            f = NULL;
            if (!ok) goto cube_done;
        } else h = old;
        if (h.nc > INT32_MAX || m->ncharts > (size_t)INT32_MAX - h.nc) { ok = 0; goto cube_done; }
        m->cube[i].header = h;
        m->cube[i].base = (int32_t)m->ncharts;
        m->ncharts += (size_t)h.nc;
        m->model_hash = sm_hash(m->model_hash, &h, sizeof h);
        m->model_hash = sm_hash(m->model_hash, m->pile[i].cube_id, strlen(m->pile[i].cube_id));
        fprintf(table, "%s\t%llu\t%llu\t%d\t%llu\t%016llx\n", m->pile[i].cube_id,
                (unsigned long long)h.nv, (unsigned long long)h.nf, m->cube[i].base,
                (unsigned long long)h.nc, (unsigned long long)h.source_hash);
        total_nv += h.nv; total_nf += h.nf;
        ok = 1;
cube_done:
        if (f != NULL) fclose(f);
        Arena_dispose(&scratch);
        if (!ok) { fclose(table); return -1; }
        if (i % 50 == 0 || i+1 == m->npile)
            fprintf(stderr, "[scroll-model] source %zu/%zu: %zu charts, %.1fs%s\n",
                    i+1, m->npile, m->ncharts, ves_clock_sec()-start, cached ? " (validated cache)" : "");
    }
    if (fclose(table) != 0) return -1;
    fprintf(stderr, "[scroll-model] bounded source pass: %llu vertices / %llu faces -> %zu chart unknowns, %.2fs\n",
            (unsigned long long)total_nv, (unsigned long long)total_nf, m->ncharts, ves_clock_sec()-start);
    return 0;
}

/* An evidence region is fixed in the full source lattice. Repeated source
 * cubes always use the SAME local branch and global chart id; overlaps do
 * not create new unknowns or parent-boundary constraints. */
static int sm_collect_region(Arena_T arena, ScrollModel_T m,
                             const size_t *cubes, size_t ncubes,
                             WindingRegisterRelation **out, size_t *nout)
{
    size_t nv = 0, nc = 0, at = 0, cat = 0;
    float *v = NULL, *normal = NULL, *boundary_direction=NULL;
    size_t invalid_source_edges=0, irregular_boundary_vertices=0;
    double *z = NULL, *r = NULL, *theta = NULL, *q = NULL;
    int32_t *chart = NULL, *size = NULL, *global = NULL;
    WindingRegisterStats stats = {0};
    for (size_t i = 0; i < ncubes; i++) {
        nv += (size_t)m->cube[cubes[i]].header.nv;
        nc += (size_t)m->cube[cubes[i]].header.nc;
    }
    *out = NULL; *nout = 0;
    if (nv == 0 || nc == 0) return 0;
    if (nc > INT32_MAX) return -1;
    if (nv > SCROLL_MODEL_REGION_MAX_VERTICES) {
        long lo[3]={LONG_MAX,LONG_MAX,LONG_MAX}, hi[3]={LONG_MIN,LONG_MIN,LONG_MIN};
        int axis=-1;
        long span=0, cut=0;
        WindingRegisterRelation *part[2]={NULL,NULL};
        size_t counts[2]={0,0};
        for (size_t i=0; i<ncubes; i++) {
            const MeshPileEntry *p=m->pile+cubes[i];
            long xyz[3]={p->oz,p->oy,p->ox};
            for (int k=0; k<3; k++) { if (xyz[k]<lo[k]) lo[k]=xyz[k]; if (xyz[k]>hi[k]) hi[k]=xyz[k]; }
        }
        for (int k=0; k<3; k++) if (hi[k]-lo[k]>span) { axis=k; span=hi[k]-lo[k]; }
        if (axis<0 || span<2*m->options.chunk) {
            fprintf(stderr,"[scroll-model] minimal halo region exceeds bounded working set (%zu vertices); refused\n",nv);
            return -1;
        }
        cut=lo[axis]+(span/m->options.chunk/2)*m->options.chunk;
        fprintf(stderr,"[scroll-model] subdivide dense evidence region (%zu vertices), axis %d at %ld; one complete cube of overlap\n",nv,axis,cut);
        for (int side=0; side<2; side++) {
            Arena_T scratch=Arena_new();
            size_t subset[100], nsubset=0;
            WindingRegisterRelation *collected=NULL;
            int ok=0;
            for (size_t i=0; i<ncubes; i++) {
                const MeshPileEntry *p=m->pile+cubes[i];
                long xyz[3]={p->oz,p->oy,p->ox};
                if ((side==0 && xyz[axis]<=cut) || (side==1 && xyz[axis]>=cut)) subset[nsubset++]=cubes[i];
            }
            ok=sm_collect_region(scratch,m,subset,nsubset,&collected,&counts[side])==0;
            if (ok) {
                part[side]=(WindingRegisterRelation *)ARENA_ALLOC(arena,counts[side]*sizeof **part);
                if (counts[side]) memcpy(part[side],collected,counts[side]*sizeof **part);
            }
            Arena_dispose(&scratch);
            if (!ok) return -1;
        }
        *nout=counts[0]+counts[1];
        *out=(WindingRegisterRelation *)ARENA_ALLOC(arena,*nout*sizeof **out);
        if (counts[0]) memcpy(*out,part[0],counts[0]*sizeof **out);
        if (counts[1]) memcpy(*out+counts[0],part[1],counts[1]*sizeof **out);
        return 0;
    }
    v = (float *)ARENA_ALLOC(arena, 3*nv*sizeof(float));
    normal = (float *)ARENA_CALLOC(arena, 3*nv, sizeof(float));
    boundary_direction=(float *)ARENA_CALLOC(arena,3*nv,sizeof(float));
    z = (double *)ARENA_ALLOC(arena, nv*sizeof(double));
    r = (double *)ARENA_ALLOC(arena, nv*sizeof(double));
    theta = (double *)ARENA_ALLOC(arena, nv*sizeof(double));
    q = (double *)ARENA_ALLOC(arena, nv*sizeof(double));
    chart = (int32_t *)ARENA_ALLOC(arena, nv*sizeof(int32_t));
    size = (int32_t *)ARENA_ALLOC(arena, nc*sizeof(int32_t));
    global = (int32_t *)ARENA_ALLOC(arena, nc*sizeof(int32_t));
    for (size_t i = 0; i < ncubes; i++) {
        Arena_Mark mark = Arena_save(arena);
        size_t cube = cubes[i], local_nc = (size_t)m->cube[cube].header.nc;
        MeshBinData mesh = {0};
        MeshBinData framed={0};
        ScrollSourceBoundaryReport boundary_report={0};
        float *straight=NULL,*directions=NULL;
        int32_t *c = NULL, *b = NULL, *sizes = NULL, *map = NULL, *vertex = NULL;
        uint8_t *core = NULL;
        if (MeshBin_read_arena(arena, m->pile[cube].path, &mesh) != 0 ||
            sm_read_cube(arena, m, cube, &c, &b, &sizes, &core) != 0) return -1;
        map = (int32_t *)ARENA_ALLOC(arena, local_nc*sizeof(int32_t));
        vertex = (int32_t *)ARENA_ALLOC(arena, mesh.nv*sizeof(int32_t));
        straight=(float *)ARENA_ALLOC(arena,3*mesh.nv*sizeof(float));
        for (size_t j = 0; j < local_nc; j++) {
            map[j] = -1;
            if (core[j]) continue;
            map[j] = (int32_t)cat;
            global[cat] = m->cube[cube].base + (int32_t)j;
            size[cat++] = sizes[j];
        }
        for (size_t j = 0; j < mesh.nv; j++) {
            double metric[3],dy,dx;
            vertex[j]=-1;
            if (sm_coordinate_point(m,mesh.verts+3*j,metric)!=0) return -1;
            dy=metric[1]-m->options.axis_y; dx=metric[2]-m->options.axis_x;
            straight[3*j]=(float)metric[0];
            straight[3*j+1]=(float)(dy+m->options.axis_y);
            straight[3*j+2]=(float)(dx+m->options.axis_x);
            if (map[c[j]] < 0) continue;
            vertex[j]=(int32_t)at;
            z[at] = metric[0]; r[at] = hypot(dy, dx); theta[at] = atan2(dy, dx);
            q[at] = (double)m->options.sense * theta[at] / SM_TAU + b[j];
            v[3*at] = (float)metric[0];
            v[3*at+1] = (float)(dy + m->options.axis_y);
            v[3*at+2] = (float)(dx + m->options.axis_x);
            chart[at++] = map[c[j]];
        }
        framed=mesh; framed.verts=straight;
        if (ScrollSource_boundary_directions(arena,&framed,&directions,&boundary_report)!=0) return -1;
        invalid_source_edges+=boundary_report.invalid_edges;
        irregular_boundary_vertices+=boundary_report.irregular_boundary_vertices;
        for (size_t j=0;j<mesh.nv;j++) if (vertex[j]>=0)
            memcpy(boundary_direction+3*(size_t)vertex[j],directions+3*j,3*sizeof(float));
        for (size_t f=0; f<mesh.nf; f++) {
            int ids[3]={vertex[mesh.faces[3*f]],vertex[mesh.faces[3*f+1]],vertex[mesh.faces[3*f+2]]};
            double e[3]={0},g[3]={0},n[3]={0};
            if (ids[0]<0 || ids[1]<0 || ids[2]<0) continue;
            for (int k=0; k<3; k++) {
                e[k]=(double)v[3*(size_t)ids[1]+k]-v[3*(size_t)ids[0]+k];
                g[k]=(double)v[3*(size_t)ids[2]+k]-v[3*(size_t)ids[0]+k];
            }
            n[0]=e[1]*g[2]-e[2]*g[1]; n[1]=e[2]*g[0]-e[0]*g[2]; n[2]=e[0]*g[1]-e[1]*g[0];
            for (int j=0; j<3; j++) for (int k=0; k<3; k++) normal[3*(size_t)ids[j]+k]+=(float)n[k];
        }
        Arena_restore(arena, mark);
    }
    for (size_t i=0; i<at; i++) {
        float *n=normal+3*i;
        double len=sqrt((double)n[0]*n[0]+(double)n[1]*n[1]+(double)n[2]*n[2]);
        if (len>0) for (int k=0; k<3; k++) n[k]=(float)(n[k]/len);
    }
    if (at == 0 || cat < 2) return 0;
    if (WindingRegister_collect(arena, v, at, z, r, theta, q, chart, (int32_t)cat,
                               size, 0.0, m->options.pitch, m->options.sense,
                               NULL, normal, boundary_direction, out, nout, &stats) != 0) return -1;
    fprintf(stderr,"[scroll-model] original open fronts: %zu/%zu vertices; %zu non-facing proximity pairs; source invalid edges=%zu irregular boundary vertices=%zu\n",
            stats.continuation_front_vertices,at,stats.continuation_front_rejected,
            invalid_source_edges,irregular_boundary_vertices);
    if (stats.observations_dropped) {
        fprintf(stderr, "[scroll-model] evidence overflow: %zu observations dropped; model refused\n", stats.observations_dropped);
        return -1;
    }
    for (size_t i = 0; i < *nout; i++) {
        (*out)[i].a = global[(*out)[i].a];
        (*out)[i].b = global[(*out)[i].b];
    }
    return 0;
}

static int sm_relation_compare(const void *pa, const void *pb)
{
    const WindingRegisterRelation *a = (const WindingRegisterRelation *)pa;
    const WindingRegisterRelation *b = (const WindingRegisterRelation *)pb;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    if (a->b != b->b) return a->b < b->b ? -1 : 1;
    if (a->kind != b->kind) return a->kind < b->kind ? -1 : 1;
    if (a->target != b->target) return a->target < b->target ? -1 : 1;
    if (a->weight != b->weight) return a->weight > b->weight ? -1 : 1;
    if (a->mode_observations != b->mode_observations) return a->mode_observations > b->mode_observations ? -1 : 1;
    return 0;
}

static int sm_evidence(ScrollModel_T m, WindingRegisterRelation **out, size_t *nout)
{
    long lo[3] = {LONG_MAX, LONG_MAX, LONG_MAX}, hi[3] = {LONG_MIN, LONG_MIN, LONG_MIN};
    const long stride[3] = {3,4,4};
    WindingRegisterRelation *all = NULL;
    size_t nall = 0, capacity = 0, region = 0;
    double start = ves_clock_sec();
    for (size_t i = 0; i < m->npile; i++) {
        long p[3] = {m->pile[i].oz, m->pile[i].oy, m->pile[i].ox};
        for (int k = 0; k < 3; k++) { if (p[k] < lo[k]) lo[k] = p[k]; if (p[k] > hi[k]) hi[k] = p[k]; }
    }
    for (long iz = lo[0]; iz <= hi[0]; iz += stride[0]*m->options.chunk)
    for (long iy = lo[1]; iy <= hi[1]; iy += stride[1]*m->options.chunk)
    for (long ix = lo[2]; ix <= hi[2]; ix += stride[2]*m->options.chunk) {
        Arena_T scratch = Arena_new();
        size_t cubes[100], ncubes = 0, nr = 0;
        WindingRegisterRelation *r = NULL;
        char path[1200];
        FILE *f = NULL;
        uint64_t header[3] = {SM_MAGIC, m->model_hash, 0}, previous[3] = {0};
        int cached = 0, ok = 1;
        for (size_t i = 0; i < m->npile; i++) {
            const MeshPileEntry *p = m->pile+i;
            if (p->oz < iz || p->oz > iz+3*m->options.chunk ||
                p->oy < iy || p->oy > iy+4*m->options.chunk ||
                p->ox < ix || p->ox > ix+4*m->options.chunk) continue;
            if (ncubes == 100) { ok = 0; break; }
            cubes[ncubes++] = i;
        }
        if (!ncubes) { Arena_dispose(&scratch); continue; }
        snprintf(path, sizeof path, "%s/regions/z%ld_y%ld_x%ld.factors", m->directory, iz, iy, ix);
        f = fopen(path, "rb");
        if (f != NULL) {
            if (fread(previous, sizeof previous, 1, f) == 1 && previous[0] == header[0] &&
                previous[1] == header[1] && previous[2] < SCROLL_MODEL_MAX_FACTORS) {
                nr = (size_t)previous[2];
                r = (WindingRegisterRelation *)ARENA_ALLOC(scratch, nr*sizeof *r);
                cached = fread(r, sizeof *r, nr, f) == nr && fgetc(f) == EOF;
            }
            fclose(f); f = NULL;
        }
        if (ok && !cached) {
            ok = sm_collect_region(scratch, m, cubes, ncubes, &r, &nr) == 0;
            if (ok) {
                f = fopen(path, "wb");
                header[2] = nr;
                ok = f != NULL && fwrite(header, sizeof header, 1, f) == 1 && fwrite(r, sizeof *r, nr, f) == nr;
                if (f != NULL && fclose(f) != 0) ok = 0;
                f = NULL;
            }
        }
        if (ok && nr > SCROLL_MODEL_MAX_FACTORS-nall) ok = 0;
        if (ok && nall+nr > capacity) {
            size_t wanted = 2*(nall+nr) + 1024;
            WindingRegisterRelation *grown = (WindingRegisterRelation *)realloc(all, wanted*sizeof *all);
            if (grown == NULL) ok = 0;
            else { all = grown; capacity = wanted; }
        }
        if (ok && nr) { memcpy(all+nall, r, nr*sizeof *all); nall += nr; }
        Arena_dispose(&scratch);
        if (!ok) { free(all); return -1; }
        fprintf(stderr, "[scroll-model] region %zu (%ld,%ld,%ld): %zu cubes, %zu factors, %.1fs%s\n",
                ++region, iz, iy, ix, ncubes, nr, ves_clock_sec()-start, cached ? " (cache)" : "");
    }
    qsort(all, nall, sizeof *all, sm_relation_compare);
    /* The same source pair witnessed in two overlapping regions is not two
     * independent measurements. Keep the strongest factor for each target,
     * preserving contradictory targets and distinct evidence families. */
    {
        size_t w = 0;
        for (size_t i = 0; i < nall; ) {
            size_t j = i+1;
            while (j < nall && all[j].a == all[i].a && all[j].b == all[i].b &&
                   all[j].kind == all[i].kind && all[j].target == all[i].target) j++;
            all[w++] = all[i];
            i = j;
        }
        nall = w;
    }
    *out = all; *nout = nall;
    return 0;
}

static int sm_solve(Arena_T arena, ScrollModel_T m,
                    const WindingRegisterRelation *relations, size_t nr)
{
    ScrollCoordinateEdge *edges = (ScrollCoordinateEdge *)ARENA_ALLOC(arena, nr*sizeof *edges);
    ScrollCoordinateReport report = {0};
    UnionFind graph = UF_new(arena, (int32_t)m->ncharts);
    UnionFind material = UF_new(arena, (int32_t)m->ncharts);
    size_t ne = 0, bad_cont = 0, bad_order = 0;
    double energy = 0.0, start = ves_clock_sec();
    char path[1200];
    FILE *f = NULL;
    for (size_t i = 0; i < nr; i++) {
        const WindingRegisterRelation *r = relations+i;
        double w = r->kind == WINDING_REL_CONTINUATION
                 ? 16.0*fmax(r->weight-1000.0, .25) : fmax(r->weight, .25);
        if (!r->eligible) continue;
        edges[ne].a = r->a; edges[ne].b = r->b; edges[ne].delta = r->target;
        edges[ne].weight = w; edges[ne].exact = 0;
        uf_union(&graph, r->a, r->b);
        ne++;
    }
    if (ScrollCoordinate_solve(arena, m->ncharts, edges, ne, &m->winding, &report) != 0) return -1;
    fprintf(stderr, "[scroll-model] global all-edge least squares: %zu charts, %zu factors, %zu islands, %zu residuals, %.2fs\n",
            m->ncharts, ne, report.gauge_components, report.unsatisfied_soft_edges, ves_clock_sec()-start);
    /* Bounded label corrections refine the coarse real-valued solve. Every
     * factor survives every round; no spanning-tree UV is ever published. */
    for (int round = 0; round < SCROLL_MODEL_MRF_ROUNDS; round++) {
        Arena_Mark mark = Arena_save(arena);
        WindingMRFSite *sites = (WindingMRFSite *)ARENA_CALLOC(arena, m->ncharts, sizeof *sites);
        WindingMRFEdge *factors = (WindingMRFEdge *)ARENA_ALLOC(arena, ne*sizeof *factors);
        int32_t *labels = NULL;
        float *confidence = NULL;
        WindingMRFOptions options;
        WindingMRFStats stats = {0};
        WindingMRF_default_options(&options);
        options.label_min = -SCROLL_MODEL_MRF_WINDOW;
        options.label_max = SCROLL_MODEL_MRF_WINDOW;
        for (size_t i = 0; i < m->ncharts; i++) sites[i].fixed = uf_find(&graph, (int32_t)i) == (int32_t)i;
        for (size_t i = 0; i < ne; i++) {
            int64_t target = (int64_t)edges[i].delta - (m->winding[edges[i].b] - m->winding[edges[i].a]);
            if (target <= INT32_MIN || target > INT32_MAX) return -1;
            factors[i].a = edges[i].a; factors[i].b = edges[i].b;
            factors[i].target = (int32_t)target; factors[i].weight = edges[i].weight;
        }
        if (WindingMRF_solve(arena, sites, m->ncharts, factors, ne, &options,
                            &labels, &confidence, &stats) != 0) return -1;
        if (stats.energy_after > stats.energy_before + 1e-8*fmax(1.0, stats.energy_before)) return -1;
        for (size_t i = 0; i < m->ncharts; i++) m->winding[i] += labels[i];
        fprintf(stderr, "[scroll-model] global integer refinement %d: %zu changes, E %.9g -> %.9g\n",
                round+1, stats.changed_labels, stats.energy_before, stats.energy_after);
        Arena_restore(arena, mark);
        if (stats.changed_labels == 0) break;
    }
    snprintf(path, sizeof path, "%s/relation_residuals.tsv", m->directory);
    f = fopen(path, "wb");
    if (f == NULL) return -1;
    fprintf(f, "a\tb\tkind\ttarget\tsolved\teligible\tmode_observations\tagreement\n");
    for (size_t i = 0; i < nr; i++) {
        const WindingRegisterRelation *r = relations+i;
        int64_t solved = m->winding[r->b] - m->winding[r->a];
        fprintf(f, "%d\t%d\t%d\t%d\t%lld\t%d\t%zu\t%.9g\n", r->a, r->b, r->kind,
                r->target, (long long)solved, r->eligible, r->mode_observations, r->agreement);
        if (!r->eligible) continue;
        if (solved != r->target) {
            if (r->kind == WINDING_REL_CONTINUATION) bad_cont++; else bad_order++;
        } else if (r->kind == WINDING_REL_CONTINUATION) uf_union(&material, r->a, r->b);
        energy += fabs((double)solved-r->target);
    }
    if (fclose(f) != 0) return -1;
    m->material = (int32_t *)ARENA_ALLOC(arena, m->ncharts*sizeof(int32_t));
    for (size_t i = 0; i < m->ncharts; i++) m->material[i] = uf_find(&material, (int32_t)i);
    fprintf(stderr, "[scroll-model] PHYSICAL EVIDENCE (not a pass): residual continuation=%zu order=%zu, unweighted L1=%.0f; %.2fs\n",
            bad_cont, bad_order, energy, ves_clock_sec()-start);
    return 0;
}

/* Fixed full-source metric profile. Summation and integration have a stable
 * source/vertex order. Queries cannot alter bin edges, samples, or origin. */
static int sm_metric(Arena_T arena, ScrollModel_T m)
{
    int64_t low = INT64_MAX, high = INT64_MIN;
    uint64_t *count = NULL;
    for (size_t i = 0; i < m->ncharts; i++) {
        if (m->winding[i] < low) low = m->winding[i];
        if (m->winding[i] > high) high = m->winding[i];
    }
    low = (low-2)*SCROLL_MODEL_PHASE_BINS;
    high = (high+2)*SCROLL_MODEL_PHASE_BINS;
    if (low <= INT32_MIN || high > INT32_MAX || high-low > SCROLL_MODEL_MAX_PROFILE_BINS) return -1;
    m->profile_begin = (int32_t)low;
    m->profile_count = (size_t)(high-low+1);
    m->profile_radius = (double *)ARENA_CALLOC(arena, m->profile_count, sizeof(double));
    m->profile_u = (double *)ARENA_CALLOC(arena, m->profile_count, sizeof(double));
    count = (uint64_t *)ARENA_CALLOC(arena, m->profile_count, sizeof(uint64_t));
    for (size_t cube = 0; cube < m->npile; cube++) {
        Arena_Mark mark = Arena_save(arena);
        MeshBinData mesh = {0};
        int32_t *c = NULL, *b = NULL, *sizes = NULL;
        uint8_t *core = NULL;
        if (MeshBin_read_arena(arena, m->pile[cube].path, &mesh) != 0 ||
            sm_read_cube(arena, m, cube, &c, &b, &sizes, &core) != 0) return -1;
        for (size_t i = 0; i < mesh.nv; i++) {
            double metric[3],q=0.,radius=0.;
            int64_t index = 0;
            if (core[c[i]]) continue;
            if (sm_coordinate_point(m,mesh.verts+3*i,metric)!=0) return -1;
            radius = hypot(metric[1]-m->options.axis_y,metric[2]-m->options.axis_x);
            q = m->options.sense*atan2(metric[1]-m->options.axis_y,metric[2]-m->options.axis_x)/SM_TAU +
                b[i] + (double)m->winding[m->cube[cube].base+c[i]];
            index = (int64_t)floor(q*SCROLL_MODEL_PHASE_BINS)-m->profile_begin;
            if (index < 0 || (size_t)index >= m->profile_count) return -1;
            m->profile_radius[index] += radius;
            count[index]++;
        }
        Arena_restore(arena, mark);
    }
    {
        size_t first = m->profile_count;
        double previous = 0.0;
        for (size_t i = 0; i < m->profile_count; i++) if (count[i]) {
            m->profile_radius[i] /= (double)count[i];
            if (first == m->profile_count) first = i;
        }
        if (first == m->profile_count) return -1;
        previous = m->profile_radius[first];
        for (size_t i = 0; i < m->profile_count; i++) {
            if (count[i]) previous = m->profile_radius[i]; else m->profile_radius[i] = previous;
            if (i) m->profile_u[i] = m->profile_u[i-1] + SM_TAU*m->profile_radius[i-1]/SCROLL_MODEL_PHASE_BINS;
        }
    }
    return 0;
}

/* Versioned coefficient snapshot. All disk scalars have explicit widths;
 * records are zero-initialized, never copies of padded in-memory pointers.
 * This is a model artifact, not a parent certificate or a query cache. */
typedef struct {
    uint64_t magic,version,endian,npile,ncharts,nbins,naxis,model_hash,options_hash,checksum;
    int64_t profile_begin,sense,chunk;
    double axis_y,axis_x,pitch,core_radius;
} SmStoreHeader;

/* A physical-coordinate snapshot has an explicitly distinct version and an
 * additional checksummed payload. Legacy snapshot bytes remain unchanged. */
typedef struct {
    uint64_t version,has_normal;
    double source_voxel[3],source_origin[3],metric_voxel,radius,jacobian,ambiguity;
    double normal[3],reference_y,reference_x;
} SmPhysicalFrame;

typedef struct {
    char path[MESH_PILE_MAX_PATH],id[24];
    int64_t oz,oy,ox,base;
    SmCubeHeader source;
} SmStoreCube;

static int sm_io(FILE *f, void *data, size_t bytes, int writing, uint64_t *hash)
{
    if ((writing ? fwrite(data,1,bytes,f) : fread(data,1,bytes,f))!=bytes) return -1;
    *hash=sm_hash(*hash,data,bytes);
    return 0;
}

static int sm_store(Arena_T arena, ScrollModel_T m, int writing)
{
    SmStoreHeader h={0};
    SmPhysicalFrame frame={0};
    const AxisWarp *axis=m->options.axis;
    uint64_t hash=UINT64_C(14695981039346656037),expected=0;
    char path[1200];
    FILE *f=NULL;
    int rc=-1;
    if (writing) {
        h.magic=UINT64_C(0x314c45444f4d5356); h.version=SCROLL_MODEL_VERSION; h.endian=UINT64_C(0x01020304);
        if (axis && axis->physical) {
            const AxisWarpPhysical *c=&axis->physical_config;
            h.version=SCROLL_MODEL_VERSION+1000;
            frame.version=1; frame.has_normal=c->has_initial_normal;
            memcpy(frame.source_voxel,c->source_voxel_um_zyx,sizeof frame.source_voxel);
            memcpy(frame.source_origin,c->source_origin_um_zyx,sizeof frame.source_origin);
            memcpy(frame.normal,c->initial_normal_zyx,sizeof frame.normal);
            frame.metric_voxel=c->metric_voxel_um; frame.radius=c->maximum_radius_um;
            frame.jacobian=c->minimum_jacobian; frame.ambiguity=c->ambiguity_um;
            frame.reference_y=axis->reference_y; frame.reference_x=axis->reference_x;
        }
        h.npile=m->npile; h.ncharts=m->ncharts; h.nbins=m->profile_count;
        h.naxis=axis && AxisWarp_valid(axis) ? axis->n : 0;
        h.model_hash=m->model_hash; h.options_hash=m->options_hash;
        h.profile_begin=m->profile_begin; h.sense=m->options.sense; h.chunk=m->options.chunk;
        h.axis_y=m->options.axis_y; h.axis_x=m->options.axis_x;
        h.pitch=m->options.pitch; h.core_radius=m->options.core_radius;
    }
    snprintf(path,sizeof path,"%s/model.bin%s",m->directory,writing ? ".tmp" : "");
    f=fopen(path,writing ? "wb" : "rb");
    if (f==NULL) return -1;
    if (writing) { if (fwrite(&h,sizeof h,1,f)!=1) goto done; }
    else {
        AxisWarp *loaded=NULL;
        if (fread(&h,sizeof h,1,f)!=1 || h.magic!=UINT64_C(0x314c45444f4d5356) ||
            (h.version!=SCROLL_MODEL_VERSION && h.version!=SCROLL_MODEL_VERSION+1000) || h.endian!=UINT64_C(0x01020304) ||
            h.npile==0 || h.npile>SCROLL_MODEL_MAX_CUBES || h.ncharts==0 || h.ncharts>INT32_MAX ||
            h.nbins<2 || h.nbins>SCROLL_MODEL_MAX_PROFILE_BINS || h.naxis>1048576 ||
            h.profile_begin<INT32_MIN || h.profile_begin>INT32_MAX ||
            (h.sense!=-1 && h.sense!=1) || h.chunk<=0 || h.chunk>LONG_MAX ||
            !isfinite(h.axis_y) || !isfinite(h.axis_x) || !isfinite(h.pitch) || h.pitch<=0 ||
            !isfinite(h.core_radius) || h.core_radius<0) goto done;
        m->npile=(size_t)h.npile; m->ncharts=(size_t)h.ncharts; m->profile_count=(size_t)h.nbins;
        m->model_hash=h.model_hash; m->options_hash=h.options_hash; m->profile_begin=(int32_t)h.profile_begin;
        m->options.axis_y=h.axis_y; m->options.axis_x=h.axis_x;
        m->options.pitch=h.pitch; m->options.core_radius=h.core_radius;
        m->options.sense=(int)h.sense; m->options.chunk=(long)h.chunk;
        m->pile=(MeshPileEntry *)ARENA_CALLOC(arena,m->npile,sizeof *m->pile);
        m->cube=(SmCube *)ARENA_CALLOC(arena,m->npile,sizeof *m->cube);
        m->winding=(int64_t *)ARENA_ALLOC(arena,m->ncharts*sizeof(int64_t));
        m->material=(int32_t *)ARENA_ALLOC(arena,m->ncharts*sizeof(int32_t));
        m->profile_u=(double *)ARENA_ALLOC(arena,m->profile_count*sizeof(double));
        m->profile_radius=(double *)ARENA_ALLOC(arena,m->profile_count*sizeof(double));
        if (h.naxis) {
            loaded=(AxisWarp *)ARENA_CALLOC(arena,1,sizeof *loaded);
            loaded->n=(size_t)h.naxis; loaded->reference_y=h.axis_y; loaded->reference_x=h.axis_x;
            loaded->z=(double *)ARENA_ALLOC(arena,loaded->n*sizeof(double));
            loaded->y=(double *)ARENA_ALLOC(arena,loaded->n*sizeof(double));
            loaded->x=(double *)ARENA_ALLOC(arena,loaded->n*sizeof(double));
            m->options.axis=loaded; axis=loaded;
        }
        expected=h.checksum; h.checksum=0;
    }
    hash=sm_hash(hash,&h,sizeof h);
    for (size_t i=0,base=0; i<m->npile; i++) {
        SmStoreCube cube={0};
        MeshPileEntry *p=m->pile+i;
        if (writing) {
            memcpy(cube.path,p->path,strlen(p->path)+1);
            memcpy(cube.id,p->cube_id,strlen(p->cube_id)+1);
            cube.oz=p->oz; cube.oy=p->oy; cube.ox=p->ox;
            cube.base=m->cube[i].base; cube.source=m->cube[i].header;
        }
        if (sm_io(f,&cube,sizeof cube,writing,&hash)!=0) goto done;
        if (!writing) {
            long z=0,y=0,x=0;
            char id[24]={0};
            if (memchr(cube.path,0,sizeof cube.path)==NULL || memchr(cube.id,0,sizeof cube.id)==NULL ||
                cube.path[0]==0 || MeshPile_parse_cube_id(cube.id,id,&z,&y,&x)!=0 ||
                z!=cube.oz || y!=cube.oy || x!=cube.ox || cube.base!=(int64_t)base ||
                cube.source.magic!=SM_MAGIC || cube.source.version!=SCROLL_SOURCE_VERSION ||
                cube.source.options_hash!=h.options_hash || cube.source.nc>cube.source.nv ||
                cube.source.nv>INT32_MAX || cube.source.nf>INT32_MAX ||
                cube.source.nc>m->ncharts-base) goto done;
            memcpy(p->path,cube.path,sizeof p->path); memcpy(p->cube_id,cube.id,sizeof p->cube_id);
            p->oz=z; p->oy=y; p->ox=x; p->has_id=1;
            m->cube[i].base=(int32_t)base; m->cube[i].header=cube.source;
            base+=(size_t)cube.source.nc;
            if (i+1==m->npile && base!=m->ncharts) goto done;
        }
    }
    if (sm_io(f,m->winding,m->ncharts*sizeof(int64_t),writing,&hash)!=0 ||
        sm_io(f,m->material,m->ncharts*sizeof(int32_t),writing,&hash)!=0 ||
        sm_io(f,m->profile_u,m->profile_count*sizeof(double),writing,&hash)!=0 ||
        sm_io(f,m->profile_radius,m->profile_count*sizeof(double),writing,&hash)!=0) goto done;
    if (h.naxis && (sm_io(f,axis->z,(size_t)h.naxis*sizeof(double),writing,&hash)!=0 ||
                    sm_io(f,axis->y,(size_t)h.naxis*sizeof(double),writing,&hash)!=0 ||
                    sm_io(f,axis->x,(size_t)h.naxis*sizeof(double),writing,&hash)!=0)) goto done;
    if (h.version==SCROLL_MODEL_VERSION+1000 &&
        (h.naxis<3 || sm_io(f,&frame,sizeof frame,writing,&hash)!=0)) goto done;
    if (writing) {
        h.checksum=hash;
        if (fseek(f,0,SEEK_SET)!=0 || fwrite(&h,sizeof h,1,f)!=1) goto done;
    } else {
        if (hash!=expected || fgetc(f)!=EOF) goto done;
        for (size_t i=0; i<m->profile_count; i++)
            if (!isfinite(m->profile_u[i]) || !isfinite(m->profile_radius[i]) ||
                m->profile_radius[i]<=0 || (i && m->profile_u[i]<=m->profile_u[i-1])) goto done;
        for (size_t i=0; i<m->ncharts; i++)
            if (m->material[i]<0 || (size_t)m->material[i]>=m->ncharts ||
                m->winding[i]<INT32_MIN || m->winding[i]>INT32_MAX) goto done;
        if (h.naxis && !AxisWarp_valid(axis)) goto done;
        if (h.version==SCROLL_MODEL_VERSION+1000) {
            AxisWarpPhysical c={0}; AxisWarp fresh;
            AxisWarp_init(&fresh);
            if (frame.version!=1 || frame.has_normal>1) goto done;
            memcpy(c.source_voxel_um_zyx,frame.source_voxel,sizeof frame.source_voxel);
            memcpy(c.source_origin_um_zyx,frame.source_origin,sizeof frame.source_origin);
            memcpy(c.initial_normal_zyx,frame.normal,sizeof frame.normal);
            c.has_initial_normal=(int)frame.has_normal; c.metric_voxel_um=frame.metric_voxel;
            c.maximum_radius_um=frame.radius; c.minimum_jacobian=frame.jacobian; c.ambiguity_um=frame.ambiguity;
            double *points=(double*)ARENA_ALLOC(arena,(size_t)h.naxis*3*sizeof(double));
            for (size_t i=0;i<(size_t)h.naxis;i++) { points[3*i]=axis->z[i]; points[3*i+1]=axis->y[i]; points[3*i+2]=axis->x[i]; }
            if (AxisWarp_create_physical(&fresh,points,(size_t)h.naxis,&c,frame.reference_y,frame.reference_x)!=0) goto done;
            AxisWarp *owned=(AxisWarp*)ARENA_ALLOC(arena,sizeof *owned);
            *owned=fresh; m->options.axis=owned; m->owned_axis=owned;
        }
    }
    rc=0;
done:
    if (fclose(f)!=0) rc=-1;
    if (!writing && rc!=0) ScrollModel_release(m);
    if (writing && rc==0) {
        char destination[1200];
        snprintf(destination,sizeof destination,"%s/model.bin",m->directory);
#ifdef _WIN32
        if (!MoveFileExA(path,destination,MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH)) rc=-1;
#else
        if (rename(path,destination)!=0) rc=-1;
#endif
    }
    return rc;
}

int ScrollModel_load(Arena_T arena, const char *directory, ScrollModel_T *out)
{
    ScrollModel_T m=NULL;
    if (arena==NULL || directory==NULL || out==NULL || strlen(directory)>=1024) return -1;
    *out=NULL;
    m=(ScrollModel_T)ARENA_CALLOC(arena,1,sizeof *m);
    memcpy(m->directory,directory,strlen(directory)+1);
    if (sm_store(arena,m,0)!=0) return -1;
    *out=m;
    return 0;
}


int ScrollModel_build(Arena_T arena, const char *source_dir,
                       const char *model_dir, const ScrollModelOptions *options,
                       ScrollModel_T *out)
{
    ScrollModel_T m = NULL;
    WindingRegisterRelation *relations = NULL;
    size_t nr = 0;
    char path[1200];
    FILE *f = NULL;
    double start = ves_clock_sec();
    if (arena == NULL || source_dir == NULL || model_dir == NULL || options == NULL || out == NULL ||
        options->chunk <= 0 || !isfinite(options->pitch) || options->pitch <= 0.0 ||
        strlen(model_dir) >= 1024) return -1;
    *out = NULL;
    m = (ScrollModel_T)ARENA_CALLOC(arena, 1, sizeof *m);
    m->options = *options;
    memcpy(m->directory, model_dir, strlen(model_dir)+1);
    m->options_hash = UINT64_C(14695981039346656037);
    m->options_hash = sm_hash(m->options_hash, &options->axis_y, sizeof(double));
    m->options_hash = sm_hash(m->options_hash, &options->axis_x, sizeof(double));
    m->options_hash = sm_hash(m->options_hash, &options->pitch, sizeof(double));
    m->options_hash = sm_hash(m->options_hash, &options->core_radius, sizeof(double));
    m->options_hash = sm_hash(m->options_hash, &options->sense, sizeof(int));
    m->options_hash = sm_hash(m->options_hash, &options->chunk, sizeof(long));
    if (options->axis != NULL && AxisWarp_valid(options->axis)) {
        m->options_hash = sm_hash(m->options_hash, options->axis->z, options->axis->n*sizeof(double));
        m->options_hash = sm_hash(m->options_hash, options->axis->y, options->axis->n*sizeof(double));
        m->options_hash = sm_hash(m->options_hash, options->axis->x, options->axis->n*sizeof(double));
        if (options->axis->physical) {
            uint64_t coordinates=AxisWarp_fingerprint(options->axis);
            m->options_hash=sm_hash(m->options_hash,&coordinates,sizeof coordinates);
        }
    }
    {
        const uint64_t version=SCROLL_MODEL_VERSION;
        const uint64_t evidence_revision=SCROLL_MODEL_EVIDENCE_REVISION;
        m->model_hash=sm_hash(m->options_hash,&version,sizeof version);
        m->model_hash=sm_hash(m->model_hash,&evidence_revision,sizeof evidence_revision);
    }
    m->pile = (MeshPileEntry *)ARENA_ALLOC(arena, SCROLL_MODEL_MAX_CUBES*sizeof *m->pile);
    if (MeshPile_scan(source_dir, m->pile, SCROLL_MODEL_MAX_CUBES, &m->npile, NULL, NULL) != 0 || m->npile == 0) return -1;
    for (size_t i = 0; i < m->npile; i++) {
        if (!m->pile[i].has_id) return -1;
#ifdef _WIN32
        {
            char absolute[MESH_PILE_MAX_PATH];
            if (_fullpath(absolute,m->pile[i].path,sizeof absolute)==NULL) return -1;
            memcpy(m->pile[i].path,absolute,strlen(absolute)+1);
        }
#else
        {
            char *absolute=realpath(m->pile[i].path,NULL);
            if (absolute==NULL) return -1;
            if (strlen(absolute)>=sizeof m->pile[i].path) { free(absolute); return -1; }
            memcpy(m->pile[i].path,absolute,strlen(absolute)+1); free(absolute);
        }
#endif
    }
    m->cube = (SmCube *)ARENA_CALLOC(arena, m->npile, sizeof *m->cube);
    snprintf(path, sizeof path, "%s/source", model_dir);
    if (sm_directory(path) != 0) return -1;
    snprintf(path, sizeof path, "%s/regions", model_dir);
    if (sm_directory(path) != 0 || sm_prepare_source(m) != 0 || m->ncharts == 0) return -1;
    if (sm_evidence(m, &relations, &nr) != 0) return -1;
    if (sm_solve(arena, m, relations, nr) != 0) { free(relations); return -1; }
    free(relations);
    if (sm_metric(arena, m) != 0 || sm_store(arena,m,1)!=0) return -1;
    snprintf(path, sizeof path, "%s/model.json", model_dir);
    f = fopen(path, "wb");
    if (f == NULL) return -1;
    fprintf(f, "{\n  \"schema\": \"scroll-model-v%u\",\n  \"source_cubes\": %zu,\n  \"source_charts\": %zu,\n  \"model_fingerprint\": \"%016llx\",\n  \"profile_bins\": %zu,\n  \"elapsed_seconds\": %.6f,\n  \"physical_qualification\": \"NOT YET QUALIFIED; inspect relation residuals and geometry\"\n}\n",
            (unsigned)SCROLL_MODEL_VERSION, m->npile, m->ncharts, (unsigned long long)m->model_hash,
            m->profile_count, ves_clock_sec()-start);
    if (fclose(f) != 0) return -1;
    *out = m;
    return 0;
}

size_t ScrollModel_cube_count(ScrollModel_T model) { return model ? model->npile : 0; }
int ScrollModel_validate_catalog(Arena_T arena,ScrollModel_T m,const char *source_dir)
{
    Arena_Mark mark;
    MeshPileEntry *pile=NULL;
    size_t count=0;
    int ok=0;
    if (!arena || !m || !source_dir) return -1;
    mark=Arena_save(arena);
    pile=(MeshPileEntry *)ARENA_ALLOC(arena,SCROLL_MODEL_MAX_CUBES*sizeof *pile);
    ok=MeshPile_scan(source_dir,pile,SCROLL_MODEL_MAX_CUBES,&count,NULL,NULL)==0 && count==m->npile;
    for (size_t i=0;ok && i<count;i++) {
        if (!pile[i].has_id || strcmp(pile[i].cube_id,m->pile[i].cube_id)!=0) { ok=0; break; }
#ifdef _WIN32
        {
            char absolute[MESH_PILE_MAX_PATH];
            ok=_fullpath(absolute,pile[i].path,sizeof absolute)!=NULL &&
               _stricmp(absolute,m->pile[i].path)==0;
        }
#else
        {
            char *absolute=realpath(pile[i].path,NULL);
            ok=absolute!=NULL && strcmp(absolute,m->pile[i].path)==0;
            free(absolute);
        }
#endif
    }
    Arena_restore(arena,mark);
    if (!ok) fprintf(stderr,"[scroll model] input catalog does not match the fixed model; refused\n");
    return ok ? 0 : -1;
}
const MeshPileEntry *ScrollModel_cube(ScrollModel_T model, size_t cube)
{ return model != NULL && cube < model->npile ? model->pile+cube : NULL; }

int ScrollModel_evaluate_cube(Arena_T arena, ScrollModel_T m,
                              size_t cube, const MeshBinData *source_mesh,
                              float **out_uv, float **out_q,
                              int32_t **out_material)
{
    int32_t *component = NULL, *branch = NULL, *sizes = NULL;
    uint8_t *core = NULL;
    float *uv = NULL, *qout = NULL;
    int32_t *material = NULL;
    size_t nv=source_mesh ? source_mesh->nv : 0;
    const float *world=source_mesh ? source_mesh->verts : NULL;
    uint64_t hash=UINT64_C(14695981039346656037);
    if (arena==NULL || m == NULL || cube >= m->npile || source_mesh==NULL || nv != m->cube[cube].header.nv ||
        source_mesh->nf!=m->cube[cube].header.nf || (source_mesh->nf && source_mesh->faces==NULL) ||
        out_uv == NULL || out_q == NULL || out_material == NULL || (nv && world == NULL)) return -1;
    hash=sm_hash(hash,world,nv*3*sizeof(float));
    hash=sm_hash(hash,source_mesh->faces,source_mesh->nf*3*sizeof(int32_t));
    if (hash!=m->cube[cube].header.source_hash) return -1;
    uv = (float *)ARENA_ALLOC(arena, 2*nv*sizeof(float));
    qout = (float *)ARENA_ALLOC(arena, nv*sizeof(float));
    material = (int32_t *)ARENA_ALLOC(arena, nv*sizeof(int32_t));
    if (sm_read_cube(arena, m, cube, &component, &branch, &sizes, &core) != 0) return -1;
    for (size_t i = 0; i < nv; i++) {
        double metric[3],q=0.,t=0.;
        int32_t chart = m->cube[cube].base+component[i];
        int64_t bin = 0;
        if (sm_coordinate_point(m,world+3*i,metric)!=0) return -1;
        q = m->options.sense*atan2(metric[1]-m->options.axis_y,metric[2]-m->options.axis_x)/SM_TAU + branch[i] + (double)m->winding[chart];
        t = q*SCROLL_MODEL_PHASE_BINS - m->profile_begin;
        bin = (int64_t)floor(t);
        if (bin < 0 || (size_t)bin+1 >= m->profile_count) return -1;
        uv[2*i] = (float)(m->profile_u[bin] + (t-(double)bin)*(m->profile_u[bin+1]-m->profile_u[bin]));
        uv[2*i+1] = (float)metric[0];
        qout[i] = (float)q;
        material[i] = core[component[i]] ? -(int32_t)core[component[i]] : m->material[chart];
    }
    *out_uv = uv; *out_q = qout; *out_material = material;
    return 0;
}

static int sm_in_area(const MeshPileEntry *p, const long *b)
{
    return b==NULL || (p->oz>=b[0] && p->oz<=b[1] && p->oy>=b[2] &&
                      p->oy<=b[3] && p->ox>=b[4] && p->ox<=b[5]);
}

typedef struct { ScrollModel_T model; const char *directory; uint64_t fingerprint; } SmQuadContext;
static uint64_t sm_quad_fingerprint(ScrollModel_T m)
{
    uint64_t h=sm_hash(m->model_hash,m->winding,m->ncharts*sizeof(int64_t));
    h=sm_hash(h,m->material,m->ncharts*sizeof(int32_t));
    return sm_hash(h,m->profile_u,m->profile_count*sizeof(double));
}
static int sm_quad_reader(void *context,Arena_T arena,size_t chunk,MeshBinData *mesh,int32_t **chart)
{
    ScrollModel_T m=(ScrollModel_T)context;
    float *uv=NULL,*q=NULL;
    if (chunk>=m->npile || MeshBin_read_arena(arena,m->pile[chunk].path,mesh)!=0 ||
        ScrollModel_evaluate_cube(arena,m,chunk,mesh,&uv,&q,chart)!=0) return -1;
    mesh->uv=uv;
    return 0;
}
static int sm_quad_checkpoint(void *context,int level,QuadField_T field,const QuadFieldReport *report)
{
    const SmQuadContext *c=(const SmQuadContext *)context;
    Arena_T scratch=Arena_new();
    MeshBinData mesh={0};
    char dir[1200],path[1400];
    FILE *f=NULL;
    int rc=-1;
    snprintf(dir,sizeof dir,"%s/level_%02d",c->directory,level);
    if (sm_directory(dir)!=0 || QuadField_mesh(scratch,field,&mesh)!=0) goto done;
    snprintf(path,sizeof path,"%s/field.bin",dir);
    if (QuadField_save(field,path,c->fingerprint)!=0) goto done;
    snprintf(path,sizeof path,"%s/mesh_world.vmesh",dir);
    if (MeshBin_write(path,mesh.verts,mesh.nv,mesh.faces,mesh.nf,mesh.uv)!=0) goto done;
    snprintf(path,sizeof path,"%s/quad_field_report.json",dir);
    f=fopen(path,"wb"); if (!f) goto done;
    fprintf(f,"{\"schema\":\"global-quad-field-v1\",\"level\":%d,\"model_coefficients\":\"%016llx\","
              "\"cells\":%zu,\"nodes\":%zu,\"true_unknowns\":%zu,\"hanging_nodes\":%zu,"
              "\"samples\":%zu,\"sample_rms_vox\":%.9g,\"sample_max_vox\":%.9g,\"linear_relative_residual\":%.9g,"
              "\"query_independent\":true,\"dense_raster\":false,\"physical_qualification\":"
              "\"NOT QUALIFIED: observation-fit error is not a continuous surface error or topology certificate\"}\n",
            level,(unsigned long long)c->fingerprint,report->cells,report->nodes,report->unknowns,
            report->hanging_nodes,report->samples,report->rms_error,report->maximum_error,report->relative_residual);
    if (fclose(f)!=0) { f=NULL; goto done; } f=NULL; rc=0;
done:
    if (f) fclose(f); Arena_dispose(&scratch);
    return rc;
}
int ScrollModel_fit_quads(Arena_T arena,ScrollModel_T m,const char *directory)
{
    SmQuadContext context={0};
    QuadField_T final=NULL;
    if (!arena || !m || !directory) return -1;
    context.model=m; context.directory=directory;
    context.fingerprint=sm_quad_fingerprint(m);
    return QuadField_fit_stream(arena,m->npile,sm_quad_reader,m,SCROLL_QUAD_ROOT_CELLS,
                                SCROLL_QUAD_MIN_STEP,SCROLL_QUAD_FIT_ERROR,
                                sm_quad_checkpoint,&context,&final);
}

int ScrollModel_write_quad_area(Arena_T arena,ScrollModel_T m,const char *field_path,
                                 const long *bounds,const char *directory)
{
    QuadField_T field=NULL;
    MeshBinStream_T supported=NULL;
    uint8_t *selected=NULL;
    int32_t *node=NULL;
    MeshBinData mesh={0};
    char dir[1200],path[1400];
    FILE *f=NULL,*ids_file=NULL;
    size_t cells=0,samples=0,faces=0,reversed=0,collapsed=0;
    double started=ves_clock_sec(),maximum=0,square=0;
    int rc=-1;
    if (!arena || !m || !field_path || !directory ||
        (bounds && (bounds[0]>bounds[1] || bounds[2]>bounds[3] || bounds[4]>bounds[5]))) return -1;
    if (QuadField_load(arena,field_path,sm_quad_fingerprint(m),&field)!=0) return -1;
    selected=(uint8_t *)ARENA_CALLOC(arena,QuadField_cell_count(field),1);
    for (size_t chunk=0;chunk<m->npile;chunk++) if (sm_in_area(m->pile+chunk,bounds)) {
        Arena_T scratch=Arena_new();
        MeshBinData source={0};
        int32_t *chart=NULL;
        int ok=sm_quad_reader(m,scratch,chunk,&source,&chart)==0;
        for (size_t i=0;ok && i<source.nv;i++) if (chart[i]>=0) {
            size_t cell=0;
            ok=QuadField_locate(field,chart[i],source.uv[2*i],source.uv[2*i+1],&cell)==0;
            if (ok) { if (!selected[cell]) cells++; selected[cell]=1; samples++; }
        }
        for (size_t i=0;ok && i<source.nf;i++) {
            const int32_t *t=source.faces+3*i;
            for (int k=0;k<3;k++) if (t[k]<0 || (size_t)t[k]>=source.nv) ok=0;
            if (!ok) break;
            if (chart[t[0]]!=chart[t[1]] || chart[t[0]]!=chart[t[2]]) { ok=0; break; }
            if (chart[t[0]]>=0) faces++;
        }
        Arena_dispose(&scratch);
        if (!ok) return -1;
    }
    if (QuadField_mesh_selection(arena,field,selected,&mesh,&node)!=0) return -1;
    snprintf(dir,sizeof dir,"%s/stage3_fit",directory);
    if (sm_directory(dir)!=0) return -1;
    /* Keep control-domain visualization separate from the physical support.
     * Whole leaves are useful for solver inspection but are NOT a sheet. */
    snprintf(path,sizeof path,"%s/quad_domain_debug.vmesh",dir);
    if (MeshBin_write(path,mesh.verts,mesh.nv,mesh.faces,mesh.nf,mesh.uv)!=0) return -1;
    snprintf(path,sizeof path,"%s/quad_domain_global_node.i32",dir);
    f=fopen(path,"wb"); if (!f) return -1;
    {
        int ok=fwrite(node,sizeof(int32_t),mesh.nv,f)==mesh.nv;
        if (fclose(f)!=0) ok=0; f=NULL;
        if (!ok) return -1;
    }
    snprintf(path,sizeof path,"%s/fit_ribbon_world.vmesh",dir);
    if (MeshBin_stream_open(arena,path,samples,faces,1,&supported)!=0) goto done;
    snprintf(path,sizeof path,"%s/fit_source_id.u64",dir);
    ids_file=fopen(path,"wb"); if (!ids_file) goto done;
    for (size_t chunk=0;chunk<m->npile;chunk++) if (sm_in_area(m->pile+chunk,bounds)) {
        Arena_T scratch=Arena_new();
        MeshBinData source={0},part={0};
        QuadFieldProjectionReport report={0};
        int32_t *chart=NULL,*original=NULL;
        int ok=sm_quad_reader(m,scratch,chunk,&source,&chart)==0 &&
               QuadField_project_mesh(scratch,field,&source,chart,&part,&original,&report)==0;
        if (ok) ok=MeshBin_stream_append(supported,&part)==0;
        for (size_t begin=0;ok && begin<part.nv;) {
            uint64_t ids[4096]={0};
            size_t count=part.nv-begin;
            if (count>4096) count=4096;
            for (size_t i=0;i<count;i++) ids[i]=((uint64_t)chunk<<32)|(uint64_t)original[begin+i];
            ok=fwrite(ids,sizeof(uint64_t),count,ids_file)==count; begin+=count;
        }
        if (report.maximum_error>maximum) maximum=report.maximum_error;
        square+=report.rms_error*report.rms_error*report.vertices;
        reversed+=report.reversed_faces; collapsed+=report.collapsed_faces;
        Arena_dispose(&scratch);
        if (!ok) goto done;
    }
    if (fclose(ids_file)!=0) { ids_file=NULL; goto done; } ids_file=NULL;
    if (MeshBin_stream_close(&supported,1)!=0) goto done;
    snprintf(path,sizeof path,"%s/quad_query.json",dir);
    f=fopen(path,"wb"); if (!f) goto done;
    fprintf(f,"{\"schema\":\"global-quad-query-v2\",\"model_coefficients\":\"%016llx\","
              "\"cells\":%zu,\"vertices\":%zu,\"faces\":%zu,\"source_samples\":%zu,\"query_solve\":false,"
              "\"sample_max_vox\":%.9g,\"unweighted_sample_rms_vox\":%.9g,\"reversed_faces\":%zu,\"collapsed_faces\":%zu,"
              "\"elapsed_seconds\":%.6f,\"restriction\":\"fixed field evaluated at original source UVs, exact source connectivity; extras retained in stage2\","
              "\"surface_tessellation\":\"dense source PL; NOT yet a coarsened trimmed quad mesh\","
              "\"control_domain_debug\":{\"vertices\":%zu,\"faces\":%zu,\"untrimmed\":true},"
              "\"physical_qualification\":\"NOT QUALIFIED: source support preserved, but fitted embedding and alignment require validation\"}\n",
            (unsigned long long)sm_quad_fingerprint(m),cells,samples,faces,samples,maximum,
            samples ? sqrt(square/samples) : 0,reversed,collapsed,ves_clock_sec()-started,mesh.nv,mesh.nf);
    if (fclose(f)!=0) { f=NULL; goto done; } f=NULL;
    fprintf(stderr,"[quad query] %zu fixed leaves; supported surface %zu faces / %zu samples, rms %.6g max %.6g, reversed %zu collapsed %zu; no solve; %.3fs\n",
            cells,faces,samples,samples ? sqrt(square/samples) : 0,maximum,reversed,collapsed,ves_clock_sec()-started);
    rc=0;
done:
    if (supported) MeshBin_stream_close(&supported,0);
    if (f) fclose(f); if (ids_file) fclose(ids_file);
    return rc;
}

/* Core curls and unresolved source charts remain explicit geometry. They
 * must not be fed to a sheet bake with an arbitrary turn-zero placeholder:
 * that polluted the canvas bounds and created enormous blank sheets. */
static int sm_query_kind(int32_t material)
{ return material>=0 ? 0 : material==-1 ? 1 : material==-2 ? 2 : -1; }

static int sm_query_subset(Arena_T arena,const MeshBinData *input,
                            const int32_t *material,int kind,
                            MeshBinData *out,int32_t **original)
{
    int32_t *map=(int32_t *)ARENA_ALLOC(arena,input->nv*sizeof(int32_t));
    memset(out,0,sizeof *out);
    for (size_t i=0;i<input->nv;i++) {
        map[i]=-1;
        if (sm_query_kind(material[i])==kind) out->nv++;
    }
    out->verts=(float *)ARENA_ALLOC(arena,3*out->nv*sizeof(float));
    out->uv=(float *)ARENA_ALLOC(arena,2*out->nv*sizeof(float));
    *original=(int32_t *)ARENA_ALLOC(arena,out->nv*sizeof(int32_t));
    for (size_t i=0,at=0;i<input->nv;i++) if (sm_query_kind(material[i])==kind) {
        map[i]=(int32_t)at; (*original)[at]=(int32_t)i;
        memcpy(out->verts+3*at,input->verts+3*i,3*sizeof(float));
        memcpy(out->uv+2*at,input->uv+2*i,2*sizeof(float)); at++;
    }
    for (size_t i=0;i<input->nf;i++) {
        int32_t a=input->faces[3*i],b=input->faces[3*i+1],c=input->faces[3*i+2];
        if (sm_query_kind(material[a])!=sm_query_kind(material[b]) ||
            sm_query_kind(material[a])!=sm_query_kind(material[c])) return -1;
        if (map[a]>=0) out->nf++;
    }
    out->faces=(int32_t *)ARENA_ALLOC(arena,3*out->nf*sizeof(int32_t));
    for (size_t i=0,at=0;i<input->nf;i++) if (map[input->faces[3*i]]>=0) {
        for (int k=0;k<3;k++) out->faces[3*at+k]=map[input->faces[3*i+k]];
        at++;
    }
    return 0;
}

int ScrollModel_write_area(Arena_T arena, ScrollModel_T m,
                            const long *bounds, const char *directory)
{
    MeshBinStream_T source=NULL, world=NULL,extras[2]={NULL,NULL};
    FILE *qfile=NULL,*materialfile=NULL,*idfile=NULL,*catalog=NULL,*report=NULL;
    FILE *extra_ids[2]={NULL,NULL};
    const char *extra_name[2]={"core_curl_extras","unresolved_extras"};
    char path[2048];
    size_t nv=0,nf=0,ncubes=0,at_v=0,at_f=0,core_vertices=0,unresolved_vertices=0;
    size_t kind_nv[3]={0},kind_nf[3]={0};
    double started=ves_clock_sec();
    int rc=-1;
    if (arena==NULL || m==NULL || directory==NULL || strlen(directory)>1500) return -1;
    if (bounds && (bounds[0]>bounds[1] || bounds[2]>bounds[3] || bounds[4]>bounds[5])) return -1;
    for (size_t i=0; i<m->npile; i++) if (sm_in_area(m->pile+i,bounds)) {
        if (m->cube[i].header.nv>INT32_MAX-nv || m->cube[i].header.nf>SIZE_MAX-nf) return -1;
        nv+=(size_t)m->cube[i].header.nv; nf+=(size_t)m->cube[i].header.nf; ncubes++;
        {
            Arena_T scratch=Arena_new();
            MeshBinData mesh={0};
            int32_t *component=NULL,*branch=NULL,*sizes=NULL;
            uint8_t *core=NULL;
            int ok=MeshBin_read_arena(scratch,m->pile[i].path,&mesh)==0 &&
                   sm_read_cube(scratch,m,i,&component,&branch,&sizes,&core)==0;
            if (ok && (mesh.nv!=m->cube[i].header.nv || mesh.nf!=m->cube[i].header.nf)) ok=0;
            for (size_t j=0;ok && j<mesh.nv;j++) {
                int kind=core[component[j]];
                if (kind>2) { ok=0; break; }
                kind_nv[kind]++;
            }
            for (size_t j=0;ok && j<mesh.nf;j++) {
                int kind=core[component[mesh.faces[3*j]]];
                if (kind!=core[component[mesh.faces[3*j+1]]] || kind!=core[component[mesh.faces[3*j+2]]]) { ok=0; break; }
                kind_nf[kind]++;
            }
            Arena_dispose(&scratch);
            if (!ok) return -1;
        }
    }
    if (!nv || !nf) return -1;
    snprintf(path,sizeof path,"%s/stage1_mesh",directory);
    if (sm_directory(path)!=0) return -1;
    snprintf(path,sizeof path,"%s/stage2_unwrap",directory);
    if (sm_directory(path)!=0) return -1;
    snprintf(path,sizeof path,"%s/stage1_mesh/mesh_world.vmesh",directory);
    if (MeshBin_stream_open(arena,path,nv,nf,0,&source)!=0) goto done;
    snprintf(path,sizeof path,"%s/stage2_unwrap/unwrap_world.vmesh",directory);
    if (MeshBin_stream_open(arena,path,kind_nv[0],kind_nf[0],1,&world)!=0) goto done;
    for (int k=0;k<2;k++) if (kind_nv[k+1]) {
        snprintf(path,sizeof path,"%s/stage2_unwrap/%s.vmesh",directory,extra_name[k]);
        if (MeshBin_stream_open(arena,path,kind_nv[k+1],kind_nf[k+1],0,&extras[k])!=0) goto done;
        snprintf(path,sizeof path,"%s/stage2_unwrap/%s_source_id.u64",directory,extra_name[k]);
        extra_ids[k]=fopen(path,"wb");
        if (!extra_ids[k]) goto done;
    }
    snprintf(path,sizeof path,"%s/stage2_unwrap/unwrap_q.f32",directory); qfile=fopen(path,"wb");
    snprintf(path,sizeof path,"%s/stage2_unwrap/unwrap_material.i32",directory); materialfile=fopen(path,"wb");
    snprintf(path,sizeof path,"%s/stage2_unwrap/unwrap_source_id.u64",directory); idfile=fopen(path,"wb");
    snprintf(path,sizeof path,"%s/stage2_unwrap/source_catalog.tsv",directory); catalog=fopen(path,"wb");
    if (!qfile || !materialfile || !idfile || !catalog) goto done;
    fprintf(catalog,"cube_index\tcube\tfirst_vertex\tvertices\tfirst_face\tfaces\n");
    for (size_t i=0; i<m->npile; i++) if (sm_in_area(m->pile+i,bounds)) {
        Arena_T scratch=Arena_new();
        MeshBinData mesh={0};
        float *uv=NULL,*q=NULL;
        int32_t *material=NULL;
        uint64_t ids[4096]={0};
        int ok=MeshBin_read_arena(scratch,m->pile[i].path,&mesh)==0;
        if (ok) ok=ScrollModel_evaluate_cube(scratch,m,i,&mesh,&uv,&q,&material)==0;
        mesh.uv=uv;
        if (ok) ok=MeshBin_stream_append(source,&mesh)==0;
        for (int kind=0;ok && kind<3;kind++) {
            Arena_Mark mark=Arena_save(scratch);
            MeshBinData part={0};
            int32_t *original=NULL;
            ok=sm_query_subset(scratch,&mesh,material,kind,&part,&original)==0;
            if (ok && part.nv) ok=MeshBin_stream_append(kind==0 ? world : extras[kind-1],&part)==0;
            for (size_t j=0;ok && j<part.nv;) {
                float qq[4096]={0};
                int32_t mm[4096]={0};
                size_t n=part.nv-j;
                if (n>4096) n=4096;
                for (size_t k=0;k<n;k++) {
                    size_t vi=(size_t)original[j+k];
                    ids[k]=((uint64_t)i<<32)|(uint64_t)vi;
                    qq[k]=q[vi]; mm[k]=material[vi];
                }
                ok=fwrite(ids,sizeof(uint64_t),n,kind==0 ? idfile : extra_ids[kind-1])==n;
                if (ok && kind==0) ok=fwrite(qq,sizeof(float),n,qfile)==n && fwrite(mm,sizeof(int32_t),n,materialfile)==n;
                j+=n;
            }
            if (ok && kind==0) {
                fprintf(catalog,"%zu\t%s\t%zu\t%zu\t%zu\t%zu\n",i,m->pile[i].cube_id,at_v,part.nv,at_f,part.nf);
                at_v+=part.nv; at_f+=part.nf;
            }
            Arena_restore(scratch,mark);
        }
        Arena_dispose(&scratch);
        if (!ok) goto done;
    }
    if (ferror(qfile) || ferror(materialfile) || ferror(idfile) || ferror(catalog)) goto done;
    if (fclose(qfile)!=0) { qfile=NULL; goto done; } qfile=NULL;
    if (fclose(materialfile)!=0) { materialfile=NULL; goto done; } materialfile=NULL;
    if (fclose(idfile)!=0) { idfile=NULL; goto done; } idfile=NULL;
    if (fclose(catalog)!=0) { catalog=NULL; goto done; } catalog=NULL;
    if (MeshBin_stream_close(&source,1)!=0 || MeshBin_stream_close(&world,1)!=0) goto done;
    for (int k=0;k<2;k++) {
        if (extras[k] && MeshBin_stream_close(&extras[k],1)!=0) goto done;
        if (extra_ids[k]) {
            int result=fclose(extra_ids[k]); extra_ids[k]=NULL;
            if (result!=0) goto done;
        }
    }
    core_vertices=kind_nv[1]; unresolved_vertices=kind_nv[2];
    snprintf(path,sizeof path,"%s/stage2_unwrap/scroll_model_query.json",directory); report=fopen(path,"wb");
    if (!report) goto done;
    fprintf(report,"{\"schema\":\"scroll-model-query-v2\",\"model_fingerprint\":\"%016llx\","
                   "\"source_cubes\":%zu,\"vertices\":%zu,\"faces\":%zu,\"core_vertices\":%zu,"
                   "\"unresolved_source_vertices\":%zu,\"primary_vertices\":%zu,\"primary_faces\":%zu,"
                   "\"core_faces\":%zu,\"unresolved_faces\":%zu,\"query_solve\":false,\"coordinate_frame\":\"world_zyx\","
                   "\"source_identity\":\"(fixed_cube_index << 32) | vertex_index\",\"elapsed_seconds\":%.6f,"
                   "\"physical_qualification\":\"NOT QUALIFIED; source-coordinate inspection only\"}\n",
            (unsigned long long)m->model_hash,ncubes,nv,nf,core_vertices,unresolved_vertices,
            kind_nv[0],kind_nf[0],kind_nf[1],kind_nf[2],ves_clock_sec()-started);
    if (fclose(report)!=0) { report=NULL; goto done; } report=NULL;
    fprintf(stderr,"[scroll model query] %zu cubes, %zu vertices, %zu faces; one-cube working geometry, no solve; %.3fs\n",
            ncubes,nv,nf,ves_clock_sec()-started);
    rc=0;
done:
    if (source) MeshBin_stream_close(&source,0);
    if (world) MeshBin_stream_close(&world,0);
    for (int k=0;k<2;k++) {
        if (extras[k]) MeshBin_stream_close(&extras[k],0);
        if (extra_ids[k]) fclose(extra_ids[k]);
    }
    if (qfile) fclose(qfile); if (materialfile) fclose(materialfile);
    if (idfile) fclose(idfile); if (catalog) fclose(catalog); if (report) fclose(report);
    return rc;
}
