#include "scroll_source.h"
#include "../common/csr.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SS_TAU 6.283185307179586476925286766559

typedef struct { int32_t a, b, opposite; } SsBoundaryEdge;

static int ss_boundary_compare(const void *lhs, const void *rhs)
{
    const SsBoundaryEdge *a=lhs, *b=rhs;
    int32_t alo=a->a<a->b ? a->a : a->b, ahi=a->a<a->b ? a->b : a->a;
    int32_t blo=b->a<b->b ? b->a : b->b, bhi=b->a<b->b ? b->b : b->a;
    if (alo!=blo) return alo<blo ? -1 : 1;
    if (ahi!=bhi) return ahi<bhi ? -1 : 1;
    if (a->a!=b->a) return a->a<b->a ? -1 : 1;
    if (a->opposite!=b->opposite) return a->opposite<b->opposite ? -1 : 1;
    return 0;
}

static int ss_boundary_same(const SsBoundaryEdge *a, const SsBoundaryEdge *b)
{
    return (a->a==b->a && a->b==b->b) || (a->a==b->b && a->b==b->a);
}

int ScrollSource_boundary_directions(Arena_T arena, const MeshBinData *mesh,
                                      float **out_direction,
                                      ScrollSourceBoundaryReport *report)
{
    ScrollSourceBoundaryReport local={0};
    size_t nv=mesh ? mesh->nv : 0, nf=mesh ? mesh->nf : 0;
    float *answer=NULL;
    if (out_direction==NULL) return -1;
    *out_direction=NULL;
    if (report==NULL) report=&local;
    memset(report,0,sizeof *report);
    if (arena==NULL || mesh==NULL || nv>INT32_MAX ||
        nv>SIZE_MAX/(3*sizeof(double)) || nf>SIZE_MAX/(3*sizeof(SsBoundaryEdge)) ||
        (nv && mesh->verts==NULL) || (nf && mesh->faces==NULL)) return -1;
    if (!nv) return nf ? -1 : 0;
    for (size_t i=0;i<3*nv;i++) if (!isfinite(mesh->verts[i])) return -1;
    for (size_t i=0;i<3*nf;i++)
        if (mesh->faces[i]<0 || (size_t)mesh->faces[i]>=nv) return -1;
    answer=(float *)ARENA_CALLOC(arena,3*nv,sizeof(float));
    if (!nf) { *out_direction=answer; return 0; }
    Arena_Mark mark=Arena_save(arena);
    double *sum=(double *)ARENA_CALLOC(arena,3*nv,sizeof(double));
    uint8_t *invalid=(uint8_t *)ARENA_CALLOC(arena,nv,sizeof(uint8_t));
    uint8_t *degree=(uint8_t *)ARENA_CALLOC(arena,nv,sizeof(uint8_t));
    SsBoundaryEdge *edges=(SsBoundaryEdge *)ARENA_ALLOC(arena,3*nf*sizeof *edges);
    for (size_t f=0;f<nf;f++) for (int k=0;k<3;k++) {
        SsBoundaryEdge *e=edges+3*f+(size_t)k;
        e->a=mesh->faces[3*f+(size_t)k];
        e->b=mesh->faces[3*f+(size_t)((k+1)%3)];
        e->opposite=mesh->faces[3*f+(size_t)((k+2)%3)];
    }
    qsort(edges,3*nf,sizeof *edges,ss_boundary_compare);
    for (size_t first=0;first<3*nf;) {
        size_t end=first+1;
        const SsBoundaryEdge *e=edges+first;
        while (end<3*nf && ss_boundary_same(e,edges+end)) end++;
        if (e->a==e->b || end-first>2 ||
            (end-first==2 && edges[first+1].a==e->a)) {
            invalid[e->a]=invalid[e->b]=1;
            report->invalid_edges++;
        } else if (end-first==1) {
            double u[3]={0},v[3]={0},n[3]={0},c[3]={0},length=0;
            for (int k=0;k<3;k++) {
                u[k]=(double)mesh->verts[3*(size_t)e->b+k]-mesh->verts[3*(size_t)e->a+k];
                v[k]=(double)mesh->verts[3*(size_t)e->opposite+k]-mesh->verts[3*(size_t)e->a+k];
            }
            n[0]=u[1]*v[2]-u[2]*v[1]; n[1]=u[2]*v[0]-u[0]*v[2]; n[2]=u[0]*v[1]-u[1]*v[0];
            length=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
            report->boundary_edges++;
            if (degree[e->a]<3) degree[e->a]++;
            if (degree[e->b]<3) degree[e->b]++;
            if (!(length>0) || !isfinite(length)) {
                invalid[e->a]=invalid[e->b]=1;
                report->invalid_edges++;
            } else {
                for (int k=0;k<3;k++) n[k]/=length;
                /* edge x face-normal points away from the incident face. */
                c[0]=u[1]*n[2]-u[2]*n[1]; c[1]=u[2]*n[0]-u[0]*n[2]; c[2]=u[0]*n[1]-u[1]*n[0];
                for (int k=0;k<3;k++) { sum[3*(size_t)e->a+k]+=c[k]; sum[3*(size_t)e->b+k]+=c[k]; }
            }
        }
        first=end;
    }
    for (size_t i=0;i<nv;i++) {
        const double *s=sum+3*i;
        double length=sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]);
        if (degree[i] && degree[i]!=2) {
            report->irregular_boundary_vertices++;
            continue;
        }
        if (invalid[i] || !(length>0) || !isfinite(length)) continue;
        for (int k=0;k<3;k++) answer[3*i+k]=(float)(s[k]/length);
        report->boundary_vertices++;
    }
    Arena_restore(arena,mark);
    *out_direction=answer;
    return 0;
}

int ScrollSource_build(Arena_T arena, const MeshBinData *mesh,
                        const AxisWarp *axis, double axis_y, double axis_x,
                        int sense, double core_radius, ScrollSource *out)
{
    CSR_T adj = NULL;
    const int32_t *off = NULL, *tgt = NULL;
    int32_t *queue = NULL;
    size_t nv = 0, nf = 0;
    if (arena == NULL || mesh == NULL || out == NULL ||
        (sense != 1 && sense != -1) || !isfinite(axis_y) ||
        !isfinite(axis_x) || !isfinite(core_radius) || core_radius < 0.0)
        return -1;
    memset(out, 0, sizeof *out);
    nv = mesh->nv;
    nf = mesh->nf;
    if (nv == 0 && nf == 0) return 0;
    if (nv > INT32_MAX || nf > INT32_MAX / 6 || mesh->verts == NULL ||
        (nf && mesh->faces == NULL)) return -1;
    out->axial = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    out->radius = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    out->theta = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    out->q = (double *)ARENA_ALLOC(arena, nv * sizeof(double));
    out->straight = (float *)ARENA_ALLOC(arena, nv * 3 * sizeof(float));
    out->normal = (float *)ARENA_CALLOC(arena, nv * 3, sizeof(float));
    out->component = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    out->branch = (int32_t *)ARENA_CALLOC(arena, nv, sizeof(int32_t));
    out->component_size = (int32_t *)ARENA_CALLOC(arena, nv, sizeof(int32_t));
    out->core = (uint8_t *)ARENA_ALLOC(arena, nv);
    queue = (int32_t *)ARENA_ALLOC(arena, nv * sizeof(int32_t));
    for (size_t i = 0; i < nv; i++) {
        double ay = axis_y, ax = axis_x;
        const float *v = mesh->verts + 3*i;
        double y = 0.0, x = 0.0;
        if (!isfinite(v[0]) || !isfinite(v[1]) || !isfinite(v[2])) return -1;
        double axial=v[0];
        if (axis && axis->physical) {
            double world[3]={v[0],v[1],v[2]},metric[3]; uint32_t flags=0;
            if (AxisWarp_to_metric(axis,world,metric,&flags)!=0) return -1;
            axial=metric[0]; ay=v[1]-(metric[1]-axis->reference_y); ax=v[2]-(metric[2]-axis->reference_x);
        } else if (axis != NULL && AxisWarp_valid(axis)) AxisWarp_eval(axis, v[0], &ay, &ax);
        y = (double)v[1] - ay;
        x = (double)v[2] - ax;
        out->straight[3*i] = (float)axial;
        out->straight[3*i+1] = (float)(y + axis_y);
        out->straight[3*i+2] = (float)(x + axis_x);
        out->axial[i] = axial;
        out->radius[i] = hypot(y, x);
        out->theta[i] = atan2(y, x);
        out->q[i] = (double)sense * out->theta[i] / SS_TAU;
        out->component[i] = -1;
    }
    for (size_t f = 0; f < nf; f++) {
        int32_t a = mesh->faces[3*f], b = mesh->faces[3*f+1], c = mesh->faces[3*f+2];
        double u[3] = {0}, v[3] = {0}, n[3] = {0};
        if (a < 0 || b < 0 || c < 0 || (size_t)a >= nv ||
            (size_t)b >= nv || (size_t)c >= nv) return -1;
        for (int k = 0; k < 3; k++) {
            u[k] = (double)out->straight[3*(size_t)b+k] - out->straight[3*(size_t)a+k];
            v[k] = (double)out->straight[3*(size_t)c+k] - out->straight[3*(size_t)a+k];
        }
        n[0] = u[1]*v[2] - u[2]*v[1];
        n[1] = u[2]*v[0] - u[0]*v[2];
        n[2] = u[0]*v[1] - u[1]*v[0];
        for (int k = 0; k < 3; k++) {
            out->normal[3*(size_t)a+k] += (float)n[k];
            out->normal[3*(size_t)b+k] += (float)n[k];
            out->normal[3*(size_t)c+k] += (float)n[k];
        }
    }
    for (size_t i = 0; i < nv; i++) {
        float *n = out->normal + 3*i;
        double len = sqrt((double)n[0]*n[0] + (double)n[1]*n[1] + (double)n[2]*n[2]);
        if (len > 0.0) for (int k = 0; k < 3; k++) n[k] = (float)(n[k]/len);
    }
    adj = CSR_from_faces(arena, mesh->faces, nf, nv);
    off = CSR_offset(adj);
    tgt = CSR_target(adj);
    for (size_t s = 0; s < nv; s++) {
        size_t head = 0, tail = 0, contradictions = 0;
        int32_t chart = (int32_t)out->ncomponents;
        double max_radius = 0.0;
        if (out->component[s] >= 0) continue;
        out->component[s] = chart;
        queue[tail++] = (int32_t)s;
        while (head < tail) {
            int32_t a = queue[head++];
            if (out->radius[a] > max_radius) max_radius = out->radius[a];
            for (int32_t e = off[a]; e < off[a+1]; e++) {
                int32_t b = tgt[e];
                int delta = -(int)floor(out->q[b] - out->q[a] + 0.5);
                int64_t branch = (int64_t)out->branch[a] + delta;
                if (branch <= INT32_MIN || branch > INT32_MAX) return -1;
                if (out->component[b] < 0) {
                    out->component[b] = chart;
                    out->branch[b] = (int32_t)branch;
                    queue[tail++] = b;
                } else if (b > a && out->branch[b] != branch) {
                    contradictions++;
                }
            }
        }
        out->component_size[chart] = (int32_t)tail;
        out->core[chart] = max_radius < core_radius ? 1 : 0;
        if (!out->core[chart] && contradictions) {
            out->contradictory_edges += contradictions;
            out->core[chart] = 2;
        }
        if (contradictions)
            fprintf(stderr, "[scroll source] chart %d: %zu vertices, rmax=%.3f, %zu cycle residuals%s\n",
                    chart, tail, max_radius, contradictions,
                    out->core[chart] == 1 ? " (core curl)" : " (unresolved cyclic chart)");
        out->ncomponents++;
    }
    for (size_t i = 0; i < nv; i++) out->q[i] += out->branch[i];
    return 0;
}

int ScrollSource_selftest(void)
{
    Arena_T arena = Arena_new();
    ScrollSource s = {0}, repeat = {0};
    float v[] = {0, .1f, -10, 0, -.1f, -10, 1, .1f, -10,
                 1, -.1f, -10, 0, 0, 1};
    int32_t f[] = {0,1,2, 1,3,2};
    MeshBinData m = {v, NULL, f, 5, 2};
    int fail = ScrollSource_build(arena, &m, NULL, 0, 0, 1, 2, &s) != 0;
    if (!fail) fail += !(s.ncomponents == 2 && s.contradictory_edges == 0 &&
                        s.branch[1] == 1 && s.branch[3] == 1 &&
                        s.component_size[0] == 4 && s.core[1] && !s.core[0]);
    fail += ScrollSource_build(arena, &m, NULL, 0, 0, 1, 2, &repeat) != 0;
    if (!fail) fail += memcmp(s.q, repeat.q, 5*sizeof(double)) != 0;
    m.nv = m.nf = 0;
    fail += ScrollSource_build(arena, &m, NULL, 0, 0, 1, 2, &s) != 0;
    {
        float xyz[]={0,0,10, 4,0,10, 4,4,10, 0,4,10, 2,2,10};
        int32_t faces[]={0,1,4, 1,2,4, 2,3,4, 3,0,4};
        MeshBinData square={xyz,NULL,faces,5,4};
        ScrollSourceBoundaryReport report={0};
        float *direction=NULL,*flipped=NULL;
        fail+=ScrollSource_boundary_directions(arena,&square,&direction,&report)!=0;
        if (direction) {
            fail+=!(report.boundary_edges==4 && report.boundary_vertices==4 && !report.invalid_edges);
            fail+=!(direction[0]<-.7f && direction[1]<-.7f && direction[2]==0);
            fail+=!(direction[12]==0 && direction[13]==0 && direction[14]==0);
        }
        for (int i=0;i<4;i++) { int32_t t=faces[3*i]; faces[3*i]=faces[3*i+1]; faces[3*i+1]=t; }
        fail+=ScrollSource_boundary_directions(arena,&square,&flipped,&report)!=0;
        if (direction && flipped) fail+=memcmp(direction,flipped,15*sizeof(float))!=0;
        square.nv=square.nf=0;
        fail+=ScrollSource_boundary_directions(arena,&square,&flipped,&report)!=0 || flipped!=NULL;
        square.nv=5; square.nf=4; faces[0]=5;
        fail+=ScrollSource_boundary_directions(arena,&square,&flipped,&report)==0 || flipped!=NULL;
    }
    {
        /* Coincident coordinates do not make distinct source vertices or
         * source fronts identical. No coordinate welding is allowed here. */
        float xyz[]={0,0,10, 4,0,10, 4,4,10, 0,4,10,
                     0,0,10, 4,0,10, 4,4,10, 0,4,10};
        int32_t faces[]={0,1,2, 0,2,3, 4,5,6, 4,6,7};
        MeshBinData squares={xyz,NULL,faces,8,4};
        ScrollSourceBoundaryReport report={0};
        float *direction=NULL;
        fail+=ScrollSource_boundary_directions(arena,&squares,&direction,&report)!=0;
        if (direction) {
            fail+=!(report.boundary_edges==8 && report.boundary_vertices==8 &&
                    !report.invalid_edges && !report.irregular_boundary_vertices);
            fail+=memcmp(direction,direction+12,12*sizeof(float))!=0;
        }
        /* Pinching two otherwise disjoint fronts at one original vertex
         * must not manufacture an averaged continuation direction there. */
        faces[6]=0; faces[9]=0;
        fail+=ScrollSource_boundary_directions(arena,&squares,&direction,&report)!=0;
        if (direction) {
            fail+=!(report.boundary_edges==8 && report.boundary_vertices==6 &&
                    !report.invalid_edges && report.irregular_boundary_vertices==1);
            fail+=!(direction[0]==0 && direction[1]==0 && direction[2]==0 &&
                    direction[12]==0 && direction[13]==0 && direction[14]==0);
        }
    }
    {
        /* An inconsistently oriented shared edge cannot be used to propose
         * a material continuation through either of its endpoints. */
        float xyz[]={0,0,10, 4,0,10, 4,4,10, 0,4,10};
        int32_t faces[]={0,1,2, 0,3,2};
        MeshBinData square={xyz,NULL,faces,4,2};
        ScrollSourceBoundaryReport report={0};
        float *direction=NULL;
        fail+=ScrollSource_boundary_directions(arena,&square,&direction,&report)!=0;
        if (direction) {
            fail+=!(report.boundary_edges==4 && report.boundary_vertices==2 &&
                    report.invalid_edges==1 && !report.irregular_boundary_vertices);
            fail+=!(direction[0]==0 && direction[1]==0 && direction[2]==0 &&
                    direction[6]==0 && direction[7]==0 && direction[8]==0);
        }
    }
    fprintf(stderr, "[selftest] scroll source %s (%d failures)\n", fail ? "FAIL" : "PASS", fail);
    Arena_dispose(&arena);
    return fail ? -1 : 0;
}
