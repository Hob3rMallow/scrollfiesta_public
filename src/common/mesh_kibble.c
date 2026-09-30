#include "mesh_kibble.h"
#include "union_find.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

int MeshKibble_filter(Arena_T arena, const MeshBinData *in, double min_area,
                       MeshBinData *out, MeshKibbleStats *stats)
{
    UnionFind uf;
    double *area;
    int32_t *map;
    size_t v, f;
    if (!arena || !in || !out || !stats || !isfinite(min_area) || min_area < 0 ||
        in->nv > INT32_MAX || (in->nv && !in->verts) || (in->nf && !in->faces))
        return -1;
    memset(out, 0, sizeof *out);
    memset(stats, 0, sizeof *stats);
    stats->input_vertices = in->nv; stats->input_faces = in->nf;
    if (min_area == 0) {
        *out = *in;
        stats->kept_vertices = in->nv; stats->kept_faces = in->nf;
        return 0;
    }
    if (!in->nv) return in->nf ? -1 : 0;
    uf = UF_new(arena, (int32_t)in->nv);
    area = ARENA_ALLOC(arena, in->nv * sizeof *area);
    map = ARENA_ALLOC(arena, in->nv * sizeof *map);
    for (v = 0; v < in->nv; ++v) {
        area[v] = -1;
        for (int k = 0; k < 3; ++k)
            if (!isfinite(in->verts[3*v+k])) return -1;
    }
    for (f = 0; f < in->nf; ++f) {
        for (int k = 0; k < 3; ++k)
            if (in->faces[3*f+k] < 0 || (size_t)in->faces[3*f+k] >= in->nv)
                return -1;
        uf_union(&uf, in->faces[3*f], in->faces[3*f+1]);
        uf_union(&uf, in->faces[3*f], in->faces[3*f+2]);
    }
    for (f = 0; f < in->nf; ++f) {
        const float *a = in->verts + 3*(size_t)in->faces[3*f];
        const float *b = in->verts + 3*(size_t)in->faces[3*f+1];
        const float *c = in->verts + 3*(size_t)in->faces[3*f+2];
        double e[3], d[3], cross[3], ar;
        int32_t r = uf_find(&uf, in->faces[3*f]);
        for (int k = 0; k < 3; ++k) { e[k] = (double)b[k]-a[k]; d[k] = (double)c[k]-a[k]; }
        cross[0] = e[1]*d[2]-e[2]*d[1];
        cross[1] = e[2]*d[0]-e[0]*d[2];
        cross[2] = e[0]*d[1]-e[1]*d[0];
        ar = 0.5*sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
        if (area[r] < 0) area[r] = 0;
        area[r] += ar;
    }
    for (v = 0; v < in->nv; ++v) {
        if (area[v] >= 0) {
            ++stats->components; stats->input_area += area[v];
            if (area[v] < min_area) {
                ++stats->removed_components; stats->removed_area += area[v];
            }
        }
        int32_t r = uf_find(&uf, (int32_t)v);
        if (area[r] < 0) ++stats->orphan_vertices;
        map[v] = area[r] >= min_area ? (int32_t)out->nv++ : -1;
    }
    for (f = 0; f < in->nf; ++f) if (map[in->faces[3*f]] >= 0) ++out->nf;
    stats->kept_vertices = out->nv; stats->kept_faces = out->nf;
    if (!out->nf) return 0;
    out->verts = ARENA_ALLOC(arena, out->nv * 3 * sizeof *out->verts);
    out->faces = ARENA_ALLOC(arena, out->nf * 3 * sizeof *out->faces);
    if (in->uv) out->uv = ARENA_ALLOC(arena, out->nv * 2 * sizeof *out->uv);
    for (v = 0; v < in->nv; ++v) if (map[v] >= 0) {
        memcpy(out->verts + 3*(size_t)map[v], in->verts + 3*v, 3*sizeof(float));
        if (in->uv) memcpy(out->uv + 2*(size_t)map[v], in->uv + 2*v, 2*sizeof(float));
    }
    size_t at = 0;
    for (f = 0; f < in->nf; ++f) if (map[in->faces[3*f]] >= 0) {
        for (int k = 0; k < 3; ++k) out->faces[3*at+k] = map[in->faces[3*f+k]];
        ++at;
    }
    return 0;
}

int MeshKibble_selftest(void)
{
    /* Two real disconnected sheets, an interleaved dust triangle, an orphan.
     * Face counts are equal: only physical area distinguishes dust. */
    float xyz[] = {0,0,0, 0,0,100, 1,1,1, 0,100,0, 1,1,2,
                   10,0,0, 10,0,80, 1,2,1, 10,80,0, 9999,9999,9999};
    float uv[20];
    int32_t faces[] = {2,4,7, 0,1,3, 5,6,8};
    int32_t expected[] = {0,1,2, 3,4,5};
    MeshBinData in = {xyz, uv, faces, 10, 3}, out, empty = {0};
    MeshKibbleStats st;
    Arena_T a = Arena_new();
    int fail = 0;
    if (!a) return 1;
    for (int i=0; i<20; ++i) uv[i] = (float)i;
    if (MeshKibble_filter(a,&in,1024,&out,&st) || out.nv!=6 || out.nf!=2 ||
        st.components!=3 || st.removed_components!=1 || st.orphan_vertices!=1 ||
        st.removed_area!=0.5 || memcmp(out.faces,expected,sizeof expected) ||
        memcmp(out.verts+6,xyz+9,3*sizeof(float)) || out.uv[4]!=uv[6]) ++fail;
    if (MeshKibble_filter(a,&in,3200,&out,&st) || out.nf!=2) ++fail;
    if (MeshKibble_filter(a,&in,5001,&out,&st) || out.nv || out.nf || st.removed_components!=3) ++fail;
    if (MeshKibble_filter(a,&in,0,&out,&st) || out.verts!=xyz || out.faces!=faces || out.nv!=10) ++fail;
    if (MeshKibble_filter(a,&empty,1024,&out,&st) || out.nv || out.nf) ++fail;
    if (!MeshKibble_filter(a,&in,-1,&out,&st)) ++fail;
    xyz[0] = (float)NAN;
    if (!MeshKibble_filter(a,&in,1024,&out,&st)) ++fail;
    Arena_dispose(&a);
    fprintf(stderr,"[selftest] mesh kibble %s (%d failures)\n",fail?"FAIL":"PASS",fail);
    return fail;
}
