/*
 * component_cull.c -- connectivity pass + surface-area filter. See header.
 *
 * After hole-fill/QEM a cube can carry small disconnected "kibble" (stray BPA
 * islands, cut-zone crumbs). This splits every mesh into its connectivity-
 * components and drops any whose area is below a fraction of the total meshed
 * area, leaving only the real sheets.
 */
#include "component_cull.h"

#include "../common/union_find.h"
#include "../common/pca.h"

#include <string.h>
#include <math.h>
#include <assert.h>

static double tri_area(const float *v, int32_t a, int32_t b, int32_t c)
{
    double e1x = (double)v[b*3+0]-v[a*3+0], e1y = (double)v[b*3+1]-v[a*3+1], e1z = (double)v[b*3+2]-v[a*3+2];
    double e2x = (double)v[c*3+0]-v[a*3+0], e2y = (double)v[c*3+1]-v[a*3+1], e2z = (double)v[c*3+2]-v[a*3+2];
    double cx = e1y*e2z - e1z*e2y, cy = e1z*e2x - e1x*e2z, cz = e1x*e2y - e1y*e2x;
    return 0.5 * sqrt(cx*cx + cy*cy + cz*cz);
}

/* Build a ComponentMesh from the faces of `src` whose verts are in UF root
 * `root`; also return its triangle-area sum. Returns 0 on success. */
static int extract_cc(Arena_T arena, const ComponentMesh *src,
                      UnionFind *uf, int32_t root,
                      ComponentMesh *out, double *out_area)
{
    size_t nf = 0;
    for (size_t f = 0; f < src->nf; f++)
        if (uf_find(uf, src->faces[f*3]) == root) nf++;
    if (nf == 0) return -1;

    int32_t *faces = (int32_t *)ARENA_ALLOC(arena, (long)(nf*3*sizeof(int32_t)));
    size_t k = 0;
    for (size_t f = 0; f < src->nf; f++)
        if (uf_find(uf, src->faces[f*3]) == root) {
            faces[k*3+0]=src->faces[f*3+0];
            faces[k*3+1]=src->faces[f*3+1];
            faces[k*3+2]=src->faces[f*3+2];
            k++;
        }

    int32_t *remap = (int32_t *)ARENA_CALLOC(arena, (long)src->nv, (long)sizeof(int32_t));
    for (size_t i = 0; i < nf*3; i++) remap[faces[i]] = 1;
    int32_t cnt = 0;
    for (size_t i = 0; i < src->nv; i++) remap[i] = remap[i] ? cnt++ : -1;
    if (cnt < 3) return -1;
    float *v = (float *)ARENA_ALLOC(arena, (long)((size_t)cnt*3*sizeof(float)));
    for (size_t i = 0; i < src->nv; i++) if (remap[i] >= 0) {
        v[remap[i]*3+0]=src->verts[i*3+0];
        v[remap[i]*3+1]=src->verts[i*3+1];
        v[remap[i]*3+2]=src->verts[i*3+2];
    }
    for (size_t i = 0; i < nf*3; i++) faces[i] = remap[faces[i]];

    double area = 0.0;
    for (size_t f = 0; f < nf; f++)
        area += tri_area(v, faces[f*3+0], faces[f*3+1], faces[f*3+2]);

    float pca[3] = {0,0,1}, cen[3] = {0,0,0};
    PCA_normal(v, (size_t)cnt, pca, cen);
    memset(out, 0, sizeof(*out));
    out->verts = v; out->faces = faces; out->nv = (size_t)cnt; out->nf = nf;
    out->pca_normal[0]=pca[0]; out->pca_normal[1]=pca[1]; out->pca_normal[2]=pca[2];
    out->centroid[0]=cen[0]; out->centroid[1]=cen[1]; out->centroid[2]=cen[2];
    out->self = out;
    *out_area = area;
    return 0;
}

int ComponentCull_by_area(Arena_T arena,
                          const ComponentMesh *in, size_t n_in,
                          float min_frac,
                          ComponentMesh **out, size_t *n_out,
                          size_t *n_culled)
{
    assert(arena && out && n_out);
    *out = NULL; *n_out = 0;
    if (n_culled) *n_culled = 0;
    if (n_in == 0) return 0;

    /* split every input into connectivity-components; collect each with its area */
    size_t cap = n_in + 8, cnt = 0;
    ComponentMesh *ccs   = (ComponentMesh *)ARENA_ALLOC(arena, (long)(cap*sizeof(ComponentMesh)));
    double        *areas = (double *)ARENA_ALLOC(arena, (long)(cap*sizeof(double)));
    size_t        *parent = (size_t *)ARENA_ALLOC(arena, (long)(cap*sizeof(size_t)));
    double *parent_total = (double *)ARENA_CALLOC(arena, (long)n_in,
                                                   (long)sizeof(double));

    for (size_t i = 0; i < n_in; i++) {
        const ComponentMesh *cm = &in[i];
        if (cm->nf == 0 || cm->nv == 0) continue;
        UnionFind uf = UF_new(arena, (int32_t)cm->nv);
        for (size_t f = 0; f < cm->nf; f++) {
            uf_union(&uf, cm->faces[f*3+0], cm->faces[f*3+1]);
            uf_union(&uf, cm->faces[f*3+1], cm->faces[f*3+2]);
            uf_union(&uf, cm->faces[f*3+0], cm->faces[f*3+2]);
        }
        /* roots that carry at least one face */
        uint8_t *has_face = (uint8_t *)ARENA_CALLOC(arena, (long)cm->nv, (long)sizeof(uint8_t));
        for (size_t f = 0; f < cm->nf; f++)
            has_face[uf_find(&uf, cm->faces[f*3])] = 1;
        for (size_t r = 0; r < cm->nv; r++) {
            if (!has_face[r]) continue;
            if (cnt >= cap) {
                size_t nc = cap*2;
                ComponentMesh *na = (ComponentMesh *)ARENA_ALLOC(arena, (long)(nc*sizeof(ComponentMesh)));
                memcpy(na, ccs, cnt*sizeof(ComponentMesh));
                for (size_t _i=0;_i<cnt;_i++) na[_i].self = &na[_i];
                double *naa = (double *)ARENA_ALLOC(arena, (long)(nc*sizeof(double)));
                memcpy(naa, areas, cnt*sizeof(double));
                size_t *npa = (size_t *)ARENA_ALLOC(arena, (long)(nc*sizeof(size_t)));
                memcpy(npa, parent, cnt*sizeof(size_t));
                ccs = na; areas = naa; parent = npa; cap = nc;
            }
            double a = 0.0;
            if (extract_cc(arena, cm, &uf, (int32_t)r, &ccs[cnt], &a) == 0) {
                areas[cnt] = a; parent[cnt] = i; parent_total[i] += a; cnt++;
            }
        }
    }
    if (cnt == 0) return 0;

    /* Cull relative to the component's UPSTREAM input sheet, not the whole cube.
     * A wrapped cube can legitimately contain dozens of similarly-sized sheets;
     * comparing each with the cube total deletes all of them once the sheet count
     * exceeds 1/min_frac. The input array already carries the semantic split, so
     * this pass should only remove disconnected crumbs within each input. */
    uint8_t *take = (uint8_t *)ARENA_CALLOC(arena, (long)cnt,
                                             (long)sizeof(uint8_t));
    uint8_t *parent_has = (uint8_t *)ARENA_CALLOC(arena, (long)n_in,
                                                   (long)sizeof(uint8_t));
    size_t *best = (size_t *)ARENA_ALLOC(arena, (long)(n_in*sizeof(size_t)));
    double *best_area = (double *)ARENA_CALLOC(arena, (long)n_in,
                                                (long)sizeof(double));
    for (size_t i = 0; i < n_in; i++) best[i] = (size_t)-1;
    for (size_t c = 0; c < cnt; c++) {
        size_t p = parent[c];
        double thresh = (double)min_frac * parent_total[p];
        if (areas[c] >= thresh) { take[c] = 1; parent_has[p] = 1; }
        if (best[p] == (size_t)-1 || areas[c] > best_area[p]) {
            best[p] = c; best_area[p] = areas[c];
        }
    }
    /* Never erase an upstream sheet entirely, even if it fragmented into more
     * than 1/min_frac equal pieces. Retain that parent's largest component. */
    for (size_t i = 0; i < n_in; i++)
        if (!parent_has[i] && best[i] != (size_t)-1) take[best[i]] = 1;

    ComponentMesh *keep = (ComponentMesh *)ARENA_ALLOC(arena, (long)(cnt*sizeof(ComponentMesh)));
    size_t nk = 0, dropped = 0;
    for (size_t c = 0; c < cnt; c++) {
        if (take[c]) {
            keep[nk] = ccs[c];
            keep[nk].self = &keep[nk];
            keep[nk].comp_id = (int)(nk + 1);
            nk++;
        } else {
            dropped++;
        }
    }

    *out = keep; *n_out = nk;
    if (n_culled) *n_culled = dropped;
    return 0;
}
