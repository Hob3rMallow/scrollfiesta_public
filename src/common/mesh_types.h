#ifndef MESH_TYPES_INCLUDED
#define MESH_TYPES_INCLUDED

#include <stdint.h>
#include <stddef.h>

typedef struct ComponentMesh {
    float   *verts;          /* [nv * 3], float32, (z,y,x) order */
    int32_t *faces;          /* [nf * 3], int32, 0-based indices */
    uint8_t *pin_mask;       /* [nv], 1 = halo vert (must not move), 0 = owned.
                              * NULL if extracted without halo (halo=0). */
    float   *vert_normals;   /* [nv * 3], float32, (z,y,x) order. Per-vertex
                              * unit normal from MLS-midpoint projection in
                              * Step 0. NULL if the mesh predates MLS or if
                              * vert_normals were not computed. Downstream
                              * code that needs a normal at a point should
                              * prefer this over the mesh-level `pca_normal`. */
    size_t   nv;             /* vertex count */
    size_t   nf;             /* face (triangle) count */
    int      comp_id;        /* 1-based, assigned by Step 4 */
    float    pca_normal[3];  /* unit normal — vertex-area-weighted mean of
                              * vert_normals when MLS is used; PCA
                              * eigenvector on the raw envelope when not. */
    float    centroid[3];    /* mean vertex position from Step 0 */
    size_t   nv_pre_fill;    /* verts before hole fill (0 = none) */
    void    *self;           /* validation sentinel: self == &this */
} ComponentMesh;

static inline int ComponentMesh_valid(const ComponentMesh *cm)
{
    return cm && cm->self == (void *)cm && cm->nv > 0 && cm->verts
           && cm->faces && cm->nf > 0;
}

#endif
