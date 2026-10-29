#ifndef RIBBON_DOMAINS_INCLUDED
#define RIBBON_DOMAINS_INCLUDED

#include "../common/mesh_bin.h"

/* Work batches contain WHOLE connected components of the original indexed
 * mesh. Coincident geometry, UVs and material labels never join components.
 * Every source vertex (including isolated points) and face has exactly one
 * owner. Component order uses its smallest original vertex, independent of
 * union-find root selection; indices inside a batch preserve source order.
 * This partitions work, not topology, and does not imply physical clearance. */
typedef struct RibbonDomains_T *RibbonDomains_T;
typedef struct {
    size_t vertices, faces, components;
    size_t input_bytes; /* XYZ/UV/faces, including arena alignment; no solver storage */
    const int32_t *source_vertices, *source_faces;
    int oversized; /* a single whole component exceeds the target batch size */
} RibbonDomainPart;
typedef struct {
    size_t components, parts, oversized_parts;
    size_t largest_component_vertices, largest_component_faces;
    size_t largest_part_vertices, largest_part_faces;
} RibbonDomainsReport;

/* Targets bound ordinary batches; an oversized connected component remains
 * whole and is reported explicitly, never split or dropped. All index storage
 * is arena-owned. Source arrays are read-only and must outlive the plan. */
int RibbonDomains_new(Arena_T arena, const MeshBinData *source,
    size_t target_faces, size_t target_vertices, RibbonDomains_T *out);
int RibbonDomains_report(RibbonDomains_T domains, RibbonDomainsReport *out);
int RibbonDomains_part(RibbonDomains_T domains, size_t part, RibbonDomainPart *out);
/* Dense read-only owner arrays are indexed by original vertex/face. */
const int32_t *RibbonDomains_vertex_owners(RibbonDomains_T domains);
const int32_t *RibbonDomains_face_owners(RibbonDomains_T domains);
/* Extract just one batch in caller scratch storage. Enforces the actual
 * XYZ/UV/face byte limit BEFORE allocating; a large component is not permission
 * to exceed it. Original IDs are supplied by part(), not inferred spatially. */
int RibbonDomains_extract(Arena_T scratch, RibbonDomains_T domains, size_t part,
    size_t maximum_bytes, MeshBinData *out);
int RibbonDomains_selftest(void);

#endif
