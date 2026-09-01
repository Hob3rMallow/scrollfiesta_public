/*
 * gco_wrap.h -- C-callable wrapper for GCO v3.0 graph-cut optimization.
 *
 * Wraps GCoptimizationGeneralGraph for multi-label MRF optimization.
 * Uses extern "C" linkage so the C pipeline can call it directly.
 *
 * All internal memory is managed by GCO (new/delete).
 * Caller is responsible for the lifetime of arrays passed in.
 */
#ifndef GCO_WRAP_INCLUDED
#define GCO_WRAP_INCLUDED

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GCO_Opaque *GCO_Handle;

/* Site-aware pairwise cost.  This is required for winding constraints whose
 * preferred label difference depends on the observed edge. */
typedef int32_t (*GCO_SmoothCostCallback)(
    int site1, int site2, int label1, int label2, void *userdata);

/* Create a general graph with num_sites nodes and num_labels labels. */
GCO_Handle GCO_create(int num_sites, int num_labels);

/* Set unary (data) costs: array[site * num_labels + label] = cost (int32). */
void GCO_set_data_cost(GCO_Handle gc, const int32_t *data);

/* Set pairwise (smooth) costs: array[l1 * num_labels + l2] (Potts matrix). */
void GCO_set_smooth_cost(GCO_Handle gc, const int32_t *smooth);

/* Set a site-aware pairwise cost callback.  userdata must remain alive until
 * optimization has finished. */
void GCO_set_smooth_cost_callback(
    GCO_Handle gc, GCO_SmoothCostCallback callback, void *userdata);

/* Add an undirected edge between s1 and s2 with given weight. */
void GCO_set_neighbor(GCO_Handle gc, int s1, int s2, int weight);

/* Run alpha-expansion for max_iter iterations. Returns final energy. */
long long GCO_expansion(GCO_Handle gc, int max_iter);

/* Run alpha-beta swap for max_iter iterations.  Shifted label-difference
 * costs used by winding registration are swap-regular but are not, in
 * general, alpha-expansion regular. */
long long GCO_swap(GCO_Handle gc, int max_iter);

/* Retrieve the label assigned to a single site. */
int GCO_what_label(GCO_Handle gc, int site);

/* Bulk retrieve labels for all sites into out_labels[num_sites]. */
void GCO_get_labels(GCO_Handle gc, int *out_labels, int num_sites);

/* Set initial label for a site (call before expansion). */
void GCO_set_label(GCO_Handle gc, int site, int label);

/* Destroy the graph and free all internal memory. */
void GCO_destroy(GCO_Handle gc);

#ifdef __cplusplus
}
#endif

#endif /* GCO_WRAP_INCLUDED */
