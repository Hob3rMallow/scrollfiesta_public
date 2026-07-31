/*
 * gco_wrap.cpp -- C++ wrapper calling GCO v3.0 alpha-expansion.
 *
 * Compiles as C++17.  Linked into the C pipeline via extern "C".
 */
#include "gco_wrap.h"
#include "GCoptimization.h"

#include <cstdio>

struct GCO_Opaque {
    GCoptimizationGeneralGraph *gc;
};

extern "C" {

GCO_Handle GCO_create(int num_sites, int num_labels)
{
    if (num_sites < 1 || num_labels < 2) return nullptr;
    try {
        auto *h = new GCO_Opaque;
        h->gc = new GCoptimizationGeneralGraph(num_sites, num_labels);
        h->gc->setVerbosity(0);
        return h;
    } catch (GCException &e) {
        fprintf(stderr, "GCO_create: %s\n", e.message);
        return nullptr;
    }
}

void GCO_set_data_cost(GCO_Handle gc, const int32_t *data)
{
    if (!gc || !data) return;
    try {
        gc->gc->setDataCost(const_cast<int32_t *>(data));
    } catch (GCException &e) {
        fprintf(stderr, "GCO_set_data_cost: %s\n", e.message);
    }
}

void GCO_set_smooth_cost(GCO_Handle gc, const int32_t *smooth)
{
    if (!gc || !smooth) return;
    try {
        gc->gc->setSmoothCost(const_cast<int32_t *>(smooth));
    } catch (GCException &e) {
        fprintf(stderr, "GCO_set_smooth_cost: %s\n", e.message);
    }
}

void GCO_set_neighbor(GCO_Handle gc, int s1, int s2, int weight)
{
    if (!gc) return;
    try {
        gc->gc->setNeighbors(s1, s2, weight);
    } catch (GCException &e) {
        fprintf(stderr, "GCO_set_neighbor: %s\n", e.message);
    }
}

long long GCO_expansion(GCO_Handle gc, int max_iter)
{
    if (!gc) return -1;
    try {
        return gc->gc->expansion(max_iter);
    } catch (GCException &e) {
        fprintf(stderr, "GCO_expansion: %s\n", e.message);
        return -1;
    }
}

int GCO_what_label(GCO_Handle gc, int site)
{
    if (!gc) return -1;
    return gc->gc->whatLabel(site);
}

void GCO_get_labels(GCO_Handle gc, int *out_labels, int num_sites)
{
    if (!gc || !out_labels) return;
    gc->gc->whatLabel(0, num_sites, out_labels);
}

void GCO_set_label(GCO_Handle gc, int site, int label)
{
    if (!gc) return;
    gc->gc->setLabel(site, label);
}

void GCO_destroy(GCO_Handle gc)
{
    if (!gc) return;
    delete gc->gc;
    delete gc;
}

} /* extern "C" */
