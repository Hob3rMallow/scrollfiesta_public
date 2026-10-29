#include "winding_mrf.h"

#include "../common/gco_wrap.h"
#include "../common/union_find.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WM_COST_MAX (INT32_MAX / 32)

typedef struct {
    int32_t a, b;
    int32_t target;
    int kind;
    double weight;
    int32_t graph_weight;
    int64_t graph_baseline;
} WmFactor;

typedef struct {
    int32_t a, b;
    size_t first, count;
} WmGroup;

typedef struct {
    const WmFactor *factors;
    const WmGroup *groups;
    size_t ngroups;
} WmSmoothContext;

static int wm_factor_compare(const void *lhs, const void *rhs)
{
    const WmFactor *a = (const WmFactor *)lhs;
    const WmFactor *b = (const WmFactor *)rhs;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    if (a->b != b->b) return a->b < b->b ? -1 : 1;
    if (a->target != b->target) return a->target < b->target ? -1 : 1;
    if (a->kind != b->kind) return a->kind < b->kind ? -1 : 1;
    return 0;
}

static int64_t wm_violation(const WmFactor *factor, int64_t difference)
{
    int64_t residual=difference-factor->target;
    if (factor->kind==WINDING_MRF_AT_LEAST) return residual<0 ? -residual : 0;
    if (factor->kind==WINDING_MRF_AT_MOST) return residual>0 ? residual : 0;
    return residual<0 ? -residual : residual;
}

static const WmGroup *wm_find_group(
    const WmGroup *groups, size_t ngroups, int32_t a, int32_t b)
{
    size_t lo = 0, hi = ngroups;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        const WmGroup *group = &groups[mid];
        if (group->a < a || (group->a == a && group->b < b)) lo = mid + 1;
        else hi = mid;
    }
    if (lo < ngroups && groups[lo].a == a && groups[lo].b == b)
        return &groups[lo];
    return NULL;
}

static int32_t wm_round_cost(double cost)
{
    if (!(cost > 0.0)) return 0;
    if (!isfinite(cost) || cost >= (double)WM_COST_MAX) return WM_COST_MAX;
    return (int32_t)floor(cost + 0.5);
}

static double wm_huber(double z, double delta)
{
    z = fabs(z);
    if (z <= delta) return 0.5 * z * z;
    return delta * (z - 0.5 * delta);
}

static double wm_unary_energy(
    const WindingMRFSite *site, int32_t label,
    const WindingMRFOptions *options)
{
    if (site->fixed)
        return label == site->initial_label ? 0.0 : (double)WM_COST_MAX;
    if (!(site->weight > 0.0) || !(site->sigma > 0.0) ||
        !isfinite(site->center) || !isfinite(site->sigma) ||
        !isfinite(site->weight))
        return 0.0;
    double z = ((double)label - site->center) / site->sigma;
    return site->weight * wm_huber(z, options->unary_huber_delta);
}

static int32_t wm_smooth_cost(
    int site1, int site2, int label1, int label2, void *userdata)
{
    WmSmoothContext *context = (WmSmoothContext *)userdata;
    int32_t a = site1 < site2 ? site1 : site2;
    int32_t b = site1 < site2 ? site2 : site1;
    const WmGroup *group = wm_find_group(
        context->groups, context->ngroups, a, b);
    if (group == NULL) return 0;
    int difference = site1 < site2 ? label2 - label1 : label1 - label2;
    int64_t energy = 0;
    for (size_t i = 0; i < group->count; i++) {
        const WmFactor *factor = &context->factors[group->first + i];
        int64_t residual = wm_violation(factor,difference)-factor->graph_baseline;
        energy += (int64_t)factor->graph_weight * residual;
        if (energy >= WM_COST_MAX) return WM_COST_MAX;
    }
    return (int32_t)energy;
}

static int32_t wm_initial_label(
    const WindingMRFSite *site, const WindingMRFOptions *options)
{
    int64_t label;
    if (site->initial_label != INT32_MIN) label = site->initial_label;
    else if (isfinite(site->center)) label = (int64_t)llround(site->center);
    else label = ((int64_t)options->label_min + options->label_max) / 2;
    if (label < options->label_min) label = options->label_min;
    if (label > options->label_max) label = options->label_max;
    return (int32_t)label;
}

static double wm_pair_energy(
    const WmFactor *factors, size_t nfactors, const int32_t *labels)
{
    double energy = 0.0;
    for (size_t i = 0; i < nfactors; i++) {
        const WmFactor *factor = &factors[i];
        int64_t difference = (int64_t)labels[factor->b] - labels[factor->a];
        energy += factor->weight * (double)wm_violation(factor,difference);
    }
    return energy;
}

static double wm_natural_energy(
    const WindingMRFSite *sites, size_t nsites,
    const WmFactor *factors, size_t nfactors,
    const WindingMRFOptions *options, const int32_t *labels,
    const double *extra_unary, int nlabels)
{
    double energy = wm_pair_energy(factors, nfactors, labels);
    for (size_t i = 0; i < nsites; i++) {
        energy += wm_unary_energy(&sites[i], labels[i], options);
        if (extra_unary != NULL) {
            int k = labels[i] - options->label_min;
            if (k >= 0 && k < nlabels)
                energy += extra_unary[i * (size_t)nlabels + (size_t)k];
        }
    }
    return energy;
}

static int64_t wm_quantized_energy(
    const WindingMRFSite *sites, size_t nsites,
    const WmGroup *groups, size_t ngroups,
    const WmSmoothContext *context, const WindingMRFOptions *options,
    const int32_t *labels, const double *extra_unary, int nlabels)
{
    int64_t energy = 0;
    for (size_t i = 0; i < nsites; i++) {
        double unary = wm_unary_energy(&sites[i], labels[i], options);
        if (extra_unary != NULL) {
            int k = labels[i] - options->label_min;
            if (k >= 0 && k < nlabels)
                unary += extra_unary[i * (size_t)nlabels + (size_t)k];
        }
        energy += wm_round_cost(options->cost_scale * unary);
    }
    for (size_t i = 0; i < ngroups; i++) {
        const WmGroup *group = &groups[i];
        int l1 = labels[group->a] - options->label_min;
        int l2 = labels[group->b] - options->label_min;
        energy += wm_smooth_cost(group->a, group->b, l1, l2,
                                 (void *)context);
    }
    return energy;
}

typedef struct {
    double natural, compensation, magnitude;
    int64_t quantized;
    size_t sites;
    int reject;
} WmMoveDelta;

static void wm_move_add(WmMoveDelta *move,double value)
{
    double corrected=value-move->compensation;
    double sum=move->natural+corrected;
    move->compensation=(sum-move->natural)-corrected;
    move->natural=sum;
    move->magnitude+=fabs(value);
}

/* Changed-site components in the ORIGINAL factor graph have independent
 * energy deltas: any edge between two changed sites puts them in the same
 * component; all remaining endpoints are unchanged boundary conditions.
 * Accepting one whole component cannot change another's energy. This keeps
 * improving moves while refusing integer-cost ties which worsen weak real
 * evidence, rather than letting four leaves drift or veto a whole scroll. */
static void wm_accept_components(
    Arena_T arena,const WindingMRFSite *sites,size_t nsites,
    const WmSmoothContext *context,size_t nfactors,const WindingMRFOptions *options,
    const double *extra_unary,int nlabels,const int32_t *initial,int32_t *labels,
    size_t *rejected_components,size_t *rejected_sites)
{
    Arena_Mark mark=Arena_save(arena);
    UnionFind graph=UF_new(arena,(int32_t)nsites);
    WmMoveDelta *delta=ARENA_CALLOC(arena,nsites,sizeof *delta);
    *rejected_components=0; *rejected_sites=0;
    for (size_t i=0;i<nfactors;i++) {
        const WmFactor *f=context->factors+i;
        if (labels[f->a]!=initial[f->a] && labels[f->b]!=initial[f->b])
            uf_union(&graph,f->a,f->b);
    }
    for (size_t i=0;i<nsites;i++) if (labels[i]!=initial[i]) {
        WmMoveDelta *move=delta+uf_find(&graph,(int32_t)i);
        double before=wm_unary_energy(sites+i,initial[i],options);
        double after=wm_unary_energy(sites+i,labels[i],options);
        if (extra_unary) {
            before+=extra_unary[i*(size_t)nlabels+(size_t)(initial[i]-options->label_min)];
            after+=extra_unary[i*(size_t)nlabels+(size_t)(labels[i]-options->label_min)];
        }
        wm_move_add(move,after-before);
        move->quantized+=(int64_t)wm_round_cost(options->cost_scale*after)-
                                  wm_round_cost(options->cost_scale*before);
        move->sites++;
        if (sites[i].fixed) move->reject=1;
    }
    for (size_t i=0;i<nfactors;i++) {
        const WmFactor *f=context->factors+i;
        int32_t changed=labels[f->a]!=initial[f->a] ? f->a :
                        labels[f->b]!=initial[f->b] ? f->b : -1;
        if (changed<0) continue;
        int64_t before=(int64_t)initial[f->b]-initial[f->a];
        int64_t after=(int64_t)labels[f->b]-labels[f->a];
        /* Subtract integer violations BEFORE multiplying by the real weight;
         * a far target's unavoidable constant cannot erase a small delta. */
        wm_move_add(delta+uf_find(&graph,changed),
                    f->weight*(double)(wm_violation(f,after)-wm_violation(f,before)));
    }
    for (size_t i=0;i<context->ngroups;i++) {
        const WmGroup *g=context->groups+i;
        int32_t changed=labels[g->a]!=initial[g->a] ? g->a :
                        labels[g->b]!=initial[g->b] ? g->b : -1;
        if (changed<0) continue;
        int32_t before=wm_smooth_cost(g->a,g->b,initial[g->a]-options->label_min,
                                      initial[g->b]-options->label_min,(void *)context);
        int32_t after=wm_smooth_cost(g->a,g->b,labels[g->a]-options->label_min,
                                     labels[g->b]-options->label_min,(void *)context);
        delta[uf_find(&graph,changed)].quantized+=(int64_t)after-before;
    }
    for (size_t i=0;i<nsites;i++) if (delta[i].sites) {
        WmMoveDelta *move=delta+i;
        if (!isfinite(move->natural) || !isfinite(move->magnitude) ||
            !(move->natural < -64*DBL_EPSILON*move->magnitude) || move->quantized>0)
            move->reject=1;
        if (move->reject) {
            (*rejected_components)++; *rejected_sites+=move->sites;
        }
    }
    for (size_t i=0;i<nsites;i++)
        if (labels[i]!=initial[i] && delta[uf_find(&graph,(int32_t)i)].reject)
            labels[i]=initial[i];
    Arena_restore(arena,mark);
}

static double wm_site_conditional_energy(
    size_t site_index, int32_t candidate,
    const WindingMRFSite *sites,
    const WmFactor *factors,
    const size_t *incident_offset, const int32_t *incident_factor,
    const WindingMRFOptions *options, const int32_t *labels,
    const double *extra_unary, int nlabels)
{
    double energy = wm_unary_energy(&sites[site_index], candidate, options);
    if (extra_unary != NULL) {
        int k = candidate - options->label_min;
        if (k >= 0 && k < nlabels)
            energy += extra_unary[site_index * (size_t)nlabels + (size_t)k];
    }
    for (size_t p = incident_offset[site_index];
         p < incident_offset[site_index + 1]; p++) {
        const WmFactor *factor = &factors[incident_factor[p]];
        if ((size_t)factor->a == site_index) {
            int64_t difference = (int64_t)labels[factor->b] - candidate;
            energy += factor->weight *
                      (double)wm_violation(factor,difference);
        } else if ((size_t)factor->b == site_index) {
            int64_t difference = (int64_t)candidate - labels[factor->a];
            energy += factor->weight *
                      (double)wm_violation(factor,difference);
        }
    }
    return energy;
}

/* The graph cut cannot locate a site whose every pairwise weight rounds to
 * zero. Solve that site's finite conditional objective with the ORIGINAL
 * weights instead. Neighbor labels are held at the current proposal, and
 * candidates may not increase the already quantized unary. Thus this exact
 * scalar solve cannot raise the graph-cut objective or inflate weak evidence.
 * A sweep is not a global optimum for mutually coupled invisible sites. */
static size_t wm_refine_invisible_sites(
    const WindingMRFSite *sites,size_t nsites,const WmFactor *factors,
    const size_t *incident_offset,const int32_t *incident_factor,
    const WindingMRFOptions *options,const double *extra_unary,int nlabels,
    const int32_t *data,const int32_t *initial,int32_t *labels)
{
    size_t refined=0;
    for (size_t i=0;i<nsites;i++) {
        int visible=sites[i].fixed;
        for (size_t p=incident_offset[i];!visible && p<incident_offset[i+1];p++)
            visible=factors[incident_factor[p]].graph_weight!=0;
        if (visible) continue;
        int32_t original=labels[i];
        int32_t budget=data[i*(size_t)nlabels+(size_t)(original-options->label_min)];
        int32_t best=data[i*(size_t)nlabels+(size_t)(initial[i]-options->label_min)]<=budget
                     ? initial[i] : original;
        double energy=wm_site_conditional_energy(i,best,sites,factors,incident_offset,
                                                 incident_factor,options,labels,extra_unary,nlabels);
        for (int k=0;k<nlabels;k++) {
            int32_t candidate=options->label_min+k;
            if (data[i*(size_t)nlabels+(size_t)k]>budget) continue;
            double value=wm_site_conditional_energy(i,candidate,sites,factors,incident_offset,
                                                    incident_factor,options,labels,extra_unary,nlabels);
            if (value<energy) { energy=value; best=candidate; }
        }
        labels[i]=best;
        refined+=best!=original;
    }
    return refined;
}

void WindingMRF_default_options(WindingMRFOptions *options)
{
    if (options == NULL) return;
    options->label_min = -32;
    options->label_max = 32;
    options->max_iterations = 8;
    options->cost_scale = 256.0;
    options->unary_huber_delta = 2.5;
    options->confidence_temperature = 1.0;
    options->abstain_threshold = 0.35;
}

static int wm_solve(
    Arena_T arena,
    const WindingMRFSite *sites, size_t nsites,
    const WindingMRFEdge *edges, size_t nedges,
    const int *kinds,
    const WindingMRFUnaryPenalty *penalties, size_t npenalties,
    const WindingMRFOptions *options_arg,
    int32_t **out_labels, float **out_confidence,
    WindingMRFStats *stats)
{
    if (out_labels!=NULL) *out_labels=NULL;
    if (out_confidence!=NULL) *out_confidence=NULL;
    WindingMRFOptions defaults;
    WindingMRF_default_options(&defaults);
    const WindingMRFOptions *options = options_arg != NULL
                                     ? options_arg : &defaults;
    if (arena == NULL || sites == NULL || nsites == 0 || out_labels == NULL ||
        nsites > (size_t)INT_MAX || nedges > (size_t)INT_MAX ||
        options->label_max <= options->label_min ||
        options->max_iterations < 1 || !(options->cost_scale > 0.0) ||
        !isfinite(options->cost_scale) ||
        !(options->unary_huber_delta > 0.0) ||
        !isfinite(options->unary_huber_delta) ||
        !(options->confidence_temperature > 0.0) ||
        !isfinite(options->confidence_temperature) ||
        !(options->abstain_threshold >= 0.0) ||
        !(options->abstain_threshold <= 1.0) ||
        (nedges > 0 && edges == NULL) ||
        (npenalties > 0 && penalties == NULL))
        return -1;
    int64_t label_count64 = (int64_t)options->label_max -
                            options->label_min + 1;
    if (label_count64 < 2 || label_count64 > INT_MAX ||
        nsites > SIZE_MAX / (size_t)label_count64)
        return -1;
    int nlabels = (int)label_count64;
    for (size_t i = 0; i < nsites; i++) {
        const WindingMRFSite *site = &sites[i];
        if (site->weight < 0.0 || !isfinite(site->weight)) return -1;
        if (site->fixed && (site->initial_label < options->label_min ||
                            site->initial_label > options->label_max))
            return -1;
    }
    for (size_t i = 0; i < nedges; i++) {
        const WindingMRFEdge *edge = &edges[i];
        if (edge->a < 0 || edge->b < 0 || edge->a == edge->b ||
            (size_t)edge->a >= nsites || (size_t)edge->b >= nsites ||
            edge->target==INT32_MIN ||
            (kinds && (kinds[i]<WINDING_MRF_EQUAL || kinds[i]>WINDING_MRF_AT_MOST)) ||
            edge->weight < 0.0 || !isfinite(edge->weight))
            return -1;
    }
    for (size_t i = 0; i < npenalties; i++) {
        const WindingMRFUnaryPenalty *penalty = &penalties[i];
        if (penalty->site < 0 || (size_t)penalty->site >= nsites ||
            penalty->label < options->label_min ||
            penalty->label > options->label_max ||
            penalty->weight < 0.0 || !isfinite(penalty->weight))
            return -1;
    }

    int32_t *labels = ARENA_ALLOC(arena, nsites * sizeof *labels);
    float *confidence = out_confidence != NULL
                      ? ARENA_ALLOC(arena, nsites * sizeof *confidence) : NULL;
    *out_labels = labels;
    if (out_confidence != NULL) *out_confidence = confidence;
    Arena_Mark scratch = Arena_save(arena);

    WmFactor *factors = nedges > 0
                      ? ARENA_ALLOC(arena, nedges * sizeof *factors) : NULL;
    size_t nfactors = 0;
    for (size_t i = 0; i < nedges; i++) {
        const WindingMRFEdge *edge = &edges[i];
        if (!(edge->weight > 0.0)) continue;
        WmFactor factor;
        factor.kind=kinds ? kinds[i] : WINDING_MRF_EQUAL;
        if (edge->a < edge->b) {
            factor.a = edge->a;
            factor.b = edge->b;
            factor.target = edge->target;
        } else {
            factor.a = edge->b;
            factor.b = edge->a;
            factor.target = -edge->target;
            if (factor.kind==WINDING_MRF_AT_LEAST) factor.kind=WINDING_MRF_AT_MOST;
            else if (factor.kind==WINDING_MRF_AT_MOST) factor.kind=WINDING_MRF_AT_LEAST;
        }
        factor.weight = edge->weight;
        factor.graph_weight = wm_round_cost(
            options->cost_scale * edge->weight);
        /* Remove only the unavoidable constant within this label window.
         * A distant, weak outlier otherwise overflows GCO even though its
         * preference between allowed labels is small. Natural energy and
         * residual diagnostics retain the complete, unshifted observation. */
        int64_t width=label_count64-1, best=factor.target;
        if (best < -width) best=-width;
        if (best > width) best=width;
        factor.graph_baseline=wm_violation(&factor,best);
        /* Do not promote arbitrarily weak evidence to one graph-cost unit:
         * one source face split over many factors would regain unbounded
         * influence through that minimum. Zero-quantized terms remain in
         * natural-energy/confidence accounting, but do not force a cut. */
        factors[nfactors++] = factor;
    }
    if (nfactors > 1)
        qsort(factors, nfactors, sizeof *factors, wm_factor_compare);
    WmGroup *groups = nfactors > 0
                    ? ARENA_ALLOC(arena, nfactors * sizeof *groups) : NULL;
    size_t ngroups = 0;
    for (size_t i = 0; i < nfactors;) {
        size_t end = i + 1;
        while (end < nfactors && factors[end].a == factors[i].a &&
               factors[end].b == factors[i].b)
            end++;
        groups[ngroups].a = factors[i].a;
        groups[ngroups].b = factors[i].b;
        groups[ngroups].first = i;
        groups[ngroups].count = end - i;
        ngroups++;
        i = end;
    }

    /* Confidence is a site-conditional posterior.  The old implementation
     * found a site's incident factors by scanning the entire factor array for
     * every candidate label, twice.  That is O(S*L*E) and made a fragmented
     * 10x10x10 (97k sites) spend minutes per conflict round on a nominally
     * diagnostic statistic.  Build a compact CSR incidence index once so all
     * conditional-energy work is O(L*(S+E)). */
    size_t *incident_offset = ARENA_CALLOC(
        arena, nsites + 1, sizeof *incident_offset);
    for (size_t i = 0; i < nfactors; i++) {
        incident_offset[(size_t)factors[i].a + 1]++;
        incident_offset[(size_t)factors[i].b + 1]++;
    }
    for (size_t i = 1; i <= nsites; i++) {
        if (incident_offset[i] > SIZE_MAX - incident_offset[i - 1]) {
            Arena_restore(arena, scratch);
            return -1;
        }
        incident_offset[i] += incident_offset[i - 1];
    }
    if (nfactors > SIZE_MAX / (2 * sizeof(int32_t))) {
        Arena_restore(arena, scratch);
        return -1;
    }
    int32_t *incident_factor = nfactors > 0
        ? ARENA_ALLOC(arena, 2 * nfactors * sizeof *incident_factor) : NULL;
    size_t *incident_cursor = ARENA_ALLOC(
        arena, nsites * sizeof *incident_cursor);
    memcpy(incident_cursor, incident_offset,
           nsites * sizeof *incident_cursor);
    for (size_t i = 0; i < nfactors; i++) {
        incident_factor[incident_cursor[factors[i].a]++] = (int32_t)i;
        incident_factor[incident_cursor[factors[i].b]++] = (int32_t)i;
    }

    size_t data_count = nsites * (size_t)nlabels;
    int32_t *data = ARENA_ALLOC(arena, data_count * sizeof *data);
    double *extra_unary = npenalties > 0
        ? ARENA_CALLOC(arena, data_count, sizeof *extra_unary) : NULL;
    if (extra_unary != NULL) {
        double cap = (double)WM_COST_MAX / options->cost_scale;
        for (size_t i = 0; i < npenalties; i++) {
            const WindingMRFUnaryPenalty *penalty = &penalties[i];
            size_t k = (size_t)(penalty->label - options->label_min);
            size_t index = (size_t)penalty->site * (size_t)nlabels + k;
            double value = extra_unary[index] + penalty->weight;
            extra_unary[index] = value < cap ? value : cap;
        }
    }
    int *gco_labels = ARENA_ALLOC(arena, nsites * sizeof *gco_labels);
    for (size_t i = 0; i < nsites; i++) {
        labels[i] = wm_initial_label(&sites[i], options);
        for (int k = 0; k < nlabels; k++) {
            int32_t label = options->label_min + k;
            double unary = wm_unary_energy(&sites[i], label, options);
            if (extra_unary != NULL)
                unary += extra_unary[i * (size_t)nlabels + (size_t)k];
            data[i * (size_t)nlabels + (size_t)k] = wm_round_cost(
                options->cost_scale * unary);
        }
    }

    WmSmoothContext context;
    context.factors = factors;
    context.groups = groups;
    context.ngroups = ngroups;
    int64_t quantized_before = wm_quantized_energy(
        sites, nsites, groups, ngroups, &context, options, labels,
        extra_unary, nlabels);
    double energy_before = wm_natural_energy(
        sites, nsites, factors, nfactors, options, labels,
        extra_unary, nlabels);

    GCO_Handle graph = GCO_create((int)nsites, nlabels);
    if (graph == NULL) {
        Arena_restore(arena, scratch);
        return -1;
    }
    GCO_set_data_cost(graph, data);
    GCO_set_smooth_cost_callback(graph, wm_smooth_cost, &context);
    for (size_t i = 0; i < ngroups; i++)
        GCO_set_neighbor(graph, groups[i].a, groups[i].b, 1);
    for (size_t i = 0; i < nsites; i++)
        GCO_set_label(graph, (int)i, labels[i] - options->label_min);
    long long optimizer_energy = GCO_swap(graph, options->max_iterations);
    if (optimizer_energy < 0) {
        GCO_destroy(graph);
        Arena_restore(arena, scratch);
        return -1;
    }
    GCO_get_labels(graph, gco_labels, (int)nsites);
    GCO_destroy(graph);

    size_t changed = 0, proposed = 0, rejected_components = 0, rejected_sites = 0;
    int32_t *initial = ARENA_ALLOC(arena, nsites * sizeof *initial);
    memcpy(initial, labels, nsites * sizeof *initial);
    for (size_t i = 0; i < nsites; i++) {
        labels[i] = options->label_min + gco_labels[i];
        proposed += labels[i] != initial[i];
    }
    size_t natural_refined = wm_refine_invisible_sites(
        sites, nsites, factors, incident_offset, incident_factor, options,
        extra_unary, nlabels, data, initial, labels);
    wm_accept_components(arena, sites, nsites, &context, nfactors, options,
                          extra_unary, nlabels, initial, labels,
                          &rejected_components, &rejected_sites);
    int32_t solution_min = options->label_max;
    int32_t solution_max = options->label_min;
    for (size_t i = 0; i < nsites; i++) {
        int32_t solved = labels[i];
        if (solved != initial[i]) changed++;
        if (solved < solution_min) solution_min = solved;
        if (solved > solution_max) solution_max = solved;
    }
    int64_t quantized_after = wm_quantized_energy(
        sites, nsites, groups, ngroups, &context, options, labels,
        extra_unary, nlabels);
    double energy_after = wm_natural_energy(
        sites, nsites, factors, nfactors, options, labels,
        extra_unary, nlabels);

    size_t abstained = 0;
    double confidence_sum = 0.0;
    double confidence_min = DBL_MAX;
    for (size_t i = 0; i < nsites; i++) {
        double minimum = DBL_MAX;
        for (int k = 0; k < nlabels; k++) {
            int32_t label = options->label_min + k;
            double value = wm_site_conditional_energy(
                i, label, sites, factors, incident_offset, incident_factor,
                options, labels,
                extra_unary, nlabels);
            if (value < minimum) minimum = value;
        }
        double sum = 0.0, maximum_probability_numerator = 0.0;
        for (int k = 0; k < nlabels; k++) {
            int32_t label = options->label_min + k;
            double value = wm_site_conditional_energy(
                i, label, sites, factors, incident_offset, incident_factor,
                options, labels,
                extra_unary, nlabels);
            double probability_numerator = exp(
                -(value - minimum) / options->confidence_temperature);
            sum += probability_numerator;
            if (label == labels[i])
                maximum_probability_numerator = probability_numerator;
        }
        double probability = sum > 0.0
                           ? maximum_probability_numerator / sum : 0.0;
        double flat = 1.0 / nlabels;
        double sharpness = (probability - flat) / (1.0 - flat);
        if (sharpness < 0.0) sharpness = 0.0;
        if (sharpness > 1.0) sharpness = 1.0;
        if (confidence != NULL) confidence[i] = (float)sharpness;
        confidence_sum += sharpness;
        if (sharpness < confidence_min) confidence_min = sharpness;
        if (sharpness < options->abstain_threshold) abstained++;
    }

    if (stats != NULL) {
        memset(stats, 0, sizeof *stats);
        stats->sites = nsites;
        stats->input_edges = nedges;
        stats->active_edges = nfactors;
        stats->unique_pairs = ngroups;
        stats->labels = (size_t)nlabels;
        stats->changed_labels = changed;
        stats->abstained_sites = abstained;
        stats->solution_min = solution_min;
        stats->solution_max = solution_max;
        stats->quantized_energy_before = quantized_before;
        stats->quantized_energy_after = quantized_after;
        stats->energy_before = energy_before;
        stats->energy_after = energy_after;
        stats->mean_confidence = confidence_sum / nsites;
        stats->min_confidence = confidence_min;
        stats->proposed_changes = proposed;
        stats->rejected_move_components = rejected_components;
        stats->rejected_move_sites = rejected_sites;
        stats->natural_refined_sites = natural_refined;
    }
    Arena_restore(arena, scratch);
    return 0;
}

int WindingMRF_solve_with_penalties(
    Arena_T arena, const WindingMRFSite *sites, size_t nsites,
    const WindingMRFEdge *edges, size_t nedges,
    const WindingMRFUnaryPenalty *penalties, size_t npenalties,
    const WindingMRFOptions *options,
    int32_t **out_labels, float **out_confidence, WindingMRFStats *stats)
{
    return wm_solve(arena,sites,nsites,edges,nedges,NULL,penalties,npenalties,
                    options,out_labels,out_confidence,stats);
}

int WindingMRF_solve_evidence(
    Arena_T arena, const WindingMRFSite *sites, size_t nsites,
    const WindingMRFEdge *edges, const int *kinds, size_t nedges,
    const WindingMRFOptions *options,
    int32_t **out_labels, float **out_confidence, WindingMRFStats *stats)
{
    return wm_solve(arena,sites,nsites,edges,nedges,kinds,NULL,0,
                    options,out_labels,out_confidence,stats);
}

int WindingMRF_solve(
    Arena_T arena,
    const WindingMRFSite *sites, size_t nsites,
    const WindingMRFEdge *edges, size_t nedges,
    const WindingMRFOptions *options,
    int32_t **out_labels, float **out_confidence,
    WindingMRFStats *stats)
{
    return WindingMRF_solve_with_penalties(
        arena, sites, nsites, edges, nedges, NULL, 0, options,
        out_labels, out_confidence, stats);
}

static void wm_selftest_check(int condition, const char *message, int *fails)
{
    if (condition) return;
    fprintf(stderr, "[winding MRF selftest] FAIL: %s\n", message);
    (*fails)++;
}

/* Seven-site extraction of the PHerc0139 4x5x5 Artifact C cycle at z=4608.
 * The numbers below are the natural-energy factors dumped by the production
 * 4x5x5 solve after adding the half-shifted axial order lattice.  Components
 * 158--204--227 form the inner ply, while 260--229--216--179 form the outer
 * ply.  The old collector saw only the continuation cycle and therefore
 * aliased both plies.  The shifted lattice contributes two independently
 * supported +1 order factors, 158--229 and 204--216.  Their combined evidence
 * should cut the weak 227--260 continuation (weight 25.57), not translate an
 * arbitrarily large part of the scroll.
 *
 * Keeping the measured cycle here makes GCO answer the same 2187-assignment
 * problem as the independent exact oracle in
 * python/scripts/winding_artifact_c_micro_mrf.py. */
static void wm_artifact_c_case(
    const char *name,
    int add_shifted_order,
    double gaussian_prior_weight,
    const int32_t expected[7],
    int *fails)
{
    WindingMRFSite sites[7];
    WindingMRFEdge edges[8];
    size_t nedges = 0;
    memset(sites, 0, sizeof sites);
    memset(edges, 0, sizeof edges);
    for (int i = 0; i < 7; i++) {
        sites[i].center = 0.0;
        sites[i].sigma = 1.0;
        sites[i].weight = gaussian_prior_weight;
        sites[i].initial_label = 0;
    }
    /* Site order: 158, 204, 227, 260, 229, 216, 179. */
    edges[nedges++] = (WindingMRFEdge){0, 1, 0, 42.681835867693736};
    edges[nedges++] = (WindingMRFEdge){1, 2, 0, 49.970723294303752};
    edges[nedges++] = (WindingMRFEdge){2, 3, 0, 25.571182009636686};
    edges[nedges++] = (WindingMRFEdge){4, 3, 0, 43.248565647354552};
    edges[nedges++] = (WindingMRFEdge){5, 4, 0, 42.802768553792704};
    edges[nedges++] = (WindingMRFEdge){6, 5, 0, 45.856226849047744};
    if (add_shifted_order) {
        edges[nedges++] = (WindingMRFEdge){1, 5, 1, 78.722653490961207};
        edges[nedges++] = (WindingMRFEdge){0, 4, 1, 65.95062021231567};
    }

    WindingMRFOptions options;
    WindingMRF_default_options(&options);
    options.label_min = -1;
    options.label_max = 1;
    options.unary_huber_delta = 100.0; /* exact Gaussian NLL in this range */
    Arena_T arena = Arena_new();
    int32_t *labels = NULL;
    WindingMRFStats stats;
    char message[160];
    snprintf(message, sizeof message, "Artifact C %s arena creates", name);
    wm_selftest_check(arena != NULL, message, fails);
    if (arena == NULL) return;
    int rc = WindingMRF_solve(
        arena, sites, 7, edges, nedges, &options, &labels, NULL, &stats);
    snprintf(message, sizeof message, "Artifact C %s solve runs", name);
    wm_selftest_check(rc == 0, message, fails);
    if (rc == 0 && labels != NULL) {
        int match = 1;
        for (int i = 0; i < 7; i++)
            if (labels[i] != expected[i]) match = 0;
        fprintf(stderr,
                 "[winding MRF Artifact C] %-24s labels="
                 "158:%d 204:%d 227:%d 260:%d 229:%d 216:%d 179:%d "
                 "E=%.9g margin=n/a\n",
                 name, labels[0], labels[1], labels[2], labels[3], labels[4],
                 labels[5], labels[6],
                 stats.energy_after);
        snprintf(message, sizeof message,
                 "Artifact C %s matches exact oracle", name);
        wm_selftest_check(match, message, fails);
    }
    Arena_dispose(&arena);
}

int WindingMRF_selftest(void)
{
    WindingMRFSite sites[6];
    WindingMRFEdge edges[7];
    memset(sites, 0, sizeof sites);
    memset(edges, 0, sizeof edges);
    for (int i = 0; i < 5; i++) {
        sites[i].center = i == 3 ? -0.5 : (double)i + 0.15;
        sites[i].sigma = 0.7;
        sites[i].weight = i == 3 ? 0.1 : 0.7;
        sites[i].initial_label = INT32_MIN;
    }
    sites[0].initial_label = 0;
    sites[0].fixed = 1;
    for (int i = 0; i < 4; i++) {
        edges[i].a = i;
        edges[i].b = i + 1;
        edges[i].target = 1;
        edges[i].weight = 3.0;
    }
    edges[4] = (WindingMRFEdge){0, 4, 4, 2.0};
    edges[5] = (WindingMRFEdge){4, 0, -4, 1.0}; /* reversed duplicate */
    edges[6] = (WindingMRFEdge){0, 4, 3, 0.15}; /* imperfect closer */
    sites[5].center = NAN;
    sites[5].sigma = 0.0;
    sites[5].weight = 0.0;
    sites[5].initial_label = 1;

    WindingMRFOptions options;
    WindingMRF_default_options(&options);
    options.label_min = -2;
    options.label_max = 6;
    Arena_T arena = Arena_new();
    int32_t *labels = NULL;
    float *confidence = NULL;
    WindingMRFStats stats;
    int fails = 0;
    wm_selftest_check(arena != NULL, "arena creates", &fails);
    if (arena != NULL) {
        wm_selftest_check(WindingMRF_solve(
            arena, sites, 6, edges, 7, &options,
            &labels, &confidence, &stats) == 0,
            "noisy loop solves", &fails);
        if (labels != NULL) {
            for (int i = 0; i < 5; i++)
                wm_selftest_check(labels[i] == i,
                                  "relative winding sequence recovers", &fails);
            wm_selftest_check(stats.energy_after <= stats.energy_before,
                              "optimization does not raise energy", &fails);
            wm_selftest_check(confidence[0] > 0.99f,
                              "fixed anchor is sharp", &fails);
            wm_selftest_check(confidence[5] < 0.01f,
                              "unobserved site abstains", &fails);
            wm_selftest_check(stats.abstained_sites >= 1,
                              "abstention is reported", &fails);
        }
        Arena_dispose(&arena);
    }
    arena = Arena_new();
    if (arena != NULL) {
        WindingMRFSite exclusion_sites[2];
        WindingMRFEdge exclusion_edge = {0, 1, 0, 1.0};
        WindingMRFUnaryPenalty exclusion = {1, 0, 100.0};
        memset(exclusion_sites, 0, sizeof exclusion_sites);
        exclusion_sites[0].center = 0.0;
        exclusion_sites[0].sigma = 1.0;
        exclusion_sites[0].weight = 1.0;
        exclusion_sites[0].initial_label = 0;
        exclusion_sites[0].fixed = 1;
        exclusion_sites[1].center = 0.0;
        exclusion_sites[1].sigma = 1.0;
        exclusion_sites[1].weight = 1.0;
        exclusion_sites[1].initial_label = 0;
        options.label_min = -2;
        options.label_max = 2;
        labels = NULL;
        wm_selftest_check(WindingMRF_solve_with_penalties(
            arena, exclusion_sites, 2, &exclusion_edge, 1,
            &exclusion, 1, &options, &labels, NULL, &stats) == 0,
            "sparse exclusion solve runs", &fails);
        if (labels != NULL)
            wm_selftest_check(labels[0] == 0 && labels[1] != 0,
                              "large unary excludes one occupied label",
                              &fails);
        Arena_dispose(&arena);
    }
    {
        const int32_t aliased[7] = {0, 0, 0, 0, 0, 0, 0};
        const int32_t separated[7] = {-1, -1, -1, 0, 0, 0, 0};
        wm_artifact_c_case(
            "old-continuation-only", 0, 32.02, aliased, &fails);
        wm_artifact_c_case(
            "shifted-order-flat", 1, 0.02, separated, &fails);
        wm_artifact_c_case(
            "shifted-order-prior32", 1, 32.02, separated, &fails);
    }
    {
        Arena_T scratch=Arena_new();
        WindingMRFSite order_sites[2]={{0}};
        WindingMRFEdge edge={0,1,1,100};
        WindingMRFOptions order_options;
        WindingMRFStats order_stats={0};
        int kind=WINDING_MRF_AT_LEAST;
        int32_t *order_labels=NULL;
        WindingMRF_default_options(&order_options);
        order_options.label_min=-4; order_options.label_max=4;
        order_sites[0].fixed=1;
        order_sites[1].center=3; order_sites[1].sigma=1; order_sites[1].weight=1;
        wm_selftest_check(WindingMRF_solve_evidence(scratch,order_sites,2,&edge,&kind,1,
            &order_options,&order_labels,NULL,&order_stats)==0 && order_labels && order_labels[1]==3,
            "one-sided order does not pull a valid separation to one turn",&fails);
        edge=(WindingMRFEdge){1,0,-1,100}; kind=WINDING_MRF_AT_MOST;
        wm_selftest_check(WindingMRF_solve_evidence(scratch,order_sites,2,&edge,&kind,1,
            &order_options,&order_labels,NULL,&order_stats)==0 && order_labels && order_labels[1]==3,
            "reversing endpoints reverses inequality direction",&fails);
        kind=WINDING_MRF_EQUAL;
        wm_selftest_check(WindingMRF_solve_evidence(scratch,order_sites,2,&edge,&kind,1,
            &order_options,&order_labels,NULL,&order_stats)==0 && order_labels && order_labels[1]==1,
            "equality remains a distinct legacy objective",&fails);
        edge=(WindingMRFEdge){0,1,2,100}; kind=WINDING_MRF_AT_LEAST;
        order_sites[1].center=-2;
        wm_selftest_check(WindingMRF_solve_evidence(scratch,order_sites,2,&edge,&kind,1,
            &order_options,&order_labels,NULL,&order_stats)==0 && order_labels && order_labels[1]==2,
            "violated order incurs the expected hinge penalty",&fails);
        kind=99;
        wm_selftest_check(WindingMRF_solve_evidence(scratch,order_sites,2,&edge,&kind,1,
            &order_options,&order_labels,NULL,&order_stats)!=0 && order_labels==NULL,
            "unknown factor family fails closed",&fails);
        Arena_dispose(&scratch);
    }
    {
        /* PHerc0139 21^3 round-7 witness: four degree-one charts, each
         * connected by one real factor below half a graph-cost quantum.
         * GCO moved all four by -4 every round, accumulating -24 before a
         * natural-energy refusal. An unrelated improving component must
         * survive rejection of those harmful zero-quantized tie moves. */
        const double weights[4]={.0017088167924956337,.0011167488587801287,
                                  .0019162146074757124,.0009111515882936122};
        for (int strong=0;strong<2;strong++) for (int drift=0;drift<2;drift++) {
            Arena_T scratch=Arena_new();
            WindingMRFSite leaf_sites[10]={{0}};
            WindingMRFEdge leaf_edges[5]={{0}};
            WindingMRFOptions leaf_options;
            WindingMRFStats leaf_stats={0};
            int32_t *leaf_labels=NULL;
            WindingMRF_default_options(&leaf_options);
            leaf_options.label_min=-4; leaf_options.label_max=4;
            for (int i=0;i<4;i++) {
                leaf_sites[2*i].fixed=1;
                leaf_edges[i]=(WindingMRFEdge){2*i,2*i+1,drift ? 24 : 0,weights[i]};
            }
            leaf_sites[8].fixed=1;
            leaf_edges[4]=(WindingMRFEdge){8,9,1,2};
            int rc=WindingMRF_solve_evidence(scratch,leaf_sites,strong ? 10 : 8,
                    leaf_edges,NULL,strong ? 5 : 4,&leaf_options,&leaf_labels,NULL,&leaf_stats);
            int stable=rc==0 && leaf_labels!=NULL;
            for (int i=0;stable && i<8;i++)
                if (leaf_labels[i]!=(drift && i%2 ? 4 : 0)) stable=0;
            if (strong && stable && (leaf_labels[8]!=0 || leaf_labels[9]!=1)) stable=0;
            fprintf(stderr,"[winding MRF selftest] quantized leaf drift=%d strong=%d: natural %.17g -> %.17g, %s\n",
                    drift,strong,leaf_stats.energy_before,leaf_stats.energy_after,stable ? "stable" : "FAIL");
            wm_selftest_check(stable,"original-weight conditional solves prevent or reverse invisible leaf drift without vetoing a strong improvement",&fails);
            Arena_dispose(&scratch);
        }
    }
    {
        /* v28's independent full-scroll residual audit found that accepting
         * changed components alone still lets a weak leaf ride along with a
         * temporarily moving strong neighbor. Its exact conditional must
         * follow that neighbor in EITHER direction, not GCO's lowest tie. */
        for (int sign=-1;sign<=1;sign+=2) {
            Arena_T scratch=Arena_new();
            WindingMRFSite moving_sites[3]={{0}};
            WindingMRFEdge moving_edges[2]={{0,1,0,2},{1,2,0,.0009111515882936122}};
            WindingMRFOptions moving_options;
            WindingMRFStats moving_stats={0};
            int32_t *moving_labels=NULL;
            WindingMRF_default_options(&moving_options);
            moving_options.label_min=-4; moving_options.label_max=4;
            moving_sites[0].fixed=1; moving_edges[0].target=sign;
            int rc=WindingMRF_solve_evidence(scratch,moving_sites,3,moving_edges,NULL,2,
                    &moving_options,&moving_labels,NULL,&moving_stats);
            wm_selftest_check(rc==0 && moving_labels && moving_labels[0]==0 &&
                moving_labels[1]==sign && moving_labels[2]==sign &&
                moving_stats.natural_refined_sites>0 && moving_stats.energy_after==0,
                "invisible weak leaf follows a moving supported neighbor without accumulated drift",&fails);
            Arena_dispose(&scratch);
        }
    }
    fprintf(stderr, "[winding MRF selftest] %s (%d failure(s))\n",
            fails == 0 ? "PASS" : "FAIL", fails);
    return fails;
}
