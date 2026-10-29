/* asm_run.c -- growable arena tables of the run. */
#include <string.h>

#include "asm_types.h"

#define ASM_GROW(run, arr, n, cap, T, init)                                   \
    do {                                                                      \
        if ((n) == (cap)) {                                                   \
            size_t ncap = (cap) ? (cap) * 2 : (init);                         \
            T *na = ARENA_ALLOC((run)->arena, ncap * sizeof(T));              \
            if ((n)) memcpy(na, (arr), (n) * sizeof(T));                      \
            (arr) = na; (cap) = ncap;                                         \
        }                                                                     \
    } while (0)

size_t AsmRun_push_chart(AsmRun *run, const AsmChart *c)
{
    ASM_GROW(run, run->charts, run->n_charts, run->cap_charts, AsmChart, 1024);
    run->charts[run->n_charts] = *c;
    return run->n_charts++;
}

size_t AsmRun_push_rel(AsmRun *run, const AsmRelation *r)
{
    ASM_GROW(run, run->rels, run->n_rels, run->cap_rels, AsmRelation, 4096);
    run->rels[run->n_rels] = *r;
    return run->n_rels++;
}

size_t AsmRun_push_corr(AsmRun *run, const AsmCorr *c)
{
    ASM_GROW(run, run->corr, run->n_corr, run->cap_corr, AsmCorr, 65536);
    run->corr[run->n_corr] = *c;
    return run->n_corr++;
}

size_t AsmRun_push_layer(AsmRun *run, const AsmLayerPair *l)
{
    ASM_GROW(run, run->layers, run->n_layers, run->cap_layers, AsmLayerPair, 4096);
    run->layers[run->n_layers] = *l;
    return run->n_layers++;
}

size_t AsmRun_push_layer_hit(AsmRun *run, const AsmLayerHit *h)
{
    ASM_GROW(run, run->layer_hits, run->n_layer_hits, run->cap_layer_hits, AsmLayerHit, 16384);
    run->layer_hits[run->n_layer_hits] = *h;
    return run->n_layer_hits++;
}
