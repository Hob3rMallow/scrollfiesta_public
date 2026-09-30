#include "material_evidence.h"
#include "../common/pipeline_constants.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct {
    double position, quality;
    size_t interval;
    int start;
} MeEvent;

typedef struct {
    int32_t a, b, target;
    int kind;
    double length, support;
    uint64_t face_a, face_b;
    uint64_t region_a, region_b;
    double budget;
} MeMass;

typedef struct { int32_t chart; uint64_t face; size_t mass; } MeFaceUse;
typedef struct { size_t factor; int side; uint64_t region; } MeRegionUse;

static int me_face_compare(const void *lhs, const void *rhs)
{
    const MeFaceUse *a=lhs, *b=rhs;
    if (a->chart!=b->chart) return a->chart<b->chart ? -1 : 1;
    if (a->face!=b->face) return a->face<b->face ? -1 : 1;
    return a->mass<b->mass ? -1 : a->mass>b->mass;
}

static int me_region_compare(const void *lhs,const void *rhs)
{
    const MeRegionUse *a=lhs,*b=rhs;
    if (a->factor!=b->factor) return a->factor<b->factor ? -1 : 1;
    if (a->side!=b->side) return a->side<b->side ? -1 : 1;
    return a->region<b->region ? -1 : a->region>b->region;
}

/* Apply a shared influence budget across all hypotheses and all partners on
 * BOTH sources. Repeating this for correlated physical blocks can only reduce
 * a face's already bounded influence. No local contradiction becomes a veto. */
static void me_cap(Arena_T arena,MeMass *mass,size_t nm,int regions)
{
    Arena_Mark mark=Arena_save(arena);
    MeFaceUse *uses=ARENA_ALLOC(arena,2*nm*sizeof *uses);
    for (size_t i=0;i<nm;i++) {
        mass[i].budget=1;
        uses[2*i]=(MeFaceUse){mass[i].a,regions ? mass[i].region_a : mass[i].face_a,i};
        uses[2*i+1]=(MeFaceUse){mass[i].b,regions ? mass[i].region_b : mass[i].face_b,i};
    }
    qsort(uses,2*nm,sizeof *uses,me_face_compare);
    for (size_t i=0;i<2*nm;) {
        size_t j=i;
        double total=0;
        while (j<2*nm && uses[j].chart==uses[i].chart && uses[j].face==uses[i].face) {
            total+=mass[uses[j].mass].support; j++;
        }
        for (size_t k=i;k<j;k++)
            mass[uses[k].mass].budget=fmax(mass[uses[k].mass].budget,total);
        i=j;
    }
    for (size_t i=0;i<nm;i++) mass[i].support/=mass[i].budget;
    Arena_restore(arena,mark);
}

static int me_interval_compare(const void *lhs, const void *rhs)
{
    const MaterialEvidenceInterval *a=lhs, *b=rhs;
    if (a->unit!=b->unit) return a->unit<b->unit ? -1 : 1;
    if (a->a!=b->a) return a->a<b->a ? -1 : 1;
    if (a->b!=b->b) return a->b<b->b ? -1 : 1;
    if (a->kind!=b->kind) return a->kind<b->kind ? -1 : 1;
    if (a->target!=b->target) return a->target<b->target ? -1 : 1;
    if (a->begin!=b->begin) return a->begin<b->begin ? -1 : 1;
    if (a->end!=b->end) return a->end<b->end ? -1 : 1;
    if (a->quality!=b->quality) return a->quality>b->quality ? -1 : 1;
    return 0;
}

static int me_event_compare(const void *lhs, const void *rhs)
{
    const MeEvent *a=lhs, *b=rhs;
    if (a->position!=b->position) return a->position<b->position ? -1 : 1;
    if (a->start!=b->start) return a->start<b->start ? -1 : 1;
    if (a->quality!=b->quality) return a->quality>b->quality ? -1 : 1;
    return a->interval<b->interval ? -1 : a->interval>b->interval;
}

static int me_mass_compare(const void *lhs, const void *rhs)
{
    const MeMass *a=lhs, *b=rhs;
    if (a->a!=b->a) return a->a<b->a ? -1 : 1;
    if (a->b!=b->b) return a->b<b->b ? -1 : 1;
    if (a->kind!=b->kind) return a->kind<b->kind ? -1 : 1;
    if (a->target!=b->target) return a->target<b->target ? -1 : 1;
    /* Stable numerical sum under input permutation. */
    if (a->support!=b->support) return a->support<b->support ? -1 : 1;
    if (a->length!=b->length) return a->length<b->length ? -1 : 1;
    return 0;
}

/* Sweep the quality upper envelope using a max heap and lazy deletion. */
static void me_integrate(Arena_T arena, const MaterialEvidenceInterval *v,
                         size_t n, double *length, double *quality_length)
{
    Arena_Mark mark=Arena_save(arena);
    MeEvent *event=ARENA_ALLOC(arena,2*n*sizeof *event);
    size_t *heap=ARENA_ALLOC(arena,n*sizeof *heap);
    uint8_t *active=ARENA_CALLOC(arena,n,sizeof *active);
    size_t nh=0;
    double previous=v[0].begin;
    *length=0; *quality_length=0;
    for (size_t i=0;i<n;i++) {
        event[2*i]=(MeEvent){v[i].begin,v[i].quality,i,1};
        event[2*i+1]=(MeEvent){v[i].end,v[i].quality,i,0};
    }
    qsort(event,2*n,sizeof *event,me_event_compare);
    for (size_t i=0;i<2*n;) {
        double position=event[i].position;
        if (nh) {
            double delta=position-previous;
            *length+=delta;
            *quality_length+=delta*v[heap[0]].quality;
        }
        while (i<2*n && event[i].position==position) {
            size_t k=event[i].interval;
            active[k]=(uint8_t)event[i].start;
            if (event[i].start) {
                size_t p=nh++;
                while (p && v[heap[(p-1)/2]].quality<v[k].quality) {
                    heap[p]=heap[(p-1)/2]; p=(p-1)/2;
                }
                heap[p]=k;
            }
            i++;
        }
        while (nh && !active[heap[0]]) {
            size_t k=heap[--nh], p=0;
            if (!nh) break;
            while (2*p+1<nh) {
                size_t c=2*p+1;
                if (c+1<nh && v[heap[c+1]].quality>v[heap[c]].quality) c++;
                if (v[k].quality>=v[heap[c]].quality) break;
                heap[p]=heap[c]; p=c;
            }
            heap[p]=k;
        }
        previous=position;
    }
    Arena_restore(arena,mark);
}

int MaterialEvidence_reduce(Arena_T arena,
                            const MaterialEvidenceInterval *intervals, size_t n,
                            double length_scale, MaterialEvidenceFactor **out,
                            size_t *count, MaterialEvidenceReport *report)
{
    MaterialEvidenceReport local={0};
    MaterialEvidenceInterval *v=NULL;
    MeMass *mass=NULL;
    size_t nm=0, nf=0;
    if (!out || !count) return -1;
    *out=NULL; *count=0;
    if (!report) report=&local;
    memset(report,0,sizeof *report);
    if (!arena || (n && !intervals) || !isfinite(length_scale) || length_scale<=0 ||
        n>SIZE_MAX/(2*sizeof(MeEvent)) || n>SIZE_MAX/sizeof *v ||
        n>SIZE_MAX/sizeof(MaterialEvidenceFactor) || n>SIZE_MAX/sizeof(MeMass) ||
        n>SIZE_MAX/(2*sizeof(MeFaceUse)) || n>SIZE_MAX/(2*sizeof(MeRegionUse))) return -1;
    report->intervals=n;
    if (!n) return 0;
    v=ARENA_ALLOC(arena,n*sizeof *v);
    mass=ARENA_ALLOC(arena,n*sizeof *mass);
    memcpy(v,intervals,n*sizeof *v);
    for (size_t i=0;i<n;i++) {
        MaterialEvidenceInterval *s=v+i;
        if (s->a<0 || s->b<0 || s->a==s->b || s->target==INT32_MIN ||
            s->kind<WINDING_MRF_EQUAL || s->kind>WINDING_MRF_AT_MOST ||
            !isfinite(s->begin) || !isfinite(s->end) || !(s->end>s->begin) ||
            !isfinite(s->end-s->begin) || !isfinite(s->quality) ||
            s->quality<0 || s->quality>1 ||
            !isfinite(s->begin_b) || !isfinite(s->end_b)) return -1;
        if (s->a>s->b) {
            int32_t t=s->a; s->a=s->b; s->b=t; s->target=-s->target;
            uint64_t face=s->source_face_a; s->source_face_a=s->source_face_b; s->source_face_b=face;
            uint64_t edge=s->source_edge_a; s->source_edge_a=s->source_edge_b; s->source_edge_b=edge;
            uint64_t region=s->region_a; s->region_a=s->region_b; s->region_b=region;
            if (s->kind==WINDING_MRF_AT_LEAST) s->kind=WINDING_MRF_AT_MOST;
            else if (s->kind==WINDING_MRF_AT_MOST) s->kind=WINDING_MRF_AT_LEAST;
        }
    }
    qsort(v,n,sizeof *v,me_interval_compare);
    for (size_t first=0;first<n;) {
        size_t last=first+1, begin_mass=nm;
        double total=0;
        while (last<n && v[last].unit==v[first].unit) {
            if (v[last].a!=v[first].a || v[last].b!=v[first].b ||
                v[last].kind!=v[first].kind ||
                v[last].source_face_a!=v[first].source_face_a ||
                v[last].source_face_b!=v[first].source_face_b ||
                v[last].source_edge_a!=v[first].source_edge_a || v[last].source_edge_b!=v[first].source_edge_b ||
                v[last].region_a!=v[first].region_a || v[last].region_b!=v[first].region_b) return -1;
            last++;
        }
        report->units++;
        for (size_t i=first;i<last;) {
            size_t j=i+1;
            double length=0, weighted=0;
            while (j<last && v[j].target==v[i].target) j++;
            me_integrate(arena,v+i,j-i,&length,&weighted);
            if (!isfinite(length) || !isfinite(weighted)) return -1;
            if (weighted>0) {
                double support=fmin(weighted/length_scale,1.0);
                mass[nm++]=(MeMass){v[i].a,v[i].b,v[i].target,v[i].kind,length,support,
                                    v[i].source_face_a,v[i].source_face_b,v[i].region_a,v[i].region_b,1.0};
                total+=support;
            }
            i=j;
        }
        /* Densifying or duplicating conflicting hypotheses cannot multiply a
         * unit's influence. Do not discard minority evidence or hard-veto. */
        for (size_t i=begin_mass;i<nm;i++) mass[i].support/=fmax(total,1.0);
        first=last;
    }
    me_cap(arena,mass,nm,0);
    me_cap(arena,mass,nm,1);
    qsort(mass,nm,sizeof *mass,me_mass_compare);
    MaterialEvidenceFactor *result=ARENA_CALLOC(arena,nm,sizeof *result);
    MeRegionUse *regions=ARENA_ALLOC(arena,2*nm*sizeof *regions);
    for (size_t i=0;i<nm;) {
        size_t j=i;
        MaterialEvidenceFactor *f=result+nf++;
        f->edge.a=mass[i].a; f->edge.b=mass[i].b;
        f->edge.target=mass[i].target; f->kind=mass[i].kind;
        while (j<nm && mass[j].a==mass[i].a && mass[j].b==mass[i].b &&
               mass[j].kind==mass[i].kind && mass[j].target==mass[i].target) {
            f->covered_length+=mass[j].length;
            f->effective_support+=mass[j].support;
            regions[2*j]=(MeRegionUse){nf-1,0,mass[j].region_a};
            regions[2*j+1]=(MeRegionUse){nf-1,1,mass[j].region_b};
            j++;
        }
        if (!isfinite(f->covered_length) || !isfinite(f->effective_support)) return -1;
        f->edge.weight=f->effective_support;
        report->covered_length+=f->covered_length;
        report->effective_support+=f->effective_support;
        i=j;
    }
    qsort(regions,2*nm,sizeof *regions,me_region_compare);
    for (size_t i=0;i<2*nm;i++) {
        const MeRegionUse *r=regions+i;
        if (i && !me_region_compare(r,regions+i-1)) continue;
        if (r->side==0) result[r->factor].regions_a++;
        else result[r->factor].regions_b++;
    }
    /* Equality factors are sorted by pair, kind, target. Their summed L1 score
     * is convex; its two neighboring-label differences give the exact local
     * preference margin. A removed unit changes either difference by <=1; an
     * adversarial replacement by <=2. No assumption of independent trials is
     * needed for this bounded-influence statement on the reduced scores. */
    for (size_t first=0;first<nf;) {
        size_t last=first+1;
        double total=result[first].effective_support, below=0;
        while (last<nf && result[last].edge.a==result[first].edge.a &&
               result[last].edge.b==result[first].edge.b && result[last].kind==result[first].kind) {
            total+=result[last].effective_support; last++;
        }
        if (result[first].kind==WINDING_MRF_EQUAL) for (size_t i=first;i<last;i++) {
            double mass_i=result[i].effective_support, above=total-below-mass_i;
            result[i].pair_margin=mass_i-fabs(below-above);
            below+=mass_i;
        }
        first=last;
    }
    for (size_t i=0;i<nf;i++) {
        MaterialEvidenceFactor *f=result+i;
        f->units=f->regions_a<f->regions_b ? f->regions_a : f->regions_b;
        f->significant=f->units>=MATERIAL_EVIDENCE_MIN_UNITS &&
                       f->effective_support>=MATERIAL_EVIDENCE_MIN_SUPPORT &&
                       (f->kind!=WINDING_MRF_EQUAL ||
                        f->pair_margin>2*MATERIAL_EVIDENCE_CONTAMINATION_BUDGET);
        report->significant_factors+=(size_t)f->significant;
    }
    report->factors=nf;
    *out=result; *count=nf;
    return 0;
}

int MaterialEvidence_selftest(void)
{
    Arena_T arena=Arena_new();
    MaterialEvidenceInterval v[80]={0};
    MaterialEvidenceFactor *f=NULL;
    MaterialEvidenceReport report={0};
    size_t n=0;
    int fail=0;
    /* Six independent pieces of paper support target 1. One bad source
     * interval at target 0 is copied 64 times: it still has only one vote. */
    for (int i=0;i<6;i++) v[i]=(MaterialEvidenceInterval){0,1,1,0,(uint64_t)i,0,8,1};
    for (int i=6;i<70;i++) v[i]=(MaterialEvidenceInterval){0,1,0,0,6,0,8,1};
    for (int i=0;i<70;i++) {
        v[i].source_face_a=(uint64_t)(i<6 ? i : 6);
        v[i].source_face_b=(uint64_t)(i<6 ? i : 6);
        v[i].region_a=v[i].region_b=(uint64_t)(i<6 ? i : 6);
    }
    fail+=MaterialEvidence_reduce(arena,v,70,8,&f,&n,&report)!=0;
    if (n!=2) fail++;
    else {
        fail+=fabs(f[0].edge.weight-1)>1e-12 || f[0].significant;
        fail+=fabs(f[1].edge.weight-6)>1e-12 || !f[1].significant;
        WindingMRFSite sites[2]={{0}};
        WindingMRFEdge edges[2]={f[0].edge,f[1].edge};
        int kinds[2]={f[0].kind,f[1].kind};
        int32_t *labels=NULL;
        WindingMRFOptions options;
        WindingMRF_default_options(&options);
        options.label_min=-2; options.label_max=2; sites[0].fixed=1;
        fail+=WindingMRF_solve_evidence(arena,sites,2,edges,kinds,2,&options,&labels,NULL,NULL)!=0;
        fail+=!labels || labels[1]!=1;
        edges[0].target=-100000;
        WindingMRFStats extreme_stats={0};
        fail+=WindingMRF_solve_evidence(arena,sites,2,edges,kinds,2,&options,&labels,NULL,&extreme_stats)!=0;
        fail+=!labels || labels[1]!=1;
        fail+=fabs(extreme_stats.energy_before-100006)>1e-9 || fabs(extreme_stats.energy_after-100001)>1e-9;
        fail+=extreme_stats.quantized_energy_before-extreme_stats.quantized_energy_after!=1280;
        edges[0].target=100000; kinds[0]=WINDING_MRF_AT_LEAST;
        fail+=WindingMRF_solve_evidence(arena,sites,2,edges,kinds,2,&options,&labels,NULL,&extreme_stats)!=0;
        fail+=!labels || labels[1]!=1 || fabs(extreme_stats.energy_after-99999)>1e-9;
    }
    /* Arbitrarily subdividing a support interval also leaves its mass alone. */
    for (int i=6;i<70;i++) {
        v[i].begin=(double)(i-6)/8; v[i].end=(double)(i-5)/8;
    }
    fail+=MaterialEvidence_reduce(arena,v,70,8,&f,&n,&report)!=0;
    fail+=n!=2 || fabs(f[0].edge.weight-1)>1e-12 || fabs(f[1].edge.weight-6)>1e-12;
    /* Many DISTINCT original faces still describe one correlated short
     * contact. Their physical block has at most one unit, on either side. */
    for (int i=6;i<70;i++) {
        v[i].unit=(uint64_t)i; v[i].begin=0; v[i].end=8;
        v[i].source_face_a=(uint64_t)i;
        v[i].source_face_b=(uint64_t)i;
        v[i].region_b=(uint64_t)i;
    }
    fail+=MaterialEvidence_reduce(arena,v,70,8,&f,&n,&report)!=0;
    fail+=n!=2 || fabs(f[0].edge.weight-1)>1e-12 || f[0].units!=1 || f[0].significant;
    for (int i=6;i<70;i++) { v[i].region_a=(uint64_t)i; v[i].region_b=6; }
    fail+=MaterialEvidence_reduce(arena,v,70,8,&f,&n,&report)!=0;
    fail+=n!=2 || fabs(f[0].edge.weight-1)>1e-12 || f[0].units!=1 || f[0].significant;
    /* A single corrupted face spanning distinct physical blocks also stays
     * bounded. Neither face IDs nor block IDs alone are sufficient. */
    for (int i=6;i<70;i++) { v[i].source_face_a=6; v[i].region_b=(uint64_t)i; }
    fail+=MaterialEvidence_reduce(arena,v,70,8,&f,&n,&report)!=0;
    fail+=n!=2 || fabs(f[0].edge.weight-1)>1e-12 || f[0].significant;
    /* Even assigning each copied sample a distinct interval cannot give its
     * one ORIGINAL bad face more than one unit of total MRF influence. */
    for (int i=6;i<70;i++) {
        v[i].unit=(uint64_t)i; v[i].begin=0; v[i].end=8;
        v[i].source_face_b=(uint64_t)i;
    }
    fail+=MaterialEvidence_reduce(arena,v,70,8,&f,&n,&report)!=0;
    fail+=n!=2 || fabs(f[0].edge.weight-1)>1e-12 || fabs(f[1].edge.weight-6)>1e-12;
    /* Conflicting targets in the SAME support unit share a budget. */
    v[0]=(MaterialEvidenceInterval){0,1,0,0,0,0,8,1};
    v[1]=(MaterialEvidenceInterval){0,1,1,0,0,0,8,1};
    fail+=MaterialEvidence_reduce(arena,v,2,8,&f,&n,&report)!=0;
    fail+=n!=2 || fabs(f[0].edge.weight-.5)>1e-12 || fabs(f[1].edge.weight-.5)>1e-12;
    /* Partial overlaps use max quality, not their sum: integral = 6. */
    v[0]=(MaterialEvidenceInterval){0,1,0,0,0,0,8,.5};
    v[1]=(MaterialEvidenceInterval){0,1,0,0,0,4,8,1};
    fail+=MaterialEvidence_reduce(arena,v,2,8,&f,&n,&report)!=0;
    fail+=n!=1 || fabs(f[0].edge.weight-.75)>1e-12 || fabs(f[0].covered_length-8)>1e-12;
    /* Substantial material support on BOTH incompatible targets is ambiguity,
     * not two significant accepted links. A two-unit margin is not robust to
     * an adversarial replacement of one unit. Minority factors still survive. */
    for (int i=0;i<10;i++) {
        v[i]=(MaterialEvidenceInterval){0,1,i<4 ? 0 : 1,0,(uint64_t)i,0,8,1};
        v[i].source_face_a=v[i].source_face_b=v[i].region_a=v[i].region_b=(uint64_t)i;
    }
    fail+=MaterialEvidence_reduce(arena,v,10,8,&f,&n,&report)!=0;
    fail+=n!=2 || f[0].significant || f[1].significant || fabs(f[1].pair_margin-2)>1e-12;
    fail+=MaterialEvidence_reduce(arena,v+1,9,8,&f,&n,&report)!=0;
    fail+=n!=2 || f[0].significant || !f[1].significant || fabs(f[1].pair_margin-3)>1e-12;
    /* Competing PARTNERS share the same source-region budget too, even if
     * every partner touches a different original triangle in that region. */
    for (int i=0;i<8;i++) {
        v[i]=(MaterialEvidenceInterval){0,i+1,0,0,(uint64_t)i,0,8,1};
        v[i].source_face_a=(uint64_t)i; v[i].source_face_b=(uint64_t)(100+i);
        v[i].region_a=0; v[i].region_b=(uint64_t)i;
    }
    fail+=MaterialEvidence_reduce(arena,v,8,8,&f,&n,&report)!=0;
    fail+=n!=8 || fabs(report.effective_support-1)>1e-12 || report.significant_factors!=0;
    v[0].quality=NAN;
    fail+=MaterialEvidence_reduce(arena,v,1,8,&f,&n,&report)==0 || f!=NULL || n!=0;
    fail+=MaterialEvidence_reduce(arena,NULL,0,8,&f,&n,&report)!=0 || f!=NULL || n!=0;
    fprintf(stderr,"[selftest] material evidence %s (%d failures)\n",fail ? "FAIL" : "PASS",fail);
    Arena_dispose(&arena);
    return fail ? -1 : 0;
}
