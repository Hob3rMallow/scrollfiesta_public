#include "material_observation.h"
#include "material_front.h"
#include "scroll_coordinate.h"
#include "../common/union_find.h"
#include "../common/ves_platform.h"
#include "../common/pipeline_constants.h"

#include <limits.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    MaterialFront *front;
    size_t nfront, charts, vertices, faces;
    int32_t base;
    uint64_t hash;
} MoCube;

typedef struct { size_t begin,end,next; } MoRange;

static int mo_directory(const char *directory)
{
    char path[1600];
    size_t n=strlen(directory);
    if (!n || n>=sizeof path) return -1;
    memcpy(path,directory,n+1);
    for (size_t i=1;i<=n;i++) {
        char c=path[i];
        if (c!='/' && c!='\\' && c!=0) continue;
        if (path[i-1]==':') continue;
        path[i]=0;
        if (ves_mkdir(path)!=0 && errno!=EEXIST) return -1;
        path[i]=c;
    }
    return 0;
}

static int mo_filter(void *context,long z,long y,long x)
{
    const long *b=context;
    return !b || (z>=b[0] && z<=b[1] && y>=b[2] && y<=b[3] && x>=b[4] && x<=b[5]);
}

static uint64_t mo_hash(uint64_t hash,const void *data,size_t bytes)
{
    const uint8_t *p=data;
    for (size_t i=0;i<bytes;i++) { hash^=p[i]; hash*=UINT64_C(1099511628211); }
    return hash;
}

/* Preserve the exact proposed graph step before any repair. These files are
 * failure diagnostics, never a completed evidence report or accepted phase. */
static int mo_phase_failure(const char *directory,int round,size_t nc,
                             const WindingMRFSite *sites,const int64_t *phase,
                             const int32_t *labels,size_t nf,const WindingMRFEdge *factors,
                             const int *kinds,const WindingMRFOptions *options,
                             const WindingMRFStats *stats,int solver_rc)
{
    char path[1600];
    FILE *file=NULL;
    int rc=-1;
    snprintf(path,sizeof path,"%s/material_phase_failure_sites.tsv",directory);
    file=fopen(path,"wb");
    if (!file) return -1;
    fprintf(file,"site\tphase\tfixed\tproposed_increment\n");
    for (size_t i=0;i<nc;i++)
        fprintf(file,"%zu\t%lld\t%d\t%d\n",i,(long long)phase[i],sites[i].fixed,labels ? labels[i] : 0);
    rc=ferror(file) ? -1 : 0;
    if (fclose(file)!=0) rc=-1;
    if (rc!=0) return -1;
    snprintf(path,sizeof path,"%s/material_phase_failure_factors.tsv",directory);
    file=fopen(path,"wb");
    if (!file) return -1;
    fprintf(file,"a\tb\ttarget\tweight\tkind\n");
    for (size_t i=0;i<nf;i++)
        fprintf(file,"%d\t%d\t%d\t%.17g\t%d\n",factors[i].a,factors[i].b,
                factors[i].target,factors[i].weight,kinds[i]);
    rc=ferror(file) ? -1 : 0;
    if (fclose(file)!=0) rc=-1;
    if (rc!=0) return -1;
    snprintf(path,sizeof path,"%s/material_phase_failure.json",directory);
    file=fopen(path,"wb");
    if (!file) return -1;
    fprintf(file,"{\"schema\":\"material-phase-failure-v1\",\"complete\":false,\"round\":%d,"
                 "\"solver_rc\":%d,\"sites\":%zu,\"factors\":%zu,\"label_min\":%d,\"label_max\":%d,"
                 "\"cost_scale\":%.17g,\"energy_before\":%.17g,\"energy_after\":%.17g,"
                 "\"quantized_before\":%lld,\"quantized_after\":%lld,\"changed_labels\":%zu}\n",
            round,solver_rc,nc,nf,options->label_min,options->label_max,options->cost_scale,
            stats->energy_before,stats->energy_after,(long long)stats->quantized_energy_before,
            (long long)stats->quantized_energy_after,stats->changed_labels);
    rc=ferror(file) ? -1 : 0;
    if (fclose(file)!=0) rc=-1;
    return rc;
}

static int mo_phase_graph(Arena_T arena,const char *directory,size_t nc,const MaterialEvidenceFactor *f,size_t nf,
                          int64_t **out)
{
    if (nc>INT32_MAX) return -1;
    ScrollCoordinateEdge *edges=ARENA_ALLOC(arena,nf*sizeof *edges);
    ScrollCoordinateReport report={0};
    UnionFind graph=UF_new(arena,(int32_t)nc);
    for (size_t i=0;i<nf;i++) {
        if (f[i].kind!=WINDING_MRF_EQUAL) return -1;
        edges[i]=(ScrollCoordinateEdge){f[i].edge.a,f[i].edge.b,f[i].edge.target,f[i].edge.weight,0};
        uf_union(&graph,edges[i].a,edges[i].b);
    }
    if (ScrollCoordinate_solve(arena,nc,edges,nf,out,&report)!=0) return -1;
    fprintf(stderr,"[material observation] phase initialization only: %zu charts, %zu factors, %zu islands\n",
            nc,nf,report.gauge_components);
    for (int round=0;round<SCROLL_MODEL_MRF_ROUNDS;round++) {
        Arena_Mark mark=Arena_save(arena);
        WindingMRFSite *sites=ARENA_CALLOC(arena,nc,sizeof *sites);
        WindingMRFEdge *factors=ARENA_ALLOC(arena,nf*sizeof *factors);
        int *kinds=ARENA_ALLOC(arena,nf*sizeof *kinds);
        WindingMRFOptions options;
        WindingMRFStats stats={0};
        int32_t *labels=NULL;
        WindingMRF_default_options(&options);
        options.label_min=-SCROLL_MODEL_MRF_WINDOW; options.label_max=SCROLL_MODEL_MRF_WINDOW;
        for (size_t i=0;i<nc;i++) sites[i].fixed=uf_find(&graph,(int32_t)i)==(int32_t)i;
        for (size_t i=0;i<nf;i++) {
            int64_t target=(int64_t)f[i].edge.target-((*out)[f[i].edge.b]-(*out)[f[i].edge.a]);
            if (target<=INT32_MIN || target>INT32_MAX) return -1;
            factors[i]=f[i].edge; factors[i].target=(int32_t)target; kinds[i]=f[i].kind;
        }
        int solver_rc=WindingMRF_solve_evidence(arena,sites,nc,factors,kinds,nf,&options,&labels,NULL,&stats);
        if (solver_rc!=0 || stats.energy_after>stats.energy_before+1e-8*fmax(1,stats.energy_before)) {
            fprintf(stderr,"[material observation] MRF %d refused: solver_rc=%d; natural %.17g -> %.17g; quantized %lld -> %lld; %zu proposed changes\n",
                    round+1,solver_rc,stats.energy_before,stats.energy_after,
                    (long long)stats.quantized_energy_before,(long long)stats.quantized_energy_after,stats.changed_labels);
            if (mo_phase_failure(directory,round+1,nc,sites,*out,solver_rc==0 ? labels : NULL,
                                  nf,factors,kinds,&options,&stats,solver_rc)!=0)
                fprintf(stderr,"[material observation] could not complete MRF failure capture\n");
            return -1;
        }
        for (size_t i=0;i<nc;i++) (*out)[i]+=labels[i];
        fprintf(stderr,"[material observation] MRF %d: %zu changes; energy %.9g -> %.9g; natural-refined %zu invisible sites; rejected %zu component(s)/%zu proposed site moves\n",
                round+1,stats.changed_labels,stats.energy_before,stats.energy_after,
                stats.natural_refined_sites,
                stats.rejected_move_components,stats.rejected_move_sites);
        Arena_restore(arena,mark);
        if (!stats.changed_labels) return 0;
    }
    fprintf(stderr,"[material observation] phase MRF stopped at iteration cap; diagnostic only\n");
    return 0;
}

int MaterialObservation_build(Arena_T arena,const char *source_dir,const char *output_dir,
                               const long *bounds,const ScrollModelOptions *options)
{
    MeshPileEntry *pile=NULL;
    MoCube *cube=NULL;
    size_t npile=0,nfront=0,nchart=0,ninterval=0,nfactor=0;
    uint64_t total_vertices=0,total_faces=0;
    MaterialFront *fronts=NULL;
    MaterialEvidenceInterval *intervals=NULL;
    MaterialEvidenceFactor *factors=NULL;
    MaterialFrontReport front_report={0};
    MaterialEvidenceReport evidence_report={0};
    int64_t *phase=NULL;
    char path[1600];
    FILE *ledger=NULL,*table=NULL;
    double started=ves_clock_sec();
    int rc=-1;
    if (!arena || !source_dir || !output_dir || !options || strlen(output_dir)>1300 ||
        (bounds && (bounds[0]>bounds[1] || bounds[2]>bounds[3] || bounds[4]>bounds[5]))) return -1;
    if (mo_directory(output_dir)!=0) {
        fprintf(stderr,"[material observation] cannot create output directory %s\n",output_dir);
        return -1;
    }
    pile=ARENA_ALLOC(arena,SCROLL_MODEL_MAX_CUBES*sizeof *pile);
    if (MeshPile_scan(source_dir,pile,SCROLL_MODEL_MAX_CUBES,&npile,mo_filter,(void *)bounds)!=0 || !npile) return -1;
    cube=ARENA_CALLOC(arena,npile,sizeof *cube);
    snprintf(path,sizeof path,"%s/material_source_ledger.json",output_dir);
    ledger=fopen(path,"wb");
    if (!ledger) return -1;
    fprintf(ledger,"{\"schema\":\"material-source-ledger-v1\",\"stage\":\"source-and-soft-evidence-only\","
                   "\"physical_qualification\":\"INCOMPLETE\",\"halo_cubes\":0,\"cubes\":[\n");
    for (size_t i=0;i<npile;i++) {
        Arena_T scratch=Arena_new();
        MeshBinData mesh={0};
        ScrollSource source={0};
        MaterialFront *part=NULL;
        size_t npart=0;
        int ok=0;
        if (!pile[i].has_id || MeshBin_read_arena(scratch,pile[i].path,&mesh)!=0 ||
            ScrollSource_build(scratch,&mesh,options->axis,options->axis_y,options->axis_x,
                               options->sense,options->core_radius,&source)!=0 ||
            source.ncomponents>(size_t)INT32_MAX-nchart) goto cube_done;
        cube[i].base=(int32_t)nchart; cube[i].charts=source.ncomponents;
        cube[i].vertices=mesh.nv; cube[i].faces=mesh.nf;
        cube[i].hash=mo_hash(UINT64_C(14695981039346656037),mesh.verts,3*mesh.nv*sizeof(float));
        cube[i].hash=mo_hash(cube[i].hash,mesh.faces,3*mesh.nf*sizeof(int32_t));
        if (MaterialFront_extract(scratch,&mesh,&source,(uint32_t)i,(int32_t)nchart,&part,&npart)!=0 ||
            nfront>SIZE_MAX-npart) goto cube_done;
        cube[i].front=ARENA_ALLOC(arena,npart*sizeof *part); cube[i].nfront=npart;
        memcpy(cube[i].front,part,npart*sizeof *part);
        fprintf(ledger,"%s{\"id\":\"%s\",\"geometry_hash\":\"%016llx\",\"vertices\":%zu,\"faces\":%zu,"
                       "\"chart_base\":%zu,\"fronts\":%zu,\"patches\":[",
                i ? ",\n" : "",pile[i].cube_id,(unsigned long long)cube[i].hash,mesh.nv,mesh.nf,nchart,npart);
        double *area=ARENA_CALLOC(scratch,source.ncomponents,sizeof *area);
        size_t *counts=ARENA_CALLOC(scratch,source.ncomponents,sizeof *counts);
        size_t *head=ARENA_ALLOC(scratch,source.ncomponents*sizeof *head);
        size_t *tail=ARENA_ALLOC(scratch,source.ncomponents*sizeof *tail);
        MoRange *ranges=ARENA_ALLOC(scratch,mesh.nf*sizeof *ranges);
        size_t nranges=0;
        int32_t *seed=ARENA_ALLOC(scratch,source.ncomponents*sizeof *seed);
        for (size_t k=0;k<source.ncomponents;k++) { seed[k]=INT32_MAX; head[k]=tail[k]=SIZE_MAX; }
        for (size_t k=0;k<mesh.nv;k++) if ((int32_t)k<seed[source.component[k]]) seed[source.component[k]]=(int32_t)k;
        for (size_t k=0;k<mesh.nf;k++) {
            int32_t chart=source.component[mesh.faces[3*k]];
            double a[3]={0},b[3]={0},cross[3]={0};
            for (int d=0;d<3;d++) {
                a[d]=(double)mesh.verts[3*(size_t)mesh.faces[3*k+1]+d]-mesh.verts[3*(size_t)mesh.faces[3*k]+d];
                b[d]=(double)mesh.verts[3*(size_t)mesh.faces[3*k+2]+d]-mesh.verts[3*(size_t)mesh.faces[3*k]+d];
            }
            cross[0]=a[1]*b[2]-a[2]*b[1]; cross[1]=a[2]*b[0]-a[0]*b[2]; cross[2]=a[0]*b[1]-a[1]*b[0];
            area[chart]+=.5*sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]); counts[chart]++;
            if (tail[chart]!=SIZE_MAX && ranges[tail[chart]].end==k) ranges[tail[chart]].end=k+1;
            else {
                ranges[nranges]=(MoRange){k,k+1,SIZE_MAX};
                if (tail[chart]==SIZE_MAX) head[chart]=nranges;
                else ranges[tail[chart]].next=nranges;
                tail[chart]=nranges++;
            }
        }
        for (size_t k=0;k<source.ncomponents;k++) {
            fprintf(ledger,"%s{\"seed_vertex\":%d,\"source_state\":%u,\"faces\":%zu,\"area\":%.17g,\"face_ranges\":[",
                    k ? "," : "",seed[k],(unsigned)source.core[k],counts[k],area[k]);
            for (size_t r=head[k];r!=SIZE_MAX;r=ranges[r].next)
                fprintf(ledger,"%s[%zu,%zu]",r==head[k] ? "" : ",",ranges[r].begin,ranges[r].end);
            fprintf(ledger,"]}");
        }
        fprintf(ledger,"]}");
        snprintf(path,sizeof path,"%s/%s_source_component.i32",output_dir,pile[i].cube_id);
        table=fopen(path,"wb");
        if (!table) goto cube_done;
        ok=fwrite(source.component,sizeof(int32_t),mesh.nv,table)==mesh.nv;
        if (fclose(table)!=0) ok=0; table=NULL;
        nchart+=source.ncomponents; nfront+=npart;
        total_vertices+=mesh.nv; total_faces+=mesh.nf;
cube_done:
        Arena_dispose(&scratch);
        if (!ok) goto done;
        if (i%25==0 || i+1==npile)
            fprintf(stderr,"[material observation] source %zu/%zu: %zu charts, %zu boundary segments, %.2fs\n",
                    i+1,npile,nchart,nfront,ves_clock_sec()-started);
    }
    fprintf(ledger,"],\"vertices\":%llu,\"faces\":%llu}\n",
            (unsigned long long)total_vertices,(unsigned long long)total_faces);
    if (fclose(ledger)!=0) { ledger=NULL; goto done; } ledger=NULL;
    fronts=ARENA_ALLOC(arena,nfront*sizeof *fronts);
    for (size_t i=0,at=0;i<npile;i++) { memcpy(fronts+at,cube[i].front,cube[i].nfront*sizeof *fronts); at+=cube[i].nfront; }
    snprintf(path,sizeof path,"%s/material_ports.tsv",output_dir); table=fopen(path,"wb");
    if (!table) goto done;
    fprintf(table,"edge\tchart\tvertex_a\tvertex_b\tlength\tcurve\tarc\tregion_base\n");
    for (size_t i=0;i<nfront;i++) {
        const MaterialFront *f=fronts+i;
        double length=0;
        for (int k=0;k<3;k++) { double d=(double)f->xyz[k+3]-f->xyz[k]; length+=d*d; }
        fprintf(table,"%llu\t%d\t%d\t%d\t%.17g\t%llu\t%.17g\t%llu\n",(unsigned long long)f->source_edge,
                f->chart,f->source_vertices[0],f->source_vertices[1],sqrt(length),
                (unsigned long long)f->curve,f->arc,(unsigned long long)f->region_base);
    }
    if (fclose(table)!=0) { table=NULL; goto done; } table=NULL;
    if (MaterialFront_collect(arena,fronts,nfront,&intervals,&ninterval,&front_report)!=0 ||
        MaterialEvidence_reduce(arena,intervals,ninterval,MATERIAL_FRONT_LENGTH_SCALE,
                                &factors,&nfactor,&evidence_report)!=0 ||
        mo_phase_graph(arena,output_dir,nchart,factors,nfactor,&phase)!=0) goto done;
    snprintf(path,sizeof path,"%s/material_phase_hint.i64",output_dir); table=fopen(path,"wb");
    if (!table) goto done;
    if (fwrite(phase,sizeof(int64_t),nchart,table)!=nchart) goto done;
    if (fclose(table)!=0) { table=NULL; goto done; } table=NULL;
    snprintf(path,sizeof path,"%s/material_factors.tsv",output_dir); table=fopen(path,"wb");
    if (!table) goto done;
    fprintf(table,"a\tb\ttarget\tsolved\tkind\tunits\teffective_support\tcovered_length\tsignificant\tregions_a\tregions_b\tpair_margin\n");
    for (size_t i=0;i<nfactor;i++) {
        const MaterialEvidenceFactor *f=factors+i;
        fprintf(table,"%d\t%d\t%d\t%lld\t%d\t%zu\t%.17g\t%.17g\t%d\t%zu\t%zu\t%.17g\n",
                f->edge.a,f->edge.b,f->edge.target,(long long)(phase[f->edge.b]-phase[f->edge.a]),
                f->kind,f->units,f->effective_support,f->covered_length,f->significant,
                f->regions_a,f->regions_b,f->pair_margin);
    }
    if (fclose(table)!=0) { table=NULL; goto done; } table=NULL;
    snprintf(path,sizeof path,"%s/material_intervals.tsv",output_dir); table=fopen(path,"wb");
    if (!table) goto done;
    fprintf(table,"a\tb\tface_a\tface_b\tunit\ttarget\tbegin\tend\tquality\tedge_a\tedge_b\tregion_a\tregion_b\tbegin_b\tend_b\n");
    for (size_t i=0;i<ninterval;i++) {
        const MaterialEvidenceInterval *v=intervals+i;
        fprintf(table,"%d\t%d\t%llu\t%llu\t%llu\t%d\t%.17g\t%.17g\t%.17g\t%llu\t%llu\t%llu\t%llu\t%.17g\t%.17g\n",
                v->a,v->b,(unsigned long long)v->source_face_a,(unsigned long long)v->source_face_b,
                (unsigned long long)v->unit,v->target,v->begin,v->end,v->quality,
                (unsigned long long)v->source_edge_a,(unsigned long long)v->source_edge_b,
                (unsigned long long)v->region_a,(unsigned long long)v->region_b,v->begin_b,v->end_b);
    }
    if (fclose(table)!=0) { table=NULL; goto done; } table=NULL;
    snprintf(path,sizeof path,"%s/material_evidence_report.json",output_dir); table=fopen(path,"wb");
    if (!table) goto done;
    fprintf(table,"{\"schema\":\"material-evidence-report-v1\",\"source_cubes\":%zu,\"charts\":%zu,"
                  "\"boundary_segments\":%zu,\"candidate_pairs\":%zu,\"roundoff_only_overlaps\":%zu,\"intervals\":%zu,\"factors\":%zu,"
                  "\"significant_factors\":%zu,\"effective_support\":%.17g,\"elapsed_seconds\":%.6f,"
                  "\"correlation_length_voxels\":%.17g,\"support_units\":\"bounded source-boundary blocks on both sides; not independent trials\","
                  "\"equality_significance\":\"material support and local pair margin greater than twice one bounded unit; not global confidence\","
                  "\"order_evidence_integrated\":false,\"accepted_interfaces\":0,"
                  "\"physical_qualification\":\"INCOMPLETE: evidence and phase initialization only; no intrinsic UV or accepted interfaces\"}\n",
            npile,nchart,nfront,front_report.candidates,front_report.degenerate_intervals,
            ninterval,nfactor,evidence_report.significant_factors,
            evidence_report.effective_support,ves_clock_sec()-started,MATERIAL_FRONT_CORRELATION_LENGTH);
    if (fclose(table)!=0) { table=NULL; goto done; } table=NULL;
    fprintf(stderr,"[material observation] %zu intervals -> %zu soft factors (%zu significant), %.2fs; no physical qualification\n",
            ninterval,nfactor,evidence_report.significant_factors,ves_clock_sec()-started);
    rc=0;
done:
    if (rc!=0) fprintf(stderr,"[material observation] FAILED; retained files are incomplete diagnostics\n");
    if (ledger) fclose(ledger);
    if (table) fclose(table);
    return rc;
}
