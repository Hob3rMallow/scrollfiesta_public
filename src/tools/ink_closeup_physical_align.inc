/* Complete-component physical alignment. All original vertices, including
 * isolated samples, survive in their original order. Component extraction is
 * temporary; no chart packing, scale change or boundary pin is introduced. */
#include "../flatten/closeup_align.h"
#ifdef VES_CLOSEUP_LOCAL_GLOBAL
#include "../flatten/closeup_local_global.h"
#endif
#include "ink_closeup_global_contact.inc"

typedef struct IcPhysicalAlignStats {
    size_t components, proposed, unchanged, unknown, orphans, trials;
    size_t reference_objective_refusals;
    size_t protected_bad25_faces_before,protected_bad25_faces_after;
    size_t rigid_proposed, rigid_trials;
    double rigid_seconds;
    size_t global_contact_calls, global_contact_refusals, global_contact_unknown;
    double area, measured_area, confidence_mass;
    double bad10_area_before, bad10_area_after;
    double bad25_area_before, bad25_area_after;
    double protected_area, protected_bad10_area_before, protected_bad10_area_after;
    double protected_bad25_area_before, protected_bad25_area_after;
    double min_sigma, max_sigma;
} IcPhysicalAlignStats;

typedef struct IcAlignPart {
    size_t vertices, faces, vertex_begin, face_begin;
} IcAlignPart;

static int32_t ic_align_root(int32_t *parent, int32_t vertex)
{
    while (parent[vertex] != vertex) {
        parent[vertex] = parent[parent[vertex]];
        vertex = parent[vertex];
    }
    return vertex;
}

static void ic_align_json_number(char text[64], double value)
{
    if (isfinite(value)) snprintf(text,64,"%.17g",value);
    else memcpy(text,"null",5);
}

/* A local allowance must reserve the bad area on every stationary component.
 * In particular, components already above 5% cannot lend imaginary capacity. */
static double ic_align_bad10_limit(double sheet_area,double sheet_bad,
                                   double local_area,double local_bad)
{
    return fmax(0,fmin(.05*local_area,.05*sheet_area-(sheet_bad-local_bad)));
}

static int ic_align_components(Arena_T arena, const MeshBinData *mesh,
                                const SlimRefineDirection *directions,
                                const uint8_t *attributes,
                                const CloseupAlignOpts *options, double *output,
                                SlimRefineStats *stats, IcPhysicalAlignStats *detail,
                                FILE *ledger, const int32_t *aliases,
                                int global_guard, Arena_T contact_arena, int rigid_only)
{
    size_t nv = mesh->nv, nf = mesh->nf, nc = 0, incident = 0;
    int32_t *parent = ARENA_ALLOC(arena,nv*sizeof *parent);
    int32_t *map = ARENA_ALLOC(arena,nv*sizeof *map);
    double total_angle_before = 0, total_angle_after = 0;
    for (size_t v = 0; v < nv; v++) {
        parent[v] = (int32_t)v; map[v] = -1;
        for (size_t k = 0; k < 3; k++) if (!isfinite(mesh->verts[3*v+k])) return -1;
        for (size_t k = 0; k < 2; k++) {
            double value = MeshBin_uv(mesh,2*v+k);
            if (!isfinite(value)) return -1;
            output[2*v+k] = value;
        }
    }
    for (size_t f = 0; f < nf; f++) {
        const int32_t *face = mesh->faces+3*f;
        for (size_t k = 0; k < 3; k++) {
            if (face[k] < 0 || (size_t)face[k] >= nv) return -1;
            map[face[k]] = 0;
        }
        for (size_t k = 1; k < 3; k++) {
            int32_t a = ic_align_root(parent,face[0]), b = ic_align_root(parent,face[k]);
            if (a != b) parent[a > b ? a : b] = a < b ? a : b;
        }
    }
    for (size_t v = 0; v < nv; v++) if (map[v] >= 0) parent[v] = ic_align_root(parent,(int32_t)v);
    /* Every parent is now a root. Negative root entries encode component ids;
     * this reuses the union array rather than allocating a second nv lookup. */
    for (size_t v = 0; v < nv; v++) if (map[v] >= 0 && parent[v] == (int32_t)v)
        parent[v] = -(int32_t)(++nc);
    for (size_t v = 0; v < nv; v++) if (map[v] >= 0 && parent[v] >= 0) parent[v] = parent[parent[v]];
    if (!nc) return -1;
    IcAlignPart *part = ARENA_CALLOC(arena,nc,sizeof *part);
    for (size_t v = 0; v < nv; v++) if (map[v] >= 0) {
        part[(size_t)(-parent[v]-1)].vertices++; incident++;
    }
    for (size_t f = 0; f < nf; f++) part[(size_t)(-parent[mesh->faces[3*f]]-1)].faces++;
    size_t vertex_at = 0, face_at = 0;
    for (size_t c = 0; c < nc; c++) {
        part[c].vertex_begin = vertex_at; part[c].face_begin = face_at;
        vertex_at += part[c].vertices; face_at += part[c].faces;
        part[c].vertices = part[c].faces = 0;
    }
    int32_t *vertices = ARENA_ALLOC(arena,incident*sizeof *vertices);
    int32_t *faces = ARENA_ALLOC(arena,nf*sizeof *faces);
    for (size_t v = 0; v < nv; v++) if (map[v] >= 0) {
        IcAlignPart *p = part+(size_t)(-parent[v]-1);
        vertices[p->vertex_begin+p->vertices++] = (int32_t)v;
    }
    for (size_t f = 0; f < nf; f++) {
        IcAlignPart *p = part+(size_t)(-parent[mesh->faces[3*f]]-1);
        faces[p->face_begin+p->faces++] = (int32_t)f;
    }
    stats->components = detail->components = nc;
    detail->orphans = nv-incident;
    detail->min_sigma = stats->min_det_before = stats->min_det_after = DBL_MAX;
    IcAlignGlobalContact *global = global_guard ? ic_align_global_new(contact_arena,mesh,output,aliases) : NULL;
    if (global_guard && !global) return -1;
    if (global) {
        global->extra_context = options->candidate_context;
        global->extra_guard = options->candidate_guard;
    }
    double sheet_mass=0,sheet_angle=0,sheet_area[2]={0},sheet_bad[2]={0};
#ifdef VES_CLOSEUP_LOCAL_GLOBAL
    if(directions&&!rigid_only) {
        CloseupAlignMetric metric={0};
        uint8_t *protected_all=ARENA_ALLOC(arena,nf);
        for(size_t f=0;f<nf;f++) {
            int original=1;
            if(attributes)for(size_t k=0;k<3;k++) {
                uint8_t value=attributes[mesh->faces[3*f+k]];original&=value==1||value==3;
            }
            protected_all[f]=(uint8_t)original;
        }
        if(CloseupAlign_measure(arena,mesh->verts,nv,mesh->faces,nf,output,directions,protected_all,1,&metric))return -1;
        sheet_mass=metric.confidence_mass;
        sheet_area[0]=metric.area;sheet_area[1]=metric.protected_area;
        sheet_bad[0]=metric.bad10_area;sheet_bad[1]=metric.protected_bad10_area;
        if(sheet_mass>0) {
            double radians=metric.fiber_rms_degrees/57.2957795130823208768;
            sheet_angle=sheet_mass*radians*radians;
        }
    }
#endif
    for (size_t c = 0; c < nc; c++) {
        Arena_Mark mark = Arena_save(arena);
        const IcAlignPart *p = part+c;
        float *xyz = ARENA_ALLOC(arena,3*p->vertices*sizeof *xyz);
        double *seed = ARENA_ALLOC(arena,2*p->vertices*sizeof *seed);
        double *candidate = ARENA_ALLOC(arena,2*p->vertices*sizeof *candidate);
        int32_t *triangles = ARENA_ALLOC(arena,3*p->faces*sizeof *triangles);
        SlimRefineDirection *arms = directions ? ARENA_ALLOC(arena,p->faces*sizeof *arms) : NULL;
        uint8_t *protected_faces = attributes ? ARENA_ALLOC(arena,p->faces) : NULL;
        CloseupAlignStats report = {0};
        for (size_t k = 0; k < p->vertices; k++) {
            size_t v = (size_t)vertices[p->vertex_begin+k];
            map[v] = (int32_t)k;
            memcpy(xyz+3*k,mesh->verts+3*v,3*sizeof *xyz);
            seed[2*k] = output[2*v]; seed[2*k+1] = output[2*v+1];
        }
        for (size_t k = 0; k < p->faces; k++) {
            size_t f = (size_t)faces[p->face_begin+k];
            for (size_t a = 0; a < 3; a++) triangles[3*k+a] = map[mesh->faces[3*f+a]];
            if (arms) memcpy(arms+k,directions+f,sizeof *arms);
            if (protected_faces) {
                int original = 1;
                for (size_t a = 0; a < 3; a++) {
                    uint8_t value = attributes[mesh->faces[3*f+a]];
                    original &= value == 1 || value == 3;
                }
                protected_faces[k] = (uint8_t)original;
            }
        }
        CloseupAlignOpts local_options = *options;
        double local_limits[2]={0};
#ifdef VES_CLOSEUP_LOCAL_GLOBAL
        if(!rigid_only&&sheet_mass>0) {
            CloseupAlignMetric local={0};
            if(CloseupAlign_measure(arena,xyz,p->vertices,triangles,p->faces,seed,arms,protected_faces,1,&local))return -1;
            local_limits[0]=ic_align_bad10_limit(sheet_area[0],sheet_bad[0],local.area,local.bad10_area);
            local_limits[1]=ic_align_bad10_limit(sheet_area[1],sheet_bad[1],local.protected_area,local.protected_bad10_area);
            local_options.local_global_bad10_limits=local_limits;
            double target=options->target_rms_degrees/57.2957795130823208768;
            double budget=sheet_mass*target*target;
            local_options.local_global_iterations=sheet_angle>=budget?1000:0;
            if(local.confidence_mass>0) {
                double radians=local.fiber_rms_degrees/57.2957795130823208768;
                double remaining=budget-(sheet_angle-local.confidence_mass*radians*radians);
                /* First spend nonlinear work on components exceeding the
                 * requested RMS. A component already below it must not be
                 * driven far lower merely to subsidize an unvisited one. */
                if(remaining>0)local_options.target_rms_degrees=fmax(options->target_rms_degrees,fmin(45,
                    sqrt(remaining/local.confidence_mass)*57.2957795130823208768)*(1-1e-12));
            }
        }
#endif
        size_t guard_calls = global ? global->calls : 0, guard_refused = global ? global->refused : 0;
        if (global) {
            global->vertices = vertices+p->vertex_begin; global->local_vertices = p->vertices;
            memset(global->moving,0,nf);
            for (size_t k = 0; k < p->faces; k++) global->moving[faces[p->face_begin+k]] = 1;
            AsmContacts_active_faces(global->query,global->moving);
            local_options.candidate_context = global;
            local_options.candidate_guard = ic_align_global_candidate;
        }
        int (*align)(Arena_T,const float *,size_t,const int32_t *,size_t,const double *,
                     const SlimRefineDirection *,const uint8_t *,const CloseupAlignOpts *,
                     double *,CloseupAlignStats *) = rigid_only ? CloseupAlign_rotate : CloseupAlign_propose;
        int status = align(arena,xyz,p->vertices,triangles,p->faces,
                                          seed,arms,protected_faces,&local_options,candidate,&report);
        if (status < 0 || report.before.invalid_evidence_faces || report.after.invalid_evidence_faces) {
            fprintf(stderr,"ERROR: physical alignment component %zu (first face %d) refused: %d, invalid evidence %zu\n",
                    c,faces[p->face_begin],status,report.before.invalid_evidence_faces+report.after.invalid_evidence_faces);
            Arena_restore(arena,mark); return -1;
        }
        const double weight = CLOSEUP_FIBER_DIRECTION_WEIGHT;
        double radians_before=report.before.fiber_rms_degrees/57.2957795130823208768;
        double radians_after=report.after.fiber_rms_degrees/57.2957795130823208768;
        double before = report.before.sd_mean*report.before.area +
            (report.before.confidence_mass>0?32*report.before.confidence_mass*radians_before*radians_before:0);
        double after = report.after.sd_mean*report.after.area +
            (report.after.confidence_mass>0?32*report.after.confidence_mass*radians_after*radians_after:0);
        int reference_ok = isfinite(before) && isfinite(after) && after <= before+1e-8;
        /* Adaptive solve weights and the fixed audit objective are distinct.
         * An otherwise useful proposal may fail this existing source contract.
         * Retain its exact original UV, with the refusal visible in the ledger. */
        if (!reference_ok) {
            detail->reference_objective_refusals++;
            memcpy(candidate,seed,2*p->vertices*sizeof *candidate);
            report.after = report.before;
            status = CLOSEUP_ALIGN_UNCHANGED;
        }
        if(sheet_mass>0&&report.before.confidence_mass>0) {
            double a=report.before.fiber_rms_degrees/57.2957795130823208768;
            double b=report.after.fiber_rms_degrees/57.2957795130823208768;
            sheet_angle+=report.before.confidence_mass*(b*b-a*a);
        }
        sheet_bad[0]+=report.after.bad10_area-report.before.bad10_area;
        sheet_bad[1]+=report.after.protected_bad10_area-report.before.protected_bad10_area;
        for (size_t k = 0; k < p->vertices; k++) {
            size_t v = (size_t)vertices[p->vertex_begin+k];
            output[2*v] = candidate[2*k]; output[2*v+1] = candidate[2*k+1];
        }
        if (global) {
            if (ic_align_global_commit(global)) {
                fprintf(stderr,"ERROR: complete-material incumbent refit failed after component %zu\n",c);
                Arena_restore(arena,mark); return -1;
            }
            detail->global_contact_calls = global->calls;
            detail->global_contact_refusals = global->refused;
            detail->global_contact_unknown = global->unknown;
            if (options->slim.verbose)
                fprintf(stderr,"[global component guard] component=%zu calls=%zu refused=%zu total-unknown=%zu\n",
                        c,global->calls-guard_calls,global->refused-guard_refused,global->unknown);
        }
        detail->proposed += status == CLOSEUP_ALIGN_PROPOSED;
        detail->unchanged += status == CLOSEUP_ALIGN_UNCHANGED;
        detail->unknown += status == CLOSEUP_ALIGN_NO_EVIDENCE;
        detail->trials += report.trials;
        for (size_t k = 0; k < report.trials; k++) {
            stats->iterations_run += report.trial[k].iterations_run;
            stats->accepted_steps += report.trial[k].accepted_steps;
        }
        detail->area += report.after.area;
        detail->measured_area += report.after.measured_area;
        detail->confidence_mass += report.after.confidence_mass;
        detail->bad10_area_before += report.before.bad10_area;
        detail->bad10_area_after += report.after.bad10_area;
        detail->bad25_area_before += report.before.bad25_area;
        detail->bad25_area_after += report.after.bad25_area;
        detail->protected_area += report.after.protected_area;
        detail->protected_bad10_area_before += report.before.protected_bad10_area;
        detail->protected_bad10_area_after += report.after.protected_bad10_area;
        detail->protected_bad25_area_before += report.before.protected_bad25_area;
        detail->protected_bad25_area_after += report.after.protected_bad25_area;
        detail->protected_bad25_faces_before += report.before.protected_bad25_faces;
        detail->protected_bad25_faces_after += report.after.protected_bad25_faces;
        detail->min_sigma = fmin(detail->min_sigma,report.after.min_sigma);
        detail->max_sigma = fmax(detail->max_sigma,report.after.max_sigma);
        stats->energy_before += report.before.sd_mean*report.before.area;
        stats->energy_after += report.after.sd_mean*report.after.area;
        stats->direction_energy_before += weight*report.before.direction_energy;
        stats->direction_energy_after += weight*report.after.direction_energy;
        stats->direction_faces += report.after.measured_faces;
        stats->min_det_before = fmin(stats->min_det_before,report.before.min_det);
        stats->min_det_after = fmin(stats->min_det_after,report.after.min_det);
        stats->flips_before += report.before.inverted_faces;
        stats->flips_after += report.after.inverted_faces;
        if (report.after.confidence_mass > 0) {
            double a = report.before.fiber_rms_degrees/57.2957795130823208768;
            double b = report.after.fiber_rms_degrees/57.2957795130823208768;
            total_angle_before += report.before.confidence_mass*a*a;
            total_angle_after += report.after.confidence_mass*b*b;
        }
        if (ledger) {
            if(report.local_global_attempted&&fprintf(ledger,
                "{\"phase\":\"local_global\",\"component\":%zu,\"status\":%d,\"termination\":%d,"
                "\"loss_converged\":%s,\"initial_loss\":%.17g,\"final_loss\":%.17g,\"selected\":%s}\n",
                c,report.local_global_status,report.local_global_termination,
                report.local_global_loss_converged?"true":"false",report.local_global_initial_loss,
                report.local_global_final_loss,report.selected_seed==3?"true":"false")<0)return -1;
            char rms_before[64], rms_after[64];
            ic_align_json_number(rms_before,report.before.fiber_rms_degrees);
            ic_align_json_number(rms_after,report.after.fiber_rms_degrees);
            if (fprintf(ledger,"{\"phase\":\"%s\",\"component\":%zu,\"first_source_face\":%d,\"vertices\":%zu,\"faces\":%zu,"
                        "\"status\":%d,\"trials\":%zu,\"selected_seed\":%d,\"selected_direction_weight\":%.17g,"
                        "\"reference_objective_accepted\":%s,\"rms_before\":%s,\"rms_after\":%s,"
                        "\"global_guard_calls\":%zu,\"global_guard_refusals\":%zu}\n",
                        rigid_only?"rigid":"slim",c,faces[p->face_begin],p->vertices,p->faces,status,report.trials,report.selected_seed,
                        report.selected_direction_weight,reference_ok?"true":"false",rms_before,rms_after,
                        global?global->calls-guard_calls:0,global?global->refused-guard_refused:0) < 0) {
                Arena_restore(arena,mark); return -1;
            }
        }
        if (options->slim.verbose) fprintf(stderr,"[ink_closeup] physical alignment component %zu/%zu: %zu faces, status %d, RMS %.8g -> %.8g\n",
                c+1,nc,p->faces,status,report.before.fiber_rms_degrees,report.after.fiber_rms_degrees);
        Arena_restore(arena,mark);
    }
    stats->objective_before = stats->energy_before+32*total_angle_before;
    stats->objective_after = stats->energy_after+32*total_angle_after;
    stats->direction_area_fraction = detail->area > 0 ? detail->measured_area/detail->area : 0;
    stats->direction_angle_before = detail->confidence_mass > 0 ? total_angle_before/detail->confidence_mass : NAN;
    stats->direction_angle_after = detail->confidence_mass > 0 ? total_angle_after/detail->confidence_mass : NAN;
    if (!isfinite(stats->objective_before) || !isfinite(stats->objective_after) ||
        stats->objective_after > stats->objective_before+1e-8 || stats->flips_after) return -1;
    return 0;
}

/* UV and summary outputs commit together after every component succeeds.
 * The optional JSONL ledger is diagnostic and may record an incomplete attempt. */
static int ic_reparameterize_physical_material(const MeshBinData *mesh,
                                      const SlimRefineDirection *directions,
                                      const uint8_t *attributes,
                                      const CloseupAlignOpts *options, double *uv,
                                      SlimRefineStats *stats, IcPhysicalAlignStats *detail,
                                      FILE *ledger, const int32_t *aliases, int global_guard)
{
    if (!mesh || !options || !uv || !stats || !detail || !mesh->verts || !mesh->faces ||
        (!mesh->uv && !mesh->uv64) || mesh->nv < 3 || !mesh->nf ||
        mesh->nv > (size_t)INT_MAX/2 || mesh->nf > (size_t)INT_MAX/3 ||
        mesh->nv > SIZE_MAX/(8*sizeof(double)) || mesh->nf > SIZE_MAX/sizeof(SlimRefineDirection)) return -1;
    Arena_T volatile arena = NULL, contact_arena = NULL;
    volatile int rc = -1;
    double started = ves_clock_sec();
    TRY {
        arena = Arena_new();
        if (global_guard) contact_arena = Arena_new();
        SlimRefineStats measured = {0};
        IcPhysicalAlignStats information = {0};
        double *candidate = ARENA_ALLOC(arena,2*mesh->nv*sizeof *candidate);
        Arena_Mark work = Arena_save(arena);
        rc = ic_align_components(arena,mesh,directions,attributes,options,candidate,&measured,&information,
                                 ledger,aliases,global_guard,contact_arena,0);
        Arena_restore(arena,work);
        if (!rc && measured.direction_faces) {
            /* First retain every SLIM decision in the original packed context.
             * Rotations then preserve that entire result's intrinsic metric. */
            double rotate_started = ves_clock_sec();
            if (contact_arena) Arena_free(contact_arena);
            MeshBinData flattened = *mesh;
            flattened.uv64 = candidate; flattened.uv = NULL;
            double *rotated = ARENA_ALLOC(arena,2*mesh->nv*sizeof *rotated);
            SlimRefineStats after = {0}; IcPhysicalAlignStats turns = {0};
            rc = ic_align_components(arena,&flattened,directions,attributes,options,rotated,&after,&turns,
                                     ledger,aliases,global_guard,contact_arena,1);
            if (!rc) {
                measured.energy_after = after.energy_after;
                measured.objective_after = after.objective_after;
                measured.direction_energy_after = after.direction_energy_after;
                measured.direction_angle_after = after.direction_angle_after;
                measured.min_det_after = after.min_det_after;
                measured.flips_after = after.flips_after;
                measured.boundary_intersections_after = after.boundary_intersections_after;
                information.bad10_area_after = turns.bad10_area_after;
                information.bad25_area_after = turns.bad25_area_after;
                information.protected_bad10_area_after = turns.protected_bad10_area_after;
                information.protected_bad25_area_after = turns.protected_bad25_area_after;
                information.protected_bad25_faces_after = turns.protected_bad25_faces_after;
                information.min_sigma = turns.min_sigma; information.max_sigma = turns.max_sigma;
                information.global_contact_calls += turns.global_contact_calls;
                information.global_contact_refusals += turns.global_contact_refusals;
                information.global_contact_unknown += turns.global_contact_unknown;
                information.reference_objective_refusals += turns.reference_objective_refusals;
                information.rigid_proposed = turns.proposed; information.rigid_trials = turns.trials;
                information.rigid_seconds = ves_clock_sec()-rotate_started;
                candidate = rotated;
            }
        }
        if (!rc) {
            measured.total_seconds = ves_clock_sec()-started;
            memcpy(uv,candidate,2*mesh->nv*sizeof *uv);
            *stats = measured; *detail = information;
        }
    } EXCEPT(Arena_Failed) {
        fprintf(stderr,"ERROR: physical alignment allocation failed; original surface retained\n");
        rc = -1;
    } FINALLY {
        Arena_T owned_contacts = contact_arena, owned_scratch = arena;
        if (owned_contacts) Arena_dispose(&owned_contacts);
        if (owned_scratch) Arena_dispose(&owned_scratch);
    } END_TRY;
    return rc;
}

/* Identity-topology callers also use the whole-material guard. Only native
 * regression probes invoke the lower-level function without the guard. */
static int ic_reparameterize_physical(const MeshBinData *mesh,
                                      const SlimRefineDirection *directions,
                                      const uint8_t *attributes,
                                      const CloseupAlignOpts *options, double *uv,
                                      SlimRefineStats *stats, IcPhysicalAlignStats *detail,
                                      FILE *ledger)
{
    return ic_reparameterize_physical_material(mesh,directions,attributes,options,uv,stats,detail,ledger,NULL,1);
}

static void ic_physical_alignment_report(FILE *fp, const IcPhysicalAlignStats *s)
{
    fprintf(fp," \"physical_alignment\": {\"schema\":\"closeup-physical-align-v1\","
            "\"component_scope\":\"all_incident_faces\",\"protected_faces\":\"all_three_DATA_or_adopted_vertices\","
            "\"boundary_free\":true,\"component_centroids_preserved\":true,"
            "\"components\":%zu,\"proposed\":%zu,\"unchanged\":%zu,\"unknown\":%zu,"
            "\"rigid_proposed\":%zu,\"rigid_trials\":%zu,\"rigid_seconds\":%.9g,"
            "\"orphan_vertices_preserved\":%zu,\"solver_trials\":%zu,\"reference_objective_refusals\":%zu,"
            "\"global_guard_calls\":%zu,\"global_guard_refusals\":%zu,\"global_guard_unknown\":%zu,"
            "\"area\":%.17g,\"measured_area\":%.17g,\"confidence_mass\":%.17g,"
            "\"bad10_area_before\":%.17g,\"bad10_area_after\":%.17g,"
            "\"bad25_area_before\":%.17g,\"bad25_area_after\":%.17g,"
            "\"protected_area\":%.17g,\"protected_bad10_area_before\":%.17g,\"protected_bad10_area_after\":%.17g,"
            "\"protected_bad25_area_before\":%.17g,\"protected_bad25_area_after\":%.17g,"
            "\"min_sigma\":%.17g,\"max_sigma\":%.17g},\n",
            s->components,s->proposed,s->unchanged,s->unknown,
            s->rigid_proposed,s->rigid_trials,s->rigid_seconds,s->orphans,s->trials,
            s->reference_objective_refusals,
            s->global_contact_calls,s->global_contact_refusals,s->global_contact_unknown,
            s->area,s->measured_area,s->confidence_mass,
            s->bad10_area_before,s->bad10_area_after,s->bad25_area_before,s->bad25_area_after,
            s->protected_area,s->protected_bad10_area_before,s->protected_bad10_area_after,
            s->protected_bad25_area_before,s->protected_bad25_area_after,
            s->min_sigma,s->max_sigma);
}
