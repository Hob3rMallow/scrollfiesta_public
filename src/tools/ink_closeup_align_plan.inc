/* Cuts are a recovery proposal, not a mandatory replacement of an already
 * useful uncut embedding. Both alternatives preserve the complete material. */
typedef struct IcAlignPlanStats {
    int uncut_status, uncut_valid, cut_attempted, cut_status, cut_valid, selected_cut;
    int reference_status, reference_remeasured;
    size_t proposed_cones, proposed_edges, proposed_fiber_cones;
    size_t cut_material_pairs;
    const char *reason;
} IcAlignPlanStats;

typedef struct IcAlignAttempt {
    double *uv;
    SlimRefineStats stats;
    IcPhysicalAlignStats physical;
    AsmContactStats contact;
    int status, valid;
} IcAlignAttempt;

static int ic_align_native_targets(const IcAlignAttempt *a)
{
    if (!a->valid || !(a->physical.area > 0)) return 0;
    if (a->physical.bad10_area_after > CLOSEUP_ALIGN_BAD_AREA_FRACTION*a->physical.area ||
        a->physical.protected_bad10_area_after > CLOSEUP_ALIGN_BAD_AREA_FRACTION*a->physical.protected_area) return 0;
    /* Missing physical evidence cannot be repaired by cutting the UV domain.
     * It still fails the independent fiber-support qualification gate. */
    if (!a->stats.direction_faces) return 1;
    double rms = sqrt(a->stats.direction_angle_after)*57.2957795130823208768;
    return isfinite(rms) && rms <= CLOSEUP_FIBER_RMS_DEG;
}

/* An embedding refusal does not make the original physical metric disappear.
 * Measure it without topology or boundary assumptions when the uncut solver
 * could not provide its ordinary before-stats. Orientation-independent sigma
 * and cross metrics are sufficient here; this is not an embedding certificate. */
static int ic_align_original_measure(const MeshBinData *mesh, const uint8_t *attributes,
                                      const SlimRefineDirection *directions,
                                      SlimRefineStats *stats, IcPhysicalAlignStats *physical)
{
    Arena_T arena = Arena_new();
    if (!arena) return -1;
    volatile int rc = -1;
    TRY {
        double *seed = ARENA_ALLOC(arena,2*mesh->nv*sizeof *seed);
        uint8_t *protected_faces = attributes ? ARENA_ALLOC(arena,mesh->nf) : NULL;
        CloseupAlignMetric metric = {0};
        for (size_t k = 0; k < 2*mesh->nv; k++) seed[k] = MeshBin_uv(mesh,k);
        if (protected_faces) for (size_t f = 0; f < mesh->nf; f++) {
            int original = 1;
            for (size_t k = 0; k < 3; k++) {
                int32_t v = mesh->faces[3*f+k];
                if (v < 0 || (size_t)v >= mesh->nv) { original = 0; break; }
                original &= attributes[v] == 1 || attributes[v] == 3;
            }
            protected_faces[f] = (uint8_t)original;
        }
        if (!CloseupAlign_measure(arena,mesh->verts,mesh->nv,mesh->faces,mesh->nf,
                                   seed,directions,protected_faces,1,&metric) &&
            !metric.invalid_evidence_faces && isfinite(metric.sd_mean) && metric.area > 0) {
            SlimRefineStats measured = {0};
            IcPhysicalAlignStats detail = {0};
            measured.energy_before = metric.area*metric.sd_mean;
            measured.direction_energy_before = CLOSEUP_FIBER_DIRECTION_WEIGHT*metric.direction_energy;
            measured.direction_faces = metric.measured_faces;
            measured.direction_area_fraction = metric.measured_area/metric.area;
            double angle = metric.fiber_rms_degrees/57.2957795130823208768;
            measured.direction_angle_before = angle*angle;
            measured.objective_before = measured.energy_before+
                (metric.confidence_mass>0?32*metric.confidence_mass*angle*angle:0);
            detail.area = metric.area; detail.protected_area = metric.protected_area;
            detail.bad10_area_before = metric.bad10_area;
            detail.bad25_area_before = metric.bad25_area;
            detail.protected_bad10_area_before = metric.protected_bad10_area;
            detail.protected_bad25_area_before = metric.protected_bad25_area;
            detail.protected_bad25_faces_before = metric.protected_bad25_faces;
            if (isfinite(measured.objective_before)) {
                *stats = measured; *physical = detail; rc = 0;
            }
        }
    } EXCEPT(Arena_Failed) {
        rc = -1;
    } FINALLY { Arena_dispose(&arena); } END_TRY;
    return rc;
}

static int ic_align_prefer_cut(const IcAlignAttempt *uncut, const IcAlignAttempt *cut,
                                const IcAlignAttempt *reference)
{
    if (!cut->valid || reference->status || !(reference->physical.area > 0) ||
        !isfinite(reference->stats.energy_before) || !isfinite(reference->stats.objective_before)) return 0;
    double tolerance = 1e-10*fmax(1,reference->physical.area);
    /* The unchanged original embedding is the reference for both routes. */
    if (cut->stats.objective_after > reference->stats.objective_before+1e-8 ||
        cut->physical.bad10_area_after > CLOSEUP_ALIGN_BAD_AREA_FRACTION*reference->physical.area+tolerance ||
        cut->physical.bad25_area_after > reference->physical.bad25_area_before+tolerance ||
        cut->physical.protected_bad10_area_after > CLOSEUP_ALIGN_BAD_AREA_FRACTION*reference->physical.protected_area+tolerance ||
        cut->physical.protected_bad25_area_after > reference->physical.protected_bad25_area_before+tolerance ||
        cut->physical.protected_bad25_faces_after > reference->physical.protected_bad25_faces_before) return 0;
    if (!uncut->valid) return 1;
    int uncut_met = ic_align_native_targets(uncut), cut_met = ic_align_native_targets(cut);
    if (uncut_met != cut_met) return cut_met;
    /* Both routes use the same physical budgets and joint angular objective. */
    if (cut->stats.objective_after > uncut->stats.objective_after+1e-8) return 0;
    if (!uncut->stats.direction_faces) return cut->stats.energy_after < uncut->stats.energy_after-1e-8;
    return isfinite(cut->stats.direction_angle_after) &&
           cut->stats.direction_angle_after < uncut->stats.direction_angle_after-1e-12;
}

static int ic_align_clone(const MeshBinData *source, const double *seed, MeshBinData *copy)
{
    copy->nv = source->nv; copy->nf = source->nf;
    copy->verts = malloc(3*source->nv*sizeof *copy->verts);
    copy->faces = malloc(3*source->nf*sizeof *copy->faces);
    copy->uv64 = malloc(2*source->nv*sizeof *copy->uv64);
    if (!copy->verts || !copy->faces || !copy->uv64) return -1;
    memcpy(copy->verts,source->verts,3*source->nv*sizeof *copy->verts);
    memcpy(copy->faces,source->faces,3*source->nf*sizeof *copy->faces);
    for (size_t k = 0; k < 2*source->nv; k++) copy->uv64[k] = seed ? seed[k] : MeshBin_uv(source,k);
    return 0;
}

static void ic_align_attempt(const IcConfig *config, const char *name, const MeshBinData *mesh,
                              const uint8_t *attributes, const SlimRefineDirection *directions,
                              const int32_t *aliases, const CloseupAlignOpts *options, IcAlignAttempt *attempt)
{
    char path[2600];
    attempt->status = -1;
    attempt->uv = malloc(2*mesh->nv*sizeof *attempt->uv);
    if (!attempt->uv) return;
    snprintf(path,sizeof path,"%s/%s_alignment_%s_components.jsonl",config->out_dir,config->id,name);
    ves_ensure_parent_dir(path);
    FILE *ledger = fopen(path,"wb");
    if (!ledger) return;
    attempt->status = ic_reparameterize_physical_material(mesh,directions,attributes,options,attempt->uv,
                                                         &attempt->stats,&attempt->physical,ledger,aliases,1);
    if (fclose(ledger)) attempt->status = -1;
    if (!attempt->status) {
        (void)ic_material_measure(mesh,attempt->uv,aliases,&attempt->contact);
        attempt->valid = ic_material_qualified(&attempt->contact);
    }
}

static int ic_reparameterize_physical_plan(const IcConfig *config, MeshBinData *mesh,
                                            const uint8_t *attributes, const SlimRefineDirection *directions,
                                            const int32_t *original_aliases,
                                            int allow_cuts, double **uv, int32_t **source_map,
                                            SlimRefineStats *stats, IcPhysicalAlignStats *physical,
                                            AsmContactStats *contact, size_t *cones, size_t *edges,
                                            size_t *fiber_cones, IcAlignPlanStats *plan)
{
    IcAlignAttempt uncut = {0}, cut = {0};
    MeshBinData cut_mesh = {0};
    int32_t *cut_map = NULL,*cut_aliases=NULL;
    uint8_t *cut_attributes = NULL;
    CloseupAlignOpts options;
    IcProjectionLinear projection = {0};
    IcAlignPlanStats result = {0};
    double started = ves_clock_sec();
    int rc = -1;
    result.cut_status = -1; result.reason = "no_feasible_proposal";
    CloseupAlign_defaults(&options);
    options.slim.iterations = CLOSEUP_FLAT_ITERS;
    options.slim.energy_tolerance = CLOSEUP_FLAT_TOL;
    ic_projection_options(&options,&projection);
    options.slim.verbose = 1;
    ic_align_attempt(config,"uncut",mesh,attributes,directions,original_aliases,&options,&uncut);
    IcAlignAttempt reference = uncut;
    result.reference_status = reference.status;
    result.uncut_status = uncut.status; result.uncut_valid = uncut.valid;
    if (!allow_cuts) result.reason = "uncut_command";
    else if (ic_align_native_targets(&uncut)) result.reason = "uncut_native_targets_met";
    else {
        result.cut_attempted = 1;
        if (ic_align_clone(mesh,uncut.valid?uncut.uv:NULL,&cut_mesh)) {
            result.reason = "cut_clone_allocation_failed";
        } else if (ic_cone_cut_directed(&cut_mesh,directions,&cut_map,&result.proposed_cones,
                                         &result.proposed_edges,&result.proposed_fiber_cones)) {
            result.reason = "cut_construction_refused";
        } else if (!result.proposed_edges) {
            result.cut_status = 0; result.reason = "no_additional_cone_cuts";
        } else {
            cut_attributes = malloc(cut_mesh.nv);
            if (!cut_attributes) result.reason = "cut_attributes_allocation_failed";
            else {
                for (size_t v = 0; v < cut_mesh.nv; v++)
                    cut_attributes[v] = attributes[cut_map?(size_t)cut_map[v]:v];
                if(ic_material_compose_aliases(mesh->nv,original_aliases,cut_mesh.nv,cut_map,&cut_aliases))goto done;
                ic_align_attempt(config,"cut",&cut_mesh,cut_attributes,directions,cut_aliases,&options,&cut);
                result.cut_status = cut.status; result.cut_valid = cut.valid;
                result.cut_material_pairs = cut.contact.pairs;
                result.reason = "cut_trial_did_not_improve";
                if (cut.valid && reference.status) {
                    result.reference_remeasured = 1;
                    reference.status = ic_align_original_measure(mesh,attributes,directions,
                                                                 &reference.stats,&reference.physical);
                    result.reference_status = reference.status;
                    if (reference.status) result.reason = "original_metric_unmeasured";
                }
                if (ic_align_prefer_cut(&uncut,&cut,&reference)) {
                    result.selected_cut = 1; result.reason = "cut_trial_improved";
                }
            }
        }
    }
    IcAlignAttempt *chosen = result.selected_cut ? &cut : &uncut;
    if (!chosen->valid) goto done;
    *stats = chosen->stats; *physical = chosen->physical; *contact = chosen->contact;
    if (result.selected_cut) {
        stats->energy_before = reference.stats.energy_before;
        stats->objective_before = reference.stats.objective_before;
        stats->direction_energy_before = reference.stats.direction_energy_before;
        stats->direction_angle_before = reference.stats.direction_angle_before;
        physical->bad10_area_before = reference.physical.bad10_area_before;
        physical->bad25_area_before = reference.physical.bad25_area_before;
        physical->protected_bad10_area_before = reference.physical.protected_bad10_area_before;
        physical->protected_bad25_area_before = reference.physical.protected_bad25_area_before;
        physical->protected_bad25_faces_before = reference.physical.protected_bad25_faces_before;
    }
    stats->iterations_run = uncut.stats.iterations_run+cut.stats.iterations_run;
    stats->accepted_steps = uncut.stats.accepted_steps+cut.stats.accepted_steps;
    stats->total_seconds = ves_clock_sec()-started;
    *uv = chosen->uv; chosen->uv = NULL;
    if (result.selected_cut) {
        MeshBin_dispose(mesh); *mesh = cut_mesh; memset(&cut_mesh,0,sizeof cut_mesh);
        *source_map = cut_map; cut_map = NULL;
        *cones = result.proposed_cones; *edges = result.proposed_edges; *fiber_cones = result.proposed_fiber_cones;
    }
    rc = 0;
done:
    *plan = result;
    fprintf(stderr,"[ink_closeup] physical cut plan: %s; uncut-valid=%d cut-valid=%d selected-cut=%d\n",
            result.reason,result.uncut_valid,result.cut_valid,result.selected_cut);
    free(uncut.uv); free(cut.uv); free(cut_map); free(cut_aliases); free(cut_attributes); MeshBin_dispose(&cut_mesh);
    return rc;
}

static void ic_alignment_plan_report(FILE *fp, const IcAlignPlanStats *s)
{
    fprintf(fp," \"cut_planning\": {\"policy\":\"uncut_control_then_required_cone_proposal_v1\","
            "\"reason\":\"%s\",\"uncut_status\":%d,\"uncut_valid\":%s,\"cut_attempted\":%s,"
            "\"reference_status\":%d,\"reference_remeasured\":%s,"
            "\"cut_status\":%d,\"cut_valid\":%s,\"selected_cut\":%s,"
            "\"proposed_cones\":%zu,\"proposed_edges\":%zu,\"proposed_fiber_cones\":%zu,"
            "\"cut_material_pairs\":%zu},\n",s->reason,s->uncut_status,s->uncut_valid?"true":"false",
            s->cut_attempted?"true":"false",s->reference_status,s->reference_remeasured?"true":"false",
            s->cut_status,s->cut_valid?"true":"false",s->selected_cut?"true":"false",
            s->proposed_cones,s->proposed_edges,s->proposed_fiber_cones,s->cut_material_pairs);
}
