/* Repair the original cached seed before any source registration. This is
 * private proposal work: qualification here admits no chart and weakens no
 * source, seam or contact obligation in the subsequent transaction. */
static int api_recover(AsmField *seed, AsmChart *chart, Arena_T arena,
                       const AsmRepairOpts *opts, const char **reason, FILE *trace)
{
    double start = ves_clock_sec(); int okay = 0; AmrResult result = {0};
    AmfStats source = {0}; const char *stop = "intrinsic_seed_unavailable";
    if (!seed) goto done;
    result.before = amr_metric(seed,0,seed->uv); result.after = result.before;
    if (apr_expired(opts)) { stop = "intrinsic_recovery_deadline"; goto done; }
    source = amf_chart(chart,NULL);
    if (!source.io) { stop = "intrinsic_source_invalid"; goto done; }
    if (source.orientation_obstructions) { stop = "intrinsic_source_needs_cut"; goto done; }
    double deadline = opts->deadline > 0 ? opts->deadline : start+opts->patch_seconds;
    if (opts->admission_solver_deadline > 0) {
#pragma omp critical(asm_region_deformation)
        {
            okay = amb_refine_region_budget(seed,0,NULL,NULL,
                fmin(deadline,opts->admission_solver_deadline),1,&result,opts,1);
        }
    } else okay = amb_refine_region_budget(seed,0,NULL,NULL,deadline,1,&result,opts,1);
    result.after = amr_metric(seed,0,seed->uv);
    if (!okay) { stop = result.reason; goto done; }
    /* The stored intrinsic format is float. Certify those exact rounded
     * coordinates, in the field's mirrored frame, before keeping the seed. */
    float *uv = ARENA_ALLOC(arena,2*chart->nv*sizeof(float));
    double sign = chart->flags & ASM_CHART_MIRROR ? -1 : 1;
    for (size_t v = 0; v < chart->nv; v++) {
        uv[2*v] = (float)(sign*seed->uv[2*v]); uv[2*v+1] = (float)seed->uv[2*v+1];
        seed->uv[2*v] = sign*(double)uv[2*v]; seed->uv[2*v+1] = uv[2*v+1];
    }
    result.after = amr_metric(seed,0,seed->uv); okay = amr_qualified(&result.after);
    if (!okay) { stop = "intrinsic_float_metric_failed"; goto done; }
    if (apr_expired(opts) || !apr_workspace_check(opts,0)) {
        okay = 0; stop = opts->workspace && opts->workspace->exceeded ? "workspace_memory_budget" : "intrinsic_recovery_deadline"; goto done;
    }
    chart->uv = uv; chart->placed_uv = NULL; chart->pose_theta = chart->pose_x = chart->pose_y = 0;
    stop = "intrinsic_metric_qualified";
done:
    if (!okay) *reason = stop;
    if (trace) fprintf(trace,"{\"phase\":\"intrinsic_recovery\",\"chart\":%d,\"reason\":\"%s\",\"qualified\":%s,"
        "\"source_obstructions\":%zu,\"single_turn_obstructions\":%zu,\"obstruction_area\":%.17g,\"before_bad\":%zu,\"after_bad\":%zu,\"after_invalid\":%zu,"
        "\"before_within10\":%.17g,\"after_within10\":%.17g,\"after_lo\":%.17g,\"after_hi\":%.17g,"
        "\"free_vertices\":%zu,\"active_faces\":%zu,\"steps\":%zu,\"factor_seconds\":%.17g,\"seconds\":%.17g}\n",
        chart->id,stop,okay?"true":"false",source.orientation_obstructions,source.obstruction_stars,source.affected_area,result.before.bad,result.after.bad,result.after.invalid,
        result.before.area>0?result.before.within10/result.before.area:0,result.after.area>0?result.after.within10/result.after.area:0,
        result.after.lo,result.after.hi,result.free_vertices,result.active_faces,result.refine_steps,result.factor_seconds,ves_clock_sec()-start);
    return okay;
}
