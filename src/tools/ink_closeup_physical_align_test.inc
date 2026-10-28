static int ic_align_check(int condition, const char *name)
{
    fprintf(stderr,"[closeup physical integration] %s %s\n",condition?"PASS":"FAIL",name);
    return !condition;
}

#include "ink_closeup_local_global_test.inc"

static int ic_align_original_reference_selftest(void)
{
    float xyz[12] = {0,0,0, 1,0,0, 1,1,0, 0,1,0};
    int32_t faces[6] = {0,1,2, 0,2,3};
    uint8_t attributes[4] = {1,3,1,2};
    double seed[8], good[8] = {0,0, 1,0, 1,1, 0,1};
    double worse[8] = {0,0, 1.4,0, 1.4,1.4, 0,1.4};
    SlimRefineDirection arms[2] = {0};
    MeshBinData mesh = {0};
    IcAlignAttempt reference = {0}, uncut = {0}, cut = {0};
    SlimRefineStats saved_stats;
    IcPhysicalAlignStats saved_physical;
    Arena_T arena = Arena_new();
    CloseupAlignMetric metric = {0};
    uint8_t protected_faces[2] = {1,0};
    int failures = 0;
    mesh.nv = 4; mesh.nf = 2; mesh.verts = xyz; mesh.faces = faces; mesh.uv64 = seed;
    for (size_t f = 0; f < 2; f++) {
        arms[f].h[0] = 1; arms[f].v[1] = 1; arms[f].weight = 1;
    }
    for (size_t v = 0; v < 4; v++) {
        double x = 1.06*good[2*v], y = .97*good[2*v+1];
        seed[2*v] = cos(.3)*x-sin(.3)*y;
        seed[2*v+1] = sin(.3)*x+cos(.3)*y;
    }
    reference.status = ic_align_original_measure(&mesh,attributes,arms,&reference.stats,&reference.physical);
    double expected_sd = 1.06*1.06+.97*.97+1/(1.06*1.06)+1/(.97*.97);
    failures += ic_align_check(!reference.status && fabs(reference.stats.energy_before-expected_sd)<1e-12 &&
                               fabs(reference.stats.direction_angle_before-.09)<1e-12 &&
                               reference.physical.area==1 && reference.physical.protected_area==.5 &&
                               !reference.physical.bad10_area_before,
                               "independent original metric retains area, adopted source mask and measured angle");
    /* Reproduce the planner branch whose uncut embedding failed its contact
     * audit. The valid cut proposal must still satisfy the original metric. */
    uncut = reference; uncut.valid = 0; cut.valid = 1;
    for (size_t trial = 0; trial < 2; trial++) {
        int rc = CloseupAlign_measure(arena,xyz,4,faces,2,trial?good:worse,arms,protected_faces,1,&metric);
        cut.stats.energy_after = metric.area*metric.sd_mean;
        cut.stats.direction_energy_after = CLOSEUP_FIBER_DIRECTION_WEIGHT*metric.direction_energy;
        cut.stats.objective_after = cut.stats.energy_after+cut.stats.direction_energy_after;
        cut.stats.direction_angle_after = 0;
        cut.stats.direction_faces = metric.measured_faces;
        cut.physical.area = metric.area; cut.physical.protected_area = metric.protected_area;
        cut.physical.bad10_area_after = metric.bad10_area;
        cut.physical.bad25_area_after = metric.bad25_area;
        cut.physical.protected_bad10_area_after = metric.protected_bad10_area;
        cut.physical.protected_bad25_area_after = metric.protected_bad25_area;
        failures += ic_align_check(!rc && ic_align_prefer_cut(&uncut,&cut,&reference)==(int)trial,
                                   trial ? "contact-invalid uncut still permits a measured improving cut" :
                                           "contact-invalid uncut cannot admit an aligned but stretched cut");
    }
    uncut.status = -1;
    failures += ic_align_check(ic_align_prefer_cut(&uncut,&cut,&reference),
                               "failed uncut solver can recover using independently measured original metric");
    reference.status = -1;
    failures += ic_align_check(!ic_align_prefer_cut(&uncut,&cut,&reference),
                               "unknown original metric cannot qualify a cut");
    reference.status = ic_align_original_measure(&mesh,attributes,NULL,&reference.stats,&reference.physical);
    failures += ic_align_check(!reference.status && !reference.stats.direction_faces &&
                               isnan(reference.stats.direction_angle_before) &&
                               !reference.stats.direction_energy_before,
                               "original metric without physical evidence makes no angular claim");
    memcpy(&saved_stats,&reference.stats,sizeof saved_stats);
    memcpy(&saved_physical,&reference.physical,sizeof saved_physical);
    for (size_t v = 0; v < 4; v++) seed[2*v+1] = 0;
    int rc = ic_align_original_measure(&mesh,attributes,arms,&reference.stats,&reference.physical);
    failures += ic_align_check(rc && !memcmp(&saved_stats,&reference.stats,sizeof saved_stats) &&
                               !memcmp(&saved_physical,&reference.physical,sizeof saved_physical),
                               "degenerate original metric refuses without modifying measured outputs");
    Arena_dispose(&arena);
    return failures;
}

static int ic_rigid_contact_selftest(void)
{
    float xyz[24];double uv[16],out[16];
    int32_t faces[12]={0,1,2,0,2,3,4,5,6,4,6,7};uint8_t attr[8]={1,1,1,1,1,1,1,1};
    SlimRefineDirection arms[4]={{0},{0},{0},{0}};
    for (int c=0;c<2;c++) for (int k=0;k<4;k++) {
        int v=4*c+k;double u=(k==1||k==2?32:-32)+66*c,w=k>=2?32:-32;
        xyz[3*v]=(float)u;xyz[3*v+1]=(float)w;xyz[3*v+2]=0;uv[2*v]=u;uv[2*v+1]=w;
    }
    for (int f=0;f<4;f++) {double a=f<2?.6:0;arms[f].h[0]=cos(a);arms[f].h[1]=sin(a);arms[f].v[0]=-sin(a);arms[f].v[1]=cos(a);arms[f].weight=1;}
    MeshBinData mesh={0};mesh.nv=8;mesh.nf=4;mesh.verts=xyz;mesh.faces=faces;mesh.uv64=uv;
    CloseupAlignOpts opts;CloseupAlign_defaults(&opts);IcProjectionLinear cache={0};ic_projection_options(&opts,&cache);
    opts.slim.threads=1;opts.slim.verbose=0;
    SlimRefineStats stats={0};IcPhysicalAlignStats detail={0};AsmContactStats contacts={0};
    Arena_T arena=Arena_new_sized(4096), contact_arena=Arena_new_sized(4096);
    int rc=ic_align_components(arena,&mesh,arms,attr,&opts,out,&stats,&detail,NULL,NULL,1,contact_arena,1);
    int failed=rc || !detail.global_contact_refusals || !detail.proposed || !(stats.direction_angle_after<stats.direction_angle_before) ||
        ic_material_measure(&mesh,out,NULL,&contacts) || !ic_material_qualified(&contacts);
    arms[1].weight=arms[3].weight=0;memset(&stats,0,sizeof stats);memset(&detail,0,sizeof detail);
    Arena_free(arena);Arena_free(contact_arena);
    rc=ic_align_components(arena,&mesh,arms,attr,&opts,out,&stats,&detail,NULL,NULL,1,contact_arena,1);
    failed|=rc || ic_fiber_qualified(&stats) || fabs(stats.direction_area_fraction-.5)>1e-12;
    fprintf(stderr,"[rigid material selftest] %s; overlap refused, smaller safe turn accepted, and partial evidence retained\n",failed?"FAIL":"ok");
    Arena_dispose(&contact_arena);Arena_dispose(&arena);
    return failed;
}

static int ic_physical_alignment_selftest(void)
{
    enum { NX=4, NY=3, PART_NV=NX*NY, PART_NF=2*(NX-1)*(NY-1), NV=2*PART_NV+4, NF=2*PART_NF+1 };
    float xyz[3*NV], original_xyz[3*NV];
    int32_t faces[3*NF], original_faces[3*NF];
    double uv[2*NV], original_uv[2*NV], result[2*NV], sentinel[2*NV];
    SlimRefineDirection arms[NF], original_arms[NF];
    uint8_t attributes[NV];
    SlimRefineStats stats = {0}, saved_stats;
    IcPhysicalAlignStats information = {0}, saved_information;
    CloseupAlignOpts options;
    MeshBinData mesh = {0};
    int failures = ic_align_original_reference_selftest()+ic_rigid_contact_selftest();
    memset(arms,0,sizeof arms); memset(attributes,1,sizeof attributes);
    CloseupAlign_defaults(&options);
    options.slim.linear_solve = ic_slim_linear;
    options.slim.iterations = 20; options.slim.verbose = 0; options.slim.threads = 2;
    for (int c = 0; c < 2; c++) {
        double scale = c ? 2 : 1, angle = c ? .47 : .32, sign = c ? -1 : 1;
        double origin[2] = {c?100:10,c?-50:20};
        for (int y = 0; y < NY; y++) for (int x = 0; x < NX; x++) {
            int v = c*PART_NV+y*NX+x;
            double u = scale*1.06*x, w = scale*.97*y;
            xyz[3*v] = (float)(scale*x); xyz[3*v+1] = (float)(scale*y); xyz[3*v+2] = (float)(10*c);
            uv[2*v] = origin[0]+cos(angle)*u-sin(angle)*w;
            uv[2*v+1] = origin[1]+sign*(sin(angle)*u+cos(angle)*w);
        }
        for (int y = 0, k = 0; y < NY-1; y++) for (int x = 0; x < NX-1; x++) {
            int v = c*PART_NV+y*NX+x, first = 2*k+c;
            k++;
            int second = 2*k+c;
            k++;
            faces[3*first] = v; faces[3*first+1] = v+1; faces[3*first+2] = v+NX+1;
            faces[3*second] = v; faces[3*second+1] = v+NX+1; faces[3*second+2] = v+NX;
            arms[first].h[0] = arms[second].h[0] = 1;
            arms[first].v[1] = arms[second].v[1] = 1;
            arms[first].weight = arms[second].weight = c?.9:.2;
        }
    }
    for (int k = 0; k < 3; k++) {
        int v = 2*PART_NV+k;
        xyz[3*v] = k==1?1.f:0.f; xyz[3*v+1] = k==2?1.f:0.f; xyz[3*v+2] = 30;
        uv[2*v] = 300+xyz[3*v]; uv[2*v+1] = 400+xyz[3*v+1];
        faces[3*(NF-1)+k] = v;
    }
    xyz[3*(NV-1)] = 999; xyz[3*(NV-1)+1] = 888; xyz[3*(NV-1)+2] = 777;
    uv[2*(NV-1)] = -0.0; uv[2*(NV-1)+1] = 123456789.125;
    attributes[NV-1] = 0;
    memcpy(original_xyz,xyz,sizeof xyz); memcpy(original_faces,faces,sizeof faces);
    memcpy(original_uv,uv,sizeof uv); memcpy(original_arms,arms,sizeof arms);
    for (int k = 0; k < 2*NV; k++) result[k] = sentinel[k] = -12345;
    mesh.nv = NV; mesh.nf = NF; mesh.verts = xyz; mesh.faces = faces; mesh.uv64 = uv;
    int rc = ic_reparameterize_physical(&mesh,arms,attributes,&options,result,&stats,&information,NULL);
    failures += ic_align_check(!rc && information.components==3 && information.orphans==1 &&
                               information.proposed==2 && information.unknown==1 && ic_fiber_qualified(&stats),
                               "interleaved complete components align with partial evidence and one exact orphan");
    failures += ic_align_check(!memcmp(xyz,original_xyz,sizeof xyz) && !memcmp(faces,original_faces,sizeof faces) &&
                               !memcmp(uv,original_uv,sizeof uv) && !memcmp(arms,original_arms,sizeof arms) &&
                               !memcmp(result+4*PART_NV,uv+4*PART_NV,8*sizeof(double)),
                               "XYZ, faces, physical arms, unknown component and signed-zero orphan remain exact");
    double expected_mass = 6*.2+24*.9;
    double expected_angle = (6*.2*.32*.32+24*.9*.47*.47)/expected_mass;
    failures += ic_align_check(fabs(information.area-30.5)<1e-12 &&
                               fabs(information.confidence_mass-expected_mass)<1e-12 &&
                               fabs(stats.direction_area_fraction-30/30.5)<1e-12 &&
                               fabs(stats.direction_angle_before-expected_angle)<1e-12,
                               "RMS uses evidence times physical area rather than an average of component RMS values");
    SlimRefineStats original_stats = {0};
    IcPhysicalAlignStats original_physical = {0};
    rc = ic_align_original_measure(&mesh,attributes,arms,&original_stats,&original_physical);
    failures += ic_align_check(!rc && fabs(original_stats.energy_before-stats.energy_before)<1e-9 &&
                               fabs(original_stats.objective_before-stats.objective_before)<1e-9 &&
                               fabs(original_stats.direction_angle_before-stats.direction_angle_before)<1e-12 &&
                               original_physical.area==information.area &&
                               original_physical.protected_area==information.protected_area,
                               "original reference measures disconnected mixed-handedness components without recertifying them");
    for (int c = 0; c < 2; c++) {
        double before[2] = {0}, after[2] = {0};
        int orientation_ok = 1, moved = 0;
        for (int k = 0; k < PART_NF; k++) {
            const int32_t *f = faces+3*(2*k+c);
            for (int a = 0; a < 2; a++) for (int j = 0; j < 3; j++) {
                before[a] += uv[2*f[j]+a]/(3.0*PART_NF);
                after[a] += result[2*f[j]+a]/(3.0*PART_NF);
            }
            double det = (result[2*f[1]]-result[2*f[0]])*(result[2*f[2]+1]-result[2*f[0]+1])-
                         (result[2*f[1]+1]-result[2*f[0]+1])*(result[2*f[2]]-result[2*f[0]]);
            orientation_ok &= c ? det<0 : det>0;
        }
        const int corners[4] = {0,NX-1,NX*(NY-1),PART_NV-1};
        for (int j = 0; j < 4; j++) {
            int v = c*PART_NV+corners[j];
            moved += hypot(result[2*v]-uv[2*v],result[2*v+1]-uv[2*v+1])>1e-3;
        }
        failures += ic_align_check(orientation_ok && moved==4 && hypot(before[0]-after[0],before[1]-after[1])<1e-10,
                                   "component handedness and centroid survive while all boundary corners move");
    }
    Arena_T arena = Arena_new();
    CloseupAlignMetric metric = {0};
    int measured = CloseupAlign_measure(arena,xyz,NV,faces,NF,result,arms,NULL,1,&metric);
    failures += ic_align_check(!measured && fabs(stats.energy_after-metric.area*metric.sd_mean)<1e-9 &&
                               fabs(stats.direction_energy_after-16*metric.direction_energy)<1e-9 &&
                               fabs(sqrt(stats.direction_angle_after)*57.2957795130823208768-metric.fiber_rms_degrees)<1e-9,
                               "scattered output independently matches fixed-weight objective and joint angular measurement");
    Arena_dispose(&arena);
    AsmContactStats contact = {0};
    (void)ic_material_measure(&mesh,result,NULL,&contact);
    failures += ic_align_check(ic_material_qualified(&contact),"assembled output passes the complete material contact audit");
    {
        double cached_result[2*NV];
        SlimRefineStats cached_stats = {0}; IcPhysicalAlignStats cached_info = {0};
        IcProjectionLinear cache = {0}; CloseupAlignOpts cached = options;
        ic_projection_options(&cached,&cache);
        int status = ic_reparameterize_physical(&mesh,arms,attributes,&cached,cached_result,&cached_stats,&cached_info,NULL);
        if (status || memcmp(cached_result,result,sizeof result)) {
            double difference = 0;
            for (int k = 0; k < 2*NV; k++) difference = fmax(difference,fabs(cached_result[k]-result[k]));
            fprintf(stderr,"[closeup projection integration] status=%d max-UV-difference=%.17g "
                    "phases=%zu factors=%zu solves=%zu energy=%.17g/%.17g angle=%.17g/%.17g\n",
                    status,difference,cache.phases,cache.factorizations,cache.solves,
                    stats.energy_after,cached_stats.energy_after,stats.direction_angle_after,cached_stats.direction_angle_after);
        }
        failures += ic_align_check(!status && !memcmp(cached_result,result,sizeof result) &&
            cached_info.proposed == information.proposed && cached_info.unknown == information.unknown &&
            cached_stats.energy_after == stats.energy_after &&
            cached_stats.direction_angle_after == stats.direction_angle_after &&
            cache.phases == 2 && cache.factorizations == 16 && cache.solves == 16 && !cache.factor && !cache.storage,
            "small-system backend preserves exact assembled UV with independent component scopes");
    }
    rc = ic_reparameterize_physical(&mesh,NULL,attributes,&options,result,&stats,&information,NULL);
    failures += ic_align_check(!rc && information.unknown==3 && !memcmp(result,uv,sizeof uv) &&
                               isnan(stats.direction_angle_after) && !ic_fiber_qualified(&stats),
                               "all-unknown evidence preserves the whole mesh without an alignment claim");
    /* The final component fails after earlier components have produced usable
     * proposals. A failed whole transaction must not leak those partial writes. */
    uv[4*PART_NV+5] = uv[4*PART_NV+1];
    memcpy(result,sentinel,sizeof result);
    memset(&stats,0x5a,sizeof stats); memcpy(&saved_stats,&stats,sizeof stats);
    memset(&information,0x5a,sizeof information); memcpy(&saved_information,&information,sizeof information);
    rc = ic_reparameterize_physical(&mesh,arms,attributes,&options,result,&stats,&information,NULL);
    failures += ic_align_check(rc && !memcmp(result,sentinel,sizeof result) &&
                               !memcmp(&stats,&saved_stats,sizeof stats) && !memcmp(&information,&saved_information,sizeof information),
                               "late infeasible component leaves all output UV and summaries unchanged");
    memcpy(uv,original_uv,sizeof uv);
    faces[3*(NF-1)+2] = NV;
    rc = ic_reparameterize_physical(&mesh,arms,attributes,&options,result,&stats,&information,NULL);
    failures += ic_align_check(rc && !memcmp(result,sentinel,sizeof result),"out-of-range source topology refuses atomically");
    failures += ic_local_global_integration_test();
    fprintf(stderr,"[closeup physical integration] %s: %d failure(s)\n",failures?"FAIL":"PASS",failures);
    return failures;
}
