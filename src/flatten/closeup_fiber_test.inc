static int cfiber_test_volume(const char *root, int missing, int bent)
{
    enum { S = 256, D = 64 };
    char path[1800];
    uint8_t *volume = malloc((size_t)S*S*D);
    if (!volume) return -1;
    double angle = 20*3.14159265358979323846/180;
    for (int y = 0; y < S; y++) for (int x = 0; x < S; x++) {
        double px=bent ? ((double)x-128)/cos(15*3.14159265358979323846/180) : (double)x;
        double a = px*cos(angle)+y*sin(angle), b = -px*sin(angle)+y*cos(angle);
        double intensity = 128+50*cos(2*3.14159265358979323846*a/12)+25*cos(2*3.14159265358979323846*b/12);
        for (int z = 0; z < D; z++) volume[((size_t)z*S+(size_t)y)*S+(size_t)x] = missing && x<128 && y<128 ? 0 : (uint8_t)intensity;
    }
    ves_mkdir(root); snprintf(path,sizeof path,"%s/0",root); ves_mkdir(path);
    snprintf(path,sizeof path,"%s/0/.zarray",root);
    FILE *fp = fopen(path,"wb");
    if (!fp) { free(volume); return -1; }
    fprintf(fp,"{\"shape\":[%d,%d,%d],\"chunks\":[%d,%d,%d],\"dtype\":\"|u1\",\"fill_value\":0,"
               "\"order\":\"C\",\"filters\":null,\"dimension_separator\":\"/\",\"compressor\":null,\"zarr_format\":2}",D,S,S,D,S,S);
    int rc = fclose(fp) != 0;
    snprintf(path,sizeof path,"%s/0/0",root); ves_mkdir(path);
    snprintf(path,sizeof path,"%s/0/0/0",root); ves_mkdir(path);
    snprintf(path,sizeof path,"%s/0/0/0/0",root); fp = fopen(path,"wb");
    if (!fp) { free(volume); return -1; }
    if (fwrite(volume,1,(size_t)S*S*D,fp) != (size_t)S*S*D) rc = -1;
    if (fclose(fp)) rc = -1;
    free(volume); return rc;
}

static int cfiber_test_borrowed(const char *raw,const float *xyz,const float *normal,
                                const double *uv,size_t nv,const int32_t *faces,size_t nf,
                                const SlimRefineDirection *expected,const CloseupFiberStats *reference)
{
    ZarrU8 *volume=calloc(1,sizeof *volume);
    SlimRefineDirection *actual=malloc(nf*sizeof *actual);
    int rc=-1,opened=0;
    if (!volume || !actual || ZarrU8_open(volume,raw)) goto done;
    opened=1;volume->slot_limit=1;
    uint8_t sample=ZarrU8_get(volume,8,16,16);
    ZarrU8CacheStats before={0},after={0};
    ZarrU8_cache_stats(volume,&before);
    for (int pass=0;pass<2;pass++) {
        size_t misses=volume->misses,missing=volume->missing;
        CloseupFiberStats measured={0};
        if (CloseupFiber_measure_volume(volume,xyz,normal,uv,nv,faces,nf,actual,&measured) ||
            memcmp(actual,expected,nf*sizeof *actual) ||
            measured.patches!=reference->patches || measured.measured_patches!=reference->measured_patches ||
            measured.faces!=reference->faces || measured.samples!=reference->samples ||
            memcmp(&measured.area,&reference->area,sizeof(double)) ||
            memcmp(&measured.measured_area,&reference->measured_area,sizeof(double)) ||
            measured.chunk_misses!=volume->misses-misses ||
            measured.chunk_missing!=volume->missing-missing || measured.chunk_misses ||
            volume->slot_limit!=1 || ZarrU8_get(volume,8,16,16)!=sample) goto done;
        ZarrU8_cache_stats(volume,&after);
        if (after.handles!=before.handles || after.chunks!=before.chunks || after.bytes!=before.bytes) goto done;
    }
    {
        /* A failed or empty measurement must not consume the borrowed handle. */
        int32_t invalid[3]={(int32_t)nv,0,1};
        CloseupFiberStats measured={0};
        if (!CloseupFiber_measure_volume(volume,xyz,normal,uv,nv,invalid,1,actual,&measured) ||
            !CloseupFiber_measure_volume(NULL,xyz,normal,uv,nv,faces,1,actual,&measured) ||
            CloseupFiber_measure_volume(volume,xyz,normal,uv,nv,faces,0,actual,&measured) ||
            measured.patches || measured.faces || measured.samples ||
            volume->slot_limit!=1 || ZarrU8_get(volume,8,16,16)!=sample) goto done;
        ZarrU8_cache_stats(volume,&after);
        if (after.handles!=before.handles || after.chunks!=before.chunks || after.bytes!=before.bytes) goto done;
    }
    rc=0;
done:
    if (opened) ZarrU8_close(volume);
    free(volume);free(actual);
    fprintf(stderr,"  physical CT fibers: borrowed/path parity and caller cache lifetime %s\n",rc?"FAILED":"ok");
    return rc;
}

int CloseupFiber_selftest(void)
{
    enum { SIDE = 105, NV = SIDE*SIDE, NF = 2*(SIDE-1)*(SIDE-1) };
    char dir[1400], raw[1600];
    float *xyz = malloc((size_t)3*NV*sizeof(float)), *normal = calloc((size_t)3*NV,sizeof(float));
    double *uv = malloc((size_t)2*NV*sizeof(double));
    int32_t *faces = malloc((size_t)3*NF*sizeof(int32_t));
    SlimRefineDirection *directions = malloc((size_t)NF*sizeof *directions);
    CloseupFiberStats stats = {0}; int rc = 1;
    if (ves_selftest_dir(dir,sizeof dir,"closeup_fiber")) goto done;
    snprintf(raw,sizeof raw,"%s/raw.zarr",dir);
    if (!xyz || !normal || !uv || !faces || !directions || cfiber_test_volume(raw,0,0)) goto done;
    {
        /* Inverting a tangent projection must recover the true 45-degree
         * material arm even when the adjoint would rotate it to 75.96 deg. */
        double angle=60*3.14159265358979323846/180;
        const double patch_normal[3]={0,0,1};
        double face_normal[3]={sin(angle),0,cos(angle)};
        double true_arm[3]={cos(angle)/sqrt(2.0),1/sqrt(2.0),-sin(angle)/sqrt(2.0)};
        double projected[3]={true_arm[0],true_arm[1],0},lifted[3];
        if (cfiber_lift(projected,patch_normal,face_normal,lifted) ||
            cfiber_dot(true_arm,lifted)<1-1e-12 ||
            !(cfiber_dot(patch_normal,face_normal)<CLOSEUP_FIBER_PATCH_NORMAL)) goto done;
    }
    for (int y = 0; y < SIDE; y++) for (int x = 0; x < SIDE; x++) {
        size_t at = (size_t)y*SIDE+(size_t)x;
        xyz[3*at] = 8; xyz[3*at+1] = (float)(16+2*y); xyz[3*at+2] = (float)(16+2*x);
        normal[3*at] = 1;
        uv[2*at] = 1.3*xyz[3*at+2]+.2*xyz[3*at+1]; uv[2*at+1] = .8*xyz[3*at+1];
    }
    size_t nf = 0;
    for (int y = 0; y+1 < SIDE; y++) for (int x = 0; x+1 < SIDE; x++) {
        int32_t a = y*SIDE+x, b = a+1, c = a+SIDE, d = c+1;
        faces[3*nf] = a; faces[3*nf+1] = b; faces[3*nf+2] = d; nf++;
        faces[3*nf] = a; faces[3*nf+1] = d; faces[3*nf+2] = c; nf++;
    }
    for (int mirror = 0; mirror < 2; mirror++) {
        if (CloseupFiber_measure(xyz,normal,uv,NV,faces,NF,raw,directions,&stats) ||
                stats.measured_area < .7*stats.area || !stats.measured_patches) goto done;
        if (!mirror && cfiber_test_borrowed(raw,xyz,normal,uv,NV,faces,NF,directions,&stats)) goto done;
        double a = 20*3.14159265358979323846/180, target[3] = {0,sin(a),cos(a)};
        double maximum = 0;
        for (size_t f = 0; f < NF; f++) if (directions[f].weight > 0) {
            double agreement = fmax(fabs(cfiber_dot(directions[f].h,target)),fabs(cfiber_dot(directions[f].v,target)));
            maximum = fmax(maximum,acos(fmin(1.0,agreement))*180/3.14159265358979323846);
        }
        fprintf(stderr,"  physical CT fibers: mirror=%d measured=%.1f%% max-angle=%.4f deg\n",
                mirror,100*stats.measured_area/stats.area,maximum);
        if (maximum > .3) goto done;
        for (size_t k = 0; k < NV; k++) uv[2*k] = -uv[2*k];
    }
    {
        /* A folded projection cannot choose whichever material face happens
         * to rasterize first, even if both read the same synthetic texture. */
        Arena_T arena=Arena_new();ZarrU8 *volume=calloc(1,sizeof *volume);
        float points[24],coords[16];const int32_t triangles[12]={0,1,2,0,2,3,4,5,6,4,6,7};
        const double square[8]={0,0,80,0,80,80,0,80};
        for (size_t v=0;v<8;v++) {
            points[3*v]=(float)(v<4?8:12);points[3*v+1]=40+(float)square[2*(v%4)+1];points[3*v+2]=40+(float)square[2*(v%4)];
            coords[2*v]=(float)square[2*(v%4)];coords[2*v+1]=(float)square[2*(v%4)+1];
        }
        AsmChart chart={0};chart.xyz=points;chart.uv=coords;chart.faces=(int32_t *)triangles;chart.nv=8;chart.nf=4;
        AsmFibreMoments moments={0};size_t samples[4]={0};int bad=!volume || ZarrU8_open(volume,raw);
        if (!bad) {
            AsmFibre_measure_patch(volume,arena,&chart,&moments,samples);
            if (moments.px || samples[0] || samples[1] || samples[2] || samples[3]) bad=1;
            /* A single sheet remains measurable, shared diagonal included. */
            chart.nv=4;chart.nf=2;AsmFibre_measure_patch(volume,arena,&chart,&moments,samples);
            if (moments.px<2000 || !samples[0] || !samples[1]) bad=1;
            ZarrU8_close(volume);
        }
        free(volume);Arena_dispose(&arena);
        if (bad) { fprintf(stderr,"  physical CT fibers: ambiguous projection support failed\n");goto done; }
    }
    {
        /* Compressed indexing UV makes a sufficiently large physical halo to
         * measure a cross while the primary cell's CT is entirely missing.
         * Those primary faces must not inherit measured coverage. */
        snprintf(raw,sizeof raw,"%s/missing.zarr",dir);
        if (cfiber_test_volume(raw,1,0)) goto done;
        for (size_t v=0;v<NV;v++) { uv[2*v]=.5*xyz[3*v+2];uv[2*v+1]=.5*xyz[3*v+1]; }
        if (CloseupFiber_measure(xyz,normal,uv,NV,faces,NF,raw,directions,&stats) || !stats.measured_patches) goto done;
        size_t missing_faces=0,measured_faces=0;
        for (size_t f=0;f<NF;f++) {
            double x=0,y=0;
            for (size_t k=0;k<3;k++) { x+=xyz[3*(size_t)faces[3*f+k]+2]/3.0;y+=xyz[3*(size_t)faces[3*f+k]+1]/3.0; }
            if (x<126 && y<126) { missing_faces++;if (directions[f].weight>0) goto done; }
            else if (directions[f].weight>0) measured_faces++;
        }
        if (!missing_faces || !measured_faces) goto done;
        fprintf(stderr,"  physical CT fibers: missing primary faces=%zu unsupported, measured halo faces=%zu\n",missing_faces,measured_faces);
    }
    {
        /* Physical CT on a gently bent sheet. Both sides project into one
         * patch and need distinct inverse tangent lifts. A sharper bend must
         * be left unmeasured by the contributor and target angular guards. */
        snprintf(raw,sizeof raw,"%s/bent.zarr",dir);
        if (cfiber_test_volume(raw,0,1)) goto done;
        for (int sharp=0;sharp<2;sharp++) {
            double tilt=(sharp?45:15)*3.14159265358979323846/180;
            for (size_t v=0;v<NV;v++) {
                double x=xyz[3*v+2],slope=(x<128?-1:1)*tan(tilt);
                xyz[3*v]=(float)(8+tan(tilt)*fabs(x-128));
                normal[3*v]=(float)cos(tilt);normal[3*v+1]=0;normal[3*v+2]=(float)(-slope*cos(tilt));
                uv[2*v]=.25*(x-128)+32;uv[2*v+1]=.25*(xyz[3*v+1]-128)+32;
            }
            if (CloseupFiber_measure(xyz,normal,uv,NV,faces,NF,raw,directions,&stats)) goto done;
            double maximum=0,a=20*3.14159265358979323846/180;
            for (size_t f=0;f<NF;f++) if (directions[f].weight>0) {
                double x=0;for (size_t k=0;k<3;k++) x+=xyz[3*(size_t)faces[3*f+k]+2]/3.0;
                double target[3]={(x<128?-1:1)*sin(tilt)*cos(a),sin(a),cos(tilt)*cos(a)};
                double agreement=fmax(fabs(cfiber_dot(directions[f].h,target)),fabs(cfiber_dot(directions[f].v,target)));
                maximum=fmax(maximum,acos(fmin(1.0,agreement))*180/3.14159265358979323846);
            }
            fprintf(stderr,"  physical CT fibers: tilt=%d measured=%.1f%% max-angle=%.4f deg\n",sharp?45:15,
                    100*stats.measured_area/stats.area,maximum);
            if (sharp ? stats.faces!=0 : (stats.measured_area<.6*stats.area || maximum>1.5)) goto done;
        }
    }
    rc = 0;
done:
    free(xyz); free(normal); free(uv); free(faces); free(directions);
    if (!rc) ves_selftest_clean(dir);
    fprintf(stderr,"[closeup fiber selftest] %s\n",rc ? "FAILED" : "ok");
    return rc;
}
