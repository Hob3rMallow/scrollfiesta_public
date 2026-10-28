static uint8_t cfv_test_encode(double value)
{
    double encoded = floor(127*value+128+.5);
    return (uint8_t)fmax(0,fmin(255,encoded));
}

static int cfv_test_field(const char *root)
{
    enum { SIDE = 48 };
    const size_t count = (size_t)SIDE*SIDE*SIDE;
    uint8_t *data = calloc(3*count,1);
    char path[1800];
    const char *names[3] = {"presence","nx","ny"};
    double angle = 20*3.14159265358979323846/180, tilt = 30*3.14159265358979323846/180;
    int rc = -1;
    if (!data) return -1;
    for (int z = 0; z < SIDE; z++) for (int y = 0; y < SIDE; y++) for (int x = 0; x < SIDE; x++) {
        size_t at = ((size_t)z*SIDE+(size_t)y)*SIDE+(size_t)x;
        double first[3] = {0,0,1}, second[3] = {0,1,0}, direction[3] = {0};
        double phase = angle;
        if (x < 24 && y >= 20) continue; /* Missing support, including every depth. */
        if (x >= 24 && y >= 20) {
            first[0] = sin(tilt); first[2] = cos(tilt);
            if (-cos(tilt)*(z-24)+sin(tilt)*(x-32) >= 0) phase += 3.14159265358979323846/2;
        } else if (x >= 24) {
            if (z == 16) continue;
            phase = z < 16 ? 0 : 3.14159265358979323846/4; /* Equal cancelling cross phases. */
        } else if (z >= 16) phase += 3.14159265358979323846/2;
        for (int k = 0; k < 3; k++) direction[k] = cos(phase)*first[k]+sin(phase)*second[k];
        if (direction[0] < 0) for (int k = 0; k < 3; k++) direction[k] = -direction[k];
        data[at] = 255; data[count+at] = cfv_test_encode(direction[2]);
        data[2*count+at] = cfv_test_encode(direction[1]);
    }
    for (int channel = 0; channel < 3; channel++) {
        snprintf(path,sizeof path,"%s/%s/0/.zarray",root,names[channel]); ves_ensure_parent_dir(path);
        FILE *fp = fopen(path,"wb");
        if (!fp) goto done;
        fprintf(fp,"{\"shape\":[%d,%d,%d],\"chunks\":[%d,%d,%d],\"dtype\":\"|u1\",\"fill_value\":0,"
                   "\"order\":\"C\",\"filters\":null,\"dimension_separator\":\"/\",\"compressor\":null,\"zarr_format\":2}",
                   SIDE,SIDE,SIDE,SIDE,SIDE,SIDE);
        if (fclose(fp)) goto done;
        snprintf(path,sizeof path,"%s/%s/0/0/0/0",root,names[channel]); ves_ensure_parent_dir(path);
        fp = fopen(path,"wb");
        if (!fp) goto done;
        int failed = fwrite(data+(size_t)channel*count,1,count,fp) != count;
        if (fclose(fp)) failed = 1;
        if (failed) goto done;
    }
    rc = 0;
done:
    free(data); return rc;
}

int CloseupFiberVolume_selftest(void)
{
    char dir[1400] = {0}, field[1600];
    float xyz[36] = {0}, normals[36] = {0};
    int32_t faces[12] = {0,1,2,3,4,5,6,7,8,9,10,11};
    SlimRefineDirection directions[4] = {0};
    CloseupFiberStats stats = {0};
    double tilt = 30*3.14159265358979323846/180, angle = 20*3.14159265358979323846/180;
    double centers[4][3] = {{16,8,8},{24,32,32},{16,32,8},{16,8,32}};
    double expected[2][3] = {{0,sin(angle),cos(angle)},
                             {cos(angle)*sin(tilt),sin(angle),cos(angle)*cos(tilt)}};
    int fails = 0;
    if (ves_selftest_dir(dir,sizeof dir,"closeup_fiber_volume")) return 1;
    snprintf(field,sizeof field,"%s/field",dir);
    if (cfv_test_field(field)) { fails++; goto done; }
    for (int f = 0; f < 4; f++) {
        double first[3] = {0,0,1}, second[3] = {0,1,0}, normal[3] = {0};
        if (f == 1) { first[0] = sin(tilt); first[2] = cos(tilt); }
        cfv_cross(first,second,normal);
        for (int corner = 0; corner < 3; corner++) for (int k = 0; k < 3; k++) {
            size_t at = 3*(size_t)(3*f+corner)+(size_t)k;
            xyz[at] = (float)(centers[f][k]+(corner == 1 ? 6 : -3)*first[k]+(corner == 2 ? 6 : -3)*second[k]);
            normals[at] = f == 0 ? 0 : (float)normal[k]; /* Exact face-normal fallback. */
        }
    }
    float original_xyz[36], original_normals[36]; int32_t original_faces[12];
    memcpy(original_xyz,xyz,sizeof xyz); memcpy(original_normals,normals,sizeof normals);
    memcpy(original_faces,faces,sizeof faces);
    if (CloseupFiberVolume_measure(xyz,normals,12,faces,4,field,directions,&stats) ||
            stats.faces != 2 || stats.measured_patches != 2 || stats.patches != 4 ||
            stats.samples != 2*(size_t)(2*FIBER_FLATTEN_HALF_DEPTH+1)+(size_t)(2*FIBER_FLATTEN_HALF_DEPTH) ||
            !(stats.measured_area > 0 && stats.measured_area < stats.area)) { fails++; goto done; }
    for (int f = 0; f < 2; f++) {
        double agreement = fmax(fabs(cfv_dot(directions[f].h,expected[f])),fabs(cfv_dot(directions[f].v,expected[f])));
        double error = acos(fmin(1.0,agreement))*180/3.14159265358979323846;
        fprintf(stderr,"  inferred physical fibers: face=%d error=%.6f deg evidence=%.6f\n",f,error,directions[f].weight);
        if (error > 1 || directions[f].weight < .95 || fabs(cfv_dot(directions[f].h,directions[f].v)) > 1e-12) fails++;
    }
    SlimRefineDirection zero = {0};
    if (memcmp(directions+2,&zero,sizeof zero) || memcmp(directions+3,&zero,sizeof zero) ||
            memcmp(original_xyz,xyz,sizeof xyz) || memcmp(original_normals,normals,sizeof normals) ||
            memcmp(original_faces,faces,sizeof faces)) fails++;
    faces[0] = 12;
    if (!CloseupFiberVolume_measure(xyz,normals,12,faces,4,field,directions,&stats)) fails++;
    faces[0] = 0; xyz[0] = NAN;
    if (!CloseupFiberVolume_measure(xyz,normals,12,faces,4,field,directions,&stats)) fails++;
    if (CloseupFiberVolume_measure(NULL,NULL,0,NULL,0,NULL,NULL,&stats) || stats.faces || stats.samples) fails++;
done:
    if (!fails) ves_selftest_clean(dir);
    fprintf(stderr,"[closeup fiber volume selftest] %s (%d failures)\n",fails ? "FAILED" : "ok",fails);
    return fails ? 1 : 0;
}
