static double fr_test_score(const float *xyz, const int32_t *faces, size_t nf,
                            const double *uv, const SlimRefineDirection *directions)
{
    double sum = 0, mass = 0;
    for (size_t f = 0; f < nf; f++) {
        double area = 0, mean = 0, delta = 0, score = 0;
        if (fr_project(xyz, faces+3*f, uv, directions+f, &area, &mean, &delta, &score)) return NAN;
        sum += area*directions[f].weight*score;
        mass += area*directions[f].weight;
    }
    return sqrt(sum/mass)*fr_degrees;
}

int FiberRigid_selftest(void)
{
    Arena_T arena = Arena_new();
    int failed = 0, checks = 0;
#define FR_CHECK(test) do { checks++; if (!(test)) { failed++; \
    fprintf(stderr, "FAIL fiber rigid line %d: %s\n", __LINE__, #test); } } while (0)
    float xyz[21] = {0,0,0, 2,0,0, 0,3,0, 4,0,0, 6,0,0, 4,3,0, 99,99,99};
    int32_t faces[6] = {0,1,2, 3,4,5};
    SlimRefineDirection directions[2] = {{{1,0,0},{0,1,0},1}, {{1,0,0},{0,1,0},1}};
    double uv[14], out[14], trial[14];
    FiberRigidStats stats;
    for (int degrees = -44; degrees <= 44; degrees += 4) {
        double angle = (double)degrees/fr_degrees, co = cos(angle), si = sin(angle);
        for (size_t v = 0; v < 7; v++) {
            uv[2*v] = 2e6+co*xyz[3*v]-si*xyz[3*v+1];
            uv[2*v+1] = -3e6+si*xyz[3*v]+co*xyz[3*v+1];
        }
        FR_CHECK(FiberRigid_fit(arena, xyz, 7, faces, 2, uv, directions, &stats) == 0);
        FR_CHECK(fabs(stats.radians+angle) < 1e-9 && stats.rms_after_degrees < 1e-7);
        FR_CHECK(FiberRigid_apply(uv, 7, &stats, out) == 0);
        FR_CHECK(fr_test_score(xyz, faces, 2, out, directions) < 1e-7);
        double old = hypot(uv[12]-uv[0], uv[13]-uv[1]);
        FR_CHECK(fabs(hypot(out[12]-out[0], out[13]-out[1])-old) < 2e-9);
        memcpy(trial, uv, sizeof uv);
        FR_CHECK(FiberRigid_apply(trial, 7, &stats, trial) == 0 && memcmp(trial, out, sizeof out) == 0);
    }
    /* Independent complementary-axis dense search: mixed local directions,
     * unequal confidence, shear, reflections and arbitrary family signs. */
    for (int c = 0; c < 32; c++) {
        for (size_t v = 0; v < 7; v++) {
            uv[2*v] = xyz[3*v]+(.13+.01*c)*xyz[3*v+1];
            uv[2*v+1] = (c%2 ? -1 : 1)*(.8+.01*c)*xyz[3*v+1];
        }
        for (size_t f = 0; f < 2; f++) {
            double a = (double)(c*37+(int)f*73)/fr_degrees;
            directions[f].h[0] = cos(a); directions[f].h[1] = sin(a);
            directions[f].v[0] = -sin(a); directions[f].v[1] = cos(a);
            directions[f].weight = f ? .13 : .91;
            if (c%3) for (size_t k = 0; k < 3; k++) directions[f].v[k] *= -1;
        }
        FR_CHECK(FiberRigid_fit(arena, xyz, 7, faces, 2, uv, directions, &stats) == 0);
        FR_CHECK(FiberRigid_apply(uv, 7, &stats, out) == 0);
        double actual = fr_test_score(xyz, faces, 2, out, directions);
        FR_CHECK(fabs(actual-stats.rms_after_degrees) < 1e-9);
        int grid_pass = 1;
        for (int k = 0; k <= 720; k++) {
            FiberRigidStats grid = stats;
            grid.radians = -.5*fr_period+fr_period*(double)k/720;
            if (FiberRigid_apply(uv, 7, &grid, trial) ||
                    actual > fr_test_score(xyz, faces, 2, trial, directions)+1e-9) grid_pass = 0;
        }
        FR_CHECK(grid_pass);
    }
    directions[0].weight = directions[1].weight = 0;
    FR_CHECK(FiberRigid_fit(arena, xyz, 7, faces, 2, uv, directions, &stats) != 0);
    directions[0].weight = 1; directions[1].weight = NAN;
    FR_CHECK(FiberRigid_fit(arena, xyz, 7, faces, 2, uv, directions, &stats) != 0);
    directions[1].weight = 0;
    FR_CHECK(FiberRigid_fit(arena, xyz, 7, faces, 2, uv, directions, &stats) == 0 && stats.supported_faces == 1);
    faces[0] = -1;
    FR_CHECK(FiberRigid_fit(arena, xyz, 7, faces, 2, uv, directions, &stats) != 0);
    faces[0] = 0;
    uv[0] = NAN;
    FR_CHECK(FiberRigid_fit(arena, xyz, 7, faces, 2, uv, directions, &stats) != 0);
    /* Point fans have no bearing on a common rigid transform. */
    uv[0] = 0; faces[3] = 0;
    FR_CHECK(FiberRigid_fit(arena, xyz, 7, faces, 2, uv, directions, &stats) == 0);
    Arena_dispose(&arena);
    fprintf(stderr, "FiberRigid: %d checks, %d failures\n", checks, failed);
#undef FR_CHECK
    return failed ? 1 : 0;
}
