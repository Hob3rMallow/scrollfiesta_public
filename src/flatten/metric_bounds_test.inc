/* Native analytic boundary and filter-vs-exact controls. */
extern "C" int MetricBounds_selftest(void)
{
    size_t cases = 0, failures = 0, exact = 0, filtered = 0;
    const double bounds[4] = {1/1.1,1.1,1/1.25,1.25};
    const float unit[9] = {0,0,0,1,0,0,0,1,0};
    for (int axis = 0; axis < 2; axis++) for (int k = 0; k < 4; k++)
        for (int offset = -1; offset <= 1; offset++) {
            double scale = bounds[k];
            if (offset) scale = std::nextafter(scale,offset < 0 ? 0 : 2);
            double uv[6] = {0,0,1,0,0,1}; uv[axis ? 5 : 2] = scale;
            MetricBoundsResult r; int rc = MetricBounds_classify(unit,uv,bounds,&r);
            unsigned expected = 0;
            for (int b = 0; b < 4; b++)
                if (b & 1 ? scale > bounds[b] : scale < bounds[b]) expected |= 1u << b;
            failures += rc || r.violated != expected; cases++;
        }
    /* Integer translations and power-of-two scales preserve these stored
     * triangles exactly, including badly conditioned physical triangles. */
    for (int e = -40; e <= 40; e += 10) for (int shift = 0; shift < 3; shift++)
        for (int sample = 0; sample < 80; sample++) {
            float h = std::ldexp(1.0f,-(sample%21));
            float x[9] = {0,0,0,1,0,0,1,h,0};
            double q = std::ldexp(1.0,e), t = shift == 0 ? 0 : std::ldexp(1.0,shift+2);
            double u[6] = {t,-t,t+q*(.85+.005*sample),-t+.02*q,
                t+q*(.85+.005*sample)+.03*q*h,-t+.02*q+q*h};
            for (int i = 0; i < 9; i++) x[i] = (float)(x[i]*q);
            MetricBoundsResult r; unsigned expected = 0;
            int rc = MetricBounds_classify(x,u,bounds,&r);
            int erc = mb_classify<BSNumber>(x,u,bounds,&expected);
            failures += rc != erc || (!rc && r.violated != expected); cases++;
            if (!rc) { exact += r.exact_fallback != 0; filtered += r.exact_fallback == 0; }
        }
    /* Huge, subnormal, reflected and zero-area UVs still get exact signs. */
    for (int e : {-1074,-1022,-500,500,1000}) {
        double scale = std::ldexp(1.0,e), u[6] = {0,0,scale,0,0,-scale};
        MetricBoundsResult r; unsigned expected = 0;
        int rc = MetricBounds_classify(unit,u,bounds,&r);
        int erc = mb_classify<BSNumber>(unit,u,bounds,&expected);
        failures += rc != erc || (!rc && r.violated != expected); cases++;
    }
    MetricBoundsResult sentinel = {123,456}, r = sentinel;
    float flat[9] = {0,0,0,1,0,0,2,0,0}; double uv[6] = {0,0,1,0,0,1};
    failures += MetricBounds_classify(flat,uv,bounds,&r) != -1 || std::memcmp(&r,&sentinel,sizeof r); cases++;
    uv[1] = std::numeric_limits<double>::quiet_NaN();
    failures += MetricBounds_classify(unit,uv,bounds,&r) != -1 || std::memcmp(&r,&sentinel,sizeof r); cases++;
    failures += !exact || !filtered; cases++;
    /* Alternate rounding and DAZ/FTZ cannot change classification of the
     * stored coordinates. The physical subnormal below is half an edge,
     * so dropping it during float-to-double promotion changes the metric. */
    uint32_t encoded[9]={0,0,0,0x00800000,0,0,0x00400000,0x00800000,0};
    float tiny[9];std::memcpy(tiny,encoded,sizeof tiny);
    double tiny_uv[6]={0,0,std::ldexp(1.1,-126),0,std::ldexp(1.1,-127),std::ldexp(1.0,-126)};
    double small_bounds[4]={std::ldexp(1.0,-1074),std::ldexp(1.0,-1022),.8,1.25};
    double small_uv[6]={0,0,std::ldexp(1.0,-1030),0,0,std::ldexp(1.0,-1030)};
    int saved_round=std::fegetround();
#if MB_HAS_SSE
    unsigned saved_control=_mm_getcsr();
#endif
    for(int rounding : {FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO})
        for(int flush=0;flush<2;flush++) {
            failures+=std::fesetround(rounding)!=0;
#if MB_HAS_SSE
            unsigned control=_mm_getcsr();_mm_setcsr((control&~0x8040u)|(flush?0x8040u:0));
#endif
            for(int which=0;which<2;which++) {
                const float *px=which?unit:tiny;const double *pu=which?small_uv:tiny_uv;
                const double *pb=which?small_bounds:bounds;
                MetricBoundsResult got;unsigned expected=0;
                int rc=MetricBounds_classify(px,pu,pb,&got),erc=mb_classify<BSNumber>(px,pu,pb,&expected);
                failures+=rc!=erc||(!rc&&got.violated!=expected);cases++;
            }
        }
    std::fesetround(saved_round);
#if MB_HAS_SSE
    _mm_setcsr(saved_control);
#endif
    std::fprintf(stderr,"[metric bounds test] %zu cases, %zu failures, %zu filtered, %zu exact\n",cases,failures,filtered,exact);
    return failures ? -1 : 0;
}
