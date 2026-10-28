    {
        /* Differential controls retain the independent general path predicate.
         * Include two-root crossings, translated/thin frames, reversals,
         * invalid data and refreshes of the same cache storage. */
        const int32_t f[3] = {0,1,2};
        const double base[6] = {0,0,1,0,0,1};
        const double two_roots[6] = {0,0,-2,0,0,-3};
        double uv[6], target[6]; SRGrowthReference reference;
        uint32_t rng = 0x18c9542bu;
        size_t cases = 0, mismatch = 0, safe = 0, unsafe = 0;
        sr_growth_reference_set(&reference,base,f);
        if (sr_growth_reference_safe(&reference,two_roots,f) ||
            UvGuard_interval(base,two_roots,f,1)) mismatch++;
        for (int frame = 0; frame < 12; frame++) {
            double scale = frame < 4 ? 1 : frame < 8 ? 1e-6 : 1e6;
            double shift = (frame & 1) ? 1e8 : 0;
            for (size_t k = 0; k < 6; k++) uv[k] = shift+scale*base[k];
            if (frame%4 == 2) uv[5] = shift+scale*1e-10;
            if (frame%4 == 3) uv[5] = shift-scale;
            sr_growth_reference_set(&reference,uv,f);
            for (int sample = 0; sample < 8192; sample++) {
                for (size_t k = 0; k < 6; k++) {
                    rng = 1664525u*rng+1013904223u;
                    target[k] = uv[k]+scale*((double)(rng & 65535)/8192-4);
                }
                if (!sample) memcpy(target,uv,sizeof uv);
                int expected = UvGuard_interval(uv,target,f,1);
                int actual = sr_growth_reference_safe(&reference,target,f);
                mismatch += actual != expected; cases++;
                safe += expected != 0; unsafe += expected == 0;
            }
        }
        for (int bad = 0; bad < 3; bad++) {
            memcpy(uv,base,sizeof uv); memcpy(target,base,sizeof target);
            uv[2] = bad == 0 ? NAN : bad == 1 ? INFINITY : -INFINITY;
            sr_growth_reference_set(&reference,uv,f);
            mismatch += sr_growth_reference_safe(&reference,target,f) != UvGuard_interval(uv,target,f,1);
            cases++;
            memcpy(uv,base,sizeof uv); target[2] = bad == 0 ? NAN : bad == 1 ? INFINITY : -INFINITY;
            sr_growth_reference_set(&reference,uv,f);
            mismatch += sr_growth_reference_safe(&reference,target,f) != UvGuard_interval(uv,target,f,1);
            cases++;
        }
        sr_growth_reference_set(&reference,base,f);
        mismatch += !sr_growth_reference_safe(&reference,base,f);
        if (mismatch || !safe || !unsafe) failures++;
        fprintf(stderr,"[slim selftest] growth reference %s: %zu cases, %zu safe, %zu unsafe, %zu mismatches\n",
            mismatch || !safe || !unsafe ? "FAIL" : "ok",cases,safe,unsafe,mismatch);
    }
