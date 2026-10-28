    {
        /* Preserve the original three-corner arithmetic, including signed zero,
         * cancellation, wide scales and nonfinite classification. Inputs and
         * result never overlap at any internal Jacobian call site. */
        SRFace face={0}; const int32_t tri[3]={1,3,4};
        double uv[10],actual[4],expected[4];
        uint32_t random=0x6a19f277u;
        const int power[7]={0,20,-20,400,-400,700,-700};
        size_t cases=0,mismatch=0,finite=0,nonfinite=0;
        for (int scale=0;scale<49;scale++) for (int sample=0;sample<1024;sample++) {
            for (int k=0;k<3;k++) {
                random=1664525u*random+1013904223u;
                face.gx[k]=ldexp((double)(random&65535)/32768-1,power[scale/7]);
                random=1664525u*random+1013904223u;
                face.gy[k]=ldexp((double)(random&65535)/32768-1,power[scale/7]);
            }
            if (sample&1) {face.gx[0]=-face.gx[1]-face.gx[2];face.gy[0]=-face.gy[1]-face.gy[2];}
            for (int k=0;k<10;k++) {
                random=1664525u*random+1013904223u;
                uv[k]=ldexp((double)(random&65535)/32768-1,power[scale%7]);
                if (sample<2) uv[k]=(sample^(k&1))?-0.0:0.0;
            }
            if (sample==2) uv[2]=NAN;
            if (sample==3) uv[2]=INFINITY;
            if (sample==4) uv[2]=-INFINITY;
            sr_face_jacobian(&face,tri,uv,actual);
            sr_jacobian_original(&face,tri,uv,expected);
            for (int k=0;k<4;k++) {
                mismatch+=!(isnan(actual[k]) && isnan(expected[k])) && memcmp(actual+k,expected+k,sizeof *actual)!=0;
                finite+=isfinite(expected[k])!=0;nonfinite+=isfinite(expected[k])==0;
            }
            cases++;
        }
        if (mismatch || !finite || !nonfinite) failures++;
        fprintf(stderr,"[slim selftest] local Jacobian %s: %zu cases, %zu finite values, %zu nonfinite values, %zu mismatches\n",
            mismatch || !finite || !nonfinite ? "FAIL":"ok",cases,finite,nonfinite,mismatch);
    }
