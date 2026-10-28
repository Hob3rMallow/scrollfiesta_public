static void smbu_test_map(double *uv,int face,double low,double high,double angle)
{
    double c=cos(angle),s=sin(angle),x=3*face;
    const double raw[6]={x,0,x+low,0,x,high};
    for(int k=0;k<3;k++){uv[6*face+2*k]=c*raw[2*k]-s*raw[2*k+1];uv[6*face+2*k+1]=s*raw[2*k]+c*raw[2*k+1];}
}

static int smbu_selftest(void)
{
    enum {NF=64,NV=3*NF,N=2*NV};int failures=0,cases=0;double derivative_error=0;
    for(int sample=0;sample<32;sample++) {
        Arena_T arena=Arena_new();float xyz[3*NV];int32_t faces[3*NF],map[NV];
        double reference[N],current[N],trial[N],delta[N],invalid[N];uint8_t masks[2*NF],pending[NF]={0};
        SmMetric metric={0};SmBudget budget={0};SmBudgetCurve curve={0};double angle=.043*sample;
        for(int f=0;f<NF;f++) {
            smbu_test_map(reference,f,1,1,0);smbu_test_map(current,f,1,1,angle);
            masks[f]=1;masks[NF+f]=(uint8_t)(f<32);
            for(int k=0;k<3;k++) {
                int v=3*f+k;faces[3*f+k]=map[v]=v;
                xyz[3*v]=(float)reference[2*v];xyz[3*v+1]=(float)reference[2*v+1];xyz[3*v+2]=0;
            }
        }
        smbu_test_map(current,0,.89,1,angle);
        smbu_test_map(current,1,1,1.1-1e-6,angle);smbu_test_map(current,40,1,1.1-1e-6,angle);
        int rc=smm_setup(arena,xyz,NV,faces,NF,reference,&metric);
        if(!rc)rc=smb_setup(arena,&metric,2,masks,&budget);
        if(!rc)rc=smbc_setup(arena,NV,&budget,current,4,&curve);
        failures+=rc;cases++;if(rc){Arena_dispose(&arena);continue;}curve.use_upper=1;
        SmBound *rows=NULL;size_t count=0;
        rc=smb_bounds(arena,&budget,map,current,&rows,&count);
        if(!rc)rc=smbc_begin(&curve,current,rows,count);
        failures+=rc||count!=1||curve.active[1]!=2||curve.active[0]||curve.active[40];cases++;
        memcpy(trial,current,sizeof trial);smbc_test_shear(trial,1,angle);smbc_test_shear(trial,40,angle);
        pending[0]=pending[1]=pending[40]=2;size_t added=0;
        rc=smbc_learn(&budget,current,trial,pending,&added);
        failures+=rc||added!=1||curve.state.count!=1||curve.state.planes[0].face!=1||curve.state.planes[0].side!=1;cases++;
        SmCurvePlane saved=curve.state.planes[0];
        failures+=smbc_learn(&budget,current,trial,pending,&added)||added||curve.state.count!=1||memcmp(&saved,curve.state.planes,sizeof saved);cases++;
        SmBound *extra=NULL;size_t extra_count=0;rc=smbc_bounds(arena,&budget,map,current,&extra,&extra_count);
        for(int k=0;k<N;k++)delta[k]=trial[k]-current[k];
        failures+=rc||extra_count!=1||extra[0].face!=1||extra[0].upper!=4||extra[0].rhs-smm_dot(extra,delta)<=1e-4;cases++;
        /* The cut is the spectral norm derivative at the refused trial. */
        double jt[4],st[2];smm_jacobian(metric.rest+1,faces+3,trial,jt);rc=smm_singular(jt,st,NULL);
        for(int k=0;k<4;k++) {
            double plus[4],minus[4],sp[2],sm[2];memcpy(plus,jt,sizeof plus);memcpy(minus,jt,sizeof minus);
            plus[k]+=1e-6;minus[k]-=1e-6;
            rc|=smm_singular(plus,sp,NULL)||smm_singular(minus,sm,NULL);
            derivative_error=fmax(derivative_error,fabs((sp[1]-sm[1])/2e-6-saved.gradient[k]));
        }
        failures+=rc||fabs(smc_dot4(saved.gradient,jt)-st[1])>1e-12;cases++;
        /* Current good maps are on the admitted side before the safety skin. */
        double jc[4];smm_jacobian(metric.rest+1,faces+3,current,jc);
        failures+=smc_dot4(saved.gradient,jc)>1.1+1e-12;cases++;
        /* An opposite-side witness must not reactivate a retained upper cut. */
        SmBound opposite={0};opposite.face=1;opposite.upper=3;
        rc=smbc_begin(&curve,current,&opposite,1);extra=NULL;extra_count=0;
        if(!rc)rc=smbc_bounds(arena,&budget,map,current,&extra,&extra_count);
        failures+=rc||extra_count||extra!=NULL||curve.state.count!=1;cases++;
        rc=smbc_begin(&curve,current,rows,count);
        SmBound badrow={0};badrow.upper=4;badrow.face=NF;
        failures+=rc||smbc_begin(&curve,current,&badrow,1)==0||curve.active[1]!=2||curve.state.count!=1;cases++;
        memcpy(invalid,trial,sizeof invalid);invalid[0]=NAN;
        failures+=smbc_learn(&budget,current,invalid,pending,&added)==0||curve.state.count!=1||memcmp(&saved,curve.state.planes,sizeof saved);cases++;
        memcpy(invalid,current,sizeof invalid);invalid[0]+=.01;
        failures+=smbc_learn(&budget,invalid,trial,pending,&added)==0||curve.state.count!=1;cases++;
        /* Recovery releases the group capacity and all upper witnesses. */
        smbu_test_map(current,0,1,1,angle);rows=NULL;count=0;
        rc=smb_bounds(arena,&budget,map,current,&rows,&count);
        if(!rc)rc=smbc_begin(&curve,current,rows,count);
        if(!rc)rc=smbc_learn(&budget,current,trial,pending,&added);
        failures+=rc||count||curve.active[1]||curve.state.count||added;cases++;
        /* A single face can supply two independent side directions. */
        smbu_test_map(current,2,1/1.1+1e-6,1.1-1e-6,angle);rows=NULL;count=0;
        rc=smb_bounds(arena,&budget,map,current,&rows,&count);
        if(!rc)rc=smbc_begin(&curve,current,rows,count);
        failures+=rc||count!=3||curve.active[1]!=2||curve.active[2]!=3||curve.active[40];cases++;
        memcpy(trial,current,sizeof trial);smbu_test_map(trial,1,1,1.11,angle);smbu_test_map(trial,2,.90,1.11,angle);
        memset(pending,0,sizeof pending);pending[1]=2;pending[2]=3;
        failures+=smbc_learn(&budget,current,trial,pending,&added)||added||curve.state.count;cases++;
        memcpy(trial,current,sizeof trial);smbc_test_shear(trial,1,angle);smbc_test_shear(trial,2,angle);
        size_t capacity=curve.state.capacity;curve.state.capacity=2;
        failures+=smbc_learn(&budget,current,trial,pending,&added)==0||curve.state.count||added;cases++;
        curve.state.capacity=capacity;rc=smbc_learn(&budget,current,trial,pending,&added);
        failures+=rc||added!=3||curve.state.count!=3||curve.state.planes[1].side!=0||curve.state.planes[2].side!=1;cases++;
        extra=NULL;extra_count=0;rc=smbc_bounds(arena,&budget,map,current,&extra,&extra_count);
        for(int k=0;k<N;k++)delta[k]=trial[k]-current[k];
        int bad=rc||extra_count!=3;
        if(!bad)for(size_t i=0;i<extra_count;i++)bad|=extra[i].rhs-smm_dot(extra+i,delta)<=1e-4;
        failures+=bad;cases++;
        /* A zero pending mask cannot create new cuts. */
        memset(pending,0,sizeof pending);
        failures+=smbc_learn(&budget,current,trial,pending,&added)||added||curve.state.count!=3;cases++;
        smbu_test_map(current,2,1,1,angle);rows=NULL;count=0;
        rc=smb_bounds(arena,&budget,map,current,&rows,&count);
        if(!rc)rc=smbc_begin(&curve,current,rows,count);
        failures+=rc||count||curve.state.count;cases++;
        Arena_dispose(&arena);
    }
    failures+=derivative_error>2e-8;cases++;
    fprintf(stderr,"[source upper budget curvature selftest] %d cases, %d failures, derivative_error=%.9g\n",cases,failures,derivative_error);
    return failures?1:0;
}
