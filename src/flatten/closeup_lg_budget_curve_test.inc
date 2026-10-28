static void smbc_test_shear(double *uv,int face,double angle)
{
    double c=cos(angle),s=sin(angle),amount=.002;
    uv[6*face+2]-=s*amount;uv[6*face+3]+=c*amount;
    uv[6*face+4]+=c*amount;uv[6*face+5]+=s*amount;
}

static int smbc_selftest(void)
{
    enum {NF=64,NV=3*NF,N=2*NV};int failures=0,cases=0;
    for(int sample=0;sample<32;sample++) {
        Arena_T arena=Arena_new();float xyz[3*NV];int32_t faces[3*NF],map[NV];
        double reference[N],current[N],trial[N],delta[N],invalid[N];uint8_t masks[2*NF],pending[NF]={0};
        SmMetric metric={0};SmBudget budget={0};SmBudgetCurve curve={0};double angle=.043*sample;
        for(int f=0;f<NF;f++) {
            smb_test_map(reference,f,1,0);smb_test_map(current,f,1,angle);
            masks[f]=1;masks[NF+f]=(uint8_t)(f<32);
            for(int k=0;k<3;k++) {
                int v=3*f+k;faces[3*f+k]=map[v]=v;
                xyz[3*v]=(float)reference[2*v];xyz[3*v+1]=(float)reference[2*v+1];xyz[3*v+2]=0;
            }
        }
        smb_test_map(current,0,.89,angle);smb_test_map(current,1,1/1.1+1e-6,angle);
        smb_test_map(current,40,1/1.1+1e-6,angle);
        int rc=smm_setup(arena,xyz,NV,faces,NF,reference,&metric);
        if(!rc)rc=smb_setup(arena,&metric,2,masks,&budget);
        if(!rc)rc=smbc_setup(arena,NV,&budget,current,4,&curve);
        failures+=rc;cases++;if(rc){Arena_dispose(&arena);continue;}
        SmBound *rows=NULL;size_t count=0;
        rc=smb_bounds(arena,&budget,map,current,&rows,&count);
        if(!rc)rc=smbc_begin(&curve,current,rows,count);
        failures+=rc||count!=1||!curve.active[1]||curve.active[0]||curve.active[40];cases++;
        memcpy(trial,current,sizeof trial);smbc_test_shear(trial,1,angle);smbc_test_shear(trial,40,angle);
        pending[0]=pending[1]=pending[40]=1;size_t added=0;
        rc=smbc_learn(&budget,current,trial,pending,&added);
        failures+=rc||added!=1||curve.state.count!=1||curve.state.planes[0].face!=1;cases++;
        SmCurvePlane saved=curve.state.planes[0];
        failures+=smbc_learn(&budget,current,trial,pending,&added)||added||curve.state.count!=1||
            memcmp(&saved,curve.state.planes,sizeof saved);cases++;
        SmBound *extra=NULL;size_t extra_count=0;rc=smbc_bounds(arena,&budget,map,current,&extra,&extra_count);
        for(int k=0;k<N;k++)delta[k]=trial[k]-current[k];
        failures+=rc||extra_count!=1||extra[0].face!=1||extra[0].upper!=3||extra[0].rhs-smm_dot(extra,delta)<=1e-4;cases++;
        /* A direction cannot survive removal of its current row witness. */
        rc=smbc_begin(&curve,current,NULL,0);extra=NULL;extra_count=0;
        if(!rc)rc=smbc_bounds(arena,&budget,map,current,&extra,&extra_count);
        failures+=rc||extra_count||extra!=NULL||curve.state.count!=1;cases++;
        rc=smbc_begin(&curve,current,rows,count);
        SmBound badrow={0};badrow.upper=3;badrow.face=NF;
        failures+=rc||smbc_begin(&curve,current,&badrow,1)==0||!curve.active[1]||curve.state.count!=1;cases++;
        memcpy(invalid,trial,sizeof invalid);invalid[0]=NAN;
        failures+=smbc_learn(&budget,current,invalid,pending,&added)==0||curve.state.count!=1||
            memcmp(&saved,curve.state.planes,sizeof saved);cases++;
        memcpy(invalid,current,sizeof invalid);invalid[0]+=.01;
        failures+=smbc_learn(&budget,invalid,trial,pending,&added)==0||curve.state.count!=1;cases++;
        /* Real recovery frees the local area's budget and releases the row. */
        smb_test_map(current,0,1,angle);rows=NULL;count=0;
        rc=smb_bounds(arena,&budget,map,current,&rows,&count);
        if(!rc)rc=smbc_begin(&curve,current,rows,count);
        if(!rc)rc=smbc_learn(&budget,current,trial,pending,&added);
        failures+=rc||count||curve.active[1]||curve.state.count||added;cases++;
        /* Two local frontier faces exceed capacity even with no used area. */
        smb_test_map(current,2,1/1.1+1e-6,angle);rows=NULL;count=0;
        rc=smb_bounds(arena,&budget,map,current,&rows,&count);
        if(!rc)rc=smbc_begin(&curve,current,rows,count);
        failures+=rc||count!=2||!curve.active[1]||!curve.active[2]||curve.active[40];cases++;
        memcpy(trial,current,sizeof trial);smbc_test_shear(trial,1,angle);smbc_test_shear(trial,2,angle);
        memset(pending,0,sizeof pending);pending[1]=pending[2]=1;
        size_t capacity=curve.state.capacity;curve.state.capacity=1;
        failures+=smbc_learn(&budget,current,trial,pending,&added)==0||curve.state.count||added;cases++;
        curve.state.capacity=capacity;rc=smbc_learn(&budget,current,trial,pending,&added);
        failures+=rc||added!=2||curve.state.count!=2;cases++;
        pending[1]=pending[2]=2;
        failures+=smbc_learn(&budget,current,trial,pending,&added)||added||curve.state.count!=2;cases++;
        Arena_dispose(&arena);
    }
    fprintf(stderr,"[source budget curvature selftest] %d cases, %d failures\n",cases,failures);
    return failures?1:0;
}
