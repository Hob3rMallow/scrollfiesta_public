static void smc_test_matrix(double a,double b,double lo,double hi,double *j)
{
    double ca=cos(a),sa=sin(a),cb=cos(b),sb=sin(b);
    j[0]=ca*lo*cb-sa*hi*sb;j[1]=-ca*lo*sb-sa*hi*cb;
    j[2]=sa*lo*cb+ca*hi*sb;j[3]=-sa*lo*sb+ca*hi*cb;
}

static double smc_test_inner(const double *current,const double *j)
{
    double length=hypot(current[0]+current[3],current[2]-current[1]);
    double c=(current[0]+current[3])/length,s=(current[2]-current[1])/length;
    double trace=c*(j[0]+j[3])+s*(j[2]-j[1]);
    double difference=c*(j[0]-j[3])+s*(j[2]+j[1]);
    double twice_off=c*(j[1]+j[2])+s*(j[3]-j[0]);
    return .5*(trace-hypot(difference,twice_off));
}

static int smc_selftest(void)
{
    int failures=0,cases=0;double derivative_error=0;
    for(int sample=0;sample<128;sample++) {
        double current[4],trial[4],other[4],sc[2],st[2],so[2];
        smc_test_matrix(.07*sample,-.03*sample,.80001+.001*sample,1.20+.0003*sample,current);
        smc_test_matrix(-.09*sample,.13*sample,.7+.002*sample,1.3+.001*sample,trial);
        smc_test_matrix(.05*sample,.19*sample,.5+.003*sample,1.4+.001*sample,other);
        int rc=smm_singular(current,sc,NULL)||smm_singular(trial,st,NULL)||smm_singular(other,so,NULL);
        failures+=rc||fabs(smc_test_inner(current,current)-sc[0])>2e-12||
            smc_test_inner(current,trial)>st[0]+2e-12;cases++;
        for(int side=0;side<2;side++) {
            double g[4]={0},value=side?st[1]:smc_test_inner(current,trial);
            rc=smc_direction(current,trial,side,g);
            double support=smc_dot4(g,other),expected=side?so[1]:smc_test_inner(current,other);
            failures+=rc||fabs(smc_dot4(g,trial)-value)>2e-12||
                (side?support>expected+2e-12:support<expected-2e-12)||fabs(smc_dot4(g,g)-1)>2e-12;cases++;
            for(int k=0;k<4;k++) {
                double plus[4],minus[4],vp=0,vm=0;
                memcpy(plus,trial,sizeof plus);memcpy(minus,trial,sizeof minus);
                plus[k]+=1e-6;minus[k]-=1e-6;
                if(side) {
                    double sp[2],sm[2];rc|=smm_singular(plus,sp,NULL)||smm_singular(minus,sm,NULL);
                    vp=sp[1];vm=sm[1];
                } else {vp=smc_test_inner(current,plus);vm=smc_test_inner(current,minus);}
                derivative_error=fmax(derivative_error,fabs((vp-vm)/2e-6-g[k]));
            }
            failures+=rc;cases++;
        }
    }
    failures+=derivative_error>2e-8;cases++;
    for(int sample=0;sample<32;sample++) {
        Arena_T arena=Arena_new();SmMetric m={0};SmCurvature curve={0};
        const float xyz[18]={0,0,0,1,0,0,0,1,0,0,0,1,1,0,1,0,1,1};
        const int32_t faces[6]={0,1,2,3,4,5},map[6]={0,1,2,3,4,5};
        const double reference[12]={0,0,1,0,0,1,7,9,8,9,7,10};
        double current[12]={0,0,.800001,0,0,1,7,9,8,9,7,10.249999};
        double trial[12],delta[12],invalid[12];memcpy(trial,current,sizeof trial);
        trial[3]=trial[4]=.002;trial[9]+=.002;trial[10]+=.002;
        double angle=.17*sample,c=cos(angle),s=sin(angle);
        for(int v=0;v<6;v++) {
            double x=current[2*v],y=current[2*v+1];current[2*v]=c*x-s*y;current[2*v+1]=s*x+c*y;
            x=trial[2*v];y=trial[2*v+1];trial[2*v]=c*x-s*y;trial[2*v+1]=s*x+c*y;
        }
        for(int k=0;k<12;k++)delta[k]=trial[k]-current[k];
        uint8_t pending[2]={1,2};size_t added=0;
        int rc=smm_setup(arena,xyz,6,faces,2,reference,&m);
        if(!rc)rc=smc_setup(arena,6,&m,current,4,&curve);
        failures+=rc;cases++;
        if(rc){Arena_dispose(&arena);continue;}
        size_t capacity=curve.capacity;curve.capacity=1;
        failures+=smc_learn(&m,current,trial,pending,&added)==0||curve.count||curve.learned||added;cases++;
        curve.capacity=capacity;
        rc=smc_learn(&m,current,trial,pending,&added);
        failures+=rc||added!=2||curve.count!=2||curve.learned!=2;cases++;
        SmCurvePlane saved[2];memcpy(saved,curve.planes,sizeof saved);
        failures+=smc_learn(&m,current,trial,pending,&added)||added||curve.count!=2||memcmp(saved,curve.planes,sizeof saved);cases++;
        SmBound *rows=NULL;size_t count=0;rc=smc_bounds(arena,&m,map,current,&rows,&count);
        int bad=rc||count!=2;
        if(!bad)for(size_t i=0;i<count;i++)bad|=rows[i].rhs-smm_dot(rows+i,delta)<=1e-4;
        failures+=bad;cases++;
        memcpy(invalid,current,sizeof invalid);invalid[2]=NAN;
        failures+=smc_begin(&curve,invalid)==0||curve.count!=2||memcmp(curve.current,current,sizeof current);cases++;
        memcpy(invalid,trial,sizeof invalid);invalid[2]=NAN;
        failures+=smc_learn(&m,current,invalid,pending,&added)==0||curve.count!=2||memcmp(saved,curve.planes,sizeof saved);cases++;
        memcpy(invalid,current,sizeof invalid);invalid[0]+=.1;
        failures+=smc_learn(&m,invalid,trial,pending,&added)==0||curve.count!=2||memcmp(saved,curve.planes,sizeof saved);cases++;
        failures+=smc_begin(&curve,invalid)||curve.count||curve.resets!=1||memcmp(curve.current,invalid,sizeof invalid);cases++;
        failures+=smc_begin(&curve,current)||curve.resets!=2;cases++;
        m.good[0]=m.good[1]=0;
        failures+=smc_learn(&m,current,trial,pending,&added)||added||curve.count;cases++;
        m.good[0]=m.good[1]=1;
        /* A fresh face is selected by the ordinary feedback path first. */
        current[2]*=1.05;current[11]*=.99;memcpy(trial,current,sizeof trial);
        failures+=smc_begin(&curve,current)||smc_learn(&m,current,trial,pending,&added)||added||curve.count;cases++;
        Arena_dispose(&arena);
    }
    {
        double current[4]={1,0,0,1},trial[4]={1.3,0,0,1.3},g[4]={4,5,6,7};
        failures+=smc_direction(current,trial,1,g)||fabs(smc_dot4(g,trial)-1.3)>1e-12;cases++;
        trial[0]=trial[3]=.7;
        failures+=smc_direction(current,trial,0,g)||fabs(smc_dot4(g,trial)-.7)>1e-12;cases++;
        double saved[4];memcpy(saved,g,sizeof saved);trial[0]=NAN;
        failures+=smc_direction(current,trial,0,g)==0||memcmp(saved,g,sizeof saved);cases++;
        trial[0]=-1;
        failures+=smc_direction(current,trial,1,g)==0||memcmp(saved,g,sizeof saved);cases++;
    }
    fprintf(stderr,"[source curvature selftest] %d cases, %d failures, derivative_error=%.9g\n",cases,failures,derivative_error);
    return failures?1:0;
}
