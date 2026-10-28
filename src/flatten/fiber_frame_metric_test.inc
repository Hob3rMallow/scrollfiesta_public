/* The same physical plane in several sheared/stretched canvas coordinates.
 * Expected output coordinates are dot products with known WORLD fibre axes,
 * independent of the solver's target-gradient construction. */
static int ff_metric_selftest(void)
{
    enum { H=13, W=21, N=H*W, CELL=8 };
    const double cases[5][4]={{1,0,1,0},{1,0,1,.23},{1,.4,1,0},{2,.4,.5,.23},{2,0,.5,.23}};
    int fails=0;
    for(size_t test=0;test<5;test++){
        Arena_T arena=Arena_new();FiberFrame ff={0};FiberFrameStats st={0};
        double a=cases[test][0],b=cases[test][1],d=cases[test][2],angle=cases[test][3];
        double co=cos(angle),si=sin(angle),hv=si/d,hu=(co-b*hv)/a;
        double error=0,lo=0,hi=0;
        ff.gh=H;ff.gw=W;ff.cell=CELL;
        ff.theta=ARENA_ALLOC(arena,N*sizeof *ff.theta);ff.weight=ARENA_ALLOC(arena,N*sizeof *ff.weight);
        ff.real=ARENA_ALLOC(arena,N);ff.closed=ARENA_ALLOC(arena,N);
        ff.pu=ARENA_CALLOC(arena,3*N,sizeof *ff.pu);ff.pv=ARENA_CALLOC(arena,3*N,sizeof *ff.pv);
        ff.U=ARENA_ALLOC(arena,N*sizeof *ff.U);ff.V=ARENA_ALLOC(arena,N*sizeof *ff.V);
        for(size_t i=0;i<N;i++){
            ff.theta[i]=(float)atan2(hv,hu);ff.weight[i]=1;ff.real[i]=1;
            /* E1=(.6,0,.8), E2=(0,1,0): a tilted world plane. */
            ff.pu[3*i]=.6*a;ff.pu[3*i+2]=.8*a;
            ff.pv[3*i]=.6*b;ff.pv[3*i+1]=d;ff.pv[3*i+2]=.8*b;
        }
        /* Both perpendicular physical fibre layers must give the same
         * fourfold phase, even though their canvas directions are oblique.
         * A normal component must not change the in-plane observation. */
        double e1[3],e2[3],metric[3];
        int frame_bad=ff_physical_frame(ff.pu,ff.pv,e1,e2,metric);
        if(!frame_bad){
            for(int arm=0;arm<2;arm++){
                double x=arm?-si:co,y=arm?co:si;
                double direction[3]={.6*x-.8*.2,y,.8*x+.6*.2};
                double phase=ff_physical_angle(direction,e1,e2);
                if(fabs(cos(4*phase)-cos(4*angle))>1e-10 || fabs(sin(4*phase)-sin(4*angle))>1e-10)frame_bad=1;
            }
            if(fabs(remainder(ff_canvas_angle(angle,metric)-atan2(hv,hu),2*M_PI))>1e-10)frame_bad=1;
        }
        int rc=FiberFrame_solve(arena,&ff,&st);
        if(!rc){
            for(int y=0;y<H;y++)for(int x=0;x<W;x++){
                size_t i=(size_t)y*W+(size_t)x;
                double world_x=CELL*(a*x+b*y),world_y=CELL*d*y;
                double eu=co*world_x+si*world_y,ev=-si*world_x+co*world_y;
                error=fmax(error,hypot(ff.U[i]-ff.U[0]-eu,ff.V[i]-ff.V[0]-ev));
            }
            size_t i=(H/2)*W+W/2;
            double ux=(ff.U[i+1]-ff.U[i])/CELL,uy=(ff.U[i+W]-ff.U[i])/CELL;
            double vx=(ff.V[i+1]-ff.V[i])/CELL,vy=(ff.V[i+W]-ff.V[i])/CELL;
            double A=ux/a,C=vx/a,B=(uy-b*A)/d,D=(vy-b*C)/d;
            hi=.5*(hypot(A+D,C-B)+hypot(A-D,C+B));lo=fabs(A*D-B*C)/hi;
        }
        int bad=frame_bad || rc || error>1e-3 || fabs(lo-1)>1e-5 || fabs(hi-1)>1e-5;
        fails+=bad;
        fprintf(stderr,"[selftest] fiber physical coordinates %zu: max position error %.9g, physical sigma %.9g..%.9g; %s\n",
                test,error,lo,hi,bad?"FAIL":"PASS");
        Arena_dispose(&arena);
    }
    fprintf(stderr,"[selftest] fiber metric coordinate test: %s (%d failures)\n",fails?"FAIL":"PASS",fails);
    return fails;
}
