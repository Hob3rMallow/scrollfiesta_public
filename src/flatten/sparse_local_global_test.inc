static int slg_expect(int condition,int *cases,int *failures)
{
    (*cases)++;if(!condition)(*failures)++;return condition;
}

int Sparse_local_global_selftest(void)
{
    Arena_T arena=Arena_new();int cases=0,failures=0;
    for(int t=0;t<24;t++) {
        Arena_Mark mark=Arena_save(arena);
        double d0=pow(10,(double)(t%7-3)),d1=1/d0,rscale=pow(10,(double)(t%9-4));
        int row[3]={0,1,1},col[3]={0,0,1};
        double h[3]={2*d0*d0,-.5*d0*d1,3*d1*d1},b[2]={-.45*d0,.6*d1};
        double known[2]={.2/d0,.3/d1},dual[4]={.7/rscale,.2*rscale,0,0};
        SparseLocalRow c[4]={0};
        c[0].value[0]=d0*rscale;c[0].lower=.2*rscale;
        c[1].index[0]=1;c[1].value[0]=d1/rscale;c[1].lower=.3/rscale;
        c[2].index[1]=1;c[2].value[0]=d0;c[2].value[1]=d1;c[2].lower=.25;
        c[3].value[0]=.25*d0;c[3].value[2]=.75*d0;c[3].lower=.1;
        for(int mode=0;mode<3;mode++) {
            double start[2]={known[0]*(mode==1?1:.9),known[1]*(mode==1?1:.9)};
            double out[2]={987,654},lambda[4]={0};SparseLocalGlobalStats s={0};
            int rc=Sparse_local_global(arena,2,3,row,col,h,b,4,c,mode?start:NULL,mode?dual:NULL,10000,out,lambda,&s);
            slg_expect(!rc&&s.converged&&fabs(out[0]*d0-.2)<1e-7&&fabs(out[1]*d1-.3)<1e-7,&cases,&failures);
            slg_expect(mode!=1||s.iterations==0,&cases,&failures);
            slg_expect(s.factors<=s.solves&&s.stationarity_residual<=1e-9&&s.primal_residual<=1e-9&&s.complementarity<=1e-10,&cases,&failures);
        }
        double out[2]={987,654},lambda[4]={7,8,9,10};SparseLocalGlobalStats s={0};
        slg_expect(Sparse_local_global(arena,2,3,row,col,h,b,4,c,NULL,NULL,0,out,lambda,&s)==1&&out[0]==987&&lambda[0]==7&&!s.converged,&cases,&failures);
        c[0].lower=NAN;
        slg_expect(Sparse_local_global(arena,2,3,row,col,h,b,4,c,NULL,NULL,20,out,lambda,&s)==-1&&out[0]==987&&lambda[0]==7,&cases,&failures);
        Arena_restore(arena,mark);
    }
    { /* No rows: the ordinary global least-squares solution. */
        int row[1]={0},col[1]={0};double h[1]={2},b[1]={6},out[1]={91};SparseLocalGlobalStats s={0};
        slg_expect(!Sparse_local_global(arena,1,1,row,col,h,b,0,NULL,NULL,NULL,3,out,NULL,&s)&&fabs(out[0]-3)<1e-12,&cases,&failures);
    }
    { /* Nearly parallel binding planes: multiplier progress must not wait
       * for excessive accuracy on each intermediate penalty problem. */
        for(int t=0;t<4;t++) {
            Arena_Mark mark=Arena_save(arena);
            int row[2]={0,1},col[2]={0,1};double h[2]={2,3},b[2]={.4,-1};
            double epsilon=t%2?1e-3:1e-2,rscale=t<2?1:100;
            SparseLocalRow c[2]={0};
            c[0].index[1]=1;c[0].value[0]=rscale;c[0].value[1]=rscale*epsilon;
            c[1].index[1]=1;c[1].value[0]=-1/rscale;c[1].value[1]=epsilon/rscale;
            double initial[2]={0,.01},dual[2]={(.5/epsilon-.2)/rscale,(.5/epsilon+.2)*rscale};
            double out[2]={987,654},lambda[2]={0};SparseLocalGlobalStats s={0};
            int rc=Sparse_local_global(arena,2,2,row,col,h,b,2,c,initial,dual,10000,out,lambda,&s);
            slg_expect(!rc&&s.converged&&fabs(out[0])<1e-6&&fabs(out[1])<1e-6,&cases,&failures);
            slg_expect(s.primal_residual<=1e-9&&s.stationarity_residual<=1e-9&&s.complementarity<=1e-10,&cases,&failures);
            Arena_restore(arena,mark);
        }
    }
    { /* Inconsistent local copies must not publish an approximate result. */
        int row[1]={0},col[1]={0};double h[1]={2},b[1]={0},out[1]={91},lambda[2]={4,5};
        SparseLocalRow c[2]={0};c[0].value[0]=1;c[0].lower=1;c[1].value[0]=-1;c[1].lower=1;
        SparseLocalGlobalStats s={0};
        slg_expect(Sparse_local_global(arena,1,1,row,col,h,b,2,c,NULL,NULL,100,out,lambda,&s)!=0&&out[0]==91&&lambda[0]==4&&!s.converged,&cases,&failures);
        slg_expect(s.iterations<=100,&cases,&failures);
    }
    Arena_dispose(&arena);
    fprintf(stderr,"local/global least squares: %d cases, %d failures\n",cases,failures);
    return failures?1:0;
}
