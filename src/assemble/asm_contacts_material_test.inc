static double ac_test_objective_weight(void *context,size_t a,size_t b,double area)
{(void)a;(void)b;(void)area;return *(const double *)context;}

static int ac_material_policy_control(void)
{
    int fails=fabs(AsmContacts_material_tolerance()-4)>1e-12;
    for(int pose=0;pose<2;pose++)for(int width=0;width<4;width++)for(int refined=0;refined<2;refined++){
        Arena_T arena=Arena_new();int n=refined?20:1;
        size_t vertices=(size_t)(n+1)*(n+1),faces=(size_t)2*n*n,nv=2*vertices,nf=2*faces;
        double *uv=ARENA_ALLOC(arena,2*nv*sizeof *uv),*before=ARENA_ALLOC(arena,2*nv*sizeof *before);
        int32_t *tri=ARENA_ALLOC(arena,3*nf*sizeof *tri);
        size_t *external=ARENA_ALLOC(arena,faces*sizeof *external);
        uint8_t *exclude=ARENA_CALLOC(arena,nf,1);
        double depth=width==0?3.5:width==1?4:width==2?4.0001:6,angle=pose?.37:0,cs=cos(angle),sn=sin(angle);
        for(int side=0;side<2;side++){
            for(int j=0;j<=n;j++)for(int i=0;i<=n;i++){
                size_t v=(size_t)side*vertices+(size_t)j*(n+1)+i;double x=20.*i/n+side*(20-depth),y=20.*j/n;
                uv[2*v]=100000+cs*x-sn*y;uv[2*v+1]=-27000+sn*x+cs*y;
            }
            for(int j=0;j<n;j++)for(int i=0;i<n;i++){
                size_t t=(size_t)side*faces+2*((size_t)j*n+i);
                int32_t v=(int32_t)((size_t)side*vertices+(size_t)j*(n+1)+i);
                tri[3*t]=v;tri[3*t+1]=v+1;tri[3*t+2]=v+n+2;
                tri[3*t+3]=v;tri[3*t+4]=v+n+2;tri[3*t+5]=v+n+1;
            }
        }
        for(size_t i=0;i<faces;i++){external[i]=faces+i;exclude[faces+i]=1;}
        AsmContacts *q=AsmContacts_new(arena,nv,tri,nf,uv,AsmContacts_material_tolerance());
        AsmContactStats serial={0},threaded={0},probe={0};int bad=0;
        if(!q){fails++;Arena_dispose(&arena);continue;}
        if(q->material->components!=2 || q->material->edges!=(size_t)8*n || q->material->strict[0] || q->material->strict[1])bad=1;
        if(AsmContacts_measure(q,uv,NULL,NULL,NULL,&serial) || !serial.complete || !!serial.pairs!=(width>=2))bad=1;
        AsmContacts_threads(q,4);
        if(AsmContacts_measure(q,uv,NULL,NULL,NULL,&threaded) || !threaded.complete ||
           threaded.pairs!=serial.pairs || threaded.area!=serial.area)bad=1;
        if(AsmContacts_probe_mesh_prepared(q,q,uv,faces,external,exclude,NULL,NULL,&probe) ||
           !probe.complete || probe.pairs!=serial.pairs || fabs(probe.area-serial.area)>1e-7)bad=1;
        /* Translating one square along its local x axis changes the true
         * intersection area at exactly -20 vox^2/vox at every subdivision.
         * The material-filtered gradient can cancel to zero on fine meshes;
         * the optimization objective must retain the true boundary force. */
        {
            double *gradient=ARENA_ALLOC(arena,2*nv*sizeof(double)),*parallel=ARENA_ALLOC(arena,2*nv*sizeof(double));
            double raw=0,raw_parallel=0,force=0;AsmContactStats objective={0};
            AsmContacts_threads(q,1);
            if(AsmContacts_objective(q,uv,gradient,NULL,NULL,&objective,&raw) || objective.pairs!=serial.pairs || objective.area!=serial.area)bad=1;
            for(size_t v=vertices;v<nv;v++)force+=cs*gradient[2*v]+sn*gradient[2*v+1];
            if(fabs(raw-20*depth)>1e-6 || fabs(force+20)>1e-6)bad=1;
            if(refined){
                size_t left=vertices+(size_t)(n/2)*(n+1),inside=left+1;
                if(fabs(cs*gradient[2*left]+sn*gradient[2*left+1]+20./n)>1e-6 ||
                   hypot(gradient[2*inside],gradient[2*inside+1])>1e-6)bad=1;
            }
            AsmContacts_threads(q,4);
            if(AsmContacts_objective(q,uv,parallel,NULL,NULL,&objective,&raw_parallel) || raw!=raw_parallel ||
               memcmp(gradient,parallel,2*nv*sizeof(double)))bad=1;
            /* Soft weighting changes force and energy, never the material
             * certificate, including boundary forces on subdivided charts. */
            for(int zero=0;zero<2;zero++){
                double weight=zero?0:.5,weighted=0;
                if(AsmContacts_weighted_objective(q,uv,parallel,NULL,NULL,ac_test_objective_weight,&weight,&objective,&weighted) ||
                   !objective.complete || objective.pairs!=serial.pairs || objective.area!=serial.area ||
                   fabs(weighted-weight*raw)>1e-7)bad=1;
                for(size_t v=0;v<2*nv;v++)if(fabs(parallel[v]-weight*gradient[v])>1e-7)bad=1;
            }
            double invalid=-1;
            if(!AsmContacts_weighted_objective(q,uv,NULL,NULL,NULL,ac_test_objective_weight,&invalid,&objective,&raw_parallel) || objective.complete)bad=1;
            if(AsmContacts_objective(q,uv,parallel,NULL,NULL,&objective,&raw_parallel) || raw!=raw_parallel ||
               memcmp(gradient,parallel,2*nv*sizeof(double)))bad=1;
            if(refined && width==3 && !pose){
                double retained_force=0;
                if(AsmContacts_measure(q,uv,gradient,NULL,NULL,&objective))bad=1;
                for(size_t v=vertices;v<nv;v++)retained_force+=gradient[2*v];
                fprintf(stderr,"  material translation derivative: full %.9g (analytic -20), filtered %.9g\n",force,retained_force);
            }
        }
        /* A single selected interior triangle still belongs to its complete
         * original patch. Fine triangles cannot supply an invented boundary. */
        if(refined && width==3){
            size_t selected=faces+2*((size_t)(n/2)*n+2);
            if(AsmContacts_probe_mesh_prepared(q,q,uv,1,&selected,exclude,NULL,NULL,&probe) || !probe.complete || !probe.pairs)bad=1;
        }
        memcpy(before,uv,2*nv*sizeof *before);
        double old_area=AsmContacts_pair_faces(q,before,faces/2+2*(n-1),faces+faces/2,NULL,NULL);
        for(size_t v=vertices;v<nv;v++){uv[2*v]+=100*cs;uv[2*v+1]+=100*sn;}
        if(AsmContacts_refit_changed(q,uv,exclude) || AsmContacts_measure(q,uv,NULL,NULL,NULL,&probe) || probe.pairs)bad=1;
        if(AsmContacts_pair_faces(q,before,faces/2+2*(n-1),faces+faces/2,NULL,NULL)!=old_area)bad=1;
        /* Rebuilding the initial geometry after a refit must restore both
         * indices even when its partition can be reused from the build. */
        memcpy(uv,before,2*nv*sizeof *uv);
        if(AsmContacts_rebuild(q,uv) || AsmContacts_probe_mesh_prepared(q,q,uv,faces,external,exclude,NULL,NULL,&probe) ||
           probe.pairs!=serial.pairs || fabs(probe.area-serial.area)>1e-7)bad=1;
        if(bad){fails++;fprintf(stderr,"  material width %.7g pose %d subdivision %d: %zu pairs, probe %zu\n",depth,pose,n,serial.pairs,probe.pairs);}
        Arena_dispose(&arena);
    }
    fprintf(stderr,"  fixed 4-voxel original-material query policy: %s (%d failures)\n",fails?"FAIL":"ok",fails);
    return fails;
}
