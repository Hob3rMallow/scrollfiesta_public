/* A flat narrow strip admits a gentle planar bend with bounded strain. The
 * varying measured directions need that bend: a common rotation cannot align
 * them. This exercises the public whole-sheet transaction and its loss audit. */
static int ic_local_global_integration_test(void)
{
#ifdef VES_CLOSEUP_LOCAL_GLOBAL
    enum {NX=41,NY=3,NV=NX*NY,NF=2*(NX-1)*(NY-1)};
    float xyz[3*NV];double uv[2*NV],out[2*NV];int32_t faces[3*NF];uint8_t attr[NV];
    SlimRefineDirection arms[NF];memset(arms,0,sizeof arms);
    for(int y=0;y<NY;y++)for(int x=0;x<NX;x++) {
        int v=y*NX+x;xyz[3*v]=(float)x;xyz[3*v+1]=(float)y;xyz[3*v+2]=0;
        uv[2*v]=x;uv[2*v+1]=y;attr[v]=1;
    }
    size_t at=0;
    for(int y=0;y<NY-1;y++)for(int x=0;x<NX-1;x++) {
        int a=y*NX+x,b=a+1,c=a+NX,d=c+1;
        int32_t pair[6]={a,b,d,a,d,c};memcpy(faces+3*at,pair,sizeof pair);
        for(int t=0;t<2;t++) {
            double center=(x+(t==0?2.0/3:1.0/3));
            double theta=.4*sin(6.28318530717958647692*center/40);
            arms[at+t].h[0]=cos(theta);arms[at+t].h[1]=-sin(theta);
            arms[at+t].v[0]=sin(theta);arms[at+t].v[1]=cos(theta);arms[at+t].weight=1;
        }
        at+=2;
    }
    MeshBinData mesh={0};mesh.verts=xyz;mesh.faces=faces;mesh.uv64=uv;mesh.nv=NV;mesh.nf=NF;
    CloseupAlignOpts options;CloseupAlign_defaults(&options);options.slim.linear_solve=ic_slim_linear;
    options.slim.verbose=0;SlimRefineStats stats={0};IcPhysicalAlignStats physical={0};
    int rc=ic_reparameterize_physical(&mesh,arms,attr,&options,out,&stats,&physical,NULL);
    double rms=sqrt(stats.direction_angle_after)*57.2957795130823208768;
    int failures=ic_align_check(!rc&&rms<10&&stats.energy_after>stats.energy_before+1e-5&&
        stats.objective_after<stats.objective_before&&physical.bad10_area_after<=.05*physical.area+1e-10&&
        !physical.protected_bad25_faces_after&&!stats.flips_after,
        "public local/global accepts bounded-strain nonrigid angular improvement under the correct loss");
    failures+=ic_align_check(fabs(ic_align_bad10_limit(100,4.9,80,1.0)-1.1)<1e-12&&
        fabs(ic_align_bad10_limit(100,4.9,20,3.9)-1.0)<1e-12&&
        ic_align_bad10_limit(100,8,20,1)==0&&ic_align_bad10_limit(0,0,0,0)==0,
        "component metric budgets reserve existing bad area elsewhere on the sheet");
    return failures;
#else
    return 0;
#endif
}
