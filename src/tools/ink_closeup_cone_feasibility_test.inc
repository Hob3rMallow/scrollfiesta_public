#include "ink_closeup_cone_feasibility_fixture.inc"

static int ic_cone_feasibility_selftest(void)
{
    int failures=ic_cone_search_selftest();
    for (int fixture=0;fixture<4;fixture++) {
        Arena_T arena=Arena_new();
        float xyz[27]={0},saved_xyz[27];double uv[18]={0},saved_uv[18];
        int32_t faces[24],saved_faces[24];uint8_t terminals[9]={2};
        SlimRefineDirection directions[8],saved_directions[8];
        for (int k=0;k<8;k++) {
            double angle=6.283185307179586*k/8;
            xyz[3*(k+1)]=(float)cos(angle);xyz[3*(k+1)+1]=(float)sin(angle);
            xyz[3*(k+1)+2]=(k&1?1.f:-1.f)*(fixture==1?.4f:.1f);
            uv[2*(k+1)]=xyz[3*(k+1)];uv[2*(k+1)+1]=xyz[3*(k+1)+1];
            faces[3*k]=0;faces[3*k+1]=k+1;faces[3*k+2]=(k+1)%8+1;
            double direction=fixture==2 && (k&1)?.7853981633974483:0;
            SlimRefineDirection d={{cos(direction),sin(direction),0},{-sin(direction),cos(direction),0},1};
            directions[k]=d;
        }
        if (fixture==3) directions[3].weight=0;
        memcpy(saved_xyz,xyz,sizeof xyz);memcpy(saved_uv,uv,sizeof uv);
        memcpy(saved_faces,faces,sizeof faces);memcpy(saved_directions,directions,sizeof directions);
        MeshBinData source={0};source.verts=xyz;source.uv64=uv;source.faces=faces;source.nv=9;source.nf=8;
        IcConeFanContext check={0};check.arena=arena;check.mesh=&source;check.directions=directions;
        SeamCutResult legacy={0},cut={0};
        int bad=SeamCut_to_boundary(arena,xyz,9,faces,8,terminals,&legacy);
        bad|=SeamCut_to_boundary_checked(arena,xyz,9,faces,8,terminals,ic_cone_fan_feasible,&check,&cut);
        if (cut.nf!=8 || legacy.seam_edges!=2 || cut.seam_edges!=(fixture?2:1) ||
            check.reasons[IC_FAN_PASS]!=(fixture?0u:1u)) bad=1;
        if (fixture==1 && !check.reasons[IC_FAN_METRIC] && !check.reasons[IC_FAN_GUARD]) bad=1;
        if (fixture==2 && !check.reasons[IC_FAN_METRIC] && !check.reasons[IC_FAN_DIRECTION] && !check.reasons[IC_FAN_GUARD]) bad=1;
        if (fixture==3 && !check.reasons[IC_FAN_SUPPORT]) bad=1;
        if (!bad) {
            for (size_t k=0;k<3*cut.nf;k++) if (cut.vmap[cut.faces[k]]!=faces[k]) bad=1;
            for (size_t v=0;v<cut.nv;v++) if (memcmp(cut.verts+3*v,xyz+3*(size_t)cut.vmap[v],3*sizeof(float))) bad=1;
        }
        if (memcmp(saved_xyz,xyz,sizeof xyz) || memcmp(saved_uv,uv,sizeof uv) ||
            memcmp(saved_faces,faces,sizeof faces) || memcmp(saved_directions,directions,sizeof directions)) bad=1;
        ic_cone_fan_report(&check);
        fprintf(stderr,"  [cone fan feasibility %s] seam=%ld legacy=%ld exact-source -> %s\n",
                fixture==0?"mild-saddle":fixture==1?"severe-saddle":fixture==2?"incompatible-arms":"missing-arms",
                cut.seam_edges,legacy.seam_edges,bad?"FAIL":"ok");
        failures+=bad!=0;Arena_dispose(&arena);
    }
    return failures;
}
