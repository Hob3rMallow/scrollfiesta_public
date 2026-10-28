typedef struct {
    const int32_t *faces;
    int response,calls,bad;
} SeamCutFanTest;

static int seamcut_test_fan_check(void *context,int32_t vertex,
                                  const int32_t *corners,const int32_t *rim,size_t nf)
{
    SeamCutFanTest *test=context;
    test->calls++;
    if (vertex!=0 || nf!=8 || rim[0]!=rim[nf]) test->bad=1;
    uint8_t seen[8]={0};
    for (size_t f=0;f<nf;f++) {
        int32_t cc=corners[f];
        if (cc<0 || cc>=24 || test->faces[cc]!=vertex) { test->bad=1;continue; }
        int32_t a=test->faces[NEXT(cc)],b=test->faces[PREV(cc)];
        if (!((a==rim[f] && b==rim[f+1]) || (b==rim[f] && a==rim[f+1])) || seen[cc/3]) test->bad=1;
        seen[cc/3]=1;
    }
    return test->response;
}

static int seamcut_fan_callback_selftest(Arena_T arena)
{
    int fails=0;
    for (int mixed=0;mixed<2;mixed++) for (int response=-1;response<=1;response++) {
        Arena_Mark mark=Arena_save(arena);
        float vertices[27]={0};int32_t faces[24];uint8_t terminals[9]={2};
        for (int k=0;k<8;k++) {
            double t=2*M_PI*k/8;
            vertices[3*(k+1)]=(float)cos(t);vertices[3*(k+1)+1]=(float)sin(t);
            vertices[3*(k+1)+2]=k&1?.1f:-.1f;
            faces[3*k]=0;faces[3*k+1]=k+1;faces[3*k+2]=(k+1)%8+1;
            if (mixed && k%3==0) { int32_t q=faces[3*k];faces[3*k]=faces[3*k+2];faces[3*k+2]=q; }
        }
        SeamCutFanTest test={faces,response,0,0};SeamCutResult cut={0};
        int rc=SeamCut_to_boundary_checked(arena,vertices,9,faces,8,terminals,seamcut_test_fan_check,&test,&cut);
        int bad=test.bad || test.calls!=1 || (response<0?rc!=-1:rc!=0);
        if (response>=0 && !rc) {
            if (cut.nf!=8 || cut.seam_edges!=(response?1:2)) bad=1;
            for (size_t k=0;k<3*cut.nf;k++) if (cut.vmap[cut.faces[k]]!=faces[k]) bad=1;
            for (size_t v=0;v<cut.nv;v++) if (memcmp(cut.verts+3*v,vertices+3*(size_t)cut.vmap[v],3*sizeof(float))) bad=1;
        }
        fprintf(stderr,"  [seamcut fan callback] mixed-winding=%d response=%d calls=%d -> %s\n",mixed,response,test.calls,bad?"FAIL":"ok");
        fails+=bad!=0;Arena_restore(arena,mark);
    }
    return fails;
}
