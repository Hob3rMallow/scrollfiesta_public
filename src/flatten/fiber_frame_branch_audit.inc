/* Diagnostic only: compare selected fibre arms in transported world frames. */
#include <limits.h>
static double ffb_dot(const double a[3], const double b[3])
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

static void ffb_cross(const double a[3], const double b[3], double c[3])
{
    for(size_t k=0;k<3;k++)c[k]=a[(k+1)%3]*b[(k+2)%3]-a[(k+2)%3]*b[(k+1)%3];
}

static int ffb_axes(const FiberFrame *ff, size_t i, double h[3], double n[3])
{
    double e1[3],e2[3],metric[3],t=ff->theta[i];
    if(!isfinite(t) || ff_physical_frame(ff->pu+3*i,ff->pv+3*i,e1,e2,metric))return -1;
    double a=metric[0]*cos(t)+metric[1]*sin(t),b=metric[2]*sin(t),length=hypot(a,b);
    if(!(length>0) || !isfinite(length))return -1;
    for(size_t k=0;k<3;k++)h[k]=(a*e1[k]+b*e2[k])/length;
    ffb_cross(e1,e2,n);return 0;
}

/* A numerical singularity, not a curvature threshold. The minimal rotation
 * is not uniquely defined for exactly opposing normals. */
static int ffb_transport(const double h[3], const double a[3], const double b[3], double out[3])
{
    double c=fmax(-1.0,fmin(1.0,ffb_dot(a,b))),v[3],first[3],second[3];
    if(1+c<=64*DBL_EPSILON)return -1;
    ffb_cross(a,b,v);ffb_cross(v,h,first);ffb_cross(v,first,second);
    for(size_t k=0;k<3;k++)out[k]=h[k]+first[k]+second[k]/(1+c);
    return 0;
}

static int ffb_path(const float *X, const uint8_t *valid, int H, int W,
                    int ax, int ay, int bx, int by, double *maximum_step)
{
    if(ax<0 || ay<0 || bx<0 || by<0 || ax>=W || bx>=W || ay>=H || by>=H || (ax!=bx && ay!=by))return -1;
    int dx=(bx>ax)-(bx<ax),dy=(by>ay)-(by<ay),x=ax,y=ay;
    size_t previous=0;int first=1;
    *maximum_step=0;
    for(;;){
        size_t i=(size_t)y*(size_t)W+(size_t)x;
        if(!valid[i])return -1;
        double d2=0;
        for(size_t k=0;k<3;k++){
            if(!isfinite(X[3*i+k]))return -1;
            if(!first){double d=(double)X[3*i+k]-X[3*previous+k];d2+=d*d;}
        }
        *maximum_step=fmax(*maximum_step,sqrt(d2));
        if(x==bx && y==by)break;
        previous=i;first=0;x+=dx;y+=dy;
    }
    return 0;
}

static void ffb_example(FiberFrameBranchStats *st, const FiberFrameBranchExample *e)
{
    size_t at=st->examples;
    if(at<8)st->examples++;
    else {
        at=0;
        for(size_t k=1;k<8;k++)if(st->example[k].coherence<st->example[at].coherence)at=k;
        if(e->coherence<=st->example[at].coherence)return;
    }
    st->example[at]=*e;
}

int FiberFrame_audit_branches(Arena_T arena, const FiberFrame *ff,
                              const float *X, const uint8_t *valid, int H, int W,
                              FiberFrameBranchStats *st)
{
    if(!arena || !ff || !X || !valid || !st || H<=0 || W<=0 || ff->gh<=0 || ff->gw<=0 ||
       ff->cell<=0 || !ff->real || !ff->weight || !ff->coherence || !ff->theta || !ff->pu || !ff->pv ||
       (size_t)ff->gh>SIZE_MAX/(size_t)ff->gw/(2*sizeof(double)) ||
       (size_t)H>SIZE_MAX/(size_t)W/3 || ff->gx0<0 || ff->gy0<0 ||
       (int64_t)ff->gx0+(int64_t)(ff->gw-1)*ff->cell>INT_MAX ||
       (int64_t)ff->gy0+(int64_t)(ff->gh-1)*ff->cell>INT_MAX)return -1;
    memset(st,0,sizeof *st);
    Arena_Mark mark=Arena_save(arena);
    size_t nn=(size_t)ff->gh*(size_t)ff->gw,nw=0;
    double *weights=ARENA_ALLOC(arena,nn*sizeof *weights);
    double *selected=ARENA_ALLOC(arena,2*nn*sizeof *selected),*residual=ARENA_ALLOC(arena,2*nn*sizeof *residual);
    double *normal=ARENA_ALLOC(arena,2*nn*sizeof *normal),*coherence=ARENA_ALLOC(arena,2*nn*sizeof *coherence);
    double *step=ARENA_ALLOC(arena,2*nn*sizeof *step);
    uint8_t *evidence=ARENA_CALLOC(arena,nn,1);
    for(size_t i=0;i<nn;i++)if(ff->real[i] && isfinite(ff->weight[i]))weights[nw++]=ff->weight[i];
    qsort(weights,nw,sizeof *weights,ff_cmp_double);
    double median=ff_pct_sorted(weights,nw,50);
    for(size_t i=0;i<nn;i++){
        evidence[i]=(uint8_t)(ff->real[i] && isfinite(ff->weight[i]) &&
            (double)ff->weight[i]>FIBER_FLATTEN_EVIDENCE_REL*median);
        st->evidence_nodes+=evidence[i];
    }
    for(int row=0;row<ff->gh;row++)for(int col=0;col<ff->gw;col++){
        size_t a=(size_t)row*(size_t)ff->gw+(size_t)col;
        if(!evidence[a])continue;
        for(int axis=0;axis<2;axis++){
            int r=row+axis,c=col+1-axis;
            if(r>=ff->gh || c>=ff->gw)continue;
            size_t b=(size_t)r*(size_t)ff->gw+(size_t)c;
            if(!evidence[b])continue;
            st->neighbor_pairs++;
            FiberFrameBranchExample e={0};
            e.col_a=ff->gx0+col*ff->cell;e.row_a=ff->gy0+row*ff->cell;
            e.col_b=ff->gx0+c*ff->cell;e.row_b=ff->gy0+r*ff->cell;
            if(ffb_path(X,valid,H,W,e.col_a,e.row_a,e.col_b,e.row_b,&e.maximum_world_step)){
                st->unsupported_paths++;continue;
            }
            double ha[3],hb[3],na[3],nb[3],transported[3],vb[3];
            if(ffb_axes(ff,a,ha,na) || ffb_axes(ff,b,hb,nb) || !isfinite(ff->coherence[a]) ||
               !isfinite(ff->coherence[b]) || ff->coherence[a]<0 || ff->coherence[a]>1 ||
               ff->coherence[b]<0 || ff->coherence[b]>1){st->invalid_frames++;continue;}
            if(ffb_transport(ha,na,nb,transported)){st->antipodal_normals++;continue;}
            ffb_cross(nb,hb,vb);
            double angle=atan2(ffb_dot(transported,vb),ffb_dot(transported,hb));
            int quarter=(int)floor(angle/(.5*M_PI)+.5);
            e.quarter=(quarter%4+4)%4;
            e.selected_deg=angle*180/M_PI;e.residual_deg=remainder(angle,.5*M_PI)*180/M_PI;
            e.normal_deg=acos(fmax(-1.0,fmin(1.0,ffb_dot(na,nb))))*180/M_PI;
            e.coherence=fmin(ff->coherence[a],ff->coherence[b]);
            size_t p=st->compared_pairs++;
            st->nearest_quarter[e.quarter]++;
            selected[p]=fabs(e.selected_deg);residual[p]=fabs(e.residual_deg);
            if(e.quarter){
                size_t s=st->switches++;
                normal[s]=e.normal_deg;coherence[s]=e.coherence;step[s]=e.maximum_world_step;
                ffb_example(st,&e);
            }
        }
    }
    st->selected_deg=ff_quant(selected,st->compared_pairs);st->residual_deg=ff_quant(residual,st->compared_pairs);
    st->switch_normal_deg=ff_quant(normal,st->switches);st->switch_coherence=ff_quant(coherence,st->switches);
    st->switch_maximum_world_step=ff_quant(step,st->switches);
    Arena_restore(arena,mark);return 0;
}

static int ff_branch_selftest(void)
{
    Arena_T arena=Arena_new();FiberFrame ff={0};FiberFrameBranchStats st={0};
    float X[45]={0},theta[2]={0},weight[2]={1,1},coherence[2]={1,1};
    uint8_t valid[15],real[2]={1,1};double pu[6]={0},pv[6]={0};int fails=0;
    memset(valid,1,sizeof valid);
    for(size_t i=0;i<15;i++){X[3*i]=(float)(i%5);X[3*i+1]=(float)(i/5);}
    ff.gh=1;ff.gw=2;ff.cell=2;ff.gx0=ff.gy0=1;ff.real=real;ff.theta=theta;
    ff.weight=weight;ff.coherence=coherence;ff.pu=pu;ff.pv=pv;
    for(int trial=0;trial<7;trial++){
        memset(pu,0,sizeof pu);memset(pv,0,sizeof pv);memset(valid,1,sizeof valid);
        for(size_t i=0;i<2;i++){
            double basis=(44+2*(double)i)*M_PI/180,psi=-basis;
            pu[3*i]=cos(basis);pu[3*i+1]=sin(basis);pv[3*i]=-sin(basis);pv[3*i+1]=cos(basis);
            theta[i]=(float)(atan2(sin(4*psi),cos(4*psi))/4);
            if(trial)theta[i]=(float)psi;
        }
        if(trial==2){ /* Known quarter-turn of the surface normal, same transported H. */
            for(size_t i=0;i<2;i++){
                pu[3*i]=1;pu[3*i+1]=0;pv[3*i]=0;pv[3*i+1]=i?0:1;pv[3*i+2]=i?1:0;theta[i]=0;
            }
        }
        if(trial==3)theta[1]+=(float)M_PI;
        if(trial==4)for(size_t k=0;k<3;k++)pv[3+k]=-pv[3+k];
        if(trial==5)valid[7]=0;
        if(trial==6)for(size_t k=0;k<3;k++)pu[3+k]=0;
        int rc=FiberFrame_audit_branches(arena,&ff,X,valid,3,5,&st),bad=rc;
        if(trial==0)bad|=st.nearest_quarter[3]!=1 || st.switches!=1 || fabs(st.selected_deg.p50-90)>1e-4 || st.residual_deg.p50>1e-4;
        if(trial==1 || trial==2)bad|=st.nearest_quarter[0]!=1 || st.switches!=0 || st.selected_deg.p50>1e-4;
        if(trial==3)bad|=st.nearest_quarter[2]!=1 || st.switches!=1 || st.residual_deg.p50>1e-4;
        if(trial==4)bad|=st.antipodal_normals!=1 || st.compared_pairs!=0;
        if(trial==5)bad|=st.unsupported_paths!=1 || st.compared_pairs!=0;
        if(trial==6)bad|=st.invalid_frames!=1 || st.compared_pairs!=0;
        fails+=bad!=0;
        fprintf(stderr,"[selftest] fiber branch control %d: %s\n",trial,bad?"FAIL":"PASS");
    }
    Arena_dispose(&arena);
    fprintf(stderr,"[selftest] fiber transported branch audit: %s (%d failures)\n",fails?"FAIL":"PASS",fails);
    return fails;
}
