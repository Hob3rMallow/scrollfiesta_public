/* The restoration initializer includes touching source triangles. Keep this
 * query distinct from positive-area overlap certification. */
static void act_pair(AsmContacts *q,size_t a,size_t b)
{
    if(q->frames && q->frames[a]!=q->frames[b])return;
    if(q->active_faces && !q->active_faces[a] && !q->active_faces[b])return;
    double x[6],y[6],first[6],second[6];int order[3];
    for(int k=0;k<3;k++)for(int d=0;d<2;d++){
        x[2*k+d]=q->uv[2*(size_t)q->faces[3*a+k]+d];y[2*k+d]=q->uv[2*(size_t)q->faces[3*b+k]+d];
    }
    ac_positive(x,first,order,x);ac_positive(y,second,order,x);
    if(!ac_overlap(first,second,-1e-10))return;
    double area=ac_clip_polygon(first,second,NULL,NULL);if(!isfinite(area)){q->failed=1;return;}
    q->stats.pairs++;q->stats.area+=area;if(q->visit)q->visit(q->context,a,b,area);
}
static void act_walk(AsmContacts *q,int32_t a,int32_t b)
{
    if(q->over_budget)return;
    const AcNode *one=q->tree+a,*two=q->tree+b;
    if(!one->active && !two->active)return;
    for(int d=0;d<2;d++)if(one->lo[d]>two->hi[d]+1e-10 || two->lo[d]>one->hi[d]+1e-10)return;
    if(one->left<0 && two->left<0){
        if(ac_leaf_pair_budget(q))return;
        for(size_t i=one->begin;i<one->end;i++)for(size_t j=a==b?i+1:two->begin;j<two->end;j++)act_pair(q,q->order[i],q->order[j]);
    }else if(a==b){act_walk(q,one->left,one->left);act_walk(q,one->left,one->right);act_walk(q,one->right,one->right);}
    else if(two->left<0 || (one->left>=0 && one->end-one->begin>=two->end-two->begin)){
        act_walk(q,one->left,b);act_walk(q,one->right,b);
    }else{act_walk(q,a,two->left);act_walk(q,a,two->right);}
}
int AsmContacts_touching(AsmContacts *q,const double *uv,AsmContactVisit visit,void *context,AsmContactStats *out)
{
    if(out)memset(out,0,sizeof *out);
    if(!out || !ac_refit(q,uv))return -1;
    q->gradient=NULL;q->visit=visit;q->context=context;memset(&q->stats,0,sizeof q->stats);q->failed=0;q->over_budget=0;
    act_walk(q,0,0);*out=q->stats;out->complete=!q->over_budget && !q->failed;
    if(q->over_budget)return -2;
    return q->failed || !isfinite(out->area)?-1:0;
}
static int act_control(void)
{
    int fails=0;const int32_t faces[]={0,1,2,3,4,5};
    for(int frame=0;frame<2;frame++){
        Arena_T arena=Arena_new();double uv[12],angle=frame?.37:0,cs=cos(angle),sn=sin(angle);
        const double source[]={0,0,2,0,0,2,2,0,3,0,2,1};
        for(int v=0;v<6;v++){
            uv[2*v]=cs*source[2*v]-sn*source[2*v+1]+10000*frame;
            uv[2*v+1]=sn*source[2*v]+cs*source[2*v+1]-3000*frame;
        }
        AsmContacts *q=AsmContacts_new(arena,6,faces,2,uv,1e-7);if(!q){fails++;Arena_dispose(&arena);continue;}
        for(int trial=0;trial<4;trial++){
            double offset=trial==0?0:trial==1?5e-11:trial==2?5e-9:-.1;
            for(int v=3;v<6;v++){
                uv[2*v]=cs*(source[2*v]+offset)-sn*source[2*v+1]+10000*frame;
                uv[2*v+1]=sn*(source[2*v]+offset)+cs*source[2*v+1]-3000*frame;
            }
            AsmContactStats touch={0},overlap={0};
            if(AsmContacts_touching(q,uv,NULL,NULL,&touch) || touch.pairs!=(size_t)(trial!=2))fails++;
            if(AsmContacts_measure(q,uv,NULL,NULL,NULL,&overlap) || overlap.pairs!=(size_t)(trial==3))fails++;
            if(trial==3 && fabs(overlap.area-.005)>1e-11)fails++;
        }
        Arena_dispose(&arena);
    }
    fprintf(stderr,"  initializer touching query and unchanged overlap certificate: %s (%d failures)\n",fails?"FAIL":"ok",fails);
    return fails;
}
