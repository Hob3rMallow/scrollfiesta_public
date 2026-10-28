static int smf_selftest(void)
{
    enum {NF=40,NV=3*NF,N=2*NV};int failures=0,cases=0;
    for(int scope=0;scope<2;scope++)for(int sample=0;sample<16;sample++)for(int scenario=0;scenario<4;scenario++) {
        Arena_T arena=Arena_new();float xyz[3*NV];int32_t faces[3*NF],map[NV];
        double reference[N],current[N],trial[N],accepted[N];uint8_t masks[2*NF],guard[2]={1,1};
        SmMetric metric={0};SmBudget budget={0};SmFeedback feedback={0};double angle=.071*sample;
        for(int f=0;f<NF;f++) {
            smb_test_map(reference,f,1,0);masks[f]=1;masks[NF+f]=(uint8_t)(scenario>=2||f<20);
            for(int k=0;k<3;k++) {
                int v=3*f+k;faces[3*f+k]=map[v]=v;
                xyz[3*v]=(float)reference[2*v];xyz[3*v+1]=(float)reference[2*v+1];xyz[3*v+2]=0;
            }
            smb_test_map(current,f,1,angle);
        }
        if(scenario>=2)smb_test_map(reference,0,.75,0);
        smb_test_map(current,0,scenario<2?.89:.75,angle);
        const double initial[]={.95,1.07,.82,1.23},rejected[]={.90,1.11,.79,1.26},safe[]={.93,1.08,.81,1.24};
        smb_test_map(current,1,initial[scenario],angle);memcpy(trial,current,sizeof current);memcpy(accepted,current,sizeof current);
        smb_test_map(trial,1,rejected[scenario],angle);smb_test_map(accepted,1,safe[scenario],angle);
        if(scenario>=2)smb_test_map(trial,0,.74,angle); /* Already bad in the original reference; never a new row. */
        int setup=smm_setup(arena,xyz,NV,faces,NF,reference,&metric);
        if(!setup)setup=smb_setup(arena,&metric,2,masks,&budget);
        if(!setup)setup=smf_setup(arena,NV,&metric,&budget,guard,current,&feedback);
        failures+=setup!=0;cases++;if(setup){Arena_dispose(&arena);continue;}
        PMAMetric measured[2]={{0}};measured[0].area=20;measured[1].area=scenario<2?10:20;
        measured[0].bad10_area=measured[1].bad10_area=1;
        if(scenario>=2){measured[0].bad25_area=measured[1].bad25_area=1;measured[0].bad25_faces=measured[1].bad25_faces=2;}
        feedback.accumulate=scope;
        PMAEvent e={0};e.iteration=1;e.line=4;e.first_failure=PMA_PHYSICAL;e.physical_group=scenario<2?1:0;
        int rc=smf_observe(&feedback,trial,&e,measured,2);int bits=scenario%2?2:1;
        uint8_t *pending=scenario<2?feedback.pending10:feedback.pending25;
        failures+=rc||!feedback.pending||pending[0]||pending[1]!=bits;cases++;
        for(int f=2;f<NF;f++)failures+=pending[f]!=0;
        e.accepted=1;e.line=5;e.first_failure=PMA_OK;
        rc=smf_observe(&feedback,accepted,&e,measured,2);uint8_t *selected=scenario<2?budget.feedback:metric.feedback;
        failures+=rc||selected[1]!=bits||feedback.commits!=1||memcmp(feedback.current,accepted,sizeof accepted);cases++;
        SmBound *rows=NULL;size_t nr=0;
        if(scenario<2)rc=smb_bounds(arena,&budget,map,accepted,&rows,&nr);
        else rc=smm_bounds(arena,&metric,map,accepted,&rows,&nr);
        failures+=rc||nr!=1||(nr&&(rows[0].face!=1||rows[0].upper!=(scenario<2?3+scenario:scenario-2)));cases++;
        if(scenario<2) {
            double recovered[N];memcpy(recovered,accepted,sizeof recovered);smb_test_map(recovered,0,1,angle);
            rows=NULL;nr=0;rc=smb_bounds(arena,&budget,map,recovered,&rows,&nr);
            failures+=rc||nr!=0;cases++; /* Capacity recovered: feedback alone does not spend it twice. */
        }
        e.accepted=0;e.iteration=2;e.line=2;e.first_failure=PMA_PHYSICAL;
        rc=smf_observe(&feedback,trial,&e,measured,2);
        e.line=3;e.first_failure=PMA_RIM_PATH;rc|=smf_observe(&feedback,trial,&e,measured,2);
        e.accepted=1;e.line=4;e.first_failure=PMA_OK;rc|=smf_observe(&feedback,accepted,&e,measured,2);
        failures+=rc||selected[1]!=(scope?bits:0)||feedback.commits!=1;cases++;
        e.accepted=0;e.iteration=3;e.line=2;e.first_failure=PMA_PHYSICAL;
        rc=smf_observe(&feedback,trial,&e,measured,2);
        e.accepted=1;e.iteration=4;e.line=0;e.first_failure=PMA_OK;rc|=smf_observe(&feedback,accepted,&e,measured,2);
        failures+=rc||selected[1]!=(scope?bits:0)||feedback.commits!=1;cases++;
        double invalid[N];memcpy(invalid,accepted,sizeof invalid);invalid[8]=NAN;
        failures+=smf_observe(&feedback,invalid,&e,measured,2)==0||memcmp(feedback.current,accepted,sizeof accepted)||selected[1]!=(scope?bits:0);cases++;
        e.accepted=0;e.iteration=5;e.line=2;e.first_failure=PMA_PHYSICAL;e.physical_group=2;
        failures+=smf_observe(&feedback,trial,&e,measured,2)==0||memcmp(feedback.current,accepted,sizeof accepted);cases++;
        e.physical_group=scenario<2?1:0;
        failures+=smf_observe(&feedback,invalid,&e,measured,2)==0||feedback.pending||memcmp(feedback.current,accepted,sizeof accepted);cases++;
        /* Refinement must add new rows without moving current, terminate on
         * duplicates, and refuse stale/current-mismatched requests. */
        memset(metric.feedback,0,NF);memset(budget.feedback,0,NF);feedback.accumulate=1;
        e.accepted=0;e.iteration=6;e.line=4;e.refinement=2;e.first_failure=PMA_PHYSICAL;
        rc=smf_observe(&feedback,trial,&e,measured,2);e.accepted=1;e.line=5;e.first_failure=PMA_OK;
        size_t added10=0,added25=0;
        int retry=smf_refine(&feedback,accepted,&e,&added10,&added25);
        failures+=rc||retry!=1||added10!=(scenario<2?1u:0u)||added25!=(scenario>=2?1u:0u)||
            selected[1]!=bits||memcmp(feedback.current,accepted,sizeof accepted);cases++;
        failures+=smf_refine(&feedback,accepted,&e,&added10,&added25)!=0||added10||added25;cases++;
        memcpy(invalid,accepted,sizeof invalid);invalid[0]+=.01;
        failures+=smf_refine(&feedback,invalid,&e,&added10,&added25)>=0||selected[1]!=bits;cases++;
        memset(metric.feedback,0,NF);memset(budget.feedback,0,NF);e.refinement=3;
        failures+=smf_refine(&feedback,accepted,&e,&added10,&added25)!=0||added10||added25||selected[1];cases++;
        Arena_dispose(&arena);
    }
    fprintf(stderr,"[source refusal feedback controls] %d cases, %d failures\n",cases,failures);
    return failures?1:0;
}
