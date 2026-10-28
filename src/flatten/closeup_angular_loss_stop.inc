/* Progress is normalized by physical loss above the SD minimum, not by its
 * constant 4*area offset. Small accepted movement alone is not convergence:
 * the touching model's directional slope must also be small, and an unresolved
 * refinement/work limit or numerical failure never supplies a certificate. */
static int pmal_begin(PMAStats *s,double objective,double area)
{
    if(!s||!isfinite(objective)||!isfinite(area)||!(area>0))return -1;
    double loss=objective-4*area;
    if(!isfinite(loss)||loss < -1e-12*fmax(1,fabs(objective)))return -1;
    s->loss_floor=4*area;s->initial_loss=s->final_loss=fmax(0,loss);
    s->termination=PMA_STOP_WORK_LIMIT;return 0;
}

static int pmal_update(PMAStats *s,const PMAOptions *o,double objective,double derivative,
                        int complete_step)
{
    if(!s||!o||o->loss_window<=0||!(o->loss_relative_tolerance>0)||
       !isfinite(o->loss_relative_tolerance)||!isfinite(objective)||!isfinite(derivative))return -1;
    double loss=objective-s->loss_floor,scale=fmax(1,s->final_loss);
    if(!isfinite(loss)||loss < -1e-12*fmax(1,fabs(objective)))return -1;
    loss=fmax(0,loss);double decrease=s->final_loss-loss;
    if(decrease < -1e-12*scale)return -1;
    s->last_loss_decrease=decrease;s->last_relative_loss_decrease=decrease/scale;
    s->last_relative_model_slope=fabs(derivative)/scale;s->final_loss=loss;s->loss_checks++;
    int small=complete_step&&decrease>=0&&s->last_relative_loss_decrease<=o->loss_relative_tolerance&&
        s->last_relative_model_slope<=o->loss_relative_tolerance;
    s->loss_stable_steps=small?s->loss_stable_steps+1:0;
    if(s->loss_stable_steps>=o->loss_window){s->loss_converged=1;s->termination=PMA_STOP_LOSS_TOLERANCE;return 1;}
    return 0;
}
