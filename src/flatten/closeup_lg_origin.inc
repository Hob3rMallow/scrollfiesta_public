/* Change only the artificial 1e-4 inward model skin. Policy 1 retains the
 * existing margin whenever the current point can satisfy it, otherwise keeps
 * the current margin. Policy 2 uses the original scale bound. Neither policy
 * weakens the original bound; boundary rows and nonlinear gates are untouched. */
static int smo_bounds(SmMetric *m,SmBound *rows,size_t count)
{
    if(!m||m->origin_policy<0||m->origin_policy>2||(count&&!rows))return -1;
    if(!m->origin_policy)return 0;
    for(size_t r=0;r<count;r++)if(!isfinite(rows[r].rhs))return -1;
    for(size_t r=0;r<count;r++) {
        SmBound *row=rows+r;
        if(row->face==SIZE_MAX||!(row->upper==0||row->upper==1||row->upper==3||row->upper==4))continue;
        double physical=row->rhs-1e-4;
        if(m->origin_policy==2)row->rhs=physical;
        else if(row->rhs>0)row->rhs=fmax(0,physical);
    }
    return 0;
}
