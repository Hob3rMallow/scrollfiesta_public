static int apr_workspace_selftest(void)
{
    AsmRepairOpts opts;AsmRepair_defaults(&opts);opts.compact=1;int failures=0;
    AprWorkspace work={0};opts.workspace=&work;
    size_t ordinary=apr_workspace_estimate(2000,4000,200,8);
    size_t dense=apr_workspace_estimate(160000,320000,2000,22);
    if(!(ordinary<dense) || !apr_workspace_check(&opts,dense) || work.exceeded || work.estimate!=dense)failures++;
    /* Equal face counts can have radically different witness/solver storage. */
    size_t constrained=apr_workspace_estimate(2000,4000,4000000,8);
    if(!(constrained>opts.patch_bytes) || apr_workspace_check(&opts,constrained) || !work.exceeded)failures++;
    work=(AprWorkspace){0};opts.patch_bytes=ordinary-1;
    if(apr_workspace_check(&opts,ordinary) || !work.exceeded)failures++;
    work=(AprWorkspace){0};opts.patch_bytes=ordinary;work.peak=ordinary+1;
    if(apr_workspace_check(&opts,0) || !work.exceeded)failures++;
    if(apr_workspace_estimate(SIZE_MAX,1,1,1)!=SIZE_MAX || apr_workspace_estimate(1,SIZE_MAX,1,1)!=SIZE_MAX)failures++;
    size_t bytes=0;
    if(!apr_process_bytes(&bytes) || !bytes)failures++;
    fprintf(stderr,"  workspace bytes, dense eligible mesh, witness cost, observed growth and overflow: %s\n",failures?"FAIL":"ok");
    return failures;
}
