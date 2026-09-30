/* Resource eligibility is measured in bytes and elapsed seconds. Original
 * tessellation counts are diagnostics, never evidence of geometric failure.
 * Process sampling includes native sparse factors and contact buffers that
 * do not belong to the field arena. Checks are cooperative, not an OS quota. */
#ifdef _WIN32
#include <psapi.h>
#endif
typedef struct AprWorkspace {
    size_t entry, peak, estimate;
    int measured, exceeded;
} AprWorkspace;

static int apr_process_bytes(size_t *bytes)
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX usage={0};usage.cb=sizeof usage;
    if(!GetProcessMemoryInfo(GetCurrentProcess(),(PROCESS_MEMORY_COUNTERS *)&usage,sizeof usage))return 0;
    *bytes=(size_t)usage.PrivateUsage;return 1;
#else
    FILE *fp=fopen("/proc/self/statm","rb");unsigned long long pages=0,resident=0;
    if(!fp)return 0;
    int okay=fscanf(fp,"%llu %llu",&pages,&resident)==2;fclose(fp);
    long page=sysconf(_SC_PAGESIZE);
    if(!okay || page<=0 || resident>SIZE_MAX/(size_t)page)return 0;
    *bytes=(size_t)resident*(size_t)page;return 1;
#endif
}

static int apr_bytes_add(size_t *bytes,size_t count,size_t width)
{
    if(width && count>(SIZE_MAX-*bytes)/width){*bytes=SIZE_MAX;return 0;}
    *bytes+=count*width;return 1;
}

/* Preflight the field, material-boundary BVH, solver vectors and original
 * sparse assembly. Three boundary edges per face is an upper bound; sparse
 * factor fill and growing contact buffers are monitored in the process.
 * Two field copies allow obstacle expansion while its entry field lives. */
static size_t apr_workspace_estimate(size_t vertices,size_t faces,size_t witnesses,size_t charts)
{
    size_t bytes=0;
    apr_bytes_add(&bytes,charts,2*sizeof(AsmChart)+8*sizeof(size_t)+16*sizeof(double));
    apr_bytes_add(&bytes,vertices,2*(3*sizeof(float)+sizeof(int32_t)+5*sizeof(double))+24*sizeof(double)+sizeof(int));
    apr_bytes_add(&bytes,faces,2*(4*sizeof(int32_t)+sizeof(AsmMetricFace)+sizeof(size_t)+
        sizeof(int32_t)+16+3*(2*sizeof(int32_t)+2*sizeof(size_t)+128))+3*sizeof(double)+sizeof(size_t));
    apr_bytes_add(&bytes,faces,6*(2*sizeof(int)+sizeof(double)));
    apr_bytes_add(&bytes,witnesses,2*sizeof(AsmWitness)+2*sizeof(double)+20*(2*sizeof(int)+sizeof(double)));
    return bytes;
}

static size_t apr_member_workspace(const AsmRun *run,const uint8_t *member)
{
    size_t vertices=0,faces=0,witnesses=0;
    for(size_t c=0;c<run->n_charts;c++)if(member[c]){
        if(!apr_bytes_add(&vertices,run->charts[c].nv,1) || !apr_bytes_add(&faces,run->charts[c].nf,1))return SIZE_MAX;
    }
    for(size_t i=0;i<run->n_rels;i++){
        const AsmRelation *r=run->rels+i;
        if((member[r->a] || member[r->b]) && !apr_bytes_add(&witnesses,r->corr_count,2))return SIZE_MAX;
    }
    return apr_workspace_estimate(vertices,faces,witnesses,run->n_charts);
}

static int apr_workspace_check(const AsmRepairOpts *opts,size_t estimate)
{
    if(!opts->compact)return 1;
    AprWorkspace *work=opts->workspace;size_t current=0;
    if(work){
        if(estimate>work->estimate)work->estimate=estimate;
        if(work->measured && apr_process_bytes(&current)){
            size_t used=current>work->entry?current-work->entry:0;
            if(used>work->peak)work->peak=used;
        }
        if(work->peak>opts->patch_bytes || estimate>opts->patch_bytes)work->exceeded=1;
        return !work->exceeded;
    }
    return estimate<=opts->patch_bytes;
}
