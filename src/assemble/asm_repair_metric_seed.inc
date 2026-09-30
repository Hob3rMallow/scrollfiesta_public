/* Bounded endpoint reinitialization of an already represented chart. The
 * original intrinsic map and a unit-scale source projection are proposals,
 * never replacement metrics. The optional local correction shares a deadline
 * and never calls the cold flattener.
 * Each proposal must qualify every original face and preserve passing seams
 * and every existing contact pair's area before the parent field is changed. */
typedef struct AmrMetric {
    size_t bad, invalid;
    double area, within10, lo, hi;
} AmrMetric;

typedef struct AmrResult {
    AmrMetric before, after;
    size_t field_faces, passing_before, passing_after, local_pairs, external_pairs, regressions;
    size_t seed_bad, seed_invalid, free_vertices, active_faces, refine_steps;
    double seconds, factor_seconds, solve_seconds;
    int kept, expired, projection;
    int trace_written;
    int local_checked, exterior_checked;
    int clearance_attempted;
    size_t clearance_steps;
    double clearance_seconds, clearance_area_before, clearance_area_after;
    size_t polish_steps;
    double polish_seconds;
    size_t vertex_polish_steps;
    double vertex_polish_seconds;
    const char *reason;
} AmrResult;

static AmrMetric amr_metric(const AsmField *f, size_t chart, const double *uv)
{
    AmrMetric m = {0}; m.lo = DBL_MAX;
    size_t first = f->face_offset[chart];
    for (size_t k = 0; k < f->run->charts[chart].nf; k++) {
        size_t t = first + k; double j[4], lo, hi, det;
        AsmField_face(f, t, uv, j, &lo, &hi, &det);
        int valid = UvGuard_positive(uv, f->faces + 3*t) && lo > 0 && isfinite(lo) && isfinite(hi);
        m.invalid += !valid; m.bad += !valid || !AsmMetric_within25(lo,hi);
        m.area += f->metric[t].area;
        if (valid) {
            m.lo = fmin(m.lo, lo); m.hi = fmax(m.hi, hi);
            if (AsmMetric_within10(lo,hi)) m.within10 += f->metric[t].area;
        }
    }
    if (m.lo == DBL_MAX) m.lo = 0;
    return m;
}

static int amr_qualified(const AmrMetric *m)
{
    return m->area > 0 && !m->bad && m->within10 >= .95*m->area - 32*DBL_EPSILON*m->area;
}

#include "asm_repair_metric_band.inc"

typedef struct AmrTraceContacts { const AsmField *field; const double *before; FILE *csv; int okay; } AmrTraceContacts;
static void amr_trace_contact(void *context, size_t a, size_t b, double area)
{
    AmrTraceContacts *s=context;const AsmField *f=s->field;
    double old=AsmContacts_pair_faces(f->contacts,s->before,a,b,NULL,NULL);
    int32_t ca = f->face_chart[a], cb = f->face_chart[b];
    s->okay &= fprintf(s->csv,"%d,%zu,%d,%zu,%.17g,%.17g\n",ca,a-f->face_offset[ca],cb,b-f->face_offset[cb],old,area) >= 0;
}

static int amr_trace(AsmField *f, size_t chart, const double *before, const char *dir)
{
    char path[2048]; ves_mkdir(dir); size_t first = f->vertex_offset[chart], nf = f->run->charts[chart].nf;
    const AsmChart *ch = f->run->charts+chart; int okay = 1;
    snprintf(path,sizeof path,"%s/candidate_vertices.csv",dir); FILE *csv = fopen(path,"wb"); if (!csv) return 0;
    okay &= fprintf(csv,"vertex,source_vertex,source_0,source_1,source_2,before_u,before_v,candidate_u,candidate_v\n") >= 0;
    for (size_t v = 0; v < ch->nv && okay; v++) okay &= fprintf(csv,"%zu,%d,%.9g,%.9g,%.9g,%.17g,%.17g,%.17g,%.17g\n",
        v,ch->vid?ch->vid[v]:-1,ch->xyz[3*v],ch->xyz[3*v+1],ch->xyz[3*v+2],before[2*(first+v)],before[2*(first+v)+1],f->uv[2*(first+v)],f->uv[2*(first+v)+1]) >= 0;
    if (fclose(csv)) okay = 0;
    snprintf(path,sizeof path,"%s/candidate_faces.csv",dir); csv = fopen(path,"wb"); if (!csv) return 0;
    okay &= fprintf(csv,"original_face,v0,v1,v2,sigma_lo,sigma_hi,det,source_area\n") >= 0;
    for (size_t k = 0; k < nf && okay; k++) {
        size_t t = f->face_offset[chart]+k; double j[4], lo, hi, det; AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
        okay &= fprintf(csv,"%zu,%zu,%zu,%zu,%.17g,%.17g,%.17g,%.17g\n",k,(size_t)f->faces[3*t]-first,
            (size_t)f->faces[3*t+1]-first,(size_t)f->faces[3*t+2]-first,lo,hi,det,f->metric[t].area) >= 0;
    }
    if (fclose(csv)) okay = 0;
    if (f->contacts) {
        snprintf(path,sizeof path,"%s/candidate_contacts.csv",dir); csv = fopen(path,"wb"); if (!csv) return 0;
        okay &= fprintf(csv,"chart_a,original_face_a,chart_b,original_face_b,before_area,candidate_area\n") >= 0;
        AmrTraceContacts context = {f,before,csv,1}; AsmContactStats stats = {0};
        if (AsmContacts_measure(f->contacts,f->uv,NULL,amr_trace_contact,&context,&stats)) okay = 0;
        okay &= context.okay; if (fclose(csv)) okay = 0;
    }
    return okay;
}

typedef struct AmrPolishVertex { size_t vertex; double gradient; } AmrPolishVertex;

static int amr_polish_order(const void *a,const void *b)
{
    const AmrPolishVertex *x=a,*y=b;
    if(x->gradient!=y->gradient)return x->gradient>y->gradient?-1:1;
    return x->vertex<y->vertex?-1:x->vertex>y->vertex;
}

typedef struct AmrPolishInterface {
    AmbContactSeeds seeds;
    uint8_t *opponents;
} AmrPolishInterface;

static void amr_polish_interface(void *context,size_t a,size_t b,double area)
{
    AmrPolishInterface *p=context;
    size_t ca=(size_t)p->seeds.field->face_chart[a],cb=(size_t)p->seeds.field->face_chart[b];
    amb_contact_seed(&p->seeds,a,b,area);
    if(ca==p->seeds.chart){p->opponents[cb]=1;}
    if(cb==p->seeds.chart){p->opponents[ca]=1;}
}

/* Retain the WHOLE chart-pair intersection, including tolerated boundary
 * triangles and zero-area touching derivatives. Filtering individual
 * certified triangles would again erase the material-boundary force. */
static double amr_polish_weight(void *context,size_t a,size_t b,double raw_area)
{
    (void)raw_area;const AmrPolishInterface *p=context;
    size_t ca=(size_t)p->seeds.field->face_chart[a],cb=(size_t)p->seeds.field->face_chart[b];
    return (ca==p->seeds.chart && p->opponents[cb]) || (cb==p->seeds.chart && p->opponents[ca]);
}

/* Contact residue can stall a simultaneous constrained solve. Probe one
 * contact vertex along the exact area gradient, checking the whole chart,
 * every passing seam and the complete local contact tree at each endpoint.
 * Spend the deadline on physical progress, not tessellation-count gates. */
static size_t amr_vertex_polish_report(AsmField *f,size_t chart,double deadline,AmrPolishStats *stats)
{
    AmrPolishStats local={0};if(!stats)stats=&local;memset(stats,0,sizeof *stats);stats->stop="deadline";
    if(ves_clock_sec()>=deadline)return 0;
    Arena_T arena = Arena_new(); size_t steps = 0;int focused=0;
    size_t first=f->vertex_offset[chart],nv=f->run->charts[chart].nv;
    double *gradient = ARENA_ALLOC(arena,2*f->nv*sizeof(double));
    AmrPolishVertex *order=ARENA_ALLOC(arena,nv*sizeof *order);
    uint8_t *vertices = ARENA_CALLOC(arena,f->nv,1), *passing = ARENA_CALLOC(arena,f->n_seams,1);
    uint8_t *opponents=ARENA_CALLOC(arena,f->run->n_charts,1);
    for (size_t i = 0; i < f->n_seams; i++) { AsmSeamMeasure m; AsmField_seam(f,i,f->uv,&m); passing[i] = (uint8_t)m.pass; }
    while(ves_clock_sec()<deadline){
        stats->rounds++;
        memset(vertices,0,f->nv);memset(opponents,0,f->run->n_charts);
        AmrPolishInterface interface={{f,chart,vertices},opponents};AsmContactStats old = {0};
        double overlap_objective=0;
        int rc;
        if(focused){
            stats->focused_rounds++;
            rc=AsmContacts_measure(f->contacts,f->uv,NULL,amr_polish_interface,&interface,&old);
            if(!rc && old.complete){
                rc=AsmContacts_weighted_objective(f->contacts,f->uv,gradient,NULL,NULL,
                    amr_polish_weight,&interface,&old,&overlap_objective);
            }
        }else{
            rc=AsmContacts_objective(f->contacts,f->uv,gradient,amr_polish_interface,&interface,&old,&overlap_objective);
        }
        if(rc || !old.complete){
            stats->query_failed++;stats->stop="contact_query_failed";break;
        }
        if(stats->rounds==1){stats->pairs_before=old.pairs;stats->area_before=old.area;}
        stats->pairs_after=old.pairs;stats->area_after=old.area;
        if(!old.pairs){stats->stop="clear";break;}
        size_t count=0;
        for(size_t v=first;v<first+nv;v++)if(focused || vertices[v]){
            double length=hypot(gradient[2*v],gradient[2*v+1]);
            if(!(length>1e-12) || !isfinite(length)){stats->zero_gradient++;continue;}
            stats->boundary_candidates+=!vertices[v];
            order[count++]=(AmrPolishVertex){v,length};
        }
        stats->candidates+=count;qsort(order,count,sizeof *order,amr_polish_order);
        size_t best_vertex=SIZE_MAX,best_pairs=old.pairs;
        double best_point[2]={0},best_area=old.area,best_objective=overlap_objective;int zero=0;
        for(size_t k=0;k<count && !zero && ves_clock_sec()<deadline;k++){
            size_t v=order[k].vertex;double length=order[k].gradient;stats->vertices++;
            double saved[2] = {f->uv[2*v],f->uv[2*v+1]}, base = fmin(.02,2*old.area/length+1e-4);
            for (int scale = 0; scale < 4 && ves_clock_sec() < deadline; scale++) {
                stats->trials++;
                double move = fmin(.05,base*ldexp(1.,scale-1));
                for (int d = 0; d < 2; d++) f->uv[2*v+d] = saved[d]-move*gradient[2*v+d]/length;
                AmrMetric metric = amr_metric(f,chart,f->uv); int valid = amr_qualified(&metric);
                stats->metric_refused+=!valid;
                for (size_t i = 0; i < f->n_seams && valid; i++) if (passing[i]) {
                    AsmSeamMeasure m; AsmField_seam(f,i,f->uv,&m);valid=m.pass;stats->seam_refused+=!valid;
                }
                AsmContactStats next = {0};
                if(valid){
                    stats->contact_queries++;
                    double next_objective=0;
                    rc=focused?AsmContacts_weighted_objective(f->contacts,f->uv,NULL,NULL,NULL,
                        amr_polish_weight,&interface,&next,&next_objective):
                        AsmContacts_measure(f->contacts,f->uv,NULL,NULL,NULL,&next);
                    double tolerance=1e-12*fmax(1,best_area);
                    int plateau=focused && next.area<=best_area && fabs(next.area-best_area)<=tolerance &&
                        next_objective<best_objective-1e-12*fmax(1,best_objective);
                    if(rc || !next.complete)stats->query_failed++;
                    else if(!next.pairs || next.area<best_area-tolerance || plateau){
                        best_vertex=v;memcpy(best_point,f->uv+2*v,sizeof best_point);best_area=next.area;best_pairs=next.pairs;
                        best_objective=next_objective;zero=!next.pairs;
                    }else stats->contact_refused++;
                }
                memcpy(f->uv+2*v,saved,sizeof saved);
                if (zero) break;
            }
        }
        if(best_vertex==SIZE_MAX){
            /* A certified triangle may be entirely inside the overlap while
             * its boundary is in the tolerated band. Follow that complete
             * interface and permit raw-area descent across a certificate
             * plateau. Physical contact area may not increase; the caller's
             * full admission certificate remains the acceptance gate. */
            if(!focused){focused=1;continue;}
            stats->stop="no_improving_vertex";break;
        }
        memcpy(f->uv+2*best_vertex,best_point,sizeof best_point); steps++;
        stats->pair_increasing_steps+=best_pairs>old.pairs;
        stats->plateau_steps+=focused && fabs(best_area-old.area)<=1e-12*fmax(1,old.area);
        stats->pairs_after=best_pairs;stats->area_after=best_area;
    }
    if(ves_clock_sec()>=deadline)stats->stop="deadline";
    Arena_dispose(&arena); return steps;
}

static size_t amr_vertex_polish(AsmField *f,size_t chart,double deadline)
{return amr_vertex_polish_report(f,chart,deadline,NULL);}
