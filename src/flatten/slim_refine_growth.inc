/* Shrinking selected vertex steps can propagate a ring of newly flipped
 * triangles without ever reaching a feasible proposal. Instead start at the
 * certified map and grow scalar steps along the SAME global SLIM destination.
 * A coordinate update retains every incident original-to-trial path. For a
 * fixed path time the determinant is affine in that vertex's scalar step, so
 * the feasible interval containing the old scalar is convex. No vertex energy
 * gate is imposed: coherent motion can require locally uphill contributions.
 * The complete proposal still needs strict exact objective descent. */
/* This predicate runs for every incident face of each growth proposal.
 * Inline its unchanged arithmetic to avoid repeated call/register traffic. */
#if defined(_MSC_VER)
#define SR_GROWTH_GUARD_INLINE __forceinline
#elif defined(__GNUC__)
#define SR_GROWTH_GUARD_INLINE inline __attribute__((always_inline))
#else
#define SR_GROWTH_GUARD_INLINE inline
#endif
static SR_GROWTH_GUARD_INLINE int sr_growth_face_safe(const SRSystem *s,const double *current,
                                const double *trial,size_t f)
{
    (void)current;
    const int32_t *face=s->faces+3*f;
    double j[4];
    sr_face_jacobian(s->face+f,face,trial,j);
    double det=j[0]*j[3]-j[1]*j[2];
    return isfinite(det) && det>SR_AREA_EPS &&
           sr_growth_reference_safe(s->growth_reference+f,trial,face);
}
#undef SR_GROWTH_GUARD_INLINE

static size_t sr_growth_candidate(const SRSystem *s,const double *current,
                                   const double *destination,size_t component,
                                   const SlimRefineOpts *opts,const double *limits,
                                   double scale,double *trial,double *vertex_alpha,
                                   int *out_sweeps,size_t *out_face_blocks,
                                   size_t *out_rim_blocks)
{
    size_t begin=s->component_vertex_offset[component];
    size_t end=s->component_vertex_offset[component+1];
    size_t moves=0,face_blocks=0,rim_blocks=0;
    int sweeps=0;
    for (size_t p=begin;p<end;p++) {
        size_t v=(size_t)s->component_vertex[p];
        vertex_alpha[v]=0;
        trial[2*v]=current[2*v];trial[2*v+1]=current[2*v+1];
    }
    for (int sweep=0;sweep<opts->local_growth_sweeps;sweep++) {
        size_t sweep_moves=0;
        sweeps++;
        for (size_t k=0;k<end-begin;k++) {
            size_t p=(sweep&1)?end-1-k:begin+k;
            size_t v=(size_t)s->component_vertex[p];
            double cap=limits?scale*limits[v]:1;
            if (vertex_alpha[v]>=cap || (opts->fixed_vertices && opts->fixed_vertices[v])) continue;
            double old[2]={trial[2*v],trial[2*v+1]},previous=vertex_alpha[v];
            double bound=1;
            /* The endpoint determinant is linear in this coordinate update.
             * Use the same 80%-of-collapse margin as the uniform and ordinary
             * vertex searches; the full path certificate below is additional. */
            for (size_t q=s->incident_offset[v];q<s->incident_offset[v+1];q++) {
                size_t f=s->incident_face[q];double a[4],b[4];
                sr_face_jacobian(s->face+f,s->faces+3*f,trial,a);
                trial[2*v]=current[2*v]+cap*(destination[2*v]-current[2*v]);
                trial[2*v+1]=current[2*v+1]+cap*(destination[2*v+1]-current[2*v+1]);
                sr_face_jacobian(s->face+f,s->faces+3*f,trial,b);
                trial[2*v]=old[0];trial[2*v+1]=old[1];
                double before=a[0]*a[3]-a[1]*a[2],after=b[0]*b[3]-b[1]*b[2];
                if (after<before) bound=fmin(bound,.8*before/(before-after));
            }
            double alpha=previous+(cap-previous)*bound;
            if (!(alpha>previous)) continue;
            for (int local_line=0;local_line<opts->line_search;local_line++) {
                trial[2*v]=current[2*v]+alpha*(destination[2*v]-current[2*v]);
                trial[2*v+1]=current[2*v+1]+alpha*(destination[2*v+1]-current[2*v+1]);
                int safe=1;
                for (size_t q=s->incident_offset[v];q<s->incident_offset[v+1];q++) {
                    if (!sr_growth_face_safe(s,current,trial,s->incident_face[q])) {
                        safe=0;face_blocks++;break;
                    }
                }
                if (safe && opts->guard_boundary && s->boundary_vertex[v] &&
                    !sr_relax_rim_safe(s,component,trial,v,old)) {
                    safe=0;rim_blocks++;
                }
                if (safe) {
                    vertex_alpha[v]=alpha;
                    if (trial[2*v]!=old[0] || trial[2*v+1]!=old[1]) sweep_moves++;
                    break;
                }
                trial[2*v]=old[0];trial[2*v+1]=old[1];
                alpha=.5*(previous+alpha);
            }
        }
        moves+=sweep_moves;
        if (!sweep_moves) break;
    }
    *out_sweeps=sweeps;*out_face_blocks=face_blocks;*out_rim_blocks=rim_blocks;
    return moves;
}

static int sr_component_growth_step(const SRSystem *s,const double *current,
                                     const double *destination,size_t component,
                                     const SlimRefineOpts *opts,double old_energy,
                                     double *trial,double *vertex_alpha,SREdgeBox *boxes,
                                     SRMeasure *candidate,size_t *out_trials,
                                     size_t *out_moves,double *out_alpha_min,
                                     double *out_alpha_mean,double *out_alpha_max)
{
    size_t begin=s->component_vertex_offset[component];
    size_t end=s->component_vertex_offset[component+1];
    size_t face_begin=s->component_face_offset[component];
    size_t face_end=s->component_face_offset[component+1];
    size_t trials=0,total_moves=0,total_face_blocks=0,total_rim_blocks=0;
    int finite=1;
    *out_trials=0;*out_moves=0;
    for (size_t p=begin;p<end;p++) {
        size_t v=(size_t)s->component_vertex[p];
        finite=finite && isfinite(destination[2*v]) && isfinite(destination[2*v+1]);
    }
    if (!finite || end==begin) return 0;
    for (size_t p=face_begin;p<face_end;p++) {
        size_t f=s->component_face[p];
        sr_growth_reference_set(s->growth_reference+f,current,s->faces+3*f);
    }
    if (opts->guard_boundary) sr_growth_rim_index(s,component,current,destination);
    /* Backtrack the first feasible alpha field's upper bounds BEFORE
     * constructing each subsequent feasible proposal. Keeping those bounds
     * preserves the first sweep's damping of a pathological proxy vertex.
     * Uniformly scaling a completed alpha field follows a different boundary
     * path than its certified sequential construction: two slit banks can
     * cross at every small scale although the full endpoint is clear. Start
     * again at current so each candidate follows its own certified rim path. */
    double scale=1;
    for (int line=0;line<opts->line_search;line++,scale*=.5) {
        size_t face_blocks=0,rim_blocks=0;
        int sweeps=0;
        size_t moves=sr_growth_candidate(s,current,destination,component,opts,
            line?s->growth_limit:NULL,scale,trial,vertex_alpha,
            &sweeps,&face_blocks,&rim_blocks);
        if (!line) for (size_t p=begin;p<end;p++) {
            size_t v=(size_t)s->component_vertex[p];
            s->growth_limit[v]=vertex_alpha[v];
        }
        total_moves+=moves;total_face_blocks+=face_blocks;total_rim_blocks+=rim_blocks;
        /* Retain the exact original-to-trial face interval and determinant floor,
         * plus the same final boundary predicate as every other SLIM proposal. */
        int safe=1;
        trials++;
        for (size_t p=face_begin;p<face_end;p++)
            if (!sr_growth_face_safe(s,current,trial,s->component_face[p])) {safe=0;break;}
        size_t contacts=safe && opts->guard_boundary?
            sr_component_boundary_intersections(s,component,trial,boxes,NULL):0;
        int measured=safe && !contacts && !sr_measure_component(s,trial,component,candidate);
        if (measured &&
            !candidate->flips && candidate->energy<old_energy) {
            double minimum=1,maximum=0,sum=0;
            for (size_t p=begin;p<end;p++) {
                size_t v=(size_t)s->component_vertex[p];
                double alpha=vertex_alpha[v];
                minimum=fmin(minimum,alpha);maximum=fmax(maximum,alpha);sum+=alpha;
            }
            *out_trials=trials;*out_moves=moves;*out_alpha_min=minimum;
            *out_alpha_mean=sum/(double)(end-begin);*out_alpha_max=maximum;
            if (opts->verbose) fprintf(stderr,
                "      [slim growth] c%zu accepted sweeps=%d moves=%zu face-blocks=%zu rim-blocks=%zu alpha=%.3e/%.3e/%.3e cap-scale=%.3e E=%.9g->%.9g\n",
                component,sweeps,moves,face_blocks,rim_blocks,
                minimum,*out_alpha_mean,maximum,scale,old_energy,candidate->energy);
            return 1;
        }
        if (opts->verbose && (line<3 || line+1==opts->line_search))
            fprintf(stderr,"      [slim growth] c%zu trial=%d face-safe=%d contacts=%zu E=%.9g->%.9g\n",
                    component,line+1,safe,contacts,old_energy,measured?candidate->energy:INFINITY);
        /* Preserve exact coordinates if the final audit ever disagrees with
         * the sequential path certificate; diagnostics do not alter guards. */
        if (opts->verbose && contacts && (line==1 || line+1==opts->line_search)) {
            fprintf(stderr,"      [slim growth rim] c%zu trial=%d cap-scale=%.17g\n",
                    component,line+1,scale);
            sr_relax_rim_diagnose(s,component,current,trial,boxes);
            size_t current_contacts=sr_component_boundary_intersections(
                s,component,current,boxes,NULL);
            fprintf(stderr,"      [slim growth rim] c%zu current-contacts=%zu\n",
                    component,current_contacts);
            if (current_contacts) sr_relax_rim_diagnose(s,component,current,current,boxes);
        }
    }
    if (opts->verbose) fprintf(stderr,
        "      [slim growth] c%zu refused moves=%zu face-blocks=%zu rim-blocks=%zu\n",
        component,total_moves,total_face_blocks,total_rim_blocks);
    *out_trials=trials;*out_moves=0;
    return 0;
}
