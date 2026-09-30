#ifndef ASM_REPAIR_INCLUDED
#define ASM_REPAIR_INCLUDED
#include "asm_field.h"
#include "asm_place.h"
struct ApaFailureCache;
struct ApaFeedback;
struct ApaContinuation;
struct AprWorkspace;
typedef struct AsmRepairOpts {
    int iterations;
    int pose_iterations; /* optional bounded chart-pose phase before vertex deformation */
    double pose_seconds;
    int regional, admission, clearance_iterations, source_trials, admission_rounds;
    int closure_iterations;
    int context_hops;     /* required-seam hops; 0 = legacy whole core, compact scheduler default of 2 */
    int compact;          /* bounded compact transactions; 0 retains the historical scheduler */
    int connectivity_first; /* prioritize required seams by represented area they reconnect */
    size_t patch_faces, patch_field_faces, patch_coordinates;
    size_t patch_bytes; /* local workspace budget; counts above are representational ceilings */
    struct AprWorkspace *workspace; /* measured transaction-local process-memory growth */
    double budget_seconds, patch_seconds;
    double deadline;      /* internal absolute ves_clock_sec deadline, shared by nested solvers */
    double admission_solver_deadline; /* private-region wall cap; serialize legacy solver when positive */
    double admission_metric_margin; /* optional interior reserve before final rigid packing */
    int clearance_method; /* 0 normal; 1 local, 2 supporting translation, 3 material-area diagnostics */
    int clearance_first; /* internal admission classification: warm metric/seam pass, contact failure */
    int finish_clearance; /* stop fitting once original metric and represented seams pass; spend remaining work on contact */
    int stall_steps;      /* 0 off: a fitting or clearance phase ends after this many accepted steps without material progress */
    double stall_ratio;   /* material progress: a defect count falls, or a continuous defect falls by this fraction */
    int restoration_only; /* 0 normal; diagnostics: 1 level, 2 ladder, 3 stretch, 4 separation, 5 forest, 6 forest+ladder */
    double restoration_grid_step; /* diagnostic coordinate level; zero uses original vertices */
    double metric_weight, seam_weight, contact_weight, overlap_weight;
    double contact_length, proximal, proximal_length;
    const uint8_t *fixed_charts; /* optional immutable halo, indexed by run chart */
    const uint8_t *admission_targets; /* optional scheduling scope; exterior charts remain collision context */
    const uint8_t *admission_moving; /* registered neighbors freed with a new chart; replacement overlay keeps current UVs */
    int admission_incident_sources; /* diagnostic: also close eligible deferred source runs incident to the core */
    const char *admission_trace_dir; /* optional original-face endpoint export for a bounded diagnostic */
    struct ApaFailureCache *admission_failure_cache; /* sequence-local completed local refusals only */
    struct ApaFeedback *admission_feedback; /* measured blockers for the automatic neighborhood scheduler */
    struct ApaContinuation *admission_continuation; /* private, uncommitted candidate between bounded attempts */
    int (*contact_merit)(const AsmField *field,const double *uv,double *blocking_area,void *context);
    void *contact_merit_context; /* internal original-material admission obligations */
    int (*contact_objective)(const AsmField *field,const double *uv,double *gradient,
        AsmContactVisit visit,void *visit_context,AsmContactStats *contacts,double *energy,void *context);
    int (*contact_copy_active)(const AsmField *field,size_t a,size_t b,void *context);
    int (*external_guard)(const AsmField *field,const double *uv,void *context);
    void *external_context; /* immutable incumbent for compact line-search guards */
    const AsmPlaceOpts *layer_place; /* post-repair retry with the original placement configuration */
    const AsmConflictOpts *layer_conflict;
    int (*checkpoint)(AsmField *field, void *context);
    void *checkpoint_context;
    int threads;                 /* contact-measurement workers; <= 1 is the serial walk */
} AsmRepairOpts;
typedef struct AsmRepairStats {
    size_t vertices, faces, charts, obligations;
    size_t accepted_rigid, accepted_deformation, accepted_contact, trials, iterations;
    size_t active_coordinates;
    double general_seconds,clearance_seconds,factor_seconds;
    size_t pose_phase_steps;
    double pose_phase_seconds;
    size_t pre_fit_clearance_steps;
    double pre_fit_clearance_seconds;
    size_t contact_merit_evaluations;
    double contact_merit_seconds;
    double energy_before, energy_after, overlap_before, overlap_after;
    size_t contacts_before, contacts_after;
    size_t regions, constrained_steps, clearance_steps, admissions, admitted_charts;
    size_t source_trials, recovered_connections;
    size_t closure_steps, post_closure_clearance_steps;
    size_t restoration_steps, restoration_trials, restoration_skipped_contacts;
    size_t restoration_levels, restoration_max_coordinates;
    size_t stretch_levels, stretch_steps, stretch_pressure_increases;
    size_t initialization_trials, initialization_steps, initialization_constraints, initialization_separated;
    size_t forest_trials, forest_initialized, forest_entry_restarts;
    size_t rejected_fixed, rejected_orientation, rejected_metric, rejected_chart;
    size_t rejected_seam, rejected_contacts, rejected_energy;
    int32_t first_metric_chart;
    size_t first_metric_face;
    int stalled;
    size_t stall_stops;   /* phases ended by stall_steps */
    /* the overlap pre-check (ASM_REPAIR_MAX_OVERLAP_FRACTION) and the contact walk budget */
    int repairable;                  /* 0 = a stack: source-repair trials and regions were skipped */
    size_t grid_cells, stacked_cells;
    double stacked_fraction, precheck_sec;
    size_t regions_refused_budget, trials_refused_budget;
    size_t patches, patches_kept, patches_deferred, patches_oversize;
    size_t source_candidates, source_unresolved, source_unattempted;
    int deadline_reached;
} AsmRepairStats;
void AsmRepair_defaults(AsmRepairOpts *opts);
/* Whole placed field, every original face and all frozen physical witnesses.
 * A bounded solve may stall; zero means a valid non-regressing field was
 * retained. The following independent audit decides geometric qualification. */
int AsmRepair_run(AsmRun *run, const AsmRepairOpts *opts, const char *out_dir, AsmRepairStats *stats);
int AsmRepair_selftest(void);
/* Read-only chart bounds and source boundary pins for the full-resolution viewer. */
int AsmRepair_sheet_review_case(const char *input,const char *out_dir);
/* Both read-only indices from one checkpoint load. */
int AsmRepair_reading_review_case(const char *input,const char *out_dir,const char *axis_path);
int AsmRepair_winding_ribbon_case(const char *input,const char *axis_path,const char *out_dir);
#endif
