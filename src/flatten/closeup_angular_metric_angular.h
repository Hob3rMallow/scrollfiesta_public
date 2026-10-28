#ifndef CLOSEUP_ANGULAR_INTERNAL_INCLUDED
#define CLOSEUP_ANGULAR_INTERNAL_INCLUDED
#include "slim_refine.h"

/* Internal angular SLIM model. Public callers use closeup_local_global.h. */
enum { PMA_MAPPED_LENGTH_V1=0, PMA_JOINT_ANGLE_V1=1 };
enum { PMA_OK=0, PMA_NONFINITE=1, PMA_FACE_PATH=2, PMA_RIM_PATH=3,
       PMA_NATIVE_BOUNDARY=4, PMA_OBJECTIVE=5, PMA_PHYSICAL=6,
       PMA_MATERIAL=7, PMA_MATERIAL_UNKNOWN=8 };
enum { PMA_STOP_UNSPECIFIED=0, PMA_STOP_WORK_LIMIT=1, PMA_STOP_LOSS_TOLERANCE=2,
       PMA_STOP_NO_SAFE_STEP=3, PMA_STOP_SOLVER_FAILURE=4, PMA_STOP_TARGET=5 };
typedef struct PMAMetric {
    size_t faces,flips,bad10_faces,bad25_faces,measured_faces;
    double area,sd,mapped16,angular16,mass,measured_area;
    double sigma_min,sigma_max,bad10_area,bad25_area,scale_error;
} PMAMetric;
typedef struct PMAEvent {
    int iteration,line,accepted,first_failure,refinement,deferred;
    size_t face_path_failures,rim_path_failures,native_boundary_contacts;
    size_t material_pairs,physical_group;
    double alpha,objective,previous_objective,proxy_derivative,fd_derivative;
} PMAEvent;
typedef struct PMAOptions {
    SlimRefineOpts slim;
    const float *exact_metric_xyz; /* NULL keeps the established diagnostic mode */
    int model,iterations,refinement_limit;
    int final_refinement_limit; /* diagnostic: only extend the final global step */
    int loss_window; /* zero preserves the fixed-work comparison mode */
    double loss_relative_tolerance;
    double target_rms_degrees; /* zero disables; strictly below the target */
    double refinement_progress_alpha; /* zero disables; progress restart is never convergence */
    size_t groups;
    const uint8_t *group_masks; /* group-major nf; immutable source/tile groups */
    const uint8_t *group_count_guard; /* protected bad25 count budget */
    const double *group_bad10_limits; /* optional absolute area caps, bounded by the normal 5% limit */
    const double *reference_uv; /* EXACT same world/faces; original low-distortion UV */
    void *context;
    /* A fresh complete exact-output material query, not per-component. */
    int (*material_guard)(void *,const double *,size_t *pairs);
    int (*observe)(void *,const double *,const PMAEvent *,const PMAMetric *,size_t);
    int (*refine)(void *,const double *,const double *,const PMAEvent *,const PMAMetric *,size_t);
} PMAOptions;
typedef struct PMAStats {
    int status,stop_reason,iterations,accepted_steps;
    size_t trials;
    size_t refinement_solves,refinement_retries,refinement_failures,refinement_limit_hits,refinement_prior_best;
    double initial_objective,final_objective;
    int termination,loss_converged,loss_stable_steps;
    size_t loss_checks,refinement_progress_stops;
    double loss_floor,initial_loss,final_loss,last_loss_decrease,last_relative_loss_decrease;
    double last_proxy_derivative,last_relative_model_slope;
} PMAStats;
int CloseupAngular_run(Arena_T,const float *,size_t,const int32_t *,size_t,
    const double *,const PMAOptions *,double *,PMAStats *);
int CloseupAngular_selftest(void);
int CloseupAngularRefineStep_selftest(void);
int CloseupAngularLossStop_selftest(void);
int CloseupAngularRefineBudget_selftest(void);
int CloseupAngularGlobalRefine_selftest(void);
int CloseupAngularProgressRefine_selftest(void);
int CloseupAngular_trace(Arena_T,const float *,size_t,const int32_t *,size_t,
    const double *,const double *,const int32_t *,const char *);
int CloseupAngular_collect_rim(Arena_T,const float *,size_t,const int32_t *,size_t,
    const double *,const double *,int32_t **,size_t *);
#endif
