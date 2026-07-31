#ifndef ATLAS_REGISTER_INCLUDED
#define ATLAS_REGISTER_INCLUDED

#include <stddef.h>
#include <stdint.h>

#include "../common/arena.h"

/*
 * Rigid registration of already-parameterized mesh components.
 *
 * The local strip/FEM solve determines shape.  This stage is intentionally
 * unable to change that shape: its only unknown is one additive u shift per
 * physical input component.  Every equation has the form
 *
 *     shift[component1] - shift[component0] = target.
 *
 * Weld equations ask adjacent cube fragments to meet.  Radial equations ask
 * neighbouring wraps to differ by a local circumference.  They are soft
 * evidence, never topology-changing unions.
 */
typedef enum {
    ATLAS_REGISTER_WELD = 0,
    ATLAS_REGISTER_RADIAL = 1,
    ATLAS_REGISTER_LOCAL = 2
} AtlasRegisterEdgeKind;

typedef struct {
    int32_t component0;
    int32_t component1;
    double target;
    double weight;
    /* Huber transition in u units.  Non-positive means quadratic. */
    double robust_scale;
    AtlasRegisterEdgeKind kind;
} AtlasRegisterEdge;

typedef struct {
    size_t ncomponents;
    const AtlasRegisterEdge *edges;
    size_t nedges;

    /* Optional geometric gauge seed.  These values do not enter the Hessian:
     * after relative registration they choose only the one removed additive
     * constant per evidence island.  NaN means unavailable. */
    const double *prior_shift;
    const double *prior_weight;
} AtlasRegisterProblem;

enum { ATLAS_REGISTER_MAX_IRLS = 16 };

typedef struct {
    int iteration;                 /* one-based accepted IRLS iteration */
    double frozen_energy_before;
    double frozen_energy_after;
    double robust_energy_before;
    double robust_energy_after;
    double lambda_min;
    double lambda_max;
    double condition_estimate;     /* of the Jacobi-equilibrated SPD system */
    double raw_diagonal_min;
    double raw_diagonal_max;
    double max_relative_step;
    size_t downweighted_edges;
    double weld_residual_rms;
    double radial_residual_rms;
    double local_residual_rms;
} AtlasRegisterIterationStats;

typedef int (*AtlasRegisterIterationFn)(
    void *context, const AtlasRegisterProblem *problem,
    const double *absolute_shift,
    const AtlasRegisterIterationStats *iteration);

typedef struct {
    double lambda_weld;
    double lambda_radial;
    double lambda_local;
    int irls_rounds;
    int eigen_power_iterations;
    int eigen_inverse_iterations;
    double max_condition;
    double energy_relative_tolerance;
    double convergence_tolerance;
    AtlasRegisterIterationFn iteration_fn;
    void *iteration_context;
} AtlasRegisterOptions;

void AtlasRegisterOptions_default(AtlasRegisterOptions *options);

typedef struct {
    size_t components;
    size_t usable_edges;
    size_t rejected_edges;
    size_t weld_edges;
    size_t radial_edges;
    size_t local_edges;
    size_t evidence_islands;
    size_t variables;
    size_t trusted_priors;
    int irls_rounds_run;
    size_t downweighted_edges;

    double robust_energy_initial;
    double robust_energy_final;
    double weld_residual_rms_before;
    double weld_residual_rms_after;
    double radial_residual_rms_before;
    double radial_residual_rms_after;
    double local_residual_rms_before;
    double local_residual_rms_after;
    double condition_min;
    double condition_max;
    double raw_diagonal_min;
    double raw_diagonal_max;
    double shift_rms;
    double shift_max;

    int solve_failed;
    int condition_failed;
    int energy_failed;
    int trace_failed;

    AtlasRegisterIterationStats iteration[ATLAS_REGISTER_MAX_IRLS];
} AtlasRegisterStats;

/* out_shift may not alias any problem array. */
int AtlasRegister_solve(Arena_T arena,
                        const AtlasRegisterProblem *problem,
                        const AtlasRegisterOptions *options,
                        double *out_shift,
                        AtlasRegisterStats *stats);

int AtlasRegister_selftest(void);

#endif
