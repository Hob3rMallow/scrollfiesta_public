#ifndef ASM_FIELD_INCLUDED
#define ASM_FIELD_INCLUDED
#include "asm_types.h"
#include "asm_metric.h"
#include "asm_contacts.h"

/* Immutable original faces, metric and measured physical seam rows surrounding
 * a mutable, float64 global field. An omitted/unplaced endpoint remains in the
 * obligation table and cannot become a passing seam by disappearing. */
typedef struct AsmWitness {
    int32_t vertex[2][4];
    double coefficient[2][4];
    double weight; /* one physical arc length; each directional row gets half */
    int valid, required;
} AsmWitness;
typedef struct AsmFieldSeam {
    size_t relation, first, count;
    int32_t a, b;
    int complete, parity, source_run, required, original_cut;
} AsmFieldSeam;
typedef struct AsmSeamMeasure {
    double mass, support, rms, energy;
    size_t observations, invalid;
    int pass;
} AsmSeamMeasure;
typedef struct AsmField {
    Arena_T arena;
    const AsmRun *run;
    size_t nv, nf, charts;
    size_t *vertex_offset, *face_offset;
    float *xyz;
    int32_t *faces, *face_chart, *vertex_chart;
    int32_t *material_vertex; /* validated source-cut aliases, NULL without cuts */
    double *uv, *mass;
    AsmMetricFace *metric;
    AsmFieldSeam *seams;
    AsmWitness *witness;
    size_t n_seams, n_witness;
    AsmContacts *contacts;
} AsmField;

/* Workers every subsequently built field may use to measure its contacts.
 * <= 1 keeps the serial walk; the measured answer is identical either way.
 * Returns the previous setting so a case can restore it. */
int AsmField_threads(int threads);

/* Exact required-obligation predicate used by fields and their audit. A
 * selected requirement survives historical hypothesis flags; unselected
 * vetoed hypotheses do not become requirements merely by carrying source
 * observations. */
int AsmField_relation_required(const AsmRelation *relation);
/* 1 for a recorded tear (ASM_REL_TORN) of a relation that would otherwise be
 * required. A tear is no obligation; the audit ledgers it. */
int AsmField_relation_torn(const AsmRelation *relation);

AsmField *AsmField_new(Arena_T arena, const AsmRun *run);
/* A LOCAL field: only the run's registered charts marked in member[] (NULL =
 * every one of them, i.e. AsmField_new without its contact tree).  A relation
 * with both endpoints absent is dropped; one with a single absent endpoint
 * stays an incomplete obligation, so it can neither pass nor regress.  The
 * coordinates come from `parent` (a whole field over the same run) where it
 * holds the chart, else from the chart itself; the per-face metric is
 * recomputed, which is bit-identical because it is a function of xyz alone.
 * The contact tree is NOT built: an admission transports its candidate first
 * and then calls AsmField_build_contacts, because a tree built before the
 * transport has leaves that span it (2026-09-17: 2^34 leaf pairs). */
AsmField *AsmField_new_local(Arena_T arena, const AsmRun *run, const uint8_t *member, const AsmField *parent);
int AsmField_build_contacts(AsmField *field);
/* Original-material vertex map for a repacked subset. SIZE_MAX offsets omit
 * charts; invalid cut evidence fails instead of exposing artificial borders. */
int AsmField_material_map(Arena_T arena, const AsmRun *run, const size_t *offset,
                          size_t nv, int32_t **material);
/* Commits pointers to the arena-owned global field, preserving source arrays,
 * clean UVs, chart membership, correspondence data and pose parity. */
int AsmField_commit(AsmField *field, AsmRun *run);
/* Source-order full-field checkpoint and exact independent readback. Original
 * winding is retained on disk; source vertex IDs and parity are explicit. */
int AsmField_write(const AsmField *field, const char *out_dir);
int AsmField_verify(const AsmField *field, const char *out_dir);
void AsmField_face(const AsmField *field, size_t face, const double *uv,
                    double j[4], double *lo, double *hi, double *det);
void AsmField_seam(const AsmField *field, size_t seam, const double *uv, AsmSeamMeasure *out);
/* Read-only complete seam certificates, indexed by original relation.
 * Uses the same source quadrature, cut validation and residuals as the full
 * audit, without constructing face metrics or a material-contact tree.
 * out must hold run->n_rels entries; non-obligations remain unresolved. */
int AsmField_measure_source_seams(const AsmRun *run, AsmSeamMeasure *out);
void AsmField_residual(const AsmWitness *w, int direction, const double *uv, double r[2]);
int AsmField_selftest(void);
#endif
