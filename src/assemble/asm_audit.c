#include "asm_audit.h"
#include "../common/mesh_bin.h"
#include "../common/sha256.h"
#include "../common/uv_guard.h"
#include "../common/ves_platform.h"
#include "asm_continuity.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../common/ves_omp.h"

typedef struct AaFace { int32_t v[3]; size_t face; } AaFace;
static void aa_key(const int32_t in[3], int32_t out[3])
{
    int at = in[1] < in[0] ? 1 : 0; if (in[2] < in[at]) at = 2;
    for (int k = 0; k < 3; k++) out[k] = in[(at+k)%3]; /* keep original winding */
}
static int aa_compare(const void *one, const void *two)
{
    const AaFace *a = one, *b = two;
    for (int k = 0; k < 3; k++) if (a->v[k] != b->v[k]) return a->v[k] < b->v[k] ? -1 : 1;
    return 0;
}
static double aa_area(const float *x, const int32_t *f)
{
    double a[3], b[3], c[3];
    for (int d = 0; d < 3; d++) { a[d] = (double)x[3*(size_t)f[1]+d]-x[3*(size_t)f[0]+d]; b[d] = (double)x[3*(size_t)f[2]+d]-x[3*(size_t)f[0]+d]; }
    for (int d = 0; d < 3; d++) c[d] = a[(d+1)%3]*b[(d+2)%3]-a[(d+2)%3]*b[(d+1)%3];
    return .5*hypot(hypot(c[0],c[1]),c[2]);
}
static void aa_digest(const void *data, size_t n, char out[65])
{ VesSha256 sha; VesSha256_init(&sha); VesSha256_update(&sha,data,n); VesSha256_hex(&sha,out); }

/* One source cube's ownership check, independent of every other cube: its
 * decoded geometry digest, and per source face the owner state (0 absent,
 * 1 excluded, 2 unplaced, 3 represented) and area.  The charts of the cube
 * claim faces in chart order, exactly as the serial walk did. */
typedef struct AaCube {
    char digest[65];
    size_t nf, altered, foreign, multiply;
    uint8_t *state;   /* malloc [nf] */
    double *area;     /* malloc [nf] */
    int failed;
} AaCube;

static void aa_cube(const AsmRun *run, size_t cube, const int32_t *charts, size_t n_charts, AaCube *out)
{
    memset(out,0,sizeof *out); out->failed = 1;
    MeshBinData mesh = {0};
    if (MeshBin_read_precise_malloc(run->pile[cube].path,&mesh)) return;
    uint64_t sizes[2] = {(uint64_t)mesh.nv,(uint64_t)mesh.nf};
    VesSha256 sha; VesSha256_init(&sha);
    VesSha256_update(&sha,sizes,sizeof sizes); VesSha256_update(&sha,mesh.verts,3*mesh.nv*sizeof(float));
    VesSha256_update(&sha,mesh.faces,3*mesh.nf*sizeof(int32_t)); VesSha256_hex(&sha,out->digest);
    AaFace *index = malloc((mesh.nf ? mesh.nf : 1)*sizeof(AaFace));
    uint8_t *state = calloc(mesh.nf ? mesh.nf : 1,1);
    double *area = malloc((mesh.nf ? mesh.nf : 1)*sizeof(double));
    if (!index || !state || !area) goto fail;
    for (size_t i = 0; i < mesh.nf; i++) { aa_key(mesh.faces+3*i,index[i].v); index[i].face = i; }
    qsort(index,mesh.nf,sizeof *index,aa_compare);
    for (size_t k = 0; k < n_charts; k++) {
        const AsmChart *c = run->charts+charts[k];
        if (!c->vid || !c->xyz || !c->faces) goto fail;
        int state_value = !AsmChart_in_layout(c) ? 1 : !AsmChart_registered(c) ? 2 : 3;
        for (size_t v = 0; v < c->nv; v++) {
            int32_t original = c->vid[v];
            if (original < 0 || (size_t)original >= mesh.nv ||
                memcmp(c->xyz+3*v,mesh.verts+3*(size_t)original,3*sizeof(float))) out->altered++;
        }
        for (size_t t = 0; t < c->nf; t++) {
            int32_t v[3]; int valid = 1;
            for (int j = 0; j < 3; j++) {
                int32_t local = c->faces[3*t+j];
                if (local < 0 || (size_t)local >= c->nv) { valid = 0; break; }
                v[j] = c->vid[local];
            }
            if (!valid) { out->foreign++; continue; }
            AaFace key; aa_key(v,key.v); size_t lo = 0, hi = mesh.nf;
            while (lo < hi) { size_t mid = lo+(hi-lo)/2; if (aa_compare(index+mid,&key) < 0) lo = mid+1; else hi = mid; }
            if (lo == mesh.nf || aa_compare(index+lo,&key)) { out->foreign++; continue; }
            size_t match = lo;
            while (match < mesh.nf && !aa_compare(index+match,&key) && state[index[match].face]) match++;
            if (match == mesh.nf || aa_compare(index+match,&key)) { out->multiply++; continue; }
            state[index[match].face] = (uint8_t)state_value;
        }
    }
    for (size_t t = 0; t < mesh.nf; t++) {
        area[t] = aa_area(mesh.verts,mesh.faces+3*t);
        if (!isfinite(area[t])) goto fail;
    }
    out->nf = mesh.nf; out->state = state; out->area = area; out->failed = 0;
    free(index); MeshBin_dispose(&mesh);   /* dispose clears mesh: nf is taken first */
    return;
fail:
    free(index); free(state); free(area); MeshBin_dispose(&mesh);
}

/* Stream each source cube once. SHA and face ownership use the SAME decoded
 * immutable arrays; a later pathname hash cannot rebind an old checkpoint.
 * Cubes are checked in parallel, a block at a time; every sum, digest and log
 * line is then taken in cube order and, within a cube, in face order, so the
 * result is the serial walk's to the bit (the 21x21x21's 9,115 cubes took
 * 472 s one at a time). */
static int aa_source(Arena_T arena, const AsmRun *run, const char *dir, AsmAuditStats *st)
{
    char path[2048]; snprintf(path,sizeof path,"%s/stage7_source.csv",dir);
    FILE *log = fopen(path,"wb"); if (!log) return -1;
    int io = fprintf(log,"cube,source_faces,accounted_faces,represented_faces,source_area,decoded_geometry_sha256\n") >= 0;
    uint8_t *has_child = ARENA_CALLOC(arena,run->n_charts ? run->n_charts : 1,1);
    for (size_t i = 0; i < run->n_charts; i++) if (run->charts[i].cube < 0 || (size_t)run->charts[i].cube >= run->n_cubes)
        st->foreign_faces += run->charts[i].nf;
    for (size_t i = 0; i < run->n_charts; i++) if (run->charts[i].parent >= 0) {
        if ((size_t)run->charts[i].parent >= run->n_charts || run->charts[i].parent == (int32_t)i) { fclose(log); return -1; }
        has_child[run->charts[i].parent] = 1;
    }
    if (run->n_charts > (size_t)INT32_MAX || run->n_cubes > (size_t)INT32_MAX) { fclose(log); return -1; }
    /* each cube's leaf charts in chart order */
    size_t *first = ARENA_CALLOC(arena,run->n_cubes+1,sizeof(size_t));
    int32_t *owned = ARENA_ALLOC(arena,(run->n_charts ? run->n_charts : 1)*sizeof(int32_t));
    for (size_t i = 0; i < run->n_charts; i++) {
        int32_t cube = run->charts[i].cube;
        if (cube >= 0 && (size_t)cube < run->n_cubes && !has_child[i]) first[(size_t)cube+1]++;
    }
    for (size_t c = 0; c < run->n_cubes; c++) first[c+1] += first[c];
    {
        size_t *fill = ARENA_CALLOC(arena,run->n_cubes ? run->n_cubes : 1,sizeof(size_t));
        for (size_t i = 0; i < run->n_charts; i++) {
            int32_t cube = run->charts[i].cube;
            if (cube >= 0 && (size_t)cube < run->n_cubes && !has_child[i]) owned[first[cube]+fill[cube]++] = (int32_t)i;
        }
    }
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads(); if (threads < 1) threads = 1;
#endif
    size_t block = (size_t)threads*4;
    AaCube *part = ARENA_CALLOC(arena,block,sizeof(AaCube));
    VesSha256 all; VesSha256_init(&all);
    int failed = 0;
    for (size_t start = 0; start < run->n_cubes && !failed; start += block) {
        size_t count = run->n_cubes-start < block ? run->n_cubes-start : block;
        int n = (int)count, k = 0;
#pragma omp parallel for schedule(dynamic,1) num_threads(threads)
        for (k = 0; k < n; k++) {
            size_t cube = start+(size_t)k;
            aa_cube(run,cube,owned+first[cube],first[cube+1]-first[cube],part+k);
        }
        for (size_t k2 = 0; k2 < count; k2++) {
            AaCube *r = part+k2; size_t cube = start+k2;
            if (r->failed || failed) { failed = 1; free(r->state); free(r->area); continue; }
            uint64_t cube_id = (uint64_t)cube; VesSha256_update(&all,&cube_id,sizeof cube_id); VesSha256_update(&all,r->digest,64);
            st->altered_vertices += r->altered; st->foreign_faces += r->foreign; st->multiply_owned_faces += r->multiply;
            size_t accounted = 0, represented = 0; double area = 0;
            for (size_t t = 0; t < r->nf; t++) {
                double a = r->area[t];
                area += a;
                if (r->state[t]) accounted++; else { st->absent_faces++; st->missing_area += a; }
                if (r->state[t] == 1) st->excluded_area += a;
                if (r->state[t] == 2) st->unplaced_area += a;
                if (r->state[t] == 3) { represented++; st->represented_area += a; }
            }
            st->source_faces += r->nf; st->accounted_faces += accounted; st->represented_faces += represented; st->source_area += area;
            if (fprintf(log,"%zu,%zu,%zu,%zu,%.17g,%s\n",cube,r->nf,accounted,represented,area,r->digest) < 0) io = 0;
            free(r->state); free(r->area); r->state = NULL; r->area = NULL;
        }
    }
    if (failed) { fclose(log); return -1; }
    VesSha256_hex(&all,st->source_digest);
    st->source_preserved = run->n_cubes > 0 && !st->absent_faces && !st->multiply_owned_faces && !st->foreign_faces && !st->altered_vertices;
    if (fclose(log)) io = 0; return io ? 0 : -1;
}

typedef struct AaContacts { const AsmField *field; AsmAuditStats *stats; FILE *log; int io; } AaContacts;
static void aa_contact(void *context, size_t a, size_t b, double area)
{
    AaContacts *c = context; const AsmField *f = c->field;
    int32_t ca = f->face_chart[a], cb = f->face_chart[b];
    if (ca == cb) c->stats->self_overlapping_pairs++;
    if (fprintf(c->log,"%d,%zu,%d,%zu,%.17g\n",ca,a-f->face_offset[ca],cb,b-f->face_offset[cb],area) < 0) c->io = 0;
}
static size_t aa_root(size_t *parent, size_t c)
{ while (c != parent[c]) { parent[c] = parent[parent[c]]; c = parent[c]; } return c; }

static int aa_edge_key_order(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return (x > y) - (x < y);
}

/* Border of the flat sheet (AsmAuditStats.border_length). A chart's boundary
 * edges are its edges used by one face; their UV length is its boundary. */
static int aa_border(Arena_T arena, const AsmRun *run, const AsmField *f, const double *stitched, AsmAuditStats *st, FILE *ledger)
{
    if (ledger && fputs("chart,boundary_length,stitched_length,border_length\n",ledger) < 0) return -1;
    for (size_t c = 0; c < run->n_charts; c++) {
        size_t at = f->vertex_offset[c]; const AsmChart *ch = run->charts+c;
        if (at == SIZE_MAX || !ch->nf) continue;
        Arena_Mark mark = Arena_save(arena);
        uint64_t *edge = ARENA_ALLOC(arena,3*ch->nf*sizeof *edge);
        for (size_t t = 0; t < ch->nf; t++) for (int k = 0; k < 3; k++) {
            uint64_t a = (uint32_t)ch->faces[3*t+(size_t)k], b = (uint32_t)ch->faces[3*t+(size_t)(k+1)%3];
            edge[3*t+(size_t)k] = a < b ? (a << 32 | b) : (b << 32 | a);
        }
        qsort(edge,3*ch->nf,sizeof *edge,aa_edge_key_order);
        double boundary = 0;
        for (size_t i = 0; i < 3*ch->nf;) {
            size_t j = i; while (j < 3*ch->nf && edge[j] == edge[i]) j++;
            if (j-i == 1) {
                size_t a = at+(size_t)(edge[i] >> 32), b = at+(size_t)(edge[i] & 0xffffffffu);
                boundary += hypot(f->uv[2*a]-f->uv[2*b],f->uv[2*a+1]-f->uv[2*b+1]);
            }
            i = j;
        }
        Arena_restore(arena,mark);
        if (!isfinite(boundary)) return -1;
        st->boundary_length += boundary; st->stitched_length += fmin(boundary,stitched[c]);
        st->border_length += fmax(0,boundary-stitched[c]);
        if (ledger && fprintf(ledger,"%zu,%.9g,%.9g,%.9g\n",c,boundary,fmin(boundary,stitched[c]),fmax(0,boundary-stitched[c])) < 0) return -1;
    }
    return 0;
}

static int aa_run(const AsmRun *run, const char *dir, AsmAuditStats *st, int require_readback)
{
    memset(st,0,sizeof *st); st->sigma_min = DBL_MAX;
    char path[2048]; snprintf(path,sizeof path,"%s/stage7_audit.json",dir);
    FILE *summary = fopen(path,"wb"); if (!summary) return -1;
    fputs("{\"audit_complete\":false,\"geometry_qualified\":false,\"whole_scroll_qualified\":false}\n",summary);
    if (fclose(summary)) return -1;
    Arena_T arena = Arena_new(); int status = -1;
    FILE *metric = NULL, *contacts = NULL, *seams = NULL, *charts = NULL;
    /* phase clock, printed with the completion line: which part of the
     * audit a larger rung pays for */
    double tp[6] = {0}, tp0 = ves_clock_sec();
    fprintf(stderr,"[assemble audit] checking original source inventory\n");
    if (aa_source(arena,run,dir,st)) goto done;
    tp[0] = ves_clock_sec() - tp0; tp0 = ves_clock_sec();
    fprintf(stderr,"[assemble audit] source checked in %.1f s; building field and contact tree\n",tp[0]);
    AsmField *f = AsmField_new(arena,run); if (!f) goto done;
    tp[1] = ves_clock_sec() - tp0; tp0 = ves_clock_sec();
    fprintf(stderr,"[assemble audit] field built in %.1f s; checking saved geometry\n",tp[1]);
    int repaired = require_readback;
    for (size_t c = 0; c < run->n_charts; c++) if (run->charts[c].placed_uv) repaired = 1;
    if (repaired && AsmField_verify(f,dir)) goto done;
    tp[2] = ves_clock_sec() - tp0; tp0 = ves_clock_sec();
    fprintf(stderr,"[assemble audit] %s in %.1f s; measuring all original faces\n",
            repaired?"saved geometry checked":"in-memory affine field ready",tp[2]);
    st->faces = f->nf; st->placed_charts = f->charts;
    if (st->faces != st->represented_faces) st->source_preserved = 0;
    snprintf(path,sizeof path,"%s/stage7_metric_faces.csv",dir); metric = fopen(path,"wb");
    snprintf(path,sizeof path,"%s/stage7_contacts.csv",dir); contacts = fopen(path,"wb");
    snprintf(path,sizeof path,"%s/stage7_seams.csv",dir); seams = fopen(path,"wb");
    snprintf(path,sizeof path,"%s/stage7_charts.csv",dir); charts = fopen(path,"wb");
    if (!metric || !contacts || !seams || !charts) goto done;
    int io = fprintf(metric,"chart,original_chart_face,sigma_min,sigma_max,determinant,source_area\n") >= 0;
    io &= fprintf(contacts,"chart_a,original_face_a,chart_b,original_face_b,intersection_area\n") >= 0;
    io &= fprintf(seams,"relation,chart_a,chart_b,declared_observations,evaluated_observations,invalid,mass,supported_mass,rms,complete,parity,source_run,pass,repair_required,reconstructed,relation_flags,original_mesh_cut\n") >= 0;
    io &= fprintf(charts,"chart,cube,flags,placed,faces,area,within10_area,within25_area,invalid,strict_bad\n") >= 0;
    double *area = ARENA_CALLOC(arena,run->n_charts,sizeof(double));
    double *original_chart_area = ARENA_CALLOC(arena,run->n_charts,sizeof(double));
    double *within10 = ARENA_CALLOC(arena,run->n_charts,sizeof(double)), *within25 = ARENA_CALLOC(arena,run->n_charts,sizeof(double));
    size_t *bad = ARENA_CALLOC(arena,run->n_charts,sizeof(size_t)), *invalid = ARENA_CALLOC(arena,run->n_charts,sizeof(size_t));
    for (size_t t = 0; t < f->nf; t++) {
        double j[4], lo, hi, det; AsmField_face(f,t,f->uv,j,&lo,&hi,&det);
        int32_t chart = f->face_chart[t]; double a = f->metric[t].area; area[chart] += a;
        int valid = UvGuard_positive(f->uv,f->faces+3*t) && det > 0 && lo > 0 && isfinite(lo) && isfinite(hi);
        if (!valid) { st->invalid_faces++; invalid[chart]++; }
        else {
            st->sigma_min = fmin(st->sigma_min,lo); st->sigma_max = fmax(st->sigma_max,hi);
            if (AsmMetric_within10(lo,hi)) { st->within10_area += a; within10[chart] += a; }
            if (AsmMetric_within25(lo,hi)) { st->within25_area += a; within25[chart] += a; }
        }
        if (!valid || !AsmMetric_within25(lo,hi)) {
            st->strict_bad_faces++; bad[chart]++;
            if (fprintf(metric,"%d,%zu,%.17g,%.17g,%.17g,%.17g\n",chart,t-f->face_offset[chart],lo,hi,det,a) < 0) io = 0;
        }
    }
    for (size_t c = 0; c < run->n_charts; c++) {
        const AsmChart *ch = run->charts+c;
        double original_area = 0;
        for (size_t t = 0; t < ch->nf; t++) original_area += aa_area(ch->xyz,ch->faces+3*t);
        original_chart_area[c] = original_area;
        if (!AsmChart_in_layout(ch)) st->excluded_charts++; else if (!AsmChart_registered(ch)) st->unplaced_charts++;
        double roundoff=32*DBL_EPSILON*area[c];
        if (f->vertex_offset[c] != SIZE_MAX && (invalid[c] || within10[c]+roundoff < .95*area[c] || within25[c]+roundoff < .99*area[c])) st->bad_metric_charts++;
        if (fprintf(charts,"%zu,%d,%u,%d,%zu,%.17g,%.17g,%.17g,%zu,%zu\n",c,ch->cube,ch->flags,AsmChart_registered(ch),ch->nf,original_area,within10[c],within25[c],invalid[c],bad[c]) < 0) io = 0;
    }
    tp[3] = ves_clock_sec() - tp0; tp0 = ves_clock_sec();
    fprintf(stderr,"[assemble audit] face metrics checked in %.1f s; walking material contacts\n",tp[3]);
    st->contacts_complete = 1;
    if (f->nf) {
        AaContacts visit = {f,st,contacts,1}; AsmContactStats cs = {0};
        int contact_rc = AsmContacts_measure(f->contacts,f->uv,NULL,aa_contact,&visit,&cs);
        if (!visit.io || (contact_rc && contact_rc != -2)) goto done;
        st->contact_leaf_pairs = cs.leaf_pairs; st->contact_leaf_pair_budget = AsmContacts_budget(f->contacts);
        if (contact_rc == -2) {
            /* a stack: the walk cannot finish, the overlaps are UNMEASURED (not zero) and the
             * sheet cannot qualify; the audit still completes so the sheet exports with its gates */
            st->contacts_complete = 0; st->overlapping_pairs = 0; st->overlap_area = 0; st->self_overlapping_pairs = 0;
        } else { st->overlapping_pairs = cs.pairs; st->overlap_area = cs.area; }
    }
    tp[4] = ves_clock_sec() - tp0; tp0 = ves_clock_sec();
    fprintf(stderr,"[assemble audit] contacts checked in %.1f s; measuring seams and hashing the field\n",tp[4]);
    size_t *parent = ARENA_ALLOC(arena,run->n_charts*sizeof(size_t));
    size_t *source_parent = ARENA_ALLOC(arena,run->n_charts*sizeof(size_t));
    double *stitched = ARENA_CALLOC(arena,run->n_charts ? run->n_charts : 1,sizeof(double));
    for (size_t c = 0; c < run->n_charts; c++) parent[c] = source_parent[c] = c;
    VesSha256 obligation_hash; VesSha256_init(&obligation_hash);
    for (size_t s = 0; s < f->n_seams; s++) {
        const AsmFieldSeam *seam = f->seams+s; AsmSeamMeasure m; AsmField_seam(f,s,f->uv,&m);
        st->obligations++; if (m.pass) st->passing_seams++; else st->unresolved_seams++;
        if (m.pass) { stitched[seam->a] += m.mass; stitched[seam->b] += m.mass; }
        if(seam->required){
            st->repair_obligations++;
            source_parent[aa_root(source_parent,seam->a)] = aa_root(source_parent,seam->b);
            if(m.pass){st->passing_repair_seams++;parent[aa_root(parent,seam->a)]=aa_root(parent,seam->b);}
            else {
                st->unresolved_repair_seams++;
                /* Missing material belongs in coverage; an unmeasurable
                 * seam between represented charts remains a real failure. */
                if(f->vertex_offset[seam->a]!=SIZE_MAX && f->vertex_offset[seam->b]!=SIZE_MAX)
                    st->unresolved_present_repair_seams++;
                else st->unresolved_missing_repair_seams++;
            }
        }
        if (fprintf(seams,"%zu,%d,%d,%zu,%zu,%zu,%.17g,%.17g,%.17g,%d,%d,%d,%d,%d,%d,%u,%d\n",seam->relation,seam->a,seam->b,seam->count,m.observations,m.invalid,m.mass,m.support,m.rms,seam->complete,seam->parity,seam->source_run,m.pass,seam->required,!!(run->rels[seam->relation].continuity&ASM_CONT_RECONSTRUCTED),run->rels[seam->relation].flags,seam->original_cut) < 0) io = 0;
        /* Hash explicit fields, excluding C struct padding. Frozen rows bind
         * the original physical witness data used by this exact audit. */
        uint64_t ids[4] = {(uint64_t)seam->relation,(uint64_t)seam->a,(uint64_t)seam->b,(uint64_t)seam->count};
        VesSha256_update(&obligation_hash,ids,sizeof ids);
        const AsmRelation *original = run->rels+seam->relation;
        VesSha256_update(&obligation_hash,&original->flags,sizeof original->flags);
        VesSha256_update(&obligation_hash,&original->continuity,sizeof original->continuity);
        VesSha256_update(&obligation_hash,&seam->required,sizeof seam->required);
        for (size_t k = 0; k < seam->count; k++) {
            const AsmCorr *observation = run->corr+(size_t)original->corr_first+k;
            VesSha256_update(&obligation_hash,&observation->va,sizeof observation->va);
            VesSha256_update(&obligation_hash,&observation->vb,sizeof observation->vb);
            VesSha256_update(&obligation_hash,observation->gap_a,sizeof observation->gap_a);
            VesSha256_update(&obligation_hash,observation->gap_b,sizeof observation->gap_b);
            VesSha256_update(&obligation_hash,&observation->valid,sizeof observation->valid);
            VesSha256_update(&obligation_hash,&observation->run,sizeof observation->run);
            const AsmTrim *trim[2] = {&observation->trim_a,&observation->trim_b};
            for (int side = 0; side < 2; side++) {
                VesSha256_update(&obligation_hash,&trim[side]->valid,sizeof trim[side]->valid);
                VesSha256_update(&obligation_hash,&trim[side]->face,sizeof trim[side]->face);
                VesSha256_update(&obligation_hash,trim[side]->inverse,sizeof trim[side]->inverse);
            }
            const AsmWitness *w = f->witness+seam->first+k;
            VesSha256_update(&obligation_hash,w->vertex,sizeof w->vertex); VesSha256_update(&obligation_hash,w->coefficient,sizeof w->coefficient);
            VesSha256_update(&obligation_hash,&w->weight,sizeof w->weight); VesSha256_update(&obligation_hash,&w->valid,sizeof w->valid);
        }
    }
    VesSha256_hex(&obligation_hash,st->obligations_digest);
    double *connected = ARENA_CALLOC(arena,run->n_charts,sizeof(double));
    for (size_t c = 0; c < run->n_charts; c++) if (f->vertex_offset[c] != SIZE_MAX) connected[aa_root(parent,c)] += area[c];
    for (size_t c = 0; c < run->n_charts; c++) if (connected[c] > 0) { st->continuity_components++; st->largest_connected_area = fmax(st->largest_connected_area,connected[c]); }
    /* A crop can contain unrelated arcs of the same sheet. Measure the
     * retained continuity inside each recorded source region, including
     * breaks across missing charts, rather than demanding nonexistent joins
     * between arcs. The complete original source ledger remains unchanged. */
    double *possible = ARENA_CALLOC(arena,run->n_charts,sizeof(double));
    double *region_largest = ARENA_CALLOC(arena,run->n_charts,sizeof(double));
    size_t *region_parts = ARENA_CALLOC(arena,run->n_charts,sizeof(size_t));
    for(size_t c=0;c<run->n_charts;c++) {
        size_t root=aa_root(source_parent,c);
        possible[root]+=original_chart_area[c];
        if(connected[c]>0){region_parts[root]++;region_largest[root]=fmax(region_largest[root],connected[c]);}
    }
    for(size_t c=0;c<run->n_charts;c++){
        st->source_regions+=possible[c]>0;
        st->largest_possible_connected_area=fmax(st->largest_possible_connected_area,possible[c]);
        st->represented_source_regions+=region_parts[c]>0;
        st->fractured_source_regions+=region_parts[c]>1;
        st->source_region_connected_area+=region_largest[c];
    }
    /* Recorded tears (ASM_REL_TORN) are no obligations and the regions above
     * are measured across them. Ledger every tear, and measure the source
     * regions again with the torn seams joined, so what the tears split stays
     * visible beside the qualification. */
    {
        snprintf(path,sizeof path,"%s/stage7_tears.csv",dir); FILE *tears = fopen(path,"wb"); if (!tears) goto done;
        io &= fprintf(tears,"relation,chart_a,chart_b,kind,seam_length,chart_a_placed,chart_b_placed\n") >= 0;
        size_t *untorn = ARENA_ALLOC(arena,run->n_charts*sizeof(size_t));
        for (size_t c = 0; c < run->n_charts; c++) untorn[c] = aa_root(source_parent,c);
        for (size_t i = 0; i < run->n_rels; i++) {
            const AsmRelation *r = run->rels+i; if (!AsmField_relation_torn(r)) continue;
            size_t a = (size_t)r->a, b = (size_t)r->b; if (a >= run->n_charts || b >= run->n_charts) continue;
            int pa = f->vertex_offset[a] != SIZE_MAX, pb = f->vertex_offset[b] != SIZE_MAX;
            st->torn_seams++; st->torn_length += r->seam_len;
            if (a == b) st->torn_cuts++; else if (pa && pb) st->torn_between_pieces++;
            untorn[aa_root(untorn,a)] = aa_root(untorn,b);
            io &= fprintf(tears,"%zu,%d,%d,%s,%.17g,%d,%d\n",i,r->a,r->b,a == b ? "intrachart_cut" : "between_charts",r->seam_len,pa,pb) >= 0;
        }
        if (fclose(tears)) io = 0;
        double *largest = ARENA_CALLOC(arena,run->n_charts,sizeof(double)), *whole = ARENA_CALLOC(arena,run->n_charts,sizeof(double));
        size_t *parts = ARENA_CALLOC(arena,run->n_charts,sizeof(size_t));
        for (size_t c = 0; c < run->n_charts; c++) {
            size_t root = aa_root(untorn,c); whole[root] += original_chart_area[c];
            if (connected[c] > 0) { parts[root]++; largest[root] = fmax(largest[root],connected[c]); }
        }
        for (size_t c = 0; c < run->n_charts; c++) {
            st->untorn_source_regions += whole[c] > 0; st->untorn_fractured_source_regions += parts[c] > 1;
            st->untorn_source_region_connected_area += largest[c];
        }
    }
    {
        snprintf(path,sizeof path,"%s/stage7_border.csv",dir); FILE *border = fopen(path,"wb");
        int bad = !border || aa_border(arena,run,f,stitched,st,border);
        if (border && fclose(border)) bad = 1;
        if (bad) goto done;
    }
    VesSha256 hash; VesSha256_init(&hash);
    uint64_t dimensions[2] = {(uint64_t)f->nv,(uint64_t)f->nf}; VesSha256_update(&hash,dimensions,sizeof dimensions);
    if (f->nv) { VesSha256_update(&hash,f->xyz,3*f->nv*sizeof(float)); VesSha256_update(&hash,f->uv,2*f->nv*sizeof(double)); }
    if (f->nf) { VesSha256_update(&hash,f->faces,3*f->nf*sizeof(int32_t)); VesSha256_update(&hash,f->face_chart,f->nf*sizeof(int32_t)); }
    VesSha256_hex(&hash,st->field_digest);
    if (!io) goto done;
    if (fclose(metric)) { metric = NULL; goto done; } metric = NULL;
    if (fclose(contacts)) { contacts = NULL; goto done; } contacts = NULL;
    if (fclose(seams)) { seams = NULL; goto done; } seams = NULL;
    if (fclose(charts)) { charts = NULL; goto done; } charts = NULL;
    st->complete = 1;
    st->legacy_geometry_qualified = st->source_preserved && st->faces > 0 && !st->invalid_faces && !st->strict_bad_faces &&
        !st->bad_metric_charts && st->contacts_complete && !st->overlapping_pairs && !st->unresolved_seams &&
        st->represented_area >= .90*st->source_area && st->largest_connected_area >= .99*st->represented_area;
    st->geometry_qualified = st->source_preserved && st->faces > 0 && !st->invalid_faces && !st->strict_bad_faces &&
        !st->bad_metric_charts && st->contacts_complete && !st->overlapping_pairs && !st->unresolved_present_repair_seams &&
        st->represented_area >= .85*st->source_area && st->source_region_connected_area >= .99*st->represented_area;
    if (st->sigma_min == DBL_MAX) st->sigma_min = 0;
    snprintf(path,sizeof path,"%s/stage7_audit.json",dir); summary = fopen(path,"wb"); if (!summary) goto done;
    int wrote = fprintf(summary,
        "{\n  \"schema\":\"native-original-field-v2\",\n  \"audit_complete\":true,\n  \"geometry_qualified\":%s,\n"
        "  \"whole_scroll_qualified\":false,\n  \"material_identity_review\":\"Geometry alone does not establish CT material identity\",\n"
        "  \"source_preserved\":%s,\n  \"source_faces\":%zu,\n  \"accounted_source_faces\":%zu,\n  \"represented_source_faces\":%zu,\n"
        "  \"absent_source_faces\":%zu,\n  \"multiply_owned_source_faces\":%zu,\n  \"foreign_faces\":%zu,\n  \"altered_source_vertices\":%zu,\n"
        "  \"source_area\":%.17g,\n  \"represented_area\":%.17g,\n  \"excluded_area\":%.17g,\n  \"unplaced_area\":%.17g,\n  \"missing_area\":%.17g,\n"
        "  \"coverage\":%.17g,\n  \"placed_charts\":%zu,\n  \"unplaced_charts\":%zu,\n  \"excluded_charts\":%zu,\n  \"audited_faces\":%zu,\n"
        "  \"invalid_faces\":%zu,\n  \"strict_bad_faces\":%zu,\n  \"bad_metric_charts\":%zu,\n  \"sigma_min\":%.17g,\n  \"sigma_max\":%.17g,\n"
        "  \"within10_area\":%.17g,\n  \"within25_area\":%.17g,\n  \"contact_walk_complete\":%s,\n  \"contact_leaf_pairs\":%zu,\n  \"contact_leaf_pair_budget\":%zu,\n"
        "  \"overlapping_triangle_pairs\":%zu,\n  \"self_overlapping_pairs\":%zu,\n  \"overlap_area_sum\":%.17g,\n"
        "  \"obligations\":%zu,\n  \"passing_seams\":%zu,\n  \"unresolved_seams\":%zu,\n  \"continuity_components\":%zu,\n  \"largest_connected_area\":%.17g,\n"
        "  \"repair_obligations\":%zu,\n  \"passing_repair_seams\":%zu,\n  \"unresolved_repair_seams\":%zu,\n  \"continuity_graph\":\"Passing selected repair obligations; deferred proposals remain separately audited\",\n"
        "  \"contact_tolerance_vox\":%.17g,\n  \"contact_policy\":\"%s\",\n  \"metric_singular_value_roundoff\":%.17g,\n  \"seam_rms_limit_vox\":%.17g,\n  \"seam_support_limit_vox\":%.17g,\n  \"minimum_supported_seam_mass_fraction\":0.9,\n"
        "  \"original_mesh_cut_requires_exact_source_edge\":true,\n  \"original_mesh_cut_endpoint_count\":2,\n  \"original_mesh_cut_supported_mass_fraction\":1,\n"
        "  \"qualification_contract\":\"85%% original area, continuity within recorded source regions, all represented required seams\",\n"
        "  \"minimum_coverage\":0.85,\n  \"minimum_source_region_coherence\":0.99,\n"
        "  \"legacy_geometry_qualified\":%s,\n  \"source_regions\":%zu,\n  \"represented_source_regions\":%zu,\n  \"fractured_source_regions\":%zu,\n"
        "  \"largest_possible_connected_area\":%.17g,\n  \"source_region_connected_area\":%.17g,\n  \"source_region_coherence\":%.17g,\n"
        "  \"unresolved_present_repair_seams\":%zu,\n  \"unresolved_missing_repair_seams\":%zu,\n"
        "  \"tear_policy\":\"recorded tears (stage7_tears.csv) are no obligations; source regions are measured across them and without them\",\n"
        "  \"torn_seams\":%zu,\n  \"torn_intrachart_cuts\":%zu,\n  \"torn_seams_between_pieces\":%zu,\n  \"torn_seam_length\":%.17g,\n"
        "  \"untorn_source_regions\":%zu,\n  \"untorn_fractured_source_regions\":%zu,\n  \"untorn_source_region_coherence\":%.17g,\n"
        "  \"border_policy\":\"placed charts' UV boundary length less the arc length of passing seams on each side\",\n"
        "  \"boundary_length_vox\":%.17g,\n  \"stitched_boundary_length_vox\":%.17g,\n  \"border_length_vox\":%.17g,\n  \"border_per_sqrt_area\":%.17g,\n"
        "  \"source_decoded_sha256\":\"%s\",\n  \"field_sha256\":\"%s\",\n  \"obligations_sha256\":\"%s\"\n}\n",
        st->geometry_qualified ? "true":"false",st->source_preserved ? "true":"false",st->source_faces,st->accounted_faces,st->represented_faces,
        st->absent_faces,st->multiply_owned_faces,st->foreign_faces,st->altered_vertices,st->source_area,st->represented_area,st->excluded_area,st->unplaced_area,st->missing_area,
        st->source_area > 0 ? st->represented_area/st->source_area : 0,st->placed_charts,st->unplaced_charts,st->excluded_charts,st->faces,
        st->invalid_faces,st->strict_bad_faces,st->bad_metric_charts,st->sigma_min,st->sigma_max,st->within10_area,st->within25_area,
        st->contacts_complete ? "true":"false",st->contact_leaf_pairs,st->contact_leaf_pair_budget,
        st->overlapping_pairs,st->self_overlapping_pairs,st->overlap_area,st->obligations,st->passing_seams,st->unresolved_seams,st->continuity_components,st->largest_connected_area,
        st->repair_obligations,st->passing_repair_seams,st->unresolved_repair_seams,AsmContacts_material_tolerance(),AsmContacts_material_policy(),(double)ASM_METRIC_ROUNDOFF,
        (double)ASM_SEAM_TOLERANCE_VOX,(double)ASM_SEAM_TOLERANCE_VOX,
        st->legacy_geometry_qualified ? "true":"false",st->source_regions,st->represented_source_regions,st->fractured_source_regions,
        st->largest_possible_connected_area,st->source_region_connected_area,
        st->represented_area>0 ? st->source_region_connected_area/st->represented_area : 0,
        st->unresolved_present_repair_seams,st->unresolved_missing_repair_seams,
        st->torn_seams,st->torn_cuts,st->torn_between_pieces,st->torn_length,
        st->untorn_source_regions,st->untorn_fractured_source_regions,
        st->represented_area>0 ? st->untorn_source_region_connected_area/st->represented_area : 0,
        st->boundary_length,st->stitched_length,st->border_length,
        st->represented_area>0 ? st->border_length/sqrt(st->represented_area) : 0,
        st->source_digest,st->field_digest,st->obligations_digest) >= 0;
    if (fclose(summary)) wrote = 0;
    if (!wrote) { st->complete = 0; goto done; }
    char overlaps[96];
    if (st->contacts_complete) snprintf(overlaps,sizeof overlaps,"%zu overlaps",st->overlapping_pairs);
    else snprintf(overlaps,sizeof overlaps,"overlaps UNMEASURED (contact walk stopped at %zu leaf pairs, budget %zu)",st->contact_leaf_pairs,st->contact_leaf_pair_budget);
    fprintf(stderr,"[assemble audit] complete: %zu/%zu source faces represented (%.3f%% area), %zu invalid, %zu strict metric failures, %s, seams %zu/%zu, %zu continuity components; geometry %s\n",
        st->represented_faces,st->source_faces,st->source_area > 0 ? 100*st->represented_area/st->source_area : 0,st->invalid_faces,st->strict_bad_faces,
        overlaps,st->passing_seams,st->obligations,st->continuity_components,st->geometry_qualified ? "PASS":"UNQUALIFIED");
    fprintf(stderr,"[assemble audit] border %.0f vox (%.2f per sqrt area): boundary %.0f vox, %.0f stitched by passing seams\n",
        st->border_length,st->represented_area>0 ? st->border_length/sqrt(st->represented_area) : 0,st->boundary_length,st->stitched_length);
    if (st->torn_seams)
        fprintf(stderr,"[assemble audit] recorded tears: %zu (%zu inside charts, %zu between placed pieces, %.1f vox); %zu of %zu source regions fractured without them\n",
            st->torn_seams,st->torn_cuts,st->torn_between_pieces,st->torn_length,st->untorn_fractured_source_regions,st->untorn_source_regions);
    tp[5] = ves_clock_sec() - tp0;
    fprintf(stderr,"[assemble audit] seconds: source %.1f, field and contact tree %.1f, export readback %.1f, face metrics %.1f, contact walk %.1f, seams and report %.1f\n",
        tp[0],tp[1],tp[2],tp[3],tp[4],tp[5]);
    status = 0;
done:
    if (metric) fclose(metric); if (contacts) fclose(contacts); if (seams) fclose(seams); if (charts) fclose(charts);
    Arena_dispose(&arena); return status;
}

int AsmAudit_run(const AsmRun *run,const char *dir,AsmAuditStats *out)
{ return aa_run(run,dir,out,0); }

int AsmAudit_checkpoint(const AsmRun *run,const char *dir,AsmAuditStats *out)
{
    int repaired=0;
    for(size_t c=0;c<run->n_charts;c++)if(run->charts[c].placed_uv)repaired=1;
    if(repaired){
        /* Export does not query contacts; the following independent audit
         * constructs and walks the full tree exactly once. */
        Arena_T arena=Arena_new();AsmField *field=AsmField_new_local(arena,run,NULL,NULL);
        int rc=field?AsmField_write(field,dir):-1;
        Arena_dispose(&arena);if(rc)return rc;
    }
    return AsmAudit_run(run,dir,out);
}

static int aa_sidecar(const char *dir, const char *stem, const char *suffix, const int32_t *expected, size_t n)
{
    char path[2048]; snprintf(path,sizeof path,"%s/%s_%s.i32",dir,stem,suffix);
    FILE *fp = fopen(path,"rb"); if (!fp) return -1; int ok = 1; int32_t data[8192];
    for (size_t at = 0; at < n;) {
        size_t take = n-at > 8192 ? 8192 : n-at;
        if (fread(data,sizeof(int32_t),take,fp) != take || memcmp(data,expected+at,take*sizeof(int32_t))) { ok = 0; break; }
        at += take;
    }
    if (fgetc(fp) != EOF || ferror(fp)) ok = 0; if (fclose(fp)) ok = 0; return ok ? 0 : -1;
}

int AsmAudit_encoded(Arena_T arena, const char *dir, const char *stem,
                      const float *xyz, const double *uv, size_t nv,
                      const int32_t *faces, size_t nf, const int32_t *chart,
                      const int32_t *component, const int32_t *cube)
{
    Arena_Mark mark = Arena_save(arena); char path[2048]; MeshBinData actual = {0}; int rc = -1;
    snprintf(path,sizeof path,"%s/%s_roundtrip.json",dir,stem); FILE *report = fopen(path,"wb"); if (!report) return -1;
    snprintf(path,sizeof path,"%s/%s.vmesh",dir,stem);
    int exact = MeshBin_read_precise_arena(arena,path,&actual) == 0 && actual.nv == nv && actual.nf == nf && actual.uv64 &&
        !memcmp(xyz,actual.verts,3*nv*sizeof(float)) && !memcmp(uv,actual.uv64,2*nv*sizeof(double)) && !memcmp(faces,actual.faces,3*nf*sizeof(int32_t));
    exact = exact && !aa_sidecar(dir,stem,"chart",chart,nv) && !aa_sidecar(dir,stem,"component",component,nv) &&
        !aa_sidecar(dir,stem,"material_identity",component,nv) && !aa_sidecar(dir,stem,"cube",cube,nv);
    AsmContactStats contacts = {0}; char digest[65] = ""; int walk_complete = 0;
    if (exact) {
        AsmContacts *q = AsmContacts_new(arena,nv,actual.faces,nf,actual.uv64,AsmContacts_material_tolerance());
        if (q && strcmp(stem,"sheet_extras") == 0) {
            int32_t *frames = ARENA_ALLOC(arena,nf*sizeof(int32_t));
            for (size_t t = 0; t < nf; t++) {
                frames[t] = component[actual.faces[3*t]];
                if (frames[t] != component[actual.faces[3*t+1]] || frames[t] != component[actual.faces[3*t+2]]) { exact = 0; break; }
            }
            AsmContacts_separate_frames(q,frames);
        }
        if (exact && q) {
            /* the exact readback is what this audit certifies; a contact walk that stops at its
             * budget (a stacked sheet) leaves the overlaps UNMEASURED, recorded below */
            int contact_rc = AsmContacts_measure(q,actual.uv64,NULL,NULL,NULL,&contacts);
            if (contact_rc == 0 || contact_rc == -2) { aa_digest(actual.uv64,2*nv*sizeof(double),digest); rc = 0; walk_complete = contact_rc == 0; }
            if (contact_rc == -2) { contacts.pairs = 0; contacts.area = 0; }
        }
    }
    if (fprintf(report,"{\"audit_complete\":%s,\"exact_mesh_and_provenance\":%s,\"vertices\":%zu,\"faces\":%zu,\"contact_walk_complete\":%s,\"overlapping_triangle_pairs\":%zu,\"overlap_area_sum\":%.17g,\"encoded_uv_sha256\":\"%s\",\"coordinate_scope\":\"%s\"}\n",
        rc == 0 ? "true":"false",exact ? "true":"false",nv,nf,walk_complete ? "true":"false",contacts.pairs,contacts.area,digest,
        strcmp(stem,"sheet_extras") == 0 ? "independent unregistered component frames":"shared global sheet frame") < 0) rc = -1;
    if (fclose(report)) rc = -1; Arena_restore(arena,mark); return rc;
}

int AsmAudit_selftest(void)
{
    int fails = 0; char digest[65];
    aa_digest("abc",3,digest); if (strcmp(digest,"ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")) fails++;
    aa_digest("",0,digest); if (strcmp(digest,"e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")) fails++;
    int32_t a[] = {2,7,5}, b[] = {7,5,2}, c[] = {7,2,5}; AaFace ka,kb,kc;
    aa_key(a,ka.v); aa_key(b,kb.v); aa_key(c,kc.v);
    if (aa_compare(&ka,&kb) || !aa_compare(&ka,&kc)) fails++;
    /* Real writer/readback plus original-source ledger negative controls. */
    char dir[1024], path[1200]; snprintf(dir,sizeof dir,"output/asm_native_audit_selftest_%d",ves_getpid()); ves_mkdir(dir);
    snprintf(path,sizeof path,"%s/source.vmesh",dir);
    float xyz[] = {0,0,0,1,0,0,0,1,0, 3,0,0,4,0,0,3,1,0};
    float local_uv[] = {0,0,1,0,0,1}; int32_t faces[] = {0,1,2,3,4,5}, local_faces[] = {0,1,2}; int32_t vid[] = {0,1,2,3,4,5};
    /* Duplicate ownership must also exceed the fixed physical contact band. */
    for(size_t k=0;k<sizeof xyz/sizeof *xyz;k++)xyz[k]*=8;
    for(size_t k=0;k<sizeof local_uv/sizeof *local_uv;k++)local_uv[k]*=8;
    Arena_T arena = Arena_new(); AsmRun run = {0}; AsmChart ch[3] = {{0}}; MeshPileEntry pile = {0};
    snprintf(pile.path,sizeof pile.path,"%s",path); run.arena=arena;run.charts=ch;run.n_charts=2;run.pile=&pile;run.n_cubes=1;
    for (int side=0;side<2;side++){ch[side].id=side;ch[side].parent=-1;ch[side].cube=0;ch[side].nv=3;ch[side].nf=1;ch[side].faces=local_faces;
        ch[side].xyz=xyz+9*side;ch[side].vid=vid+3*side;ch[side].uv=local_uv;ch[side].pose_x=24*side;ch[side].placed=1;ch[side].placement_state=ASM_PLACE_ROOT;}
    AsmAuditStats st;
    if (MeshBin_write(path,xyz,6,faces,2,NULL) || AsmAudit_run(&run,dir,&st) || !st.source_preserved || st.represented_faces!=2 ||
        st.invalid_faces || st.overlapping_pairs || st.continuity_components!=2 || !st.geometry_qualified ||
        st.legacy_geometry_qualified || st.source_regions!=2 || st.source_region_connected_area!=st.represented_area) fails++;
    /* Border: two unstitched triangles of legs 8 have their whole perimeter
     * as border; stitched arc comes off each chart's boundary, never below 0. */
    double perimeter = 16+8*sqrt(2.);
    if (fabs(st.boundary_length-2*perimeter) > 1e-9 || st.border_length != st.boundary_length || st.stitched_length != 0) fails++;
    {
        Arena_T border_arena = Arena_new(); AsmField *bf = AsmField_new(border_arena,&run);
        double stitch[2] = {10,1000}; AsmAuditStats bs; memset(&bs,0,sizeof bs);
        if (!bf || aa_border(border_arena,&run,bf,stitch,&bs,NULL) || fabs(bs.border_length-(perimeter-10)) > 1e-9 ||
            fabs(bs.stitched_length-(10+perimeter)) > 1e-9) fails++;
        Arena_dispose(&border_arena);
    }
    AsmRelation layout = {0}; layout.a=0; layout.b=1; layout.flags=ASM_REL_LAYOUT; layout.corr_first=-1;
    run.rels=&layout; run.n_rels=1;
    if (AsmAudit_run(&run,dir,&st) || !st.complete || !st.source_preserved || st.represented_faces!=2 ||
        st.obligations!=1 || st.unresolved_seams!=1 || st.passing_seams || st.continuity_components!=2 || !st.geometry_qualified ||
        st.unresolved_present_repair_seams || st.legacy_geometry_qualified) fails++;
    /* A required connection with no usable witnesses cannot be ignored just
     * because its numerical seam is incomplete. Both charts are present. */
    layout.continuity=ASM_CONT_SOURCE|ASM_CONT_REQUIRED;
    if(AsmAudit_run(&run,dir,&st) || st.geometry_qualified || st.source_regions!=1 ||
       st.fractured_source_regions!=1 || st.unresolved_present_repair_seams!=1 || st.unresolved_missing_repair_seams ||
       st.source_region_connected_area!=.5*st.represented_area)fails++;
    /* The same seam recorded as a tear: no obligation, two regions across it,
     * ledgered, and still one fractured region without it. */
    layout.flags=ASM_REL_LAYOUT|ASM_REL_TORN;
    if(AsmAudit_run(&run,dir,&st) || !st.geometry_qualified || st.obligations || st.torn_seams!=1 || st.torn_between_pieces!=1 ||
       st.torn_cuts || st.source_regions!=2 || st.fractured_source_regions || st.untorn_source_regions!=1 ||
       st.untorn_fractured_source_regions!=1 || st.untorn_source_region_connected_area!=.5*st.represented_area)fails++;
    layout.flags=ASM_REL_LAYOUT;
    ch[1].placement_state=ASM_PLACE_NONE;
    if(AsmAudit_run(&run,dir,&st) || st.geometry_qualified || st.unresolved_present_repair_seams ||
       st.unresolved_missing_repair_seams!=1 || st.represented_faces!=1)fails++;
    ch[1].placement_state=ASM_PLACE_ROOT;
    run.rels=NULL; run.n_rels=0;
    AsmField *checkpoint=AsmField_new(arena,&run);
    if(!checkpoint || AsmField_write(checkpoint,dir) || AsmField_verify(checkpoint,dir))fails++;
    else {
        checkpoint->uv[0]+=.125;if(!AsmField_verify(checkpoint,dir))fails++;checkpoint->uv[0]-=.125;
        snprintf(path,sizeof path,"%s/stage6_repaired_source_vertex.i32",dir);FILE *bad=fopen(path,"r+b");
        if(!bad)fails++;else{int32_t wrong=999;fwrite(&wrong,sizeof wrong,1,bad);fclose(bad);if(!AsmField_verify(checkpoint,dir))fails++;}
        /* A checkpoint has authoritative double coordinates but no required
         * pre-existing export beside it. Its diagnostic explicitly exports
         * that state, while the ordinary published-file audit still rejects
         * the corrupted sidecar. */
        if(AsmField_commit(checkpoint,&run) || !AsmAudit_run(&run,dir,&st) ||
           AsmAudit_checkpoint(&run,dir,&st) || !st.complete || !st.source_preserved ||
           st.represented_faces!=2 || st.invalid_faces || st.overlapping_pairs)fails++;
        for(size_t i=0;i<run.n_charts;i++)ch[i].placed_uv=NULL;
    }
    run.n_charts=1;
    if(AsmAudit_run(&run,dir,&st) || st.source_faces!=2 || st.absent_faces!=1 || st.source_preserved || st.geometry_qualified)fails++;
    run.n_charts=3;ch[2]=ch[0];ch[2].id=2;
    if(AsmAudit_run(&run,dir,&st) || st.multiply_owned_faces!=1 || st.source_preserved || !st.overlapping_pairs)fails++;
    run.n_charts=2;ch[1].flags=ASM_CHART_MIRROR;
    if(AsmAudit_run(&run,dir,&st) || st.invalid_faces || st.strict_bad_faces)fails++;
    ch[1].flags=0;float altered[9];memcpy(altered,xyz,sizeof altered);altered[0]+=.1f;ch[0].xyz=altered;
    if(AsmAudit_run(&run,dir,&st) || st.altered_vertices!=1 || st.source_preserved)fails++;
    ch[0].xyz=xyz;ch[0].placement_state=ch[1].placement_state=ASM_PLACE_NONE;
    AsmRelation relation={0};relation.a=0;relation.b=1;relation.continuity=ASM_CONT_SOURCE;run.rels=&relation;run.n_rels=1;
    if(AsmAudit_run(&run,dir,&st) || !st.complete || st.faces || st.represented_faces || st.obligations!=1 || st.unresolved_seams!=1 || st.geometry_qualified)fails++;
    /* Both charts still have local poses. The emitted ledger must agree
     * with the global field's empty registration, not the legacy pose flag. */
    snprintf(path,sizeof path,"%s/stage7_charts.csv",dir); FILE *ledger=fopen(path,"rb");
    if(!ledger)fails++;
    else {
        char row[1024]; int chart,cube,placed; unsigned flags;
        if(!fgets(row,sizeof row,ledger))fails++;
        for(int i=0;i<2;i++)
            if(!fgets(row,sizeof row,ledger) || sscanf(row,"%d,%d,%u,%d,",&chart,&cube,&flags,&placed)!=4 || chart!=i || placed)fails++;
        if(fgets(row,sizeof row,ledger))fails++;
        fclose(ledger);
    }
    /* Coverage uses original physical area, not face count or represented
     * area as its denominator. One of two source faces carries exactly 85%
     * of the area; missing the other face is disclosed without losing it. */
    {
        float source[]={0,0,0,17,0,0,0,10,0,30,0,0,33,0,0,30,10,0};
        float cover_uv[]={0,0,17,0,0,10};
        run.n_rels=0;run.rels=NULL;ch[0].xyz=source;ch[1].xyz=source+9;
        ch[0].uv=cover_uv;ch[0].placement_state=ASM_PLACE_ROOT;
        if(MeshBin_write(pile.path,source,6,faces,2,NULL) || AsmAudit_run(&run,dir,&st) ||
           !st.geometry_qualified || st.legacy_geometry_qualified || !st.source_preserved ||
           st.source_faces!=2 || st.represented_faces!=1 || st.source_area!=100 || st.represented_area!=85)fails++;
        source[3]=16; /* 80/95 area, while the UV stretch remains below 10%. */
        if(MeshBin_write(pile.path,source,6,faces,2,NULL) || AsmAudit_run(&run,dir,&st) ||
           st.geometry_qualified || !st.source_preserved || st.invalid_faces || st.strict_bad_faces ||
           st.source_area!=95 || st.represented_area!=80)fails++;
    }
    Arena_dispose(&arena);
    fprintf(stderr,"  independent native audit primitives: %s (%d failures)\n",fails ? "FAIL":"ok",fails); return fails;
}
