#include "asm_field.h"
#include "asm_boundary.h"
#include "asm_continuity.h"
#include "../common/mesh_bin.h"
#include "../common/pipeline_constants.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int aff_obligation(const AsmRelation *r)
{
    if (r->flags & ASM_REL_TORN) return 0; /* a recorded tear */
    if ((r->flags & ASM_REL_CONTACT) && !(r->continuity & ASM_CONT_INTRACHART_CUT)) return 0;
    return AsmRel_is_join(r) || (r->continuity & (ASM_CONT_SOURCE | ASM_CONT_CUT | ASM_CONT_VERIFIED | ASM_CONT_INCUMBENT | ASM_CONT_REQUIRED | ASM_CONT_INTRACHART_CUT));
}

int AsmField_relation_required(const AsmRelation *r)
{
    if (!r || !aff_obligation(r)) return 0;
    return !!(r->continuity & (ASM_CONT_REQUIRED | ASM_CONT_CUT | ASM_CONT_INTRACHART_CUT)) ||
           (AsmRel_is_join(r) && !!(r->continuity & (ASM_CONT_SOURCE | ASM_CONT_CUT | ASM_CONT_VERIFIED)));
}

int AsmField_relation_torn(const AsmRelation *r)
{
    if (!r || !(r->flags & ASM_REL_TORN)) return 0;
    AsmRelation untorn = *r; untorn.flags &= ~(uint32_t)ASM_REL_TORN;
    return AsmField_relation_required(&untorn);
}

#include "asm_field_cut.inc"

static int aff_trim(const AsmChart *c, int32_t vertex, const AsmTrim *trim,
                     const float gap[2], double sign, size_t offset, int32_t ids[3], double w[3])
{
    if (!trim->valid || trim->face < 0 || (size_t)trim->face >= c->nf ||
        vertex < 0 || (size_t)vertex >= c->nv) return 0;
    const int32_t *f = c->faces+3*(size_t)trim->face;
    double a = sign*(trim->inverse[0]*gap[0]+trim->inverse[1]*gap[1]);
    double b = sign*(trim->inverse[2]*gap[0]+trim->inverse[3]*gap[1]);
    w[0] = -a-b; w[1] = a; w[2] = b;
    int found = 0;
    for (int k = 0; k < 3; k++) {
        if (f[k] < 0 || (size_t)f[k] >= c->nv || !isfinite(w[k])) return 0;
        ids[k] = (int32_t)(offset+(size_t)f[k]);
        if (f[k] == vertex) { w[k] += 1; found++; }
    }
    return found == 1 && fabs(w[0]+w[1]+w[2]-1) <= 1e-10;
}

int AsmField_material_map(Arena_T arena, const AsmRun *run, const size_t *offset, size_t nv, int32_t **material)
{
    *material = NULL; if (nv > INT32_MAX) return -1;
    int32_t *parent = NULL;
    for (size_t i = 0; i < run->n_rels; i++) {
        const AsmRelation *r = run->rels+i;
        if (!(r->continuity & ASM_CONT_INTRACHART_CUT)) continue;
        if (r->a < 0 || r->a != r->b || (size_t)r->a >= run->n_charts) return -1;
        size_t base = offset[r->a]; if (base == SIZE_MAX) continue;
        const AsmChart *ch = run->charts+r->a;
        if (base > nv || ch->nv > nv-base || r->corr_first < 0 || r->corr_count != 2 ||
            (size_t)r->corr_first > run->n_corr || 2 > run->n_corr-(size_t)r->corr_first) return -1;
        const AsmCorr *corr = run->corr+r->corr_first; double weights[2];
        if (!aff_cut_weights(ch,r,corr,weights)) return -1;
        if (!parent) { parent = ARENA_ALLOC(arena,nv*sizeof(int32_t)); for (size_t v = 0; v < nv; v++) parent[v] = (int32_t)v; }
        for (int k = 0; k < 2; k++) {
            int32_t ids[3]; double coefficients[3];
            if (!aff_trim(ch,corr[k].va,&corr[k].trim_a,corr[k].gap_a,1,base,ids,coefficients) ||
                !aff_trim(ch,corr[k].vb,&corr[k].trim_b,corr[k].gap_b,-1,base,ids,coefficients)) return -1;
            int32_t a = aff_material_root(parent,(int32_t)(base+(size_t)corr[k].va));
            int32_t b = aff_material_root(parent,(int32_t)(base+(size_t)corr[k].vb)); parent[a>b?a:b] = a<b?a:b;
        }
    }
    if (parent) for (size_t v = 0; v < nv; v++) parent[v] = aff_material_root(parent,(int32_t)v);
    *material = parent; return 0;
}

/* How many workers a field's contact measurement may use.  Every field built
 * on this path wants the same answer -- the repair's admission builds one per
 * attempt -- and the walk's result does not depend on the count, so this is a
 * configured default rather than per-field state. */
static int af_threads = 1;

int AsmField_threads(int threads) { int previous = af_threads; af_threads = threads > 1 ? threads : 1; return previous; }

static AsmField *af_build(Arena_T arena, const AsmRun *run, const uint8_t *member, const AsmField *parent, int build_contacts, int seams_only)
{
    if (!arena || !run || !run->n_charts || run->n_charts > INT32_MAX) return NULL;
    AsmField *f = ARENA_CALLOC(arena,1,sizeof *f); f->arena = arena; f->run = run;
    f->vertex_offset = ARENA_ALLOC(arena,run->n_charts*sizeof(size_t));
    f->face_offset = ARENA_ALLOC(arena,run->n_charts*sizeof(size_t));
    for (size_t i = 0; i < run->n_charts; i++) {
        const AsmChart *c = run->charts+i;
        f->vertex_offset[i] = f->face_offset[i] = SIZE_MAX;
        if (!AsmChart_registered(c) || (member && !member[i])) continue;
        if (c->id != (int32_t)i || !c->nv || !c->nf || c->nv > INT32_MAX-f->nv || c->nf > SIZE_MAX/96-f->nf) return NULL;
        f->vertex_offset[i] = f->nv; f->face_offset[i] = f->nf;
        f->nv += c->nv; f->nf += c->nf; f->charts++;
    }
    size_t alloc_nv = f->nv ? f->nv : 1, alloc_nf = f->nf ? f->nf : 1;
    f->uv = ARENA_ALLOC(arena,2*alloc_nv*sizeof(double));
    if (!seams_only) {
        f->xyz = ARENA_ALLOC(arena,3*alloc_nv*sizeof(float));
        f->faces = ARENA_ALLOC(arena,3*alloc_nf*sizeof(int32_t));
        f->face_chart = ARENA_ALLOC(arena,alloc_nf*sizeof(int32_t));
        f->vertex_chart = ARENA_ALLOC(arena,alloc_nv*sizeof(int32_t));
        f->mass = ARENA_CALLOC(arena,alloc_nv,sizeof(double));
        f->metric = ARENA_ALLOC(arena,alloc_nf*sizeof(AsmMetricFace));
    }
    /* Every chart owns disjoint vertex and face ranges and its faces index
     * only its own vertices, so charts are copied, measured and weighed in
     * parallel with the serial per-chart order of every sum (a whole
     * 21x21x21 field is 350 M faces; one thread took ~150 s). */
    int nc_int = (int)run->n_charts, ci = 0, bad = 0;
#pragma omp parallel for schedule(dynamic,16) num_threads(af_threads) reduction(|:bad) if(f->nf >= 1000000)
    for (ci = 0; ci < nc_int; ci++) {
        size_t i = (size_t)ci;
        if (f->vertex_offset[i] == SIZE_MAX) continue;
        const AsmChart *c = run->charts+i; size_t v0 = f->vertex_offset[i], f0 = f->face_offset[i];
        int chart_bad = 0;
        if (!seams_only) memcpy(f->xyz+3*v0,c->xyz,3*c->nv*sizeof(float));
        /* the parent field's CURRENT coordinates where it holds this chart (a repair in progress
         * owns them; AsmChart_point would read whatever was last committed) */
        const double *from = parent && parent->vertex_offset[i] != SIZE_MAX ? parent->uv+2*parent->vertex_offset[i] : NULL;
        if (from) memcpy(f->uv+2*v0,from,2*c->nv*sizeof(double));
        for (size_t v = 0; v < c->nv; v++) {
            if (!seams_only) f->vertex_chart[v0+v] = (int32_t)i;
            if (!from) AsmChart_point(c,v,f->uv+2*(v0+v));
            if (!isfinite(f->uv[2*(v0+v)]) || !isfinite(f->uv[2*(v0+v)+1])) chart_bad = 1;
        }
        for (size_t t = 0; !seams_only && t < c->nf; t++) {
            f->face_chart[f0+t] = (int32_t)i;
            for (int k = 0; k < 3; k++) {
                int source = (c->flags & ASM_CHART_MIRROR) && k ? 3-k : k;
                int32_t v = c->faces[3*t+source];
                if (v < 0 || (size_t)v >= c->nv) { chart_bad = 1; v = 0; }
                f->faces[3*(f0+t)+k] = (int32_t)(v0+(size_t)v);
            }
        }
        if (!chart_bad && !seams_only) {
            AsmChart view = {0}; view.nv = f->nv; view.nf = c->nf; view.xyz = f->xyz; view.faces = f->faces+3*f0;
            if (AsmMetric_prepare(&view,f->metric+f0)) chart_bad = 1;
            else for (size_t t = f0; t < f0+c->nf; t++) for (int k = 0; k < 3; k++)
                f->mass[f->faces[3*t+k]] += f->metric[t].area/3;
        }
        bad |= chart_bad;
    }
    if (bad) return NULL;
    /* In a LOCAL field a relation whose two endpoints are both absent is not an obligation of
     * this field at all; on the whole-field path the rule is unchanged (an obligation between
     * two unregistered charts is still counted, and the audit's digest depends on that). */
    #define AFF_LOCAL_DROP(r) (member && (size_t)(r)->a < run->n_charts && (size_t)(r)->b < run->n_charts && \
                               f->vertex_offset[(r)->a] == SIZE_MAX && f->vertex_offset[(r)->b] == SIZE_MAX)
    for (size_t i = 0; i < run->n_rels; i++) if (aff_obligation(run->rels+i) && !AFF_LOCAL_DROP(run->rels+i)) {
        const AsmRelation *r = run->rels+i;
        /* Layout proposals have no measured boundary samples and use -1 as
         * their empty range's sentinel. Retain the unresolved obligation;
         * only a nonempty range needs a valid offset and correspondence data. */
        if (r->corr_count < 0 || (r->corr_count > 0 &&
            (!run->corr || r->corr_first < 0 || (size_t)r->corr_first > run->n_corr ||
             (size_t)r->corr_count > run->n_corr-(size_t)r->corr_first))) return NULL;
        if ((size_t)r->corr_count > SIZE_MAX/sizeof(AsmWitness)-f->n_witness) return NULL;
        f->n_seams++; f->n_witness += (size_t)r->corr_count;
    }
    f->seams = ARENA_CALLOC(arena,f->n_seams ? f->n_seams : 1,sizeof(AsmFieldSeam));
    f->witness = ARENA_CALLOC(arena,f->n_witness ? f->n_witness : 1,sizeof(AsmWitness));
    const AsmChart *source_charts = aff_source_charts(f);
    size_t ns = 0, nw = 0;
    /* Each chart's boundary loops are built once for all its seams. */
    Arena_T loops_arena = Arena_new_sized((size_t)1<<20);
    AsmBoundaryLoops **loops = ARENA_CALLOC(loops_arena,run->n_charts ? run->n_charts : 1,sizeof *loops);
    for (size_t i = 0; i < run->n_rels; i++) if (aff_obligation(run->rels+i) && !AFF_LOCAL_DROP(run->rels+i)) {
        const AsmRelation *r = run->rels+i;
        AsmFieldSeam *s = f->seams+ns++;
        s->relation = i; s->a = r->a; s->b = r->b; s->first = nw; s->count = (size_t)r->corr_count; nw += s->count;
        s->original_cut = !!(r->continuity & ASM_CONT_INTRACHART_CUT);
        if (r->a < 0 || r->b < 0 || (size_t)r->a >= run->n_charts || (size_t)r->b >= run->n_charts ||
            (r->a == r->b && !s->original_cut)) { Arena_dispose(&loops_arena); return NULL; }
        const AsmChart *a = run->charts+r->a, *b = run->charts+r->b;
        size_t oa = f->vertex_offset[r->a], ob = f->vertex_offset[r->b];
        s->parity = !!((a->flags ^ b->flags) & ASM_CHART_MIRROR) == !!(r->flags & ASM_REL_PARITY);
        s->source_run = !!(r->continuity & (ASM_CONT_SOURCE | ASM_CONT_CUT | ASM_CONT_VERIFIED));
        s->required = AsmField_relation_required(r);
        if (!s->count || oa == SIZE_MAX || ob == SIZE_MAX || !s->parity) continue;
        Arena_Mark mark = Arena_save(arena);
        double *weights = ARENA_ALLOC(arena,s->count*sizeof(double));
        const AsmCorr *corr = run->corr+r->corr_first;
        s->complete = s->original_cut ? aff_cut_weights(a,r,corr,weights) : aff_boundary_weights(f,source_charts,loops_arena,loops,r,corr,s->count,weights) == 0;
        if (s->complete) for (size_t k = 0; k < s->count; k++) {
            const AsmCorr *c = corr+k; AsmWitness *w = f->witness+s->first+k;
            w->weight = weights[k];
            if (!aff_trim(a,c->va,&c->trim_a,c->gap_a,1,oa,w->vertex[0],w->coefficient[0]) ||
                !aff_trim(b,c->vb,&c->trim_b,c->gap_b,-1,ob,w->vertex[1],w->coefficient[1])) { s->complete = 0; continue; }
            w->vertex[0][3] = (int32_t)(ob+(size_t)c->vb);
            w->vertex[1][3] = (int32_t)(oa+(size_t)c->va);
            for (int d = 0; d < 2; d++) w->coefficient[d][3] = -(w->coefficient[d][0]+w->coefficient[d][1]+w->coefficient[d][2]);
            w->valid = 1;
            w->required = s->required;
        }
        Arena_restore(arena,mark);
    }
    Arena_dispose(&loops_arena);
    #undef AFF_LOCAL_DROP
    if (!build_contacts) return f;
    return AsmField_build_contacts(f) ? NULL : f;
}

AsmField *AsmField_new(Arena_T arena, const AsmRun *run)
{ return af_build(arena,run,NULL,NULL,1,0); }

AsmField *AsmField_new_local(Arena_T arena, const AsmRun *run, const uint8_t *member, const AsmField *parent)
{ return af_build(arena,run,member,parent,0,0); }

int AsmField_build_contacts(AsmField *f)
{
    if (!f) return -1;
    if (!f->nf) return 0;
    int32_t *material = NULL;
    if (!aff_material_vertices(f,&material)) return -1;
    f->contacts = AsmContacts_new_source(f->arena,f->nv,f->faces,f->nf,f->uv,AsmContacts_material_tolerance(),af_threads,material);
    if (!f->contacts) return -1;
    AsmContacts_threads(f->contacts,af_threads);
    return 0;
}

int AsmField_commit(AsmField *field, AsmRun *run)
{
    if (!field || field->run != run) return -1;
    for (size_t c = 0; c < run->n_charts; c++) if (field->vertex_offset[c] != SIZE_MAX)
        run->charts[c].placed_uv = field->uv+2*field->vertex_offset[c];
    return 0;
}

static int aff_labels(const AsmField *f, const char *dir, int verify)
{
    const char *names[] = {"chart","cube","source_vertex","parity"};
    for (int type = 0; type < 4; type++) {
        char path[2048]; snprintf(path,sizeof path,"%s/stage6_repaired_%s.i32",dir,names[type]);
        FILE *fp = fopen(path,verify ? "rb":"wb"); if (!fp) return -1; int ok = 1;
        int32_t expected[8192], actual[8192];
        for (size_t c = 0; c < f->run->n_charts && ok; c++) if (f->vertex_offset[c] != SIZE_MAX) {
            const AsmChart *chart = f->run->charts+c; if (!chart->vid) { ok = 0; break; }
            for (size_t start = 0; start < chart->nv && ok;) {
                size_t take = chart->nv-start > 8192 ? 8192 : chart->nv-start;
                for (size_t k = 0; k < take; k++) expected[k] = type == 0 ? chart->id : type == 1 ? chart->cube :
                    type == 2 ? chart->vid[start+k] : (chart->flags & ASM_CHART_MIRROR) ? -1 : 1;
                if (verify) ok = fread(actual,sizeof(int32_t),take,fp) == take && !memcmp(actual,expected,take*sizeof(int32_t));
                else ok = fwrite(expected,sizeof(int32_t),take,fp) == take;
                start += take;
            }
        }
        if (verify && (fgetc(fp) != EOF || ferror(fp))) ok = 0;
        if (fclose(fp)) ok = 0; if (!ok) return -1;
    }
    return 0;
}

int AsmField_write(const AsmField *f, const char *dir)
{
    if (!f->nv) return 0;
    Arena_T arena = Arena_new(); char path[2048]; snprintf(path,sizeof path,"%s/stage6_repaired.vmesh",dir);
    MeshBinStream_T writer = NULL; int rc = MeshBin_stream_open_uv64(arena,path,f->nv,f->nf,&writer);
    if (rc) { Arena_dispose(&arena); return -1; }
    for (size_t c = 0; c < f->run->n_charts && !rc; c++) if (f->vertex_offset[c] != SIZE_MAX) {
        const AsmChart *chart = f->run->charts+c; size_t offset = f->vertex_offset[c];
        MeshBinData chunk = {0}; chunk.nv = chart->nv; chunk.nf = chart->nf; chunk.verts = f->xyz+3*offset;
        chunk.uv64 = f->uv+2*offset; chunk.faces = chart->faces; rc = MeshBin_stream_append(writer,&chunk);
    }
    if (MeshBin_stream_close(&writer,rc == 0)) rc = -1;
    if (!rc) rc = aff_labels(f,dir,0);
    Arena_dispose(&arena); return rc;
}

int AsmField_verify(const AsmField *f, const char *dir)
{
    if (!f->nv) return 0;
    Arena_T arena = Arena_new(); MeshBinData mesh = {0}; char path[2048]; int rc = -1;
    snprintf(path,sizeof path,"%s/stage6_repaired.vmesh",dir);
    if (MeshBin_read_precise_arena(arena,path,&mesh) || mesh.nv != f->nv || mesh.nf != f->nf || !mesh.uv64 ||
        memcmp(mesh.verts,f->xyz,3*f->nv*sizeof(float)) || memcmp(mesh.uv64,f->uv,2*f->nv*sizeof(double))) goto done;
    for (size_t c = 0; c < f->run->n_charts; c++) if (f->vertex_offset[c] != SIZE_MAX) {
        const AsmChart *chart = f->run->charts+c; size_t offset = f->vertex_offset[c], first = f->face_offset[c];
        for (size_t k = 0; k < 3*chart->nf; k++) if (mesh.faces[3*first+k] != chart->faces[k]+(int32_t)offset) goto done;
    }
    rc = aff_labels(f,dir,1);
done:
    Arena_dispose(&arena); return rc;
}

void AsmField_face(const AsmField *f, size_t t, const double *uv,
                    double j[4], double *lo, double *hi, double *det)
{
    const int32_t *v = f->faces+3*t; const AsmMetricFace *r = f->metric+t;
    double u = uv[2*(size_t)v[0]], w = uv[2*(size_t)v[0]+1];
    j[0] = (uv[2*(size_t)v[1]]-u)/r->length;
    j[2] = (uv[2*(size_t)v[1]+1]-w)/r->length;
    j[1] = (uv[2*(size_t)v[2]]-u-r->along*j[0])/r->height;
    j[3] = (uv[2*(size_t)v[2]+1]-w-r->along*j[2])/r->height;
    *det = j[0]*j[3]-j[1]*j[2];
    *hi = .5*(hypot(j[0]+j[3],j[2]-j[1])+hypot(j[0]-j[3],j[2]+j[1]));
    *lo = *hi > 0 ? fabs(*det)/ *hi : 0;
}

void AsmField_residual(const AsmWitness *w, int direction, const double *uv, double r[2])
{
    r[0] = r[1] = 0;
    /* Subtract a common point before multiplying large extrapolation weights.
     * Rows are translation invariant to a strict absolute sum tolerance. */
    size_t origin = (size_t)w->vertex[direction][3];
    for (int k = 0; k < 3; k++) for (int d = 0; d < 2; d++)
        r[d] += w->coefficient[direction][k]*(uv[2*(size_t)w->vertex[direction][k]+d]-uv[2*origin+d]);
}

void AsmField_seam(const AsmField *f, size_t seam, const double *uv, AsmSeamMeasure *out)
{
    memset(out,0,sizeof *out); const AsmFieldSeam *s = f->seams+seam;
    for (size_t k = 0; k < s->count; k++) {
        const AsmWitness *w = f->witness+s->first+k;
        if (!w->valid || !(w->weight > 0)) { out->invalid++; continue; }
        double a[2], b[2]; AsmField_residual(w,0,uv,a); AsmField_residual(w,1,uv,b);
        double one = a[0]*a[0]+a[1]*a[1], two = b[0]*b[0]+b[1]*b[1], worst = fmax(one,two);
        if (!isfinite(one) || !isfinite(two)) { out->invalid++; continue; }
        out->mass += w->weight; out->rms += w->weight*worst;
        out->energy += w->weight*.5*(one+two); out->observations++;
        if (worst <= ASM_SEAM_TOLERANCE_VOX*ASM_SEAM_TOLERANCE_VOX) out->support += w->weight;
    }
    out->rms = out->mass > 0 && !out->invalid ? sqrt(out->rms/out->mass) : INFINITY;
    if (s->original_cut)
        out->pass = s->complete && s->parity && s->source_run && !out->invalid && out->observations == 2 &&
                    out->mass > 0 && out->support == out->mass && out->rms <= ASM_SEAM_TOLERANCE_VOX;
    else
        out->pass = s->complete && s->parity && s->source_run && !out->invalid && out->observations >= 6 &&
                    out->mass >= 8 && out->support >= .9*out->mass && out->rms <= ASM_SEAM_TOLERANCE_VOX;
}

int AsmField_measure_source_seams(const AsmRun *run, AsmSeamMeasure *out)
{
    if (!run || (run->n_rels && !out)) return -1;
    if (!run->n_rels) return 0;
    memset(out,0,run->n_rels*sizeof *out);
    Arena_T arena=Arena_new();
    AsmField *field=af_build(arena,run,NULL,NULL,0,1);int rc=-1;
    if (field) {
        for (size_t s=0;s<field->n_seams;s++)
            AsmField_seam(field,s,field->uv,out+field->seams[s].relation);
        rc=0;
    }
    Arena_dispose(&arena);return rc;
}

#include "asm_field_cut_test.inc"

static int aff_seam_policy_selftest(void)
{
    int failures = 0;
    for (int test = 0; test < 11; test++) {
        AsmWitness witnesses[10] = {{0}}; double uv[22] = {0};
        AsmFieldSeam seam = {0}; AsmField field = {0}; AsmSeamMeasure measure = {0};
        seam.count = 10; seam.complete = seam.parity = seam.source_run = 1;
        field.seams = &seam; field.n_seams = 1; field.witness = witnesses; field.n_witness = 10;
        for (int k = 0; k < 10; k++) {
            AsmWitness *w = witnesses+k; w->valid = 1; w->weight = 1;
            uv[2*k] = test == 0 ? 3.5 : test == 1 ? 4 : test == 2 ? 4.01 : 0;
            for (int side = 0; side < 2; side++) {
                w->vertex[side][0] = k; w->vertex[side][3] = 10;
                w->coefficient[side][0] = 1; w->coefficient[side][3] = -1;
            }
        }
        /* Ninety percent support alone cannot excuse excessive RMS, and
         * low RMS alone cannot excuse insufficient boundary support. */
        if (test == 3) uv[0] = 12;
        if (test == 4) uv[0] = 13;
        if (test == 5) uv[0] = uv[2] = 5;
        if (test == 6) seam.parity = 0;
        if (test == 7) seam.source_run = 0;
        if (test == 8) seam.complete = 0;
        if (test == 9) witnesses[0].valid = 0;
        if (test == 10) for (int k = 0; k < 10; k++) witnesses[k].weight = .5;
        AsmField_seam(&field,0,uv,&measure);
        int expected = test == 0 || test == 1 || test == 3;
        int okay = measure.pass == expected;
        fprintf(stderr,"  native four-voxel seam policy %d: %s\n",test,okay?"ok":"FAIL");
        failures += !okay;
    }
    return failures;
}

int AsmField_selftest(void)
{
    int fails = aff_cut_selftest()+aff_seam_policy_selftest();
    AsmWitness w = {0}; w.vertex[0][0] = 0; w.vertex[0][1] = 1; w.vertex[0][2] = 2; w.vertex[0][3] = 3;
    w.coefficient[0][0] = -1; w.coefficient[0][1] = 1; w.coefficient[0][2] = 1; w.coefficient[0][3] = -1;
    double uv[] = {0,0,1,0,0,1,1,1}, p[2]; AsmField_residual(&w,0,uv,p);
    if (hypot(p[0],p[1]) != 0) fails++;
    for (int i = 0; i < 8; i++) uv[i] += 1e9;
    AsmField_residual(&w,0,uv,p); if (hypot(p[0],p[1]) != 0) fails++;
    {
        /* A raster join must survive the handoff to repair without gaining
         * a physical seam certificate. Malformed nonempty ranges still fail. */
        Arena_T arena = Arena_new();
        AsmRun run = {0}; AsmChart charts[2] = {{0}}; AsmRelation relation = {0};
        run.arena = arena; run.charts = charts; run.n_charts = 2;
        run.rels = &relation; run.n_rels = 1;
        relation.a = 0; relation.b = 1; relation.flags = ASM_REL_LAYOUT; relation.corr_first = -1;
        AsmField *f = AsmField_new(arena,&run); AsmSeamMeasure seam = {0};
        if (!f || f->n_seams != 1 || f->n_witness || f->seams[0].complete || f->seams[0].source_run) fails++;
        else {
            AsmField_seam(f,0,f->uv,&seam);
            if (seam.pass || seam.observations || seam.mass) fails++;
        }
        relation.corr_count = 1;
        if (AsmField_new(arena,&run)) fails++;
        AsmCorr corr = {0}; run.corr = &corr; run.n_corr = 1;
        if (AsmField_new(arena,&run)) fails++;
        relation.corr_first = 1;
        if (AsmField_new(arena,&run)) fails++;
        relation.corr_first = 0; relation.corr_count = -1;
        if (AsmField_new(arena,&run)) fails++;
        Arena_dispose(&arena);
    }
    {
        /* A LOCAL field over every registered chart must be the whole field, field for field; a
         * member mask must drop the charts it leaves out, keep a half-present obligation
         * incomplete, and start from the parent's CURRENT coordinates. */
        Arena_T arena = Arena_new();
        enum { NC = 3, NV = 4, NF = 2 };
        AsmRun run = {0}; run.arena = arena;
        AsmChart *charts = ARENA_CALLOC(arena,NC,sizeof(AsmChart));
        float *xyz = ARENA_ALLOC(arena,NC*3*NV*sizeof(float));
        float *uvc = ARENA_ALLOC(arena,NC*2*NV*sizeof(float));
        int32_t *faces = ARENA_ALLOC(arena,NC*3*NF*sizeof(int32_t));
        for (int c = 0; c < NC; c++) {
            const float quad[8] = {0,0, 10,0, 0,10, 10,10};
            for (int v = 0; v < NV; v++) {
                xyz[c*3*NV+3*v] = quad[2*v]; xyz[c*3*NV+3*v+1] = quad[2*v+1]; xyz[c*3*NV+3*v+2] = (float)(40*c);
                uvc[c*2*NV+2*v] = quad[2*v]; uvc[c*2*NV+2*v+1] = quad[2*v+1];
            }
            const int32_t tri[6] = {0,1,3, 0,3,2};
            memcpy(faces+c*3*NF,tri,sizeof tri);
            charts[c].id = c; charts[c].nv = NV; charts[c].nf = NF; charts[c].placed = 1;
            charts[c].placement_state = ASM_PLACE_SEAM; charts[c].area3d = 100.0;
            charts[c].xyz = xyz+c*3*NV; charts[c].uv = uvc+c*2*NV; charts[c].faces = faces+c*3*NF;
            charts[c].pose_x = 40.0*c;
        }
        run.charts = charts; run.n_charts = NC;
        AsmRelation *rels = ARENA_CALLOC(arena,2,sizeof(AsmRelation));
        rels[0].a = 0; rels[0].b = 1; rels[0].flags = ASM_REL_LAYOUT; rels[0].corr_first = -1;
        rels[1].a = 1; rels[1].b = 2; rels[1].flags = ASM_REL_LAYOUT; rels[1].corr_first = -1;
        run.rels = rels; run.n_rels = 2;
        AsmField *whole = AsmField_new(arena,&run);
        uint8_t all[NC] = {1,1,1}, two[NC] = {1,1,0};
        AsmField *same = AsmField_new_local(arena,&run,all,NULL);
        int ok = whole && same && same->nv == whole->nv && same->nf == whole->nf && same->charts == whole->charts &&
                 same->n_seams == whole->n_seams && same->n_witness == whole->n_witness && !same->contacts &&
                 !memcmp(same->vertex_offset,whole->vertex_offset,NC*sizeof(size_t)) &&
                 !memcmp(same->face_offset,whole->face_offset,NC*sizeof(size_t)) &&
                 !memcmp(same->faces,whole->faces,3*whole->nf*sizeof(int32_t)) &&
                 !memcmp(same->uv,whole->uv,2*whole->nv*sizeof(double)) &&
                 !memcmp(same->xyz,whole->xyz,3*whole->nv*sizeof(float)) &&
                 !memcmp(same->metric,whole->metric,whole->nf*sizeof(AsmMetricFace)) &&
                 !memcmp(same->mass,whole->mass,whole->nv*sizeof(double));
        if (!ok) fails++;
        if (same && AsmField_build_contacts(same)) fails++;
        else if (same && !same->contacts) fails++;
        /* the parent's coordinates, not the charts': move one vertex in the whole field */
        if (whole) whole->uv[0] = -123.0;
        AsmField *local = AsmField_new_local(arena,&run,two,whole);
        int ok2 = local && local->charts == 2 && local->nv == 2*NV && local->vertex_offset[2] == SIZE_MAX &&
                  local->n_seams == 2 && local->uv[0] == -123.0;
        /* relation 1-2 keeps its obligation (one endpoint present) and cannot pass; a field
         * without either endpoint drops it */
        uint8_t one[NC] = {1,0,0};
        AsmField *single = AsmField_new_local(arena,&run,one,whole);
        if (!ok2 || !single || single->n_seams != 1) fails++;
        fprintf(stderr,"  local seam field over a chart subset: %s\n",fails ? "FAIL" : "ok");
        Arena_dispose(&arena);
    }
    /* The fixed policy reaches newly built local/global fields: shallow
     * overlap is accepted and penetration beyond the band is still found. */
    {
        Arena_T arena=Arena_new();
        int32_t tri[6]={0,1,2,3,4,5};double uv[12]={0,0,20,0,0,20,18.5,0,38.5,0,18.5,20};
        AsmField field={0};field.arena=arena;field.nv=6;field.nf=2;field.faces=tri;field.uv=uv;
        AsmContactStats got={0};
        if (AsmField_build_contacts(&field) ||
            AsmContacts_measure(field.contacts,uv,NULL,NULL,NULL,&got) || !got.complete || got.pairs) fails++;
        /* Overlap is a right triangle of leg 9.5: its largest contained
         * disk has diameter (2-sqrt(2))*9.5 > 4 voxels. */
        for (int k=0;k<3;k++) uv[6+2*k]-=8;
        field.contacts=NULL;
        if (AsmField_build_contacts(&field) || AsmContacts_measure(field.contacts,uv,NULL,NULL,NULL,&got) ||
            !got.complete || got.pairs!=1) fails++;
        Arena_dispose(&arena);
    }
    fprintf(stderr,"  immutable native seam field: %s (%d failures)\n",fails ? "FAIL" : "ok",fails);
    return fails;
}
