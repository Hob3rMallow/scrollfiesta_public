/*
 * fiber_flatten -- fibre-aligned re-flattening of a tifxyz canvas (src/flatten/fiber_frame.c; 2026-10-15, the user
 * after PHerc0211 line 11's close-up page measured a 2,089 px arc by the fibres: "Do all the things").
 *
 *   fiber_flatten <job.json>
 *   fiber_flatten --audit-map <job.json>
 *   fiber_flatten --audit-branches <job.json>
 *   fiber_flatten --selftest
 *
 * job.json (paths and the window only; tuning lives in pipeline_constants.h FIBER_FLATTEN_*):
 *   "seg"    the input tifxyz (x.tif y.tif z.tif: float32 planes, -1 = invalid; 1 node per voxel)
 *   "fibers" the fibre field: an uncompressed uint8 zarr group with presence/, nx/, ny/ (fiber_field_infer.py)
 *   "out"    the output directory
 *   "rect"   [x0, y0, x1, y1] the canvas window (optional; default the whole canvas)
 *   "uuid"   the new segment's name (optional; default <seg name>_fiberaligned)
 * Threads: OMP_NUM_THREADS.
 *
 * Writes <out>/<uuid>/{x,y,z,mask}.tif + meta.json (the re-flattened tifxyz, 1 px per voxel of fibre length),
 * <out>/fiber_map_u.tif and fiber_map_v.tif (the forward map at the node grid: float32, NaN off the solved nodes;
 * node (i, j) = canvas row gy0 + i * cell, column gx0 + j * cell; new pixel = (u' - u_min, v' - v_min)),
 * <out>/fiber_profile.csv (a fibre's height across the page, old and new), <out>/fiber_flatten.json (inputs,
 * constants, the audit, timings).
 */
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../common/arena.h"
#include "../common/json_read.h"
#include "../common/pipeline_constants.h"
#include "../common/tiff_io.h"
#include "../common/ves_platform.h"
#include "../common/zarr_u8.h"
#include "../flatten/fiber_frame.h"
#include "../flatten/sparse_solve.h"

static double now_s(void)
{
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

static void stage(const char *name)
{
    printf("[orch] stage %s\n", name);
    fflush(stdout);
}

static int threads_env(void)
{
    const char *s = getenv("OMP_NUM_THREADS");
    int t = s != NULL ? atoi(s) : 1;
    return t > 0 ? t : 1;
}

/* a JSON string literal: quotes, backslashes (Windows paths) and control characters escaped */
static void json_str(FILE *f, const char *s)
{
    fputc('"', f);
    for (; s != NULL && *s != '\0'; s++) {
        unsigned char ch = (unsigned char)*s;
        if (ch == '"' || ch == '\\') {
            fputc('\\', f);
            fputc((int)ch, f);
        } else if (ch < 0x20) {
            fprintf(f, "\\u%04x", (unsigned)ch);
        } else {
            fputc((int)ch, f);
        }
    }
    fputc('"', f);
}

/* the three coordinate planes into one interleaved [H*W*3] array; valid = all three >= 0 and finite */
static int load_tifxyz(const char *dir, float **Xo, uint8_t **vo, int *Ho, int *Wo)
{
    const char *names[3] = { "x.tif", "y.tif", "z.tif" };
    Arena_T io = Arena_new();
    float *pl[3] = { NULL, NULL, NULL }, *X = NULL;
    uint8_t *v = NULL;
    int W = 0, H = 0;
    for (int k = 0; k < 3; k++) {
        char path[1200];
        int w = 0, h = 0;
        snprintf(path, sizeof path, "%s/%s", dir, names[k]);
        if (TiffIO_load_float2d(io, path, &pl[k], &w, &h) != 0 || (k > 0 && (w != W || h != H))) {
            fprintf(stderr, "[fiber_flatten] cannot read %s\n", path);
            Arena_dispose(&io);
            return -1;
        }
        W = w;
        H = h;
    }
    X = (float *)malloc((size_t)W * (size_t)H * 3 * sizeof *X);
    v = (uint8_t *)malloc((size_t)W * (size_t)H);
    if (X == NULL || v == NULL) { free(X); free(v); Arena_dispose(&io); return -1; }
    for (size_t i = 0; i < (size_t)W * (size_t)H; i++) {
        float x = pl[0][i], y = pl[1][i], z = pl[2][i];
        X[i * 3] = x; X[i * 3 + 1] = y; X[i * 3 + 2] = z;
        v[i] = (uint8_t)(x >= 0.0f && y >= 0.0f && z >= 0.0f && isfinite(x) && isfinite(y) && isfinite(z));
    }
    Arena_dispose(&io);
    *Xo = X; *vo = v; *Ho = H; *Wo = W;
    return 0;
}

/* a tifxyz segment from interleaved coordinates (-1 where invalid) and a 0/1 mask */
static int write_tifxyz(const char *dir, const char *uuid, const float *X, const uint8_t *vn, int H, int W)
{
    char path[1200];
    size_t n = (size_t)H * (size_t)W;
    float *pl = (float *)malloc(n * sizeof *pl);
    uint8_t *m = (uint8_t *)malloc(n);
    double lo[3] = { INFINITY, INFINITY, INFINITY }, hi[3] = { -INFINITY, -INFINITY, -INFINITY };
    FILE *f = NULL;
    const char *names[3] = { "x.tif", "y.tif", "z.tif" };
    if (pl == NULL || m == NULL) { free(pl); free(m); return -1; }
    ves_mkdir(dir);
    for (size_t i = 0; i < n; i++) {
        m[i] = vn[i] ? 255 : 0;
        if (vn[i])
            for (int k = 0; k < 3; k++) {
                if (X[i * 3 + (size_t)k] < lo[k]) lo[k] = X[i * 3 + (size_t)k];
                if (X[i * 3 + (size_t)k] > hi[k]) hi[k] = X[i * 3 + (size_t)k];
            }
    }
    for (int k = 0; k < 3; k++) {
        for (size_t i = 0; i < n; i++) pl[i] = vn[i] ? X[i * 3 + (size_t)k] : -1.0f;
        snprintf(path, sizeof path, "%s/%s", dir, names[k]);
        if (TiffIO_save_float2d(path, pl, W, H) != 0) { free(pl); free(m); return -1; }
    }
    snprintf(path, sizeof path, "%s/mask.tif", dir);
    if (TiffIO_save(path, m, 1, H, W) != 0) { free(pl); free(m); return -1; }
    snprintf(path, sizeof path, "%s/meta.json", dir);
    f = fopen(path, "w");
    if (f == NULL) { free(pl); free(m); return -1; }
    fprintf(f, "{\n \"bbox\": [[%.4f, %.4f, %.4f], [%.4f, %.4f, %.4f]],\n \"type\": \"seg\",\n \"uuid\": ",
            isfinite(lo[0]) ? lo[0] : 0.0, isfinite(lo[1]) ? lo[1] : 0.0, isfinite(lo[2]) ? lo[2] : 0.0,
            isfinite(hi[0]) ? hi[0] : 0.0, isfinite(hi[1]) ? hi[1] : 0.0, isfinite(hi[2]) ? hi[2] : 0.0);
    json_str(f, uuid);
    fprintf(f, ",\n \"format\": \"tifxyz\",\n \"scale\": [1.0, 1.0]\n}\n");
    fclose(f);
    free(pl);
    free(m);
    return 0;
}

static int write_map(const char *path, const double *G, int gh, int gw)
{
    size_t n = (size_t)gh * (size_t)gw;
    float *f = (float *)malloc(n * sizeof *f);
    int rc = 0;
    if (f == NULL) return -1;
    for (size_t i = 0; i < n; i++) f[i] = (float)G[i];
    rc = TiffIO_save_float2d(path, f, gw, gh);
    free(f);
    return rc;
}

static void json_number(FILE *f, double value, int digits)
{
    if (isfinite(value)) fprintf(f, "%.*f", digits, value);
    else fputs("null", f);  /* small selections may have no full audit patches */
}

static void jq_values(FILE *f, FiberQuant q, int digits)
{
    fputs("\"p50\": ", f); json_number(f, q.p50, digits);
    fputs(", \"p90\": ", f); json_number(f, q.p90, digits);
    fputs(", \"p99\": ", f); json_number(f, q.p99, digits);
}

static void jq(FILE *f, const char *key, FiberQuant q, const char *tail)
{
    fprintf(f, "  \"%s\": {", key); jq_values(f, q, 4); fprintf(f, "}%s\n", tail);
}

static void jclimb(FILE *f, double range, double bow, double end, double width, double p10, double p90)
{
    fputs("{\"range_px\": ", f); json_number(f, range, 1);
    fputs(", \"bow_after_best_line_px\": ", f); json_number(f, bow, 1);
    fputs(", \"end_to_end_px\": ", f); json_number(f, end, 1);
    fputs(", \"width_px\": ", f); json_number(f, width, 1);
    fputs(", \"column_angle_p10_deg\": ", f); json_number(f, p10, 2);
    fputs(", \"column_angle_p90_deg\": ", f); json_number(f, p90, 2);
    fputc('}', f);
}

/* Audit inputs are either a saved float32 frame or the float64 map written
 * beside an ink_closeup --reuv surface report. This only reads maps. */
static int load_audit_map(Arena_T arena, const JsonValue *spec, FiberMapResult *map)
{
    const char *frame = Json_as_string(Json_object_get(spec,"frame"));
    const char *binary = Json_as_string(Json_object_get(spec,"map_file"));
    const char *report = Json_as_string(Json_object_get(spec,"report")), *error = NULL;
    char path[1400];
    if ((frame != NULL) == (binary != NULL)) return -1;
    if (frame) { snprintf(path,sizeof path,"%s/fiber_flatten.json",frame); report=path; }
    if (!report) return -1;
    const JsonValue *meta = Json_parse_file(arena,report,&error);
    const JsonValue *grid = Json_object_get(meta,frame?"grid":"guarded_map");
    if (!grid) return -1;
    long rows=Json_as_long(Json_object_get(grid,"rows"),0), cols=Json_as_long(Json_object_get(grid,"cols"),0);
    double cell=Json_as_double(Json_object_get(grid,"cell"),NAN);
    double row0=Json_as_double(Json_object_get(grid,"row0"),NAN), col0=Json_as_double(Json_object_get(grid,"col0"),NAN);
    if (rows<2 || cols<2 || rows>INT_MAX || cols>INT_MAX || !(cell>0) ||
        !isfinite(cell) || !isfinite(row0) || !isfinite(col0) ||
        (size_t)rows>SIZE_MAX/(size_t)cols/2/sizeof(double)) return -1;
    size_t count=(size_t)rows*(size_t)cols;
    *map=(FiberMapResult){(int)rows,(int)cols,row0,col0,cell,ARENA_ALLOC(arena,2*count*sizeof(double))};
    if (frame) {
        float *u=NULL,*v=NULL; int uw=0,uh=0,vw=0,vh=0;
        snprintf(path,sizeof path,"%s/fiber_map_u.tif",frame);
        if (TiffIO_load_float2d(arena,path,&u,&uw,&uh)) return -1;
        snprintf(path,sizeof path,"%s/fiber_map_v.tif",frame);
        if (TiffIO_load_float2d(arena,path,&v,&vw,&vh) || uw!=cols || vw!=cols || uh!=rows || vh!=rows) return -1;
        for (size_t i=0;i<count;i++) { map->uv[2*i]=u[i]; map->uv[2*i+1]=v[i]; }
    } else {
        const char *storage=Json_as_string(Json_object_get(grid,"storage"));
        if (!storage || strcmp(storage,"interleaved_float64_uv")) return -1;
        FILE *stream=fopen(binary,"rb");
        if (!stream) return -1;
        int bad=fread(map->uv,sizeof(double),2*count,stream)!=2*count;
        bad |= fgetc(stream)!=EOF || ferror(stream);
        if (fclose(stream)) bad=1;
        if (bad) return -1;
    }
    return 0;
}

static int audit_retained_maps(Arena_T arena, const JsonValue *job, const FiberFrame *ff,
                               const char *out, const char *seg, const char *fibers)
{
    const JsonValue *maps=Json_object_get(job,"audit_maps");
    size_t count=Json_array_len(maps);
    if (!count || count>SIZE_MAX/sizeof(FiberFrameMapStats)) return -1;
    FiberFrameMapStats *stats=ARENA_CALLOC(arena,count,sizeof *stats);
    for (size_t k=0;k<count;k++) {
        const JsonValue *spec=Json_array_get(maps,k);
        FiberMapResult map={0};
        Arena_Mark mark=Arena_save(arena);
        if (!Json_as_string(Json_object_get(spec,"name")) || load_audit_map(arena,spec,&map) ||
            FiberFrame_audit_map(arena,ff,&map,stats+k)) {
            fprintf(stderr,"[fiber_flatten] retained-map audit failed at input %zu\n",k);
            Arena_restore(arena,mark); return -1;
        }
        Arena_restore(arena,mark);
    }
    char path[1400];
    snprintf(path,sizeof path,"%s/fiber_map_quality.json",out);
    FILE *stream=fopen(path,"w");
    if (!stream) return -1;
    fputs("{\n \"schema\": \"fiber-retained-map-quality-v1\",\n \"seg\": ",stream); json_str(stream,seg);
    fputs(",\n \"fibers\": ",stream); json_str(stream,fibers);
    fprintf(stream,",\n \"field_grid\": {\"rows\": %d, \"cols\": %d, \"cell\": %d, \"row0\": %d, \"col0\": %d},\n",
            ff->gh,ff->gw,ff->cell,ff->gy0,ff->gx0);
    fputs(" \"method\": \"Exact bilinear derivatives at measured field nodes; worst incident-cell angle; physical perpendicular arms; line angles modulo 180 degrees; paired before/after quantiles\",\n"
          " \"scope\": \"Sampled fibre alignment only; no whole-mesh orientation, metric, contact or anatomical certificate\",\n \"maps\": [\n",stream);
    for (size_t k=0;k<count;k++) {
        const FiberFrameMapStats *s=stats+k;
        fputs(" {\"name\": ",stream); json_str(stream,Json_as_string(Json_object_get(Json_array_get(maps,k),"name")));
        fprintf(stream,", \"real_nodes\": %zu, \"evidence_nodes\": %zu, \"supported_nodes\": %zu, "
                "\"missing_nodes\": %zu, \"invalid_nodes\": %zu, \"folded_nodes\": %zu, \"cell_sides\": %zu,\n",
                s->real_nodes,s->evidence_nodes,s->supported_nodes,s->missing_nodes,s->invalid_nodes,s->folded_nodes,s->cell_sides);
        jq(stream,"H_old_abs_deg",s->H_old_abs,","); jq(stream,"H_new_abs_deg",s->H_new_abs,",");
        jq(stream,"V_old_abs_deg",s->V_old_abs,","); jq(stream,"V_new_abs_deg",s->V_new_abs,",");
        fputs("  \"determinant_min\": ",stream); json_number(stream,s->determinant_min,10);
        fputs(", \"determinant_max\": ",stream); json_number(stream,s->determinant_max,10);
        fprintf(stream,"}%s\n",k+1<count?",":"");
        fprintf(stderr,"[fiber_flatten] retained map %zu: H p90 %.4f -> %.4f deg; V p90 %.4f -> %.4f deg; "
                "%zu/%zu supported, %zu missing, %zu invalid, %zu folded\n",k,s->H_old_abs.p90,s->H_new_abs.p90,
                s->V_old_abs.p90,s->V_new_abs.p90,s->supported_nodes,s->evidence_nodes,s->missing_nodes,s->invalid_nodes,s->folded_nodes);
    }
    fputs(" ]\n}\n",stream);
    int bad=ferror(stream);
    if (fclose(stream)) bad=1;
    return bad?-1:0;
}

static int audit_branches(Arena_T arena, const FiberFrame *ff, const float *X,
                           const uint8_t *valid, int H, int W, const char *out)
{
    FiberFrameBranchStats st={0};char path[1400];
    if(FiberFrame_audit_branches(arena,ff,X,valid,H,W,&st))return -1;
    snprintf(path,sizeof path,"%s/fiber_branch_audit.json",out);
    FILE *stream=fopen(path,"w");if(!stream)return -1;
    fputs("{\n \"schema\": \"fiber-transported-branch-audit-v1\",\n"
          " \"scope\": \"Supported canvas-grid paths; shortest normal transport; not original-mesh adjacency, same-leaf identity or a repair certificate\",\n",stream);
    fprintf(stream," \"evidence_nodes\": %zu, \"neighbor_pairs\": %zu, \"unsupported_paths\": %zu, \"invalid_frames\": %zu,\n"
            " \"antipodal_normals\": %zu, \"compared_pairs\": %zu, \"nearest_quarter\": [%zu,%zu,%zu,%zu], \"switches\": %zu,\n",
            st.evidence_nodes,st.neighbor_pairs,st.unsupported_paths,st.invalid_frames,st.antipodal_normals,st.compared_pairs,
            st.nearest_quarter[0],st.nearest_quarter[1],st.nearest_quarter[2],st.nearest_quarter[3],st.switches);
    jq(stream,"selected_abs_deg",st.selected_deg,",");jq(stream,"nearest_arm_residual_abs_deg",st.residual_deg,",");
    jq(stream,"switch_normal_rotation_deg",st.switch_normal_deg,",");jq(stream,"switch_coherence",st.switch_coherence,",");
    jq(stream,"switch_maximum_world_step",st.switch_maximum_world_step,",");
    fputs(" \"highest_coherence_switch_examples\": [\n",stream);
    for(size_t k=0;k<st.examples;k++){
        const FiberFrameBranchExample *e=st.example+k;
        fprintf(stream,"  {\"a\": [%d,%d], \"b\": [%d,%d], \"quarter\": %d, \"selected_deg\": %.10g, "
                "\"residual_deg\": %.10g, \"normal_deg\": %.10g, \"coherence\": %.10g, \"maximum_world_step\": %.10g}%s\n",
                e->col_a,e->row_a,e->col_b,e->row_b,e->quarter,e->selected_deg,e->residual_deg,e->normal_deg,e->coherence,
                e->maximum_world_step,k+1<st.examples?",":"");
    }
    fputs(" ]\n}\n",stream);int bad=ferror(stream);if(fclose(stream))bad=1;
    fprintf(stderr,"[fiber_flatten] transported branches: %zu switches / %zu compared; %zu unsupported paths, %zu invalid frames, %zu opposing normals\n",
            st.switches,st.compared_pairs,st.unsupported_paths,st.invalid_frames,st.antipodal_normals);
    return bad?-1:0;
}

static int audit_map_reader_selftest(const char *dir)
{
    char report[1400], binary[1400], spec_path[1400], path[1400];
    const double u[4]={1,2,NAN,4}, v[4]={5,6,7,8};
    const double uv[8]={1,5,2,6,3,7,4,8};
    Arena_T arena=Arena_new(); FiberMapResult map={0}; const char *error=NULL;
    int fails=0;
    snprintf(report,sizeof report,"%s/fiber_flatten.json",dir);
    snprintf(binary,sizeof binary,"%s/audit_map.f64",dir);
    snprintf(spec_path,sizeof spec_path,"%s/audit_map_spec.json",dir);
    FILE *f=fopen(report,"w");
    if (!f) { Arena_dispose(&arena); return 1; }
    fputs("{\"grid\":{\"rows\":2,\"cols\":2,\"cell\":4,\"row0\":12,\"col0\":24},"
          "\"guarded_map\":{\"rows\":2,\"cols\":2,\"cell\":4,\"row0\":12,\"col0\":24,"
          "\"storage\":\"interleaved_float64_uv\"}}",f); fclose(f);
    snprintf(path,sizeof path,"%s/fiber_map_u.tif",dir); fails+=write_map(path,u,2,2)!=0;
    snprintf(path,sizeof path,"%s/fiber_map_v.tif",dir); fails+=write_map(path,v,2,2)!=0;
    f=fopen(spec_path,"w");
    if (!f) { Arena_dispose(&arena); return 1; }
    fputs("{\"frame\":",f); json_str(f,dir); fputs("}",f); fclose(f);
    const JsonValue *spec=Json_parse_file(arena,spec_path,&error);
    int rc=load_audit_map(arena,spec,&map);
    fails+=rc || map.rows!=2 || map.cols!=2 || fabs(map.col0-24)>1e-12 ||
           (!rc && (!isnan(map.uv[4]) || fabs(map.uv[7]-8)>1e-12));
    f=fopen(spec_path,"w");
    if (!f) { Arena_dispose(&arena); return 1; }
    fputs("{\"map_file\":",f); json_str(f,binary); fputs(",\"report\":",f); json_str(f,report); fputs("}",f); fclose(f);
    spec=Json_parse_file(arena,spec_path,&error);
    for (size_t n=7;n<=9;n++) {
        f=fopen(binary,"wb");
        if (!f) { fails++; continue; }
        fwrite(uv,sizeof *uv,n<8?n:8,f);
        if (n==9) fputc(0,f);
        fclose(f);
        rc=load_audit_map(arena,spec,&map);
        if (n==8) fails+=rc || (!rc && memcmp(map.uv,uv,sizeof uv)!=0);
        else fails+=rc==0;
    }
    Arena_dispose(&arena);
    fprintf(stderr,"[selftest] fiber retained map reader: %s (%d failures)\n",fails?"FAIL":"PASS",fails);
    return fails;
}

static int selftest(void)
{
    int fails = 0;
    char dir[512], sub[700];
    fails += ZarrU8_selftest();
    fails += FiberFrame_selftest();
    fails += Sparse_selftest() != 0;
    /* the tifxyz writer round trip (3 x 4 nodes, one invalid) */
    if (ves_selftest_dir(dir, sizeof dir, "fiber_flatten") == 0) {
        fails += audit_map_reader_selftest(dir);
        float X[12 * 3], *Y = NULL;
        uint8_t v[12], *vb = NULL;
        int H = 0, W = 0, bad = 0;
        for (int i = 0; i < 12; i++) {
            X[i * 3] = (float)(10 + i); X[i * 3 + 1] = (float)(20 + i); X[i * 3 + 2] = (float)(30 + i);
            v[i] = (uint8_t)(i != 5);
        }
        snprintf(sub, sizeof sub, "%s/seg", dir);
        if (write_tifxyz(sub, "seg", X, v, 3, 4) != 0 || load_tifxyz(sub, &Y, &vb, &H, &W) != 0 || H != 3 || W != 4) {
            fprintf(stderr, "[selftest] fiber_flatten FAIL: tifxyz round trip\n");
            fails++;
        } else {
            for (int i = 0; i < 12; i++) {
                if (vb[i] != v[i]) bad++;
                if (v[i] && (Y[i * 3] != X[i * 3] || Y[i * 3 + 2] != X[i * 3 + 2])) bad++;
            }
            if (bad) { fprintf(stderr, "[selftest] fiber_flatten FAIL: tifxyz values\n"); fails++; }
        }
        free(Y);
        free(vb);
        /* a Windows-style name with a backslash and a quote must come back intact from the JSON reader */
        snprintf(sub, sizeof sub, "%s/seg2", dir);
        if (write_tifxyz(sub, "C:\\a \"b\"\\c", X, v, 3, 4) == 0) {
            Arena_T ja = Arena_new();
            const char *jerr = NULL;
            char mp[800];
            const JsonValue *meta = NULL;
            snprintf(mp, sizeof mp, "%s/meta.json", sub);
            meta = Json_parse_file(ja, mp, &jerr);
            if (meta == NULL || Json_as_string(Json_object_get(meta, "uuid")) == NULL ||
                strcmp(Json_as_string(Json_object_get(meta, "uuid")), "C:\\a \"b\"\\c") != 0) {
                fprintf(stderr, "[selftest] fiber_flatten FAIL: JSON string escaping\n");
                fails++;
            }
            Arena_dispose(&ja);
        } else {
            fprintf(stderr, "[selftest] fiber_flatten FAIL: tifxyz write (escaping)\n");
            fails++;
        }
        if (fails == 0) ves_selftest_clean(dir);
    }
    fprintf(stderr, "[selftest] fiber_flatten %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails ? 1 : 0;
}

int main(int argc, char **argv)
{
    Arena_T arena = NULL;
    const JsonValue *job = NULL, *rect = NULL;
    const char *err = NULL, *seg = NULL, *fibers = NULL, *out = NULL, *uuid_in = NULL;
    char uuid[256], path[1400];
    float *X = NULL, *orc = NULL, *Xn = NULL;
    uint8_t *valid = NULL, *vn = NULL;
    int H = 0, W = 0, x0 = 0, y0 = 0, x1 = 0, y1 = 0, Hn = 0, Wn = 0, threads = threads_env();
    double u_min = 0.0, v_min = 0.0, t0 = now_s(), t_field = 0.0, t_solve = 0.0, t_audit = 0.0, t_inv = 0.0, t_write = 0.0;
    size_t n_valid = 0, old_valid = 0;
    FiberVolume_T fv = NULL;
    FiberFrame ff;
    FiberFrameStats st;
    FILE *f = NULL;
    int branch_only = argc == 3 && strcmp(argv[1],"--audit-branches") == 0;
    int audit_only = branch_only || (argc == 3 && strcmp(argv[1],"--audit-map") == 0);
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc != 2 && !audit_only) {
        fprintf(stderr, "usage: fiber_flatten <job.json>\n       fiber_flatten --audit-map <job.json>\n       fiber_flatten --audit-branches <job.json>\n       fiber_flatten --selftest\n");
        return 2;
    }
    arena = Arena_new();
    const char *job_path=argv[audit_only?2:1];
    job = Json_parse_file(arena, job_path, &err);
    if (job == NULL) { fprintf(stderr, "[fiber_flatten] %s: %s\n", job_path, err ? err : "parse error"); return 1; }
    seg = Json_as_string(Json_object_get(job, "seg"));
    fibers = Json_as_string(Json_object_get(job, "fibers"));
    out = Json_as_string(Json_object_get(job, "out"));
    uuid_in = Json_as_string(Json_object_get(job, "uuid"));
    if (seg == NULL || fibers == NULL || out == NULL) {
        fprintf(stderr, "[fiber_flatten] the job needs \"seg\", \"fibers\" and \"out\"\n");
        return 1;
    }
    if (uuid_in != NULL) snprintf(uuid, sizeof uuid, "%s", uuid_in);
    else {
        const char *b = strrchr(seg, '/'), *b2 = strrchr(seg, '\\');
        if (b2 != NULL && (b == NULL || b2 > b)) b = b2;
        snprintf(uuid, sizeof uuid, "%s_fiberaligned", b != NULL ? b + 1 : seg);
    }
    stage("load");
    if (load_tifxyz(seg, &X, &valid, &H, &W) != 0) return 1;
    x0 = 0; y0 = 0; x1 = W; y1 = H;
    rect = Json_object_get(job, "rect");
    if (rect != NULL && Json_array_len(rect) == 4) {
        x0 = (int)Json_as_long(Json_array_get(rect, 0), 0);
        y0 = (int)Json_as_long(Json_array_get(rect, 1), 0);
        x1 = (int)Json_as_long(Json_array_get(rect, 2), W);
        y1 = (int)Json_as_long(Json_array_get(rect, 3), H);
    }
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 > W) x1 = W;
    if (y1 > H) y1 = H;
    for (int r = y0; r < y1; r++)
        for (int c = x0; c < x1; c++) old_valid += valid[(size_t)r * (size_t)W + (size_t)c];
    fprintf(stderr, "[fiber_flatten] %s: canvas %d x %d, window [%d, %d) x [%d, %d), %zu valid nodes; %d threads\n",
            seg, H, W, x0, x1, y0, y1, old_valid, threads);
    if (FiberVolume_open(fibers, &fv) != 0) { fprintf(stderr, "[fiber_flatten] cannot open the fibre field %s\n", fibers); return 1; }
    ves_mkdir(out);
    stage("field");
    if (FiberFrame_field(arena, X, valid, H, W, x0, y0, x1, y1, fv, &ff) != 0) { fprintf(stderr, "[fiber_flatten] field failed\n"); return 1; }
    FiberVolume_close(&fv);
    t_field = now_s();
    fprintf(stderr, "[fiber_flatten] field: %d x %d nodes, %zu samples, %zu with fibres (%.0f s)\n", ff.gh, ff.gw, ff.samples,
            ff.samples_with_fibres, t_field - t0);
    if(branch_only){
        stage("audit_transported_branches");
        int rc=audit_branches(arena,&ff,X,valid,H,W,out);
        free(X);free(valid);Arena_dispose(&arena);
        if(!rc)stage("done");return rc?1:0;
    }
    snprintf(path,sizeof path,"%s/fiber_directions.csv",out);
    if (FiberFrame_write_directions(arena,&ff,path)) {
        fprintf(stderr,"[fiber_flatten] cannot save measured fibre directions\n"); return 1;
    }
    if (audit_only) {
        stage("audit_retained_maps");
        int rc=audit_retained_maps(arena,job,&ff,out,seg,fibers);
        free(X); free(valid); Arena_dispose(&arena);
        if (!rc) stage("done");
        return rc?1:0;
    }
    stage("solve");
    if (FiberFrame_solve(arena, &ff, &st) != 0) { fprintf(stderr, "[fiber_flatten] solve failed\n"); return 1; }
    t_solve = now_s();
    fprintf(stderr, "[fiber_flatten] solve: %zu real + %zu gap nodes, %zu edges, rms residual %.3f / %.3f px per edge (%.0f s)\n",
            st.nodes_real, st.nodes_gap, st.edges, st.rms_residual_u, st.rms_residual_v, t_solve - t_field);
    stage("audit");
    snprintf(path, sizeof path, "%s/fiber_profile.csv", out);
    if (FiberFrame_audit(arena, &ff, &st, path) != 0) { fprintf(stderr, "[fiber_flatten] audit failed\n"); return 1; }
    t_audit = now_s();
    fprintf(stderr, "[fiber_flatten] fibre climb %.0f -> %.0f px over %.0f px; H |p50| %.2f -> %.2f deg; turn %.2f deg; "
                    "non-rigid p50 %.1f px; patches turn p50 %.1f deg; area p1..p99 %.3f..%.3f; %zu folded\n",
            st.climb_old_range, st.climb_new_range, st.climb_old_width, st.H_old_abs.p50, st.H_new_abs.p50, st.rotation_deg,
            st.nonrigid_px.p50, st.tile_rotation_deg.p50, st.area_p1, st.area_p99, st.flipped);
    stage("inverse");
    if (FiberFrame_inverse(arena, &ff, &st, threads, &Hn, &Wn, &u_min, &v_min, &orc) != 0) {
        fprintf(stderr, "[fiber_flatten] inverse failed\n");
        return 1;
    }
    if (FiberFrame_resample(arena, X, valid, H, W, x0, y0, x1, y1, orc, Hn, Wn, &Xn, &vn, &n_valid) != 0) {
        fprintf(stderr, "[fiber_flatten] resample failed\n");
        return 1;
    }
    t_inv = now_s();
    fprintf(stderr, "[fiber_flatten] new page %d x %d, %zu valid px (%.4f of the old window's)\n", Hn, Wn, n_valid,
            old_valid ? (double)n_valid / (double)old_valid : 0.0);
    stage("write");
    snprintf(path, sizeof path, "%s/%s", out, uuid);
    if (write_tifxyz(path, uuid, Xn, vn, Hn, Wn) != 0) { fprintf(stderr, "[fiber_flatten] cannot write %s\n", path); return 1; }
    snprintf(path, sizeof path, "%s/fiber_map_u.tif", out);
    write_map(path, ff.U, ff.gh, ff.gw);
    snprintf(path, sizeof path, "%s/fiber_map_v.tif", out);
    write_map(path, ff.V, ff.gh, ff.gw);
    t_write = now_s();
    snprintf(path, sizeof path, "%s/fiber_flatten.json", out);
    f = fopen(path, "w");
    if (f == NULL) return 1;
    fprintf(f, "{\n  \"schema\": \"fiber-flatten-v1\",\n  \"seg\": ");
    json_str(f, seg);
    fprintf(f, ",\n  \"fibers\": ");
    json_str(f, fibers);
    fprintf(f, ",\n  \"out\": ");
    json_str(f, out);
    fprintf(f, ",\n  \"uuid\": ");
    json_str(f, uuid);
    fprintf(f, ",\n");
    fprintf(f, "  \"rect\": [%d, %d, %d, %d],\n  \"canvas\": [%d, %d],\n  \"threads\": %d,\n", x0, y0, x1, y1, H, W, threads);
    fprintf(f, "  \"constants\": {\"sample_px\": %d, \"half_depth\": %d, \"cell_px\": %d, \"sigma_px\": %.1f, \"evidence_rel\": %.3f, "
               "\"eps\": %g, \"gap_cells\": %d, \"gap_weight\": %.3f, \"inv_iters\": %d, \"inv_step\": %g, \"inv_tol\": %g, "
               "\"pad_cells\": %d, \"tile_px\": %d, \"profile_min\": %d},\n",
            FIBER_FLATTEN_SAMPLE_PX, FIBER_FLATTEN_HALF_DEPTH, FIBER_FLATTEN_CELL_PX, FIBER_FLATTEN_SIGMA_PX,
            FIBER_FLATTEN_EVIDENCE_REL, FIBER_FLATTEN_EPS, FIBER_FLATTEN_GAP_CELLS, FIBER_FLATTEN_GAP_WEIGHT,
            FIBER_FLATTEN_INV_ITERS, FIBER_FLATTEN_INV_STEP, FIBER_FLATTEN_INV_TOL, FIBER_FLATTEN_PAD_CELLS,
            FIBER_FLATTEN_TILE_PX, FIBER_FLATTEN_PROFILE_MIN);
    fprintf(f, "  \"grid\": {\"rows\": %d, \"cols\": %d, \"cell\": %d, \"row0\": %d, \"col0\": %d},\n", ff.gh, ff.gw, ff.cell, ff.gy0, ff.gx0);
    fprintf(f, "  \"field\": {\"samples\": %zu, \"with_fibres\": %zu, \"evidence_share\": %.4f},\n", ff.samples,
            ff.samples_with_fibres, st.evidence_share);
    fprintf(f, "  \"solve\": {\"real_nodes\": %zu, \"gap_nodes\": %zu, \"edges\": %zu, \"rms_residual_u\": %.4f, "
               "\"rms_residual_v\": %.4f, \"L_h_p50\": %.4f, \"L_v_p50\": %.4f},\n",
            st.nodes_real, st.nodes_gap, st.edges, st.rms_residual_u, st.rms_residual_v, st.L_h_p50, st.L_v_p50);
    fputs("  \"H_from_rows_old_signed_mean_deg\": ", f); json_number(f, st.H_old_signed_mean, 3);
    fputs(",\n  \"H_from_rows_new_signed_mean_deg\": ", f); json_number(f, st.H_new_signed_mean, 3); fputs(",\n", f);
    jq(f, "H_from_rows_old_abs_deg", st.H_old_abs, ",");
    jq(f, "H_from_rows_new_abs_deg", st.H_new_abs, ",");
    fputs("  \"fibre_climb\": {\"old\": ", f);
    jclimb(f, st.climb_old_range, st.climb_old_bow, st.climb_old_end, st.climb_old_width, st.col_angle_old_p10, st.col_angle_old_p90);
    fputs(", \"new\": ", f);
    jclimb(f, st.climb_new_range, st.climb_new_bow, st.climb_new_end, st.climb_new_width, st.col_angle_new_p10, st.col_angle_new_p90);
    fputs("},\n  \"map\": {\"rotation_deg\": ", f); json_number(f, st.rotation_deg, 3);
    fputs(", \"nonrigid_px\": {", f); jq_values(f, st.nonrigid_px, 2);
    fputs(", \"max\": ", f); json_number(f, st.nonrigid_max, 2);
    fprintf(f, "}, \"patches\": %zu, \"patch_nonrigid_px\": {", st.tiles); jq_values(f, st.tile_nonrigid_px, 2);
    fputs("}, \"patch_rotation_deg\": {", f); jq_values(f, st.tile_rotation_deg, 2);
    fprintf(f, "}, \"folded_nodes\": %zu, \"area_ratio\": {\"p1\": %.4f, \"p10\": %.4f, \"p50\": %.4f, \"p90\": %.4f, \"p99\": %.4f}, "
               "\"affine\": [[%.6f, %.6f, %.3f], [%.6f, %.6f, %.3f]]},\n",
            st.flipped, st.area_p1, st.area_p10, st.area_p50,
            st.area_p90, st.area_p99, st.A[0], st.A[1], st.b[0], st.A[2], st.A[3], st.b[1]);
    fprintf(f, "  \"new_page\": {\"rows\": %d, \"cols\": %d, \"u_min\": %.1f, \"v_min\": %.1f, \"valid_px\": %zu, "
               "\"old_window_valid_px\": %zu, \"kept\": %.4f},\n",
            Hn, Wn, u_min, v_min, n_valid, old_valid, old_valid ? (double)n_valid / (double)old_valid : 0.0);
    fprintf(f, "  \"seconds\": {\"field\": %.1f, \"solve\": %.1f, \"audit\": %.1f, \"inverse\": %.1f, \"write\": %.1f}\n}\n",
            t_field - t0, t_solve - t_field, t_audit - t_solve, t_inv - t_audit, t_write - t_inv);
    fclose(f);
    fprintf(stderr, "[fiber_flatten] wrote %s/%s and fiber_flatten.json\n", out, uuid);
    stage("done");
    free(X);
    free(valid);
    Arena_dispose(&arena);
    return 0;
}
