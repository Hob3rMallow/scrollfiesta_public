/* Close-up windows are presentation crops, never resampling transforms. A median
 * translation alone can cut off a small piece beside a much larger one. Retain
 * every supported required-piece sample, with a sampling apron, and give the
 * expanded rectangle free space below stationary neighbours when necessary. */
#include "closeup_window.h"
#include "arena.h"
#include "json_read.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

typedef struct CwPlan {
    int window[4], page[4], added[4], shift_y, expanded;
    size_t samples;
} CwPlan;

static int cw_int(const JsonValue *v, int *out)
{
    double x = Json_as_double(v, NAN);
    if (!isfinite(x) || floor(x) != x || fabs(x) > INT_MAX / 4) return -1;
    *out = (int)x;
    return 0;
}

static int cw_array(const JsonValue *v, int *out, size_t n)
{
    if (Json_type(v) != JSON_ARRAY || Json_array_len(v) != n) return -1;
    for (size_t i = 0; i < n; ++i)
        if (cw_int(Json_array_get(v, i), &out[i])) return -1;
    return 0;
}

static int cw_rect(const int r[4]) { return r[0] < r[2] && r[1] < r[3]; }
static int cw_min(int a, int b) { return a < b ? a : b; }
static int cw_max(int a, int b) { return a > b ? a : b; }
static int cw_meet(const int a[4], const int b[4])
{
    return a[0] < b[2] && b[0] < a[2] && a[1] < b[3] && b[1] < a[3];
}

static int cw_plan(Arena_T arena, const JsonValue *j, CwPlan *out)
{
    int old[4], page[4], origin[2], canvas[2], pad, outside = 0;
    const JsonValue *required = Json_object_get(j, "required");
    const JsonValue *matches = Json_object_get(j, "matches");
    const JsonValue *obstacles = Json_object_get(j, "obstacles");
    size_t nr = Json_array_len(required), nm = Json_array_len(matches), no = Json_array_len(obstacles);
    if (cw_array(Json_object_get(j, "window"), old, 4) || !cw_rect(old) ||
        cw_array(Json_object_get(j, "page_rect"), page, 4) || !cw_rect(page) ||
        cw_array(Json_object_get(j, "page_origin"), origin, 2) ||
        cw_array(Json_object_get(j, "canvas"), canvas, 2) || canvas[0] <= 0 || canvas[1] <= 0 ||
        cw_int(Json_object_get(j, "pad"), &pad) || pad < 0 || pad > 4096 ||
        page[2] - page[0] != old[2] - old[0] || page[3] - page[1] != old[3] - old[1] ||
        Json_type(required) != JSON_ARRAY || Json_type(matches) != JSON_ARRAY || Json_type(obstacles) != JSON_ARRAY)
        return -1;
    memset(out, 0, sizeof *out);
    memcpy(out->window, old, sizeof old);
    memcpy(out->page, page, sizeof page);
    int *pieces = ARENA_CALLOC(arena, nr ? nr : 1, sizeof *pieces);
    size_t *hits = ARENA_CALLOC(arena, nr ? nr : 1, sizeof *hits);
    int (*blocks)[4] = ARENA_CALLOC(arena, no ? no : 1, sizeof *blocks);
    for (size_t i = 0; i < nr; ++i) {
        if (cw_int(Json_array_get(required, i), &pieces[i]) || pieces[i] < 0) return -1;
        for (size_t k = 0; k < i; ++k) if (pieces[k] == pieces[i]) return -1;
    }
    for (size_t i = 0; i < no; ++i)
        if (cw_array(Json_array_get(obstacles, i), blocks[i], 4) || !cw_rect(blocks[i])) return -1;
    for (size_t i = 0; i < nm; ++i) {
        int m[3];
        if (cw_array(Json_array_get(matches, i), m, 3) || m[1] < 0 || m[1] >= canvas[0] ||
            m[2] < 0 || m[2] >= canvas[1]) return -1;
        for (size_t k = 0; k < nr; ++k) if (m[0] == pieces[k]) {
            hits[k]++;
            out->samples++;
            outside |= m[1] < old[0] || m[1] >= old[2] || m[2] < old[1] || m[2] >= old[3];
            out->window[0] = cw_min(out->window[0], cw_max(0, m[1] - pad));
            out->window[1] = cw_min(out->window[1], cw_max(0, m[2] - pad));
            out->window[2] = cw_max(out->window[2], cw_min(canvas[0], m[1] + pad + 1));
            out->window[3] = cw_max(out->window[3], cw_min(canvas[1], m[2] + pad + 1));
            break;
        }
    }
    for (size_t i = 0; i < nr; ++i) if (!hits[i]) return -1;
    /* New material has no correspondence in the old piece. Its measured cell
     * bounds are a separate obligation, never fabricated required-piece hits. */
    const JsonValue *material = Json_object_get(j, "material_bounds");
    if (material) {
        int bounds[4];
        if (cw_array(material, bounds, 4) || !cw_rect(bounds) || bounds[0] < 0 || bounds[1] < 0 ||
            bounds[2] > canvas[0] || bounds[3] > canvas[1]) return -1;
        outside |= bounds[0] < old[0] || bounds[1] < old[1] || bounds[2] > old[2] || bounds[3] > old[3];
        out->window[0] = cw_min(out->window[0], bounds[0]);
        out->window[1] = cw_min(out->window[1], bounds[1]);
        out->window[2] = cw_max(out->window[2], bounds[2]);
        out->window[3] = cw_max(out->window[3], bounds[3]);
    }
    if (!outside) {
        memcpy(out->window, old, sizeof old);
        return 0; /* an apron alone is not a reason to move a fitting crop */
    }
    out->added[0] = old[0] - out->window[0];
    out->added[1] = old[1] - out->window[1];
    out->added[2] = out->window[2] - old[2];
    out->added[3] = out->window[3] - old[3];
    out->expanded = memcmp(old, out->window, sizeof old) != 0;
    if (!out->expanded) return 0; /* preserve existing placements that already fit */
    int width = out->window[2] - out->window[0], height = out->window[3] - out->window[1];
    out->page[0] = cw_max(origin[0], page[0] - out->added[0]);
    out->page[1] = cw_max(origin[1], page[1] - out->added[1]);
    out->page[2] = out->page[0] + width;
    out->page[3] = out->page[1] + height;
    int start_y = out->page[1];
    /* Each move passes at least one finite obstacle bottom. Re-scan because a
     * move may meet a lower fragment. Touching edges do not overlap. */
    for (size_t pass = 0; pass <= no; ++pass) {
        int y = out->page[1];
        for (size_t i = 0; i < no; ++i)
            if (cw_meet(out->page, blocks[i])) y = cw_max(y, blocks[i][3]);
        if (y == out->page[1]) {
            out->shift_y = y - start_y;
            return 0;
        }
        out->page[1] = y;
        out->page[3] = y + height;
    }
    return -1;
}

int CloseupWindow_run(const char *input)
{
    Arena_T arena = Arena_new();
    const char *error = NULL;
    const JsonValue *j = Json_parse_file(arena, input, &error);
    CwPlan p = {0};
    int rc = j == NULL ? -1 : cw_plan(arena, j, &p);
    if (rc) fprintf(stderr, "closeup window: invalid plan input or a required piece has no supported samples%s%s\n",
                    error ? ": " : "", error ? error : "");
    else if (printf("{\"schema\":\"closeup-window-v1\",\"window\":[%d,%d,%d,%d],"
                    "\"page_rect\":[%d,%d,%d,%d],\"added_px\":[%d,%d,%d,%d],"
                    "\"expanded\":%s,\"neighbor_shift_y\":%d,\"required_samples\":%zu}\n",
                    p.window[0],p.window[1],p.window[2],p.window[3],p.page[0],p.page[1],p.page[2],p.page[3],
                    p.added[0],p.added[1],p.added[2],p.added[3],p.expanded ? "true" : "false",p.shift_y,p.samples) < 0)
        rc = -1;
    Arena_dispose(&arena);
    return rc ? 1 : 0;
}

static int cw_case(const char *json, const int *window, const int *page, int bad)
{
    Arena_T arena = Arena_new();
    const JsonValue *j = Json_parse(arena, json, NULL);
    CwPlan p = {0};
    int rc = cw_plan(arena, j, &p);
    int fail = bad ? rc == 0 : rc != 0 || memcmp(p.window, window, sizeof p.window) || memcmp(p.page, page, sizeof p.page);
    Arena_dispose(&arena);
    return fail != 0;
}

int CloseupWindow_selftest(void)
{
    int fails = 0;
    /* All old-piece matches fit; a newly grown lobe still needs room. A
     * stationary neighbour forces relocation without clipping or resampling. */
    fails += cw_case("{\"window\":[20,20,80,80],\"page_rect\":[20,20,80,80],\"page_origin\":[0,0],"
        "\"canvas\":[120,100],\"pad\":0,\"required\":[1],\"matches\":[[1,40,40]],"
        "\"material_bounds\":[10,10,120,100],\"obstacles\":[[80,0,140,90]]}",
        (int[4]){10,10,120,100}, (int[4]){10,90,120,180}, 0);
    /* A partial selection can grow without removing any entire piece. */
    fails += cw_case("{\"window\":[20,20,80,80],\"page_rect\":[20,20,80,80],\"page_origin\":[0,0],"
        "\"canvas\":[120,100],\"pad\":0,\"required\":[],\"matches\":[],"
        "\"material_bounds\":[10,10,120,100],\"obstacles\":[]}",
        (int[4]){10,10,120,100}, (int[4]){10,10,120,100}, 0);
    fails += cw_case("{\"window\":[0,0,10,10],\"page_rect\":[0,0,10,10],\"page_origin\":[0,0],"
        "\"canvas\":[10,10],\"pad\":0,\"required\":[1],\"matches\":[],"
        "\"material_bounds\":[0,0,10,10],\"obstacles\":[]}", NULL,NULL,1);
    fails += cw_case("{\"window\":[0,0,10,10],\"page_rect\":[0,0,10,10],\"page_origin\":[0,0],"
        "\"canvas\":[10,10],\"pad\":0,\"required\":[],\"matches\":[],"
        "\"material_bounds\":[0,0,11,10],\"obstacles\":[]}", NULL,NULL,1);
    /* A small displaced piece must survive even though the large piece sets the
     * median. Neither the neighbour nor the pixels' scale may change. */
    fails += cw_case("{\"window\":[100,0,900,200],\"page_rect\":[1000,0,1800,200],\"page_origin\":[0,0],"
        "\"canvas\":[1000,400],\"pad\":16,\"required\":[1,2],\"matches\":[[1,50,40],[2,850,150]],"
        "\"obstacles\":[[800,0,1000,80],[900,100,1000,140]]}",
        (int[4]){34,0,900,200}, (int[4]){934,140,1800,340}, 0);
    /* A distant incidental piece cannot inflate the crop. An already fitting
     * crop keeps its established placement, including intentional box overlap. */
    fails += cw_case("{\"window\":[100,0,900,200],\"page_rect\":[1000,0,1800,200],\"page_origin\":[0,0],"
        "\"canvas\":[1200,400],\"pad\":16,\"required\":[1],\"matches\":[[1,200,40],[9,1100,350]],"
        "\"obstacles\":[[900,0,1100,100]]}",
        (int[4]){100,0,900,200}, (int[4]){1000,0,1800,200}, 0);
    /* Growth on every side; clamp the sampling apron to actual canvas cells. */
    fails += cw_case("{\"window\":[20,20,80,80],\"page_rect\":[0,0,60,60],\"page_origin\":[0,0],"
        "\"canvas\":[100,100],\"pad\":16,\"required\":[1],\"matches\":[[1,0,0],[1,99,99]],\"obstacles\":[]}",
        (int[4]){0,0,100,100}, (int[4]){0,0,100,100}, 0);
    fails += cw_case("{\"window\":[20,20,80,80],\"page_rect\":[0,0,60,60],\"page_origin\":[0,0],"
        "\"canvas\":[100,100],\"pad\":16,\"required\":[1],\"matches\":[[1,20,20],[1,79,79]],\"obstacles\":[]}",
        (int[4]){20,20,80,80}, (int[4]){0,0,60,60}, 0);
    fails += cw_case("{\"window\":[0,0,10,10],\"page_rect\":[0,0,10,10],\"page_origin\":[0,0],"
        "\"canvas\":[10,10],\"pad\":1,\"required\":[1,2],\"matches\":[[1,2,2]],\"obstacles\":[]}", NULL,NULL,1);
    fails += cw_case("{\"window\":[0,0,10,10],\"page_rect\":[0,0,10,10],\"page_origin\":[0,0],"
        "\"canvas\":[10,10],\"pad\":1,\"required\":[1],\"matches\":[[1,10,2]],\"obstacles\":[]}", NULL,NULL,1);
    fails += cw_case("{\"window\":[0,0,10,10],\"page_rect\":[0,0,9,10],\"page_origin\":[0,0],"
        "\"canvas\":[10,10],\"pad\":1,\"required\":[],\"matches\":[],\"obstacles\":[]}", NULL,NULL,1);
    fprintf(stderr, "[selftest] closeup window: %s (%d failures)\n", fails ? "FAIL" : "PASS", fails);
    return fails != 0;
}
