/*
 * welded_winding_lift.c -- lift the shipped per-chart winding atlas through an
 * already-certified grid_weld mesh.
 *
 * This deliberately does not re-parameterize a welded disk. Surviving source
 * vertices retain the UV written by scroll_whole, seam-refinement vertices
 * inherit exact midpoint UV through grid_weld lineage, and only the small
 * parentless residue is harmonically interpolated.  Chart gauges are then
 * translated from either the registered winding field or one topology-lifted
 * winding coordinate. Final sheet components are moved rigidly into one
 * horizontal strip.
 */

#include "../common/arena.h"
#include "../common/mesh_bin.h"
#include "../common/obj_io.h"
#include "../common/pipeline_constants.h"
#include "../common/ves_platform.h"
#include "../flatten/ribbon.h"
#include "../holefill/hole_fill.h"
#include "../unroll/piece_set.h"
#include "../unroll/scaffold.h"
#include "../whole/cube_register.h"

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WWL_PI 3.1415926535897932384626433832795
#define WWL_2PI (2.0 * WWL_PI)

typedef struct {
    size_t *offset;
    int32_t *neighbor;
    size_t nv;
} WwlGraph;

typedef struct {
    int32_t *parent;
    uint8_t *rank;
    size_t n;
} WwlDsu;

typedef struct {
    int32_t chart;
    double source_u, source_v;
    double target_u, target_v;
} WwlChartObservation;

typedef struct {
    size_t observations;
    size_t reflected_charts;
    size_t point_gauge_charts;
    double residual_median;
    double residual_maximum;
} WwlChartGaugeStats;

typedef struct {
    size_t charts;
    size_t anchors;
    int32_t minimum_turn;
    int32_t maximum_turn;
    double residual_median_turns;
    double residual_maximum_turns;
} WwlTopologyGaugeStats;

typedef struct {
    int32_t component;
    size_t faces;
    double umin, umax, vmin, vmax;
} WwlPackRecord;

typedef struct {
    int32_t lo, hi;
    int32_t base, moving;
    double source_u, source_v;
    double target_u, target_v;
    double maximum_cross_edge;
    double maximum_physical_winding_delta;
    double maximum_raw_phi_delta_turns;
    double maximum_registered_phi_delta_turns;
} WwlBridgeObservation;

typedef struct {
    int32_t lo, hi;
    int32_t base, moving;
    double cosine, sine;
    double translate_u, translate_v;
    size_t observations;
    double median_residual;
    double maximum_residual;
    size_t pair_observations, lo_observations, hi_observations;
    double translation_u, translation_v;
    double translation_median_residual;
    double translation_maximum_residual;
    double cross_edge_median, cross_edge_maximum;
    double physical_winding_median, physical_winding_maximum;
    size_t raw_phi_observations, registered_phi_observations;
    double raw_phi_median_turns, raw_phi_maximum_turns;
    double registered_phi_median_turns, registered_phi_maximum_turns;
} WwlBridgeGaugeEdge;

typedef struct {
    const char *path;
    const size_t *chart_faces;
    const int32_t *component;
    const size_t *component_faces;
    size_t ncomponent;
    const int32_t *anchor_source;
    const PieceSet *source;
    const ScaffoldCalib *calib;
    const float *registered_phi;
} WwlBridgeAudit;

typedef struct {
    size_t bridge_faces;
    size_t same_chart_faces;
    size_t two_chart_faces;
    size_t three_chart_faces;
    size_t invalid_faces;
    size_t observations;
    size_t gauge_edges;
    size_t gauge_components;
    size_t gauge_cycles;
    size_t orientation_normalized_charts;
    size_t translation_selected_edges;
    double observation_residual_median;
    double observation_residual_maximum;
    double cycle_residual_maximum;
} WwlBridgeGaugeStats;

typedef struct {
    uint32_t xyz[3];
    int32_t source;
    uint8_t occupied;
} WwlPositionSlot;

typedef struct {
    int32_t *parent0;
    int32_t *parent1;
    char (*cube_ids)[48];
    size_t nv;
    size_t source_faces;
    size_t n_cubes;
    size_t source_anchors;
    uint32_t version;
} WwlLineage;

#define WWL_CHART_LAYOUT_V3_KNOTS 9
typedef struct {
    char cube_id[48];
    int32_t local_chart;
    int32_t representative_vertex;
    double shift_u;
    double shift_u_per_v;
    double v_reference;
    double v_step;
    double shift_v[WWL_CHART_LAYOUT_V3_KNOTS];
    int32_t winding_turn;
    int32_t relation_component;
} WwlChartLayoutNode;

typedef struct {
    uint32_t version;
    double shift_u;
    double shift_u_per_v;
    double v_reference;
    double v_step;
    double shift_v[WWL_CHART_LAYOUT_V3_KNOTS];
    int32_t winding_turn;
    int32_t relation_component;
    int32_t layout_node;
} WwlChartLayoutCorrection;

typedef struct {
    uint32_t a;
    uint32_t b;
    double support;
} WwlChartLayoutEdge;

typedef struct {
    int32_t a;
    int32_t b;
    double support;
} WwlMappedChartRelation;

static int wwl_mapped_chart_relation_compare(const void *va,const void *vb)
{
    const WwlMappedChartRelation *a=(const WwlMappedChartRelation*)va;
    const WwlMappedChartRelation *b=(const WwlMappedChartRelation*)vb;
    if(a->a!=b->a)return a->a<b->a?-1:1;
    if(a->b!=b->b)return a->b<b->b?-1:1;
    if(a->support!=b->support)return a->support>b->support?-1:1;
    return 0;
}

static int wwl_chart_layout_edge_compare(const void *va,const void *vb)
{
    const WwlChartLayoutEdge *a=(const WwlChartLayoutEdge*)va;
    const WwlChartLayoutEdge *b=(const WwlChartLayoutEdge*)vb;
    if(a->a!=b->a)return a->a<b->a?-1:1;
    if(a->b!=b->b)return a->b<b->b?-1:1;
    if(a->support!=b->support)return a->support>b->support?-1:1;
    return 0;
}

typedef struct {
    WwlChartLayoutNode *node;
    WwlChartLayoutEdge *edge;
    size_t nnode;
    size_t relation_edges;
    size_t cycle_edges;
    size_t relation_components;
    uint32_t version;
} WwlChartLayout;

static void wwl_usage(const char *program)
{
    fprintf(stderr,
        "usage:\n"
        "  %s --selftest\n"
        "  %s <placed_dir> <welded.obj> <out.obj> [--source-faces N] "
        "[--sweeps N] [--gap F] [--no-pack] [--raw-uv] [--dump-obj] "
        "[--carry-registered|"
        "--preserve-chart-gauges|--assemble-bridge-gauges|"
        "--scaffold-oriented-bridge-gauges|"
        "--registered-chart-gauges|--topology-winding-gauges|"
        "--topology-winding-scaffold] "
        "[--pin-v-to-z] [--fit-ribbon] [--atlas-v-rebase] "
        "[--lineage <grid_weld.lineage.bin>] "
        "[--chart-layout <chart_layout.gwlayout>] "
        "[--bridge-edge-audit <audit.csv>] "
        "[--placed-chart-diagnostic|--source-chart-diagnostic|"
        "--final-source-only]\n",
        program, program);
}

static void wwl_chart_layout_dispose(WwlChartLayout *layout)
{
    free(layout->node);
    free(layout->edge);
    memset(layout, 0, sizeof(*layout));
}

static int wwl_chart_layout_load(const char *path,WwlChartLayout *out)
{
    static const char expected[8]={'G','W','L','A','Y','T','1','\n'};
    char magic[8];uint32_t version=0,flags=0;
    uint64_t nnode64=0,nedge64=0,ncycle64=0;
    FILE *file=NULL;int tail;int32_t maximum_component=-1;
    if(!path||!out)return -1;
    memset(out,0,sizeof(*out));
    file=fopen(path,"rb");if(!file)return -1;
    if(fread(magic,1,sizeof(magic),file)!=sizeof(magic)||
       fread(&version,sizeof(version),1,file)!=1||
       fread(&flags,sizeof(flags),1,file)!=1||
       fread(&nnode64,sizeof(nnode64),1,file)!=1||
       fread(&nedge64,sizeof(nedge64),1,file)!=1||
       fread(&ncycle64,sizeof(ncycle64),1,file)!=1||
       memcmp(magic,expected,sizeof(magic))!=0||
       (version!=1&&version!=2&&version!=3&&version!=4)||flags!=0||
       nnode64==0||nnode64>SIZE_MAX||nedge64>SIZE_MAX||
       ncycle64>nedge64||nnode64>SIZE_MAX/sizeof(*out->node))
        goto fail;
    out->node=(WwlChartLayoutNode*)calloc(
        (size_t)nnode64,sizeof(*out->node));
    if(!out->node)goto fail;
    out->nnode=(size_t)nnode64;
    out->relation_edges=(size_t)nedge64;
    out->cycle_edges=(size_t)ncycle64;
    out->version=version;
    for(size_t i=0;i<out->nnode;i++){
        WwlChartLayoutNode *node=&out->node[i];
        if(fread(node->cube_id,1,sizeof(node->cube_id),file)!=
               sizeof(node->cube_id)||
           fread(&node->local_chart,sizeof(node->local_chart),1,file)!=1||
            fread(&node->representative_vertex,
                  sizeof(node->representative_vertex),1,file)!=1)
             goto fail;
        if(version<=2&&
           fread(&node->shift_u,sizeof(node->shift_u),1,file)!=1)
             goto fail;
        if(version==2&&
           (fread(&node->shift_u_per_v,
                  sizeof(node->shift_u_per_v),1,file)!=1||
            fread(&node->v_reference,
                  sizeof(node->v_reference),1,file)!=1))
             goto fail;
        if(version>=3&&
           (fread(&node->v_reference,
                  sizeof(node->v_reference),1,file)!=1||
            fread(&node->v_step,sizeof(node->v_step),1,file)!=1||
            fread(node->shift_v,sizeof(node->shift_v[0]),
                  WWL_CHART_LAYOUT_V3_KNOTS,file)!=
                WWL_CHART_LAYOUT_V3_KNOTS))
             goto fail;
        if(fread(&node->winding_turn,
                 sizeof(node->winding_turn),1,file)!=1||
            fread(&node->relation_component,
                  sizeof(node->relation_component),1,file)!=1||
            memchr(node->cube_id,'\0',sizeof(node->cube_id))==NULL||
            node->local_chart<0||node->representative_vertex<0||
            node->relation_component<0||
            (uint64_t)node->relation_component>=nnode64||
            !isfinite(node->shift_u)||!isfinite(node->shift_u_per_v)||
            !isfinite(node->v_reference)||
            (version>=3&&(!isfinite(node->v_step)||node->v_step<=0.0)))
             goto fail;
        if(node->relation_component>maximum_component)
            maximum_component=node->relation_component;
        if(version>=3)
            for(size_t k=0;k<WWL_CHART_LAYOUT_V3_KNOTS;k++)
                if(!isfinite(node->shift_v[k]))goto fail;
    }
    if(version>=4&&out->relation_edges!=0){
        if(out->relation_edges>SIZE_MAX/sizeof(*out->edge))goto fail;
        out->edge=(WwlChartLayoutEdge*)malloc(
            out->relation_edges*sizeof(*out->edge));
        if(!out->edge)goto fail;
        for(size_t e=0;e<out->relation_edges;e++){
            WwlChartLayoutEdge *edge=&out->edge[e];
            if(fread(&edge->a,sizeof(edge->a),1,file)!=1||
               fread(&edge->b,sizeof(edge->b),1,file)!=1||
               fread(&edge->support,sizeof(edge->support),1,file)!=1||
               edge->a>=out->nnode||edge->b>=out->nnode||
               edge->a==edge->b||!isfinite(edge->support)||
               edge->support<=0.0||
               out->node[edge->a].relation_component!=
                   out->node[edge->b].relation_component||
               out->node[edge->a].winding_turn!=
                   out->node[edge->b].winding_turn)
                goto fail;
        }
        qsort(out->edge,out->relation_edges,sizeof(*out->edge),
              wwl_chart_layout_edge_compare);
        for(size_t e=1;e<out->relation_edges;e++)
            if(out->edge[e-1].a==out->edge[e].a&&
               out->edge[e-1].b==out->edge[e].b)goto fail;
    }
    tail=fgetc(file);fclose(file);file=NULL;
    if(tail!=EOF){wwl_chart_layout_dispose(out);return -1;}
    {
        size_t ncomponent=(size_t)maximum_component+1;
        uint8_t *component_seen=(uint8_t*)calloc(ncomponent,1);
        if(!component_seen){wwl_chart_layout_dispose(out);return -1;}
        for(size_t i=0;i<out->nnode;i++)
            component_seen[out->node[i].relation_component]=1;
        for(size_t c=0;c<ncomponent;c++)if(!component_seen[c]){
            free(component_seen);wwl_chart_layout_dispose(out);return -1;
        }
        free(component_seen);
        out->relation_components=ncomponent;
    }
    return 0;
fail:
    if(file)fclose(file);
    wwl_chart_layout_dispose(out);
    return -1;
}

static void wwl_lineage_dispose(WwlLineage *lineage)
{
    free(lineage->parent0);
    free(lineage->parent1);
    free(lineage->cube_ids);
    memset(lineage, 0, sizeof(*lineage));
}

static int wwl_lineage_load(const char *path, size_t expected_nv,
                            const float *verts, WwlLineage *out,
                            size_t *out_midpoints)
{
    const char magic_v1[8] = {'G','W','L','I','N','1','\r','\n'};
    const char magic_v2[8] = {'G','W','L','I','N','2','\r','\n'};
    char magic[8];
    uint32_t version = 0, reserved = 0;
    uint64_t nv64 = 0, source_faces64 = 0, cube_count64 = 0;
    FILE *fp = NULL;
    size_t midpoint_count = 0, source_anchor_count = 0, mismatch_count = 0;
    memset(out, 0, sizeof(*out));
    fp = fopen(path, "rb");
    if (fp == NULL) return -1;
    if (fread(magic, 1, sizeof(magic), fp) != sizeof(magic) ||
        fread(&version, sizeof(version), 1, fp) != 1 ||
        fread(&reserved, sizeof(reserved), 1, fp) != 1 ||
        fread(&nv64, sizeof(nv64), 1, fp) != 1 ||
        fread(&source_faces64, sizeof(source_faces64), 1, fp) != 1 ||
        reserved != 0 || nv64 != (uint64_t)expected_nv ||
        source_faces64 == 0 || source_faces64 > (uint64_t)SIZE_MAX)
        goto fail;
    if (memcmp(magic, magic_v1, sizeof(magic)) == 0 && version == 1) {
        cube_count64 = 0;
    } else if (memcmp(magic, magic_v2, sizeof(magic)) == 0 && version == 2) {
        if (fread(&cube_count64, sizeof(cube_count64), 1, fp) != 1 ||
            cube_count64 == 0 || cube_count64 > (uint64_t)INT32_MAX ||
            cube_count64 > (uint64_t)SIZE_MAX)
            goto fail;
        out->cube_ids = (char (*)[48])malloc(
            (size_t)cube_count64 * sizeof(*out->cube_ids));
        if (out->cube_ids == NULL ||
            fread(out->cube_ids, sizeof(*out->cube_ids),
                  (size_t)cube_count64, fp) != (size_t)cube_count64)
            goto fail;
        for (size_t c = 0; c < (size_t)cube_count64; c++)
            if (memchr(out->cube_ids[c], '\0',
                       sizeof(out->cube_ids[c])) == NULL)
                goto fail;
    } else {
        goto fail;
    }
    out->parent0 = (int32_t *)malloc(
        (expected_nv ? expected_nv : 1) * sizeof(*out->parent0));
    out->parent1 = (int32_t *)malloc(
        (expected_nv ? expected_nv : 1) * sizeof(*out->parent1));
    if (out->parent0 == NULL || out->parent1 == NULL ||
        fread(out->parent0, sizeof(*out->parent0), expected_nv, fp) != expected_nv ||
        fread(out->parent1, sizeof(*out->parent1), expected_nv, fp) != expected_nv ||
        fgetc(fp) != EOF)
        goto fail;
    fclose(fp);
    fp = NULL;
    out->nv = expected_nv;
    out->source_faces = (size_t)source_faces64;
    out->n_cubes = (size_t)cube_count64;
    out->version = version;

    /* A lineage captured from a different zipper result is acceptable only
     * when its entire ordered vertex stream is still exactly the same.  Check
     * every recorded midpoint against the mesh before using a single parent. */
    for (size_t v = 0; v < expected_nv; v++) {
        int32_t p0 = out->parent0[v], p1 = out->parent1[v];
        if (version == 2 && p0 <= -2) {
            int64_t cube = -(int64_t)p0 - 2;
            if (p1 < 0 || cube < 0 || (uint64_t)cube >= cube_count64)
                mismatch_count++;
            else
                source_anchor_count++;
            continue;
        }
        if (p0 == -1 && p1 == -1) continue;
        if (p0 < 0 || p1 < 0 || (size_t)p0 >= v || (size_t)p1 >= v) {
            mismatch_count++;
            continue;
        }
        midpoint_count++;
        for (int axis = 0; axis < 3; axis++) {
            float midpoint = (float)(0.5 *
                ((double)verts[(size_t)p0 * 3 + (size_t)axis] +
                 (double)verts[(size_t)p1 * 3 + (size_t)axis]));
            if (memcmp(&midpoint,
                       &verts[v * 3 + (size_t)axis], sizeof(midpoint)) != 0) {
                mismatch_count++;
                break;
            }
        }
    }
    out->source_anchors = source_anchor_count;
    if (mismatch_count != 0) {
        fprintf(stderr,
            "welded_winding_lift: lineage rejected: %zu ancestry/geometry "
            "mismatch(es)\n", mismatch_count);
        goto fail;
    }
    *out_midpoints = midpoint_count;
    return 0;
fail:
    if (fp != NULL) fclose(fp);
    wwl_lineage_dispose(out);
    return -1;
}

static int wwl_parse_size(const char *text, size_t *out)
{
    char *end = NULL;
    errno = 0;
    unsigned long long value = strtoull(text, &end, 10);
    if (errno || end == text || *end != '\0' || value > SIZE_MAX) return -1;
    *out = (size_t)value;
    return 0;
}

static int wwl_parse_double(const char *text, double *out)
{
    char *end = NULL;
    errno = 0;
    double value = strtod(text, &end);
    if (errno || end == text || *end != '\0' || !isfinite(value)) return -1;
    *out = value;
    return 0;
}

static int wwl_compare_double(const void *left, const void *right)
{
    double a = *(const double *)left;
    double b = *(const double *)right;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static double wwl_median(double *value, size_t n)
{
    if (n == 0) return NAN;
    qsort(value, n, sizeof(*value), wwl_compare_double);
    if ((n & 1u) != 0) return value[n / 2];
    return 0.5 * (value[n / 2 - 1] + value[n / 2]);
}

static double wwl_wrap_pi(double angle)
{
    angle = fmod(angle + WWL_PI, WWL_2PI);
    if (angle < 0.0) angle += WWL_2PI;
    return angle - WWL_PI;
}

static int wwl_dsu_init(WwlDsu *dsu, size_t n)
{
    memset(dsu, 0, sizeof(*dsu));
    if (n > (size_t)INT32_MAX) return -1;
    dsu->parent = (int32_t *)malloc((n ? n : 1) * sizeof(*dsu->parent));
    dsu->rank = (uint8_t *)calloc(n ? n : 1, sizeof(*dsu->rank));
    if (dsu->parent == NULL || dsu->rank == NULL) return -1;
    dsu->n = n;
    for (size_t i = 0; i < n; i++) dsu->parent[i] = (int32_t)i;
    return 0;
}

static void wwl_dsu_dispose(WwlDsu *dsu)
{
    free(dsu->parent);
    free(dsu->rank);
    memset(dsu, 0, sizeof(*dsu));
}

static int32_t wwl_dsu_find(WwlDsu *dsu, int32_t value)
{
    int32_t root = value;
    while (dsu->parent[root] != root) root = dsu->parent[root];
    while (dsu->parent[value] != value) {
        int32_t next = dsu->parent[value];
        dsu->parent[value] = root;
        value = next;
    }
    return root;
}

static void wwl_dsu_union(WwlDsu *dsu, int32_t a, int32_t b)
{
    int32_t ra = wwl_dsu_find(dsu, a);
    int32_t rb = wwl_dsu_find(dsu, b);
    if (ra == rb) return;
    if (dsu->rank[ra] < dsu->rank[rb]) {
        int32_t temporary = ra; ra = rb; rb = temporary;
    }
    dsu->parent[rb] = ra;
    if (dsu->rank[ra] == dsu->rank[rb]) dsu->rank[ra]++;
}

static int wwl_label_components(
    size_t nv, const int32_t *faces, size_t nf,
    int32_t **out_vertex_component, size_t *out_ncomponent,
    size_t **out_face_count)
{
    WwlDsu dsu;
    uint8_t *used = NULL;
    int32_t *root_component = NULL;
    int32_t *component = NULL;
    size_t *face_count = NULL;
    size_t ncomponent = 0;
    if (wwl_dsu_init(&dsu, nv) != 0) return -1;
    used = (uint8_t *)calloc(nv ? nv : 1, 1);
    root_component = (int32_t *)malloc((nv ? nv : 1) * sizeof(*root_component));
    component = (int32_t *)malloc((nv ? nv : 1) * sizeof(*component));
    if (used == NULL || root_component == NULL || component == NULL) goto fail;
    for (size_t i = 0; i < nv; i++) root_component[i] = component[i] = -1;
    for (size_t f = 0; f < nf; f++) {
        int32_t a = faces[f * 3], b = faces[f * 3 + 1], c = faces[f * 3 + 2];
        if (a < 0 || b < 0 || c < 0 ||
            (size_t)a >= nv || (size_t)b >= nv || (size_t)c >= nv)
            goto fail;
        used[a] = used[b] = used[c] = 1;
        wwl_dsu_union(&dsu, a, b);
        wwl_dsu_union(&dsu, b, c);
    }
    for (size_t v = 0; v < nv; v++) {
        if (!used[v]) continue;
        int32_t root = wwl_dsu_find(&dsu, (int32_t)v);
        if (root_component[root] < 0)
            root_component[root] = (int32_t)ncomponent++;
        component[v] = root_component[root];
    }
    face_count = (size_t *)calloc(ncomponent ? ncomponent : 1,
                                  sizeof(*face_count));
    if (face_count == NULL) goto fail;
    for (size_t f = 0; f < nf; f++) {
        int32_t c = component[faces[f * 3]];
        if (c < 0 || (size_t)c >= ncomponent) goto fail;
        face_count[c]++;
    }
    free(used);
    free(root_component);
    wwl_dsu_dispose(&dsu);
    *out_vertex_component = component;
    *out_ncomponent = ncomponent;
    if (out_face_count != NULL) *out_face_count = face_count;
    else free(face_count);
    return 0;
fail:
    free(used);
    free(root_component);
    free(component);
    free(face_count);
    wwl_dsu_dispose(&dsu);
    return -1;
}

static int wwl_graph_build(
    size_t nv, const int32_t *faces, size_t nf, WwlGraph *out)
{
    size_t *degree = NULL, *cursor = NULL;
    memset(out, 0, sizeof(*out));
    if (nf > SIZE_MAX / 6) return -1;
    degree = (size_t *)calloc(nv ? nv : 1, sizeof(*degree));
    if (degree == NULL) return -1;
    for (size_t f = 0; f < nf; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t a = faces[f * 3 + (size_t)k];
            int32_t b = faces[f * 3 + (size_t)((k + 1) % 3)];
            if (a < 0 || b < 0 || (size_t)a >= nv || (size_t)b >= nv)
                goto fail;
            degree[a]++;
            degree[b]++;
        }
    }
    out->offset = (size_t *)malloc((nv + 1) * sizeof(*out->offset));
    if (out->offset == NULL) goto fail;
    out->offset[0] = 0;
    for (size_t v = 0; v < nv; v++) {
        if (out->offset[v] > SIZE_MAX - degree[v]) goto fail;
        out->offset[v + 1] = out->offset[v] + degree[v];
    }
    out->neighbor = (int32_t *)malloc(
        (out->offset[nv] ? out->offset[nv] : 1) * sizeof(*out->neighbor));
    cursor = (size_t *)malloc((nv ? nv : 1) * sizeof(*cursor));
    if (out->neighbor == NULL || cursor == NULL) goto fail;
    memcpy(cursor, out->offset, nv * sizeof(*cursor));
    for (size_t f = 0; f < nf; f++) {
        for (int k = 0; k < 3; k++) {
            int32_t a = faces[f * 3 + (size_t)k];
            int32_t b = faces[f * 3 + (size_t)((k + 1) % 3)];
            out->neighbor[cursor[a]++] = b;
            out->neighbor[cursor[b]++] = a;
        }
    }
    out->nv = nv;
    free(degree);
    free(cursor);
    return 0;
fail:
    free(degree);
    free(cursor);
    free(out->offset);
    free(out->neighbor);
    memset(out, 0, sizeof(*out));
    return -1;
}

static void wwl_graph_dispose(WwlGraph *graph)
{
    free(graph->offset);
    free(graph->neighbor);
    memset(graph, 0, sizeof(*graph));
}

static int wwl_same_xyz(const float *a, const float *b)
{
    return memcmp(a, b, 3 * sizeof(float)) == 0;
}

static uint64_t wwl_position_hash(const float *point)
{
    uint32_t bits[3];
    memcpy(bits, point, sizeof(bits));
    uint64_t hash = 1469598103934665603ull;
    for (int axis = 0; axis < 3; axis++) {
        hash ^= bits[axis];
        hash *= 1099511628211ull;
    }
    hash ^= hash >> 33;
    hash *= 0xff51afd7ed558ccdull;
    hash ^= hash >> 33;
    return hash;
}

/* grid_weld appends fill/refinement vertices and its final compactor preserves
 * relative order, but a placed directory need not enumerate cube files in that
 * same order. Join the original block by exact float position instead. The
 * production geometry guarantees adjacent cubes do not share bitwise seam
 * vertices; any true duplicate is reported as ambiguous and never anchored. */
static int wwl_map_source_exact(
    const PieceSet *source, const float *welded, size_t welded_nv,
    int32_t *anchor_source, size_t *out_matched, size_t *out_missing,
    size_t *out_ambiguous)
{
    size_t capacity = 1;
    while (capacity < source->nv * 2 + 1) {
        if (capacity > SIZE_MAX / 2) return -1;
        capacity *= 2;
    }
    WwlPositionSlot *slot = (WwlPositionSlot *)calloc(
        capacity, sizeof(*slot));
    if (slot == NULL) return -1;
    for (size_t v = 0; v < welded_nv; v++) anchor_source[v] = -1;
    size_t source_ambiguous = 0;
    for (size_t s = 0; s < source->nv; s++) {
        const float *point = &source->verts[s * 3];
        uint32_t bits[3];
        memcpy(bits, point, sizeof(bits));
        size_t at = (size_t)wwl_position_hash(point) & (capacity - 1);
        while (slot[at].occupied &&
               memcmp(slot[at].xyz, bits, sizeof(bits)) != 0)
            at = (at + 1) & (capacity - 1);
        if (!slot[at].occupied) {
            slot[at].occupied = 1;
            memcpy(slot[at].xyz, bits, sizeof(bits));
            slot[at].source = (int32_t)s;
        } else if (slot[at].source >= 0) {
            slot[at].source = -2;
            source_ambiguous++;
        }
    }
    size_t matched = 0, final_ambiguous = 0;
    size_t scan = source->nv < welded_nv ? source->nv : welded_nv;
    for (size_t v = 0; v < scan; v++) {
        const float *point = &welded[v * 3];
        uint32_t bits[3];
        memcpy(bits, point, sizeof(bits));
        size_t at = (size_t)wwl_position_hash(point) & (capacity - 1);
        while (slot[at].occupied &&
               memcmp(slot[at].xyz, bits, sizeof(bits)) != 0)
            at = (at + 1) & (capacity - 1);
        if (slot[at].occupied && slot[at].source >= 0) {
            anchor_source[v] = slot[at].source;
            matched++;
        } else if (slot[at].occupied) {
            final_ambiguous++;
        }
    }
    free(slot);
    size_t missing = source->nv - matched;
    *out_matched = matched;
    *out_missing = missing;
    *out_ambiguous = final_ambiguous;
    if (source_ambiguous > 0)
        fprintf(stderr,
            "welded_winding_lift: %zu duplicate source position(s), "
            "%zu surviving vertices left unanchored\n",
            source_ambiguous, final_ambiguous);
    return source->nv > 0 && matched * 1000 >= source->nv * 999 ? 0 : -1;
}

/* GWLIN2 records source identity as (cube id, local vertex), so coincident
 * vertices from different charts remain distinct.  Every anchor is still
 * checked bit-for-bit against the welded geometry before its UV is accepted. */
static int wwl_map_source_lineage(
    const PieceSet *source, const float *welded, size_t welded_nv,
    const WwlLineage *lineage, int32_t *anchor_source,
    size_t *out_matched, size_t *out_missing, size_t *out_ambiguous)
{
    int32_t *cube_map = NULL;
    uint8_t *seen = NULL;
    size_t matched = 0, mismatch = 0, shared = 0, distinct = 0, adjusted = 0;
    double maximum_displacement = 0.0;
    const double allowed_source_move =
        (double)HOLEFILL_CHART_PINCH_MAX_STEP_VOX + 1.0e-3;
    if (lineage == NULL || lineage->version < 2 ||
        lineage->n_cubes == 0 || lineage->source_anchors == 0)
        return -1;
    cube_map = (int32_t *)malloc(
        lineage->n_cubes * sizeof(*cube_map));
    seen = (uint8_t *)calloc(source->nv ? source->nv : 1, 1);
    if (cube_map == NULL || seen == NULL) goto fail;
    for (size_t c = 0; c < lineage->n_cubes; c++) {
        cube_map[c] = -1;
        for (size_t s = 0; s < source->n_cubes; s++) {
            if (strcmp(lineage->cube_ids[c], source->ids[s]) != 0) continue;
            if (cube_map[c] >= 0) goto fail;
            cube_map[c] = (int32_t)s;
        }
        if (cube_map[c] < 0) {
            fprintf(stderr,
                "welded_winding_lift: lineage cube %s is absent from placed "
                "source\n", lineage->cube_ids[c]);
            goto fail;
        }
    }
    for (size_t v = 0; v < welded_nv; v++) anchor_source[v] = -1;
    for (size_t v = 0; v < welded_nv; v++) {
        int32_t p0 = lineage->parent0[v], local = lineage->parent1[v];
        if (p0 > -2) continue;
        int64_t lineage_cube = -(int64_t)p0 - 2;
        if (lineage_cube < 0 ||
            (uint64_t)lineage_cube >= (uint64_t)lineage->n_cubes ||
            local < 0) {
            mismatch++;
            continue;
        }
        size_t source_cube =
            (size_t)cube_map[(size_t)lineage_cube];
        size_t begin = source->cube_voff[source_cube];
        size_t end = source->cube_voff[source_cube + 1];
        if ((size_t)local >= end - begin) {
            mismatch++;
            continue;
        }
        size_t source_vertex = begin + (size_t)local;
        if (!wwl_same_xyz(&source->verts[source_vertex * 3],
                          &welded[v * 3])) {
            double delta[3];
            for (int axis = 0; axis < 3; axis++)
                delta[axis] =
                    (double)welded[v * 3 + (size_t)axis] -
                    (double)source->verts[source_vertex * 3 + (size_t)axis];
            double displacement = sqrt(
                delta[0] * delta[0] + delta[1] * delta[1] +
                delta[2] * delta[2]);
            if (displacement > maximum_displacement)
                maximum_displacement = displacement;
            if (displacement > allowed_source_move && mismatch < 16)
                fprintf(stderr,
                    "welded_winding_lift: moved source anchor final=%zu "
                    "cube=%s local=%d displacement=%.9g "
                    "source=(%.9g %.9g %.9g) final=(%.9g %.9g %.9g)\n",
                    v, lineage->cube_ids[(size_t)lineage_cube], local,
                    displacement,
                    (double)source->verts[source_vertex * 3],
                    (double)source->verts[source_vertex * 3 + 1],
                    (double)source->verts[source_vertex * 3 + 2],
                    (double)welded[v * 3],
                    (double)welded[v * 3 + 1],
                    (double)welded[v * 3 + 2]);
            if (displacement > allowed_source_move) {
                mismatch++;
                continue;
            }
            adjusted++;
        }
        /* grid_weld splits a source vertex whenever it peels a bowtie, cuts a
         * source chart to a disk, or resolves a non-manifold edge.  Both halves
         * carry the same (cube id, local vertex) provenance and, as the
         * bit-exact check above has just confirmed, sit at the very same
         * position -- so both must inherit that source vertex's UV.  Only a
         * geometry disagreement is an error; a shared parent is not. */
        if (seen[source_vertex]) {
            shared++;
        } else {
            seen[source_vertex] = 1;
            distinct++;
        }
        anchor_source[v] = (int32_t)source_vertex;
        matched++;
    }
    if (matched != lineage->source_anchors || mismatch != 0) {
        fprintf(stderr,
            "welded_winding_lift: explicit source provenance rejected "
            "(anchors=%zu/%zu geometry/id mismatch=%zu shared=%zu "
            "max_displacement=%.9g)\n",
            matched, lineage->source_anchors, mismatch, shared,
            maximum_displacement);
        goto fail;
    }
    if (shared > 0)
        fprintf(stderr,
            "welded_winding_lift: %zu welded vertex/vertices share %zu "
            "source parent(s) (weld-time vertex splits); each inherits "
            "its parent UV\n", shared, matched - distinct);
    if (adjusted > 0)
        fprintf(stderr,
            "welded_winding_lift: accepted %zu explicit source anchor(s) "
            "with bounded chart-pinch adjustment (max %.9g vox, cap %.9g)\n",
            adjusted, maximum_displacement, allowed_source_move);
    *out_matched = matched;
    *out_missing = source->nv > distinct ? source->nv - distinct : 0;
    *out_ambiguous = shared;
    free(cube_map);
    free(seen);
    return 0;
fail:
    free(cube_map);
    free(seen);
    return -1;
}

static int wwl_source_order_diagnostic(
    const PieceSet *source, const float *welded, size_t welded_nv,
    size_t *out_ordered)
{
    size_t next = 0, matched = 0;
    for (size_t s = 0; s < source->nv; s++) {
        if (next < welded_nv &&
            wwl_same_xyz(&source->verts[s * 3], &welded[next * 3])) {
            next++;
            matched++;
        } else {
            continue;
        }
    }
    *out_ordered = matched;
    return 0;
}

/* out_phi (nullable): carry source->phi through the SAME anchor / midpoint /
 * BFS / harmonic machinery as uv.  phi here is the lifted continuous phase
 * (2pi*k already applied by registration), so plain averaging is correct --
 * there is no wrap boundary inside a chart. */
static int wwl_interpolate_chart_uv(
    const float *verts, size_t nv, const WwlGraph *graph,
    const int32_t *chart, const int32_t *anchor_source,
    const PieceSet *source, const WwlLineage *lineage,
    size_t sweeps, double *uv, float *out_phi, size_t *out_unknown,
    size_t *out_midpoint_fixed, size_t *out_residual,
    double *out_last_delta)
{
    uint8_t *fixed = (uint8_t *)calloc(nv ? nv : 1, 1);
    int32_t *queue = (int32_t *)malloc((nv ? nv : 1) * sizeof(*queue));
    if (fixed == NULL || queue == NULL) { free(fixed); free(queue); return -1; }
    size_t head = 0, tail = 0, unknown = 0;
    for (size_t v = 0; v < nv; v++) {
        if (chart[v] < 0) continue;
        if (anchor_source[v] >= 0) {
            size_t s = (size_t)anchor_source[v];
            uv[v * 2] = source->uv[s * 2];
            uv[v * 2 + 1] = source->uv[s * 2 + 1];
            if (out_phi != NULL) out_phi[v] = source->phi[s];
            fixed[v] = 1;
            queue[tail++] = (int32_t)v;
        } else {
            uv[v * 2] = uv[v * 2 + 1] = NAN;
            unknown++;
        }
    }
    size_t midpoint_fixed = 0;
    if (lineage != NULL) {
        for (size_t v = 0; v < nv; v++) {
            int32_t p0 = lineage->parent0[v], p1 = lineage->parent1[v];
            if (chart[v] < 0 || fixed[v] || p0 < 0 || p1 < 0) continue;
            if (!fixed[(size_t)p0] || !fixed[(size_t)p1] ||
                chart[(size_t)p0] != chart[v] || chart[(size_t)p1] != chart[v])
                continue;
            uv[v * 2] = 0.5 *
                (uv[(size_t)p0 * 2] + uv[(size_t)p1 * 2]);
            uv[v * 2 + 1] = 0.5 *
                (uv[(size_t)p0 * 2 + 1] + uv[(size_t)p1 * 2 + 1]);
            if (out_phi != NULL)
                out_phi[v] = 0.5f * (out_phi[(size_t)p0]
                                     + out_phi[(size_t)p1]);
            fixed[v] = 1;
            queue[tail++] = (int32_t)v;
            midpoint_fixed++;
        }
    }
    size_t residual = unknown - midpoint_fixed;
    /* Nearest-anchor initialization keeps deep refinement chains finite before
     * the deterministic harmonic sweeps begin.  With exact lineage this is
     * restricted to original duplicate positions and chart-fill vertices; the
     * safe winding subdivision itself remains fixed bit-for-bit in UV. */
    while (head < tail) {
        int32_t from = queue[head++];
        for (size_t at = graph->offset[from]; at < graph->offset[from + 1]; at++) {
            int32_t to = graph->neighbor[at];
            if (chart[to] != chart[from] || isfinite(uv[(size_t)to * 2])) continue;
            uv[(size_t)to * 2] = uv[(size_t)from * 2];
            uv[(size_t)to * 2 + 1] = uv[(size_t)from * 2 + 1];
            if (out_phi != NULL) out_phi[to] = out_phi[from];
            queue[tail++] = to;
        }
    }
    for (size_t v = 0; v < nv; v++)
        if (chart[v] >= 0 && !isfinite(uv[v * 2])) {
            free(fixed); free(queue); return -1;
        }
    /* Sweep a worklist of the FREE vertices only.  The old loops scanned all
     * nv per sweep; at 124M welded vertices with 231k free ones and the fixed
     * 160-sweep cap (the residual plateaus above the 1e-7 exit), that was
     * 160 full-mesh scans ~= 10+ minutes of pre-fit time on the 21x21x21.
     * The list is built in ascending v order, so the relaxation order and the
     * result are bit-identical to the full scan. */
    size_t nfree = 0;
    int32_t *free_list = (int32_t *)malloc(
        (residual ? residual : 1) * sizeof(*free_list));
    if (free_list == NULL) { free(fixed); free(queue); return -1; }
    for (size_t v = 0; v < nv; v++)
        if (chart[v] >= 0 && !fixed[v]) {
            if (nfree < (residual ? residual : 1))
                free_list[nfree] = (int32_t)v;
            nfree++;
        }
    if (nfree > residual) {   /* defensive: never overrun the list */
        free(free_list);
        free(fixed); free(queue);
        return -1;
    }
    double last_delta = 0.0;
    for (size_t sweep = 0; sweep < sweeps; sweep++) {
        double max_delta = 0.0;
        for (size_t t = 0; t < nfree; t++) {
            size_t v = (size_t)free_list[t];
            double sum_u = 0.0, sum_v = 0.0, sum_p = 0.0, sum_w = 0.0;
            for (size_t at = graph->offset[v]; at < graph->offset[v + 1]; at++) {
                int32_t n = graph->neighbor[at];
                if (chart[n] != chart[v]) continue;
                double dz = (double)verts[v * 3] - verts[(size_t)n * 3];
                double dy = (double)verts[v * 3 + 1] - verts[(size_t)n * 3 + 1];
                double dx = (double)verts[v * 3 + 2] - verts[(size_t)n * 3 + 2];
                double length = sqrt(dz * dz + dy * dy + dx * dx);
                double weight = 1.0 / fmax(length, 1.0e-6);
                sum_u += weight * uv[(size_t)n * 2];
                sum_v += weight * uv[(size_t)n * 2 + 1];
                if (out_phi != NULL) sum_p += weight * (double)out_phi[n];
                sum_w += weight;
            }
            if (sum_w > 0.0) {
                double next_u = sum_u / sum_w;
                double next_v = sum_v / sum_w;
                double delta = hypot(next_u - uv[v * 2], next_v - uv[v * 2 + 1]);
                if (delta > max_delta) max_delta = delta;
                uv[v * 2] = next_u;
                uv[v * 2 + 1] = next_v;
                if (out_phi != NULL) out_phi[v] = (float)(sum_p / sum_w);
            }
        }
        last_delta = max_delta;
        if (max_delta < 1.0e-7) break;
    }
    free(free_list);
    free(fixed);
    free(queue);
    *out_unknown = unknown;
    *out_midpoint_fixed = midpoint_fixed;
    *out_residual = residual;
    *out_last_delta = last_delta;
    return 0;
}

/* --carry-registered support: the registered per-cube field is authoritative,
 * so anchors that disagree wildly with their own chart are provenance noise
 * (the 0.02% u-tail that stretched the 4x21x21 span 1.196M -> 1.515M and
 * inflated every raster).  Demote anchors farther than max_dev from their
 * chart's median anchor u; the median element itself always survives, and
 * charts with fewer than 8 anchors are left alone. Returns demoted count. */
typedef struct {
    double phi, u;
    int32_t vert;
} WwlPhiAnchor;

static int wwl_compare_phi_anchor(const void *pa, const void *pb)
{
    const WwlPhiAnchor *a = (const WwlPhiAnchor *)pa;
    const WwlPhiAnchor *b = (const WwlPhiAnchor *)pb;
    if (a->phi != b->phi) return a->phi < b->phi ? -1 : 1;
    if (a->u != b->u) return a->u < b->u ? -1 : 1;
    return a->vert < b->vert ? -1 : a->vert > b->vert ? 1 : 0;
}

static size_t wwl_carry_winsorize_anchors(
    size_t nv, const int32_t *chart, size_t nchart,
    int32_t *anchor_source, const PieceSet *source, double max_dev)
{
    size_t *count = (size_t *)calloc(nchart ? nchart : 1, sizeof(*count));
    size_t *offset = (size_t *)calloc(nchart + 1, sizeof(*offset));
    size_t demoted = 0;
    if (count == NULL || offset == NULL) { free(count); free(offset); return 0; }
    for (size_t v = 0; v < nv; v++)
        if (chart[v] >= 0 && anchor_source[v] >= 0)
            count[chart[v]]++;
    for (size_t c = 0; c < nchart; c++)
        offset[c + 1] = offset[c] + count[c];
    size_t total = offset[nchart];
    double *au = (double *)malloc((total ? total : 1) * sizeof(*au));
    double *ap = (double *)malloc((total ? total : 1) * sizeof(*ap));
    int32_t *av = (int32_t *)malloc((total ? total : 1) * sizeof(*av));
    size_t *cursor = (size_t *)malloc((nchart ? nchart : 1) * sizeof(*cursor));
    double *med = (double *)malloc((nchart ? nchart : 1) * sizeof(*med));
    double *pmed = (double *)malloc((nchart ? nchart : 1) * sizeof(*pmed));
    double *beta = (double *)malloc((nchart ? nchart : 1) * sizeof(*beta));
    if (au == NULL || ap == NULL || av == NULL || cursor == NULL ||
        med == NULL || pmed == NULL || beta == NULL) {
        free(count); free(offset); free(au); free(ap); free(av);
        free(cursor); free(med); free(pmed); free(beta);
        return 0;
    }
    memcpy(cursor, offset, nchart * sizeof(*cursor));
    for (size_t v = 0; v < nv; v++) {
        if (chart[v] < 0 || anchor_source[v] < 0) continue;
        size_t c = (size_t)chart[v];
        au[cursor[c]] = source->uv[(size_t)anchor_source[v] * 2];
        ap[cursor[c]] = source->phi[(size_t)anchor_source[v]];
        av[cursor[c]] = (int32_t)v;
        cursor[c]++;
    }
    /* Winsorize per (chart, PHI-BIN).  Neither a chart median nor any
     * per-chart u(phi) FIT survives the real data: the welded mega-sheet
     * winds 159 turns (58k anchors, u span 751k vox), where arc length is
     * strongly quadratic in phi -- a raw median gate demoted everything
     * off one wrap, and a robust linear fit still shed 10.5k true anchors
     * to curvature.  Registered phi separates wraps by construction, so
     * within a small phi window a correct chart occupies ONE tight u
     * range regardless of how many turns the chart makes overall; genuine
     * weld artifacts remain far from their bin's median.  Bin width
     * adapts so sparse fragments still fill bins (~8 anchors), capped so
     * the in-bin arc stays under the gate. */
    size_t *cdem = (size_t *)calloc(nchart ? nchart : 1, sizeof(*cdem));
    for (size_t c = 0; c < nchart; c++) {
        size_t n = count[c];
        med[c] = 0.0;
        pmed[c] = 0.0;
        beta[c] = 0.0;
        if (n < 8) continue;
        WwlPhiAnchor *slice = (WwlPhiAnchor *)malloc(
            n * sizeof(*slice));
        double *vals = (double *)malloc(n * sizeof(*vals));
        if (slice == NULL || vals == NULL) {
            free(slice); free(vals);
            continue;
        }
        for (size_t k = 0; k < n; k++) {
            slice[k].phi = ap[offset[c] + k];
            slice[k].u = au[offset[c] + k];
            slice[k].vert = av[offset[c] + k];
        }
        qsort(slice, n, sizeof(*slice), wwl_compare_phi_anchor);
        double span = slice[n - 1].phi - slice[0].phi;
        double width = span * 8.0 / (double)n;
        if (width < 0.5) width = 0.5;
        if (width > 4.0) width = 4.0;
        size_t b0 = 0;
        while (b0 < n) {
            size_t b1 = b0;
            while (b1 < n && slice[b1].phi < slice[b0].phi + width) b1++;
            size_t bn = b1 - b0;
            if (bn >= 5) {
                for (size_t k = 0; k < bn; k++) vals[k] = slice[b0 + k].u;
                qsort(vals, bn, sizeof(*vals), wwl_compare_double);
                double bmed = vals[bn / 2];
                for (size_t k = b0; k < b1; k++) {
                    if (fabs(slice[k].u - bmed) > max_dev) {
                        anchor_source[slice[k].vert] = -1;
                        demoted++;
                        if (cdem != NULL) cdem[c]++;
                    }
                }
            }
            b0 = b1;
        }
        free(slice);
        free(vals);
    }
    /* diagnostic: the charts that lose the most anchors tell us WHAT the
     * demoted population is (fusion branches, scattered lineage, failed
     * beta fits ...) */
    if (cdem != NULL && demoted > 0) {
        for (int top = 0; top < 12; top++) {
            size_t best = 0, bc = nchart;
            for (size_t c = 0; c < nchart; c++)
                if (cdem[c] > best) { best = cdem[c]; bc = c; }
            if (bc == nchart || best == 0) break;
            double ulo = 1e300, uhi = -1e300, plo = 1e300, phi_hi = -1e300;
            for (size_t k = offset[bc]; k < offset[bc + 1]; k++) {
                if (au[k] < ulo) ulo = au[k];
                if (au[k] > uhi) uhi = au[k];
                if (ap[k] < plo) plo = ap[k];
                if (ap[k] > phi_hi) phi_hi = ap[k];
            }
            fprintf(stderr, "  [winsorize] chart %zu: anchors=%zu "
                    "demoted=%zu u=[%.0f,%.0f] "
                    "phi=[%.2f,%.2f] (%.2f turns)\n",
                    bc, count[bc], cdem[bc], ulo, uhi,
                    plo, phi_hi, (phi_hi - plo) / 6.283185307);
            cdem[bc] = 0;
        }
    }
    free(cdem);
    free(count); free(offset); free(au); free(ap); free(av);
    free(cursor); free(med); free(pmed); free(beta);
    return demoted;
}

/* --carry-registered support: a chart made entirely of fill/refinement
 * vertices has no registered anchor and would fail interpolation.  Seed it by
 * copying uv from anchored neighbours across the FULL graph (bridge faces
 * included), iterating so chains of anchorless charts drain outward.  Charts
 * still dry after that are counted in *out_dry (interpolation will then fail
 * loudly rather than silently). */
static void wwl_carry_seed_anchorless(
    size_t nv, const WwlGraph *full_graph, const int32_t *chart,
    size_t nchart, int32_t *anchor_source, size_t *out_seeded_charts,
    size_t *out_dry_charts)
{
    uint8_t *has_anchor = (uint8_t *)calloc(nchart ? nchart : 1, 1);
    *out_seeded_charts = 0;
    *out_dry_charts = 0;
    if (has_anchor == NULL) return;
    for (size_t v = 0; v < nv; v++)
        if (chart[v] >= 0 && anchor_source[v] >= 0)
            has_anchor[chart[v]] = 1;
    for (int round = 0; round < 8; round++) {
        size_t new_seeds = 0;
        for (size_t v = 0; v < nv; v++) {
            if (chart[v] < 0 || anchor_source[v] < 0) continue;
            for (size_t at = full_graph->offset[v];
                 at < full_graph->offset[v + 1]; at++) {
                int32_t n = full_graph->neighbor[at];
                if (chart[n] < 0 || has_anchor[chart[n]] ||
                    anchor_source[n] >= 0)
                    continue;
                anchor_source[n] = anchor_source[v];
                new_seeds++;
            }
        }
        if (new_seeds == 0) break;
        for (size_t v = 0; v < nv; v++)
            if (chart[v] >= 0 && anchor_source[v] >= 0 &&
                !has_anchor[chart[v]]) {
                has_anchor[chart[v]] = 1;
                (*out_seeded_charts)++;
            }
    }
    for (size_t c = 0; c < nchart; c++)
        if (!has_anchor[c]) (*out_dry_charts)++;
    free(has_anchor);
}

static int wwl_compare_observation(const void *left, const void *right)
{
    const WwlChartObservation *a = (const WwlChartObservation *)left;
    const WwlChartObservation *b = (const WwlChartObservation *)right;
    return a->chart < b->chart ? -1 : (a->chart > b->chart ? 1 : 0);
}

static int wwl_fit_chart_isometry_candidate(
    const WwlChartObservation *observation, size_t first, size_t last,
    int reflected, double *scratch, double matrix[4], double translate[2],
    double *out_median, double *out_maximum, int *out_point_gauge)
{
    size_t count = last - first;
    if (count == 0) return -1;
    double sx = 0.0, sy = 0.0, tx = 0.0, ty = 0.0;
    for (size_t i = first; i < last; i++) {
        sx += observation[i].source_u;
        sy += observation[i].source_v;
        tx += observation[i].target_u;
        ty += observation[i].target_v;
    }
    sx /= (double)count; sy /= (double)count;
    tx /= (double)count; ty /= (double)count;
    double a = 0.0, b = 0.0, spread = 0.0;
    for (size_t i = first; i < last; i++) {
        double x = observation[i].source_u - sx;
        double y = observation[i].source_v - sy;
        double u = observation[i].target_u - tx;
        double v = observation[i].target_v - ty;
        if (reflected) {
            a += x * u - y * v;
            b += y * u + x * v;
        } else {
            a += x * u + y * v;
            b += x * v - y * u;
        }
        spread += x * x + y * y;
    }
    int point_gauge = count < 2 || spread <= 1.0e-12 ||
                      hypot(a, b) <= 1.0e-12;
    double cosine = 1.0, sine = 0.0;
    if (!point_gauge) {
        double norm = hypot(a, b);
        cosine = a / norm;
        sine = b / norm;
    } else if (reflected) {
        return -1;
    }
    if (reflected) {
        matrix[0] = cosine; matrix[1] = sine;
        matrix[2] = sine;   matrix[3] = -cosine;
    } else {
        matrix[0] = cosine; matrix[1] = -sine;
        matrix[2] = sine;   matrix[3] = cosine;
    }
    for (size_t i = first; i < last; i++)
        scratch[i - first] = observation[i].target_u -
            (matrix[0] * observation[i].source_u +
             matrix[1] * observation[i].source_v);
    translate[0] = wwl_median(scratch, count);
    for (size_t i = first; i < last; i++)
        scratch[i - first] = observation[i].target_v -
            (matrix[2] * observation[i].source_u +
             matrix[3] * observation[i].source_v);
    translate[1] = wwl_median(scratch, count);
    double maximum = 0.0;
    for (size_t i = first; i < last; i++) {
        double u = matrix[0] * observation[i].source_u +
                   matrix[1] * observation[i].source_v + translate[0];
        double v = matrix[2] * observation[i].source_u +
                   matrix[3] * observation[i].source_v + translate[1];
        double residual = hypot(u - observation[i].target_u,
                                v - observation[i].target_v);
        scratch[i - first] = residual;
        if (residual > maximum) maximum = residual;
    }
    *out_median = wwl_median(scratch, count);
    *out_maximum = maximum;
    *out_point_gauge = point_gauge;
    return 0;
}

/* Fit one Euclidean isometry per raw source chart to a globally meaningful
 * target coordinate.  Rotations and reflections are gauges, not deformation:
 * every raw edge length, triangle area, and injectivity certificate survives
 * bit-for-bit up to floating-point evaluation. */
static int wwl_apply_target_chart_isometries(
    size_t nv, const int32_t *chart, size_t nchart,
    const int32_t *anchor_source, const double *target_uv, double *uv,
    WwlChartGaugeStats *stats)
{
    WwlChartObservation *observation = (WwlChartObservation *)malloc(
        (nv ? nv : 1) * sizeof(*observation));
    double *matrix = (double *)malloc(
        (nchart ? nchart : 1) * 4 * sizeof(*matrix));
    double *translate = (double *)malloc(
        (nchart ? nchart : 1) * 2 * sizeof(*translate));
    double *scratch = (double *)malloc((nv ? nv : 1) * sizeof(*scratch));
    double *all_residual = (double *)malloc(
        (nv ? nv : 1) * sizeof(*all_residual));
    if (observation == NULL || matrix == NULL || translate == NULL ||
        scratch == NULL || all_residual == NULL) goto fail;
    if (stats != NULL) memset(stats, 0, sizeof(*stats));

    size_t n = 0;
    for (size_t v = 0; v < nv; v++) {
        int32_t c = chart[v];
        if (c < 0 || anchor_source[v] < 0) continue;
        if ((size_t)c >= nchart || !isfinite(target_uv[v * 2]) ||
            !isfinite(target_uv[v * 2 + 1]) ||
            !isfinite(uv[v * 2]) || !isfinite(uv[v * 2 + 1]))
            goto fail;
        observation[n].chart = c;
        observation[n].source_u = uv[v * 2];
        observation[n].source_v = uv[v * 2 + 1];
        observation[n].target_u = target_uv[v * 2];
        observation[n].target_v = target_uv[v * 2 + 1];
        n++;
    }
    if (n == 0) goto fail;
    qsort(observation, n, sizeof(*observation), wwl_compare_observation);
    for (size_t c = 0; c < nchart; c++) matrix[c * 4] = NAN;

    size_t reflected_charts = 0, point_gauge_charts = 0;
    for (size_t first = 0; first < n;) {
        size_t last = first + 1;
        while (last < n && observation[last].chart == observation[first].chart)
            last++;
        int32_t c = observation[first].chart;
        double proper_m[4], proper_t[2], proper_med, proper_max;
        double reflect_m[4], reflect_t[2], reflect_med = DBL_MAX;
        double reflect_max = DBL_MAX;
        int proper_point = 0, reflect_point = 0;
        if (wwl_fit_chart_isometry_candidate(
                observation, first, last, 0, scratch, proper_m, proper_t,
                &proper_med, &proper_max, &proper_point) != 0)
            goto fail;
        int have_reflection = wwl_fit_chart_isometry_candidate(
                observation, first, last, 1, scratch, reflect_m, reflect_t,
                &reflect_med, &reflect_max, &reflect_point) == 0;
        int use_reflection = have_reflection &&
            (reflect_med + 1.0e-9 < proper_med ||
             (fabs(reflect_med - proper_med) <= 1.0e-9 &&
              reflect_max < proper_max));
        const double *chosen_m = use_reflection ? reflect_m : proper_m;
        const double *chosen_t = use_reflection ? reflect_t : proper_t;
        memcpy(&matrix[(size_t)c * 4], chosen_m, 4 * sizeof(double));
        memcpy(&translate[(size_t)c * 2], chosen_t, 2 * sizeof(double));
        if (use_reflection) reflected_charts++;
        if (proper_point) point_gauge_charts++;
        first = last;
    }
    for (size_t c = 0; c < nchart; c++)
        if (!isfinite(matrix[c * 4])) goto fail;
    for (size_t v = 0; v < nv; v++) {
        int32_t c = chart[v];
        if (c < 0 || (size_t)c >= nchart) goto fail;
        const double *m = &matrix[(size_t)c * 4];
        const double *t = &translate[(size_t)c * 2];
        double u = uv[v * 2], w = uv[v * 2 + 1];
        uv[v * 2] = m[0] * u + m[1] * w + t[0];
        uv[v * 2 + 1] = m[2] * u + m[3] * w + t[1];
    }
    double maximum = 0.0;
    for (size_t i = 0; i < n; i++) {
        int32_t c = observation[i].chart;
        const double *m = &matrix[(size_t)c * 4];
        const double *t = &translate[(size_t)c * 2];
        double u = m[0] * observation[i].source_u +
                   m[1] * observation[i].source_v + t[0];
        double v = m[2] * observation[i].source_u +
                   m[3] * observation[i].source_v + t[1];
        double residual = hypot(u - observation[i].target_u,
                                v - observation[i].target_v);
        all_residual[i] = residual;
        if (residual > maximum) maximum = residual;
    }
    if (stats != NULL) {
        stats->observations = n;
        stats->reflected_charts = reflected_charts;
        stats->point_gauge_charts = point_gauge_charts;
        stats->residual_median = wwl_median(all_residual, n);
        stats->residual_maximum = maximum;
    }
    free(observation); free(matrix); free(translate); free(scratch);
    free(all_residual);
    return 0;

fail:
    free(observation); free(matrix); free(translate); free(scratch);
    free(all_residual);
    return -1;
}

/* Transfer only the rigid gauge of the registered winding solution onto the
 * safe pre-registration charts. Reregistration may contain a vertex-wise warp;
 * copying it directly can fold a chart, while a chart isometry cannot. */
static int wwl_apply_registered_chart_gauges(
    size_t nv, size_t source_nv, const int32_t *chart, size_t nchart,
    const int32_t *anchor_source, const float *registered_uv, double *uv,
    WwlChartGaugeStats *stats)
{
    double *target = (double *)malloc((nv ? nv : 1) * 2 * sizeof(*target));
    if (target == NULL) return -1;
    for (size_t v = 0; v < nv; v++) {
        int32_t s = anchor_source[v];
        target[v * 2] = target[v * 2 + 1] = NAN;
        if (s < 0) continue;
        if ((size_t)s >= source_nv) { free(target); return -1; }
        target[v * 2] = registered_uv[(size_t)s * 2];
        target[v * 2 + 1] = registered_uv[(size_t)s * 2 + 1];
    }
    int rc = wwl_apply_target_chart_isometries(
        nv, chart, nchart, anchor_source, target, uv, stats);
    free(target);
    return rc;
}

static int wwl_compare_pack(const void *left, const void *right)
{
    const WwlPackRecord *a = (const WwlPackRecord *)left;
    const WwlPackRecord *b = (const WwlPackRecord *)right;
    if (a->faces != b->faces) return a->faces > b->faces ? -1 : 1;
    return a->component < b->component ? -1 :
           (a->component > b->component ? 1 : 0);
}

static int wwl_compare_pack_umin(const void *left, const void *right)
{
    const WwlPackRecord *a = (const WwlPackRecord *)left;
    const WwlPackRecord *b = (const WwlPackRecord *)right;
    if (a->umin != b->umin) return a->umin < b->umin ? -1 : 1;
    if (a->umax != b->umax) return a->umax < b->umax ? -1 : 1;
    return a->component < b->component ? -1 :
           (a->component > b->component ? 1 : 0);
}

static void wwl_face_orientation_stats(
    const double *uv, const int32_t *faces, size_t nf,
    size_t *positive, size_t *negative, size_t *degenerate);
static double wwl_area2_points(const double a[2], const double b[2],
                               const double c[2]);

static int wwl_assemble_bridge_gauges(
    const float *verts, size_t nv, const int32_t *faces,
    size_t source_nf, size_t nf, const int32_t *chart, size_t nchart,
    size_t expected_components, int allow_rotation, int apply_translation,
    int axial_only, double *uv,
    WwlBridgeGaugeStats *stats, const WwlBridgeAudit *audit);

/* Select one integer turn per certified source chart from the topology-lifted
 * polar phase, then apply the exact nonlinear winding correction used by the
 * shipped scroll_whole workflow.  Unlike a chart isometry, this keeps winding
 * U as winding U: no chart may rotate its axial direction into the scroll
 * direction merely to fit a coarse global target. */
static int wwl_apply_topology_chart_turns(
    size_t nv, const int32_t *chart, size_t nchart,
    const int32_t *anchor_source, const PieceSet *source,
    const ScaffoldCalib *calib, const double *lift, double *uv,
    WwlTopologyGaugeStats *stats)
{
    int result = -1;
    int32_t *turn = NULL;
    double *scratch = NULL;
    size_t residual_count = 0;
    WwlTopologyGaugeStats local;
    memset(&local, 0, sizeof(local));
    local.minimum_turn = INT32_MAX;
    local.maximum_turn = INT32_MIN;
    if (nv == 0 || nchart == 0 || chart == NULL || anchor_source == NULL ||
        source == NULL || source->phi == NULL || calib == NULL ||
        lift == NULL || uv == NULL)
        return -1;
    turn = (int32_t *)malloc(nchart * sizeof(*turn));
    scratch = (double *)malloc(nv * sizeof(*scratch));
    if (turn == NULL || scratch == NULL) goto cleanup;

    for (size_t c = 0; c < nchart; c++) {
        size_t count = 0;
        for (size_t v = 0; v < nv; v++) {
            int32_t s = anchor_source[v];
            if (chart[v] != (int32_t)c || s < 0) continue;
            if ((size_t)s >= source->nv || !isfinite(lift[v])) goto cleanup;
            scratch[count++] =
                (lift[v] - (double)source->phi[(size_t)s]) / WWL_2PI;
        }
        if (count == 0) goto cleanup;
        double median_turn = wwl_median(scratch, count);
        long long rounded = llround(median_turn);
        if (rounded < INT32_MIN || rounded > INT32_MAX) goto cleanup;
        turn[c] = (int32_t)rounded;
        if (turn[c] < local.minimum_turn) local.minimum_turn = turn[c];
        if (turn[c] > local.maximum_turn) local.maximum_turn = turn[c];
        local.anchors += count;
    }

    for (size_t v = 0; v < nv; v++) {
        int32_t c = chart[v], s = anchor_source[v];
        if (c < 0 || (size_t)c >= nchart || !isfinite(lift[v])) goto cleanup;
        if (s >= 0) {
            if ((size_t)s >= source->nv) goto cleanup;
            double delta_turn =
                (lift[v] - (double)source->phi[(size_t)s]) / WWL_2PI;
            double residual = fabs(delta_turn - (double)turn[c]);
            scratch[residual_count++] = residual;
            if (residual > local.residual_maximum_turns)
                local.residual_maximum_turns = residual;
        }
    }
    if (residual_count != local.anchors || residual_count == 0) goto cleanup;
    local.residual_median_turns = wwl_median(scratch, residual_count);

    for (size_t v = 0; v < nv; v++) {
        int32_t c = chart[v], s = anchor_source[v];
        double raw_phi = s >= 0
            ? (double)source->phi[(size_t)s]
            : lift[v] - WWL_2PI * (double)turn[c];
        uv[v * 2] += CubeReg_deltaU(
            calib->spiral_a, calib->spiral_b, raw_phi, turn[c]);
    }
    local.charts = nchart;
    if (stats != NULL) *stats = local;
    result = 0;

cleanup:
    free(turn);
    free(scratch);
    return result;
}

static int wwl_lift_components(
    const float *verts, size_t nv, const int32_t *faces,
    size_t source_nf, size_t nf, const WwlGraph *full_graph,
    const int32_t *component, size_t ncomponent,
    const int32_t *chart, size_t nchart,
    const int32_t *anchor_source, const PieceSet *source,
    const ScaffoldCalib *calib, double *uv,
    int topology_mode,
    int *out_sign, double *out_offset, double *out_cycle_error,
    WwlChartGaugeStats *out_gauge, WwlTopologyGaugeStats *out_topology,
    WwlBridgeGaugeStats *out_bridge, float *out_phase)
{
    double cp = 0.0, sp = 0.0, cm = 0.0, sm = 0.0;
    size_t phase_samples = 0;
    double *theta = (double *)malloc((nv ? nv : 1) * sizeof(*theta));
    double *lift = (double *)malloc((nv ? nv : 1) * sizeof(*lift));
    int32_t *queue = (int32_t *)malloc((nv ? nv : 1) * sizeof(*queue));
    uint8_t *seen = (uint8_t *)calloc(nv ? nv : 1, 1);
    double *scratch = (double *)malloc((nv ? nv : 1) * sizeof(*scratch));
    size_t *component_count = (size_t *)calloc(
        ncomponent ? ncomponent : 1, sizeof(*component_count));
    size_t *component_offset = (size_t *)malloc(
        (ncomponent + 1) * sizeof(*component_offset));
    size_t *component_cursor = (size_t *)malloc(
        (ncomponent ? ncomponent : 1) * sizeof(*component_cursor));
    double *component_turn = (double *)malloc(
        (ncomponent ? ncomponent : 1) * sizeof(*component_turn));
    double *target_u = (double *)malloc((nv ? nv : 1) * sizeof(*target_u));
    double *target_v = (double *)malloc((nv ? nv : 1) * sizeof(*target_v));
    if (theta == NULL || lift == NULL || queue == NULL || seen == NULL ||
        scratch == NULL || component_count == NULL ||
        component_offset == NULL || component_cursor == NULL ||
        component_turn == NULL || target_u == NULL || target_v == NULL)
        goto fail;

    for (size_t v = 0; v < nv; v++) {
        double raw = atan2((double)verts[v * 3 + 1] - calib->axis_point[1],
                           (double)verts[v * 3 + 2] - calib->axis_point[2]);
        theta[v] = raw;
        if (anchor_source[v] >= 0) {
            double phi = source->phi[(size_t)anchor_source[v]];
            cp += cos(phi - raw); sp += sin(phi - raw);
            cm += cos(phi + raw); sm += sin(phi + raw);
            phase_samples++;
        }
    }
    if (phase_samples == 0) goto fail;
    double coherence_plus = hypot(cp, sp) / (double)phase_samples;
    double coherence_minus = hypot(cm, sm) / (double)phase_samples;
    int sign = coherence_plus >= coherence_minus ? 1 : -1;
    double offset = sign > 0 ? atan2(sp, cp) : atan2(sm, cm);
    for (size_t v = 0; v < nv; v++) {
        theta[v] = (double)sign * theta[v] + offset;
        lift[v] = NAN;
    }

    double maximum_cycle_error = 0.0;
    for (size_t seed = 0; seed < nv; seed++) {
        if (component[seed] < 0 || seen[seed]) continue;
        size_t head = 0, tail = 0;
        seen[seed] = 1;
        lift[seed] = theta[seed];
        queue[tail++] = (int32_t)seed;
        while (head < tail) {
            int32_t from = queue[head++];
            for (size_t at = full_graph->offset[from];
                 at < full_graph->offset[from + 1]; at++) {
                int32_t to = full_graph->neighbor[at];
                if (component[to] != component[from]) continue;
                double candidate = lift[from] +
                                   wwl_wrap_pi(theta[to] - theta[from]);
                if (!seen[to]) {
                    seen[to] = 1;
                    lift[to] = candidate;
                    queue[tail++] = to;
                } else {
                    double error = fabs(candidate - lift[to]);
                    double turns = nearbyint(error / WWL_2PI);
                    error = fabs(error - turns * WWL_2PI);
                    if (error > maximum_cycle_error)
                        maximum_cycle_error = error;
                }
            }
        }
    }

    /* One whole-turn constant per welded sheet. Group anchors in two linear
     * passes instead of rescanning every vertex once per component; the latter
     * became the dominant cost on full-scroll grids with hundreds of sheets. */
    for (size_t v = 0; v < nv; v++) {
        int32_t c = component[v];
        if (c < 0 || (size_t)c >= ncomponent) goto fail;
        if (anchor_source[v] >= 0) component_count[(size_t)c]++;
    }
    component_offset[0] = 0;
    for (size_t c = 0; c < ncomponent; c++) {
        if (component_count[c] == 0 ||
            component_offset[c] > nv - component_count[c])
            goto fail;
        component_offset[c + 1] = component_offset[c] + component_count[c];
        component_cursor[c] = component_offset[c];
    }
    for (size_t v = 0; v < nv; v++) {
        int32_t c = component[v], s = anchor_source[v];
        if (s >= 0)
            scratch[component_cursor[(size_t)c]++] =
                (source->phi[(size_t)s] - lift[v]) / WWL_2PI;
    }
    for (size_t c = 0; c < ncomponent; c++)
        component_turn[c] = nearbyint(wwl_median(
            &scratch[component_offset[c]], component_count[c]));
    for (size_t v = 0; v < nv; v++)
        lift[v] += component_turn[(size_t)component[v]] * WWL_2PI;
    if (out_phase != NULL)
        for (size_t v = 0; v < nv; v++) out_phase[v] = (float)lift[v];

    if (topology_mode == 1) {
        if (wwl_apply_topology_chart_turns(
                nv, chart, nchart, anchor_source, source, calib, lift, uv,
                out_topology) != 0 ||
            wwl_assemble_bridge_gauges(
                verts, nv, faces, source_nf, nf, chart, nchart, ncomponent,
                0, 0, 0, uv, out_bridge, NULL) != 0)
            goto fail;
        free(theta); free(lift); free(queue); free(seen); free(scratch);
        free(component_count); free(component_offset);
        free(component_cursor); free(component_turn);
        free(target_u); free(target_v);
        *out_sign = sign;
        *out_offset = offset;
        *out_cycle_error = maximum_cycle_error;
        return 0;
    }

    for (size_t v = 0; v < nv; v++) {
        double phi = lift[v];
        target_u[v] = calib->spiral_a * phi +
            calib->spiral_b * phi * phi / (2.0 * WWL_2PI);
        target_v[v] =
            (double)verts[v * 3] * calib->axis_dir[0] +
            (double)verts[v * 3 + 1] * calib->axis_dir[1] +
            (double)verts[v * 3 + 2] * calib->axis_dir[2];
    }

    if (topology_mode == 2) {
        for (size_t v = 0; v < nv; v++) {
            uv[v * 2] = target_u[v];
            uv[v * 2 + 1] = target_v[v];
        }
        free(theta); free(lift); free(queue); free(seen); free(scratch);
        free(component_count); free(component_offset);
        free(component_cursor); free(component_turn);
        free(target_u); free(target_v);
        *out_sign = sign;
        *out_offset = offset;
        *out_cycle_error = maximum_cycle_error;
        return 0;
    }

    /* The independently safe winding charts do not share an orientation
     * gauge. Fit a full rigid isometry, not merely a translation, to the
     * topology-lifted developable scaffold. */
    double *target = (double *)malloc((nv ? nv : 1) * 2 * sizeof(*target));
    if (target == NULL) goto fail;
    for (size_t v = 0; v < nv; v++) {
        target[v * 2] = target_u[v];
        target[v * 2 + 1] = target_v[v];
    }
    int gauge_rc = wwl_apply_target_chart_isometries(
        nv, chart, nchart, anchor_source, target, uv, out_gauge);
    free(target);
    if (gauge_rc != 0) goto fail;

    free(theta); free(lift); free(queue); free(seen); free(scratch);
    free(component_count); free(component_offset);
    free(component_cursor); free(component_turn);
    free(target_u); free(target_v);
    *out_sign = sign;
    *out_offset = offset;
    *out_cycle_error = maximum_cycle_error;
    return 0;
fail:
    free(theta); free(lift); free(queue); free(seen); free(scratch);
    free(component_count); free(component_offset);
    free(component_cursor); free(component_turn);
    free(target_u); free(target_v);
    return -1;
}

/* A source-face prefix may contain a small chart made by certified hole/seam
 * construction, or a placed chart whose sampled winding field is rank one.
 * Neither supplies a usable local 2-D map.  Seed only charts with no source
 * anchor or no strict orientation majority from the global topology-winding /
 * axial scaffold.  Every already-valid winding chart is left bit-for-bit
 * alone.  A normal atlas export still requires a strong local orientation
 * majority.  The fitted-ribbon path uses this target only as a winding/bridge
 * gauge before rebuilding geometry on physical slices, so mixed projected
 * orientation is admissible there and is audited by the final grid gates. */
static int wwl_seed_invalid_chart_targets(
    size_t nv, const int32_t *faces, size_t source_nf,
    const int32_t *chart, size_t nchart, const int32_t *anchor_source,
    const double *target_uv, double *uv, int gauge_only,
    size_t *out_seeded)
{
    typedef struct WwlSeedChart {
        size_t anchors;
        size_t raw_positive, raw_negative;
        size_t target_positive, target_negative, target_degenerate;
        uint8_t seed;
    } WwlSeedChart;
    WwlSeedChart *state = NULL;
    size_t seeded = 0;
    size_t rejected = 0;
    size_t strict_seeded = 0, majority_seeded = 0, mixed_seeded = 0;
    if (faces == NULL || chart == NULL || anchor_source == NULL ||
        target_uv == NULL || uv == NULL || nchart == 0)
        return -1;
    state = (WwlSeedChart *)calloc(nchart, sizeof(*state));
    if (state == NULL) return -1;
    for (size_t v = 0; v < nv; v++) {
        int32_t c = chart[v];
        if (c < 0 || (size_t)c >= nchart) { free(state); return -1; }
        if (anchor_source[v] >= 0) state[(size_t)c].anchors++;
    }
    /* Accumulate both orientation audits in one face pass.  The previous
     * chart-by-face nested loops were O(nchart * source_nf). */
    for (size_t f = 0; f < source_nf; f++) {
        int32_t a = faces[f * 3], b = faces[f * 3 + 1];
        int32_t d = faces[f * 3 + 2];
        if (a < 0 || b < 0 || d < 0 || (size_t)a >= nv ||
            (size_t)b >= nv || (size_t)d >= nv) {
            free(state); return -1;
        }
        int32_t c = chart[(size_t)a];
        if (c < 0 || (size_t)c >= nchart || chart[(size_t)b] != c ||
            chart[(size_t)d] != c) {
            free(state); return -1;
        }
        WwlSeedChart *s = &state[(size_t)c];
        double raw_area = wwl_area2_points(&uv[(size_t)a * 2],
                                           &uv[(size_t)b * 2],
                                           &uv[(size_t)d * 2]);
        double target_area = wwl_area2_points(&target_uv[(size_t)a * 2],
                                              &target_uv[(size_t)b * 2],
                                              &target_uv[(size_t)d * 2]);
        if (raw_area > 1.0e-12) s->raw_positive++;
        else if (raw_area < -1.0e-12) s->raw_negative++;
        if (target_area > 1.0e-12) s->target_positive++;
        else if (target_area < -1.0e-12) s->target_negative++;
        else s->target_degenerate++;
    }
    for (size_t c = 0; c < nchart; c++) {
        WwlSeedChart *s = &state[c];
        size_t total = s->target_positive + s->target_negative;
        size_t minority = s->target_positive < s->target_negative
                        ? s->target_positive : s->target_negative;
        int strict = (s->target_positive == 0) != (s->target_negative == 0);
        int strong_majority = total >= 10 && minority <= total / 10;
        if (s->anchors != 0 && s->raw_positive != s->raw_negative) continue;
        if ((gauge_only && total == 0) ||
            (!gauge_only &&
             (s->target_degenerate != 0 || (!strict && !strong_majority)))) {
            fprintf(stderr,
                "welded_winding_lift: topology scaffold rejected for "
                "anchorless chart %zu (+%zu -%zu zero=%zu)\n",
                c, s->target_positive, s->target_negative,
                s->target_degenerate);
            rejected++;
            continue;
        }
        s->seed = 1;
        if (strict) {
            strict_seeded++;
        } else if (strong_majority) {
            majority_seeded++;
        } else {
            mixed_seeded++;
        }
        if (!strict)
            fprintf(stderr,
                "welded_winding_lift: seeded invalid chart %zu from %s "
                "topology scaffold (anchors=%zu raw=+%zu/-%zu "
                "target=+%zu/-%zu/zero%zu)\n",
                c, strong_majority ? "90%-majority" : "gauge-only mixed",
                s->anchors, s->raw_positive, s->raw_negative,
                s->target_positive, s->target_negative,
                s->target_degenerate);
        seeded++;
    }
    if (rejected != 0) {
        fprintf(stderr,
            "welded_winding_lift: %zu invalid chart scaffold(s) rejected\n",
            rejected);
        free(state);
        return -1;
    }
    for (size_t v = 0; v < nv; v++)
        if (state[(size_t)chart[v]].seed) {
            uv[v * 2] = target_uv[v * 2];
            uv[v * 2 + 1] = target_uv[v * 2 + 1];
        }
    if (seeded != 0)
        fprintf(stderr,
            "welded_winding_lift: invalid-chart scaffold classes: "
            "strict=%zu 90%%-majority=%zu gauge-only-mixed=%zu\n",
            strict_seeded, majority_seeded, mixed_seeded);
    free(state);
    if (out_seeded != NULL) *out_seeded = seeded;
    return 0;
}

static int wwl_pack_components(
    size_t nv, double *uv, const int32_t *component,
    size_t ncomponent, const size_t *face_count, double gap)
{
    WwlPackRecord *record = (WwlPackRecord *)malloc(
        (ncomponent ? ncomponent : 1) * sizeof(*record));
    double *shift_u = (double *)malloc(
        (ncomponent ? ncomponent : 1) * sizeof(*shift_u));
    double *shift_v = (double *)malloc(
        (ncomponent ? ncomponent : 1) * sizeof(*shift_v));
    if (record == NULL || shift_u == NULL || shift_v == NULL) {
        free(record); free(shift_u); free(shift_v); return -1;
    }
    double global_vmin = DBL_MAX;
    for (size_t c = 0; c < ncomponent; c++) {
        record[c].component = (int32_t)c;
        record[c].faces = face_count[c];
        record[c].umin = record[c].vmin = DBL_MAX;
        record[c].umax = record[c].vmax = -DBL_MAX;
    }
    for (size_t v = 0; v < nv; v++) {
        int32_t c = component[v];
        if (c < 0) continue;
        WwlPackRecord *r = &record[c];
        if (uv[v * 2] < r->umin) r->umin = uv[v * 2];
        if (uv[v * 2] > r->umax) r->umax = uv[v * 2];
        if (uv[v * 2 + 1] < r->vmin) r->vmin = uv[v * 2 + 1];
        if (uv[v * 2 + 1] > r->vmax) r->vmax = uv[v * 2 + 1];
        if (uv[v * 2 + 1] < global_vmin) global_vmin = uv[v * 2 + 1];
    }
    if (!isfinite(global_vmin)) {
        free(record); free(shift_u); free(shift_v); return -1;
    }
    qsort(record, ncomponent, sizeof(*record), wwl_compare_pack);
    double cursor = 0.0;
    for (size_t order = 0; order < ncomponent; order++) {
        WwlPackRecord *r = &record[order];
        shift_u[r->component] = cursor - r->umin;
        /* This is a horizontal atlas pack.  V is the shared axial/world-Z
         * coordinate, so independently rebasing every component's V minimum
         * destroys their physical row alignment.  A single global rebase keeps
         * the output raster-friendly without changing any relative V value. */
        shift_v[r->component] = -global_vmin;
        cursor += (r->umax - r->umin) + gap;
    }
    for (size_t v = 0; v < nv; v++) {
        int32_t c = component[v];
        if (c < 0) continue;
        uv[v * 2] += shift_u[c];
        uv[v * 2 + 1] += shift_v[c];
    }
    free(record); free(shift_u); free(shift_v);
    return 0;
}

/* A stitched welded component already owns its intrinsic metric U.  The
 * topology scaffold is used only to recover that component's additive cover
 * gauge; averaging one constant cannot shear, scale, or rotate the component. */
static int wwl_align_component_u_to_target(
    size_t nv, double *uv, const double *target_uv,
    const int32_t *component, size_t ncomponent,
    double *out_mean_residual, double *out_maximum_residual)
{
    double *sum = NULL, *compensation = NULL, *shift = NULL;
    size_t *count = NULL;
    if (uv == NULL || target_uv == NULL || component == NULL ||
        ncomponent == 0)
        return -1;
    sum = (double *)calloc(ncomponent, sizeof(*sum));
    compensation = (double *)calloc(ncomponent, sizeof(*compensation));
    shift = (double *)malloc(ncomponent * sizeof(*shift));
    count = (size_t *)calloc(ncomponent, sizeof(*count));
    if (sum == NULL || compensation == NULL || shift == NULL || count == NULL)
        goto fail;
    for (size_t v = 0; v < nv; v++) {
        int32_t c = component[v];
        if (c < 0 || (size_t)c >= ncomponent ||
            !isfinite(uv[v * 2]) || !isfinite(target_uv[v * 2]))
            goto fail;
        double value = target_uv[v * 2] - uv[v * 2];
        double corrected = value - compensation[c];
        double next = sum[c] + corrected;
        compensation[c] = (next - sum[c]) - corrected;
        sum[c] = next;
        count[c]++;
    }
    for (size_t c = 0; c < ncomponent; c++) {
        if (count[c] == 0) goto fail;
        shift[c] = sum[c] / (double)count[c];
    }
    double residual_sum = 0.0, residual_maximum = 0.0;
    for (size_t v = 0; v < nv; v++) {
        int32_t c = component[v];
        uv[v * 2] += shift[c];
        double residual = fabs(uv[v * 2] - target_uv[v * 2]);
        residual_sum += residual;
        if (residual > residual_maximum) residual_maximum = residual;
    }
    if (out_mean_residual != NULL)
        *out_mean_residual = nv ? residual_sum / (double)nv : 0.0;
    if (out_maximum_residual != NULL)
        *out_maximum_residual = residual_maximum;
    free(sum); free(compensation); free(shift); free(count);
    return 0;
fail:
    free(sum); free(compensation); free(shift); free(count);
    return -1;
}

/* Preserve cover placement wherever component U intervals overlap, and remove
 * only the empty, therefore unobservable, gaps between interval islands.  This
 * is the direct-mesh analogue of Ribbon's metric-atlas island pack. */
static int wwl_pack_component_intervals(
    size_t nv, double *uv, const int32_t *component,
    size_t ncomponent, const size_t *face_count, double gap,
    size_t *out_islands, double *out_cover_span, double *out_atlas_span)
{
    WwlPackRecord *record = NULL;
    double *shift_u = NULL;
    double global_vmin = DBL_MAX, cover_umin = DBL_MAX, cover_umax = -DBL_MAX;
    if (uv == NULL || component == NULL || face_count == NULL ||
        ncomponent == 0)
        return -1;
    record = (WwlPackRecord *)malloc(ncomponent * sizeof(*record));
    shift_u = (double *)malloc(ncomponent * sizeof(*shift_u));
    if (record == NULL || shift_u == NULL) goto fail;
    for (size_t c = 0; c < ncomponent; c++) {
        record[c].component = (int32_t)c;
        record[c].faces = face_count[c];
        record[c].umin = record[c].vmin = DBL_MAX;
        record[c].umax = record[c].vmax = -DBL_MAX;
    }
    for (size_t v = 0; v < nv; v++) {
        int32_t c = component[v];
        if (c < 0 || (size_t)c >= ncomponent) goto fail;
        WwlPackRecord *r = &record[c];
        double u = uv[v * 2], w = uv[v * 2 + 1];
        if (!isfinite(u) || !isfinite(w)) goto fail;
        if (u < r->umin) r->umin = u;
        if (u > r->umax) r->umax = u;
        if (w < r->vmin) r->vmin = w;
        if (w > r->vmax) r->vmax = w;
        if (u < cover_umin) cover_umin = u;
        if (u > cover_umax) cover_umax = u;
        if (w < global_vmin) global_vmin = w;
    }
    for (size_t c = 0; c < ncomponent; c++)
        if (!(record[c].umin <= record[c].umax)) goto fail;
    qsort(record, ncomponent, sizeof(*record), wwl_compare_pack_umin);

    size_t islands = 0;
    double cursor = 0.0, island_end = -DBL_MAX, island_shift = 0.0;
    for (size_t order = 0; order < ncomponent; order++) {
        WwlPackRecord *r = &record[order];
        if (order == 0 || r->umin > island_end + 1.0e-9) {
            if (order != 0) cursor += island_end + island_shift + gap - cursor;
            island_shift = cursor - r->umin;
            island_end = r->umax;
            islands++;
        } else if (r->umax > island_end) {
            island_end = r->umax;
        }
        shift_u[r->component] = island_shift;
    }
    double atlas_span = island_end + island_shift;
    for (size_t v = 0; v < nv; v++) {
        int32_t c = component[v];
        uv[v * 2] += shift_u[c];
        uv[v * 2 + 1] -= global_vmin;
    }
    if (out_islands != NULL) *out_islands = islands;
    if (out_cover_span != NULL) *out_cover_span = cover_umax - cover_umin;
    if (out_atlas_span != NULL) *out_atlas_span = atlas_span;
    free(record); free(shift_u);
    return 0;
fail:
    free(record); free(shift_u);
    return -1;
}

static int wwl_load_raw_uv(const char *placed_dir, PieceSet *source)
{
    for (size_t cube = 0; cube < source->n_cubes; cube++) {
        char path[2048];
        size_t first = source->cube_voff[cube];
        size_t count = source->cube_voff[cube + 1] - first;
        if (snprintf(path, sizeof(path), "%s/%s_uvphi_raw.f32",
                     placed_dir, source->ids[cube]) < 0)
            return -1;
        FILE *file = fopen(path, "rb");
        if (file == NULL) return -1;
        float *triples = (float *)malloc(
            (count ? count : 1) * 3 * sizeof(*triples));
        if (triples == NULL ||
            fread(triples, 3 * sizeof(*triples), count, file) != count ||
            fgetc(file) != EOF) {
            free(triples);
            fclose(file);
            return -1;
        }
        fclose(file);
        for (size_t i = 0; i < count; i++) {
            source->uv[(first + i) * 2] = triples[i * 3];
            source->uv[(first + i) * 2 + 1] = triples[i * 3 + 1];
            source->phi[first + i] = triples[i * 3 + 2];
        }
        free(triples);
    }
    return 0;
}

static int wwl_chart_layout_map(
    const WwlChartLayout *layout,const PieceSet *source,size_t welded_nv,
    const int32_t *anchor_source,const int32_t *chart,size_t nchart,
    WwlChartLayoutCorrection **out_correction,size_t *out_mapped,
    double *out_minimum,double *out_maximum)
{
    WwlChartLayoutCorrection *correction=NULL;uint8_t *seen=NULL;
    int32_t *source_chart=NULL;size_t source_nv=0,last_cube=(size_t)-1;
    size_t mapped=0;double minimum=DBL_MAX,maximum=-DBL_MAX;
    if(!layout||!source||!anchor_source||!chart||!out_correction||
       layout->nnode==0||nchart==0)return -1;
    *out_correction=NULL;
    correction=(WwlChartLayoutCorrection*)malloc(
        nchart*sizeof(*correction));
    seen=(uint8_t*)calloc(nchart,1);
    if(!correction||!seen)goto fail;
    for(size_t c=0;c<nchart;c++){
        correction[c].version=0;
        correction[c].shift_u=NAN;
        correction[c].shift_u_per_v=NAN;
        correction[c].v_reference=NAN;
        correction[c].v_step=NAN;
        correction[c].winding_turn=0;
        correction[c].relation_component=-1;
        correction[c].layout_node=-1;
        for(size_t k=0;k<WWL_CHART_LAYOUT_V3_KNOTS;k++)
            correction[c].shift_v[k]=NAN;
    }
    /* Invert the anchor map ONCE.  Scanning every welded vertex per layout
       node is quadratic: 97,077 nodes against 124M vertices on the
       21x21x21 is about 1.2e13 probes and tens of terabytes of memory
       traffic, and it simply never finished.  One pass over the welded
       vertices builds the reverse map, after which each node is a lookup.
       -1 means no welded vertex is anchored there; -2 means several are
       and they disagree about the chart -- both are failures, but only
       when a node actually asks, which is what the per-node scan did. */
    source_nv=source->cube_voff[source->n_cubes];
    source_chart=(int32_t*)malloc((source_nv?source_nv:1)*
                                  sizeof(*source_chart));
    if(!source_chart)goto fail;
    for(size_t s=0;s<source_nv;s++)source_chart[s]=-1;
    for(size_t v=0;v<welded_nv;v++){
        int32_t s=anchor_source[v],c;
        if(s<0)continue;
        if((size_t)s>=source_nv)goto fail;
        c=chart[v];
        if(c<0||(size_t)c>=nchart)source_chart[s]=-2;
        else if(source_chart[s]==-1)source_chart[s]=c;
        else if(source_chart[s]!=c)source_chart[s]=-2;
    }
    for(size_t i=0;i<layout->nnode;i++){
        const WwlChartLayoutNode *node=&layout->node[i];
        size_t cube=source->n_cubes,source_vertex;
        int32_t mapped_chart=-1;
        /* Charts are emitted grouped by cube, so the previous hit answers
           almost every lookup and the linear scan runs ~once per cube
           rather than once per chart. */
        if(last_cube<source->n_cubes&&
           strcmp(source->ids[last_cube],node->cube_id)==0){
            cube=last_cube;
        }else{
            for(size_t c=0;c<source->n_cubes;c++)
                if(strcmp(source->ids[c],node->cube_id)==0){cube=c;break;}
            last_cube=cube;
        }
        if(cube==source->n_cubes||
           (size_t)node->representative_vertex>=
               source->cube_voff[cube+1]-source->cube_voff[cube])
            goto fail;
        source_vertex=source->cube_voff[cube]+
                      (size_t)node->representative_vertex;
        mapped_chart=source_chart[source_vertex];
        if(mapped_chart<0||seen[mapped_chart])goto fail;
        correction[mapped_chart].version=layout->version;
        correction[mapped_chart].shift_u=node->shift_u;
        correction[mapped_chart].shift_u_per_v=node->shift_u_per_v;
        correction[mapped_chart].v_reference=node->v_reference;
        correction[mapped_chart].v_step=node->v_step;
        correction[mapped_chart].winding_turn=node->winding_turn;
        correction[mapped_chart].relation_component=
            node->relation_component;
        correction[mapped_chart].layout_node=(int32_t)i;
        memcpy(correction[mapped_chart].shift_v,node->shift_v,
               sizeof(node->shift_v));
        seen[mapped_chart]=1;mapped++;
        if(layout->version>=3){
            for(size_t k=0;k<WWL_CHART_LAYOUT_V3_KNOTS;k++){
                if(node->shift_v[k]<minimum)minimum=node->shift_v[k];
                if(node->shift_v[k]>maximum)maximum=node->shift_v[k];
            }
        }else{
            if(node->shift_u<minimum)minimum=node->shift_u;
            if(node->shift_u>maximum)maximum=node->shift_u;
        }
    }
    if(mapped!=layout->nnode)goto fail;
    /* A chart with no match at all is never a graph node, so the layout
       cannot speak for it.  Requiring a perfect 1:1 rejected the whole
       parameterization over three such charts out of 97,080.  Give them an
       IDENTITY correction -- their carried u is used exactly as it stands
       -- and an island of their own, so the reference-island check sees a
       component that spans nothing and the fit leaves them alone.  They
       are isolated by construction: having no relation, nothing else
       depends on where they sit. */
    if(mapped<nchart){
        size_t orphan=0;
        for(size_t c=0;c<nchart;c++){
            if(seen[c])continue;
            correction[c].version=layout->version;
            correction[c].shift_u=0.0;
            correction[c].shift_u_per_v=0.0;
            correction[c].v_reference=0.0;
            correction[c].v_step=1.0;
            correction[c].winding_turn=0;
            for(size_t k=0;k<WWL_CHART_LAYOUT_V3_KNOTS;k++)
                correction[c].shift_v[k]=0.0;
            correction[c].relation_component=
                (int32_t)(layout->relation_components+orphan);
            correction[c].layout_node=-1;
            orphan++;
        }
        fprintf(stderr,
            "welded_winding_lift: %zu chart(s) have no layout node "
            "(unmatched); carried u used as-is, each in its own island\n",
            orphan);
    }
    free(seen);free(source_chart);*out_correction=correction;
    if(out_mapped)*out_mapped=mapped;
    if(out_minimum)*out_minimum=minimum;
    if(out_maximum)*out_maximum=maximum;
    return 0;
fail:
    free(seen);free(source_chart);free(correction);return -1;
}

static int wwl_chart_layout_map_relations(
    const WwlChartLayout *layout,
    const WwlChartLayoutCorrection *correction,size_t nchart,
    int32_t **out_a,int32_t **out_b,double **out_support)
{
    int32_t *node_chart=NULL,*a=NULL,*b=NULL;
    double *support=NULL;
    WwlMappedChartRelation *mapped=NULL;
    if(!layout||!correction||!out_a||!out_b||!out_support||
       layout->version<4)return -1;
    *out_a=NULL;*out_b=NULL;*out_support=NULL;
    node_chart=(int32_t*)malloc(
        (layout->nnode?layout->nnode:1)*sizeof(*node_chart));
    a=(int32_t*)malloc((layout->relation_edges?layout->relation_edges:1)*
                      sizeof(*a));
    b=(int32_t*)malloc((layout->relation_edges?layout->relation_edges:1)*
                      sizeof(*b));
    support=(double*)malloc(
        (layout->relation_edges?layout->relation_edges:1)*sizeof(*support));
    mapped=(WwlMappedChartRelation*)malloc(
        (layout->relation_edges?layout->relation_edges:1)*sizeof(*mapped));
    if(!node_chart||!a||!b||!support||!mapped)goto fail;
    for(size_t i=0;i<layout->nnode;i++)node_chart[i]=-1;
    for(size_t c=0;c<nchart;c++){
        int32_t node=correction[c].layout_node;
        /* Unmatched charts carry no node; the check that matters is that
           every NODE is claimed exactly once, which the loop below makes. */
        if(node<0)continue;
        if((size_t)node>=layout->nnode||node_chart[node]>=0)goto fail;
        node_chart[node]=(int32_t)c;
    }
    for(size_t i=0;i<layout->nnode;i++)if(node_chart[i]<0)goto fail;
    for(size_t e=0;e<layout->relation_edges;e++){
        mapped[e].a=node_chart[layout->edge[e].a];
        mapped[e].b=node_chart[layout->edge[e].b];
        mapped[e].support=layout->edge[e].support;
        if(mapped[e].a<0||mapped[e].b<0||mapped[e].a==mapped[e].b)
            goto fail;
        if(mapped[e].b<mapped[e].a){
            int32_t swap=mapped[e].a;
            mapped[e].a=mapped[e].b;mapped[e].b=swap;
        }
    }
    qsort(mapped,layout->relation_edges,sizeof(*mapped),
          wwl_mapped_chart_relation_compare);
    for(size_t e=0;e<layout->relation_edges;e++){
        a[e]=mapped[e].a;b[e]=mapped[e].b;support[e]=mapped[e].support;
    }
    free(node_chart);free(mapped);
    *out_a=a;*out_b=b;*out_support=support;
    return 0;
fail:
    free(node_chart);free(a);free(b);free(support);free(mapped);return -1;
}

/* Vertex island identity is the chart's relation component; every SOURCE
 * face therefore lies in exactly one island by construction and is checked
 * strictly.  A welded mesh component or a BRIDGE face MAY legitimately span
 * two islands: that is precisely a physically-welded junction whose carried-U
 * constraint the layout deferred for winding synchronization.  Those spans
 * are counted and reported, never failed -- the chart-first fit consumes
 * source faces only, so a deferred junction can never leak into emission. */
static int wwl_chart_layout_reference_islands(
    size_t nv,const int32_t *chart,size_t nchart,
    const WwlChartLayoutCorrection *correction,
    const int32_t *mesh_component,size_t nmesh_component,
    const int32_t *faces,size_t nf,size_t source_faces,size_t nisland,
    size_t *out_spanning_components,size_t *out_spanning_bridge_faces,
    int32_t **out_vertex_island,size_t **out_island_faces)
{
    int32_t *vertex_island=NULL,*mesh_island=NULL;
    size_t *island_faces=NULL;uint8_t *island_seen=NULL,*mesh_spans=NULL;
    size_t spanning_components=0,spanning_bridge_faces=0;
    if(!chart||!correction||!mesh_component||!faces||
       !out_vertex_island||!out_island_faces||nchart==0||
       nmesh_component==0||nisland==0||source_faces>nf)return -1;
    *out_vertex_island=NULL;*out_island_faces=NULL;
    if(out_spanning_components)*out_spanning_components=0;
    if(out_spanning_bridge_faces)*out_spanning_bridge_faces=0;
    vertex_island=(int32_t*)malloc((nv?nv:1)*sizeof(*vertex_island));
    mesh_island=(int32_t*)malloc(
        nmesh_component*sizeof(*mesh_island));
    island_faces=(size_t*)calloc(nisland,sizeof(*island_faces));
    island_seen=(uint8_t*)calloc(nisland,1);
    mesh_spans=(uint8_t*)calloc(nmesh_component,1);
    if(!vertex_island||!mesh_island||!island_faces||!island_seen||
       !mesh_spans)goto fail;
    for(size_t c=0;c<nmesh_component;c++)mesh_island[c]=-1;
    for(size_t v=0;v<nv;v++){
        int32_t source_chart=chart[v],mc=mesh_component[v],island;
        if(source_chart<0||(size_t)source_chart>=nchart||
           mc<0||(size_t)mc>=nmesh_component)goto fail;
        island=correction[source_chart].relation_component;
        if(island<0||(size_t)island>=nisland)goto fail;
        vertex_island[v]=island;island_seen[island]=1;
        if(mesh_island[mc]<0)mesh_island[mc]=island;
        else if(mesh_island[mc]!=island)mesh_spans[mc]=1;
    }
    for(size_t c=0;c<nmesh_component;c++)
        spanning_components+=mesh_spans[c]?1u:0u;
    for(size_t f=0;f<nf;f++){
        int32_t a=faces[f*3],b=faces[f*3+1],c=faces[f*3+2];
        int same=0;
        if(a<0||b<0||c<0||(size_t)a>=nv||(size_t)b>=nv||
           (size_t)c>=nv)goto fail;
        same=vertex_island[a]==vertex_island[b]&&
             vertex_island[a]==vertex_island[c];
        if(!same){
            if(f<source_faces)goto fail; /* charts never span islands */
            spanning_bridge_faces++;
            continue;
        }
        island_faces[vertex_island[a]]++;
    }
    for(size_t island=0;island<nisland;island++)
        if(!island_seen[island]||island_faces[island]==0)goto fail;
    free(mesh_island);free(island_seen);free(mesh_spans);
    if(out_spanning_components)
        *out_spanning_components=spanning_components;
    if(out_spanning_bridge_faces)
        *out_spanning_bridge_faces=spanning_bridge_faces;
    *out_vertex_island=vertex_island;
    *out_island_faces=island_faces;
    return 0;
fail:
    free(vertex_island);free(mesh_island);free(island_faces);
    free(island_seen);free(mesh_spans);return -1;
}

static int wwl_chart_layout_apply(size_t nv,const float *verts,double *uv,
                                  const int32_t *chart,size_t nchart,
                                  const WwlChartLayoutCorrection *correction)
{
    if(!verts||!uv||!chart||!correction)return -1;
    for(size_t v=0;v<nv;v++){
        int32_t c=chart[v];
        double shift;
        if(c<0||(size_t)c>=nchart)return -1;
        if(correction[c].version>=3){
            if(!isfinite(correction[c].v_reference)||
               !isfinite(correction[c].v_step)||
               correction[c].v_step<=0.0)return -1;
            double t=((double)verts[v*3]-correction[c].v_reference)/
                     correction[c].v_step;
            if(t<=0.0)shift=correction[c].shift_v[0];
            else if(t>=(double)(WWL_CHART_LAYOUT_V3_KNOTS-1))
                shift=correction[c].shift_v[WWL_CHART_LAYOUT_V3_KNOTS-1];
            else{
                size_t k=(size_t)floor(t);double f=t-(double)k;
                shift=(1.0-f)*correction[c].shift_v[k]+
                      f*correction[c].shift_v[k+1];
            }
        }else{
            if((correction[c].version!=1&&correction[c].version!=2)||
               !isfinite(correction[c].shift_u)||
               !isfinite(correction[c].shift_u_per_v)||
               !isfinite(correction[c].v_reference))return -1;
            shift=correction[c].shift_u+
                  correction[c].shift_u_per_v*
                  ((double)verts[v*3]-correction[c].v_reference);
        }
        if(!isfinite(shift))return -1;
        uv[v*2]+=shift;
    }
    return 0;
}

static int wwl_write_placed_chart_diagnostic(
    const char *output_path, const PieceSet *source, double gap)
{
    int32_t *chart = NULL;
    size_t nchart = 0, *face_count = NULL;
    double *uv = NULL;
    int result = -1;
    if (wwl_label_components(source->nv, source->faces, source->nf,
                             &chart, &nchart, &face_count) != 0)
        goto cleanup;
    uv = (double *)malloc(source->nv * 2 * sizeof(*uv));
    if (uv == NULL) goto cleanup;
    for (size_t v = 0; v < source->nv; v++) {
        uv[v * 2] = source->uv[v * 2];
        uv[v * 2 + 1] = source->uv[v * 2 + 1];
    }
    if (wwl_pack_components(source->nv, uv, chart, nchart,
                            face_count, gap) != 0)
        goto cleanup;
    size_t positive = 0, negative = 0, degenerate = 0;
    wwl_face_orientation_stats(uv, source->faces, source->nf,
                               &positive, &negative, &degenerate);
    fprintf(stderr,
        "welded_winding_lift: placed diagnostic %zu charts; "
        "UV signs +%zu -%zu zero=%zu\n",
        nchart, positive, negative, degenerate);
    if (ves_ensure_parent_dir(output_path) != 0 ||
        ObjIO_write_uv_double(output_path, source->verts, source->nv,
                              source->faces, source->nf, uv) != 0)
        goto cleanup;
    result = 0;
cleanup:
    free(chart);
    free(face_count);
    free(uv);
    return result;
}

static void wwl_face_orientation_stats(
    const double *uv, const int32_t *faces, size_t nf,
    size_t *positive, size_t *negative, size_t *degenerate)
{
    *positive = *negative = *degenerate = 0;
    for (size_t f = 0; f < nf; f++) {
        size_t a = (size_t)faces[f * 3];
        size_t b = (size_t)faces[f * 3 + 1];
        size_t c = (size_t)faces[f * 3 + 2];
        double area = (uv[b * 2] - uv[a * 2]) *
                      (uv[c * 2 + 1] - uv[a * 2 + 1]) -
                      (uv[b * 2 + 1] - uv[a * 2 + 1]) *
                      (uv[c * 2] - uv[a * 2]);
        if (area > 1.0e-12) (*positive)++;
        else if (area < -1.0e-12) (*negative)++;
        else (*degenerate)++;
    }
}

static double wwl_area2_points(const double a[2], const double b[2],
                               const double c[2])
{
    return (b[0] - a[0]) * (c[1] - a[1]) -
           (b[1] - a[1]) * (c[0] - a[0]);
}

static double wwl_distance3(const float *a, const float *b)
{
    double d0 = (double)b[0] - a[0];
    double d1 = (double)b[1] - a[1];
    double d2 = (double)b[2] - a[2];
    return sqrt(d0 * d0 + d1 * d1 + d2 * d2);
}

/* The same branch-cut-free local winding invariant used by ChartZipper's
 * hard weld gate.  The production scroll axis is stored in (z,y,x), with the
 * polar phase atan2(y,x); retaining that exact convention makes this audit
 * directly comparable to the transaction that admitted the bridge. */
static double wwl_physical_winding_delta(
    const ScaffoldCalib *calib, const float *a, const float *b)
{
    if (calib == NULL || !(calib->pitch > 0.0)) return NAN;
    double ay = (double)a[1] - calib->axis_point[1];
    double ax = (double)a[2] - calib->axis_point[2];
    double by = (double)b[1] - calib->axis_point[1];
    double bx = (double)b[2] - calib->axis_point[2];
    double ar = hypot(ay, ax), br = hypot(by, bx);
    if (!(ar > 1.0e-12) || !(br > 1.0e-12)) return NAN;
    double dtheta = wwl_wrap_pi(atan2(by, bx) - atan2(ay, ax));
    return fabs((br - ar) / calib->pitch - dtheta / WWL_2PI);
}

static int32_t wwl_source_cube_index(const PieceSet *source, int32_t vertex)
{
    if (source == NULL || vertex < 0 || (size_t)vertex >= source->nv)
        return -1;
    size_t lo = 0, hi = source->n_cubes;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if ((size_t)vertex < source->cube_voff[mid]) hi = mid;
        else if ((size_t)vertex >= source->cube_voff[mid + 1]) lo = mid + 1;
        else return (int32_t)mid;
    }
    return -1;
}

static int wwl_compare_bridge_observation(const void *left, const void *right)
{
    const WwlBridgeObservation *a = (const WwlBridgeObservation *)left;
    const WwlBridgeObservation *b = (const WwlBridgeObservation *)right;
    if (a->lo != b->lo) return a->lo < b->lo ? -1 : 1;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    return 0;
}

/* Fit target = R*source+t for one directed side of a chart pair.  R is a
 * proper 2-D rotation because handedness has already been normalized. */
static int wwl_fit_bridge_isometry(
    const WwlBridgeObservation *observation, size_t first, size_t last,
    int32_t base, int allow_rotation, WwlBridgeGaugeEdge *out,
    double *scratch)
{
    size_t count = 0;
    double sx = 0.0, sy = 0.0, tx = 0.0, ty = 0.0;
    for (size_t i = first; i < last; i++) {
        if (observation[i].base != base) continue;
        sx += observation[i].source_u;
        sy += observation[i].source_v;
        tx += observation[i].target_u;
        ty += observation[i].target_v;
        count++;
    }
    if (count == 0) return -1;
    sx /= (double)count; sy /= (double)count;
    tx /= (double)count; ty /= (double)count;
    double dot = 0.0, cross = 0.0, source_spread = 0.0;
    for (size_t i = first; i < last; i++) {
        if (observation[i].base != base) continue;
        double x = observation[i].source_u - sx;
        double y = observation[i].source_v - sy;
        double a = observation[i].target_u - tx;
        double b = observation[i].target_v - ty;
        dot += x * a + y * b;
        cross += x * b - y * a;
        source_spread += x * x + y * y;
    }
    double cosine = 1.0, sine = 0.0;
    if (allow_rotation && count >= 2 && source_spread > 1.0e-12 &&
        hypot(dot, cross) > 1.0e-12) {
        double norm = hypot(dot, cross);
        cosine = dot / norm;
        sine = cross / norm;
    }
    double translate_u, translate_v;
    if (allow_rotation) {
        translate_u = tx - (cosine * sx - sine * sy);
        translate_v = ty - (sine * sx + cosine * sy);
    } else {
        size_t at = 0;
        for (size_t i = first; i < last; i++)
            if (observation[i].base == base)
                scratch[at++] = observation[i].target_u -
                                observation[i].source_u;
        translate_u = wwl_median(scratch, at);
        at = 0;
        for (size_t i = first; i < last; i++)
            if (observation[i].base == base)
                scratch[at++] = observation[i].target_v -
                                observation[i].source_v;
        translate_v = wwl_median(scratch, at);
    }
    size_t at = 0;
    double maximum = 0.0;
    for (size_t i = first; i < last; i++) {
        if (observation[i].base != base) continue;
        double px = cosine * observation[i].source_u -
                    sine * observation[i].source_v + translate_u;
        double py = sine * observation[i].source_u +
                    cosine * observation[i].source_v + translate_v;
        double residual = hypot(px - observation[i].target_u,
                                py - observation[i].target_v);
        scratch[at++] = residual;
        if (residual > maximum) maximum = residual;
    }
    out->lo = observation[first].lo;
    out->hi = observation[first].hi;
    out->base = base;
    out->moving = base == out->lo ? out->hi : out->lo;
    out->cosine = cosine;
    out->sine = sine;
    out->translate_u = translate_u;
    out->translate_v = translate_v;
    out->observations = count;
    out->median_residual = wwl_median(scratch, count);
    out->maximum_residual = maximum;
    return 0;
}

/* Every edge of the certified chart forest is a topological generator: remove
 * it and exactly one welded sheet splits in two.  Report that leverage beside
 * the bridge's independent geometric evidence so a suspect join can be
 * corrected as one whole chart transaction, never as a scatter of faces or
 * duplicated vertices. */
static int wwl_write_bridge_edge_audit(
    const WwlBridgeAudit *audit, size_t nv, const int32_t *chart,
    size_t nchart, const WwlBridgeGaugeEdge *edge, size_t nedge,
    const size_t *offset, const int32_t *adj_edge)
{
    int result = -1;
    FILE *file = NULL;
    int32_t *chart_component = NULL, *chart_cube = NULL, *queue = NULL;
    size_t *component_charts = NULL, *component_source_faces = NULL;
    uint8_t *seen = NULL;
    if (audit == NULL || audit->path == NULL || audit->chart_faces == NULL ||
        audit->component == NULL || audit->component_faces == NULL ||
        audit->ncomponent == 0 || audit->anchor_source == NULL ||
        audit->source == NULL)
        return -1;

    chart_component = (int32_t *)malloc(nchart * sizeof(*chart_component));
    chart_cube = (int32_t *)malloc(nchart * sizeof(*chart_cube));
    queue = (int32_t *)malloc(nchart * sizeof(*queue));
    component_charts = (size_t *)calloc(
        audit->ncomponent, sizeof(*component_charts));
    component_source_faces = (size_t *)calloc(
        audit->ncomponent, sizeof(*component_source_faces));
    seen = (uint8_t *)malloc(nchart);
    if (chart_component == NULL || chart_cube == NULL || queue == NULL ||
        component_charts == NULL || component_source_faces == NULL ||
        seen == NULL)
        goto cleanup;
    for (size_t c = 0; c < nchart; c++) {
        chart_component[c] = -1;
        chart_cube[c] = -1;
    }
    for (size_t v = 0; v < nv; v++) {
        int32_t c = chart[v], component = audit->component[v];
        if (c < 0 || (size_t)c >= nchart || component < 0 ||
            (size_t)component >= audit->ncomponent)
            goto cleanup;
        if (chart_component[c] < 0) chart_component[c] = component;
        else if (chart_component[c] != component) goto cleanup;
        if (audit->anchor_source[v] >= 0) {
            int32_t cube = wwl_source_cube_index(
                audit->source, audit->anchor_source[v]);
            if (cube < 0) goto cleanup;
            if (chart_cube[c] == -1) chart_cube[c] = cube;
            else if (chart_cube[c] >= 0 && chart_cube[c] != cube)
                chart_cube[c] = -2;
        }
    }
    for (size_t c = 0; c < nchart; c++) {
        int32_t component = chart_component[c];
        if (component < 0 || (size_t)component >= audit->ncomponent)
            goto cleanup;
        component_charts[component]++;
        component_source_faces[component] += audit->chart_faces[c];
    }

    if (ves_ensure_parent_dir(audit->path) != 0) goto cleanup;
    file = fopen(audit->path, "wb");
    if (file == NULL) goto cleanup;
    fprintf(file,
        "edge,component,chart_lo,chart_hi,cube_lo,cube_hi,"
        "component_charts,component_source_faces,component_final_faces,"
        "chart_lo_faces,chart_hi_faces,side_lo_charts,side_lo_source_faces,"
        "side_hi_charts,side_hi_source_faces,small_side_charts,"
        "small_side_source_faces,pair_observations,lo_observations,"
        "hi_observations,fit_base,fit_moving,rotation_degrees,"
        "rigid_residual_median,rigid_residual_maximum,translation_u,"
        "translation_v,translation_residual_median,"
        "translation_residual_maximum,cross_edge_median,"
        "cross_edge_maximum,physical_winding_median_turns,"
        "physical_winding_maximum_turns,raw_phi_observations,"
        "raw_phi_median_turns,raw_phi_maximum_turns,"
        "registered_phi_observations,registered_phi_median_turns,"
        "registered_phi_maximum_turns,incompatibility_score\n");

    for (size_t e = 0; e < nedge; e++) {
        int32_t component = chart_component[edge[e].lo];
        if (component < 0 || chart_component[edge[e].hi] != component)
            goto cleanup;
        memset(seen, 0, nchart);
        size_t head = 0, tail = 0, side_lo_charts = 0;
        size_t side_lo_faces = 0;
        seen[edge[e].lo] = 1;
        queue[tail++] = edge[e].lo;
        while (head < tail) {
            int32_t here = queue[head++];
            side_lo_charts++;
            side_lo_faces += audit->chart_faces[here];
            for (size_t at = offset[here]; at < offset[here + 1]; at++) {
                int32_t edge_id = adj_edge[at];
                if (edge_id < 0 || (size_t)edge_id >= nedge) goto cleanup;
                if ((size_t)edge_id == e) continue;
                int32_t other = edge[edge_id].lo == here
                    ? edge[edge_id].hi : edge[edge_id].lo;
                if (!seen[other]) {
                    seen[other] = 1;
                    queue[tail++] = other;
                }
            }
        }
        if (seen[edge[e].hi]) goto cleanup;
        size_t total_charts = component_charts[component];
        size_t total_faces = component_source_faces[component];
        if (side_lo_charts >= total_charts || side_lo_faces > total_faces)
            goto cleanup;
        size_t side_hi_charts = total_charts - side_lo_charts;
        size_t side_hi_faces = total_faces - side_lo_faces;
        size_t small_charts = side_lo_charts < side_hi_charts
            ? side_lo_charts : side_hi_charts;
        size_t small_faces = side_lo_faces < side_hi_faces
            ? side_lo_faces : side_hi_faces;
        const char *cube_lo = chart_cube[edge[e].lo] >= 0
            ? audit->source->ids[chart_cube[edge[e].lo]]
            : (chart_cube[edge[e].lo] == -2 ? "MULTI" : "UNANCHORED");
        const char *cube_hi = chart_cube[edge[e].hi] >= 0
            ? audit->source->ids[chart_cube[edge[e].hi]]
            : (chart_cube[edge[e].hi] == -2 ? "MULTI" : "UNANCHORED");
        double score = edge[e].maximum_residual /
            fmax(1.0, edge[e].cross_edge_median);
        if (isfinite(edge[e].physical_winding_maximum))
            score = fmax(score, edge[e].physical_winding_maximum /
                                SEAM_WIND_HARD_TOL_DEFAULT_TURNS);
        if (isfinite(edge[e].registered_phi_maximum_turns))
            score = fmax(score, edge[e].registered_phi_maximum_turns /
                                SEAM_WIND_HARD_TOL_DEFAULT_TURNS);

        fprintf(file,
            "%zu,%d,%d,%d,%s,%s,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,"
            "%zu,%zu,%zu,%zu,%zu,%d,%d",
            e, component, edge[e].lo, edge[e].hi, cube_lo, cube_hi,
            total_charts, total_faces, audit->component_faces[component],
            audit->chart_faces[edge[e].lo], audit->chart_faces[edge[e].hi],
            side_lo_charts, side_lo_faces, side_hi_charts, side_hi_faces,
            small_charts, small_faces, edge[e].pair_observations,
            edge[e].lo_observations, edge[e].hi_observations,
            edge[e].base, edge[e].moving);
        fprintf(file,
            ",%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
            "%.17g,%.17g,%zu,%.17g,%.17g,%zu,%.17g,%.17g,%.17g\n",
            atan2(edge[e].sine, edge[e].cosine) * 180.0 / WWL_PI,
            edge[e].median_residual, edge[e].maximum_residual,
            edge[e].translation_u, edge[e].translation_v,
            edge[e].translation_median_residual,
            edge[e].translation_maximum_residual,
            edge[e].cross_edge_median, edge[e].cross_edge_maximum,
            edge[e].physical_winding_median,
            edge[e].physical_winding_maximum,
            edge[e].raw_phi_observations, edge[e].raw_phi_median_turns,
            edge[e].raw_phi_maximum_turns,
            edge[e].registered_phi_observations,
            edge[e].registered_phi_median_turns,
            edge[e].registered_phi_maximum_turns, score);
    }
    if (ferror(file) || fclose(file) != 0) {
        file = NULL;
        goto cleanup;
    }
    file = NULL;
    result = 0;

cleanup:
    if (file != NULL) fclose(file);
    free(chart_component); free(chart_cube); free(queue);
    free(component_charts); free(component_source_faces); free(seen);
    return result;
}

/* Recover one rigid translation per source chart directly from the accepted
 * bridge forest.  A bridge triangle with two vertices on chart A and one on
 * chart B is developed isometrically from A's already-safe UV edge.  This
 * gives an observation of t_B-t_A without changing either chart's intrinsic
 * parameterization.  The robust per-pair medians are then integrated over the
 * chart forest.  Because grid_weld certifies that forest, no least-squares
 * compromise or topology-changing chart split is required. */
static int wwl_assemble_bridge_gauges(
    const float *verts, size_t nv, const int32_t *faces,
    size_t source_nf, size_t nf, const int32_t *chart, size_t nchart,
    size_t expected_components, int allow_rotation, int apply_translation,
    int axial_only, double *uv,
    WwlBridgeGaugeStats *stats, const WwlBridgeAudit *audit)
{
    int result = -1, have_dsu = 0;
    int64_t *orientation_balance = NULL;
    int8_t *orientation_sign = NULL;
    WwlBridgeObservation *observation = NULL;
    WwlBridgeGaugeEdge *edge = NULL;
    double *scratch = NULL, *all_residual = NULL;
    size_t *degree = NULL, *offset = NULL, *cursor = NULL;
    int32_t *adj_edge = NULL, *queue = NULL;
    double *rotate_c = NULL, *rotate_s = NULL;
    double *shift_u = NULL, *shift_v = NULL;
    WwlDsu dsu;
    size_t nobs = 0, nedge = 0, nresidual = 0;
    memset(&dsu, 0, sizeof dsu);
    if (stats != NULL) memset(stats, 0, sizeof(*stats));
    if (verts == NULL || faces == NULL || chart == NULL || uv == NULL ||
        source_nf > nf || nchart == 0 || nchart > (size_t)INT32_MAX)
        return -1;
    if (axial_only && allow_rotation) return -1;

    orientation_balance = (int64_t *)calloc(
        nchart, sizeof(*orientation_balance));
    orientation_sign = (int8_t *)malloc(nchart * sizeof(*orientation_sign));
    observation = (WwlBridgeObservation *)malloc(
        ((nf - source_nf) ? (nf - source_nf) : 1) * sizeof(*observation));
    if (orientation_balance == NULL || orientation_sign == NULL ||
        observation == NULL)
        goto cleanup;

    for (size_t f = 0; f < source_nf; f++) {
        int32_t a = faces[f * 3], b = faces[f * 3 + 1], c = faces[f * 3 + 2];
        if (a < 0 || b < 0 || c < 0 || (size_t)a >= nv ||
            (size_t)b >= nv || (size_t)c >= nv || chart[a] < 0 ||
            chart[a] != chart[b] || chart[a] != chart[c] ||
            (size_t)chart[a] >= nchart) {
            fprintf(stderr,
                "welded_winding_lift: source-prefix chart invariant failed "
                "at face %zu vertices=%d/%d/%d charts=%d/%d/%d\n",
                f, a, b, c,
                (a >= 0 && (size_t)a < nv) ? chart[a] : -2,
                (b >= 0 && (size_t)b < nv) ? chart[b] : -2,
                (c >= 0 && (size_t)c < nv) ? chart[c] : -2);
            goto cleanup;
        }
        double area = wwl_area2_points(&uv[(size_t)a * 2],
                                       &uv[(size_t)b * 2],
                                       &uv[(size_t)c * 2]);
        if (area > 1.0e-12) orientation_balance[chart[a]]++;
        else if (area < -1.0e-12) orientation_balance[chart[a]]--;
    }
    for (size_t c = 0; c < nchart; c++) {
        if (orientation_balance[c] == 0) {
            /* In the fitted axial path U already comes from the global
             * topology scaffold and V is pinned to world Z.  A tied projected
             * face count therefore does not leave a reflection gauge to solve;
             * preserve the scaffold's U sign. */
            if (axial_only) {
                orientation_sign[c] = 1;
                continue;
            }
            size_t positive = 0, negative = 0, degenerate = 0;
            size_t vertices = 0, anchors = 0, bridge_faces = 0;
            double umin = INFINITY, umax = -INFINITY;
            double vmin = INFINITY, vmax = -INFINITY;
            for (size_t f = 0; f < source_nf; f++) {
                int32_t a = faces[f * 3];
                if (chart[a] != (int32_t)c) continue;
                int32_t b = faces[f * 3 + 1], d = faces[f * 3 + 2];
                double area = wwl_area2_points(&uv[(size_t)a * 2],
                                               &uv[(size_t)b * 2],
                                               &uv[(size_t)d * 2]);
                if (area > 1.0e-12) positive++;
                else if (area < -1.0e-12) negative++;
                else degenerate++;
            }
            for (size_t v = 0; v < nv; v++) {
                if (chart[v] != (int32_t)c) continue;
                vertices++;
                if (audit != NULL && audit->anchor_source != NULL &&
                    audit->anchor_source[v] >= 0)
                    anchors++;
                if (uv[v * 2] < umin) umin = uv[v * 2];
                if (uv[v * 2] > umax) umax = uv[v * 2];
                if (uv[v * 2 + 1] < vmin) vmin = uv[v * 2 + 1];
                if (uv[v * 2 + 1] > vmax) vmax = uv[v * 2 + 1];
            }
            for (size_t f = source_nf; f < nf; f++) {
                int touches = 0;
                for (int k = 0; k < 3; k++)
                    if (chart[faces[f * 3 + (size_t)k]] == (int32_t)c)
                        touches = 1;
                bridge_faces += (size_t)touches;
            }
            fprintf(stderr,
                "welded_winding_lift: source chart %zu has no orientation "
                "majority (+%zu -%zu zero=%zu); vertices=%zu anchors=%zu "
                "bridge_faces=%zu uv=[%.9g,%.9g]x[%.9g,%.9g]\n",
                c, positive, negative, degenerate, vertices, anchors,
                bridge_faces, umin, umax, vmin, vmax);
            goto cleanup;
        }
        orientation_sign[c] = orientation_balance[c] > 0 ? 1 : -1;
    }
    /* The independently parameterized charts may choose opposite UV
     * handedness.  In the unconstrained diagnostic modes either reflection is
     * a valid isometry.  In the production axial mode V is world Z, however,
     * so V is not a gauge: normalize handedness with (u,v)->(-u,v), and later
     * recover only an additive U gauge. */
    for (size_t c = 0; c < nchart; c++) {
        if (orientation_sign[c] < 0) {
            if (stats != NULL) stats->orientation_normalized_charts++;
            orientation_sign[c] = 1;
        }
    }
    for (size_t v = 0; v < nv; v++) {
        int32_t c = chart[v];
        if (c < 0 || (size_t)c >= nchart) goto cleanup;
        if (orientation_balance[c] < 0) {
            if (axial_only) uv[v * 2] = -uv[v * 2];
            else uv[v * 2 + 1] = -uv[v * 2 + 1];
        }
    }

    for (size_t f = source_nf; f < nf; f++) {
        int32_t vertex[3] = {
            faces[f * 3], faces[f * 3 + 1], faces[f * 3 + 2]
        };
        int32_t fc[3];
        if (stats != NULL) stats->bridge_faces++;
        for (int k = 0; k < 3; k++) {
            if (vertex[k] < 0 || (size_t)vertex[k] >= nv ||
                chart[vertex[k]] < 0 || (size_t)chart[vertex[k]] >= nchart) {
                if (stats != NULL) stats->invalid_faces++;
                goto cleanup;
            }
            fc[k] = chart[vertex[k]];
        }
        if (fc[0] == fc[1] && fc[1] == fc[2]) {
            if (stats != NULL) stats->same_chart_faces++;
            continue;
        }
        if (fc[0] != fc[1] && fc[0] != fc[2] && fc[1] != fc[2]) {
            if (stats != NULL) stats->three_chart_faces++;
            continue;
        }

        int base0, base1, single;
        if (fc[0] == fc[1]) { base0 = 0; base1 = 1; single = 2; }
        else if (fc[1] == fc[2]) { base0 = 1; base1 = 2; single = 0; }
        else { base0 = 2; base1 = 0; single = 1; }
        int32_t base_chart = fc[base0], single_chart = fc[single];
        if (orientation_sign[base_chart] != orientation_sign[single_chart]) {
            if (stats != NULL) stats->invalid_faces++;
            continue;
        }
        const double *u0 = &uv[(size_t)vertex[base0] * 2];
        const double *u1 = &uv[(size_t)vertex[base1] * 2];
        const double *us = &uv[(size_t)vertex[single] * 2];
        double ex = u1[0] - u0[0], ey = u1[1] - u0[1];
        double uv_base = hypot(ex, ey);
        double xyz_base = wwl_distance3(&verts[(size_t)vertex[base0] * 3],
                                        &verts[(size_t)vertex[base1] * 3]);
        double xyz0 = wwl_distance3(&verts[(size_t)vertex[base0] * 3],
                                    &verts[(size_t)vertex[single] * 3]);
        double xyz1 = wwl_distance3(&verts[(size_t)vertex[base1] * 3],
                                    &verts[(size_t)vertex[single] * 3]);
        if (!(uv_base > 1.0e-9) || !(xyz_base > 1.0e-9) ||
            !(xyz0 > 1.0e-9) || !(xyz1 > 1.0e-9)) {
            if (stats != NULL) stats->invalid_faces++;
            continue;
        }
        double scale = uv_base / xyz_base;
        double d0 = xyz0 * scale, d1 = xyz1 * scale;
        double along = (d0 * d0 - d1 * d1 + uv_base * uv_base) /
                       (2.0 * uv_base);
        double height2 = d0 * d0 - along * along;
        double tolerance = 1.0e-8 *
            (d0 * d0 + d1 * d1 + uv_base * uv_base + 1.0);
        if (height2 < -tolerance) {
            if (stats != NULL) stats->invalid_faces++;
            continue;
        }
        double height = sqrt(fmax(0.0, height2));
        double dx = ex / uv_base, dy = ey / uv_base;
        double candidate[2][2] = {
            {u0[0] + along * dx - height * dy,
             u0[1] + along * dy + height * dx},
            {u0[0] + along * dx + height * dy,
             u0[1] + along * dy - height * dx}
        };
        int selected = -1;
        for (int side = 0; side < 2; side++) {
            double point[3][2];
            for (int k = 0; k < 3; k++) {
                point[k][0] = uv[(size_t)vertex[k] * 2];
                point[k][1] = uv[(size_t)vertex[k] * 2 + 1];
            }
            point[single][0] = candidate[side][0];
            point[single][1] = candidate[side][1];
            double area = wwl_area2_points(point[0], point[1], point[2]);
            int sign = area > 1.0e-12 ? 1 : (area < -1.0e-12 ? -1 : 0);
            if (sign == orientation_sign[base_chart]) {
                if (selected >= 0) { selected = -2; break; }
                selected = side;
            }
        }
        if (selected < 0) {
            if (stats != NULL) stats->invalid_faces++;
            continue;
        }
        int32_t lo = base_chart < single_chart ? base_chart : single_chart;
        int32_t hi = base_chart < single_chart ? single_chart : base_chart;
        observation[nobs].lo = lo;
        observation[nobs].hi = hi;
        observation[nobs].base = base_chart;
        observation[nobs].moving = single_chart;
        observation[nobs].source_u = us[0];
        observation[nobs].source_v = us[1];
        observation[nobs].target_u = candidate[selected][0];
        observation[nobs].target_v = candidate[selected][1];
        observation[nobs].maximum_cross_edge = 0.0;
        observation[nobs].maximum_physical_winding_delta = NAN;
        observation[nobs].maximum_raw_phi_delta_turns = NAN;
        observation[nobs].maximum_registered_phi_delta_turns = NAN;
        {
            const int base_index[2] = {base0, base1};
            for (int pair = 0; pair < 2; pair++) {
                int32_t bv = vertex[base_index[pair]], sv = vertex[single];
                double distance = wwl_distance3(
                    &verts[(size_t)bv * 3], &verts[(size_t)sv * 3]);
                if (distance > observation[nobs].maximum_cross_edge)
                    observation[nobs].maximum_cross_edge = distance;
                if (audit != NULL && audit->calib != NULL) {
                    double winding = wwl_physical_winding_delta(
                        audit->calib, &verts[(size_t)bv * 3],
                        &verts[(size_t)sv * 3]);
                    if (isfinite(winding) &&
                        (!isfinite(observation[nobs].maximum_physical_winding_delta) ||
                         winding > observation[nobs].maximum_physical_winding_delta))
                        observation[nobs].maximum_physical_winding_delta = winding;
                }
                if (audit != NULL && audit->anchor_source != NULL &&
                    audit->source != NULL) {
                    int32_t bs = audit->anchor_source[bv];
                    int32_t ss = audit->anchor_source[sv];
                    if (bs >= 0 && ss >= 0 &&
                        (size_t)bs < audit->source->nv &&
                        (size_t)ss < audit->source->nv) {
                        double raw_delta = fabs(
                            (double)audit->source->phi[(size_t)bs] -
                            (double)audit->source->phi[(size_t)ss]) / WWL_2PI;
                        if (isfinite(raw_delta) &&
                            (!isfinite(observation[nobs].maximum_raw_phi_delta_turns) ||
                             raw_delta > observation[nobs].maximum_raw_phi_delta_turns))
                            observation[nobs].maximum_raw_phi_delta_turns = raw_delta;
                        if (audit->registered_phi != NULL) {
                            double registered_delta = fabs(
                                (double)audit->registered_phi[(size_t)bs] -
                                (double)audit->registered_phi[(size_t)ss]) /
                                WWL_2PI;
                            if (isfinite(registered_delta) &&
                                (!isfinite(observation[nobs].maximum_registered_phi_delta_turns) ||
                                 registered_delta > observation[nobs].maximum_registered_phi_delta_turns))
                                observation[nobs].maximum_registered_phi_delta_turns =
                                    registered_delta;
                        }
                    }
                }
            }
        }
        nobs++;
        if (stats != NULL) stats->two_chart_faces++;
    }
    if (stats != NULL) stats->observations = nobs;
    if (nobs == 0 || (stats != NULL && stats->three_chart_faces != 0))
        goto cleanup;

    qsort(observation, nobs, sizeof(*observation),
          wwl_compare_bridge_observation);
    edge = (WwlBridgeGaugeEdge *)malloc(nobs * sizeof(*edge));
    scratch = (double *)malloc(nobs * sizeof(*scratch));
    all_residual = (double *)malloc(nobs * sizeof(*all_residual));
    if (edge == NULL || scratch == NULL || all_residual == NULL) goto cleanup;
    for (size_t first = 0; first < nobs;) {
        size_t last = first + 1;
        while (last < nobs && observation[last].lo == observation[first].lo &&
               observation[last].hi == observation[first].hi)
            last++;
        WwlBridgeGaugeEdge fit_lo, fit_hi;
        int have_lo = wwl_fit_bridge_isometry(
            observation, first, last, observation[first].lo, allow_rotation,
            &fit_lo, scratch) == 0;
        int have_hi = wwl_fit_bridge_isometry(
            observation, first, last, observation[first].hi, allow_rotation,
            &fit_hi, scratch) == 0;
        if (!have_lo && !have_hi) goto cleanup;
        WwlBridgeGaugeEdge chosen;
        if (!have_hi || (have_lo &&
            (fit_lo.observations > fit_hi.observations ||
             (fit_lo.observations == fit_hi.observations &&
              fit_lo.median_residual <= fit_hi.median_residual))))
            chosen = fit_lo;
        else
            chosen = fit_hi;
        {
            WwlBridgeGaugeEdge translation_fit;
            if (wwl_fit_bridge_isometry(
                    observation, first, last, chosen.base, 0,
                    &translation_fit, scratch) != 0)
                goto cleanup;
            chosen.pair_observations = last - first;
            chosen.lo_observations = chosen.hi_observations = 0;
            for (size_t i = first; i < last; i++) {
                if (observation[i].base == chosen.lo)
                    chosen.lo_observations++;
                else if (observation[i].base == chosen.hi)
                    chosen.hi_observations++;
            }
            chosen.translation_u = translation_fit.translate_u;
            chosen.translation_v = translation_fit.translate_v;
            chosen.translation_median_residual =
                translation_fit.median_residual;
            chosen.translation_maximum_residual =
                translation_fit.maximum_residual;
            /* The unconstrained Procrustes angle is an L2 fit and can be
             * captured by a handful of branch/outlier observations.  The
             * identity rotation plus robust median translation is part of
             * the same proper-isometry family and preserves winding/axial
             * semantics.  Apply whichever has the lower robust residual;
             * retain the translation diagnostics either way. */
            if (translation_fit.median_residual < chosen.median_residual ||
                (translation_fit.median_residual == chosen.median_residual &&
                 translation_fit.maximum_residual < chosen.maximum_residual)) {
                chosen.cosine = translation_fit.cosine;
                chosen.sine = translation_fit.sine;
                chosen.translate_u = translation_fit.translate_u;
                chosen.translate_v = translation_fit.translate_v;
                chosen.median_residual = translation_fit.median_residual;
                chosen.maximum_residual = translation_fit.maximum_residual;
                if (stats != NULL) stats->translation_selected_edges++;
            }

            size_t count = 0;
            chosen.cross_edge_maximum = 0.0;
            for (size_t i = first; i < last; i++) {
                double value = observation[i].maximum_cross_edge;
                scratch[count++] = value;
                if (value > chosen.cross_edge_maximum)
                    chosen.cross_edge_maximum = value;
            }
            chosen.cross_edge_median = wwl_median(scratch, count);

            count = 0;
            chosen.physical_winding_maximum = NAN;
            for (size_t i = first; i < last; i++) {
                double value = observation[i].maximum_physical_winding_delta;
                if (!isfinite(value)) continue;
                scratch[count++] = value;
                if (!isfinite(chosen.physical_winding_maximum) ||
                    value > chosen.physical_winding_maximum)
                    chosen.physical_winding_maximum = value;
            }
            chosen.physical_winding_median =
                count ? wwl_median(scratch, count) : NAN;

            count = 0;
            chosen.raw_phi_maximum_turns = NAN;
            for (size_t i = first; i < last; i++) {
                double value = observation[i].maximum_raw_phi_delta_turns;
                if (!isfinite(value)) continue;
                scratch[count++] = value;
                if (!isfinite(chosen.raw_phi_maximum_turns) ||
                    value > chosen.raw_phi_maximum_turns)
                    chosen.raw_phi_maximum_turns = value;
            }
            chosen.raw_phi_observations = count;
            chosen.raw_phi_median_turns =
                count ? wwl_median(scratch, count) : NAN;

            count = 0;
            chosen.registered_phi_maximum_turns = NAN;
            for (size_t i = first; i < last; i++) {
                double value =
                    observation[i].maximum_registered_phi_delta_turns;
                if (!isfinite(value)) continue;
                scratch[count++] = value;
                if (!isfinite(chosen.registered_phi_maximum_turns) ||
                    value > chosen.registered_phi_maximum_turns)
                    chosen.registered_phi_maximum_turns = value;
            }
            chosen.registered_phi_observations = count;
            chosen.registered_phi_median_turns =
                count ? wwl_median(scratch, count) : NAN;
        }
        edge[nedge] = chosen;
        for (size_t i = first; i < last; i++) {
            if (observation[i].base != chosen.base) continue;
            double px = chosen.cosine * observation[i].source_u -
                        chosen.sine * observation[i].source_v +
                        chosen.translate_u;
            double py = chosen.sine * observation[i].source_u +
                        chosen.cosine * observation[i].source_v +
                        chosen.translate_v;
            all_residual[nresidual++] =
                hypot(px - observation[i].target_u,
                      py - observation[i].target_v);
        }
        nedge++;
        first = last;
    }
    if (stats != NULL) stats->gauge_edges = nedge;

    if (wwl_dsu_init(&dsu, nchart) != 0) goto cleanup;
    have_dsu = 1;
    size_t cycles = 0;
    for (size_t e = 0; e < nedge; e++) {
        if (wwl_dsu_find(&dsu, edge[e].lo) ==
            wwl_dsu_find(&dsu, edge[e].hi))
            cycles++;
        else
            wwl_dsu_union(&dsu, edge[e].lo, edge[e].hi);
    }
    size_t graph_components = 0;
    for (size_t c = 0; c < nchart; c++)
        if (wwl_dsu_find(&dsu, (int32_t)c) == (int32_t)c)
            graph_components++;
    if (stats != NULL) {
        stats->gauge_components = graph_components;
        stats->gauge_cycles = cycles;
    }
    if (cycles != 0 || graph_components != expected_components ||
        nedge + graph_components != nchart) {
        fprintf(stderr,
            "welded_winding_lift: bridge chart graph certificate failed: "
            "charts=%zu edges=%zu components=%zu expected=%zu cycles=%zu\n",
            nchart, nedge, graph_components, expected_components, cycles);
        goto cleanup;
    }

    degree = (size_t *)calloc(nchart, sizeof(*degree));
    offset = (size_t *)malloc((nchart + 1) * sizeof(*offset));
    cursor = (size_t *)malloc(nchart * sizeof(*cursor));
    adj_edge = (int32_t *)malloc(2 * nedge * sizeof(*adj_edge));
    queue = (int32_t *)malloc(nchart * sizeof(*queue));
    rotate_c = (double *)malloc(nchart * sizeof(*rotate_c));
    rotate_s = (double *)malloc(nchart * sizeof(*rotate_s));
    shift_u = (double *)malloc(nchart * sizeof(*shift_u));
    shift_v = (double *)malloc(nchart * sizeof(*shift_v));
    if (degree == NULL || offset == NULL || cursor == NULL ||
        adj_edge == NULL || queue == NULL || rotate_c == NULL ||
        rotate_s == NULL || shift_u == NULL || shift_v == NULL)
        goto cleanup;
    for (size_t e = 0; e < nedge; e++) {
        degree[edge[e].lo]++;
        degree[edge[e].hi]++;
    }
    offset[0] = 0;
    for (size_t c = 0; c < nchart; c++) offset[c + 1] = offset[c] + degree[c];
    memcpy(cursor, offset, nchart * sizeof(*cursor));
    for (size_t e = 0; e < nedge; e++) {
        adj_edge[cursor[edge[e].lo]++] = (int32_t)e;
        adj_edge[cursor[edge[e].hi]++] = (int32_t)e;
    }
    if (audit != NULL && audit->path != NULL) {
        if (wwl_write_bridge_edge_audit(
                audit, nv, chart, nchart, edge, nedge, offset,
                adj_edge) != 0)
            goto cleanup;
        fprintf(stderr,
            "welded_winding_lift: wrote %zu chart-level bridge generator "
            "record(s) to %s\n", nedge, audit->path);
    }
    for (size_t c = 0; c < nchart; c++) {
        rotate_c[c] = rotate_s[c] = NAN;
        shift_u[c] = shift_v[c] = NAN;
    }
    double maximum_cycle_residual = 0.0;
    size_t propagated_components = 0;
    for (size_t seed = 0; seed < nchart; seed++) {
        if (isfinite(shift_u[seed])) continue;
        propagated_components++;
        size_t head = 0, tail = 0;
        rotate_c[seed] = 1.0; rotate_s[seed] = 0.0;
        shift_u[seed] = shift_v[seed] = 0.0;
        queue[tail++] = (int32_t)seed;
        while (head < tail) {
            int32_t here = queue[head++];
            for (size_t at = offset[here]; at < offset[here + 1]; at++) {
                WwlBridgeGaugeEdge *ge = &edge[adj_edge[at]];
                int32_t other;
                double candidate_c, candidate_s, candidate_u, candidate_v;
                if (here == ge->base) {
                    other = ge->moving;
                    candidate_c = rotate_c[here] * ge->cosine -
                                  rotate_s[here] * ge->sine;
                    candidate_s = rotate_s[here] * ge->cosine +
                                  rotate_c[here] * ge->sine;
                    candidate_u = rotate_c[here] * ge->translate_u -
                                  rotate_s[here] * ge->translate_v +
                                  shift_u[here];
                    candidate_v = rotate_s[here] * ge->translate_u +
                                  rotate_c[here] * ge->translate_v +
                                  shift_v[here];
                } else {
                    other = ge->base;
                    candidate_c = rotate_c[here] * ge->cosine +
                                  rotate_s[here] * ge->sine;
                    candidate_s = rotate_s[here] * ge->cosine -
                                  rotate_c[here] * ge->sine;
                    candidate_u = shift_u[here] -
                        (candidate_c * ge->translate_u -
                         candidate_s * ge->translate_v);
                    candidate_v = shift_v[here] -
                        (candidate_s * ge->translate_u +
                         candidate_c * ge->translate_v);
                }
                if (!isfinite(shift_u[other])) {
                    rotate_c[other] = candidate_c;
                    rotate_s[other] = candidate_s;
                    shift_u[other] = candidate_u;
                    shift_v[other] = candidate_v;
                    queue[tail++] = other;
                } else {
                    double translation_residual =
                        hypot(candidate_u - shift_u[other],
                              candidate_v - shift_v[other]);
                    double rotation_residual = fabs(atan2(
                        candidate_s * rotate_c[other] -
                            candidate_c * rotate_s[other],
                        candidate_c * rotate_c[other] +
                            candidate_s * rotate_s[other]));
                    double residual = translation_residual + rotation_residual;
                    if (residual > maximum_cycle_residual)
                        maximum_cycle_residual = residual;
                }
            }
        }
    }
    if (propagated_components != expected_components) goto cleanup;
    if (apply_translation)
        for (size_t v = 0; v < nv; v++) {
            int32_t c = chart[v];
            if (c < 0 || (size_t)c >= nchart) goto cleanup;
            double u = uv[v * 2], w = uv[v * 2 + 1];
            if (axial_only) {
                uv[v * 2] = u + shift_u[c];
            } else {
                uv[v * 2] = rotate_c[c] * u - rotate_s[c] * w + shift_u[c];
                uv[v * 2 + 1] = rotate_s[c] * u + rotate_c[c] * w + shift_v[c];
            }
        }

    if (stats != NULL) {
        stats->observation_residual_median =
            wwl_median(all_residual, nresidual);
        stats->observation_residual_maximum = 0.0;
        for (size_t i = 0; i < nresidual; i++)
            if (all_residual[i] > stats->observation_residual_maximum)
                stats->observation_residual_maximum = all_residual[i];
        stats->cycle_residual_maximum = maximum_cycle_residual;
    }
    result = 0;

cleanup:
    if (have_dsu) wwl_dsu_dispose(&dsu);
    free(orientation_balance); free(orientation_sign); free(observation);
    free(edge); free(scratch); free(all_residual); free(degree); free(offset);
    free(cursor); free(adj_edge); free(queue); free(rotate_c); free(rotate_s);
    free(shift_u); free(shift_v);
    return result;
}

static int wwl_selftest(void)
{
    WwlDsu dsu;
    if (wwl_dsu_init(&dsu, 4) != 0) return 1;
    wwl_dsu_union(&dsu, 0, 1);
    wwl_dsu_union(&dsu, 2, 3);
    int ok = wwl_dsu_find(&dsu, 0) == wwl_dsu_find(&dsu, 1) &&
             wwl_dsu_find(&dsu, 0) != wwl_dsu_find(&dsu, 2) &&
             fabs(wwl_wrap_pi(3.0 * WWL_PI) + WWL_PI) < 1.0e-12;
    wwl_dsu_dispose(&dsu);
    {
        /* Position-only lookup cannot distinguish these coincident vertices.
         * GWLIN2 cube/local provenance must recover both exact source UV
         * anchors without choosing or averaging. */
        float source_verts[2 * 3] = {1,2,3, 1,2,3};
        float welded_verts[2 * 3] = {1,2,3, 1,2,3};
        char source_ids[2][48] = {{0}};
        char lineage_ids[2][48] = {{0}};
        size_t cube_voff[3] = {0,1,2};
        int32_t parent0[2] = {-2,-3};
        int32_t parent1[2] = {0,0};
        int32_t anchor[2] = {-1,-1};
        size_t matched = 0, missing = 0, ambiguous = 0;
        PieceSet source;
        WwlLineage lineage;
        memset(&source, 0, sizeof source);
        memset(&lineage, 0, sizeof lineage);
        snprintf(source_ids[0], sizeof(source_ids[0]), "cube-a");
        snprintf(source_ids[1], sizeof(source_ids[1]), "cube-b");
        memcpy(lineage_ids, source_ids, sizeof source_ids);
        source.verts = source_verts;
        source.nv = 2;
        source.ids = source_ids;
        source.cube_voff = cube_voff;
        source.n_cubes = 2;
        lineage.parent0 = parent0;
        lineage.parent1 = parent1;
        lineage.cube_ids = lineage_ids;
        lineage.nv = 2;
        lineage.n_cubes = 2;
        lineage.source_anchors = 2;
        lineage.version = 2;
        ok = ok && wwl_map_source_lineage(
                &source, welded_verts, 2, &lineage, anchor,
                &matched, &missing, &ambiguous) == 0 &&
             matched == 2 && missing == 0 && ambiguous == 0 &&
             anchor[0] == 0 && anchor[1] == 1;
    }
    {
        /* grid_weld splits one source vertex into two welded vertices whenever
         * it peels a bowtie or cuts a source chart to a disk.  Both halves keep
         * the same (cube, local) provenance and the same position, so both must
         * inherit that source vertex's UV.  Rejecting the whole mapping over a
         * shared parent threw away 28.8M good anchors on the 4x21x21 grid. */
        float source_verts[1 * 3] = {1,2,3};
        float welded_verts[2 * 3] = {1,2,3, 1,2,3};
        char source_ids[1][48] = {{0}};
        char lineage_ids[1][48] = {{0}};
        size_t cube_voff[2] = {0,1};
        int32_t parent0[2] = {-2,-2};
        int32_t parent1[2] = {0,0};
        int32_t anchor[2] = {-1,-1};
        size_t matched = 0, missing = 0, ambiguous = 0;
        PieceSet source;
        WwlLineage lineage;
        memset(&source, 0, sizeof source);
        memset(&lineage, 0, sizeof lineage);
        snprintf(source_ids[0], sizeof(source_ids[0]), "cube-a");
        memcpy(lineage_ids, source_ids, sizeof source_ids);
        source.verts = source_verts;
        source.nv = 1;
        source.ids = source_ids;
        source.cube_voff = cube_voff;
        source.n_cubes = 1;
        lineage.parent0 = parent0;
        lineage.parent1 = parent1;
        lineage.cube_ids = lineage_ids;
        lineage.nv = 2;
        lineage.n_cubes = 1;
        lineage.source_anchors = 2;
        lineage.version = 2;
        ok = ok && wwl_map_source_lineage(
                &source, welded_verts, 2, &lineage, anchor,
                &matched, &missing, &ambiguous) == 0 &&
             matched == 2 && missing == 0 && ambiguous == 1 &&
             anchor[0] == 0 && anchor[1] == 0;
        fprintf(stderr,
            "[selftest] weld-split shared parent matched=%zu missing=%zu "
            "shared=%zu anchors=%d/%d\n",
            matched, missing, ambiguous, anchor[0], anchor[1]);
    }
    {
        /* Registered per-vertex polish is deliberately reduced to one rigid
         * isometry per chart; raw triangle geometry must remain exact. */
        const int32_t chart[6] = {0,0,0,1,1,1};
        const int32_t anchor[6] = {0,1,2,3,4,5};
        const float registered[12] = {
            10.1f,20.0f, 10.9f,20.0f, 10.0f,21.0f,
            -5.0f,3.0f, -6.0f,3.0f, -5.0f,4.0f
        };
        double uv[12] = {
            0,0, 1,0, 0,1,
            0,0, 1,0, 0,1
        };
        WwlChartGaugeStats gauge;
        int rc = wwl_apply_registered_chart_gauges(
            6, 6, chart, 2, anchor, registered, uv, &gauge);
        fprintf(stderr,
            "[selftest] registered isometry rc=%d obs=%zu reflected=%zu "
            "point=%zu residual=%.6g/%.6g edge=%.12g/%.12g\n",
            rc, gauge.observations, gauge.reflected_charts,
            gauge.point_gauge_charts, gauge.residual_median,
            gauge.residual_maximum,
            hypot(uv[0] - uv[2], uv[1] - uv[3]),
            hypot(uv[6] - uv[8], uv[7] - uv[9]));
        ok = ok && rc == 0 && gauge.observations == 6 &&
             gauge.reflected_charts == 1 &&
             fabs(hypot(uv[0] - uv[2], uv[1] - uv[3]) - 1.0) < 1.0e-12 &&
             fabs(hypot(uv[6] - uv[8], uv[7] - uv[9]) - 1.0) < 1.0e-12 &&
              gauge.residual_maximum < 0.25;
    }
    {
        /* A topology phase chooses only an integer chart gauge.  The actual U
         * update must remain the shipped nonlinear winding correction, while
         * axial V is untouched. */
        const int32_t chart[4] = {0,0,1,1};
        const int32_t anchor[4] = {0,1,2,3};
        float phi[4] = {0.1f,0.2f,-0.3f,-0.2f};
        double lift[4] = {
            (double)phi[0] + WWL_2PI, (double)phi[1] + WWL_2PI,
            (double)phi[2] - 2.0 * WWL_2PI,
            (double)phi[3] - 2.0 * WWL_2PI
        };
        double uv[8] = {1,10, 2,11, 3,12, 4,13};
        const double original_u[4] = {1,2,3,4};
        PieceSet source;
        ScaffoldCalib calib;
        WwlTopologyGaugeStats gauge;
        memset(&source, 0, sizeof source);
        Scaffold_calib_default(&calib);
        source.phi = phi;
        source.nv = 4;
        int rc = wwl_apply_topology_chart_turns(
            4, chart, 2, anchor, &source, &calib, lift, uv, &gauge);
        ok = ok && rc == 0 && gauge.charts == 2 && gauge.anchors == 4 &&
             gauge.minimum_turn == -2 && gauge.maximum_turn == 1 &&
             gauge.residual_maximum_turns < 1.0e-12 &&
             fabs(uv[0] - (original_u[0] + CubeReg_deltaU(
                 calib.spiral_a, calib.spiral_b, (double)phi[0], 1))) < 1.0e-10 &&
             fabs(uv[4] - (original_u[2] + CubeReg_deltaU(
                 calib.spiral_a, calib.spiral_b, (double)phi[2], -2))) < 1.0e-10 &&
             uv[1] == 10.0 && uv[3] == 11.0 &&
             uv[5] == 12.0 && uv[7] == 13.0;
    }
    {
        /* Two independently translated safe charts joined by one bridge
         * triangle.  Intrinsic development of the bridge must recover the
         * exact rigid chart translation without deforming either chart. */
        const float verts[18] = {
            0,0,0, 1,0,0, 0,1,0,
            1,1,0, 2,1,0, 1,2,0
        };
        const int32_t faces[9] = {0,1,2, 3,4,5, 1,3,2};
        const int32_t chart[6] = {0,0,0, 1,1,1};
        double uv[12] = {
            0,0, 1,0, 0,1,
            101,1, 102,1, 101,2
        };
        WwlBridgeGaugeStats gauge;
        int rc = wwl_assemble_bridge_gauges(
            verts, 6, faces, 2, 3, chart, 2, 1, 0, 1, 1, uv, &gauge, NULL);
        ok = ok && rc == 0 && gauge.observations == 1 &&
             gauge.gauge_edges == 1 && gauge.gauge_components == 1 &&
             gauge.gauge_cycles == 0 &&
             fabs(uv[6] - 1.0) < 1.0e-12 &&
             fabs(uv[7] - 1.0) < 1.0e-12 &&
              fabs(uv[10] - 1.0) < 1.0e-12 &&
              fabs(uv[11] - 2.0) < 1.0e-12;
    }
    {
        /* A generated chart with no source anchor receives the strict
         * topology target, while an anchored winding chart is untouched. */
        const int32_t faces[6] = {0,1,2, 3,4,5};
        const int32_t chart[6] = {0,0,0, 1,1,1};
        const int32_t anchor[6] = {0,1,2, -1,-1,-1};
        const double target[12] = {
            10,10, 11,10, 10,11,
            20,20, 21,20, 20,21
        };
        double uv[12] = {
            0,0, 1,0, 0,1,
            0,7, 0,8, 0,9
        };
        size_t seeded = 0;
        int rc = wwl_seed_invalid_chart_targets(
            6, faces, 2, chart, 2, anchor, target, uv, 0, &seeded);
        ok = ok && rc == 0 && seeded == 1 &&
             uv[0] == 0.0 && uv[2] == 1.0 &&
             uv[6] == 20.0 && uv[7] == 20.0 &&
             uv[10] == 20.0 && uv[11] == 21.0;
    }
    {
        /* carry mode: phi rides the same anchor / lineage-midpoint / harmonic
         * machinery as uv.  v3 is an exact lineage midpoint of v0,v1; v4 is a
         * free vertex equidistant from fixed v1,v2, so one Jacobi sweep lands
         * it on their mean. */
        const float verts[15] = {
            0,0,0, 0,0,2, 0,2,0, 0,0,1, 0,2,2
        };
        const int32_t faces[9] = {0,3,1, 0,1,2, 1,4,2};
        const int32_t chart[5] = {0,0,0,0,0};
        int32_t anchor[5] = {0,1,2,-1,-1};
        float source_uv[6] = {0,0, 2,0, 0,2};
        float source_phi[3] = {1.0f, 3.0f, 5.0f};
        int32_t parent0[5] = {-1,-1,-1, 0,-1};
        int32_t parent1[5] = {-1,-1,-1, 1,-1};
        double uv[10];
        float phi[5] = {0};
        PieceSet source;
        WwlLineage lineage;
        WwlGraph graph;
        size_t unknown = 0, mid = 0, resid = 0;
        double delta = 0.0;
        memset(&source, 0, sizeof source);
        memset(&lineage, 0, sizeof lineage);
        memset(&graph, 0, sizeof graph);
        source.uv = source_uv;
        source.phi = source_phi;
        source.nv = 3;
        lineage.parent0 = parent0;
        lineage.parent1 = parent1;
        lineage.nv = 5;
        int rc = wwl_graph_build(5, faces, 3, &graph) != 0 ||
                 wwl_interpolate_chart_uv(verts, 5, &graph, chart, anchor,
                                          &source, &lineage, 64, uv, phi,
                                          &unknown, &mid, &resid, &delta) != 0;
        ok = ok && rc == 0 && mid == 1 &&
             fabs(uv[6] - 1.0) < 1.0e-12 && fabs(uv[7] - 0.0) < 1.0e-12 &&
             fabs((double)phi[3] - 2.0) < 1.0e-6 &&
             fabs(uv[8] - 1.0) < 1.0e-5 && fabs(uv[9] - 1.0) < 1.0e-5 &&
             fabs((double)phi[4] - 4.0) < 1.0e-5;
        fprintf(stderr,
            "[selftest] carry phi: mid=%zu uv3=(%.3f,%.3f) phi3=%.3f "
            "uv4=(%.3f,%.3f) phi4=%.3f -> %s\n",
            mid, uv[6], uv[7], (double)phi[3], uv[8], uv[9],
            (double)phi[4], ok ? "ok" : "FAIL");
        wwl_graph_dispose(&graph);
    }
    {
        /* carry winsorization: one 9-anchor chart with a 5000-vox u-tail
         * anchor loses exactly that anchor; a 3-anchor chart keeps its
         * outlier (below the 8-anchor floor). */
        int32_t chart[12] = {0,0,0,0,0,0,0,0,0, 1,1,1};
        int32_t anchor[12] = {0,1,2,3,4,5,6,7,8, 9,10,11};
        float source_uv[24] = {
            0,0, 0.5f,0, 1,0, 1.5f,0, 2,0, 2.5f,0, 3,0, 3.5f,0, 5000,0,
            0,0, 1,0, 9000,0
        };
        float source_phi[12] = {0};
        PieceSet source;
        memset(&source, 0, sizeof source);
        source.uv = source_uv;
        source.phi = source_phi;
        source.nv = 12;
        size_t demoted = wwl_carry_winsorize_anchors(
            12, chart, 2, anchor, &source, 1000.0);
        ok = ok && demoted == 1 && anchor[8] == -1 &&
             anchor[0] == 0 && anchor[7] == 7 && anchor[11] == 11;
        fprintf(stderr,
            "[selftest] carry winsorize: demoted=%zu tail=%d small-chart=%d "
            "-> %s\n", demoted, anchor[8], anchor[11], ok ? "ok" : "FAIL");
    }
    {
        /* carry winsorization, MULTI-WRAP chart: u tracks phi over three
         * turns (span ~5700 vox >> the 1000-vox gate); the phi-linear
         * residual must keep every true anchor and demote only the one
         * genuine 8000-vox artifact.  The old raw |u - median| gate
         * demoted the whole outer/inner windings of such charts and let
         * the harmonic fill stack the carried field (z6400: 29k demotions,
         * 27.9% fill). */
        enum { MW = 40 };
        int32_t chart[MW];
        int32_t anchor[MW];
        float source_uv[MW * 2];
        float source_phi[MW];
        PieceSet source;
        for (int i = 0; i < MW; i++) {
            double ph = (double)i * (18.84955592153876 / (MW - 1)); /* 6pi */
            chart[i] = 0;
            anchor[i] = i;
            source_phi[i] = (float)ph;
            source_uv[i * 2] = (float)(300.0 * ph);
            source_uv[i * 2 + 1] = 0.0f;
        }
        source_uv[(MW / 2) * 2] += 8000.0f;   /* the one true artifact */
        memset(&source, 0, sizeof source);
        source.uv = source_uv;
        source.phi = source_phi;
        source.nv = MW;
        size_t demoted = wwl_carry_winsorize_anchors(
            MW, chart, 1, anchor, &source, 1000.0);
        int kept = 1;
        for (int i = 0; i < MW; i++)
            if (i != MW / 2 && anchor[i] != i) kept = 0;
        ok = ok && demoted == 1 && anchor[MW / 2] == -1 && kept;
        fprintf(stderr,
            "[selftest] carry winsorize multi-wrap: demoted=%zu "
            "artifact=%d kept=%d -> %s\n", demoted, anchor[MW / 2],
            kept, ok ? "ok" : "FAIL");
    }
    {
        /* carry anchorless seeding: chart 1 (no anchors) borrows provenance
         * from chart 0 across the bridge face; the isolated chart 2 stays
         * dry and is reported, not silently zeroed. */
        const int32_t faces[12] = {0,1,2, 3,4,5, 1,3,4, 6,7,8};
        int32_t chart[9] = {0,0,0, 1,1,1, 2,2,2};
        int32_t anchor[9] = {0,1,2, -1,-1,-1, -1,-1,-1};
        WwlGraph graph;
        size_t seeded = 0, dry = 0;
        memset(&graph, 0, sizeof graph);
        int rc = wwl_graph_build(9, faces, 4, &graph);
        if (rc == 0)
            wwl_carry_seed_anchorless(9, &graph, chart, 3, anchor,
                                      &seeded, &dry);
        ok = ok && rc == 0 && seeded == 1 && dry == 1 &&
             anchor[3] == 1 && anchor[4] == 1 && anchor[6] == -1;
        fprintf(stderr,
            "[selftest] carry anchorless: seeded=%zu dry=%zu a3=%d a4=%d "
            "-> %s\n", seeded, dry, anchor[3], anchor[4],
            ok ? "ok" : "FAIL");
        wwl_graph_dispose(&graph);
    }
    {
        const char *path="welded_winding_lift_chart_layout_selftest.bin";
        const char magic[8]={'G','W','L','A','Y','T','1','\n'};
        uint32_t version=1,flags=0;uint64_t nodes=2,edges=3,cycles=1;
        char ids[2][48]={{0}};int32_t local[2]={0,0},rep[2]={0,0};
        double delta[2]={2.5,-1.5};int32_t turn[2]={0,0},component_id[2]={0,0};
        char source_ids[2][48]={{0}};size_t cube_voff[3]={0,1,2};
        int32_t anchor[2]={0,1},chart[2]={0,1};double uv[4]={10,0,20,0};
        float test_verts[6]={0,0,0,0,0,0};
        PieceSet source;WwlChartLayout layout;
        WwlChartLayoutCorrection *correction=NULL;
        size_t mapped=0;double minimum=0,maximum=0;FILE *file=NULL;int rc=-1;
        snprintf(ids[0],sizeof(ids[0]),"cube-a");
        snprintf(ids[1],sizeof(ids[1]),"cube-b");
        memcpy(source_ids,ids,sizeof(ids));memset(&source,0,sizeof(source));
        memset(&layout,0,sizeof(layout));source.ids=source_ids;
        source.cube_voff=cube_voff;source.n_cubes=2;source.nv=2;
        file=fopen(path,"wb");
        if(file&&fwrite(magic,1,sizeof(magic),file)==sizeof(magic)&&
           fwrite(&version,sizeof(version),1,file)==1&&
           fwrite(&flags,sizeof(flags),1,file)==1&&
           fwrite(&nodes,sizeof(nodes),1,file)==1&&
           fwrite(&edges,sizeof(edges),1,file)==1&&
           fwrite(&cycles,sizeof(cycles),1,file)==1){
            rc=0;
            for(size_t i=0;i<2&&rc==0;i++)
                if(fwrite(ids[i],1,sizeof(ids[i]),file)!=sizeof(ids[i])||
                   fwrite(&local[i],sizeof(local[i]),1,file)!=1||
                   fwrite(&rep[i],sizeof(rep[i]),1,file)!=1||
                   fwrite(&delta[i],sizeof(delta[i]),1,file)!=1||
                   fwrite(&turn[i],sizeof(turn[i]),1,file)!=1||
                   fwrite(&component_id[i],sizeof(component_id[i]),1,file)!=1)
                    rc=-1;
        }
        if(file&&fclose(file)!=0)rc=-1;
        if(rc==0)rc=wwl_chart_layout_load(path,&layout);
        if(rc==0)rc=wwl_chart_layout_map(
            &layout,&source,2,anchor,chart,2,&correction,&mapped,
            &minimum,&maximum);
        if(rc==0)rc=wwl_chart_layout_apply(
            2,test_verts,uv,chart,2,correction);
        ok=ok&&rc==0&&mapped==2&&layout.relation_edges==3&&
           layout.cycle_edges==1&&fabs(minimum+1.5)<1e-12&&
           fabs(maximum-2.5)<1e-12&&fabs(uv[0]-12.5)<1e-12&&
           fabs(uv[2]-18.5)<1e-12;
        fprintf(stderr,
            "[selftest] chart graph layout load/map/apply: nodes=%zu "
            "edges=%zu cycles=%zu -> %s\n",mapped,layout.relation_edges,
            layout.cycle_edges,ok?"ok":"FAIL");
        remove(path);free(correction);wwl_chart_layout_dispose(&layout);
    }
    {
        const char *path="welded_winding_lift_chart_layout_v2_selftest.bin";
        const char magic[8]={'G','W','L','A','Y','T','1','\n'};
        uint32_t version=2,flags=0;uint64_t nodes=1,edges=0,cycles=0;
        char id[48]={0},source_ids[1][48]={{0}};
        int32_t local=0,rep=0,turn=0,component_id=0;
        double shift=2.0,slope=0.25,reference=10.0;
        size_t cube_voff[2]={0,1};int32_t anchor[1]={0},chart[1]={0};
        float test_verts[3]={14,0,0};double uv[2]={7,0};
        PieceSet source;WwlChartLayout layout;
        WwlChartLayoutCorrection *correction=NULL;
        size_t mapped=0;double minimum=0,maximum=0;FILE *file=NULL;int rc=-1;
        snprintf(id,sizeof(id),"cube-v");memcpy(source_ids[0],id,sizeof(id));
        memset(&source,0,sizeof(source));memset(&layout,0,sizeof(layout));
        source.ids=source_ids;source.cube_voff=cube_voff;
        source.n_cubes=1;source.nv=1;
        file=fopen(path,"wb");
        if(file&&fwrite(magic,1,sizeof(magic),file)==sizeof(magic)&&
           fwrite(&version,sizeof(version),1,file)==1&&
           fwrite(&flags,sizeof(flags),1,file)==1&&
           fwrite(&nodes,sizeof(nodes),1,file)==1&&
           fwrite(&edges,sizeof(edges),1,file)==1&&
           fwrite(&cycles,sizeof(cycles),1,file)==1&&
           fwrite(id,1,sizeof(id),file)==sizeof(id)&&
           fwrite(&local,sizeof(local),1,file)==1&&
           fwrite(&rep,sizeof(rep),1,file)==1&&
           fwrite(&shift,sizeof(shift),1,file)==1&&
           fwrite(&slope,sizeof(slope),1,file)==1&&
           fwrite(&reference,sizeof(reference),1,file)==1&&
           fwrite(&turn,sizeof(turn),1,file)==1&&
           fwrite(&component_id,sizeof(component_id),1,file)==1)rc=0;
        if(file&&fclose(file)!=0)rc=-1;
        if(rc==0)rc=wwl_chart_layout_load(path,&layout);
        if(rc==0)rc=wwl_chart_layout_map(
            &layout,&source,1,anchor,chart,1,&correction,&mapped,
            &minimum,&maximum);
        if(rc==0)rc=wwl_chart_layout_apply(
            1,test_verts,uv,chart,1,correction);
        ok=ok&&rc==0&&mapped==1&&layout.version==2&&
           fabs(uv[0]-10.0)<1e-12;
        fprintf(stderr,
            "[selftest] chart graph V-dependent layout: nodes=%zu -> %s\n",
            mapped,ok?"ok":"FAIL");
        remove(path);free(correction);wwl_chart_layout_dispose(&layout);
    }
    {
        const char *path="welded_winding_lift_chart_layout_v3_selftest.bin";
        const char magic[8]={'G','W','L','A','Y','T','1','\n'};
        uint32_t version=3,flags=0;uint64_t nodes=1,edges=8,cycles=2;
        char id[48]={0},source_ids[1][48]={{0}};
        int32_t local=0,rep=0,turn=0,component_id=0;
        double reference=10.0,step=2.0;
        double shifts[WWL_CHART_LAYOUT_V3_KNOTS]={0,1,2,3,4,5,6,7,8};
        size_t cube_voff[2]={0,1};int32_t anchor[1]={0},chart[1]={0};
        float test_verts[3]={15,0,0};double uv[2]={7,0};
        PieceSet source;WwlChartLayout layout;
        WwlChartLayoutCorrection *correction=NULL;
        size_t mapped=0;double minimum=0,maximum=0;FILE *file=NULL;int rc=-1;
        snprintf(id,sizeof(id),"cube-v3");memcpy(source_ids[0],id,sizeof(id));
        memset(&source,0,sizeof(source));memset(&layout,0,sizeof(layout));
        source.ids=source_ids;source.cube_voff=cube_voff;
        source.n_cubes=1;source.nv=1;
        file=fopen(path,"wb");
        if(file&&fwrite(magic,1,sizeof(magic),file)==sizeof(magic)&&
           fwrite(&version,sizeof(version),1,file)==1&&
           fwrite(&flags,sizeof(flags),1,file)==1&&
           fwrite(&nodes,sizeof(nodes),1,file)==1&&
           fwrite(&edges,sizeof(edges),1,file)==1&&
           fwrite(&cycles,sizeof(cycles),1,file)==1&&
           fwrite(id,1,sizeof(id),file)==sizeof(id)&&
           fwrite(&local,sizeof(local),1,file)==1&&
           fwrite(&rep,sizeof(rep),1,file)==1&&
           fwrite(&reference,sizeof(reference),1,file)==1&&
           fwrite(&step,sizeof(step),1,file)==1&&
           fwrite(shifts,sizeof(shifts[0]),WWL_CHART_LAYOUT_V3_KNOTS,file)==
               WWL_CHART_LAYOUT_V3_KNOTS&&
           fwrite(&turn,sizeof(turn),1,file)==1&&
           fwrite(&component_id,sizeof(component_id),1,file)==1)rc=0;
        if(file&&fclose(file)!=0)rc=-1;
        if(rc==0)rc=wwl_chart_layout_load(path,&layout);
        if(rc==0)rc=wwl_chart_layout_map(
            &layout,&source,1,anchor,chart,1,&correction,&mapped,
            &minimum,&maximum);
        if(rc==0)rc=wwl_chart_layout_apply(
            1,test_verts,uv,chart,1,correction);
        ok=ok&&rc==0&&mapped==1&&layout.version==3&&
           layout.relation_edges==8&&layout.cycle_edges==2&&
           fabs(minimum)<1e-12&&fabs(maximum-8.0)<1e-12&&
           fabs(uv[0]-9.5)<1e-12;
        fprintf(stderr,
            "[selftest] chart graph piecewise-V layout: nodes=%zu -> %s\n",
            mapped,ok?"ok":"FAIL");
        remove(path);free(correction);wwl_chart_layout_dispose(&layout);
    }
    {
        const char *path="welded_winding_lift_chart_layout_v4_selftest.bin";
        const char magic[8]={'G','W','L','A','Y','T','1','\n'};
        uint32_t version=4,flags=0;uint64_t nodes=2,edges=1,cycles=0;
        char ids[2][48]={{0}},source_ids[2][48]={{0}};
        int32_t local[2]={0,0},rep[2]={0,0},turn[2]={0,0};
        int32_t component_id[2]={0,0};
        double reference[2]={0,0},step[2]={16,16};
        double shifts[2][WWL_CHART_LAYOUT_V3_KNOTS]={{0}};
        uint32_t edge_a=0,edge_b=1;double edge_support=23.0;
        size_t cube_voff[3]={0,1,2};int32_t anchor[2]={0,1};
        int32_t chart[2]={1,0};
        PieceSet source;WwlChartLayout layout;
        WwlChartLayoutCorrection *correction=NULL;
        int32_t *mapped_a=NULL,*mapped_b=NULL;double *mapped_support=NULL;
        size_t mapped=0;FILE *file=NULL;int rc=-1;
        snprintf(ids[0],sizeof(ids[0]),"cube-v4-a");
        snprintf(ids[1],sizeof(ids[1]),"cube-v4-b");
        memcpy(source_ids,ids,sizeof(ids));memset(&source,0,sizeof(source));
        memset(&layout,0,sizeof(layout));source.ids=source_ids;
        source.cube_voff=cube_voff;source.n_cubes=2;source.nv=2;
        file=fopen(path,"wb");
        if(file&&fwrite(magic,1,sizeof(magic),file)==sizeof(magic)&&
           fwrite(&version,sizeof(version),1,file)==1&&
           fwrite(&flags,sizeof(flags),1,file)==1&&
           fwrite(&nodes,sizeof(nodes),1,file)==1&&
           fwrite(&edges,sizeof(edges),1,file)==1&&
           fwrite(&cycles,sizeof(cycles),1,file)==1){
            rc=0;
            for(size_t i=0;i<2&&rc==0;i++)
                if(fwrite(ids[i],1,sizeof(ids[i]),file)!=sizeof(ids[i])||
                   fwrite(&local[i],sizeof(local[i]),1,file)!=1||
                   fwrite(&rep[i],sizeof(rep[i]),1,file)!=1||
                   fwrite(&reference[i],sizeof(reference[i]),1,file)!=1||
                   fwrite(&step[i],sizeof(step[i]),1,file)!=1||
                   fwrite(shifts[i],sizeof(shifts[i][0]),
                          WWL_CHART_LAYOUT_V3_KNOTS,file)!=
                       WWL_CHART_LAYOUT_V3_KNOTS||
                   fwrite(&turn[i],sizeof(turn[i]),1,file)!=1||
                   fwrite(&component_id[i],sizeof(component_id[i]),1,file)!=1)
                    rc=-1;
        }
        if(rc==0&&
           (fwrite(&edge_a,sizeof(edge_a),1,file)!=1||
            fwrite(&edge_b,sizeof(edge_b),1,file)!=1||
            fwrite(&edge_support,sizeof(edge_support),1,file)!=1))rc=-1;
        if(file&&fclose(file)!=0)rc=-1;
        if(rc==0)rc=wwl_chart_layout_load(path,&layout);
        if(rc==0)rc=wwl_chart_layout_map(
            &layout,&source,2,anchor,chart,2,&correction,&mapped,NULL,NULL);
        if(rc==0)rc=wwl_chart_layout_map_relations(
            &layout,correction,2,&mapped_a,&mapped_b,&mapped_support);
        ok=ok&&rc==0&&mapped==2&&layout.version==4&&layout.edge!=NULL&&
           mapped_a[0]==0&&mapped_b[0]==1&&fabs(mapped_support[0]-23.0)<1e-12;
        fprintf(stderr,
            "[selftest] chart graph direct-edge handoff: edge=%d-%d "
            "support=%.3f -> %s\n",mapped_a?mapped_a[0]:-1,
            mapped_b?mapped_b[0]:-1,mapped_support?mapped_support[0]:-1.0,
            ok?"ok":"FAIL");
        remove(path);free(correction);free(mapped_a);free(mapped_b);
        free(mapped_support);wwl_chart_layout_dispose(&layout);
    }
    fprintf(stderr, "[selftest] welded winding lift primitives/provenance -> %s\n",
             ok ? "ok" : "FAIL");
    return ok ? 0 : 1;
}

int main(int argc, char **argv)
{
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return wwl_selftest();
    if (argc < 6) { wwl_usage(argv[0]); return 2; }

    const char *placed_dir = argv[1];
    const char *welded_path = argv[2];
    const char *output_path = argv[3];
    const char *lineage_path = NULL;
    const char *chart_layout_path = NULL;
    const char *bridge_edge_audit_path = NULL;
    size_t source_faces = 0, sweeps = 160;
    double gap = 20.0;
    int gap_set = 0;
    int pack = 1, source_chart_diagnostic = 0, final_source_only = 0;
    int fit_ribbon_output = 0;
    int placed_chart_diagnostic = 0, raw_uv = 0;
    int carry_registered = 0;
    int preserve_chart_gauges = 0, assemble_bridge_gauges = 0;
    int scaffold_oriented_bridge_gauges = 0;
    int registered_chart_gauges = 0, topology_winding_gauges = 0;
    int topology_winding_scaffold = 0, pin_v_to_z = 0;
    int pack_cover_intervals = 0;
    /* The authoritative output is always the .vmesh companion of <out.obj>;
     * the text OBJ itself is a viewer dump (57s of fprintf on the 21x21x21)
     * and only written on request. */
    int dump_obj = 0;
    /* --atlas-v-rebase: drop each atlas column run to its own v minimum at
     * emission (atlas height = tallest run; the 21x21x21 atlas was 80%
     * reserved-but-empty rows).  <out>_runs.json records the offsets. */
    int atlas_v_rebase = 0;
    /* --v-strips N: axial strip height in slice planes for the parallel
     * fitted-grid solve.  Default/0 = auto (~2 strips per thread; a fixed
     * 64-plane strip left slab-height grids on 5 of 32 threads); N > 0
     * pins the height (64 = one 128-vox cube of 2-vox planes). */
    int v_strips = -1;
    for (int i = 4; i < argc; i++) {
        if (strcmp(argv[i], "--source-faces") == 0 && i + 1 < argc) {
            if (wwl_parse_size(argv[++i], &source_faces) != 0) return 2;
        } else if (strcmp(argv[i], "--sweeps") == 0 && i + 1 < argc) {
            if (wwl_parse_size(argv[++i], &sweeps) != 0) return 2;
        } else if (strcmp(argv[i], "--gap") == 0 && i + 1 < argc) {
            if (wwl_parse_double(argv[++i], &gap) != 0 || gap < 0.0) return 2;
            gap_set = 1;
        } else if (strcmp(argv[i], "--no-pack") == 0) {
            pack = 0;
        } else if (strcmp(argv[i], "--dump-obj") == 0) {
            dump_obj = 1;
        } else if (strcmp(argv[i], "--atlas-v-rebase") == 0) {
            atlas_v_rebase = 1;
        } else if (strcmp(argv[i], "--v-strips") == 0 && i + 1 < argc) {
            size_t vs = 0;
            if (wwl_parse_size(argv[++i], &vs) != 0) return 2;
            v_strips = vs == 0 ? -1 : (int)vs;   /* 0 = auto */
        } else if (strcmp(argv[i], "--pin-v-to-z") == 0) {
            pin_v_to_z = 1;
        } else if (strcmp(argv[i], "--fit-ribbon") == 0) {
            fit_ribbon_output = 1;
        } else if (strcmp(argv[i], "--source-chart-diagnostic") == 0) {
            source_chart_diagnostic = 1;
        } else if (strcmp(argv[i], "--final-source-only") == 0) {
            final_source_only = 1;
        } else if (strcmp(argv[i], "--placed-chart-diagnostic") == 0) {
            placed_chart_diagnostic = 1;
        } else if (strcmp(argv[i], "--raw-uv") == 0) {
            raw_uv = 1;
        } else if (strcmp(argv[i], "--carry-registered") == 0) {
            carry_registered = 1;
        } else if (strcmp(argv[i], "--preserve-chart-gauges") == 0) {
            preserve_chart_gauges = 1;
        } else if (strcmp(argv[i], "--assemble-bridge-gauges") == 0) {
            assemble_bridge_gauges = 1;
        } else if (strcmp(argv[i],
                          "--scaffold-oriented-bridge-gauges") == 0) {
            scaffold_oriented_bridge_gauges = 1;
        } else if (strcmp(argv[i], "--registered-chart-gauges") == 0) {
            registered_chart_gauges = 1;
        } else if (strcmp(argv[i], "--topology-winding-gauges") == 0) {
            topology_winding_gauges = 1;
        } else if (strcmp(argv[i], "--topology-winding-scaffold") == 0) {
            topology_winding_scaffold = 1;
        } else if (strcmp(argv[i], "--lineage") == 0 && i + 1 < argc) {
            lineage_path = argv[++i];
        } else if (strcmp(argv[i], "--chart-layout") == 0 && i + 1 < argc) {
            chart_layout_path = argv[++i];
        } else if (strcmp(argv[i], "--bridge-edge-audit") == 0 &&
                   i + 1 < argc) {
            bridge_edge_audit_path = argv[++i];
        } else {
            fprintf(stderr, "welded_winding_lift: unknown option %s\n", argv[i]);
            return 2;
        }
    }
    if (source_faces == 0 && lineage_path == NULL) {
        fprintf(stderr,
            "welded_winding_lift: --source-faces or --lineage is required\n");
        return 2;
    }
    int gauge_modes = carry_registered + preserve_chart_gauges +
                      assemble_bridge_gauges +
                      scaffold_oriented_bridge_gauges +
                      registered_chart_gauges + topology_winding_gauges +
                      topology_winding_scaffold;
    if (gauge_modes > 1 ||
        ((assemble_bridge_gauges || scaffold_oriented_bridge_gauges ||
          registered_chart_gauges ||
          topology_winding_gauges || topology_winding_scaffold) && !raw_uv)) {
        fprintf(stderr,
            "welded_winding_lift: bridge/registered/topology winding gauge "
            "assembly requires --raw-uv and gauge modes are mutually "
            "exclusive\n");
        return 2;
    }
    if (carry_registered && raw_uv) {
        fprintf(stderr,
            "welded_winding_lift: --carry-registered carries the "
            "RE-REGISTERED field; it is incompatible with --raw-uv\n");
        return 2;
    }
    if(chart_layout_path!=NULL&&!carry_registered){
        fprintf(stderr,
            "welded_winding_lift: --chart-layout currently requires "
            "--carry-registered\n");
        return 2;
    }
    /* The registered atlas is contiguous; the default 20-vox island gap was
     * tuned for re-derived gauges and paints 10-px black bars between every
     * island at du=2. */
    if (carry_registered && !gap_set)
        gap = 4.0;
    if (bridge_edge_audit_path != NULL && !assemble_bridge_gauges &&
        !scaffold_oriented_bridge_gauges) {
        fprintf(stderr,
            "welded_winding_lift: --bridge-edge-audit requires "
            "--raw-uv --assemble-bridge-gauges\n");
        return 2;
    }
    if (fit_ribbon_output && carry_registered && !pin_v_to_z) {
        fprintf(stderr,
            "welded_winding_lift: --fit-ribbon --carry-registered "
            "requires --pin-v-to-z\n");
        return 2;
    }
    if (fit_ribbon_output && !carry_registered &&
        (!raw_uv ||
         (!assemble_bridge_gauges && !topology_winding_gauges &&
          !topology_winding_scaffold && gauge_modes != 0) ||
         !pin_v_to_z)) {
        fprintf(stderr,
            "welded_winding_lift: --fit-ribbon requires --carry-registered "
            "--pin-v-to-z, or the legacy --raw-uv "
            "the default topology lift, or --assemble-bridge-gauges/"
            "--topology-winding-gauges/--topology-winding-scaffold "
            "--pin-v-to-z\n");
        return 2;
    }

    Arena_T arena = Arena_new();
    PieceSet source;
    ScaffoldCalib calib;
    WwlGraph bridge_graph;
    float *verts = NULL;
    int32_t *faces = NULL;
    size_t nv = 0, nf = 0;
    int32_t *anchor_source = NULL, *chart = NULL, *component = NULL;
    size_t *chart_faces = NULL, *component_faces = NULL;
    size_t nchart = 0, ncomponent = 0;
    WwlGraph source_graph, full_graph;
    WwlLineage lineage;
    WwlChartLayout chart_layout;
    double *uv = NULL;
    float *registered_uv = NULL;
    float *registered_phi = NULL;
    float *fit_phi = NULL;
    double *fit_atlas_shift = NULL;
    double *fit_reference_uv = NULL;
    double *fit_reference_u = NULL;
    int32_t *fit_reference_island = NULL;
    size_t *fit_reference_island_faces = NULL;
    /* Relation components PLUS one island per unmatched chart.  Every
       consumer of fit_reference_island must use this, not
       chart_layout.relation_components, or the orphan islands sit past the
       end of the array. */
    size_t fit_reference_island_count = 0;
    int32_t *fit_relation_a = NULL;
    int32_t *fit_relation_b = NULL;
    double *fit_relation_support = NULL;
    WwlChartLayoutCorrection *chart_layout_correction = NULL;
    int result = 1;
    memset(&source, 0, sizeof(source));
    memset(&source_graph, 0, sizeof(source_graph));
    memset(&full_graph, 0, sizeof(full_graph));
    memset(&bridge_graph, 0, sizeof(bridge_graph));
    memset(&lineage, 0, sizeof(lineage));
    memset(&chart_layout,0,sizeof(chart_layout));

    if (PieceSet_build(arena, placed_dir, &source) != 0) {
        fprintf(stderr, "welded_winding_lift: cannot load placed dir %s\n",
                placed_dir);
        goto cleanup;
    }
    if(chart_layout_path!=NULL&&
       wwl_chart_layout_load(chart_layout_path,&chart_layout)!=0){
        fprintf(stderr,"welded_winding_lift: cannot load chart layout %s\n",
                chart_layout_path);
        goto cleanup;
    }
    if (bridge_edge_audit_path != NULL) {
        registered_phi = (float *)malloc(
            (source.nv ? source.nv : 1) * sizeof(*registered_phi));
        if (registered_phi == NULL) goto cleanup;
        memcpy(registered_phi, source.phi,
               source.nv * sizeof(*registered_phi));
    }
    if (registered_chart_gauges) {
        registered_uv = (float *)malloc(
            (source.nv ? source.nv : 1) * 2 * sizeof(*registered_uv));
        if (registered_uv == NULL) goto cleanup;
        memcpy(registered_uv, source.uv,
               source.nv * 2 * sizeof(*registered_uv));
    }
    if (raw_uv && wwl_load_raw_uv(placed_dir, &source) != 0) {
        fprintf(stderr,
            "welded_winding_lift: cannot load complete raw UV field from %s\n",
            placed_dir);
        goto cleanup;
    }
    if (raw_uv)
        fprintf(stderr,
            "welded_winding_lift: using pre-reregister _uvphi_raw field\n");
    if (placed_chart_diagnostic) {
        if (wwl_write_placed_chart_diagnostic(output_path, &source, gap) != 0) {
            fprintf(stderr,
                "welded_winding_lift: placed chart diagnostic failed\n");
            goto cleanup;
        }
        fprintf(stderr,
            "welded_winding_lift: wrote placed chart diagnostic %s\n",
            output_path);
        result = 0;
        goto cleanup;
    }
    Scaffold_calib_default(&calib);
    if (Scaffold_read_calib(placed_dir, &calib) != 0)
        fprintf(stderr,
                "welded_winding_lift: using default scaffold calibration\n");
    {
        char binary_path[2048];
        MeshBinData mesh;
        if (MeshBin_companion_path(welded_path, binary_path,
                                   sizeof binary_path) != 0 ||
            MeshBin_read_arena(arena, binary_path, &mesh) != 0 ||
            mesh.nv == 0 || mesh.nf == 0) {
            fprintf(stderr,
                "welded_winding_lift: cannot load authoritative welded "
                "container %s (OBJ fallback is disabled)\n",
                MeshBin_companion_path(welded_path, binary_path,
                                       sizeof binary_path) == 0
                    ? binary_path : welded_path);
            goto cleanup;
        }
        verts = mesh.verts; nv = mesh.nv;
        faces = mesh.faces; nf = mesh.nf;
        fprintf(stderr,
                "welded_winding_lift: loaded %s (%zu vertices, %zu faces)\n",
                binary_path, nv, nf);
    }
    if (lineage_path != NULL) {
        size_t lineage_midpoints = 0;
        if (wwl_lineage_load(lineage_path, nv, verts, &lineage,
                             &lineage_midpoints) != 0) {
            fprintf(stderr,
                "welded_winding_lift: cannot validate lineage %s against %s\n",
                lineage_path, welded_path);
            goto cleanup;
        }
        if (source_faces != 0 && source_faces != lineage.source_faces) {
            fprintf(stderr,
                "welded_winding_lift: source-face mismatch: option=%zu "
                "lineage=%zu\n", source_faces, lineage.source_faces);
            goto cleanup;
        }
        source_faces = lineage.source_faces;
        fprintf(stderr,
            "welded_winding_lift: lineage v%u geometry PASS: %zu exact "
            "source anchors, %zu exact midpoint vertices, %zu source faces\n",
            (unsigned)lineage.version, lineage.source_anchors,
            lineage_midpoints, source_faces);
    }
    if (source_faces == 0 || source_faces > nf) {
        fprintf(stderr,
            "welded_winding_lift: invalid source-face prefix %zu/%zu\n",
            source_faces, nf);
        goto cleanup;
    }
    anchor_source = (int32_t *)malloc(nv * sizeof(*anchor_source));
    uv = (double *)malloc(nv * 2 * sizeof(*uv));
    if (anchor_source == NULL || uv == NULL) goto cleanup;

    size_t matched = 0, missing = 0, ambiguous = 0, ordered = 0;
    (void)wwl_source_order_diagnostic(&source, verts, nv, &ordered);
    int source_map_rc =
        lineage.version >= 2 && lineage.source_anchors > 0
            ? wwl_map_source_lineage(&source, verts, nv, &lineage,
                                     anchor_source, &matched, &missing,
                                     &ambiguous)
            : wwl_map_source_exact(&source, verts, nv, anchor_source,
                                   &matched, &missing, &ambiguous);
    if (source_map_rc != 0) {
        fprintf(stderr,
            "welded_winding_lift: exact source join failed "
            "(%zu/%zu matches, %zu ambiguous)\n",
            matched, source.nv, ambiguous);
        goto cleanup;
    }
    fprintf(stderr,
        "welded_winding_lift: source provenance PASS: %zu exact %s anchors, "
        "%zu absent; ordered diagnostic %zu\n",
        matched, lineage.version >= 2 ? "cube/local" : "position",
        missing, ordered);

    if (wwl_label_components(nv, faces, source_faces,
                             &chart, &nchart, &chart_faces) != 0 ||
        wwl_label_components(nv, faces, nf,
                             &component, &ncomponent, &component_faces) != 0) {
        fprintf(stderr, "welded_winding_lift: component labeling failed\n");
        goto cleanup;
    }
    fprintf(stderr,
        "welded_winding_lift: %zu source charts -> %zu welded components\n",
        nchart, ncomponent);
    if(chart_layout_path!=NULL){
        size_t mapped=0,moved_turns=0;double minimum=0.0,maximum=0.0;
        int32_t minimum_turn=INT32_MAX,maximum_turn=INT32_MIN;
        if(wwl_chart_layout_map(
                &chart_layout,&source,nv,anchor_source,chart,nchart,
                &chart_layout_correction,&mapped,&minimum,&maximum)!=0){
            fprintf(stderr,
                "welded_winding_lift: chart layout provenance map failed "
                "(%zu graph node(s), %zu source chart(s))\n",
                chart_layout.nnode,nchart);
            goto cleanup;
        }
        for(size_t i=0;i<chart_layout.nnode;i++){
            int32_t turn=chart_layout.node[i].winding_turn;
            if(turn<minimum_turn)minimum_turn=turn;
            if(turn>maximum_turn)maximum_turn=turn;
            if(turn!=0)moved_turns++;
        }
        fprintf(stderr,
            "welded_winding_lift: chart graph provenance PASS: %zu node(s), "
            "%zu relation edge(s), %zu cycle edge(s), %zu relation "
            "component(s), shift=[%.6g,%.6g], discrete turns=[%d,%d] "
            "on %zu node(s)\n",
            mapped,chart_layout.relation_edges,chart_layout.cycle_edges,
            chart_layout.relation_components,minimum,maximum,
            minimum_turn,maximum_turn,moved_turns);
        size_t island_spanning_components=0,island_spanning_bridges=0;
        fit_reference_island_count=chart_layout.relation_components+
            (nchart>mapped?nchart-mapped:0);
        if(wwl_chart_layout_reference_islands(
                nv,chart,nchart,chart_layout_correction,
                component,ncomponent,faces,nf,source_faces,
                fit_reference_island_count,
                &island_spanning_components,&island_spanning_bridges,
                &fit_reference_island,
                &fit_reference_island_faces)!=0){
            fprintf(stderr,
                "welded_winding_lift: graph relation components conflict "
                "with welded topology\n");
            goto cleanup;
        }
        if(chart_layout.version>=4&&
           wwl_chart_layout_map_relations(
               &chart_layout,chart_layout_correction,nchart,
               &fit_relation_a,&fit_relation_b,
               &fit_relation_support)!=0){
            fprintf(stderr,
                "welded_winding_lift: direct chart relation map failed\n");
            goto cleanup;
        }
        fprintf(stderr,
            "welded_winding_lift: relation-island topology PASS: %zu "
            "graph component(s) across %zu welded component(s); "
            "%zu component(s) and %zu bridge face(s) span deferred "
            "junctions (excluded from the fit)\n",
            chart_layout.relation_components,ncomponent,
            island_spanning_components,island_spanning_bridges);
    }
    /* full_graph (~2.4 GB CSR on the 21x21x21, 302M scattered writes) is
     * consumed only by wwl_lift_components -- the re-derivation paths.  The
     * carry/preserve/registered modes never read it; carry's anchorless
     * seeding needs only the BRIDGE-face edges (charts touch other charts
     * exclusively through bridge faces), a ~50x smaller graph. */
    int need_full_graph = topology_winding_scaffold || topology_winding_gauges
                        || assemble_bridge_gauges
                        || scaffold_oriented_bridge_gauges || gauge_modes == 0;
    if (wwl_graph_build(nv, faces, source_faces, &source_graph) != 0 ||
        (need_full_graph &&
         wwl_graph_build(nv, faces, nf, &full_graph) != 0) ||
        (carry_registered &&
         wwl_graph_build(nv, faces + source_faces * 3,
                         nf - source_faces, &bridge_graph) != 0)) {
        fprintf(stderr, "welded_winding_lift: adjacency build failed\n");
        goto cleanup;
    }
    if (carry_registered) {
        /* Winsorize the registered anchors (the 0.02% u-tail), then seed
         * charts that carry no provenance at all from anchored neighbours
         * across bridge faces, so the harmonic fill has support everywhere. */
        size_t demoted = wwl_carry_winsorize_anchors(
            nv, chart, nchart, anchor_source, &source, 1000.0);
        size_t seeded = 0, dry = 0;
        wwl_carry_seed_anchorless(nv, &bridge_graph, chart, nchart,
                                  anchor_source, &seeded, &dry);
        fprintf(stderr,
            "welded_winding_lift: carry prep: %zu tail anchor(s) demoted, "
            "%zu anchorless chart(s) seeded from neighbours, %zu dry\n",
            demoted, seeded, dry);
        if (fit_ribbon_output) {
            fit_phi = (float *)calloc(nv ? nv : 1, sizeof(*fit_phi));
            if (fit_phi == NULL) goto cleanup;
        }
    }
    size_t unknown = 0, midpoint_fixed = 0, residual = 0;
    double last_delta = 0.0;
    if (wwl_interpolate_chart_uv(verts, nv, &source_graph, chart,
                                 anchor_source, &source,
                                 lineage_path != NULL ? &lineage : NULL,
                                 sweeps, uv,
                                 carry_registered ? fit_phi : NULL,
                                 &unknown, &midpoint_fixed,
                                 &residual, &last_delta) != 0) {
        fprintf(stderr, "welded_winding_lift: chart interpolation failed"
                "%s\n", carry_registered
                ? " (a chart still has no carried anchor -- see 'dry' above)"
                : "");
        goto cleanup;
    }
    fprintf(stderr,
        "welded_winding_lift: retained %zu raw anchors; fixed %zu exact "
        "midpoint descendants; harmonic residue %zu/%zu "
        "(%zu sweeps, residual %.6g)\n",
        matched, midpoint_fixed, residual, unknown, sweeps, last_delta);
    if(chart_layout_correction!=NULL){
        if(wwl_chart_layout_apply(
                nv,verts,uv,chart,nchart,chart_layout_correction)!=0){
            fprintf(stderr,"welded_winding_lift: chart layout apply failed\n");
            goto cleanup;
        }
        fprintf(stderr,
            "welded_winding_lift: applied graph-synchronized U %s to %zu "
            "source chart(s)\n",chart_layout.version>=2
            ? "V-dependent correction" : "translation",nchart);
    }

    if (source_chart_diagnostic) {
        if (pack &&
            wwl_pack_components(nv, uv, chart, nchart, chart_faces, gap) != 0) {
            fprintf(stderr,
                "welded_winding_lift: source-chart diagnostic packing failed\n");
            goto cleanup;
        }
        size_t positive = 0, negative = 0, degenerate = 0;
        wwl_face_orientation_stats(uv, faces, source_faces,
                                   &positive, &negative, &degenerate);
        fprintf(stderr,
            "welded_winding_lift: source-chart UV signs +%zu -%zu zero=%zu\n",
            positive, negative, degenerate);
        if (ves_ensure_parent_dir(output_path) != 0 ||
            ObjIO_write_uv_double(output_path, verts, nv, faces,
                                  source_faces, uv) != 0) {
            fprintf(stderr, "welded_winding_lift: cannot write %s\n",
                    output_path);
            goto cleanup;
        }
        fprintf(stderr,
            "welded_winding_lift: wrote source-chart diagnostic %s\n",
            output_path);
        result = 0;
        goto cleanup;
    }

    if (topology_winding_scaffold) {
        int polar_sign = 0;
        double polar_offset = 0.0, cycle_error = 0.0;
        WwlChartGaugeStats unused_chart_gauge;
        memset(&unused_chart_gauge, 0, sizeof unused_chart_gauge);
        if (wwl_lift_components(
                verts, nv, faces, source_faces, nf, &full_graph,
                component, ncomponent, chart, nchart, anchor_source, &source,
                &calib, uv, 2, &polar_sign, &polar_offset, &cycle_error,
                &unused_chart_gauge, NULL, NULL, NULL) != 0) {
            fprintf(stderr,
                "welded_winding_lift: direct topology winding scaffold failed\n");
            goto cleanup;
        }
        fprintf(stderr,
            "welded_winding_lift: direct topology winding scaffold PASS: "
            "sign=%d offset=%.9g cycle_residual=%.6g; "
            "u=F(topology phase), v=axial\n",
            polar_sign, polar_offset, cycle_error);
    } else if (topology_winding_gauges) {
        int polar_sign = 0;
        double polar_offset = 0.0, cycle_error = 0.0;
        WwlChartGaugeStats unused_chart_gauge;
        WwlTopologyGaugeStats topology_gauge;
        WwlBridgeGaugeStats bridge_gauge;
        memset(&unused_chart_gauge, 0, sizeof unused_chart_gauge);
        memset(&topology_gauge, 0, sizeof topology_gauge);
        memset(&bridge_gauge, 0, sizeof bridge_gauge);
        if (wwl_lift_components(
                verts, nv, faces, source_faces, nf, &full_graph,
                component, ncomponent, chart, nchart, anchor_source, &source,
                &calib, uv, 1, &polar_sign, &polar_offset, &cycle_error,
                &unused_chart_gauge, &topology_gauge, &bridge_gauge, NULL) != 0) {
            fprintf(stderr,
                "welded_winding_lift: topology winding gauge assembly failed\n");
            goto cleanup;
        }
        fprintf(stderr,
            "welded_winding_lift: topology winding gauges PASS: sign=%d "
            "offset=%.9g cycle_residual=%.6g; %zu chart turn(s) [%d,%d], "
            "%zu anchor(s), turn residual median/max=%.6g/%.6g\n",
            polar_sign, polar_offset, cycle_error, topology_gauge.charts,
            topology_gauge.minimum_turn, topology_gauge.maximum_turn,
            topology_gauge.anchors, topology_gauge.residual_median_turns,
            topology_gauge.residual_maximum_turns);
        fprintf(stderr,
            "welded_winding_lift: certified bridge translation diagnostic "
            "PASS (raw chart origins preserved): "
            "%zu observation(s), %zu chart edge(s), %zu forest(s), "
            "cycles=%zu; residual median/max=%.6g/%.6g vox, cycle=%.6g; "
            "orientation-normalized=%zu\n",
            bridge_gauge.observations, bridge_gauge.gauge_edges,
            bridge_gauge.gauge_components, bridge_gauge.gauge_cycles,
            bridge_gauge.observation_residual_median,
            bridge_gauge.observation_residual_maximum,
            bridge_gauge.cycle_residual_maximum,
            bridge_gauge.orientation_normalized_charts);
    } else if (registered_chart_gauges) {
        WwlChartGaugeStats gauge;
        if (wwl_apply_registered_chart_gauges(
                nv, source.nv, chart, nchart, anchor_source, registered_uv, uv,
                &gauge) != 0) {
            fprintf(stderr,
                "welded_winding_lift: registered chart gauge transfer failed\n");
            goto cleanup;
        }
        fprintf(stderr,
            "welded_winding_lift: registered chart gauges PASS: "
            "%zu anchor observation(s), %zu reflected, %zu point-gauge; "
            "residual median/max=%.6g/%.6g vox; "
            "raw chart geometry unchanged\n",
            gauge.observations, gauge.reflected_charts,
            gauge.point_gauge_charts, gauge.residual_median,
            gauge.residual_maximum);
    } else if (assemble_bridge_gauges || scaffold_oriented_bridge_gauges) {
        WwlBridgeGaugeStats gauge;
        WwlBridgeAudit audit;
        double *scaffold_uv = NULL;
        int32_t *pose_anchor = NULL;
        int scaffold_sign = 0;
        double scaffold_offset = 0.0, scaffold_cycle_error = 0.0;
        WwlChartGaugeStats unused_scaffold_gauge;
        WwlChartGaugeStats scaffold_pose_gauge;
        size_t scaffold_seeded = 0;
        memset(&audit, 0, sizeof audit);
        memset(&unused_scaffold_gauge, 0, sizeof unused_scaffold_gauge);
        memset(&scaffold_pose_gauge, 0, sizeof scaffold_pose_gauge);
        audit.path = bridge_edge_audit_path;
        audit.chart_faces = chart_faces;
        audit.component = component;
        audit.component_faces = component_faces;
        audit.ncomponent = ncomponent;
        audit.anchor_source = anchor_source;
        audit.source = &source;
        audit.calib = &calib;
        audit.registered_phi = registered_phi;
        scaffold_uv = (double *)malloc(
            (nv ? nv : 1) * 2 * sizeof(*scaffold_uv));
        if (fit_ribbon_output)
            fit_phi = (float *)malloc(
                (nv ? nv : 1) * sizeof(*fit_phi));
        if (scaffold_uv == NULL || (fit_ribbon_output && fit_phi == NULL) ||
            wwl_lift_components(
                verts, nv, faces, source_faces, nf, &full_graph,
                component, ncomponent, chart, nchart, anchor_source, &source,
                &calib, scaffold_uv, 2, &scaffold_sign, &scaffold_offset,
                &scaffold_cycle_error, &unused_scaffold_gauge, NULL, NULL,
                fit_phi) != 0 ||
            wwl_seed_invalid_chart_targets(
                nv, faces, source_faces, chart, nchart, anchor_source,
                scaffold_uv, uv, fit_ribbon_output && pin_v_to_z,
                &scaffold_seeded) != 0) {
            free(scaffold_uv);
            fprintf(stderr,
                "welded_winding_lift: invalid-chart topology scaffold failed\n");
            goto cleanup;
        }
        if (scaffold_oriented_bridge_gauges) {
            pose_anchor = (int32_t *)malloc(
                (nv ? nv : 1) * sizeof(*pose_anchor));
            if (pose_anchor == NULL) {
                free(scaffold_uv);
                goto cleanup;
            }
            for (size_t v = 0; v < nv; v++) pose_anchor[v] = 0;
            if (wwl_apply_target_chart_isometries(
                    nv, chart, nchart, pose_anchor, scaffold_uv, uv,
                    &scaffold_pose_gauge) != 0) {
                free(pose_anchor);
                free(scaffold_uv);
                fprintf(stderr,
                    "welded_winding_lift: absolute scaffold chart pose failed\n");
                goto cleanup;
            }
            free(pose_anchor);
            pose_anchor = NULL;
            fprintf(stderr,
                "welded_winding_lift: absolute scaffold chart poses PASS: "
                "%zu vertex observation(s), %zu reflected, %zu point-gauge; "
                "residual median/max=%.6g/%.6g vox\n",
                scaffold_pose_gauge.observations,
                scaffold_pose_gauge.reflected_charts,
                scaffold_pose_gauge.point_gauge_charts,
                scaffold_pose_gauge.residual_median,
                scaffold_pose_gauge.residual_maximum);
        }
        if (scaffold_seeded != 0)
            fprintf(stderr,
                "welded_winding_lift: topology scaffold seeded %zu "
                "invalid local chart(s); sign=%d offset=%.9g "
                "cycle_residual=%.6g\n",
                scaffold_seeded, scaffold_sign, scaffold_offset,
                scaffold_cycle_error);
        if (fit_ribbon_output) {
            fit_reference_uv = (double *)malloc(
                (nv ? nv : 1) * 2 * sizeof(*fit_reference_uv));
            if (fit_reference_uv == NULL) {
                free(scaffold_uv);
                goto cleanup;
            }
            memcpy(fit_reference_uv, scaffold_uv,
                   nv * 2 * sizeof(*fit_reference_uv));
        }
        if (wwl_assemble_bridge_gauges(
                verts, nv, faces, source_faces, nf, chart, nchart,
                ncomponent,
                pin_v_to_z ? 0 : (scaffold_oriented_bridge_gauges ? 0 : 1),
                1, pin_v_to_z,
                uv, &gauge,
                bridge_edge_audit_path != NULL ? &audit : NULL) != 0) {
            fprintf(stderr,
                "welded_winding_lift: bridge gauge assembly failed "
                "(faces=%zu same=%zu two=%zu three=%zu invalid=%zu; "
                "obs=%zu edges=%zu components=%zu cycles=%zu)\n",
                gauge.bridge_faces, gauge.same_chart_faces,
                gauge.two_chart_faces, gauge.three_chart_faces,
                gauge.invalid_faces, gauge.observations, gauge.gauge_edges,
                gauge.gauge_components, gauge.gauge_cycles);
            free(scaffold_uv);
            goto cleanup;
        }
        fprintf(stderr,
            "welded_winding_lift: bridge gauge assembly PASS: "
            "%zu observation(s), %zu chart edge(s), %zu forest(s), "
            "cycles=%zu; residual median/max=%.6g/%.6g vox, "
            "cycle=%.6g; orientation-normalized=%zu; "
            "translation-selected=%zu\n",
            gauge.observations, gauge.gauge_edges,
            gauge.gauge_components, gauge.gauge_cycles,
            gauge.observation_residual_median,
            gauge.observation_residual_maximum,
            gauge.cycle_residual_maximum,
            gauge.orientation_normalized_charts,
            gauge.translation_selected_edges);
        fprintf(stderr,
            "welded_winding_lift: bridge faces same/two/three/invalid="
            "%zu/%zu/%zu/%zu\n",
            gauge.same_chart_faces, gauge.two_chart_faces,
            gauge.three_chart_faces, gauge.invalid_faces);
        if (pin_v_to_z) {
            double gauge_residual_mean = 0.0, gauge_residual_maximum = 0.0;
            if (wwl_align_component_u_to_target(
                    nv, uv, scaffold_uv, component, ncomponent,
                    &gauge_residual_mean, &gauge_residual_maximum) != 0) {
                free(scaffold_uv);
                fprintf(stderr,
                    "welded_winding_lift: component cover-gauge alignment failed\n");
                goto cleanup;
            }
            fprintf(stderr,
                "welded_winding_lift: aligned %zu stitched component U gauge(s) "
                "to topology cover (residual mean/max=%.6g/%.6g vox)\n",
                ncomponent, gauge_residual_mean, gauge_residual_maximum);
            pack_cover_intervals = 1;
        }
        free(scaffold_uv);
    } else if (carry_registered) {
        /* THE point of this mode: the re-registered field (r2=1.000, matches
         * the analytic spiral arc length to 0.2%) was interpolated into uv
         * above and is carried through untouched -- no scaffold lift, no
         * bridge-gauge assembly, no mean-shift component alignment, none of
         * the re-derivation whose alignment residual reached a mean of
         * 23,815 vox on the 21x21x21 and shredded the atlas. */
        fprintf(stderr,
            "welded_winding_lift: carrying re-registered field verbatim "
            "(no gauge re-derivation)\n");
        /* interval packing preserves relative placement INSIDE an island --
         * exactly the registration being carried; per-chart packing would
         * scramble it */
        pack_cover_intervals = 1;
        if (fit_ribbon_output) {
            fit_reference_uv = (double *)malloc(
                (nv ? nv : 1) * 2 * sizeof(*fit_reference_uv));
            if (fit_reference_uv == NULL) goto cleanup;
            memcpy(fit_reference_uv, uv,
                   nv * 2 * sizeof(*fit_reference_uv));
        }
    } else if (preserve_chart_gauges) {
        fprintf(stderr,
            "welded_winding_lift: preserving original winding-chart gauges\n");
    } else {
        int polar_sign = 0;
        double polar_offset = 0.0, cycle_error = 0.0;
        WwlChartGaugeStats gauge;
        if (wwl_lift_components(
                verts, nv, faces, source_faces, nf, &full_graph,
                component, ncomponent, chart, nchart, anchor_source, &source,
                &calib, uv, 0, &polar_sign, &polar_offset,
                &cycle_error, &gauge, NULL, NULL, NULL) != 0) {
            fprintf(stderr,
                    "welded_winding_lift: topology winding lift failed\n");
            goto cleanup;
        }
        fprintf(stderr,
            "welded_winding_lift: topology winding sign=%d offset=%.9g "
            "cycle_residual=%.6g; chart isometries=%zu reflected, "
            "%zu point-gauge, residual median/max=%.6g/%.6g vox\n",
            polar_sign, polar_offset, cycle_error, gauge.reflected_charts,
            gauge.point_gauge_charts, gauge.residual_median,
            gauge.residual_maximum);
    }

    /* V is world Z in this parameterisation -- scroll_whole writes it that way
     * and rawtex_bake audits _diagverr = |vt.v - (z - zmin)|.  Keep that
     * contract exact even for generated chart vertices and float roundoff. */
    if (pin_v_to_z) {
        double vmoved = 0.0;
        size_t nmoved = 0;
        for (size_t v = 0; v < nv; v++) {
            double want = (double)verts[v * 3 + 0];
            double delta = fabs(uv[v * 2 + 1] - want);
            if (delta > 1.0e-6) { nmoved++; if (delta > vmoved) vmoved = delta; }
            uv[v * 2 + 1] = want;
        }
        fprintf(stderr,
            "welded_winding_lift: pinned v to world z on %zu/%zu vertices "
            "(max correction %.6g vox)\n", nmoved, nv, vmoved);
    }

    if (fit_ribbon_output) {
        if (fit_reference_uv == NULL) {
            fit_reference_uv = (double *)malloc(
                (nv ? nv : 1) * 2 * sizeof(*fit_reference_uv));
            if (fit_reference_uv == NULL) goto cleanup;
            memcpy(fit_reference_uv, uv,
                   nv * 2 * sizeof(*fit_reference_uv));
        }
        RibbonOpts ribbon_opts;
        RibbonResult ribbon;
        RibbonWriteStats written;
        fit_atlas_shift = (double *)malloc(
            (nv ? nv : 1) * sizeof(*fit_atlas_shift));
        fit_reference_u = (double *)malloc(
            (nv ? nv : 1) * sizeof(*fit_reference_u));
        if (fit_atlas_shift == NULL || fit_reference_u == NULL ||
            fit_reference_uv == NULL)
            goto cleanup;
        for (size_t v = 0; v < nv; v++)
            fit_atlas_shift[v] = -fit_reference_uv[v * 2];
        if (pack && pack_cover_intervals) {
            size_t atlas_islands = 0;
            double cover_span = 0.0, atlas_span = 0.0;
            const int32_t *pack_component=fit_reference_island!=NULL
                                         ? fit_reference_island:component;
            size_t pack_count=fit_reference_island!=NULL
                            ? fit_reference_island_count:ncomponent;
            const size_t *pack_faces=fit_reference_island_faces!=NULL
                                    ? fit_reference_island_faces:component_faces;
            if (wwl_pack_component_intervals(
                    nv, fit_reference_uv, pack_component, pack_count,
                    pack_faces, gap,
                    &atlas_islands, &cover_span, &atlas_span) != 0) {
                fprintf(stderr,
                    "welded_winding_lift: pre-fit cover packing failed\n");
                goto cleanup;
            }
            fprintf(stderr,
                "welded_winding_lift: pre-fit topology-gauge pack %zu "
                "component(s) -> %zu interval island(s), span %.6g -> "
                "%.6g vox (gap %.3g)\n",
                pack_count, atlas_islands, cover_span, atlas_span, gap);
        }
        for (size_t v = 0; v < nv; v++) {
            fit_atlas_shift[v] += fit_reference_uv[v * 2];
            fit_reference_u[v] = fit_reference_uv[v * 2];
        }
        RibbonOpts_default(&ribbon_opts);
        memcpy(ribbon_opts.axis_point, calib.axis_point,
               sizeof(ribbon_opts.axis_point));
        memcpy(ribbon_opts.axis_dir, calib.axis_dir,
               sizeof(ribbon_opts.axis_dir));
        ribbon_opts.slice_h = 2.0f;
        ribbon_opts.sample_h = 2.0f;
        ribbon_opts.grid_u = 2.0f;
        ribbon_opts.wrap_spacing = (float)calib.pitch;
        ribbon_opts.reference_u = fit_reference_u;
        ribbon_opts.reference_phi = fit_phi;
        ribbon_opts.reference_atlas_shift = fit_atlas_shift;
        ribbon_opts.reference_island = fit_reference_island!=NULL
                                     ? fit_reference_island:component;
        ribbon_opts.reference_island_count = fit_reference_island!=NULL
                                           ? fit_reference_island_count
                                           : ncomponent;
        ribbon_opts.reference_chart = chart_layout.version>=4 ? chart : NULL;
        ribbon_opts.reference_chart_count = chart_layout.version>=4
                                          ? nchart : 0;
        ribbon_opts.reference_relation_a = fit_relation_a;
        ribbon_opts.reference_relation_b = fit_relation_b;
        ribbon_opts.reference_relation_support = fit_relation_support;
        ribbon_opts.reference_relation_count = chart_layout.version>=4
                                             ? chart_layout.relation_edges : 0;
        ribbon_opts.component_global = 1;
        /* Carried graph registration supplies one absolute observation for
         * every chain gauge; physical cross-row matches refine, never replace,
         * that globally synchronized frame. */
        ribbon_opts.reference_anchor_gauges = carry_registered;
        ribbon_opts.v_strip_planes = v_strips;
        ribbon_opts.radial_bridge_cut = -1;
        memset(&ribbon, 0, sizeof ribbon);
        memset(&written, 0, sizeof written);
        fprintf(stderr,
            "welded_winding_lift: fitting regular ribbon directly from "
            "optimized material U (metric solve skipped)\n");
        size_t ribbon_faces = chart_layout.version>=4 ? source_faces : nf;
        if (Ribbon_run(arena, verts, nv, faces, ribbon_faces,
                       &ribbon_opts, &ribbon) != 0 ||
            ves_ensure_parent_dir(output_path) != 0 ||
            Ribbon_write_obj(output_path, &ribbon, dump_obj,
                             atlas_v_rebase, &written) != 0) {
            fprintf(stderr,
                "welded_winding_lift: fast fitted-ribbon output failed\n");
            goto cleanup;
        }
        fprintf(stderr,
            "welded_winding_lift: wrote fitted ribbon %s: "
            "%zu vertices, %zu faces, %zu atlas columns in %zu run(s); "
            "claims=%zu replaced=%zu\n",
            output_path, written.vertices, written.faces,
            written.atlas_columns, written.atlas_runs,
            ribbon.grid_claim_conflicts, ribbon.grid_claim_replaced);
        fprintf(stderr,
            "welded_winding_lift: ownership by construction: "
            "junctions=%zu ambiguous=%zu fuse(same=%zu alternate=%zu "
            "violations=%zu max=%.2f) transitions(gap=%zu certified=%zu "
            "uncertified=%zu)\n",
            ribbon.grid_merge_junctions, ribbon.grid_merge_ambiguous,
            ribbon.grid_fuse_same_logical, ribbon.grid_fuse_alternate,
            ribbon.grid_fuse_violations, ribbon.grid_fuse_violation_max,
            ribbon.grid_transition_gap, ribbon.grid_transition_certified,
            ribbon.grid_transition_uncertified);
        result = 0;
        goto cleanup;
    }

    if (pack) {
        const int32_t *pack_component=fit_reference_island!=NULL
                                     ? fit_reference_island:component;
        size_t pack_count=fit_reference_island!=NULL
                        ? fit_reference_island_count:ncomponent;
        const size_t *pack_faces=fit_reference_island_faces!=NULL
                                ? fit_reference_island_faces:component_faces;
        if (pack_cover_intervals) {
            size_t atlas_islands = 0;
            double cover_span = 0.0, atlas_span = 0.0;
            if (wwl_pack_component_intervals(
                    nv, uv, pack_component, pack_count, pack_faces,
                    gap,
                    &atlas_islands, &cover_span, &atlas_span) != 0) {
                fprintf(stderr,
                    "welded_winding_lift: cover-interval packing failed\n");
                goto cleanup;
            }
            fprintf(stderr,
                "welded_winding_lift: metric cover pack %zu component(s) -> "
                "%zu interval island(s), span %.6g -> %.6g vox (gap %.3g)\n",
                pack_count, atlas_islands, cover_span, atlas_span, gap);
        } else if (wwl_pack_components(
                       nv, uv, pack_component, pack_count,
                       pack_faces, gap) != 0) {
            fprintf(stderr, "welded_winding_lift: component packing failed\n");
            goto cleanup;
        }
    }
    size_t positive = 0, negative = 0, degenerate = 0;
    wwl_face_orientation_stats(uv, faces, nf,
                               &positive, &negative, &degenerate);
    fprintf(stderr,
        "welded_winding_lift: UV face signs +%zu -%zu zero=%zu "
        "(component-wise audit still required)\n",
        positive, negative, degenerate);
    size_t output_faces = final_source_only ? source_faces : nf;
    if (ves_ensure_parent_dir(output_path) != 0) {
        fprintf(stderr, "welded_winding_lift: cannot write %s\n", output_path);
        goto cleanup;
    }
    {   /* authoritative UV'd VMESH companion (obj_bake_raw reads only this);
         * the text OBJ is a viewer dump behind --dump-obj */
        char binary_path[2048];
        float *uv_f = (float *)malloc((nv ? nv : 1) * 2 * sizeof(*uv_f));
        int wrc = 0;
        if (uv_f == NULL) goto cleanup;
        for (size_t v = 0; v < nv; v++) {
            uv_f[v * 2 + 0] = (float)uv[v * 2 + 0];
            uv_f[v * 2 + 1] = (float)uv[v * 2 + 1];
        }
        wrc = MeshBin_companion_path(output_path, binary_path,
                                     sizeof binary_path) != 0 ||
              strcmp(binary_path, output_path) == 0 ||
              MeshBin_write(binary_path, verts, nv, faces,
                            output_faces, uv_f) != 0;
        free(uv_f);
        if (wrc) {
            fprintf(stderr,
                "welded_winding_lift: cannot write authoritative companion "
                "for %s\n", output_path);
            goto cleanup;
        }
    }
    if (dump_obj &&
        ObjIO_write_uv_double(output_path, verts, nv, faces,
                              output_faces, uv) != 0) {
        fprintf(stderr, "welded_winding_lift: cannot write %s\n", output_path);
        goto cleanup;
    }
    fprintf(stderr, "welded_winding_lift: wrote %s (%zu/%zu faces%s%s)\n",
            output_path, output_faces, nf,
            final_source_only ? ", final source prefix only" : "",
            dump_obj ? "" : ", vmesh companion only");
    result = 0;

cleanup:
    wwl_chart_layout_dispose(&chart_layout);
    wwl_lineage_dispose(&lineage);
    wwl_graph_dispose(&source_graph);
    wwl_graph_dispose(&full_graph);
    wwl_graph_dispose(&bridge_graph);
    free(anchor_source);
    free(chart);
    free(chart_faces);
    free(component);
    free(component_faces);
    free(uv);
    free(registered_uv);
    free(registered_phi);
    free(fit_phi);
    free(fit_atlas_shift);
    free(fit_reference_uv);
    free(fit_reference_u);
    free(fit_reference_island);
    free(fit_reference_island_faces);
    free(fit_relation_a);
    free(fit_relation_b);
    free(fit_relation_support);
    free(chart_layout_correction);
    Arena_dispose(&arena);
    return result;
}
