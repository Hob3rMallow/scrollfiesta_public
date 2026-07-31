/*
 * atlas_overlap_fix_tool.c
 *
 * Move welded sheet fragments by whole wraps until they stop overlapping, and
 * park what cannot be placed.  Reads a scroll_whole placed directory and its
 * registered UV, writes before/after parameter-space and world-space OBJs plus
 * a per-group decision table.
 */

#include "../common/arena.h"
#include "../common/ves_platform.h"
#include "../unroll/piece_set.h"
#include "../unroll/scaffold.h"
#include "../common/ves_png.h"
#include "../whole/atlas_overlap_fix.h"
#include "../whole/atlas_register.h"
#include "../whole/atlas_solution.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define AOF_PATH_CAP 2048

static const char *action_name(int action)
{
    switch (action) {
    case ATLAS_OVERLAP_FIX_SETTLED: return "settled";
    case ATLAS_OVERLAP_FIX_SHIFTED: return "shifted";
    case ATLAS_OVERLAP_FIX_PARKED:  return "parked";
    case ATLAS_OVERLAP_FIX_STUCK:   return "stuck";
    default:                        return "unknown";
    }
}

static FILE *open_out(const char *dir, const char *name)
{
    char path[AOF_PATH_CAP];
    int n = snprintf(path, AOF_PATH_CAP, "%s/%s", dir, name);
    if (n < 0 || n >= AOF_PATH_CAP) return NULL;
    if (ves_ensure_parent_dir(path) != 0) return NULL;
    return fopen(path, "wb");
}

static int parse_long_full(const char *text, long *out)
{
    char *end = NULL;
    long value = strtol(text, &end, 10);
    if (end == text || *end != '\0') return -1;
    *out = value;
    return 0;
}

/* Sidecars are the block-coordinate handoff: a tiled driver commits only the
 * movable core files into its one global placed state and leaves the halo
 * files untouched.  Phi is source data; only solved u/v are replaced. */
static int write_uvphi_sidecars(const char *dir, const PieceSet *ps,
                                const double *u, const double *v)
{
    for (size_t q = 0; q < ps->n_cubes; q++) {
        size_t first = ps->cube_voff[q], last = ps->cube_voff[q + 1];
        size_t n = last - first;
        if (n == 0 || n > SIZE_MAX / (3 * sizeof(float))) return -1;
        float *uvphi = (float *)malloc(n * 3 * sizeof(float));
        if (uvphi == NULL) return -1;
        for (size_t i = 0; i < n; i++) {
            uvphi[i * 3 + 0] = (float)u[first + i];
            uvphi[i * 3 + 1] = (float)v[first + i];
            uvphi[i * 3 + 2] = ps->phi[first + i];
        }
        char name[128];
        int len = snprintf(name, sizeof name, "solved_uvphi/%s_uvphi.f32",
                           ps->ids[q]);
        FILE *fp = len > 0 && (size_t)len < sizeof name
                 ? open_out(dir, name) : NULL;
        int failed = fp == NULL ||
            fwrite(uvphi, sizeof(float), n * 3, fp) != n * 3;
        if (fp != NULL && fclose(fp) != 0) failed = 1;
        free(uvphi);
        if (failed) return -1;
    }
    return 0;
}

static void group_color(int32_t group, double rgb[3])
{
    unsigned x = (unsigned)(group >= 0 ? group : 0) * 2654435761u;
    rgb[0] = 0.25 + 0.70 * (double)((x >>  0) & 255u) / 255.0;
    rgb[1] = 0.25 + 0.70 * (double)((x >>  8) & 255u) / 255.0;
    rgb[2] = 0.25 + 0.70 * (double)((x >> 16) & 255u) / 255.0;
}

/* Parameter-space mesh, one colour per welded group. */
static int write_uv_obj(const char *dir, const char *name, const PieceSet *ps,
                        const double *u, const double *v,
                        const uint8_t *face_keep, const int32_t *vertex_group)
{
    FILE *fp = open_out(dir, name);
    if (fp == NULL) return -1;
    for (size_t i = 0; i < ps->nv; i++) {
        double rgb[3] = {0.5, 0.5, 0.5};
        if (vertex_group != NULL) group_color(vertex_group[i], rgb);
        fprintf(fp, "v %.9g %.9g 0 %.4f %.4f %.4f\n",
                isfinite(u[i]) ? u[i] : 0.0, isfinite(v[i]) ? v[i] : 0.0,
                rgb[0], rgb[1], rgb[2]);
    }
    for (size_t f = 0; f < ps->nf; f++) {
        if (face_keep != NULL && !face_keep[f]) continue;
        fprintf(fp, "f %d %d %d\n", ps->faces[f * 3] + 1,
                ps->faces[f * 3 + 1] + 1, ps->faces[f * 3 + 2] + 1);
    }
    fclose(fp);
    return 0;
}

/* World-space mesh coloured by group, so a moved fragment is identifiable in
 * the geometry it came from. */
static int write_xyz_obj(const char *dir, const char *name, const PieceSet *ps,
                         const uint8_t *face_keep, const int32_t *vertex_group)
{
    FILE *fp = open_out(dir, name);
    if (fp == NULL) return -1;
    for (size_t i = 0; i < ps->nv; i++) {
        double rgb[3] = {0.5, 0.5, 0.5};
        if (vertex_group != NULL) group_color(vertex_group[i], rgb);
        fprintf(fp, "v %.6g %.6g %.6g %.4f %.4f %.4f\n",
                (double)ps->verts[i * 3], (double)ps->verts[i * 3 + 1],
                (double)ps->verts[i * 3 + 2], rgb[0], rgb[1], rgb[2]);
    }
    for (size_t f = 0; f < ps->nf; f++) {
        if (face_keep != NULL && !face_keep[f]) continue;
        fprintf(fp, "f %d %d %d\n", ps->faces[f * 3] + 1,
                ps->faces[f * 3 + 1] + 1, ps->faces[f * 3 + 2] + 1);
    }
    fclose(fp);
    return 0;
}

static int write_groups_csv(const char *dir, const AtlasOverlapFixGroup *group,
                            size_t ngroups)
{
    FILE *fp = open_out(dir, "groups.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "group,color,action,charts,faces,vertices,mean_radius,"
                "radius_spread,circumference,wraps,shift_u,shift_v,"
                "overlap_before,overlap_after,intra_overlap,"
                "intra_wrongwind\n");
    for (size_t g = 0; g < ngroups; g++) {
        const AtlasOverlapFixGroup *e = &group[g];
        double rgb[3];
        group_color(e->group, rgb);
        fprintf(fp, "%d,#%02x%02x%02x,%s,%zu,%zu,%zu,%.6g,%.6g,%.6g,%d,%.6g,"
                    "%.6g,%zu,%zu,%zu,%zu\n",
                e->group, (unsigned)(rgb[0] * 255.0),
                (unsigned)(rgb[1] * 255.0), (unsigned)(rgb[2] * 255.0),
                action_name(e->action), e->ncharts, e->nfaces,
                e->nvertices, e->mean_radius, e->radius_spread,
                e->circumference, e->wraps, e->shift_u, e->shift_v,
                e->overlap_before, e->overlap_after, e->intra_overlap,
                e->intra_wrongwind);
    }
    fclose(fp);
    return 0;
}

/* Winding depth colour: 0 = grey, outward = warm, inward = cool, saturating
 * with |k| so a staircase of depths reads as a gradient. */
static void depth_color(int32_t k, double rgb[3])
{
    if (k == 0) { rgb[0] = rgb[1] = rgb[2] = 0.62; return; }
    double t = (double)(k < 0 ? -k : k);
    if (t > 4.0) t = 4.0;
    t /= 4.0;
    if (k > 0) {
        rgb[0] = 0.55 + 0.45 * t;
        rgb[1] = 0.45 - 0.25 * t;
        rgb[2] = 0.15;
    } else {
        rgb[0] = 0.15;
        rgb[1] = 0.45 - 0.25 * t;
        rgb[2] = 0.55 + 0.45 * t;
    }
}

/* Solved parameter-space mesh coloured by assigned winding depth. */
static int write_depth_obj(const char *dir, const PieceSet *ps,
                           const double *u, const double *v,
                           const uint8_t *face_keep,
                           const AtlasOverlapFixTabuReport *rep)
{
    FILE *fp = open_out(dir, "uv_depth.obj");
    if (fp == NULL) return -1;
    for (size_t i = 0; i < ps->nv; i++) {
        double rgb[3] = {0.3, 0.3, 0.3};
        int32_t c = rep->vertex_chart != NULL ? rep->vertex_chart[i] : -1;
        if (c >= 0 && (size_t)c < rep->nchart)
            depth_color(rep->chart[c].k, rgb);
        fprintf(fp, "v %.9g %.9g 0 %.4f %.4f %.4f\n",
                isfinite(u[i]) ? u[i] : 0.0, isfinite(v[i]) ? v[i] : 0.0,
                rgb[0], rgb[1], rgb[2]);
    }
    for (size_t f = 0; f < ps->nf; f++) {
        if (face_keep != NULL && !face_keep[f]) continue;
        fprintf(fp, "f %d %d %d\n", ps->faces[f * 3] + 1,
                ps->faces[f * 3 + 1] + 1, ps->faces[f * 3 + 2] + 1);
    }
    fclose(fp);
    return 0;
}

/* One OBJ per winding layer: only the faces of charts assigned depth k, in
 * the PRE-relayout tabu placement, original group colours.  This is the view
 * that shows whether a region moved whole (one colour lands in one layer) or
 * was shredded (one colour smeared across layers). */
static int write_layer_objs(const char *dir, int round, int span,
                            const PieceSet *ps,
                            const double *u, const double *v,
                            const uint8_t *face_keep,
                            const int32_t *vertex_group,
                            const AtlasOverlapFixTabuReport *rep)
{
    if (rep->vertex_chart == NULL || u == NULL) return 0;
    for (int k = -span; k <= span; k++) {
        int any = 0;
        for (size_t c = 0; c < rep->nchart && !any; c++)
            if (rep->chart[c].k == k) any = 1;
        if (!any) continue;
        char name[80];
        snprintf(name, sizeof(name), "round%d_wind%+d_uv.obj", round, k);
        FILE *fp = open_out(dir, name);
        if (fp == NULL) return -1;
        for (size_t i = 0; i < ps->nv; i++) {
            double rgb[3] = {0.5, 0.5, 0.5};
            if (vertex_group != NULL) group_color(vertex_group[i], rgb);
            fprintf(fp, "v %.9g %.9g 0 %.4f %.4f %.4f\n",
                    isfinite(u[i]) ? u[i] : 0.0,
                    isfinite(v[i]) ? v[i] : 0.0, rgb[0], rgb[1], rgb[2]);
        }
        for (size_t f = 0; f < ps->nf; f++) {
            if (face_keep != NULL && !face_keep[f]) continue;
            int32_t c = rep->vertex_chart[ps->faces[f * 3]];
            if (c < 0 || (size_t)c >= rep->nchart || rep->chart[c].k != k)
                continue;
            fprintf(fp, "f %d %d %d\n", ps->faces[f * 3] + 1,
                    ps->faces[f * 3 + 1] + 1, ps->faces[f * 3 + 2] + 1);
        }
        fclose(fp);
    }
    return 0;
}

static int write_charts_csv(const char *dir, const char *name,
                            const PieceSet *ps,
                            const AtlasOverlapFixTabuReport *rep)
{
    FILE *fp = open_out(dir, name);
    if (fp == NULL) return -1;
    fprintf(fp, "chart,cube,group,neighbourhood,phase_flat,phase_region,"
                "region_offset,faces,vertices,"
                "qc_max,radius_med,u_med,w_phi0,"
                "w_rad,w_phase,phase_gauge,phase_bin,wind,k,gbin,shift_u,"
                "relayout_rank,relayout_metric_shift_u,relayout_shift_u,"
                "cells_initial,cells_final\n");
    for (size_t c = 0; c < rep->nchart; c++) {
        const AtlasOverlapFixChart *e = &rep->chart[c];
        const char *cube = e->cube >= 0 && (size_t)e->cube < ps->n_cubes
                         ? ps->ids[e->cube] : "?";
        fprintf(fp, "%d,%s,%d,%d,%d,%d,%d,%zu,%zu,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%d,%d,"
                    "%.6g,%d,%d,%.6g,%zu,%.6g,%.6g,%lld,%lld\n",
                e->chart, cube, e->group, e->neighbourhood, e->phase_flat,
                e->phase_region, e->region_offset,
                e->nfaces, e->nvertices,
                e->qc_max, e->radius_med, e->u_med, e->w_phi0, e->w_rad,
                e->w_phase, e->phase_gauge, e->phase_bin, e->wind,
                e->k, e->gbin,
                e->shift_u, e->relayout_rank, e->relayout_metric_shift_u,
                e->relayout_shift_u,
                (long long)e->cells_initial,
                (long long)e->cells_final);
    }
    fclose(fp);
    return 0;
}

static int write_moves_csv(const char *dir,
                           const AtlasOverlapFixTabuReport *rep)
{
    FILE *fp = open_out(dir, "tabu_moves.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "iter,kind,id,from,to,d_overlap,d_neighbour,d_radial,"
                "energy_after\n");
    for (size_t m = 0; m < rep->nmove; m++) {
        const AtlasOverlapFixTabuMove *e = &rep->move[m];
        const char *kind = e->kind == 0 ? "chart"
                         : e->kind == 1 ? "group"
                         : e->kind == 2 ? "nbhd"
                         : e->kind == 3 ? "graded" : "region";
        fprintf(fp, "%d,%s,%d,%d,%d,%.6g,%.6g,%.6g,%.6g\n",
                e->iter, kind, e->id, e->from, e->to,
                e->d_overlap, e->d_neighbour, e->d_radial, e->energy_after);
    }
    fclose(fp);
    return 0;
}

/* One row per chart changed by an accepted move.  Unlike the compact move
 * log, this expands compound moves, making the exact onset of fragmentation
 * directly inspectable without reconstructing internal neighbourhoods. */
static int write_chart_moves_csv(const char *dir,
                                 const AtlasOverlapFixTabuReport *rep)
{
    FILE *fp = open_out(dir, "tabu_chart_moves.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "iter,kind,id,chart,group,gbin,k_before,k_after\n");
    if (rep->state_k != NULL && rep->nstate > 1) {
        size_t nm = rep->nmove < rep->nstate - 1
                  ? rep->nmove : rep->nstate - 1;
        for (size_t m = 0; m < nm; m++) {
            const AtlasOverlapFixTabuMove *e = &rep->move[m];
            const char *kind = e->kind == 0 ? "chart"
                             : e->kind == 1 ? "group"
                             : e->kind == 2 ? "nbhd"
                             : e->kind == 3 ? "graded" : "region";
            const int8_t *before = &rep->state_k[m * rep->nchart];
            const int8_t *after = &rep->state_k[(m + 1) * rep->nchart];
            for (size_t c = 0; c < rep->nchart; c++) {
                if (before[c] == after[c]) continue;
                fprintf(fp, "%d,%s,%d,%zu,%d,%d,%d,%d\n",
                        e->iter, kind, e->id, c, rep->chart[c].group,
                        rep->chart[c].gbin, (int)before[c], (int)after[c]);
            }
        }
    }
    fclose(fp);
    return 0;
}

static int write_candidates_csv(const char *dir,
                                const AtlasOverlapFixTabuReport *rep)
{
    FILE *fp = open_out(dir, "tabu_candidates.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "candidate,state,exact_cross,exact_intra,search_overlap,"
                "neighbour,radial,energy,layout_abs,exact_score,"
                "weld_mean_abs,weld_rms,weld_max,weld_jumps_25,"
                "weld_jumps_100,winner\n");
    for (size_t i = 0; i < rep->ncandidate; i++) {
        const AtlasOverlapFixTabuCandidate *c = &rep->candidate[i];
        fprintf(fp, "%zu,%zu,%zu,%zu,%.9g,%.9g,%.9g,%.9g,%zu,%.9g,"
                    "%.9g,%.9g,%.9g,%zu,%zu,%d\n",
                i, c->state, c->exact_cross, c->exact_intra,
                c->search_overlap, c->neighbour, c->radial, c->energy,
                c->layout_abs, c->exact_score, c->weld_mean_abs,
                c->weld_rms, c->weld_max, c->weld_jumps_25,
                c->weld_jumps_100, c->winner);
    }
    fclose(fp);
    return 0;
}

static int write_edges_csv(const char *dir,
                           const AtlasOverlapFixTabuReport *rep)
{
    FILE *fp = open_out(dir, "tabu_edges.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "edge,kind,chart0,group0,chart1,group1,target,final_delta,"
                "torn,family,weight,du,final_du,phase_turn,phase_nearest,"
                "phase_residual\n");
    for (size_t i = 0; i < rep->nedge; i++) {
        const AtlasOverlapFixTabuEdge *e = &rep->edge[i];
        int32_t g0 = rep->chart[e->chart0].group;
        int32_t g1 = rep->chart[e->chart1].group;
        int nearest = (int)floor(e->phase_turn + 0.5);
        fprintf(fp, "%d,%s,%d,%d,%d,%d,%d,%d,%d,%d,%.9g,%.9g,%.9g,"
                    "%.9g,%d,%.9g\n",
                e->edge, e->lateral ? "lateral" : "weld",
                e->chart0, g0, e->chart1, g1, e->target, e->final_delta,
                e->final_delta != e->target, e->family, e->weight, e->du,
                e->final_du, e->phase_turn, nearest, e->phase_residual);
    }
    fclose(fp);
    return 0;
}
static int write_overlap_pairs_csv(const char *dir,
                                   const AtlasOverlapFixTabuReport *rep)
{
    FILE *fp = open_out(dir, "tabu_overlap_pairs.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "chart0,group0,chart1,group1,face_pairs,direct_kind,edge,"
                "target,final_delta,intact,radius_delta,min_xyz\n");
    for (size_t i = 0; i < rep->noverlap; i++) {
        const AtlasOverlapFixTabuOverlap *o = &rep->overlap[i];
        int32_t g0 = rep->chart[o->chart0].group;
        int32_t g1 = rep->chart[o->chart1].group;
        const char *kind = o->edge < 0 ? "none" :
                           o->lateral ? "lateral" : "weld";
        fprintf(fp, "%d,%d,%d,%d,%u,%s,%d,%d,%d,%d,%.9g,%.9g\n",
                o->chart0, g0, o->chart1, g1, o->face_pairs, kind,
                o->edge, o->target, o->final_delta,
                o->edge >= 0 && o->final_delta == o->target,
                o->radius_delta, o->min_xyz);
    }
    fclose(fp);
    return 0;
}


/* The DIAGONAL: u across, radius up, one splat per kept vertex in its group
 * colour.  A correct arc-length atlas of a spiral is a clean rising diagonal
 * -- one wind per circumference; misplaced charts sit off the line and
 * collapsed winds stack vertically. */

static const AtlasOverlapFixTabuEdge *find_relayout_edge(
    const AtlasOverlapFixTabuReport *rep, int32_t chart0, int32_t chart1)
{
    for (size_t i = 0; i < rep->nedge; i++) {
        const AtlasOverlapFixTabuEdge *e = &rep->edge[i];
        if ((e->chart0 == chart0 && e->chart1 == chart1) ||
            (e->chart0 == chart1 && e->chart1 == chart0))
            return e;
    }
    return NULL;
}

static int write_relayout_bounds_csv(
    const char *dir, const AtlasOverlapFixTabuReport *rep)
{
    FILE *fp = open_out(dir, "relayout_collision_bounds.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "chart_lo,group_lo,rank_lo,chart_hi,group_hi,rank_hi,"
                "source_face_lo,source_face_hi,lower,final_delta,slack,"
                "initial,active,updates,topology,collision,direct_kind,edge,"
                "target,discrete_delta,metric_delta\n");
    for (size_t i = 0; i < rep->nrelayout_bound; i++) {
        const AtlasOverlapFixRelayoutBound *b = &rep->relayout_bound[i];
        if (b->chart_lo < 0 || b->chart_hi < 0 ||
            (size_t)b->chart_lo >= rep->nchart ||
            (size_t)b->chart_hi >= rep->nchart) {
            fclose(fp);
            return -1;
        }
        const AtlasOverlapFixChart *lo = &rep->chart[b->chart_lo];
        const AtlasOverlapFixChart *hi = &rep->chart[b->chart_hi];
        const AtlasOverlapFixTabuEdge *edge = find_relayout_edge(
            rep, b->chart_lo, b->chart_hi);
        const char *kind = edge == NULL ? "none" :
                           edge->lateral ? "lateral" : "weld";
        double metric_delta = hi->relayout_metric_shift_u -
                              lo->relayout_metric_shift_u;
        fprintf(fp, "%d,%d,%zu,%d,%d,%zu,%d,%d,%.9g,%.9g,%.9g,"
                    "%d,%d,%u,%d,%d,%s,%d,%d,%d,%.9g\n",
                b->chart_lo, lo->group, lo->relayout_rank,
                b->chart_hi, hi->group, hi->relayout_rank,
                b->face_lo, b->face_hi, b->lower, b->final_delta,
                b->slack, b->initial, b->active, b->updates,
                b->topology, b->collision, kind,
                edge != NULL ? edge->edge : -1,
                edge != NULL ? edge->target : 0,
                edge != NULL ? edge->final_delta : 0, metric_delta);
    }
    fclose(fp);
    return 0;
}

static int write_relayout_residual_csv(
    const char *dir, const AtlasOverlapFixTabuReport *rep)
{
    FILE *fp = open_out(dir, "relayout_residual_pairs.csv");
    if (fp == NULL) return -1;
    fprintf(fp, "chart0,group0,face0,chart1,group1,face1,allowed,class,"
                "spatial_class,face_min_xyz,radius_delta,same_chart,"
                "chart0_faces,chart1_faces,chart0_qc_max,chart1_qc_max\n");
    for (size_t i = 0; i < rep->nrelayout_residual; i++) {
        const AtlasOverlapFixRelayoutResidual *r =
            &rep->relayout_residual[i];
        if (r->chart0 < 0 || r->chart1 < 0 ||
            (size_t)r->chart0 >= rep->nchart ||
            (size_t)r->chart1 >= rep->nchart) {
            fclose(fp);
            return -1;
        }
        const AtlasOverlapFixChart *a = &rep->chart[r->chart0];
        const AtlasOverlapFixChart *b = &rep->chart[r->chart1];
        const char *classification = r->reason == 1 ? "direct" :
                                     r->reason == 2 ? "xyz_near" :
                                     r->reason == 3 ? "self" :
                                     r->reason == 4 ? "topology_limited" : "hard";
        const char *spatial = !(r->min_xyz >= 0.0) ? "unknown" :
                              r->min_xyz <= 8.0 ? "near_0_8" :
                              r->min_xyz <= 24.0 ? "local_8_24" :
                                                  "remote_gt24";
        fprintf(fp, "%d,%d,%d,%d,%d,%d,%d,%s,%s,%.9g,%.9g,%d,"
                    "%zu,%zu,%.9g,%.9g\n",
                r->chart0, a->group, r->face0,
                r->chart1, b->group, r->face1,
                r->allowed, classification, spatial,
                r->min_xyz, r->radius_delta, r->chart0 == r->chart1,
                a->nfaces, b->nfaces,
                a->qc_max, b->qc_max);
    }
    fclose(fp);
    return 0;
}

typedef struct {
    int32_t chart0;
    int32_t chart1;
    size_t total;
    size_t direct;
    size_t xyz_near;
    size_t self;
    size_t topology_limited;
    size_t hard;
    size_t spatial_near;
    size_t spatial_local;
    size_t spatial_remote;
    size_t spatial_unknown;
    double min_face_xyz;
    double max_face_xyz;
} RelayoutResidualPairSummary;

static int relayout_residual_pair_compare(const void *left, const void *right)
{
    const RelayoutResidualPairSummary *a =
        (const RelayoutResidualPairSummary *)left;
    const RelayoutResidualPairSummary *b =
        (const RelayoutResidualPairSummary *)right;
    if (a->chart0 != b->chart0) return a->chart0 < b->chart0 ? -1 : 1;
    if (a->chart1 != b->chart1) return a->chart1 < b->chart1 ? -1 : 1;
    return 0;
}

static const AtlasOverlapFixTabuOverlap *find_discrete_overlap(
    const AtlasOverlapFixTabuReport *rep, int32_t chart0, int32_t chart1)
{
    for (size_t i = 0; i < rep->noverlap; i++) {
        const AtlasOverlapFixTabuOverlap *o = &rep->overlap[i];
        if ((o->chart0 == chart0 && o->chart1 == chart1) ||
            (o->chart0 == chart1 && o->chart1 == chart0))
            return o;
    }
    return NULL;
}

static int write_relayout_residual_chart_pairs_csv(
    const char *dir, const AtlasOverlapFixTabuReport *rep)
{
    if (rep->nchart > 0 && rep->nchart > SIZE_MAX / rep->nchart) return -1;
    size_t matrix_size = rep->nchart * rep->nchart;
    int32_t *index = (int32_t *)malloc(
        (matrix_size ? matrix_size : 1) * sizeof(*index));
    RelayoutResidualPairSummary *summary =
        (RelayoutResidualPairSummary *)calloc(
            rep->nrelayout_residual ? rep->nrelayout_residual : 1,
            sizeof(*summary));
    if (index == NULL || summary == NULL) {
        free(index);
        free(summary);
        return -1;
    }
    for (size_t i = 0; i < matrix_size; i++) index[i] = -1;

    size_t nsummary = 0;
    for (size_t i = 0; i < rep->nrelayout_residual; i++) {
        const AtlasOverlapFixRelayoutResidual *r =
            &rep->relayout_residual[i];
        if (r->chart0 < 0 || r->chart1 < 0 ||
            (size_t)r->chart0 >= rep->nchart ||
            (size_t)r->chart1 >= rep->nchart) {
            free(index);
            free(summary);
            return -1;
        }
        int32_t c0 = r->chart0, c1 = r->chart1;
        if (c0 > c1) { int32_t t = c0; c0 = c1; c1 = t; }
        size_t at = (size_t)c0 * rep->nchart + (size_t)c1;
        int32_t si = index[at];
        if (si < 0) {
            if (nsummary > (size_t)INT32_MAX) {
                free(index);
                free(summary);
                return -1;
            }
            si = (int32_t)nsummary++;
            index[at] = si;
            summary[si].chart0 = c0;
            summary[si].chart1 = c1;
            summary[si].min_face_xyz = DBL_MAX;
            summary[si].max_face_xyz = -DBL_MAX;
        }
        RelayoutResidualPairSummary *s = &summary[si];
        s->total++;
        if (r->reason == 1) s->direct++;
        else if (r->reason == 2) s->xyz_near++;
        else if (r->reason == 3) s->self++;
        else if (r->reason == 4) s->topology_limited++;
        else s->hard++;
        if (!(r->min_xyz >= 0.0)) s->spatial_unknown++;
        else {
            if (r->min_xyz <= 8.0) s->spatial_near++;
            else if (r->min_xyz <= 24.0) s->spatial_local++;
            else s->spatial_remote++;
            if (r->min_xyz < s->min_face_xyz)
                s->min_face_xyz = r->min_xyz;
            if (r->min_xyz > s->max_face_xyz)
                s->max_face_xyz = r->min_xyz;
        }
    }
    free(index);
    qsort(summary, nsummary, sizeof(*summary),
          relayout_residual_pair_compare);

    FILE *fp = open_out(dir, "relayout_residual_chart_pairs.csv");
    if (fp == NULL) {
        free(summary);
        return -1;
    }
    fprintf(fp, "chart0,group0,chart1,group1,total_face_pairs,direct_pairs,"
                "xyz_near_pairs,self_pairs,topology_limited_pairs,hard_pairs,"
                "face_xyz_near_0_8,face_xyz_local_8_24,face_xyz_remote_gt24,"
                "face_xyz_unknown,pair_spatial_class,"
                "whole_chart_remote_gt24,min_face_xyz,max_face_xyz,radius_delta,"
                "initial_face_pairs,chart_min_xyz,direct_kind,intact,"
                "rank_delta\n");
    for (size_t i = 0; i < nsummary; i++) {
        const RelayoutResidualPairSummary *s = &summary[i];
        const AtlasOverlapFixChart *a = &rep->chart[s->chart0];
        const AtlasOverlapFixChart *b = &rep->chart[s->chart1];
        const AtlasOverlapFixTabuOverlap *initial = find_discrete_overlap(
            rep, s->chart0, s->chart1);
        const AtlasOverlapFixTabuEdge *edge = find_relayout_edge(
            rep, s->chart0, s->chart1);
        const char *kind = edge == NULL ? "none" :
                           edge->lateral ? "lateral" : "weld";
        int intact = edge != NULL && edge->final_delta == edge->target;
        const char *pair_spatial = s->spatial_near > 0 ? "near_contact" :
                                   s->spatial_local > 0 ? "local_or_mixed" :
                                   s->spatial_remote > 0 ? "remote_only" :
                                                          "unknown";
        int chart_remote = initial != NULL && initial->min_xyz > 24.0;
        double min_xyz = s->min_face_xyz < DBL_MAX
                       ? s->min_face_xyz : -1.0;
        double max_xyz = s->max_face_xyz > -DBL_MAX
                       ? s->max_face_xyz : -1.0;
        long long rank_delta = (long long)b->relayout_rank -
                               (long long)a->relayout_rank;
        fprintf(fp, "%d,%d,%d,%d,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu,"
                    "%zu,%s,%d,%.9g,%.9g,%.9g,%u,%.9g,%s,%d,%lld\n",
                s->chart0, a->group, s->chart1, b->group,
                s->total, s->direct, s->xyz_near, s->self,
                s->topology_limited, s->hard,
                s->spatial_near, s->spatial_local, s->spatial_remote,
                s->spatial_unknown, pair_spatial, chart_remote, min_xyz, max_xyz,
                fabs(a->radius_med - b->radius_med),
                initial != NULL ? initial->face_pairs : 0,
                initial != NULL ? initial->min_xyz : -1.0,
                kind, intact, rank_delta);
    }
    fclose(fp);
    free(summary);
    return 0;
}

static int write_diag_png(const char *dir, const char *name,
                          const PieceSet *ps, const ScaffoldCalib *cal,
                          const double *u, const int32_t *vertex_group)
{
    enum { DW = 2400, DH = 900 };
    double u_lo = 0.0, u_hi = 0.0, r_lo = 0.0, r_hi = 0.0;
    int have = 0;
    double *radius = (double *)malloc(ps->nv * sizeof(double));
    uint8_t *img = (uint8_t *)calloc((size_t)DW * DH * 3, 1);
    if (radius == NULL || img == NULL) {
        free(radius);
        free(img);
        return -1;
    }
    for (size_t i = 0; i < ps->nv; i++) {
        const float *p = &ps->verts[i * 3];
        double d0 = (double)p[0] - (double)cal->axis_point[0];
        double d1 = (double)p[1] - (double)cal->axis_point[1];
        double d2 = (double)p[2] - (double)cal->axis_point[2];
        double along = d0 * (double)cal->axis_dir[0] +
                       d1 * (double)cal->axis_dir[1] +
                       d2 * (double)cal->axis_dir[2];
        double r0 = d0 - along * (double)cal->axis_dir[0];
        double r1 = d1 - along * (double)cal->axis_dir[1];
        double r2 = d2 - along * (double)cal->axis_dir[2];
        radius[i] = sqrt(r0 * r0 + r1 * r1 + r2 * r2);
        if (vertex_group == NULL || vertex_group[i] < 0 || !isfinite(u[i]))
            continue;
        if (!have) {
            u_lo = u_hi = u[i];
            r_lo = r_hi = radius[i];
            have = 1;
        }
        if (u[i] < u_lo) u_lo = u[i];
        if (u[i] > u_hi) u_hi = u[i];
        if (radius[i] < r_lo) r_lo = radius[i];
        if (radius[i] > r_hi) r_hi = radius[i];
    }
    if (!have || !(u_hi > u_lo) || !(r_hi > r_lo)) {
        free(radius);
        free(img);
        return 0;
    }
    for (size_t i = 0; i < ps->nv; i++) {
        if (vertex_group == NULL || vertex_group[i] < 0 || !isfinite(u[i]))
            continue;
        int x = (int)((u[i] - u_lo) / (u_hi - u_lo) * (DW - 1));
        int y = (DH - 1) -
                (int)((radius[i] - r_lo) / (r_hi - r_lo) * (DH - 1));
        if (x < 0 || x >= DW || y < 0 || y >= DH) continue;
        double rgb[3];
        group_color(vertex_group[i], rgb);
        uint8_t *px = &img[((size_t)y * DW + (size_t)x) * 3];
        px[0] = (uint8_t)(rgb[0] * 255.0);
        px[1] = (uint8_t)(rgb[1] * 255.0);
        px[2] = (uint8_t)(rgb[2] * 255.0);
    }
    char path[AOF_PATH_CAP];
    int n = snprintf(path, AOF_PATH_CAP, "%s/%s", dir, name);
    int rc = -1;
    if (n > 0 && n < AOF_PATH_CAP && ves_ensure_parent_dir(path) == 0)
        rc = VesPng_write_rgb(path, img, DW, DH);
    free(radius);
    free(img);
    return rc;
}

/* World-space mesh with the solved atlas as vt, ready for obj_bake_raw:
 * bake RAW CT through this to eyeball the textured strips. */
static int write_bake_obj(const char *dir, const PieceSet *ps,
                          const double *u, const double *v,
                          const uint8_t *face_keep)
{
    FILE *fp = open_out(dir, "atlas_bake.obj");
    if (fp == NULL) return -1;
    /* The rawtex raster's origin is vt (0,0): shift the atlas into the
     * positive quadrant or every face lands outside the image. */
    double u0 = DBL_MAX, v0 = DBL_MAX;
    for (size_t f = 0; f < ps->nf; f++) {
        if (face_keep != NULL && !face_keep[f]) continue;
        for (int k = 0; k < 3; k++) {
            int32_t i = ps->faces[f * 3 + (size_t)k];
            if (!isfinite(u[i]) || !isfinite(v[i])) continue;
            if (u[i] < u0) u0 = u[i];
            if (v[i] < v0) v0 = v[i];
        }
    }
    if (u0 == DBL_MAX) { u0 = 0.0; v0 = 0.0; }
    for (size_t i = 0; i < ps->nv; i++)
        fprintf(fp, "v %.6g %.6g %.6g\n", (double)ps->verts[i * 3],
                (double)ps->verts[i * 3 + 1], (double)ps->verts[i * 3 + 2]);
    for (size_t i = 0; i < ps->nv; i++)
        fprintf(fp, "vt %.9g %.9g\n",
                isfinite(u[i]) ? u[i] - u0 : 0.0,
                isfinite(v[i]) ? v[i] - v0 : 0.0);
    for (size_t f = 0; f < ps->nf; f++) {
        if (face_keep != NULL && !face_keep[f]) continue;
        int32_t a = ps->faces[f * 3] + 1;
        int32_t b = ps->faces[f * 3 + 1] + 1;
        int32_t c = ps->faces[f * 3 + 2] + 1;
        fprintf(fp, "f %d/%d %d/%d %d/%d\n", a, a, b, b, c, c);
    }
    fclose(fp);
    return 0;
}

static int write_stats_json(const char *dir, const AtlasOverlapFixStats *s,
                            double seconds)
{
    FILE *fp = open_out(dir, "overlap_fix_stats.json");
    if (fp == NULL) return -1;
    fprintf(fp,
        "{\n"
        "  \"total_charts\": %zu,\n"
        "  \"charts_before_cull\": %zu,\n"
        "  \"total_groups\": %zu,\n"
        "  \"largest_group_charts\": %zu,\n"
        "  \"total_faces\": %zu,\n"
        "  \"faces_killed\": %zu,\n"
        "  \"faces_killed_stretch\": %zu,\n"
        "  \"faces_killed_charts\": %zu,\n"
        "  \"charts_killed_small\": %zu,\n"
        "  \"charts_killed_qc\": %zu,\n"
        "  \"chart_size_hist\": [%zu, %zu, %zu, %zu, %zu, %zu, %zu, %zu, %zu, %zu],\n"
        "  \"chart_qc_hist\": [%zu, %zu, %zu, %zu, %zu, %zu, %zu, %zu, %zu],\n"
        "  \"authentic_neighbour_edges\": %zu,\n"
        "  \"weld_phase_residual_hist\": [%zu, %zu, %zu, %zu, %zu],\n"
        "  \"weld_phase_edges_nonzero\": %zu,\n"
        "  \"overlap_pairs_before\": %zu,\n"
        "  \"overlap_pairs_after\": %zu,\n"
        "  \"intra_overlap_pairs\": %zu,\n"
        "  \"intra_radial_hist\": [%zu, %zu, %zu, %zu, %zu, %zu],\n"
        "  \"groups_overlapping_before\": %zu,\n"
        "  \"groups_overlapping_after\": %zu,\n"
        "  \"settled_groups\": %zu,\n"
        "  \"shifted_groups\": %zu,\n"
        "  \"parked_groups\": %zu,\n"
        "  \"stuck_groups\": %zu,\n"
        "  \"park_v_base\": %.9g,\n"
        "  \"tabu_iterations\": %zu,\n"
        "  \"tabu_moves_chart\": %zu,\n"
        "  \"tabu_moves_group\": %zu,\n"
        "  \"tabu_moves_nbhd\": %zu,\n"
        "  \"tabu_moves_graded\": %zu,\n"
        "  \"tabu_moves_region\": %zu,\n"
        "  \"tabu_graded_groups\": %zu,\n"
        "  \"tabu_phase_regions\": %zu,\n"
        "  \"tabu_charts_moved\": %zu,\n"
        "  \"tabu_active_charts\": %zu,\n"
        "  \"tabu_welds_torn\": %zu,\n"
        "  \"tabu_edges_excluded\": %zu,\n"
        "  \"tabu_best_iter\": %zu,\n"
        "  \"tabu_archive_states\": %zu,\n"
        "  \"tabu_leaderboard\": %zu,\n"
        "  \"tabu_winner\": %zu,\n"
        "  \"tabu_winner_iter\": %zu,\n"
        "  \"tabu_exact_score\": %.9g,\n"
        "  \"tabu_exact_queries\": %zu,\n"
        "  \"tabu_exact_aabb_rejects\": %zu,\n"
        "  \"tabu_exact_cache_hits\": %zu,\n"
        "  \"tabu_exact_cache_misses\": %zu,\n"
        "  \"tabu_exact_cache_entries\": %zu,\n"
        "  \"tabu_exact_sat_tests\": %zu,\n"
        "  \"tabu_energy_initial\": %.9g,\n"
        "  \"tabu_energy_final\": %.9g,\n"
        "  \"tabu_ladder_step\": %.17g,\n"
        "  \"tabu_ladder_u0\": %.17g,\n"
        "  \"tabu_ladder_fixed\": %d,\n"
        "  \"tabu_ladder_pairs\": %zu,\n"
        "  \"tabu_phase_gauge_islands\": %zu,\n"
        "  \"tabu_phase_flat_groups\": %zu,\n"
        "  \"tabu_phase_gauge_inconsistent\": %zu,\n"
        "  \"tabu_lateral_edges\": %zu,\n"
        "  \"tabu_phase_lateral_trusted\": %zu,\n"
        "  \"tabu_phase_lateral_nonzero\": %zu,\n"
        "  \"tabu_lateral_torn\": %zu,\n"
        "  \"tabu_family_edges\": %zu,\n"
        "  \"tabu_family_torn\": %zu,\n"
        "  \"intra_overlap_pairs_after\": %zu,\n"
        "  \"intra_radial_hist_after\": [%zu, %zu, %zu, %zu, %zu, %zu],\n"
        "  \"relayout_edges\": [%zu, %zu, %zu],\n"
        "  \"relayout_weld_rms\": [%.9g, %.9g],\n"
        "  \"relayout_radial_rms\": [%.9g, %.9g],\n"
        "  \"relayout_metric_shift_rms\": %.9g,\n"
        "  \"relayout_metric_shift_max\": %.9g,\n"
        "  \"relayout_shift_rms\": %.9g,\n"
        "  \"relayout_shift_max\": %.9g,\n"
        "  \"relayout_atlas_width\": [%.9g, %.9g],\n"
        "  \"relayout_moved_charts\": %zu,\n"
        "  \"relayout_exact_pairs\": [%zu, %zu],\n"
        "  \"relayout_cross_pairs\": [%zu, %zu],\n"
        "  \"relayout_allowed_pairs\": [%zu, %zu],\n"
        "  \"relayout_hard_pairs\": [%zu, %zu],\n"
        "  \"relayout_topology_limited_pairs\": [%zu, %zu],\n"
        "  \"relayout_total_bounds\": %zu,\n"
        "  \"relayout_topology_bounds\": %zu,\n"
        "  \"relayout_collision_bounds\": %zu,\n"
        "  \"relayout_collision_bounds_added\": %zu,\n"
        "  \"relayout_collision_bounds_rejected\": %zu,\n"
        "  \"relayout_outer_iterations\": %zu,\n"
        "  \"relayout_qp_iterations\": %d,\n"
        "  \"relayout_qp_active_bounds\": %d,\n"
        "  \"relayout_qp_objective\": %.9g,\n"
        "  \"relayout_qp_stationarity\": %.9g,\n"
        "  \"relayout_qp_max_violation\": %.9g,\n"
        "  \"relayout_failed\": %d,\n"
        "  \"relayout_topology_infeasible\": %d,\n"
        "  \"relayout_collision_failed\": %d,\n"
        "  \"relayout_reverted\": %d,\n"
        "  \"seconds\": %.3f\n"
        "}\n",
        s->total_charts, s->charts_before_cull, s->total_groups,
        s->largest_group_charts, s->total_faces, s->faces_killed,
        s->faces_killed_stretch, s->faces_killed_charts,
        s->charts_killed_small, s->charts_killed_qc,
        s->chart_size_hist[0], s->chart_size_hist[1],
        s->chart_size_hist[2], s->chart_size_hist[3],
        s->chart_size_hist[4], s->chart_size_hist[5],
        s->chart_size_hist[6], s->chart_size_hist[7],
        s->chart_size_hist[8], s->chart_size_hist[9],
        s->chart_qc_hist[0], s->chart_qc_hist[1], s->chart_qc_hist[2],
        s->chart_qc_hist[3], s->chart_qc_hist[4], s->chart_qc_hist[5],
        s->chart_qc_hist[6], s->chart_qc_hist[7], s->chart_qc_hist[8],
        s->authentic_neighbour_edges,
        s->weld_phase_residual_hist[0], s->weld_phase_residual_hist[1],
        s->weld_phase_residual_hist[2], s->weld_phase_residual_hist[3],
        s->weld_phase_residual_hist[4], s->weld_phase_edges_nonzero,
        s->overlap_pairs_before, s->overlap_pairs_after,
        s->intra_overlap_pairs,
        s->intra_radial_hist[0], s->intra_radial_hist[1],
        s->intra_radial_hist[2], s->intra_radial_hist[3],
        s->intra_radial_hist[4], s->intra_radial_hist[5],
        s->groups_overlapping_before,
        s->groups_overlapping_after, s->settled_groups, s->shifted_groups,
        s->parked_groups, s->stuck_groups, s->park_v_base,
        s->tabu_iterations, s->tabu_moves_chart, s->tabu_moves_group,
        s->tabu_moves_nbhd,
        s->tabu_moves_graded, s->tabu_moves_region,
        s->tabu_graded_groups, s->tabu_phase_regions,
        s->tabu_charts_moved, s->tabu_active_charts,
        s->tabu_welds_torn, s->tabu_edges_excluded,
        s->tabu_best_iter, s->tabu_archive_states, s->tabu_leaderboard,
        s->tabu_winner,
        s->tabu_winner_iter, s->tabu_exact_score,
        s->tabu_exact_queries, s->tabu_exact_aabb_rejects,
        s->tabu_exact_cache_hits, s->tabu_exact_cache_misses,
        s->tabu_exact_cache_entries, s->tabu_exact_sat_tests,
        s->tabu_energy_initial, s->tabu_energy_final,
        s->tabu_ladder_step, s->tabu_ladder_u0, s->tabu_ladder_fixed,
        s->tabu_ladder_pairs,
        s->tabu_phase_gauge_islands, s->tabu_phase_flat_groups,
        s->tabu_phase_gauge_inconsistent,
        s->tabu_lateral_edges, s->tabu_phase_lateral_trusted,
        s->tabu_phase_lateral_nonzero, s->tabu_lateral_torn,
        s->tabu_family_edges, s->tabu_family_torn,
        s->intra_overlap_pairs_after,
        s->intra_radial_hist_after[0], s->intra_radial_hist_after[1],
        s->intra_radial_hist_after[2], s->intra_radial_hist_after[3],
        s->intra_radial_hist_after[4], s->intra_radial_hist_after[5],
        s->relayout_edges_weld, s->relayout_edges_lateral,
        s->relayout_edges_radial,
        s->relayout_weld_rms_before, s->relayout_weld_rms_after,
        s->relayout_radial_rms_before, s->relayout_radial_rms_after,
        s->relayout_metric_shift_rms, s->relayout_metric_shift_max,
        s->relayout_shift_rms, s->relayout_shift_max,
        s->relayout_atlas_width_before, s->relayout_atlas_width_after,
        s->relayout_moved_charts,
        s->relayout_exact_pairs_before, s->relayout_exact_pairs_after,
        s->relayout_cross_pairs_before, s->relayout_cross_pairs_after,
        s->relayout_allowed_pairs_before,
        s->relayout_allowed_pairs_after,
        s->relayout_hard_pairs_before,
        s->relayout_hard_pairs_after,
        s->relayout_topology_limited_before,
        s->relayout_topology_limited_after,
        s->relayout_total_bounds, s->relayout_topology_bounds,
        s->relayout_collision_bounds, s->relayout_collision_bounds_added,
        s->relayout_collision_bounds_rejected,
        s->relayout_outer_iterations, s->relayout_qp_iterations,
        s->relayout_qp_active_bounds, s->relayout_qp_objective,
        s->relayout_qp_stationarity, s->relayout_qp_max_violation,
        s->relayout_failed, s->relayout_topology_infeasible,
        s->relayout_collision_failed,
        s->relayout_reverted,
        seconds);
    fclose(fp);
    return 0;
}

static void usage(const char *exe)
{
    fprintf(stderr,
        "Usage: %s <placed_dir> <out_dir> [options]\n"
        "       %s --selftest\n\n"
        "  --production         use the validated atlas recipe: exact #58\n"
        "                       settings on small inputs, bounded historical\n"
        "                       settings on production-size inputs\n"
        "  --aspect-kill X      UV/XYZ edge stretch above which a face dies (8)\n"
        "  --floater-faces N    kill whole charts with at most N faces before\n"
        "                       fitting (0 disables)\n"
        "  --floater-group-faces N\n"
        "                       kill welded groups with at most N total faces\n"
        "  --chart-qc-kill X    kill a chart if any face's projected Sander\n"
        "                       sigma1/sigma2 exceeds X (0 disables)\n"
        "  --max-shift N        search +/-1 .. +/-N wraps (4)\n"
        "  --neighbour-dist X   XYZ proximity for an authentic weld (6)\n"
        "  --neighbour-dot X    |dot| of contact normals to accept a weld (0.7)\n"
        "  --neighbour-shared N close vertex pairs required to accept it (3)\n"
        "  --park-margin X      blank v gap around a parked group (16)\n"
        "  --no-park            leave unplaceable groups where they are\n"
        "  --no-wind-correct    skip the pre-tabu weld wind reconciliation\n"
        "  --wind-only          globally reconcile weld gauges, emit UV\n"
        "                       sidecars, and skip every overlap scan/search\n"
        "  --phase-targets      use source phase, not current u, for relative\n"
        "                       weld depth targets (controlled ablation)\n"
        "  --graded-targets     use weld-chain winding differences as targets\n"
        "                       inside geometrically graded groups\n"
        "  --phase-graded       derive graded bins from source phase minus the\n"
        "                       winding already present in the incoming atlas\n"
        "  --phase-all-groups   apply phase grading to every trusted group\n"
        "  --phase-flat-lock    keep a trusted one-bin welded group rigid; no\n"
        "                       chart or mixed-neighbourhood move may split it\n"
        "  --phase-lateral-targets\n"
        "                       derive lateral relative gauges from source phase\n"
        "  --phase-regions      preserve consensus relative gauges with one\n"
        "                       compound move per connected group region\n"
        "  --lock-graded        after a graded move, allow only rigid or graded\n"
        "                       moves for that group\n"
        "  --postgraded-refine  protect a graded group until its first graded\n"
        "                       move, then allow exact chart/neighbourhood\n"
        "                       cleanup (implies --lock-graded)\n"
        "\n"
        "  --tabu               per-CHART winding-depth tabu search instead of\n"
        "                       the rigid greedy group placement\n"
        "  --tabu-span N        depth range +/-N wraps (4)\n"
        "  --tabu-tenure N      iterations a departed depth stays banned (40)\n"
        "  --tabu-iters N       hard iteration cap (2000)\n"
        "  --tabu-stall N       stop after N non-improving iterations (300)\n"
        "                       inputs >=400 cubes or >=750k faces default to\n"
        "                       one iteration/round; --tabu-iters opts out\n"
        "  --tabu-cell X        legacy occupancy-audit cell size, vox (4)\n"
        "  --tabu-ladder-step X reuse a global signed field-ladder scale\n"
        "  --tabu-ladder-u0 X   reuse its global u intercept; supply both\n"
        "  --active-core Z0 Z1 Y0 Y1 X0 X1\n"
        "                       move only cubes in this half-open origin box;\n"
        "                       cubes elsewhere in the loaded tile are halo\n"
        "  --write-uvphi        emit solved per-cube sidecars for core commit\n"
        "  --lambda-ov X        energy per exact overlapping face pair (1)\n"
        "  --lambda-nb X        energy per contact-weighted weld tear (4)\n"
        "  --lambda-rad X       energy per turn of radial deviation (2)\n"
        "  --lambda-gauge X     quadratic common-gauge retention (0)\n"
        "  --nb-unhappy-w X     steering weight of edges to a colliding\n"
        "                       neighbour (0.25)\n"
        "  --no-happy           disable the happiness steering entirely\n"
        "  --lateral-w X        contact-equivalent weight coupling same-wind\n"
        "                       non-welded neighbours (8; 0 disables)\n"
        "  --lateral-gap X      max |du| to call charts side-by-side (200)\n"
        "  --family-w X         boost on edges inside a happy family -- a\n"
        "                       group with zero internal overlap (8; <=1 off)\n"
        "\n"
        "  --relayout           after tabu, re-register continuous per-chart\n"
        "                       u shifts with the windings held fixed\n"
        "  --no-relayout        disable relayout selected by a preset\n"
        "  --relayout-seam-cap X\n"
        "                       max per-edge deviation from the coherent metric\n"
        "                       solution (25; 0 disables topology priority)\n"
        "  --rounds N           alternate tabu <-> relayout N times, feeding\n"
        "                       each round's field to the next; stops early\n"
        "                       when a round moves no charts (implies both\n"
        "                       --tabu and --relayout)\n",
        exe, exe);
}

static int production_size(size_t n_cubes, size_t nfaces)
{
    return n_cubes >= 400 || nfaces >= 750000;
}

/* The July 29 reference atlas (overlap_fix_ablation/58_ribbon_checkpoint).
 * Keep this as a named recipe rather than silently relying on nineteen
 * independent command-line defaults. */
static void apply_reference58_preset(AtlasOverlapFixOptions *opts)
{
    opts->wind_correct = 0;
    opts->floater_max_faces = 6;
    opts->chart_qc_kill = 50.0;
    opts->floater_group_max_faces = 512;
    opts->phase_graded = 1;
    opts->graded_targets = 1;
    opts->phase_flat_lock = 1;
    opts->phase_lateral_targets = 1;
    opts->phase_regions = 1;
    opts->lock_graded = 1;
    opts->tabu = 1;
    opts->tabu_span = 16;
    opts->tabu_max_iters = 80;
    opts->tabu_stall = 80;
    opts->lambda_gauge = 4.0;
    opts->relayout = 1;
    opts->relayout_seam_cap = 25.0;
}

/* The last legitimate 10x10x10 atlas used this deliberately bounded solve.
 * One iteration took 45.6 s on 996 cubes; applying the small-fixture 80x
 * budget (or the later seven-round wrapper) is the apparent infinite run. */
static void apply_large_production_preset(AtlasOverlapFixOptions *opts)
{
    opts->floater_max_faces = 6;
    opts->chart_qc_kill = 50.0;
    opts->tabu = 1;
    opts->tabu_span = 4;
    opts->tabu_max_iters = 1;
    opts->tabu_stall = 1;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0)
        return AtlasOverlapFix_selftest() == 0 &&
               AtlasRegister_selftest() == 0 ? 0 : 1;
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    const char *placed_dir = argv[1];
    const char *out_dir = argv[2];
    int rounds = 1;
    int tabu_iters_explicit = 0;
    int production_preset = 0;
    int tabu_ladder_step_explicit = 0;
    int tabu_ladder_u0_explicit = 0;
    int active_core_requested = 0;
    int write_uvphi = 0;
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--production") == 0) production_preset = 1;
        if (strcmp(argv[i], "--active-core") == 0) active_core_requested = 1;
    }

    Arena_T arena = Arena_new();
    PieceSet ps;
    if (PieceSet_build(arena, placed_dir, &ps) != 0) {
        fprintf(stderr, "%s: cannot load placed dir %s\n", argv[0], placed_dir);
        Arena_dispose(&arena);
        return 1;
    }

    AtlasOverlapFixOptions opts;
    AtlasOverlapFixOptions_default(&opts);
    if (production_preset) {
        if (active_core_requested || production_size(ps.n_cubes, ps.nf))
            apply_large_production_preset(&opts);
        else
            apply_reference58_preset(&opts);
    }

    for (int i = 3; i < argc; i++) {
        char *end = NULL;
        if (strcmp(argv[i], "--production") == 0) {
            continue;
        } else if (strcmp(argv[i], "--aspect-kill") == 0 && i + 1 < argc) {
            opts.aspect_ratio_kill = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--floater-faces") == 0 && i + 1 < argc) {
            opts.floater_max_faces = (size_t)strtoul(argv[++i], &end, 10);
        } else if (strcmp(argv[i], "--floater-group-faces") == 0 && i + 1 < argc) {
            opts.floater_group_max_faces =
                (size_t)strtoul(argv[++i], &end, 10);
        } else if (strcmp(argv[i], "--chart-qc-kill") == 0 && i + 1 < argc) {
            opts.chart_qc_kill = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--max-shift") == 0 && i + 1 < argc) {
            opts.max_shift_wraps = (int)strtol(argv[++i], &end, 10);
        } else if (strcmp(argv[i], "--neighbour-dist") == 0 && i + 1 < argc) {
            opts.neighbour_distance = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--neighbour-dot") == 0 && i + 1 < argc) {
            opts.neighbour_normal_dot = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--neighbour-shared") == 0 && i + 1 < argc) {
            opts.neighbour_min_shared = (size_t)strtoul(argv[++i], &end, 10);
        } else if (strcmp(argv[i], "--park-margin") == 0 && i + 1 < argc) {
            opts.park_margin = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--no-park") == 0) {
            opts.park_unplaceable = 0;
            continue;
        } else if (strcmp(argv[i], "--no-wind-correct") == 0) {
            opts.wind_correct = 0;
            continue;
        } else if (strcmp(argv[i], "--wind-only") == 0) {
            opts.wind_only = 1;
            continue;
        } else if (strcmp(argv[i], "--phase-targets") == 0) {
            opts.phase_targets = 1;
            continue;
        } else if (strcmp(argv[i], "--graded-targets") == 0) {
            opts.graded_targets = 1;
            continue;
        } else if (strcmp(argv[i], "--phase-graded") == 0) {
            opts.phase_graded = 1;
            continue;
        } else if (strcmp(argv[i], "--phase-all-groups") == 0) {
            opts.phase_all_groups = 1;
            continue;
        } else if (strcmp(argv[i], "--phase-flat-lock") == 0) {
            opts.phase_flat_lock = 1;
            continue;
        } else if (strcmp(argv[i], "--phase-lateral-targets") == 0) {
            opts.phase_lateral_targets = 1;
            continue;
        } else if (strcmp(argv[i], "--phase-regions") == 0) {
            opts.phase_regions = 1;
            continue;
        } else if (strcmp(argv[i], "--lock-graded") == 0) {
            opts.lock_graded = 1;
            continue;
        } else if (strcmp(argv[i], "--postgraded-refine") == 0) {
            opts.lock_graded = 1;
            opts.postgraded_refine = 1;
            continue;
        } else if (strcmp(argv[i], "--tabu") == 0) {
            opts.tabu = 1;
            continue;
        } else if (strcmp(argv[i], "--tabu-span") == 0 && i + 1 < argc) {
            opts.tabu_span = (int)strtol(argv[++i], &end, 10);
        } else if (strcmp(argv[i], "--tabu-tenure") == 0 && i + 1 < argc) {
            opts.tabu_tenure = (int)strtol(argv[++i], &end, 10);
        } else if (strcmp(argv[i], "--tabu-iters") == 0 && i + 1 < argc) {
            opts.tabu_max_iters = (int)strtol(argv[++i], &end, 10);
            tabu_iters_explicit = 1;
        } else if (strcmp(argv[i], "--tabu-stall") == 0 && i + 1 < argc) {
            opts.tabu_stall = (int)strtol(argv[++i], &end, 10);
        } else if (strcmp(argv[i], "--tabu-cell") == 0 && i + 1 < argc) {
            opts.tabu_cell = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--tabu-ladder-step") == 0 &&
                   i + 1 < argc) {
            opts.tabu_ladder_step = strtod(argv[++i], &end);
            tabu_ladder_step_explicit = 1;
        } else if (strcmp(argv[i], "--tabu-ladder-u0") == 0 &&
                   i + 1 < argc) {
            opts.tabu_ladder_u0 = strtod(argv[++i], &end);
            tabu_ladder_u0_explicit = 1;
        } else if (strcmp(argv[i], "--active-core") == 0 && i + 6 < argc) {
            long box[6];
            int bad = 0;
            for (int q = 0; q < 6; q++)
                if (parse_long_full(argv[++i], &box[q]) != 0) bad = 1;
            if (bad) {
                fprintf(stderr, "%s: bad --active-core bound\n", argv[0]);
                Arena_dispose(&arena);
                return 1;
            }
            opts.tabu_active_core = 1;
            opts.tabu_core_z_lo = box[0]; opts.tabu_core_z_hi = box[1];
            opts.tabu_core_y_lo = box[2]; opts.tabu_core_y_hi = box[3];
            opts.tabu_core_x_lo = box[4]; opts.tabu_core_x_hi = box[5];
            continue;
        } else if (strcmp(argv[i], "--write-uvphi") == 0) {
            write_uvphi = 1;
            continue;
        } else if (strcmp(argv[i], "--lambda-ov") == 0 && i + 1 < argc) {
            opts.lambda_ov = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--lambda-nb") == 0 && i + 1 < argc) {
            opts.lambda_nb = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--lambda-rad") == 0 && i + 1 < argc) {
            opts.lambda_rad = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--lambda-gauge") == 0 && i + 1 < argc) {
            opts.lambda_gauge = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--nb-unhappy-w") == 0 && i + 1 < argc) {
            opts.nb_unhappy_w = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--no-happy") == 0) {
            opts.nb_unhappy_w = 1.0;
            continue;
        } else if (strcmp(argv[i], "--lateral-w") == 0 && i + 1 < argc) {
            opts.lateral_w = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--lateral-gap") == 0 && i + 1 < argc) {
            opts.lateral_u_gap = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--family-w") == 0 && i + 1 < argc) {
            opts.family_w = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--relayout-seam-cap") == 0 &&
                   i + 1 < argc) {
            opts.relayout_seam_cap = strtod(argv[++i], &end);
        } else if (strcmp(argv[i], "--relayout") == 0) {
            opts.relayout = 1;
            continue;
        } else if (strcmp(argv[i], "--no-relayout") == 0) {
            opts.relayout = 0;
            continue;
        } else if (strcmp(argv[i], "--rounds") == 0 && i + 1 < argc) {
            rounds = (int)strtol(argv[++i], &end, 10);
        } else {
            fprintf(stderr, "%s: unknown or incomplete option %s\n",
                    argv[0], argv[i]);
            Arena_dispose(&arena);
            return 1;
        }
        if (end == argv[i]) {
            fprintf(stderr, "%s: bad value for %s\n", argv[0], argv[i - 1]);
            Arena_dispose(&arena);
            return 1;
        }
    }

    if (tabu_ladder_step_explicit != tabu_ladder_u0_explicit) {
        fprintf(stderr,
                "%s: --tabu-ladder-step and --tabu-ladder-u0 must be "
                "supplied together\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    opts.tabu_ladder_override = tabu_ladder_step_explicit;

    if (opts.wind_only) {
        opts.wind_correct = 1;
        opts.tabu = 0;
        opts.relayout = 0;
        rounds = 1;
        write_uvphi = 1;
    }

    if (opts.tabu_active_core && !opts.tabu) {
        fprintf(stderr, "%s: --active-core requires --tabu (or --production)\n",
                argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    if (opts.tabu_active_core && opts.wind_correct) {
        fprintf(stderr, "[overlap_fix] active core: disabling tile-local "
                        "pre-wind so the halo remains frozen\n");
        opts.wind_correct = 0;
    }
    if (opts.tabu_active_core && opts.relayout) {
        fprintf(stderr, "[overlap_fix] active core: disabling global relayout\n");
        opts.relayout = 0;
    }

    if (rounds < 1) rounds = 1;
    if (rounds > 1) {
        /* --rounds IS the alternating loop, so it turns both halves on itself
         * rather than making the caller know that the loop is tabu plus
         * relayout underneath. This is what the usage text has always
         * promised; requiring an explicit --tabu alongside it only ever
         * produced an error for a request that was already unambiguous. */
        opts.tabu = 1;
        opts.relayout = 1;
    }

    char probe[AOF_PATH_CAP];
    if (snprintf(probe, AOF_PATH_CAP, "%s/.", out_dir) >= AOF_PATH_CAP ||
        ves_ensure_parent_dir(probe) != 0) {
        fprintf(stderr, "%s: cannot create %s\n", argv[0], out_dir);
        Arena_dispose(&arena);
        return 1;
    }

    if (production_preset)
        fprintf(stderr, "[overlap_fix] production preset: %s\n",
                active_core_requested || production_size(ps.n_cubes, ps.nf)
                    ? "large-grid bounded" : "reference-58");

    /* Candidate evaluation is super-linear in chart count.  The old 10^3
     * validation deliberately used one tabu iteration and finished in under
     * a minute; asking the small-atlas 2000x7 budget of the same input looks
     * like a hang (about 55 s/iteration on PHerc0139).  Keep the exhaustive
     * small-fixture search, but use the proven bounded policy by default for
     * production-size grids.  Supplying --tabu-iters is an explicit opt-out. */
    if (opts.tabu && !tabu_iters_explicit &&
        (opts.tabu_active_core || production_size(ps.n_cubes, ps.nf))) {
        opts.tabu_max_iters = 1;
        opts.tabu_stall = 1;
        if (rounds > 1) rounds = 1;
        fprintf(stderr,
            "[overlap_fix] large-input guard: %zu cubes/%zu faces; "
            "tabu_iters=1 rounds=1 (pass --tabu-iters N to override)\n",
            ps.n_cubes, ps.nf);
    }
    if (opts.tabu_active_core)
        fprintf(stderr,
                "[overlap_fix] active core: z[%ld,%ld) y[%ld,%ld) x[%ld,%ld); "
                "loaded complement is fixed halo\n",
                opts.tabu_core_z_lo, opts.tabu_core_z_hi,
                opts.tabu_core_y_lo, opts.tabu_core_y_hi,
                opts.tabu_core_x_lo, opts.tabu_core_x_hi);
    ScaffoldCalib cal;
    if (Scaffold_read_calib(placed_dir, &cal) != 0) {
        fprintf(stderr, "%s: missing placed calibration\n", argv[0]);
        Arena_dispose(&arena);
        return 1;
    }
    float *registered_u = (float *)ARENA_ALLOC(arena, ps.nv * sizeof(float));
    for (size_t i = 0; i < ps.nv; i++) {
        registered_u[i] = ps.uv[i * 2];
        if (!isfinite(registered_u[i])) {
            fprintf(stderr, "%s: non-finite registered u at vertex %zu\n",
                    argv[0], i);
            Arena_dispose(&arena);
            return 1;
        }
    }
    fprintf(stderr, "[overlap_fix] %zu cubes, %zu vertices, %zu faces, "
                    "pitch=%.4g\n", ps.n_cubes, ps.nv, ps.nf, cal.pitch);

    double *out_u = (double *)ARENA_ALLOC(arena, ps.nv * sizeof(double));
    double *out_v = (double *)ARENA_ALLOC(arena, ps.nv * sizeof(double));
    double *before_u = (double *)ARENA_ALLOC(arena, ps.nv * sizeof(double));
    double *before_v = (double *)ARENA_ALLOC(arena, ps.nv * sizeof(double));
    uint8_t *face_keep = (uint8_t *)ARENA_ALLOC(arena, ps.nf * sizeof(uint8_t));
    int32_t *vertex_group = (int32_t *)ARENA_ALLOC(
        arena, ps.nv * sizeof(int32_t));
    AtlasOverlapFixGroup *group = NULL;
    size_t ngroups = 0;
    AtlasOverlapFixTabuReport tabu_report;
    AtlasOverlapFixStats stats;
    memset(&tabu_report, 0, sizeof(tabu_report));
    memset(&stats, 0, sizeof(stats));

    float *u_iter = (float *)ARENA_ALLOC(arena, ps.nv * sizeof(float));
    memcpy(u_iter, registered_u, ps.nv * sizeof(float));

    double t0 = ves_clock_sec();
    int rounds_run = 0;
    for (int round = 1; round <= rounds; round++) {
        if (AtlasOverlapFix_solve(arena, &ps, &cal, u_iter, &opts,
                                  out_u, out_v, face_keep, vertex_group,
                                  &group, &ngroups, &tabu_report,
                                  &stats) != 0) {
            fprintf(stderr, "%s: solve failed (round %d)\n", argv[0], round);
            Arena_dispose(&arena);
            return 1;
        }
        rounds_run = round;
        if (opts.tabu && tabu_report.nchart > 0) {
            write_layer_objs(out_dir, round, opts.tabu_span, &ps,
                             tabu_report.u_post_tabu, out_v, face_keep,
                             vertex_group, &tabu_report);
            char cname[64];
            snprintf(cname, sizeof(cname), "round%d_charts.csv", round);
            write_charts_csv(out_dir, cname, &ps, &tabu_report);
        }
        if (rounds > 1) {
            fprintf(stderr,
                "[overlap_fix] round %d: moved %zu charts, torn %zu welds / "
                "%zu lateral, SAT %zu cross + %zu intra%s\n",
                round, stats.tabu_charts_moved, stats.tabu_welds_torn,
                stats.tabu_lateral_torn, stats.overlap_pairs_after,
                stats.intra_overlap_pairs_after,
                opts.relayout && stats.relayout_failed
                    ? "  (RELAYOUT FAILED)" : "");
        }
        if (round < rounds) {
            char round_name[64];
            snprintf(round_name, sizeof(round_name), "round%d_uv.obj", round);
            write_uv_obj(out_dir, round_name, &ps, out_u, out_v, face_keep,
                         vertex_group);
            for (size_t i = 0; i < ps.nv; i++) u_iter[i] = (float)out_u[i];
            /* The wind-correct pre-pass re-measures welds each round; on a
             * field where an earlier round DELIBERATELY tore some, it would
             * drag those subtrees a whole wind back before tabu runs.  Later
             * rounds run without it: the edge targets then encode standing
             * tears and preserve prior decisions. */
            opts.wind_correct = 0;
        }
        if (opts.tabu && stats.tabu_charts_moved == 0 && round > 1) {
            fprintf(stderr, "[overlap_fix] converged after round %d "
                            "(no charts moved)\n", round);
            break;
        }
    }
    double seconds = ves_clock_sec() - t0;
    if (opts.wind_only) {
        int io = write_stats_json(out_dir, &stats, seconds);
        io |= write_uvphi_sidecars(out_dir, &ps, out_u, out_v);
        if (io != 0) {
            fprintf(stderr, "%s: cannot write wind-only outputs under %s\n",
                    argv[0], out_dir);
            Arena_dispose(&arena);
            return 1;
        }
        fprintf(stderr,
                "[overlap_fix] wind-only: corrected %zu/%zu charts across "
                "%zu gauge islands in %.2fs\n",
                stats.charts_wind_corrected, stats.total_charts,
                stats.gauge_islands, seconds);
        Arena_dispose(&arena);
        return 0;
    }
    if (rounds > 1)
        fprintf(stderr, "[overlap_fix] %d round(s), %.1fs total\n",
                rounds_run, seconds);

    /* "Before" reuses the solve's own kill mask and grouping so the two
     * pictures differ only by the placement, which is the thing under test. */
    for (size_t i = 0; i < ps.nv; i++) {
        before_u[i] = (double)registered_u[i];
        before_v[i] = out_v[i];
    }
    for (size_t g = 0; g < ngroups; g++) {
        if (group[g].shift_v == 0.0) continue;
        for (size_t i = 0; i < ps.nv; i++)
            if (vertex_group[i] == group[g].group)
                before_v[i] -= group[g].shift_v;
    }

    int io = 0;
    io |= write_uv_obj(out_dir, "before_uv.obj", &ps, before_u, before_v,
                       face_keep, vertex_group);
    io |= write_uv_obj(out_dir, "after_uv.obj", &ps, out_u, out_v,
                       face_keep, vertex_group);
    io |= write_xyz_obj(out_dir, "groups_xyz.obj", &ps, face_keep,
                        vertex_group);
    io |= write_groups_csv(out_dir, group, ngroups);
    io |= write_stats_json(out_dir, &stats, seconds);
    if (write_uvphi)
        io |= write_uvphi_sidecars(out_dir, &ps, out_u, out_v);
    if (opts.tabu && tabu_report.nchart > 0) {
        char solution_path[AOF_PATH_CAP];
        int solution_name = snprintf(solution_path, sizeof solution_path,
                                     "%s/atlas_solution.bin", out_dir);
        if (solution_name < 0 ||
            (size_t)solution_name >= sizeof solution_path)
            io = -1;
        else
            io |= AtlasSolution_write(solution_path, &ps, out_u, out_v,
                                      face_keep, &tabu_report);
        io |= write_depth_obj(out_dir, &ps, out_u, out_v, face_keep,
                              &tabu_report);
        io |= write_charts_csv(out_dir, "charts.csv", &ps, &tabu_report);
        io |= write_moves_csv(out_dir, &tabu_report);
        io |= write_chart_moves_csv(out_dir, &tabu_report);
        io |= write_candidates_csv(out_dir, &tabu_report);
        io |= write_edges_csv(out_dir, &tabu_report);
        io |= write_bake_obj(out_dir, &ps, out_u, out_v, face_keep);
        io |= write_diag_png(out_dir, "diag_before.png", &ps, &cal,
                             before_u, vertex_group);
        io |= write_overlap_pairs_csv(out_dir, &tabu_report);
        if (opts.relayout) {
            io |= write_relayout_bounds_csv(out_dir, &tabu_report);
            io |= write_relayout_residual_csv(out_dir, &tabu_report);
            io |= write_relayout_residual_chart_pairs_csv(out_dir, &tabu_report);
        }
        io |= write_diag_png(out_dir, "diag_after.png", &ps, &cal,
                             out_u, vertex_group);
    }
    if (io != 0) {
        fprintf(stderr, "%s: cannot write outputs under %s\n", argv[0],
                out_dir);
        Arena_dispose(&arena);
        return 1;
    }

    for (size_t g = 0; g < ngroups; g++) {
        if (group[g].action == ATLAS_OVERLAP_FIX_SETTLED) continue;
        fprintf(stderr,
            "[overlap_fix]   group %d: %s, %zu charts, r=%.1f (spread %.1f), "
            "wraps=%d du=%.1f dv=%.1f, overlap %zu -> %zu\n",
            group[g].group, action_name(group[g].action), group[g].ncharts,
            group[g].mean_radius, group[g].radius_spread, group[g].wraps,
            group[g].shift_u, group[g].shift_v, group[g].overlap_before,
            group[g].overlap_after);
    }

    size_t hist_total = 0;
    for (int b = 0; b < 5; b++) hist_total += stats.weld_residual_hist[b];
    double clean_fraction = hist_total > 0
        ? 100.0 * (double)stats.weld_residual_hist[0] / (double)hist_total
        : 0.0;
    double phase_clean_fraction = hist_total > 0
        ? 100.0 * (double)stats.weld_phase_residual_hist[0] / (double)hist_total
        : 0.0;

    fprintf(stderr,
        "\n=== PREFIT COMPONENT CULL ===\n"
        "  charts before / after ... %zu / %zu\n"
        "  size hist (faces) ........ 1:%zu 2-3:%zu 4-7:%zu 8-15:%zu "
        "16-31:%zu 32-63:%zu 64-127:%zu 128-255:%zu 256-511:%zu 512+:%zu\n"
        "  max Sander QC hist ....... <=2:%zu <=4:%zu <=8:%zu <=16:%zu "
        "<=32:%zu <=64:%zu <=128:%zu <=256:%zu >256:%zu\n"
        "  charts killed ............ %zu small, %zu projected-degenerate\n"
        "  welded groups before/after %zu / %zu\n"
        "  group size hist .......... <=1:%zu <=3:%zu <=7:%zu <=15:%zu "
        "<=31:%zu <=63:%zu <=127:%zu <=255:%zu <=511:%zu <=1023:%zu "
        "<=2047:%zu >2047:%zu\n"
        "  groups killed ............ %zu small\n"
        "  faces killed ............. %zu stretch + %zu chart + %zu group\n"
        "==========================\n",
        stats.charts_before_cull, stats.total_charts,
        stats.chart_size_hist[0], stats.chart_size_hist[1],
        stats.chart_size_hist[2], stats.chart_size_hist[3],
        stats.chart_size_hist[4], stats.chart_size_hist[5],
        stats.chart_size_hist[6], stats.chart_size_hist[7],
        stats.chart_size_hist[8], stats.chart_size_hist[9],
        stats.chart_qc_hist[0], stats.chart_qc_hist[1],
        stats.chart_qc_hist[2], stats.chart_qc_hist[3],
        stats.chart_qc_hist[4], stats.chart_qc_hist[5],
        stats.chart_qc_hist[6], stats.chart_qc_hist[7],
        stats.chart_qc_hist[8], stats.charts_killed_small,
        stats.charts_killed_qc, stats.groups_before_cull, stats.total_groups,
        stats.group_size_hist[0], stats.group_size_hist[1],
        stats.group_size_hist[2], stats.group_size_hist[3],
        stats.group_size_hist[4], stats.group_size_hist[5],
        stats.group_size_hist[6], stats.group_size_hist[7],
        stats.group_size_hist[8], stats.group_size_hist[9],
        stats.group_size_hist[10], stats.group_size_hist[11],
        stats.groups_killed_small, stats.faces_killed_stretch,
        stats.faces_killed_charts, stats.faces_killed_groups);

    fprintf(stderr,
        "\n=== OVERLAP FIX ===\n"
        "  charts .................. %zu\n"
        "  welded groups ........... %zu (largest %zu charts)\n"
        "  authentic weld edges .... %zu\n"
        "  faces killed (total) .... %zu / %zu\n"
        "\n  -- is the error a whole number of turns? --\n"
        "  weld |du/circ - round| .. <0.05:%zu 0.15:%zu 0.25:%zu 0.35:%zu "
        "more:%zu  (%.1f%% clean)\n"
        "  welds off by >=1 wrap ... %zu / %zu\n"
        "  phase |turn - round| ..... <0.05:%zu 0.15:%zu 0.25:%zu 0.35:%zu "
        "more:%zu  (%.1f%% clean)\n"
        "  nonzero phase targets .... %zu / %zu\n"
        "  trusted / untrusted ..... %zu / %zu (+%zu over the wrap cap)\n"
        "  wrap range .............. %d .. %d\n"
        "\n  -- gauge reconciliation --\n"
        "  charts wind-corrected ... %zu / %zu\n"
        "  gauge islands ........... %zu (tree %zu, cycle %zu)\n"
        "  cycle consistent ........ %zu\n"
        "  cycle INCONSISTENT ...... %zu  (holonomy; no gauge can fix these)\n"
        "\n  -- overlap --\n"
        "  cross-group ............. %zu raw -> %zu post-wind -> %zu final\n"
        "  intra-group ............. %zu raw -> %zu post-wind\n"
        "  intra |dr|/pitch ........ ~0:%zu ~0.5:%zu ~1:%zu ~1.5:%zu ~2:%zu "
        ">2.5:%zu\n"
        "    (~0..0.5 = folds no gauge fixes; ~1+ = wrong-wrap, re-gauge CAN fix)\n"
        "  groups overlapping ...... %zu -> %zu\n"
        "  settled / shifted ....... %zu / %zu\n"
        "  parked / stuck .......... %zu / %zu\n"
        "  time .................... %.2fs\n"
        "===================\n",
        stats.total_charts, stats.total_groups, stats.largest_group_charts,
        stats.authentic_neighbour_edges, stats.faces_killed, stats.total_faces,
        stats.weld_residual_hist[0], stats.weld_residual_hist[1],
        stats.weld_residual_hist[2], stats.weld_residual_hist[3],
        stats.weld_residual_hist[4], clean_fraction,
        stats.weld_edges_nonzero_wrap, stats.authentic_neighbour_edges,
        stats.weld_phase_residual_hist[0],
        stats.weld_phase_residual_hist[1],
        stats.weld_phase_residual_hist[2],
        stats.weld_phase_residual_hist[3],
        stats.weld_phase_residual_hist[4], phase_clean_fraction,
        stats.weld_phase_edges_nonzero, stats.authentic_neighbour_edges,
        stats.weld_edges_trusted, stats.weld_edges_untrusted,
        stats.weld_edges_oversized,
        stats.wind_wraps_min, stats.wind_wraps_max,
        stats.charts_wind_corrected, stats.total_charts,
        stats.gauge_islands, stats.gauge_tree_edges, stats.gauge_cycle_edges,
        stats.gauge_cycle_consistent, stats.gauge_cycle_inconsistent,
        stats.overlap_pairs_raw, stats.overlap_pairs_before,
        stats.overlap_pairs_after,
        stats.intra_overlap_pairs_raw, stats.intra_overlap_pairs,
        stats.intra_radial_hist[0], stats.intra_radial_hist[1],
        stats.intra_radial_hist[2], stats.intra_radial_hist[3],
        stats.intra_radial_hist[4], stats.intra_radial_hist[5],
        stats.groups_overlapping_before, stats.groups_overlapping_after,
        stats.settled_groups, stats.shifted_groups,
        stats.parked_groups, stats.stuck_groups, seconds);

    if (opts.tabu) {
        fprintf(stderr,
            "\n=== TABU RE-GAUGE ===\n"
            "  depth span .............. +/-%d wraps, tenure %d\n"
            "  field ladder ............ %+.9g x 2pi*r per wind "
            "(%zu pairs%s)\n"
            "  ladder u origin ......... %.17g\n"
            "  iterations .............. %zu (canonical best at %zu)\n"
            "  moves ................... %zu chart + %zu group + %zu nbhd "
            "+ %zu graded (of %zu eligible groups)\n"
            "  phase-flat regions ...... %zu%s\n"
            "  phase lateral targets ... %zu trusted, %zu nonzero%s\n"
            "  energy (cell units) ..... %.1f -> %.1f\n"
            "  movable charts .......... %zu / %zu\n"
            "  charts moved ............ %zu / %zu\n"
            "  welds torn .............. %zu / %zu usable (+%zu excluded)\n"
            "  lateral pairs coupled ... %zu (torn %zu)\n"
            "  happy-family bonds ...... %zu (torn %zu)\n"
            "  state archive ........... %zu kept, %zu SAT-verified\n"
            "  exact winner ............ #%zu at iter %zu, score %.1f\n"
            "  hash rebuilds ........... %zu\n"
            "  cross after (exact) ..... %zu pairs\n"
            "  intra after (exact) ..... %zu pairs\n"
            "  intra |dr|/pitch after .. ~0:%zu ~0.5:%zu ~1:%zu ~1.5:%zu "
            "~2:%zu >2.5:%zu\n"
            "=====================\n",
            opts.tabu_span, opts.tabu_tenure,
            stats.tabu_ladder_step, stats.tabu_ladder_pairs,
            stats.tabu_ladder_fixed ? "; FIXED" :
            (stats.tabu_ladder_pairs > 0 ? "" : "; SPIRAL-MAP FALLBACK"),
            stats.tabu_ladder_u0,
            stats.tabu_iterations, stats.tabu_best_iter,
            stats.tabu_moves_chart, stats.tabu_moves_group,
            stats.tabu_moves_nbhd, stats.tabu_moves_graded,
            stats.tabu_graded_groups,
            stats.tabu_phase_flat_groups,
            opts.phase_flat_lock ? " locked" : " measured only",
            stats.tabu_phase_lateral_trusted,
            stats.tabu_phase_lateral_nonzero,
            opts.phase_lateral_targets ? " active" : " measured only",
            stats.tabu_energy_initial, stats.tabu_energy_final,
            stats.tabu_active_charts, stats.total_charts,
            stats.tabu_charts_moved, stats.total_charts,
            stats.tabu_welds_torn,
            stats.authentic_neighbour_edges - stats.tabu_edges_excluded,
            stats.tabu_edges_excluded,
            stats.tabu_lateral_edges, stats.tabu_lateral_torn,
            stats.tabu_family_edges, stats.tabu_family_torn,
            stats.tabu_archive_states, stats.tabu_leaderboard,
            stats.tabu_winner, stats.tabu_winner_iter, stats.tabu_exact_score,
            stats.tabu_hash_rebuilds,
            stats.overlap_pairs_after,
            stats.intra_overlap_pairs_after,
            stats.intra_radial_hist_after[0], stats.intra_radial_hist_after[1],
            stats.intra_radial_hist_after[2], stats.intra_radial_hist_after[3],
            stats.intra_radial_hist_after[4],
            stats.intra_radial_hist_after[5]);
        fprintf(stderr,
            "  exact workload .......... %zu queries; %zu AABB rejects; "
            "%zu cache hits / %zu misses (%zu entries); %zu SAT tests\n",
            stats.tabu_exact_queries, stats.tabu_exact_aabb_rejects,
            stats.tabu_exact_cache_hits, stats.tabu_exact_cache_misses,
            stats.tabu_exact_cache_entries, stats.tabu_exact_sat_tests);
    }

    if (opts.relayout) {
        fprintf(stderr,
            "\n=== COLLISION-AWARE U REGISTRATION (windings fixed) ===\n"
            "  edges ................... %zu weld + %zu lateral + %zu radial\n"
            "  weld rms ................ %.2f -> %.2f vox\n"
            "  radial rms .............. %.1f -> %.1f vox\n"
            "  exact pairs ............. %zu total / %zu cross -> "
                                        "%zu total / %zu cross\n"
            "  allowed coverage ........ %zu -> %zu face pairs\n"
            "  hard collisions ......... %zu -> %zu face pairs\n"
            "  topology-limited ........ %zu -> %zu face pairs\n"
            "  bounds .................. %zu total (%zu topology + %zu collision)\n"
            "  collision admissions .... %zu added / %zu rejected in "
                                        "%zu audit(s)\n"
            "  unconstrained shift ..... %.1f rms / %.1f max vox\n"
            "  accepted shift .......... %zu charts; %.1f rms / %.1f max vox\n"
            "  atlas width ............. %.1f -> %.1f vox (%+.2f%%)\n"
            "  QP ...................... %d iterations; %d active; "
                                        "violation %.3g; stationarity %.3g\n"
            "%s"
            "=========================================================\n",
            stats.relayout_edges_weld, stats.relayout_edges_lateral,
            stats.relayout_edges_radial,
            stats.relayout_weld_rms_before, stats.relayout_weld_rms_after,
            stats.relayout_radial_rms_before, stats.relayout_radial_rms_after,
            stats.relayout_exact_pairs_before,
            stats.relayout_cross_pairs_before,
            stats.relayout_exact_pairs_after,
            stats.relayout_cross_pairs_after,
            stats.relayout_allowed_pairs_before,
            stats.relayout_allowed_pairs_after,
            stats.relayout_hard_pairs_before,
            stats.relayout_hard_pairs_after,
            stats.relayout_topology_limited_before,
            stats.relayout_topology_limited_after,
            stats.relayout_total_bounds, stats.relayout_topology_bounds,
            stats.relayout_collision_bounds,
            stats.relayout_collision_bounds_added,
            stats.relayout_collision_bounds_rejected,
            stats.relayout_outer_iterations,
            stats.relayout_metric_shift_rms, stats.relayout_metric_shift_max,
            stats.relayout_moved_charts,
            stats.relayout_shift_rms, stats.relayout_shift_max,
            stats.relayout_atlas_width_before,
            stats.relayout_atlas_width_after,
            stats.relayout_atlas_width_before > 0.0
                ? 100.0 * (stats.relayout_atlas_width_after /
                           stats.relayout_atlas_width_before - 1.0) : 0.0,
            stats.relayout_qp_iterations, stats.relayout_qp_active_bounds,
            stats.relayout_qp_max_violation, stats.relayout_qp_stationarity,
            stats.relayout_topology_infeasible
                ? "  TOPOLOGY FAILURE: no continuity-feasible initializer\n"
                :
            stats.relayout_collision_failed
                ? "  COLLISION FAILURE: u left as tabu placed it\n"
                : stats.relayout_failed
                ? "  SOLVE FAILED: u left as tabu placed it\n"
                : stats.relayout_reverted
                ? "  REVERTED: independent exact audit rejected result\n" : "");
    }

    Arena_dispose(&arena);
    return 0;
}
