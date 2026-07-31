#ifndef GROUP_GRAPH_INCLUDED
#define GROUP_GRAPH_INCLUDED

#include <stdint.h>
#include <stddef.h>
#include "../common/arena.h"
#include "cube_schedule.h"
#include "placed_cube.h"

/* ============================================================================
 * group_graph.h -- GLOBAL integer-turn registration over the (cube, winding-
 * group) graph. Replaces the greedy spiral chain + Gauss-Seidel (which plant
 * loop-closure branch cuts that coordinate descent provably cannot lift; at
 * 4x21x21, 594 of 3714 cube-seams end a whole turn off) with:
 *
 *   nodes  = every winding group present in a cube's RAW boundary skin;
 *            unknowns are one integer k (whole turns) + one real du each.
 *   edges  = per-(gidA,gidB) constraints measured from raw skin pairs
 *            (both-side gid bucketing -- the old CubeReg_solve pooled all
 *            neighbor groups into one median, blending different wraps):
 *                k_a - k_b  = obs_e = lround(median(phi_braw-phi_araw)/2pi)
 *                du_a - du_b = du_e (k-independent; see DeltaU composition)
 *            gated by 3D distance, |dr| <= pitch/2 (cross-wrap mispairs are
 *            geometrically impossible pairs), and near-integer dphi.
 *   solve  = max-confidence spanning forest (Kruskal; the most trustworthy
 *            equations connect everything, contaminated sparse edges land
 *            OFF-tree as checks) -> BFS propagation -> per-component gauge
 *            from the radius prior (mean over member skin verts of
 *            (r-a)/b - phi_raw/2pi) -> optional EXACT collective-shift
 *            min-cut moves closing residues a tree can't (Maxflow/KZ) ->
 *            du by tree propagation + weighted-median relaxation.
 *
 * The struct is transparent so the selftest can inject synthetic graphs and
 * the driver can serialize per-node solutions.
 * ==========================================================================*/

typedef struct {
    int32_t cube;        /* index into the caller's CubeNode array */
    int32_t gid;         /* per-cube winding-group id */
    int32_t k;           /* SOLVED whole-turn correction */
    double  du;          /* SOLVED residual u shift */
    double  prior;       /* radius prior for the ABSOLUTE turn (mean of
                            (r-a)/b - phi_raw/2pi over the group's skin) */
    int32_t n_prior;     /* skin verts behind the prior */
    int32_t comp;        /* forest component id */
} GGNode;

typedef struct {
    int32_t a, b;        /* node ids; constraint k_a - k_b = obs */
    int32_t obs;         /* observed whole-turn offset */
    double  du;          /* observed du_a - du_b */
    double  conf;        /* min(n, conf_n_cap) / (1 + mad/conf_mad0) */
    double  mad;         /* median |dphi - 2pi*obs| over the pairs, rad */
    int32_t n;           /* pairs behind the edge */
    uint8_t in_tree;
} GGEdge;

typedef struct {
    GGNode  *nodes;
    size_t   n_nodes;
    GGEdge  *edges;
    size_t   n_edges;
    int32_t *cube_node0;   /* [n_cubes+1] node range per cube (gid-sorted) */
    size_t   n_cubes;
    /* solve diagnostics */
    int32_t  n_comp;
    double   energy;        /* sum conf*|residual| over edges, after solve */
    size_t   n_frustrated;  /* edges with nonzero residual after solve */
    int32_t  moves_applied;
    double   edge_redundancy;        /* non-tree edges / forest edges */
    double   anchor_weight_effective;
    size_t   pairs_used, pairs_rej_radius, pairs_rej_frac;
    size_t   edges_rej_prior, edges_rej_dr;
} GroupGraph;

typedef struct {
    double pair_gate;        /* 3D pairing gate, vox (def 3.5) */
    double radius_gate;      /* |dr| pair gate; <=0 => 0.5*pitch */
    double frac_gate;        /* max |wrap_pi(dphi)|, rad (def 0.5) */
    double prior_gate;       /* OPTIONAL cross-edge gate on |obs - dprior|,
                                turns. Default 0 = OFF: cross-cube prior
                                differences carry the scroll's eccentric
                                wobble (up to +-1 turn between cubes), so
                                gating on them rejects legitimate edges
                                wholesale (measured: 86% at 4x21x21). The
                                collapse defense is intra_prior below. */
    double edge_dr_gate;     /* bucket-median |dr| admission, vox (def 1.75).
                                THE contact-edge discriminator: same-wrap
                                pairs cross the trim gap TANGENTIALLY (same
                                surface, |dr| ~ 0.1-0.5); fused-core contact
                                pairs join two TOUCHING SHEETS separated
                                radially by the papyrus thickness (|dr| ~
                                2-3.5). Contact edges are internally
                                consistent (low mad, high conf) yet their
                                obs bridges distinct wraps as "same turn"
                                (ribbon's collapsed winding) -- admitted,
                                they stack wraps onto one u. <=0 off */
    int    intra_prior;      /* EXPERIMENTAL intra-cube prior ladders
                                (def 0: measured noisier than theorized --
                                wobble does not cancel across a cube's
                                radial extent; 13k frustrated edges) */
    double intra_frac_max;   /* |dprior - round(dprior)| admission (0.35) */
    double intra_conf;       /* confidence scale of intra edges (128) */
    int    min_edge_pairs;   /* bucket -> edge admission (def 3) */
    int    conf_n_cap;       /* pair-count saturation (def 256) */
    double conf_mad0;        /* mad scale in the confidence (def 0.02) */
    int    prior_min_verts;  /* node joins the gauge vote at >= this (def 8) */
    int    raw_component_gauge; /* shift each solved forest component closest
                                   to Ribbon's shared raw k=0 chart */
    int    raw_du_gauge;      /* preserve Ribbon's shared raw continuous chart:
                                leave every solved du at zero */
    double anchor_weight;     /* soft per-node radius anchor in collective
                                min-cut moves; 0 disables (default) */
    double anchor_redundancy_ref; /* >0 scales anchor_weight by
                                min(1,sqrt(edge_redundancy/ref)) */
    int    do_moves;         /* collective-shift min-cut moves (def 1) */
    int    max_moves;        /* cap (def 40) */
    int    du_passes;        /* weighted-median relaxation rounds (def 4) */
    int    verbose;
} GroupGraphOpts;

void GroupGraphOpts_default(GroupGraphOpts *o);

/* Build nodes + edges from the RAW skins (skins[c]/nskin[c] indexed like
 * cnodes; NULL/0 entries = cube absent). Uses cnodes[c].nbr for seam
 * adjacency; axis_point (z,y,x, axis = +z, matching CubeSched) for radii.
 * All outputs arena-allocated. Returns 0; -1 on degenerate input. */
int GroupGraph_build(Arena_T arena, const CubeNode *cnodes, size_t n_cubes,
                     SkinVert *const *skins, const size_t *nskin,
                     const float axis_point[3],
                     double spiral_a, double spiral_b, double pitch,
                     const GroupGraphOpts *opts, GroupGraph *out);

/* Solve k (forest + gauge + moves) and du (tree + relaxation) in place. */
int GroupGraph_solve(Arena_T arena, GroupGraph *g, const GroupGraphOpts *opts);

/* Extract one cube's PlacedReg from the solved graph. Table sized to the
 * cube's max gid + 1; gids without a node fall back to the cube-level
 * (n_prior-weighted majority) correction. Arena-allocated. */
int GroupGraph_cube_reg(Arena_T arena, const GroupGraph *g, size_t cube,
                        PlacedReg *out);

/* In-process unit tests (consistent cycle, branch-cut ring only a segment
 * flip repairs, two-component gauge, build filters, du recovery, empty).
 * Returns 0 if all pass. */
int GroupGraph_selftest(void);

#endif
