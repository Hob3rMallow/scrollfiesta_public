/*
 * scaffold_audit.c -- L1 scaffold coverage/consistency audit over a
 * scroll_whole placed dir. Entry point only: parse args, load the piece set,
 * run Scaffold_kstats (E0: derived-k coverage + per-face winding coherence),
 * print the gate verdicts. Occupancy classification (E2) is added to
 * scaffold.c and surfaced here once E0 passes.
 *
 * usage: scaffold_audit <placed_dir> [--sense S] [--hist out.json]
 *        scaffold_audit --selftest
 */
#include "../common/ves_platform.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/arena.h"
#include "../unroll/piece_set.h"
#include "../unroll/scaffold.h"

int main(int argc, char **argv)
{
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) {
        int f = Scaffold_selftest();
        f += PieceSet_selftest();
        fprintf(stderr, "=== scaffold_audit selftest %s (%d failure%s) ===\n",
                f ? "FAILED" : "PASSED", f, f == 1 ? "" : "s");
        return f ? 1 : 0;
    }
    if (argc < 2) {
        fprintf(stderr,
            "usage: scaffold_audit <placed_dir> [--sense S] [--hist out.json]\n"
            "       scaffold_audit --selftest\n");
        return 1;
    }

    const char *placed = argv[1];
    int   sense_override = 0, sense_val = 0;
    const char *hist = NULL;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--sense") && i + 1 < argc) {
            sense_override = 1; sense_val = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--hist") && i + 1 < argc) {
            hist = argv[++i];
        } else {
            fprintf(stderr, "scaffold_audit: unknown arg %s\n", argv[i]);
            return 1;
        }
    }

    Arena_T ar = Arena_new();

    ScaffoldCalib c;
    if (Scaffold_read_calib(placed, &c) != 0) {
        fprintf(stderr, "scaffold_audit: no placed_index.json under %s "
                "(using 0139 defaults)\n", placed);
    }
    if (sense_override) c.sense = sense_val;
    fprintf(stderr, "calib: a=%.6f b=%.6f pitch=%.3f sense=%d "
            "umbilicus=(%.1f,%.1f,%.1f)\n",
            c.spiral_a, c.spiral_b, c.pitch, c.sense,
            (double)c.axis_point[0], (double)c.axis_point[1],
            (double)c.axis_point[2]);

    PieceSet ps;
    if (PieceSet_build(ar, placed, &ps) != 0) {
        fprintf(stderr, "scaffold_audit: PieceSet_build failed for %s\n", placed);
        Arena_dispose(&ar);
        return 1;
    }
    fprintf(stderr, "pieces: %zu cubes, %zu verts, %zu kept faces "
            "(u=[%.0f,%.0f] v=[%.0f,%.0f])\n",
            ps.n_cubes, ps.nv, ps.nf, ps.u_min, ps.u_max, ps.v_min, ps.v_max);

    ScaffoldKStats st;
    if (Scaffold_kstats(ar, &ps, &c, &st, hist) != 0) {
        fprintf(stderr, "scaffold_audit: Scaffold_kstats failed\n");
        Arena_dispose(&ar);
        return 1;
    }

    printf("=== scaffold E0 (derived-k) ===\n");
    printf("  nv_used = %zu   nf_used = %zu\n", st.nv_used, st.nf_used);
    printf("  k range [%d, %d]  extent = %d turns\n",
           st.k_min, st.k_max, st.k_extent);
    printf("  per-face unwound-phi spread: p50=%.4f  p99=%.4f rad\n",
           st.p50_face_spread, st.p99_face_spread);
    printf("  faces spread >= pi: %zu = %.4f%% of kept faces\n",
           st.face_spread_ge_pi, 100.0 * st.face_spread_frac);
    printf("  groups (cube,gid): %zu   spread >= pi: %zu = %.2f%%\n",
           st.n_groups, st.group_spread_ge_pi,
           st.n_groups ? 100.0 * (double)st.group_spread_ge_pi / (double)st.n_groups : 0.0);
    printf("  GATE face-spread < 1%%: %s\n",
           st.face_spread_frac < 0.01 ? "PASS" : "FAIL");
    if (hist) printf("  wrote k x cube-z histogram -> %s\n", hist);

    Arena_dispose(&ar);
    return 0;
}
