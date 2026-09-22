# Vesuvius-C build (Linux / GCC).
#
# This Makefile mirrors the source lists in the MSBuild projects at the repo
# root (see scrollfiesta.sln), which
# are the canonical, CI-verified build. The old "7-step" pipeline was retired
# in 2026-05; the live stages are:
#
#     extract/  ->  split/  ->  remesh/  +  holefill/   (orchestrated by pipeline/)
#
# Run from src/:
#     make            per-cube pipeline (cube_mesh), debug+ASan
#     make release    optimized
#     make grid_weld grid_pipeline scroll_whole welded_winding_lift obj_bake_raw
#     make test       unit tests (test_common)
#
# NOTE: the source lists below were reconciled against the .vcxproj on Windows
# but this Makefile has NOT been compile-tested on Linux since the reorg — give
# it a `make` on the deploy box before trusting it.

CC       = gcc
CXX      = g++
# -std=gnu99, not -std=c99: ves_platform.h calls posix_fadvise() and main.c
# calls setenv(), which strict ISO C hides. Under GCC 15 an implicit declaration
# is an error rather than a warning, so the very first object file fails.
#
# -Werror is opt-in rather than always on. GCC 15 promotes several warning
# classes on code that GCC 13 accepted -- sign-conversion, clobbered,
# misleading-indentation, format-truncation -- and with -Werror the build stops
# at the third source file, long before anything links. The warnings still
# print; only the promotion to error is gated, so CI keeps its strictness by
# building with SF_WERROR=1 and a newer compiler can still produce binaries.
#
#     make release                # warnings print, build proceeds
#     SF_WERROR=1 make release    # the strict build CI runs
WERROR   = $(if $(filter 1,$(SF_WERROR)),-Werror,)
CFLAGS   = -std=gnu99 -Wall -Wextra -Wshadow -Wformat=2 $(WERROR) -Wconversion
# -include limits: the vendored andres/graph headers use std::numeric_limits
# without including <limits>. libstdc++ used to pull it in transitively; since
# GCC 15 it does not, so greedy-additive.hxx and kernighan-lin.hxx fail to
# compile. Forcing the include here leaves the vendored tree untouched.
CXXFLAGS = -std=c++17 -O2 -Wno-register -include limits

# System libs + local deps (pre-installed on the target; see CLAUDE.md).
# The vendored libtiff that build-deps.sh stages into deps/prefix is preferred;
# where it has not been staged this falls back to -ltiff, the system library,
# exactly as before. zlib is staged into the same prefix and a static libtiff
# needs it.
TIFF_A     := $(wildcard ../deps/prefix/lib/libtiff.a)
TIFF_EXTRA := $(wildcard ../deps/prefix/lib/libz.a)
TIFF_LIB   := $(if $(TIFF_A),$(TIFF_A) $(TIFF_EXTRA),-ltiff)
LDFLAGS  = -lm $(TIFF_LIB) -fopenmp -lpthread
EXT_INC  = -I../deps/include -I../deps/src/graph/include -I../docs/research/gco-v3.0
EXT_LIB  = -L../deps/lib -ltriangle -lClipper2 -lClipper2Z -lstdc++

DEBUG_FLAGS   = -g -O0 -fsanitize=address,undefined -DVESUVIUS_DEBUG
RELEASE_FLAGS = -O2 -DNDEBUG -fopenmp
BUILD_FLAGS   = $(DEBUG_FLAGS)
TEST_INC      = -I.

# ---- Source groups (mirror cube_mesh.vcxproj) ----
# common/pitch_table.c e common/closest_tri.c sono su disco e referenziati (PitchTable_*,
# ClosestTri_*) ma non erano in nessuna lista: il link di cube_mesh li chiedeva come
# simboli non definiti.
COMMON_SRC = common/pitch_table.c common/closest_tri.c common/arena.c common/except.c common/csr.c common/kdtree.c \
             common/union_find.c common/cc_color.c common/bfs.c common/pca.c common/mls_project.c \
             common/vert_weld.c common/mesh_normals.c common/tiff_io.c common/halo_loader.c \
             common/raw_sample.c \
             common/mesh_trim.c common/ply_io.c common/obj_io.c common/mesh_bin.c common/dump_obj.c \
             common/ves_platform.c common/raw_snapback.c common/snap_cg.c \
             common/eig3.c common/dev_metric.c common/mesh_manifold.c
EXTRACT_SRC  = whole/axis_warp.c extract/mesh_extract.c extract/mesh_resplit.c extract/pred_reject.c
SPLIT_SRC    = split/overlap_quality.c split/oracle.c split/bridge_cut.c split/overlap_sep.c split/overlap_merge.c split/dev_cut.c split/depth_peel.c
REMESH_SRC   = remesh/patch_repair.c remesh/cvt_remesh.c remesh/rvd_clip.c remesh/ball_pivot.c remesh/orient_mesh.c remesh/pinhole_fill.c remesh/normal_orient.c remesh/component_cull.c remesh/manifold_guard.c remesh/dev_gate.c
HOLEFILL_SRC = holefill/hole_fill.c
# Seamster genus reduction + topology, for the per-cube short-handle sever.
FLATTEN_SRC  = flatten/mesh_topo.c flatten/topology_invariants.c flatten/seam_cut.c flatten/developability.c \
               flatten/surface_snap.c flatten/snap_darkcut.c flatten/snap_quilt.c \
               flatten/snap_solve.c flatten/recto_refine.c
PIPELINE_SRC = pipeline/pipeline_cube.c

CORE_OBJ = $(COMMON_SRC:.c=.o) $(EXTRACT_SRC:.c=.o) $(SPLIT_SRC:.c=.o) \
           $(REMESH_SRC:.c=.o) $(HOLEFILL_SRC:.c=.o) $(FLATTEN_SRC:.c=.o) \
           $(PIPELINE_SRC:.c=.o)
COMMON_OBJ = $(COMMON_SRC:.c=.o)

# C++ wrappers (Clipper2 boolean ops; lifted-multicut graph solver).
# Le tre righe di attene servono a cube_mesh perche' remesh/cvt_remesh.c chiama
# CVT_remesh -> rvd_clip -> AttenePred_*. Le regole per costruirle esistevano gia' piu'
# in basso, ma solo WELD_CPP_OBJ le usava.
CPP_OBJ = holefill/clipper2_wrap.o split/multicut_wrap.o common/gco_wrap.o \
          common/attene_predicates_wrap.o \
          ../deps/src/attene-predicates/uinteger.o \
          ../deps/src/attene-predicates/ieee_binary.o \
          ../docs/research/gco-v3.0/GCoptimization.o \
          ../docs/research/gco-v3.0/LinkedBlockList.o

# grid_weld pulls these extra remesh modules (mirror grid_weld.vcxproj).
WELD_OBJ = remesh/seam_weld.o remesh/seam_planes.o remesh/seam_refine.o remesh/seam_band_cvt.o \
           remesh/cvt_remesh.o remesh/rvd_clip.o common/closest_tri.o \
           remesh/sheet_reweld.o remesh/fold_cleanup.o remesh/orient_weld.o remesh/weld_cleanup.o \
           remesh/intersection_cleanup.o remesh/boundary_arc_weld.o \
           remesh/chart_lineage.o remesh/chart_bridge_forest.o remesh/chart_zipper.o \
           remesh/stream_weld_format.o
# seam_band_cvt pulls the CVT/RVD exact predicates (C++17), as cube_mesh does.
WELD_CPP_OBJ = common/attene_predicates_wrap.o \
               ../deps/src/attene-predicates/uinteger.o \
               ../deps/src/attene-predicates/ieee_binary.o

HDR = $(wildcard common/*.h extract/*.h split/*.h remesh/*.h holefill/*.h \
                 pipeline/*.h flatten/*.h whole/*.h unroll/*.h)

.PHONY: all debug release test dev_cut_test depth_peel_test sheet_reweld_test obj_reorient_test mesh_trim_cut_test weld_recoarsen_test clean
all: debug

debug:   BUILD_FLAGS = $(DEBUG_FLAGS)
debug:   cube_mesh
release: BUILD_FLAGS = $(RELEASE_FLAGS)
release: cube_mesh

# Generic C rule (covers common/ extract/ split/ remesh/ holefill/ pipeline/ tools/).
%.o: %.c $(HDR)
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $(EXT_INC) -c $< -o $@

main.o: main.c $(HDR)
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $(EXT_INC) -c $< -o $@

holefill/clipper2_wrap.o: holefill/clipper2_wrap.cpp
	$(CXX) $(CXXFLAGS) -I../deps/include -c $< -o $@
split/multicut_wrap.o: split/multicut_wrap.cpp
	$(CXX) $(CXXFLAGS) $(EXT_INC) -DNOMINMAX -c $< -o $@
common/gco_wrap.o: common/gco_wrap.cpp
	$(CXX) $(CXXFLAGS) $(EXT_INC) -c $< -o $@
../docs/research/gco-v3.0/%.o: ../docs/research/gco-v3.0/%.cpp
	$(CXX) $(CXXFLAGS) $(EXT_INC) -c $< -o $@

# Per-cube pipeline binary (== cube_mesh.exe).
cube_mesh: main.o $(CORE_OBJ) $(CPP_OBJ)
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) $(EXT_LIB)

# CVT-based LOD simplifier (== cvt_simplify.exe): per-component CVT decimation,
# boundary-preserving; the sole supported LOD simplifier in the pyramid.
cvt_simplify: tools/cvt_simplify.o remesh/cvt_remesh.o remesh/rvd_clip.o \
              common/closest_tri.o common/kdtree.o \
              remesh/manifold_guard.o remesh/pinhole_fill.o remesh/orient_mesh.o \
              common/mesh_manifold.o common/obj_io.o common/arena.o common/except.o \
              common/attene_predicates_wrap.o \
              ../deps/src/attene-predicates/uinteger.o ../deps/src/attene-predicates/ieee_binary.o
	$(CXX) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) -lstdc++

# Winding-coordinate topology auditor (== wind_audit.exe): flags inter-wrap
# fusions (|dw|~1 radial short-circuits) + intra-sheet splits (same-turn pieces
# left disconnected). Standalone (own OBJ parser, libm only), like obj_cc_color;
# not part of the pipeline link. Relaxed flags: not -Werror -Wconversion tested
# on Linux yet -- MSVC wind_audit.vcxproj is the CI-verified build.
wind_audit: tools/wind_audit.c
	$(CC) -std=c99 -O2 -Wall $< -o $@ -lm

# Component topology certificate + minimum-cardinality H1 generator basis
# (== topology_audit.exe).  Betti ranks, orientability, edge/vertex manifold
# checks, disk verdicts, and tree-cotree loop representatives stay in native C.
topology_audit: tools/topology_audit.o flatten/topology_invariants.o \
                flatten/disk_topology_repair.o flatten/seam_cut.o \
                flatten/mesh_topo.o common/union_find.o common/ves_platform.o \
                common/mesh_manifold.o common/obj_io.o common/arena.o \
                common/except.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ -lm

# Native tifxyz -> Villa winding evidence/PCL miner.  Confidence and acceptance
# policy intentionally stay in C; higher-level packages only launch this tool.
winding_evidence: tools/winding_evidence.o common/tiff_io.o common/kdtree.o \
                  common/arena.o common/except.o common/ves_platform.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS)

# These standalone boundary tools are warning-clean under the canonical MSVC
# Level4 builds. Keep their first Linux build diagnostic rather than failing on
# the legacy tree-wide -Wconversion policy before that platform is calibrated.
tools/winding_evidence.o: tools/winding_evidence.c $(HDR)
	$(CC) -std=c99 -Wall -Wextra $(BUILD_FLAGS) $(EXT_INC) -c $< -o $@

# Verified same-wrap tifxyz strip/patch assembler.  Native code owns lattice,
# overlap, conflict-hiding, and contested-area acceptance policy.
tifxyz_weld: tools/tifxyz_weld.o common/tiff_io.o common/arena.o \
             common/except.o common/ves_platform.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS)

tools/tifxyz_weld.o: tools/tifxyz_weld.c $(HDR)
	$(CC) -std=c99 -Wall -Wextra $(BUILD_FLAGS) $(EXT_INC) -c $< -o $@

# Cross-cube weld (== grid_weld.exe).
grid_weld: tools/grid_weld.o remesh/ball_pivot.o remesh/pinhole_fill.o remesh/orient_mesh.o \
           $(WELD_OBJ) $(WELD_CPP_OBJ) holefill/hole_fill.o holefill/clipper2_wrap.o \
           flatten/disk_topology_repair.o flatten/topology_invariants.o \
           flatten/seam_cut.o flatten/mesh_topo.o \
           common/arena.o common/except.o common/obj_io.o common/mesh_bin.o common/ves_platform.o \
           common/csr.o common/pca.o common/kdtree.o common/union_find.o \
           common/cc_color.o common/mesh_manifold.o common/pitch_table.o common/vert_weld.o \
           common/mesh_normals.o common/bfs.o \
           remesh/manifold_guard.o remesh/seam_hole_fill.o remesh/component_cull.o \
           flatten/developability.o common/eig3.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) $(EXT_LIB)

common/attene_predicates_wrap.o: common/attene_predicates_wrap.cpp
	$(CXX) $(CXXFLAGS) -I../deps/src/attene-predicates -c $< -o $@
../deps/src/attene-predicates/%.o: ../deps/src/attene-predicates/%.cpp
	$(CXX) $(CXXFLAGS) -I../deps/src/attene-predicates -c $< -o $@

# Grid orchestrator subprocess (== grid_pipeline.exe). Links the garbage-cube
# detector (extract/pred_reject) + TIFF loader for the pre-spawn reject gate.
grid_pipeline: tools/grid_pipeline.o common/ves_platform.o common/arena.o \
               common/except.o common/tiff_io.o extract/pred_reject.o common/mesh_bin.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS)

# Register and place per-cube charts in one scroll coordinate frame.
scroll_whole: tools/scroll_whole.o whole/cube_schedule.o whole/axis_warp.o \
              whole/placed_cube.o whole/cube_register.o whole/group_graph.o \
              whole/weld_track.o whole/weld_build.o whole/arc_reg.o \
              common/maxflow.o flatten/ribbon.o flatten/ribbon_grid_io.o \
              remesh/intersection_cleanup.o flatten/seam_cut.o \
              flatten/mesh_topo.o flatten/topology_invariants.o \
              common/arena.o common/except.o common/obj_io.o common/mesh_bin.o \
              common/mesh_manifold.o common/ves_platform.o common/pca.o \
              common/csr.o common/union_find.o common/kdtree.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS)

# Preserve scroll_whole's winding-safe chart UV while lifting it through an
# accepted grid_weld topology (== welded_winding_lift.exe).
welded_winding_lift: tools/welded_winding_lift.o unroll/piece_set.o \
                     unroll/scaffold.o flatten/ribbon.o \
                     flatten/ribbon_grid_io.o remesh/intersection_cleanup.o \
                     common/arena.o common/except.o common/obj_io.o common/mesh_bin.o \
                     common/ves_platform.o common/pca.o common/csr.o \
                     common/union_find.o common/kdtree.o common/mesh_normals.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS)

# Solid quad-ribbon deliverable assembler (== solid_quad_ribbon.exe): one full
# num_v x num_u quad strip per scroll component, positions sampled from the
# fitted field and Dirichlet-Laplacian filled (SnapCG) elsewhere.  Mirrors
# solid_quad_ribbon.vcxproj.  -lz is for libtiff, not ves_png (stb is self-
# contained).  json_read.o / quad_strip.o build via the generic %.o rule.
solid_quad_ribbon: tools/solid_quad_ribbon.o flatten/quad_strip.o \
                   flatten/quad_strip_snap.o common/json_read.o \
                   common/csr.o common/snap_cg.o common/eig3.o common/raw_sample.o \
                   common/mesh_normals.o common/tiff_io.o common/obj_io.o \
                   common/mesh_bin.o common/ves_png.o common/arena.o common/except.o \
                   common/ves_platform.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) -lz

# Local terminal-ARAP correction for a detector-confirmed marbled UV patch.
marble_repair: tools/marble_repair.o flatten/quad_strip.o common/csr.o \
               common/snap_cg.o common/eig3.o common/mesh_bin.o \
               common/arena.o common/except.o common/ves_platform.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS)

# Local UV-only marble correction. The detector changes trust continuously;
# the actual solve is StrokeStrip's relaxed/final variational parameterization
# on a coarse ribbon grid. XYZ, topology, and axial V never move. The canonical
# Windows build owns the TAUCS/CLAPACK source list in marble_uv_repair.vcxproj.
marble_uv_repair: tools/marble_uv_repair.o flatten/marble_strip_uv.o \
                  whole/atlas_strip.o whole/monotone_qp.o \
                  flatten/sparse_solve.o common/mesh_bin.o common/tiff_io.o \
                  common/arena.o common/except.o common/ves_platform.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) -ltiff -lz

# Its unit test is the embedded --selftest (json_read + quad_strip + mesh_bin).
solid_quad_ribbon_test: solid_quad_ribbon
	./solid_quad_ribbon --selftest

# Per-chart integer winding-depth search over the placed frame, with an exact
# triangle-SAT overlap term (== atlas_overlap_fix.exe).  Mirrors
# atlas_overlap_fix.vcxproj; CubeReg_deltaU is static inline in the header, so
# cube_register.o is deliberately absent.
atlas_overlap_fix: tools/atlas_overlap_fix_tool.o \
                   tools/atlas_overlap_fix_sparse.o \
                   whole/atlas_overlap_fix.o \
                   whole/atlas_collision_register.o \
                   whole/atlas_overlap_audit.o whole/atlas_register.o \
                   whole/atlas_solution.o whole/monotone_qp.o \
                   unroll/piece_set.o unroll/scaffold.o common/arena.o \
                   common/except.o common/union_find.o common/obj_io.o \
                   common/mesh_bin.o common/raw_sample.o \
                   common/mesh_normals.o common/tiff_io.o \
                   common/ves_platform.o common/ves_png.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) -lz

# Sample the source RAW volumes onto the authoritative parameterized VMESH.
obj_bake_raw: tools/obj_bake_raw.o flatten/rawtex_bake.o flatten/marble_diag.o common/ves_png.o \
              common/tif_strip.o common/raw_sample.o common/arena.o \
              common/mesh_bin.o common/except.o common/ves_platform.o \
              common/tiff_io.o common/mesh_normals.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) -lz

# Orientation-repair tool (== obj_reorient.exe): fixes backward-wound components
# + non-orientable knots on a welded mesh. Mirrors obj_reorient.vcxproj.
obj_reorient: tools/obj_reorient.o common/obj_io.o common/mesh_manifold.o \
              remesh/orient_mesh.o remesh/orient_weld.o remesh/pinhole_fill.o \
              common/arena.o common/except.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS)

# Its unit test is the embedded --selftest (like obj_orient_audit/manifold_check).
obj_reorient_test: obj_reorient
	./obj_reorient --selftest

# Unit tests (shared-infra suite). The full harness lives in the all_tests
# MSBuild project; this covers the common/ asserts on Linux.
test: ../scripts/test_harness.c ../scripts/test_common.c $(COMMON_OBJ)
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(TEST_INC) $^ -o test_common $(LDFLAGS)
	./test_common

# Developability-cut unit test (split/dev_cut + its Crane-energy / eig3 deps).
# Linux parity for dev_cut_test.vcxproj.
dev_cut_test: ../scripts/split/test_dev_cut.c split/dev_cut.o flatten/developability.o \
              common/eig3.o common/arena.o common/except.o common/kdtree.o \
              common/union_find.o common/pca.o
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(TEST_INC) $^ -o dev_cut_test $(LDFLAGS)
	./dev_cut_test

# Stacked-wrap depth-peel unit test (split/depth_peel + arena/kdtree/uf/pca deps).
# Linux parity for depth_peel_test.vcxproj.
depth_peel_test: ../scripts/split/test_depth_peel.c split/depth_peel.o \
              common/arena.o common/except.o common/kdtree.o \
              common/union_find.o common/pca.o
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(TEST_INC) $^ -o depth_peel_test $(LDFLAGS)
	./depth_peel_test

# Phase-2 sheet-correspondence re-weld unit test (remesh/sheet_reweld + the
# seam-bridge deps). Linux parity for sheet_reweld_test.vcxproj.
sheet_reweld_test: ../scripts/remesh/test_sheet_reweld.c remesh/sheet_reweld.o \
              remesh/seam_weld.o remesh/seam_planes.o remesh/ball_pivot.o common/kdtree.o \
              common/union_find.o common/arena.o common/except.o
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(TEST_INC) $^ -o sheet_reweld_test $(LDFLAGS)
	./sheet_reweld_test

# Cut-at-plane trim unit test. Linux parity for mesh_trim_cut_test.vcxproj.
mesh_trim_cut_test: ../scripts/extract/mesh_trim_cut_test.c common/mesh_trim.o \
              common/arena.o common/except.o
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(TEST_INC) $^ -o mesh_trim_cut_test $(LDFLAGS)
	./mesh_trim_cut_test

# Seam-band recoarsen unit test. Linux parity for weld_recoarsen_test.vcxproj.
weld_recoarsen_test: ../scripts/remesh/weld_recoarsen_test.c remesh/weld_cleanup.o \
              remesh/seam_planes.o common/arena.o common/except.o
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(TEST_INC) $^ -o weld_recoarsen_test $(LDFLAGS)
	./weld_recoarsen_test

clean:
	rm -f common/*.o extract/*.o split/*.o remesh/*.o holefill/*.o pipeline/*.o \
	      flatten/*.o whole/*.o unroll/*.o tools/*.o *.o cube_mesh grid_weld \
	      grid_pipeline scroll_whole welded_winding_lift obj_bake_raw \
	      obj_reorient topology_audit test_common dev_cut_test depth_peel_test sheet_reweld_test
