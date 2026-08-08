# Vesuvius-C build (Linux / GCC).
#
# This Makefile mirrors the source lists in the MSBuild projects at the repo
# root (cube_mesh.vcxproj / grid_weld.vcxproj / grid_pipeline.vcxproj), which
# are the canonical, CI-verified build. The old "7-step" pipeline was retired
# in 2026-05; the live stages are:
#
#     extract/  ->  split/  ->  remesh/  +  holefill/   (orchestrated by pipeline/)
#
# Run from src/:
#     make            per-cube pipeline (cube_mesh), debug+ASan
#     make release    optimized
#     make grid_weld grid_pipeline
#     make test       unit tests (test_common)
#
# NOTE: the source lists below were reconciled against the .vcxproj on Windows
# but this Makefile has NOT been compile-tested on Linux since the reorg — give
# it a `make` on the deploy box before trusting it.

CC       = gcc
CXX      = g++
# -fno-strict-aliasing: gcc 16.1.0 miscompiles correct interleaved-store mesh
# loops at -O3 under strict aliasing (see scripts/gcc16_aliasing_repro.c).
CFLAGS   = -std=c99 -Wall -Wextra -Wshadow -Wformat=2 -Werror=implicit-function-declaration -D_GNU_SOURCE -fno-strict-aliasing
CXXFLAGS = -std=c++17 -O2 -Wno-register

# System libs + local deps (pre-installed on the target; see CLAUDE.md).
LDFLAGS  = -lm -ltiff -fopenmp -lpthread
EXT_INC  = -I../deps/include -I../deps/src/graph/include -I../docs/research/gco-v3.0
EXT_LIB  = -L../deps/lib -ltriangle -lClipper2 -lstdc++

DEBUG_FLAGS   = -g -O0 -fsanitize=address,undefined -DVESUVIUS_DEBUG
RELEASE_FLAGS = -O2 -DNDEBUG -fopenmp
BUILD_FLAGS   = $(DEBUG_FLAGS)
TEST_INC      = -I.

# ---- Source groups (mirror cube_mesh.vcxproj) ----
COMMON_SRC = common/arena.c common/except.c common/csr.c common/kdtree.c \
             common/union_find.c common/cc_color.c common/bfs.c common/pca.c common/mls_project.c \
             common/vert_weld.c common/tiff_io.c common/halo_loader.c \
             common/raw_sample.c \
             common/mesh_trim.c common/obj_io.c common/dump_obj.c \
             common/ves_platform.c common/run_ctx.c common/qem.c common/raw_snapback.c common/snap_cg.c \
             common/eig3.c common/dev_metric.c common/mesh_manifold.c
EXTRACT_SRC  = extract/mesh_extract.c extract/mesh_resplit.c extract/pred_reject.c
SPLIT_SRC    = split/oracle.c split/bridge_cut.c split/overlap_sep.c split/overlap_merge.c split/dev_cut.c split/depth_peel.c
REMESH_SRC   = remesh/ball_pivot.c remesh/orient_mesh.c remesh/pinhole_fill.c remesh/normal_orient.c remesh/component_cull.c remesh/manifold_guard.c remesh/dev_gate.c
HOLEFILL_SRC = holefill/hole_fill.c
# Seamster genus reduction + topology, for the per-cube short-handle sever.
FLATTEN_SRC  = flatten/mesh_topo.c flatten/seam_cut.c flatten/developability.c \
               flatten/surface_snap.c flatten/snap_darkcut.c flatten/snap_quilt.c \
               flatten/snap_solve.c flatten/recto_refine.c
PIPELINE_SRC = pipeline/pipeline_cube.c

CORE_OBJ = $(COMMON_SRC:.c=.o) $(EXTRACT_SRC:.c=.o) $(SPLIT_SRC:.c=.o) \
           $(REMESH_SRC:.c=.o) $(HOLEFILL_SRC:.c=.o) $(FLATTEN_SRC:.c=.o) \
           $(PIPELINE_SRC:.c=.o)
COMMON_OBJ = $(COMMON_SRC:.c=.o)

# C++ wrappers (Clipper2 boolean ops; lifted-multicut graph solver).
CPP_OBJ = holefill/clipper2_wrap.o split/multicut_wrap.o common/gco_wrap.o \
          ../docs/research/gco-v3.0/GCoptimization.o \
          ../docs/research/gco-v3.0/LinkedBlockList.o

# grid_weld pulls these extra remesh modules (mirror grid_weld.vcxproj).
WELD_OBJ = remesh/seam_weld.o remesh/seam_planes.o remesh/seam_refine.o remesh/seam_band_cvt.o \
           remesh/cvt_remesh.o remesh/rvd_clip.o common/closest_tri.o \
           remesh/sheet_reweld.o remesh/fold_cleanup.o remesh/orient_weld.o remesh/weld_cleanup.o
# seam_band_cvt pulls the CVT/RVD exact predicates (C++17), as cube_mesh does.
WELD_CPP_OBJ = common/attene_predicates_wrap.o \
               ../deps/src/attene-predicates/uinteger.o \
               ../deps/src/attene-predicates/ieee_binary.o

HDR = $(wildcard common/*.h extract/*.h split/*.h remesh/*.h holefill/*.h pipeline/*.h flatten/*.h)

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
# boundary-preserving; the quality-keeping replacement for qslim in the pyramid.
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

# Cross-cube weld (== grid_weld.exe).
grid_weld: tools/grid_weld.o remesh/ball_pivot.o remesh/pinhole_fill.o remesh/orient_mesh.o \
           $(WELD_OBJ) $(WELD_CPP_OBJ) holefill/hole_fill.o holefill/clipper2_wrap.o \
           common/arena.o common/except.o common/obj_io.o common/ves_platform.o \
           common/csr.o common/pca.o common/kdtree.o common/union_find.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) $(EXT_LIB)

common/attene_predicates_wrap.o: common/attene_predicates_wrap.cpp
	$(CXX) $(CXXFLAGS) -I../deps/src/attene-predicates -c $< -o $@
../deps/src/attene-predicates/%.o: ../deps/src/attene-predicates/%.cpp
	$(CXX) $(CXXFLAGS) -I../deps/src/attene-predicates -c $< -o $@

# Grid orchestrator subprocess (== grid_pipeline.exe). Links the garbage-cube
# detector (extract/pred_reject) + TIFF loader for the pre-spawn reject gate.
grid_pipeline: tools/grid_pipeline.o common/ves_platform.o common/arena.o \
               common/except.o common/tiff_io.o extract/pred_reject.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS)


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
	      flatten/*.o tools/*.o *.o cube_mesh grid_weld grid_pipeline \
	      obj_reorient test_common dev_cut_test depth_peel_test sheet_reweld_test
