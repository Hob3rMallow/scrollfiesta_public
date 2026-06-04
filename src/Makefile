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
CFLAGS   = -std=c99 -Wall -Wextra -Wshadow -Wformat=2 -Werror -Wconversion
CXXFLAGS = -std=c++17 -O2 -Wno-register

# System libs + local deps (pre-installed on the target; see CLAUDE.md).
LDFLAGS  = -lm -ltiff -fopenmp -lpthread
EXT_INC  = -I../deps/include -I../deps/src/graph/include
EXT_LIB  = -L../deps/lib -ltriangle -lClipper2 -lClipper2Z -lstdc++

DEBUG_FLAGS   = -g -O0 -fsanitize=address,undefined -DVESUVIUS_DEBUG
RELEASE_FLAGS = -O2 -DNDEBUG -fopenmp
BUILD_FLAGS   = $(DEBUG_FLAGS)
TEST_INC      = -I.

# ---- Source groups (mirror cube_mesh.vcxproj) ----
COMMON_SRC = common/arena.c common/except.c common/csr.c common/kdtree.c \
             common/union_find.c common/bfs.c common/pca.c common/mls_project.c \
             common/vert_weld.c common/tiff_io.c common/halo_loader.c \
             common/mesh_trim.c common/ply_io.c common/obj_io.c common/dump_obj.c \
             common/ves_platform.c common/qem.c common/raw_snapback.c common/snap_cg.c
EXTRACT_SRC  = extract/mesh_extract.c extract/mesh_resplit.c extract/pred_reject.c
SPLIT_SRC    = split/oracle.c split/bridge_cut.c split/overlap_sep.c
REMESH_SRC   = remesh/ball_pivot.c remesh/orient_mesh.c remesh/pinhole_fill.c remesh/normal_orient.c remesh/component_cull.c
HOLEFILL_SRC = holefill/hole_fill.c
# Seamster genus reduction + topology, for the per-cube short-handle sever.
TOPO_SRC     = topology/mesh_topo.c topology/seam_cut.c
PIPELINE_SRC = pipeline/pipeline_cube.c

CORE_OBJ = $(COMMON_SRC:.c=.o) $(EXTRACT_SRC:.c=.o) $(SPLIT_SRC:.c=.o) \
           $(REMESH_SRC:.c=.o) $(HOLEFILL_SRC:.c=.o) $(TOPO_SRC:.c=.o) \
           $(PIPELINE_SRC:.c=.o)
COMMON_OBJ = $(COMMON_SRC:.c=.o)

# C++ wrappers (Clipper2 boolean ops; lifted-multicut graph solver).
CPP_OBJ = holefill/clipper2_wrap.o split/multicut_wrap.o

# grid_weld pulls these extra remesh modules (mirror grid_weld.vcxproj).
WELD_OBJ = remesh/seam_weld.o remesh/fold_cleanup.o remesh/orient_weld.o remesh/weld_cleanup.o

HDR = $(wildcard common/*.h extract/*.h split/*.h remesh/*.h holefill/*.h pipeline/*.h topology/*.h)

.PHONY: all debug release test clean
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

# Per-cube pipeline binary (== cube_mesh.exe).
cube_mesh: main.o $(CORE_OBJ) $(CPP_OBJ)
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) $(EXT_LIB)

# Cross-cube weld (== grid_weld.exe).
grid_weld: tools/grid_weld.o remesh/ball_pivot.o remesh/pinhole_fill.o remesh/orient_mesh.o \
           $(WELD_OBJ) holefill/hole_fill.o holefill/clipper2_wrap.o \
           common/arena.o common/except.o common/obj_io.o common/ves_platform.o \
           common/csr.o common/pca.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) $(EXT_LIB)

# Grid orchestrator subprocess (== grid_pipeline.exe). Links the garbage-cube
# detector (extract/pred_reject) + TIFF loader for the pre-spawn reject gate.
grid_pipeline: tools/grid_pipeline.o common/ves_platform.o common/arena.o \
               common/except.o common/tiff_io.o extract/pred_reject.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS)

# Unit tests (shared-infra suite). The full harness lives in the all_tests
# MSBuild project; this covers the common/ asserts on Linux.
test: ../scripts/test_harness.c ../scripts/test_common.c $(COMMON_OBJ)
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(TEST_INC) $^ -o test_common $(LDFLAGS)
	./test_common

clean:
	rm -f common/*.o extract/*.o split/*.o remesh/*.o holefill/*.o pipeline/*.o \
	      topology/*.o tools/*.o *.o cube_mesh grid_weld grid_pipeline test_common
