# Vesuvius-C build (Linux / GCC).
#
# This Makefile mirrors the source lists in the MSBuild projects at the repo
# root (cube_mesh.vcxproj / grid_weld.vcxproj / grid_pipeline.vcxproj, plus the
# obj_components / pinhole_verdict / seam_audit tool projects), which are the
# canonical, CI-verified build. The old "7-step" pipeline was retired in
# 2026-05; the live stages are:
#
#     extract/  ->  split/  ->  remesh/  +  holefill/   (orchestrated by pipeline/)
#
# Run from src/:
#     make            per-cube pipeline (cube_mesh), debug+ASan
#     make release    all six binaries, optimized (cube_mesh, grid_pipeline,
#                     grid_weld, obj_components, pinhole_verdict, seam_audit)
#     make grid_weld grid_pipeline
#     make test       unit tests (test_common)
#
# Linux prerequisites: run ../build-deps.sh once to stage the vendored static
# libs into ../deps/lib (libtriangle.a, libClipper2.a, libClipper2Z.a), and
# install the system libtiff (libtiff-dev) that -ltiff links against — its
# headers live in the multiarch include dir that GCC searches by default.
# Ported to / verified with GCC 13; POSIX/XSI symbols are exposed via -D_GNU_SOURCE.

CC       = gcc
CXX      = g++
CFLAGS   = -std=c99 -Wall -Wextra -Wshadow -Wformat=2 -Werror=implicit-function-declaration -Wconversion -D_GNU_SOURCE
CXXFLAGS = -std=c++17 -O2 -Wno-register -D_GNU_SOURCE

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
             common/mesh_manifold.c common/mesh_trim.c common/ply_io.c common/obj_io.c common/dump_obj.c \
             common/ves_platform.c common/qem.c common/raw_snapback.c common/snap_cg.c \
             common/eig3.c common/dev_metric.c
EXTRACT_SRC  = extract/mesh_extract.c extract/mesh_resplit.c extract/pred_reject.c
SPLIT_SRC    = split/oracle.c split/bridge_cut.c split/overlap_sep.c split/depth_peel.c split/dev_cut.c
REMESH_SRC   = remesh/ball_pivot.c remesh/orient_mesh.c remesh/pinhole_fill.c remesh/normal_orient.c remesh/component_cull.c remesh/manifold_guard.c remesh/dev_gate.c
HOLEFILL_SRC = holefill/hole_fill.c
# Seamster genus reduction + topology, for the per-cube short-handle sever.
TOPO_SRC     = topology/mesh_topo.c topology/seam_cut.c topology/developability.c
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
release: cube_mesh grid_pipeline grid_weld obj_components pinhole_verdict seam_audit

# Generic C rule (covers common/ extract/ split/ remesh/ holefill/ pipeline/ tools/).
%.o: %.c $(HDR)
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $(EXT_INC) -c $< -o $@

main.o: main.c $(HDR)
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $(EXT_INC) -c $< -o $@

holefill/clipper2_wrap.o: holefill/clipper2_wrap.cpp
	$(CXX) $(CXXFLAGS) -I../deps/include -c $< -o $@
# The vendored andres/graph headers use std::numeric_limits without including
# <limits>; GCC 13's libstdc++ no longer pulls it in transitively, so force it.
split/multicut_wrap.o: split/multicut_wrap.cpp
	$(CXX) $(CXXFLAGS) $(EXT_INC) -DNOMINMAX -include limits -c $< -o $@

# Per-cube pipeline binary (== cube_mesh.exe).
cube_mesh: main.o $(CORE_OBJ) $(CPP_OBJ)
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) $(EXT_LIB)

# Cross-cube weld (== grid_weld.exe).
grid_weld: tools/grid_weld.o remesh/ball_pivot.o remesh/pinhole_fill.o remesh/orient_mesh.o \
           remesh/manifold_guard.o $(WELD_OBJ) holefill/hole_fill.o holefill/clipper2_wrap.o \
           common/arena.o common/except.o common/obj_io.o common/ves_platform.o \
           common/mesh_manifold.o common/csr.o common/pca.o \
           topology/developability.o common/eig3.o common/vert_weld.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS) $(EXT_LIB)

# Grid orchestrator subprocess (== grid_pipeline.exe). Links the garbage-cube
# detector (extract/pred_reject) + TIFF loader for the pre-spawn reject gate.
grid_pipeline: tools/grid_pipeline.o common/ves_platform.o common/arena.o \
               common/except.o common/tiff_io.o extract/pred_reject.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ $(LDFLAGS)

# Diagnostic tools (== obj_components.exe / pinhole_verdict.exe / seam_audit.exe).
# Each is a single self-contained translation unit with its own main(), mirroring
# its .vcxproj (one ClCompile, no shared objects); only libm is needed.
obj_components: tools/obj_components.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ -lm
pinhole_verdict: tools/pinhole_verdict.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ -lm
seam_audit: tools/seam_audit.o
	$(CC) $(CFLAGS) $(BUILD_FLAGS) $^ -o $@ -lm

# Unit tests (shared-infra suite). On Linux this builds test_common.c standalone
# (it provides its own main() when TEST_HARNESS is undefined). test_harness.c is
# the Windows all_tests dispatcher -- it needs every *_test_main symbol, so it is
# intentionally NOT compiled here.
test: ../scripts/test_common.c $(COMMON_OBJ)
	$(CC) $(CFLAGS) $(DEBUG_FLAGS) $(TEST_INC) $^ -o test_common $(LDFLAGS)
	./test_common

clean:
	rm -f common/*.o extract/*.o split/*.o remesh/*.o holefill/*.o pipeline/*.o \
	      topology/*.o tools/*.o *.o cube_mesh grid_weld grid_pipeline \
	      obj_components pinhole_verdict seam_audit test_common
