# ScrollFiesta — CMake build.
#
# Builds the public C library (target `scrollfiesta`, alias
# `scrollfiesta::scrollfiesta`, header include/scrollfiesta.h) plus, when
# top-level, the supported pipeline CLI tools and small diagnostics. The
# legacy MSBuild and Makefile builds remain available for narrower workflows;
# CMake is the portable, CI-verified build of record.
#
# Standalone:
#     cmake -S . -B build && cmake --build build
# Embedded (e.g. VC3D via FetchContent): tools/TIFF/tests default OFF, only
# the library builds, with no dependencies beyond the vendored sources:
#     FetchContent_MakeAvailable(scrollfiesta)  ->  scrollfiesta::scrollfiesta
#
# OpenMP note: this project NEVER calls find_package(OpenMP) and no library
# TU includes <omp.h> (see ves_platform.c). Embedding hosts (VC3D) alias
# OpenMP:: targets to no-op stubs and inject a C++-only omp.h include dir;
# raw -fopenmp//openmp flags keep us independent of all of that.

cmake_minimum_required(VERSION 3.24)

project(scrollfiesta VERSION 0.10.0 LANGUAGES C CXX)

# ---------------------------------------------------------------------------
# Options (standalone builds get everything; embedded builds get the library)
# ---------------------------------------------------------------------------
option(SCROLLFIESTA_BUILD_TOOLS
       "Build the supported CLI front ends plus small diagnostic tools"
       ${PROJECT_IS_TOP_LEVEL})
option(SCROLLFIESTA_BUILD_TESTS
       "Build the ctest suite"
       ${PROJECT_IS_TOP_LEVEL})
option(SCROLLFIESTA_WITH_TIFF
       "Enable TIFF file I/O (staged vendored libtiff on Windows, find_package elsewhere)"
       ${PROJECT_IS_TOP_LEVEL})
option(SCROLLFIESTA_OPENMP
       "Compile the hot loops with OpenMP (raw flags, no find_package)"
       ON)
option(SCROLLFIESTA_INSTALL
       "Generate install/export rules"
       ${PROJECT_IS_TOP_LEVEL})

if(SCROLLFIESTA_BUILD_TOOLS AND NOT SCROLLFIESTA_WITH_TIFF)
    message(FATAL_ERROR "SCROLLFIESTA_BUILD_TOOLS requires SCROLLFIESTA_WITH_TIFF=ON "
                        "(the CLI tools read/write prediction TIFFs).")
endif()

# GCC 16.1.0 (current MSYS2 UCRT64 and winlibs releases) miscompiles correct,
# UB-free interleaved-store loops at -O3 when strict aliasing is enabled:
# scripts/gcc16_aliasing_repro.c builds a triangle grid whose face indices come
# out wrong under -O3 but right under -O3 -fno-strict-aliasing or -O2. Mesh
# code here is exactly that shape, and the MSVC baseline never exploits strict
# aliasing, so disable it for GNU project-wide (directory scope: library,
# tools, and tests) until the upstream fix lands.
add_compile_options(
    "$<$<AND:$<COMPILE_LANGUAGE:C>,$<C_COMPILER_ID:GNU>>:-fno-strict-aliasing>"
    "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:GNU>>:-fno-strict-aliasing>")

# The canonical library artifact is the shared library (consumed via dlopen +
# sf_get_api). Standalone builds default to shared unless the user says
# otherwise; embedded builds inherit the host's BUILD_SHARED_LIBS.
if(PROJECT_IS_TOP_LEVEL AND NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS ON)
endif()

if(SCROLLFIESTA_BUILD_TESTS)
    enable_testing()
endif()

set(SF_SRC  ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SF_DEPS ${CMAKE_CURRENT_SOURCE_DIR}/deps)

# ---------------------------------------------------------------------------
# Vendored dependency: Shewchuk Triangle (CDT for hole filling).
# Compiled from source with the exact defines of build-deps.sh / the .vcxproj.
# NB: Triangle has its own (non-commercial-distribution) license — see
# THIRD_PARTY_LICENSES.md.
# ---------------------------------------------------------------------------
add_library(sf_triangle STATIC ${SF_DEPS}/src/triangle/triangle.c)
target_compile_definitions(sf_triangle PRIVATE
    TRILIBRARY ANSI_DECLARATORS NO_TIMER REAL=double VOID=int)
target_include_directories(sf_triangle PUBLIC ${SF_DEPS}/include)
set_target_properties(sf_triangle PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(MSVC)
    target_compile_options(sf_triangle PRIVATE /W1 /wd4311 /wd4312)
    target_compile_definitions(sf_triangle PRIVATE _CRT_SECURE_NO_WARNINGS)
else()
    target_compile_options(sf_triangle PRIVATE -w)
endif()

# ---------------------------------------------------------------------------
# Vendored dependency: Clipper2 (+Z variant) via its own CMake.
# Always static and folded into our library, regardless of BUILD_SHARED_LIBS.
# ---------------------------------------------------------------------------
set(CLIPPER2_UTILS    OFF CACHE BOOL   "" FORCE)
set(CLIPPER2_EXAMPLES OFF CACHE BOOL   "" FORCE)
set(CLIPPER2_TESTS    OFF CACHE BOOL   "" FORCE)
set(CLIPPER2_USINGZ   "ON" CACHE STRING "" FORCE)
set(_sf_saved_build_shared ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(${SF_DEPS}/src/Clipper2/CPP clipper2 EXCLUDE_FROM_ALL)
set(BUILD_SHARED_LIBS ${_sf_saved_build_shared})
set_target_properties(Clipper2 Clipper2Z PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ---------------------------------------------------------------------------
# TIFF (only when SCROLLFIESTA_WITH_TIFF).
# Windows: prefer the static libs staged by build-deps.ps1 into
# deps/lib/win64/<Config>/ (tiff.lib + zs.lib), matching the .sln build.
# Otherwise (Linux/macOS/vcpkg): find_package(TIFF).
# ---------------------------------------------------------------------------
set(SF_TIFF_LIBS "")
if(SCROLLFIESTA_WITH_TIFF)
    if(WIN32 AND EXISTS ${SF_DEPS}/lib/win64/Release/tiff.lib)
        add_library(sf_tiff_staged STATIC IMPORTED)
        set_target_properties(sf_tiff_staged PROPERTIES
            IMPORTED_LOCATION         ${SF_DEPS}/lib/win64/Release/tiff.lib
            IMPORTED_LOCATION_RELEASE ${SF_DEPS}/lib/win64/Release/tiff.lib
            IMPORTED_LOCATION_DEBUG   ${SF_DEPS}/lib/win64/Debug/tiff.lib
            MAP_IMPORTED_CONFIG_RELWITHDEBINFO Release
            MAP_IMPORTED_CONFIG_MINSIZEREL     Release
            INTERFACE_INCLUDE_DIRECTORIES
                "${SF_DEPS}/src/tiff-4.7.1/libtiff;${SF_DEPS}/src/tiff-4.7.1/bld/libtiff")
        add_library(sf_zlib_staged STATIC IMPORTED)
        set_target_properties(sf_zlib_staged PROPERTIES
            IMPORTED_LOCATION         ${SF_DEPS}/lib/win64/Release/zs.lib
            IMPORTED_LOCATION_RELEASE ${SF_DEPS}/lib/win64/Release/zs.lib
            IMPORTED_LOCATION_DEBUG   ${SF_DEPS}/lib/win64/Debug/zs.lib
            MAP_IMPORTED_CONFIG_RELWITHDEBINFO Release
            MAP_IMPORTED_CONFIG_MINSIZEREL     Release)
        set(SF_TIFF_LIBS sf_tiff_staged sf_zlib_staged)
        message(STATUS "scrollfiesta: using staged vendored libtiff (deps/lib/win64)")
    else()
        # Linux/macOS: build-deps.sh stages the vendored zlib+libtiff into
        # deps/prefix when the host has no libtiff development package. That
        # is the Windows arrangement, one directory later: look there first,
        # then fall back to whatever the host provides.
        if(EXISTS ${SF_DEPS}/prefix/include/tiffio.h)
            list(PREPEND CMAKE_PREFIX_PATH ${SF_DEPS}/prefix)
            # The static libtiff package config names ZLIB::ZLIB and
            # CMath::CMath in its link interface without requiring them, so
            # both imported targets have to exist before it is included.
            find_package(ZLIB QUIET)
            if(NOT TARGET CMath::CMath)
                add_library(CMath::CMath INTERFACE IMPORTED)
                find_library(SF_LIBM m)
                if(SF_LIBM)
                    set_target_properties(CMath::CMath PROPERTIES
                        INTERFACE_LINK_LIBRARIES ${SF_LIBM})
                endif()
            endif()
        endif()
        find_package(TIFF REQUIRED)
        set(SF_TIFF_LIBS TIFF::TIFF)
        message(STATUS "scrollfiesta: using find_package(TIFF)")
    endif()
endif()

# ---------------------------------------------------------------------------
# scrollfiesta_core — every algorithm module, as one static library.
# Source lists mirror src/Makefile (COMMON/EXTRACT/SPLIT/REMESH/HOLEFILL/
# TOPO/PIPELINE) plus the grid_weld-side WELD modules and the C++ wrappers.
# ---------------------------------------------------------------------------
set(SF_COMMON_SRC
    ${SF_SRC}/common/arena.c
    ${SF_SRC}/common/except.c
    ${SF_SRC}/common/csr.c
    ${SF_SRC}/common/kdtree.c
    ${SF_SRC}/common/union_find.c
    ${SF_SRC}/common/cc_color.c
    ${SF_SRC}/common/bfs.c
    ${SF_SRC}/common/pca.c
    ${SF_SRC}/common/mls_project.c
    ${SF_SRC}/common/vert_weld.c
    ${SF_SRC}/common/halo_loader.c
    ${SF_SRC}/common/raw_sample.c
    ${SF_SRC}/common/closest_tri.c
    ${SF_SRC}/common/mesh_manifold.c
    ${SF_SRC}/common/mesh_trim.c

    ${SF_SRC}/common/obj_io.c
    ${SF_SRC}/common/dump_obj.c
    ${SF_SRC}/common/ves_platform.c
    ${SF_SRC}/common/raw_snapback.c
    ${SF_SRC}/common/snap_cg.c
    ${SF_SRC}/common/eig3.c
    ${SF_SRC}/common/dev_metric.c
    ${SF_SRC}/common/run_ctx.c
    ${SF_SRC}/common/mesh_bin.c
    ${SF_SRC}/common/mesh_normals.c
    ${SF_SRC}/common/pitch_table.c)
if(SCROLLFIESTA_WITH_TIFF)
    list(APPEND SF_COMMON_SRC ${SF_SRC}/common/tiff_io.c)
else()
    list(APPEND SF_COMMON_SRC ${SF_SRC}/common/tiff_io_stub.c)
endif()

set(SF_EXTRACT_SRC
    ${SF_SRC}/extract/mesh_extract.c
    ${SF_SRC}/extract/mesh_resplit.c
    ${SF_SRC}/extract/pred_reject.c)
set(SF_SPLIT_SRC
    ${SF_SRC}/split/oracle.c
    ${SF_SRC}/split/bridge_cut.c
    ${SF_SRC}/split/overlap_sep.c
    ${SF_SRC}/split/overlap_merge.c
    ${SF_SRC}/split/overlap_quality.c
    ${SF_SRC}/split/depth_peel.c
    ${SF_SRC}/split/dev_cut.c)
set(SF_REMESH_SRC
    ${SF_SRC}/remesh/ball_pivot.c
    ${SF_SRC}/remesh/patch_repair.c
    ${SF_SRC}/remesh/cvt_remesh.c
    ${SF_SRC}/remesh/rvd_clip.c
    ${SF_SRC}/remesh/orient_mesh.c
    ${SF_SRC}/remesh/pinhole_fill.c
    ${SF_SRC}/remesh/normal_orient.c
    ${SF_SRC}/remesh/component_cull.c
    ${SF_SRC}/remesh/manifold_guard.c
    ${SF_SRC}/remesh/dev_gate.c
    ${SF_SRC}/remesh/intersection_cleanup.c
    ${SF_SRC}/remesh/chart_lineage.c
    ${SF_SRC}/remesh/chart_zipper.c
    ${SF_SRC}/remesh/chart_bridge_forest.c)
set(SF_WELD_SRC
    ${SF_SRC}/remesh/seam_weld.c
    ${SF_SRC}/remesh/seam_planes.c
    ${SF_SRC}/remesh/seam_refine.c
    ${SF_SRC}/remesh/seam_band_cvt.c
    ${SF_SRC}/remesh/sheet_reweld.c
    ${SF_SRC}/remesh/seam_hole_fill.c
    ${SF_SRC}/remesh/fold_cleanup.c
    ${SF_SRC}/remesh/orient_weld.c
    ${SF_SRC}/remesh/weld_cleanup.c
    ${SF_SRC}/remesh/boundary_arc_weld.c
    ${SF_SRC}/remesh/stream_weld_format.c)
set(SF_HOLEFILL_SRC ${SF_SRC}/holefill/hole_fill.c)
set(SF_TOPO_SRC
    ${SF_SRC}/topology/mesh_topo.c
    ${SF_SRC}/topology/seam_cut.c
    ${SF_SRC}/topology/developability.c)
set(SF_FLATTEN_SRC
    ${SF_SRC}/flatten/surface_snap.c
    ${SF_SRC}/flatten/snap_darkcut.c
    ${SF_SRC}/flatten/snap_quilt.c
    ${SF_SRC}/flatten/snap_solve.c
    ${SF_SRC}/flatten/recto_refine.c
    ${SF_SRC}/flatten/disk_topology_repair.c
    ${SF_SRC}/flatten/topology_invariants.c
    ${SF_SRC}/flatten/ribbon_grid_io.c)
# Whole-scroll modules the per-cube tools also reach for (axis straightening is
# used by grid_weld and scroll_whole alike, and it warps along the fitted
# physical shaft of shaft_warp.c), so they belong to the core library rather
# than to one front end's source list.
set(SF_WHOLE_SRC
    ${SF_SRC}/whole/axis_warp.c
    ${SF_SRC}/whole/shaft_warp.c)

# ---------------------------------------------------------------------------
# TAUCS: the supernodal multifrontal Cholesky behind the ribbon/metric solves.
# scroll_whole is its only consumer, so the whole block — detection included —
# is skipped for library-only and embedded builds. A FetchContent consumer
# never pays for a BLAS it cannot reach.
# The source list mirrors marble_uv_repair.vcxproj, which owned it on the
# MSBuild side.  Only the double-precision LLT path is built, and
# taucs_stubs.c supplies the LAPACK entry points it would otherwise need, so
# no BLAS/CLAPACK dependency is pulled in.  Kept as its own static library and
# linked only where it is used, rather than folded into scrollfiesta_core.
# ---------------------------------------------------------------------------
set(SF_HAVE_LAPACK FALSE)
if(SCROLLFIESTA_BUILD_TOOLS)
set(SF_TAUCS_DIR ${SF_DEPS}/taucs)
set(SF_TAUCS_SRC
    ${SF_TAUCS_DIR}/build/wrap/taucs_complex_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_complex_D.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_logging_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_timer_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_malloc_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_ccs_base_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_ccs_base_D.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_vec_base_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_vec_base_D.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_ccs_ops_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_ccs_ops_D.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_ccs_order_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_sn_llt_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_sn_llt_D.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_ccs_factor_llt_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_ccs_factor_llt_D.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_ccs_solve_llt_GEN.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_ccs_solve_llt_D.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_iter_D.c
    ${SF_TAUCS_DIR}/build/wrap/taucs_linsolve_GEN.c
    ${SF_TAUCS_DIR}/external/src/genmmd.c
    ${SF_TAUCS_DIR}/build/taucs_stubs.c)

# TAUCS needs a BLAS/LAPACK underneath (dgemm_, dpotrf_, dtrsm_, dsyrk_,
# dnrm2_).  The MSBuild path links the vendored CLAPACK built out of
# deps/clapack/clapack.sln, which build-deps does not stage; a system LAPACK
# serves just as well.  If neither is present we skip TAUCS and the one front
# end that needs it rather than failing the whole build.
find_library(SF_CLAPACK_LIB NAMES clapack_nowrap clapack
             PATHS ${SF_DEPS}/clapack/LIB/x64 ${SF_DEPS}/lib/win64 NO_DEFAULT_PATH)
find_library(SF_BLAS_LIB    NAMES BLAS_nowrap blas
             PATHS ${SF_DEPS}/clapack/LIB/x64 ${SF_DEPS}/lib/win64 NO_DEFAULT_PATH)
find_library(SF_F2C_LIB     NAMES libf2c f2c
             PATHS ${SF_DEPS}/clapack/LIB/x64 ${SF_DEPS}/lib/win64 NO_DEFAULT_PATH)

if(SF_CLAPACK_LIB AND SF_BLAS_LIB AND SF_F2C_LIB)
    set(SF_LAPACK_LIBS ${SF_CLAPACK_LIB} ${SF_BLAS_LIB} ${SF_F2C_LIB})
    set(SF_HAVE_LAPACK TRUE)
    message(STATUS "scrollfiesta: TAUCS enabled against vendored CLAPACK")
else()
    find_package(LAPACK QUIET)
    if(LAPACK_FOUND)
        set(SF_LAPACK_LIBS ${LAPACK_LIBRARIES})
        set(SF_HAVE_LAPACK TRUE)
        message(STATUS "scrollfiesta: TAUCS enabled against system LAPACK")
    elseif(EXISTS ${SF_DEPS}/clapack/CMakeLists.txt)
        # Last resort: build the vendored CLAPACK from source. Correct
        # everywhere and needs nothing installed, but it is ~1,650 translation
        # units, so it is tried only after a staged or system LAPACK.
        add_subdirectory(${SF_DEPS}/clapack ${CMAKE_BINARY_DIR}/clapack EXCLUDE_FROM_ALL)
        set(SF_LAPACK_LIBS clapack_nowrap BLAS_nowrap libf2c)
        set(SF_HAVE_LAPACK TRUE)
        message(STATUS "scrollfiesta: TAUCS enabled against the vendored CLAPACK source")
    else()
        set(SF_HAVE_LAPACK FALSE)
        message(STATUS
            "scrollfiesta: no BLAS/LAPACK found -- skipping TAUCS and the "
            "scroll_whole front end. Install a system LAPACK, or run "
            "build-deps, to enable them.")
    endif()
endif()

if(SF_HAVE_LAPACK)
add_library(sf_taucs STATIC ${SF_TAUCS_SRC} ${SF_SRC}/flatten/sparse_solve.c)
target_link_libraries(sf_taucs PUBLIC ${SF_LAPACK_LIBS})
target_include_directories(sf_taucs PUBLIC
    ${SF_TAUCS_DIR}/src
    ${SF_TAUCS_DIR}/build
    ${SF_TAUCS_DIR}/external/src
    PRIVATE ${SF_SRC})
set_target_properties(sf_taucs PROPERTIES C_STANDARD 17 C_EXTENSIONS ON)
if(MSVC)
    # Vendored third-party C; do not hold it to the project warning level.
    target_compile_options(sf_taucs PRIVATE /W1)
    target_compile_definitions(sf_taucs PRIVATE
        _CRT_SECURE_NO_WARNINGS _CRT_NONSTDC_NO_DEPRECATE)
else()
    target_compile_options(sf_taucs PRIVATE -w)
endif()
# sparse_solve.c serializes TAUCS under an OpenMP critical section, and its
# callers run in parallel regions; compiled without the OpenMP flag (as the
# MSVC projects never are) that lock would vanish.
if(SCROLLFIESTA_OPENMP)
    if(MSVC)
        target_compile_options(sf_taucs PRIVATE /openmp)
    else()
        target_compile_options(sf_taucs PRIVATE -fopenmp)
        target_link_options(sf_taucs INTERFACE -fopenmp)
    endif()
endif()
endif()  # SF_HAVE_LAPACK
endif()  # SCROLLFIESTA_BUILD_TOOLS

set(SF_PIPELINE_SRC ${SF_SRC}/pipeline/pipeline_cube.c)
set(SF_CXX_WRAP_SRC
    ${SF_SRC}/holefill/clipper2_wrap.cpp
    ${SF_SRC}/split/multicut_wrap.cpp
    ${SF_SRC}/common/gco_wrap.cpp
    ${SF_SRC}/common/attene_predicates_wrap.cpp
    ${SF_DEPS}/src/attene-predicates/uinteger.cpp
    ${SF_DEPS}/src/attene-predicates/ieee_binary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/docs/research/gco-v3.0/GCoptimization.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/docs/research/gco-v3.0/LinkedBlockList.cpp)

add_library(scrollfiesta_core STATIC
    ${SF_COMMON_SRC} ${SF_EXTRACT_SRC} ${SF_SPLIT_SRC} ${SF_REMESH_SRC}
    ${SF_WELD_SRC} ${SF_HOLEFILL_SRC} ${SF_TOPO_SRC} ${SF_FLATTEN_SRC}
    ${SF_WHOLE_SRC} ${SF_PIPELINE_SRC} ${SF_CXX_WRAP_SRC})
target_include_directories(scrollfiesta_core PUBLIC
    ${SF_SRC}
    ${SF_DEPS}/include
    ${SF_DEPS}/src/graph/include
    PRIVATE
    ${SF_DEPS}/src/attene-predicates
    ${CMAKE_CURRENT_SOURCE_DIR}/docs/research/gco-v3.0)
target_link_libraries(scrollfiesta_core PRIVATE
    sf_triangle Clipper2 Clipper2Z ${SF_TIFF_LIBS})
set_target_properties(scrollfiesta_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_STANDARD 17
    C_EXTENSIONS ON
    CXX_STANDARD 17
    CXX_EXTENSIONS OFF)

if(WIN32)
    target_compile_definitions(scrollfiesta_core PRIVATE _CRT_SECURE_NO_WARNINGS NOMINMAX)
else()
    target_compile_definitions(scrollfiesta_core PRIVATE _GNU_SOURCE)
    target_link_libraries(scrollfiesta_core PUBLIC m)
endif()
if(NOT MSVC)
    # The folded-in C++ wrappers need the C++ runtime at the final link of C
    # consumers. MSVC objects self-describe their runtime; gcc-family linkers
    # (including MinGW, which is WIN32 but not MSVC) do not.
    target_link_libraries(scrollfiesta_core PRIVATE stdc++)
endif()

if(MSVC)
    target_compile_options(scrollfiesta_core PRIVATE /W4)
    # seam_cut.c ships at warning level 1 in the .vcxproj.
    set_source_files_properties(${SF_SRC}/topology/seam_cut.c
                                PROPERTIES COMPILE_OPTIONS /W1)
    # Keep project code at /W4 without inheriting warning noise from the
    # vendored third-party GCoptimization 3.0 sources.
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/docs/research/gco-v3.0/GCoptimization.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/docs/research/gco-v3.0/LinkedBlockList.cpp
        PROPERTIES COMPILE_OPTIONS /W1)
else()
    target_compile_options(scrollfiesta_core PRIVATE
        $<$<COMPILE_LANGUAGE:C>:-Wall;-Wextra;-Wshadow;-Wformat=2;-Werror=implicit-function-declaration>)
    # andres/graph uses std::numeric_limits without including <limits>;
    # GCC 13's libstdc++ no longer provides it transitively (see src/Makefile).
    set_source_files_properties(${SF_SRC}/split/multicut_wrap.cpp
                                PROPERTIES COMPILE_OPTIONS "-include;limits")
    # GCO v3.0 calls memset/memcpy/strcpy without <cstring> (its TU also
    # inlines graph.cpp/maxflow.cpp); MSVC's headers provide it transitively,
    # GCC 13's do not.
    set_source_files_properties(
        ${CMAKE_CURRENT_SOURCE_DIR}/docs/research/gco-v3.0/GCoptimization.cpp
        PROPERTIES COMPILE_OPTIONS "-include;cstring")
endif()

# OpenMP: raw flags only. Never find_package(OpenMP); no TU includes <omp.h>.
if(SCROLLFIESTA_OPENMP)
    if(MSVC)
        target_compile_options(scrollfiesta_core PRIVATE /openmp)
    else()
        target_compile_options(scrollfiesta_core PRIVATE -fopenmp)
        # Propagate the runtime to every final link that folds the core in.
        target_link_options(scrollfiesta_core INTERFACE -fopenmp)
    endif()
endif()

# ---------------------------------------------------------------------------
# scrollfiesta — the public library (include/scrollfiesta.h + src/api/).
# Shared or static per BUILD_SHARED_LIBS; the canonical embedding is a
# shared library consumed via dlopen + sf_get_api().
# ---------------------------------------------------------------------------
add_library(scrollfiesta
    ${SF_SRC}/api/scrollfiesta_api.c
    ${SF_SRC}/api/sf_mesh_convert.c
    ${SF_SRC}/api/sf_pipeline.c)
add_library(scrollfiesta::scrollfiesta ALIAS scrollfiesta)
target_include_directories(scrollfiesta PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(scrollfiesta PRIVATE scrollfiesta_core)
target_compile_definitions(scrollfiesta PRIVATE SCROLLFIESTA_BUILD)
set_target_properties(scrollfiesta PROPERTIES
    C_STANDARD 17
    C_EXTENSIONS ON
    POSITION_INDEPENDENT_CODE ON
    SOVERSION ${PROJECT_VERSION_MAJOR}
    VERSION ${PROJECT_VERSION})
get_target_property(_sf_lib_type scrollfiesta TYPE)
if(_sf_lib_type STREQUAL "SHARED_LIBRARY")
    target_compile_definitions(scrollfiesta PUBLIC SCROLLFIESTA_SHARED)
endif()
if(WIN32)
    target_compile_definitions(scrollfiesta PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

# ---------------------------------------------------------------------------
# CLI tools (top-level builds).
# ---------------------------------------------------------------------------
if(SCROLLFIESTA_BUILD_TOOLS)
    add_executable(cube_mesh     ${SF_SRC}/main.c)
    add_executable(grid_pipeline ${SF_SRC}/tools/grid_pipeline.c)
    add_executable(grid_weld     ${SF_SRC}/tools/grid_weld.c)
    foreach(tool cube_mesh grid_pipeline grid_weld)
        target_link_libraries(${tool} PRIVATE scrollfiesta_core)
        set_target_properties(${tool} PROPERTIES C_STANDARD 17 C_EXTENSIONS ON)
        if(WIN32)
            target_compile_definitions(${tool} PRIVATE _CRT_SECURE_NO_WARNINGS)
        else()
            target_compile_definitions(${tool} PRIVATE _GNU_SOURCE)
        endif()
    endforeach()

    # Supported whole-scroll path. These source lists mirror the corresponding
    # MSVC projects and deliberately replace the retired scroll_unwrap /
    # scroll_atlas / scroll_ribbon front ends.
    #
    # scroll_whole reaches the TAUCS Cholesky through flatten/ribbon.c, so it
    # is built only when a BLAS/LAPACK was located above.
    if(SF_HAVE_LAPACK)
    add_executable(scroll_whole
        ${SF_SRC}/tools/scroll_whole.c
        ${SF_SRC}/whole/cube_schedule.c
        ${SF_SRC}/whole/placed_cube.c
        ${SF_SRC}/whole/cube_register.c
        ${SF_SRC}/whole/group_graph.c
        ${SF_SRC}/whole/weld_track.c
        ${SF_SRC}/whole/weld_build.c
        ${SF_SRC}/whole/arc_reg.c
        ${SF_SRC}/common/maxflow.c
        ${SF_SRC}/flatten/ribbon.c
        ${SF_SRC}/flatten/winding_mrf.c
        ${SF_SRC}/flatten/lane_turn_mrf.c
        ${SF_SRC}/flatten/seam_cut.c
        ${SF_SRC}/flatten/mesh_topo.c
        ${SF_SRC}/common/arena.c
        ${SF_SRC}/common/except.c
        ${SF_SRC}/common/obj_io.c
        ${SF_SRC}/common/ves_platform.c
        ${SF_SRC}/common/pca.c
        ${SF_SRC}/common/csr.c
        ${SF_SRC}/common/union_find.c
        ${SF_SRC}/common/kdtree.c)
    target_link_libraries(scroll_whole PRIVATE sf_taucs)
    endif()  # SF_HAVE_LAPACK

    add_executable(scroll_unroll
        ${SF_SRC}/tools/scroll_unroll.c
        ${SF_SRC}/unroll/piece_set.c
        ${SF_SRC}/unroll/scaffold.c
        ${SF_SRC}/unroll/export_atlas.c
        ${SF_SRC}/unroll/scroll_raster.c
        ${SF_SRC}/unroll/strip_metrics.c
        ${SF_SRC}/unroll/strip_preview.c
        ${SF_SRC}/unroll/seam_own.c
        ${SF_SRC}/unroll/snap_grid.c
        ${SF_SRC}/unroll/tifxyz_export.c
        ${SF_SRC}/unroll/uv_weld.c
        ${SF_SRC}/unroll/relax_grid.c
        ${SF_SRC}/flatten/overlap_quilt.c
        ${SF_SRC}/flatten/snap_darkcut.c
        ${SF_SRC}/flatten/snap_quilt.c
        ${SF_SRC}/flatten/recto_refine.c
        ${SF_SRC}/flatten/ribbon_relax.c
        ${SF_SRC}/flatten/fiber_field.c
        ${SF_SRC}/flatten/rawtex_bake.c
        ${SF_SRC}/common/pca.c
        ${SF_SRC}/common/csr.c
        ${SF_SRC}/common/kdtree.c
        ${SF_SRC}/common/snap_cg.c
        ${SF_SRC}/common/gco_wrap.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/docs/research/gco-v3.0/GCoptimization.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/docs/research/gco-v3.0/LinkedBlockList.cpp
        ${SF_SRC}/split/multicut_wrap.cpp
        ${SF_SRC}/common/arena.c
        ${SF_SRC}/common/except.c
        ${SF_SRC}/common/obj_io.c
        ${SF_SRC}/common/ves_platform.c
        ${SF_SRC}/common/raw_sample.c
        ${SF_SRC}/common/tiff_io.c
        ${SF_SRC}/common/tif_strip.c
        ${SF_SRC}/common/ves_png.c)

    add_executable(atlas_overlap_fix
        ${SF_SRC}/tools/atlas_overlap_fix_tool.c
        ${SF_SRC}/tools/atlas_overlap_fix_sparse.c
        ${SF_SRC}/whole/atlas_collision_register.c
        ${SF_SRC}/whole/atlas_overlap_audit.c
        ${SF_SRC}/whole/monotone_qp.c
        ${SF_SRC}/whole/atlas_overlap_fix.c
        ${SF_SRC}/whole/atlas_register.c
        ${SF_SRC}/whole/atlas_solution.c
        ${SF_SRC}/unroll/piece_set.c
        ${SF_SRC}/unroll/scaffold.c
        ${SF_SRC}/common/arena.c
        ${SF_SRC}/common/except.c
        ${SF_SRC}/common/union_find.c
        ${SF_SRC}/common/obj_io.c
        ${SF_SRC}/common/raw_sample.c
        ${SF_SRC}/common/tiff_io.c
        ${SF_SRC}/common/ves_platform.c
        ${SF_SRC}/common/ves_png.c)

    # Terminal stage of the whole-scroll path: fits one collision-free ribbon
    # through the registered atlas that atlas_overlap_fix hands over as
    # atlas_solution.bin. Mirrors atlas_ribbon_fit.vcxproj.
    add_executable(atlas_ribbon_fit
        ${SF_SRC}/tools/atlas_ribbon_fit_tool.c
        ${SF_SRC}/tools/atlas_overlap_fix_sparse.c
        ${SF_SRC}/whole/atlas_ribbon_fit.c
        ${SF_SRC}/whole/atlas_ribbon_texture.c
        ${SF_SRC}/whole/atlas_overlap_audit.c
        ${SF_SRC}/whole/monotone_qp.c
        ${SF_SRC}/whole/atlas_solution.c
        ${SF_SRC}/unroll/piece_set.c
        ${SF_SRC}/unroll/scaffold.c
        ${SF_SRC}/common/arena.c
        ${SF_SRC}/common/except.c
        ${SF_SRC}/common/union_find.c
        ${SF_SRC}/common/obj_io.c
        ${SF_SRC}/common/raw_sample.c
        ${SF_SRC}/common/tiff_io.c
        ${SF_SRC}/common/ves_platform.c)

    set(SF_BIG_TOOLS scroll_unroll atlas_overlap_fix atlas_ribbon_fit)
    if(SF_HAVE_LAPACK)
        list(PREPEND SF_BIG_TOOLS scroll_whole)
    endif()

    foreach(tool ${SF_BIG_TOOLS})
        # These four carry standalone source lists mirroring their MSVC
        # projects, but still reach for shared core modules (axis warping, mesh
        # bin I/O, normals, intersection cleanup).  Link the core library so a
        # module added there does not have to be repeated in four places; a
        # static-library member is only pulled in when the symbol is still
        # undefined, so this cannot collide with a tool's own sources.
        target_link_libraries(${tool} PRIVATE scrollfiesta_core)
        target_include_directories(${tool} PRIVATE
            ${SF_SRC}
            ${SF_DEPS}/include
            ${SF_DEPS}/src/graph/include
            ${CMAKE_CURRENT_SOURCE_DIR}/docs/research/gco-v3.0)
        set_target_properties(${tool} PROPERTIES
            C_STANDARD 17 C_EXTENSIONS ON CXX_STANDARD 17 CXX_EXTENSIONS OFF)
        if(WIN32)
            target_compile_definitions(${tool} PRIVATE
                _CRT_SECURE_NO_WARNINGS _CRT_NONSTDC_NO_DEPRECATE NOMINMAX)
        else()
            target_compile_definitions(${tool} PRIVATE _GNU_SOURCE)
            target_link_libraries(${tool} PRIVATE m)
        endif()
    endforeach()
    target_link_libraries(scroll_unroll PRIVATE ${SF_TIFF_LIBS})
    target_link_libraries(atlas_overlap_fix PRIVATE ${SF_TIFF_LIBS})
    target_link_libraries(atlas_ribbon_fit PRIVATE ${SF_TIFF_LIBS})
    if(SCROLLFIESTA_OPENMP)
        set(SF_OPENMP_TOOLS scroll_unroll)
        if(SF_HAVE_LAPACK)
            list(APPEND SF_OPENMP_TOOLS scroll_whole)
        endif()
        foreach(tool ${SF_OPENMP_TOOLS})
            if(MSVC)
                target_compile_options(${tool} PRIVATE /openmp)
            else()
                target_compile_options(${tool} PRIVATE -fopenmp)
                target_link_options(${tool} PRIVATE -fopenmp)
            endif()
        endforeach()
    endif()

    # Diagnostic tools: single self-contained translation units (only libm).
    foreach(tool obj_components pinhole_verdict seam_audit)
        add_executable(${tool} ${SF_SRC}/tools/${tool}.c)
        target_link_libraries(${tool} PRIVATE scrollfiesta_core)
        set_target_properties(${tool} PROPERTIES C_STANDARD 17 C_EXTENSIONS ON)
        if(WIN32)
            target_compile_definitions(${tool} PRIVATE _CRT_SECURE_NO_WARNINGS)
        else()
            target_compile_definitions(${tool} PRIVATE _GNU_SOURCE)
            target_link_libraries(${tool} PRIVATE m)
        endif()
    endforeach()

    # -----------------------------------------------------------------------
    # The native sheet assembler (September 2026) and the tools that check,
    # texture and lay out its sheets. Source lists mirror the MSVC projects of
    # the same names; modules that are already in scrollfiesta_core come from
    # the library.
    # -----------------------------------------------------------------------
    function(sf_sheet_tool tool)
        target_link_libraries(${tool} PRIVATE scrollfiesta_core ${SF_TIFF_LIBS})
        target_include_directories(${tool} PRIVATE ${SF_SRC} ${SF_DEPS}/include)
        set_target_properties(${tool} PROPERTIES C_STANDARD 17 C_EXTENSIONS ON)
        if(WIN32)
            target_compile_definitions(${tool} PRIVATE
                _CRT_SECURE_NO_WARNINGS _CRT_NONSTDC_NO_DEPRECATE NOMINMAX)
        else()
            target_compile_definitions(${tool} PRIVATE _GNU_SOURCE)
            target_link_libraries(${tool} PRIVATE m)
        endif()
        if(SCROLLFIESTA_OPENMP)
            if(MSVC)
                target_compile_options(${tool} PRIVATE /openmp)
            else()
                target_compile_options(${tool} PRIVATE -fopenmp)
                target_link_options(${tool} PRIVATE -fopenmp)
            endif()
        endif()
        # The MSVC projects reserve 512 MB of stack: TAUCS recursion and the
        # SLIM/repair solves overflow the 1 MB default on large sheets. (On
        # Linux the main thread's stack is set by `ulimit -s`.)
        if(MSVC)
            target_link_options(${tool} PRIVATE /STACK:536870912)
        endif()
    endfunction()

    add_executable(ribbon_verdict ${SF_SRC}/tools/ribbon_verdict.c)
    add_executable(swirl_slim
        ${SF_SRC}/tools/swirl_slim.c
        ${SF_SRC}/flatten/slim_refine.c)
    add_executable(obj_bake_raw
        ${SF_SRC}/tools/obj_bake_raw.c
        ${SF_SRC}/flatten/rawtex_bake.c
        ${SF_SRC}/flatten/textured_obj.c
        ${SF_SRC}/flatten/marble_diag.c
        ${SF_SRC}/common/ves_png.c
        ${SF_SRC}/common/tif_strip.c)
    add_executable(sheet_strip
        ${SF_SRC}/tools/sheet_strip.c
        ${SF_SRC}/common/strip_layout.c
        ${SF_SRC}/common/json_read.c
        ${SF_SRC}/common/ves_png.c)
    add_executable(sheet_reading
        ${SF_SRC}/tools/sheet_reading.c
        ${SF_SRC}/common/strip_layout.c
        ${SF_SRC}/common/json_read.c
        ${SF_SRC}/common/ves_png.c)
    add_executable(sheet_layers
        ${SF_SRC}/tools/sheet_layers.c
        ${SF_SRC}/flatten/rawtex_layers.c
        ${SF_SRC}/flatten/rawtex_bake.c
        ${SF_SRC}/assemble/asm_axis.c
        ${SF_SRC}/common/ves_png.c
        ${SF_SRC}/common/tif_strip.c)
    add_executable(vmesh_tifxyz
        ${SF_SRC}/tools/vmesh_tifxyz.c
        ${SF_SRC}/flatten/tifxyz_write.c)
    add_executable(scroll_axis_track
        ${SF_SRC}/tools/scroll_axis_track.c
        ${SF_SRC}/whole/axis_track.c
        ${SF_SRC}/common/mesh_pile.c
        ${SF_SRC}/common/png_write.c)
    set(SF_SHEET_TOOLS ribbon_verdict swirl_slim obj_bake_raw sheet_strip
                       sheet_reading sheet_layers vmesh_tifxyz scroll_axis_track)

    # sheet_assemble solves its constrained QPs with the vendored OSQP 1.0
    # (built-in QDLDL + AMD; configuration in deps/osqp/osqp_configure.h) and
    # its linear systems with TAUCS, so it needs the BLAS/LAPACK found above.
    if(SF_HAVE_LAPACK)
        set(SF_OSQP_DIR ${SF_DEPS}/osqp)
        add_library(sf_osqp STATIC
            ${SF_OSQP_DIR}/src/algebra_libs.c
            ${SF_OSQP_DIR}/src/auxil.c
            ${SF_OSQP_DIR}/src/csc_math.c
            ${SF_OSQP_DIR}/src/csc_utils.c
            ${SF_OSQP_DIR}/src/error.c
            ${SF_OSQP_DIR}/src/kkt.c
            ${SF_OSQP_DIR}/src/matrix.c
            ${SF_OSQP_DIR}/src/osqp_api.c
            ${SF_OSQP_DIR}/src/polish.c
            ${SF_OSQP_DIR}/src/qdldl.c
            ${SF_OSQP_DIR}/src/qdldl_interface.c
            ${SF_OSQP_DIR}/src/scaling.c
            ${SF_OSQP_DIR}/src/util.c
            ${SF_OSQP_DIR}/src/vector.c
            ${SF_OSQP_DIR}/amd/src/SuiteSparse_config.c
            ${SF_OSQP_DIR}/amd/src/amd_1.c
            ${SF_OSQP_DIR}/amd/src/amd_2.c
            ${SF_OSQP_DIR}/amd/src/amd_aat.c
            ${SF_OSQP_DIR}/amd/src/amd_control.c
            ${SF_OSQP_DIR}/amd/src/amd_defaults.c
            ${SF_OSQP_DIR}/amd/src/amd_info.c
            ${SF_OSQP_DIR}/amd/src/amd_order.c
            ${SF_OSQP_DIR}/amd/src/amd_postorder.c
            ${SF_OSQP_DIR}/amd/src/amd_post_tree.c
            ${SF_OSQP_DIR}/amd/src/amd_preprocess.c
            ${SF_OSQP_DIR}/amd/src/amd_valid.c)
        target_include_directories(sf_osqp
            PUBLIC ${SF_OSQP_DIR} ${SF_OSQP_DIR}/inc/public
            PRIVATE ${SF_OSQP_DIR}/inc/private ${SF_OSQP_DIR}/amd/include)
        set_target_properties(sf_osqp PROPERTIES C_STANDARD 17 C_EXTENSIONS ON)
        if(MSVC)
            target_compile_options(sf_osqp PRIVATE /W1)
            target_compile_definitions(sf_osqp PRIVATE _CRT_SECURE_NO_WARNINGS)
        else()
            target_compile_options(sf_osqp PRIVATE -w)
            target_link_libraries(sf_osqp PUBLIC m)
        endif()

        # The assembler factors many systems at once. With MSVC the BLAS/LAPACK
        # underneath TAUCS is f2c output whose locals are static, so, as in
        # sheet_assemble.vcxproj, it links deps/clapack/TS (the routines TAUCS
        # reaches, with automatic locals) ahead of that library and calls TAUCS
        # without sparse_solve's lock. f2c's `integer` is `long`, which matches
        # the `int` TAUCS passes only on LLP64, so other toolchains keep the
        # ordinary sf_taucs, whose TAUCS calls stay serialized.
        if(MSVC)
            add_library(sf_taucs_mt STATIC ${SF_TAUCS_SRC} ${SF_SRC}/flatten/sparse_solve.c)
            target_compile_definitions(sf_taucs_mt PUBLIC SPARSE_TAUCS_THREAD_SAFE
                PRIVATE _CRT_SECURE_NO_WARNINGS _CRT_NONSTDC_NO_DEPRECATE)
            target_link_libraries(sf_taucs_mt PUBLIC ${SF_LAPACK_LIBS})
            target_include_directories(sf_taucs_mt PUBLIC
                ${SF_TAUCS_DIR}/src
                ${SF_TAUCS_DIR}/build
                ${SF_TAUCS_DIR}/external/src
                PRIVATE ${SF_SRC})
            set_target_properties(sf_taucs_mt PROPERTIES C_STANDARD 17 C_EXTENSIONS ON)
            target_compile_options(sf_taucs_mt PRIVATE /W1)
            if(SCROLLFIESTA_OPENMP)
                target_compile_options(sf_taucs_mt PRIVATE /openmp)
            endif()
            set(SF_ASM_TAUCS sf_taucs_mt)
            set(SF_CLAPACK_TS_SRC
                ${SF_DEPS}/clapack/TS/dgemm.c
                ${SF_DEPS}/clapack/TS/dsyrk.c
                ${SF_DEPS}/clapack/TS/dtrsm.c
                ${SF_DEPS}/clapack/TS/ddot.c
                ${SF_DEPS}/clapack/TS/dgemv.c
                ${SF_DEPS}/clapack/TS/dscal.c
                ${SF_DEPS}/clapack/TS/lsame.c
                ${SF_DEPS}/clapack/TS/dpotrf.c
                ${SF_DEPS}/clapack/TS/dpotf2.c
                ${SF_DEPS}/clapack/TS/ilaenv.c
                ${SF_DEPS}/clapack/TS/ieeeck.c
                ${SF_DEPS}/clapack/TS/iparmq.c)
            set_source_files_properties(${SF_CLAPACK_TS_SRC} PROPERTIES
                INCLUDE_DIRECTORIES ${SF_DEPS}/clapack/INCLUDE
                COMPILE_DEFINITIONS NO_BLAS_WRAP
                COMPILE_OPTIONS /W0)
        else()
            set(SF_ASM_TAUCS sf_taucs)
            set(SF_CLAPACK_TS_SRC "")
        endif()

        add_executable(sheet_assemble
            ${SF_SRC}/tools/sheet_assemble.c
            ${SF_SRC}/assemble/asm_run.c
            ${SF_SRC}/assemble/asm_flatten.c
            ${SF_SRC}/assemble/asm_clean.c
            ${SF_SRC}/assemble/asm_store.c
            ${SF_SRC}/assemble/asm_report.c
            ${SF_SRC}/assemble/asm_relate.c
            ${SF_SRC}/assemble/asm_pose.c
            ${SF_SRC}/assemble/asm_conflict.c
            ${SF_SRC}/assemble/asm_layer_cut.c
            ${SF_SRC}/assemble/asm_emit.c
            ${SF_SRC}/assemble/asm_verdict.c
            ${SF_SRC}/assemble/asm_place.c
            ${SF_SRC}/assemble/asm_continuity.c
            ${SF_SRC}/assemble/asm_boundary.c
            ${SF_SRC}/assemble/asm_metric.c
            ${SF_SRC}/assemble/asm_contacts.c
            ${SF_SRC}/assemble/asm_field.c
            ${SF_SRC}/assemble/asm_repair.c
            ${SF_SRC}/assemble/asm_discover.c
            ${SF_SRC}/assemble/asm_audit.c
            ${SF_SRC}/assemble/asm_axis.c
            ${SF_SRC}/assemble/asm_axis_derive.c
            ${SF_SRC}/assemble/asm_reading_order.c
            ${SF_SRC}/whole/atlas_winding_sync.c
            ${SF_SRC}/common/mesh_pile.c
            ${SF_SRC}/common/maxflow.c
            ${SF_SRC}/common/ves_png.c
            ${SF_SRC}/common/json_read.c
            ${SF_SRC}/flatten/convex_qp.c
            ${SF_SRC}/flatten/active_set_qp.c
            ${SF_SRC}/flatten/slim_refine.c
            ${SF_CLAPACK_TS_SRC})
        target_link_libraries(sheet_assemble PRIVATE sf_osqp ${SF_ASM_TAUCS})
        list(APPEND SF_SHEET_TOOLS sheet_assemble)
    endif()  # SF_HAVE_LAPACK

    foreach(tool ${SF_SHEET_TOOLS})
        sf_sheet_tool(${tool})
    endforeach()

    # Test programs of the maintained per-cube modules (each has its own main).
    add_executable(holefill_interior_test ${CMAKE_CURRENT_SOURCE_DIR}/scripts/holefill/holefill_interior_test.c)
    add_executable(pinhole_fill_test      ${CMAKE_CURRENT_SOURCE_DIR}/scripts/extract/pinhole_fill_test.c)
    add_executable(mesh_resplit_test      ${CMAKE_CURRENT_SOURCE_DIR}/scripts/extract/mesh_resplit_test.c)
    set(SF_MODULE_TESTS holefill_interior_test pinhole_fill_test mesh_resplit_test)
    foreach(tool ${SF_MODULE_TESTS})
        sf_sheet_tool(${tool})
    endforeach()
endif()

# ---------------------------------------------------------------------------
# Tests.
# ---------------------------------------------------------------------------
if(SCROLLFIESTA_BUILD_TESTS)
    # Existing shared-infrastructure suite (provides its own main()).
    add_executable(sf_test_common ${CMAKE_CURRENT_SOURCE_DIR}/scripts/test_common.c)
    target_link_libraries(sf_test_common PRIVATE scrollfiesta_core)
    set_target_properties(sf_test_common PROPERTIES C_STANDARD 17 C_EXTENSIONS ON)
    # TIFF-backed cases self-skip when the library is built with the stub.
    target_compile_definitions(sf_test_common PRIVATE
        SF_TEST_WITH_TIFF=$<BOOL:${SCROLLFIESTA_WITH_TIFF}>)
    if(WIN32)
        target_compile_definitions(sf_test_common PRIVATE _CRT_SECURE_NO_WARNINGS)
    else()
        target_compile_definitions(sf_test_common PRIVATE _GNU_SOURCE)
    endif()
    add_test(NAME sf_test_common COMMAND sf_test_common)

    # Direct-link smoke test of the public API (broad op coverage).
    add_executable(sf_api_smoke ${CMAKE_CURRENT_SOURCE_DIR}/tests/api_smoke.c)
    target_link_libraries(sf_api_smoke PRIVATE scrollfiesta)
    set_target_properties(sf_api_smoke PROPERTIES C_STANDARD 17 C_EXTENSIONS ON)
    if(WIN32)
        target_compile_definitions(sf_api_smoke PRIVATE _CRT_SECURE_NO_WARNINGS)
    else()
        target_compile_definitions(sf_api_smoke PRIVATE _GNU_SOURCE)
        target_link_libraries(sf_api_smoke PRIVATE m)
    endif()
    add_test(NAME sf_api_smoke COMMAND sf_api_smoke)

    # The canonical consumer test: links NOTHING from scrollfiesta, loads the
    # shared library at runtime and drives it through sf_get_api — exactly
    # what an embedding host (VC3D) does. Only meaningful for shared builds.
    if(_sf_lib_type STREQUAL "SHARED_LIBRARY")
        add_executable(sf_api_dlopen ${CMAKE_CURRENT_SOURCE_DIR}/tests/api_dlopen.c)
        target_include_directories(sf_api_dlopen PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/include)
        set_target_properties(sf_api_dlopen PROPERTIES C_STANDARD 17 C_EXTENSIONS ON)
        if(WIN32)
            target_compile_definitions(sf_api_dlopen PRIVATE _CRT_SECURE_NO_WARNINGS)
        else()
            target_link_libraries(sf_api_dlopen PRIVATE ${CMAKE_DL_LIBS})
        endif()
        add_test(NAME sf_api_dlopen
                 COMMAND sf_api_dlopen $<TARGET_FILE:scrollfiesta>)
    endif()

    # Every sheet tool carries its module selftests (--selftest). They write
    # scratch files into the working directory, so each gets its own. Leak
    # checking stays with the library tests above: under the sanitizer CI job
    # these whole-tool suites should fail on memory errors, not on process-
    # lifetime allocations still held at exit.
    if(SCROLLFIESTA_BUILD_TOOLS)
        foreach(tool ${SF_SHEET_TOOLS})
            set(_sf_wd ${CMAKE_CURRENT_BINARY_DIR}/selftest/${tool})
            file(MAKE_DIRECTORY ${_sf_wd})
            add_test(NAME ${tool}_selftest COMMAND ${tool} --selftest
                     WORKING_DIRECTORY ${_sf_wd})
            set_tests_properties(${tool}_selftest PROPERTIES
                ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
        endforeach()
        foreach(tool ${SF_MODULE_TESTS})
            set(_sf_wd ${CMAKE_CURRENT_BINARY_DIR}/selftest/${tool})
            file(MAKE_DIRECTORY ${_sf_wd})
            add_test(NAME ${tool} COMMAND ${tool} WORKING_DIRECTORY ${_sf_wd})
            set_tests_properties(${tool} PROPERTIES
                ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
        endforeach()
        if(TARGET sheet_assemble)
            add_test(NAME sheet_assemble_cli
                     COMMAND ${CMAKE_COMMAND}
                         -DSHEET_ASSEMBLE=$<TARGET_FILE:sheet_assemble>
                         -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                         -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/selftest/sheet_cli
                         -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/sheet_cli.cmake)
            set_tests_properties(sheet_assemble_cli PROPERTIES
                ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
            # The assembler's suite runs whole miniature assemblies (about 100 s
            # on 8 threads).
            set_tests_properties(sheet_assemble_selftest PROPERTIES TIMEOUT 3600)
        endif()
    endif()
endif()

# ---------------------------------------------------------------------------
# Install / export.
# ---------------------------------------------------------------------------
if(SCROLLFIESTA_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)

    install(FILES include/scrollfiesta.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    install(FILES LICENSE THIRD_PARTY_LICENSES.md
            DESTINATION ${CMAKE_INSTALL_DATADIR}/scrollfiesta)

    if(_sf_lib_type STREQUAL "SHARED_LIBRARY")
        # Shared: PRIVATE-linked internals are folded into the binary, so the
        # export set is self-contained and find_package(scrollfiesta) works.
        install(TARGETS scrollfiesta
                EXPORT scrollfiestaTargets
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
        install(EXPORT scrollfiestaTargets
                NAMESPACE scrollfiesta::
                DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scrollfiesta)
        write_basic_package_version_file(
            ${CMAKE_CURRENT_BINARY_DIR}/scrollfiesta-config-version.cmake
            COMPATIBILITY SameMajorVersion)
        file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/scrollfiesta-config.cmake
             "include(\${CMAKE_CURRENT_LIST_DIR}/scrollfiestaTargets.cmake)\n")
        install(FILES
                ${CMAKE_CURRENT_BINARY_DIR}/scrollfiesta-config.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/scrollfiesta-config-version.cmake
                DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/scrollfiesta)
    else()
        # Static: the archive would need every vendored dep exported alongside
        # it; install just the archive + header, no CMake package. (Static
        # embedders should use add_subdirectory/FetchContent instead.)
        install(TARGETS scrollfiesta ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
        message(STATUS "scrollfiesta: static build — installing archive+header only, no CMake package export")
    endif()

    if(SCROLLFIESTA_BUILD_TOOLS)
        install(FILES configs/default.json configs/axis.csv
                DESTINATION ${CMAKE_INSTALL_DATADIR}/scrollfiesta/configs)
        set(SF_INSTALL_TOOLS cube_mesh grid_pipeline grid_weld scroll_unroll
                             atlas_overlap_fix atlas_ribbon_fit
                             obj_components pinhole_verdict seam_audit
                             ${SF_SHEET_TOOLS})
        if(SF_HAVE_LAPACK)
            list(APPEND SF_INSTALL_TOOLS scroll_whole)
        endif()
        install(TARGETS ${SF_INSTALL_TOOLS}
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()
endif()
