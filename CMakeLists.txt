# ScrollFiesta — CMake build.
#
# Builds the public C library (target `scrollfiesta`, alias
# `scrollfiesta::scrollfiesta`, header include/scrollfiesta.h) plus, when
# top-level, the six CLI tools that the .sln / src/Makefile also build.
# The legacy builds remain supported; this file mirrors their source lists
# and flags (see src/Makefile and cube_mesh.vcxproj, the CI-verified pair).
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

project(scrollfiesta VERSION 0.9.0 LANGUAGES C CXX)

# ---------------------------------------------------------------------------
# Options (standalone builds get everything; embedded builds get the library)
# ---------------------------------------------------------------------------
option(SCROLLFIESTA_BUILD_TOOLS
       "Build cube_mesh/grid_pipeline/grid_weld + diagnostic tools"
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
    ${SF_SRC}/common/bfs.c
    ${SF_SRC}/common/pca.c
    ${SF_SRC}/common/mls_project.c
    ${SF_SRC}/common/vert_weld.c
    ${SF_SRC}/common/halo_loader.c
    ${SF_SRC}/common/mesh_manifold.c
    ${SF_SRC}/common/mesh_trim.c
    ${SF_SRC}/common/ply_io.c
    ${SF_SRC}/common/obj_io.c
    ${SF_SRC}/common/dump_obj.c
    ${SF_SRC}/common/ves_platform.c
    ${SF_SRC}/common/qem.c
    ${SF_SRC}/common/raw_snapback.c
    ${SF_SRC}/common/snap_cg.c
    ${SF_SRC}/common/eig3.c
    ${SF_SRC}/common/dev_metric.c
    ${SF_SRC}/common/run_ctx.c)
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
    ${SF_SRC}/split/depth_peel.c
    ${SF_SRC}/split/dev_cut.c)
set(SF_REMESH_SRC
    ${SF_SRC}/remesh/ball_pivot.c
    ${SF_SRC}/remesh/orient_mesh.c
    ${SF_SRC}/remesh/pinhole_fill.c
    ${SF_SRC}/remesh/normal_orient.c
    ${SF_SRC}/remesh/component_cull.c
    ${SF_SRC}/remesh/manifold_guard.c
    ${SF_SRC}/remesh/dev_gate.c)
set(SF_WELD_SRC
    ${SF_SRC}/remesh/seam_weld.c
    ${SF_SRC}/remesh/seam_hole_fill.c
    ${SF_SRC}/remesh/fold_cleanup.c
    ${SF_SRC}/remesh/orient_weld.c
    ${SF_SRC}/remesh/weld_cleanup.c)
set(SF_HOLEFILL_SRC ${SF_SRC}/holefill/hole_fill.c)
set(SF_TOPO_SRC
    ${SF_SRC}/topology/mesh_topo.c
    ${SF_SRC}/topology/seam_cut.c
    ${SF_SRC}/topology/developability.c)
set(SF_PIPELINE_SRC ${SF_SRC}/pipeline/pipeline_cube.c)
set(SF_CXX_WRAP_SRC
    ${SF_SRC}/holefill/clipper2_wrap.cpp
    ${SF_SRC}/split/multicut_wrap.cpp)

add_library(scrollfiesta_core STATIC
    ${SF_COMMON_SRC} ${SF_EXTRACT_SRC} ${SF_SPLIT_SRC} ${SF_REMESH_SRC}
    ${SF_WELD_SRC} ${SF_HOLEFILL_SRC} ${SF_TOPO_SRC} ${SF_PIPELINE_SRC}
    ${SF_CXX_WRAP_SRC})
target_include_directories(scrollfiesta_core PUBLIC
    ${SF_SRC}
    ${SF_DEPS}/include
    ${SF_DEPS}/src/graph/include)
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
    # The two C++ wrappers are folded in; consumers of the static core need
    # the C++ runtime at final link.
    target_link_libraries(scrollfiesta_core PRIVATE stdc++)
endif()

if(MSVC)
    target_compile_options(scrollfiesta_core PRIVATE /W4)
    # seam_cut.c ships at warning level 1 in the .vcxproj.
    set_source_files_properties(${SF_SRC}/topology/seam_cut.c
                                PROPERTIES COMPILE_OPTIONS /W1)
else()
    target_compile_options(scrollfiesta_core PRIVATE
        $<$<COMPILE_LANGUAGE:C>:-Wall;-Wextra;-Wshadow;-Wformat=2;-Werror=implicit-function-declaration>)
    # andres/graph uses std::numeric_limits without including <limits>;
    # GCC 13's libstdc++ no longer provides it transitively (see src/Makefile).
    set_source_files_properties(${SF_SRC}/split/multicut_wrap.cpp
                                PROPERTIES COMPILE_OPTIONS "-include;limits")
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

    # Diagnostic tools: single self-contained translation units (only libm).
    foreach(tool obj_components pinhole_verdict seam_audit)
        add_executable(${tool} ${SF_SRC}/tools/${tool}.c)
        set_target_properties(${tool} PROPERTIES C_STANDARD 17 C_EXTENSIONS ON)
        if(WIN32)
            target_compile_definitions(${tool} PRIVATE _CRT_SECURE_NO_WARNINGS)
        else()
            target_compile_definitions(${tool} PRIVATE _GNU_SOURCE)
            target_link_libraries(${tool} PRIVATE m)
        endif()
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
    if(WIN32)
        target_compile_definitions(sf_test_common PRIVATE _CRT_SECURE_NO_WARNINGS)
    else()
        target_compile_definitions(sf_test_common PRIVATE _GNU_SOURCE)
    endif()
    add_test(NAME sf_test_common COMMAND sf_test_common)

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
        install(TARGETS cube_mesh grid_pipeline grid_weld
                        obj_components pinhole_verdict seam_audit
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()
endif()
