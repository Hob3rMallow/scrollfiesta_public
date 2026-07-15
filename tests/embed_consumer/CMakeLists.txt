# Synthetic embedding host that simulates VC3D's build environment:
# no-op OpenMP:: alias targets and a C++-only stub omp.h on the include path
# are installed BEFORE scrollfiesta is added, exactly like VC3D does on
# MSVC/Clang. If scrollfiesta ever calls find_package(OpenMP) or a library TU
# includes <omp.h> by name, this project stops building — which is the point.
#
# Configure with -DSCROLLFIESTA_CHECKOUT=<path to the scrollfiesta repo>.
cmake_minimum_required(VERSION 3.24)
project(scrollfiesta_embed_consumer LANGUAGES C CXX)

if(NOT SCROLLFIESTA_CHECKOUT)
    message(FATAL_ERROR "pass -DSCROLLFIESTA_CHECKOUT=<scrollfiesta repo root>")
endif()

# --- the hostile-host simulation -------------------------------------------
add_library(embed_openmp_stub INTERFACE)
add_library(OpenMP::OpenMP_CXX ALIAS embed_openmp_stub)
add_library(OpenMP::OpenMP_C   ALIAS embed_openmp_stub)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/omp_stub)

# --- consume scrollfiesta the way VC3D does --------------------------------
set(BUILD_SHARED_LIBS ON)
set(SCROLLFIESTA_BUILD_TOOLS OFF)
set(SCROLLFIESTA_BUILD_TESTS OFF)
set(SCROLLFIESTA_WITH_TIFF   OFF)
set(SCROLLFIESTA_INSTALL     OFF)

include(FetchContent)
FetchContent_Declare(scrollfiesta SOURCE_DIR ${SCROLLFIESTA_CHECKOUT})
FetchContent_MakeAvailable(scrollfiesta)

# The consumer binary links NOTHING from scrollfiesta — pure dlopen.
add_executable(embed_consumer ${SCROLLFIESTA_CHECKOUT}/tests/api_dlopen.c)
target_include_directories(embed_consumer PRIVATE
    ${SCROLLFIESTA_CHECKOUT}/include)
if(NOT WIN32)
    target_link_libraries(embed_consumer PRIVATE ${CMAKE_DL_LIBS})
endif()

enable_testing()
add_test(NAME embed_consumer
         COMMAND embed_consumer $<TARGET_FILE:scrollfiesta>)
