# ---------------------------------------------------------------------------
# Vendored CLAPACK (f2c-translated LAPACK + reference BLAS + the f2c runtime).
#
# TAUCS needs a BLAS/LAPACK underneath its supernodal Cholesky (dgemm_,
# dpotrf_, dtrsm_, dsyrk_, dnrm2_).  Upstream ships clapack.sln, but that is
# a VS2008 .vcproj solution modern MSBuild cannot consume without an
# interactive upgrade, so the build is expressed here instead.  The result is
# what the historical "without wrap" configuration produced:
#
#     BLAS_nowrap  clapack_nowrap  libf2c
#
# "Without wrap" means NO_BLAS_WRAP is defined, so the BLAS entry points keep
# their plain f2c names (dgemm_) rather than being macro-renamed to f2c_dgemm.
# TAUCS calls the plain names, so the wrap must stay off.
#
# Consumed two ways: add_subdirectory from the top-level build, or built
# standalone by build-deps.ps1, which stages the three archives into
# deps/lib/win64 where the top-level find_library picks them up.
#
# These are third-party sources; warnings are silenced rather than fixed.
# ---------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.24)
project(sf_clapack C)

set(CL ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SF_BLAS_SRC   ${CL}/BLAS/SRC/*.c)
file(GLOB SF_LAPACK_SRC ${CL}/SRC/*.c)
file(GLOB SF_F2C_SRC    ${CL}/F2CLIBS/libf2c/*.c)

# Trim the glob to what upstream's own libf2c.vcproj builds, plus main.c:
#   main.c      the Fortran program entry point; in a library it would fight
#               the application's own main.
#   arithchk.c  build-time generator for arith.h, which ships pre-generated.
#   the rest    the quad/64-bit integer variants, which are not in upstream's
#               Windows build and do not compile there ('longint' is undefined,
#               so the declarations parse as K&R and MSVC rejects them).
list(REMOVE_ITEM SF_F2C_SRC
     ${CL}/F2CLIBS/libf2c/main.c
     ${CL}/F2CLIBS/libf2c/arithchk.c
     ${CL}/F2CLIBS/libf2c/ftell64_.c
     ${CL}/F2CLIBS/libf2c/pow_qq.c
     ${CL}/F2CLIBS/libf2c/qbitbits.c
     ${CL}/F2CLIBS/libf2c/qbitshft.c)

if(NOT SF_BLAS_SRC OR NOT SF_LAPACK_SRC OR NOT SF_F2C_SRC)
    message(FATAL_ERROR "sf_clapack: vendored CLAPACK sources not found under ${CL}")
endif()

add_library(libf2c        STATIC ${SF_F2C_SRC})
add_library(BLAS_nowrap   STATIC ${SF_BLAS_SRC})
add_library(clapack_nowrap STATIC ${SF_LAPACK_SRC})

foreach(target libf2c BLAS_nowrap clapack_nowrap)
    # Only INCLUDE goes on the search path. It carries f2c.h (byte-identical
    # to the libf2c copy) alongside clapack.h and blaswrap.h.
    #
    # F2CLIBS/libf2c must NOT be added: it ships a VC4.2-era math.h shim whose
    # body is "#include <math.h>", which resolves back to itself the moment its
    # own directory is searched, and the compiler dies at include depth 1024.
    # libf2c's genuinely private headers (arith.h, fio.h, signal1.h, sysdep1.h)
    # are all quoted includes, so they resolve next to the including source
    # without a -I.
    target_include_directories(${target} PUBLIC ${CL}/INCLUDE)
    target_compile_definitions(${target} PUBLIC NO_BLAS_WRAP)
    set_target_properties(${target} PROPERTIES
        C_STANDARD 11 C_EXTENSIONS ON
        POSITION_INDEPENDENT_CODE ON)
    if(MSVC)
        target_compile_options(${target} PRIVATE /w)
        # MSDOS/USE_CLOCK/NO_ONEXIT are exactly what upstream's Windows
        # project defines: they steer the timing routines onto clock() instead
        # of sys/times.h, and off the Unix-only unistd.h and onexit paths.
        #
        # _CRT_DECLARE_NONSTDC_NAMES=0 keeps the UCRT from declaring the
        # POSIX-spelled names. libf2c declares isatty() itself (fio.h,
        # s_paus.c), and the UCRT's dllimport-qualified declaration in
        # corecrt_io.h is a different linkage, which is a hard C2375 error.
        target_compile_definitions(${target} PRIVATE
            MSDOS USE_CLOCK NO_ONEXIT
            _CRT_SECURE_NO_WARNINGS _CRT_DECLARE_NONSTDC_NAMES=0)
    else()
        # f2c'd code trips a lot of modern diagnostics; none are actionable
        # here, and several are errors by default from GCC 14 on.
        target_compile_options(${target} PRIVATE
            -w -fno-strict-aliasing
            -Wno-implicit-function-declaration -Wno-int-conversion
            -Wno-incompatible-pointer-types -Wno-return-mismatch)
    endif()
endforeach()

target_link_libraries(clapack_nowrap PUBLIC BLAS_nowrap libf2c)
target_link_libraries(BLAS_nowrap PUBLIC libf2c)

# The f2c runtime wants the math library on Unix.
if(UNIX)
    target_link_libraries(libf2c PUBLIC m)
endif()
