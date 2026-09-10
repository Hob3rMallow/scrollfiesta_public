#!/usr/bin/env bash
#
# build-deps.sh -- Build the vendored C/C++ dependencies for ScrollFiesta (Linux).
#
# Linux analog of build-deps.ps1. Produces the three static libraries that
# src/Makefile links from deps/lib/:
#
#     deps/lib/libtriangle.a    (Shewchuk's Triangle CDT, compiled from source)
#     deps/lib/libClipper2.a    (Clipper2 polygon boolean ops)
#     deps/lib/libClipper2Z.a   (Clipper2 with the Z coordinate, USINGZ)
#
# Unlike the Windows build, zlib and libtiff are NOT built here: the Makefile's
# `-ltiff` links the *system* libtiff (libtiff-dev; headers under the multiarch
# /usr/include/<triple>/ that GCC searches by default). Install it with:
#
#     sudo apt-get install libtiff-dev
#
# Run this ONCE before building the binaries:
#
#     ./build-deps.sh
#     cd src && make release -j"$(nproc)"
#
# Requires: gcc, g++, cmake (>= 3.24), make, ar.

set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
deps="$root/deps/src"
lib="$root/deps/lib"
jobs="$(nproc 2>/dev/null || echo 4)"

mkdir -p "$lib"

cyan() { printf '\033[36m==> %s\033[0m\n' "$*"; }
green() { printf '\033[32m%s\033[0m\n' "$*"; }

# --- Triangle ---------------------------------------------------------------
# Shewchuk's Triangle. The Windows build (cube_mesh.vcxproj / grid_weld.vcxproj)
# compiles triangle.c directly with these exact defines at WarningLevel1:
#   TRILIBRARY  -- build as a callable library (no main())
#   ANSI_DECLARATORS, NO_TIMER, REAL=double, VOID=int
# It is third-party numerical code, so suppress warnings (-w) the same way the
# .vcxproj drops it to the lowest warning level. NO_TIMER avoids the <sys/time.h>
# timing path; we deliberately do NOT define LINUX (no x87 FPU control word).
cyan "libtriangle.a"
tri_obj="$(mktemp -d)/triangle.o"
gcc -O2 -w \
    -DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER -DREAL=double -DVOID=int \
    -c "$deps/triangle/triangle.c" -o "$tri_obj"
ar rcs "$lib/libtriangle.a" "$tri_obj"
rm -rf "$(dirname "$tri_obj")"
echo "    staged libtriangle.a"

# --- Clipper2 ---------------------------------------------------------------
# CMake build of the static Clipper2 / Clipper2Z libraries. Tests git-clone
# googletest, so keep them off (as build-deps.ps1 does). Clipper2's CMakeLists
# hard-codes -Werror for non-MSVC; it is warning-clean under GCC 13, but pass
# -DCMAKE_COMPILE_WARNING_AS_ERROR=OFF as a cheap safety net for newer compilers.
cyan "libClipper2.a + libClipper2Z.a"
clipper_src="$deps/Clipper2/CPP"
clipper_bld="$clipper_src/build"
cmake -S "$clipper_src" -B "$clipper_bld" \
    -DCLIPPER2_TESTS=OFF -DCLIPPER2_EXAMPLES=OFF -DCLIPPER2_UTILS=OFF \
    -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_COMPILE_WARNING_AS_ERROR=OFF
cmake --build "$clipper_bld" --config Release -j"$jobs"

# Stage the two .a's (CMake names them libClipper2.a / libClipper2Z.a on Linux).
for name in libClipper2.a libClipper2Z.a; do
    found="$(find "$clipper_bld" -name "$name" -print -quit)"
    if [ -z "$found" ]; then
        echo "ERROR: $name not produced by the Clipper2 build" >&2
        exit 1
    fi
    cp -f "$found" "$lib/"
    echo "    staged $name"
done

# --- zlib + libtiff (only when the host does not provide them) --------------
# build-deps.ps1 builds both from deps/src on Windows. On Linux this script
# built neither, and its comment sent the user to `sudo apt-get install
# libtiff-dev`, so a host without that package (or without root) could not
# configure at all even though the sources are vendored right here. Build them
# into deps/prefix, which CMakeLists.txt picks up the way it already picks up
# deps/lib/win64 on Windows.
prefix="$root/deps/prefix"
if printf '#include <tiffio.h>\nint main(void){return 0;}\n' \
        | cc -x c - -o /dev/null -ltiff >/dev/null 2>&1; then
    echo "    system libtiff found; not building the vendored copy"
else
    cyan "vendored zlib + libtiff -> deps/prefix"
    # The two build directories are pure output. Starting them clean keeps this
    # step idempotent: zlib caches ZLIB_CONF_WRITTEN, so a half-finished build
    # directory makes the next configure fail on a missing zconf.h.cmakein.
    rm -rf "$deps/zlib/build-host" "$deps/tiff-4.7.1/bld-host"
    [ -f "$deps/zlib/zconf.h" ] || cp "$deps/zlib/zconf.h.in" "$deps/zlib/zconf.h"
    cmake -S "$deps/zlib" -B "$deps/zlib/build-host" \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DZLIB_BUILD_SHARED=OFF -DZLIB_BUILD_TESTING=OFF -DZLIB_INSTALL=ON \
        -DCMAKE_INSTALL_PREFIX="$prefix"
    cmake --build "$deps/zlib/build-host" -j"$jobs"
    cmake --install "$deps/zlib/build-host"
    cmake -S "$deps/tiff-4.7.1" -B "$deps/tiff-4.7.1/bld-host" \
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DBUILD_SHARED_LIBS=OFF -Dtiff-tools=OFF -Dtiff-tests=OFF \
        -Dtiff-contrib=OFF -Dtiff-docs=OFF \
        -DCMAKE_PREFIX_PATH="$prefix" -DCMAKE_INSTALL_PREFIX="$prefix"
    cmake --build "$deps/tiff-4.7.1/bld-host" -j"$jobs"
    cmake --install "$deps/tiff-4.7.1/bld-host"
    echo "    staged libz.a and libtiff.a into deps/prefix"
fi

echo
green "All dependencies staged into $lib"
ls -l "$lib"
