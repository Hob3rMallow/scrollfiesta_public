#!/usr/bin/env bash
# Build scrollfiesta's cube_mesh + grid_weld natively on macOS (Apple Silicon).
#
# The repo Makefile targets Linux/GCC; this script builds the vendored deps
# (Triangle, Clipper2) into deps/lib and invokes make with clang + Homebrew
# libomp/libtiff, relaxing -Werror. Idempotent.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC="$REPO/src"
DEPS="$REPO/deps"
LIBDIR="$DEPS/lib"
mkdir -p "$LIBDIR"

LIBOMP="$(brew --prefix libomp)"
LIBTIFF="$(brew --prefix libtiff)"
echo "REPO=$REPO"
echo "libomp=$LIBOMP  libtiff=$LIBTIFF"

# ---- 1. Triangle -> libtriangle.a -----------------------------------------
if [ ! -f "$LIBDIR/libtriangle.a" ]; then
  echo "== building Triangle =="
  clang -O2 -w -DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER \
    -c "$DEPS/src/Triangle/triangle.c" -o "$LIBDIR/triangle.o"
  ar rcs "$LIBDIR/libtriangle.a" "$LIBDIR/triangle.o"
  rm -f "$LIBDIR/triangle.o"
fi

# ---- 2. Clipper2 -> libClipper2.a + libClipper2Z.a ------------------------
if [ ! -f "$LIBDIR/libClipper2.a" ] || [ ! -f "$LIBDIR/libClipper2Z.a" ]; then
  echo "== building Clipper2 =="
  CB="$DEPS/build-clipper2"
  cmake -S "$DEPS/src/Clipper2/CPP" -B "$CB" -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=OFF -DCLIPPER2_TESTS=OFF -DCLIPPER2_EXAMPLES=OFF \
    -DCLIPPER2_UTILS=OFF -DCLIPPER2_USINGZ=ON >/dev/null
  cmake --build "$CB" --config Release -j >/dev/null
  find "$CB" -name 'libClipper2*.a' -exec cp {} "$LIBDIR/" \;
fi
ls -1 "$LIBDIR"

# ---- 3. cube_mesh + grid_weld ---------------------------------------------
echo "== building cube_mesh + grid_weld =="
CFLAGS="-std=c99 -O2 -DNDEBUG -fno-common -Wall -Wno-error \
-Wno-implicit-function-declaration -Xpreprocessor -fopenmp \
-I$LIBOMP/include -I$LIBTIFF/include"
CXXFLAGS="-std=c++17 -O2 -w"
LDFLAGS="-lm -L$LIBTIFF/lib -ltiff -L$LIBOMP/lib -lomp -lpthread \
-Wl,-rpath,$LIBOMP/lib -Wl,-rpath,$LIBTIFF/lib"
EXT_LIB="-L$DEPS/lib -ltriangle -lClipper2 -lClipper2Z -lc++"

make -C "$SRC" cube_mesh grid_weld \
  CC=clang CXX=clang++ BUILD_FLAGS=" " \
  CFLAGS="$CFLAGS" CXXFLAGS="$CXXFLAGS" LDFLAGS="$LDFLAGS" EXT_LIB="$EXT_LIB"

echo "== done =="
ls -l "$SRC/cube_mesh" "$SRC/grid_weld"
