<#
  build-deps.ps1 -- Build the vendored C/C++ dependencies for ScrollFiesta.

  Builds zlib, Clipper2 and libtiff as static x64 libraries and stages them
  into deps\lib\win64\<Config>\ , which is where scrollfiesta.sln expects to
  find tiff.lib, zs.lib, Clipper2.lib and Clipper2Z.lib.

  Run this ONCE before building the solution. Requires CMake and Visual
  Studio 2022 (toolset v143). Afterwards, build the solution itself:

      msbuild scrollfiesta.sln -p:Configuration=Release -p:Platform=x64 -m

  Usage:
      pwsh -File build-deps.ps1                 # Release (default)
      pwsh -File build-deps.ps1 -Config Debug
#>
param(
    [ValidateSet('Release', 'Debug')] [string]$Config = 'Release',
    [string]$Generator = 'Visual Studio 17 2022'
)

$ErrorActionPreference = 'Stop'
$root  = $PSScriptRoot
$deps  = Join-Path $root 'deps\src'
$stage = Join-Path $root "deps\lib\win64\$Config"
New-Item -ItemType Directory -Force $stage | Out-Null

function Build-Dep {
    param([string]$Name, [string]$Src, [string]$Bld, [string[]]$CMakeArgs, [string]$LibGlob)
    Write-Host "==> $Name" -ForegroundColor Cyan
    cmake -S $Src -B $Bld -G $Generator -A x64 @CMakeArgs
    if ($LASTEXITCODE) { throw "$Name : configure failed" }
    cmake --build $Bld --config $Config
    if ($LASTEXITCODE) { throw "$Name : build failed" }
    Get-ChildItem -Recurse $Bld -Filter $LibGlob | ForEach-Object {
        Copy-Item $_.FullName $stage -Force
        Write-Host "    staged $($_.Name)"
    }
}

# --- zlib -------------------------------------------------------------------
# This repo's customized zlib CMakeLists reads zconf.h as a build input, so
# make sure it exists (restore it from the .in template if a clean ever removed
# it). Static lib comes out as "zs.lib" (OUTPUT_NAME suffix "s").
$zconf = Join-Path $deps 'zlib\zconf.h'
if (-not (Test-Path $zconf)) { Copy-Item (Join-Path $deps 'zlib\zconf.h.in') $zconf }
Build-Dep 'zlib' (Join-Path $deps 'zlib') (Join-Path $deps 'zlib\build') `
    @('-DZLIB_BUILD_SHARED=OFF', '-DZLIB_BUILD_TESTING=OFF', '-DZLIB_INSTALL=OFF') 'zs.lib'

# --- Clipper2 ---------------------------------------------------------------
# Disable tests (they git-clone googletest), examples and utils.
Build-Dep 'Clipper2' (Join-Path $deps 'Clipper2\CPP') (Join-Path $deps 'Clipper2\CPP\build') `
    @('-DCLIPPER2_TESTS=OFF', '-DCLIPPER2_EXAMPLES=OFF', '-DCLIPPER2_UTILS=OFF', '-DBUILD_SHARED_LIBS=OFF') 'Clipper2*.lib'

# --- libtiff ----------------------------------------------------------------
# Static, linked against the zlib we just built. The build directory MUST be
# named "bld": scrollfiesta.sln's include path hardcodes
# deps\src\tiff-4.7.1\bld\libtiff for the generated tif_config.h / tiffvers.h.
$zinc = Join-Path $deps 'zlib'
$zlib = Join-Path $deps "zlib\build\$Config\zs.lib"
Build-Dep 'libtiff' (Join-Path $deps 'tiff-4.7.1') (Join-Path $deps 'tiff-4.7.1\bld') `
    @('-DBUILD_SHARED_LIBS=OFF', '-Dtiff-tools=OFF', '-Dtiff-tests=OFF', '-Dtiff-contrib=OFF',
      '-Dtiff-docs=OFF', "-DZLIB_INCLUDE_DIR=$zinc", "-DZLIB_LIBRARY=$zlib") 'tiff.lib'

Write-Host "`nAll dependencies staged into $stage" -ForegroundColor Green
Get-ChildItem $stage -Filter *.lib | Select-Object Name, @{n = 'KB'; e = { [int]($_.Length / 1KB) } }
