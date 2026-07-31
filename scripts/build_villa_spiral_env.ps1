[CmdletBinding()]
param(
    [string]$VillaRoot = 'D:\work\villa\volume-cartographer',
    [string]$VcpkgRoot = 'D:\vcpkg',
    [string]$EnvironmentDir,
    [string]$BuildDir,
    [switch]$Offline,
    [switch]$RestoreVcpkg
)

$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$villa = (Resolve-Path $VillaRoot).Path
$vcpkg = (Resolve-Path $VcpkgRoot).Path
if (-not $EnvironmentDir) {
    $EnvironmentDir = Join-Path $repoRoot '.toolchains\villa-spiral-py314'
}
if (-not $BuildDir) {
    $BuildDir = Join-Path $repoRoot 'build\villa-python-surface-index'
}
$EnvironmentDir = [IO.Path]::GetFullPath($EnvironmentDir)
$BuildDir = [IO.Path]::GetFullPath($BuildDir)

$triplet = 'x64-windows-rel'
$installedRoot = Join-Path $villa 'build\windows-msvc\vcpkg_installed'
$triplets = Join-Path $villa 'triplets'
$prefix = Join-Path $installedRoot $triplet
$ceresConfig = Join-Path $prefix 'share\ceres\CeresConfig.cmake'
$vcpkgExe = Join-Path $vcpkg 'vcpkg.exe'
$toolchain = Join-Path $vcpkg 'scripts\buildsystems\vcpkg.cmake'

if ($RestoreVcpkg) {
    & $vcpkgExe install `
        --triplet $triplet `
        --vcpkg-root $vcpkg `
        --x-wait-for-lock `
        "--x-manifest-root=$villa" `
        "--x-install-root=$installedRoot" `
        "--overlay-triplets=$triplets"
    if ($LASTEXITCODE -ne 0) {
        throw "Villa vcpkg restore failed with exit code $LASTEXITCODE"
    }
}
if (-not (Test-Path -LiteralPath $ceresConfig)) {
    throw "Villa's $triplet dependency tree is incomplete; rerun with -RestoreVcpkg"
}

$env:UV_PROJECT_ENVIRONMENT = $EnvironmentDir
$syncArgs = @('sync', '--locked', '--project', $villa, '--package', 'spiral',
              '--no-install-workspace')
if ($Offline) {
    $syncArgs += '--offline'
}
& uv @syncArgs
if ($LASTEXITCODE -ne 0) {
    throw "Villa spiral dependency sync failed with exit code $LASTEXITCODE"
}

$python = Join-Path $EnvironmentDir 'Scripts\python.exe'
& uv pip install --python $python 'nanobind>=2.0'
if ($LASTEXITCODE -ne 0) {
    throw "nanobind install failed with exit code $LASTEXITCODE"
}

$cmakeArgs = @(
    '-S', $villa,
    '-B', $BuildDir,
    '-G', 'Visual Studio 17 2022',
    '-A', 'x64',
    "-DCMAKE_TOOLCHAIN_FILE=$toolchain",
    '-DVCPKG_MANIFEST_MODE=OFF',
    "-DVCPKG_INSTALLED_DIR=$installedRoot",
    "-DVCPKG_TARGET_TRIPLET=$triplet",
    "-DVCPKG_OVERLAY_TRIPLETS=$triplets",
    "-DCMAKE_PREFIX_PATH=$prefix",
    '-DVC_BUILD_PYTHON=ON',
    '-DVC_BUILD_APPS=OFF',
    '-DVC_BUILD_UI_TRACER=OFF',
    '-DVC_BUILD_FLATBOI=OFF',
    '-DVC_TESTING=OFF',
    '-DVC_WITH_SCROLLFIESTA=OFF',
    '-DVC_ENABLE_AMGX=OFF',
    "-DPython_EXECUTABLE=$python"
)
& cmake @cmakeArgs
if ($LASTEXITCODE -ne 0) {
    throw "Villa Python CMake configure failed with exit code $LASTEXITCODE"
}
& cmake --build $BuildDir --target vc_surface_index vc_volume --config Release -j 8
if ($LASTEXITCODE -ne 0) {
    throw "Villa Python binding build failed with exit code $LASTEXITCODE"
}

$sitePackages = (& $python -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])').Trim()
$target = Join-Path $sitePackages 'vc'
$moduleDir = Join-Path $BuildDir 'python\vc\Release'
$runtimeDir = Join-Path $BuildDir 'bin\Release'
New-Item -ItemType Directory -Path $target -Force | Out-Null
Copy-Item -LiteralPath (Join-Path $villa 'python\vc\__init__.py') -Destination $target -Force
foreach ($name in @('surface_index', 'volume')) {
    $module = Get-ChildItem -LiteralPath $moduleDir -Filter "$name*.pyd" |
        Select-Object -First 1
    if (-not $module) {
        throw "missing built Villa module $name in $moduleDir"
    }
    Copy-Item -LiteralPath $module.FullName -Destination $target -Force
}
foreach ($dir in @($runtimeDir, $moduleDir)) {
    Get-ChildItem -LiteralPath $dir -Filter '*.dll' | ForEach-Object {
        Copy-Item -LiteralPath $_.FullName -Destination $target -Force
    }
}

& $python -c 'import vc; from vc import surface_index; print("VILLA_PYTHON_BINDINGS OK:", surface_index.__file__)'
if ($LASTEXITCODE -ne 0) {
    throw "Villa Python binding import test failed with exit code $LASTEXITCODE"
}
Write-Output "Villa spiral environment ready: $EnvironmentDir"
