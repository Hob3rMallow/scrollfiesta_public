[CmdletBinding()]
param(
    [ValidateSet("cuda", "hip", "cpu", "rust-cpu")]
    [string]$Backend = "cuda",
    [switch]$SkipTests
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$crate = Join-Path $repo "deps\scrollfiesta-mls-cubecl"
$cargoHome = Join-Path $repo ".toolchains\cargo"
$rustupHome = Join-Path $repo ".toolchains\rustup"
$cargo = Join-Path $cargoHome "bin\cargo.exe"
if (-not (Test-Path -LiteralPath $cargo)) {
    $cargoCommand = Get-Command cargo -ErrorAction SilentlyContinue
    if (-not $cargoCommand) {
        throw "cargo not found; run scripts\install_rust_local.ps1 first"
    }
    $cargo = $cargoCommand.Source
} else {
    $env:CARGO_HOME = $cargoHome
    $env:RUSTUP_HOME = $rustupHome
}

$features = switch ($Backend) {
    "cuda" { "rust-cpu,cubecl-cuda" }
    "hip" { "rust-cpu,cubecl-hip" }
    "cpu" { "rust-cpu,cubecl-cpu" }
    default { "rust-cpu" }
}

if ($Backend -eq "cuda") {
    $localCuda = Join-Path $repo ".toolchains\cuda-13.2-wheels\nvidia\cu13"
    if (-not $env:CUDA_PATH) {
        if (Test-Path -LiteralPath (Join-Path $localCuda "include\cccl")) {
            $env:CUDA_PATH = $localCuda
            Write-Host "Using project-local CUDA payload at $localCuda"
        } else {
            Write-Warning "CUDA_PATH is unset. Run scripts\install_cuda_local.ps1 before the CUDA smoke test."
        }
    }
    if ($env:CUDA_PATH) {
        if (-not (Test-Path -LiteralPath (Join-Path $env:CUDA_PATH "include\cccl"))) {
            throw "CUDA_PATH does not contain the CCCL headers expected by CubeCL: $env:CUDA_PATH"
        }
        $cudaBin = Join-Path $env:CUDA_PATH "bin"
        $cudaDllBin = Join-Path $cudaBin "x86_64"
        $env:PATH = "$cudaDllBin;$cudaBin;$env:PATH"
    }
}

Push-Location $crate
try {
    & $cargo fmt --all -- --check
    if ($LASTEXITCODE -ne 0) { throw "cargo fmt check failed" }
    & $cargo build --release --locked --no-default-features --features $features
    if ($LASTEXITCODE -ne 0) { throw "cargo build failed" }
    if (-not $SkipTests) {
        & $cargo test --release --locked --no-default-features --features $features
        if ($LASTEXITCODE -ne 0) { throw "cargo test failed" }
    }
} finally {
    Pop-Location
}

$target = Join-Path $crate "target\release"
$dll = Join-Path $target "herculaneum_mls_cubecl.dll"
$import = Join-Path $target "herculaneum_mls_cubecl.dll.lib"
if (-not (Test-Path -LiteralPath $dll)) { throw "Cargo did not produce $dll" }
if (-not (Test-Path -LiteralPath $import)) {
    throw "Cargo did not produce the expected MSVC import library $import"
}

$dest = Join-Path $repo "deps\lib\win64\Release"
New-Item -ItemType Directory -Force -Path $dest | Out-Null
# The import library embeds this original DLL basename. Keep it unchanged;
# renaming only the file causes Windows loader error 0xC0000135.
Copy-Item -Force -LiteralPath $dll -Destination (Join-Path $dest "herculaneum_mls_cubecl.dll")
Copy-Item -Force -LiteralPath $import -Destination (Join-Path $dest "mls_cubecl_import.lib")

$static = Join-Path $target "herculaneum_mls_cubecl.lib"
if (Test-Path -LiteralPath $static) {
    Copy-Item -Force -LiteralPath $static -Destination (Join-Path $dest "mls_cubecl_static.lib")
}
if ($Backend -eq "cuda" -and $env:CUDA_PATH) {
    # NVRTC is loaded dynamically by cudarc, so PATH alone is too fragile for
    # Python/Start-Process children on Windows (which can carry both Path and
    # PATH entries). Stage the official runtime DLLs beside cube_mesh instead.
    $cudaRuntimeSource = Join-Path $env:CUDA_PATH "bin\x86_64"
    $cudaRuntimeDest = Join-Path $dest "cubecl-cuda-runtime"
    $cudaRuntimeDlls = Get-ChildItem -LiteralPath $cudaRuntimeSource -Filter "*.dll" -File
    if (-not ($cudaRuntimeDlls.Name -like "nvrtc*.dll")) {
        throw "CUDA runtime payload has no NVRTC DLL under $cudaRuntimeSource"
    }
    New-Item -ItemType Directory -Force -Path $cudaRuntimeDest | Out-Null
    $cudaRuntimeDlls | Copy-Item -Force -Destination $cudaRuntimeDest
    Write-Host "Staged $($cudaRuntimeDlls.Count) CUDA runtime DLLs for child processes"
}
Write-Host "CubeCL MLS $Backend artifacts staged in $dest"
