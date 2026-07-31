[CmdletBinding()]
param(
    [string]$Version = "13.2.86"
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
$target = Join-Path $repo ".toolchains\cuda-13.2-wheels"
$uv = Get-Command uv -ErrorAction SilentlyContinue
if (-not $uv) {
    throw "uv not found; install uv or place a supported CUDA toolkit in CUDA_PATH"
}

New-Item -ItemType Directory -Force -Path $target | Out-Null
$packages = @(
    "nvidia-cuda-runtime==$Version",
    "nvidia-cuda-nvrtc==$Version",
    "nvidia-cuda-nvcc==$Version",
    "nvidia-cuda-cccl==$Version"
)
& $uv.Source pip install --target $target --upgrade @packages
if ($LASTEXITCODE -ne 0) {
    throw "CUDA wheel installation failed with exit code $LASTEXITCODE"
}

$cudaRoot = Join-Path $target "nvidia\cu13"
$required = @(
    "bin\nvcc.exe",
    "bin\x86_64\nvrtc64_130_0.dll",
    "include\cuda.h",
    "include\cccl"
)
foreach ($relative in $required) {
    if (-not (Test-Path -LiteralPath (Join-Path $cudaRoot $relative))) {
        throw "CUDA payload is incomplete: missing $relative under $cudaRoot"
    }
}

Write-Host "Project-local CUDA build/runtime payload installed at $cudaRoot"
Write-Host "scripts\build_mls_cubecl.ps1 will discover it automatically."
