#requires -Version 7.0
<#
.SYNOPSIS
  Thin wrapper around grid_pipeline.exe.

.DESCRIPTION
  After the May-2026 cleanup the per-cube + weld orchestration lives in
  src/tools/grid_pipeline.c. This script preserves the same invocation
  for legacy callers but delegates to the C binary.

.PARAMETER GridDir          Grid directory containing cubes_PRED/*.tif.
.PARAMETER OutputDir        Output directory.
.PARAMETER Halo             Halo voxels (default 8).
.PARAMETER MaxCubes         Stop after N cubes (default all).
.PARAMETER ThreadsPerCube   OpenMP threads per cube (default 4).
.PARAMETER MaxConcurrent    Max concurrent subprocesses (default cores/tpc).
.PARAMETER Qem              Enable QEM decimation (default).
.PARAMETER SkipWeld         Skip running grid_weld at the end.
.PARAMETER CheckDeterminism Run each cube twice and cmp OBJs.
#>

[CmdletBinding()]
param(
    [string]$GridDir          = 'PHerc0139-4x5x5',
    [string]$OutputDir        = 'output/halo_grid',
    [int]   $Halo             = 13,
    [int]   $MaxCubes         = 0,
    [int]   $ThreadsPerCube   = 4,
    [int]   $MaxConcurrent    = 0,
    [string]$ExePath,
    [string]$WeldExePath,
    [string]$GridPipelinePath,
    [switch]$Qem,
    [switch]$SkipWeld,
    [switch]$CheckDeterminism
)

$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path "$PSScriptRoot/..").Path
if (-not $GridPipelinePath) {
    $GridPipelinePath = Join-Path $repoRoot 'build/Release/grid_pipeline.exe'
}
if (-not $ExePath)     { $ExePath     = Join-Path $repoRoot 'build/Release/cube_mesh.exe' }
if (-not $WeldExePath) { $WeldExePath = Join-Path $repoRoot 'build/Release/grid_weld.exe' }

foreach ($p in $GridPipelinePath, $ExePath, $WeldExePath) {
    if (-not (Test-Path $p)) { throw "Missing path: $p" }
}

$gridDir   = (Resolve-Path (Join-Path $repoRoot $GridDir)).Path
$outputDir = Join-Path $repoRoot $OutputDir
New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

$argList = @(
    $gridDir, $outputDir,
    '--halo', $Halo,
    '--threads-per-cube', $ThreadsPerCube,
    '--exe', $ExePath,
    '--weld', $WeldExePath
)
if ($MaxCubes -gt 0)       { $argList += @('--max-cubes', $MaxCubes) }
if ($MaxConcurrent -gt 0)  { $argList += @('--max-concurrent', $MaxConcurrent) }
if (-not $Qem)             { $argList += '--no-qem' }
if ($SkipWeld)             { $argList += '--skip-weld' }
if ($CheckDeterminism)     { $argList += '--check-determinism' }

Write-Host "grid_pipeline:" ($argList -join ' ')
& $GridPipelinePath @argList
exit $LASTEXITCODE
