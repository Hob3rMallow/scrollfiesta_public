<#
  run_hierarchical.ps1 -- build + run the hierarchical LOD-pyramid orchestrator.

  Welds a grid of already-meshed leaf cubes up a fanout^3 pyramid, decimating a
  lot between levels (see src/tools/hierarchical_weld.c). Emits one LOD tier per
  level under <OutDir>/level<N>/, with per-level _summary.csv + _manifold.log.

  Prereqs: leaf cubes already meshed under <LeafDump>/<cube>/<cube>_step12_final/.
  Produce them with:
    build/Release/grid_pipeline.exe <grid> <OutDir> --halo 13 --skip-weld --skip-existing

  Usage:
    pwsh scripts/run_hierarchical.ps1                     # 4x21x21 defaults
    pwsh scripts/run_hierarchical.ps1 -KeepRatio 0.125 -MaxConcurrent 64
    pwsh scripts/run_hierarchical.ps1 -DryRun
#>
[CmdletBinding()]
param(
  [string]$LeafDump      = "output/grid_4x21x21_lod/dump",
  [string]$OutDir        = "output/grid_4x21x21_lod",
  [string]$LeafStage     = "step12_final",
  [double]$KeepRatio     = 0.25,
  [int]$Fanout           = 2,
  [int]$MaxConcurrent    = 48,
  [int]$StopNodes        = 1,
  [int]$MaxLevels        = 16,
  [double]$BlockTimeout  = 900,
  [double]$UmbY          = 3405,
  [double]$UmbX          = 2878,
  [double]$WrapPitch     = 9.5,
  [switch]$SkipBuild,
  [switch]$SkipExisting,
  [switch]$SimplifyTop,
  [switch]$DryRun
)
$ErrorActionPreference = "Stop"
Set-Location (Split-Path $PSScriptRoot -Parent)   # repo root

$exe     = "build/Release/hierarchical_weld.exe"
$msbuild = "C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe"

if (-not $SkipBuild) {
  Write-Host "== building hierarchical_weld.vcxproj (Release x64) =="
  & $msbuild hierarchical_weld.vcxproj /p:Configuration=Release /p:Platform=x64 /m /v:minimal
  if ($LASTEXITCODE -ne 0) { throw "build failed" }
}
foreach ($dep in "build/Release/grid_weld.exe","build/Release/qslim_obj.exe","build/Release/manifold_check.exe") {
  if (-not (Test-Path $dep)) { Write-Warning "$dep missing -- build it (msbuild <name>.vcxproj)" }
}
if (-not (Test-Path $LeafDump)) { throw "leaf dump not found: $LeafDump" }

$args = @(
  $LeafDump, $OutDir,
  "--leaf-stage", $LeafStage,
  "--fanout", $Fanout,
  "--keep-ratio", $KeepRatio,
  "--stop-nodes", $StopNodes,
  "--max-levels", $MaxLevels,
  "--max-concurrent", $MaxConcurrent,
  "--block-timeout", $BlockTimeout,
  "--umbilicus-y", $UmbY, "--umbilicus-x", $UmbX, "--wrap-pitch", $WrapPitch
)
if ($SkipExisting) { $args += "--skip-existing" }
if ($SimplifyTop)  { $args += "--simplify-top" }
if ($DryRun)       { $args += "--dry-run" }

New-Item -ItemType Directory -Force $OutDir | Out-Null
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$log   = Join-Path $OutDir "hierarchical_$stamp.log"
Write-Host "== hierarchical_weld $($args -join ' ') =="
Write-Host "== logging to $log =="

$sw = [System.Diagnostics.Stopwatch]::StartNew()
& $exe @args 2>&1 | Tee-Object -FilePath $log
$rc = $LASTEXITCODE
$sw.Stop()

if ($DryRun) { exit $rc }

# ---- end-of-run summary: one row per level from its _summary.csv ----
Write-Host ""
Write-Host "===== PYRAMID SUMMARY ($([int]$sw.Elapsed.TotalMinutes) min, tool exit $rc) ====="
"{0,-7} {1,6} {2,5} {3,5} {4,5} {5,14} {6,10}" -f "level","nodes","ok","warn","fail","total_faces","worst_nm" | Write-Host
Get-ChildItem "$OutDir/level*" -Directory -ErrorAction SilentlyContinue |
  Sort-Object { [int]($_.Name -replace 'level','') } |
  ForEach-Object {
    $csv = Join-Path $_.FullName "_summary.csv"
    if (Test-Path $csv) {
      $rows = Import-Csv $csv
      $ok   = ($rows | Where-Object status -eq 'ok').Count
      $warn = ($rows | Where-Object status -eq 'warn').Count
      $fail = ($rows | Where-Object status -eq 'fail').Count
      $tf   = ($rows | Measure-Object final_faces -Sum).Sum
      $wnm  = ($rows | Measure-Object non_manifold -Maximum).Maximum
      "{0,-7} {1,6} {2,5} {3,5} {4,5} {5,14} {6,10}" -f $_.Name,$rows.Count,$ok,$warn,$fail,$tf,$wnm | Write-Host
    }
  }
$top = Get-ChildItem "$OutDir/level*" -Directory -ErrorAction SilentlyContinue | Sort-Object { [int]($_.Name -replace 'level','') } | Select-Object -Last 1
if ($top) { Write-Host "Top LOD tier: $($top.FullName)" }
exit $rc
