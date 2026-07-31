<#
  run_pred_reject_scan.ps1 -- build + run the garbage-prediction-cube detector
  over a grid's cubes_PRED, plus the labelled calibration examples.

  Runner only (the analysis is the C tool src/tools/pred_reject.c). Writes the
  reject list + per-cube stats CSV under output/tools/pred_reject/ and prints
  the verdicts for the known-garbage and known-real reference cubes.

  Usage:
    pwsh scripts/extract/run_pred_reject_scan.ps1 [-Grid PHerc0139-4x21x21-cubes] [-NoBuild]
#>
param(
  [string]$Grid = "PHerc0139-4x21x21-cubes",
  [switch]$NoBuild
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)   # repo root
Set-Location $root

$msbuild = "C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe"
$exe = "build\Release\pred_reject.exe"
if (-not $NoBuild) {
  & $msbuild pred_reject.vcxproj /p:Configuration=Release /p:Platform=x64 /m /v:minimal
  if ($LASTEXITCODE -ne 0) { throw "build failed" }
}

$outDir = "output\tools\pred_reject"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null
$rejects = Join-Path $outDir "rejects.txt"
$stats   = Join-Path $outDir "stats.csv"
$pred    = Join-Path $Grid "cubes_PRED"

Write-Host "== scanning $pred ==" -ForegroundColor Cyan
& $exe --scan $pred --out $rejects --csv $stats

Write-Host "`n== labelled reference cubes ==" -ForegroundColor Cyan
# 4 known-garbage rectangles + the dense slab originally mistaken for 'real'.
$refs = @(
  "z04352_y02048_x01536","z04352_y02176_x01536","z04608_y04608_x03200",
  "z04608_y04608_x04096","z04736_y04608_x03968"
)
foreach ($c in $refs) {
  $tif = Join-Path $pred "$c.tif"
  if (Test-Path $tif) { & $exe $tif }
}

$nRej = (Get-Content $rejects | Measure-Object -Line).Lines
Write-Host "`nrejects: $nRej  (list: $rejects, stats: $stats)" -ForegroundColor Green
