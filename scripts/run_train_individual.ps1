#requires -Version 7.0
<#
.SYNOPSIS
  Run a bunch of Kaggle training cubes (nnunet-preds/*.tif) through the per-cube
  pipeline INDIVIDUALLY -- no cross-cube welding. Cubes run CONCURRENTLY (each is
  serial-bound at ~1 core, so running them one-at-a-time wastes the box); threads
  per cube are capped so the parallel regions don't oversubscribe. Full per-cube
  logging + a summary table at the end.

.DESCRIPTION
  Each cube is an isolated 128^3 nnUNet prediction with an arbitrary integer ID
  (NOT the z#####_y#####_x##### grid naming), so --halo 0 (no neighbour load) is
  the right call -- there are no spatial neighbours to load from nnunet-preds.

  Concurrency model (mirrors grid_pipeline's tpc=1 + many-concurrent pattern,
  scaled for a handful of cubes): up to -MaxConcurrent cubes in flight at once,
  each pinned to VESUVIUS_THREADS=-Threads. Default 7 cubes x 4 threads = 28 <= 32
  logical cores. Most of the per-cube pipeline (BPA front growth, Hoppe orient,
  hole fill) is serial, so the real win is cube-level concurrency, not threads.

  This is a RUNNER only (PS as a launcher); any mesh analysis belongs in a C tool.

  Output lands under output/<Experiment>/:
    <id>/<id>_<stage>/...   per-cube OBJ dumps (cube id is in the path)
    logs/<id>.log / .err    per-cube stdout / stderr
    summary.txt             the end-of-run summary table

.PARAMETER Cubes         Cube IDs to run. Default = the 7 canonical eval cubes.
.PARAMETER Experiment    Output subdir under output/ (default: train_individual).
.PARAMETER PredDir       Dir holding <id>.tif predictions (default: nnunet-preds).
.PARAMETER Exe           cube_mesh.exe path (default build/Release/cube_mesh.exe).
.PARAMETER Halo          Halo voxels (default 0 = isolated cube).
.PARAMETER Threads       OpenMP threads per cube (VESUVIUS_THREADS). Default 4.
.PARAMETER MaxConcurrent Max cubes in flight at once. Default 8.
.PARAMETER NoQem         Disable QEM decimation (default: QEM on).
.PARAMETER Force         Overwrite an existing output/<Experiment> dir.
#>
[CmdletBinding()]
param(
    [string[]]$Cubes = @('86701140','118632705','1006462223','1013184726',
                         '3290306825','3294954456','3394433588'),
    [string]  $Experiment    = 'train_individual',
    [string]  $PredDir       = 'nnunet-preds',
    [string]  $Exe           = 'build/Release/cube_mesh.exe',
    [int]     $Halo          = 0,
    [int]     $Threads       = 4,
    [int]     $MaxConcurrent = 8,
    [switch]  $NoQem,
    [switch]  $Force
)

$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path "$PSScriptRoot/..").Path
Set-Location $repoRoot

if (-not (Test-Path $Exe)) { throw "Missing exe: $Exe (build it first)" }
$exeAbs = (Resolve-Path $Exe).Path

$outDir = Join-Path $repoRoot "output/$Experiment"
if (Test-Path $outDir) {
    if (-not $Force) { throw "output/$Experiment already exists. Use -Force or pick a new -Experiment." }
    Remove-Item -Recurse -Force $outDir
}
$logDir = Join-Path $outDir 'logs'
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$env:VESUVIUS_THREADS = "$Threads"
Write-Host ("==== run_train_individual: {0} cube(s), halo={1}, threads/cube={2}, maxconc={3}, qem={4} ====" -f `
    $Cubes.Count, $Halo, $Threads, $MaxConcurrent, (-not $NoQem)) -ForegroundColor Cyan

$pending = [System.Collections.Generic.Queue[string]]::new()
$Cubes | ForEach-Object { $pending.Enqueue($_) }
$running = @{}        # procId -> tracking object
$results = @()
$started = Get-Date

while ($pending.Count -gt 0 -or $running.Count -gt 0) {
    # fill the pool
    while ($running.Count -lt $MaxConcurrent -and $pending.Count -gt 0) {
        $id  = $pending.Dequeue()
        $tif = Join-Path $PredDir "$id.tif"
        $log = Join-Path $logDir "$id.log"
        $err = Join-Path $logDir "$id.err"
        if (-not (Test-Path $tif)) {
            Write-Host "  SKIP $id : no $tif" -ForegroundColor Yellow
            $results += [pscustomobject]@{ cube=$id; exit='no-tif'; sec=0; comps=''; log=$log }
            continue
        }
        $argList = @($tif, (Join-Path $outDir "$id.tif"), '--dump-obj', $outDir, '--halo', "$Halo")
        if ($NoQem) { $argList += '--no-qem' }
        $proc = Start-Process -FilePath $exeAbs -ArgumentList $argList `
                    -RedirectStandardOutput $log -RedirectStandardError $err `
                    -NoNewWindow -PassThru
        $running[$proc.Id] = [pscustomobject]@{ cube=$id; proc=$proc; t0=(Get-Date); log=$log; err=$err }
        Write-Host ("  launched {0} (pid {1})  [{2} running]" -f $id, $proc.Id, $running.Count)
    }
    # reap finished
    Start-Sleep -Milliseconds 500
    foreach ($procId in @($running.Keys)) {
        $r = $running[$procId]
        if ($r.proc.HasExited) {
            $sec  = [math]::Round(((Get-Date) - $r.t0).TotalSeconds, 1)
            $code = $r.proc.ExitCode
            $comps = ''
            foreach ($f in @($r.log, $r.err)) {
                if (-not (Test-Path $f)) { continue }
                $m = Select-String -Path $f -Pattern '(\d+)\s+component' -EA SilentlyContinue | Select-Object -Last 1
                if (-not $m) { $m = Select-String -Path $f -Pattern 'comp(?:onent)?s?\s*[:=]\s*(\d+)' -EA SilentlyContinue | Select-Object -Last 1 }
                if ($m) { $comps = $m.Matches[0].Groups[1].Value }
            }
            $color = if ($code -eq 0) { 'Green' } else { 'Red' }
            Write-Host ("  done {0} (pid {1}): exit={2} {3:n1}s" -f $r.cube, $procId, $code, $sec) -ForegroundColor $color
            $results += [pscustomobject]@{ cube=$r.cube; exit=$code; sec=$sec; comps=$comps; log=$r.log }
            $running.Remove($procId)
        }
    }
}

$wall = [math]::Round(((Get-Date) - $started).TotalSeconds, 1)

# ---- SUMMARY ----
$results = $results | Sort-Object { [array]::IndexOf($Cubes, $_.cube) }
$summaryPath = Join-Path $outDir 'summary.txt'
$ok     = ($results | Where-Object { $_.exit -eq 0 }).Count
$tot    = $results.Count
$cpuSec = [math]::Round((($results | Measure-Object -Property sec -Sum).Sum), 1)
$tbl = $results | Format-Table cube, exit, sec, comps -AutoSize | Out-String
$summary = @"
==== train_individual summary ($Experiment) ====
cubes: $tot   ok: $ok   failed: $($tot-$ok)
wall: ${wall}s   sum-of-per-cube: ${cpuSec}s   (concurrency win = sum/wall)
halo=$Halo  threads/cube=$Threads  maxconc=$MaxConcurrent  qem=$(-not $NoQem)  exe=$Exe
$tbl
artifacts: $outDir  (per-cube dumps under <id>/, logs under logs/)
"@
$summary | Tee-Object -FilePath $summaryPath
if ($ok -ne $tot) { exit 1 }
exit 0
