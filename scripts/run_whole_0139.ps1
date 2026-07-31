# run_whole_0139.ps1 -- carve + mesh the WHOLE PHerc0139 scroll in the background.
#
# Phase 1: carve cubes_PRED for the whole scroll from the local m7 zarr mirror
#          (PRED only -- meshing does not need RAW; that is carved per-slab later).
# Phase 2: grid_pipeline per-cube mesh (NO weld -- we want per-cube meshes for the
#          whole-scroll unroll), capped cores, with a healing loop that re-runs
#          --skip-existing until the done-count stops rising (robust to killed
#          workers). Everything tee'd to output/rebuild_0139_full/run.log.
#
# CPU: -MaxConcurrent cubes x -Threads per cube = total threads. Default 8x2 = 16
#      threads, leaving the rest of the 64c/128t box free for other work.
# Resume: re-run with -SkipCarve after an interruption; --skip-existing only
#         redoes incomplete cubes.
param(
    [int]$MaxConcurrent = 8,
    [int]$Threads = 2,
    [switch]$SkipCarve
)
$ErrorActionPreference = 'Continue'
$repo = "D:\work\vesuvius-c"
Set-Location $repo
$zarr = "$repo\PHerc0139-full\20250728140407-surface-20260413222639-surface-m7-L0-th0.2.zarr"
$grid = "$repo\PHerc0139-full-grid"                 # cubes_PRED lands here (gitignored)
$meshOut = "$repo\output\rebuild_0139_full"
New-Item -ItemType Directory -Force $meshOut | Out-Null
$log = "$meshOut\run.log"
function Log($m) { $s = "[$(Get-Date -Format 'MM-dd HH:mm:ss')] $m"; Write-Output $s; Add-Content $log $s }

Log "=== WHOLE PHerc0139 carve+mesh START  (max-concurrent $MaxConcurrent, threads/cube $Threads) ==="
$t0 = Get-Date

# ---- Phase 1: carve cubes_PRED (whole scroll, PRED only) ----
if (-not $SkipCarve -and -not (Test-Path "$grid\cubes_PRED\present.json")) {
    Log "Phase 1: carving cubes_PRED (bbox 0 20864 0 6528 0 6528) -> $grid"
    py -m uv run --project python python python/scripts/carve_grid_tifs.py `
        --pred-zarr $zarr --bbox 0 20864 0 6528 0 6528 `
        --umbilicus 3405 2878 --out $grid *>> $log
    Log "Phase 1: carve exit $LASTEXITCODE"
} else {
    Log "Phase 1: carve skipped (present.json exists or -SkipCarve)"
}
$nPred = @(Get-ChildItem "$grid\cubes_PRED\*.tif" -ErrorAction SilentlyContinue).Count
Log "cubes_PRED present: $nPred"
if ($nPred -lt 1) { Log "ABORT: no PRED cubes carved"; exit 1 }

# ---- Phase 2: per-cube mesh, no weld, healing loop ----
$gp = "$repo\build\Release\grid_pipeline.exe"
$cube = "$repo\build\Release\cube_mesh.exe"
$weld = "$repo\build\Release\grid_weld.exe"
$env:VESUVIUS_THREADS = "$Threads"; $env:OMP_NUM_THREADS = "$Threads"
function Count-Done {
    @(Get-ChildItem "$meshOut\dump" -Directory -ErrorAction SilentlyContinue | Where-Object {
        Test-Path "$($_.FullName)\$($_.Name)_step12_final\$($_.Name)_step12_final_all.obj"
    }).Count
}
$prev = -1
for ($heal = 0; $heal -lt 6; $heal++) {
    $done = Count-Done
    Log "Phase 2 heal iter $heal : $done / ~$nPred cubes done (prev $prev)"
    if ($heal -gt 0 -and $done -eq $prev) { Log "Phase 2 converged"; break }
    $prev = $done
    # NOTE (2026-07-21): --qem only sets skip_qem=0 since the CVT default flip
    # (51bcf5a); it does NOT select the engine. The ~13.4k cubes meshed before
    # the flip are QEM; use --simplify qem explicitly so a resume stays a
    # single-engine scroll. (Switching this run to CVT means re-meshing from
    # scratch into a fresh out dir -- decide deliberately, don't mix.)
    & $gp $grid $meshOut --halo 13 --threads-per-cube $Threads --max-concurrent $MaxConcurrent `
        --exe $cube --weld $weld --simplify qem --trim-inset 0 --skip-weld --skip-existing *>> $log
    Log "Phase 2 grid_pipeline pass $heal exit $LASTEXITCODE"
}

$done = Count-Done
$mins = [int]((Get-Date) - $t0).TotalMinutes
Log "=== DONE: $done / ~$nPred cubes meshed in $mins min -> $meshOut\dump ==="
Log "Next: global scroll_whole registration, then per-z-slab scroll_unroll --z-range + tifxyz_weld."
