#requires -Version 7.0
<#
.SYNOPSIS
  One-stop shop: rebuild the pipeline binaries, regenerate the full per-cube
  grid from scratch, seam-weld it, and summarize the result.

.DESCRIPTION
  This is the canonical "rebuild + regenerate + test on real cube data" entry
  point. It does NOT reimplement the grid orchestration -- the per-cube run and
  the seam weld already live in grid_pipeline.exe (also wrapped by the legacy
  scripts/run_grid_halo.ps1). This script just adds the two things that loop was
  missing for a post-major-change workflow: a from-scratch BUILD, and a VERIFY/
  summary pass (with a world-space guard for the stacked-cube bug). The actual
  run is delegated to the same grid_pipeline.exe.

  Stages:
    1. BUILD   cube_mesh.exe, grid_weld.exe, grid_pipeline.exe (Release x64).
    2. RUN     grid_pipeline.exe (the existing engine): per-cube extract+dump in
               world coords, then grid_weld --seam-weld (peel + BPA bridge).
    3. VERIFY  read welded.obj.weld_report.json, compute the welded bbox, and
               HARD-FAIL if the cubes did not land in world space (the
               "all cubes stacked at [0,cube]" bug) -- the exact failure that
               motivated this script.

  All output lands under output/<Experiment>/ (one directory per experiment).

.PARAMETER Experiment   Output subdir under output/ (default: grid_seamweld).
.PARAMETER Grid         Grid dir with cubes_PRED/*.tif (default PHerc0139-4x5x5).
.PARAMETER Cubes        Comma-separated cube-id tokens to RESTRICT the run to a
                        subset (e.g. a 2x1x1 seam test). Each token is matched as
                        a substring against the TIFF basenames; matching TIFFs are
                        linked into a filtered cubes_PRED and only those are run.
                        Empty = whole grid. Origins come from the filename, so no
                        manifest is needed for the subset.
.PARAMETER Halo         Halo voxels (default 13).
.PARAMETER MaxCubes     Stop after N cubes (default 0 = all). Quick smoke; note
                        FS order is not guaranteed adjacent -- prefer -Cubes for a
                        spatially meaningful subset.
.PARAMETER NoQem        Disable QEM decimation (QEM is on by default).
.PARAMETER Clean        Clean each project before building.
.PARAMETER SkipBuild    Skip the rebuild; run with whatever is in build/Release.
.PARAMETER Force        Overwrite an existing output/<Experiment> dir.

.EXAMPLE
  pwsh scripts/run_pipeline_grid.ps1                       # full 100-cube grid
.EXAMPLE
  pwsh scripts/run_pipeline_grid.ps1 -Experiment seam_2x1x1 `
       -Cubes z04352_y03072_x02560,z04352_y03072_x02688   # one x-seam, fast
#>
[CmdletBinding()]
param(
    [string]  $Experiment = 'grid_seamweld',
    [string]  $Grid       = 'PHerc0139-4x5x5',
    [string[]]$Cubes      = @(),
    [int]     $Halo       = 13,
    [int]   $MaxCubes   = 0,
    [switch]$NoQem,
    [switch]$Clean,
    [switch]$SkipBuild,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path "$PSScriptRoot/..").Path
Set-Location $repoRoot   # keep relative dependency and output paths stable
$msbuild  = 'C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe'
$cubeSize = 128.0   # source-voxel cube edge; seam planes sit at multiples of this

function Section($t) { Write-Host ""; Write-Host "==== $t ====" -ForegroundColor Cyan }

# ---- 1. BUILD ------------------------------------------------------------
$projects = 'cube_mesh.vcxproj', 'grid_weld.vcxproj', 'grid_pipeline.vcxproj'
if (-not $SkipBuild) {
    Section "BUILD (Release x64)"
    foreach ($proj in $projects) {
        $path = Join-Path $repoRoot $proj
        if ($Clean) {
            & $msbuild $path /t:Clean /p:Configuration=Release /p:Platform=x64 /v:quiet | Out-Null
        }
        Write-Host "  building $proj ..."
        & $msbuild $path /p:Configuration=Release /p:Platform=x64 /m /v:minimal /nologo
        if ($LASTEXITCODE -ne 0) { throw "BUILD FAILED: $proj" }
    }
    Write-Host "  all binaries built OK" -ForegroundColor Green
} else {
    Section "BUILD skipped (-SkipBuild)"
}

$exe   = Join-Path $repoRoot 'build/Release/cube_mesh.exe'
$weld  = Join-Path $repoRoot 'build/Release/grid_weld.exe'
$gpipe = Join-Path $repoRoot 'build/Release/grid_pipeline.exe'
foreach ($p in $exe, $weld, $gpipe) { if (-not (Test-Path $p)) { throw "Missing binary: $p" } }

# ---- 2. RUN --------------------------------------------------------------
$fullGrid = (Resolve-Path (Join-Path $repoRoot $Grid)).Path
$outDir   = Join-Path $repoRoot "output/$Experiment"
if (Test-Path $outDir) {
    if (-not $Force) { throw "output/$Experiment already exists. Use -Force to overwrite, or pick a new -Experiment." }
    Write-Host "  -Force: removing existing $outDir"
    Remove-Item -Recurse -Force $outDir
}
New-Item -ItemType Directory -Force -Path $outDir | Out-Null
$pipeLog = Join-Path $outDir 'pipeline.log'

# Subset selection: when -Cubes is given, build a filtered cubes_PRED holding
# only the matching TIFFs and point grid_pipeline at it. grid_pipeline only reads
# cubes_PRED/*.tif and derives each origin from the filename, so a subset dir with
# just the TIFFs is sufficient (no manifest needed).
$gridDir = $fullGrid
$nSelected = 0
$tokens = @($Cubes | ForEach-Object { $_ -split ',' } | ForEach-Object { $_.Trim() } | Where-Object { $_ })
if ($tokens.Count -gt 0) {
    $srcPred = Join-Path $fullGrid 'cubes_PRED'
    $subPred = Join-Path $outDir '_griddir/cubes_PRED'
    New-Item -ItemType Directory -Force -Path $subPred | Out-Null
    foreach ($tif in Get-ChildItem (Join-Path $srcPred '*.tif')) {
        if ($tokens | Where-Object { $tif.Name -like "*$_*" }) {
            Copy-Item $tif.FullName (Join-Path $subPred $tif.Name)
            $nSelected++
        }
    }
    if ($nSelected -eq 0) { throw "-Cubes '$($tokens -join ',')' matched no TIFFs in $srcPred" }
    $gridDir = Join-Path $outDir '_griddir'
    Write-Host "  subset: $nSelected cube(s) matched -Cubes '$($tokens -join ',')'" -ForegroundColor Yellow
}

Section "RUN grid_pipeline -> output/$Experiment"
$gpArgs = @($gridDir, $outDir, '--halo', $Halo, '--exe', $exe, '--weld', $weld)
if ($MaxCubes -gt 0) { $gpArgs += @('--max-cubes', $MaxCubes) }
if ($NoQem)          { $gpArgs += '--no-qem' } else { $gpArgs += '--qem' }
Write-Host "  grid_pipeline $($gpArgs -join ' ')"
$started = Get-Date
& $gpipe @gpArgs 2>&1 | Tee-Object -FilePath $pipeLog
$gpExit = $LASTEXITCODE
$elapsed = (Get-Date) - $started

# ---- 3. VERIFY + SUMMARIZE ----------------------------------------------
Section "SUMMARY"
$weldObj    = Join-Path $outDir 'welded.obj'
$reportPath = "$weldObj.weld_report.json"
$weldLog    = Join-Path $outDir 'grid_weld.log'
$fail = $false
$cubesProcessed = 0

Write-Host ("  grid_pipeline exit: {0}   elapsed: {1:n1} min" -f $gpExit, $elapsed.TotalMinutes)

if (Test-Path $reportPath) {
    $r = Get-Content $reportPath -Raw | ConvertFrom-Json
    $a = $r.manifold_audit
    $cubesProcessed = [int]$r.cubes_processed
    Write-Host ("  cubes processed : {0}" -f $r.cubes_processed)
    Write-Host ("  unique verts    : {0:n0}" -f $r.total_unique_verts)
    Write-Host ("  unique faces    : {0:n0}" -f $r.total_unique_faces)
    Write-Host  "  manifold audit  :"
    Write-Host ("     manifold pairs : {0:n0}" -f $a.manifold_pairs)
    Write-Host ("     unpaired edges : {0:n0}  (open boundary; seam-weld should drive this DOWN)" -f $a.unpaired)
    Write-Host ("     non-manifold   : {0}" -f $a.non_manifold) -ForegroundColor ($(if ($a.non_manifold -gt 0) {'Yellow'} else {'Gray'}))
    Write-Host ("     same-dir pairs : {0}  (foldovers; must be 0)" -f $a.same_dir_pairs) -ForegroundColor ($(if ($a.same_dir_pairs -gt 0) {'Red'} else {'Gray'}))
    if ($a.same_dir_pairs -gt 0) { $fail = $true }
} else {
    Write-Host "  NO weld report found ($reportPath)" -ForegroundColor Red
    $fail = $true
}

# Seam-weld bridge stats from the grid_weld log.
if (Test-Path $weldLog) {
    $seam = Select-String -Path $weldLog -Pattern '\[seam\] planes=|\[seam\] merged:|Seam-weld:|bridge:'
    if ($seam) {
        Write-Host "  seam-weld:"
        $seam | ForEach-Object { Write-Host ("     " + $_.Line.Trim()) }
    } else {
        Write-Host "  seam-weld: NO [seam] activity logged (detect_planes found 0 planes?)" -ForegroundColor Yellow
    }
}

# World-space sanity check: stream the welded.obj bbox. The stacked-cube bug
# (grid_weld not applying cube origins) collapses every cube into [0,cubeSize];
# a correct multi-cube weld must span several cube widths on at least one axis.
if (Test-Path $weldObj) {
    $mnx = [double]::PositiveInfinity; $mxx = [double]::NegativeInfinity
    $mny = [double]::PositiveInfinity; $mxy = [double]::NegativeInfinity
    foreach ($line in [System.IO.File]::ReadLines($weldObj)) {
        if ($line.Length -lt 2 -or $line[0] -ne 'v' -or $line[1] -ne ' ') { continue }
        $p = $line.Split(' ')
        # welded.obj is "v z y x"
        $y = [double]$p[2]; $x = [double]$p[3]
        if ($x -lt $mnx) {$mnx=$x}; if ($x -gt $mxx) {$mxx=$x}
        if ($y -lt $mny) {$mny=$y}; if ($y -gt $mxy) {$mxy=$y}
    }
    $spanX = $mxx - $mnx; $spanY = $mxy - $mny
    Write-Host ("  welded bbox     : x[{0:n1},{1:n1}] (span {2:n1})  y[{3:n1},{4:n1}] (span {5:n1})" -f $mnx,$mxx,$spanX,$mny,$mxy,$spanY)
    if ($cubesProcessed -le 1) {
        Write-Host "  world-space check: skipped (single cube; span ~one cube is correct)" -ForegroundColor Gray
    } elseif ($spanX -le ($cubeSize * 1.5) -and $spanY -le ($cubeSize * 1.5)) {
        Write-Host "  WORLD-SPACE CHECK FAILED: $cubesProcessed cubes collapsed into ~one cube -- origins not applied (stacked-cube bug)." -ForegroundColor Red
        $fail = $true
    } else {
        Write-Host "  world-space check: OK (spans multiple cubes)" -ForegroundColor Green
    }
}

# ---- 4. GATHER pre-weld per-cube meshes into one flat dir for MeshLab QC ----
# Each cube's final pre-weld mesh (world-coord, all components) is the step12_final
# stage grid_weld consumes: dump/<cube>/<cube>_step12_final/<cube>_step12_final_all.obj.
# Copy them all into output/<Experiment>/final_cubes/<cube_id>.obj so they can be
# loaded into MeshLab in one shot (they're in world coords, so they assemble in
# place as the un-welded grid) and stepped through layer-by-layer for QC.
Section "GATHER pre-weld per-cube meshes -> output/$Experiment/final_cubes"
$finalDir = Join-Path $outDir 'final_cubes'
New-Item -ItemType Directory -Force -Path $finalDir | Out-Null
$dumpDir = Join-Path $outDir 'dump'
$nGathered = 0
if (Test-Path $dumpDir) {
    foreach ($cubeDir in Get-ChildItem $dumpDir -Directory) {
        $cubeId = $cubeDir.Name
        $src = Join-Path $cubeDir.FullName "${cubeId}_step12_final/${cubeId}_step12_final_all.obj"
        if (Test-Path $src) {
            Copy-Item $src (Join-Path $finalDir "$cubeId.obj") -Force
            $nGathered++
        }
    }
}
Write-Host ("  gathered {0} per-cube mesh(es) -> {1}" -f $nGathered, $finalDir)
Write-Host  "  MeshLab QC: File > Import Mesh, select every *.obj in final_cubes (world-coord, pre-weld)"
if ($nGathered -eq 0) { Write-Host "  (none found -- did grid_pipeline dump step12_final? check output/$Experiment/dump)" -ForegroundColor Yellow }

Write-Host ""
Write-Host ("  artifacts: {0}" -f $outDir)
Write-Host ("    welded.obj, welded.obj.weld_report.json, grid_weld.log, pipeline.log")
Write-Host ("    final_cubes/  ({0} pre-weld per-cube OBJs for MeshLab QC)" -f $nGathered)
Write-Host ""
if ($fail -or $gpExit -ne 0) {
    Write-Host "RESULT: FAIL" -ForegroundColor Red
    exit 1
}
Write-Host "RESULT: OK" -ForegroundColor Green
exit 0
