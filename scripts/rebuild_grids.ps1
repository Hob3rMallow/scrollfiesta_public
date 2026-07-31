#requires -Version 7.0
<#
.SYNOPSIS
  Rebuild the 4x5x5 and 4x21x21 grids from scratch: re-mesh every cube, then
  weld each grid into one giant mesh.

.DESCRIPTION
  Two things, nothing else -- for each grid:
    1. RUN  every cube through cube_mesh (per-cube surface mesh), and
    2. WELD the per-cube meshes into one combined mesh.

  The two grids weld DIFFERENTLY, because they are different scales:

    4x5x5   (100 cubes)   FLAT weld. grid_pipeline.exe meshes every cube then
                          runs grid_weld once over the whole dump -> a single
                          full-resolution welded.obj. 100 cubes fit in memory
                          (~1.9 GB), so no decimation is needed.

    4x21x21 (1753 cubes)  HIERARCHICAL weld. A flat monolith weld of 1753 cubes
                          is memory-infeasible, so we use the LOD-pyramid welder
                          (hierarchical_weld.exe): partition the cubes into
                          fanout^3 (2x2x2) blocks, weld each block, decimate it,
                          and treat each decimated block as a node for the next
                          level -- repeating until one node remains. Every weld
                          is bounded to <= 8 inputs and each level is decimated
                          before the next, so a full-res monolith is never
                          materialized. The giant cube is the top LOD tier.

    -BigWeld flat         forces the 4x21x21 down the flat path too (will very
                          likely OOM on 1753 cubes; provided only as an escape).

  Pipeline per grid:
    flat  : grid_pipeline <grid> <out> [--no-qem]                     (mesh+weld)
    hier  : grid_pipeline <grid> <out> --skip-weld                    (mesh only)
            hierarchical_weld <out>/dump <out> --fanout 2 --keep-ratio R ... (weld)

  Grids (fixed):
    4x5x5    PHerc0139-4x5x5          -> output/rebuild_4x5x5
    4x21x21  PHerc0139-4x21x21-cubes  -> output/rebuild_4x21x21

  The 4x21x21 run is a big job (per-cube meshing is a few hours on the 7950X).
  Leaf meshes are written incrementally, so if the pyramid weld dies you can
  re-run with -Resume (skips already-meshed cubes AND already-welded blocks).

.PARAMETER Grids          Which grids: both (default), 4x5x5, 4x21x21.
.PARAMETER BigWeld        4x21x21 weld path: hierarchical (default) or flat.
.PARAMETER Halo           Halo voxels per cube (default 13 = MLS R+1).
.PARAMETER ThreadsPerCube OpenMP threads per cube_mesh subprocess (default 1).
.PARAMETER MaxConcurrent  Max concurrent cube subprocesses (default 0 = cores/tpc).
.PARAMETER MaxCubes       Stop after N cubes per grid (default 0 = all). Smoke test.
.PARAMETER NoQem          Disable QEM decimation in cube_mesh (on by default).
.PARAMETER KeepRatio      Hierarchical: keep this fraction of faces per level (0.25).
.PARAMETER Fanout         Hierarchical: block edge in nodes (default 2 = 2x2x2).
.PARAMETER BlockConcurrent Hierarchical: concurrent block welds (default 48).
.PARAMETER SimplifyTop    Hierarchical: also decimate + re-orient the top tier.
.PARAMETER SkipWeld       Re-mesh cubes but do NOT weld (defer the giant cube).
.PARAMETER Resume         Skip already-meshed cubes / already-welded blocks.
.PARAMETER SkipBuild      Use binaries already in build/Release, don't rebuild.
.PARAMETER Clean          Clean each project before building.
.PARAMETER Force          Overwrite an existing output/rebuild_<grid> dir.

.EXAMPLE
  pwsh scripts/rebuild_grids.ps1                          # both grids, correct weld each
.EXAMPLE
  pwsh scripts/rebuild_grids.ps1 -Grids 4x21x21 -Resume   # continue the big pyramid
.EXAMPLE
  pwsh scripts/rebuild_grids.ps1 -Grids 4x5x5 -MaxCubes 2 -SkipBuild -Force   # smoke
#>
[CmdletBinding()]
param(
    [ValidateSet('both','4x5x5','4x21x21')]
    [string]$Grids          = 'both',
    [ValidateSet('hierarchical','flat')]
    [string]$BigWeld        = 'hierarchical',
    [int]   $Halo           = 13,
    [int]   $ThreadsPerCube = 1,
    [int]   $MaxConcurrent  = 0,
    [int]   $MaxCubes       = 0,
    [switch]$NoQem,
    [double]$KeepRatio      = 0.25,
    [int]   $Fanout         = 2,
    [int]   $BlockConcurrent= 48,
    [switch]$SimplifyTop,
    [switch]$SkipWeld,
    [switch]$Resume,
    [switch]$SkipBuild,
    [switch]$Clean,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path "$PSScriptRoot/..").Path
Set-Location $repoRoot   # keep relative dependency and output paths stable
$msbuild  = 'C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe'
$cubeSize = 128.0        # source-voxel cube edge; seam planes sit at multiples of this

# Weld-only winding/phase anchors (full-scroll source-space voxels). grid_pipeline
# deliberately does not use these during per-cube BPA: doing so fragmented the
# byte-identical z04480_y03328_x02816 input from 5 coherent sheets into 24 pieces.
# hierarchical_weld uses the same geometry after leaf meshing.
$umbY = 3405.0; $umbX = 2878.0; $wrapPitch = 9.5

function Section($t) { Write-Host ""; Write-Host "==== $t ====" -ForegroundColor Cyan }

# Fixed grid table: id -> (source grid dir, output dir, default weld mode).
$gridTable = [ordered]@{
    '4x5x5'   = @{ Dir = 'PHerc0139-4x5x5';         Out = 'output/rebuild_4x5x5';   Weld = 'flat' }
    '4x21x21' = @{ Dir = 'PHerc0139-4x21x21-cubes'; Out = 'output/rebuild_4x21x21'; Weld = $BigWeld }
}
$selected = if ($Grids -eq 'both') { @('4x5x5','4x21x21') } else { @($Grids) }
$anyHier  = $selected | Where-Object { $gridTable[$_].Weld -eq 'hierarchical' }

# ---- Binaries ------------------------------------------------------------
$exe        = Join-Path $repoRoot 'build/Release/cube_mesh.exe'
$weld       = Join-Path $repoRoot 'build/Release/grid_weld.exe'
$gpipe      = Join-Path $repoRoot 'build/Release/grid_pipeline.exe'
$hierExe    = Join-Path $repoRoot 'build/Release/hierarchical_weld.exe'
$qslimExe   = Join-Path $repoRoot 'build/Release/qslim_obj.exe'
$manifoldEx = Join-Path $repoRoot 'build/Release/manifold_check.exe'

# ---- 1. BUILD (once) -----------------------------------------------------
$projects = [System.Collections.Generic.List[string]]::new()
$projects.AddRange([string[]]@('cube_mesh.vcxproj','grid_weld.vcxproj','grid_pipeline.vcxproj'))
if ($anyHier) {
    # The pyramid welder shells out to grid_weld + qslim_obj + manifold_check.
    $projects.AddRange([string[]]@('hierarchical_weld.vcxproj','qslim_obj.vcxproj','manifold_check.vcxproj'))
}

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

$needBins = @($exe, $weld, $gpipe)
if ($anyHier) { $needBins += @($hierExe, $qslimExe, $manifoldEx) }
foreach ($p in $needBins) { if (-not (Test-Path $p)) { throw "Missing binary: $p" } }

# Sanity: the orchestrators' own unit tests must pass before we drive them.
if (-not $SkipBuild) {
    foreach ($t in @(@{N='grid_pipeline'; P=$gpipe}) + $(if ($anyHier) { @(@{N='hierarchical_weld'; P=$hierExe}) } else { @() })) {
        Write-Host "  $($t.N) --selftest ..."
        & $t.P --selftest | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "$($t.N) --selftest FAILED" }
    }
    Write-Host "  selftests OK" -ForegroundColor Green
}

# ---- SUMMARY: flat weld (single welded.obj) ------------------------------
function Write-FlatSummary {
    param([string]$Id, [string]$OutDir, [int]$GpExit, [double]$Minutes, [bool]$Welded)
    $result = [pscustomobject]@{ Id = $Id; Ok = $true; Out = '' }
    Section "SUMMARY  $Id  (flat weld)"
    Write-Host ("  grid_pipeline exit: {0}   elapsed: {1:n1} min" -f $GpExit, $Minutes)
    if ($GpExit -ne 0) { $result.Ok = $false }

    $summaryCsv = Join-Path $OutDir 'pipeline_summary.csv'
    if (Test-Path $summaryCsv) {
        $nOk = 0; $nReject = 0; $nFail = 0
        foreach ($row in Import-Csv $summaryCsv) {
            switch ([int]$row.exit_code) { 0 { $nOk++ } -2 { $nReject++ } default { $nFail++ } }
        }
        Write-Host ("  cubes meshed    : {0} ok, {1} rejected (solid-slab/empty), {2} failed" -f $nOk, $nReject, $nFail)
        if ($nFail -gt 0) { $result.Ok = $false }
    } else {
        Write-Host "  pipeline_summary.csv MISSING" -ForegroundColor Red; $result.Ok = $false
    }

    if (-not $Welded) { Write-Host "  weld: skipped (-SkipWeld)" -ForegroundColor Yellow; return $result }

    $weldObj    = Join-Path $OutDir 'welded.obj'
    $reportPath = "$weldObj.weld_report.json"
    $result.Out = $weldObj
    $cubesProcessed = 0
    if (Test-Path $reportPath) {
        $r = Get-Content $reportPath -Raw | ConvertFrom-Json
        $a = $r.manifold_audit
        $cubesProcessed = [int]$r.cubes_processed
        Write-Host ("  welded cubes    : {0}" -f $r.cubes_processed)
        Write-Host ("  unique verts    : {0:n0}" -f $r.total_unique_verts)
        Write-Host ("  unique faces    : {0:n0}" -f $r.total_unique_faces)
        Write-Host ("  manifold pairs  : {0:n0}" -f $a.manifold_pairs)
        Write-Host ("  unpaired edges  : {0:n0}  (open boundary)" -f $a.unpaired)
        Write-Host ("  non-manifold    : {0}" -f $a.non_manifold) -ForegroundColor ($(if ($a.non_manifold -gt 0) {'Yellow'} else {'Gray'}))
        Write-Host ("  same-dir pairs  : {0}  (foldovers; must be 0)" -f $a.same_dir_pairs) -ForegroundColor ($(if ($a.same_dir_pairs -gt 0) {'Red'} else {'Gray'}))
        if ($a.same_dir_pairs -gt 0) { $result.Ok = $false }
    } else {
        Write-Host "  NO weld report found ($reportPath)" -ForegroundColor Red; $result.Ok = $false
    }

    if (Test-Path $weldObj) {
        $mnx = [double]::PositiveInfinity; $mxx = [double]::NegativeInfinity
        $mny = [double]::PositiveInfinity; $mxy = [double]::NegativeInfinity
        foreach ($line in [System.IO.File]::ReadLines($weldObj)) {
            if ($line.Length -lt 2 -or $line[0] -ne 'v' -or $line[1] -ne ' ') { continue }
            $p = $line.Split(' ')            # welded.obj is "v z y x"
            $y = [double]$p[2]; $x = [double]$p[3]
            if ($x -lt $mnx) {$mnx=$x}; if ($x -gt $mxx) {$mxx=$x}
            if ($y -lt $mny) {$mny=$y}; if ($y -gt $mxy) {$mxy=$y}
        }
        $spanX = $mxx - $mnx; $spanY = $mxy - $mny
        Write-Host ("  welded bbox     : x[{0:n1},{1:n1}] (span {2:n1})  y[{3:n1},{4:n1}] (span {5:n1})" -f $mnx,$mxx,$spanX,$mny,$mxy,$spanY)
        if ($cubesProcessed -le 1) {
            Write-Host "  world-space     : n/a (single cube)" -ForegroundColor Gray
        } elseif ($spanX -le ($cubeSize * 1.5) -and $spanY -le ($cubeSize * 1.5)) {
            Write-Host "  WORLD-SPACE CHECK FAILED: cubes collapsed into ~one cube (origins not applied)." -ForegroundColor Red
            $result.Ok = $false
        } else {
            Write-Host "  world-space     : OK (spans multiple cubes)" -ForegroundColor Green
        }
    } else {
        Write-Host "  welded.obj MISSING" -ForegroundColor Red; $result.Ok = $false
    }
    return $result
}

# ---- SUMMARY: hierarchical weld (LOD pyramid) ----------------------------
function Write-PyramidSummary {
    param([string]$Id, [string]$OutDir, [int]$MeshExit, [int]$HierExit, [double]$Minutes, [bool]$Welded)
    $result = [pscustomobject]@{ Id = $Id; Ok = $true; Out = '' }
    Section "SUMMARY  $Id  (hierarchical LOD pyramid)"
    Write-Host ("  mesh exit: {0}   pyramid exit: {1}   elapsed: {2:n1} min" -f $MeshExit, $HierExit, $Minutes)
    if ($MeshExit -ne 0) { $result.Ok = $false }

    $summaryCsv = Join-Path $OutDir 'pipeline_summary.csv'
    if (Test-Path $summaryCsv) {
        $nOk = 0; $nReject = 0; $nFail = 0
        foreach ($row in Import-Csv $summaryCsv) {
            switch ([int]$row.exit_code) { 0 { $nOk++ } -2 { $nReject++ } default { $nFail++ } }
        }
        Write-Host ("  leaf cubes      : {0} ok, {1} rejected (solid-slab/empty), {2} failed" -f $nOk, $nReject, $nFail)
        if ($nFail -gt 0) { $result.Ok = $false }
    } else {
        Write-Host "  pipeline_summary.csv MISSING" -ForegroundColor Red; $result.Ok = $false
    }

    if (-not $Welded) { Write-Host "  pyramid weld: skipped (-SkipWeld)" -ForegroundColor Yellow; return $result }
    if ($HierExit -ne 0) { $result.Ok = $false }

    # Per-level pyramid table from each level<N>/_summary.csv.
    $levels = Get-ChildItem "$OutDir/level*" -Directory -ErrorAction SilentlyContinue |
              Sort-Object { [int]($_.Name -replace 'level','') }
    if (-not $levels) {
        Write-Host "  NO level* tiers produced" -ForegroundColor Red; $result.Ok = $false; return $result
    }
    "  {0,-7} {1,6} {2,5} {3,5} {4,5} {5,14} {6,9}" -f "level","nodes","ok","warn","fail","total_faces","worst_nm" | Write-Host
    $topLevel = $null
    foreach ($lvl in $levels) {
        $csv = Join-Path $lvl.FullName '_summary.csv'
        if (-not (Test-Path $csv)) { continue }
        $rows = Import-Csv $csv
        $ok   = ($rows | Where-Object status -eq 'ok').Count
        $warn = ($rows | Where-Object status -eq 'warn').Count
        $fail = ($rows | Where-Object status -eq 'fail').Count
        $tf   = ($rows | Measure-Object final_faces -Sum).Sum
        $wnm  = ($rows | Measure-Object non_manifold -Maximum).Maximum
        "  {0,-7} {1,6} {2,5} {3,5} {4,5} {5,14:n0} {6,9}" -f $lvl.Name,$rows.Count,$ok,$warn,$fail,$tf,$wnm | Write-Host
        if ($fail -gt 0) { $result.Ok = $false }
        $topLevel = $lvl
    }

    # The giant cube = the single node OBJ at the top tier.
    if ($topLevel) {
        $topNode = Get-ChildItem (Join-Path $topLevel.FullName '*/*_all.obj') -ErrorAction SilentlyContinue |
                   Sort-Object Length -Descending | Select-Object -First 1
        $nTop = @(Get-ChildItem $topLevel.FullName -Directory -ErrorAction SilentlyContinue).Count
        if ($topNode) {
            Write-Host ("  giant cube      : {0}" -f $topNode.FullName) -ForegroundColor Green
            $result.Out = $topNode.FullName
            if ($nTop -gt 1) {
                Write-Host ("  NOTE: top tier has {0} nodes (did not converge to 1)" -f $nTop) -ForegroundColor Yellow
            }
        } else {
            Write-Host "  top-tier node OBJ MISSING" -ForegroundColor Red; $result.Ok = $false
        }
    }
    return $result
}

# ---- 2. RUN + WELD each grid --------------------------------------------
$results = @()
foreach ($id in $selected) {
    $g       = $gridTable[$id]
    $mode    = $g.Weld
    $gridDir = Join-Path $repoRoot $g.Dir
    $outDir  = Join-Path $repoRoot $g.Out
    if (-not (Test-Path (Join-Path $gridDir 'cubes_PRED'))) {
        throw "Grid '$id': $gridDir/cubes_PRED not found"
    }

    if ((Test-Path $outDir) -and -not $Resume) {
        if (-not $Force) {
            throw "$($g.Out) already exists. Use -Force to overwrite or -Resume to continue."
        }
        Write-Host "  -Force: removing existing $outDir"
        Remove-Item -Recurse -Force $outDir
    }
    New-Item -ItemType Directory -Force -Path $outDir | Out-Null

    # -- Stage A: mesh every cube (grid_pipeline). Hierarchical meshes leaves
    #    only (--skip-weld); flat lets grid_pipeline weld inline unless -SkipWeld.
    $gpArgs = @(
        $gridDir, $outDir,
        '--halo', $Halo, '--threads-per-cube', $ThreadsPerCube,
        '--exe', $exe, '--weld', $weld,
        # Seam winding/phase-gate anchors for the FLAT weld (grid_pipeline scopes
        # these to grid_weld only). The hierarchical path gets
        # the same values via hierarchical_weld's own flags below. Without them
        # the weld runs gates-off: Jul-12..17 builds shipped 661 seam-band gaps
        # on the 4x5x5 instead of 16 (the upstream weld-hole complaint).
        '--umb-y', $umbY, '--umb-x', $umbX, '--wrap-pitch', $wrapPitch
    )
    if ($MaxConcurrent -gt 0) { $gpArgs += @('--max-concurrent', $MaxConcurrent) }
    if ($MaxCubes -gt 0)      { $gpArgs += @('--max-cubes', $MaxCubes) }
    if ($NoQem)               { $gpArgs += '--no-qem' } else { $gpArgs += '--qem' }
    if ($Resume)              { $gpArgs += '--skip-existing' }
    if ($mode -eq 'hierarchical' -or $SkipWeld) { $gpArgs += '--skip-weld' }

    Section "REBUILD $id  ($($g.Dir) -> $($g.Out))  weld=$mode"
    $log = Join-Path $outDir 'rebuild.log'
    Write-Host "  [mesh] grid_pipeline $($gpArgs -join ' ')"
    $started = Get-Date
    & $gpipe @gpArgs 2>&1 | Tee-Object -FilePath $log
    $meshExit = $LASTEXITCODE

    if ($mode -eq 'flat') {
        $minutes = ((Get-Date) - $started).TotalMinutes
        $results += (Write-FlatSummary -Id $id -OutDir $outDir -GpExit $meshExit -Minutes $minutes -Welded (-not $SkipWeld))
        continue
    }

    # -- Stage B: hierarchical 2x2x2 LOD-pyramid weld over the leaf dump.
    $hierExit = 0
    if ($SkipWeld) {
        Write-Host "  [weld] skipped (-SkipWeld); leaves meshed only" -ForegroundColor Yellow
    } elseif ($meshExit -ne 0) {
        Write-Host "  [weld] SKIPPED: leaf meshing exited $meshExit" -ForegroundColor Red
        $hierExit = $meshExit
    } else {
        $leafDump = Join-Path $outDir 'dump'
        $hierArgs = @(
            $leafDump, $outDir,
            '--leaf-stage', 'step12_final',
            '--fanout', $Fanout,
            '--keep-ratio', $KeepRatio,
            '--stop-nodes', 1,
            '--max-levels', 16,
            '--max-concurrent', $BlockConcurrent,
            '--block-timeout', 900,
            '--umbilicus-y', $umbY, '--umbilicus-x', $umbX, '--wrap-pitch', $wrapPitch,
            '--grid-weld', $weld, '--qslim', $qslimExe, '--manifold', $manifoldEx
        )
        if ($Resume)      { $hierArgs += '--skip-existing' }
        if ($SimplifyTop) { $hierArgs += '--simplify-top' }
        Write-Host "  [weld] hierarchical_weld $($hierArgs -join ' ')"
        & $hierExe @hierArgs 2>&1 | Tee-Object -FilePath $log -Append
        $hierExit = $LASTEXITCODE
    }

    $minutes = ((Get-Date) - $started).TotalMinutes
    $results += (Write-PyramidSummary -Id $id -OutDir $outDir -MeshExit $meshExit -HierExit $hierExit -Minutes $minutes -Welded (-not $SkipWeld))
}

# ---- 3. FINAL SUMMARY ----------------------------------------------------
Section "DONE"
$allOk = $true
foreach ($res in $results) {
    $tag = if ($res.Ok) { 'OK' } else { 'FAIL' }
    $col = if ($res.Ok) { 'Green' } else { 'Red' }
    if (-not $res.Ok) { $allOk = $false }
    $where = if ($res.Out) { $res.Out } else { '(no weld output)' }
    Write-Host ("  {0,-8} {1,-5} -> {2}" -f $res.Id, $tag, $where) -ForegroundColor $col
}
Write-Host ""
if ($allOk) { Write-Host "RESULT: OK" -ForegroundColor Green; exit 0 }
Write-Host "RESULT: FAIL" -ForegroundColor Red; exit 1
