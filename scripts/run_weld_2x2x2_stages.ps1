#requires -Version 7.0
<#
.SYNOPSIS
  Weld ONE 2x2x2 (configurable) cube block from the cached per-cube grid dump,
  dumping every weld stage -- and every interior hole -- as its own OBJ, with a
  wall-clock timeout so a hang is reported instead of spinning forever.

.DESCRIPTION
  This is the focused debug counterpart to run_weld_debug_2x2x2.ps1 (which tiles
  the whole grid). It re-welds a SINGLE block of cached per-cube OBJ dumps
  (no expensive per-cube re-extraction) and runs grid_weld with --dump-stages,
  so you get:

    <block>/stages/<block>_00_concat.obj      plain concat + dedup
    <block>/stages/<block>_01_refine.obj      after seam-band refinement
    <block>/stages/<block>_02_bridge.obj      raw BPA seam bridge
    <block>/stages/<block>_03_preorient.obj   after pre-fold orientation
    <block>/stages/<block>_04_foldcleanup.obj after fold-flap removal
    <block>/stages/<block>_05_pinhole.obj     after pinhole reclose
    <block>/stages/<block>_06_cull.obj        after tiny-component cull
    <block>/stages/<block>_07_cleanup.obj     after flip+collapse cleanup
    <block>/stages/<block>_08_holefill.obj    after interior-hole fill
    <block>/stages/<block>_09_recoarsen.obj   after seam-band recoarsen
    <block>/stages/<block>_10_bandcvt.obj     after seam-band CVT beautification
    <block>/stages/<block>_11_final.obj       final (== welded.obj)
    <block>/stages/holes/<block>_hole####_loop#_##v_in.obj   each hole's INPUT
    <block>/stages/holes/<block>_hole####_loop#_filled.obj   each hole's PATCH

  Each stage is written as it completes, and each hole's *_in.obj is flushed
  BEFORE the fill, so if grid_weld hangs the last stage OBJ + the last *_in.obj
  pinpoint exactly where (the full 100-cube run hangs in hole-fill on a
  pathological CDT polygon -- this lets you catch the equivalent locally).

  grid_weld runs under -TimeoutSec; on timeout the process tree is killed and
  the tail of grid_weld.log (the hung hole) is printed.

  SAFETY: members are linked into <block>/_dump/ via DIRECTORY JUNCTIONS. Cleanup
  deletes those links NON-recursively before removing the block dir, so we can
  never recurse through a junction into the real (large) dump tree.

.PARAMETER Anchor   Min-corner cube id of the block, "z#####_y#####_x#####"
                    (default z04352_y03072_x02560 -- the grid's first corner).
.PARAMETER Size     Block edge in cubes per axis (default 2 -> 2x2x2).
.PARAMETER DumpDir  Source per-cube dump dir (default output/grid_seamweld/dump).
.PARAMETER OutDir   Output root (default output/weld_2x2x2_stages).
.PARAMETER Weld     grid_weld.exe path (default build/Release/grid_weld.exe).
.PARAMETER TimeoutSec  Kill grid_weld after this many seconds (default 600).
.PARAMETER NoHolefill  Pass --no-holefill (bisect: stop before the hanging stage).
.PARAMETER Force    Re-run even if welded.obj already exists.

.EXAMPLE
  pwsh scripts/run_weld_2x2x2_stages.ps1
.EXAMPLE
  pwsh scripts/run_weld_2x2x2_stages.ps1 -Anchor z04480_y03328_x02816 -TimeoutSec 300
.EXAMPLE
  # bisect: get a clean 00..06 + 08 with no hole-fill (skips the hanging stage)
  pwsh scripts/run_weld_2x2x2_stages.ps1 -NoHolefill
#>
[CmdletBinding()]
param(
    [string] $Anchor     = 'z04352_y03072_x02560',
    [int]    $Size       = 2,
    [string] $DumpDir    = 'output/grid_seamweld/dump',
    [string] $OutDir     = 'output/weld_2x2x2_stages',
    [string] $Weld       = 'build/Release/grid_weld.exe',
    [int]    $TimeoutSec = 600,
    [switch] $NoHolefill,
    [switch] $Force
)

$ErrorActionPreference = 'Stop'
$repoRoot = (Resolve-Path "$PSScriptRoot/..").Path
Set-Location $repoRoot
$cubeStep = 128   # source-voxel cube edge

$DumpDir = (Resolve-Path $DumpDir).Path
$OutDir  = Join-Path $repoRoot $OutDir
$Weld    = (Resolve-Path $Weld).Path
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

function Section($t) { Write-Host ""; Write-Host "==== $t ====" -ForegroundColor Cyan }
function Id([int]$z,[int]$y,[int]$x) { 'z{0:D5}_y{1:D5}_x{2:D5}' -f $z,$y,$x }

# --- Safe removal: strip junctions (link-only) before recursing ----------
function Remove-BlockDirSafe([string]$dir) {
    if (-not (Test-Path -LiteralPath $dir)) { return }
    $sub = Join-Path $dir '_dump'
    if (Test-Path -LiteralPath $sub) {
        foreach ($child in Get-ChildItem -LiteralPath $sub -Force -Directory -ErrorAction SilentlyContinue) {
            if ($child.Attributes -band [IO.FileAttributes]::ReparsePoint) {
                [IO.Directory]::Delete($child.FullName, $false)   # link only
            }
        }
    }
    Remove-Item -LiteralPath $dir -Recurse -Force
}

# --- Parse anchor cube id ------------------------------------------------
if ($Anchor -notmatch '^z(\d{5})_y(\d{5})_x(\d{5})$') {
    throw "Bad -Anchor '$Anchor' (need z#####_y#####_x#####)"
}
$z0 = [int]$Matches[1]; $y0 = [int]$Matches[2]; $x0 = [int]$Matches[3]

# --- Build the Size^3 member set, keeping only cubes present on disk -----
$members = @()
$missing = @()
for ($dz = 0; $dz -lt $Size; $dz++) {
  for ($dy = 0; $dy -lt $Size; $dy++) {
    for ($dx = 0; $dx -lt $Size; $dx++) {
        $id = Id ($z0 + $dz*$cubeStep) ($y0 + $dy*$cubeStep) ($x0 + $dx*$cubeStep)
        if (Test-Path -LiteralPath (Join-Path $DumpDir $id)) { $members += $id }
        else { $missing += $id }
    }
  }
}

Write-Host ("Block anchor {0}  Size {1}^3  ->  {2}/{3} cubes present" -f `
    $Anchor, $Size, $members.Count, ($Size*$Size*$Size))
if ($missing.Count -gt 0) {
    Write-Host ("  missing (not in dump): {0}" -f ($missing -join ', ')) -ForegroundColor Yellow
}
if ($members.Count -lt 2) { throw "Need >=2 present cubes for a seam; got $($members.Count)." }

# --- Assemble the block dir ----------------------------------------------
$blockName = 'blk_{0}_{1}x{1}x{1}' -f $Anchor, $Size
$bdir = Join-Path $OutDir $blockName
$wobj = Join-Path $bdir 'welded.obj'
if ((Test-Path $wobj) -and -not $Force) {
    Write-Host "  ${blockName}: welded.obj exists, skip (-Force to redo)" -ForegroundColor DarkGray
    return
}
Remove-BlockDirSafe $bdir
$sub       = Join-Path $bdir '_dump'
$stagesDir = Join-Path $bdir 'stages'
New-Item -ItemType Directory -Force -Path $sub       | Out-Null
New-Item -ItemType Directory -Force -Path $stagesDir | Out-Null
New-Item -ItemType Directory -Force -Path (Join-Path $stagesDir 'holes') | Out-Null
foreach ($id in $members) {
    New-Item -ItemType Junction -Path (Join-Path $sub $id) -Target (Join-Path $DumpDir $id) | Out-Null
}
@{ anchor=$Anchor; size=$Size; cubes=$members; missing=$missing } |
    ConvertTo-Json -Depth 5 | Set-Content (Join-Path $bdir 'block.json')

# --- Run grid_weld with stage dumps, under a wall-clock timeout ----------
Section ("WELD $blockName  ($($members.Count) cubes, timeout ${TimeoutSec}s)")
$weldArgs = @($sub, $wobj, '--dump-stages', $stagesDir, '--stage-prefix', $blockName)
if ($NoHolefill) { $weldArgs += '--no-holefill' }

$log    = Join-Path $bdir 'grid_weld.log'
$outLog = Join-Path $bdir 'grid_weld.out'
Write-Host ("  {0} {1}" -f $Weld, ($weldArgs -join ' '))
Write-Host ("  log -> {0}" -f $log)

$sw = [System.Diagnostics.Stopwatch]::StartNew()
$proc = Start-Process -FilePath $Weld -ArgumentList $weldArgs -NoNewWindow -PassThru `
        -RedirectStandardError $log -RedirectStandardOutput $outLog
$exited = $proc.WaitForExit($TimeoutSec * 1000)
$sw.Stop()
$elapsed = [math]::Round($sw.Elapsed.TotalSeconds, 1)

$hung = $false
if (-not $exited) {
    $hung = $true
    Write-Host ("  *** HUNG: no exit after ${TimeoutSec}s -- killing process tree ***") -ForegroundColor Red
    try { $proc.Kill($true) } catch { try { $proc.Kill() } catch {} }
    Start-Sleep -Milliseconds 500
    $rc = -1
} else {
    $rc = $proc.ExitCode
}

# --- Report --------------------------------------------------------------
Section "RESULT"
Write-Host ("  elapsed: {0}s   exit: {1}   {2}" -f $elapsed, $rc, ($(if($hung){'HUNG'}elseif($rc -eq 0){'ok'}else{'FAILED'})))

# Stage OBJs written (size = a quick "did this stage produce geometry" check)
Write-Host ""
Write-Host "  Stage OBJs:" -ForegroundColor Cyan
$stageObjs = Get-ChildItem -LiteralPath $stagesDir -Filter '*.obj' -File -ErrorAction SilentlyContinue |
             Sort-Object Name
if ($stageObjs) {
    foreach ($f in $stageObjs) {
        Write-Host ("    {0,-40} {1,10:N0} KB" -f $f.Name, ($f.Length/1KB))
    }
} else { Write-Host "    (none)" -ForegroundColor Yellow }

# Per-hole dumps: count *_in.obj (attempted) vs *_filled.obj (succeeded). A hole
# with an *_in.obj but no *_filled.obj is the one that hung / failed.
$holesDir = Join-Path $stagesDir 'holes'
$inObjs     = @(Get-ChildItem -LiteralPath $holesDir -Filter '*_in.obj'     -File -ErrorAction SilentlyContinue)
$filledObjs = @(Get-ChildItem -LiteralPath $holesDir -Filter '*_filled.obj' -File -ErrorAction SilentlyContinue)
Write-Host ""
Write-Host ("  Per-hole dumps: {0} attempted, {1} filled" -f $inObjs.Count, $filledObjs.Count) -ForegroundColor Cyan
if ($inObjs.Count -gt 0) {
    $lastIn = $inObjs | Sort-Object LastWriteTime | Select-Object -Last 1
    $filledStems = $filledObjs | ForEach-Object { ($_.Name -replace '_filled\.obj$','') }
    $lastStem = ($lastIn.Name -replace '_\d+v_in\.obj$','')
    $lastFilled = $filledStems -contains $lastStem
    if ($hung -or -not $lastFilled) {
        Write-Host ("    LAST hole attempted (likely culprit): {0}" -f $lastIn.Name) -ForegroundColor Red
        Write-Host ("      -> open this *_in.obj to see the polygon that hung/failed.") -ForegroundColor Red
    }
}

# weld_report.json (only present if grid_weld reached the audit/write step)
$rep = "$wobj.weld_report.json"
if (Test-Path $rep) {
    $j = Get-Content $rep -Raw | ConvertFrom-Json
    $a = $j.manifold_audit
    Write-Host ""
    Write-Host ("  Manifold audit: unpaired={0} non_manifold={1} same_dir={2} manifold_pairs={3}" -f `
        $a.unpaired, $a.non_manifold, $a.same_dir_pairs, $a.manifold_pairs) -ForegroundColor Cyan
    Write-Host ("  Verts {0:N0}  Faces {1:N0}" -f $j.total_unique_verts, $j.total_unique_faces)
} else {
    Write-Host ""
    Write-Host "  (no weld_report.json -- grid_weld did not reach the write step)" -ForegroundColor Yellow
}

# Tail of the log -- on a hang this shows the hole that hung.
Write-Host ""
Write-Host "  Last 25 log lines ($log):" -ForegroundColor Cyan
if (Test-Path $log) {
    Get-Content $log -Tail 25 | ForEach-Object { Write-Host "    $_" }
}
Write-Host ""
Write-Host ("  Open in MeshLab: {0}" -f $stagesDir)
if ($hung) { exit 2 } elseif ($rc -ne 0) { exit 1 } else { exit 0 }
