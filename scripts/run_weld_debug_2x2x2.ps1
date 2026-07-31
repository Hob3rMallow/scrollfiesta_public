#requires -Version 7.0
<#
.SYNOPSIS
  Split a cached per-cube grid dump into small 2x2x2 (configurable) subregions
  and run grid_weld on each, so cross-cube seam-weld bugs can be inspected one
  block at a time.

.DESCRIPTION
  This does NOT re-run the (expensive) per-cube pipeline. It re-welds the
  per-cube OBJ dumps that grid_pipeline already produced -- grid_weld only ever
  reads <cube_id>/<cube_id>_step12_final/<cube_id>_step12_final_all.obj per cube, and
  those were extracted with full halos in the original full-grid run. So a
  subset weld reproduces the EXACT per-seam behaviour of the full weld for every
  seam interior to the block (validated: foldover #2 reproduces as same_dir=1 in
  its 2x2x2 block). The block's outer faces become open boundary (unpaired) --
  expected, they are just the cut edge of the subregion.

  Each member cube is linked into the block's _dump/ via a DIRECTORY JUNCTION
  (instant, no copy). grid_weld --seam-weld runs over that _dump/ and writes
  welded.obj + welded.obj.weld_report.json + welded.obj.bad_edges.obj (and
  .nonmanifold.obj if any) into the block dir. All grid_weld stderr (per-plane
  seam detail) is captured to grid_weld.log.

  A top-level summary.csv / summary.md aggregates every block's manifold audit,
  sorted worst-first (same_dir desc, then non_manifold, then unpaired), so the
  blocks that contain foldovers float to the top.

  SAFETY: _dump/ holds junctions pointing at the real 33 GB dump tree. Cleanup
  removes those junctions NON-RECURSIVELY (link only) before deleting any block
  dir, so we can never recurse through a junction and delete the source dumps.

.PARAMETER DumpDir   Source per-cube dump dir (default output/grid_seamweld/dump).
.PARAMETER OutDir    Test dir for the subregion welds (default output/weld_debug_2x2x2).
.PARAMETER Weld      grid_weld.exe path (default build/Release/grid_weld.exe).
.PARAMETER Size      Window edge in cubes per axis (default 2 -> 2x2x2).
.PARAMETER Stride    Anchor step in cubes. Default 2 (disjoint tiling). Set 1 for
                     overlapping windows so every internal seam is interior to a
                     block (more blocks, some redundant).
.PARAMETER Blocks    Optional explicit min-corner anchors "z,y,x;z,y,x;..." in
                     SOURCE-VOXEL ORIGIN coords. When given, REPLACES auto-tiling;
                     each anchor expands to a Size^3 window clamped to present
                     cubes. Use to add targeted windows around known foldovers.
                     Runs append to OutDir (idempotent) so you can tile first,
                     then add -Blocks for the seams disjoint tiling missed.
.PARAMETER Force     Re-weld blocks that already have a welded.obj (default skips
                     them, so the run is resumable / append-only).
.PARAMETER SummaryOnly  Skip welding; just re-aggregate summary.{csv,md} from the
                        weld_report.json files already on disk.

.EXAMPLE
  pwsh scripts/run_weld_debug_2x2x2.ps1                       # 18 disjoint blocks
.EXAMPLE
  # add the 2 foldover seams that fall on disjoint-tile boundaries
  pwsh scripts/run_weld_debug_2x2x2.ps1 -Blocks "4352,3200,2560;4480,3584,2688"
.EXAMPLE
  pwsh scripts/run_weld_debug_2x2x2.ps1 -SummaryOnly          # rebuild the table
#>
[CmdletBinding()]
param(
    [string] $DumpDir = 'output/grid_seamweld/dump',
    [string] $OutDir  = 'output/weld_debug_2x2x2',
    [string] $Weld    = 'build/Release/grid_weld.exe',
    [int]    $Size    = 2,
    [int]    $Stride  = 2,
    [string] $Blocks  = '',
    [switch] $Force,
    [switch] $SummaryOnly
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

# --- Safe removal: strip junctions (link-only) before recursing ----------
function Remove-BlockDirSafe([string]$dir) {
    if (-not (Test-Path -LiteralPath $dir)) { return }
    $sub = Join-Path $dir '_dump'
    if (Test-Path -LiteralPath $sub) {
        # _dump children are junctions -> delete each link NON-recursively so
        # we never traverse into the real dump tree.
        foreach ($child in Get-ChildItem -LiteralPath $sub -Force -Directory -ErrorAction SilentlyContinue) {
            if ($child.Attributes -band [IO.FileAttributes]::ReparsePoint) {
                [IO.Directory]::Delete($child.FullName, $false)   # link only
            }
        }
    }
    Remove-Item -LiteralPath $dir -Recurse -Force   # now junction-free
}

# --- Enumerate cached cubes ----------------------------------------------
$cubeRe = '^z(\d{5})_y(\d{5})_x(\d{5})$'
$present = @{}            # cube_id -> $true
$zs = [System.Collections.Generic.SortedSet[int]]::new()
$ys = [System.Collections.Generic.SortedSet[int]]::new()
$xsv = [System.Collections.Generic.SortedSet[int]]::new()
foreach ($d in Get-ChildItem -LiteralPath $DumpDir -Directory) {
    if ($d.Name -match $cubeRe) {
        $present[$d.Name] = $true
        [void]$zs.Add([int]$Matches[1]); [void]$ys.Add([int]$Matches[2]); [void]$xsv.Add([int]$Matches[3])
    }
}
if ($present.Count -eq 0) { throw "No z#####_y#####_x##### cube dirs under $DumpDir" }
$zA = @($zs); $yA = @($ys); $xA = @($xsv)
Write-Host ("Cubes present: {0}   z[{1}..{2}] y[{3}..{4}] x[{5}..{6}]" -f `
    $present.Count, $zA[0], $zA[-1], $yA[0], $yA[-1], $xA[0], $xA[-1])

function Id([int]$z,[int]$y,[int]$x) { 'z{0:D5}_y{1:D5}_x{2:D5}' -f $z,$y,$x }

# axis anchors: indices into the sorted axis value array
function AxisAnchors([int]$len) {
    if ($len -le 0) { return @() }
    $last = if ($Stride -ge $Size) { $len - 1 } else { [Math]::Max(0, $len - $Size) }
    $a = @(); $i = 0
    while ($i -le $last) { $a += $i; $i += $Stride }
    return $a
}

# --- Build the window list -----------------------------------------------
# Each window = list of axis-value arrays {zvals, yvals, xvals}; we emit a block
# from the cubes present in the cartesian product.
$windows = @()
if ($Blocks.Trim()) {
    # explicit min-corner anchors -> Size^3 windows clamped to present values
    foreach ($tok in ($Blocks -split ';' | Where-Object { $_.Trim() })) {
        $p = $tok -split ',' | ForEach-Object { [int]$_.Trim() }
        if ($p.Count -ne 3) { throw "Bad -Blocks token '$tok' (need z,y,x)" }
        $zi = [Array]::IndexOf($zA, $p[0]); $yi = [Array]::IndexOf($yA, $p[1]); $xi = [Array]::IndexOf($xA, $p[2])
        if ($zi -lt 0 -or $yi -lt 0 -or $xi -lt 0) { Write-Host "  -Blocks anchor $tok not in grid -- skip" -ForegroundColor Yellow; continue }
        $windows += ,@{
            z = @($zA[$zi..([Math]::Min($zi+$Size-1, $zA.Count-1))])
            y = @($yA[$yi..([Math]::Min($yi+$Size-1, $yA.Count-1))])
            x = @($xA[$xi..([Math]::Min($xi+$Size-1, $xA.Count-1))])
        }
    }
} else {
    foreach ($zi in AxisAnchors $zA.Count) {
      foreach ($yi in AxisAnchors $yA.Count) {
        foreach ($xi in AxisAnchors $xA.Count) {
            $windows += ,@{
                z = @($zA[$zi..([Math]::Min($zi+$Size-1, $zA.Count-1))])
                y = @($yA[$yi..([Math]::Min($yi+$Size-1, $yA.Count-1))])
                x = @($xA[$xi..([Math]::Min($xi+$Size-1, $xA.Count-1))])
            }
        }
      }
    }
}

# --- Weld each block ------------------------------------------------------
if (-not $SummaryOnly) {
    Section ("WELD {0} candidate window(s)  (Size={1} Stride={2})" -f $windows.Count, $Size, $Stride)
    $bi = 0
    foreach ($w in $windows) {
        $bi++
        $members = @()
        foreach ($z in $w.z) { foreach ($y in $w.y) { foreach ($x in $w.x) {
            $id = Id $z $y $x
            if ($present[$id]) { $members += $id }
        }}}
        if ($members.Count -lt 2) { continue }   # need >=2 cubes for a seam
        $zmin = ($w.z | Measure-Object -Minimum).Minimum
        $ymin = ($w.y | Measure-Object -Minimum).Minimum
        $xmin = ($w.x | Measure-Object -Minimum).Minimum
        $name = 'blk_{0}__{1}x{2}x{3}' -f (Id $zmin $ymin $xmin), $w.z.Count, $w.y.Count, $w.x.Count
        $bdir = Join-Path $OutDir $name
        $wobj = Join-Path $bdir 'welded.obj'

        if ((Test-Path $wobj) -and -not $Force) {
            Write-Host ("  [{0}/{1}] {2}: exists, skip (-Force to redo)" -f $bi, $windows.Count, $name) -ForegroundColor DarkGray
            continue
        }
        Remove-BlockDirSafe $bdir
        $sub = Join-Path $bdir '_dump'
        New-Item -ItemType Directory -Force -Path $sub | Out-Null
        foreach ($id in $members) {
            New-Item -ItemType Junction -Path (Join-Path $sub $id) -Target (Join-Path $DumpDir $id) | Out-Null
        }
        # member manifest
        @{ name=$name; window=@{z=$w.z; y=$w.y; x=$w.x}; cubes=$members } |
            ConvertTo-Json -Depth 5 | Set-Content (Join-Path $bdir 'block.json')

        $log = Join-Path $bdir 'grid_weld.log'
        & $Weld $sub $wobj --seam-weld 2>&1 | Tee-Object -FilePath $log | Out-Null
        $rc = $LASTEXITCODE
        $rep = "$wobj.weld_report.json"
        $sd = '?'; $nm = '?'
        if (Test-Path $rep) {
            $j = Get-Content $rep -Raw | ConvertFrom-Json
            $sd = $j.manifold_audit.same_dir_pairs; $nm = $j.manifold_audit.non_manifold
        }
        $tag = if ($rc -gt 255) { 'CRASH' } elseif ($sd -gt 0 -or $nm -gt 0) { 'BAD' } else { 'ok' }
        $col = switch ($tag) { 'CRASH' {'Red'} 'BAD' {'Yellow'} default {'Green'} }
        Write-Host ("  [{0}/{1}] {2}  {3} cubes  same_dir={4} non_manifold={5}  [{6}]" -f `
            $bi, $windows.Count, $name, $members.Count, $sd, $nm, $tag) -ForegroundColor $col
    }
}

# --- Aggregate summary (scan all blocks on disk) -------------------------
Section "SUMMARY"
$rows = @()
foreach ($bdir in Get-ChildItem -LiteralPath $OutDir -Directory -Filter 'blk_*') {
    $rep = Join-Path $bdir.FullName 'welded.obj.weld_report.json'
    if (-not (Test-Path $rep)) { continue }
    $j = Get-Content $rep -Raw | ConvertFrom-Json
    $a = $j.manifold_audit
    $rows += [pscustomobject]@{
        block          = $bdir.Name
        cubes          = $j.cubes_processed
        same_dir       = [int]$a.same_dir_pairs
        non_manifold   = [int]$a.non_manifold
        unpaired       = [int]$a.unpaired
        manifold_pairs = [int]$a.manifold_pairs
        verts          = [int]$j.total_unique_verts
        faces          = [int]$j.total_unique_faces
    }
}
$rows = $rows | Sort-Object -Property @{e='same_dir';Descending=$true}, @{e='non_manifold';Descending=$true}, @{e='unpaired';Descending=$true}

$csv = Join-Path $OutDir 'summary.csv'
$rows | Export-Csv -NoTypeInformation -Path $csv
$md = Join-Path $OutDir 'summary.md'
$bad = @($rows | Where-Object { $_.same_dir -gt 0 -or $_.non_manifold -gt 0 })
$lines = @()
$lines += "# 2x2x2 weld-debug summary ($($rows.Count) blocks, $($bad.Count) with foldovers/non-manifold)"
$lines += ""
$lines += "| block | cubes | same_dir | non_manifold | unpaired | manifold_pairs |"
$lines += "|---|---|---|---|---|---|"
foreach ($r in $rows) {
    $lines += "| $($r.block) | $($r.cubes) | $($r.same_dir) | $($r.non_manifold) | $($r.unpaired) | $($r.manifold_pairs) |"
}
$lines | Set-Content $md

Write-Host ("  {0} blocks welded; {1} with same_dir>0 or non_manifold>0" -f $rows.Count, $bad.Count)
Write-Host ("  summary: {0}" -f $csv)
Write-Host ("           {0}" -f $md)
Write-Host ""
if ($bad.Count -gt 0) {
    Write-Host "  Blocks needing attention (worst first):" -ForegroundColor Yellow
    $bad | Format-Table block, cubes, same_dir, non_manifold, unpaired -AutoSize | Out-String | Write-Host
    Write-Host "  Inspect a block: open <block>/welded.obj (green=weld verts) +"
    Write-Host "  <block>/welded.obj.bad_edges.obj (yellow=same_dir foldover edges) in MeshLab."
} else {
    Write-Host "  No foldovers / non-manifold edges in any block." -ForegroundColor Green
}
