#requires -Version 7.0
<#
.SYNOPSIS
  Run the T&L chain-stitch Y-axis zipper on one cube pair.

.DESCRIPTION
  Loads two adjacent cubes' raw_snap meshes, strips wall caps, extracts
  ordered wall-boundary chains, matches A chains to B chains, T&L
  zipper-stitches matched pairs with bridge triangles, and centroid-fan
  caps any unmatched chains. Cube meshes are never modified — bridge and
  cap triangles are pure additions.

.PARAMETER CubeA
  Cube ID of cube A (Y-adjacent: vy_A + 128 = vy_B).
.PARAMETER CubeB
  Cube ID of cube B.
.PARAMETER GridDir
  Grid extract directory (used only for path parity).
.PARAMETER OutputDir
  Where OBJ dumps live: raw_snap meshes read from
  <OutputDir>/obj/<cube_id>/<cube_id>_raw_snap/...
.PARAMETER ExePath
  Path to zipper.exe. Defaults to build/Release/zipper.exe.
#>

[CmdletBinding()]
param(
    [string]$CubeA = 'z04480_y03200_x02816',
    [string]$CubeB = 'z04480_y03328_x02816',
    [string]$GridDir = 'PHerc0139-4x5x5',
    [string]$OutputDir = 'output/PHerc0139-4x5x5',
    [string]$ExePath
)

$ErrorActionPreference = 'Stop'

$repoRoot = (Resolve-Path "$PSScriptRoot/..").Path
if (-not $ExePath) { $ExePath = Join-Path $repoRoot 'build/Release/zipper.exe' }
if (-not (Test-Path $ExePath)) { throw "zipper.exe not found at $ExePath" }

$outputDir = Join-Path $repoRoot $OutputDir
$objDir    = Join-Path $outputDir 'obj'

function Parse-Origin([string]$id) {
    if ($id -notmatch '^z(\d+)_y(\d+)_x(\d+)$') { throw "Bad cube id: $id" }
    return [pscustomobject]@{ z=[int]$Matches[1]; y=[int]$Matches[2]; x=[int]$Matches[3] }
}
$oA = Parse-Origin $CubeA
$oB = Parse-Origin $CubeB
if ($oA.z -ne $oB.z -or $oA.x -ne $oB.x -or ($oB.y - $oA.y) -ne 128) {
    throw "Cubes are not Y-adjacent. A=$($oA), B=$($oB). Need vyB - vyA = 128, same vz/vx."
}

$aObj = Join-Path $objDir "$CubeA/$($CubeA)_raw_snap/$($CubeA)_raw_snap_all.obj"
$bObj = Join-Path $objDir "$CubeB/$($CubeB)_raw_snap/$($CubeB)_raw_snap_all.obj"
foreach ($p in $aObj, $bObj) {
    if (-not (Test-Path $p)) { throw "Missing: $p" }
}

$outObj = Join-Path $outputDir "zipper_$($CubeA)_$($CubeB).obj"
Write-Host "Zipper input:"
Write-Host "  A obj: $aObj"
Write-Host "  B obj: $bObj"
Write-Host "  Out  : $outObj"
Write-Host "  Origin (cube A): Z=$($oA.z) Y=$($oA.y) X=$($oA.x)"
Write-Host ''

& $ExePath `
    $aObj $bObj $outObj `
    --origin-z $oA.z `
    --origin-y $oA.y `
    --origin-x $oA.x
$rc = $LASTEXITCODE
if ($rc -ne 0) { throw "zipper.exe failed with rc=$rc" }

Write-Host ''
Write-Host "Wrote $outObj"
Get-ChildItem $outObj | Format-List FullName, Length, LastWriteTime
