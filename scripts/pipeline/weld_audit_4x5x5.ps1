<#
  weld_audit_4x5x5.ps1 -- Directional pairwise-weld audit of the 4x5x5 grid.

  For every adjacent cube pair in each direction (x, y, z), weld JUST that pair
  (grid_weld --subgrid over the tight 2-cube bbox, with the production winding
  gate armed), then measure single-triangle / small gaps AT the shared seam with
  seam_gap, and record the manifold report. One row per pair -> results.csv.

  PowerShell is a RUNNER here only; all geometry analysis is in the C tools
  (grid_weld, seam_gap). See feedback_c_tools_not_powershell.

  The pairwise weld isolates ONE seam per run (all other faces are outer
  perimeter), so seam_gap's near-seam band cleanly separates weld defects from
  each cube's own trim boundary. Caveat: a 2-cube bridge derives its adaptive rho
  from just those two cubes' boundary edges, so a flagged defect should be
  re-checked against the full-grid weld before it is called production-real.
#>
[CmdletBinding()]
param(
  [string]$Dump     = "output/grid_4x5x5_20260709/dump",
  [string]$OutDir   = "output/weld_audit_4x5x5",
  [int]   $Parallel = 10,
  [double]$Band     = 6.0
)
$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
Set-Location $repo

$exeWeld = Join-Path $repo "build\Release\grid_weld.exe"
$exeGap  = Join-Path $repo "build\Release\seam_gap.exe"
foreach ($e in @($exeWeld,$exeGap)) { if (-not (Test-Path $e)) { throw "missing exe: $e" } }
$DumpFull = Join-Path $repo $Dump
if (-not (Test-Path $DumpFull)) { throw "missing dump dir: $DumpFull" }

$weldDir = Join-Path $repo "$OutDir\welds"
$gapDir  = Join-Path $repo "$OutDir\gaps"
New-Item -ItemType Directory -Force -Path $weldDir,$gapDir | Out-Null

# Grid axes (source-voxel origins), step 128. NB PowerShell variables are
# case-insensitive, so the axis arrays MUST NOT be named $Z/$Y/$X -- the loop
# scalars $z/$y/$x would clobber them (a $x=$X[$i] assignment turns $X into a
# scalar and kills the loop after one iteration).
$AxZ = 4352,4480,4608,4736
$AxY = 3072,3200,3328,3456,3584
$AxX = 2560,2688,2816,2944,3072

$pairs = [System.Collections.Generic.List[object]]::new()
# x-direction: seam plane = X = x+128 (axis 2)
foreach ($z in $AxZ) { foreach ($y in $AxY) { for ($i=0; $i -lt $AxX.Count-1; $i++) {
  $x=$AxX[$i]; $xn=$AxX[$i+1]
  $pairs.Add([pscustomobject]@{ dir='x'; z0=$z;z1=$z; y0=$y;y1=$y; x0=$x;x1=$xn; axis=2; coord=$xn
    id=("x_z{0:D5}_y{1:D5}_x{2:D5}" -f $z,$y,$x) }) } } }
# y-direction: seam plane = Y = y+128 (axis 1)
foreach ($z in $AxZ) { foreach ($x in $AxX) { for ($i=0; $i -lt $AxY.Count-1; $i++) {
  $y=$AxY[$i]; $yn=$AxY[$i+1]
  $pairs.Add([pscustomobject]@{ dir='y'; z0=$z;z1=$z; y0=$y;y1=$yn; x0=$x;x1=$x; axis=1; coord=$yn
    id=("y_z{0:D5}_y{1:D5}_x{2:D5}" -f $z,$y,$x) }) } } }
# z-direction: seam plane = Z = z+128 (axis 0)
foreach ($y in $AxY) { foreach ($x in $AxX) { for ($i=0; $i -lt $AxZ.Count-1; $i++) {
  $z=$AxZ[$i]; $zn=$AxZ[$i+1]
  $pairs.Add([pscustomobject]@{ dir='z'; z0=$z;z1=$zn; y0=$y;y1=$y; x0=$x;x1=$x; axis=0; coord=$zn
    id=("z_z{0:D5}_y{1:D5}_x{2:D5}" -f $z,$y,$x) }) } } }

Write-Host ("weld_audit: {0} pairs ({1} x, {2} y, {3} z) at {4}-way" -f `
  $pairs.Count, ($pairs|?{$_.dir -eq 'x'}).Count, ($pairs|?{$_.dir -eq 'y'}).Count, `
  ($pairs|?{$_.dir -eq 'z'}).Count, $Parallel)

$sw = [System.Diagnostics.Stopwatch]::StartNew()
$results = $pairs | ForEach-Object -ThrottleLimit $Parallel -Parallel {
  $p = $_
  $exeWeld = $using:exeWeld; $exeGap = $using:exeGap
  $DumpFull = $using:DumpFull; $weldDir = $using:weldDir; $gapDir = $using:gapDir
  $Band = $using:Band
  # production winding gate (arms only with these three set)
  $env:SEAM_UMBILICUS_Y = '3405'; $env:SEAM_UMBILICUS_X = '2878'; $env:SEAM_WRAP_PITCH = '9.5'

  $obj = Join-Path $weldDir "$($p.id).obj"
  $log = Join-Path $weldDir "$($p.id).log"
  $gapObj = Join-Path $gapDir "$($p.id)_gaps.obj"

  $wArgs = @($DumpFull, $obj, '--subgrid', $p.z0,$p.z1,$p.y0,$p.y1,$p.x0,$p.x1)
  & $exeWeld @wArgs *> $log
  $weldRc = $LASTEXITCODE

  $rep = $null
  $rj = "$obj.weld_report.json"
  if (Test-Path $rj) { try { $rep = Get-Content $rj -Raw | ConvertFrom-Json } catch {} }

  $gj = $null; $gapRc = $null
  if (Test-Path $obj) {
    $gout = & $exeGap $obj '--plane' $p.axis $p.coord '--band' $Band '--out' $gapObj 2>&1
    $gapRc = $LASTEXITCODE
    $line = ($gout | Select-String '^SEAMGAP_JSON: ' | Select-Object -First 1)
    if ($line) { try { $gj = ($line.ToString() -replace '^SEAMGAP_JSON: ','') | ConvertFrom-Json } catch {} }
  }

  [pscustomobject]@{
    id           = $p.id
    dir          = $p.dir
    seam_axis    = $p.axis
    seam_coord   = $p.coord
    weld_rc      = $weldRc
    verts        = if ($rep) { $rep.total_unique_verts } else { $null }
    faces        = if ($rep) { $rep.total_unique_faces } else { $null }
    unpaired     = if ($rep) { $rep.manifold_audit.unpaired } else { $null }
    non_manifold = if ($rep) { $rep.manifold_audit.non_manifold } else { $null }
    same_dir     = if ($rep) { $rep.manifold_audit.same_dir_pairs } else { $null }
    pinch        = if ($rep) { $rep.manifold_audit.pinch_verts } else { $null }
    tri_seam     = if ($gj) { $gj.tri_gaps_seam } else { $null }
    tri_straddle = if ($gj) { $gj.tri_straddle } else { $null }
    small_seam   = if ($gj) { $gj.small_gaps_seam } else { $null }
    big_seam     = if ($gj) { $gj.big_seam_loops } else { $null }
    seam_loops   = if ($gj) { $gj.seam_loops } else { $null }
    nm_edges_gap = if ($gj) { $gj.nm_edges } else { $null }
    max_loop     = if ($gj) { $gj.max_loop } else { $null }
    pinch_steps  = if ($gj) { $gj.pinch_steps } else { $null }
  }
}
$sw.Stop()

$csv = Join-Path $repo "$OutDir\results.csv"
$results | Sort-Object -Property @{e='tri_seam';Descending=$true}, @{e='small_seam';Descending=$true}, `
  @{e='big_seam';Descending=$true}, @{e='non_manifold';Descending=$true} |
  Export-Csv -Path $csv -NoTypeInformation
Write-Host ("DONE in {0:N0}s -> {1}" -f $sw.Elapsed.TotalSeconds, $csv)

# quick console summary
$flag = $results | Where-Object { $_.tri_seam -gt 0 -or $_.small_seam -gt 0 -or $_.big_seam -gt 0 `
  -or $_.non_manifold -gt 0 -or $_.pinch -gt 0 -or $_.same_dir -gt 0 }
Write-Host ("flagged pairs (any seam gap / manifold defect): {0} / {1}" -f $flag.Count, $results.Count)
$results | Measure-Object -Property tri_seam,small_seam,big_seam -Sum | ForEach-Object {
  Write-Host ("  total {0} = {1}" -f $_.Property, $_.Sum) }
