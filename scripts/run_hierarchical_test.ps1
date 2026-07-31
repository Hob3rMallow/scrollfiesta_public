<#
  run_hierarchical_test.ps1 -- live end-to-end test for hierarchical_weld.exe.

  Validates the LOD-pyramid orchestrator on a small, all-present corner of an
  existing leaf dump (no re-meshing):
    1. --selftest passes (pure partitioner math incl. the {242,72,12,6,4,1} counts).
    2. A 2x2x4 leaf corner -> Level 1 = 2 blocks, each weld+simplify, both MANIFOLD.
    3. Level 2 welds the 2 L1 blocks into 1 node; its re-orientation HEALS the
       qslim winding blemish -> L2 same_dir < sum(L1 same_dir).

  Junctions (mklink /J) the source cube dirs so nothing is copied; they are
  removed with `rmdir` (never Remove-Item -Recurse, which would follow the link
  into the real dump).

  Usage: pwsh scripts/run_hierarchical_test.ps1 [-SrcDump <dir>] [-KeepRatio 0.25] [-SkipBuild]
#>
[CmdletBinding()]
param(
  [string]$SrcDump   = "output/grid_4x21x21_foldfix/dump",
  [double]$KeepRatio = 0.25,
  [switch]$SkipBuild
)
$ErrorActionPreference = "Stop"
Set-Location (Split-Path $PSScriptRoot -Parent)   # repo root

$exe      = "build/Release/hierarchical_weld.exe"
$msbuild  = "C:\Program Files\Microsoft Visual Studio\2022\Professional\MSBuild\Current\Bin\MSBuild.exe"
$testDump = "output/_hw_test_dump"
$testOut  = "output/_hw_test_out"
$fails    = 0

function Fail($m) { Write-Host "FAIL: $m" -ForegroundColor Red; $script:fails++ }
function Pass($m) { Write-Host "PASS: $m" -ForegroundColor Green }

# ---- 1. build ----
if (-not $SkipBuild) {
  Write-Host "== building hierarchical_weld.vcxproj (Release x64) =="
  & $msbuild hierarchical_weld.vcxproj /p:Configuration=Release /p:Platform=x64 /m /v:minimal
  if ($LASTEXITCODE -ne 0) { throw "build failed" }
}
if (-not (Test-Path $exe)) { throw "$exe not found" }

# ---- 2. selftest ----
Write-Host "== --selftest =="
& $exe --selftest
if ($LASTEXITCODE -eq 0) { Pass "selftest" } else { Fail "selftest exit=$LASTEXITCODE" }

# ---- 3. junction a 2x2x4 all-present corner ----
function Remove-JunctionDir($d) {
  if (Test-Path $d) {
    Get-ChildItem $d -Directory -ErrorAction SilentlyContinue | ForEach-Object { cmd /c rmdir "$($_.FullName)" 2>$null }
    Remove-Item $d -Force -Recurse -ErrorAction SilentlyContinue
  }
}
Remove-JunctionDir $testDump
Remove-Item $testOut -Recurse -Force -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force $testDump | Out-Null
$n = 0
foreach ($z in 4352,4480) { foreach ($y in 2048,2176) { foreach ($x in 2048,2176,2304,2432) {
  $id = "z{0:D5}_y{1:D5}_x{2:D5}" -f $z,$y,$x
  $s  = Join-Path $SrcDump $id
  if (Test-Path $s) { cmd /c mklink /J (Join-Path $testDump $id) (Resolve-Path $s).Path | Out-Null; $n++ }
  else { Write-Host "  (skip missing $id)" }
}}}
Write-Host "== junctioned $n leaf cubes =="
if ($n -lt 12) { Fail "too few present cubes ($n) for a meaningful test" }

try {
  # ---- 4. run L1 + L2 ----
  Write-Host "== running pyramid (max-levels 2, keep $KeepRatio) =="
  & $exe $testDump $testOut --max-levels 2 --max-concurrent 4 --keep-ratio $KeepRatio
  if ($LASTEXITCODE -ne 0) { Fail "tool exit=$LASTEXITCODE" }

  function Parse-Manifold($log) {
    $rows = @()
    if (-not (Test-Path $log)) { return $rows }
    foreach ($ln in Get-Content $log) {
      if ($ln -match 'same_dir=(\d+)\s+VERDICT=(\w+)') {
        $rows += [pscustomobject]@{ same_dir = [int]$Matches[1]; verdict = $Matches[2] }
      }
    }
    return ,$rows
  }
  $l1 = Parse-Manifold "$testOut/level1/_manifold.log"
  $l2 = Parse-Manifold "$testOut/level2/_manifold.log"

  if ($l1.Count -eq 2) { Pass "L1 produced 2 nodes" } else { Fail "L1 node count = $($l1.Count), expected 2" }
  if (($l1 | Where-Object verdict -ne 'MANIFOLD').Count -eq 0 -and $l1.Count -gt 0) { Pass "L1 nodes all MANIFOLD" } else { Fail "L1 has non-manifold node(s)" }
  if ($l2.Count -eq 1) { Pass "L2 produced 1 node" } else { Fail "L2 node count = $($l2.Count), expected 1" }
  if (($l2 | Where-Object verdict -ne 'MANIFOLD').Count -eq 0 -and $l2.Count -gt 0) { Pass "L2 node MANIFOLD" } else { Fail "L2 node non-manifold" }

  $l1sd = ($l1 | Measure-Object same_dir -Sum).Sum
  $l2sd = ($l2 | Measure-Object same_dir -Sum).Sum
  if ($l2.Count -gt 0 -and $l2sd -lt $l1sd) { Pass "winding self-heal: L2 same_dir=$l2sd < L1 sum=$l1sd" }
  else { Fail "no winding self-heal (L1 sum=$l1sd, L2=$l2sd)" }
}
finally {
  Remove-JunctionDir $testDump
  Remove-Item $testOut -Recurse -Force -ErrorAction SilentlyContinue
}

Write-Host ""
if ($fails -eq 0) { Write-Host "hierarchical_weld TEST: ALL PASS" -ForegroundColor Green; exit 0 }
else { Write-Host "hierarchical_weld TEST: $fails FAILURE(S)" -ForegroundColor Red; exit 1 }
