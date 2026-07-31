<#
  wind_audit_regression.ps1 -- functional progress-tracker for scroll-mesh
  topology (inter-wrap fusions + intra-sheet splits).

  Analysis is done by the C tool build/Release/wind_audit.exe (winding-coordinate
  auditor); this script is only a RUNNER: it invokes wind_audit --json on a target
  mesh (and, by default, on the known-good QEM reference), then prints a
  side-by-side and gates PASS/FAIL on the reliable metrics so any pipeline change
  can be scored against the clean baseline.

  Reliable gate metrics (see src/tools/wind_audit.c):
    merge_ft_per_outer_edge  full-turn (|dw| in [0.7,1.3]) fusions / outer edge.
                             Clean QEM ref ~= 3e-6; coarse CVT weld ~= 1.5e-2.
    split_comp_pairs         # distinct component-pairs that are same-turn but
                             disconnected across a small gap (intra-sheet splits).
    split_pairs_seam_frac    fraction of split pairs sitting on a 128-vox cube
                             seam -> high => the WELD is failing to bridge.

  Usage:
    scripts/remesh/wind_audit_regression.ps1 -Mesh <top.obj> [-Baseline <ref.obj>]
        [-UmbY 3405] [-UmbX 2878] [-Pitch 9.5] [-OutDir output/wind_audit]
        [-MaxFtPerEdge 5e-4] [-MaxSplitCompPairs 40]
#>
param(
  [Parameter(Mandatory=$true)][string]$Mesh,
  [string]$Baseline = "output/rebuild_4x5x5/welded.obj",
  [double]$UmbY = 3405,
  [double]$UmbX = 2878,
  [double]$Pitch = 9.5,
  [string]$OutDir = "output/wind_audit",
  [double]$MaxFtPerEdge = 5e-4,       # 100x the clean baseline; well below coarse-weld 1.5e-2
  [int]$MaxSplitCompPairs = 40
)

$ErrorActionPreference = "Stop"
$exe = "build/Release/wind_audit.exe"
if (-not (Test-Path $exe)) { throw "wind_audit.exe not found at $exe -- build wind_audit.vcxproj first" }
if (-not (Test-Path $Mesh)) { throw "target mesh not found: $Mesh" }
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

function Run-Audit([string]$obj, [string]$tag) {
  $json = Join-Path $OutDir "$tag.json"
  $log  = Join-Path $OutDir "$tag.log"
  Write-Host "== wind_audit: $tag ($obj) =="
  & $exe $obj --umb-y $UmbY --umb-x $UmbX --pitch $Pitch --seam-pitch 128 --json $json --top 6 2>&1 |
    Tee-Object -FilePath $log | Select-String "MERGERS|FULL-TURN|histogram|SPLITS|total candidate" | ForEach-Object { "   $_" }
  if (-not (Test-Path $json)) { throw "wind_audit produced no json for $tag" }
  return (Get-Content $json -Raw | ConvertFrom-Json)
}

$T = Run-Audit $Mesh "target"
$B = $null
if (Test-Path $Baseline) { $B = Run-Audit $Baseline "baseline" }
else { Write-Host "WARN: baseline mesh not found ($Baseline) -- gating on absolute thresholds only" }

Write-Host ""
Write-Host "===================== SUMMARY ====================="
$fmt = "{0,-28} {1,16} {2,16}"
Write-Host ($fmt -f "metric", "target", "baseline(QEM)")
Write-Host ($fmt -f "components", $T.components, ($(if($B){$B.components}else{"-"})))
Write-Host ($fmt -f "full-turn fusions", $T.merge_ft, ($(if($B){$B.merge_ft}else{"-"})))
Write-Host ($fmt -f "  per outer edge", ("{0:e2}" -f $T.merge_ft_per_outer_edge), ($(if($B){"{0:e2}" -f $B.merge_ft_per_outer_edge}else{"-"})))
Write-Host ($fmt -f "  seam fraction", ("{0:p1}" -f $T.merge_ft_seam_frac), ($(if($B){"{0:p1}" -f $B.merge_ft_seam_frac}else{"-"})))
Write-Host ($fmt -f "split comp-pairs", $T.split_comp_pairs, ($(if($B){$B.split_comp_pairs}else{"-"})))
Write-Host ($fmt -f "  pairs @ seam", ("{0:p1}" -f $T.split_pairs_seam_frac), ($(if($B){"{0:p1}" -f $B.split_pairs_seam_frac}else{"-"})))
Write-Host "---------------------------------------------------"

$fail = @()
if ($T.merge_ft_per_outer_edge -gt $MaxFtPerEdge) { $fail += "FUSION: merge_ft_per_outer_edge $("{0:e2}" -f $T.merge_ft_per_outer_edge) > $("{0:e2}" -f $MaxFtPerEdge)" }
if ($T.split_comp_pairs -gt $MaxSplitCompPairs)   { $fail += "SPLIT: split_comp_pairs $($T.split_comp_pairs) > $MaxSplitCompPairs" }

if ($fail.Count -eq 0) {
  Write-Host "RESULT: PASS" -ForegroundColor Green
  exit 0
} else {
  foreach ($f in $fail) { Write-Host "RESULT: FAIL -- $f" -ForegroundColor Red }
  exit 1
}
