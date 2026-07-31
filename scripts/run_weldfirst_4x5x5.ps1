# run_weldfirst_4x5x5.ps1 -- "weld dense, then RVD-decimate" pipeline.
#
# Reorders the mesh pipeline so welding is decoupled from coarseness: mesh cubes
# DENSE (--no-qem) so their ~1-vox seams bridge reliably, weld hierarchically to
# bound memory, then RVD-remesh the SEAMLESS welded result per wrap-component as
# coarse as wanted. Because the remesh runs AFTER welding, there is no seam-density
# constraint and the interior is freely light; because it runs per-component it
# cannot fuse wraps. Validated on PHerc0139-4x5x5: complete weld (seam_audit clean),
# manifold, ~107k faces at DECIM=10 (vs graded 1.32M / uniform-0.012 904k).
#
#   pwsh scripts/run_weldfirst_4x5x5.ps1 [-Grid PHerc0139-4x5x5] [-Out output/cvt_weldfirst_4x5x5]
#        [-Decim 10]     # per-component target = component_faces / (2*Decim) sites; smaller => finer
#        [-Conc 16] [-SkipLeaves] [-SkipLod]
param(
  [string]$Grid = "PHerc0139-4x5x5",
  [string]$Out  = "output\cvt_weldfirst_4x5x5",
  [int]$Decim   = 10,
  [int]$Conc    = 16,
  # dir holding cvt_remesh_obj.exe, mesh_clean.exe, mesh_quality.exe. These three
  # currently build via scratchpad cl scripts (TODO: wire vcxprojs); point -Tools at
  # wherever they are (default assumes they've been copied into build\Release).
  [string]$Tools = "build\Release",
  [switch]$SkipLeaves,
  [switch]$SkipLod
)
$ErrorActionPreference = "Stop"
$B = "build\Release"
$S = $Tools
New-Item -ItemType Directory -Force $Out | Out-Null
Remove-Item Env:\VES_CVT_RATIO -ErrorAction SilentlyContinue
function Log($m){ Write-Host ("[{0}] {1}" -f (Get-Date -Format HH:mm:ss), $m) }

# 1. DENSE trimmed leaves (no per-cube simplify) -- seams stay ~1 vox, weld reliably.
if (-not $SkipLeaves) {
  Log "dense leaves (--no-qem)"
  & "$B\grid_pipeline.exe" $Grid $Out --no-qem --halo 13 --skip-weld `
      --max-concurrent $Conc --threads-per-cube 1 *> "$Out\leaves.err"
}

# 2. Memory-bounded seamless weld: hierarchical weld + qslim (pins open boundaries so
#    each level still welds; winding gate armed by PHerc0139 umbilicus/pitch defaults).
if (-not $SkipLod) {
  Log "hierarchical weld + decimate"
  & "$B\hierarchical_weld.exe" "$Out\dump" "$Out\lod" --leaf-stage step12_final `
      --fanout 2 --keep-ratio 0.25 --max-concurrent $Conc --no-remesh *> "$Out\lod.err"
}
$top = Get-ChildItem "$Out\lod\level*" -Recurse -Filter *_all.obj |
       Sort-Object { [int]($_.Directory.Parent.Name -replace '\D','0') } |
       Select-Object -Last 1
if (-not $top) { throw "no LOD top OBJ found under $Out\lod" }
Log ("LOD top: {0}" -f $top.FullName)

# 3-5. Split the seamless top into wrap-components, RVD-remesh each coarse (per-component
#      => no fusion; standalone tool => no size gate), merge.
Log "split + per-component CVT remesh (Decim=$Decim)"
New-Item -ItemType Directory -Force "$Out\topcomps","$Out\topcvt" | Out-Null
Get-ChildItem "$Out\topcomps\c_*.obj" -ErrorAction SilentlyContinue | Remove-Item
& "$B\obj_components.exe" $top.FullName --dump-sep "$Out\topcomps\c" --brief *> "$Out\components.log"
$cvtOut = @()
foreach ($f in Get-ChildItem "$Out\topcomps\c_*.obj") {
  if ($f.Name -match '_(\d+)f\.obj$') { $nf = [int]$Matches[1] } else { continue }
  if ($nf -lt 200) { continue }                          # drop tiny fragments
  $tv = [Math]::Max(50, [int]($nf / $Decim))
  $o  = "$Out\topcvt\$($f.BaseName)_cvt.obj"
  & "$S\cvt_remesh_obj.exe" $f.FullName $o $tv --iters 10 *> $null
  if (Test-Path $o) { $cvtOut += $o }
}
& "$B\obj_merge.exe" "$Out\weldfirst_merged.obj" @cvtOut *> $null

# 6. Manifold cleanup (split the few CVT-dual bowties, re-assert winding).
Log "mesh_clean (--reorient)"
& "$S\mesh_clean.exe" "$Out\weldfirst_merged.obj" "$Out\weldfirst_final.obj" --reorient

# Report.
Log "validate"
& "$B\manifold_check.exe" "$Out\weldfirst_final.obj" 2>&1 | Select-Object -Last 4
& "$B\seam_audit.exe"     "$Out\weldfirst_final.obj" --cube 128 2>&1 | Select-Object -Last 3
& "$S\mesh_quality.exe"   "$Out\weldfirst_final.obj" 2>&1 | Select-Object -First 6
& "$B\mesh_render.exe"    "$Out\weldfirst_final.obj" "$Out\weldfirst_final.png" --views --size 800 --no-edges *> $null
Log ("DONE -> {0}\weldfirst_final.obj" -f $Out)
