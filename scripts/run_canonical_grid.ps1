# run_canonical_grid.ps1 -- produce a canonical worked run for one cube grid.
#
# Runs the whole-scroll path end to end with the validated ribbon settings and
# publishes review images for the CT-snapped and final-relaxed surfaces before
# baking the finished ribbon:
#
#   grid_pipeline -> scroll_whole -> --reregister -> --audit
#                 -> scroll_unroll (join/ownership/snap/final relax reviews)
#                 -> atlas_overlap_fix --production -> atlas texture
#                 -> atlas_ribbon_fit -> ribbon texture
#
# Everything runs out of this repository's own build/Release, so what is
# validated here is what the repository ships.
#
#   pwsh scripts/run_canonical_grid.ps1 -Grid ..\vesuvius-c\PHerc0139-4x5x5 `
#        -Out output\pherc0139_4x5x5_run
param(
  [Parameter(Mandatory)][string]$Grid,
  [Parameter(Mandatory)][string]$Out,
  [string]$Tree   = "D:\work\scrollfiesta_public",
  # Optional real CT source override. Some archived test grids retain valid
  # predictions beside placeholder/zero RAW cubes; never bake from those just
  # because they share the grid directory.
  [string]$RawDir = "",
  [int]$Conc      = 16,
  [int]$Tpc       = 2,
  # Zero means use the validated adaptive production preset. Positive values
  # are explicit diagnostic overrides.
  [int]$Rounds    = 0,
  [int]$TabuSpan  = 0,
  # Radial pitch, vox/turn. Left EMPTY on purpose: scroll_whole estimates it
  # from the geometry, and the estimate is what produced our best atlas. Pin it
  # only when you know the scroll's pitch and the estimate is wrong -- a wrong
  # pin mis-measures the distance between wraps and fragments the atlas.
  [string]$WrapSpacing = "",
  # Geometry for grid_weld's winding/phase gates. These values are deliberately
  # weld-scoped: they must not alter per-cube BPA or fragment coherent sheets.
  # All three must be set together; left empty, the weld gates stay disarmed.
  # PHerc0139: -UmbY 3405 -UmbX 2878 -WrapPitch 9.5
  [string]$UmbY = "",
  [string]$UmbX = "",
  [string]$WrapPitch = "",
  # Production meshing uses the GPU MLS implementation. The backend is passed
  # explicitly so a stale parent shell cannot silently select Rust CPU/auto.
  [ValidateSet("cubecl-cuda", "rust-cpu", "cubecl-cpu", "cubecl-hip")]
  [string]$MlsBackend = "cubecl-cuda",
  # CUDA_PATH is needed by NVRTC for headers. Prefer an explicit path, then an
  # existing environment value, then this repository's local wheel toolkit.
  [string]$CudaPath = "",
  # The atlas path consumes per-cube dumps directly; a flat whole-grid weld is
  # wasted work here and is memory-infeasible on grids such as 4x21x21.
  # Hierarchical LOD welding is a separate workflow. Keep the old SkipWeld
  # switch for compatibility and require an explicit opt-in for a flat weld.
  [switch]$SkipWeld,
  [switch]$RunFlatWeld,
  [switch]$SkipMesh,
  [switch]$AtlasOnly,
  # Grids wider than five cubes on any axis use bounded Gauss-Seidel tiles by
  # default. This switch is an explicit diagnostic escape hatch only.
  [switch]$MonolithicAtlas,
  [int]$AtlasCoreSize = 3,
  [int]$AtlasHalo = 1,
  [int]$AtlasSweeps = 1,
  # Four voxels along the very long u axis keeps a 4x21x21 atlas under the
  # baker's 2^26-pixel safety budget while retaining three-voxel v detail.
  [double]$AtlasRasterDu = 4.0,
  [double]$AtlasRasterDv = 3.0
)
$ErrorActionPreference = "Stop"
$B    = Join-Path $Tree "build\Release"
$logs = Join-Path $Out "logs"
$Raw = if ($RawDir -ne "") {
  (Resolve-Path -LiteralPath $RawDir).Path
} else {
  (Resolve-Path -LiteralPath (Join-Path $Grid "cubes_RAW")).Path
}
New-Item -ItemType Directory -Force $Out, $logs | Out-Null
$summary = Join-Path $logs "SUMMARY.txt"
Set-Content -Path $summary -Value "canonical run -- grid=$Grid tree=$Tree started $(Get-Date -Format s)"

function Log($m) { Write-Host ("[{0}] {1}" -f (Get-Date -Format HH:mm:ss), $m) }

# A nonzero exit from --audit is a VERDICT about the atlas, not a crash: record
# it and carry on, so a bad grid still produces its (honest) texture.
function Stage($name, $exe, $argList, [switch]$allowFail) {
  $log = Join-Path $logs "$name.log"
  Log ("  {0,-18} -> {1}" -f $name, (Split-Path $log -Leaf))
  $sw = [Diagnostics.Stopwatch]::StartNew()
  & $exe @argList *> $log
  $code = $LASTEXITCODE
  $sw.Stop()
  Add-Content -Path $summary -Value ("{0}: exit={1} t={2:N1}s" -f $name, $code, $sw.Elapsed.TotalSeconds)
  if ($code -ne 0 -and -not $allowFail) { throw "STAGE FAILED -- $name exit=$code (see $log)" }
  Log ("  {0,-18}    exit={1} t={2:N1}s" -f $name, $code, $sw.Elapsed.TotalSeconds)
}

function RequireArtifact($label, $path, $notOlderThan = $null) {
  if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
    throw "CANONICAL ARTIFACT MISSING -- $label ($path)"
  }
  $item = Get-Item -LiteralPath $path
  if ($item.Length -le 0) {
    throw "CANONICAL ARTIFACT EMPTY -- $label ($path)"
  }
  if ($null -ne $notOlderThan -and $item.LastWriteTimeUtc -lt $notOlderThan) {
    throw "CANONICAL ARTIFACT STALE -- $label ($path)"
  }
}

if (-not $SkipMesh) {
  $env:MLS_BACKEND = $MlsBackend
  if ($MlsBackend -eq "cubecl-cuda") {
    $cudaCandidates = @()
    if ($CudaPath -ne "") { $cudaCandidates += $CudaPath }
    if ($env:CUDA_PATH) { $cudaCandidates += $env:CUDA_PATH }
    $cudaCandidates += (Join-Path $Tree ".toolchains\cuda-13.2-wheels\nvidia\cu13")
    $cuda = $cudaCandidates | Where-Object {
      Test-Path (Join-Path $_ "include\cuda.h")
    } | Select-Object -First 1
    if (-not $cuda) {
      throw "CubeCL CUDA requested, but no CUDA_PATH with include\cuda.h was found. Pass -CudaPath explicitly."
    }
    $cuda = (Resolve-Path $cuda).Path
    $env:CUDA_PATH = $cuda
    $cudaBin = Join-Path $cuda "bin"
    $cudaArchBin = Join-Path $cudaBin "x86_64"
    $env:Path = "$B;$cudaArchBin;$cudaBin;$($env:Path)"
    foreach ($runtime in @("herculaneum_mls_cubecl.dll", "cudart64_13.dll", "nvrtc64_130_0.dll")) {
      if (-not (Test-Path (Join-Path $B $runtime))) {
        throw "CubeCL CUDA runtime is incomplete: missing $B\$runtime"
      }
    }
    Add-Content -Path $summary -Value "mls_backend=$MlsBackend cuda_path=$cuda"
    Log "Cube MLS backend: $MlsBackend ($cuda)"
  } else {
    Add-Content -Path $summary -Value "mls_backend=$MlsBackend"
    Log "Cube MLS backend: $MlsBackend"
  }
}

if (-not $SkipMesh) {
  $gp = @($Grid, "$Out\mesh", "--halo", "13",
          "--trim-inset", "0",
          "--max-concurrent", "$Conc", "--threads-per-cube", "$Tpc",
          "--exe", "$B\cube_mesh.exe", "--weld", "$B\grid_weld.exe")
  if ($UmbY -ne "" -and $UmbX -ne "" -and $WrapPitch -ne "") {
    $gp += @("--umb-y", $UmbY, "--umb-x", $UmbX, "--wrap-pitch", $WrapPitch)
  }
  if ($SkipWeld -or -not $RunFlatWeld) {
    $gp += "--skip-weld"
    Log "Flat grid weld skipped (atlas consumes per-cube dumps)"
  }
  Stage "grid_pipeline" "$B\grid_pipeline.exe" $gp
}

$sw = @("$Out\mesh\dump", "$Out\placed",
        "--axis-point", "0", "3405", "2878",
        "--pair-gate", "3.5", "--skin", "4.0", "--max-concurrent", "$Conc")
if ($WrapSpacing -ne "") { $sw += @("--wrap-spacing", $WrapSpacing) }
Stage "scroll_whole" "$B\scroll_whole.exe" $sw

Stage "reregister" "$B\scroll_whole.exe" @("$Out\placed", "--reregister")
Stage "audit"      "$B\scroll_whole.exe" @("$Out\placed", "--audit") -allowFail

$coords = @(Get-ChildItem -LiteralPath "$Out\placed" -Filter '*_uvphi.f32' -File | ForEach-Object {
  if ($_.Name -match '^z(\d+)_y(\d+)_x(\d+)_uvphi\.f32$') {
    [pscustomobject]@{ Z=$Matches[1]; Y=$Matches[2]; X=$Matches[3] }
  }
})
if ($coords.Count -eq 0) { throw "No placed UV sidecars under $Out\placed" }
$zCount = @($coords.Z | Sort-Object -Unique).Count
$yCount = @($coords.Y | Sort-Object -Unique).Count
$xCount = @($coords.X | Sort-Object -Unique).Count
$useTiledAtlas = !$MonolithicAtlas -and ($zCount -gt 5 -or $yCount -gt 5 -or $xCount -gt 5)
$atlasSolution = "$Out\atlas\atlas_solution.bin"

# The atlas and ribbon paths do not perform the CT surface snap or its final
# light UV relaxation. Run the canonical five-stage review path as a sibling
# product so those two states remain directly inspectable. Coarser raster
# spacing keeps the review bounded on the long 4x21x21 grid; it changes only
# the review raster, not the snapped/relaxed geometry.
$reviewDir = "$Out\snap_relax"
$reviewId = "canonical"
$reviewDu = $AtlasRasterDu.ToString("R", [Globalization.CultureInfo]::InvariantCulture)
$reviewDv = $AtlasRasterDv.ToString("R", [Globalization.CultureInfo]::InvariantCulture)
$reviewStartedUtc = [DateTime]::UtcNow.AddSeconds(-2)
Stage "snap_relax" "$B\scroll_unroll.exe" @(
  "$Out\placed", $reviewDir, "--raw", $Raw,
  "--steps", "12345", "--id", $reviewId,
  "--du", $reviewDu, "--dv", $reviewDv, "--threads", "$Conc")

$snapPreview = Join-Path $reviewDir "${reviewId}_step4_snap_rawtex_preview.png"
$snapStrip = Join-Path $reviewDir "${reviewId}_step4_snap_rawtex_strip.png"
$snapTif = Join-Path $reviewDir "${reviewId}_step4_snap_rawtex.tif"
$snapStats = Join-Path $reviewDir "${reviewId}_step4_snap_stats.json"
$relaxPreview = Join-Path $reviewDir "${reviewId}_step5_relax_rawtex_preview.png"
$relaxStrip = Join-Path $reviewDir "${reviewId}_step5_relax_rawtex_strip.png"
$relaxTif = Join-Path $reviewDir "${reviewId}_step5_relax_rawtex.tif"
$relaxStats = Join-Path $reviewDir "${reviewId}_step5_relax_stats.json"
RequireArtifact "snapped CT preview" $snapPreview $reviewStartedUtc
RequireArtifact "snapped CT readable strip" $snapStrip $reviewStartedUtc
RequireArtifact "snapped CT full-resolution raster" $snapTif $reviewStartedUtc
RequireArtifact "snapped CT metrics" $snapStats $reviewStartedUtc
RequireArtifact "final-relaxed CT preview" $relaxPreview $reviewStartedUtc
RequireArtifact "final-relaxed CT readable strip" $relaxStrip $reviewStartedUtc
RequireArtifact "final-relaxed CT full-resolution raster" $relaxTif $reviewStartedUtc
RequireArtifact "final-relaxed metrics" $relaxStats $reviewStartedUtc
try {
  $snapReport = Get-Content -LiteralPath $snapStats -Raw | ConvertFrom-Json
  $relaxReport = Get-Content -LiteralPath $relaxStats -Raw | ConvertFrom-Json
} catch {
  throw "CANONICAL ARTIFACT INVALID -- snap/relax metrics are not valid JSON: $_"
}
if ($snapReport.tool -ne "scroll_unroll step4 snap_grid") {
  throw "CANONICAL ARTIFACT INVALID -- unexpected snap metrics tool '$($snapReport.tool)'"
}
if ($relaxReport.tool -ne "scroll_unroll step5 relax_grid") {
  throw "CANONICAL ARTIFACT INVALID -- unexpected relax metrics tool '$($relaxReport.tool)'"
}
$reviewIndex = Join-Path $Out "CANONICAL_OUTPUTS.md"
@"
# Canonical review outputs

| State | Preview | Readable strip | Full resolution | Metrics |
|---|---|---|---|---|
| CT-snapped surface | [preview](snap_relax/${reviewId}_step4_snap_rawtex_preview.png) | [PNG](snap_relax/${reviewId}_step4_snap_rawtex_strip.png) | [TIF](snap_relax/${reviewId}_step4_snap_rawtex.tif) | [JSON](snap_relax/${reviewId}_step4_snap_stats.json) |
| Final-relaxed surface | [preview](snap_relax/${reviewId}_step5_relax_rawtex_preview.png) | [PNG](snap_relax/${reviewId}_step5_relax_rawtex_strip.png) | [TIF](snap_relax/${reviewId}_step5_relax_rawtex.tif) | [JSON](snap_relax/${reviewId}_step5_relax_stats.json) |
"@ | Set-Content -LiteralPath $reviewIndex -Encoding utf8
RequireArtifact "canonical review index" $reviewIndex
Add-Content -Path $summary -Value "snap_preview=$snapPreview"
Add-Content -Path $summary -Value "snap_strip=$snapStrip"
Add-Content -Path $summary -Value "snap_tif=$snapTif"
Add-Content -Path $summary -Value "snap_stats=$snapStats"
Add-Content -Path $summary -Value "relax_preview=$relaxPreview"
Add-Content -Path $summary -Value "relax_strip=$relaxStrip"
Add-Content -Path $summary -Value "relax_tif=$relaxTif"
Add-Content -Path $summary -Value "relax_stats=$relaxStats"
Add-Content -Path $summary -Value "review_index=$reviewIndex"
Log "  snapped review     -> $snapPreview"
Log "  relaxed review     -> $relaxPreview"
Log "  review index       -> $reviewIndex"

if ($useTiledAtlas) {
  Log "Atlas mode: bounded tiles for ${zCount}x${yCount}x${xCount} grid"
  $tiledOut = "$Out\tiled_atlas"
  Stage "tiled_atlas" "C:\Program Files\PowerShell\7\pwsh.exe" @(
    "-NoProfile", "-File", (Join-Path $Tree "scripts\run_tiled_atlas.ps1"),
    "-Placed", "$Out\placed", "-Out", $tiledOut,
    "-CoreSize", "$AtlasCoreSize", "-Halo", "$AtlasHalo",
    "-Sweeps", "$AtlasSweeps", "-WindSeed",
    "-Raw", $Raw,
    "-RasterDu", "$AtlasRasterDu", "-RasterDv", "$AtlasRasterDv")
  $atlasSolution = "$tiledOut\atlas\atlas_solution.bin"
} else {
  $atlasArgs = @("$Out\placed", "$Out\atlas", "--production")
  if ($Rounds -gt 0)   { $atlasArgs += @("--rounds", "$Rounds") }
  if ($TabuSpan -gt 0) { $atlasArgs += @("--tabu-span", "$TabuSpan") }
  Stage "atlas_overlap_fix" "$B\atlas_overlap_fix.exe" $atlasArgs

  # Always make the atlas inspectable before the ribbon is allowed to obscure
  # where a failure entered the pipeline.
  Stage "obj_bake_atlas" "$B\obj_bake_raw.exe" @(
    "$Out\atlas\atlas_bake.obj", $Raw,
    "$Out\atlas_texture", "--id", "atlas", "--require-contrast")
}
if ($AtlasOnly) {
  Log "=== atlas-only complete ==="
  Get-Content $summary
  exit 0
}

Stage "atlas_ribbon_fit"  "$B\atlas_ribbon_fit.exe"  @(
  "$Out\placed", $atlasSolution, "$Out\ribbon",
  "--collision-rounds", "24", "--collision-polish", "4",
  "--collision-sweeps", "20", "--register-sweeps", "20",
  "--lambda-register-v", "1024")

# One fresh volume sample per texel through the finished collision-registered
# ribbon. obj_bake_raw writes both the full-length texture and row-cut strip.
Stage "obj_bake_raw" "$B\obj_bake_raw.exe" @(
  "$Out\ribbon\ribbon.obj", $Raw, "$Out\texture", "--id", "ribbon",
  "--require-contrast")

Log "=== complete ==="
Get-Content $summary
