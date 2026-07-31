#requires -Version 5.1
<#
.SYNOPSIS
  Run the vesuvius-c pipeline on every cube in a grid extract and stitch
  the per-cube meshes into one global OBJ.

.DESCRIPTION
  Reads <GridDir>/manifest.json + <GridDir>/cubes_PRED/present.json,
  runs build/Release/cube_mesh.exe on each cube TIFF, then walks each
  cube's per-cube combined OBJ (default stage 'step12_final') and emits one
  global OBJ with vertices offset by the cube origin parsed from the
  filename (cube_id format: z{vz}_y{vy}_x{vx}, vertex order: Z Y X).

.PARAMETER GridDir
  Directory containing manifest.json and cubes_PRED/.
  Example: PHerc0139-4x5x5

.PARAMETER OutputDir
  Output root. Created if missing.
  Layout:
    <OutputDir>/cubes/<cube_id>.tif     (per-cube cleaned binary TIFF)
    <OutputDir>/obj/<cube_id>/...       (per-cube --dump-obj tree)
    <OutputDir>/logs/<cube_id>.log
    <OutputDir>/results.csv
    <OutputDir>/combined.obj            (the stitched mesh)

.PARAMETER ExePath
  Optional. Defaults to <repo>/build/Release/cube_mesh.exe.

.PARAMETER Stage
  Which dump (mesh) stage to stitch. Default 'step12_final' (the weld input).
  Valid: step1_bpa, step4_bpa, step5_cc, step7_cc_bpa, step8_holefill,
  step10_qem, step11_kibble, step12_final.

.PARAMETER AllStages
  Stitch every cached stage in one pass, writing combined_<stage>.obj for
  each. Iterates cubes once per stage; faster than running the script N
  times. Ignores -Stage.

.PARAMETER MaxCubes
  If >0, only process the first N cubes (testing).

.PARAMETER Only
  Optional list of cube_ids to process (overrides MaxCubes).

.PARAMETER CenterCrop
  Take only the center Nz x Ny x Nx cubes (e.g. "2x3x3"). Computed from
  the grid shape in manifest.json. Cheaper than -Only when you want a
  centered subset of a larger extract.

.PARAMETER OutputName
  Base name for the combined OBJ (default: 'combined'). Used to avoid
  clobbering an earlier full-grid stitch when stitching a subset.

.PARAMETER SkipRun
  Skip per-cube pipeline invocations. Use to re-stitch existing dumps.
  Combine with -AllStages to debug seam-gluing at every stage without
  rerunning the pipeline.

.PARAMETER SkipStitch
  Run pipeline only, no final stitch.

.EXAMPLE
  pwsh scripts/run_grid.ps1 -GridDir PHerc0139-4x5x5 -OutputDir output/PHerc0139-4x5x5

.EXAMPLE
  # Re-stitch all stages from existing dumps (no pipeline rerun)
  pwsh scripts/run_grid.ps1 -GridDir PHerc0139-4x5x5 -OutputDir output/PHerc0139-4x5x5 -SkipRun -AllStages
#>

[CmdletBinding()]
param(
    [Parameter(Mandatory=$true)][string]$GridDir,
    [Parameter(Mandatory=$true)][string]$OutputDir,
    [string]$ExePath,
    [ValidateSet('step1_bpa','step4_bpa','step5_cc','step7_cc_bpa','step8_holefill','step10_qem','step11_kibble','step12_final')]
    [string]$Stage = 'step12_final',
    [switch]$AllStages,
    [int]$MaxCubes = 0,
    [string[]]$Only,
    [string]$CenterCrop,
    [string]$OutputName = 'combined',
    [switch]$SkipRun,
    [switch]$SkipStitch
)

$ALL_STAGES = @('step1_bpa','step4_bpa','step5_cc','step7_cc_bpa','step8_holefill','step10_qem','step11_kibble','step12_final')

$ErrorActionPreference = 'Stop'
$script:WarningPreference = 'Continue'

# ============================================================
# Setup
# ============================================================

$repoRoot = (Resolve-Path "$PSScriptRoot/..").Path
$GridDir  = (Resolve-Path $GridDir).Path

if (-not $ExePath) {
    $ExePath = Join-Path $repoRoot 'build/Release/cube_mesh.exe'
}
if (-not (Test-Path $ExePath)) {
    throw "cube_mesh.exe not found at $ExePath. Build with: /vs build Release"
}

$manifestPath = Join-Path $GridDir 'manifest.json'
$presentPath  = Join-Path $GridDir 'cubes_PRED/present.json'
$cubesDir     = Join-Path $GridDir 'cubes_PRED'

foreach ($p in $manifestPath, $presentPath, $cubesDir) {
    if (-not (Test-Path $p)) { throw "Missing: $p" }
}

if (-not (Test-Path $OutputDir)) { New-Item -ItemType Directory -Path $OutputDir | Out-Null }
$OutputDir = (Resolve-Path $OutputDir).Path
$outCubes  = Join-Path $OutputDir 'cubes'
$outObj    = Join-Path $OutputDir 'obj'
$outLogs   = Join-Path $OutputDir 'logs'
foreach ($d in $outCubes, $outObj, $outLogs) {
    if (-not (Test-Path $d)) { New-Item -ItemType Directory -Path $d | Out-Null }
}

# ============================================================
# Load manifest + present list
# ============================================================

$manifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
$present  = Get-Content -LiteralPath $presentPath  -Raw | ConvertFrom-Json

$cubeIds = @($present.emitted)
if ($Only -and $Only.Count -gt 0) {
    $set = [System.Collections.Generic.HashSet[string]]::new([string[]]$Only)
    $cubeIds = @($cubeIds | Where-Object { $set.Contains($_) })
}

if ($CenterCrop) {
    if ($CenterCrop -notmatch '^(\d+)x(\d+)x(\d+)$') {
        throw "CenterCrop must be 'NzxNyxNx', got: $CenterCrop"
    }
    $cropZ = [int]$Matches[1]; $cropY = [int]$Matches[2]; $cropX = [int]$Matches[3]
    $gz, $gy, $gx = $manifest.n_chunks  # (Z,Y,X)
    foreach ($p in @(@($cropZ,$gz,'Z'), @($cropY,$gy,'Y'), @($cropX,$gx,'X'))) {
        if ($p[0] -gt $p[1]) { throw "CenterCrop $($p[2])=$($p[0]) exceeds grid $($p[2])=$($p[1])" }
    }
    # Unique sorted origins per axis from cube IDs
    $allOrigins = $present.emitted | ForEach-Object {
        if ($_ -match '^z(\d+)_y(\d+)_x(\d+)$') {
            [pscustomobject]@{ z=[int]$Matches[1]; y=[int]$Matches[2]; x=[int]$Matches[3] }
        }
    }
    $zSet = @($allOrigins.z | Sort-Object -Unique)
    $ySet = @($allOrigins.y | Sort-Object -Unique)
    $xSet = @($allOrigins.x | Sort-Object -Unique)
    $zStart = [Math]::Floor(($zSet.Count - $cropZ) / 2)
    $yStart = [Math]::Floor(($ySet.Count - $cropY) / 2)
    $xStart = [Math]::Floor(($xSet.Count - $cropX) / 2)
    $zKeep = $zSet[$zStart..($zStart + $cropZ - 1)]
    $yKeep = $ySet[$yStart..($yStart + $cropY - 1)]
    $xKeep = $xSet[$xStart..($xStart + $cropX - 1)]
    $zKeepSet = [System.Collections.Generic.HashSet[int]]::new([int[]]$zKeep)
    $yKeepSet = [System.Collections.Generic.HashSet[int]]::new([int[]]$yKeep)
    $xKeepSet = [System.Collections.Generic.HashSet[int]]::new([int[]]$xKeep)
    $cubeIds = @($cubeIds | Where-Object {
        if ($_ -match '^z(\d+)_y(\d+)_x(\d+)$') {
            $zKeepSet.Contains([int]$Matches[1]) -and
                $yKeepSet.Contains([int]$Matches[2]) -and
                $xKeepSet.Contains([int]$Matches[3])
        } else { $false }
    })
    Write-Host "CenterCrop ${cropZ}x${cropY}x${cropX}: Z=$($zKeep -join ',') Y=$($yKeep -join ',') X=$($xKeep -join ',')"
}

if ($MaxCubes -gt 0 -and $cubeIds.Count -gt $MaxCubes) {
    $cubeIds = $cubeIds[0..($MaxCubes - 1)]
}

Write-Host "Grid:        $GridDir"
Write-Host "Output:      $OutputDir"
Write-Host "Exe:         $ExePath"
Write-Host "Stage:       $Stage"
Write-Host "Cubes to do: $($cubeIds.Count) of $($present.emitted.Count) total"
Write-Host "Chunk size:  $($manifest.chunk_size)"
Write-Host "Grid shape:  $($manifest.n_chunks -join 'x') (Z,Y,X)"

# Copy manifest into output so the stitch is self-describing
Copy-Item -LiteralPath $manifestPath -Destination (Join-Path $OutputDir 'manifest_used.json') -Force

# ============================================================
# Helpers
# ============================================================

function Parse-CubeOrigin {
    param([string]$cubeId)
    # Format: z{vz:05d}_y{vy:05d}_x{vx:05d}
    if ($cubeId -match '^z(\d+)_y(\d+)_x(\d+)$') {
        return [pscustomobject]@{
            cube_id = $cubeId
            vz = [int]$Matches[1]
            vy = [int]$Matches[2]
            vx = [int]$Matches[3]
        }
    }
    throw "Cube ID does not match z#####_y#####_x##### : $cubeId"
}

function Get-StageObj {
    param([string]$cubeId, [string]$stage)
    $stageDir = Join-Path $outObj $cubeId
    $stageDir = Join-Path $stageDir "${cubeId}_$stage"
    return (Join-Path $stageDir "${cubeId}_${stage}_all.obj")
}

# ============================================================
# Run the pipeline on each cube
# ============================================================

$csvPath = Join-Path $OutputDir 'results.csv'
$csvSb   = [System.Text.StringBuilder]::new()
[void]$csvSb.AppendLine('cube_id,total_s,s0,s1,s2,s3,s4,s5,s6,rawsnap,s7,status,obj_path,n_verts,n_faces')

$tStart = Get-Date
$nOk = 0
$nFail = 0
$idx = 0

foreach ($cubeId in $cubeIds) {
    $idx++
    $inputTif  = Join-Path $cubesDir "$cubeId.tif"
    $outputTif = Join-Path $outCubes "$cubeId.tif"
    $logPath   = Join-Path $outLogs "$cubeId.log"
    $objPath   = Get-StageObj -cubeId $cubeId -stage $Stage

    Write-Host ("[{0,3}/{1}] {2}" -f $idx, $cubeIds.Count, $cubeId) -NoNewline

    if (-not (Test-Path $inputTif)) {
        Write-Host "  MISSING INPUT"
        [void]$csvSb.AppendLine("$cubeId,,,,,,,,,,,MISSING_INPUT,,")
        $nFail++
        continue
    }

    $tCube = Get-Date
    $status = 'OK'

    if (-not $SkipRun) {
        # Stream log file so we can tail in real time
        & $ExePath $inputTif $outputTif --dump-obj $outObj 2>$logPath
        $rc = $LASTEXITCODE
        if ($rc -ne 0) {
            $status = "FAILED_rc$rc"
            $nFail++
        }
    } elseif (-not (Test-Path $objPath)) {
        $status = 'MISSING_OBJ'
        $nFail++
    }

    # Parse step timings from log (whether we ran it now or it already existed)
    $s0=$s1=$s2=$s3=$s4=$s5=$s6=$s7=$rawsnap=$totalS = ''
    if (Test-Path $logPath) {
        $logText = Get-Content -LiteralPath $logPath -Raw
        if ($logText -match 'Step timings:\s*denoise=([\d.]+)\s+s0=([\d.]+)\s+s1=([\d.]+)\s+s2=([\d.]+)\s+s3=([\d.]+)\s+s4=([\d.]+)\s+s5=([\d.]+)\s+s6=([\d.]+)\s+s6\.5=([\d.]+)\s+rawsnap=([\d.]+)\s+s7=([\d.]+)') {
            $s0 = $Matches[2]; $s1 = $Matches[3]; $s2 = $Matches[4]
            $s3 = $Matches[5]; $s4 = $Matches[6]; $s5 = $Matches[7]
            $s6 = $Matches[8]; $rawsnap = $Matches[10]; $s7 = $Matches[11]
        }
        if ($logText -match 'Total:\s+([\d.]+)s\s+Status:\s+(OK|FAILED)') {
            $totalS = $Matches[1]
            if ($Matches[2] -eq 'FAILED' -and $status -eq 'OK') {
                $status = 'FAILED_pipeline'
                $nFail++
            }
        }
    }

    $nv = 0; $nf = 0
    if (Test-Path $objPath) {
        # Cheap counts via streaming
        $reader = [System.IO.File]::OpenText($objPath)
        try {
            while (-not $reader.EndOfStream) {
                $ln = $reader.ReadLine()
                if ($ln.Length -ge 2 -and $ln[0] -eq 'v' -and $ln[1] -eq ' ') { $nv++ }
                elseif ($ln.Length -ge 2 -and $ln[0] -eq 'f' -and $ln[1] -eq ' ') { $nf++ }
            }
        } finally { $reader.Dispose() }
    } else {
        if ($status -eq 'OK') {
            $status = 'NO_OBJ'
            $nFail++
        }
    }

    if ($status -eq 'OK') { $nOk++ }

    $dt = ((Get-Date) - $tCube).TotalSeconds
    Write-Host ("  {0,7:F2}s  {1}  nv={2} nf={3}" -f $dt, $status, $nv, $nf)

    $objRel = if (Test-Path $objPath) {
        [System.IO.Path]::GetRelativePath($OutputDir, $objPath)
    } else { '' }
    [void]$csvSb.AppendLine("$cubeId,$totalS,$s0,$s1,$s2,$s3,$s4,$s5,$s6,$rawsnap,$s7,$status,$objRel,$nv,$nf")
}

Set-Content -LiteralPath $csvPath -Value $csvSb.ToString() -Encoding ASCII

$wallSec = ((Get-Date) - $tStart).TotalSeconds
Write-Host ''
Write-Host '====================================='
Write-Host ('Pipeline run: {0} OK / {1} FAILED of {2}, {3:F1}s wall' -f $nOk, $nFail, $cubeIds.Count, $wallSec)
Write-Host "Results CSV: $csvPath"
Write-Host '====================================='

if ($SkipStitch) { return }

# ============================================================
# Stitch per-cube _all.obj into one global OBJ with origin offsets.
# Streaming: read each input line by line, rewrite v / f lines.
# ============================================================

function Invoke-StitchStage {
    param(
        [Parameter(Mandatory=$true)][string]$StageName,
        [Parameter(Mandatory=$true)][string]$OutputObj,
        [Parameter(Mandatory=$true)][string[]]$CubeIds,
        [Parameter(Mandatory=$true)]$Manifest,
        [Parameter(Mandatory=$true)][string]$GridDescr
    )

    Write-Host ''
    Write-Host "Stitching $StageName -> $OutputObj"

    $tStitch = Get-Date
    $writer = [System.IO.StreamWriter]::new($OutputObj, $false, [System.Text.UTF8Encoding]::new($false))
    $cubesIncluded = 0
    $missingCubes  = 0
    $totalV        = [long]0
    try {
        $writer.WriteLine("# vesuvius-c grid stitch")
        $writer.WriteLine("# grid: $GridDescr")
        $writer.WriteLine("# stage: $StageName")
        $writer.WriteLine("# vertex order: Z Y X  (source-space voxels)")
        $writer.WriteLine("# grid range: z=$($Manifest.z_range_l0 -join '..') y=$($Manifest.y_range_l0 -join '..') x=$($Manifest.x_range_l0 -join '..')")

        $vertOffset = [long]0
        $invariant = [System.Globalization.CultureInfo]::InvariantCulture

        foreach ($cubeId in $CubeIds) {
            $obj = Get-StageObj -cubeId $cubeId -stage $StageName
            if (-not (Test-Path $obj)) {
                $missingCubes++
                continue
            }
            $origin = Parse-CubeOrigin -cubeId $cubeId
            $oz = [double]$origin.vz
            $oy = [double]$origin.vy
            $ox = [double]$origin.vx

            $writer.WriteLine("# --- cube $cubeId  origin Z=$($origin.vz) Y=$($origin.vy) X=$($origin.vx) ---")
            $vertsInCube = [long]0

            $reader = [System.IO.StreamReader]::new($obj)
            try {
                while (-not $reader.EndOfStream) {
                    $ln = $reader.ReadLine()
                    if ($ln.Length -lt 2) { continue }
                    $c0 = $ln[0]; $c1 = $ln[1]

                    if ($c0 -eq 'v' -and $c1 -eq ' ') {
                        $parts = $ln.Split(' ', [StringSplitOptions]::RemoveEmptyEntries)
                        if ($parts.Length -lt 4) { continue }
                        $vz = [double]::Parse($parts[1], $invariant) + $oz
                        $vy = [double]::Parse($parts[2], $invariant) + $oy
                        $vx = [double]::Parse($parts[3], $invariant) + $ox
                        if ($parts.Length -ge 7) {
                            $writer.WriteLine([string]::Format($invariant,
                                "v {0:0.000000} {1:0.000000} {2:0.000000} {3} {4} {5}",
                                $vz, $vy, $vx, $parts[4], $parts[5], $parts[6]))
                        } else {
                            $writer.WriteLine([string]::Format($invariant,
                                "v {0:0.000000} {1:0.000000} {2:0.000000}",
                                $vz, $vy, $vx))
                        }
                        $vertsInCube++
                    }
                    elseif ($c0 -eq 'f' -and $c1 -eq ' ') {
                        $parts = $ln.Split(' ', [StringSplitOptions]::RemoveEmptyEntries)
                        if ($parts.Length -lt 4) { continue }
                        $idxs = New-Object int[] ($parts.Length - 1)
                        for ($i = 1; $i -lt $parts.Length; $i++) {
                            $tok = $parts[$i].Split('/')[0]
                            $idxs[$i-1] = [int]$tok + [int]$vertOffset
                        }
                        $writer.WriteLine("f " + ($idxs -join ' '))
                    }
                    elseif ($c0 -eq 'o' -and $c1 -eq ' ') {
                        $writer.WriteLine("o ${cubeId}_$($ln.Substring(2))")
                    }
                    # comments / blanks / unknown lines silently dropped
                }
            } finally { $reader.Dispose() }

            $vertOffset += $vertsInCube
            $totalV += $vertsInCube
            $cubesIncluded++
        }
    } finally { $writer.Dispose() }

    $dt = ((Get-Date) - $tStitch).TotalSeconds
    Write-Host ("  {0}: {1} cubes, {2} verts, {3:F1}s ({4} cubes missing this stage)" -f `
        $StageName, $cubesIncluded, $totalV, $dt, $missingCubes)
}

$gridDescr = $GridDir

if ($AllStages) {
    foreach ($s in $ALL_STAGES) {
        $obj = Join-Path $OutputDir "${OutputName}_$s.obj"
        Invoke-StitchStage -StageName $s -OutputObj $obj `
            -CubeIds $cubeIds -Manifest $manifest -GridDescr $gridDescr
    }
} else {
    $obj = Join-Path $OutputDir "$OutputName.obj"
    Invoke-StitchStage -StageName $Stage -OutputObj $obj `
        -CubeIds $cubeIds -Manifest $manifest -GridDescr $gridDescr
}
