[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$Placed,
    [Parameter(Mandatory = $true)][string]$Out,
    [int]$CoreSize = 3,
    [int]$Halo = 1,
    [int]$Sweeps = 1,
    [ValidateSet('z','y','x')][string[]]$ForceSplitAxes = @(),
    [switch]$WindSeed,
    [switch]$SkipAudit,
    [string]$Raw = '',
    [string]$AtlasExe = '',
    [double]$RasterDu = 1.0,
    [double]$RasterDv = 1.0
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

if ($CoreSize -lt 1 -or $CoreSize + 2 * $Halo -gt 5) {
    throw "CoreSize + 2*Halo must be <= 5 (got $CoreSize + 2*$Halo)"
}
if ($Sweeps -lt 1) { throw 'Sweeps must be positive' }

$repo = Split-Path -Parent $PSScriptRoot
if ($AtlasExe -eq '') {
    $AtlasExe = Join-Path $repo 'build\Release\atlas_overlap_fix.exe'
}
$AtlasExe = (Resolve-Path -LiteralPath $AtlasExe).Path
$Placed = (Resolve-Path -LiteralPath $Placed).Path
$Out = [IO.Path]::GetFullPath($Out)
if (Test-Path -LiteralPath $Out) {
    throw "Output already exists; refusing to overwrite: $Out"
}
New-Item -ItemType Directory -Path $Out | Out-Null

function New-HardLink([string]$Path, [string]$Target) {
    try {
        New-Item -ItemType HardLink -Path $Path -Target $Target -ErrorAction Stop | Out-Null
    } catch {
        # Sandboxed validation may allow reading a reference tree without
        # granting the filesystem right needed to hard-link to it. A copy is
        # semantically identical; normal in-workspace production runs retain
        # the cheap hard-link path.
        if (!(Test-Path -LiteralPath $Target) -or (Test-Path -LiteralPath $Path)) {
            throw
        }
        Copy-Item -LiteralPath $Target -Destination $Path
    }
}

function Axis-Step([long[]]$Values) {
    if ($Values.Count -lt 2) { return 128L }
    $groups = @(for ($i = 1; $i -lt $Values.Count; $i++) {
        $Values[$i] - $Values[$i - 1]
    }) | Group-Object | Sort-Object Count -Descending
    return [long]$groups[0].Name
}

function Axis-Ranges([long[]]$Values, [string]$Axis) {
    $split = $Values.Count -gt 5 -or $ForceSplitAxes -contains $Axis
    $stride = if ($split) { $CoreSize } else { $Values.Count }
    $ranges = @()
    for ($first = 0; $first -lt $Values.Count; $first += $stride) {
        $last = [Math]::Min($Values.Count, $first + $stride)
        $ranges += [pscustomobject]@{ First = $first; Last = $last }
    }
    return $ranges
}

$cube = @()
foreach ($file in Get-ChildItem -LiteralPath $Placed -Filter '*_uvphi.f32' -File) {
    if ($file.Name -notmatch '^(z(\d+)_y(\d+)_x(\d+))_uvphi\.f32$') { continue }
    $cube += [pscustomobject]@{
        Id = $Matches[1]
        Z = [long]$Matches[2]
        Y = [long]$Matches[3]
        X = [long]$Matches[4]
    }
}
if ($cube.Count -eq 0) { throw "No placed UV sidecars under $Placed" }

[long[]]$zv = @($cube.Z | Sort-Object -Unique)
[long[]]$yv = @($cube.Y | Sort-Object -Unique)
[long[]]$xv = @($cube.X | Sort-Object -Unique)
$zstep = Axis-Step $zv; $ystep = Axis-Step $yv; $xstep = Axis-Step $xv
$zindex = @{}; $yindex = @{}; $xindex = @{}
for ($i = 0; $i -lt $zv.Count; $i++) { $zindex[$zv[$i]] = $i }
for ($i = 0; $i -lt $yv.Count; $i++) { $yindex[$yv[$i]] = $i }
for ($i = 0; $i -lt $xv.Count; $i++) { $xindex[$xv[$i]] = $i }
foreach ($c in $cube) {
    Add-Member -InputObject $c -NotePropertyName IZ -NotePropertyValue $zindex[$c.Z]
    Add-Member -InputObject $c -NotePropertyName IY -NotePropertyValue $yindex[$c.Y]
    Add-Member -InputObject $c -NotePropertyName IX -NotePropertyValue $xindex[$c.X]
}

$zr = @(Axis-Ranges $zv 'z'); $yr = @(Axis-Ranges $yv 'y'); $xr = @(Axis-Ranges $xv 'x')
$tiles = @()
foreach ($rz in $zr) { foreach ($ry in $yr) { foreach ($rx in $xr) {
    $tiles += [pscustomobject]@{ Z = $rz; Y = $ry; X = $rx }
}}}
$singleFullGridTile = $tiles.Count -eq 1
Write-Host (("[tiled_atlas] grid={0}x{1}x{2}, cubes={3}, tiles={4}, " +
            "core={5}, halo={6}, sweeps={7}") -f
            $zv.Count, $yv.Count, $xv.Count, $cube.Count, $tiles.Count,
            $CoreSize, $Halo, $Sweeps)

$seedUv = $Placed
$seedStats = $null
if ($WindSeed) {
    $seedOut = Join-Path $Out 'global_wind_seed'
    Write-Host '[tiled_atlas] global wind + ladder calibration seed'
    # A no-move global production pass performs the same winding
    # reconciliation as the historical large solve and also measures one
    # signed field ladder for every tile. A corner tile can legitimately have
    # zero calibration pairs (spiral-map fallback), so it must never define
    # the shared ladder on its own.
    & $AtlasExe $Placed $seedOut --production --no-relayout `
        --tabu-iters 0 --tabu-stall 1 --write-uvphi
    if ($LASTEXITCODE -ne 0) { throw "global wind/calibration seed failed: $LASTEXITCODE" }
    $seedUv = Join-Path $seedOut 'solved_uvphi'
    $seedStats = Get-Content -LiteralPath (Join-Path $seedOut 'overlap_fix_stats.json') -Raw | ConvertFrom-Json
}

$state = Join-Path $Out 'state_uvphi'
New-Item -ItemType Directory -Path $state | Out-Null
foreach ($c in $cube) {
    $name = "$($c.Id)_uvphi.f32"
    Copy-Item -LiteralPath (Join-Path $seedUv $name) -Destination (Join-Path $state $name)
}

$placedIndex = Join-Path $Placed 'placed_index.json'
if (!(Test-Path -LiteralPath $placedIndex)) { throw "Missing $placedIndex" }
$ladderStep = if ($null -ne $seedStats -and
    [Math]::Abs([double]$seedStats.tabu_ladder_step) -gt 1.0e-12) {
    [double]$seedStats.tabu_ladder_step
} else { $null }
$ladderU0 = if ($null -ne $ladderStep) { [double]$seedStats.tabu_ladder_u0 } else { $null }
if ($null -ne $ladderStep) {
    Write-Host ("[tiled_atlas] global calibration step={0:R}, u0={1:R}" -f $ladderStep,$ladderU0)
}
$singleTileAtlas = $null
$summary = [Collections.Generic.List[object]]::new()

for ($sweep = 0; $sweep -lt $Sweeps; $sweep++) {
    $ordered = if (($sweep % 2) -eq 0) { @($tiles) } else { @($tiles[($tiles.Count-1)..0]) }
    $tileNo = 0
    foreach ($tile in $ordered) {
        $tileNo++
        $zl = [Math]::Max(0, $tile.Z.First - $Halo)
        $zh = [Math]::Min($zv.Count, $tile.Z.Last + $Halo)
        $yl = [Math]::Max(0, $tile.Y.First - $Halo)
        $yh = [Math]::Min($yv.Count, $tile.Y.Last + $Halo)
        $xl = [Math]::Max(0, $tile.X.First - $Halo)
        $xh = [Math]::Min($xv.Count, $tile.X.Last + $Halo)
        if ($zh-$zl -gt 5 -or $yh-$yl -gt 5 -or $xh-$xl -gt 5) {
            throw 'Internal error: loaded tile exceeds 5 cubes on an axis'
        }

        $tag = 'z{0:D2}_y{1:D2}_x{2:D2}' -f $tile.Z.First,$tile.Y.First,$tile.X.First
        $root = Join-Path $Out ('tiles\sweep_{0:D2}\{1}' -f ($sweep+1),$tag)
        $input = Join-Path $root 'input'
        $atlas = Join-Path $root 'atlas'
        New-Item -ItemType Directory -Path $input -Force | Out-Null
        New-HardLink (Join-Path $input 'placed_index.json') $placedIndex

        $loaded = @($cube | Where-Object {
            $_.IZ -ge $zl -and $_.IZ -lt $zh -and
            $_.IY -ge $yl -and $_.IY -lt $yh -and
            $_.IX -ge $xl -and $_.IX -lt $xh
        })
        foreach ($c in $loaded) {
            foreach ($file in Get-ChildItem -LiteralPath $Placed -Filter "$($c.Id)_*" -File) {
                $target = Join-Path $input $file.Name
                if ($file.Name -like '*_uvphi.f32') {
                    New-HardLink $target (Join-Path $state $file.Name)
                } else {
                    New-HardLink $target $file.FullName
                }
            }
        }

        $z0 = $zv[$tile.Z.First]
        $z1 = if ($tile.Z.Last -lt $zv.Count) { $zv[$tile.Z.Last] } else { $zv[-1] + $zstep }
        $y0 = $yv[$tile.Y.First]
        $y1 = if ($tile.Y.Last -lt $yv.Count) { $yv[$tile.Y.Last] } else { $yv[-1] + $ystep }
        $x0 = $xv[$tile.X.First]
        $x1 = if ($tile.X.Last -lt $xv.Count) { $xv[$tile.X.Last] } else { $xv[-1] + $xstep }
        $args = @($input,$atlas,'--production')
        if (!$singleFullGridTile) {
            $args += @('--active-core',
                       "$z0","$z1","$y0","$y1","$x0","$x1")
        }
        # A grid that already fits in one <=5x5x5 tile must use the exact #58
        # production preset. Adding --active-core there would incorrectly
        # select the one-iteration large-grid preset even though there is no
        # halo to freeze and no Gauss-Seidel boundary to protect.
        $args += '--write-uvphi'
        if ($null -ne $ladderStep) {
            $args += @('--tabu-ladder-step',([double]$ladderStep).ToString('R',[Globalization.CultureInfo]::InvariantCulture),
                       '--tabu-ladder-u0',([double]$ladderU0).ToString('R',[Globalization.CultureInfo]::InvariantCulture))
        }
        Write-Host ("[tiled_atlas] sweep {0}/{1}, tile {2}/{3} {4}, loaded={5}" -f
                    ($sweep+1),$Sweeps,$tileNo,$tiles.Count,$tag,$loaded.Count)
        & $AtlasExe @args
        if ($LASTEXITCODE -ne 0) { throw "tile failed: $tag ($LASTEXITCODE)" }
        if ($singleFullGridTile) { $singleTileAtlas = $atlas }
        $stats = Get-Content -LiteralPath (Join-Path $atlas 'overlap_fix_stats.json') -Raw | ConvertFrom-Json
        if ($null -eq $ladderStep -and
            [Math]::Abs([double]$stats.tabu_ladder_step) -gt 1.0e-12 -and
            [int64]$stats.tabu_ladder_pairs -gt 0) {
            $ladderStep = [double]$stats.tabu_ladder_step
            $ladderU0 = [double]$stats.tabu_ladder_u0
            Write-Host ("[tiled_atlas] fixed calibration step={0:R}, u0={1:R}" -f $ladderStep,$ladderU0)
        } elseif ($null -eq $ladderStep) {
            Write-Host "[tiled_atlas] tile $tag has no ladder pairs; leaving calibration unfixed"
        }

        $coreCubes = @($cube | Where-Object {
            $_.IZ -ge $tile.Z.First -and $_.IZ -lt $tile.Z.Last -and
            $_.IY -ge $tile.Y.First -and $_.IY -lt $tile.Y.Last -and
            $_.IX -ge $tile.X.First -and $_.IX -lt $tile.X.Last
        })
        foreach ($c in $coreCubes) {
            $name = "$($c.Id)_uvphi.f32"
            Copy-Item -LiteralPath (Join-Path $atlas "solved_uvphi\$name") `
                      -Destination (Join-Path $state $name) -Force
        }
        $summary.Add([pscustomobject]@{
            sweep=$sweep+1; tile=$tag; loaded_cubes=$loaded.Count
            core_cubes=$coreCubes.Count; active_charts=$stats.tabu_active_charts
            charts_moved=$stats.tabu_charts_moved
            overlap_before=$stats.overlap_pairs_before
            overlap_after=$stats.overlap_pairs_after
            welds_torn=$stats.tabu_welds_torn; seconds=$stats.seconds
        })
    }
}
$summary | Export-Csv -LiteralPath (Join-Path $Out 'tile_summary.csv') -NoTypeInformation

$finalPlaced = Join-Path $Out 'placed'
New-Item -ItemType Directory -Path $finalPlaced | Out-Null
foreach ($file in Get-ChildItem -LiteralPath $Placed -File) {
    $target = Join-Path $finalPlaced $file.Name
    if ($file.Name -like '*_uvphi.f32') {
        New-HardLink $target (Join-Path $state $file.Name)
    } else {
        New-HardLink $target $file.FullName
    }
}

$calibration = [pscustomobject]@{
    ladder_step = $ladderStep; ladder_u0 = $ladderU0
    grid = @($zv.Count,$yv.Count,$xv.Count); cubes = $cube.Count
    tiles = $tiles.Count; sweeps = $Sweeps; wind_seed = [bool]$WindSeed
}
$calibration | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath (Join-Path $Out 'tiled_atlas.json')

if (!$SkipAudit) {
    $atlasFinal = Join-Path $Out 'atlas'
    if ($singleFullGridTile) {
        # Preserve the exact #58 result. Re-loading its float sidecars for a
        # nominal no-move export introduces a needless float round-trip and
        # changes otherwise byte-identical OBJ output.
        Write-Host '[tiled_atlas] one-tile result is the final atlas (exact, no sidecar round-trip)'
        New-Item -ItemType Directory -Path $atlasFinal | Out-Null
        foreach ($file in Get-ChildItem -LiteralPath $singleTileAtlas -File) {
            New-HardLink (Join-Path $atlasFinal $file.Name) $file.FullName
        }
    } else {
        Write-Host '[tiled_atlas] final no-move global audit/export'
        $finalArgs = @($finalPlaced,$atlasFinal,'--production','--no-wind-correct',
                       '--no-relayout','--tabu-iters','0','--tabu-stall','1')
        if ($null -ne $ladderStep) {
            $stepText = ([double]$ladderStep).ToString('R',[Globalization.CultureInfo]::InvariantCulture)
            $u0Text = ([double]$ladderU0).ToString('R',[Globalization.CultureInfo]::InvariantCulture)
            $finalArgs += @('--tabu-ladder-step',$stepText,'--tabu-ladder-u0',$u0Text)
        }
        & $AtlasExe @finalArgs
        if ($LASTEXITCODE -ne 0) { throw "final audit failed: $LASTEXITCODE" }
    }
    if ($Raw -ne '') {
        $bake = Join-Path (Split-Path -Parent $AtlasExe) 'obj_bake_raw.exe'
        $Raw = (Resolve-Path -LiteralPath $Raw).Path
        Write-Host '[tiled_atlas] final RAW texture'
        $duText = $RasterDu.ToString('R',[Globalization.CultureInfo]::InvariantCulture)
        $dvText = $RasterDv.ToString('R',[Globalization.CultureInfo]::InvariantCulture)
        & $bake (Join-Path $atlasFinal 'atlas_bake.obj') $Raw (Join-Path $Out 'texture') `
            --id tiled --raster-du $duText --raster-dv $dvText --require-contrast
        if ($LASTEXITCODE -ne 0) { throw "texture bake failed: $LASTEXITCODE" }
    }
}

Write-Host "[tiled_atlas] complete: $Out"
