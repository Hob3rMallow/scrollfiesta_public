[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$SourceDump,

    [Parameter(Mandatory = $true)]
    [string]$OutDump,

    [ValidateRange(0.001, 0.999)]
    [double]$KeepRatio = 0.05,

    [ValidateRange(1, 64)]
    [int]$Jobs = 16,

    [string]$QslimExe
)

$ErrorActionPreference = 'Stop'

if ($PSVersionTable.PSVersion.Major -lt 7) {
    throw 'make_atlas_proxy.ps1 requires PowerShell 7 for bounded parallel execution.'
}

$repoRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
if (-not $QslimExe) {
    $QslimExe = Join-Path $repoRoot 'build\Release\qslim_obj.exe'
}

$SourceDump = [IO.Path]::GetFullPath($SourceDump)
$OutDump = [IO.Path]::GetFullPath($OutDump)
$QslimExe = [IO.Path]::GetFullPath($QslimExe)

if (-not (Test-Path -LiteralPath $SourceDump -PathType Container)) {
    throw "Source dump does not exist: $SourceDump"
}
if (-not (Test-Path -LiteralPath $QslimExe -PathType Leaf)) {
    throw "qslim_obj is not built: $QslimExe"
}

[IO.Directory]::CreateDirectory($OutDump) | Out-Null

$items = @(
    Get-ChildItem -LiteralPath $SourceDump -Directory |
        Sort-Object Name |
        ForEach-Object {
            $id = $_.Name
            $leaf = "${id}_step12_final_all.obj"
            $input = Join-Path $_.FullName "${id}_step12_final\$leaf"
            if (-not (Test-Path -LiteralPath $input -PathType Leaf)) {
                throw "Missing step12 final OBJ for ${id}: $input"
            }
            [pscustomobject]@{
                Id = $id
                Input = $input
                Output = Join-Path $OutDump "$id\${id}_step12_final\$leaf"
            }
        }
)

if ($items.Count -eq 0) {
    throw "No cube directories found under $SourceDump"
}

$todo = @(
    $items | Where-Object {
        -not (Test-Path -LiteralPath $_.Output -PathType Leaf) -or
        (Get-Item -LiteralPath $_.Output).Length -eq 0
    }
)

Write-Host ("Atlas proxy: {0} cubes, {1} to simplify, keep ratio {2:P2}, jobs {3}" -f
    $items.Count, $todo.Count, $KeepRatio, $Jobs)

$started = Get-Date
$results = @(
    $todo | ForEach-Object -ThrottleLimit $Jobs -Parallel {
        $item = $_
        $ratio = $using:KeepRatio
        [IO.Directory]::CreateDirectory([IO.Path]::GetDirectoryName($item.Output)) | Out-Null

        $psi = [Diagnostics.ProcessStartInfo]::new()
        $psi.FileName = $using:QslimExe
        $psi.UseShellExecute = $false
        $psi.CreateNoWindow = $true
        $psi.RedirectStandardOutput = $true
        $psi.RedirectStandardError = $true
        $psi.ArgumentList.Add($item.Input)
        $psi.ArgumentList.Add($item.Output)
        $psi.ArgumentList.Add($ratio.ToString('R', [Globalization.CultureInfo]::InvariantCulture))

        $process = [Diagnostics.Process]::Start($psi)
        $stdout = $process.StandardOutput.ReadToEndAsync()
        $stderr = $process.StandardError.ReadToEndAsync()
        $process.WaitForExit()
        [Threading.Tasks.Task]::WaitAll(@($stdout, $stderr))

        [pscustomobject]@{
            Id = $item.Id
            ExitCode = $process.ExitCode
            Stdout = $stdout.Result.Trim()
            Stderr = $stderr.Result.Trim()
        }
    }
)

$failed = @($results | Where-Object ExitCode -ne 0)
foreach ($failure in $failed) {
    Write-Error ("qslim_obj failed for {0} (exit {1})`n{2}`n{3}" -f
        $failure.Id, $failure.ExitCode, $failure.Stdout, $failure.Stderr) -ErrorAction Continue
}
if ($failed.Count -gt 0) {
    throw "$($failed.Count) cube simplification job(s) failed; successful outputs are resumable."
}

$manifest = @(
    foreach ($item in $items) {
        if (-not (Test-Path -LiteralPath $item.Output -PathType Leaf)) {
            throw "Proxy output is missing after simplification: $($item.Output)"
        }
        $sourceInfo = Get-Item -LiteralPath $item.Input
        $proxyInfo = Get-Item -LiteralPath $item.Output
        [pscustomobject]@{
            cube_id = $item.Id
            source_bytes = $sourceInfo.Length
            proxy_bytes = $proxyInfo.Length
            byte_ratio = if ($sourceInfo.Length) {
                [Math]::Round($proxyInfo.Length / $sourceInfo.Length, 6)
            } else { 0 }
            source_obj = $item.Input
            proxy_obj = $item.Output
        }
    }
)

$manifestPath = Join-Path $OutDump 'proxy_manifest.csv'
$manifest | Export-Csv -LiteralPath $manifestPath -NoTypeInformation
$elapsed = (Get-Date) - $started
$sourceBytes = ($manifest | Measure-Object source_bytes -Sum).Sum
$proxyBytes = ($manifest | Measure-Object proxy_bytes -Sum).Sum

Write-Host ("Atlas proxy complete: {0} cubes, {1:N2} GiB -> {2:N2} GiB ({3:P2}), elapsed {4}" -f
    $manifest.Count,
    ($sourceBytes / 1GB),
    ($proxyBytes / 1GB),
    ($proxyBytes / $sourceBytes),
    $elapsed)
Write-Host "Manifest: $manifestPath"
