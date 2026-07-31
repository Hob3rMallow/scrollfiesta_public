<#
  scan_components.ps1 -- sweep obj_components --brief over every per-cube
  component OBJ in a dump tree to find BPA front-splits (one sheet left as two
  disconnected islands a sub-voxel gap apart). Runner only; all analysis is in
  src/tools/obj_components.c. Logs every line, then summarises by class and
  ranks the SEAM-SPLIT / *-DOUBLE candidates (the "should-be-one-sheet" cases).

  Usage:
    pwsh scripts/step0/scan_components.ps1 [-Stage pre_simplify] [-Gap 1.5] [-Top 30]
#>
param(
    [string]$DumpRoot = "D:\work\vesuvius-c\output\grid_seamweld\dump",
    [string]$Exe      = "D:\work\vesuvius-c\build\Release\obj_components.exe",
    [string]$Stage    = "pre_simplify",     # matches *_<Stage>_comp*.obj
    [string]$LogDir   = "D:\work\vesuvius-c\output\tools",
    [double]$Gap      = 1.5,
    [int]$Top         = 30
)
$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
$log = Join-Path $LogDir ("scan_components_{0}.log" -f $Stage)

$pattern = "*_{0}_comp*.obj" -f $Stage
$files = Get-ChildItem -Path $DumpRoot -Recurse -Filter $pattern -File
Write-Host ("scanning {0} files (stage={1}, gap={2}) ..." -f $files.Count, $Stage, $Gap)

$rows = New-Object System.Collections.Generic.List[object]
$lines = New-Object System.Collections.Generic.List[string]
$rx = [regex]'comps=(?<c>\d+)\s+(?:c(?<a>\d+)-c(?<b>\d+)\s+gap=(?<g>[\d.]+)\s+near=(?<n>\d+)\s+cover=(?<ca>[\d.]+)/(?<cb>[\d.]+)%\s+)?(?<tag>\S+)\s*$'
foreach ($f in $files) {
    $line = (& $Exe $f.FullName --brief --gap $Gap) -join ""
    $lines.Add($line)
    $m = $rx.Match($line)
    if ($m.Success) {
        $rows.Add([PSCustomObject]@{
            File  = $f.Name
            Comps = [int]$m.Groups['c'].Value
            Tag   = $m.Groups['tag'].Value
            Gap   = if ($m.Groups['g'].Success)  { [double]$m.Groups['g'].Value } else { [double]::NaN }
            Near  = if ($m.Groups['n'].Success)  { [int]$m.Groups['n'].Value }   else { 0 }
            CovA  = if ($m.Groups['ca'].Success) { [double]$m.Groups['ca'].Value } else { 0 }
            CovB  = if ($m.Groups['cb'].Success) { [double]$m.Groups['cb'].Value } else { 0 }
        })
    }
}
$lines | Out-File -FilePath $log -Encoding utf8
Write-Host ("full log -> {0}" -f $log)

Write-Host "`n=== class counts ==="
$rows | Group-Object Tag | Sort-Object Count -Descending |
    ForEach-Object { "{0,-16} {1}" -f $_.Name, $_.Count }

$split = $rows | Where-Object { $_.Tag -in @('SEAM-SPLIT','PARTIAL-DOUBLE','DOUBLE') }
Write-Host ("`n=== {0} split candidates, ranked by near-pair count ===" -f $split.Count)
$split | Sort-Object Near -Descending | Select-Object -First $Top |
    Format-Table File, Comps, Tag, @{N='Gap';E={'{0:N3}' -f $_.Gap}}, Near,
        @{N='Cov%';E={'{0:N0}/{1:N0}' -f $_.CovA, $_.CovB}} -AutoSize

Write-Host "=== cleanest 2-component splits (no main body to confuse the picture) ==="
$split | Where-Object { $_.Comps -eq 2 } | Sort-Object Near -Descending | Select-Object -First $Top |
    Format-Table File, Tag, @{N='Gap';E={'{0:N3}' -f $_.Gap}}, Near,
        @{N='Cov%';E={'{0:N0}/{1:N0}' -f $_.CovA, $_.CovB}} -AutoSize
