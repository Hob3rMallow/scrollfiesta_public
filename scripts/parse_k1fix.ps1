$files = Get-ChildItem 'D:\work\vesuvius-c\output\run_20cubes_k1fix\*_log.txt'
$results = @()
foreach ($f in $files) {
    $cube = $f.BaseName -replace '_log',''
    $adj = Select-String -Path $f.FullName -Pattern '\(adj\)' | Select-Object -Last 1
    if ($adj) {
        $line = $adj.Line
        $lb = if ($line -match 'LB=([0-9.]+)') { [double]$Matches[1] } else { 0 }
        $topo = if ($line -match 'Topo=([0-9.]+)') { [double]$Matches[1] } else { 0 }
        $sd = if ($line -match 'SDice=([0-9.]+)') { [double]$Matches[1] } else { 0 }
        $voi = if ($line -match 'VOI=([0-9.]+)') { [double]$Matches[1] } else { 0 }

        $k1line = Select-String -Path $f.FullName -Pattern 'k1 adjust:' | Select-Object -Last 1
        $k1adj = 'N/A'
        $bl_k1 = 'N/A'
        $pipe_k1 = 'N/A'
        if ($k1line -and $k1line.Line -match 'pipeline_F1_k1=([0-9.]+)\s+baseline_F1_k1=([0-9.]+)\s+\|delta\|=([0-9.]+)\s+adjusted=([0-9.]+)') {
            $pipe_k1 = $Matches[1]
            $bl_k1 = $Matches[2]
            $k1adj = $Matches[4]
        }

        $results += [PSCustomObject]@{
            Cube=$cube; LB=$lb; Topo=$topo; K1Adj=$k1adj; BL_K1=$bl_k1; Pipe_K1=$pipe_k1; SDice=$sd; VOI=$voi
        }
    }
}

Write-Host 'Cube            LB      Topo    K1adj   BL_K1   Pipe_K1 SDice   VOI'
Write-Host '------------ ------- ------- ------- ------- ------- ------- -------'
$sumLB=0; $sumTopo=0; $sumSD=0; $sumVOI=0; $n=0
foreach ($r in $results | Sort-Object Cube) {
    Write-Host ('{0,-12} {1,7:F4} {2,7:F4} {3,7} {4,7} {5,7} {6,7:F4} {7,7:F4}' -f $r.Cube, $r.LB, $r.Topo, $r.K1Adj, $r.BL_K1, $r.Pipe_K1, $r.SDice, $r.VOI)
    $sumLB += $r.LB; $sumTopo += $r.Topo; $sumSD += $r.SDice; $sumVOI += $r.VOI; $n++
}
Write-Host '------------ ------- ------- ------- ------- ------- ------- -------'
Write-Host ('MEAN ({0,2})    {1,7:F4} {2,7:F4}                         {3,7:F4} {4,7:F4}' -f $n, ($sumLB/$n), ($sumTopo/$n), ($sumSD/$n), ($sumVOI/$n))
Write-Host ''
Write-Host 'Previous run (no k1 adjustment): Mean LB=0.5385, Topo=0.2721'
