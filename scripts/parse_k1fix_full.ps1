$files = Get-ChildItem 'D:\work\vesuvius-c\output\run_20cubes_k1fix\*_log.txt'
$results = @()
foreach ($f in $files) {
    $cube = $f.BaseName -replace '_log',''
    $adj = Select-String -Path $f.FullName -Pattern '\(adj\)' | Select-Object -Last 1
    if (-not $adj) { continue }
    $line = $adj.Line
    $lb = if ($line -match 'LB=([0-9.]+)') { [double]$Matches[1] } else { 0 }
    $topo = if ($line -match 'Topo=([0-9.]+)') { [double]$Matches[1] } else { 0 }
    $sd = if ($line -match 'SDice=([0-9.]+)') { [double]$Matches[1] } else { 0 }
    $voi = if ($line -match 'VOI=([0-9.]+)') { [double]$Matches[1] } else { 0 }
    $split = if ($line -match 'split=([0-9.]+)') { [double]$Matches[1] } else { 0 }
    $merge = if ($line -match 'merge=([0-9.]+)') { [double]$Matches[1] } else { 0 }

    # Get raw (unadjusted) line for raw topo/LB
    $raw = Select-String -Path $f.FullName -Pattern '\(raw\)' | Select-Object -Last 1
    $rawLB = 0; $rawTopo = 0
    if ($raw) {
        if ($raw.Line -match 'LB=([0-9.]+)') { $rawLB = [double]$Matches[1] }
        if ($raw.Line -match 'Topo=([0-9.]+)') { $rawTopo = [double]$Matches[1] }
    }

    # Get Betti F1 line
    $betti = Select-String -Path $f.FullName -Pattern 'Betti F1:' | Select-Object -Last 1
    $k0 = 'nan'; $k0detail = ''
    $k1raw = 'nan'; $k1detail = ''
    $k2 = 'nan'
    if ($betti) {
        $bl = $betti.Line
        if ($bl -match 'k0=([0-9.]+)\s+\(([^)]+)\)') { $k0 = $Matches[1]; $k0detail = $Matches[2] }
        if ($bl -match 'k1=([0-9.]+)\s+\(([^)]+)\)') { $k1raw = $Matches[1]; $k1detail = $Matches[2] }
        if ($bl -match 'k2=([a-z0-9.]+)') { $k2 = $Matches[1] }
    }

    # Get k1 adjust details
    $k1line = Select-String -Path $f.FullName -Pattern 'k1 adjust:' | Select-Object -Last 1
    $k1adj = 'N/A'; $bl_k1 = 'N/A'; $pipe_k1 = 'N/A'; $delta = 'N/A'
    if ($k1line -and $k1line.Line -match 'pipeline_F1_k1=([0-9.]+)\s+baseline_F1_k1=([0-9.]+)\s+\|delta\|=([0-9.]+)\s+adjusted=([0-9.]+)') {
        $pipe_k1 = $Matches[1]
        $bl_k1 = $Matches[2]
        $delta = $Matches[3]
        $k1adj = $Matches[4]
    }

    # Get baseline nnpred line
    $blline = Select-String -Path $f.FullName -Pattern 'Baseline \(raw nnPred\)\s' | Select-Object -Last 1
    $blLB = ''; $blTopo = ''; $blK0 = ''; $blK1 = ''; $blSD = ''; $blVOI = ''
    if ($blline) {
        $bll = $blline.Line
        if ($bll -match 'LB=([0-9.]+)') { $blLB = $Matches[1] }
        if ($bll -match 'Topo=([0-9.]+)') { $blTopo = $Matches[1] }
        if ($bll -match 'F1_0=([0-9.]+)') { $blK0 = $Matches[1] }
        if ($bll -match 'F1_1=([0-9.]+)') { $blK1 = $Matches[1] }
        if ($bll -match 'SDice=([0-9.]+)') { $blSD = $Matches[1] }
        if ($bll -match 'VOI=([0-9.]+)') { $blVOI = $Matches[1] }
    }

    # Get timing
    $tline = Select-String -Path $f.FullName -Pattern 'Total:.*Status:' | Select-Object -Last 1
    $time = 'N/A'; $status = 'N/A'
    if ($tline -and $tline.Line -match 'Total:\s+([0-9.]+)s\s+Status:\s+(\w+)') {
        $time = $Matches[1]
        $status = $Matches[2]
    }

    # Get component count
    $ncomp = 0
    $complines = Select-String -Path $f.FullName -Pattern 'comp_id='
    if ($complines) { $ncomp = $complines.Count }

    $results += [PSCustomObject]@{
        Cube=$cube; LB=$lb; Topo=$topo; K0=[double]$k0; K1Adj=[double]$k1adj; K2=$k2
        SDice=$sd; VOI=$voi; Split=$split; Merge=$merge
        RawLB=$rawLB; RawTopo=$rawTopo; K1Raw=$k1raw; K1Detail=$k1detail; K0Detail=$k0detail
        BL_K1=$bl_k1; Pipe_K1=$pipe_k1; Delta=$delta
        BlLB=$blLB; BlTopo=$blTopo; BlK0=$blK0; BlK1=$blK1; BlSD=$blSD; BlVOI=$blVOI
        Time=$time; Status=$status; NComp=$ncomp
    }
}

$sorted = $results | Sort-Object Cube

Write-Host '===================================================================================='
Write-Host '  VESUVIUS-C  20-CUBE BENCHMARK  (k1 baseline adjustment, bay fill OFF)'
Write-Host '===================================================================================='
Write-Host ''
Write-Host '--- ADJUSTED SCORES (k1 = 1.0 - |pipeline_F1_k1 - baseline_F1_k1|) ---'
Write-Host ''
Write-Host ('  {0,-12} {1,7} {2,7} {3,7} {4,7} {5,7} {6,7} {7,7} {8,7} {9,7} {10,5} {11,6}' -f 'Cube','LB','Topo','F1_k0','F1_k1','F1_k2','SDice','VOI','Split','Merge','Time','#Comp')
Write-Host ('  {0,-12} {1,7} {2,7} {3,7} {4,7} {5,7} {6,7} {7,7} {8,7} {9,7} {10,5} {11,6}' -f '------------','-------','-------','-------','-------','-------','-------','-------','-------','-------','-----','------')

$sumLB=0; $sumTopo=0; $sumK0=0; $sumK1=0; $sumSD=0; $sumVOI=0; $sumSplit=0; $sumMerge=0; $n=0; $nk1=0
foreach ($r in $sorted) {
    $k2str = if ($r.K2 -eq 'nan') { '  nan' } else { '{0,7:F4}' -f [double]$r.K2 }
    $k1str = '{0,7:F4}' -f $r.K1Adj
    Write-Host ('  {0,-12} {1,7:F4} {2,7:F4} {3,7:F4} {4} {5} {6,7:F4} {7,7:F4} {8,7:F3} {9,7:F3} {10,5}s {11,5}' -f $r.Cube, $r.LB, $r.Topo, $r.K0, $k1str, $k2str, $r.SDice, $r.VOI, $r.Split, $r.Merge, $r.Time, $r.NComp)
    $sumLB += $r.LB; $sumTopo += $r.Topo; $sumK0 += $r.K0; $sumSD += $r.SDice; $sumVOI += $r.VOI
    $sumSplit += $r.Split; $sumMerge += $r.Merge
    if ($r.K1Adj -ne 'N/A' -and $r.K1Adj -ge 0) { $sumK1 += $r.K1Adj; $nk1++ }
    $n++
}
Write-Host ('  {0,-12} {1,7} {2,7} {3,7} {4,7} {5,7} {6,7} {7,7} {8,7} {9,7}' -f '------------','-------','-------','-------','-------','-------','-------','-------','-------','-------')
Write-Host ('  MEAN ({0,2})   {1,7:F4} {2,7:F4} {3,7:F4} {4,7:F4}         {5,7:F4} {6,7:F4} {7,7:F3} {8,7:F3}' -f $n, ($sumLB/$n), ($sumTopo/$n), ($sumK0/$n), ($sumK1/$nk1), ($sumSD/$n), ($sumVOI/$n), ($sumSplit/$n), ($sumMerge/$n))

Write-Host ''
Write-Host '--- K1 ADJUSTMENT DETAILS ---'
Write-Host ''
Write-Host ('  {0,-12} {1,10} {2,10} {3,10} {4,10}' -f 'Cube','BL_F1_k1','Pipe_F1k1','|delta|','Adjusted')
Write-Host ('  {0,-12} {1,10} {2,10} {3,10} {4,10}' -f '------------','----------','----------','----------','----------')
foreach ($r in $sorted) {
    Write-Host ('  {0,-12} {1,10} {2,10} {3,10} {4,10}' -f $r.Cube, $r.BL_K1, $r.Pipe_K1, $r.Delta, $r.K1Adj)
}

Write-Host ''
Write-Host '--- RAW BETTI COUNTS (before adjustment) ---'
Write-Host ''
Write-Host ('  {0,-12} {1,10} {2,30} {3,7} {4,7}' -f 'Cube','F1_k0','k0 (matched/pred+gt)','raw_k1','k1 detail')
Write-Host ('  {0,-12} {1,10} {2,30} {3,7} {4,7}' -f '------------','----------','------------------------------','-------','----------')
foreach ($r in $sorted) {
    Write-Host ('  {0,-12} {1,10:F4} {2,30} {3,7} {4}' -f $r.Cube, $r.K0, $r.K0Detail, $r.K1Raw, $r.K1Detail)
}

Write-Host ''
Write-Host '--- BASELINE (raw nnPred vs GT, no pipeline) ---'
Write-Host ''
Write-Host ('  {0,-12} {1,7} {2,7} {3,7} {4,7} {5,7} {6,7}' -f 'Cube','LB','Topo','F1_k0','F1_k1','SDice','VOI')
Write-Host ('  {0,-12} {1,7} {2,7} {3,7} {4,7} {5,7} {6,7}' -f '------------','-------','-------','-------','-------','-------','-------')
foreach ($r in $sorted) {
    Write-Host ('  {0,-12} {1,7} {2,7} {3,7} {4,7} {5,7} {6,7}' -f $r.Cube, $r.BlLB, $r.BlTopo, $r.BlK0, $r.BlK1, $r.BlSD, $r.BlVOI)
}
