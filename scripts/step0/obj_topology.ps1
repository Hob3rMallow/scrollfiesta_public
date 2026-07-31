<#
  obj_topology.ps1 — mesh topology analyzer for the genus-0 debugging effort.

  Reports, for an OBJ surface mesh:
    V, F, E, Euler characteristic chi = V - E + F
    connected components (vertex connectivity)
    boundary edges (1 incident face), non-manifold edges (>2 faces)
    boundary loops: count + sorted lengths  (largest = outer; rest = holes)
    total genus G = (2*C - B - chi) / 2     (C=components, B=boundary loops)

  A clean genus-0 sheet (topological disk) is:  C=1, B=1, chi=1, G=0,
  non-manifold edges = 0.

  Usage:
    pwsh scripts/step0/obj_topology.ps1 -Obj <path.obj> [-DumpBoundary <out.obj>]

  -DumpBoundary writes the boundary loops as colored OBJ line segments
  (outer loop green, internal holes red) for MeshLab inspection.
#>
param(
    [Parameter(Mandatory=$true)][string]$Obj,
    [string]$DumpBoundary = "",
    [switch]$Quiet
)

$ErrorActionPreference = "Stop"

Add-Type -TypeDefinition @'
using System;
using System.Collections.Generic;
using System.IO;
using System.Globalization;
using System.Text;

public static class ObjTopo {
    static long Key(int a, int b){ int lo = a<b?a:b; int hi = a<b?b:a; return ((long)lo<<32) | (uint)hi; }

    public class Report {
        public int V, F, E, boundaryE, nonManifoldE, components;
        public int chi;
        public long genus;
        public List<int> loopLens = new List<int>();
        public int pinchVerts;   // boundary verts with out-degree != in-degree (branchy boundary)
    }

    public static Report Analyze(string path, string boundaryOut) {
        var vx = new List<float>(); var vy = new List<float>(); var vz = new List<float>();
        var faces = new List<int>();
        char[] sep = null; // split on whitespace
        foreach (var raw in File.ReadLines(path)) {
            if (raw.Length < 2) continue;
            char c0 = raw[0];
            if (c0=='v' && raw[1]==' ') {
                var t = raw.Split(sep, StringSplitOptions.RemoveEmptyEntries);
                vx.Add(float.Parse(t[1], CultureInfo.InvariantCulture));
                vy.Add(float.Parse(t[2], CultureInfo.InvariantCulture));
                vz.Add(float.Parse(t[3], CultureInfo.InvariantCulture));
            } else if (c0=='f' && raw[1]==' ') {
                var t = raw.Split(sep, StringSplitOptions.RemoveEmptyEntries);
                int n = t.Length - 1;
                if (n < 3) continue;
                int[] idx = new int[n];
                for (int i=0;i<n;i++){
                    string tok = t[i+1];
                    int sl = tok.IndexOf('/');
                    if (sl>=0) tok = tok.Substring(0, sl);
                    int vi = int.Parse(tok, CultureInfo.InvariantCulture);
                    if (vi<0) vi = vx.Count + vi; else vi = vi - 1;
                    idx[i] = vi;
                }
                for (int i=1;i+1<n;i++){ faces.Add(idx[0]); faces.Add(idx[i]); faces.Add(idx[i+1]); }
            }
        }
        int Vtot = vx.Count;
        int F = faces.Count/3;

        var edgeCount = new Dictionary<long,int>(F*3);
        for (int f=0; f<F; f++){
            int a=faces[f*3], b=faces[f*3+1], c=faces[f*3+2];
            long k0=Key(a,b), k1=Key(b,c), k2=Key(c,a);
            int cur;
            edgeCount.TryGetValue(k0, out cur); edgeCount[k0]=cur+1;
            edgeCount.TryGetValue(k1, out cur); edgeCount[k1]=cur+1;
            edgeCount.TryGetValue(k2, out cur); edgeCount[k2]=cur+1;
        }
        int E = edgeCount.Count;
        int boundaryE=0, nonManifoldE=0;
        foreach (var kv in edgeCount){ if (kv.Value==1) boundaryE++; else if (kv.Value>2) nonManifoldE++; }

        // Union-find over referenced vertices (component count)
        var parent = new Dictionary<int,int>();
        Func<int,int> find = null;
        find = delegate(int x){ int r=x; while(parent[r]!=r) r=parent[r]; while(parent[x]!=r){ int nx=parent[x]; parent[x]=r; x=nx;} return r; };
        Action<int> add = delegate(int x){ if(!parent.ContainsKey(x)) parent[x]=x; };
        for (int f=0; f<F; f++){
            int a=faces[f*3], b=faces[f*3+1], c=faces[f*3+2];
            add(a); add(b); add(c);
            int ra=find(a), rb=find(b), rc=find(c);
            if(ra!=rb) parent[rb]=ra;
            ra=find(a); rc=find(c);
            if(ra!=rc) parent[rc]=ra;
        }
        var roots = new HashSet<int>();
        foreach (var v in parent.Keys) roots.Add(find(v));
        int components = roots.Count;
        int Vref = parent.Count;

        // Directed boundary half-edges (interior on the left). Build next[] map.
        var nextOf = new Dictionary<int,int>();
        var indeg  = new Dictionary<int,int>();
        var outdeg = new Dictionary<int,int>();
        for (int f=0; f<F; f++){
            int a=faces[f*3], b=faces[f*3+1], c=faces[f*3+2];
            int[] hu = {a,b,c}; int[] hv = {b,c,a};
            for(int e=0;e<3;e++){
                int u=hu[e], v=hv[e];
                if (edgeCount[Key(u,v)]==1){
                    nextOf[u]=v;
                    int d; outdeg.TryGetValue(u,out d); outdeg[u]=d+1;
                    indeg.TryGetValue(v,out d); indeg[v]=d+1;
                }
            }
        }
        int pinch=0;
        var allBV = new HashSet<int>();
        foreach(var u in outdeg.Keys) allBV.Add(u);
        foreach(var u in indeg.Keys) allBV.Add(u);
        foreach(var u in allBV){ int o,i2; outdeg.TryGetValue(u,out o); indeg.TryGetValue(u,out i2); if(o!=1||i2!=1) pinch++; }

        // Trace boundary loops
        var visited = new HashSet<int>();
        var loopLens = new List<int>();
        var loops = new List<List<int>>();
        foreach (var start in nextOf.Keys){
            if (visited.Contains(start)) continue;
            int cur=start, len=0; var seq=new List<int>();
            while (!visited.Contains(cur) && nextOf.ContainsKey(cur)){
                visited.Add(cur); seq.Add(cur); cur=nextOf[cur]; len++;
                if (len > Vtot+5) break;
            }
            loopLens.Add(len); loops.Add(seq);
        }
        loopLens.Sort(); loopLens.Reverse();

        var rep = new Report();
        rep.V=Vref; rep.F=F; rep.E=E; rep.boundaryE=boundaryE; rep.nonManifoldE=nonManifoldE;
        rep.components=components; rep.chi=Vref - E + F; rep.loopLens=loopLens; rep.pinchVerts=pinch;
        rep.genus = (2L*components - loopLens.Count - rep.chi)/2;

        if (!string.IsNullOrEmpty(boundaryOut)){
            // largest loop = outer (green), rest = holes (red)
            int maxLen = -1, maxIdx = -1;
            for (int i=0;i<loops.Count;i++) if (loops[i].Count>maxLen){ maxLen=loops[i].Count; maxIdx=i; }
            var sb = new StringBuilder();
            var remap = new Dictionary<int,int>(); int nextId=1;
            for (int i=0;i<loops.Count;i++){
                bool outer = (i==maxIdx);
                float r = outer?0.1f:1.0f, g = outer?1.0f:0.1f, bl = 0.1f;
                foreach (int v in loops[i]){
                    if (!remap.ContainsKey(v)){
                        remap[v]=nextId++;
                        sb.Append("v ").Append(vx[v].ToString(CultureInfo.InvariantCulture)).Append(' ')
                          .Append(vy[v].ToString(CultureInfo.InvariantCulture)).Append(' ')
                          .Append(vz[v].ToString(CultureInfo.InvariantCulture)).Append(' ')
                          .Append(r).Append(' ').Append(g).Append(' ').Append(bl).Append('\n');
                    }
                }
            }
            for (int i=0;i<loops.Count;i++){
                var seq=loops[i];
                for (int j=0;j<seq.Count;j++){
                    int a=remap[seq[j]], b=remap[seq[(j+1)%seq.Count]];
                    sb.Append("l ").Append(a).Append(' ').Append(b).Append('\n');
                }
            }
            File.WriteAllText(boundaryOut, sb.ToString());
        }
        return rep;
    }
}
'@ -Language CSharp

$rep = [ObjTopo]::Analyze($Obj, $DumpBoundary)

if ($Quiet) {
    $holes   = [Math]::Max(0, $rep.loopLens.Count - 1)
    $maxhole = if ($rep.loopLens.Count -gt 1) { $rep.loopLens[1] } else { 0 }
    Write-Output ("{0,-46} F={1,-7} comps={2,-2} nonmanifold={3,-3} loops={4,-3} holes={5,-3} maxhole={6,-4} genus={7}" -f `
        (Split-Path $Obj -Leaf), $rep.F, $rep.components, $rep.nonManifoldE, $rep.loopLens.Count, $holes, $maxhole, $rep.genus)
    return
}

$name = Split-Path $Obj -Leaf
Write-Host ("=== topology: {0} ===" -f $name)
Write-Host ("  V (referenced) : {0}" -f $rep.V)
Write-Host ("  F (triangles)  : {0}" -f $rep.F)
Write-Host ("  E (edges)      : {0}" -f $rep.E)
Write-Host ("  Euler chi      : {0}   (V - E + F)" -f $rep.chi)
Write-Host ("  components     : {0}" -f $rep.components)
Write-Host ("  non-manifold E : {0}   (>2 faces; must be 0)" -f $rep.nonManifoldE)
Write-Host ("  boundary E     : {0}" -f $rep.boundaryE)
Write-Host ("  boundary loops : {0}" -f $rep.loopLens.Count)
if ($rep.pinchVerts -gt 0) {
    Write-Host ("  PINCH verts    : {0}   (boundary branch points; loop trace approximate)" -f $rep.pinchVerts)
}
$top = $rep.loopLens | Select-Object -First 12
Write-Host ("  loop lengths   : {0}" -f ($top -join ", "))
Write-Host ("  total genus G  : {0}   (0 = clean)" -f $rep.genus)
$holes = [Math]::Max(0, $rep.loopLens.Count - 1)
$verdict = if ($rep.components -eq 1 -and $rep.loopLens.Count -eq 1 -and $rep.nonManifoldE -eq 0 -and $rep.genus -eq 0) { "CLEAN genus-0 disk" } else { "NOT clean: $($rep.components) comp(s), $($rep.loopLens.Count) loop(s) => $holes internal hole(s), genus $($rep.genus), $($rep.nonManifoldE) non-manifold" }
Write-Host ("  VERDICT        : {0}" -f $verdict)
if ($DumpBoundary -ne "") { Write-Host ("  boundary OBJ   : {0}" -f $DumpBoundary) }
