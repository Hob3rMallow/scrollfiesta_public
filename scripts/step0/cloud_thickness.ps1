<#
  cloud_thickness.ps1 — single-vs-double-layer test for a LOP point cloud.

  For each point, gather neighbours within -Radius, do PCA on the local
  neighbourhood, and measure the spread along the SMALLEST eigenvector
  (the through-sheet direction):
      thickness = max-min of neighbour projections onto the normal.
  A single median surface reads ~0 (sub-voxel noise/curvature); a double
  envelope (two faces a gap apart) reads ~the gap (1-2 vox), because the
  neighbourhood spans both faces.

  Also reports in-plane NN spacing for scale. Prints a thickness histogram
  + medians + a verdict.

  Usage: pwsh scripts/step0/cloud_thickness.ps1 -Obj <cloud.obj> [-Radius 1.5]
#>
param(
    [Parameter(Mandatory=$true)][string]$Obj,
    [double]$Radius = 1.5
)
$ErrorActionPreference = "Stop"

Add-Type -TypeDefinition @'
using System;
using System.Collections.Generic;
using System.IO;
using System.Globalization;

public static class CloudThk {
    // Cyclic Jacobi for a symmetric 3x3; returns eigenvalues ascending in ev[]
    // and corresponding eigenvectors as columns of V[][].
    static void Jacobi(double[][] a, double[] ev, double[][] V) {
        for (int i=0;i<3;i++){ for(int j=0;j<3;j++) V[i][j] = (i==j)?1.0:0.0; }
        for (int sweep=0; sweep<24; sweep++) {
            double off = Math.Abs(a[0][1])+Math.Abs(a[0][2])+Math.Abs(a[1][2]);
            if (off < 1e-18) break;
            for (int p=0;p<2;p++) for (int q=p+1;q<3;q++) {
                if (Math.Abs(a[p][q]) < 1e-20) continue;
                double th = (a[q][q]-a[p][p]) / (2.0*a[p][q]);
                double t = Math.Sign(th) / (Math.Abs(th)+Math.Sqrt(th*th+1.0));
                if (th==0) t = 1.0;
                double c = 1.0/Math.Sqrt(t*t+1.0), s = t*c;
                double app=a[p][p], aqq=a[q][q], apq=a[p][q];
                a[p][p]=c*c*app-2*s*c*apq+s*s*aqq;
                a[q][q]=s*s*app+2*s*c*apq+c*c*aqq;
                a[p][q]=0; a[q][p]=0;
                int r=3-p-q;
                double arp=a[r][p], arq=a[r][q];
                a[r][p]=c*arp-s*arq; a[p][r]=a[r][p];
                a[r][q]=s*arp+c*arq; a[q][r]=a[r][q];
                for (int k=0;k<3;k++){ double vkp=V[k][p], vkq=V[k][q]; V[k][p]=c*vkp-s*vkq; V[k][q]=s*vkp+c*vkq; }
            }
        }
        for (int i=0;i<3;i++) ev[i]=a[i][i];
        // sort ascending, reorder eigenvector columns
        for (int i=0;i<2;i++) for (int j=i+1;j<3;j++) if (ev[j]<ev[i]) {
            double tmp=ev[i]; ev[i]=ev[j]; ev[j]=tmp;
            for (int k=0;k<3;k++){ double t2=V[k][i]; V[k][i]=V[k][j]; V[k][j]=t2; }
        }
    }

    public static string Go(string path, double radius) {
        var P = new List<double[]>();
        char[] sep = null;
        foreach (var raw in File.ReadLines(path)) {
            if (raw.Length < 2 || raw[0] != 'v' || raw[1] != ' ') continue;
            var t = raw.Split(sep, StringSplitOptions.RemoveEmptyEntries);
            P.Add(new double[]{ double.Parse(t[1],CultureInfo.InvariantCulture),
                                double.Parse(t[2],CultureInfo.InvariantCulture),
                                double.Parse(t[3],CultureInfo.InvariantCulture) });
        }
        int n = P.Count;
        // uniform grid, cell = radius
        double xmin=double.MaxValue,ymin=double.MaxValue,zmin=double.MaxValue;
        for (int i=0;i<n;i++){ if(P[i][0]<xmin)xmin=P[i][0]; if(P[i][1]<ymin)ymin=P[i][1]; if(P[i][2]<zmin)zmin=P[i][2]; }
        double cs = radius;
        var grid = new Dictionary<long,List<int>>();
        Func<int,int,int,long> key = (a,b,c)=>((long)(a&0x1FFFFF)<<42)|((long)(b&0x1FFFFF)<<21)|(long)(c&0x1FFFFF);
        for (int i=0;i<n;i++){
            int cx=(int)((P[i][0]-xmin)/cs), cy=(int)((P[i][1]-ymin)/cs), cz=(int)((P[i][2]-zmin)/cs);
            long k=key(cx,cy,cz); List<int> l; if(!grid.TryGetValue(k,out l)){ l=new List<int>(); grid[k]=l; } l.Add(i);
        }

        var thick = new List<double>();
        var nnsp  = new List<double>();
        double r2 = radius*radius;
        for (int i=0;i<n;i++){
            int cx=(int)((P[i][0]-xmin)/cs), cy=(int)((P[i][1]-ymin)/cs), cz=(int)((P[i][2]-zmin)/cs);
            var nb = new List<int>();
            double nn = double.MaxValue;
            for (int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++)for(int dz=-1;dz<=1;dz++){
                List<int> l; if(!grid.TryGetValue(key(cx+dx,cy+dy,cz+dz),out l))continue;
                foreach (int j in l){
                    double d0=P[j][0]-P[i][0],d1=P[j][1]-P[i][1],d2=P[j][2]-P[i][2];
                    double d=d0*d0+d1*d1+d2*d2;
                    if (j!=i && d<nn) nn=d;
                    if (d<=r2) nb.Add(j);
                }
            }
            if (nn<double.MaxValue) nnsp.Add(Math.Sqrt(nn));
            if (nb.Count < 6) continue;
            double mx=0,my=0,mz=0;
            foreach(int j in nb){ mx+=P[j][0]; my+=P[j][1]; mz+=P[j][2]; }
            mx/=nb.Count; my/=nb.Count; mz/=nb.Count;
            double[][] cov = new double[][]{ new double[3], new double[3], new double[3] };
            foreach(int j in nb){
                double e0=P[j][0]-mx,e1=P[j][1]-my,e2=P[j][2]-mz;
                cov[0][0]+=e0*e0; cov[0][1]+=e0*e1; cov[0][2]+=e0*e2;
                cov[1][1]+=e1*e1; cov[1][2]+=e1*e2; cov[2][2]+=e2*e2;
            }
            cov[1][0]=cov[0][1]; cov[2][0]=cov[0][2]; cov[2][1]=cov[1][2];
            double[] ev=new double[3]; double[][] V=new double[][]{ new double[3],new double[3],new double[3] };
            Jacobi(cov, ev, V);
            double nx=V[0][0], ny=V[1][0], nz=V[2][0];  // smallest-eigenvalue eigenvector = normal
            double pmin=double.MaxValue, pmax=-double.MaxValue;
            foreach(int j in nb){
                double pr=(P[j][0]-mx)*nx+(P[j][1]-my)*ny+(P[j][2]-mz)*nz;
                if(pr<pmin)pmin=pr; if(pr>pmax)pmax=pr;
            }
            thick.Add(pmax-pmin);
        }

        thick.Sort(); nnsp.Sort();
        Func<List<double>,double,double> pct = (L,p)=> L.Count==0?0:L[(int)Math.Min(L.Count-1,Math.Max(0,p*L.Count))];
        double[] edges = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,2.0,3.0,1e9};
        int[] bins = new int[edges.Length-1];
        foreach(double v in thick){ for(int b=0;b<bins.Length;b++){ if(v>=edges[b]&&v<edges[b+1]){bins[b]++;break;} } }
        var sb = new System.Text.StringBuilder();
        sb.AppendFormat("points={0}  measured={1}  radius={2}\n", n, thick.Count, radius);
        sb.AppendFormat("in-plane NN spacing: median={0:F3}  p10={1:F3}  p90={2:F3} vox\n",
            pct(nnsp,0.5), pct(nnsp,0.1), pct(nnsp,0.9));
        sb.AppendFormat("local THICKNESS (spread along normal): median={0:F3}  p90={1:F3}  p99={2:F3} vox\n",
            pct(thick,0.5), pct(thick,0.9), pct(thick,0.99));
        sb.Append("thickness histogram (vox):\n");
        string[] lbl = {"0.00-0.25","0.25-0.50","0.50-0.75","0.75-1.00","1.00-1.25","1.25-1.50","1.50-2.00","2.00-3.00","3.00+   "};
        for(int b=0;b<bins.Length;b++){
            double frac = thick.Count>0 ? 100.0*bins[b]/thick.Count : 0;
            int bar = (int)Math.Round(frac/2.0);
            sb.AppendFormat("  {0}  {1,7} {2,5:F1}%  {3}\n", lbl[b], bins[b], frac, new string('#',bar));
        }
        double thickFrac = 0; foreach(double v in thick) if(v>0.9) thickFrac++; thickFrac = thick.Count>0?100.0*thickFrac/thick.Count:0;
        string verdict = (pct(thick,0.5) < 0.4 && thickFrac < 10)
            ? "SINGLE layer (thin/planar — Hoppe is a clean drop-in)"
            : ("DOUBLE/thick: " + thickFrac.ToString("F0") + "% of points have >0.9 vox local thickness (Hoppe needs per-layer / thin-first)");
        sb.AppendFormat("VERDICT: {0}\n", verdict);
        return sb.ToString();
    }
}
'@ -Language CSharp

Write-Host ([CloudThk]::Go($Obj, $Radius))
