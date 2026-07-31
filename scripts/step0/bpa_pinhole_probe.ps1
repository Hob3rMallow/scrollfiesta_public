<#
  bpa_pinhole_probe.ps1 — for each small (single-triangle) hole in a BPA mesh,
  reconstruct the ball-pivot closing ball and decide WHY BPA left it open:

    * circumradius R_c of the hole triangle vs rho:
        R_c > rho  => the radius-rho ball physically cannot touch all 3 verts.
                      The hole is "too big for rho" (a radius problem, not float).
    * if R_c <= rho, compute both candidate ball centers O+ / O- and find the
      nearest OTHER vertex. empty-ball margin = dist(O, nearest) - (rho - 1e-4):
        margin >> 0            => ball empty; BPA *should* have closed it.
                                  Rejection was orientation/winding/theta, NOT
                                  empty-ball. (front-collision artifact)
        0 < margin < ~5e-4     => BORDERLINE: a point sits within float noise of
                                  the ball surface. THIS is the precision case.
        margin < 0 (blocker)   => a 4th point genuinely inside the ball; the
                                  empty-ball rejection is correct (dense/double
                                  layer), not a bug.

  Reports the distribution so we know which mechanism dominates.

  NOTE: probes the dumped (world float32) coords. Geometry (distances, radii) is
  scale-invariant; only the borderline (~1e-4) bucket is sensitive to the world-
  vs-local quantization difference, and those are flagged as the precision cases.

  Usage: pwsh scripts/step0/bpa_pinhole_probe.ps1 -Obj <path.obj> [-Rho 1.5] [-MaxLoop 3]
#>
param(
    [Parameter(Mandatory=$true)][string]$Obj,
    [double]$Rho = 1.5,
    [int]$MaxLoop = 3
)
$ErrorActionPreference = "Stop"

Add-Type -TypeDefinition @'
using System;
using System.Collections.Generic;
using System.IO;
using System.Globalization;

public static class Pinhole {
    static long Key(int a,int b){ int lo=a<b?a:b,hi=a<b?b:a; return ((long)lo<<32)|(uint)hi; }

    public static string Probe(string path, double rho, int maxLoop){
        var X=new List<double>(); var Y=new List<double>(); var Z=new List<double>();
        var faces=new List<int>();
        char[] sep=null;
        foreach(var raw in File.ReadLines(path)){
            if(raw.Length<2) continue;
            if(raw[0]=='v'&&raw[1]==' '){
                var t=raw.Split(sep,StringSplitOptions.RemoveEmptyEntries);
                X.Add(double.Parse(t[1],CultureInfo.InvariantCulture));
                Y.Add(double.Parse(t[2],CultureInfo.InvariantCulture));
                Z.Add(double.Parse(t[3],CultureInfo.InvariantCulture));
            } else if(raw[0]=='f'&&raw[1]==' '){
                var t=raw.Split(sep,StringSplitOptions.RemoveEmptyEntries);
                int n=t.Length-1; if(n<3) continue;
                int[] id=new int[n];
                for(int i=0;i<n;i++){ string tk=t[i+1]; int s=tk.IndexOf('/'); if(s>=0) tk=tk.Substring(0,s);
                    int vi=int.Parse(tk,CultureInfo.InvariantCulture); id[i]= vi<0? X.Count+vi : vi-1; }
                for(int i=1;i+1<n;i++){ faces.Add(id[0]); faces.Add(id[i]); faces.Add(id[i+1]); }
            }
        }
        int F=faces.Count/3, V=X.Count;
        var ec=new Dictionary<long,int>(F*3);
        for(int f=0;f<F;f++){ int a=faces[f*3],b=faces[f*3+1],c=faces[f*3+2];
            long[] ks={Key(a,b),Key(b,c),Key(c,a)}; foreach(var k in ks){ int cu; ec.TryGetValue(k,out cu); ec[k]=cu+1; } }
        // directed boundary next[]
        var nextOf=new Dictionary<int,int>();
        for(int f=0;f<F;f++){ int a=faces[f*3],b=faces[f*3+1],c=faces[f*3+2];
            int[] hu={a,b,c}; int[] hv={b,c,a};
            for(int e=0;e<3;e++){ int u=hu[e],v=hv[e]; if(ec[Key(u,v)]==1) nextOf[u]=v; } }
        // trace loops
        var vis=new HashSet<int>(); var loops=new List<List<int>>();
        foreach(var st in nextOf.Keys){ if(vis.Contains(st))continue; int cur=st; var seq=new List<int>();
            while(!vis.Contains(cur)&&nextOf.ContainsKey(cur)){ vis.Add(cur); seq.Add(cur); cur=nextOf[cur]; if(seq.Count>V+5)break; }
            loops.Add(seq); }

        // spatial grid for nearest-other-point queries (cell = 2*rho)
        double cs=2.0*rho;
        double xmin=double.MaxValue,ymin=double.MaxValue,zmin=double.MaxValue;
        for(int i=0;i<V;i++){ if(X[i]<xmin)xmin=X[i]; if(Y[i]<ymin)ymin=Y[i]; if(Z[i]<zmin)zmin=Z[i]; }
        var grid=new Dictionary<long,List<int>>();
        Func<int,int,int,long> ck=(a,b,c)=>((long)(a&0x1FFFFF)<<42)|((long)(b&0x1FFFFF)<<21)|(long)(c&0x1FFFFF);
        for(int i=0;i<V;i++){ int cx=(int)((X[i]-xmin)/cs),cy=(int)((Y[i]-ymin)/cs),cz=(int)((Z[i]-zmin)/cs);
            long k=ck(cx,cy,cz); List<int> l; if(!grid.TryGetValue(k,out l)){ l=new List<int>(); grid[k]=l; } l.Add(i); }

        Func<double,double,double,int,int,int,double> nearestOther=(ox,oy,oz,e0,e1,e2)=>{
            int cx=(int)((ox-xmin)/cs),cy=(int)((oy-ymin)/cs),cz=(int)((oz-zmin)/cs);
            double best=double.MaxValue;
            for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++)for(int dz=-1;dz<=1;dz++){
                List<int> l; if(!grid.TryGetValue(ck(cx+dx,cy+dy,cz+dz),out l))continue;
                foreach(int p in l){ if(p==e0||p==e1||p==e2)continue;
                    double d=(X[p]-ox)*(X[p]-ox)+(Y[p]-oy)*(Y[p]-oy)+(Z[p]-oz)*(Z[p]-oz);
                    if(d<best)best=d; } }
            return best==double.MaxValue? -1.0 : Math.Sqrt(best);
        };

        // sphere_center port; returns false if no radius-rho ball touches a,b,c
        Func<double[],double[],double[],int,double[],bool> sc=(a,b,c,sign,O)=>{
            double[] ab={b[0]-a[0],b[1]-a[1],b[2]-a[2]};
            double[] ac={c[0]-a[0],c[1]-a[1],c[2]-a[2]};
            double[] nn={ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0]};
            double nl2=nn[0]*nn[0]+nn[1]*nn[1]+nn[2]*nn[2]; if(nl2<1e-30)return false;
            double abd=ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2];
            double acd=ac[0]*ac[0]+ac[1]*ac[1]+ac[2]*ac[2];
            double aba=ab[0]*ac[0]+ab[1]*ac[1]+ab[2]*ac[2];
            double D=abd*acd-aba*aba; if(D<1e-30)return false;
            double al=acd*(abd-aba)/(2*D), be=abd*(acd-aba)/(2*D);
            double[] co={a[0]+al*ab[0]+be*ac[0], a[1]+al*ab[1]+be*ac[1], a[2]+al*ab[2]+be*ac[2]};
            double rc2=(co[0]-a[0])*(co[0]-a[0])+(co[1]-a[1])*(co[1]-a[1])+(co[2]-a[2])*(co[2]-a[2]);
            if(rc2>rho*rho)return false;
            double h=Math.Sqrt(rho*rho-rc2), nl=Math.Sqrt(nl2);
            double[] nh={nn[0]/nl,nn[1]/nl,nn[2]/nl}; if(sign<0){nh[0]=-nh[0];nh[1]=-nh[1];nh[2]=-nh[2];}
            O[0]=co[0]+nh[0]*h; O[1]=co[1]+nh[1]*h; O[2]=co[2]+nh[2]*h; return true;
        };

        int nTri=0, tooBigRadius=0, blocked=0, borderline=0, empty=0;
        var sb=new System.Text.StringBuilder();
        int shown=0;
        foreach(var seq in loops){
            if(seq.Count<3||seq.Count>maxLoop) continue;
            if(seq.Count!=3) continue; // only single-triangle holes here
            nTri++;
            int p=seq[0],q=seq[1],r=seq[2];
            double[] A={X[p],Y[p],Z[p]}, B={X[q],Y[q],Z[q]}, C={X[r],Y[r],Z[r]};
            double e0=Math.Sqrt((A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]));
            double e1=Math.Sqrt((B[0]-C[0])*(B[0]-C[0])+(B[1]-C[1])*(B[1]-C[1])+(B[2]-C[2])*(B[2]-C[2]));
            double e2=Math.Sqrt((C[0]-A[0])*(C[0]-A[0])+(C[1]-A[1])*(C[1]-A[1])+(C[2]-A[2])*(C[2]-A[2]));
            double emax=Math.Max(e0,Math.Max(e1,e2));
            double[] Op=new double[3]; double[] Om=new double[3];
            bool hp=sc(A,B,C,1,Op), hm=sc(A,B,C,-1,Om);
            if(!hp&&!hm){ tooBigRadius++;
                if(shown<14){ sb.AppendFormat("  [tri {0,3}] emax={1:F3} circumR>rho => NO radius-{2} ball touches all 3\n",nTri,emax,rho); shown++; }
                continue; }
            // best (most-empty) side
            double bestMargin=-1e9; double bestDmin=-1; string side="";
            if(hp){ double dmin=nearestOther(Op[0],Op[1],Op[2],p,q,r); double m=(dmin<0? 1e9 : dmin-(rho-1e-4)); if(m>bestMargin){bestMargin=m;bestDmin=dmin;side="+";} }
            if(hm){ double dmin=nearestOther(Om[0],Om[1],Om[2],p,q,r); double m=(dmin<0? 1e9 : dmin-(rho-1e-4)); if(m>bestMargin){bestMargin=m;bestDmin=dmin;side="-";} }
            string cls;
            if(bestMargin<0){ blocked++; cls="BLOCKED"; }
            else if(bestMargin<5e-4){ borderline++; cls="BORDERLINE(precision)"; }
            else { empty++; cls="EMPTY(should-close: orient/theta)"; }
            if(shown<14){ sb.AppendFormat("  [tri {0,3}] emax={1:F3} bestSide={2} dmin={3:F4} margin={4:F5} -> {5}\n",
                nTri,emax,side,bestDmin,bestMargin,cls); shown++; }
        }
        var outp=new System.Text.StringBuilder();
        outp.AppendFormat("single-triangle holes probed: {0}  (rho={1})\n",nTri,rho);
        outp.AppendFormat("  too-big-for-rho (circumR>rho) : {0}\n",tooBigRadius);
        outp.AppendFormat("  genuine BLOCKER inside ball    : {0}\n",blocked);
        outp.AppendFormat("  BORDERLINE (float precision)   : {0}\n",borderline);
        outp.AppendFormat("  EMPTY ball (orient/theta rej)  : {0}\n",empty);
        outp.Append("samples:\n").Append(sb.ToString());
        return outp.ToString();
    }
}
'@ -Language CSharp

Write-Host ([Pinhole]::Probe($Obj, $Rho, $MaxLoop))
