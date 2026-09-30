/* Continuous material-width predicate. An overlap is within the tolerance
 * iff every point is within half that tolerance of an ORIGINAL boundary.
 * Internal triangulation edges never enter this query. A segment's band is
 * its rectangle plus two endpoint disks. Subtract the rectangles exactly;
 * maxima of distance to the remaining disk centers occur at polygon
 * vertices, edge/bisector intersections or Voronoi vertices. No pixel grid,
 * finite probe sampling, or polygonal approximation of circles certifies it.
 * Storage/work exhaustion returns -1, never a false clear result. */
#define ACB_VERTICES 32
#define ACB_PIECES 128
#define ACB_SEGMENTS 128
#define ACB_POINTS (2*ACB_SEGMENTS)
#define ACB_CANDIDATES 65536u

typedef struct AcBandPoly { double xy[2*ACB_VERTICES]; int n; } AcBandPoly;
typedef struct AcBandWork {
    AcBandPoly first[ACB_PIECES],second[ACB_PIECES],temp[3];
    double point[2*ACB_POINTS];
    size_t candidates;
} AcBandWork;

static double acb_segment_distance2(const double p[2],const double segment[4])
{
    double x=segment[2]-segment[0],y=segment[3]-segment[1];
    double a=p[0]-segment[0],b=p[1]-segment[1],length=x*x+y*y;
    double t=length>0?ac_max(0,ac_min(1,(a*x+b*y)/length)):0;
    a-=t*x;b-=t*y;return a*a+b*b;
}

static int acb_uncovered(const double p[2],const double *segments,size_t count,double radius2)
{
    for(size_t i=0;i<count;i++)if(acb_segment_distance2(p,segments+4*i)<=radius2)return 0;
    return 1;
}

static double acb_area(const AcBandPoly *p)
{
    double area=0;
    for(int i=1;i+1<p->n;i++)area+=ac_cross(p->xy[2*i]-p->xy[0],p->xy[2*i+1]-p->xy[1],
        p->xy[2*i+2]-p->xy[0],p->xy[2*i+3]-p->xy[1]);
    return .5*fabs(area);
}

static int acb_append(AcBandPoly *p,double x,double y)
{
    if(p->n && p->xy[2*p->n-2]==x && p->xy[2*p->n-1]==y)return 1;
    if(p->n==ACB_VERTICES)return 0;
    p->xy[2*p->n]=x;p->xy[2*p->n+1]=y;p->n++;return 1;
}

/* Keep nx*x+ny*y >= height. */
static int acb_clip(const AcBandPoly *input,double nx,double ny,double height,AcBandPoly *out)
{
    out->n=0;
    for(int i=0;i<input->n;i++){
        int j=(i+1)%input->n;const double *a=input->xy+2*i,*b=input->xy+2*j;
        double da=nx*a[0]+ny*a[1]-height,db=nx*b[0]+ny*b[1]-height;
        if(da>=0 && !acb_append(out,a[0],a[1]))return 0;
        if((da<0 && db>0) || (da>0 && db<0)){
            double t=da/(da-db);
            if(!acb_append(out,a[0]+t*(b[0]-a[0]),a[1]+t*(b[1]-a[1])))return 0;
        }
    }
    if(out->n>1 && out->xy[0]==out->xy[2*out->n-2] && out->xy[1]==out->xy[2*out->n-1])out->n--;
    return 1;
}

static int acb_inside(const AcBandPoly *p,const double at[2],double roundoff)
{
    for(int i=0;i<p->n;i++){
        int j=(i+1)%p->n;double x=p->xy[2*j]-p->xy[2*i],y=p->xy[2*j+1]-p->xy[2*i+1];
        if(ac_cross(x,y,at[0]-p->xy[2*i],at[1]-p->xy[2*i+1]) < -roundoff*hypot(x,y))return 0;
    }
    return 1;
}

static int acb_witness(AcBandWork *work,const AcBandPoly *poly,const double p[2],int points,
                       const double *segments,size_t count,double radius2,double roundoff,double witness[2])
{
    if(++work->candidates>ACB_CANDIDATES)return -1;
    if(!isfinite(p[0]) || !isfinite(p[1]))return -1;
    if(!acb_inside(poly,p,roundoff))return 0;
    double nearest=DBL_MAX;
    for(int i=0;i<points;i++){
        double x=p[0]-work->point[2*i],y=p[1]-work->point[2*i+1];nearest=ac_min(nearest,x*x+y*y);
    }
    if(nearest<=radius2)return 0;
    if(acb_uncovered(p,segments,count,radius2)){if(witness)memcpy(witness,p,2*sizeof(double));return 1;}
    /* Rectangle subtraction leaves closed polygon pieces. A maximizing
     * vertex can lie exactly ON a removed rectangle even though its open
     * neighborhood is uncovered. Move toward this piece's interior by less
     * than half the endpoint-disk clearance, then validate against the exact
     * original segments. Testing only the vertex would silently miss the
     * eroded interior of a triangle cut by an internal triangulation edge. */
    double center[2]={0};
    for(int i=0;i<poly->n;i++){center[0]+=poly->xy[2*i]/poly->n;center[1]+=poly->xy[2*i+1]/poly->n;}
    double length=hypot(center[0]-p[0],center[1]-p[1]);
    double step=points && length>0?ac_min(.5,.25*(sqrt(nearest)-sqrt(radius2))/length):.5;
    double inside[2]={p[0]+step*(center[0]-p[0]),p[1]+step*(center[1]-p[1])};
    if(!acb_uncovered(inside,segments,count,radius2))return -1;
    if(witness)memcpy(witness,inside,sizeof inside);return 1;
}

/* Complete complement of the endpoint disks inside one remaining convex
 * piece. Collinear centers need only the polygon-edge candidates. */
static int acb_disk_complement(AcBandWork *work,const AcBandPoly *poly,int points,
                                const double *segments,size_t count,double radius2,double roundoff,double witness[2])
{
    for(int i=0;i<poly->n;i++){
        int found=acb_witness(work,poly,poly->xy+2*i,points,segments,count,radius2,roundoff,witness);
        if(found)return found;
    }
    /* One disk containing every vertex contains this complete convex piece. */
    for(int i=0;i<points;i++){
        int inside=1;const double *p=work->point+2*i;
        for(int j=0;j<poly->n && inside;j++){
            double x=poly->xy[2*j]-p[0],y=poly->xy[2*j+1]-p[1];inside=x*x+y*y<=radius2;
        }
        if(inside)return 0;
    }
    for(int i=0;i<points;i++)for(int j=i+1;j<points;j++){
        const double *a=work->point+2*i,*b=work->point+2*j;
        double nx=b[0]-a[0],ny=b[1]-a[1],height=.5*(nx*nx+ny*ny);
        if(!(height>0))continue;
        for(int k=0;k<poly->n;k++){
            const double *start=poly->xy+2*k,*end=poly->xy+2*((k+1)%poly->n);
            double denominator=nx*(end[0]-start[0])+ny*(end[1]-start[1]);
            if(denominator==0)continue;
            double t=(height-nx*(start[0]-a[0])-ny*(start[1]-a[1]))/denominator;
            if(t<0 || t>1)continue;
            double p[2]={start[0]+t*(end[0]-start[0]),start[1]+t*(end[1]-start[1])};
            int found=acb_witness(work,poly,p,points,segments,count,radius2,roundoff,witness);if(found)return found;
        }
        for(int k=j+1;k<points;k++){
            const double *c=work->point+2*k;
            double x=c[0]-a[0],y=c[1]-a[1],denominator=2*ac_cross(nx,ny,x,y);
            if(denominator==0)continue;
            double ab=nx*nx+ny*ny,ac=x*x+y*y;
            double p[2]={a[0]+(ab*y-ac*ny)/denominator,a[1]+(ac*nx-ab*x)/denominator};
            int found=acb_witness(work,poly,p,points,segments,count,radius2,roundoff,witness);if(found)return found;
        }
    }
    return 0;
}

/* polygon and segments are recentered in the same frame. The caller supplies
 * a scale-aware floating-point allowance, not a new physical tolerance. */
static int acb_contact(AcBandWork *work,const double *polygon,int vertices,
                       const double *segments,size_t count,double tolerance,double roundoff,double witness[2])
{
    if(!work || !polygon || vertices<3 || vertices>ACB_VERTICES || (count && !segments) || count>ACB_SEGMENTS ||
       !isfinite(tolerance) || tolerance<0 || !isfinite(roundoff) || roundoff<0)return -1;
    work->candidates=0;AcBandPoly *a=work->first,*b=work->second;
    a[0].n=vertices;memcpy(a[0].xy,polygon,2*(size_t)vertices*sizeof(double));
    for(int k=0;k<2*vertices;k++)if(!isfinite(polygon[k]))return -1;
    for(size_t k=0;k<4*count;k++)if(!isfinite(segments[k]))return -1;
    if(!(acb_area(a)>0))return 0;
    double radius=.5*tolerance+roundoff,radius2=radius*radius,center[2]={0};
    if(!isfinite(radius2))return -1;
    for(int k=0;k<vertices;k++){center[0]+=polygon[2*k]/vertices;center[1]+=polygon[2*k+1]/vertices;}
    if(acb_uncovered(center,segments,count,radius2)){if(witness)memcpy(witness,center,sizeof center);return 1;}
    for(int k=0;k<vertices;k++)if(acb_uncovered(polygon+2*k,segments,count,radius2)){
        if(witness)memcpy(witness,polygon+2*k,2*sizeof(double));return 1;
    }
    for(size_t i=0;i<count;i++){
        int inside=1;
        for(int k=0;k<vertices && inside;k++)inside=acb_segment_distance2(polygon+2*k,segments+4*i)<=radius2;
        if(inside)return 0;
    }
    int pieces=1;
    for(size_t i=0;i<count && pieces;i++){
        const double *s=segments+4*i;double nx=s[2]-s[0],ny=s[3]-s[1],length=hypot(nx,ny);
        if(!(length>0))continue;nx/=length;ny/=length;
        double ax=nx*s[0]+ny*s[1],ay=-ny*s[0]+nx*s[1];int next=0;
        const double planes[12]={nx,ny,ax,-nx,-ny,-ax-length,-ny,nx,ay-radius,ny,-nx,-ay-radius};
        for(int j=0;j<pieces;j++){
            AcBandPoly *current=work->temp,*inside=work->temp+1,*outside=work->temp+2;*current=a[j];
            for(int k=0;k<4 && current->n>=3;k++){
                const double *h=planes+3*k;
                if(!acb_clip(current,-h[0],-h[1],-h[2],outside) || !acb_clip(current,h[0],h[1],h[2],inside))return -1;
                if(outside->n>=3 && acb_area(outside)>0){if(next==ACB_PIECES)return -1;b[next++]=*outside;}
                *current=*inside;
            }
        }
        pieces=next;AcBandPoly *swap=a;a=b;b=swap;
    }
    if(!pieces)return 0;
    double lo[2]={DBL_MAX,DBL_MAX},hi[2]={-DBL_MAX,-DBL_MAX};
    for(int i=0;i<vertices;i++)for(int d=0;d<2;d++){lo[d]=ac_min(lo[d],polygon[2*i+d]);hi[d]=ac_max(hi[d],polygon[2*i+d]);}
    int points=0;
    for(size_t i=0;i<count;i++)for(int end=0;end<2;end++){
        const double *p=segments+4*i+2*end;int exists=0;
        double dx=ac_max(0,ac_max(lo[0]-p[0],p[0]-hi[0])),dy=ac_max(0,ac_max(lo[1]-p[1],p[1]-hi[1]));
        if(dx*dx+dy*dy>radius2)continue;
        for(int j=0;j<points;j++)if(work->point[2*j]==p[0] && work->point[2*j+1]==p[1]){exists=1;break;}
        if(!exists){if(points==ACB_POINTS)return -1;memcpy(work->point+2*points++,p,2*sizeof(double));}
    }
    for(int i=0;i<pieces;i++){
        int found=acb_disk_complement(work,a+i,points,segments,count,radius2,roundoff,witness);if(found)return found;
    }
    return 0;
}

static int acb_selftest(void)
{
    Arena_T arena=Arena_new();AcBandWork *work=ARENA_ALLOC(arena,sizeof *work);int failures=0;
    /* A material triangle of leg 8 has inscribed diameter 4.686. Cutting
     * its overlap along x=y must retain that contact; eroded-piece vertices
     * lie exactly on rectangle boundaries and need an interior witness. */
    {
        double half[6]={8,8,16,8,12,12},segments[12]={8,8,16,8,16,8,8,16,8,16,8,8},witness[2];
        if(acb_contact(work,half,3,segments,3,4,1e-10,witness)!=1 || !acb_uncovered(witness,segments,3,4))failures++;
    }
    /* Two squares overlap in strips of known physical width. Subdividing
     * either mesh changes none of these original boundary segments. */
    for(int pose=0;pose<2;pose++)for(int width=0;width<4;width++)for(int fine=0;fine<2;fine++){
        double depth=width==0?3.5:width==1?4:width==2?4.0001:6;
        double angle=pose?.37:0,cs=cos(angle),sn=sin(angle),segments[32],base[32]={
            0,0,20,0,20,0,20,20,20,20,0,20,0,20,0,0,
            20-depth,0,40-depth,0,40-depth,0,40-depth,20,40-depth,20,20-depth,20,20-depth,20,20-depth,0};
        for(int i=0;i<16;i++){segments[2*i]=cs*base[2*i]-sn*base[2*i+1];segments[2*i+1]=sn*base[2*i]+cs*base[2*i+1];}
        int hit=0,complete=1,cells=fine?20:1;
        for(int j=0;j<cells;j++)for(int i=0;i<cells;i++){
            double lo=20-depth+depth*i/cells,hi=20-depth+depth*(i+1)/cells,y0=20.*j/cells,y1=20.*(j+1)/cells;
            double box[8]={lo,y0,hi,y0,hi,y1,lo,y1},poly[8],witness[2];
            for(int k=0;k<4;k++){poly[2*k]=cs*box[2*k]-sn*box[2*k+1];poly[2*k+1]=sn*box[2*k]+cs*box[2*k+1];}
            int result=acb_contact(work,poly,4,segments,8,4,1e-10,witness);complete&=result>=0;hit|=result>0;
            if(result>0 && !acb_uncovered(witness,segments,8,(2+1e-10)*(2+1e-10)))failures++;
        }
        if(!complete || hit!=(width>=2)){failures++;fprintf(stderr,"  boundary band pose %d width %.7g fine %d: hit %d complete %d\n",pose,depth,fine,hit,complete);}
    }
    /* Vertices and centroid are covered, while internal Voronoi vertices
     * contain holes. Increasing the radius covers the entire square using
     * a union of disks; no individual disk covers it. */
    {
        double poly[8]={0,0,6,0,6,6,0,6},segments[36],witness[2];
        for(int j=0;j<3;j++)for(int i=0;i<3;i++){
            size_t k=(size_t)j*3+i;segments[4*k]=segments[4*k+2]=3*i;segments[4*k+1]=segments[4*k+3]=3*j;
        }
        if(acb_contact(work,poly,4,segments,9,4.1,1e-10,witness)!=1 ||
           !acb_uncovered(witness,segments,9,2.05*2.05) || witness[0]<=0 || witness[0]>=6 || witness[1]<=0 || witness[1]>=6)failures++;
        if(acb_contact(work,poly,4,segments,9,4.25,1e-10,NULL)!=0)failures++;
    }
    {
        double poly[8]={-3,-2,3,-2,3,2,-3,2},segments[16]={-3,-2,-3,-2,3,-2,3,-2,3,2,3,2,-3,2,-3,2},witness[2];
        if(acb_contact(work,poly,4,segments,4,4,1e-10,witness)!=1)failures++;
        double small[8]={-.8,-.8,.8,-.8,.8,.8,-.8,.8},disk[4]={0,0,0,0};
        if(acb_contact(work,small,4,disk,1,4,1e-10,NULL)!=0)failures++;
        if(acb_contact(work,small,4,NULL,0,4,1e-10,NULL)!=1 || acb_contact(work,small,4,disk,ACB_SEGMENTS+1,4,0,NULL)!=-1)failures++;
    }
    fprintf(stderr,"  original-boundary material band: %s (%d failures)\n",failures?"FAIL":"ok",failures);
    Arena_dispose(&arena);return failures;
}
