/* Narrow phase on original material boundaries. The retained pair's full
 * intersection area supplies the existing solver objective and derivative;
 * it is not the area of the eroded intersection. The optional raw objective
 * also retains allowed pairs, so internal face derivatives cancel correctly
 * and leave forces on the actual material boundary. */
static AcBandWork *ac_work(const AsmContacts *q)
{
    /* A private query may itself be called by an outer repair worker.
     * Only our own leaf-pair walk selects the numbered worker workspace. */
    return q->n_band_work?q->band_work:NULL;
}

static void ac_face_xy(const AsmContacts *q,const double *uv,size_t face,double xy[6])
{
    for(int k=0;k<3;k++)for(int d=0;d<2;d++)xy[2*k+d]=uv[2*(size_t)q->faces[3*face+k]+d];
}

static void ac_pair_gradient(const double a[6],const double b[6],const int ia[3],const int ib[3],double *ga,double *gb)
{
    double g[6];
    if(ga){ac_boundary_gradient(a,b,g);for(int k=0;k<3;k++)for(int d=0;d<2;d++)ga[2*ia[k]+d]=g[2*k+d];}
    if(gb){ac_boundary_gradient(b,a,g);for(int k=0;k<3;k++)for(int d=0;d<2;d++)gb[2*ib[k]+d]=g[2*k+d];}
}

#if defined(_MSC_VER)
#define AC_NOINLINE __declspec(noinline)
#elif defined(__GNUC__)
#define AC_NOINLINE __attribute__((noinline))
#else
#define AC_NOINLINE
#endif

/* The material-boundary test of a positive-area pair, kept out of line: its
 * segment buffer (4*ACB_SEGMENTS doubles) would otherwise make every pair
 * test probe a large stack frame. Returns the band contact verdict, or -1
 * when the boundary gather is incomplete. */
static AC_NOINLINE int ac_material_contact(AcBandWork *work,
    const AcMaterial *one,size_t face_a,const double *uv_a,const double a[6],
    const AcMaterial *two,size_t face_b,const double *uv_b,const double b[6],
    const double origin[2],const double *polygon,int vertices,double tolerance,double roundoff)
{
    double segments[4*ACB_SEGMENTS],radius=.5*tolerance+roundoff;
    AcMaterialGather gather={segments,{DBL_MAX,DBL_MAX,-DBL_MAX,-DBL_MAX},{origin[0],origin[1]},NULL,0,1};
    for(int k=0;k<vertices;k++)for(int d=0;d<2;d++){
        gather.box[d]=ac_min(gather.box[d],polygon[2*k+d]-radius);
        gather.box[2+d]=ac_max(gather.box[2+d],polygon[2*k+d]+radius);
    }
    for(int side=0;side<2;side++){
        const AcMaterial *m=side?two:one;const double *uv=side?uv_b:uv_a,*tri=side?b:a;
        if(m)acm_gather(m,side?face_b:face_a,uv,&gather);
        else for(int k=0;k<3;k++){
            if(gather.count==ACB_SEGMENTS){gather.complete=0;break;}
            double *segment=segments+4*gather.count++;int j=(k+1)%3;
            segment[0]=tri[2*k];segment[1]=tri[2*k+1];segment[2]=tri[2*j];segment[3]=tri[2*j+1];
        }
    }
    if(!gather.complete)return -1;
    return acb_contact(work,polygon,vertices,segments,gather.count,tolerance,roundoff,NULL);
}

static double ac_material_pair(AcBandWork *work,
    const AcMaterial *one,size_t face_a,const double *uv_a,const double in_a[6],
    const AcMaterial *two,size_t face_b,const double *uv_b,const double in_b[6],
    double tolerance,double ga[6],double gb[6],double *raw_area)
{
    if(raw_area)*raw_area=0;
    if(ga)memset(ga,0,6*sizeof(double));if(gb)memset(gb,0,6*sizeof(double));
    double scale=1;
    for(int k=0;k<6;k++){
        if(!isfinite(in_a[k]) || !isfinite(in_b[k]))return NAN;
        scale=ac_max(scale,ac_max(fabs(in_a[k]),fabs(in_b[k])));
    }
    double roundoff=64*DBL_EPSILON*scale;
    double a[6],b[6],origin[2]={in_a[0],in_a[1]},polygon[24];int ia[3],ib[3],vertices=0;
    ac_positive(in_a,a,ia,origin);ac_positive(in_b,b,ib,origin);
    if(!ac_overlap(a,b,roundoff)){
        /* Coincident internal mesh edges have zero-area pairs on one side.
         * Their symmetric boundary derivative completes the force from the
         * positive pairs on the other side. Omitting it makes an aligned
         * fine grid appear stationary despite a changing overlap width. */
        if(raw_area && (ga || gb) && ac_overlap(a,b,-roundoff))ac_pair_gradient(a,b,ia,ib,ga,gb);
        return 0;
    }
    double area=ac_clip_polygon(a,b,polygon,&vertices);
    if(!(area>0) || !isfinite(area))return area;
    if(raw_area)*raw_area=area;
    if(ga || gb)ac_pair_gradient(a,b,ia,ib,ga,gb);
    if(tolerance>0 && !(one && one->strict[one->component[face_a]]) &&
       !(two && two->strict[two->component[face_b]])){
        int contact=ac_material_contact(work,one,face_a,uv_a,a,two,face_b,uv_b,b,origin,polygon,vertices,tolerance,roundoff);
        if(contact<0)return NAN;
        if(!contact){
            if(!raw_area){if(ga)memset(ga,0,6*sizeof(double));if(gb)memset(gb,0,6*sizeof(double));}
            return 0;
        }
    }
    return area;
}

double AsmContacts_pair_faces(const AsmContacts *q,const double *uv,size_t a,size_t b,double ga[6],double gb[6])
{
    if(!q || !uv || a>=q->nf || b>=q->nf)return NAN;
    double x[6],y[6];ac_face_xy(q,uv,a,x);ac_face_xy(q,uv,b,y);
    return ac_material_pair(ac_work(q),q->material,a,uv,x,q->material,b,uv,y,q->tolerance,ga,gb,NULL);
}
