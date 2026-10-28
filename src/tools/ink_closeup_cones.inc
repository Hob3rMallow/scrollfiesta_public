/* Conformal initialization is only a start for SLIM. The final map still
 * receives SLIM's orientation and boundary checks and the publication audit. */
static int ic_reparameterize_mode = 0;

static int ic_cone_key_compare(const void *a,const void *b)
{
    uint64_t x=*(const uint64_t *)a,y=*(const uint64_t *)b;
    return x<y?-1:x>y;
}

static int ic_cone_key_has(const uint64_t *keys,size_t count,uint64_t key)
{
    size_t lo=0,hi=count;
    while (lo<hi) { size_t mid=lo+(hi-lo)/2;if (keys[mid]<key) lo=mid+1;else hi=mid; }
    return lo<count && keys[lo]==key;
}

#include "ink_closeup_cone_feasibility.inc"

static int32_t ic_uv_root(int32_t *p, int32_t v)
{
    while (p[v]!=v) { p[v]=p[p[v]]; v=p[v]; }
    return v;
}

static int ic_conformal_seed(const MeshBinData *m, double *uv)
{
    size_t nv=m->nv,nf=m->nf,cap=21*nf+2*nv,at=0;
    if (nv>INT_MAX/2 || cap>INT_MAX) return -1;
    int32_t *parent=malloc(nv*sizeof *parent), *pins=malloc(2*nv*sizeof *pins);
    int *rows=malloc(cap*sizeof *rows),*cols=malloc(cap*sizeof *cols);
    double *values=malloc(cap*sizeof *values),*rhs=calloc(2*nv,sizeof *rhs);
    double *sum=calloc(9*nv,sizeof *sum);
    uint8_t *fixed=calloc(nv,1);
    int8_t *orientation=calloc(nv,1);
    int rc=-1;
    if (!parent || !pins || !rows || !cols || !values || !rhs || !sum || !fixed || !orientation) goto done;
    for (size_t v=0;v<nv;v++) { parent[v]=(int32_t)v; pins[2*v]=pins[2*v+1]=-1; }
    for (size_t f=0;f<nf;f++) for (int k=1;k<3;k++) {
        int32_t a=ic_uv_root(parent,m->faces[3*f]),b=ic_uv_root(parent,m->faces[3*f+(size_t)k]);
        if (a!=b) parent[b]=a;
    }
    for (size_t v=0;v<nv;v++) {
        int32_t p=ic_uv_root(parent,(int32_t)v); parent[v]=p;
        int32_t a=pins[2*(size_t)p],b=pins[2*(size_t)p+1];
        if (a<0 || uv[2*v]<uv[2*(size_t)a] || (uv[2*v]==uv[2*(size_t)a] && uv[2*v+1]<uv[2*(size_t)a+1])) pins[2*(size_t)p]=(int32_t)v;
        if (b<0 || uv[2*v]>uv[2*(size_t)b] || (uv[2*v]==uv[2*(size_t)b] && uv[2*v+1]>uv[2*(size_t)b+1])) pins[2*(size_t)p+1]=(int32_t)v;
    }
    for (size_t v=0;v<nv;v++) if (parent[v]==(int32_t)v) {
        if (pins[2*v]<0 || pins[2*v]==pins[2*v+1]) goto done;
        fixed[pins[2*v]]=fixed[pins[2*v+1]]=1;
    }
    for (size_t f=0;f<nf;f++) {
        size_t a=m->faces[3*f],b=m->faces[3*f+1],c=m->faces[3*f+2],p=(size_t)parent[a];
        double det=(uv[2*b]-uv[2*a])*(uv[2*c+1]-uv[2*a+1])-(uv[2*b+1]-uv[2*a+1])*(uv[2*c]-uv[2*a]);
        int8_t sign=det<0?-1:det>0?1:0;
        if (!sign || (orientation[p] && orientation[p]!=sign)) goto done;
        orientation[p]=sign;
    }
    for (size_t v=0;v<nv;v++) if (orientation[parent[v]]<0) uv[2*v+1]=-uv[2*v+1];
    for (size_t f=0;f<nf;f++) {
        const int32_t *face=m->faces+3*f;
        double e[3],d[3],cr[3],l2=0,dot=0;
        for (int k=0;k<3;k++) { e[k]=(double)m->verts[3*(size_t)face[1]+k]-m->verts[3*(size_t)face[0]+k]; d[k]=(double)m->verts[3*(size_t)face[2]+k]-m->verts[3*(size_t)face[0]+k]; l2+=e[k]*e[k]; dot+=e[k]*d[k]; }
        for (int k=0;k<3;k++) cr[k]=e[(k+1)%3]*d[(k+2)%3]-e[(k+2)%3]*d[(k+1)%3];
        double len=sqrt(l2),a2=sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
        if (!(len>0 && a2>0)) goto done;
        double along=dot/len,height=a2/len;
        double gx[3]={-1/len,1/len,0},gy[3]={(along/len-1)/height,-along/(len*height),1/height};
        double grad[6][2]; int index[6];
        for (int k=0;k<3;k++) {
            index[2*k]=2*face[k];index[2*k+1]=2*face[k]+1;
            grad[2*k][0]=gx[k];grad[2*k][1]=gy[k];grad[2*k+1][0]=-gy[k];grad[2*k+1][1]=gx[k];
        }
        sum[9*(size_t)parent[face[0]]]+=a2/2;
        for (int i=0;i<6;i++) for (int j=0;j<=i;j++) {
            int a=index[i],b=index[j]; double value=.5*a2*(grad[i][0]*grad[j][0]+grad[i][1]*grad[j][1]);
            if (fixed[a/2] || fixed[b/2]) {
                if (!fixed[a/2]) rhs[a]-=value*uv[b];
                else if (!fixed[b/2]) rhs[b]-=value*uv[a];
                continue;
            }
            rows[at]=a>b?a:b;cols[at]=a>b?b:a;values[at++]=value;
        }
    }
    for (size_t v=0;v<nv;v++) for (int k=0;k<2;k++) {
        /* Pinch-connected fans can have a conformal rotation nullspace. This
         * tiny initialization prior retains their position; SLIM subsequently
         * optimizes the actual, unmodified 3-D metric. */
        double weight=fixed[v]?1.:1e-9;
        rows[at]=cols[at]=(int)(2*v)+k;values[at++]=weight;rhs[2*v+k]+=weight*uv[2*v+k];
    }
    fprintf(stderr,"[ink_closeup] conformal start: %zu vertices / %zu triangles\n",nv,nf);
    double *start=malloc(2*nv*sizeof *start);
    if (!start) goto done;
    memcpy(start,uv,2*nv*sizeof *start);
    rc=Sparse_solve_sym((int)(2*nv),(int)at,rows,cols,values,rhs,start,SPARSE_SPD);
    if (rc) { free(start); goto done; }
    for (size_t f=0;f<nf;f++) {
        size_t a=m->faces[3*f],b=m->faces[3*f+1],c=m->faces[3*f+2],p=(size_t)parent[a];
        sum[9*p+1]+=.5*fabs((start[2*b]-start[2*a])*(start[2*c+1]-start[2*a+1])-(start[2*b+1]-start[2*a+1])*(start[2*c]-start[2*a]));
    }
    for (size_t v=0;v<nv;v++) {
        size_t p=(size_t)parent[v]; double *s=sum+9*p;
        s[2]+=start[2*v];s[3]+=start[2*v+1];s[4]+=uv[2*v];s[5]+=uv[2*v+1];s[6]++;
    }
    for (size_t v=0;v<nv;v++) if (parent[v]==(int32_t)v) for (int k=2;k<6;k++) sum[9*v+k]/=sum[9*v+6];
    for (size_t v=0;v<nv;v++) {
        double *s=sum+9*(size_t)parent[v];double x=start[2*v]-s[2],y=start[2*v+1]-s[3],u=uv[2*v]-s[4],w=uv[2*v+1]-s[5];
        s[7]+=x*u+y*w;s[8]+=x*w-y*u;
    }
    for (size_t v=0;v<nv;v++) {
        double *s=sum+9*(size_t)parent[v],norm=hypot(s[7],s[8]);
        if (!(s[1]>0 && norm>0)) { rc=-1; break; }
        double scale=sqrt(s[0]/s[1]),co=s[7]/norm,si=s[8]/norm,x=start[2*v]-s[2],y=start[2*v+1]-s[3];
        uv[2*v]=s[4]+scale*(co*x-si*y);uv[2*v+1]=s[5]+scale*(si*x+co*y);
    }
    free(start);
done:
    free(parent);free(pins);free(rows);free(cols);free(values);free(rhs);free(sum);free(fixed);free(orientation);
    return rc;
}

/* Use the existing 0.7-radian cone criterion in either curvature direction.
 * Distributed cones include adjacent same-sign vertices above 0.1 radian.
 * Negative cones request an extra fan branch unless the exact single-slit fan
 * receives a strict physical metric/direction certificate. No triangles are
 * excised. Curvature-only legacy callers retain unconditional fan branches. */
static int ic_cone_cut_directed(MeshBinData *m,const SlimRefineDirection *directions,
                                 int32_t **source_map,size_t *cones,size_t *edges,size_t *fiber_cones)
{
    Arena_T arena=Arena_new();
    double *angle=calloc(m->nv,sizeof *angle);
    uint8_t *term=calloc(m->nv,1),*boundary=calloc(m->nv,1);
    uint64_t *keys=malloc(3*m->nf*sizeof *keys);
    int rc=-1;
    *source_map=NULL;*cones=*edges=0;
    if (fiber_cones) *fiber_cones=0;
    if (!arena || !angle || !term || !boundary || !keys) goto done;
    for (size_t f=0;f<m->nf;f++) for (int k=0;k<3;k++) {
        size_t a=m->faces[3*f+k],b=m->faces[3*f+(k+1)%3],c=m->faces[3*f+(k+2)%3];
        keys[3*f+k]=((uint64_t)(a<b?a:b)<<32)|(uint64_t)(a<b?b:a);
        double e[3],d[3],l=0,h=0,dot=0;
        for (int q=0;q<3;q++) { e[q]=(double)m->verts[3*b+q]-m->verts[3*a+q];d[q]=(double)m->verts[3*c+q]-m->verts[3*a+q];l+=e[q]*e[q];h+=d[q]*d[q];dot+=e[q]*d[q]; }
        if (!(l>0 && h>0)) goto done;
        angle[a]+=acos(fmin(1.,fmax(-1.,dot/sqrt(l*h))));
    }
    /* Native edge-key comparator is local to this module. */
    qsort(keys,3*m->nf,sizeof *keys,ic_cone_key_compare);
    for (size_t k=0;k<3*m->nf;) {
        size_t end=k+1;while (end<3*m->nf && keys[end]==keys[k]) end++;
        if (end-k==1) { boundary[keys[k]>>32]=1;boundary[keys[k]&UINT32_MAX]=1; }
        if (end-k>2) goto done;
        k=end;
    }
    for (size_t v=0;v<m->nv;v++) if (angle[v]>0) {
        double deficit=6.283185307179586-angle[v];
        if (!boundary[v] && deficit>=ASM_SOLID_CONE_RAD) term[v]=1;
        if (deficit<=-ASM_SOLID_CONE_RAD) term[v]=2;
    }
    /* Curvature of the same sign in a cone's immediate fan opens with it. */
    for (size_t f=0;f<m->nf;f++) {
        int positive=0,negative=0;
        for (int k=0;k<3;k++) {
            size_t v=m->faces[3*f+k];double deficit=6.283185307179586-angle[v];
            positive|=term[v]==1 && deficit>=ASM_SOLID_CONE_RAD;
            negative|=term[v]==2 && deficit<=-ASM_SOLID_CONE_RAD;
        }
        for (int k=0;k<3;k++) {
            size_t v=m->faces[3*f+k];double deficit=6.283185307179586-angle[v];
            if (!boundary[v] && positive && deficit>=ASM_SOLID_CONE_NODE_RAD) term[v]=1;
            if (!boundary[v] && negative && deficit<=-ASM_SOLID_CONE_NODE_RAD) term[v]=2;
        }
    }
    size_t field_count=0;
    if (directions && SeamCut_direction_terminals(arena,m->verts,m->nv,m->faces,m->nf,
                                                  directions,term,&field_count)) goto done;
    if (fiber_cones) *fiber_cones=field_count;
    /* A field singularity may expose a mildly negative-curvature fan below
     * the curvature-only threshold. Test its physical overfull sector before
     * deciding whether it needs a second branch; flat/positive sectors stay. */
    if (directions) for (size_t v=0;v<m->nv;v++)
        if (term[v] && angle[v]>6.283185307179586) term[v]=2;
    for (size_t v=0;v<m->nv;v++) *cones+=term[v]!=0;
    if (!*cones) { rc=0;goto done; }
    SeamCutResult cut={0};
    IcConeFanContext check={0};check.arena=arena;check.mesh=m;check.directions=directions;
    int cutting=SeamCut_to_boundary_checked(arena,m->verts,m->nv,m->faces,m->nf,term,
                                           directions?ic_cone_fan_feasible:NULL,&check,&cut);
    if (directions) ic_cone_fan_report(&check);
    if (cutting) goto done;
    *edges=(size_t)cut.seam_edges;
    size_t cap=m->nv+cut.nv,nv=m->nv;
    int32_t *first=malloc(m->nv*sizeof *first),*new_id=malloc(cut.nv*sizeof *new_id),*mapping=malloc(cap*sizeof *mapping);
    MeshBinData out={0};out.nf=m->nf;
    out.verts=malloc(3*cap*sizeof(float));out.uv64=malloc(2*cap*sizeof(double));out.faces=malloc(3*m->nf*sizeof(int32_t));
    if (!first || !new_id || !mapping || !out.verts || !out.uv64 || !out.faces) {
        free(first);free(new_id);free(mapping);MeshBin_dispose(&out);goto done;
    }
    for (size_t v=0;v<m->nv;v++) { first[v]=-1;mapping[v]=(int32_t)v; }
    for (size_t v=0;v<cut.nv;v++) {
        int32_t source=cut.vmap[v];
        if (first[source]<0) first[source]=new_id[v]=source;
        else { new_id[v]=(int32_t)nv;mapping[nv++]=source; }
    }
    out.nv=nv;
    for (size_t v=0;v<nv;v++) {
        size_t s=(size_t)mapping[v];memcpy(out.verts+3*v,m->verts+3*s,3*sizeof(float));
        out.uv64[2*v]=MeshBin_uv(m,2*s);out.uv64[2*v+1]=MeshBin_uv(m,2*s+1);
    }
    for (size_t k=0;k<3*m->nf;k++) out.faces[k]=new_id[cut.faces[k]];
    /* Open each bank into the half-planes of its own boundary faces. A fan's
     * centroid can lie outside one of those half-planes near a reflex corner,
     * making newly duplicated slit banks cross even for tiny displacements. */
    double alpha=0;
    int opening=ic_cone_open_banks(arena,&out,mapping,m->nv,&alpha);
    if (opening) { free(first);free(new_id);free(mapping);MeshBin_dispose(&out);goto done; }
    fprintf(stderr,"[ink_closeup] certified slit-bank opening alpha=%.9g\n",alpha);
    fprintf(stderr,"[ink_closeup] %zu cone terminals / %zu seam edges / %zu duplicate vertices; all %zu triangles retained\n",*cones,*edges,nv-m->nv,m->nf);
    MeshBin_dispose(m);*m=out;*source_map=mapping;free(first);free(new_id);rc=0;
done:
    free(angle);free(term);free(boundary);free(keys);Arena_dispose(&arena);return rc;
}

static int ic_cone_cut(MeshBinData *m,int32_t **source_map,size_t *cones,size_t *edges)
{
    return ic_cone_cut_directed(m,NULL,source_map,cones,edges,NULL);
}
