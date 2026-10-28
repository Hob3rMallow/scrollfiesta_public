static int sc_direction_phase(const SlimRefineDirection *d,const double normal[3],
                               const double ray[3],double *phase)
{
    if (!(d->weight>0) || !isfinite(d->weight)) return 0;
    double across[3],real=0,imaginary=0;
    for (int k=0;k<3;k++) across[k]=normal[(k+1)%3]*ray[(k+2)%3]-normal[(k+2)%3]*ray[(k+1)%3];
    for (int arm=0;arm<2;arm++) {
        const double *w=arm?d->v:d->h;
        double x=vdot(w,ray),y=vdot(w,across),length=hypot(x,y),norm=vlen(w);
        if (!isfinite(norm) || !(length>1e-10*norm)) return 0;
        double a=4*atan2(y,x);real+=cos(a);imaginary+=sin(a);
    }
    if (!(hypot(real,imaginary)>1e-8)) return 0;
    *phase=.25*atan2(imaginary,real);return 1;
}

/* Walk one complete, oriented interior fan. Disconnected pinched fans,
 * inconsistent winding, missing evidence and ambiguous quarter-turn jumps
 * cannot establish an interior cross-field index. */
static int sc_direction_index(const Conn *c,const double *normal,
                               const SlimRefineDirection *directions,size_t vertex,
                               long *quarter_index)
{
    int32_t begin=c->voff[vertex],end=c->voff[vertex+1];
    if (begin==end) return 0;
    int32_t start=c->vcor[begin],corner=start;double rotation=0,angles=0;
    size_t count=0,degree=(size_t)(end-begin);
    do {
        int32_t twin=c->twin[corner];
        if (twin<0 || c->faces[twin]!=c->faces[NEXT(corner)] ||
            c->faces[NEXT(twin)]!=(int32_t)vertex) return 0;
        double ray[3],a=0,b=0;
        vsub(c->verts+3*(size_t)c->faces[NEXT(corner)],c->verts+3*vertex,ray);
        double length=vlen(ray);if (!(length>0)) return 0;
        for (int k=0;k<3;k++) ray[k]/=length;
        if (!sc_direction_phase(directions+corner/3,normal+3*(size_t)(corner/3),ray,&a) ||
            !sc_direction_phase(directions+twin/3,normal+3*(size_t)(twin/3),ray,&b)) return 0;
        double delta=.25*atan2(sin(4*(a-b)),cos(4*(a-b)));
        if (fabs(fabs(delta)-.25*M_PI)<1e-9) return 0;
        rotation+=delta;
        angles+=corner_angle(c->verts+3*vertex,c->verts+3*(size_t)c->faces[NEXT(corner)],
                             c->verts+3*(size_t)c->faces[PREV(corner)]);
        corner=NEXT(twin);count++;
    } while (corner!=start && count<degree);
    if (corner!=start || count!=degree) return 0;
    double quarters=(rotation+2*M_PI-angles)/(.5*M_PI);
    if (!isfinite(quarters) || fabs(quarters)>INT_MAX) return 0;
    long index=lround(quarters);
    if (fabs(quarters-(double)index)>1e-7*(double)degree) return 0;
    *quarter_index=index;return 1;
}

int SeamCut_direction_terminals(Arena_T arena,const float *verts,size_t nv,
                                 const int32_t *faces,size_t nf,
                                 const SlimRefineDirection *directions,
                                 uint8_t *terminals,size_t *count)
{
    if (count) *count=0;
    if (!arena || !verts || !faces || !directions || !terminals ||
        nv>(size_t)INT_MAX || nf>(size_t)INT_MAX/3) return -1;
    if (!nv || !nf) return 0;
    for (size_t k=0;k<3*nf;k++) if (faces[k]<0 || (size_t)faces[k]>=nv) return -1;
    Arena_Mark mark=Arena_save(arena);Conn c;conn_build(arena,verts,nv,faces,nf,&c);
    double *normal=ARENA_ALLOC(arena,3*nf*sizeof *normal);int rc=-1;
    size_t eligible=0,positive=0,negative=0;
    if (c.nonmanifold) goto done;
    for (size_t f=0;f<nf;f++) {
        double e[3],g[3];vsub(verts+3*(size_t)faces[3*f+1],verts+3*(size_t)faces[3*f],e);
        vsub(verts+3*(size_t)faces[3*f+2],verts+3*(size_t)faces[3*f],g);
        for (int k=0;k<3;k++) normal[3*f+k]=e[(k+1)%3]*g[(k+2)%3]-e[(k+2)%3]*g[(k+1)%3];
        double length=vlen(normal+3*f);if (!(length>0) || !isfinite(length)) goto done;
        for (int k=0;k<3;k++) normal[3*f+k]/=length;
    }
    for (size_t v=0;v<nv;v++) {
        long index=0;
        if (!sc_direction_index(&c,normal,directions,v,&index)) continue;
        eligible++;
        if (index) { if (!terminals[v]) terminals[v]=1;positive+=index>0;negative+=index<0; }
    }
    if (count) *count=positive+negative;
    fprintf(stderr,"  [seamcut fiber] supported interior vertices=%zu positive=%zu negative=%zu\n",eligible,positive,negative);
    rc=0;
done:
    Arena_restore(arena,mark);return rc;
}
