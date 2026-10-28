/* Compact a disjoint chunk on disk without allocating packed attribute or
 * face arrays. Retained vertices must be numbered in original order. */
int MeshBin_stream_append_mapped(MeshBinStream_T w, const MeshBinData *m,
                                 const int32_t *map, size_t nv)
{
    union { float xyz[3072]; float uv[2048]; double uv64[2048]; int32_t face[3072]; } buffer;
    size_t used = 0;
    if (!w || !m || !w->f || w->failed) return -1;
    const void *uv = w->uv_size == sizeof(double) ? (const void *)m->uv64 : (const void *)m->uv;
    if (m->nv > (size_t)INT32_MAX || m->nf > SIZE_MAX/3 || nv > m->nv ||
        nv > w->nv-w->at_v || m->nf > w->nf-w->at_f || (m->nv && !map) ||
        (nv && (!m->verts || (w->uvoff && !uv))) || (m->nf && !m->faces)) goto fail;
    for (size_t v = 0; v < m->nv; v++) {
        if (map[v] == -1) continue;
        if (map[v] < 0 || (size_t)map[v] != used) goto fail;
        used++;
    }
    if (used != nv) goto fail;
    for (size_t i = 0; i < 3*m->nf; i++)
        if (m->faces[i] < 0 || (size_t)m->faces[i] >= m->nv || map[m->faces[i]] < 0) goto fail;
    /* Most charts have no orphan vertices: borrow their original arrays. */
    if (nv == m->nv) return MeshBin_stream_append(w,m);
    if (!w->ff && mb_seek(w->f,MB_HEADER_BYTES+w->at_v*12)) goto fail;
    used = 0;
    for (size_t v = 0; v < m->nv; v++) if (map[v] >= 0) {
        memcpy(buffer.xyz+3*used,m->verts+3*v,3*sizeof(float));
        if (++used == 1024) {
            if (!mb_write_bytes(w->f,buffer.xyz,used*3*sizeof(float))) goto fail;
            used = 0;
        }
    }
    if (!mb_write_bytes(w->f,buffer.xyz,used*3*sizeof(float))) goto fail;
    if (w->uvoff) {
        FILE *fu = w->fu ? w->fu : w->f;
        if (!w->fu && mb_seek(fu,w->uvoff+w->at_v*2*w->uv_size)) goto fail;
        used = 0;
        for (size_t v = 0; v < m->nv; v++) if (map[v] >= 0) {
            if (w->uv_size == sizeof(double)) memcpy(buffer.uv64+2*used,m->uv64+2*v,2*sizeof(double));
            else memcpy(buffer.uv+2*used,m->uv+2*v,2*sizeof(float));
            if (++used == 1024) {
                if (!mb_write_bytes(fu,&buffer,used*2*w->uv_size)) goto fail;
                used = 0;
            }
        }
        if (!mb_write_bytes(fu,&buffer,used*2*w->uv_size)) goto fail;
    }
    FILE *ff = w->ff ? w->ff : w->f;
    if (!w->ff && mb_seek(ff,w->foff+w->at_f*12)) goto fail;
    for (size_t i = 0; i < 3*m->nf; ) {
        size_t n = 3*m->nf-i;
        if (n > 3072) n = 3072;
        for (size_t k = 0; k < n; k++) buffer.face[k] = map[m->faces[i+k]]+(int32_t)w->at_v;
        if (!mb_write_bytes(ff,buffer.face,n*sizeof(int32_t))) goto fail;
        i += n;
    }
    w->at_v += nv; w->at_f += m->nf;
    return 0;
fail:
    w->failed = 1;
    return -1;
}

static int mb_mapped_files_equal(const char *a, const char *b)
{
    FILE *fa = fopen(a,"rb"), *fb = fopen(b,"rb");
    int same = fa && fb;
    unsigned char x[4096], y[4096];
    while (same) {
        size_t nx = fread(x,1,sizeof x,fa), ny = fread(y,1,sizeof y,fb);
        if (nx != ny || memcmp(x,y,nx) || ferror(fa) || ferror(fb)) { same = 0; break; }
        if (!nx) break;
    }
    if (fa) fclose(fa);
    if (fb) fclose(fb);
    return same;
}

static int mb_mapped_selftest(void)
{
    const char *path = "mesh_bin_mapped_selftest.vmesh", *expected = "mesh_bin_mapped_expected.vmesh";
    float xyz[15] = {99,99,99, 1,2,3, 4,5,6, 99,99,99, 7,8,9};
    float compact[18] = {1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9};
    float uv[10] = {99,99, -0.0f,2, 3,4, 99,99, 5,6}, cuv[12] = {-0.0f,2,3,4,5,6, -0.0f,2,3,4,5,6};
    double uv64[10] = {99,99, -0.0,1e12+0.125, 3,4, 99,99, 5,6};
    double cuv64[12] = {-0.0,1e12+0.125,3,4,5,6, -0.0,1e12+0.125,3,4,5,6};
    int32_t faces[6] = {1,2,4, 4,2,1}, cf[12] = {0,1,2, 2,1,0, 3,4,5, 5,4,3};
    int32_t map[5] = {-1,0,1,-1,2};
    int fails = 0;
    Arena_T arena = Arena_new_sized(16384);
    for (int mode = 0; mode < 3; mode++) {
        MeshBinStream_T w = NULL;
        MeshBinData m = {xyz,uv,faces,5,2,uv64}, first = {compact,cuv,cf,3,2,cuv64}, empty = {0};
        int rc = mode == 2 ? MeshBin_write_uv64(expected,compact,6,cf,4,cuv64) :
                            MeshBin_write(expected,compact,6,cf,4,mode ? cuv : NULL);
        rc |= mode == 2 ? MeshBin_stream_open_uv64(arena,path,6,4,&w) : MeshBin_stream_open(arena,path,6,4,mode,&w);
        if (!rc) rc = MeshBin_stream_append_mapped(w,&empty,NULL,0) || MeshBin_stream_append(w,&first) ||
                       MeshBin_stream_append_mapped(w,&m,map,3);
        if (w && MeshBin_stream_close(&w,!rc)) rc = -1;
        if (rc || !mb_mapped_files_equal(path,expected)) fails++;
        for (int bad = 0; bad < 7; bad++) {
            size_t nv = 3;
            m.nv = 5; m.nf = 2; map[2] = 1; map[4] = 2; faces[0] = 1;
            if (bad == 0) map[2] = 2;               /* out of order / duplicated */
            if (bad == 1) { map[2] = -1; map[4] = 1; nv = 2; } /* referenced vertex dropped */
            if (bad == 2) nv = 2;                   /* declared count differs */
            if (bad == 3) faces[0] = 5;
            if (bad == 4) faces[0] = -1;
            if (bad == 5) m.nv = SIZE_MAX;
            if (bad == 6) m.nf = SIZE_MAX;
            rc = mode == 2 ? MeshBin_stream_open_uv64(arena,path,6,4,&w) : MeshBin_stream_open(arena,path,6,4,mode,&w);
            if (rc || MeshBin_stream_append_mapped(w,&m,map,nv) == 0) fails++;
            if (w) MeshBin_stream_close(&w,0);
            if (!mb_mapped_files_equal(path,expected)) fails++;
        }
        map[2] = 1; map[4] = 2; faces[0] = 1;
        Arena_free(arena);
    }
    /* Cross every fixed-buffer boundary, including an exact full face block
     * and a final partial attribute block. Expected arrays are contiguous. */
    {
        const size_t nv = 4099, kept = 2050, nf = 2048;
        float *big = ARENA_ALLOC(arena,3*nv*sizeof *big), *packed = ARENA_ALLOC(arena,3*kept*sizeof *packed);
        double *bu = ARENA_ALLOC(arena,2*nv*sizeof *bu), *pu = ARENA_ALLOC(arena,2*kept*sizeof *pu);
        float *bf = ARENA_ALLOC(arena,2*nv*sizeof *bf), *pf = ARENA_ALLOC(arena,2*kept*sizeof *pf);
        int32_t *bm = ARENA_ALLOC(arena,nv*sizeof *bm), *tri = ARENA_ALLOC(arena,3*nf*sizeof *tri), *pt = ARENA_ALLOC(arena,3*nf*sizeof *pt);
        for (size_t v = 0; v < nv; v++) {
            bm[v] = v%2 ? -1 : (int32_t)(v/2);
            for (size_t d = 0; d < 3; d++) big[3*v+d] = (float)(v*3+d)*0.25f;
            for (size_t d = 0; d < 2; d++) { bu[2*v+d] = 1e12+(double)(2*v+d)*0.125; bf[2*v+d] = (float)(2*v+d)*0.125f; }
            if (v%2 == 0) {
                memcpy(packed+3*(v/2),big+3*v,3*sizeof *big);
                memcpy(pu+2*(v/2),bu+2*v,2*sizeof *bu); memcpy(pf+2*(v/2),bf+2*v,2*sizeof *bf);
            }
        }
        for (size_t t = 0; t < nf; t++) for (size_t d = 0; d < 3; d++) {
            pt[3*t+d] = (int32_t)(t+d); tri[3*t+d] = 2*pt[3*t+d];
        }
        MeshBinData m = {big,bf,tri,nv,nf,bu};
        for (int mode = 0; mode < 3; mode++) {
            MeshBinStream_T w = NULL;
            int rc = mode == 2 ? MeshBin_write_uv64(expected,packed,kept,pt,nf,pu) : MeshBin_write(expected,packed,kept,pt,nf,mode ? pf : NULL);
            rc |= mode == 2 ? MeshBin_stream_open_uv64(arena,path,kept,nf,&w) : MeshBin_stream_open(arena,path,kept,nf,mode,&w);
            if (!rc) rc = MeshBin_stream_append_mapped(w,&m,bm,kept);
            if (w && MeshBin_stream_close(&w,!rc)) rc = -1;
            if (rc || !mb_mapped_files_equal(path,expected)) fails++;
        }
    }
    FILE *temp = fopen("mesh_bin_mapped_selftest.vmesh.tmp","rb");
    if (temp) { fails++; fclose(temp); }
    Arena_dispose(&arena);
    if (!fails) { remove(path); remove(expected); }
    fprintf(stderr,"MeshBin mapped stream: exact geometry/UV32/UV64 bytes, rebasing, invalid-map refusal; %d failure(s)\n",fails);
    return fails;
}
