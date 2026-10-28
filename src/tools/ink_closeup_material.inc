/* A per-component simple boundary does not certify placement among different
 * components. Query every actual final triangle under the existing material
 * overlap-width rule. Source aliases remove artificial cone-cut boundaries
 * from the tolerance exemption; they never remove triangles from the query. */
static int ic_material_measure_visit(const MeshBinData *m, const double *uv,
                                      const int32_t *material_vertex, AsmContactStats *stats,
                                      AsmContactVisit visit, void *context)
{
    Arena_T arena = Arena_new();
    AsmContacts *query = NULL;
    int rc = -1, threads = 1;
    memset(stats,0,sizeof *stats);
    if (!arena || !m || !uv || !m->nv || !m->nf) goto done;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    query = AsmContacts_new_source(arena,m->nv,m->faces,m->nf,uv,
                                   AsmContacts_material_tolerance(),threads,material_vertex);
    if (!query) goto done;
    rc = AsmContacts_measure(query,uv,NULL,visit,context,stats);
    if (!stats->complete || stats->unmeasured || !isfinite(stats->area) || stats->area < 0) rc = -1;
    fprintf(stderr,"[ink_closeup] final material contact: policy=%s width=%.9g complete=%d "
            "pairs=%zu area=%.12g unmeasured=%zu\n",AsmContacts_material_policy(),
            AsmContacts_material_tolerance(),stats->complete,stats->pairs,stats->area,stats->unmeasured);
done:
    if (rc) stats->complete = 0;
    Arena_dispose(&arena);
    return rc;
}

static int ic_material_measure(const MeshBinData *m, const double *uv,
                                const int32_t *material_vertex, AsmContactStats *stats)
{
    return ic_material_measure_visit(m,uv,material_vertex,stats,NULL,NULL);
}

static int ic_material_qualified(const AsmContactStats *stats)
{
    return stats->complete && !stats->unmeasured && !stats->pairs && isfinite(stats->area) && stats->area >= 0;
}

static void ic_material_report(FILE *report, const AsmContactStats *stats)
{
    char area[64] = "null";
    if (isfinite(stats->area) && stats->area >= 0) snprintf(area,sizeof area,"%.12g",stats->area);
    fprintf(report," \"material_contact\": {\"policy\": \"%s\", \"tolerance\": %.12g, "
            "\"complete\": %s, \"pairs\": %zu, \"area\": %s, \"unmeasured\": %zu, "
            "\"qualified\": %s},\n",AsmContacts_material_policy(),AsmContacts_material_tolerance(),
            stats->complete?"true":"false",stats->pairs,area,stats->unmeasured,
            ic_material_qualified(stats)?"true":"false");
}

static int ic_material_read_aliases(const char *prefix,size_t nv,int32_t **out)
{
    char path[2800];*out=NULL;
    snprintf(path,sizeof path,"%s_surface_material_vertex.i32",prefix);
    FILE *fp=fopen(path,"rb");
    if(!fp)return errno==ENOENT?0:-1;
    if(nv>INT32_MAX||nv>SIZE_MAX/sizeof(int32_t)){fclose(fp);return -1;}
    int32_t *data=malloc(nv*sizeof *data);int rc=-1;
    if(data&&fread(data,sizeof *data,nv,fp)==nv&&fgetc(fp)==EOF&&!ferror(fp)) {
        rc=0;for(size_t v=0;v<nv;v++)if(data[v]<0||(size_t)data[v]>=nv){rc=-1;break;}
    }
    if(fclose(fp))rc=-1;
    if(rc)free(data);else *out=data;
    return rc;
}

static int ic_material_compose_aliases(size_t original_nv,const int32_t *original,
                                       size_t nv,const int32_t *source,int32_t **out)
{
    *out=NULL;if(original_nv>INT32_MAX||nv>INT32_MAX||(!source&&nv!=original_nv))return -1;
    int32_t *data=malloc(nv*sizeof *data);if(!data)return -1;
    for(size_t v=0;v<nv;v++) {
        int32_t p=source?source[v]:(int32_t)v;
        if(p<0||(size_t)p>=original_nv){free(data);return -1;}
        data[v]=original?original[p]:p;
        if(data[v]<0||(size_t)data[v]>=original_nv){free(data);return -1;}
    }
    *out=data;return 0;
}

static int ic_material_write_aliases(const char *target, size_t nv, const int32_t *material_vertex)
{
    char path[2800];
    int32_t *identity = NULL;
    if (!material_vertex) {
        if (nv > INT32_MAX) return -1;
        identity = malloc(nv*sizeof *identity);
        if (!identity) return -1;
        for (size_t v = 0; v < nv; v++) identity[v] = (int32_t)v;
        material_vertex = identity;
    }
    snprintf(path,sizeof path,"%s_surface_material_vertex.i32",target);
    int rc = ic_write_file(path,material_vertex,nv*sizeof *material_vertex);
    free(identity);
    return rc;
}

static int ic_material_selftest(void)
{
    /* Each square is individually injective. Only a complete all-component
     * audit sees their broad overlap. The established four-voxel band stays
     * allowed, while a six-voxel band must be refused. */
    const int32_t faces[12] = {0,1,2,0,2,3,4,5,6,4,6,7};
    double uv[16] = {0,0,20,0,20,20,0,20,14,0,34,0,34,20,14,20};
    MeshBinData mesh = {0}; AsmContactStats stats = {0};
    mesh.nv = 8; mesh.nf = 4; mesh.faces = (int32_t *)faces;
    int fails = ic_material_measure(&mesh,uv,NULL,&stats) || !stats.pairs || ic_material_qualified(&stats);
    for (int v = 4; v < 8; v++) uv[2*v] += 2;
    fails += ic_material_measure(&mesh,uv,NULL,&stats) || !ic_material_qualified(&stats);
    for (int v = 4; v < 8; v++) uv[2*v] += 8;
    fails += ic_material_measure(&mesh,uv,NULL,&stats) || !ic_material_qualified(&stats);
    int32_t bad_alias[8] = {0,1,2,3,4,5,6,8};
    fails += !ic_material_measure(&mesh,uv,bad_alias,&stats) || ic_material_qualified(&stats);
    {
        const int32_t original[4]={0,1,2,0},source[6]={0,1,2,3,1,3},expected[6]={0,1,2,0,1,0};
        int32_t *composed=NULL;
        fails+=ic_material_compose_aliases(4,original,6,source,&composed)!=0;
        if(composed){fails+=memcmp(composed,expected,sizeof expected)!=0;free(composed);}
        fails+=ic_material_compose_aliases(4,original,4,NULL,&composed)!=0;
        if(composed){fails+=memcmp(composed,original,sizeof original)!=0;free(composed);}
        const int32_t invalid[1]={4};
        fails+=ic_material_compose_aliases(4,original,1,invalid,&composed)==0||composed!=NULL;
    }
    fprintf(stderr,"[ink_closeup] all-component material contact: %s\n",fails?"FAIL":"ok");
    return fails != 0;
}

/* Read-only final-UV diagnostic for an immutable prior candidate. A literal
 * '-' explicitly requests identity topology; cut candidates need their alias
 * proof. This writes only the requested report, never the input mesh. */
static void ic_material_pair_csv(void *context, size_t a, size_t b, double area)
{
    fprintf((FILE *)context,"%zu,%zu,%.17g\n",a,b,area);
}

static int ic_material_audit(const char *mesh_path, const char *alias_path, const char *report_path,
                              const char *pairs_path)
{
    MeshBinData mesh = {0}; AsmContactStats stats = {0};
    int32_t *aliases = NULL; double *uv = NULL;
    FILE *report = NULL, *pairs = NULL; size_t got = 0; int rc = 1;
    if (MeshBin_read_precise_malloc(mesh_path,&mesh)) goto done;
    uv = malloc(2*mesh.nv*sizeof *uv);
    if (!uv) goto done;
    for (size_t k = 0; k < 2*mesh.nv; k++) uv[k] = MeshBin_uv(&mesh,k);
    if (strcmp(alias_path,"-")) {
        aliases = ic_read_file(alias_path,mesh.nv*sizeof *aliases,&got);
        if (!aliases || got != mesh.nv*sizeof *aliases) goto done;
    }
    if (pairs_path) {
        ves_ensure_parent_dir(pairs_path); pairs = fopen(pairs_path,"wb");
        if (!pairs) goto done;
        fprintf(pairs,"face_a,face_b,area\n");
    }
    (void)ic_material_measure_visit(&mesh,uv,aliases,&stats,pairs?ic_material_pair_csv:NULL,pairs);
    if (pairs) { int bad = ferror(pairs); bad |= fclose(pairs) != 0; pairs = NULL; if (bad) goto done; }
    ves_ensure_parent_dir(report_path); report = fopen(report_path,"wb");
    if (!report) goto done;
    fprintf(report,"{\n"); ic_material_report(report,&stats);
    fprintf(report," \"schema\": \"closeup-material-contact-audit-v1\", \"vertices\": %zu, \"faces\": %zu, "
            "\"source_aliases\": %s, \"source_modified\": false\n}\n",mesh.nv,mesh.nf,aliases?"true":"false");
    rc = fclose(report) != 0 || !ic_material_qualified(&stats); report = NULL;
done:
    if (report) fclose(report);
    if (pairs) fclose(pairs);
    free(aliases); free(uv); MeshBin_dispose(&mesh);
    return rc;
}
