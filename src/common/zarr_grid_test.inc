/* Native coordinate, codec, readback and refusal checks. Included privately. */
static uint8_t zgt_value(long z, long y, long x)
{
    return (uint8_t)((37*z + 11*y + 3*x) % 251);
}

static int zgt_metadata(const char *root, int compressed, const char *order)
{
    char path[ZG_PATH];
    FILE *f;
    int bad;
    snprintf(path,sizeof path,"%s/.zarray",root);
    f = fopen(path,"w");
    if (f == NULL) return -1;
    bad = fprintf(f,"{\"zarr_format\":2,\"shape\":[9,10,11],"
        "\"chunks\":[3,5,4],\"dtype\":\"|u1\",\"order\":\"%s\","
        "\"dimension_separator\":\"/\",\"fill_value\":0,\"filters\":null,"
        "\"compressor\":%s}",order,
        compressed?"{\"id\":\"blosc\",\"cname\":\"zstd\",\"clevel\":1,\"shuffle\":1}":"null") < 0;
    if (fclose(f)) bad = 1;
    return bad ? -1 : 0;
}

static int zgt_fixture(const char *root, int compressed)
{
    char path[ZG_PATH];
    uint8_t data[60], encoded[60+BLOSC2_MAX_OVERHEAD];
    if (ves_mkdir(root) || zgt_metadata(root,compressed,"C")) return -1;
    for (long z = 0; z < 3; z++) {
        snprintf(path,sizeof path,"%s/%ld",root,z);
        if (ves_mkdir(path)) return -1;
        for (long y = 0; y < 2; y++) {
            snprintf(path,sizeof path,"%s/%ld/%ld",root,z,y);
            if (ves_mkdir(path)) return -1;
            for (long x = 0; x < 3; x++) {
                const uint8_t *bytes = data;
                size_t n = sizeof data;
                FILE *f;
                int bad;
                for (long dz = 0; dz < 3; dz++)
                for (long dy = 0; dy < 5; dy++)
                for (long dx = 0; dx < 4; dx++)
                    data[(dz*5+dy)*4+dx] = zgt_value(z*3+dz,y*5+dy,x*4+dx);
                if (compressed) {
                    int nb = blosc1_compress(1,BLOSC_SHUFFLE,1,sizeof data,
                                            data,encoded,sizeof encoded);
                    if (nb <= 0) return -1;
                    bytes = encoded; n = (size_t)nb;
                }
                snprintf(path,sizeof path,"%s/%ld/%ld/%ld",root,z,y,x);
                f = fopen(path,"wb");
                if (f == NULL) return -1;
                bad = fwrite(bytes,1,n,f) != n;
                if (fclose(f)) bad = 1;
                if (bad) return -1;
            }
        }
    }
    return 0;
}

int ZarrGrid_selftest(void)
{
    char base[ZG_PATH], root[ZG_PATH], out[ZG_PATH], path[ZG_PATH];
    int fails = 0;
    snprintf(base,sizeof base,"zarr_grid_selftest_%d",ves_getpid());
    if (ves_mkdir(base)) return 1;
#define ZGT_CHECK(expr,what) do { if (!(expr)) { \
    fprintf(stderr,"[zarr grid selftest] FAIL %s\n",what); fails++; } } while(0)
    for (int compressed = 0; compressed <= 1; compressed++) {
        Arena_T arena = Arena_new();
        ZgArray a;
        uint8_t cube[64];
        const long origin[3] = {5,6,7}, zero[3] = {0,0,0};
        const long bad_origin[3] = {6,6,7};
        const long bbox[6] = {0,8,0,8,0,8};
        const long bad_bbox[6] = {0,12,0,8,0,8};
        snprintf(root,sizeof root,"%s/%s",base,compressed?"blosc":"raw");
        ZGT_CHECK(zgt_fixture(root,compressed) == 0,"fixture creation");
        if (zg_open(arena,root,&a)) {
            ZGT_CHECK(0,"valid metadata"); Arena_dispose(&arena); continue;
        }
        ZGT_CHECK(zg_cube(&a,origin,4,0,cube) == 0,"edge and chunk-straddling read");
        for (int z = 0; z < 4; z++)
        for (int y = 0; y < 4; y++)
        for (int x = 0; x < 4; x++)
            ZGT_CHECK(cube[(z*4+y)*4+x] == zgt_value(z+5,y+6,x+7),"ZYX byte identity");
        ZGT_CHECK(zg_cube(&a,zero,4,1,cube) == 0,"binary read");
        for (int z = 0; z < 4; z++)
        for (int y = 0; y < 4; y++)
        for (int x = 0; x < 4; x++)
            ZGT_CHECK(cube[(z*4+y)*4+x] == (zgt_value(z,y,x)?255:0),"binary preserves all nonzero voxels");
        ZGT_CHECK(zg_cube(&a,bad_origin,4,0,cube) != 0,"out-of-volume refused");
        snprintf(out,sizeof out,"%s/export",root);
        ZGT_CHECK(zg_export(&a,out,bad_bbox,4,0) != 0,"out-of-volume export refused");
        ZGT_CHECK(zg_export(&a,out,bbox,4,0) == 0,"TIFF export and complete readback");
        ZGT_CHECK(zg_export(&a,out,bbox,4,0) != 0,"existing output not overwritten");
        snprintf(path,sizeof path,"%s/z00004_y00004_x00004.tif",out);
        {
            uint8_t *v = NULL;
            int d = 0,h = 0,w = 0;
            int ok = TiffIO_load(arena,path,&v,&d,&h,&w) == 0 && d == 4 && h == 4 && w == 4;
            ZGT_CHECK(ok,"TIFF world-coordinate name and dimensions");
            if (ok) for (int z = 0; z < 4; z++)
            for (int y = 0; y < 4; y++)
            for (int x = 0; x < 4; x++)
                ZGT_CHECK(v[(z*4+y)*4+x] == zgt_value(z+4,y+4,x+4),"TIFF independent expected coordinates");
        }
        /* Break our own fixture, then open a fresh reader so cached chunks
         * cannot conceal a corrupt or missing file. */
        snprintf(path,sizeof path,"%s/0/0/0",root);
        {
            FILE *f = fopen(path,"wb");
            if (f) { fputc(0,f); fclose(f); }
            ZGT_CHECK(f != NULL,"corruption fixture opened");
        }
        ZGT_CHECK(zg_open(arena,root,&a) == 0,"fresh corrupt-chunk reader");
        ZGT_CHECK(zg_cube(&a,zero,4,0,cube) != 0,"truncated source refused");
        ZGT_CHECK(ves_unlink(path) == 0,"remove owned corruption fixture");
        ZGT_CHECK(zg_open(arena,root,&a) == 0,"fresh missing-chunk reader");
        ZGT_CHECK(zg_cube(&a,zero,4,0,cube) != 0,"missing source is not zero-filled");
        ZGT_CHECK(zgt_metadata(root,compressed,"F") == 0,"unsupported metadata fixture");
        ZGT_CHECK(zg_open(arena,root,&a) != 0,"Fortran-order input refused");
        Arena_dispose(&arena);
    }
#undef ZGT_CHECK
    fprintf(stderr,"[zarr grid selftest] %s (%d failures)\n",fails?"FAIL":"PASS",fails);
    return fails;
}
