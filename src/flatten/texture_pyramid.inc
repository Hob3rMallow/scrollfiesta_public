/* Native display pyramid. Level 1 copies every baked pixel exactly; coarser
 * levels are box averages for navigation only. No mesh or CT data changes. */
static void tp_reduce(const uint8_t *source,int w,int h,int step,int x0,int y0,int tw,int th,uint8_t *tile)
{
    for(int y=0;y<th;y++)for(int x=0;x<tw;x++){
        size_t sum=0,n=0;int sx=(x0+x)*step,sy=(y0+y)*step;
        int ex=(int)fmin(w,sx+step),ey=(int)fmin(h,sy+step);
        for(int b=sy;b<ey;b++)for(int a=sx;a<ex;a++){sum+=source[(size_t)b*w+a];n++;}
        tile[(size_t)y*tw+x]=(uint8_t)((sum+n/2)/n);
    }
}

int TexturedObj_pyramid(const char *input,const char *prefix)
{
    if(!input || !prefix || strlen(prefix)>3800)return -1;
    Arena_T arena=Arena_new();TiffRowReader reader={0};FILE *fp=NULL;int rc=-1;char path[4096];
    if(TiffIO_rows_open(arena,input,&reader))goto done;
    int w=reader.W,h=reader.H;if(w<=0 || h<=0 || (size_t)w>SIZE_MAX/(size_t)h)goto done;
    uint8_t *source=ARENA_ALLOC(arena,(size_t)w*h),*tile=ARENA_ALLOC(arena,(size_t)1024*1024);
    for(int y=0;y<h;y++)if(TiffIO_rows_read(&reader,y,0,w,source+(size_t)y*w))goto done;
    snprintf(path,sizeof path,"%s_pyramid.json",prefix);fp=fopen(path,"wb");if(!fp)goto done;
    fprintf(fp,"{\"schema\":\"native-texture-pyramid-v1\",\"width\":%d,\"height\":%d,\"tile_size\":1024,\"finest_level_pixel_exact\":true,\"levels\":[",w,h);
    int steps[]={1,4,16,64};
    for(size_t level=0;level<sizeof steps/sizeof *steps;level++){
        int step=steps[level],lw=(w-1)/step+1,lh=(h-1)/step+1;
        fprintf(fp,"%s{\"step\":%d,\"width\":%d,\"height\":%d}",level?",":"",step,lw,lh);
        for(int y=0;y<lh;y+=1024)for(int x=0;x<lw;x+=1024){
            int tw=(int)fmin(1024,lw-x),th=(int)fmin(1024,lh-y);
            tp_reduce(source,w,h,step,x,y,tw,th,tile);
            snprintf(path,sizeof path,"%s_s%d_y%d_x%d.png",prefix,step,y/1024,x/1024);
            if(VesPng_write_gray(path,tile,tw,th))goto done;
        }
    }
    fputs("]}\n",fp);if(ferror(fp))goto done;rc=0;
done:
    if(fp && fclose(fp))rc=-1;TiffIO_rows_close(&reader);Arena_dispose(&arena);
    fprintf(stderr,"[texture pyramid] %s: %s\n",prefix,rc?"FAILED":"complete; level 1 pixel-exact");return rc;
}

static int tp_selftest(void)
{
    uint8_t source[35],copy[35]={0},coarse[4]={0};for(int i=0;i<35;i++)source[i]=(uint8_t)i;
    tp_reduce(source,7,5,1,0,0,7,5,copy);int fail=memcmp(source,copy,sizeof source)!=0;
    tp_reduce(source,7,5,4,0,0,2,2,coarse);
    fail|=coarse[0]!=12 || coarse[1]!=16 || coarse[2]!=30 || coarse[3]!=33;
    uint8_t edge[6]={0};tp_reduce(source,7,5,1,5,2,2,3,edge);
    fail|=edge[0]!=19 || edge[1]!=20 || edge[4]!=33 || edge[5]!=34;
    fprintf(stderr,"[selftest] native exact texture tiles and partial box averages %s\n",fail?"FAIL":"PASS");return fail;
}
