    for (int quarter=-1;quarter<=1;quarter++) {
        float vertices[27]={0};int32_t faces[24];SlimRefineDirection directions[8];
        uint8_t terminals[9]={0};size_t count=0;
        for (int k=0;k<8;k++) {
            double t=2*M_PI*(double)k/8,phase=(double)quarter*2*M_PI*((double)k+.5)/32;
            vertices[3*(k+1)]=(float)cos(t);vertices[3*(k+1)+1]=(float)sin(t);
            faces[3*k]=0;faces[3*k+1]=k+1;faces[3*k+2]=(k+1)%8+1;
            SlimRefineDirection d={{cos(phase),sin(phase),0},{-sin(phase),cos(phase),0},1};directions[k]=d;
        }
        int bad=SeamCut_direction_terminals(arena,vertices,9,faces,8,directions,terminals,&count);
        if (count!=(quarter?1u:0u) || terminals[0]!=(quarter?1:0)) bad=1;
        if (quarter && !bad) {
            SeamCutResult cut={0};bad=SeamCut_to_boundary(arena,vertices,9,faces,8,terminals,&cut);
            if (!bad) {
                Conn before,after;conn_build(arena,vertices,9,faces,8,&before);
                conn_build(arena,cut.verts,cut.nv,cut.faces,cut.nf,&after);
                int center_on_boundary=0;
                for (size_t c=0;c<after.nc;c++) if (after.twin[c]<0 &&
                    (cut.vmap[after.faces[c]]==0 || cut.vmap[after.faces[NEXT(c)]]==0)) center_on_boundary=1;
                if (cut.nf!=8 || !cut.seam_edges || !center_on_boundary || fabs(distortion_r0(&before,0))>1e-10) bad=1;
                for (size_t c=0;c<3*cut.nf;c++) if (cut.vmap[cut.faces[c]]!=faces[c]) bad=1;
                for (size_t v=0;v<cut.nv;v++) if (memcmp(cut.verts+3*v,vertices+3*(size_t)cut.vmap[v],3*sizeof(float))) bad=1;
            }
            /* One absent face breaks the evidence loop; it cannot create a
             * terminal, even though all other faces still measure a cross. */
            directions[3].weight=0;memset(terminals,0,sizeof terminals);
            bad|=SeamCut_direction_terminals(arena,vertices,9,faces,8,directions,terminals,&count);
            if (count || terminals[0]) bad=1;
        }
        fprintf(stderr,"  [seamcut fiber quarter %d] low-curvature cut/XYZ/faces/support -> %s\n",quarter,bad?"FAIL":"ok");
        fails+=bad!=0;
    }
    {
        /* Two closed fans sharing only their center have no single interior
         * loop. Summing their angle deficits would invent a full-turn index. */
        float vertices[51]={0};int32_t faces[48];SlimRefineDirection directions[16];uint8_t terminals[17]={0};size_t count=0;
        for (int fan=0;fan<2;fan++) for (int k=0;k<8;k++) {
            int v=1+8*fan+k,f=8*fan+k;double t=2*M_PI*k/8;
            vertices[3*v]=(float)cos(t);vertices[3*v+(fan?2:1)]=(float)sin(t);
            faces[3*f]=0;faces[3*f+1]=v;faces[3*f+2]=1+8*fan+(k+1)%8;
            SlimRefineDirection d={{1,0,0},{0,fan?0:1,fan?1:0},1};directions[f]=d;
        }
        int bad=SeamCut_direction_terminals(arena,vertices,17,faces,16,directions,terminals,&count);
        if (bad || count || terminals[0]) { fprintf(stderr,"  [seamcut fiber] pinched fans invented an index\n");fails++; }
    }
