    /* A contained segment must not tunnel through the moving rim even when
     * neither endpoint intersects the final edge. */
    {
        double uv_test[8]={0,0,0,2,.5,.5,.7,.5},old[2]={0,0};
        SRBoundaryEdge rim[2]={{0,1,0},{2,3,0}};
        size_t comp_offsets[2]={0,2},comp_edges[2]={0,1};
        size_t rim_offsets[5]={0,1,2,3,4},rim_edges[4]={0,0,1,1};
        SREdgeBox boxes_test[2];double prefix_test[2];
        SRSystem system={0};system.boundary=rim;
        system.component_boundary_offset=comp_offsets;system.component_boundary=comp_edges;
        system.rim_offset=rim_offsets;system.rim_edge=rim_edges;
        system.relax_boxes=boxes_test;system.relax_prefix=prefix_test;
        sr_relax_rim_index(&system,0,uv_test);
        uv_test[0]=2;
        int bad=sr_relax_rim_safe(&system,0,uv_test,0,old);
        uv_test[0]=-.2;
        bad|=!sr_relax_rim_safe(&system,0,uv_test,0,old);
        if (bad) { fprintf(stderr,"[slim selftest] swept-rim collision guard failed\n");failures++; }
    }
    /* A nearly closed slit exposes the endpoint-dependent collinearity
     * tolerance: the swept now->fixed side is clear, but the certificate's
     * stored fixed->now edge makes contact. Catch it during the local move so
     * one tiny slit cannot discard progress across the whole component. */
    {
        double uv_test[8]={0,0,1,.999,.9,.900000000002,.9,1.2};
        double old[2]={1,.999};
        SRBoundaryEdge rim[2]={{0,1,0},{2,3,0}};
        size_t comp_offsets[2]={0,2},comp_edges[2]={0,1};
        size_t rim_offsets[5]={0,1,2,3,4},rim_edges[4]={0,0,1,1};
        SREdgeBox boxes_test[2],audit_boxes[2];double prefix_test[2];
        SRSystem system={0};system.boundary=rim;
        system.component_boundary_offset=comp_offsets;system.component_boundary=comp_edges;
        system.rim_offset=rim_offsets;system.rim_edge=rim_edges;
        system.relax_boxes=boxes_test;system.relax_prefix=prefix_test;
        int bad=sr_component_boundary_intersections(&system,0,uv_test,audit_boxes,NULL)!=0;
        sr_relax_rim_index(&system,0,uv_test);
        uv_test[3]=1;
        bad|=sr_segments_intersect(uv_test+2,uv_test,uv_test+4,uv_test+6)!=0;
        bad|=sr_component_boundary_intersections(&system,0,uv_test,audit_boxes,NULL)!=1;
        bad|=sr_relax_rim_safe(&system,0,uv_test,1,old)!=0;
        uv_test[3]=.998;
        bad|=!sr_relax_rim_safe(&system,0,uv_test,1,old);
        if (bad) { fprintf(stderr,"[slim selftest] stored-order rim contact guard failed\n");failures++; }
    }
    /* All vertices of this stretched quad lie on the rim. Boundary relaxation
     * must reduce its true metric energy without creating a fold or crossing. */
    {
        float x_test[12]={0,0,0,1,0,0,1,1,0,0,1,0};
        int32_t f_test[6]={0,1,2,0,2,3};
        double work_test[8]={0,0,1.8,0,1.8,.7,0,.7};
        Arena_T scratch=Arena_new();SRSystem system={0};SlimRefineOpts options;
        SRMeasure before={0},after={0};SREdgeBox boxes_test[4];
        SlimRefine_defaults(&options);options.local_relaxation=4;options.local_relax_boundary=1;
        size_t moves=0;int bad=sr_build_system(scratch,x_test,4,f_test,2,&options,1,&system);
        if (!bad) {
            sr_relax_setup(scratch,&system);
            bad=sr_measure_component(&system,work_test,0,&before);
            for (int k=0;!bad && k<10;k++) {
                bad=sr_assemble(&system,work_test,&options);
                if (!bad) moves+=sr_component_relax(&system,&options,0,work_test);
            }
            bad|=sr_measure_component(&system,work_test,0,&after);
            bad|=sr_component_boundary_intersections(&system,0,work_test,boxes_test,NULL)!=0;
        }
        if (bad || !moves || after.flips || !(after.energy<.8*before.energy)) {
            fprintf(stderr,"[slim selftest] guarded-rim relaxation failed\n");failures++;
        }
        Arena_dispose(&scratch);
    }
    /* A slit-opening proposal must backtrack through the same face and rim
     * certificates as optimization. A reflex-corner-sized move can otherwise
     * invert a thin incident face before SLIM is allowed to start. */
    {
        float xyz_test[12]={0,0,0,1,0,0,1,1,0,0,1,0};
        int32_t faces_test[6]={0,1,2,0,2,3};
        double uv_test[8]={0,0,1,0,1,1,0,1};
        double displacement[8]={0,0,0,0,-2,0,0,0},opened[8]={0},alpha=0;
        Arena_T scratch=Arena_new();
        int bad=SlimRefine_open_seam_seed(scratch,xyz_test,4,faces_test,2,
                                         uv_test,displacement,opened,&alpha);
        bad|=fabs(alpha-.25)>1e-12 || fabs(opened[4]-.5)>1e-12 || !UvGuard_mesh(opened,faces_test,2);
        if (bad) { fprintf(stderr,"[slim selftest] certified slit opening failed\n");failures++; }
        Arena_dispose(&scratch);
    }
