/* Private angular model controls. No scratch artifacts or production API. */
int PrivateJointAngle_selftest(void)
{
    int fails=0;double maximum_fd=0,maximum_proxy=0,maximum_symmetry=0;
    const int32_t face[3]={0,1,2};
    for(int test=0;test<160;test++) {
        double t=.037*(double)test,c=cos(t),s=sin(t);
        SRDirection d={{-c-s,c,s},{s-c,-s,c},1.7};
        double a=.7+.4*(1+sin(.17*test)),b=.8+.3*(1+cos(.21*test));
        double turn=.091*test,ct=cos(turn),st=sin(turn),shear=.18*sin(.13*test);
        double uv[6]={.23,-.17,.23+a*ct,-.17+a*st,.23+shear*ct-b*st,-.17+shear*st+b*ct};
        double r[2],g[2][6];if(sra_model(&d,face,uv,r,g)){fails++;continue;}
        double analytic[6]={0},matrix[36]={0};
        for(int arm=0;arm<2;arm++)for(int i=0;i<6;i++) {
            analytic[i]+=2*d.weight*r[arm]*g[arm][i];
            for(int j=0;j<6;j++)matrix[6*i+j]+=2*d.weight*g[arm][i]*g[arm][j];
        }
        for(int i=0;i<6;i++) {
            double step=1e-6,probe[6];memcpy(probe,uv,sizeof probe);probe[i]+=step;
            double plus=sra_energy(&d,face,probe);probe[i]-=2*step;
            double minus=sra_energy(&d,face,probe),fd=(plus-minus)/(2*step);
            maximum_fd=fmax(maximum_fd,fabs(fd-analytic[i]));
            if(fabs(fd-analytic[i])>2e-7*fmax(1,fabs(analytic[i])))fails++;
        }
        double radial=0,psd=0,delta[6];
        for(int i=0;i<6;i++){radial+=analytic[i]*uv[i];delta[i]=sin(.37*(test+1)*(i+1));}
        for(int i=0;i<6;i++)for(int j=0;j<6;j++)psd+=delta[i]*matrix[6*i+j]*delta[j];
        if(fabs(radial)>1e-12 || psd< -1e-12)fails++;
        double transformed[6],energy=sra_energy(&d,face,uv);
        for(int i=0;i<6;i++)transformed[i]=2.7*uv[i];
        if(fabs(sra_energy(&d,face,transformed)-energy)>1e-12)fails++;
        for(int k=0;k<3;k++){transformed[2*k]=-uv[2*k+1];transformed[2*k+1]=uv[2*k];}
        if(fabs(sra_energy(&d,face,transformed)-energy)>1e-12)fails++;
        for(int mode=0;mode<4;mode++) {
            SRDirection alternate=d;
            for(int k=0;k<3;k++) {
                alternate.h[k]=mode&1?d.v[k]:-d.h[k];
                alternate.v[k]=mode&1?d.h[k]:d.v[k];
                if(mode&2)alternate.v[k]=-alternate.v[k];
            }
            double ar[2],ag[2][6];if(sra_model(&alternate,face,uv,ar,ag)){fails++;continue;}
            if(fabs(sra_energy(&alternate,face,uv)-energy)>1e-12)fails++;
            for(int i=0;i<6;i++) {
                double value=2*d.weight*(ar[0]*ag[0][i]+ar[1]*ag[1][i]);
                maximum_symmetry=fmax(maximum_symmetry,fabs(value-analytic[i]));
                if(fabs(value-analytic[i])>1e-11)fails++;
            }
        }
        /* Exercise the actual sparse assembly, including off-diagonal XY
         * blocks and RHS, against the exact local residual gradient. */
        SRSystem system={0};int32_t edges[3]={0,1,2};
        double diagonal[9]={0},blocks[9]={0},rhs[6]={0};
        system.faces=face;system.direction=&d;system.face_edge=edges;
        system.diagonal=diagonal;system.edge_block=blocks;system.rhs=rhs;
        if(sra_assemble_face(&system,0,uv)){fails++;continue;}
        double product[6]={0};
        for(int k=0;k<3;k++) {
            product[2*k]=diagonal[3*k]*uv[2*k]+diagonal[3*k+1]*uv[2*k+1]-rhs[2*k];
            product[2*k+1]=diagonal[3*k+1]*uv[2*k]+diagonal[3*k+2]*uv[2*k+1]-rhs[2*k+1];
        }
        for(int k=0;k<3;k++) {
            int q=(k+1)%3;double *bb=blocks+3*k;
            product[2*k]+=bb[0]*uv[2*q]+bb[1]*uv[2*q+1];
            product[2*k+1]+=bb[1]*uv[2*q]+bb[2]*uv[2*q+1];
            product[2*q]+=bb[0]*uv[2*k]+bb[1]*uv[2*k+1];
            product[2*q+1]+=bb[1]*uv[2*k]+bb[2]*uv[2*k+1];
        }
        for(int i=0;i<6;i++) {
            maximum_proxy=fmax(maximum_proxy,fabs(2*product[i]-analytic[i]));
            if(fabs(2*product[i]-analytic[i])>1e-11)fails++;
        }
    }
    {
        /* Exact equal-cost quarter boundary: energy uses the true minimum;
         * its right derivative is the smaller branch derivative, not a
         * fictitious smooth derivative across the branch switch. */
        SRDirection d={{-1,1,0},{-1,0,1},1};
        double q=sqrt(.5),uv[6]={0,0,q,q,-q,q},probe[6],r[2],g[2][6];
        if(sra_model(&d,face,uv,r,g))fails++;
        double e=sra_energy(&d,face,uv),step=1e-7;
        for(int sign=-1;sign<=1;sign+=2) {
            double angle=sign*step,ct=cos(angle),st=sin(angle);
            for(int k=0;k<3;k++){probe[2*k]=ct*uv[2*k]-st*uv[2*k+1];probe[2*k+1]=st*uv[2*k]+ct*uv[2*k+1];}
            if(fabs((sra_energy(&d,face,probe)-e)/step+3.14159265358979323846)>1e-5)fails++;
        }
        /* Missing support produces exactly zero without inspecting bad arms. */
        d.weight=0;for(int i=0;i<6;i++)probe[i]=NAN;
        if(sra_energy(&d,face,probe)!=0 || sra_model(&d,face,probe,r,g))fails++;
        for(int a=0;a<2;a++){if(r[a]!=0)fails++;for(int i=0;i<6;i++)if(g[a][i]!=0)fails++;}
    }
    fprintf(stderr,"[angular GN controls] gradient FD max=%.17g sparse proxy max=%.17g arm-equivariance max=%.17g failures=%d\n",maximum_fd,maximum_proxy,maximum_symmetry,fails);
    return fails;
}
