/** Declaration of the output files */
FILE * fDrop;
FILE * fCA;
FILE * fPA;
FILE * falpha;
FILE * fUf;
FILE * fUd;

FILE * Dist_x;
char name_x[] = "Dist_x.csv";
char name_y[] = "Dist_y.csv";
FILE * Dist_y;
#if dimension == 3
FILE * Dist_z;
char name_z[] = "Dist_z.csv";
#endif

// FILE * Para;

/** definition of some usfull macros */

#if dimension == 2
#define POS {x,y} 
#else 
#define POS {x,y,z} 
#endif

/** Function used to initialise the drops on a square grid so the maximum 
 * possible volume fraction $\phi=\frac{\pi}{4}\approx0.785$*/

double geometry( double x, double y, double z){
  double inout=0.,circ=0.;
  srand(3458829);
  for(double i=(-Ls/2.+DP/2.); i<= Ls/2; i+=DP){
    for( double j=(-Ls/2.+DP/2.); j<= Ls/2.; j+=DP){
      #if dimension == 2
      int A=rand();
      int B=rand();
      circ= - sq(x-i+RND*(A*2./RAND_MAX-1)) 
            - sq(y-j+RND*(B*2./RAND_MAX-1)) 
            + sq(D/2);
      if(circ>0.) inout = 1.;
      #elif dimension == 3
      for( double k=(-Ls/2.+DP/2.); k<= Ls/2.; k+=DP){
        int A = rand();
        int B = rand(); 
        int C = rand();
        circ= - sq(x - i + RND*(A*2./RAND_MAX - 1)) 
              - sq(y - j + RND*(B*2./RAND_MAX - 1)) 
              - sq(z - k + RND*(C*2./RAND_MAX - 1)) 
              + sq(D/2);
        if(circ>0.) inout = 1.;
      }
      #endif
    }
  }
  return(inout);
}

/**Function that return true if the current process is of rank 0 */
int main_process(){
  #if _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  return rank == 0;
  #else
  return 1;
  #endif
}

/** Print header files */ 
void init_and_print_header(){
  if(main_process()){
    // files where the particules stats are stored 
    fDrop = fopen("fDrop.csv","w");
    fprintf(fDrop,"i,t,tagmax");
    fprintf(fDrop,",vol");
    fprintf(fDrop,",S");
    fprintf(fDrop,",ed");
    fprintf(fDrop,",Gnorm");
    fprintf(fDrop,",Pnorm");
    fprintf(fDrop,",p");
    fprintf(fDrop,",W2");
    #if dimension ==2
    fprintf(fDrop,",Y_x");
    fprintf(fDrop,",Y_y");
    fprintf(fDrop,",U_x");
    fprintf(fDrop,",U_y");
    fprintf(fDrop,",Uc_y");
    fprintf(fDrop,",Uc_x");
    fprintf(fDrop,",D_nbr_x");
    fprintf(fDrop,",D_nbr_y");
    fprintf(fDrop,",Fh_x");
    fprintf(fDrop,",Fh_y");
    // I J
    fprintf(fDrop,",G_x_x");
    fprintf(fDrop,",G_x_y");
    fprintf(fDrop,",G_y_x");
    fprintf(fDrop,",G_y_y");

    fprintf(fDrop,",P_x_x");
    fprintf(fDrop,",P_x_y");
    fprintf(fDrop,",P_y_x");
    fprintf(fDrop,",P_y_y");
    
    fprintf(fDrop,",WW_x_x");
    fprintf(fDrop,",WW_x_y");
    fprintf(fDrop,",WW_y_x");
    fprintf(fDrop,",WW_y_y");
    
    fprintf(fDrop,",Sig_x_x");
    fprintf(fDrop,",Sig_x_y");
    fprintf(fDrop,",Sig_y_x");
    fprintf(fDrop,",Sig_y_y");

    fprintf(fDrop,",Mh_x_x");
    fprintf(fDrop,",Mh_x_y");
    fprintf(fDrop,",Mh_y_x");
    fprintf(fDrop,",Mh_y_y");
    
    fprintf(fDrop,",Ms_x_x");
    fprintf(fDrop,",Ms_x_y");
    fprintf(fDrop,",Ms_y_x");
    fprintf(fDrop,",Ms_y_y");
    
    #elif dimension ==3
    // X
    fprintf(fDrop,",Y_x");
    fprintf(fDrop,",Y_y");
    fprintf(fDrop,",Y_z");
    fprintf(fDrop,",U_x");
    fprintf(fDrop,",U_y");
    fprintf(fDrop,",U_z");
    fprintf(fDrop,",Uc_x");
    fprintf(fDrop,",Uc_y");
    fprintf(fDrop,",Uc_z");
    fprintf(fDrop,",D_nbr_x");
    fprintf(fDrop,",D_nbr_y");
    fprintf(fDrop,",D_nbr_z");
    fprintf(fDrop,",Fh_x");
    fprintf(fDrop,",Fh_y");
    fprintf(fDrop,",Fh_z");

    // I J
    fprintf(fDrop,",G_x_x");
    fprintf(fDrop,",G_x_y");
    fprintf(fDrop,",G_x_z");
    fprintf(fDrop,",G_y_x");
    fprintf(fDrop,",G_y_y");
    fprintf(fDrop,",G_y_z");
    fprintf(fDrop,",G_z_x");
    fprintf(fDrop,",G_z_y");
    fprintf(fDrop,",G_z_z");

    fprintf(fDrop,",P_x_x");
    fprintf(fDrop,",P_x_y");
    fprintf(fDrop,",P_x_z");
    fprintf(fDrop,",P_y_x");
    fprintf(fDrop,",P_y_y");
    fprintf(fDrop,",P_y_z");
    fprintf(fDrop,",P_z_x");
    fprintf(fDrop,",P_z_y");
    fprintf(fDrop,",P_z_z");

    fprintf(fDrop,",WW_x_x");
    fprintf(fDrop,",WW_x_y");
    fprintf(fDrop,",WW_x_z");
    fprintf(fDrop,",WW_y_x");
    fprintf(fDrop,",WW_y_y");
    fprintf(fDrop,",WW_y_z");
    fprintf(fDrop,",WW_z_x");
    fprintf(fDrop,",WW_z_y");
    fprintf(fDrop,",WW_z_z");

    fprintf(fDrop,",Sig_x_x");
    fprintf(fDrop,",Sig_x_y");
    fprintf(fDrop,",Sig_x_z");
    fprintf(fDrop,",Sig_y_x");
    fprintf(fDrop,",Sig_y_y");
    fprintf(fDrop,",Sig_y_z");
    fprintf(fDrop,",Sig_z_x");
    fprintf(fDrop,",Sig_z_y");
    fprintf(fDrop,",Sig_z_z");

    fprintf(fDrop,",Mh_x_x");
    fprintf(fDrop,",Mh_x_y");
    fprintf(fDrop,",Mh_x_z");
    fprintf(fDrop,",Mh_y_x");
    fprintf(fDrop,",Mh_y_y");
    fprintf(fDrop,",Mh_y_z");
    fprintf(fDrop,",Mh_z_x");
    fprintf(fDrop,",Mh_z_y");
    fprintf(fDrop,",Mh_z_z");

    fprintf(fDrop,",Ms_x_x");
    fprintf(fDrop,",Ms_x_y");
    fprintf(fDrop,",Ms_x_z");
    fprintf(fDrop,",Ms_y_x");
    fprintf(fDrop,",Ms_y_y");
    fprintf(fDrop,",Ms_y_z");
    fprintf(fDrop,",Ms_z_x");
    fprintf(fDrop,",Ms_z_y");
    fprintf(fDrop,",Ms_z_z");

    #endif
    fprintf(fDrop,",ct");
    fprintf(fDrop,",tagmin");
    fprintf(fDrop,",tag");
    fprintf(fDrop,"\n");
    fclose(fDrop);
    
    // files where the fluid and solid phase stats are stored 
    fCA = fopen("fCA.csv","w");
    fprintf(fCA,"i,t,NVOF,PHI,S");
    fprintf(fCA,",Vd,Vf,dissd,dissf");
    #if dimension == 2
    fprintf(fCA,",U_x");
    fprintf(fCA,",U_y");
    fprintf(fCA,",Ud_x");
    fprintf(fCA,",Ud_y");
    fprintf(fCA,",Uf_x");
    fprintf(fCA,",Uf_y");
    fprintf(fCA,",RhoU_x");
    fprintf(fCA,",RhoU_y");
    fprintf(fCA,",A_x");
    fprintf(fCA,",A_y");

    fprintf(fCA,",UpUpf_x_x");
    fprintf(fCA,",UpUpf_x_y");
    fprintf(fCA,",UpUpf_y_x");
    fprintf(fCA,",UpUpf_y_y");
    
    fprintf(fCA,",UpUpd_x_x");
    fprintf(fCA,",UpUpd_x_y");
    fprintf(fCA,",UpUpd_y_x");
    fprintf(fCA,",UpUpd_y_y");

    fprintf(fCA,",UpUpUpf_x");
    fprintf(fCA,",UpUpUpf_x");
    
    fprintf(fCA,",UpUpUpd_x");
    fprintf(fCA,",UpUpUpd_y");
    #elif dimension == 3
    fprintf(fCA,",U_x");
    fprintf(fCA,",U_y");
    fprintf(fCA,",U_z");
    fprintf(fCA,",Ud_x");
    fprintf(fCA,",Ud_y");
    fprintf(fCA,",Ud_z");
    fprintf(fCA,",Uf_x");
    fprintf(fCA,",Uf_y");
    fprintf(fCA,",Uf_z");
    fprintf(fCA,",RhoU_x");
    fprintf(fCA,",RhoU_y");
    fprintf(fCA,",RhoU_z");
    fprintf(fCA,",A_x");
    fprintf(fCA,",A_y");
    fprintf(fCA,",A_z");

    fprintf(fCA,",UpUpf_x_x");
    fprintf(fCA,",UpUpf_x_y");
    fprintf(fCA,",UpUpf_x_z");
    fprintf(fCA,",UpUpf_y_x");
    fprintf(fCA,",UpUpf_y_y");
    fprintf(fCA,",UpUpf_y_z");
    fprintf(fCA,",UpUpf_z_x");
    fprintf(fCA,",UpUpf_z_y");
    fprintf(fCA,",UpUpf_z_z");

    fprintf(fCA,",UpUpd_x_x");
    fprintf(fCA,",UpUpd_x_y");
    fprintf(fCA,",UpUpd_x_z");
    fprintf(fCA,",UpUpd_y_x");
    fprintf(fCA,",UpUpd_y_y");
    fprintf(fCA,",UpUpd_y_z");
    fprintf(fCA,",UpUpd_z_x");
    fprintf(fCA,",UpUpd_z_y");
    fprintf(fCA,",UpUpd_z_z");

    fprintf(fCA,",UpUpUpf_x");
    fprintf(fCA,",UpUpUpf_y");
    fprintf(fCA,",UpUpUpf_z");

    fprintf(fCA,",UpUpUpd_x");
    fprintf(fCA,",UpUpUpd_y");
    fprintf(fCA,",UpUpUpd_z");

    #endif
    fprintf(fCA,"\n");
    fclose(fCA);
  

    // files where the particular averaged quantities are stored 
    fPA = fopen("fPA.csv","w");
    fprintf(fPA,"i,t,V,S,ed");
    fprintf(fPA,",Gnorm");
    fprintf(fPA,",Pnorm");
    fprintf(fPA,",W2");
    #if dimension == 2
    fprintf(fPA,",U_x");
    fprintf(fPA,",U_y");
    fprintf(fPA,",Uc_x");
    fprintf(fPA,",Uc_y");

    fprintf(fPA,",G_x_x");
    fprintf(fPA,",G_x_y");
    fprintf(fPA,",G_y_x");
    fprintf(fPA,",G_y_y");

    fprintf(fPA,",P_x_x");
    fprintf(fPA,",P_x_y");
    fprintf(fPA,",P_y_x");
    fprintf(fPA,",P_y_y");

    fprintf(fPA,",VpUp_x");
    fprintf(fPA,",VpUp_y");
    
    fprintf(fPA,",VpUpUp_x_x");
    fprintf(fPA,",VpUpUp_x_y");
    fprintf(fPA,",VpUpUp_y_x");
    fprintf(fPA,",VpUpUp_y_y");    

    fprintf(fPA,",UpUp_x_x");
    fprintf(fPA,",UpUp_x_y");
    fprintf(fPA,",UpUp_y_x");
    fprintf(fPA,",UpUp_y_y");

    fprintf(fPA,",UpUpUp_x");
    fprintf(fPA,",UpUpUp_y");
    
    fprintf(fPA,",Sig_x_x");
    fprintf(fPA,",Sig_x_y");
    fprintf(fPA,",Sig_y_x");
    fprintf(fPA,",Sig_y_y");
    
    fprintf(fPA,",Mh_x_x");
    fprintf(fPA,",Mh_x_y");
    fprintf(fPA,",Mh_y_x");
    fprintf(fPA,",Mh_y_y");

    fprintf(fPA,",Ms_x_x");
    fprintf(fPA,",Ms_x_y");
    fprintf(fPA,",Ms_y_x");
    fprintf(fPA,",Ms_y_y");

    fprintf(fPA,",WW_x_x");
    fprintf(fPA,",WW_x_y");
    fprintf(fPA,",WW_y_x");
    fprintf(fPA,",WW_y_y");
    fprintf(fPA,",PFP_x_x");
    fprintf(fPA,",PFP_x_y");
    fprintf(fPA,",PFP_y_x");
    fprintf(fPA,",PFP_y_y");
    #elif dimension == 3
    fprintf(fPA,",U_x");
    fprintf(fPA,",U_y");
    fprintf(fPA,",U_z");
    fprintf(fPA,",Uc_x");
    fprintf(fPA,",Uc_y");
    fprintf(fPA,",Uc_z");

    fprintf(fPA,",G_x_x");
    fprintf(fPA,",G_x_y");
    fprintf(fPA,",G_x_z");
    fprintf(fPA,",G_y_x");
    fprintf(fPA,",G_y_y");
    fprintf(fPA,",G_y_z");
    fprintf(fPA,",G_z_x");
    fprintf(fPA,",G_z_y");
    fprintf(fPA,",G_z_z");

    fprintf(fPA,",P_x_x");
    fprintf(fPA,",P_x_y");
    fprintf(fPA,",P_x_z");
    fprintf(fPA,",P_y_x");
    fprintf(fPA,",P_y_y");
    fprintf(fPA,",P_y_z");
    fprintf(fPA,",P_z_x");
    fprintf(fPA,",P_z_y");
    fprintf(fPA,",P_z_z");

    fprintf(fPA,",VpUp_x");
    fprintf(fPA,",VpUp_y");
    fprintf(fPA,",VpUp_z");

    fprintf(fPA,",VpUpUp_x_x");
    fprintf(fPA,",VpUpUp_x_y");
    fprintf(fPA,",VpUpUp_x_z");
    fprintf(fPA,",VpUpUp_y_x");
    fprintf(fPA,",VpUpUp_y_y");    
    fprintf(fPA,",VpUpUp_y_z");    
    fprintf(fPA,",VpUpUp_z_x");
    fprintf(fPA,",VpUpUp_z_y");    
    fprintf(fPA,",VpUpUp_z_z");    

    fprintf(fPA,",UpUp_x_x");
    fprintf(fPA,",UpUp_x_y");
    fprintf(fPA,",UpUp_x_z");
    fprintf(fPA,",UpUp_y_x");
    fprintf(fPA,",UpUp_y_y");
    fprintf(fPA,",UpUp_y_z");
    fprintf(fPA,",UpUp_z_x");
    fprintf(fPA,",UpUp_z_y");
    fprintf(fPA,",UpUp_z_z");

    fprintf(fPA,",UpUpUp_x");
    fprintf(fPA,",UpUpUp_y");
    fprintf(fPA,",UpUpUp_z");
        
    fprintf(fPA,",Sig_x_x");
    fprintf(fPA,",Sig_x_y");
    fprintf(fPA,",Sig_x_z");
    fprintf(fPA,",Sig_y_x");
    fprintf(fPA,",Sig_y_y");
    fprintf(fPA,",Sig_y_z");
    fprintf(fPA,",Sig_z_x");
    fprintf(fPA,",Sig_z_y");
    fprintf(fPA,",Sig_z_z");
    
    fprintf(fPA,",Mh_x_x");
    fprintf(fPA,",Mh_x_y");
    fprintf(fPA,",Mh_x_z");
    fprintf(fPA,",Mh_y_x");
    fprintf(fPA,",Mh_y_y");
    fprintf(fPA,",Mh_y_z");
    fprintf(fPA,",Mh_z_x");
    fprintf(fPA,",Mh_z_y");
    fprintf(fPA,",Mh_z_z");

    fprintf(fPA,",Ms_x_x");
    fprintf(fPA,",Ms_x_y");
    fprintf(fPA,",Ms_x_z");
    fprintf(fPA,",Ms_y_x");
    fprintf(fPA,",Ms_y_y");
    fprintf(fPA,",Ms_y_z");
    fprintf(fPA,",Ms_z_x");
    fprintf(fPA,",Ms_z_y");
    fprintf(fPA,",Ms_z_z");

    fprintf(fPA,",WW_x_x");
    fprintf(fPA,",WW_x_y");
    fprintf(fPA,",WW_x_z");
    fprintf(fPA,",WW_y_x");
    fprintf(fPA,",WW_y_y");
    fprintf(fPA,",WW_y_z");
    fprintf(fPA,",WW_z_x");
    fprintf(fPA,",WW_z_y");
    fprintf(fPA,",WW_z_z");

    fprintf(fPA,",PFP_x_x");
    fprintf(fPA,",PFP_x_y");
    fprintf(fPA,",PFP_x_z");
    fprintf(fPA,",PFP_y_x");
    fprintf(fPA,",PFP_y_y");
    fprintf(fPA,",PFP_y_z");
    fprintf(fPA,",PFP_z_x");
    fprintf(fPA,",PFP_z_y");
    fprintf(fPA,",PFP_z_z");
    #endif
    fprintf(fPA,"\n");
    fclose(fPA);


    // Continuity wave file header 
    falpha = fopen("falpha.csv","w");
    int nL = round(Ls/(D/LDn));
    fprintf(falpha,"t");
    for (int i = 0; i < nL; i++)
      fprintf(falpha,",L%d",i);
    fprintf(falpha,"\n");
    fclose(falpha);

    fUf = fopen("fUf.csv","w");
    fprintf(fUf,"t");
    for (int i = 0; i < nL; i++)
      fprintf(fUf,",L%d",i);
    fprintf(fUf,"\n");
    fclose(fUf);

    fUd = fopen("fUd.csv","w");
    fprintf(fUd,"t");
    for (int i = 0; i < nL; i++)
      fprintf(fUd,",L%d",i);
    fprintf(fUd,"\n");
    fclose(fUd);

    foreach_dimension(){
        Dist_x = fopen(name_x,"w");
        fprintf(Dist_x,"i,t");
        for (int j = 0; j < Nb; j++)
          for (int k=j+1; k < Nb; k++)
            fprintf(Dist_x,",b%db%d",j,k);
      fprintf(Dist_x,"\n");
      fclose(Dist_x);
    }
  }
}

void copy_file(char * sname,char * tname, FILE * target){
  char ch;
  FILE * Source = fopen(sname,"r");
  target = fopen(tname,"w");
  while( ( ch = fgetc(Source) ) != EOF )
    fputc(ch, target); 
  fclose(Source);
  fclose(target);
}

void restore_old_csv(){
    copy_file("../dump/fDrop.csv","fDrop.csv",fDrop);
    copy_file("../dump/fCA.csv","fCA.csv",fCA);
    copy_file("../dump/fPA.csv","fPA.csv",fPA);
    copy_file("../dump/falpha.csv","falpha.csv",falpha);
    copy_file("../dump/fUf.csv","fUf.csv",fUf);
    copy_file("../dump/fUd.csv","fUd.csv",fUd);
    copy_file("../dump/Dist_x.csv","Dist_x.csv",Dist_x);
    copy_file("../dump/Dist_y.csv","Dist_y.csv",Dist_y);
    #if dimension == 3
    copy_file("../dump/Dist_z.csv","Dist_z.csv",Dist_z);
    #endif
    system("cp ../dump/*mp4*  .");
    fprintf(stdout,"RESTORED!");
}




/** compute the averaged density */

double avg_rho(){
  double rhv=0.,vol=0;
  foreach(reduction(+:rhv)reduction(+:vol)){
    vol += dv();
    rhv += dv()*(f[]*(rho1-rho2)+rho2);
  }
  return rhv/vol;
}
/** compute the bulk velocity  */
coord avg_U(){
  coord U={0};
  double vol=0;
  foreach(reduction(+:U)reduction(+:vol))
    foreach_dimension(){
      U.x += u.x[]*dv();
      vol += dv();
    }
  return div_coord(U,vol);
}

/** Macro that define the basic operator 
*/

#if dimension == 2
#define vec_grad(u,G)\
  {\
    G.x.x =  ( u.x[1,0] - u.x[-1,0] )/(2.*Delta);\
    G.x.y =  ( u.x[0,1] - u.x[0,-1] )/(2.*Delta);\
    G.y.x =  ( u.y[1,0] - u.y[-1,0] )/(2.*Delta);\
    G.y.y =  ( u.y[0,1] - u.y[0,-1] )/(2.*Delta);\
  }

#define vec_lap(u,L)\
  {\
      L.x = (u.x[1,0] + u.x[-1,0] - 2*u.x[]                  \
              + u.x[0,1] + u.x[0,-1] - 2*u.x[] )/(Delta*Delta);\
      L.y = (u.y[1,0] + u.y[-1,0] - 2*u.y[]                  \
              + u.y[0,1] + u.y[0,-1] - 2*u.y[] )/(Delta*Delta);\
  }

#define scalar_grad(p,G)\
  {\
      G.x = (p[1,0] - p[-1,0] )/(2.*Delta);\
      G.y = (p[0,1] - p[0,-1] )/(2.*Delta);\
  }

#elif dimension == 3
#define vec_grad(u,G)\
  {\
    G.x.x = (u.x[1]     - u.x[-1]    )/(2.*Delta);\
    G.x.y = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);\
    G.x.z = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);\
    G.y.x = (u.y[1]     - u.y[-1]    )/(2.*Delta);\
    G.y.y = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);\
    G.y.z = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);\
    G.z.x = (u.z[1]     - u.z[-1]    )/(2.*Delta);\
    G.z.y = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);\
    G.z.z = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);\
  }

#define vec_lap(u,L)\
  {\
      L.x = (u.x[1] + u.x[-1] - 2*u.x[] \
              + u.x[0,1] + u.x[0,-1] - 2*u.x[] \
              + u.x[0,0,1] + u.x[0,0,-1] - 2*u.x[] )/(Delta*Delta);\
      L.y = (u.y[1] + u.y[-1] - 2*u.y[] \
              + u.y[0,1] + u.y[0,-1] - 2*u.y[] \
              + u.y[0,0,1] + u.y[0,0,-1] - 2*u.y[] )/(Delta*Delta);\
      L.z = (u.z[1] + u.z[-1] - 2*u.z[] \
              + u.z[0,1] + u.z[0,-1] - 2*u.z[] \
              + u.z[0,0,1] + u.z[0,0,-1] - 2*u.z[] )/(Delta*Delta);\
  }

#define scalar_grad(p,G)\
  {\
      G.x = (p[1,0,0] - p[-1,0,0] )/(2.*Delta);\
      G.y = (p[0,1,0] - p[0,-1,0] )/(2.*Delta);\
      G.z = (p[0,0,1] - p[0,0,-1] )/(2.*Delta);\
  }
#endif