/** Declaration of the output files */
FILE * fDrop;
FILE * fCA;
FILE * fPA;
#if BI_DISPERSE
FILE * fPA1;
FILE * fPA2;
#endif
FILE * falpha;
FILE * fCAnst;
// FILE * Para;

/** definition of some usfull macros */

#if dimension == 2
#define POS {x,y} 
#else 
#define POS {x,y,z} 
#endif
      // if (s[] < 1. - EPS){
        // coord  normal = mycs(point, s), pvol;
        // double alpha = plane_alpha(s[], normal);
        // // Compute volume of fragment with centroid pcell
        // plane_center (normal, alpha, s[], &pvol);
        // // The centroid of the interface volume is corrected (not equal to cell centroid)
        // pos = add_coord(pos,mult_coord(pvol,Delta));
      // }
/** Function used to initialise the drops on a square grid so the maximum 
 * possible volume fraction $\phi=\frac{\pi}{4}\approx0.785$
 * This funciton initialize an ordered array of spheres or circle with a
 * small random shift for each sphere. 
 * */

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

/**
 * This funciton initialize a random array of non touching spheres
 * so that the simulation reach a random state faster. 
 */

int check_all_dist(const coord * Position,const  double * Radii, int j){
  for (int k = 0; k < j; k++)
    if(dist_perio(Position[j],Position[k]) < ((Radii[j]+Radii[k]) * 1.2))
      return 0; 
  return 1;
}

void generate_drops_pos(coord * Position, double * Radii){
  srand(time(NULL));
  // srand(1716990199);
  int j =0;
  while (j < Nb)
  {
    double Rm = (double) RAND_MAX + 1.0;
    Position[j].x = ((double) rand()/Rm) * Ls - Ls/2;
    Position[j].y = ((double) rand()/Rm) * Ls - Ls/2;
    Position[j].z = ((double) rand()/Rm) * Ls - Ls/2;
    if(check_all_dist(Position, Radii,j))
      j++;
  }
}

int geometry_rand(double x, double y, double z, const coord * Position, const double * Radii){
  for (int j = 0; j < Nb; j++)
  {
    // compute the periodic dist to the center of the particules
    coord Pos = POS;
    coord r = dist_perio_coord(Pos,Position[j]);
    #if dimension == 2
    double circ = - sq(r.x) - sq(r.y) + sq(Radii[j]);
    #else
    double circ = - sq(r.x) - sq(r.y) - sq(r.z) + sq(Radii[j]);
    #endif
    if(circ > 0.)
      return 1;
  }
  return 0;
}

double refine_frac(coord center, coord * dpos, double * dradius, double step, int level,int maxlevel){
  int Rpd=0, Rmd=0;
  for (int j = 0; j < Nb; j++){
      coord r = dist_perio_coord(center,dpos[j]);
      Rmd =  (sq(r.x) + sq(r.y) + sq(r.z) 
              < sq(dradius[j] - step *0.5* sqrt(dimension))) ;
      if(Rmd) return 1;
      
      Rpd +=  (sq(r.x) + sq(r.y) + sq(r.z) 
              < sq(dradius[j] + step *0.5* sqrt(dimension))) ;
  }
  if(Rpd == 0) return 0;

  double floc = 0;
  if(level == maxlevel) { // If the cell is mixeleveld and the max level of refinement is reached, just return the centered value
    for (int j = 0; j < Nb; j++){
      coord r = dist_perio_coord(center,dpos[j]);
      floc += (- sq(r.x) - sq(r.y) - sq(r.z) + sq(dradius[j]) > 0);
    }
    return floc;
  }else{ // Cell is mixed -> Subdivide it into 8 subcells
    coord subcenter; 
    for (int icx = 0; icx < 2; icx++){
      subcenter.x = center.x + 0.5 * (icx - 0.5)*step;
      for (int icy = 0; icy < 2; icy++){
        subcenter.y = center.y + 0.5 * (icy - 0.5)*step;
        #if dimension == 2
        floc += refine_frac(subcenter, dpos, dradius, step * 0.5, level+1,maxlevel);
        #else 
        for (int icz = 0; icz < 2; icz++){
          subcenter.z = center.z + 0.5 * (icz - 0.5)*step;
          // Recursive call until level max or empty/full is reached
          floc += refine_frac(subcenter, dpos, dradius, step *0.5, level+1,maxlevel);
        }
        #endif
      }
    }
    return floc/pow(2,dimension);
  }
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
    print_header_drop(fDrop);
    fclose(fDrop);
    // files where the fluid and solid phase stats are stored 
    fCA = fopen("fCA.csv","w");
    print_header_CA(fCA);
    fclose(fCA);
    
    fCAnst = fopen("fCAnst.csv","w");
    print_header_CAnst(fCAnst);
    fclose(fCAnst);
  

    // files where the particular averaged quantities are stored 
    fPA = fopen("fPA.csv","w");
    print_header_PA(fPA);
    fclose(fPA);
    #if BI_DISPERSE
    fPA1 = fopen("fPA1.csv","w");
    print_header_PA(fPA1);
    fclose(fPA1);
    fPA2 = fopen("fPA2.csv","w");
    print_header_PA(fPA2);
    fclose(fPA2);
    #endif


    // Continuity wave file header 
    falpha = fopen("falpha.csv","w");
    print_header_layer_files(falpha);
    fclose(falpha);
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
    copy_file("../dump/fCAnst.csv","fCAnst.csv",fCAnst);
    copy_file("../dump/fPA.csv","fPA.csv",fPA);
    copy_file("../dump/falpha.csv","falpha.csv",falpha);
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
  foreach(reduction(+:U) reduction(+:vol))
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
      L.x = (u.x[1,0] + u.x[-1,0]                   \
              + u.x[0,1] + u.x[0,-1] - 4*u.x[] )/(Delta*Delta);\
      L.y = (u.y[1,0] + u.y[-1,0]   \
              + u.y[0,1] + u.y[0,-1] - 4*u.y[] )/(Delta*Delta);\
  }

#define scalar_grad(p,G)\
  {\
      G.x = (p[1,0] - p[-1,0] )/(2.*Delta);\
      G.y = (p[0,1] - p[0,-1] )/(2.*Delta);\
  }

#elif dimension == 3
// gradient centred
#define vec_grad(u,Gv)\
   {\
     Gv.x.x = (u.x[1]     - u.x[-1]    )/(2.*Delta);\
     Gv.x.y = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);\
     Gv.x.z = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);\
     Gv.y.x = (u.y[1]     - u.y[-1]    )/(2.*Delta);\
     Gv.y.y = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);\
     Gv.y.z = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);\
     Gv.z.x = (u.z[1]     - u.z[-1]    )/(2.*Delta);\
     Gv.z.y = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);\
     Gv.z.z = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);\
   }
// laplacien  centred
#define vec_lap(u,L)\
  {\
      L.x = (u.x[1] + u.x[-1]  \
            + u.x[0,1] + u.x[0,-1]  \
            + u.x[0,0,1] + u.x[0,0,-1] - 6*u.x[] )/(Delta*Delta);\
      L.y = (u.y[1] + u.y[-1]  \
            + u.y[0,1] + u.y[0,-1]  \
            + u.y[0,0,1] + u.y[0,0,-1] - 6*u.y[] )/(Delta*Delta);\
      L.z = (u.z[1] + u.z[-1]  \
            + u.z[0,1] + u.z[0,-1]  \
            + u.z[0,0,1] + u.z[0,0,-1] - 6*u.z[] )/(Delta*Delta);\
  }

// grad p centred
#define scalar_grad(p,Gv)\
  {\
      Gv.x = (p[1,0,0] - p[-1,0,0] )/(2.*Delta);\
      Gv.y = (p[0,1,0] - p[0,-1,0] )/(2.*Delta);\
      Gv.z = (p[0,0,1] - p[0,0,-1] )/(2.*Delta);\
  }
#endif