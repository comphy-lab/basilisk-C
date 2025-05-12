/**
##HEADER FILE - O5NS_Flux.h
*/

foreach_dimension(){

static double weno5_left_x (Point point, scalar X, double gradL){
  return ( (1./30.)*(X[-1,0]-2.*Delta*gradL) - (13./60.)*X[-2,0] + (47./60.)*X[-1,0] + (9./20.)*X[0,0] - (1./20.)*X[1,0] );
}

static double weno5_right_x (Point point, scalar X, double gradR){
  return ( (-1./20.)*X[-2,0] + (27./60.)*X[-1,0] + (47./60.)*X[0,0] - (13./60.)*X[1,0] + (1./30.)*(X[0,0]+2.*Delta*gradR) );
}

#if dimension > 1
static double * quadrature3p_2D_x (Point point, face vector WENO_Avg ){

  static double QuadValues[3];

  QuadValues[0] = ( (-9.-22.*sqrt(15.))*WENO_Avg.x[0,2] + (116.+164.*sqrt(15.))*WENO_Avg.x[0,1] + 2186.*WENO_Avg.x[] + (116.-164.*sqrt(15.))*WENO_Avg.x[0,-1] + (-9.+22.*sqrt(15.))*WENO_Avg.x[0,-2] )/2400.;

  QuadValues[1] = (9*WENO_Avg.x[0,2] - 116.*WENO_Avg.x[0,1] + 2134.*WENO_Avg.x[] - 116.*WENO_Avg.x[0,-1] + 9.*WENO_Avg.x[0,-2])/1920.;

  QuadValues[2] = ( (-9.+22.*sqrt(15.))*WENO_Avg.x[0,2] + (116.-164.*sqrt(15.))*WENO_Avg.x[0,1] + 2186.*WENO_Avg.x[] + (116.+164.*sqrt(15.))*WENO_Avg.x[0,-1] + (-9.-22.*sqrt(15.))*WENO_Avg.x[0,-2] )/2400.;

  return(QuadValues);

}
#endif

}

void Tracer_fluxes (scalar f, face vector u, face vector flux_value){
  
  vector gradf[];
  foreach()
    foreach_dimension()
       gradf.x[] = (f[1]-f[-1])/(2.*Delta);
  boundary((scalar *){gradf});

  #if dimension == 1

      foreach_face(){
          if(u.x[] >= 0)
             flux_value.x[] = u.x[]*weno5_left_x  (point,f,gradf.x[-2]);
          else
             flux_value.x[] = u.x[]*weno5_right_x (point,f,gradf.x[ 1]);
      }
     
  #elif dimension == 2

      face vector lineavg[];
      foreach_face(){
          if(u.x[]>=0)
             lineavg.x[] = weno5_left_x  (point,f,gradf.x[-2]);
          else
             lineavg.x[] = weno5_right_x (point,f,gradf.x[ 1]);
      }
      boundary((scalar *){lineavg});

      double *temp;
      double fq[3],uq[3];
      foreach_face(){
         temp = quadrature3p_2D_x (point,lineavg);
         fq[0] = *(temp);
         fq[1] = *(temp+1);
         fq[2] = *(temp+2);
         temp = quadrature3p_2D_x (point,u);
         uq[0] = *(temp);
         uq[1] = *(temp+1);
         uq[2] = *(temp+2);   
         flux_value.x[] = (5.*fq[0]*uq[0] + 8.*fq[1]*uq[1] + 5.*fq[2]*uq[2])/18.;
      }

  #endif

  boundary_flux({flux_value});
}

struct Advection {
  scalar * tracers;
  face vector u;
  double dt;
};

void advection (struct Advection p){
 
  for (scalar f in p.tracers) {
      face vector flux[];
      Tracer_fluxes(f,p.u,flux);
      foreach()
         foreach_dimension()
            f[] += p.dt*(flux.x[]-flux.x[1])/Delta;
   }
   boundary(p.tracers);
}
