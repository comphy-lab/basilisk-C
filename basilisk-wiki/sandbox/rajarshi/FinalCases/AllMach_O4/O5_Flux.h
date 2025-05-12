/**
# Solver for computing weno fluxes.
*/

foreach_dimension(){

double WENO5_Reconstruction_Left_x (Point point, scalar X, double gradL){
  return ( (1./30.)*(X[-1,0]-2.*Delta*gradL) - (13./60.)*X[-2,0] + (47./60.)*X[-1,0] + (9./20.)*X[0,0] - (1./20.)*X[1,0] );
}

double WENO5_Reconstruction_Right_x (Point point, scalar X, double gradR){
  return ( (-1./20.)*X[-2,0] + (27./60.)*X[-1,0] + (47./60.)*X[0,0] - (13./60.)*X[1,0] + (1./30.)*(X[0,0]+2.*Delta*gradR) );
}

double * QuadratureValues2D_x (Point point, face vector WENO_Avg ){

  static double QuadValues[2];
  QuadValues[0] = 2.*(((-1.+70.*sqrt(3.))/8640.)*WENO_Avg.x[0,-2] + ((4.-500.*sqrt(3.))/8640.)*WENO_Avg.x[0,-1] + (4314./8640.)*WENO_Avg.x[0,0] + ((4.+500.*sqrt(3.))/8640.)*WENO_Avg.x[0,1] + ((-1.-70.*sqrt(3.))/8640.)*WENO_Avg.x[0,2]);
  QuadValues[1] = 2.*(((-1.-70.*sqrt(3.))/8640.)*WENO_Avg.x[0,-2] + ((4.+500.*sqrt(3.))/8640.)*WENO_Avg.x[0,-1] + (4314./8640.)*WENO_Avg.x[0,0] + ((4.-500.*sqrt(3.))/8640.)*WENO_Avg.x[0,1] + ((-1.+70.*sqrt(3.))/8640.)*WENO_Avg.x[0,2]);
  return(QuadValues);
}

double Transverse_Derivative_x (Point point, face vector WENO_Avg ){
  return( (WENO_Avg.y[-2,0] - 8.*WENO_Avg.y[-1,0] + 8.*WENO_Avg.y[1,0] - WENO_Avg.y[2,0])/(48.*Delta) );
}

double Normal_Derivative_x (Point point, scalar V){
  return( (V[-2,0] - 15.*V[-1,0] + 15.*V[] - V[1,0] )/(24.*Delta) );
}

}

void tracer_fluxes (scalar f, face vector u, face vector flux_value){

  vector grad_f[];
  foreach()
    foreach_dimension()
        grad_f.x[] = (f[1]-f[-1])/(2.*Delta);
  boundary((scalar *){grad_f}); 
  
  #if dimension == 1

      foreach_face(){
          if(u.x[] >= 0)
              flux_value.x[] = u.x[]*WENO5_Reconstruction_Left_x  (point,f,grad_f.x[-2,0]);
          else
              flux_value.x[] = u.x[]*WENO5_Reconstruction_Right_x (point,f,grad_f.x[ 1,0]);
       }
 
  #elif dimension == 2

      face vector Weno_AvgLine[];
      foreach_face(){  
          if(u.x[] >= 0)
              Weno_AvgLine.x[] = WENO5_Reconstruction_Left_x  (point,f,grad_f.x[-2,0]);
          else
              Weno_AvgLine.x[] = WENO5_Reconstruction_Right_x (point,f,grad_f.x[ 1,0]);
       }
      boundary((scalar *){Weno_AvgLine});

      double *temp;
      double fq[2],uq[2];
      foreach_face(){
         temp = QuadratureValues2D_x (point,Weno_AvgLine);
         fq[0] = *(temp);
         fq[1] = *(temp+1);
         temp = QuadratureValues2D_x (point,u);
         uq[0] = *(temp);
         uq[1] = *(temp+1); 
         flux_value.x[] = (fq[0]*uq[0] + fq[1]*uq[1])/2.;
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
      tracer_fluxes(f,p.u,flux);
      foreach()
         foreach_dimension()
            f[] += p.dt*(flux.x[]-flux.x[1])/Delta;
   }
   boundary(p.tracers);
}