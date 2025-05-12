/**
This is not a bug: periodic boundaries are indeed not "real" boundaries.

Here we try to compute the flowrate at the right boundary's section of
the tri-periodic (cubic) domain of a flow driven by a pressure drop.

Problem: when using the foreach_boundary(right) syntax (i.e with the
compute_flowrate_at_boundary_right function defined bellow), there is
no cell transversal and thus the flowrate remains 0.

When using a manual detection of the cells neighboring the right
boundary, (i.e with compute_flowrate_basic_right function defined
bellow) one can see that there is indeed a flowrate but this solution
is not satisfactory as the location of the plane is not the
boundary's. */

#define LEVEL 6
#include "grid/octree.h"
#include "navier-stokes/centered.h"

#define rhoval 1000.
#define muval 60.

#define mydt 0.00375
#define maxiter (20)

double deltau;
scalar un[];

int main() {		
  
  L0 = 1; 
  
  DT = mydt;
  
  init_grid(1 << (LEVEL));


  /* Tri-periodic flow driven by body force */ 
  
  foreach_dimension()
    periodic(top);

  const face vector dp[] = {10./L0/rhoval, 0, 0};
  a = dp;
  
  run();
}


event init (i = 0) {
  origin (0, 0, 0);
  
/* set dynamic viscosity */
  const face vector muc[] = {muval, muval, muval};
  mu = muc;

/* set density of the flow */ 
  const scalar rhoc[] = rhoval;
  rho = rhoc;
  
/* The flow is at rest initially. */
  foreach(){
    u.x[] = 0.;
    un[] = u.x[];
  }
}


void compute_flowrate_at_boundary_right(const double Vol, coord * vel) {

  foreach_dimension()
    vel->x = 0.;

  boundary(all);
  int k =0;
  foreach_boundary(right){
  k ++;
    if (k == 1) printf("thread %d, x = %g",pid(),x);
    foreach_dimension()
      vel->x += sq(Delta)*u.x[];
  }
  
  foreach_dimension()
    vel->x /= Vol;
  
#if _MPI
  foreach_dimension(){
    mpi_all_reduce(vel->x, MPI_DOUBLE, MPI_SUM);
  }
#endif


}


void compute_flowrate_basic_right(const double Vol,  coord * supvelo){

  Cache intDomain = {0};
  Point lpoint;
  double xboundary = L0;

 
  // allocate cache for the integration domain (one slice of x, 2D domain)
   
  foreach(){
    if (x > (xboundary - Delta)) 
      { 
	lpoint = locate(x, y, z);
	cache_append(&intDomain, lpoint, 0);

      }
  }

  cache_shrink(&intDomain); 
  printf("thread %d, Number of point in integration domain %d\n",pid(),intDomain.n);

  foreach_dimension()
    supvelo->x = 0.;

  int k = 0;

  foreach_cache(intDomain){
    k ++;
    if (k == 1) printf("thread %d, x = %g",pid(),x);
    foreach_dimension()
      supvelo->x += sq(Delta)*u.x[];
  }
 
  foreach_dimension()
    supvelo->x /= Vol;

  free(intDomain.p);

#if _MPI
  foreach_dimension(){
    mpi_all_reduce(supvelo->x, MPI_DOUBLE, MPI_SUM);
  }
#endif
   
}

event flowrate(i++;i<maxiter){

  coord flowrate;
  double Volume = pow(L0,3);

  compute_flowrate_at_boundary_right(Volume, &flowrate);
  /* compute_flowrate_basic_right(Volume, &flowrate); */

  fprintf(ferr,"%g %g %g %g %g\n",t, flowrate.x, flowrate.y, flowrate.z, Volume);
  
}

/** 
~~~gnuplot flowrate and time
plot 'log' u 1:2 w l
~~~
*/
