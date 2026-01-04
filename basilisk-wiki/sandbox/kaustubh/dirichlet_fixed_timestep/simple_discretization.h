/**
#Advection solvers

These functions are only redefinition of other functions (fluxes and advection) 
in order to ignore the presence of an embedded boundary.

*/
#include "alex_functions.h"
/**

## Timestep

We need to define a timestep on the mesh without taking into accound the
embedded boundary.
We only take into account cells that will be updated, i.e. such that 
$|dist[]| < maxd$.
We also don't want the previous variable to be initialized at 0 for simple test
cases, we want the timestep to be constant throughout the simulation.
*/
double timestep_LS (vector u, double dtmax,
  scalar dist, double maxd)
{
  dtmax /= CFL;
  foreach(reduction(min:dtmax))
    foreach_dimension(){
      if (u.x[] != 0. && fabs(dist[]) < maxd) {
        double dt = Delta/fabs(u.x[]);
        if (dt < dtmax) dtmax = dt;
      }
    }

  return dtmax *= CFL;
}


/**
## Advection

We use a simple advection scheme with co-located variables.
*/

#include "weno2.h"

void FE(scalar f, scalar fi, vector ndist, double dt, 
  double maxd){
/**
Simple Forward Euler scheme.
*/

  foreach(){
    if(fabs(f[]) < maxd)
      foreach_dimension(){
        double graplus  = WENOdiff_x(point, fi,-1);
        double graminus = WENOdiff_x(point, fi,1);
        f[] -= dt * ( max(ndist.x[],0.)*graplus
          + min(ndist.x[],0.)*graminus);
    }
  }
  boundary({f});
  restriction({f});
}
/**
Runge Kutta 2 scheme
*/
void RK2(scalar f, vector ndist, double dt, double maxd){
  scalar f1[],f2[];

  foreach(){
    f1[] = f[];
  }
  boundary({f1});
  restriction({f1});

  FE(f1, f, ndist, dt, maxd);
  
  foreach(){
    f2[] = f1[];
  }
  boundary({f2});
  restriction({f2});

  FE(f2, f1, ndist, dt, maxd);
  foreach(){
    f[] = (f2[] + f[])*0.5;
  }
  boundary({f});
  restriction({f});
}


/**
Runge Kutta 3 advection scheme
*/
void RK3(scalar f, vector ndist, double dt, double maxd){
  static const double coeff[2][2] = {
    {0.75, 0.25},
    {1./3.,2./3.}
  };
  scalar f1[],f2[];

  foreach(){
    f1[] = f[];
  }
  boundary({f1});
  restriction({f1});

  FE(f1, f, ndist, dt, maxd);
  foreach(){
    f2[] = f1[];
  }
  boundary({f2});
  restriction({f2});

  FE(f2, f1, ndist, dt, maxd);
  foreach(){
    f1[] = coeff[0][0]*f[]+coeff[0][1]*f2[];
    f2[] = f1[];
  }
  boundary({f1,f2});
  restriction({f1,f2});
  FE(f2, f1, ndist, dt, maxd);
  foreach(){
    f[] = coeff[1][0]*f[] + coeff[1][1]*f2[];
  }
  boundary({f});
  restriction({f});

}

/**
Runge Kutta 3 advection scheme
*/

#include "weno5.h"
void RK3_WENO5(scalar f, vector ndist, double dt, double maxd){
  static const double coeff[2][2] = {
    {0.75, 0.25},
    {1./3.,2./3.}
  };
  scalar f1[],f2[];

  foreach(){
    f1[] = f[];
  }
  boundary({f1});
  restriction({f1});

  FE_WENO5(f1, f, ndist, dt, maxd);
  foreach(){
    f2[] = f1[];
  }
  boundary({f2});
  restriction({f2});

  FE_WENO5(f2, f1, ndist, dt, maxd);
  foreach(){
    f1[] = coeff[0][0]*f[]+coeff[0][1]*f2[];
    f2[] = f1[];
  }
  boundary({f1,f2});
  restriction({f1,f2});
  FE_WENO5(f2, f1, ndist, dt, maxd);
  foreach(){
    f[] = coeff[1][0]*f[] + coeff[1][1]*f2[];
  }
  boundary({f});
  restriction({f});

}

void RK3_WENO3(scalar f, vector ndist, double dt, double maxd){
  static const double coeff[2][2] = {
    {0.75, 0.25},
    {1./3.,2./3.}
  };
  scalar f1[],f2[];

  foreach(){
    f1[] = f[];
  }
  boundary({f1});
  restriction({f1});

  FE_WENO3(f1, f, ndist, dt, maxd);
  foreach(){
    f2[] = f1[];
  }
  boundary({f2});
  restriction({f2});

  FE_WENO3(f2, f1, ndist, dt, maxd);
  foreach(){
    f1[] = coeff[0][0]*f[]+coeff[0][1]*f2[];
    f2[] = f1[];
  }
  boundary({f1,f2});
  restriction({f1,f2});
  FE_WENO3(f2, f1, ndist, dt, maxd);
  foreach(){
    f[] = coeff[1][0]*f[] + coeff[1][1]*f2[];
  }
  boundary({f});
  restriction({f});

}
