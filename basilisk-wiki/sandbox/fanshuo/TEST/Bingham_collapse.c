/**
# Propose
Simulate the collpase of a rectangle mass having Bingham rheology over inclined plane, using modified multi-layer solver (hydroMT.h). */
/**
# (in process) Analytical solution..

*/
/**
# Code

*/
#define ML 1
#define HYDRO 1
#define MUI 1
#define BINGHAM 1

#include "grid/multigrid1D.h"
# include "hydroMT.h"
# define phi q
# if !HYDRO
#   include "layered/nh.h"
# endif
# include "layered/remap.h"
//# include "layered/perfs.h"



const double HR = 1., NU = 0.1, T0 = 1;
scalar uold[];


/*- - - - - - - - - - - - - - - - - - - - - - - - - */
double Ubgm(Point point)
{
  // fluide released from static
  return 0.;
}


FILE *file1;
/*- - - - - - - - - - - - - - - - - - - - - - - - - */
int main()
{

  X0 = -.5;
  L0 = 1.5;
  G = 1.;
  N = 256*4; 
  nu = NU;
  nl = 64;
  
#if !HYDRO  
  NITERMIN = 2;
#endif
  slope = 1;
/**
Bingham fluide parameters, here $\tau_y=B$, to be consistent with literature  
*/
  tauy =1.25 ;
  mu = 1.;
  
  run();
}


/**
We initialise the topography and the initial thickness of each layer *h*. 
**here $z_b$ non-null is considered**
*/

event init (i = 0)
{
  foreach() {

    zb[] = -(x-X0)*slope;
    
    foreach_layer()
      h[] = (HR)/nl*(fabs(x)<0.25)+dry;

  }
  
//We initialize the velocity 
  foreach (){ 
    foreach_layer()
      u.x[] = Ubgm(point);
    uold[] = 0; 
  }

  u.x[left] = 0.;
  u.x[right] = neumann(0.);
/**
In the non-hydrostatic case, a boundary condition is required for
the non-hydrostatic pressure $\phi$ of each layer. */
  
#if !HYDRO && ML
  phi[right] = dirichlet(0.);
#endif

		file1 = fopen("profil.dat","w");
}
/*
event acc(i++){
  foreach ()
    foreach_layer()
      u.x[] = u.x[] + G*sin(slope)*dt;
}
*/
/** We stop when the difference between two time steps is small enough */
event logfile (t += 0.1; t<=T0) {

  double du = change (u.x,uold);
    if (i > 0 && du < 1e-6)
    return 1; /* stop */
}



/**
For the hydrostatic case, we compute a diagnostic vertical velocity
field `w`. Note that this needs to be done within this event because
it relies on the fluxes `hu` and face heights `hf`, which are only
defined temporarily in the [multilayer solver](hydro.h#update_eta). */

#if HYDRO
scalar w = {-1};

event update_eta (i++)
{
  if (w.i < 0)
    w = new scalar[nl];
  vertical_velocity (w, hu, hf);

  /**
  The layer interface values are averaged at the center of each
  layer. */
  
  foreach() {
    double wm = 0.;
    foreach_layer() {
      double w1 = w[];
      w[] = (w1 + wm)/2.;
      wm = w1;
    }
  }
}
#endif // HYDRO

/**
save mass shape $h(x,t)$ */

event output (t += 0.2;t<=T0)
{
    
  double zz = 0;
  foreach () {
      zz = 0;
      foreach_layer() {  
        zz += h[];
    }

      fprintf (stdout,"%g %g \n", x, zz);
    }
  fprintf (stdout,"\n \n");   
 }
  
/**
# evolution of shape $h(x,t)$
~~~gnuplot 
set xlabel "x"
set ylabel "h(x,t)"
p "out" u 1:2 w l t"0;0.2;0.4;...1"
~~~
*/























