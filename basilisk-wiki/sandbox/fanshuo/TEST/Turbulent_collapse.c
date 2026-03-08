
#define ML 1
#define HYDRO 1
#define MUI 1
#define TURBULENT  1

#include "grid/multigrid1D.h"
#include "hydroMT.h"

# define phi q
# if !HYDRO
#   include "layered/nh.h"
# endif
# include "layered/remap.h"
//# include "layered/perfs.h"

const double NU = 0.1, T0 = 3, HR = 1.;
scalar uold[];

double Ut(Point point)
{
  return 0.;
}
/*- - - - - - - - - - - - - - - - - - - - - - - - - */
int main()
{
 
  slope = 1;  
  L0 = 2;
  X0 = -0.5;
  G = 1.;
  N = 256*4; 
  nu = 1;
  nl =64; 
  
#if ML
#if NOMETRIC
  max_slope = 0.;
#endif
#if !HYDRO  
  NITERMIN = 2;
#endif
#endif

  rugo = 0.1;
  run();
}

/**
We initialise the topography and the initial thickness of each layer *h*. */

event init (i = 0)
{
  foreach() {
#if !ML
    h[] = HR - zb[];
#else
    foreach_layer()
      h[] = (HR)/nl*(fabs(x)<0.05)+dry;
#endif
  }

//Mettre une vitesse initiale dans le domaine 
  foreach () {
    foreach_layer()
      u.x[] = Ut(point);
    uold[] = 0;
  }
		u.x[right] = neumann(0.);

		u.x[left] = neumann(0.);
}

event acc(i++){
  foreach () 
    foreach_layer()
      u.x[] = u.x[] + G*sin(slope)*dt;
}


// /**exit criteria  */
event logfile (t += 1; t<=T0) {

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
We also save the velocity and non-hydrostatic pressure fields. */
event output (t += 0.5;t<=T0 )
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

~~~gnuplot 
 set xlabel "y"
 set ylabel "h(x,t)"
 p "out" u 1:2 w l t"0,0.2,0.4,...,4"
~~~
*/ 























