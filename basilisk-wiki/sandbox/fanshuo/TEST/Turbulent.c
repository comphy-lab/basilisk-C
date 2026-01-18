//To be completed,,, 2D flow velocity profil with a turbulent viscosity 

#define ML 1
#define HYDRO 1
#define MUI 1
#define TURBULENT  1

#include "grid/multigrid1D.h"
#if !ML
# include "saint-venant.h"
#else // ML
#include "hydroMT.h"

# define phi q
# if !HYDRO
#   include "layered/nh.h"
# endif
# include "layered/remap.h"
#endif // ML


const double NU = 0.1, T0 = 10000, HR = 1.;
scalar uold[];

double Ub(Point point)
{
  double Ia = 0.; 
  double zc = 0.;
  for (int l = - point.l; l < nl - point.l; l++) {
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;

  Ia = (tan(slope)-mu0)*I0/(mu0+deltamu-tan(slope));
 
  return 3./3.*Ia/dg*sqrt(G*pow(HR,3)*cos(slope))*(1-pow((1-zc/HR),1.5));
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - */
int main()
{
  periodic(right);
 
  slope = 0.404;  // ~23 degrees
  L0 = 1;
  G = 1.;
  N = 32; 
  
  
  nu = NU;
  nl = 256; 


  run();
}

/**
We initialise the topography and the initial thickness of each layer *h*. */

event init (i = 0)
{
  foreach() {
    zb[] = 0.;
#if !ML
    h[] = HR - zb[];
#else
    foreach_layer()
      h[] = (HR)/nl;
#endif
 	eta[] = HR;  
  }

//Mettre une vitesse initiale dans le domaine 
  foreach () {
    foreach_layer()
      u.x[] = Ub(point);
    uold[] = 0;
  }
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
## Outputs

We save profiles at regular intervals. */
#if 0
event profiles (t += 5;t<=T0)
{
  foreach (serial) {
#if !ML
    double H = h[];
#else
    double H = 0.;
    foreach_layer()
      H += h[];
#endif
    fprintf (stderr, "%g %g %g\n", x, zb[] + H, zb[]);
  }
  fprintf (stderr, "\n\n");
}
#endif

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

event output (t = end)
{
 
    foreach () {
      double z = zb[];
    foreach_layer() {
      z += h[];
      fprintf (stdout, "%g %g %g \n", x, z, u.x[]);
    }
  }
#if HYDRO
  delete ({w});
#endif
 
}



/**

~~~gnuplot Velocity and shear profiles for Turbulent flow
 set xlabel "y"
 set ylabel "u"
 p "out" u 2:3 w p t'U computed', 
~~~
*/ 























