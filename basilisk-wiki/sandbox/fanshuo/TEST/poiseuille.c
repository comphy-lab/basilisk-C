#define ML 1
#define HYDRO 1



#include "grid/multigrid1D.h"
#if !ML
# include "saint-venant.h"
#else // ML
# include "layered/hydro.h"
# define phi q
# if !HYDRO
#   include "layered/nh.h"
# endif
# include "layered/remap.h"
# include "layered/perfs.h"
#endif // ML


const double QL = 1., HR = 1., NU = 0.1, BA = 0.4, T0 = 10000;
double Alpha;
FILE *file1;
int main()
{

  periodic(right);
 
  Alpha = 0.349;  
  L0 = 1.;
  G = 1;
  N = 64; 
  
  nu = NU;
  nl = 64; // going to 30 changes very little
#if ML
#if NOMETRIC
  max_slope = 0.;
#endif
#if !HYDRO  
  NITERMIN = 2;
#endif
#endif

  file1 = fopen("profil.dat", "w");
  fclose(file1);

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
 foreach (serial) {
	 foreach_layer(){
          u.x[] = y/DT;
	}
    }


  

  /**
  In the non-hydrostatic case, a boundary condition is required for
  the non-hydrostatic pressure $\phi$ of each layer. */
  
#if !HYDRO && ML
  phi[right] = dirichlet(0.);
#endif
}

event acc(i++){
  foreach (serial) {
	 foreach_layer(){
          u.x[] = u.x[] + 1*sin(Alpha)*dt;
     }
  }
}




/**
We check for convergence. */

double uold,unew;

event logfile (t=1;t += 0.1; t <= T0) {
   unew = 0.;
  foreach(){
    unew = unew + u.x[0,0,nl-1];
  }
  if (i > 1 && fabs(unew-uold) < 1e-5)
    return 1;
  else uold = unew;
}


/**
## Outputs

We save profiles at regular intervals. */

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

event output (t++)
{
  file1 = fopen ("profil.dat", "w");
  foreach (serial) {
    double z = zb[];
#if HYDRO
    printf ("%g %g %g\n", x, z, u.x[]);
    foreach_layer() {
      z += h[];
      fprintf (file1,"%g %g %g\n", x, z, u.x[]);
    }
#elif // ML
    printf ("%g %g %g %g %g\n", x, z, u.x[], w[], phi[]);
    foreach_layer() {
      z += h[];
      printf ("%g %g %g %g %g\n", x, z, u.x[], w[], phi[]);
    }
#endif // ML
    printf ("\n");
    
    fprintf (file1,"\n \n");   
  }
//#if HYDRO
//  delete ({w});
//#endif
  printf ("# end = %g\n", t);
  fclose(file1);
}
/**

~~~gnuplot Velocity and stress profiles for flow
 set xlabel "y"
 set ylabel "u"
p 'profil.dat' u 2:3 w lp t'U comp'
~~~
*/ 