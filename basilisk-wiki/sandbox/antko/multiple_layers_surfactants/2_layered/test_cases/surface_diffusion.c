/**
## MWE for surface diffusion
In this test case, the interface is stationary and the velocity field at zero. Surface concentration is initialized nonuniformly on the interface, and it diffuses along the interface to become uniform. The interface is a straight line and diffusion is purely horizontal. The exact solution is :
$$
\Gamma(x,t) = \frac{1}{2} + \frac{1}{2} \exp \left(\frac{- \xi^2}{Pe_S}\right) \cos(\xi x)
$$
*/

double initial_mass = 0.;
#include "grid/multigrid1D.h"
#include "../hydro_c.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "../surface.h"
#include "../solute.h"

/**
## Geometry and resolution*/
#define L 1.
#define h_ 0.5
#define LEVEL 6
#define layers 10
#define T_END 0.2
#define DELTA_T (T_END/500.)
#define MY_TOLERANCE 1e-11

/**
## Constants*/
#define k_  (2.*pi)

/**
## Physical parameters*/
#define g_ 0.
#define Re 1.
#define solute0 0.
#define surface0 1.

/**
##Main function
There is no advection, no adsorption, no viscosity, no gravity... but there is horizontal diffusion !
*/
int main()
{
  L0 = L;
  origin (0.);
  periodic (right);
  N = 1 << LEVEL;
  TOLERANCE = MY_TOLERANCE;
  nl = layers;
  G = g_;
  D_hor = 1.;
  run();
}

/**
## Initialisation 
A flat surface with a non-uniform surface concentration is initialised.*/
event init (i = 0)
{
  scalar h, c;
  for (h,c in hl,cl) {
    foreach() {
      double H = h_;
      h[] = H/nl;
  	  c[] = solute0;
    }
    boundary({c});
  }
  interface_area(hl, area);
  foreach()
    surface[] = (0.5 + 0.5 * cos(k_*x))*area[];
}

/**
## Outputs*/
event outputs(t += DELTA_T; t = DELTA_T; t <= T_END) {
  static FILE * fpSurf = fopen("dataSurf.txt", "w");
  static FILE * fp = fopen("surf_mass.txt", "w");
  static FILE * fpMax = fopen("surf_max.txt", "w");

  double total_solute = 0 ;
  double total_surface = 0 ;
  double Gamma = 0;
  foreach (){
    total_surface += surface[]*Delta;
    Gamma = surface[]/area[];
    fprintf (fpSurf, "%g %g %.17g\n", t, x, Gamma);
  }
  fprintf(fpSurf, "\n");
  fprintf(fpMax, "%g %g %g\n", t, Gamma - 0.5, Gamma - 0.5 - 0.5*exp(-sq(k_)*t*D_hor));
  
  initial_mass = (initial_mass == 0) ? total_solute + total_surface : initial_mass;
  fprintf (fp, "%g %g %g %g\n", t, total_solute, total_surface, total_solute + total_surface - initial_mass);
}

/**
#Results
The total mass of surfactant is well conserved. The error is the machine accuracy.

~~~gnuplot Difference between total surface mass and initial surface mass
set terminal @PNG enhanced size 640,640 font ",8"
set output 'sum.png'
set xlabel "t"
set ylabel "Mass error"
set key left
plot \
  './surf_mass.txt' u 1:4 w l
~~~

Surfactants diffuse over the surface and tend to homogenise the surface concentration. The form of the curves is supposed to be as:
$$
\Gamma(x,t) = \frac{1}{2} + \frac{1}{2} \exp \left(\frac{- \xi^2}{Pe_S}\right) \cos(\xi x)
$$

~~~gnuplot Surface concentration during diffusion over x for different time step
set output 'surface.png'
set xlabel "x"
set ylabel "Gamma"
set key left

plot \
  './dataSurf.txt' u 2:3 w l
~~~

To compare the simulation with the analytical, we draw $\Gamma - \frac{1}{2}$ and the difference with $\frac{1}{2} \exp \left(\frac{- \xi^2}{Pe_S}\right)$

~~~gnuplot Surface concentration at x = 1 over time 
set output 'rescaled.png'
set xlabel "t"
set ylabel "Gamma"
set key left

plot \
  './surf_max.txt' u 1:2 t '\sim' w l
~~~

~~~gnuplot Concentration error over time
set output 'error.png'
set xlabel "t"
set ylabel "Concentration error"
set key left

plot \
  './surf_max.txt' u 1:3 t 'Error' w l
~~~

The error at short times can be quite high and is linear with the time step.
*/
