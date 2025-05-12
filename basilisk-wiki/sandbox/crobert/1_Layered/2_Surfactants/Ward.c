/**
## MWE for Ward & Tordai
In this test case, the interface is flat and stationary, the velocity field at zero. Surface concentration is null. Bulk concentration is initialized uniformly. Surfactants adsorb on the surface to reach the equilibrium. The adsorption flux is compensed by vertical diffusion, which is limitant. The exact solution is a convolution function:
*/

#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../surfactants.h"

/**
## Geometry and resolution*/
#define L 1.
#define h_ 5.
#define LEVEL 1
#define layers 32
#define T_END 2.
#define DELTA_T (T_END/100.)
#define MY_TOLERANCE 1e-11

/**
## Physical parameters*/
#define g_ 0.
#define Re 1.
#define solute0 1.
#define surface0 0.
#define Diff 1.

/**
##Main function
There is no advection, no viscosity, no gravity... but there is adsorption and vertical diffusion !
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
#ifdef SIGMA
  M.K = SIGMA;
#else
  M.K = 1.;
#endif
  M.omega = 100./(M.K);
  D_vert = Diff;

  run();
}

/**
## Initialisation 
A flat surface with a empty surface and a uniform bulk concentration is initialised.*/
event init (i = 0)
{
  foreach() {
    double H = h_;
    foreach_layer() {
      h[] = H/nl;
	  	c[] = solute0;
    }
  }
}


/**
## Outputs*/
event profiles (i++) 
{
  static FILE * fpW = fopen("dataW.txt", "w");
  foreach()
   fprintf (fpW, "%g %g %g\n", t, M[]/area_m(point), c[0,0,nl-1]);
}

double initial_mass = 0.;
event outputs(t += DELTA_T; t = DELTA_T; t <= T_END) {
  double total_solute = 0 ;
  double total_surface = 0 ;
  foreach (){
    foreach_layer ()
      total_solute += h[]*c[]*Delta;
    total_surface += M[]*Delta;
  }
  static FILE * fp = fopen("surf_mass.txt", "w");
  initial_mass = (initial_mass == 0) ? total_solute + total_surface : initial_mass;
  fprintf (fp, "%g %g %g %g %g\n", t, total_solute, total_surface, total_solute + total_surface - initial_mass, M.K);
  
}

/**
#Results
The total mass of surfactant is well conserved. The error is the machine accuracy.

~~~gnuplot Difference between total mass and initial total mass
set terminal @PNG enhanced size 640,640 font ",8"
set output 'sum.png'
set xlabel "t"
set ylabel "Mass error"
set key left
plot \
  './surf_mass.txt' u 1:4 w l
~~~

Surface and sucsurface concentrations are drawn. During the first phase, adsorption is quick, subsurface concentration decrease quickly and then re-increase when local equilibrium is reached and a part of the diffusion flux accumulates in the cell.

~~~gnuplot Surface concentration over t
set output 'surface.png'
set xlabel "t"
set ylabel "Gamma"
set key left
plot \
  '../Ward1/dataW.txt' u 1:2 t 'c' w l, \
  '../Ward5/dataW.txt' u 1:2 t 'c' w l, \
  '../Ward10/dataW.txt' u 1:2 t 'c' w l
~~~

~~~gnuplot Subsurface concentration over t
set output 'subsurface.png'
set ylabel "c"
set key left
plot \
  '../Ward1/dataW.txt' u 1:3 t 'c' w l,\
  '../Ward5/dataW.txt' u 1:3 t 'c' w l,\
  '../Ward10/dataW.txt' u 1:3 t 'c' w l
~~~
*/
