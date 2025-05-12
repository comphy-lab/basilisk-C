/**
## MWE for surface diffusion
In this test case, the interface is stationary and the velocity field at zero. Surface concentration is initialized nonuniformly on the interface, and it diffuses along the interface to become uniform. The interface is a straight line and diffusion is purely horizontal. The exact solution is :
$$
\Gamma(x,t) = \frac{1}{2} + \frac{1}{2} \exp \left({- D_S * \xi^2 * t}\right) \cos(\xi x)
$$ 
*/
double Nt = 64;
double max_error = 0.;
double initial_mass = 0.;
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../solute.h"
#include "../nh_PhiS.h"
#include "../surface.h"

/**
## Geometry and resolution*/
#define L 1.
#define h_ 0.5
#define LEVEL 4
#define layers 10
#define T_END 1
#define DELTA_T (T_END/10.)
#define MY_TOLERANCE 1e-11

/**
## Constants*/
#define k_  (2.*pi)

/**
## Physical parameters*/
#define g_ 0.
#define Re 1.

/**
##Main function
There is no advection, no adsorption, no viscosity, no gravity... but there is horizontal diffusion !
*/
int main()
{
  L0 = L;
  origin (0.);
  periodic (right);
  TOLERANCE = MY_TOLERANCE;
  nl = layers;
  G = g_;
  M.D_surf = 1.;
  for (N = 8; N <= 128; N *= 4) {
    for (Nt = 64; Nt <= 1024; Nt *= 2) {
      initial_mass = 0;
      max_error = 0.;
      run();
    }
  }
}

/**
## Initialisation 
A flat surface with a non-uniform surface concentration is initialised.*/
event init (i = 0)
{
  foreach() {
    eta[] = zb[];
    foreach_layer() {
      h[] = h_/nl;
      c[] = (0.5 + 0.5 * cos(k_*x));
      eta[] += h[];
    }
  }
  boundary({h, eta, c});
  foreach()
    M[] = (0.5 + 0.5 * cos(k_*x))*area_m(point);
  boundary({M});
}

/**
## Outputs*/
event outputs(t += DELTA_T/Nt; t = 0.; t <= T_END) {
  char name[80];
  sprintf (name, "surf_mass-N-%d-Nt-%g", N, Nt);
  static FILE * fp = fopen (name, "w");
  sprintf (name, "dataSurf-N-%d-Nt-%g", N, Nt);
  static FILE * fpSurf = fopen(name, "w");

  double total_surface = 0 ;
  double Gamma = 0;
  foreach (){
    total_surface += M[]*Delta;
    Gamma = M[]/area_m(point);
    fprintf (fpSurf, "%g %g %g %g\n", t, x, Gamma, c[0,0,nl-1]);
  }
  fprintf(fpSurf, "\n");
  initial_mass = (initial_mass == 0) ? total_surface : initial_mass;
  fprintf(fp, "%g %g %g %g %g\n", t, total_surface - initial_mass, Gamma - 0.5, 0.5*exp(-sq(k_)*t*M.D_surf), (Gamma - 0.5 - 0.5*exp(-sq(k_)*t*M.D_surf))/0.5);
  max_error = max(max_error, (Gamma - 0.5 - 0.5*exp(-sq(k_)*t*M.D_surf))/0.5);
}

event error (t = end)
{
 static FILE * err = fopen ("error.txt", "w");
  fprintf (err, "%g %g\n", Nt/T_END, max_error);
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
  './surf_mass-N-128-Nt-1024' u 1:2 w l
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
  './dataSurf-N-128-Nt-1024' u 2:3 w l
~~~

To compare the simulation with the analytical, we draw $\Gamma - \frac{1}{2}$ and the difference with $\frac{1}{2} \exp \left(\frac{- \xi^2}{Pe_S}\right)$.
The simulation fits well the analytical model with errors of the magnitude of $10^{-3}. The maximal error is observed in the first time steps of the simulation.

~~~gnuplot Surface concentration at x = 1 over time 
set output 'rescaled.png'
set xlabel "t"
set ylabel "Gamma"
set key left

plot \
  './surf_mass-N-128-Nt-1024' u 1:3 t '\sim' w l, \
  './surf_mass-N-128-Nt-1024' u 1:4 t '\th' w l
~~~

~~~gnuplot Relative concentration error over time compared to the initial amplitude
set output 'relative error.png'
set xlabel "t"
set ylabel "Relative concentration error"
set key left

plot \
  './surf_mass-N-128-Nt-1024' u 1:5 t 'Error' w l
~~~

~~~gnuplot Relative error for the maximum height as a function of the frequency.
set output 'Error.png'
set xlabel 'Number of timesteps'
set ylabel 'Relative error'
set logscale y
set logscale x 2
set grid
plot [64:1024][1e-4:1e-1]'error.txt' every ::0::5 u 1:2 t "N = 8" w p,\
  'error.txt' every ::5::10 u 1:2 t "N = 32" w p,\
  'error.txt' every ::10::15 u 1:2 t "N = 128" w p
~~~

The maximal error (short times) is strongly dependant of the timestep and decrease linearly with it. 
On the contrary, it is independant of the resolution or of the number of layers.
*/
