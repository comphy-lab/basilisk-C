/**
## MWE for surface diffusion
In this test case, the interface is stationary and the velocity field at zero.
 Surface concentration is initialized nonuniformly on the interface, and it
 diffuses along the interface to become uniform. The interface is a straight
 line and diffusion is purely horizontal. The exact solution is :
$$
\Gamma(x,t) = 1. + \exp({-D_S*\xi^2 * t}) \cos(\xi x)
$$ 
*/
double Nt = 64;
double max_error = 0.;
double mass_init = 0.;
#include "grid/multigrid1D.h"
#include "crobert/2_Implicit/hydro-tension.h"
#include "layered/nh.h"
#include "crobert/2_Implicit/surface.h"
#include "crobert/2_Implicit/solute.h"
#include "layered/remap.h"

/**
## Geometry and resolution*/
#define h_ 0.5
#define T_END 1
#define DELTA_T (T_END/10.)
#define MY_TOLERANCE 1e-11
#define k_  (2.*pi)
#define D 1.

/**
##Main function
There is no advection, no adsorption, no viscosity, no gravity... 
 but there is horizontal diffusion !
*/
int main()
{
  L0 = 1.;
  origin (0.);
  periodic (right);
  TOLERANCE = MY_TOLERANCE;
  nl = 10;
  G = 0.;
  D_surf = D;
  for (N = 8; N <= 128; N *= 4) {
    for (Nt = 64; Nt <= 512; Nt *= 2) {
      mass_init = 0;
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
      c[] = (1. + cos(k_*x));
      eta[] += h[];
    }
  }
  boundary(all);
  foreach()
    M[] = (1. + cos(k_*x));
}

/**
## Outputs*/
event outputs(t += DELTA_T/Nt; t = 0.; t <= T_END) {
  char name[80];
  double mass_tot = 0;
  double val = 0;
  double th = 0;
  double error_step = 0;
  
  foreach (){
    mass_tot += M[]*Delta;
    val = M[] - 1;
    th = cos(k_*x) * exp(-sq(k_)*t*D);
    error_step = max(error_step, (val - th));
  }
  if (mass_init == 0)
    mass_init = mass_tot;
  max_error = max(max_error, error_step);
  
  if(N == 32 && Nt ==64) {
    sprintf (name, "dataSurf-N-%d-Nt-%g", N, Nt);
    static FILE * fpSurf = fopen(name, "w");
    foreach ()
      fprintf (fpSurf, "%g %g %g %g\n", t, x, M[], c[0,0,nl-1]);   
    fprintf(fpSurf, "\n"); 
  }

  if(N == 128 && Nt ==512) {
    sprintf (name, "surf_mass-N-%d-Nt-%g", N, Nt);
    static FILE * fp = fopen (name, "w");
    fprintf(fp, "%g %g %g %g %g\n", t, mass_tot - mass_init, val, th, error_step);
  }
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
set terminal @SVG enhanced size 640,640 font ",8"
set xlabel "t"
set ylabel "Mass error"
set key left
plot \
  './surf_mass-N-128-Nt-512' u 1:2 w l
~~~

Surfactants diffuse over the surface and tend to homogenise the surface 
concentration. The form of the curves is supposed to be as:
$$
\Gamma(x,t) = 1. + \exp \left(\frac{- \xi^2}{Pe_S}\right) \cos(\xi x)
$$

~~~gnuplot Surface concentration during diffusion over x for different time step
set xlabel "x"
set ylabel "Gamma"
set key left

plot \
  './dataSurf-N-32-Nt-64' u 2:3 w l
~~~

To compare the simulation with the analytical, we draw $\Gamma - 1$ 
and the difference with $\exp \left(\frac{- \xi^2}{Pe_S}\right)$.
The simulation fits well the analytical model with errors of the magnitude 
of $10^{-3}. The maximal error is observed in the first time steps of the 
simulation.

~~~gnuplot Surface concentration at x = 1 over time 
set xlabel "t"
set ylabel "Gamma"
set key left

plot \
  './surf_mass-N-128-Nt-512' u 1:3 t '\sim' w l, \
  './surf_mass-N-128-Nt-512' u 1:4 t '\th' w l
~~~

~~~gnuplot Relative concentration error over time compared to the initial amplitude
set xlabel "t"
set ylabel "Relative concentration error"
set key left

plot \
  './surf_mass-N-128-Nt-512' u 1:5 t 'Error' w l
~~~

~~~gnuplot Relative error for the maximum height as a function of the frequency.
set xlabel 'Number of timesteps'
set ylabel 'Relative error'
set logscale y
set logscale x 2
set grid
plot [64:512][1e-4:1e-1]'error.txt' every ::0::5 u 1:2 t "N = 8" w p ps 2,\
  'error.txt' every ::4::8 u 1:2 t "N = 32" w p ps 2,\
  'error.txt' every ::8::12 u 1:2 t "N = 128" w p ps 2
~~~

The maximal error (short times) is strongly dependant of the timestep and 
decrease linearly with it. On the contrary, it is independant of the 
resolution or of the number of layers.
*/
