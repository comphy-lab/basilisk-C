/**
# Vanishing waves 

This test looks on the behaviour of a liquid film when an obstacle on the ground
moves at uniform velocity V. We take the referential of the obstacle so the test 
simulates a uniform flow at velocity V on a perturbated ground in presence of 
gravity and surface tension.

After the obstacle a wave can be produced if it verifies the condition $c_e = V$.
However the medium is dispersive for waves with the relation : 

$$
c_e = \sqrt{(\frac{G}{k} + \gamma k)\tanh(kh)}
$$ 

For sufficiently high velocity, two waves can be created with different wave number : 
either capillary or gravity waves. However, for low velocities, no wave should appear.

In our situation, $H_0 = 2$, $G=2*, $\sigma = 1$, the minimum velocity is around 1.7

![V=1.6](vanishing_waves/movie-V1.6-sig1.mp4)
For $ V = 1.6$, no wave is observed.

![V=1.8](vanishing_waves/movie-V1.6-sig1.mp4)
For $V = 1.8$, gravity waves are clearly oberved. We are still looking for the capillary waves.
*/

/**
## Include and parameters
This test uses the multilayer solver and the surface tension additive. */
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "../nh_PhiS.h"
#include "layered/remap.h"
#include "../tension.h"

double V = 1.;
double T_END = 10.;
double sig = 1.;
/**
## Initialisation */
#define H0 2
int main()
{
  L0 = 150;
  X0 = -L0/2;
  N = 1024;
  nl = 20;
  G = 2.;

  V = 1.6;
  T_END = 2.*L0/V;
  run();

  V = 1.8;
  T_END = 3.*L0/V;
  run();
  
  V = 2.2;
  T_END = 3.*L0/V;
  run();
}

/**
The initial condition is a flat film of thickness H
and with uniform velocity V*/

event init (t = 0) {
  foreach () {
    zb[] = H0/50.*exp(-sq(x)) - H0;
    double H = - zb[];
    foreach_layer () {
      h[] = H/nl;
      u.x[] = V * H0/H;
      hu.x[] = V * h[];
    }
    sigma[] = sig;
  }
//Radiation BC seems not to work for this case (even with graviy only)
  u.n[left] = dirichlet(V);
  eta[left] = dirichlet(0.);

  u.n[right] = neumann(0.);
  eta[right] = dirichlet(0.);

  boundary(all);
}

/**
# Output
Output the final surface profile */
event logfile(t = end) {
  char name[80];
  sprintf (name, "profile-V-%g", V);
  static FILE * fpp = fopen (name, "w");
  foreach() {
    double H = 0;
    double Q = 0;
    foreach_layer() {
      H += h[];
      Q += h[]*u.x[];
    }
    fprintf (fpp, "%g %g %g\n", x, H, Q);  
  }
}

/**# Movie
Plot the evolution in time with the max speed*/
double z_visu = -0.15;
void setup (FILE * fp)
{
  fprintf (fp,
	   "set pm3d map interpolate 2,2\n"
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set cbrange [%g:%g]\n"
	   "set xlabel 'x'\n"
	   "set ylabel 'height'\n"
	   "set xrange [%g:%g]\n"
	   "set yrange [%g:%g]\n"
     "set lmargin at screen 0.1\n"
     "set rmargin at screen 0.9\n", V - 0.1, V + 0.1, X0 + + L0/4., (X0+ 3*L0/4.), z_visu, -z_visu
	   );
}

void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "sp '-' u 1:2:3\n", t);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g\n", x, max(z_visu, z), u.x[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g\n", x, max(z_visu, min(z,-z_visu)), u.x[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);  
}

event gnuplot (t += T_END/100. ; t = 0.; t <= T_END)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  fprintf (fp,
	   "set term pngcairo font \",10\" size 2000,400\n"
	   "set output 'plot-%04d.png'\n", (int) (t*100./T_END));
  if (i == 0)
    setup (fp);
  plot (fp);
}

event moviemaker (t = end) {
  char text_system[100];
  sprintf (text_system, "for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
    "ppm2mp4 movie-V-%g-sig-%g.mp4", V, sig);
  system (text_system);
}
