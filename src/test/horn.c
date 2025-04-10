/**
# Internal solitary waves

This test case is taken from section 5.6 of [Vitousek & Fringer,
2014](#vitousek2014) itself based on the laboratory experiments of
[Horn et al., 2001](#horn2001). In the experiment, two layers with
different densities separated by a (diffusive) interface are initially
tilted relative to the horizontal and released. A train of solitary
waves develops and bounces back and forth between the boundaries of
the tank. The amplitude of these waves is varied by varying the
initial slope.

The figures below give the measured and modelled wave amplitudes as a
function of time measured in the middle of the tank.

Note that both $\Delta\rho$ and $\nu$ are tuned differently from
Vitousek & Fringer, and that the results are very sensitive to this
tuning.

~~~gnuplot Model (Vitousek & Fringer, 2014) and experimental interface displacements (Figure 9 of [Vitousek & Fringer, 2014](#vitousek2014)).
set term svg enhanced size 625,800 font ",10"
set multiplot layout 5,1 scale 1,1.1
set xrange [0:400]
set yrange [-5:5]
set key top left
plot 'horn.horn' u 1:($2*100.) index 'ai = 0.01305' w l t 'Horn 2001', \
     'horn.vf' u 1:($2*100.) index 'ai = 0.01305' w l t 'V and F 2014'
unset key
plot 'horn.horn' u 1:($2*100.) index 'ai = 0.0261' w l t 'Horn 2001', \
     'horn.vf' u 1:($2*100.) index 'ai = 0.0261' w l t 'V and F 2014'
plot 'horn.horn' u 1:($2*100.) index 'ai = 0.03915' w l t 'Horn 2001', \
     'horn.vf' u 1:($2*100.) index 'ai = 0.03915' w l t 'V and F 2014'
plot 'horn.horn' u 1:($2*100.) index 'ai = 0.0522' w l t 'Horn 2001', \
     'horn.vf' u 1:($2*100.) index 'ai = 0.0522' w l t 'V and F 2014'
plot 'horn.horn' u 1:($2*100.) index 'ai = 0.0783' w l t 'Horn 2001', \
     'horn.vf' u 1:($2*100.) index 'ai = 0.0783' w l t 'V and F 2014'
unset multiplot
~~~

The results given by the multilayer solver are significantly better
than those of Vitousek & Fringer. This is mostly due to the improved
dispersion characteristics of the Keller box scheme.

~~~gnuplot Model (multilayer) and experimental interface displacement.
set term svg enhanced size 625,800 font ",10"
set multiplot layout 5,1 scale 1,1.1
set xrange [0:400]
set yrange [-5:5]
set key top left
plot 'horn.horn' u 1:($2*100.) index 'ai = 0.01305' w l t 'Horn 2001', \
     'log' u 1:($2*100.) index 'ai = 0.01305' w l t 'basilisk'
unset key
plot 'horn.horn' u 1:($2*100.) index 'ai = 0.0261' w l t 'Horn 2001', \
     'log' u 1:($2*100.) index 'ai = 0.0261' w l t 'basilisk'
plot 'horn.horn' u 1:($2*100.) index 'ai = 0.03915' w l t 'Horn 2001', \
     'log' u 1:($2*100.) index 'ai = 0.03915' w l t 'basilisk'
plot 'horn.horn' u 1:($2*100.) index 'ai = 0.0522' w l t 'Horn 2001', \
     'log' u 1:($2*100.) index 'ai = 0.0522' w l t 'basilisk'
plot 'horn.horn' u 1:($2*100.) index 'ai = 0.0783' w l t 'Horn 2001', \
     'log' u 1:($2*100.) index 'ai = 0.0783' w l t 'basilisk'
unset multiplot
~~~

## References

~~~bib
@article{vitousek2014,
  title={A nonhydrostatic, isopycnal-coordinate ocean model for internal waves},
  author={Vitousek, Sean and Fringer, Oliver B},
  journal={Ocean Modelling},
  volume={83},
  pages={118--144},
  year={2014},
  publisher={Elsevier},
  doi={10.1016/j.ocemod.2014.08.008}
}

@article{horn2001, 
  title={The degeneration of large-scale interfacial gravity waves in lakes}, 
  volume={434}, 
  doi={10.1017/S0022112001003536}, 
  journal={Journal of Fluid Mechanics}, 
  publisher={Cambridge University Press}, 
  author={Horn, D. A. and Imberger, J. and Ivey, G. N.}, 
  year={2001}, 
  pages={181–207}
}
~~~
*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/isopycnal.h"

/**
We add an optional sanity check for the free surface evolution. */

#include "layered/check_eta.h"

/**
The parameters are viscosity, height of the layer interfaces and
initial amplitude. The relative density of each layer is given by
`drho`. */

double nu_H = 2e-5;

#define H0 0.29
#define hi 0.087
double ai;

double * drho = (double []){ 0.018, 0 };

int main()
{
  size (6.);
  origin (- L0/2.);
  N = 256;
  nl = 2;
  G = 9.81;
  nu = nu_H;
  DT = 0.05;
#if 0
  lambda_b[] = {HUGE}; // free slip
#endif

  /**
  We do several initial amplitudes. */
  
  ai = 1.305e-2;
  run();
  ai = 2.61e-2;
  run();
  ai = 3.915e-2;
  run();
  ai = 5.22e-2;
  run();
  ai = 7.83e-2;
  run();
}

/**
The initial layer positions. */

event init (i = 0)
{
  foreach() {
    double zi = ai*x/(L0/2.);
    foreach_layer()
      h[] = point.l < nl/2 ? (hi - zi)/(nl/2) : (H0 - hi + zi)/(nl/2);
  }
}

/**
This adds (an approximation of) horizontal viscosity. */

event viscous_term (i++) {
#if NH
  horizontal_diffusion ({u, w}, nu_H, dt);
#else
  horizontal_diffusion ((scalar *){u}, nu_H, dt);
#endif
}

/**
We log the position of the layer in the middle of the tank. */

event logfile (i += 3; t <= 400)
{
  if (i == 0)
    fprintf (stderr, "\n\n# ai = %g\n", ai);
  double z = 0.;
  foreach_layer()
    if (_layer < nl/2)
      z += interpolate (h, 0.);
  fprintf (stderr, "%g %g %g %d\n", t, z - hi, dt,
#if NH	   
	   mgp.i
#else
	   0
#endif
	   );
}

/**
We can optionally display the evolving layers and velocity field. */

#if 0
void plot_setup (FILE * fp)
{
  fprintf (fp,
	   "set pm3d map corners2color c2\n"
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set xlabel 'x'\n"
	   "set ylabel 'depth'\n"
	   );
}

void output_fields (FILE * fp)
{
  foreach (serial) {
    double z = zb[];
#if !NH
    scalar w = u.x;
#endif
    fprintf (fp, "%g %g %g %g\n", x, z, u.x[], w[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g %g\n", x, z, u.x[], w[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");  
}

void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.1f s'\n"
	   "sp '-' u 1:2:3\n",
	   t);
  output_fields (fp);
  //  fprintf (fp, "pause .01\n");

  fflush (fp);  
}

event gnuplot (i += 100)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0) {    
    fprintf (fp, "set term x11\n");
    plot_setup (fp);
  }
  plot (fp);
  //  fprintf (fp, "pause .1\n");
}
#endif

/**
This optionally displays consistency between `res_eta` and `deta`
(corresponding to the `check_eta.h` option above). */

#if 0 // show consistency between res_eta and deta
#include "deta.h"
#endif
