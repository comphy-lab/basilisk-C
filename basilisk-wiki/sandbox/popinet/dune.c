/**
# Dune in a shear flow

The shear flow over a dune is solved using the [hydrostatic
multilayer](/src/layered/hydro.h) approximation. 

The dune moves due to [erosion/deposition](erosion.h).

We compare the [steady](erosion.h) and [unsteady](erosion-t.h)
solvers. They seem to converge toward similar but not identical
solutions.

~~~gnuplot Dune evolution $\epsilon = 0.3, a = 1$
set term @SVG size 1024,200
set xlabel 'x'
set ylabel 'z'
set yrange [0:0.16]
set key below
plot 'out' u 1:3 index 0:10 w l t 'steady', \
     '../dune-unsteady/out' u 1:3 index 0:10 w l t 'unsteady'
~~~

With the same mass but half the initial height, the dune shape
converges toward the same solution.

~~~gnuplot Dune evolution $\epsilon = 0.3, a = 0.5$
plot 'out' u 1:3 index 11:21 w l t 'steady', \
     '../dune-unsteady/out' u 1:3 index 11:21 w l t 'unsteady'
~~~

A smaller dune which still grows (this is unexpected).

~~~gnuplot Dune evolution $\epsilon = 0.05, a = 0.25$
plot 'out' u 1:3 index 22:32 w l t 'steady', \
     '../dune-unsteady/out' u 1:3 index 22:32 w l t 'unsteady'
~~~

~~~gnuplot The corresponding final frictions.
reset
set xlabel 'x'
set ylabel 'dudz'
set key top left
plot 'out' u 1:5 index 10 w l t 'eps = 0.3, a = 1', \
     'out' u 1:5 index 21 w l t 'eps = 0.3, a = 0.5', \
     'out' u 1:5 index 32 w l t 'eps = 0.05, a = 0.25'
~~~

## References

~~~bib
@article{kouakou2006,
  title={Evolution of a model dune in a shear flow},
  author={Kouakou, Kouam{\'e} Kan Jacques and Lagr{\'e}e, Pierre-Yves},
  journal={European Journal of Mechanics-B/Fluids},
  volume={25},
  number={3},
  pages={348--359},
  year={2006},
  publisher={Elsevier},
  pdf = {http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/kouakoulagree06.pdf}
}
~~~
*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/implicit.h"
#include "layered/remap.h"
#include "layered/perfs.h"

/**
We add erosion/deposition to the bathymetry. */

#if UNSTEADY
# include "erosion-t.h"
#else
# include "erosion.h"
#endif

/**
These are the dune size parameters. */

double eps = 0.3, a = 1.;

int main()
{
  periodic (right);
  L0 = 2. [1];
  X0 = - 0.4;
  G = 10;
  nu = 0.1;

  /**
  We drive the flow with $\partial_zu = 1$ on the free-surface. */

  const vector du0[] = {1};
  dut = du0;

  /**
  Numerical parameters: */
  
  N = 256;
  nl = 30;
  theta_H = 1.; // more free-surface damping
  
  /**
  These are the parameters for the dune material. */

  tau_s = 0.9 [0,-1];
#if !UNSTEADY 
  l_s = 1./64., E = 4e-3 [2];
#else

  /**
  For the [unsteady solver](erosion-t.h), we set equivalent parameters
  (with $l_a$ the additional "advection length" parameter). */
  
  double l_s = 1./64., E = 4e-3 [2];
  l_a = 8./64., l_e = E/l_s, B = l_a/l_s;

#endif
  
  /**
  The non-erodible "bedrock" is at $z_b = 0$. */
  
  z_br = zeroc;

  /**
  We run several cases to illustrate the different dynamics. */

  eps = 0.3, a = 1.;
  run();

  /**
  Same mass, but half the initial amplitude. */

  eps = 0.3, a = 0.5;
  run();

  /**
  Half the amplitude and half the mass. */

  eps = 0.05, a = 0.25;
  run();
}

/**
We initialise the topography and the initial thickness of each layer
*h*. */

event init (i = 0)
{
  fprintf (stderr, "# eps = %g, a = %g\n", eps, a);
  printf ("# eps = %g, a = %g\n", eps, a);

  double r = 0.06 [1], b = 1. [0,-1];
  foreach() {
    zb[] = 0.5*a*eps*exp(- sq(a*x/r));
    double z = 0.;
    foreach_layer() {
      h[] = (1. - zb[])/nl;
      z += h[]/2.;
      u.x[] = b*z;
      z += h[]/2.;
    }
  }
  boundary (all);
}

/**
Uncomment this part if you want on-the-fly animation. */

#if 1
event gnuplot (i += 10)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term x11\n"
	     // "set size ratio -1\n"
	     );
  fprintf (fp,
           "set title 'nl=%d, t=%f'\n"
           "set xl 'x'\nset yl 'h'\n"
           "plot "
	   //	   "'-' u 1:2 w l t 'eta'," 
	   "'-' u 1:3 w l t 'zb',"
	   "'-' u 1:4 w l t 'dz'\n",
	   nl, t); 
  foreach()
    fprintf (fp, "%g %g %g %g\n", x, eta[], zb[], dz[]);
  fprintf (fp, "e\n");
  fflush (fp);
  //  usleep (10000);
}
#endif

/**
## Outputs

We output the profiles at regular intervals. */

event profiles (t += 5; t <= 50)
{
  foreach (serial)
    printf ("%g %g %g %g %g\n",
	    x, eta[], zb[], dz[], dudz(u));
  printf ("\n\n");
  fflush (stdout);
}

event logfile (i += 10)
{
  fprintf (stderr, "%g %g %g\n", t, dt, statsf(zb).sum);
}
