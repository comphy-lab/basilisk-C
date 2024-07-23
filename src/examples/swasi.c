/**
# A Shallow Water Analogue for the Standing Accretion Shock Instability

We use the Saint-Venant solver to reproduce the experimental setup of
[Foglizzo et al. 2012](#references). We use radial coordinates and
Basilisk View to visualise the results. */

#include "grid/multigrid.h"
#include "radial.h"
#include "saint-venant.h"
#if dimension > 1
#  include "view.h"
#endif

/**
We use seven levels of refinement (i.e. $128^2$ grid points) by
default. The geometrical parameters and other inputs are those of the
experimental setup. Note that we chose centimetres and seconds as
length and time units. */

int LEVEL = 7;

const double Ri = 4. [1]; // inner radius (cm)
const double Ro = 32.;    // outer radius (cm)
const double Ho = 0.074;  // injection layer thickness (cm)
const double Q  = 1e3;    // flow rate (cm^3/s)
const double R_45 = 5.6;  // hyperbolic surface scaling (cm)

/**
We double the equivalent viscosity compared to that given in Foglizzo
et al. 2012 (0.03 cm^2/s), to obtain a more stable spiral mode. The
"tube height" is set (by trial and error using 1D runs) so that the
shock radius obtained numerically is close to 20 cm. */

const double NU = 0.06;   // equivalent viscosity (cm^2/s)
const double STEP = 6.98; // height of exit tube (cm)

/**
The injection velocity is the flow rate divided by the surface area of
the outer injection ring. */

const double Uo = Q/(2.*pi*Ro*Ho);

/**
The maximum runtime (seconds). */

const double TMAX = 150.;

/**
## Boundary conditions

The outer injection ring corresponds to the right boundary of the
radial coordinate system. 

Note that this code can also be run in one (radial) dimension by
changing the first line to "grid/multigrid1D.h". */

u.n[right] = - Uo;
#if dimension > 1
u.t[right] = dirichlet(0.);
#endif
h[right] = Ho;

/**
The inner outflow is the left boundary. */

u.n[left] = dirichlet(- Ro*Uo*Ho/(STEP*Ri));

/**
## Main program

The level of refinement can be given as a command-line argument. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  dtheta = 2.*pi;
  G = 981.; // acceleration of gravity (cm.s^-2)
  size (Ro - Ri);
  origin (Ri, 0.);
#if dimension > 1
  periodic (top);
#endif
  init_grid (1 << LEVEL);
  run();
}

/**
## Initial conditions

We setup the hyperbolic bottom profile and set an initial uniform
fluid layer and radial velocity identical to the injection
conditions. */

event init (i = 0)
{
  foreach() {
    zb[] = - sq(R_45)/r;
    h[] = Ho;
    u.x[] = - Uo;
  }
  zb[left]  = dirichlet(- sq(R_45)/r);
  zb[right] = dirichlet(- sq(R_45)/r);
}

/**
## Perturbation

At $t=40$ seconds the shock is close to being stationary, we perturb
the velocity field to trigger the instability. 

We add a systematic mode 1 perturbation of the radial
component of the velocity, tapering off at $R_i$ and
$R_o$. This gives a clean growth of the mode and allows a
simple estimate of the frequency and growth rate (see below). 

Optionally we can also add a small, constant azimuthal perturbation
which will trigger the "spiral" instability (see second movie
below). */

event perturb (t = 40) {
  foreach() {
    u.x[] += 1e-2*Uo*((r - Ri)/(Ro - Ri))*(1. - (r  - Ri)/(Ro - Ri))*sin(theta);
#if 0
    u.y[] += 1e-2*Uo;
#endif
  }
}

/**
## Friction

We use a time-implicit discretisation of the linear friction. */

event friction (i++) {
  if (NU > 0.)
    foreach() {
      double a = h[] < dry ? HUGE : 1. + NU*dt/sq(h[]);
      foreach_dimension()
	u.x[] /= a;
    }
}

/**
## Outputs

We log the evolution of the fluid height and norms of the azimuthal
(in 2D) or radial (in 1D) component of the velocity as well as the
coordinates of the center-of-mass of the fluid. */

event logfile (i += 10; t <= TMAX) {
  stats s = statsf (h);
#if dimension > 1  
  norm n = normf (u.y);
#else
  norm n = normf (u.x);
#endif
  double sumx = 0., sumy = 0.;
  foreach (reduction(+:sumx) reduction(+:sumy))
    sumx += r*cos(theta)*h[]*dv(),
    sumy += r*sin(theta)*h[]*dv();
  fprintf (stderr, "%g %d %g %g %g %g %g %g %g\n", t, i, s.min, s.max, s.sum,
	   n.rms, n.max, sumx/s.sum, sumy/s.sum);
}

/**
The evolution with time of the y-coordinate of the center of mass is
illustrated below. It clearly shows the exponential growth of the
oscillation of this position followed by non-linear saturation (see
also the animation below).

~~~gnuplot Evolution of the y-coordinate of the center-of-mass of the liquid
set xlabel 'Time (sec)'
set ylabel 'Position (cm)'
plot [0:150]'log' u 1:9 w l t ''
~~~

This can be used to estimate the period of oscillation (2.74 s) and
the initial growth rate (0.07 $s^{-1}$) as illustrated below

~~~gnuplot Same as above but using a logscale
set ylabel 'Absolute value of position (cm)'
set logscale y
set key top left
set grid
plot [40:150][1e-3:10]'log' u 1:(abs($9)) w l t 'numerical', \
     '' u 1:(exp(0.07*$1)/1050*abs(sin(2.*pi*$1/2.74+0.55))) w l t 'exp(0.07*x)*|sin(2.*pi*x/2.74)|'
~~~

While we run, we display the evolution of a vertical cross-section of
the flow. */

event profiles (i += 100) {
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp,
	     "set term x11 noraise\n"
	     "set grid\n"
	     "set xrange [%g:%g]\n", X0, X0 + L0);
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "p '-' u 1:3:2 w filledcu lc 3 t ''\n", t);
  foreach (serial)
#if dimension > 1
    if (y < Delta)
#endif
      fprintf (fp, "%g %g %g\n", x, eta[], zb[]);
  fprintf (fp, "e\n\n");
  fflush (fp);
}

/**
We save the stationary profile at $t=50$. */

event stationary (t = 50) {
  FILE * fp = fopen ("stationary", "w");
  foreach (serial)
#if dimension > 1
    if (y < Delta)
#endif
      fprintf (fp, "%g %g %g %g\n", x, eta[], zb[], u.x[]);
  fclose (fp);
}

/**
~~~gnuplot Stationary profile at $t=50$
reset
set xlabel 'radius (cm)'
set grid
p 'stationary' u 1:3:2 w filledcu lc 3 t '', \
  '' u 1:(-$4/20.) w l lt 1 lw 2 t '-u_r/20'
~~~

## Animation

We create an animation of the free surface height $\eta$. We first
define a coordinate mapping function which will be used by Basilisk
View. */

#if dimension > 1
void radial (coord * p) {
  double r = p->x, theta = p->y*dtheta/L0;
  p->x = r*cos(theta), p->y = r*sin(theta);
}

event movie (t = 60; t += 0.1) {
  view (map = radial, fov = 45, width = 600, height = 600, samples = 1);
  clear();
  squares ("eta", min = -1.4, max = -0.2, linear = true);
  save ("eta.mp4");
}
#endif

/**
![SWASI instability. Growth and non-linear saturation of mode 1.](swasi/eta.mp4)

![SWASI instability. Spiral mode produced by adding a small azimuthal
perturbation.](swasi/eta-spiral.mp4)

## References

~~~bib
@article{foglizzo2012,
  title={Shallow water analogue of the standing accretion shock instability: 
  Experimental demonstration and a two-dimensional model},
  author={Foglizzo, Thierry and Masset, Fr{\'e}d{\'e}ric and Guilet, 
  J{\'e}r{\^o}me and Durand, Gilles},
  journal={Physical Review Letters},
  volume={108},
  number={5},
  pages={051103},
  year={2012},
  publisher={APS},
  url = {https://arxiv.org/pdf/1112.3448}
}
~~~

* [Videos](http://irfu.cea.fr/Projets/SN2NS/outreach.html)
* [PRL supplementary material](http://irfu.cea.fr/Projets/SN2NS/PRL_Supp_Mat.pdf)
*/
