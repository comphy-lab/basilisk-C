/**
# Two Phase Poiseuille flow

See twophase_couette.c for details. 
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "view.h"

double dpdx = -1;
double muv = 1;
double H = 0.1;

double EPS;

double chi = 1e1;
u.n[left] = neumann(0.);
u.n[right] = neumann(0.);
/**
We should set a consistent pressure drop
 */

p[left] = dirichlet(-dpdx*L0);
p[right] = dirichlet(0.); // This is fine
int lvl;
int main()
{
  
  size (1. [0]);
  origin (-L0/2., -L0/2.);

    mu1 = 1e0; rho1 = 1e0;
  mu2 = mu1/chi; rho2 = 1e0;
  stokes = true;
  // quadratic interpolation for quadratic profile
  u.x.third = true;
  TOLERANCE = 1e-5;
  for ( lvl = 7; lvl <= 9; lvl++){ // lvl = 7 can be problematic due to small fluid-fraction cells
    N = pow(2, lvl);
    run();
  }
}

scalar un[];

event init (t = 0) {
  EPS = L0/111.;
  solid (cs, fs, intersection((H + EPS) - y, (H + EPS) + y));
  multigrid_restriction ({cs, fs});

  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);

  fraction(f, y);

}


/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-8)
    return 1; /* stop */
}

event profile (t = end)
{ 
    char file[80];
    sprintf(file, "data_poiseuille_%d", N);
    FILE * fp = fopen(file, "w");
    foreach() {
      if (point.i == N/2)
        fprintf (fp, "%g %g %g\n", y, u.x[], u.y[]);
    }
}

/**
## Results

The error is worse near the interface. This only means that the continuity of shear is not respected (numerically it is I suppose).

~~~gnuplot Horizontal velocity
set xlabel 'y'
set ylabel 'u(x)'
set xrange [-0.12:0.12]
set yrange [0:0.022]
h = (0.1 + 1./111.) 
dpdx = 1.
mu1 = 1.
mu2 = 0.1
chi = mu1/mu2
beta = (1./mu1 - 1./mu2)
A1 = -dpdx*h/(2. * mu1) * (mu1 - mu2)/(mu2 + mu1)
A2 =  -dpdx*h/(2. * mu2) * (mu1 - mu2)/(mu2 + mu1)
B1 = dpdx/(2.)*h**2 * (2./(mu1 + mu2))
B2 = B1
set key top right
pois_high(y) = -dpdx/(2.*mu1)*y**2 + A1*y + B1
pois_low(y) = -dpdx/(2.*mu2)*y**2 + A2*y + B2
p \
  'data_poiseuille_128' u 1:2 title 'N = 128', \
  'data_poiseuille_256' u 1:2 title 'N = 256', \
  'data_poiseuille_512' u 1:2 title 'N = 512', \
  "<echo '0 -0.02'" w impulse dt 2 lc rgb 'black' title 'interface',\
  pois_high(x) dt 3 lc rgb 'black' lw 2 title 'uanal top', \
  pois_low(x) dt 2 lc rgb 'black' lw 2 title 'uanal bottom'
~~~

~~~gnuplot Horizontal velocity (Zoom)
reset
set xlabel 'y'
set ylabel 'u(x)'
set xrange [-0.02:0.02]
set yrange [0.008:0.019]
h = (0.1 + 1./111.)
dpdx = 1.
mu1 = 1.
mu2 = 0.1
chi = mu1/mu2
beta = (1./mu1 - 1./mu2)
A1 = -dpdx*h/(2. * mu1) * (mu1 - mu2)/(mu2 + mu1)
A2 =  -dpdx*h/(2. * mu2) * (mu1 - mu2)/(mu2 + mu1)
B1 = dpdx/(2.)*h**2 * (2./(mu1 + mu2))
B2 = B1
set key top right
pois_high(y) = -dpdx/(2.*mu1)*y**2 + A1*y + B1
pois_low(y) = -dpdx/(2.*mu2)*y**2 + A2*y + B2
p \
  'data_poiseuille_128'  u 1:2 w l title 'N = 128', \
  'data_poiseuille_256'  u 1:2 w l title 'N = 256', \
  'data_poiseuille_512'  u 1:2 w l title 'N = 512', \
  "<echo '0 0.1'" w impulse dt 2 lc rgb 'black' title 'interface',\
  pois_high(x > 0 ? x : 1/0) dt 3 lc rgb 'black' lw 2 title 'uanal top', \
  pois_low(x < 0 ? x : 1/0) dt 2 lc rgb 'black' lw 2 title 'uanal bottom'
~~~

~~~gnuplot Error convergence 
reset 
set xlabel 'y'
set ylabel 'u(x)'
h = (0.1 + 1./111.) 
set xrange [-h:h]
dpdx = 1.
mu1 = 1.
mu2 = 0.1
chi = mu1/mu2
beta = (1./mu1 - 1./mu2)
A1 = -dpdx*h/(2. * mu1) * (mu1 - mu2)/(mu2 + mu1)
A2 =  -dpdx*h/(2. * mu2) * (mu1 - mu2)/(mu2 + mu1)
B1 = dpdx/(2.)*h**2 * (2./(mu1 + mu2))
B2 = B1
set key bottom left
pois_high(y) = -dpdx/(2.*mu1)*y**2 + A1*y + B1
pois_low(y) = -dpdx/(2.*mu2)*y**2 + A2*y + B2
set logscale y
set format y "10^{%L}"
p \
  'data_poiseuille_128'  u ($1 > 0. ? $1 : 1./0.):(abs(($2 - pois_high($1 > 0. ? $1 : 1./0.))/pois_high($1 > 0. ? $1 : 1./0.)))  lc rgb 'blue'  notitle, \
  'data_poiseuille_128'  u ($1 < 0. ? $1 : 1./0.):(abs(($2 - pois_low($1 < 0. ? $1 : 1./0.))/pois_low($1 < 0. ? $1 : 1./0.)))    lc rgb 'blue'  title 'N = 128, In channel: approx 25', \
  'data_poiseuille_256'  u ($1 > 0. ? $1 : 1./0.):(abs(($2 - pois_high($1 > 0. ? $1 : 1./0.))/pois_high($1 > 0. ? $1 : 1./0.)))  lc rgb 'green' notitle, \
  'data_poiseuille_256'  u ($1 < 0. ? $1 : 1./0.):(abs(($2 - pois_low($1 < 0. ? $1 : 1./0.))/pois_low($1 < 0. ? $1 : 1./0.)))    lc rgb 'green' title 'N = 256, In channel: approx 50', \
  'data_poiseuille_512' u ($1 > 0. ? $1 : 1./0.):(abs(($2 - pois_high($1 > 0. ? $1 : 1./0.))/pois_high($1 > 0. ? $1 : 1./0.))) lc rgb 'black'  notitle, \
  'data_poiseuille_512' u ($1 < 0. ? $1 : 1./0.):(abs(($2 - pois_low($1 < 0. ? $1 : 1./0.))/pois_low($1 < 0. ? $1 : 1./0.)))   lc rgb 'black'  title 'N = 512, In channel: approx 100'
~~~
*/

