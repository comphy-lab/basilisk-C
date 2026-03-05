/**
# Two-phase Couette flow

Consider a two-phase Couette flow in a channel of height 2H with the origin in the center at $y = 0$, which is periodic in the $x$-direction. The top and bottom boundary has velocity $U_1$ and $U_2$, respecitvely. 

The flow profile is a piecewise-linear where the slopes are related by the viscosity ratio $\chi = \mu_1/\mu_2$. 

The vertical flow profile $u_x(y)$ can be analytically computed with no-slip boundary conditions are the interface ($y=0$) and at the walls and a continuity of shear stress at the interface $\mu_1\partial_yu_1(0) = \mu_2\partial_yu_2(0)$, 

$$u_1(y) = \frac{(U_1 - U_2)}{H(1+\chi)}y + \frac{\chi U_1 + U_2}{1+\chi}, \hspace{3cm}u_2(y) = \frac{\chi(U_1 - U_2)}{H(1+\chi)}y + \frac{\chi U_1 + U_2}{1+\chi},$$

where 1 denotes the upper phase and 2 the lower phase. .*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"

double U1 = 1.;
double U2 = 0.;
double chi = 1e1;

u.t[top] = dirichlet(U1);
u.n[top] = dirichlet(0.);

u.t[bottom] = dirichlet(U2);
u.n[bottom] = dirichlet(0.);
int main()
{

  size (1. [0]);
  DT = 1. [0];
  
  origin (-L0/2., -L0/2.);

  periodic(right);
  
  stokes = true;

  mu1 = 1e0; rho1 = 1e0;
  mu2 = mu1/chi; rho2 = 1e0;

  TOLERANCE = 1e-5;
  for (int lvl = 4; lvl <= 8; lvl++){
    N = pow(2, lvl);
    run();
  }
}

scalar un[];

#define WIDTH 0.5

event init (t = 0) {
  fraction(f, y);
}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 100000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-7)
    return 1; /* stop */
}

/**
We compute error norms and display the horizontal velocity, pressure and
error fields using bview. */


#define couette_low(y, chi, U1, U2)  chi*(U1 - U2)/(L0/2. * (1 + chi)) * y + (chi*U1 + U2)/(1 + chi)
#define couette_high(y, chi, U1, U2)  (U1 - U2)/(L0/2. * (1 + chi)) * y + (chi*U1 + U2)/(1 + chi)

event profile (t = end)
{
  scalar e[];
  foreach() {
    if (y > 0. && (point.i == N/2)) {
      // When looking at the log, e[] is not the difference.
      e[] = u.x[] - couette_high(y, chi, U1, U2);
      fprintf(stderr, "high %g %g %g %g\n", y, u.x[], couette_high(y, chi, U1, U2), e[]);
    }
    if (y < 0. && (point.i == N/2)){
      e[] = u.x[] - couette_low(y, chi, U1, U2);
      fprintf(stderr, "low %g %g %g %g\n", y, u.x[], couette_low(y, chi, U1, U2), e[]);
    }

  }

  norm n = normf (e);
  fprintf (stderr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
	   N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  
    char file[80];
    sprintf(file, "data_couette_%d", N);
    FILE * fp = fopen(file, "w");
    foreach() {
      if (point.i == N/2)
        fprintf (fp, "%g %g %g %g %g\n", y, u.x[], u.y[], couette_low(y, chi, U1, U2), couette_high(y, chi, U1, U2));
    }
}

/**
## Results

# The results seem good but a discrepancy can be seen near the interface. It could be that my analytical solution is not fully correct (Grid cell off). The interface should in principle be located at y = 0 exactly. 


~~~gnuplot Horizontal velocity
set xlabel 'y'
set ylabel 'u(x)'
set xrange [-0.5:0.5]
set yrange [0:1]
U1 = 1.
U2 = 0. 
chi = 1e1; 
L0 = 1.;
couette_low(y) = chi*(U1 - U2)/(L0/2. * (1 + chi)) * y + (chi*U1 + U2)/(1 + chi)
couette_high(y) = (U1 - U2)/(L0/2. * (1 + chi)) * y + (chi*U1 + U2)/(1 + chi)
set key bottom right
p 'data_couette_32' u 1:2 title 'N = 32', \
  'data_couette_64' u 1:2 title 'N = 64', \
  'data_couette_128' u 1:2 title 'N = 128', \
  'data_couette_256' u 1:2 title 'N =256', \
  couette_low(x) dt 2 lc rgb 'black' lw 2 title 'Analytical bottom', \
  couette_high(x) dt 3 lc rgb 'black' lw 2 title 'Analytical top'
~~~


~~~gnuplot Horizontal velocity (Zoom)
reset
set xlabel 'y'
set ylabel 'u(x)'
set xrange [-0.05:0.05]
U1 = 1.
U2 = 0. 
chi = 1e1; 
L0 = 1.;
couette_low(y) = chi*(U1 - U2)/(L0/2. * (1 + chi)) * y + (chi*U1 + U2)/(1 + chi)
couette_high(y) = (U1 - U2)/(L0/2. * (1 + chi)) * y + (chi*U1 + U2)/(1 + chi)
set key bottom right
p \
  'data_couette_64' u 1:2 ps 2 pt 3 title 'N = 64', \
  'data_couette_256' u 1:2 ps 2 pt 1 title 'N =256', \
  couette_low(x) dt 2 lc rgb 'black' lw 2 title 'Analytical bottom', \
  couette_high(x) dt 3 lc rgb 'black' lw 2 title 'Analytical top'
~~~

~~~gnuplot Error convergence 
reset
set xlabel 'y'
set ylabel 'abs((usim - uanal)/uanal)'
set logscale y
set format y "10^{%L}"

U1 = 1.
U2 = 0. 
chi = 1e1; 
L0 = 1.;

couette_low(y) = chi*(U1 - U2)/(L0/2. * (1 + chi)) * y + (chi*U1 + U2)/(1 + chi)
couette_high(y) = (U1 - U2)/(L0/2. * (1 + chi)) * y + (chi*U1 + U2)/(1 + chi)

p \
  'data_couette_32'  u ($1 > 0. ? $1 : 1./0.):(abs(($2 - couette_high($1 > 0. ? $1 : 1./0.))/couette_high($1 > 0. ? $1 : 1./0.)))  lc rgb 'blue'  notitle, \
  'data_couette_32'  u ($1 < 0. ? $1 : 1./0.):(abs(($2 - couette_low($1 < 0. ? $1 : 1./0.))/couette_low($1 < 0. ? $1 : 1./0.)))    lc rgb 'blue'  title 'N = 32', \
  'data_couette_64'  u ($1 > 0. ? $1 : 1./0.):(abs(($2 - couette_high($1 > 0. ? $1 : 1./0.))/couette_high($1 > 0. ? $1 : 1./0.)))  lc rgb 'green' notitle, \
  'data_couette_64'  u ($1 < 0. ? $1 : 1./0.):(abs(($2 - couette_low($1 < 0. ? $1 : 1./0.))/couette_low($1 < 0. ? $1 : 1./0.)))    lc rgb 'green' title 'N = 64', \
  'data_couette_128' u ($1 > 0. ? $1 : 1./0.):(abs(($2 - couette_high($1 > 0. ? $1 : 1./0.))/couette_high($1 > 0. ? $1 : 1./0.))) lc rgb 'black'  notitle, \
  'data_couette_128' u ($1 < 0. ? $1 : 1./0.):(abs(($2 - couette_low($1 < 0. ? $1 : 1./0.))/couette_low($1 < 0. ? $1 : 1./0.)))   lc rgb 'black'  title 'N = 128', \
  'data_couette_256' u ($1 > 0. ? $1 : 1./0.):(abs(($2 - couette_high($1 > 0. ? $1 : 1./0.))/couette_high($1 > 0. ? $1 : 1./0.))) lc rgb 'black'  notitle, \
  'data_couette_256' u ($1 < 0. ? $1 : 1./0.):(abs(($2 - couette_low($1 < 0. ? $1 : 1./0.))/couette_low($1 < 0. ? $1 : 1./0.)))   lc rgb 'black'  title 'N = 256'
~~~
*/

