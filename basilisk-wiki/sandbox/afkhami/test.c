/**
Layered Poiseuille flow with embedded boundaries. 

- Liquid 1 is located at the bottom between *y=-H/2* and *y=-H$_2$/2*, and
- Liquid 2 is located at the top of Liquid 1 between *y=-H$_2$/2* and *y=H$_2$/2*, and
- Liquid 3 is located at the top between *y=H$_2$/2* and *y=H/2*.


$\rho_1, \rho_2, \rho_3$, and
$\mu_1, \mu_2, \mu_3$ are the densities and viscosities of each fluid.
*/

#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

/**
Dynamic viscosity:
*/

#define MU1 1.0
#define MU2 1.0
#define MU3 1.0

#define FILTERED

/**
Interface between fluids is defined using Volume-of-Fluid, defined as
the interfaces 12 and 23. We define Y$_2$ as the vertical position of the center of the middle layer, and H$_2$ as its thickness.
*/

#include "vof.h"
scalar f12[], f23[], ff12[], ff23[];
scalar * interfaces = {f12,f23};
#define Y2 0.0
#define H2 0.1

/**
The properties of the fluid ($\tilde{\rho},\tilde{\mu},\dots$) are evaluated as functions of the volume fraction of each fluid.

  $$ \tilde{\rho} = (1-f_{12})(1-f_{23})\tilde{\rho}_1 + f_{12}~(1-f_{23})~\tilde{\rho}_2 + f_{12}~f_{23}~\tilde{\rho}_3 $$
  $$ \tilde{\mu} = (1-f_{12})(1-f_{23})\tilde{\mu}_1 + f_{12}~(1-f_{23})~\tilde{\mu}_2 + f_{12}~f_{23}~\tilde{\mu}_3 $$
  $$ \vdots $$
*/

#define FLUID(f12,f23,x1,x2,x3) (((1.-(clamp(f12,0,1)))*(1.-(clamp(f23,0,1)))*x1) + ((clamp(f12,0,1))*(1.-(clamp(f23,0,1)))*x2) + ((clamp(f12,0,1))*(clamp(f23,0,1))*x3))

#define MUTILDE(f12,f23) 	FLUID(f12,f23,MU1,MU2,MU3)

u.n[left] = neumann(0.);
u.n[right] = neumann(0.);

p[left] = dirichlet(128.);
p[right] = dirichlet(0.);

double EPS;

int main()
{

  size (1. [0]);
  DT = 0.001 [0];
  
  origin (0., -L0/2.);

  stokes = true;
  
  
  TOLERANCE = 1e-6;
  
  for (N = 32; N <= 256; N *= 2)
    run();
}

scalar un[];

event init (t = 0) {

  scalar phi12[];
  foreach_vertex()
    phi12[] = y - (Y2 - H2/2.);
  fractions (phi12, f12);

  scalar phi23[];
  foreach_vertex()
    phi23[] = y - (Y2 + H2/2.);
  fractions (phi23, f23);
  
  boundary ({f12,f23});

  mu = new face vector;
  
  /**
  The viscosity is unity. 
  Space and time should be dimensionless to be able to use the ‘mu = fm’ trick.
  */
   
  //mu = fm;
  
  /**
  The geometry of the embedded boundary is defined as two plates located at $y = \pm (L0 + 2\Delta)/8$. */
  
  EPS = L0/N;
  solid (cs, fs, intersection((L0/8. + EPS/4.) - y, (L0/8. + EPS/4.) + y));


  /**
  The boundary condition is zero velocity on the embedded boundaries. */
  
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
  
  /**
  We use "third-order" face flux interpolation. */

  for (scalar s in {u})
    s.third = true;
  
  /**
  We initialize the reference velocity. */
  
  foreach()
    un[] = u.y[];
}

event properties (i++) {

  face vector muv = mu;
  
#ifdef FILTERED //2D
  foreach() {
    ff12[] = (4.*f12[] + 
	      2.*(f12[0,1] + f12[0,-1] + f12[1,0] + f12[-1,0]) +
	      f12[-1,-1] + f12[1,-1] + f12[1,1] + f12[-1,1])/16.;
    ff23[] = (4.*f23[] + 
	      2.*(f23[0,1] + f23[0,-1] + f23[1,0] + f23[-1,0]) +
	      f23[-1,-1] + f23[1,-1] + f23[1,1] + f23[-1,1])/16.;
  }  
  boundary ({ff12,ff23});
#endif

  foreach_face() {
#ifdef FILTERED
    double f12m = (ff12[] + ff12[-1])/2.;
    double f23m = (ff23[] + ff23[-1])/2.;
#else
    double f12m = (f12[] + f12[-1])/2.;
    double f23m = (f23[] + f23[-1])/2.;
#endif   
    //if (MU1 || MU2 || MU3)
    muv.x[] = fm.x[]*MUTILDE(f12m,f23m); //fm field contains the product of the metric and solid factors.
  }
}

/**
We look for a stationary solution. */

event logfile (t += 0.01; i <= 1000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-12)
    return 1; /* stop */
}

/**
We compute error norms and display the horizontal velocity, pressure and error as a function of resolution. */

#define poiseuille(y) -128./2. * (y*y - pow(L0/8. + EPS/4., 2))

event profile (t = end)
{
  scalar e[];
  foreach() {
    //if (cs[] > 0.) 
      e[] = u.x[] - poiseuille (y);
    //
    //else
      //e[] = p[] = u.x[] = u.y[] = nodata;
  }

  norm n = normf (e);
  fprintf (stderr, "%d %.3g %.3g %.3g %d %d %d %d %d\n",
           N, n.avg, n.rms, n.max, i, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  
  dump();
  
  //draw_vof ("cs", "fs", filled = -1, fc = {1, 1, 1});
  //squares ("u.x", linear = true, spread = -1);
  draw_vof ("cs", "fs");
  squares ("u.x", spread = -1);
  save ("ux.png");

  draw_vof ("cs", "fs");
  squares ("p", spread = -1);
  save ("p.png");

  draw_vof ("cs", "fs");
  squares ("e", spread = -1);
  save ("e.png");
  
  char file[80];
  sprintf(file, "data_velo_%d", N);
  FILE * fp = fopen(file, "w");
  
  foreach() {
    //if (cs[] > 0.)
      fprintf (fp, "%g %g %g %g %g\n",
          y, u.x[], u.y[], p[], e[]);
  }
}

/**
## Results

![Velocity field](test/ux.png)
Velocity field
![Pressure field](test/p.png)
Pressure field
![Error field](test/e.png)
Error field

~~~gnuplot Velocity profile (N = 256)
L0 = 1.
N = 256.
set xlabel 'y'
set ylabel 'u_x'
plate_loc = (L0/8. + (L0/N/4.))
poiseuille(y) = - 128/2 * (y*y - plate_loc * plate_loc);
set grid
# set arrow from 0.25, graph 0 to 0.25, graph 1 nohead
# set arrow from 0.5, graph 0 to 0.5, graph 1 nohead
plot [-plate_loc:plate_loc][] 'data_velo_256' u 1:2 t 'numerics', poiseuille(x) t 'theory'
~~~


~~~gnuplot Error convergence
unset arrow
set xrange [*:*]
ftitle(a,b) = sprintf("%f/x^{%4.2f}", exp(a), -b)
f(x) = a + b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x) = a2 + b2*x
fit f2(x) '' u (log($1)):(log($2)) via a2,b2
set xlabel 'Resolution'
set logscale
set xtics 8,2,1024
set ytics format "% .0e"
set grid ytics
set cbrange [1:2]
set xrange [16:512]
set ylabel 'Error'
set yrange [*:*]
set key top right
plot '' u 1:4 pt 6 t 'max', exp(f(log(x))) t ftitle(a,b), \
     '' u 1:2 t 'avg', exp(f2(log(x))) t ftitle(a2,b2)
~~~

*/