/**
# Charge relaxation in a planar cross-section

This is the same problem as in [cyl_axi.c]() but in a planar
cross-section of the column. */

#define NAVIER 1
int LEVEL;

#include "grid/multigrid.h"
#include "ehd/implicit.h"
#if NAVIER
# include "navier-stokes/centered.h"
# include "ehd/stress.h"
#endif
#include "fractions.h"

/**
Far away from the conducting column the electric potential is zero. */

phi[bottom] = dirichlet(0);
phi[top]    = dirichlet(0);
phi[right]  = dirichlet(0);
phi[left]   = dirichlet(0);
#if NAVIER
p[top]      = dirichlet(0);
u.n[top]    = neumann(0);
p[bottom]   = dirichlet(0);
u.n[bottom] = neumann(0);
p[left]     = dirichlet(0);
u.n[left]   = neumann(0);
p[right]    = dirichlet(0);
u.n[right]  = neumann(0);
#endif

#define beta 3.
#define cond 3.
#define rhoini 0.5
#define R 0.1
#define circle(x,y) (sq(R) - sq(x) - sq(y))

scalar f[];

event init (t = 0) {
  face vector s[];
  solid (f, s, circle(x, y));
  foreach()
    rhoe[] = rhoini*f[];

  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    epsilon.x[] = (ff*beta + (1. - ff))*fm.x[];
    K.x[] = cond*s.x[]*fm.x[];
  }
}

event chargesum (i++) {
  double Q = statsf(rhoe).sum;
  static double Q0;
  if (i == 0)
    Q0 = Q;
  else
    assert (fabs(Q - Q0) < 1e-7);
}

event epfield (t = 10) {
  char name[80];
  sprintf (name, "Er-%d", LEVEL);
  FILE * fp = fopen (name, "w");

  scalar ee[];
  foreach() {
    double Ex = (phi[-1,0] - phi[1,0])/(2*Delta);
    double Ey = (phi[0,-1] - phi[0,1])/(2*Delta);
    double r = sqrt(x*x + y*y);
    double En = sqrt(Ex*Ex + Ey*Ey);
    ee[] = En - (r < R ? 0. : 0.5*sq(R)*rhoini/r);
#if NAVIER
    fprintf (fp, "%g %g %g %g\n", r, En, p[], rhoe[]);
#else
    fprintf (fp, "%g %g %g\n", r, En, rhoe[]);
#endif
  }
  fclose (fp);
  
  norm n = normf (ee);
  fprintf (stderr, "%d %g %g %g\n", LEVEL, n.avg, n.rms, n.max);
}

#if TREE // fixme: this does not work
event adapt (i++) {
  double Qb = statsf(rhoe).sum;
  adapt_wavelet ({f}, (double []){1e-6},
		 maxlevel = LEVEL + 1, minlevel = LEVEL);
  double Qa = statsf(rhoe).sum;
  assert (fabs(Qa - Qb) < 1e-10);
}
#endif

int main() {

  /**
  The computational domain spans [-1:1][-1:1]. */

  X0 = Y0 = -1.;
  L0 = 2.;
  DT = 1;
  TOLERANCE = 1e-7;

  /**
  We compute the solution for different levels of refinement. */
  
  for (LEVEL = 6; LEVEL <= 8; LEVEL++) {
    N = 1 << LEVEL;
    run();
  }
}

/**
## Results

~~~gnuplot Radial electric field distribution as a function of grid refinement.
set term PNG enhanced font ",10"
set output 're.png'
set xlabel 'r'
set ylabel 'E_r'
set xrange[0:1]
set yrange[0.:0.03]
set sample 1000
R = 0.1
rhoini = 0.5
E(x) = x < R ? 0 : (R*R*rhoini/2/x)
plot 'Er-6' u 1:2 t "Level = 6",  \
     'Er-7' u 1:2 t "Level = 7",  \
     'Er-8' u 1:2 t "Level = 8",  \
     E(x) w l t "Analytical"
~~~

~~~gnuplot Pressure distribution as a function of grid refinement.
set output 'p.png'
set xlabel 'r'
set ylabel 'p'
set xrange[0:0.4]
set yrange[-0.0005:0.00005]
set sample 1000
set key right bottom
R = 0.1
rhoini = 0.5
p(x) = x > R ? 0 : -(R*R*rhoini*rhoini/8)
plot 'Er-6' u 1:3 t "Level = 6",  \
     'Er-7' u 1:3 t "Level = 7",  \
     'Er-8' u 1:3 t "Level = 8",  \
      p(x) w l t "Analytical"
~~~

~~~gnuplot Convergence of error on the norm of the electric field
reset
set term @SVG
R = 0.1
set xlabel '2R/Delta'
set ylabel 'Error norms on |E|'
set logscale
set xtics 3.2,2,102.4
set xrange [5:30]
plot 'log' u (R*2**$1):4 w lp t 'Max', \
     'log' u (R*2**$1):3 w lp t 'RMS', \
     .025/x t 'first-order'
~~~

## See also

* [Charge relaxation in an axisymmetric insulated conducting column](cyl_axi.c)
* [Same test with 
   Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/cylinder.html#planar)
*/
