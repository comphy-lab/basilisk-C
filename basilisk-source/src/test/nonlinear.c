/**
# Non-linear geostrophic adjustment

For a circular vortex defined by a tangential velocity $u_\theta(r)$, the
radial height/pressure profile is a solution of
$$
h'_0(r) = \frac{u_\theta}{g}\left(2\Omega + \frac{u_\theta}{r}\right)
$$
with $\Omega$ the angular velocity. For this test case we take
$$
u_\theta(r) = (r < 0.4)\epsilon(1 + \cos(\pi(r - 0.2)/0.2))/2
$$
The control parameters are the Froude number $U/\sqrt{gH}$ and the Rossby
number $Ro=U/\Omega L$. We set the Froude number to 0.1 and consider Rossby
numbers 0.1 and $\infty$ (no rotation). In the case without rotation the
errors reflect only the accuracy of the momentum advection terms. With
rotation, the errors also depend on the accuracy of the discretisation
of the geostrophic balance.

The figures below illustrate the evolution of the errors on free
surface height/pressure for the different solvers, with and without
rotation.

~~~gnuplot Evolution of the maximum relative error on free-surface height. $Ro = \infty$.
set xlabel 'Time'
set ylabel 'Maximum relative error'
set logscale y
plot 'log' index 'F0 = 0' u 1:2 w l t 'C grid', \
     '../nonlinear-ml/log' index 'F0 = 0' u 1:2 w l t 'multilayer'
~~~

~~~gnuplot Evolution of the maximum relative error on free-surface height. $Ro = 0.1$.
plot 'log' index 'F0 = 0.1' u 1:2 w l t 'C grid', \
     '../nonlinear-ml/log' index 'F0 = 0.1' u 1:2 w l t 'multilayer'
~~~

The total energy is exactly conserved for the C-grid scheme and very
well conserved for the multilayer scheme.

~~~gnuplot Evolution of the total energy
set ylabel 'Normalised total energy'
unset logscale
set key bottom left
plot 'log' index 'F0 = 0.1' u 1:3 w l t 'C grid (Ro = 0.1)', \
     '../nonlinear-ml/log' index 'F0 = 0.1' u 1:($3+$4) w l  \
        t 'multilayer (Ro = 0.1)',			     \
     'log' index 'F0 = 0' u 1:3 w l t 'C grid (Ro = infty)', \
     '../nonlinear-ml/log' index 'F0 = 0' u 1:($3+$4) w l    \
        t 'multilayer (Ro = infty)'
~~~

## See also

* [Same test with
  Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/nonlinear.html)
*/

#include "grid/multigrid.h"
#include <gsl/gsl_integration.h>
#pragma autolink -lgsl -lgslcblas
#if ML
# include "layered/hydro.h"
# include "layered/implicit.h"
double F0 = 0.;
# define F0() F0
# include "layered/coriolis.h"
#else
# include "atmosphere.h"
#endif

int main (int argc, char * argv[])
{
  // coordinates of lower-left corner
  origin (-0.5 + 1e-12, -0.5 + 1e-12);
  // number of grid points
  init_grid (argc > 1 ? atoi(argv[1]) : 64);
  // size of the box
  size (1.[0]);
  DT = HUGE[0];
  // acceleration of gravity
  G = 1.;
  // Coriolis parameter
  F0 = 0.1;
  // CFL number: the C-grid model is unstable for larger CFL
  CFL = 0.25;
#if !ML
  if (N > 256)
    NU = 1e-4; // need some viscosity to stabilize at high resolutions
#endif
  for (F0 = 0.; F0 <= 0.1; F0 += 0.1) {
    fprintf (stderr, "# F0 = %g\n", F0);
    run();
    fprintf (stderr, "\n\n");
  }
}

/* ---------------- Initial conditions ------------------- */

#define H0 1.
#define FROUDE 0.1

double vtheta (double r) {
  return FROUDE*(r < 0.4)*(1. + cos((r - 0.2)/0.2*M_PI))/2.;
}

double h0p (double r, void * p) {
  double vt = vtheta(r);
  return vt*(F0 + vt/r)/G;
}

double h0 (double r) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double result, error;
  gsl_function F;
  F.function = &h0p;
  gsl_integration_qags (&F, 0, r, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
  return result;
}

scalar h1[];

event init (i = 0)
{
#if ML
  CFL_H = 0.25;
#endif
  foreach (cpu) {
    zb[] = - H0;
    h1[] = h[] = (H0 + h0(sqrt (x*x + y*y)));
  }
#if ML
  foreach() {
    double r = sqrt (x*x + y*y), vt = vtheta(r);
    u.x[] = - vt*y/r;
    u.y[] =   vt*x/r;
  }
#else
  foreach_face(x) {
    double r = sqrt (x*x + y*y);
    u.x[] = - vtheta(r)*y/r;
  }
  foreach_face(y) {
    double r = sqrt (x*x + y*y);
    u.y[] =   vtheta(r)*x/r;
  }
#endif
}

/* ------------------ Output helper functions --------------- */

scalar e[];

double error()
{
  double max = 0.;
  foreach(reduction(max:max)) {
    e[] = fabs (h1[]  - h[]);
    if (e[] > max) max = e[];
  }
  return max/(normf(h1).max - H0);
}

typedef struct {
  double ke, pe;
} Energy;

Energy energy()
{
  double KE = 0., PE = 0.;
  foreach(reduction(+:KE) reduction(+:PE)) {
#if ML
    double ke = (sq(u.x[]) + sq(u.y[]))/2.;
#else
    double ke = (sq(u.x[] + u.x[1,0]) + sq(u.y[] + u.y[0,1]))/8.;
#endif
    KE += h[]*ke*dv();
    PE += G*sq(zb[] + h[])/2.*dv();
  }
  return (Energy){ KE, PE };
}

event logfile (i += 10; t <= 5.)
{
  static Energy E0 = {0., 0.};
  Energy E = energy();
  if (i == 0)
    E0 = E;
  fprintf (stderr, "%5g %5g %5g %5g\n", t, error(),
	   E.ke/(E0.ke + E0.pe), E.pe/(E0.ke + E0.pe));
}

event plots (t = end)
{
  output_ppm (e, file = "e.png", spread = -1, n = 256);
}
