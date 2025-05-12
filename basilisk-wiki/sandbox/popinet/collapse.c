/**
# Collapse of a bubble

A movie of the gradient of density is
[here](collapse/dy.mpg). Higher-resolution (*LEVEL = 11*) movies are
[here](collapse/level-11/dy.mpg). */

////////////////////////////////////////////////////////////

#include "all-mach.h"
#include "vof.h"

scalar f[], rhov[], rho1[], rho2[], * interfaces = {f};

face vector alphav[];

scalar rhoc2v[];

double CFLacous = 0.5, rhomin = 1e-6;

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  rhoc2 = rhoc2v;
  
  foreach()
    rho1[] = rho2[] = f[] = 1.;
  boundary ({rho1,rho2,f});
}

static scalar * interfaces1 = NULL;

event vof (i++) {
  vector q1 = q, q2[];
  foreach()
    foreach_dimension() {
      double u = q.x[]/rho[];
      q1.x[] = rho1[]*u;
      q2.x[] = rho2[]*u;
    }
  boundary ((scalar *){q1,q2});
  theta = 1.; // strict minmod limiter
  foreach_dimension() {
    q2.x.inverse = true;
    q1.x.gradient = q2.x.gradient = minmod2;
  }
  rho2.inverse = true;
  rho1.gradient = rho2.gradient = minmod2;
  f.tracers = {rho1,rho2,q1,q2};
  vof_advection ({f}, i);
  foreach() {
    foreach_dimension()
      q.x[] = q1.x[] + q2.x[];
    if (rho1[] < rhomin*f[])
      rho1[] = rhomin*f[];
    if (rho2[] < rhomin*(1. - f[]))
      rho2[] = rhomin*(1. - f[]);
  }
  boundary ((scalar *){q,rho1,rho2});
  interfaces1 = interfaces, interfaces = NULL;
}

event tracer_advection (i++) {
  interfaces = interfaces1;
}

double pstate (double rho, double rhoref, double pref,
	       double B, double gamma) {
  return pref*((1. + B)*pow(rho/rhoref, gamma) - B);
}

double cstate (double rho, double rhoref, double pref,
	       double B, double gamma) {
  return (1. + B)*pref/rhoref*gamma*pow(rho/rhoref, gamma - 1.);
}

event stability (i++) {
  foreach() {
    double dt = f[] > 0.5 ?
      Delta/sqrt(cstate (rho1[]/f[], 1., 1., 3000., 4.)) :
      Delta/sqrt(cstate (rho2[]/(1. - f[]), 0.001, 1., 0., 1.4));
    dt *= CFLacous;
    if (dt < dtmax)
      dtmax = dt;
  }
}

event properties (i++) {
  foreach() {
    rhov[] = rho1[] + rho2[];
    if (f[] > 0.5) {
      double r = rho1[]/f[];
      ps[] = pstate (r, 1., 1., 3000., 4.);
      rhoc2v[] = rhov[]*cstate (r, 1., 1., 3000., 4.);
    }
    else {
      double r = rho2[]/(1. - f[]);
      ps[] = pstate (r, 0.001, 1., 0., 1.4);
      rhoc2v[] = rhov[]*cstate (r, 0.001, 1., 0., 1.4);
    }
  }
  boundary ({rhov});
  
  foreach_face()
    alphav.x[] = 2./(rho[] + rho[-1]);
}

////////////////////////////////////////////////////////////

#include "tension.h"

#define LEVEL 7

int main() {
  X0 = Y0 = -L0/2.;
  f.sigma = 1.;
  run();
}

event init (i = 0) {
  if (!restore (file = "dump"))
    do {
      fraction (f, - (sq(0.25) - sq(x) - sq(y + 0.2)));
      foreach() {
	rho1[] = f[]*1.15;
	rho2[] = (1. - f[])/1000.;
	rhov[] = rho1[] + rho2[];
      }
      boundary ({rho1,rho2,rhov});
    } while (adapt_wavelet ({f,rho}, (double[]){5e-3,5e-4}, LEVEL).nf);
}

event logfile (i++) {
  stats s1 = statsf(rho1), s2 = statsf(rho2);
  scalar rhot[];
  foreach()
    rhot[] = rho1[] + rho2[];
  stats st = statsf(rhot);
  if (i == 0)
    fprintf (stderr,"t f.sum s1.min s1.sum s1.max "
	     "s2.min s2.sum s2.max st.min st.sum st.max\n");
  fprintf (stderr, "%g %.12f %.12f %.12f %.12f %.12f"
	   " %.12f %.12f %.12f %.12f %.12f\n",
	   t, statsf(f).sum,
	   s1.min, s1.sum, s1.max,
	   s2.min, s2.sum, s2.max,
	   st.min, st.sum, st.max);
}

#if 1
event gfsview (i += 10)
{
  static FILE * fp = popen ("gfsview2D -s collapse.gfv", "w");
  output_gfs (fp, t = t);
}
#endif

#if 0
event snapshot (i = 3184) {
  dump (file = "dump");
}

event snapshots (i++) {
  char name[80];
  sprintf (name, "snapshot-%d.gfs", i);
  output_gfs (file = name);
}
#endif

event movies (t += 5e-5; t <= 0.05)
{
  output_ppm (rho, linear = true, spread = 2);
  static FILE * fp1 = popen("ppm2mpeg > p.mpg", "w");
  output_ppm (p, fp1, linear = true, spread = 2);
  scalar dr[];
  foreach() {
    dr[] = 0.;
    foreach_dimension()
      dr[] += sq(rho[1] - rho[-1]);
    dr[] = sqrt(dr[])/(2.*Delta);
  }
  boundary({dr});
  static FILE * fp2 = popen("ppm2mpeg > dr.mpg", "w");
  output_ppm (dr, fp2, map = gray, min = 0, max = 2, n = 512, linear = true);

  foreach()
    dr[] = (rho[0,1] - rho[0,-1])/(2.*Delta);
  boundary({dr});
  static FILE * fp3 = popen("ppm2mpeg > dy.mpg", "w");
  output_ppm (dr, fp3, map = gray, min = -2, max = 2, n = 512, linear = true);

  foreach()
    dr[] = (rho[0,1] + rho[0,-1] - 2.*rho[])/sq(Delta);
  boundary({dr});
  static FILE * fp4 = popen("ppm2mpeg > d2y.mpg", "w");
  output_ppm (dr, fp4, map = gray, min = -10, max = 10, n = 512, linear = true);

  foreach()
    dr[] = level;
  boundary({dr});
  static FILE * fp5 = popen("ppm2mpeg > level.mpg", "w");
  output_ppm (dr, fp5, min = 0, max = LEVEL, n = 512);
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,rho}, (double[]){5e-3,5e-4}, LEVEL);
}
#endif

/**
## See also

* [Momentum-conserving incompressible two-phase flow](/src/momentum.h)
*/
