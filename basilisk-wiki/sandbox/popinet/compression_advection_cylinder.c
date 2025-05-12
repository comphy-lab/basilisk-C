/**
# Problem with the advection of a cylinder with the all-mach compressible solver 
*/

#include "grid/multigrid.h"
#include "all-mach.h"

////////////////////////////////////////////////////////////

#include "vof.h"
#include "tension.h"

scalar c[], rhov[], rho1[], rho2[], * interfaces = {c}, * interfaces1 = {c};
face vector alphav[];

scalar rhoc2v[];

double c0 = 20;
double r = 0.001; // density ratio
double uadv = 1.; 

double cstate (double gamma) {
  return c0*gamma;
}

event vof (i++) {
  vector q1 = q, q2[];
  foreach()
    foreach_dimension() {
      double u = q.x[]/rho[];
      q1.x[] = rho1[]*u;
      q2.x[] = rho2[]*u;
    }   
  boundary ((scalar *){q1,q2});
  theta = 1.; 
  foreach_dimension() {
    q2.x.inverse = true;
    q1.x.gradient = q2.x.gradient = minmod2;
  }
  rho2.inverse = true;
  rho1.gradient = rho2.gradient = minmod2;
  c.tracers = {rho1,rho2,q1,q2};
  vof_advection ({c}, i); 
  foreach()
    foreach_dimension()
      q.x[] = q1.x[] + q2.x[];
  boundary ((scalar *){q});
  interfaces = NULL;
}

event acceleration (i++) {
  interfaces = interfaces1;
}


double pstate (double rho, double rhoref, double pref,
         double B, double gamma) {
  return pref + cstate(gamma)*(rho-rhoref);
}

event properties (i++) {
  alpha = alphav;
  rho = rhov;
  rhoc2 = rhoc2v;
  
  foreach() {
    rhov[] = rho1[] + rho2[];
    if (c[] > 0.5) {
      ps[] = pstate (rho1[]/c[], 1., 1., 0., 1.4);
      rhoc2v[] = rhov[]*cstate (1.4);
      //      if (1. - c[] < 1e-6)
      //        rho2[] = 0.;
    }
    else {
      ps[] = pstate (rho2[]/(1. - c[]), r, 1., 0., 1.4);
      rhoc2v[] = rhov[]*cstate (1.8);
      //      if (c[] < 1e-6)
      //        rho1[] = 0.;
    }
  }
  boundary ({rhov,rho1,rho2});
  
  foreach_face()
    alphav.x[] = 2./(rho[] + rho[-1]);
}

int main() {
  X0 = Y0 = -L0/2.;
  c.sigma = 0.;
  periodic (right); 
  N = 32;
  TOLERANCE = 1e-3;
  CFL = 0.5;
  run();
}

event defaults (i = 0) {
  foreach()
    rho1[] = rho2[] = c[] = 1.;
  boundary ({rho1,rho2,c});
}

event init (i = 0) {
  fraction (c, - (sq(0.25) - sq(x) - sq(y)));
  foreach() {
    rho1[] = c[];
    rho2[] = (1. - c[])*r;
    q.x[] = uadv*(rho1[] + rho2[]);
    q.y[] = 0.;
    p[] = 1.;
  }
  boundary ({rho1,rho2,q});
}

event logfile (i++; t <= 2.) {
  scalar errorp[];
  foreach()
    errorp[] = fabs(p[] - 1.);
  stats st = statsf(errorp);
  
  static FILE * fp    = fopen ("err.dat", "w");

  if (i == 0) 
    fprintf (fp,"#t c.sum st.min st.avg st.max "
       "st.stddev \n");
  
  fprintf (fp, "%g %.12f %.12f %.12f %.12f %.12f\n",
     t, dtnext, st.min, st.sum/st.volume, st.max,
     st.stddev);
  fprintf (stderr, "%g %.12f %.12f %.12f %.12f %.12f\n",
     t, 0., st.min, st.sum/st.volume, st.max,
     st.stddev);
}



/**
~~~gnuplot Temporal evolution of the error
set output 'error.png'
set xlabel 'Time'
set ylabel 'Error'
p "./err.dat" u 1:5 not w l
~~~
*/
