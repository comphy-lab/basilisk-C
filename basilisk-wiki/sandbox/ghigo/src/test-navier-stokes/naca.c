/**

# 2D NACA2414 airfoil at $Re=1,000$

This test case is inspired from the Gerris test case
[starting](http://gerris.dalembert.upmc.fr/gerris/examples/examples/starting.html). 

We solve here the Navier-Stokes equations and add the NACA2414 using
an [embedded boundary](/src/embed.h). */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

double Reynolds = 1000.;
int maxlevel = 8;
face vector muv[];



/**
## Reference solution */

#define chord   (1.) // NACA2414 chord length
#define uref    (1.) // Reference velocity, uref
#define tref    ((chord)/(uref)) // Reference time, tref

/**
We also define the shape of the domain. */

#define naca00xx(x,y,a) (sq (y) - sq (5.*(a)*(0.2969*sqrt   ((x))	\
					      - 0.1260*((x))		\
					      - 0.3516*sq   ((x))	\
					      + 0.2843*cube ((x))	\
					      - 0.1036*pow  ((x), 4.)))) // -0.1015 or -0.1036


/**
The domain is four units long, centered vertically. */

int main() {
  L0 = 4;
  origin (-0.5, -L0/2.);
  N = 256;
  mu = muv;
  run(); 
}



event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]/Reynolds;
}

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The airfoil is no-slip. */

u.n[embed] = dirichlet(0.); 
u.t[embed] = dirichlet(0.);


double naca(double x, double y)
{
  // NACA2414 parameters
  double mm = 0.02;
  double pp = 0.4;
  double tt = 0.14;

if (x >= 0. && x <= (chord)) {
  // Camber line coordinates, adimensional
  double xc = x/(chord), yc = y/(chord), thetac = 0.;
  if (xc < pp) {
    yc     -= mm/sq (pp)*(2.*pp*xc - sq (xc));
    thetac = atan (2.*mm/sq (pp)*(pp - xc));
  }
  else {
    yc     -= mm/sq (1. - pp)*(1. - 2.*pp + 2.*pp*xc - sq (xc));
    thetac = atan (2.*mm/sq (1. - pp)*(pp - xc));
  }
 return naca00xx (xc, yc, tt*cos (thetac));
 }
 else
   return 1.;
}

event init (t = 0)
{
  solid (cs, fs, naca(x,y));

  /**
  We set the initial velocity field. */
  
  foreach()
    u.x[] =  1.;

  scalar omega[];
  vorticity (u, omega);
  view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.211, ty = 0.053, tz = -2.508,
      width = 818, height = 768);
  box();
  squares (color = "omega");
  draw_vof(c="cs",lw=2, lc = {1,1,1});
  cells();
  save("omega-init.png");
}

event logfile (i += 10){
  scalar omega[];
  vorticity (u,omega);
  foreach()
    {
      if(cs[]< 1.0)
	{
	  omega[]=nodata;
	}
      else
	{
	  omega[]=fabs(omega[]);
	}
    }
  stats om  = statsf(omega);
  fprintf (stderr, "%d %g %g %g\n", i, t, om.max, om.stddev);
}

/**
We produce animations of the vorticity field... */

event movies (t += 0.1; t <= 4.)
{
  scalar omega[];
  vorticity (u, omega);
  view (quat = {0.000, 0.000, 0.000, 1.000},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.211, ty = 0.053, tz = -2.508,
      width = 818, height = 768);
  box();
  squares (color = "omega");
  draw_vof(c="cs",lw=2, lc = {1,1,1});
  cells();
  save("omega.mp4");
}

/**
We adapt according to the error on the embedded geometry and velocity fields. */

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-3,1e-3,1e-3}, maxlevel, 4);
}

/**
## Results

![Vorticity](naca/omega.mp4)
*/