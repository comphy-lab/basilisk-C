/**
# Starting vortex of a 2D NACA2414 airfoil at $Re=10000$

This test case is inspired from the Gerris test case
[starting](http://gerris.dalembert.upmc.fr/gerris/examples/examples/starting.html). A
starting vortex is a vortex which forms in the air adjacent to the
trailing edge of an airfoil as it is accelerated from rest in a
fluid. It leaves the airfoil and remains (nearly) stationary in the
flow. This phenomenon does not depend for its existence on the
viscosity of the air, or the third dimension. However, it does depend
on the trailing edge being sharp, as it is for most aerofoils.

![Animation of the vorticity field. The colorscale corresponds to
$\pm$ 20 U/L.](naca2414-starting/vorticity.mp4)(autoplay width="100%")

We solve the Navier-Stokes equations and add the NACA2414 airfoil using
an [embedded boundary](/src/embed.h). */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "navier-stokes/perfs.h"
#include "view.h"

/**
## Parameters */

const double chord = 1.;            // NACA2414 chord length
const double p_theta = 6.*pi/180.;  // Incidence of the NACA2414 airfoil
const double Re = 10000.;           // Reynolds number based on the cord length
const double uref = 1.;             // Reference velocity, uref
const double tref = chord/uref;     // Reference time, tref

/**
The NACA profile. */

double naca00xx (double x, double y, double a) {
  return sq(y) - sq(5.*(a)*(0.2969*sqrt(x)
			    - 0.1260*(x)
			    - 0.3516*sq(x)
			    + 0.2843*cube(x)
			    - 0.1036*pow(x, 4.))); // -0.1015 or -0.1036
}

void airfoil (scalar c, face vector f, coord p = {0})
{
  // NACA2414 parameters
  double mm = 0.02;
  double pp = 0.4;
  double tt = 0.14;

  // Rotation parameters around the position p,
  // located at the position cc in the airfoil referential
  double theta = p_theta;
  coord cc = {0.25*chord, 0.};
  
  vertex scalar phi[];
  foreach_vertex() {

    // Coordinates with respect to the center of rotation of the airfoil p
    // where the head of the airfoil is identified as xrot = 0, yrot = 0
    double xrot = cc.x + (x - p.x)*cos (theta) - (y - p.y)*sin (theta);
    double yrot = cc.y + (x - p.x)*sin (theta) + (y - p.y)*cos (theta);

    if (xrot >= 0. && xrot <= chord) {

      // Camber line coordinates, adimensional
      double xc = xrot/chord, yc = yrot/chord, thetac = 0.;
      if (xc < pp) {
	yc     -= mm/sq (pp)*(2.*pp*xc - sq (xc));
	thetac = atan (2.*mm/sq (pp)*(pp - xc));
      }
      else {
	yc     -= mm/sq (1. - pp)*(1. - 2.*pp + 2.*pp*xc - sq (xc));
	thetac = atan (2.*mm/sq (1. - pp)*(pp - xc));
      }
      
      // Thickness
      phi[] = (naca00xx (xc, yc, tt*cos (thetac)));
    }
    else
      phi[] = 1.;
  }
  fractions (phi, c, f);
  fractions_cleanup (c, f);
}

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

const int lmax = 12; // Max mesh refinement level (l=12 is 256pt/c)
const double cmax = 1e-3*uref; // Absolute refinement criteria for the velocity field

int main()
{  
  /**
  The domain is $16\times 16$, centered on the origin. */

  size (16 [1]);
  origin (-L0/2., -L0/2.);

  mu = muv;
  run();
}

/**
## Boundary conditions 

Inflow on the left and outflow on the right. */

u.n[left] = dirichlet (uref);

u.n[right]  = neumann(0.);
p[right]    = dirichlet(0.);
pf[right]   = dirichlet(0.);

/**
## Properties

The viscosity is set to match the desired Reynolds number. */

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*uref*chord/Re;
}

/**
## Initial conditions */

event init (i = 0)
{
  
  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */
    
  for (scalar s in {u, p, pf})
    s.third = true;

  /**
  We initialize the embedded boundary. */
  
#if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */
  
  astats ss;
  int ic = 100;
  do {
    airfoil (cs, fs);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30}, maxlevel = lmax);
  } while ((ss.nf || ss.nc) && --ic);
#endif // TREE
  
  airfoil (cs, fs);

  /**
  The initial velocity in the fluid is uniform. Note that this does
  not respect the incompressibility condition (due to the presence of
  the wing) and that the first few timesteps are necessary to obtain a
  divergence-free velocity field which matches the boundary
  conditons. */
  
  foreach()
    u.x[] = cs[]*uref;
  
  /**
  The boundary condition is zero velocity on the embedded boundary. */
    
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);  
}

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2, cmax, cmax}, maxlevel = lmax);
}
#endif // TREE

/**
## Outputs */

event logfile (i++)
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq(uref)*chord);
  double CL = (Fp.y + Fmu.y)/(0.5*sq(uref)*chord);

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g\n",
	   i, (t)/tref, dt/tref,
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   CD, CL);
  fflush (stderr);
}

/**
## Animation */

scalar omega[];

void picture()
{
  vorticity (u, omega);
  view (fov = 30, near = 0.01, far = 1000,
	tx = -0.043, ty = 0.004, tz = -0.087,
	width = 1600, height = 500);
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("omega", map = blue_white_red, min = -20, max = 20, linear = true);
}

event movie (i = 25; i += 10)
{
  picture();
  save ("vorticity.mp4");
}

event ending (t = tref)
{
  picture();
  save ("vorticity.png");
}
