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

We solve here the Navier-Stokes equations and add the NACA2414 using
an [embedded boundary](/src/embed.h). */

#include "grid/quadtree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define chord   (1.) // NACA2414 chord length
#define p_theta (6.*M_PI/180.) // Incidence of the NACA2414 airfoil
#define Re      (10000.) // Reynolds number based on the cord length
#define uref    (1.) // Reference velocity, uref
#define tref    ((chord)/(uref)) // Reference time, tref

/**
We also define the shape of the domain. */

#define naca00xx(x,y,a) (sq (y) - sq (5.*(a)*(0.2969*sqrt   ((x))	\
					      - 0.1260*((x))		\
					      - 0.3516*sq   ((x))	\
					      + 0.2843*cube ((x))	\
					      - 0.1036*pow  ((x), 4.)))) // -0.1015 or -0.1036

void p_shape (scalar c, face vector f, coord p)
{
  // NACA2414 parameters
  double mm = 0.02;
  double pp = 0.4;
  double tt = 0.14;

  // Rotation parameters around the position p,
  // located at the position cc in the airfoil referential
  double theta = (p_theta);
  coord cc = {0.25*(chord), 0.};
  
  vertex scalar phi[];
  foreach_vertex() {

    // Coordinates with respect to the center of rotation of the airfoil p
    // where the head of the airfoil is identified as xrot = 0, yrot = 0
    double xrot = cc.x + (x - p.x)*cos (theta) - (y - p.y)*sin (theta);
    double yrot = cc.y + (x - p.x)*sin (theta) + (y - p.y)*cos (theta);

    if (xrot >= 0. && xrot <= (chord)) {

      // Camber line coordinates, adimensional
      double xc = xrot/(chord), yc = yrot/(chord), thetac = 0.;
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
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

#define lmin (6) // Min mesh refinement level (l=6 is 4pt/c)
#define lmax (12) // Max mesh refinement level (l=12 is 256pt/c)
#define cmax (1.e-3*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
  /**
  The domain is $16\times 16$. */

  L0 = 16.;
  size (L0);
  origin (-L0/2., -L0/2.);

  /**
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-4;
  TOLERANCE_MU = 1.e-4*(uref);  

  /**
  We initialize the grid. */
  
  N = 1 << (lmin);
  init_grid (N);

  run ();
}

/**
## Boundary conditions 

We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[bottom] = 0;
uf.n[top]    = 0;

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = (uref)*(chord)/(Re)*fm.x[];
  boundary ((scalar *) {muv});
}

/**
## Initial conditions */

event init (i = 0)
{
  /**
  We set the viscosity field in the event *properties*. */

  mu = muv;
  
  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */
    
#if ORDER2
  for (scalar s in {u, p, pf})
    s.third = false;
#else
  for (scalar s in {u, p, pf})
    s.third = true;
#endif // ORDER2

  /**
  We use a slope-limiter to reduce the errors made in small-cells. */

#if SLOPELIMITER
  for (scalar s in {u, p, pf}) {
    s.gradient = minmod2;
  }
#endif // SLOPELIMITER
  
#if TREE
  /**
  When using *TREE* and in the presence of embedded boundaries, we
  should also define the gradient of *u* at the cell center of
  cut-cells. */
#endif // TREE
  
  /**
  We initialize the embedded boundary. */
  
#if TREE
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */
  
  astats ss;
  int ic = 0;
  do {
    ic++;
    p_shape (cs, fs, p_p);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
			maxlevel = (lmax), minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // TREE
  
  p_shape (cs, fs, p_p);

  /**
  We initialize the particle's speed and accelerating. */
  
  p_u.x  = -(uref);
}

/**
## Embedded boundaries 

The particle's position is advanced to time $t + \Delta t$. */

event advection_term (i++)
{
  p_p.x += (p_u.x)*(dt);
}

/**
## Adaptive mesh refinement */

#if TREE
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
  		 maxlevel = (lmax), minlevel = (1));

  /**
  We do not need here to reset the embedded fractions to avoid
  interpolation errors on the geometry as the is already done when
  moving the embedded boundaries. It might be necessary to do this
  however if surface forces are computed around the embedded
  boundaries. */
}
#endif // TREE

/**
## Outputs */

event logfile (i++; t <= (0.5)*(tref))
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq ((uref))*(chord));
  double CL = (Fp.y + Fmu.y)/(0.5*sq ((uref))*(chord));

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g\n",
	   i, (t)/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   p_p.x/(chord),
	   CD, CL);
  fflush (stderr);
}

/**
## Animation */

event movie (i += 25)
{
  scalar omega[];
  vorticity (u, omega);

  view (fov = 1, camera = "front",
	tx = -0.5/(L0), ty = 1.e-12,
	bg = {1,1,1},
	width = 800, height = 400);
  
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("omega", map = cool_warm, min = -50, max = 50);
  save ("vorticity.mp4");
}

/**
## Results

![Vorticity](naca2414-starting/vorticity.mp4)
*/
