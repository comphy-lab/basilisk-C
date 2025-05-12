/**
# Wingtip vortex of a 3D NACA2414 wing at $Re=1000$

This test case is inspired from the Gerris test case
[wingtip](http://gerris.dalembert.upmc.fr/gerris/examples/examples/wingtip.html). Wingtip
vortices are tubes of circulating air which are left behind a wing as
it generates lift. One wingtip vortex trails from the tip of each
wing. The cores of vortices spin at very high speed and are regions of
very low pressure. Vortex-lines cannot end in fluid, therefore in
three dimensions, the starting vortex shed by the trailing edge on
take-off remains joined to the wing by counterrotating wingtip
vortices, which are associated with the escape of air from the
underside to above via the wingtips, following the negative pressure
difference from below to above which must be present if the wing is
generating lift. This vortex-line consisting of the starting vortex
and pair of wingtip vortices is notionally continued and closed by the
‘bound vortex’ inside the wing, associated with the circulation around
the wing, and proportional to its lift in accordance with the
Kutta–Joukowsky theorem. As the flow over the wing evolves from
start-up to steady-state, the starting vortex recedes and the wing-tip
vortices extend back to infinity.

We solve here the Navier-Stokes equations and add the NACA2414 using
an [embedded boundary](/src/embed.h). */

#include "grid/octree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving.h"
#include "../myperfs.h"
#include "view.h"
#include "lambda2.h"

/**
## Reference solution */

#define chord   (1.) // NACA2414 chord length
#define wing    (2.) // Length of the wing
#define p_theta (6.*M_PI/180.) // Incidence of the NACA2414 airfoil
#define Re      (1000.) // Reynolds number based on the cord length
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
  coord cc = {0.25*(chord), 0., 0.};
  
  vertex scalar phi[];
  foreach_vertex() {

    // Coordinates with respect to the center of rotation of the airfoil p
    // where the head of the airfoil is identified as xrot = 0, yrot = 0
    double xrot = cc.x + (x - p.x)*cos (theta) - (y - p.y)*sin (theta);
    double yrot = cc.y + (x - p.x)*sin (theta) + (y - p.y)*cos (theta);
    double zrot =        (z - p.z);

    if (xrot >= 0. && xrot <= (chord) && fabs (zrot) <= (wing)) {
      
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

#define lmin (5) // Min mesh refinement level (l=5 is 4pt/c)
#define lmax (9) // Max mesh refinement level (l=9 is 64pt/c)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field

int main ()
{  
  /**
  The domain is $8\times 8\times 8$. */

  L0 = 8.;
  size (L0);
  origin (-L0/2., -L0/2., -L0/2.);

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
uf.n[back]   = 0;
uf.n[front]  = 0;

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
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax),(cmax)},
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

  double CD =  (Fp.x + Fmu.x)/(0.5*sq ((uref))*pi*sq ((chord)/2.));
  double CLy = (Fp.y + Fmu.y)/(0.5*sq ((uref))*pi*sq ((chord)/2.));
  double CLz = (Fp.z + Fmu.z)/(0.5*sq ((uref))*pi*sq ((chord)/2.));

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g\n",
	   i, (t)/(tref), dt/(tref),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   p_p.x/(chord),
	   CD, CLy, CLz);
  fflush (stderr);
}

/**
## Animation */

event movie (i += 25)
{
  scalar l2[];
  lambda2 (u, l2);

  view (fov = 13.8568,
	quat = {-0.204453,0.521493,0.0481286,0.826999},
	tx = 0.0823552, ty = -0.0249205, bg = {0.3,0.4,0.6},
	width = 600, height = 600, samples = 1);
  
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.4);
  cells (n = {0,0,1}, alpha = 1.e-12);
  save ("l2.mp4", opt = "-r 10");
}

/**
## Results

![Isosurface $\lambda_2 = -0.4$](naca2414-starting3D/l2.mp4)
*/
