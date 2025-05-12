/**
# Vortex shedding behind a sphere at Reynolds = 300

![Animation of the $\lambda_2$ vortices coloured with the vorticity
 component aligned with the flow.](sphere/movie.mp4)(loop)

We solve the Navier--Stokes equations on an adaptive octree and use
embedded boundaries to define the sphere. */

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "view.h"

/**
We will use the $\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) for vortex detection. */

#include "lambda2.h"

/**
This is the maximum level of refinement i.e. an equivalent maximum
resolution of $256^3$. */

int maxlevel = 8;

/**
We need a new field to define the viscosity. */

face vector muv[];

/**
The domain size is $16^3$. We move the origin so that the center of
the unit sphere is not too close to boundaries. */

double D = 1. [1];

int main()
{
  init_grid (64);
  size (16.*D);
  origin (- 3.*D, -L0/2., -L0/2.);
  mu = muv;
  run();
}


/**
The viscosity is given by the Reynolds number, the sphere diameter and
the inflow velocity */

double U0 = 1., Re = 300.;

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*D*U0/Re;
}

/**
The boundary conditions are inflow with unit velocity on the
left-hand-side and outflow on the right-hand-side. */

u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The boundary condition is no slip on the embedded boundary. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);

event init (t = 0) {

  /**
  We initially refine only in a sphere, slightly larger than the solid
  sphere. */

  refine (x*x + y*y + z*z < sq(1.2*D/2.) && level < maxlevel);

  /**
  We define the unit sphere. */

  solid (cs, fs, x*x + y*y + z*z - sq(D/2.));

  /**
  We set the initially horizontal velocity to the inflow velocity
  everywhere (outside the sphere). */
  
  foreach()
    u.x[] = cs[] ? U0 : 0.;
}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We use Basilisk view to create the animated isosurface of $\lambda_2$
for $30 <= t <= 60$. */

event movies (t = 30; t += 0.25; t <= 60)
{

  /**
  Here we compute two new fields, $\lambda_2$ and the vorticity
  component in the $y-z$ plane. */
  
  scalar l2[], vyz[];
  foreach()
    vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);
  lambda2 (u, l2);

  view (fov = 11.44, quat = {0.072072,0.245086,0.303106,0.918076},
	tx = -0.307321, ty = 0.22653, bg = {1,1,1},
	width = 802, height = 634);
  draw_vof ("cs", "fs");
  isosurface ("l2", -0.01, color = "vyz", min = -1, max = 1,
	      linear = true, map = cool_warm);
  save ("movie.mp4");
}

/**
We set an adaptation criterion with an error threshold of 0.02 on all
velocity components and $10^{-2}$ on the geometry. */

event adapt (i++) {
  astats s = adapt_wavelet ({cs,u}, (double[]){1e-2,0.02,0.02,0.02}, // fixme: these are not seen....
			    maxlevel, 4);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
