/**
# Lid-driven cavity of a Bingham fluid at Re=0 

The double projection is necessary here for hight yield stresses. */

#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "view.h"

/**
## Mesh

We define here the level of refinement. */

#define lmax (7)

/**
## Yield stress */

double T1 = 0.;

/**
We need a face-centered field for the viscosity. */

face vector muv[];

/**
## Code */

int main()
{
  /**
  The domain is $1\times 1$. */

  L0 = 1.;
  size (L0);
  origin (0., 0.);

  /**
  We set the numerical parameters. */

  DT = 1.e-3;
  
  N = 1 << (lmax);
  init_grid (N);

  /**
  We tune the Poisson solver. */

  stokes = true;
  TOLERANCE = 1.e-5;
  
  run ();
}

/**
## Boundary conditions */

//y = L0 (slip)
u.t[top] = dirichlet(1);

//others (no-slip)
u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);

/**
We also make sure there is no flow through the walls, otherwise the
compatibility condition for the Poisson equation can be violated. */

uf.n[bottom] = dirichlet (0.);
uf.n[top]    = dirichlet (0.);
uf.n[left]   = dirichlet (0.);
uf.n[right]  = dirichlet (0.);

/**
## Initial conditions 

We define a reference velocity field *un*. */

scalar un[];

event init (i = 0)
{  
  /**
  We set the physical parameters. */

  mu = muv;
      
  /**
  We set the initial velocity of the fluid. */

  foreach() {
    foreach_dimension()
      u.x[] = 0.;
  }
}

/**
## Properties */

event properties (i++)
{
  /**
  We now set the viscosity and account for the metric. */

  foreach_face()
    muv.x[] = fm.x[];
}

/**
## Mesh adaptation */

event adapt (i++)
{
  adapt_wavelet ((scalar *) {u}, 
		 (double[]) {1.e-3, 1.e-3}, 
		 maxlevel = (lmax), minlevel = (2));
}

/**
## Post-processing and results */

event logfile (t += 0.1; i <= 1000)
{
  double du = change (u.x, un);
  fprintf (ferr, "%d %g %g %g %g %d %d\n",
	   i, t, dt,
	   T1,
	   du,
	   mgp.i, mgu.i);
  fflush (ferr);

  if (i > 0 && du < 1e-5)
    return 1; /* stop */
}

/**
We compute the streamfunction $\psi$ and the unyielded field and
output that to a file. */

event streamfunction (t = end)
{
  scalar psi[], omega[];

  psi[top]    = dirichlet(0);
  psi[bottom] = dirichlet(0);
  psi[left]   = dirichlet(0);
  psi[right]  = dirichlet(0);
  
  foreach() {
    omega[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
    psi[] = 0.;
  }
  boundary ((scalar *) {omega,psi});

  poisson (psi, omega);

  scalar yuyz[];
  foreach()
    yuyz[] = 1.;

  char name[80];
  sprintf (name, "fields-%g", T1);
  FILE * fp = fopen (name, "w");
  output_field ({psi, yuyz}, fp);
  fclose (fp);

  /**
  We also output a vertical cross-section of the velocity
  components. */
  
  sprintf (name, "xprof-%g", T1);
  fp = fopen (name, "w");

  double step = (L0)/(1 << (lmax));
  for (double y = step/2.; y <= (L0) - step/2.; y += step)
    fprintf (fp, "%g %g %g\n", y,
	     interpolate (u.x, (L0)/2., y),
	     interpolate (u.y, (L0)/2., y));
  fclose (fp);

  /**
  And finally we output the maximum of the streamfunction on standard
  error. */
  
  fprintf (fout, "%g %d %g %g\n", t, i, T1, statsf(psi).max);
  fflush (fout);
}

event movie (i += 10)
{
  clear ();
  
  view (fov = 20,
	tx = -0.5, ty = -0.5,
	bg = {1,1,1},
	width = 400, height = 400);
  
  squares ("u.x", linear = false, map = cool_warm);
  cells ();
  save ("ux.mp4");

  squares ("u.y", linear = false, map = cool_warm);
  cells ();
  save ("uy.mp4");
}

/**
## Results

See [Vola et al, 2003, Fig. 3](http://dx.doi.org/10.1016/S0021-9991(03)00118-9).

~~~gnuplot Section of horizontal velocity along the vertical mid-plane.
set key top left spacing 1.1
set xlabel 'y'
set ylabel 'u.x'
set grid
set xrange [0:1]
plot 'xprof-0' w l t 'tau_0=0', \
     '../data/Vola2003/vola2003-3' every :::0::0 lt 1 t 'Vola et al, 2003'
~~~

~~~gnuplot Streamlines and unyielded areas: $\tau_y=0$
reset
set contour base
set cntrparam levels 20
unset surface
set table 'cont.dat'
splot 'fields-0' u 1:2:3 w l
unset table

set size ratio -1
unset key
unset xtics
unset ytics
set palette gray
unset colorbox
set xrange [0:1]
set yrange [0:1]
set cbrange [-1:1]
plot 'fields-0' u 1:2:(1.-$4) w image, 'cont.dat' w l lt -1
~~~
*/
