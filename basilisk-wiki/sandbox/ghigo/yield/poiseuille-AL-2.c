/**
# 2D flow Poiseuille flow of a Bingham fluid */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "src/yield-AL-2.h"

/**
## Yield stress */

#define T1 (2.5e-1)

/**
We need a face-centered field for the viscosity and the yield stress
*T*. */

face vector muv[], Tv[];

/**
## Code */

int main()
{ 
  /**
  The domain is $0.5\times 0.5$. */

  L0 = 0.5;
  size (L0);
  origin (0., 0.);

  periodic (left);

  /**
  We set the numerical parameters. */

  DT = 1.e-2;
  
  /**
  We tune the Poisson solver. */

  stokes = true;
  TOLERANCE = 1.e-5;

  for (N = 8; N <= 128; N *= 2) {
    init_grid (N);
    run();
  }
}

/**
## Boundary conditions */

//y=0 (no-slip)
u.n[bottom]  = dirichlet (0.);
u.t[bottom]  = dirichlet (0.);

//y=L0 (default slip)
u.n[top]  = dirichlet (0.);
u.t[top]  = neumann (0.);

/**
We also make sure there is no flow through the top and bottom
boundary, otherwise the compatibility condition for the Poisson
equation can be violated. */

uf.n[bottom] = (0.);
uf.n[top]    = (0.);

/**
## Initial conditions 

We define a reference velocity field *un*. */

scalar un[];

event init (i = 0)
{
  /**
  We set the physical parameters. */

  const face vector av[] = {1., 0.};
  a = av;

  mu = muv;
  
  T = new face vector;
  Tv = T;

  r = 30.;
  
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
  We now set the viscosity and yield stress *T* and account for the
  metric. */

  foreach_face() {
    muv.x[] = fm.x[];
    Tv.x[]  = fm.x[]*(T1);
  }
  boundary ((scalar *) {Tv});
}

/**
## Post-processing and results */

event logfile (t += 0.1; i < 1000) 
{
  double du = change (u.x, un);
  /* fprintf(ferr, "%d %g %g %g\n", */
  /* 	  i, t, dt, */
  /* 	  du */
  /* 	  ); */
  /* fflush (ferr); */
  
  if (i > 0 && du < 1.e-6)
    return 1; /* stop */
}

double velocity (double y)
{
  double h = 0.5, dpdx = 1., mu = 1.;
  double U = sq((T1) - h*dpdx)/(2.*mu*dpdx);
  double Y = h - (T1)/dpdx;
  return (y < Y ? dpdx*(h*y - sq(y)/2.) - (T1)*y : U);
}

double shear_stress (double y)
{
  double h = 0.5, dpdx = 1.;
  double Y = h - (T1)/dpdx;
  return (y < Y ? Y - y : 0);
}

event profile (t = end)
{
  if (N == 16) {
    scalar shear[];
    foreach()
      shear[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    boundary ((scalar *) {shear});

    double step = (L0)/(N);
    for (double y = 0. + step/2.; y < (L0) - step/2.; y += step)
      fprintf (fout, "%g %g %g %g %g\n",
	       y,
	       interpolate (u.x, (L0)/2, y),
	       velocity (y),
	       interpolate (shear, (L0)/2, y),
	       shear_stress (y)
	       );
    fflush (fout);
  }
  
  scalar e[];
  foreach()
    e[] = u.x[] - velocity (y);
  norm n = normf (e);
  fprintf (ferr, "%d %g %g %g %g\n",
	   N, n.avg, n.rms, n.max, delta_D());
  fflush (ferr);
}

event snapshot (t = end)
{
  if (N == 16) {
    clear ();

    view (fov = 20,
	  tx = -0.5, ty = -0.5,
	  bg = {1,1,1},
	  width = 400, height = 400);
    
    squares ("u.x", min = 0, max = 1./32.,
	     linear = false, map = cool_warm);
    cells ();
    save ("ux.png");

    squares ("yuyz", min = 0, max = 1.,
	     linear = false, map = cool_warm);
    cells ();
    save ("yield.png");

    scalar e[];
    foreach()
      e[] = fabs (u.x[] - velocity (y));
      
    squares ("e", spread = -1, linear = false, map = cool_warm);
    cells ();
    save ("err.png");
  }
}

/**
## Results

![*u.x* velocity field](poiseuille-AL-2/ux.png)

![Yielded region](poiseuille-AL-2/yield.png)

![Error](poiseuille-AL-2/err.png)

#### Profiles

~~~gnuplot Velocity profile
set key top left spacing 1.1
set xlabel 'y'
set ylabel 'u_x'
set grid
set xrange[0:0.5]

plot 'out' u 1:3 w p ps 0.75 pt 6 lc rgb "black" t "Exact", \
     'out' u 1:2 w l lw 2 lc rgb "blue" t 'AL, mg'
~~~

~~~gnuplot Shear rate profile
set key top right
set ylabel 'gamma'

plot 'out' u 1:5 w p ps 0.75 pt 6 lc rgb "black" t "Exact", \
     'out' u 1:4 w l lw 2 lc rgb "blue" t 'AL, mg'
~~~

~~~gnuplot Error
reset
set xlabel 'N'
set ylabel 'error'
set xtics 8,2,128
set grid
set logscale

ftitle(a,b) = sprintf("order %4.2f", -b)
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($3)) via a2,b2
fm(x)=am+bm*x
fit fm(x) 'log' u (log($1)):(log($4)) via am,bm

plot exp (fm(log(x))) t ftitle(am,bm), \
     exp (f2(log(x))) t ftitle(a2,b2), \
     'log' u 1:4 t '|h|_{max}' ps 1.5, \
     'log' u 1:3 t '|h|_2' ps 1.5
~~~

*/
