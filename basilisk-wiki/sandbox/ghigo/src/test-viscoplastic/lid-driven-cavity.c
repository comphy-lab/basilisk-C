/**
# Lid-driven cavity of a Bingham fluid at Re=0 */

#include "grid/multigrid.h"
#include "../mycentered2.h"
#include "../myviscosity-viscoplastic.h"

/**
## Reference solution */

#define l    (1.) // Reference length
#define nu   (1.) // Viscosity
#define uref (1.) // Reference velocity
#define tref (sq (l)/(nu)) // Reference time

/**
## Set up

We need a field for the viscosity and the yield stress. It also seems
that we need a field for the density, but this may be a bug. */

face vector muv[], Tv[];
scalar rhov[];
double T1 = 0.;
int iT = 0;

/**
We also need a reference velocity field *un*. */

scalar un[];

/**
Finally, we define the level of refinement. */

#define lvl (7)

/**
## Code */

int main()
{
  /**
  The domain is $1\times 1$. */

  L0 = (l);
  size (L0);
  origin (0., 0.);

  /**
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
  
  /**
  We set the tolerance of the Poisson solver. */

  stokes       = true;
  TOLERANCE    = 1.e-5;
  TOLERANCE_MU = 1.e-5*(uref);
  NITERMAX     = 500;
  
  /**
  We initialize the grid. We start with a Newtonian fluid. */

  N = 1 << (lvl);
  init_grid (N);

  T1 = 0.;
  
  run ();
}

/**
## Boundary conditions */

u.t[top]    = dirichlet(1);
u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver and to make
sure the compatibility condition for the Poisson equation is not
violated. */

uf.n[bottom] = (0.);
uf.n[top]    = (0.);
uf.n[left]   = (0.);
uf.n[right]  = (0.);

/**
## Properties */

event properties (i++)
{
  /**
  We set the density and account for the metric. */

  foreach()
    rhov[] = cm[];
  boundary ({rhov});
  
  /**
  We now set the viscosity and yield stress *T* and account for the
  metric. */

  foreach_face() {
    muv.x[] = fm.x[];
    Tv.x[]  = fm.x[]*(T1);
  }
  boundary ((scalar *) {muv, Tv});
}

/**
## Initial conditions */

event init (i = 0)
{  
  /**
  We set the viscosity and density fields in the event
  *properties*. */

  rho = rhov;
  mu  = muv;

  /**
  We also set the yield stress and the regularization parameters. */  
  
  T = new face vector;
  Tv = T;

  /**
  We finally initialize the reference velocity field. */
  
  foreach()
    un[] = u.x[];
}

/**
## Outputs */

event logfile (i++; i <= 1000)
{
  double du = change (u.x, un);
  fprintf (ferr, "%d %g %g %g %g %g %d %d\n",
	   i, t, dt, eps,
	   T1,
	   du,
	   mgp.i, mgu.i);
  fflush (ferr);

  if (i > 0 && du < 1e-5) {

    /**
    We compute the streamfunction $\psi$ and the unyielded field and
    output that to a file. */

    scalar psi[], omega[];

    psi[top]    = dirichlet(0);
    psi[bottom] = dirichlet(0);
    psi[left]   = dirichlet(0);
    psi[right]  = dirichlet(0);
  
    foreach() {
      omega[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
      psi[] = 0.;
    }
    boundary ((scalar *) {omega, psi});

    poisson (psi, omega);
    
    yielded_region ();

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

    double step = (L0)/(N);
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

    /**
    We increase the value of the yield stress and use the previous
    solution as an initial guess. */

    T1 = pow (10., iT)/sqrt (2);

    /**
    To ensure convergence of the multigrid solver for large Bingham
    numbers, we reduce the time step and increase the regularization
    parameter. These parameters were chosen by trial and error. */

    DT  = 1.e-2*(tref)/max (1, pow (10., iT));
    eps = min (5.e-2, max (1.e-4, 5.e-4*pow (10., iT))); 
    
    event ("properties");
    boundary (all);

    if (iT == 4)
      return 1; /* stop */
    iT++;
  }
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
     'xprof-0.707107' w l t 'tau_0=1/sqrt(2)', \
     'xprof-7.07107' w l t 'tau_0=10/sqrt(2)', \
     'xprof-70.7107' w l t 'tau_0=100/sqrt(2)', \
     'xprof-707.107' w l t 'tau_0=1000/sqrt(2)', \
     '../data/Vola2003/vola2003-3' every :::0::0 lt 1 t 'Vola et al, 2003', \
     '../data/Vola2003/vola2003-3' every :::1::1 t '', \
     '../data/Vola2003/vola2003-3' every :::2::2 t '', \
     '../data/Vola2003/vola2003-3' every :::3::3 t '', \
     '../data/Vola2003/vola2003-3' every :::4::4 t ''
~~~

~~~gnuplot Vortex intensity
set key top right
set xlabel 'tau_0/sqrt(2) (Vola) or tau_0'
set ylabel 'max(psi)'
set grid
set xrange [0.1:1000]
set logscale x
plot '../data/Vola2003/vola2003-2a' u ($1/sqrt(2)):2 t 'Vola et al, 2003, fig 2a', \
     'out' u 3:4  w p t 'Basilisk'
~~~

~~~gnuplot Streamlines for $\tau_y=0$
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

~~~gnuplot Streamlines for $\tau_y=1/\sqrt{2}$
set table 'cont.dat'
splot 'fields-0.707107' u 1:2:3 w l
unset table
plot 'fields-0.707107' u 1:2:(1.-$4) w image, 'cont.dat' w l lt -1
~~~

~~~gnuplot Streamlines for $\tau_y=10/\sqrt{2}$
set table 'cont.dat'
splot 'fields-7.07107' u 1:2:3 w l
unset table
plot 'fields-7.07107' u 1:2:(1.-$4) w image, 'cont.dat' w l lt -1
~~~
*/
