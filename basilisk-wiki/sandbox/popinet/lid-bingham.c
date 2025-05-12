/**
# Lid-driven cavity of a Bingham fluid at Re=0

We use the multigrid implementation (rather than the default tree
implementation) and the yield-stress fluid solver. */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "yield.h"

double tau_y = 0;

static void run_tauy()
{
  const face vector tauc[] = {tau_y,tau_y};
  tau_0 = tauc;
  run();  
}

int main()
{ 
  init_grid (128);
  // viscosity
  double mu0 = 1.;
  const face vector muc[] = {mu0,mu0};
  mu = muc;
  // Stokes flow
  stokes = true;
  // the "r" parameter from the augmented Lagrangian formulation
  r = 10.*mu0;

  /**
  We vary the yield-stress and run a simulation for each value. */

  tau_y = 0; run_tauy();
  tau_y = 1./sqrt(2); run_tauy();
  tau_y = 10./sqrt(2); run_tauy();
  tau_y = 100./sqrt(2); run_tauy();
  tau_y = 1000./sqrt(2); run_tauy();
    
#if 0 // uncomment this to get more points for max(psi) vs tau_0
  for (tau_y = 0.1; tau_y <= 1000; tau_y *= 1.5)
    run_tauy ();
#endif
}

/**
The default boundary conditions are symmetry (i.e. slip walls). We
need no-slip on three boundaries and $u=1$ on the top
boundary i.e. */

u.t[top] = dirichlet(1);

/**
For the other no-slip boundaries this gives */

u.t[bottom] = dirichlet(0);
u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);

/**
For the colocated solver, imposing boundary conditions for the normal
components of the (face-centered) advection velocity improves the
results. Ideally, this should be done automatically by the solver. */

uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[top]    = 0;
uf.n[bottom] = 0;

/**
We want the simulation to stop when we are close to steady state. To
do this we store the `u.x` field of the previous timestep in an
auxilliary variable `un`. */

scalar un[];

/**
Every 0.1 time units we check whether $u$ has changed by more than
10^-5^. If it has not, the event returns 1 which stops the
simulation. We also output the evolution of the difference and
multigrid iterations counts on standard output. */

event logfile (t += 0.1; i <= 10000) {
  double du = change (u.x, un);
  if (i > 0 && du < 1e-5)
    return 1; /* stop */
  fprintf (stdout, "du %f %d %g %d %d\n", t, i, du, mgp.i, mgu.i);
  fflush (stdout);
}

/**
This event will happen after completion of the simulation. We compute
the streamfunction $\psi$ and the unyielded field and output that to a
file. */

event streamfunction (t = end)
{
  scalar psi[], omega[];

  psi[top] = dirichlet(0);
  psi[bottom] = dirichlet(0);
  psi[left]   = dirichlet(0);
  psi[right]  = dirichlet(0);
  
  foreach() {
    omega[] = (u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta);
    psi[] = 0.;
  }
  boundary ({omega,psi});

  poisson (psi, omega);

  scalar unyielded[];
  foreach()
    unyielded[] = (d.x.x[] == 0. && d.x.x[1] == 0. &&
		   d.y.y[] == 0. && d.y.y[0,1] == 0.);

  char name[80];
  sprintf (name, "fields-%g", tau_y);
  FILE * fp = fopen (name, "w");
  output_field ({psi, unyielded}, fp);
  fclose (fp);

  /**
  We also output a vertical cross-section of the velocity components. */
  
  sprintf (name, "xprof-%g", tau_y);
  fp = fopen (name, "w");
  for (double y = 0; y <= 1; y += 0.01)
    fprintf (fp, "%g %g %g\n", y,
	     interpolate (u.x, 0.5, y), interpolate (u.y, 0.5, y));
  fclose (fp);

  /**
  And finally we output the maximum of the streamfunction on standard
  error. */
  
  fprintf (stderr, "%g %d %g %g\n", t, i, tau_y, statsf(psi).max);
}

/**
## Results

~~~gnuplot Vortex intensity
set xlabel 'tau_0/sqrt(2) (Vola) or tau_0'
set ylabel 'max(psi)'
set logscale x
set xrange [0.1:1000]
set grid
plot '../vola2003-2a' u ($1/sqrt(2)):2 t 'Vola et al, 2003, fig 2a', \
     'log' u 3:4  w p t 'Basilisk'
~~~

~~~gnuplot Number of timesteps
set xlabel 'tau_0'
set ylabel ''
plot 'log' u 3:2  w p t ''
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

~~~gnuplot Streamlines and unyielded areas: $\tau_y=1/\sqrt{2}$
set table 'cont.dat'
splot 'fields-0.707107' u 1:2:3 w l
unset table
plot 'fields-0.707107' u 1:2:(1.-$4) w image, 'cont.dat' w l lt -1
~~~

~~~gnuplot Streamlines and unyielded areas: $\tau_y=10/\sqrt{2}$
set table 'cont.dat'
splot 'fields-7.07107' u 1:2:3 w l
unset table
plot 'fields-7.07107' u 1:2:(1.-$4) w image, 'cont.dat' w l lt -1
~~~

~~~gnuplot Streamlines and unyielded areas: $\tau_y=100/\sqrt{2}$
set table 'cont.dat'
splot 'fields-70.7107' u 1:2:3 w l
unset table
plot 'fields-70.7107' u 1:2:(1.-$4) w image, 'cont.dat' w l lt -1
~~~

~~~gnuplot Streamlines and unyielded areas: $\tau_y=1000/\sqrt{2}$
set table 'cont.dat'
splot 'fields-707.107' u 1:2:3 w l
unset table
plot 'fields-707.107' u 1:2:(1.-$4) w image, 'cont.dat' w l lt -1
~~~

See [Vola et al, 2003, Fig. 3](http://dx.doi.org/10.1016/S0021-9991(03)00118-9).

~~~gnuplot Section of horizontal velocity along the vertical mid-plane.
reset
set xlabel 'y'
set ylabel 'u.x'
set grid
set key top left
set xrange [0:1]
plot 'xprof-0' w l t 'tau_0=0', \
     'xprof-0.707107' w l t 'tau_0=1/sqrt(2)', \
     'xprof-7.07107' w l t 'tau_0=10/sqrt(2)', \
     'xprof-70.7107' w l t 'tau_0=100/sqrt(2)', \
     'xprof-707.107' w l t 'tau_0=1000/sqrt(2)' lt -1, \
     '../vola2003-3' every :::0::0 lt 1 t 'Vola et al, 2003', \
     '../vola2003-3' every :::1::1 t '', \
     '../vola2003-3' every :::2::2 t '', \
     '../vola2003-3' every :::3::3 t '', \
     '../vola2003-3' every :::4::4 t ''
~~~
*/
