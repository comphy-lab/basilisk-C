/**
## Constant wind over constant depth pool

Example of a constant wind over a square constant depth pool

A constant wind of 1 m/s blows in the x direction over a square pool of 1000 m side length with constant depth 2 m.

We assume quadratic bottom friction with coefficient 2.5e-3. Wind stress $\tau$ is defined by the formula
$$
\frac{\tau}{\rho} = C_{10} U_{10}^2,
$$
where $C_{10}$ is taken from [Wu (1982) JGRC](http://champs.cecs.ucf.edu/Library/Journal_Articles/pdfs/Wind-stress%20coefficients%20over%20sea%20surface%20from%20breeze%20to%20hurricane.pdf)
$$
C_{10}=(0.8 + 0.065 U_{10}) \times 10^{-3}.
$$
So for a $U_{10}=10 m/s$ we have a wind stress of
$$
\frac{\tau}{\rho} = 0.145.
$$ */
//#include "grid/cartesian1D.h"
#include "sandbox/Emily/storm_surge.h"

#define MAXLEVEL 8
#define MINLEVEL 4
#define ETAE 1e-8

double t1=7200.;

/**
Here we can set standard parameters in Basilisk */

int main()
{
#if QUADTREE
  // 32^2 grid points to start with
  N = 1 << MINLEVEL;
#else // Cartesian
  // 1024^2 grid points
  N = 1 << MAXLEVEL;
#endif
/**
Here we setup the domain geometry. For the moment Basilisk only
supports square domains. This case uses metres east and north. We set
the size of the box `L0` and the coordinates of the lower-left corner
`(X0,Y0)`.
In this case we are assuming a square 'pool' of length 1000 m. */
  
  // the domain is
  L0 = 1000.;
  X0 = -L0/2.;
  Y0 = -L0/2.;

/**
`G` is the acceleration of gravity required by the Saint-Venant
solver. This is the only dimensional parameter.. */

  G = 9.81;
  run();
}

/**
## Adaptation

Here we define an auxilliary function which we will use several times
in what follows. Again we have two `#if...#else` branches selecting
whether the simulation is being run on an (adaptive) quadtree or a
(static) Cartesian grid.

We want to adapt according to two criteria: an estimate of the error
on the free surface position -- to track the wave in time -- and an
estimate of the error on the maximum wave height `hmax` -- to make
sure that the final maximum wave height field is properly resolved.

We first define a temporary field (in the
[automatic variable](http://en.wikipedia.org/wiki/Automatic_variable)
`η`) which we set to $h+z_b$ but only for "wet" cells. If we used
$h+z_b$ everywhere (i.e. the default $\eta$ provided by the
Saint-Venant solver) we would also refine the dry topography, which is
not useful. */

int adapt() {
#if QUADTREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  astats s = adapt_wavelet ({eta}, (double[]){ETAE},
			    MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}

event initiate(i=0)
{
/**
The initial still water surface is at $z=0$ so that the water depth
$h$ is... */

  foreach() {
    zb[] = -2;
    h[] = max(0., - zb[]);
    ts.x[] = 0;
    ts.y[] = 0;
  }
  boundary ({h,zb}); 
}
/**
We want the simulation to stop when we are close to steady state. To
do this we store the `h` field of the previous timestep in an
auxilliary variable `hn`. */

scalar hn[];

event init_hn (i = 0) {
  foreach() {
    hn[] = h[];
  }
}

/**
Every timestep we check whether $h$ has changed by more than
10^-8^. If it has not, the event returns 1 which stops the
simulation. We also output running statistics to the
standard error. */
event logfile (i++; i <= 100000) {
  double dh = change (h, hn);
  if ( (t> t1 && dh < 1e-8) || i==100000) {
  foreach()
    fprintf (stderr, "%g %g\n", x, h[]);
    return 1; /* stop */
  }
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dh  dt\n");
  fprintf (stderr, "%12.8g %d %12.8g %12.8g %12.8g %12.8g %12.8g %12.8g %12.8g\n", t, i,
	    s.min, s.max, s.sum, n.rms, n.max, dh, dt);
}
/**
We also use a simple implicit scheme to implement quadratic bottom
friction i.e.
$$
\frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
$$
with $C_f=2.5 \times 10^{-3}$. 

Also assume that we have a constant wind blowing in the x direction*/
event source (i++) {
  double ramp = t < t1 ? t/t1 : 1.;
  foreach() {
    ts.x[] = ramp*1.7e-4;
    ts.y[] = 0;
    double a_inv = (h[] < dry ? 0. : h[]/(h[] + 2.5e-3*dt*norm(u)));
    //double a = h[] < dry ? HUGE : 1. + 2.5e-3*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] = u.x[] * a_inv;
  }
  boundary ({h,u});
}
  
/**
Every 5 minutes, the $h$, $z_b$ and `hmax` fields are interpolated
bilinearly onto a `n x n` regular grid and written on standard
output. */

event snapshots (t += 300) {
  scalar ux[], uy[];
  foreach(){
    ux[]=u.x[];
    uy[]=u.y[];
  }
  printf ("%% file: t-%g\n", t);
  output_field ({h, eta, zb, ux, uy}, stdout, n = 1 << MAXLEVEL, linear = true);
}

event movies (t+=60) {
  static FILE * fp = NULL;
  if (!fp) fp = popen ("ppm2mpeg > eta.mpg", "w");
  scalar m[], etam[];
  foreach() {
    etam[] = eta[]*(h[] > dry);
    m[] = etam[] - zb[];
  }
  boundary ({m, etam});
  output_ppm (etam, fp, min = -0.005, max = 0.005 , n = 512, linear = true);
  //  output_ppm (etam, fp, n = 512, linear = true);
#if 0
  static FILE * fp1 = NULL;
  if (!fp1) fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l = etam;
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = MINLEVEL, max = MAXLEVEL, n = 512);
#endif
}

/**
## Adaptivity

We apply our `adapt()` function at every timestep. */

event do_adapt (i++) adapt();
