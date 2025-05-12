/*
Western Boundary Current intensification

This is using Stommel's western boundary current as a test case  
Stommel, Henry (April 1948). "The Westward Intensification of Wind-Driven Ocean Currents". Transactions, American Geophysical Union 29 (2): 202–206. doi:10.1029/tr029i002p00202
See http://ido.at.fcen.uba.ar/index_archivos/Stommel_1948.pdf

A sinusoidal wind with wind stress in the x direction of $1.7 \times 10^{-4} \cos(\pi y / L0)$ m/s blows in the x direction over a square pool of L0 = 1000 km side length with constant depth 2 km.  There is linear friction with coefficient 2.5e-3.
A beta-plane approximation is made with
f=-9.5e-5+1.7e-11*y
Which corresponds to the origin being at 41 degrees South.
$$ */
#include "grid/cartesian.h"
#include "storm_surge.h"

#define MAXLEVEL 8
#define MINLEVEL 4
#define ETAE 1e-8

double t1=8*604800.;

/**
Here we can set standard parameters in Basilisk */

int main()
{
/**
Here we setup the domain geometry. For the moment Basilisk only
supports square domains. This case uses metres east and north. We set
the size of the box `L0` and the origin to be the lower-left corner
`(X0,Y0)`.
In this case we are assuming a square 'pool' of length 1000 km. */
  
  // the domain is
  L0 = 1000000.;

/**
`G` is the acceleration of gravity required by the Saint-Venant
solver. This is the only dimensional parameter.. */

  G = 9.81;

#if QUADTREE
  // 32^2 grid points to start with
  init_grid( 1 << MINLEVEL );
#else // Cartesian
  // 1024^2 grid points
  init_grid( 1 << MAXLEVEL );
#endif

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
  foreach() {
    fcor[]=-9.5e-5+1.7e-11*y;
    zb[] = -2000.;
    h[] = max(0., - zb[]);
    ts.x[] = -1.7e-4*cos(pi*y/L0);
    ts.y[] = 0.;
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
We output running statistics to the standard error. */
event logfile (i+=10; t<t1) {
  double dh = change (h, hn);
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dh  dt\n");
  fprintf (stderr, "%12.8g %d %12.8g %12.8g %12.8g %12.8g %12.8g %12.8g %12.8g\n", t, i,
	    s.min, s.max, s.sum, n.rms, n.max, dh, dt);
}
/**
We also use a simple implicit scheme to implement linear bottom
friction i.e.
$$
\frac{d\mathbf{u}}{dt} = - C_f\frac{\mathbf{u}}{h}
$$
with $C_f=2.5 \times 10^{-3}$. */
event source (i++) {
  foreach() {
    fcor[]=-9.5e-5+1.7e-11*y;
    ts.x[] = - 1.7e-4*cos(pi*y/L0);
    ts.y[] = 0.;
    double a_inv = (h[] < dry ? 0. : h[]/(h[] + 2.5e-3*dt));
    foreach_dimension()
      u.x[] = u.x[] * a_inv;
  }
  boundary ({h,u,fcor,ts});
}
  
/**
Every 12 hours, various fields are interpolated
bilinearly onto a `n x n` regular grid and written on standard
output. */

event snapshots (t += 86400.) {
  scalar ux[], uy[], tsx[], tsy[];
  int day = (int) t/86400.;
  foreach(){
    ux[]=u.x[];
    uy[]=u.y[];
    tsx[]=ts.x[];
    tsy[]=ts.y[];
  }
  printf ("%% file: t-%d\n", day);
  output_field ({h, eta, zb, ux, uy, tsx, tsy, fcor}, stdout, n = 1 << MAXLEVEL, linear = true);
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
  output_ppm (etam, fp, min = -0.01, max = 0.01 , n = 512, linear = true);
#if QUADTREE
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
