/**
## Spatially Varying Pressure
Example of a constant spatially varying pressure over a square 
constant depth water body


A pressure field of $Pa(x,y)-P_{ref}=Amp*sin(2*\pi*x/L)*cos(2*\pi*y/L)$ is 
applied to a square pool

If grid/cartesian1D.h is specified we are working only in 1D
 */
//#include "grid/cartesian1D.h"
#include "storm_surge.h"

#define MAXLEVEL 8
#define MINLEVEL 4
#define ETAE 1e-8

/**
Inverse barometer conversion hPa to metres
$$
\eta = - \frac{\Delta p}{\rho_w g}
$$
Thus for pressure in hectoPascals the coefficient is $\frac{100}{\rho_w g}$. */

double rho_inv = 1./1025.;
double amp = 200.;

double compare (scalar v, scalar vn)
{
  double maxC = 0.;
  foreach(reduction(max:maxC)) {
    double dv = fabs (v[] - vn[]);
    if (dv > maxC)
      maxC = dv;
  }
  return maxC;
}

/**
We declare variables that are used to store the pressure 
and its gradient */


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
## Boundary conditions

We set the normal velocity component on the left, right bottom and top
boundaries to a "radiation condition" with a reference sealevel of
inverse barometer. Note that the sign is important and needs to reflect the
orientation of the boundary. */

u.n[left]   = - radiation(-Pa2m*(pa[]));
u.n[right]  = + radiation(-Pa2m*(pa[]));
u.n[bottom] = - radiation(-Pa2m*(pa[]));
u.n[top]    = + radiation(-Pa2m*(pa[]));

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
`?`) which we set to $h+z_b$ but only for "wet" cells. If we used
$h+z_b$ everywhere (i.e. the default $\eta$ provided by the
Saint-Venant solver) we would also refine the dry topography, which is
not useful. */

int adapt() {
#if QUADTREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});
  
/**
$h+z_b$ everywhere (i.e. the default $\eta$ provided by the
Saint-Venant solver) we would also refine the dry topography, which is
not useful.

We can now use wavelet adaptation on the list of scalars `{?,hmax}`
with thresholds `{ETAE,HMAXE}`. The compiler is not clever enough yet
and needs to be told explicitly that this is a list of `double`s,
hence the `(double[])`
[type casting](http://en.wikipedia.org/wiki/Type_conversion) . 

The function then returns the number of cells refined. */
  
  astats s = adapt_wavelet ({eta}, (double[]){ETAE},
			    MINLEVEL, MAXLEVEL);
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

  double ramp = t < 120. ? t/120. : 1.;
  foreach() {
    zb[] = -2;
    h[] = max(0., - zb[]);
    pa[] = ramp*amp*sin(2.*pi*x/L0)*cos(2.*pi*y/L0);
  }
 boundary ({h,zb}); 
}

/**
We want the simulation to stop when we are close to steady state. To
do this we store the `h` field of the previous timestep in an
auxilliary variable `hn`. */

scalar hn[],etaref[];

event init_hn (i = 0) {
  foreach() {
    hn[] = h[];
    etaref[] = -Pa2m*pa[];
  }
}

/**
Every 10 timesteps we check whether $u$ has changed by more than
10^-8^. If it has not, the event returns 1 which stops the
simulation. We also output running statistics to the
standard error. */

event logfile (i+= 10; i <= 10000) {
  double dh = change (h, hn);
  if ( (i > 100 && dh < 1e-8 ) || i==10000) {
    double dhref = compare (eta,etaref);
    fprintf (stderr, "# Max difference is %.9g\n# Time to convergence %g",
	     dhref,t); 
    FILE * fp = fopen("diff", "w");
    fprintf(fp,"# 1:x 2:y 3:h 4:zb 5:eta 6:etaref 7:diff\n");
    for (double x = -500; x <= 500; x += 10) {
      for (double y = -500; y <= 500; y += 10)
	fprintf (fp, "%g %g %.9g %.9g %.9g %.9g %.9g\n", x, y, 
		 interpolate (h, x, y), 
		 interpolate (zb, x, y), 
		 interpolate (eta, x, y), 
		 interpolate (etaref, x, y), 
		 interpolate (etaref, x, y) - interpolate (eta, x, y));
      fprintf( fp, "\n");
    }
    fclose (fp);
    return 1; /* stop */
  }

  stats s = statsf (h);
  norm n = normf (u.x);
  double detaref = compare (eta,etaref);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dh detaref dt\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g %g %g\n", t, i, s.min, s.max, 
	   s.sum, n.rms, n.max, dh, detaref, dt);
}

/**
We also use a simple implicit scheme to implement quadratic bottom
friction i.e.
$$
\frac{d\mathbf{u}}{dt} = - C_f|\mathbf{u}|\frac{\mathbf{u}}{h}
$$
with $C_f=10^{-4}$. 
We ramp up the pressure from zero*/

event source (i++) {
  double ramp = t < 120. ? t/120. : 1.;
  foreach() {
    pa[] = ramp*amp*sin(2.*pi*x/L0)*cos(2.*pi*y/L0);
    etaref[] = -Pa2m*pa[]; 
    eta[] = h[] > dry ? h[] + zb[] : 0.;
  }
  foreach() {
    double a_inv = h[] < dry ? 0. : h[]/(h[] + 1e-4*dt*norm(u));
    //double a = h[] < dry ? HUGE : 1. + 1e-2*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] = u.x[]*a_inv;
  }
  boundary ({h,u,zb,pa,etaref});
}

/**
Every 5 minutes, the $h$, $z_b$ and `hmax` fields are interpolated
bilinearly onto a `n x n` regular grid and written on standard
output. */

event snapshots (t += 300) {
  printf ("# file: t-%g\n", t);
  output_field ({eta}, stdout, n = 1 << MAXLEVEL, linear = true);
  scalar l[];
  foreach()
    l[]=level;
  printf ("# file: level-%g\n", t);
  output_field ({l});
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
  output_ppm (etam, fp, min = -0.002, max = 0.002 , n = 512, linear = true);
  //  output_ppm (etam, fp, n = 512, linear = true);
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

/**
## Results

After running and processing by gnuplot (using [pressure_pool.plot]()) we get
the following pictures and animations.

![[Evolution of the vorticity field with time.](pressure_pool/eta.mpg)](pressure_pool/plot.png)

![Difference between calculated and reference surface level](pressure_pool/diff.png) */
