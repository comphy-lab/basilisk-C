/*
Cyclone test

A test case cyclone with a minimum pressure of XXXX
and a radius of maximum velocity of YYYY
approaches on a path in the -x direction
$$ */
//#include "grid/cartesian.h"
#include "storm_surge.h"

#define MAXLEVEL 7
#define MINLEVEL 4
#define ETAE 1.e-3
#define HINVE 1.e-1


//double t1=604800.;
double t1=10*24*3600.;
double t_ramp=86400.;

struct bathy {
  int (* iterate) (void);
};

void Init_bathy (struct bathy p) {
  int counti=0;
  do {
    foreach() {
      zb[] = (x < 200000. ? 40 - 0.0002 * x : ( x < 300000. ? 200. - 0.001 * x : 1400 - 0.005 * x)) ;
      h[] = max(0., - zb[]);
	}
  } while (p.iterate && p.iterate() && counti++<20);
}

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
  L0 = 2000000.;

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
  scalar eta[],wet[];
  foreach() {
    eta[] = h[] > dry ? h[] + zb[] : 0.;
    wet[] = h[] > dry ? 1. : 0.;
  };
  boundary ({eta,wet});
  astats s = adapt_wavelet ({eta,wet}, (double[]){ETAE,HINVE},
			    MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}

struct cyclone_param {
  double x;
  double y;
  double min_press;
  double R_maxwind;
  double maxwind;
};

struct cyclone_param cyclone_parameters(){
  struct cyclone_param Cparam;
  Cparam.x=750000.-1.*t;
  Cparam.y=1000000.;
  Cparam.min_press=95000.;  //in Pascals - divide by 100 to get hectoPascals
  Cparam.R_maxwind=50000.;
  Cparam.maxwind=50.;
  return Cparam;
};

void cyclone_fields(){
  // Get cyclone parameters for this timestep
  struct cyclone_param Cparam = cyclone_parameters();
  double x0=Cparam.x;
  double y0=Cparam.y;
  double min_press=Cparam.min_press;
  double R_maxwind=Cparam.R_maxwind;
  double maxwind=Cparam.maxwind;
  double press,windspeed,rr,phi,alpha,BB,km,KK,Pref,rhoaw;
  // Convert cyclone parameters into wind and pressure fields
  KK=2.5e-3;
  km=1.;
  BB=2.-(min_press-90000.)/16000.;
  Pref=101300.;
  rhoaw=1.216e-3;
  double ramp = t < t_ramp ? t/t_ramp : 1.;
  foreach(){
    rr=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))+0.01;
    phi=atan2(x-x0,y-y0);
    alpha= rr<3.*R_maxwind ? rr/(3.*R_maxwind)*pi/6. : pi/6.;
    press=(Pref-min_press)*exp(-pow(R_maxwind/rr,BB))+min_press;
    windspeed=km*maxwind*sqrt(pow(R_maxwind/rr,BB)*exp(1.-pow(R_maxwind/rr,BB)));
    pa[]=ramp*(press-Pref);
    //ts.x[]=ramp*rhoaw*KK*windspeed*windspeed*cos(phi+alpha);
    //ts.y[]=-ramp*rhoaw*KK*windspeed*windspeed*sin(phi+alpha);
    ts.x[]=( h[] < 1. ? 0. : dt/h[]*ramp*rhoaw*KK*windspeed*windspeed*cos(phi+alpha));
    ts.y[]=( h[] < 1. ? 0. : -dt/h[]*ramp*rhoaw*KK*windspeed*windspeed*sin(phi+alpha));
  }
}


event initiate(i=0)
{
  Init_bathy( iterate = adapt );
  foreach() {
    fcor[]=-5.5e-5+1.7e-11*y;
  }
  // Get wind and pressure fields
  cyclone_fields () ;
  boundary ({h,u,eta,fcor,zb,pa,ts}); 
}

/**
We output running statistics to the standard error. */
event logfile (i+=10; t<t1) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%12.8g %d %12.8g %12.8g %12.8g %12.8g %12.8g %12.8g\n", t, i,
	   s.min, s.max, s.sum, n.rms, n.max, dt);
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
    fcor[]=-5.5e-5+1.7e-11*y;
    // Get wind and pressure fields
    cyclone_fields () ;
    double a_inv = (h[] < dry ? 0. : h[]/(h[] + 2.5e-3*dt));
    foreach_dimension()
      u.x[] = u.x[] * a_inv;
  }
  boundary ({h,u,fcor,zb,pa,ts});
}

/**
## Boundary conditions

We set the normal velocity component on the right, top and bottom
boundaries to a "radiation condition" with a reference sealevel of
zero. The left boundary is always "dry" in this example so can be left
alone. Note that the sign is important and needs to reflect the
orientation of the boundary. */

//u.x[right]  = + radiation(-Pa2m*pa[]);
//u.y[bottom] = - radiation(-Pa2m*pa[]);
//u.y[top]    = + radiation(-Pa2m*pa[]);


  
/**
Every 12 hours, various fields are interpolated
bilinearly onto a `n x n` regular grid and written on standard
output. */

event snapshots (t += 3600.) {
  scalar ux[], uy[], tsx[], tsy[];
  int hour = (int) t/3600.;
  foreach(){
    ux[]=u.x[];
    uy[]=u.y[];
    tsx[]=ts.x[];
    tsy[]=ts.y[];
  }
  printf ("%% file: t-%d\n", hour);
  output_field ({h, eta, zb, ux, uy, pa, tsx, tsy, fcor}, stdout, n = 1 << MAXLEVEL, linear = true);
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
  output_ppm (etam, fp, min = -1., max = 1. , n = 512, linear = true);
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
