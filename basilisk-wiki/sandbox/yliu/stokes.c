/**
# Breaking wave (Forked from popinet/wave.c)

We solve the two-phase Navier--Stokes equations with surface tension
and using a momentum-conserving transport of each phase. Gravity is
taken into account using the "reduced gravity approach" and the
results are visualised using Basilisk view. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
//#include "tension.h"
#include "reduced.h"
#include "view.h"
#include "tag.h"

/**
We log some profiling information. */

#include "navier-stokes/perfs.h"
#include "profiling.h"
#include "output_htg.h"
/**
The primary parameters are the wave steepness $ak$, the Bond and
Reynolds numbers. */

double ak = 0.1;
double BO = 1000.;
double RE = 57470.;

/**
The default maximum level of refinement depends on the dimension. */

int LEVEL = dimension == 2 ? 11 : 7;

/**
The error on the components of the velocity field used for adaptive
refinement. */

double uemax = 0.005;

/**
The density and viscosity ratios are those of air and water. */

#define RATIO (1./850.)
#define MURATIO (17.4e-6/8.9e-4)

/**
Define if we want to use a Dirac viscous layer initialization. */

int DIRAC = 0;

/**
The wave number, fluid depth and acceleration of gravity are set to
these values. */

#define k_  (2.*pi)
#define h_   (2.*pi/k_)
#define g_   1.0
#define T_   (2.*pi/sqrt(g_*k_))

#define L0_  (4*2.*pi/k_)

/**
Boundary conditions. */
u.n[bottom] = 0.;

/**
The program takes optional arguments which are the level of
refinement, steepness, Bond and Reynolds numbers, and optional Dirac
initialisation. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  if (argc > 2)
    ak = atof(argv[2]);
  if (argc > 3)
    BO = atof(argv[3]);
  if (argc > 4)
    RE = atof(argv[4]);    
  if (argc > 5)
    DIRAC = atof(argv[5]);

  /**
  The domain is a cubic box centered on the origin and of length
  $L0=1$, periodic in the x- and z-directions. */
   L0 = L0_;
  origin (-L0/2, -h_, -L0/2);
  periodic (right);
#if dimension > 2
  periodic (front);
#endif

  /**
  Here we set the densities and viscosities corresponding to the
  parameters above. */
  
  rho1 = 1.225/RATIO;
  rho2 = 1.225;
  mu1 = 1.0/RE/MURATIO;
  mu2 = 1.0/RE; //using wavelength as length scale
  //f.sigma = 1./(BO*sq(k_));
  G.y = -g_;

  /**
  When we use adaptive refinement, we start with a coarse mesh which
  will be refined as required when initialising the wave. */
  
#if TREE  
  N = 32;
#else
  N = 1 << LEVEL;
#endif
  run();
}

/**
## Initial conditions

These functions return the shape of a third-order Stokes wave with the
wavenumber and steepness given by the parameters above ($ak$ and
$_k_$). */

double wave (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3 - y;
}

double eta (double x, double y) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + ak*eta2 + sq(ak)*eta3;
}

double u_x (double x, double y)
{
  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*cosh(k_*(y + h_))/cosh(k_*h_)*k_*cos(k_*x) +
    ak*3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    cosh(2.0*k_*(y + h_))*2.*k_*cos(2.0*k_*x)/cosh(2.0*k_*h_) +
    ak*ak*1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*3.*k_*cos(3.*k_*x);
}

double u_y (double x, double y)
{
  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*k_*sinh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x) +
    ak*3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    2.*k_*sinh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_) +
    ak*ak*1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    3.*k_*sinh(3.*k_*(y + h_))/cosh(3.*k_*h_)*a_*a_*k_*k_*A_*sin(3.*k_*x);
}
/**
We also calculate an approximation to a Dirac distribution on the wave
surface.  This allows us to calculate a vortex sheet on the surface to
provide a boundary layer in the air above the water surface. */

double gaus (double y, double yc, double T){
  double deltaw = sqrt(2.0/RE)/k_;
  double deltaa = sqrt(2.0/RE*MURATIO/RATIO)/k_;
  double r = y - yc;
  return 2.0/(sqrt(2.0*pi*sq(deltaa)) + sqrt(2.0*pi*sq(deltaw))) *
    (T*exp(-sq(r)/(2.0*sq(deltaw))) + (1.0 - T)*exp(-sq(r)/(2.0*sq(deltaa))));
}

/**
We either restart (if a "restart" file exists), or initialise the wave
using the third-order Stokes wave solution. */

event init (i = 0)
{
  if (!restore ("restart")) {
    do {
      fraction (f, wave(x,y));
    
      /**
      To initialise the velocity field, we first define the potential. */
        foreach(){
        u.x[] = f[] > 1e-3 ? u_x(x,y)*f[] : 0;
        u.y[] = f[] > 1e-3 ? u_y(x,y)*f[] : 0;
        }
 
      boundary ((scalar *){u});
    }

    /**
    On trees, we repeat this initialisation until mesh adaptation does
    not refine the mesh anymore. */

#if TREE  
    while (adapt_wavelet ({f,u},
			  (double[]){1e-8,uemax,uemax,uemax}, LEVEL, 5).nf);
#else
    while (0);
#endif
  }
}

/**
## Outputs

We are interested in the viscous dissipation rate in both water and air. */

int dissipation_rate (double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			     sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			     sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
    rateWater += mu1/rho[]*f[]*sqterm; //water
    rateAir   += mu2/rho[]*(1. - f[])*sqterm; //air
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}

/**
We log the evolution of the kinetic and potential energies and
dissipation rate as functions of the non-dimensional time. */

event graphs (t += T_/40.0) {
  static FILE * fpwater = fopen("budgetWater.dat", "w");
  static FILE * fpair = fopen("budgetAir.dat", "w");
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:ke) reduction(+:gpe) 
	  reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += rho[]*norm2*f[]*dv();
    keAir += rho[]*norm2*(1.0-f[])*dv();
    gpe += rho1*g_*y*f[]*dv();
    gpeAir += rho2*g_*y*(1.0-f[])*dv();
  }
  double rates[2];
  dissipation_rate(rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
    if (i == 0) {
    fprintf (fpwater, "t ke gpe dissipation\n");
    fprintf (fpair, "t ke gpe dissipation\n");
    }
  double gpe0 = -(L0_*rho1*g_*0.5*sq(h_));
   double gpeAir0 = (L0_*rho2*g_*0.5*sq(h_));
  fprintf (fpwater, "%g %g %g %g\n",
	   t/T_, ke/2., gpe - gpe0, dissWater);
  fprintf (fpair, "%g %g %g %g\n",
	   t/T_, keAir/2., gpeAir - gpeAir0, dissAir);
  fprintf (ferr, "%g %g %g %g %g\n",
	   t/T_, ke/2. + gpe - gpe0, ke/2., gpe - gpe0, dissWater);
}

/**
~~~gnuplot Evolution of kinetic, potential and total energies and dissipation rate.
set key auto col
plot 'log' u 1:2 w l, '' u 1:3 w l, '' u 1:(($2+$3)/2.) t 'total/2' w l, \
'' u 1:4 w l
~~~

## Visualisation

We use Basilisk view (and output_ppm()) to display animations of the
results. */
event movies (t += T_/40.0) {

  /**
  We first do simple movies of the volume fraction, level of
  refinement fields. In 3D, these are in a $z=0$ cross-section. */

  {
    static FILE * fp = popen ("ppm2mp4 f.mp4", "w");
    output_ppm (f, fp, min = 0, max = 1, n = 512);
  }

#if TREE
  {
    scalar l[];
    foreach()
    l[] = level;
    static FILE * fp = popen ("ppm2mp4 level.mp4", "w");
    output_ppm (l, fp, min = 5, max = LEVEL, n = 512);
  }
#endif

  /**
  ![Wave breaking. Animation of the level of refinement.](wave/level.mp4)

  We use Basilisk view differently in 2D and 3D. */
  
  scalar omega[];
  vorticity (u, omega);
#if dimension == 2
  view (width = 800, height = 600, fov = 18.8);
  clear();

  /**
  We repeat the drawing periodically in the x-direction. */
  
  for (double x = -L0; x <= L0; x += L0)
    translate (x) {
      draw_vof ("f");
      squares ("omega", linear = true);
    }

  /**
  This gives the following movie.
  ![](wave/movie.mp4)
  */

#else // dimension == 3
  /**
  In 3D, we generate a first movie seen from below. */
  
  view (width = 1600, height = 1200, theta = pi/4, phi = -pi/6, fov = 20);
  clear();
  for (double x = -2*L0; x <= L0; x += L0)
    translate (x) {
      squares ("omega", linear = true, n = {0,0,1}, alpha = -L0/2);
      for (double z = -3*L0; z <= L0; z += L0)
	translate (z = z)
	  draw_vof ("f");
    }
  save ("below.mp4");

  /**
  And a second movie, seen from above. */
  
  view (width = 1600, height = 1200, theta = pi/4, phi = pi/6, fov = 20);
  clear();

  /**
  In 3D, we are doubly-periodic (along x and z). */
  
  for (double x = -2*L0; x <= L0; x += L0)
    translate (x) {
      squares ("omega", linear = true, n = {0,0,1}, alpha = -L0/2);
      for (double z = -3*L0; z <= L0; z += L0)
	translate (z = z)
	  draw_vof ("f");
    }
#endif // dimension == 3
  save ("movie.mp4");
}


/**
Write to hypretree format. */
event write2htg (t += T_/10.0) {
  scalar vort[], p_rgh[];
  vorticity (u, vort);

  vector U[]; // to slove the compatibility issue of restore/dump
  foreach()
    foreach_dimension()
      U.x[] = u.x[];

  char fname[50];
  sprintf(fname, "result.%06d", i);
  scalar* output_scalars = {f, p, vort};
  output_htg(output_scalars,(vector *){U}, ".", fname, i, t);

  fprintf (stderr, "write output to %s\n", fname);
}

/**
## Dump/restore

To be able to restart, we dump the entire simulation at regular
intervals. */

event snapshot (i += 2000) {
  dump();
}

/**
## End 

The wave period is `k_/sqrt(g_*k_)`. We want to run up to 2
(alternatively 4) periods. */

event end (t = 100.0*T_) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}

/**
## Mesh adaptation

On trees, we adapt the mesh according to the error on volume fraction
and velocity. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-8,uemax,uemax,uemax}, LEVEL, 5);
}
#endif

/**
## Running in parallel

This file will work in 2D or 3D, either with parallel multigrid
(without adaptivity), using for example:

~~~bash
qcc -source -D_MPI=1 -grid=multigrid3D wave.c
scp _wave.c occigen.cines.fr:
~~~

and then following a recipe similar to that of the
[atomisation](/src/examples/atomisation.c#on-occigen) example.

To use adaptivity, just do something like:

~~~bash
qcc -source -D_MPI=1 -grid=octree wave.c
scp _wave.c occigen.cines.fr:
~~~
*/
