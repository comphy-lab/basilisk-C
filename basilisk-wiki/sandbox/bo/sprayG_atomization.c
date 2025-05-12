/**
## Atomisation of a gasoline surrogate jet under non-evaporative ECN spray G operating conditions. 

The Spray G conditions coorespond to the the non-reacting spray-guided gasoline injection, see [ECN Spray G](https://ecn.sandia.gov/gasoline-spray-combustion/target-condition/spray-g-operating-condition/). 

Here we follow the Argone National Lab experiment and consider a non-evaporative condtion. Therefore, the two-phase Navier--Stokes equations with surface tension is solved. We need the *tag()* function to collect droplets statistics. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "fractions.h"
#include "lambda2.h"
#include "navier-stokes/conserving.h"

/**
We define the radius and length of the inner-hole (equal to the jet radius at the inlet) as r1 and x1, the radius and length of the counterbore as r2 and x2. The values are from the ECN website. The domain size is 32 times of jet diameter at the inlet. */

#define r1 0.5
#define r2 1.123
#define x1 0.878
#define x2 2.283
#define length 0.878
#define box_length 32.
#define Uinj 1.
#define norm(v) (sqrt(sq(v.x[]) + sq(v.y[]) + sq(v.z[])))
#define magnitude(v) (sqrt(sq(v.x) + sq(v.y) + sq(v.z)))
#define xprobe2mm 11.56
#define epsilon 0.03125

/**
The default maximum refinement level is 11 and the errors threshold
for mesh adpatation for both velocity and volume fraction are 0.01. We remove droplets with a volume smaller than vol_cut. Utan is 
the tangential component of the injection velocity, which is introduced to mimic the effect of internal flow on atomization.*/

int maxlevel = 11;
double uemax = 0.01;
double femax = 0.01;
double maxruntime = 661420;
double vol_cut = 3.0518e-5;
const double PI = 3.1415926;
double Utan = 0.5;

/**
To impose boundary conditions on a short cyliner we use an auxilliary volume
fraction field *f0* which is one inside the cylinder and zero outside.*/

scalar f0[];

/** 
The scalar field nozzle[ ] is used to define the embeded solid nozzle.*/

scalar nozzle[];

/**
We invoke Dirichlet boundary condtion for velocity on the left surface and set an outflow boundary condition for the right surface.*/

u.n[left]  = dirichlet(f0[]*Uinj);
u.t[left]  = dirichlet(f0[]*Utan/1.4142);
u.r[left]  = dirichlet(f0[]*Utan/1.4142);
p[left]    = neumann(0);
f[left]    = f0[];

u.n[right] = neumann(0);
p[right]   = dirichlet(0);

/**
The program can take five optional command-line arguments: the maximum
level, the maximum run time, the tangential injection velocity,
the error threshold on velocity, and the error throshold of fraction. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    maxruntime = atof (argv[2]); 
  if (argc > 3)
    Utan = atof (argv[3]);
  if (argc > 4)
    uemax = atof (argv[4]);
  if (argc > 5)
    femax = atof (argv[5]);

  /**
  The initial domain is discretised with $64^3$ grid points. The origin is shifted to the center of the left surface. We also set the domain size. */
  
  init_grid (64);
  origin (0, -box_length/2, -box_length/2);
  size (box_length);

  /**
  We set the density and viscosity of each phase as well as the
  surface tension coefficient and start the simulation. */
  
  rho1 = 1., rho2 = 4.296e-3;
  mu1 = 7.478e-5, mu2 = 1.373e-6;
  f.sigma = 2.424e-5;

  run();
}

/**
## Initial conditions */

event init (t = 0) {
  if (!restore (file = "dump")) {

    /**
    We use a static refinement down to *maxlevel* in a cylinder 1.2 times 
    longer than the initial jet and twice the square root of radius. */

    refine (x < 1.2*length && sq(y) + sq(z) < 2.*sq(r1) && level < maxlevel);
    
    /**
    We initialise the auxilliary volume fraction field for a cylinder of
    constant radius. */
    
    fraction (f0, sq(r1) - sq(y) - sq(z));
    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0}); 
    
    /**
    We then use this to define the initial jet and its velocity. */

    foreach() {
      f[] = f0[]*(x < length);
      u.x[] = f[]*Uinj;
    }
    boundary ({f,u.x});
  }
}

/**
## Define the fraction field nozzle[] to hold the solid part of the
   simulation including the inner-hole and counterbore. */

event static_nozzle (i++) {
  fraction (nozzle, x<x1 ? sq(y)+sq(z)-sq(r1) : x<(x1+x2) ? sq(y)+sq(z)-sq(r2) : -1);

  foreach()
    foreach_dimension()
      u.x[] = (1.-nozzle[])*u.x[];
  boundary ((scalar *){u});
}

/**
## Outputs
We log some statistics on the solver. */

event logfile (i+=10) {
  double voljet = 0., umjet = 0., areajet = 0.;
  double xjet=0., xjetmax = -1e100;

  foreach(reduction(+:voljet) reduction(+:umjet) reduction(+:areajet) reduction(max:xjetmax)) {
    double dv =  f[]*dv();
    voljet += dv;
    umjet  += dv*u.x[];

    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = mycs (point, f), p;
      double alpha = plane_alpha (f[], n);
      areajet += pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    }
    
    xjet = x*f[];
    if (xjet > xjetmax)
      xjetmax = xjet;
  }

  if (i == 0)
    fprintf (ferr,
	     "i t dt vol_jet area_jet um_jet xjetmax mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  fprintf (ferr, "%d %g %g %g %g %g %g %d %d %d %ld %g %g\n",  
	   i, t, dt, voljet, areajet, umjet/voljet, xjetmax, mgp.i, mgpf.i, mgu.i, grid->tn, perf.t, perf.speed);
}

/**
 Output the dump files for postprocessing and restarting the simulation. */

event output_dump (t = 0.; t += 0.1; t <= 64.) {
  char name[80];
  sprintf (name, "dump-%3.1f",t);
  dump (name);
}

/**
  Remove small droplets every 10 steps. */

event remove_drops ( i+=10 )
{
  scalar m[];
  foreach()
    m[] = f[] > 0.;
  int n = tag (m);
  double v[n];
  int remove_flag[n];
  for (int j = 0; j < n; j++) {
    v[j] = 0.;
    remove_flag[j] = 0;
  }
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
    }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  for (int j = 0; j < n; j++)
    if ( v[j] < vol_cut )
       remove_flag[j] = 1;

  foreach()
    if (m[] > 0) {
      int j = m[] - 1;
      if ( remove_flag[j] == 1 )
         f[]=0.;
    }
}

/**
  Mesh adaptation
  We adapt the mesh according to the error on the volume fraction field
and the velocity. */

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, maxlevel);
  unrefine ( x > box_length*0.9 && level > 5);
}

/**
  Define the maxruntime so that we could restart from the specified time point. */

event runtime (i += 10) {
  mpi_all_reduce (perf.t, MPI_DOUBLE, MPI_MAX);
  if (perf.t/60 >= maxruntime) {
    dump (file = "dump"); 
    return 1; 
  }
}
