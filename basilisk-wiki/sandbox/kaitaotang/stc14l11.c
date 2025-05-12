/**
# Hole expansion on a spherical shell

We use the centered Navier--Stokes solver and log performance
statistics. */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"

/**
We have two phases e.g. air and water. For large viscosity and density
ratios, the harmonic mean for the viscosity tends to work better than
the default arithmetic mean. We "overload" the default by defining the
*mu()* macro before including the code for two phases. */

// #define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"

/**
We also need surface tension, and in 3D only we will use the
$\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) to display the vortices using
Basilisk View. */

#include "tension.h"
#if dimension == 3
# include "lambda2.h"
#endif
#include "navier-stokes/conserving.h"
#include "view.h"
#include "tag.h"

/**
We can control the maximum runtime. */

#include "maxruntime.h"



/**
The density ratio is 833 and the dynamic viscosity ratio 55. */

#define RHOR 833.
#define MUR 55.
#define eR 0.004
#define re 1.5
#define OH 0.01

#define MAXTIME 300.


/**
We choose as length unit the diameter of the bubble. The domain is
$10^3$. *X0* is the initial position of the bubble relative to the
front wall. R0 is the radius of the droplet. 
U is the speed of the ambient flow. Gravity is neglected. */

#define WIDTH 3
#define R0 1

int MAXLEV = 13, MINLEV = 5;

FILE * fk = NULL;
FILE * fe = NULL;

double diss_w = 0., diss_a = 0;
double e0 = 0., rh = 0.;
double xb = 0., yb = 0., zb = 0., sb = 0.;
double vbx = 0., vby = 0., vbz = 0.;

// Boundary conditions
u.n[left] = dirichlet(0.);
u.n[right] = dirichlet(0.);

u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);

u.n[front] = dirichlet(0.);
u.n[back] = dirichlet(0.);

scalar l2[], omegay[], pp[];



/**
The main function can take two optional parameters: the maximum level
of adaptive refinement (as well as an optional maximum runtime). */

int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    MAXLEV = atoi (argv[1]);
  
  /**
  We set the domain geometry and initial refinement. */
  size (WIDTH);
  origin (-L0/2., -L0/2., 0.);
  init_grid (2 << MINLEV);

  /**
  We set the physical parameters: densities, viscosities and surface
  tension. Note that the viscosity uses the gas density for the density
  scale. */
  
  rho1 = 1.;
  rho2 = rho1/RHOR;
  e0 = R0 * eR;
  rh = e0 * re;

  f.sigma = 1.6e-4;
  mu1 = OH * sqrt(rho1*f.sigma*e0);
  mu2 = mu1/MUR;

  /**
  We reduce the tolerance on the divergence of the flow. This is
  important to minimise mass conservation errors for these simulations
  which are very long. */
  
  TOLERANCE = 1e-5;

  char name[80];
  sprintf (name, "3Denergy-%d.dat", MAXLEV);
  fe = fopen (name, "w");

  sprintf (name, "3Dkinetics-%d.dat", MAXLEV);
  fk = fopen (name, "w");

  run();

  fclose (fe); fclose (fk);
}

/**
For the initial conditions, we first try to restore the simulation
from a previous "restart", if this fails we refine the mesh locally to
the maximum level, in a sphere of diameter 1.5 around the bubble. We
then initialise the volume fraction for a bubble initially at (0,0,0)
of diameter unity. In the fraction command, the signs are chosen to
ensure that a droplet, not a bubble is formed (i.e. f=1 inside, f=0
outside).*/

event init (t = 0) {
  if (!restore (file = "restart")) {
	refine (sq(x) + sq(y) + sq(z) - sq(R0+e0) < 0 && sq(x) + sq(y) + sq(z) - sq(R0-e0) > 0 && level < MAXLEV);
	fraction(f, (R0 + 0.5*e0 - sqrt(sq(x) + sq(y) + sq(z))) * (-R0 + 0.5*e0 + sqrt(sq(x) + sq(y) + sq(z))));
	
	boundary({f});	
  }
  else fprintf (stderr, "Dumpfile restored!\n");
}

event poke (t = 0.02) {
  foreach() if ((z > 0) && (rh - sqrt(sq(x) + sq(y)) >= 0)) f[] = 0.;
}

/**
We adapt the mesh by controlling the error on the volume fraction and
velocity field. */

event adapt (i++) {
  double femax = 1e-2, uemax = 1e-2;
  vector gf[];
  gradients ({f}, {gf});

  adapt_wavelet ({gf.x, gf.y, gf.z, u.x, u.y, u.z}, (double[]){femax, femax, femax, uemax, uemax, uemax}, MAXLEV, MINLEV);
}

/**
## Outputs

Every ten timesteps, we output the time, volume, position, and
velocity of the bubble. */

int dissipation_rate (vector u, double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1] - u.x[-1])/(2.0*Delta);
    double dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]   - u.y[-1]  )/(2.*Delta);
    double dvdy = (u.y[0,1] - u.y[0,-1])/(2.0*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.0*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.0*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			      sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			      sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz));
    
    rateWater += mu1*f[]*sqterm; //water
    rateAir   += mu2*(1. - f[])*sqterm; //air
  }

  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}



event logfile (t = 0; t <= MAXTIME; t += 0.01) {
  xb = 0., yb = 0., zb = 0., sb = 0.;
  vbx = 0., vby = 0., vbz = 0.;
  double ke_d = 0., ke_a = 0., a = 0.;

  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb) 
      reduction(+:vbx) reduction(+:vby) reduction(+:vbz) 
      reduction(+:sb) reduction(+:ke_d) reduction(+:ke_a) reduction(+:a)) {
    double dv_l = f[]*dv();
    xb += x*dv_l;
    yb += y*dv_l;
    zb += z*dv_l;
    vbx += u.x[]*dv_l;
    vby += u.y[]*dv_l;
    vbz += u.z[]*dv_l;
    sb += dv_l;

    // Kinetic energy of the droplet and ambient airflow
    ke_d += 0.5*dv()*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])) * rho1*f[];
    ke_a += 0.5*dv()*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])) * rho2*(1-f[]);
  }

  // Droplet surface area
  a = interface_area(f);



  // Viscous dissipation rate
  double rates[2];
  dissipation_rate(u, rates);
  diss_w += rates[0], diss_a += rates[1];
  
  // Droplet CM's location and velocity
  xb = xb/sb; yb = yb/sb; zb = zb/sb; vbx = vbx/sb; vby = vby/sb; vbz = vbz/sb;

  double h_max = -100.0, h_min = 100.0, t_max = -100.0, t_min = 100.0;
  foreach(reduction(max:h_max) reduction(min:h_min) reduction(max:t_max) reduction(min:t_min))
    if (f[] > 1e-6 && f[] < 1.-1e-6) {
      h_max = fmax(h_max, y); h_min = fmin(h_min, y); 
      t_max = fmax(t_max, x); t_min = fmin(t_min, x);
    }

  fprintf (stderr,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb, a, xb, yb, zb,
	   vbx, vby, vbz, h_max-h_min, t_max-t_min);
  fflush (stderr);

  fprintf (fk,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb, a, xb, yb, zb, vbx, vby, vbz, h_max-h_min, t_max-t_min);
  fflush (fk);

  fprintf (fe, 
          "%.8f %.8f %.8f %.8f %.8f %.8f\n", 
          t, ke_d, ke_a, diss_w, diss_a, a*f.sigma);
  fflush (fe);
}

/**
Every time unit, we output a full snapshot of the simulation, to be
able to restart and for visualisation. In three dimensions, we compute
the value of the $\lambda_2$ field which will be used for
visualisation of vortices, as well as the streamwise vorticity
$\omega_y = \partial_x u_z - \partial_z u_x$. */

event snapshot (t = 0; t <= MAXTIME; t += 0.2) {  
#if dimension == 3
  foreach() {
    omegay[] = (u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
    pp[] = p[];
  }
  boundary ({omegay, pp});
#endif
  
  char name[80];
  sprintf (name, "dump-%04g", t);
  dump (file = name);
}

event outputInterface(t = 0; t <= MAXTIME; t += 0.05) {
  char resultname[100];
  sprintf(resultname, "./interface_profiles/results%4.2f_%d.dat", t, pid());
  FILE * fp = fopen(resultname, "w");

  scalar xpos[], ypos[], zpos[];
  position (f, xpos, {1, 0, 0});
  position (f, ypos, {0, 1, 0});
  position (f, zpos, {0, 0, 1});
  foreach()
      if (xpos[] != nodata)
        fprintf (fp, "%g %g %g\n", xpos[], ypos[], zpos[]);
  
  fclose (fp);
}



event count_droplets (t = 0; t <= MAXTIME; t += 0.01) {
  scalar m_cd[];
  foreach() m_cd[] = f[] > 0.5;
  int n = tag (m_cd);

  /**
  Once each cell is tagged with a unique droplet index, we can easily
  compute the volume *v* and position *b* of each droplet. Note that
  we use *foreach (serial)* to avoid doing a parallel traversal when
  using OpenMP. This is because we don't have reduction operations for
  the *v* and *b* arrays (yet). */

  double v[n], ux_drop[n], uy_drop[n], uz_drop[n];
  coord b[n];
  vector u_drop[n];

  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = ux_drop[j] = uy_drop[j] = uz_drop[j] = 0.;

  foreach_leaf()
    if (m_cd[] > 0) {
      int j = m_cd[] - 1;
      v[j] += dv()*f[];

      coord p = {x,y,z};
      foreach_dimension()
	    b[j].x += dv()*f[]*p.x;

      ux_drop[j] += dv()*f[]*u.x[];
	  uy_drop[j] += dv()*f[]*u.y[];
	  uz_drop[j] += dv()*f[]*u.z[];
    }

 /**
 When using MPI we need to perform a global reduction to get the
 volumes and positions of droplets which span multiple processes. */

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, ux_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uy_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uz_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /**
  Finally we output the volume and position of each droplet to
  standard output. */

  static FILE * fdrop = fopen("droplets.dat", "w");

  if (n > 1) 
    for (int j=0; j<n; j++) {
      fprintf (fdrop, "%d %g %d %g %g %g %g %g %g %g\n", i, t, j, v[j], b[j].x/v[j], b[j].y/v[j], b[j].z/v[j], ux_drop[j]/v[j], uy_drop[j]/v[j], uz_drop[j]/v[j]);
    }
}



#define POPEN(name, mode) fopen (name ".ppm", mode)

event movie (t += 0.01)
{
  view (quat = {0.315, 0.080, 0.315, 0.891}, fov = 25, near = 0.01, far = 1000, tx = -0.025, ty = 0.017, tz = -2.549, width = 720, height = 720, bg = {1,1,1}, samples = 4);
  clear();
  draw_vof ("f");

  static FILE * fp = POPEN ("movie", "w");
  save (fp = fp);
}