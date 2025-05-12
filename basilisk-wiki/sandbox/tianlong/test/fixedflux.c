/**
The bubble evaporation and condensation with a fixed flux. This file is copied from [Edorado's sandbox](http://basilisk.fr/sandbox/ecipriano/README) and modified with EBIT.
*/
#define SEMUPC 1
#define ADV_SCHEME 2
#define USE_FIXED_MASSFLUX 2 //1:evaporation 2:condensation
#define INTGRAD_3rd 1
#define USE_DOUBLE_VEL 1
#define TEST "data/"

#if USE_DOUBLE_VEL
#include "double-evaporation.h"
#else
#include "centered-evaporation.h"
#endif //USE_DOUBLE_VEL
#include "semushin_two-phase.h"
#include "semushin-phase-change.h"
#include "mytension.h"

//out flow for the right
u.n[right] = neumann (0);
p[right] = dirichlet (0);
pf[right] = dirichlet (0);

//out flow for the top
u.n[top] = neumann (0);
p[top] = dirichlet (0);
pf[top] = dirichlet (0);

//out flow for the bottom
u.n[bottom] = neumann (0);
p[bottom] = dirichlet (0);
pf[bottom] = dirichlet (0);

//out flow for the left
u.n[left] = neumann (0);
p[left] = dirichlet (0);
pf[left] = dirichlet (0);


//use symmetry for the top and bottom (imposed in the solid setup)

int maxlevel = 6, minlevel = 3;
double betaGrowth;
const double tshift = 0.0;
#if (USE_FIXED_MASSFLUX == 1)
const double R0 = 1e-3;
const double mdot_fixed = -0.1;
#else
const double R0 = 2e-3;
const double mdot_fixed = 0.1;
#endif
const double t_stop = 0.01;

const double femax = 1e-6;
const int num_refine = 2;

int main (int argc, char * argv[]) {
  if (pid() == 0)
  {
    if (argc == 2)
    {
      maxlevel = atoi(argv[1]);
    }
    }

#if _MPI
  MPI_Bcast(&maxlevel, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  rho1 = 1000.0; rho2 = 1.0;
  mu1 = 1e-3; mu2 = 1.78e-5;

  //we still set here, though they are not needed in this case
  lambda1 = 0.679, lambda2 = 0.025;
  cp1 = 4216., cp2 = 2030.;
  dhev = 2.26e+6;

  /**
  The initial temperature and the interface
  temperature are set to the same value. */

  Tsat = 373.15;

  /**
  We change the dimension of the domain
  and the surface tension coefficient. */
  double dtlist[] = {8e-5, 4e-5, 1.4e-5, 5e-6};
  DT = dtlist[maxlevel - 5];
  TOLERANCE = 1E-5;
  L0 = 8e-3;
  f.sigma = 0.059;
  /**
  We define a list with the maximum time
  steps and the maximum levels of refinement. */
  init_grid (1 << maxlevel);
  origin(-L0/2.0, -L0/2.0);
  run();
}

ADV_SCHEME
event defaults(i = 0)
{

}

/**
We initialize the volume fraction field and the temperature
in the gas and in liquid phase. */

event init (i = 0) {
  vertex scalar phi[];
  foreach_vertex()
  {
    phi[] = sqrt(x * x + y * y) - R0;
  }

  init_markers(phi);

  foreach ()
  {
    TL[] = Tsat;
    TG[] = Tsat;
    T[] = f[] > 0.5 ? TL[] : TG[];
    foreach_dimension()
    {
      u.x[] = 0.0;
    }
  }

  boundary({T, TL, TG});

  getColorExact(color_cc);
  getMdot(color_cc, mdot);
}


void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}

event outLog(i++)
{
  if(pid() == 0 && i % 10 == 0)
  {
    printf("i = %d t = %g dt = %g\n", i, t, dt);
  }
}

/**
We refine the interface and the region where the
temperature field changes. */

#if TREE
event adapt (i++) {
  //correct the strange change of vertex color
  foreach_vertex()
  {
    double vofs = f[] + f[0, -1] + f[-1] + f[-1, -1];
    if (vofs == 0.0)
      color_pha[] = 0.0;
    else if (vofs == 4.0)
      color_pha[] = 1.0;
  }
}
#endif


double exact (double time) {
  return R0 - mdot_fixed / rho2 * t;
}

/**
### Output Files

We write the thickness of the vapor layer and the analytic
solution on a file.
*/

event movie(t += t_stop / 20.0)
{
  scalar fg[];
  foreach ()
    fg[] = 1. - f[];

  double effective_radius = pow(statsf(fg).sum / pi, 1. / 2.);

  double rsol = exact(t + tshift);
  double relerr = (rsol > 0.) ? fabs(rsol - effective_radius) / rsol : 0.;

  char name[80];
  sprintf(name, "data/OutputData-%d", maxlevel);

  static FILE *fp = fopen(name, "w");
  fprintf(fp, "%g %g %g %g\n",
          t + tshift, effective_radius, exact(t + tshift), relerr);
  fflush(fp);
}

event finalEvent(t = t_stop)
{
  char name[80];
  sprintf (name, "data/infsemu-%d", maxlevel);
  output_facets_semushin(name);

  sprintf(name, "data/dataend-%d.gfs", maxlevel);
  FILE *file = fopen(name, "w");
  output_gfs(file);
  // do nothing
}
