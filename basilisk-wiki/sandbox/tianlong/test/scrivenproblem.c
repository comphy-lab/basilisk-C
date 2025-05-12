/**
The Scriven problem setup. This file is copied from [Edorado's sandbox](http://basilisk.fr/sandbox/ecipriano/README) and modified with EBIT.
*/
#define SEMUPC 1
#define INTGRAD_3rd 1
#define USE_DOUBLE_VEL 0
#define ADV_SCHEME 2 
#define TEST "data/"

#include "axi.h"
#if USE_DOUBLE_VEL
#include "double-evaporation.h"
#else
#include "centered-evaporation.h"
#endif //USE_DOUBLE_VEL
#include "semushin_two-phase.h"
#include "semushin-phase-change.h"
#include "mytension.h"
#include "view.h"

#include <gsl/gsl_integration.h>

//out flow for the right
u.n[right] = neumann (0);
p[right] = dirichlet (0);
pf[right] = dirichlet (0);

//out flow for the top
u.n[top] = neumann (0);
p[top] = dirichlet (0);
pf[top] = dirichlet (0);

//wall for the left
//wall for the bottom
uf.n[left] = 0.0;
uf.n[bottom] = 0.0;


//use symmetry for the top and bottom (imposed in the solid setup)

int maxlevel = 6, minlevel = 6;
double Tbulk;
double betaGrowth;
const double tshift = 94.7e-6;
const double R0 = 50.0e-6;
const double t_stop = 3.0 * tshift;

const double femax = 1e-6;
const int num_refine = 2;

double intfun (double x, void * params) {
  double beta = *(double *) params;
  return exp(-sq(beta)*(pow(1. - x, -2.) - 2.*(1. - rho2/rho1)*x - 1 ));
}

double tempsol (double r, double R) {
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  double result, error;
  double beta = betaGrowth;
  gsl_function F;
  F.function = &intfun;
  F.params = &beta;
  gsl_integration_qags (&F, 1.-R/r, 1., 1.e-9, 1.e-5, 1000,
                        w, &result, &error);
  gsl_integration_workspace_free (w);
  return Tbulk - 2.*sq(beta)*(rho2*(dhev + (cp1 - cp2)*(Tbulk - Tsat))/rho1/cp1)*result;
}


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
  rho1 = 958.4; rho2 = 0.597;
  mu1 = 2.80e-4; mu2 = 1.26e-6;
  lambda1 = 0.679, lambda2 = 0.025;
  cp1 = 4216., cp2 = 2030.;
  dhev = 2.26e+6;

  /**
  The initial temperature and the interface
  temperature are set to the same value. */

  Tbulk = 375.15, Tsat = 373.15;

  /**
  We change the dimension of the domain
  and the surface tension coefficient. */

  L0 = 160e-6;
  f.sigma = 0.059;
  /**
  We define a list with the maximum time
  steps and the maximum levels of refinement. */
  TOLERANCE = 1.e-7;
  init_grid (1 << maxlevel);
  run();
}


event defaults(i = 0)
{

}

/**
We initialize the volume fraction field and the temperature
in the gas and in liquid phase. */

event init (i = 0) {

  double alpha = lambda1/rho1/cp1;
  betaGrowth = R0 / 2.0 / sqrt(alpha * tshift);

  vertex scalar phi[];
  foreach_vertex()
  {
    phi[] = sqrt(x * x + y * y) - R0;
  }

  init_markers(phi);

  foreach ()
  {
    double r = sqrt(x * x + y * y);
    TL[] = r < R0 ? Tsat : tempsol(r, R0);
    TG[] = Tsat;
    T[] = f[] > 0.5 ? TL[] : TG[];
    foreach_dimension()
    {
      u.x[] = 0.0;
    }
  }

  boundary({T, TL, TG});

  getColorExact(color_cc);
  //for the first step computation
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
  if(pid() == 0 && i % 100 == 0)
  {
    printf("i = %d t = %g dt = %g\n", i, t, dt);
  }
}

/**
We refine the interface and the region where the
temperature field changes. */

#if TREE
event adapt (i++) {
  scalar solid_refine[];
  scalar solid_refiney[];
  scalar f_refine[];

  fillRefineVOFs(f_refine, solid_refine, solid_refiney);

  adapt_wavelet({T, f_refine, solid_refine, solid_refiney}, (double[]){1.e-3, femax, femax, femax},
                maxlevel, minlevel);

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
  return 2.*betaGrowth*sqrt(lambda1/rho1/cp1*time);
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

  double effective_radius = pow(3. * statsf(fg).sum, 1. / 3.);

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
  sprintf(name, "data/Temperature-%d", maxlevel);
  Array *arrtemp = array_new();
  for (double x = 0.0; x < L0; x += L0 / (1 << maxlevel))
  {
    double val = interpolate(T, x, 0.);
    val = (val == nodata) ? 0. : val;
    array_append(arrtemp, &val, sizeof(double));
  }
  double *temps = (double *)arrtemp->p;

#if _MPI
  int size = arrtemp->len / sizeof(double);
  MPI_Allreduce(MPI_IN_PLACE, temps, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (pid() == 0)
  {
    FILE *fpp = fopen(name, "w");
    int count = 0;
    for (double x = 0; x < L0; x += L0 / (1 << maxlevel))
    {
      double r = x;
      double R = exact(t + tshift);
      double tempexact = (r >= R) ? tempsol(r, R) : Tsat;
      fprintf(fpp, "%g %g %g\n", x, temps[count], tempexact);
      count++;
    }
    fflush(fpp);
    fclose(fpp);
  }
  array_free(arrtemp);

  sprintf(name, "data/dumpend-%d", maxlevel);
  dump(name);

  sprintf(name, "data/dataend-%d.gfs", maxlevel);
  FILE *file = fopen(name, "w");
  output_gfs(file);
  // do nothing
}
