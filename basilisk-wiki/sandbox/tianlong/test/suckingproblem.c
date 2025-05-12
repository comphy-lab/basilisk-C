/**
The Sucking problem setup. This file is copied from [Edorado's sandbox](http://basilisk.fr/sandbox/ecipriano/README) and modified with EBIT.
*/
#define SEMUPC 1
#define IS_1D 1
#define INTGRAD_3rd 1
#define USE_DOUBLE_VEL 1
#define ADV_SCHEME 2
#define TEST "data/"

#if USE_DOUBLE_VEL
#include "double-evaporation.h"
#else
#include "centered-evaporation.h"
#endif //USE_DOUBLE_VEL
#include "semushin_two-phase.h"
#include "semushin-phase-change.h"
#include "mytension.h"

#if USE_MY_SOLID
#include "mysolid.h"
#endif

#include "fsolve.h"


//out flow for the right
u.n[right] = neumann (0);
p[right] = dirichlet (0);
pf[right] = dirichlet (0);

//wall for the left
u.t[left] = dirichlet (0);
uf.n[left] = 0.0;
TG[left] = dirichlet (Tsat);
TL[left] = dirichlet (Tsat);
T[left] = dirichlet (Tsat);

//use symmetry for the top and bottom (imposed in the solid setup)
uf.n[top] = 0.0;
uf.n[bottom] = 0.0;


int maxlevel = 6, minlevel = 3;
double Tbulk;
double tshift, betaGrowth;
const double delta_0 = 0.0022033868116443861;
const double t_stop = 1.0;

const double femax = 1e-6;
const int num_refine = 2;

#if USE_MY_SOLID
void solidSetup()
{
  double delta = L0 / (double)(1 << maxlevel);
  SOLID_LEN_y = ((1 << maxlevel) - 1) * delta;
}
#endif

int betafun (const gsl_vector * x, void * params, gsl_vector * f) {
  double * xdata = x->data;
  double * fdata = f->data;

  double beta = xdata[0];
  double alpha1 = lambda1/rho1/cp1;
  double alpha2 = lambda2/rho2/cp2;

  fdata[0] = exp(sq(beta))*erf(beta)*(beta -
          ( (Tbulk - Tsat)*cp2*lambda1*sqrt (alpha2)*exp
          (-sq(beta)*sq(rho2)*alpha2/sq(rho1)/alpha1) ) /
          (dhev*lambda2*sqrt(pi*alpha1)*
          erfc(beta*rho2*sqrt(alpha2)/rho1/sqrt(alpha1))));

  return GSL_SUCCESS;
}

double tempexact (double x, double beta, double t) {
  double alpha1 = lambda1/rho1/cp1;
  double alpha2 = lambda2/rho2/cp2;

  return Tbulk - ((Tbulk - Tsat)/erfc(beta*rho2*sqrt(alpha2)/rho1/sqrt(alpha1)))
      * erfc (x/2./sqrt(alpha1*t) + beta*(rho2 - rho1)/rho1*sqrt(alpha2/alpha1));
}

#if !USE_MY_SOLID
  scalar is_solid[];
#endif

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
  rho1 = 958.4, rho2 = 0.597;
  mu1 = 2.80e-4, mu2 = 1.26e-5;
  lambda1 = 0.679, lambda2 = 0.025;
  cp1 = 4216., cp2 = 2030.;
  dhev = 2.26e+6;

  /**
  The initial temperature and the interface
  temperature are set to the same value. */

  Tbulk = 378.15, Tsat = 373.15;

  /**
  We change the dimension of the domain
  and the surface tension coefficient. */

  L0 = 10e-3;
  
  /**
  We define a list with the maximum time
  steps and the maximum levels of refinement. */

  double dtlist[] = {0.001, 0.0005, 0.0001};
  DT = dtlist[maxlevel - 6];
  TOLERANCE = 1E-4;
  NITERMAX = 200;
  init_grid (1 << maxlevel);
#if USE_MY_SOLID
  solidSetup();
  origin(-SOLID_LEN_x, -SOLID_LEN_y);
#endif
  run();
}


event defaults(i = 0)
{
#if USE_MY_SOLID
  for (int ib = 0; ib < nboundary; ib++)
  {
    pf.boundary[ib] = p.boundary[ib];
  }
  setSolidFlag();
  foreach_dimension()
  {
    if(IS_SOLID_x)
    {
      f.boundarySolid_x = boundarySolidNeumman_x;
      uf.x.boundarySolid_x = boundarySolidVelF_x;
#if USE_DOUBLE_VEL
      uf2.x.boundarySolid_x = boundarySolidVelF_x;
#endif

      color_pha_cen.boundarySolid_x =  boundarySolidNeumman_x;
      s.x.boundarySolid_x = boundarySolidVectorZero_x;
      s_tmp.x.boundarySolid_x = boundarySolidVectorZero_x;
      with_marker.x.boundarySolid_x = boundarySolidVectorZero_x;
      ss_tmp.x.boundarySolid_x = boundarySolidVectorZero_x;

      color_cc.boundarySolid_x = boundarySolidNeumman_x;
      mdot.boundarySolid_x = boundarySolidNeumman_x;
      phi_dis.boundarySolid_x = boundarySolidNeumman_x;
      T.boundarySolid_x = boundarySolidNeumman_x;
      TL.boundarySolid_x = boundarySolidNeumman_x;
      TG.boundarySolid_x = boundarySolidNeumman_x;
      dTdnL.boundarySolid_x = boundarySolidNeumman_x;
      dTdnG.boundarySolid_x = boundarySolidNeumman_x;
    }
  }
#endif
}

/**
We initialize the volume fraction field and the temperature
in the gas and in liquid phase. */

event init (i = 0) {

#if !USE_MY_SOLID
  foreach()
  {
    is_solid[] = 1.0;

  }

  foreach_boundary(top)
  {
    is_solid[] = 0.0;
  }

  boundary({is_solid});
#endif

  vertex scalar phi[];
  foreach_vertex(){
    phi[] = x - delta_0;
  }
  
  init_markers(phi);

  double effective_height = 0.0;
  foreach(reduction(+:effective_height))
  {
    if(is_solid[] == 0)
    {
      effective_height += (1.0 - f[]) * Delta;
    }
  }


  Array * arrUnk = array_new();
  {
    double betafg = 0.9;
    array_append (arrUnk, &betafg, sizeof(double));
    double * unks = (double *)arrUnk->p;
    fsolve (betafun, arrUnk, NULL);
    betaGrowth = unks[0];
  }
  array_free (arrUnk);

  tshift = rho2*cp2/lambda2*sq (effective_height/2./betaGrowth);

  foreach() {
    TL[] = tempexact (x, betaGrowth, t+tshift);;
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
    printf ("%d %g %g %g %d \n", mg.i, mg.resb, mg.resa,
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
  //correct the strange change of vertex color
  foreach_vertex()
  {
    double vofs = f[] + f[0, -1] + f[-1] + f[-1, -1];
    if (vofs == 0.0)
      color_pha[] = 0.0;
    else if (vofs == 4.0)
      color_pha[] = 1.0;
#if USE_MY_SOLID
    color_pha[] *= (1.0 - is_solid_vertex[]);
#endif
  }
}
#endif


double exact (double time) {
  return 2.*betaGrowth*sqrt(lambda2/rho2/cp2*time);
}

/**
### Output Files

We write the thickness of the vapor layer and the analytic
solution on a file.
*/

event movie (t += t_stop / 20.0) {
  double effective_height = 0.;
  foreach(reduction(+:effective_height))
  {
    if(is_solid[] == 0)
    {
      effective_height += (1.0 - f[]) * Delta;
    }
  }

  double relerr = fabs (exact(t+tshift) - effective_height) / exact(t+tshift);

  char name[80];
  sprintf (name, "data/OutputData-%d", maxlevel);
  static FILE * fp = fopen (name, "w");

  fprintf (fp, "%g %g %g %g\n", t+tshift, effective_height, exact (t+tshift), relerr);
  fflush (fp);
}

event finalEvent(t = t_stop)
{
  double effective_height = 0.;
  foreach(reduction(+:effective_height))
  {
    if(is_solid[] == 0)
    {
      effective_height += (1.0 - f[]) * Delta;
    }
  }

  char name[80];
  sprintf(name, "data/Temperature-%d", maxlevel);
  Array *arrtemp = array_new();
  for (double x = 0.; x < L0; x += 0.5 * L0 / (1 << maxlevel))
  {
    double val = x > effective_height ? interpolate(TL, x, 0.) : interpolate(TG, x, 0.);
    val = (val == nodata) ? 0. : val;
    array_append(arrtemp, &val, sizeof(double));
  }
  double *temps = (double *)arrtemp->p;
#if _MPI
  int size = arrtemp->len / sizeof(double);
  MPI_Allreduce(MPI_IN_PLACE, temps, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  if(pid() == 0)
  {
    FILE *fpp = fopen(name, "w");
    int count = 0;
    for (double x = 0.; x < L0; x += 0.5 * L0 / (1 << maxlevel))
    {
      double R = exact(t + tshift);
      double temp = x > R ? tempexact (x, betaGrowth, t+tshift) : Tsat;
      fprintf(fpp, "%g %g %g\n", x, temps[count], temp);
      count++;
    }
    fflush(fpp);
    fclose(fpp);
  }
  array_free(arrtemp);
  
  sprintf(name, "data/dump-%d", maxlevel);
  dump(name);
  // do nothing
}
