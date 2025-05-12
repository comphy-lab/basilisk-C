/**
The Stefan problem setup. This file is copied from [Edorado's sandbox](http://basilisk.fr/sandbox/ecipriano/README) and modified with EBIT.
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

int maxlevel = 6, minlevel = 3;
double Twall;
double lambdaval = 0.06779249298045148;
double delta0 = 322.5e-6;
double t_stop = 10.0;
double tshift, teff;

const double femax = 1e-6;
const int num_refine = 2;

//out flow for the right
u.n[right] = neumann (0);
p[right] = dirichlet (0);
pf[right] = dirichlet (0);

//wall for the left
u.t[left] = dirichlet (0);
uf.n[left] = 0.0;
TG[left] = dirichlet (Twall);
TL[left] = dirichlet (Twall);
T[left] = dirichlet (Twall);

double tempsol (double time, double x) {
  return Twall + ((Tsat - Twall)/erf(lambdaval))*
    erf(x/2./sqrt(lambda2/rho2/cp2*time));
}

//use symmetry for the top and bottom (imposed in the solid setup)

#if USE_MY_SOLID
void solidSetup()
{
  double delta = L0 / (double)(1 << maxlevel);
  SOLID_LEN_y = ((1 << maxlevel) - 1) * delta;
}
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
  rho1 = 958., rho2 = 0.6;
  mu1 = 2.82e-4, mu2 = 1.23e-5;
  lambda1 = 0.68, lambda2 = 0.025;
  cp1 = 4216., cp2 = 2080.;
  dhev = 2.256e6,

  /**
  The initial temperature and the interface
  temperature are set to the same value. */

  Tsat = 373.15, Twall = 383.15;

  /**
  We change the dimension of the domain
  and the surface tension coefficient. */

  L0 = 10e-3;
  
  /**
  We define a list with the maximum time
  steps and the maximum levels of refinement. */

  double dtlist[] = {0.01, 0.005, 0.001};
  DT = dtlist[maxlevel - 4];
  TOLERANCE = 1.0e-6;
  NITERMAX = 300;
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

  vertex scalar phi[];
  foreach_vertex(){
    phi[] = x - delta0;
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

  tshift = sq(effective_height/2./lambdaval)*rho2*cp2/lambda2;

  foreach() {
    TL[] = Tsat;
    TG[] = x < delta0 ? tempsol(t + tshift, x) : Tsat;
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
    if(vofs == 0.0)
      color_pha[] = 0.0;
    else if(vofs == 4.0)
      color_pha[] = 1.0;
#if USE_MY_SOLID
    color_pha[] *= (1.0 - is_solid_vertex[]);
#endif
  }
}
#endif


double exact (double time) {
  return 2.*lambdaval*sqrt(lambda2/rho2/cp2*time);
}

/**
### Output Files

We write the thickness of the vapor layer and the analytic
solution on a file.
*/

event movie (t += 0.1) {
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
      double tempexact = x < R ? tempsol(t + tshift, x) : Tsat;
      fprintf(fpp, "%g %g %g\n", x, temps[count], tempexact);
      count++;
    }
    fflush(fpp);
    fclose(fpp);
  }
  array_free(arrtemp);

  int Ncell = 1 << grid->maxdepth;

  //the last facial velocity will be omitted
  double *ucx = (double *)calloc(Ncell, sizeof(double));
  double *ufx = (double *)calloc(Ncell, sizeof(double));

  for(int ii = 0; ii < Ncell; ++ii)
  {
    ucx[ii] = 0.0;
    ufx[ii] = 0.0;
  }

  foreach_boundary(top, reduction(+ : ucx[:Ncell]) reduction(+ : ufx[:Ncell]))
  {
#if !USE_DOUBLE_VEL
    ucx[point.i - 2] = u.x[];
    ufx[point.i - 2] = uf.x[];
#else
    ucx[point.i - 2] = x > effective_height ? u.x[] : u2.x[];
    ufx[point.i - 2] = (x - 0.5 * Delta) > effective_height ? uf.x[] : uf2.x[];
#endif
  }
  char name_uc[60];
  char name_uf[60];
  sprintf(name_uc, "data/uc-%d", maxlevel);
  sprintf(name_uf, "data/uf-%d", maxlevel);
  if (pid() == 0)
  {
    FILE *fpc = fopen(name_uc, "w");
    FILE *fpf = fopen(name_uf, "w");
    double delta = L0 / Ncell;
    int count = 0;
    for (double x = 0.5 * delta; x < L0; x += delta)
    {
      fprintf(fpc, "%g %g\n", x, ucx[count]);
      fprintf(fpf, "%g %g\n", x - 0.5 * delta, ufx[count]);
      ++count;
    }
    fflush(fpc);
    fclose(fpc);
    fflush(fpf);
    fclose(fpf);

    // exact velocity
    double alpha_g = lambda2 / rho2 / cp2;
    double ugamma = lambdaval * sqrt(alpha_g / (t + tshift));
    double uliq = ugamma - rho2 * ugamma / rho1;

    sprintf(name, "data/uexact");
    FILE *fpp = fopen(name, "w");
    fprintf(fpp, "%g %g\n", exact(t + tshift), uliq);
    fclose(fpp);
  }
  free(ucx);
  free(ufx);

  sprintf(name, "data/dump-%d", maxlevel);
  dump(name);
  // do nothing

}
