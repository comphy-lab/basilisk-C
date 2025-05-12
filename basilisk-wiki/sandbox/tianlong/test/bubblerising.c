/**
The setup for bubble rising in superheated liquid
*/
#define SEMUPC 1
#define INTGRAD_3rd 1
#define USE_DOUBLE_VEL 1
#define ADV_SCHEME 1
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

#if USE_MY_SOLID
#include "mysolid.h"
#endif

#include <gsl/gsl_integration.h>

//out flow for the right
u.n[right] = neumann (0);
p[right] = dirichlet (0);
pf[right] = dirichlet (0);

//wall for the left
u.t[left] = dirichlet (0);
uf.n[left] = 0.0;

//wall for the top
u.t[top] = dirichlet (0);
uf.n[top] = 0.0;

//symmetry imposed in the solid


int maxlevel = 10, minlevel = 6;
double Tbulk;
double betaGrowth;
const double tshift = 0.0056;
const double R0 = 210.0e-6;
const double t_stop = 0.08;

const double femax = 1e-6;
const int num_refine = 2;


#if USE_MY_SOLID
void solidSetup()
{
  double delta = L0 / (double)(1 << 10);
  int ncs = ceil(16.0e-3 / delta);
  double len = ncs * delta;
  SOLID_LEN_x = 0.0;
  SOLID_LEN_y = len;
}
#endif

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

  rho1 = 757.0, rho2 = 1.435;
  mu1 = 4.29e-4, mu2 = 1.04e-5;
  lambda1 = 0.154, lambda2 = 0.02;
  cp1 = 3000., cp2 = 1830;
  dhev = 9.63e5;

  /**
  The initial temperature and the interface
  temperature are set to the same value. */

  Tbulk = 354.55, Tsat = 351.45;
  f.sigma = 0.018;

  /**
  We change the dimension of the domain
  and the surface tension coefficient. */

  L0 = 20e-3;
  
  /**
  We define a list with the maximum time
  steps and the maximum levels of refinement. */
  TOLERANCE = 1E-6;
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

  double alpha = lambda1/rho1/cp1;
  betaGrowth = R0 / 2.0 / sqrt(alpha * tshift);

  vertex scalar phi[];
  foreach_vertex()
  {
    phi[] = sqrt((x - 1.0e-3) * (x - 1.0e-3) + y * y) - R0;
  }

  init_markers(phi);

  foreach ()
  {
    double r = sqrt((x - 1.0e-3) * (x - 1.0e-3) + y * y);
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


event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 9.81;
  boundary ((scalar *){av});
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
int dump_num = 0;
event movie (t += t_stop / 100.0) 
{
  char name[60];
  sprintf(name,"data/dump-%03d", dump_num);
  dump(name);

  if (dump_num != 0 && dump_num % 10 == 0)
  {
    sprintf(name, "data/data-%03d.gfs", dump_num);

    FILE *file = fopen(name, "w");
    output_gfs(file);
  }

  sprintf (name, "data/infsemu-%03d", dump_num);
  output_facets_semushin(name);

  dump_num++;
}

event finalEvent(t = t_stop)
{
  
}
