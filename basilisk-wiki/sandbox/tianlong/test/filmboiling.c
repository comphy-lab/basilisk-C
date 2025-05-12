/**
The film boiling setup. This file is copied from [Edorado's sandbox](http://basilisk.fr/sandbox/ecipriano/README) and modified.
*/
#define SEMUPC 1
#define ADV_SCHEME 2
#define INTGRAD_3rd 1
#define USE_DOUBLE_VEL 1
#define TEST "data/"

#include "axi.h"
#if USE_DOUBLE_VEL
#include "double-evaporation.h"
#else
#include "centered-evaporation.h"
#endif //USE_DOUBLE_VEL
#include "tag.h"
#include "semushin_two-phase.h"
#include "semushin-phase-change2.h" //with pinch off and outflow
#include "mytension.h"
#include "myreduced.h"
//#include "view.h"

#if USE_MY_SOLID
#include "mysolid.h"
#endif

int maxlevel = 8, minlevel = 6;
double ldnum = 0.0787;
double Twall;
const double t_stop = 1.6;
double gas_vol0;
double gas_vol_out = 0.0;

const double femax = 1e-6;
const int num_refine = 2;

//out flow for the right
u.n[right] = neumann (0);
p[right] = dirichlet (0);
pf[right] = dirichlet (0);

//top and bottom symmetry
uf.n[top] = 0.0;
uf.n[bottom] = 0.0;

//No-slip wall for the left
uf.n[left] = 0.0;
u.t[left] = dirichlet (0);
TL[left] = dirichlet (Twall);
TG[left] = dirichlet (Twall);
T[left] = dirichlet (Twall);

#if USE_MY_SOLID
void solidSetup()
{
  double delta = L0 / (1 << 8);
  double len = ceil((L0 - 0.5 * 0.078644) / delta) * delta;
  SOLID_LEN_y = len;
  ldnum = (L0 - len) * 2.0; // as the solid boundary will be shifted a little
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
  rho1 = 200.0; rho2 = 5.0;
  mu1 = 0.1; mu2 = 0.005;

  //we still set here, though they are not needed in this case
  lambda1 = 40., lambda2 = 1.;
  cp1 = 400., cp2 = 200.;
  dhev = 1e4;

  /**
  The initial temperature and the interface
  temperature are set to the same value. */

  Tsat = 1.0;
  Twall = Tsat + 5.0;
  /**
  We change the dimension of the domain
  and the surface tension coefficient. */
  f.sigma = 0.1;
  L0 = 0.23967695;
  G.x = -9.81;
  TOLERANCE = 1.0E-6;
  /**
  We define a list with the maximum time
  steps and the maximum levels of refinement. */
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
    if (IS_SOLID_x)
    {
      f.boundarySolid_x = boundarySolidNeumman_x;
      uf.x.boundarySolid_x = boundarySolidVelF_x;
#if USE_DOUBLE_VEL
      uf2.x.boundarySolid_x = boundarySolidVelF_x;
#endif
      color_pha_cen.boundarySolid_x = boundarySolidNeumman_x;
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
  foreach_vertex()
  {
    phi[] = x - ldnum / 128.0 * (4 + cos(2.0 * pi * y / ldnum));
  }

  init_markers(phi);

  foreach (reduction(+:gas_vol0))
  {
#if USE_MY_SOLID
    if((int)is_solid[] == 1)
      continue;
#endif
    gas_vol0 += (1. - f[]) * dv();
  }

  foreach ()
  {
    double pos_interf = ldnum / 128.0 * (4 + cos(2.0 * pi * y / ldnum));
    TL[] = Tsat;
    TG[] = x > pos_interf ? Tsat : Twall - 5.0 * x / pos_interf;
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

#if TREE
event adapt (i++) {
  scalar solid_refine[];
  scalar solid_refiney[];
  scalar f_refine[];

  fillRefineVOFs(f_refine, solid_refine, solid_refiney);

  scalar T_refine[];
  T_refine[left] = dirichlet(Twall);
  foreach ()
  {
    T_refine[] = T[];
    if ((int)is_solid[] == 1)
    {
      T_refine[] = Twall;
    }
  }

  adapt_wavelet({T_refine, f_refine, solid_refine, solid_refiney}, (double[]){1.e-3, femax, femax, femax},
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

/**
### Output Files

We write the thickness of the vapor layer and the analytic
solution on a file.
*/

int dump_num = 0;
event movie(t += 0.01)
{
  double gas_vol = 0.;
  foreach (reduction(+:gas_vol))
  {
#if USE_MY_SOLID
    if ((int)is_solid[] == 1)
      continue;
#endif
    gas_vol += (1. - f[]) * dv();
  }
  gas_vol += gas_vol_out;
  double ld2 = sqrt(f.sigma/(rho1 - rho2)/9.81);
  double Nu = 0.;
  foreach_boundary (left, reduction(+:Nu)) {
#if USE_MY_SOLID
    if ((int)is_solid[] == 1)
      continue;
#endif
    double Tghost = 2.*Twall - T[];
    Nu += (T[] - Tghost) * y;
  }
  Nu *= -ld2 * 2.0/(Twall - Tsat)/(0.5 * ldnum) / (0.5 * ldnum);

  char name[80];
  sprintf(name, "data/OutputData-%d", maxlevel);

  static FILE *fp = fopen(name, "w");
  fprintf(fp, "%g %g %g\n", t, gas_vol / gas_vol0, Nu);
  fflush(fp);

  sprintf (name, "data/infsemu-%d", dump_num);
  output_facets_semushin(name);

  sprintf(name, "data/data-%03d.gfs", dump_num);

  FILE *file = fopen(name, "w");
  output_gfs(file);

  sprintf(name, "data/dump-%03d", dump_num++);
  dump(name);

}


event finalEvent(t = t_stop)
{
  char name[60];
  sprintf(name, "data/infsemu-%g", t);
  output_facets_semushin(name);
}
