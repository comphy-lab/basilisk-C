#define USE_CONJUGATE_HEAT 2   //0:no conjugate heat-transfer 1:explicit coupling 2:implicit coupling

#define INT_USE_UF
#define CONSISTENTPHASE2
#define SHIFT_TO_GAS
#define INIT_TEMP


#include "axi.h"
#include "mypoisson.h"   //I don't understant why, but I has to add this line here to compile the code
#include "centered-evaporation.h"
#include "centered-doubled.h"
#include "two-phase.h"
#include "mytension.h"
#include "evaporation.h"
#include "temperature-gradient.h"
#if USE_MY_SOLID
#include "mysolid.h"
#endif

double lambda1, lambda2, cp1, cp2, dhev;
double TL0, TG0, TIntVal, Tsat, Tbulk;

#if USE_CONJUGATE_HEAT
double rhos, lambdas, cps;
#endif //USE_CONJUGATE_HEAT

double rhos2, lambdas2, cps2;
double thickness_heater;
scalar vof_heater[]; //volume fraction of heater


int maxlevel = 10, minlevel = 6;
int dump_num = 0;

const double femax = 1.0e-6;
const int num_refine = 5;

double Lsize = 7.0e-3;
const double final_time = 88.1e-3;
int out_num = 10;
double heatflux = 481;
const double tshift = 0.0;

//outlet boundary for top and right
u.n[top] = neumann(0);
p[top] = dirichlet(0);
u.n[right] = neumann(0);
p[right] = dirichlet(0);

void solidSetup()
{
  double delta = Lsize / (1024.0);
  int ncs = ceil(1.0e-3 / delta);
  double lenx = ncs * delta;
  if(pid() == 0.0)
  {
    printf("solid_len_x: %g\n",lenx);
  }
  SOLID_LEN_x = lenx;
  SOLID_LEN_y = 0.0;//0.00099951171875;
}

#include "myfunctions.h"

int Ncell = 64;

double delta_q_case(double x, double t)
{
  double grad = (0.004 - x) / 0.004 * heatflux;
  grad = max(grad * 1000.0, 0.0);

  // double pos = (1.0 + sqrt(2.0)) / 2.0 * 1.5 * 1.0e-3;
  // grad = x <= pos ? 425e3 : 0.0;
  return grad;
}

void initProperties();

int main (int argc, char * argv[]) {

  if (pid() == 0)
  {
    if (argc == 2)
    {
      maxlevel = atoi(argv[1]);
    }

    if (argc == 3)
    {
      maxlevel = atoi(argv[1]);
      heatflux = atof(argv[2]);
    }
  }

#if _MPI
  MPI_Bcast(&maxlevel, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&heatflux, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  Ncell = 1 << 8;
  init_grid(Ncell);
  solidSetup();
  size(Lsize);
  origin(-SOLID_LEN_x, -SOLID_LEN_y);
  DT = 50.0E-5;
  //  NITERMAX = 300;
  TOLERANCE = 1.0E-6;
  TOL_DIFFUSION = 1.0E-5;

  if(pid() == 0.0)
  {
    printf("heatflux = %g\n", heatflux);
  }

  initProperties();
  run();
}

void initProperties()
{
    // fluid 2
    rho2 = 0.598;
    mu2 = 1.22e-5;
    lambda2 = 0.0246;
    cp2 = 2080;

    // fluid 1
    rho1 = 958;
    mu1 = 2.82e-4;
    lambda1 = 0.677;
    cp1 = 4220;

    dhev = 2.26e6;
    Tsat = 373.15;
    TIntVal = 373.15;

    rhos = 4510.0;
    lambdas = 17.0;
    cps = 544.0;

    rhos2 = 3980.0;
    lambdas2 = 25.1;
    cps2 = 929.0;

    thickness_heater = 500e-9;
    double Rwater = 461.52277955158064;
    double omega = 0.0345;
    double coef = sqrt(Tsat * Tsat * Tsat * 2.0 * pi * Rwater) / omega / dhev / dhev / rho2;
    Rcc = 0.0;//coef;

    delta_q_volume = delta_q_case;
    //delta_q = delta_q_case;

    f.sigma = 0.0;//0.0589;
}

event defaults(i = 0)
{
    // set the related boundary condition
    for (int ib = 0; ib < nboundary; ib++)
    {
      pf.boundary[ib] = p.boundary[ib];
      pext.boundary[ib] = p.boundary[ib];
      pfext.boundary[ib] = pf.boundary[ib];

      pf.boundary_homogeneous[ib] = p.boundary_homogeneous[ib];
      pext.boundary_homogeneous[ib] = p.boundary_homogeneous[ib];
      pfext.boundary_homogeneous[ib] = pf.boundary_homogeneous[ib];

      foreach_dimension()
      {
        uext.x.boundary[ib] = u.x.boundary[ib];
        ufext.x.boundary[ib] = uf.x.boundary[ib];

        uext.x.boundary_homogeneous[ib] = u.x.boundary_homogeneous[ib];
        ufext.x.boundary_homogeneous[ib] = uf.x.boundary_homogeneous[ib];
      }
    }

#if USE_MY_SOLID
  setSolidFlag();
  foreach_dimension()
  {
    if(IS_SOLID_x)
    {
      f.boundarySolid_x = boundarySolidNeumman_x;
      fu.boundarySolid_x = boundarySolidNeumman_x;
      fuext.boundarySolid_x = boundarySolidNeumman_x;
      uf.x.boundarySolid_x = boundarySolidVelF_x;
      ufext.x.boundarySolid_x = boundarySolidVelF_x;
      uf_save.x.boundarySolid_x = boundarySolidVelF_x;

      mEvapTot.boundarySolid_x = boundarySolidNeumman_x;
      fL.boundarySolid_x = boundarySolidNeumman_x;
      fG.boundarySolid_x = boundarySolidNeumman_x;
      T.boundarySolid_x = boundarySolidNeumman_x;
      TL.boundarySolid_x = boundarySolidNeumman_x;
      TG.boundarySolid_x = boundarySolidNeumman_x;
    }
  }

#if USE_CONJUGATE_HEAT
  TL.boundarySolid_x = NULL;
  TG.boundarySolid_x = NULL;
  TS.boundarySolid_y = boundarySolidNeumman_y;
#endif //USE_CONJUGATE_HEAT

  foreach()
  {
    if(is_solid[])
    {
      f[] = 0.0;
    }
  }
  
  boundary({f});
#endif

  vof_heater.refine = vof_heater.prolongation = fraction_refine;
  foreach()
  {
    vof_heater[] = is_solid[] * (1.0 - is_solid_y[]);
    if(vof_heater[] > 0.0)
    {
      double xm = x - Delta / 2.0;
      double xp = x + Delta / 2.0;

      xm = xm + thickness_heater;
      xp = xp + thickness_heater;

      xm = max(xm, 0.0);
      xp = max(xp, 0.0);

      double frac = (xp - xm) / Delta;
      vof_heater[] = clamp(frac, 0.0, 1.0);
    }
  }
}

void fillRefineVOFs(scalar solid_refine, scalar solid_refiney)
{
  foreach ()
  {
    solid_refine[] = is_solid[];
    solid_refiney[] = is_solid_y[];
  }

  foreach ()
  {
    if (is_solid[] == 1.0 && is_solid[1] == 0.0)
    {
      solid_refine[] = 0.5;
    }
    if (is_solid[] == 0.0 && is_solid[-1] == 1.0)
    {
      solid_refine[] = 0.5;
    }
    if (is_solid_y[] == 0.0 && is_solid_y[0, -1] == 1.0)
    {
      solid_refiney[] = 0.5;
    }
    if (is_solid_y[] == 1.0 && is_solid_y[0, 1] == 0.0)
    {
      solid_refiney[] = 0.5;
    }
  }

  setJumpVarAdaptation(solid_refine, num_refine);
  setJumpVarAdaptation(solid_refiney, num_refine);
}

void refineGrid()
{
  scalar solid_refine[];
  scalar solid_refiney[];

  for(int il = 0; il < maxlevel - 8; ++il)
  {
    fillRefineVOFs(solid_refine, solid_refiney);
    adapt_wavelet({solid_refine, solid_refiney},(double[]){femax, femax}, maxlevel, minlevel);
  }

  Ncell = 1 << maxlevel;

  return;
}

event init (i = 0) {
  foreach()
  {
    f[] = 0.0;
  }

  foreach()
  {
    f[] = 1. - f[];
    fL[] = f[];
    fG[] = 1.0 - f[];
    fu[] = f[];
    fuext[] = f[];
  }

  IS_SOLID_x = false;
  IS_SOLID_y = false;

  refineGrid();
  setSolidFlag();

  foreach()
  {
    double Tlinear = Tsat;
    TL[] = Tsat;
    TG[] = Tsat;
    TS[] = Tsat;
    TL[] *= f[];
    TG[] *= (1.0 - f[]);
    T[] = TL[] + TG[];
  }

  foreach ()
  {
    vof_heater[] = is_solid[] * (1.0 - is_solid_y[]);
    if (vof_heater[] > 0.0)
    {
      double xm = x - Delta / 2.0;
      double xp = x + Delta / 2.0;

      xm = xm + thickness_heater;
      xp = xp + thickness_heater;

      xm = max(xm, 0.0);
      xp = max(xp, 0.0);

      double frac = (xp - xm) / Delta;
      vof_heater[] = clamp(frac, 0.0, 1.0);
    }
  }
}

#if TREE
event adapt (i++) {
  scalar solid_refine[];
  scalar solid_refiney[];
  foreach()
  {
    solid_refine[] = 0.0;
    solid_refiney[] = 0.0;
#if USE_MY_SOLID
    solid_refine[] = is_solid[];
    solid_refiney[] = is_solid_y[];
#endif
  }

#if USE_MY_SOLID
  foreach()
  {
    if(is_solid[] == 1.0 && is_solid[1] == 0.0)
    {
      solid_refine[] = 0.5;
    }
    if(is_solid[] == 0.0 && is_solid[-1] == 1.0)
    {
      solid_refine[] = 0.5;
    }
    if (is_solid_y[] == 0.0 && is_solid_y[0, -1] == 1.0)
    {
      solid_refiney[] = 0.5;
    }
    if (is_solid_y[] == 1.0 && is_solid_y[0, 1] == 0.0)
    {
      solid_refiney[] = 0.5;
    }
  }
#endif
  setJumpVarAdaptation(solid_refine, num_refine);
  setJumpVarAdaptation(solid_refiney, num_refine);
  adapt_wavelet({solid_refine, solid_refiney},
                (double[]){femax, femax}, maxlevel, minlevel);
#if USE_CONJUGATE_HEAT
  vector solidfacex[],solidfacey[];
  vector coefs[];

  foreach()
  {
    foreach_dimension()
    {
      solidfacex.x[] = is_solid_x_face.x[];
      solidfacey.x[] = is_solid_y_face.x[];
      coefs.x[] = lambdasf.x[];
    }
  }
#endif //USE_CONJUGATE_HEAT
}
#endif


void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
	    mg.nrelax);
}

event outLog(i++)
{
  if(pid() == 0 && i % 1 == 0)
  {
    printf("i = %d t = %g dt = %g\n", i, t, dt);
    mg_print (mgT);
    putchar ('\n');
  }
}

/**
### Movie

We write the animation with the evolution of the
temperature field and the gas-liquid interface. */
event output(t+=final_time/(double)out_num)
{
  char filename[20];
  sprintf(filename, "data/dump-%03d", dump_num++);
  dump(filename, list = {TS, is_solid, vof_heater});

  scalar Tout[];
  foreach()
  {
    Tout[] = TS[];
  }
  boundarySolidNeummanNoauto(Tout);

  scalar Tout2[];
  foreach()
  {
    Tout2[] = TS[];
  }

  foreach ()
  {
    if (is_solid[] == 0.0)
    {
      if (is_solid[-1] == 1.)
      {
        double T_sol = Tout2[-1];
        double T_liq = Tout2[];
        double tmp = getTbcS(Delta, T_sol, T_liq, 1.0, 0.0);
        double tmp2 = getTbcF(Delta, T_sol, T_liq, 1.0, 0.0);
        Tout2[] = 2.0 * tmp2 - T_sol;
      }
    }
  }

  char Tname[80];
  sprintf(Tname, "data/T_wall-%03d.dat", dump_num - 1);
  outputWallTemperature(Tname, Tout2);

}

int num_compare = 0;
event outTempCompare(t = {25.0e-3,50.0e-3,75.0e-3, 88.1e-3})
{
  char Tname[80];
  scalar Tout[];
  foreach()
  {
    Tout[] = TS[];
  }

  boundarySolidNeummanNoauto(Tout);

  scalar Tout2[];
  foreach()
  {
    Tout2[] = TS[];
  }

  foreach ()
  {
    if (is_solid[] == 0.0)
    {
      if (is_solid[-1] == 1.)
      {
        double T_sol = Tout2[-1];
        double T_liq = Tout2[];
        double tmp = getTbcS(Delta, T_sol, T_liq, 1.0, 0.0);
        double tmp2 = getTbcF(Delta, T_sol, T_liq, 1.0, 0.0);
        Tout2[] = 2.0 * tmp2 - T_sol;
      }
    }
  }

  sprintf(Tname, "data/T_wall-comp-%03d.dat", ++num_compare);
  outputWallTemperature(Tname, Tout2);
}