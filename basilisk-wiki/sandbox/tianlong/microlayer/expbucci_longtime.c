
#define USE_SMALL_ANGLE 1
#define USE_CONTACT_ANGLE 1
#define USE_CONJUGATE_HEAT 2   //0:no conjugate heat-transfer 1:explicit coupling 2:implicit coupling
#define USE_GFM 1

#if USE_GFM
//#define INT_USE_UF
#define NOSHIFTING
#else
#define INT_USE_UF
#define SHIFT_TO_GAS
#endif

#define CONSISTENTPHASE2
#define INIT_TEMP
#define FILTERED
#define INTGRAD_3rd


#include "axi.h"
#if USE_GFM
#include "velocity-jump.h"
#else
#include "centered-evaporation.h"
#include "centered-doubled.h"
//#include "velocity-potential.h"
#endif
#include "contact-evaporation.h"
#include "mytwo-phase.h"
#include "mytension.h"
#include "evaporation.h"
#include "temperature-gradient.h"
#include "mysolid.h"
#include "mytag.h"
#include "mybc.h"
#include "myfunctions.h"

vector h[];
double theta0 = 5.0;
h.t[left] = contact_angle (theta0*pi/180.);

double lambda1, lambda2, cp1, cp2, dhev;
double TL0, TG0, TIntVal, Tsat, Tbulk;

double rhos, lambdas, cps;
double rhos2, lambdas2, cps2;
double thickness_heater;
scalar vof_heater[]; //volume fraction of heater
double pos_old = 0.0;
double thick_old = 0.0;

int is_restart = 0;
int maxlevel = 11, minlevel = 6;
int dump_num = 0;
int output_period = 180;

const double uemax = 1e-2;
const double Temax = 1e-2;
const double femax = 1e-6;
const double film_height = 15e-3;   //this will not be used in this case
const double remove_height = 8.0e-3; //this will not be used in this case

const int num_refine = 1;

double Lsize = 5.5e-3;

const double R0 = 20e-6;
const double final_time = 18.0e-3;
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
  ncs = ceil(3.0e-3 / delta);
  double leny = ncs * delta;
  if(pid() == 0.0)
  {
    printf("solid_len_y: %g\n",0.0);
  }
  SOLID_LEN_x = lenx;
  SOLID_LEN_y = 0.0;//leny;
}

int Ncell = 64;

double delta_q_case(double x, double t)
{
  
#if 1
  double grad = (0.004 - x) / 0.004 * 481.0;
  grad = max(grad * 1000.0, 0.0);
#else
  double grad = 425e3;
  double pos = (1.0 + sqrt(2.0)) / 2.0 * 1.5 * 1.0e-3;
  grad = x <= pos ? 425e3 : 0.0;
#endif
  return grad;
}

double Hold = 0.0;
double theta_old = 5.0;
double omega = 0.0345;
bool is_microlayer_vanish = false;

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
      omega = atof(argv[2]);
    }

    if (argc == 4)
    {
      maxlevel = atoi(argv[1]);
      omega = atof(argv[2]);
      output_period = atoi(argv[3]);
    }

    if (argc == 5)
    {
      maxlevel = atoi(argv[1]);
      omega = atof(argv[2]);
      output_period = atoi(argv[3]);
      is_restart = 1;
      dump_num = atoi(argv[4]);
    }

    if (argc == 6)
    {
      maxlevel = atoi(argv[1]);
      omega = atof(argv[2]);
      output_period = atoi(argv[3]);
      is_restart = 1;
      dump_num = atoi(argv[4]);
      theta_old = atof(argv[5]);
    }
  }

#if _MPI
  MPI_Bcast(&maxlevel, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&omega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&is_restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&output_period, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dump_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&theta_old, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  Ncell = 1 << 8;
  init_grid(Ncell);
  solidSetup();
  size(Lsize);
  origin(-SOLID_LEN_x, -SOLID_LEN_y);
  CFL = 0.1;
  CFL_SF = 0.5;
  //DT = 2.0E-8;
  //  NITERMAX = 300;
  TOLERANCE = 1.0E-6;
  TOL_DIFFUSION = 1.0E-5;
  theta_old = theta_old * pi / 180.0;
  f.height = h;
  // foreach_dimension()
  //   h.x.nodump = true;

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
    //double omega = 0.0460;//0.0345;//0.0460;
    double coef = sqrt(Tsat * Tsat * Tsat * 2.0 * pi * Rwater) / omega / dhev / dhev / rho2;
    Rcc = coef;

    delta_q_volume = delta_q_case;
    //delta_q = delta_q_case;

    f.sigma = 0.0589;
}

event defaults(i = 0)
{
  // set the related boundary condition
  for (int ib = 0; ib < nboundary; ib++)
  {
    pf.boundary[ib] = p.boundary[ib];
    pf.boundary_homogeneous[ib] = p.boundary_homogeneous[ib];
#if VEL_POTENTIAL
    ps.boundary[ib] = p.boundary[ib];
    ps.boundary_homogeneous[ib] = p.boundary_homogeneous[ib];
#elif VELOCITY_JUMP
    //no additonal pressure boundary condition will be needed
#else
    pext.boundary[ib] = p.boundary[ib];
    pfext.boundary[ib] = pf.boundary[ib];

    pext.boundary_homogeneous[ib] = p.boundary_homogeneous[ib];
    pfext.boundary_homogeneous[ib] = pf.boundary_homogeneous[ib];
#endif

    foreach_dimension()
    {
#if VEL_POTENTIAL
      ufs.x.boundary[ib] = uf.x.boundary[ib];
      ufs.x.boundary_homogeneous[ib] = uf.x.boundary_homogeneous[ib];
#elif VELOCITY_JUMP
      u1.x.boundary[ib] = u.x.boundary[ib];
      u1.x.boundary_homogeneous[ib] = u.x.boundary_homogeneous[ib];
      u2.x.boundary[ib] = u.x.boundary[ib];
      u2.x.boundary_homogeneous[ib] = u.x.boundary_homogeneous[ib];
      uf1.x.boundary[ib] = uf.x.boundary[ib];
      uf1.x.boundary_homogeneous[ib] = uf.x.boundary_homogeneous[ib];
      uf2.x.boundary[ib] = uf.x.boundary[ib];
      uf2.x.boundary_homogeneous[ib] = uf.x.boundary_homogeneous[ib];
#else
      uext.x.boundary[ib] = u.x.boundary[ib];
      uext.x.boundary_homogeneous[ib] = u.x.boundary_homogeneous[ib];
#endif
      ufext.x.boundary[ib] = uf.x.boundary[ib];
      ufext.x.boundary_homogeneous[ib] = uf.x.boundary_homogeneous[ib];
    }
  }

  setSolidFlag();
  foreach_dimension()
  {
    if (IS_SOLID_x)
    {
      f.boundarySolid_x = boundarySolidNeumman_x;
      fu.boundarySolid_x = boundarySolidNeumman_x;
      fuext.boundarySolid_x = boundarySolidNeumman_x;
      uf.x.boundarySolid_x = boundarySolidVelF_x;
      ufext.x.boundarySolid_x = boundarySolidVelF_x;
      uf_save.x.boundarySolid_x = boundarySolidVelF_x;
#if VEL_POTENTIAL
      ufs.x.boundarySolid_x = boundarySolidVelF_x;
#endif
#if VELOCITY_JUMP
      mEvapTotE.boundarySolid_x = boundarySolidNeumman_x;
      n.x.boundarySolid_x = boundarySolidVectorNeumann_x;
      uf1.x.boundarySolid_x = boundarySolidVelF_x;
      uf2.x.boundarySolid_x = boundarySolidVelF_x;
#endif
      mEvapTot.boundarySolid_x = boundarySolidNeumman_x;
      fL.boundarySolid_x = boundarySolidNeumman_x;
      fG.boundarySolid_x = boundarySolidNeumman_x;
      T.boundarySolid_x = boundarySolidNeumman_x;
      TL.boundarySolid_x = boundarySolidNeumman_x;
      TG.boundarySolid_x = boundarySolidNeumman_x;
    }
  }

#if USE_SMALL_ANGLE
  f.boundarySolid_x = boundarySolidVOFSmallAngles_x;
#endif

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
}

#define circle(x, xc, y, yc, R) (sq(R) - sq(x - xc) - sq(y - yc))
#define lineVOF(x, film_height) (film_height - x)

void initVOF(scalar f, double xc, double yc, double R, double fheight)
{
  vertex scalar phi[];
  foreach_vertex()
  {
    double d1 = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc)) - R;
    double d2 = fheight - x;
    phi[] = min(d1, d2);
  }

  fractions(phi, f);
}

void fillRefineVOFs(scalar f_refine, scalar solid_refine, scalar solid_refiney)
{
  foreach ()
  {
    f_refine[] = f[];
    solid_refine[] = is_solid[];
    solid_refiney[] = is_solid_y[];

    f_refine[] *= (1.0 - is_solid[]);
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
}

event init (i = 0) {
  if (is_restart == 1)
  {
    char name[30];
    sprintf(name, "./data/dump-%05d", dump_num - 1);

    // mpi_boundary will be called in restore
    //  for now the solid boundary coniditons shouldn't be considered
    IS_SOLID_x = false;
    IS_SOLID_y = false;

    restore(file = name, list = all);

    // recover the solid flags
    setSolidFlag();
  }
  else
  {
    double xc = 0.0;
    double yc = 0.0;
    char dump_name[30];
    sprintf(dump_name, "dumpTini-large");
#if _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    IS_SOLID_x = false;
    IS_SOLID_y = false;

    restoreModified(dump_name);

    xc = 0.0;//R0 * cos(theta0 * pi / 180.0);
    yc = 0.0;
    initVOF(f, xc, yc, R0, film_height);

    scalar f_refine[];
    scalar solid_refine[];
    scalar solid_refiney[];

    for (int il = 0; il < (maxlevel - minlevel); ++il)
    {
      fillRefineVOFs(f_refine, solid_refine, solid_refiney);
      setJumpVarAdaptation(solid_refine, num_refine);
      setJumpVarAdaptation(solid_refiney, num_refine);
      setJumpVarAdaptation(f_refine, num_refine);
      adapt_wavelet({f_refine, solid_refine, solid_refiney}, (double[]){femax, femax, femax}, maxlevel, minlevel);
      initVOF(f, xc, yc, R0, film_height);
    }

    setSolidFlag();
    foreach ()
    {
      fL[] = f[];
      fG[] = 1.0 - f[];
      fu[] = f[];
      fuext[] = f[];
    }

    boundary({f, fL, fG, fu, fuext});

    foreach ()
    {
      TL[] = TS[];
      TG[] = Tsat;
      TL[] *= f[];
      TG[] *= (1.0 - f[]);
      T[] = TL[] + TG[];
    }

    vof_heater.refine = vof_heater.prolongation = fraction_refine;
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

  double Hnew = 0.0;
  foreach(reduction(+:Hnew))
  {
    if((int)is_solid[] == 0 && (int)is_solid[-1] == 1)
    {
      Hnew += f[] * Delta;
    }
  }

  Hold = Hnew;
}

event beforePhaseChange(i++)
{

  int num = max(0, (maxlevel - 11));
  int size = 10 * (1 << num);
  if (i > 100)
  {
    remove_droplets_phasechange(f, threshold = 1.0e-20, minsize = size);
    remove_droplets_phasechange(f, threshold = 1.0e-20, minsize = size, bubbles = true);
  }

  boundary({f});

}

#if TREE
event adapt (i++) {
  scalar solid_refine[];
  scalar solid_refiney[];
  scalar f_refine[];

  vector u_refine[];

  for (int ib = 0; ib < nboundary; ib++)
  {
    foreach_dimension()
      u_refine.x.boundary[ib] = u.x.boundary[ib];
  }

  foreach()
  {
    f_refine[] = f[] * (1.0 - is_solid[]);
  }
  
  foreach()
  {
    f[] = f_refine[];
  }

  foreach()
  {
    f_refine[] = f[];
    foreach_dimension()
    {
      u_refine.x[] = f[] > 0.5 ? u.x[] : (ufext.x[] + ufext.x[1])/(fm.x[] + fm.x[1]);
    }
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
  setJumpVarAdaptation(f_refine, num_refine);
  boundarySolidNeummanNoauto(f_refine);
  adapt_wavelet({u_refine, solid_refine, solid_refiney, f_refine},
                (double[]){uemax, uemax, femax, femax, femax}, maxlevel, minlevel);

}
#endif

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
  {
#if USE_MY_SOLID
    if((int)is_solid_face.x[] == 1)
      continue;
#endif
    double coef = 1.0;
    av.x[] -= coef * 9.81;
  }
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

  double Hnew = 0.0;
  foreach(reduction(+:Hnew))
  {
    if((int)is_solid[] == 0 && (int)is_solid[-1] == 1)
    {
      Hnew += f[] * Delta;
    }
  }

  double delta = Lsize / (1 << grid->maxdepth);

  //here we modify the contact angles
  double theta_adv = 55.0 * pi / 180.0;
  double theta_rec = 5.0 * pi / 180.0;
  double theta_change = theta_old + 2.0 * sin(theta_old) * sin(theta_old) * (Hnew - Hold) / delta;


  if(t > 1e-3 && !is_microlayer_vanish) 
  {
    //if there is no flat interface within two layers
    //we think microlayer already vanish
    int num_flat = 0;
    foreach(reduction(+:num_flat))
    {
      if((int)is_solid[] == 0 && ((int)is_solid[-1] == 1 || (int)is_solid[-2] == 1))
      {
        if(f[] > 0.0 && f[] < 1.0)
        {
          coord nc = normal(point, f);
          double tang = fabs(nc.y / nc.x);
          if(tang < 0.015)
          {
            num_flat++;
          }
        }
      }
    }

    is_microlayer_vanish = num_flat == 0;
  }

  theta_old = max(theta_rec, min(theta_adv, theta_change));
  theta0 = theta_old * 180.0 / pi;

  if(is_microlayer_vanish)
  {
    Rcc = 0.0;
  }

  if(pid() == 0 && i % 1 == 0)
  {
    printf("i = %d t = %g dt = %g theta = %g\n", i, t, dt, theta0);
  }

  Hold = Hnew;
}

/**
### Output Files

We write the bubble radius and the analytic solution on 
a file. */

event output(t += final_time / (double)output_period)
{
  scalar fg[];
  foreach()
    fg[] = 1. - f[];
  foreach()
  {
    fg[] *= (1.0 - is_solid[]);
  }

  double coef = 2.0;
  double effective_radius = pow(3. / coef *statsf(fg).sum, 1./3.);

  double rsol = 0.0;
  double relerr = (rsol > 0.) ? fabs (rsol - effective_radius) / rsol : 0.;

  char name[80];
  sprintf (name, "data/OutputData-%d", maxlevel);

  static FILE * fp = fopen (name, "w");
  fprintf (fp, "%g %g %g\n",
      t+tshift, effective_radius, theta0);
  fflush (fp);
}


event movie(t += final_time / (double)output_period)
{

  vector ufc[];
  vector ufsavec[];
  vector ufextc[];

  foreach()
  {
    T[] = (int)is_solid[] == 1 ? TS[] : T[];
    foreach_dimension()
    {
      ufextc.x[] = (ufext.x[] + ufext.x[1]) / (fm.x[] + fm.x[1])* (1.0 - is_solid[]);
      ufc.x[] = (uf.x[] + uf.x[1]) / (fm.x[] + fm.x[1])* (1.0 - is_solid[]);
      u.x[] = u.x[] * (1.0 - is_solid[]);

#if !VEL_POTENTIAL && !VELOCITY_JUMP
      uext.x[] = uext.x[] * (1.0 - is_solid[]);
#endif

#if VELOCITY_JUMP
      u1.x[] = u1.x[] * (1.0 - is_solid[]);
      u2.x[] = u2.x[] * (1.0 - is_solid[]);
      u1.x[] = f[] > 0.5 ? u1.x[] : u2.x[];

      ufextc.x[] = (uf2.x[] + uf2.x[1]) / (fm.x[] + fm.x[1]) * (1.0 - is_solid[]);
      ufc.x[] = (uf1.x[] + uf1.x[1]) / (fm.x[] + fm.x[1]) * (1.0 - is_solid[]);
      ufc.x[] = f[] > 0.5 ? ufc.x[] : ufextc.x[];
#endif
    }
  }

  foreach()
  {
    T[] = (int)is_solid[] == 1 ? TS[] : T[];
  }
  char filename[20];
  sprintf(filename, "data/dump-%05d", dump_num++);
  dump(filename);
}

event finalResult(t = final_time)
{
  dump("dump_end");
  char filename[20];
  sprintf(filename, "data/interface-%d", Ncell);
  outputFacetsParallel(filename);
}
