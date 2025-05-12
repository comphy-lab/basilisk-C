
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
int output_period = 100;

const double uemax = 1e-2;
const double Temax = 1e-2;
const double femax = 1e-6;
const double film_height = 15e-3;   //this will not be used in this case
const double remove_height = 8.0e-3; //this will not be used in this case

const int num_refine = 1;

double Lsize = 2.3e-3;

const double R0 = 20e-6;
const double final_time = 0.5e-3;
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
  double leny = ncs * delta;
  if (pid() == 0.0)
  {
    printf("solid_len_x: %g\n", lenx);
    printf("solid_len_y: %g\n", leny);
  }
  SOLID_LEN_x = lenx;
  SOLID_LEN_y = leny;
}

int Ncell = 64;

double delta_q_case(double x, double t)
{
  
  double grad = 425e3;

  // double grad = (0.004 - x) / 0.004 * 481.0;
  // grad = max(grad * 1000.0, 0.0);

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
      theta0 = atof(argv[2]);
    }

    if (argc == 4)
    {
      maxlevel = atoi(argv[1]);
      theta0 = atof(argv[2]);
      output_period = atoi(argv[3]);
    }

    if (argc == 5)
    {
      maxlevel = atoi(argv[1]);
      theta0 = atof(argv[2]);
      output_period = atoi(argv[3]);
      is_restart = 1;
      dump_num = atoi(argv[4]);
    }
  }

#if _MPI
  MPI_Bcast(&maxlevel, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&theta0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&is_restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&output_period, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dump_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
    double omega = 0.0460; // 0.0345;//0.0460;
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
    sprintf(dump_name, "dumpTini");
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
}

event beforePhaseChange(i++)
{
  remove_droplets_phasechange (f, threshold = 1.0e-20, minsize = 5);
  remove_droplets_phasechange (f, threshold = 1.0e-20, minsize = 5, bubbles = true);

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
    av.x[] -= 9.81;

  //set dynamic contact angle
  double nanoscale = 10e-9;
  double Cae = 0.0;
  double thetahock = 0.0;
  double deltat_nucl = 12.55;
  double dx = L0 / (1 << maxlevel);
  Cae = mu1 * deltat_nucl/rho1/dhev/Rcc/f.sigma;
  thetahock = pow(12*Cae*log(0.5*dx/nanoscale),0.25);

  double vel = 0.0;
  //compute contact line velocity
  foreach(reduction (max:vel))
  {
    if(is_solid[] == 0.0 && is_solid[-1] == 1.0)
    {
      bool is_cut = f[] > 0.0 && f[] < 1.0;
      if(is_cut)
      {
        vel = fabs(u.y[]);
      }
    }
  }


  double clvel = vel;

  double Cau = mu1 * max(1.0e-12, clvel) / f.sigma;
  double cangnew = min(Cae / Cau, thetahock) * 180. / pi;
  double cangle = min(Cae / Cau, thetahock);
  //theta0 = cangnew;
  if(pid() == 0)
  {
    //printf("Cae = %g thetahock = %g cvel = %g theta0 = %g Rcc = %g\n", Cae, thetahock, clvel, cangnew, Rcc);
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
  if(pid() == 0 && i % 1 == 0)
  {
    printf("i = %d t = %g dt = %g\n", i, t, dt);
  }
}

// event calInitMicroLayer(t += 0.5e-6)
// {
//   scalar mt[];
//   scalar tang[];
//   int pos_j = 0;
//   double thick = 0.0;
//   double pos = 0.0;

//   foreach()
//   {
//     mt[] = 0.0;
//     tang[] = 0.0;
//     if((int)is_solid[] == 0 && f[] > 0.0 && f[] < 1.0)
//     {
//       coord nc = normal(point, f);
//       tang[] = fabs(nc.y / nc.x);
//     }
//   }

//   for(int it = 0; it < 5; ++it)
//   {
//     foreach()
//     {
//       if((int)is_solid[] == 1)
//         continue;

//       if (tang[] == 0.0 && tang[1] != 0.0)
//       {
//         tang[] = tang[1];
//       }
//       else if (tang[] == 0.0 && tang[2] != 0.0)
//       {
//         tang[] = tang[2];
//       }
//     }
//   }

//   foreach()
//   {
//     if((int)is_solid[] == 0 && (int)is_solid[-1] == 1)
//     {
//       double threold = 0.015;//tan(3 * pi / 180.0);
//       if (mt[] == 0.0)
//       {
//         if(tang[] > 0.0 && tang[] < threold && h.x[] != nodata)
//         {
//           mt[] = (h.x[] + 0.5) * Delta;
//         }
//       }
//     }
//   }
//   double delta = Lsize / (1<<grid->maxdepth);
//   foreach(reduction(max:pos_j))
//   {
//     if((int)is_solid[] == 0 && (int)is_solid[-1] == 1)
//     {
//       if(mt[] > 0.0)
//       {
//         if(point.j > pos_j)
//         {
//           pos_j = point.j;
//         }
//       }
//     }
//   }

//   foreach (reduction(max : thick) reduction(max : pos))
//   {
//     if ((int)is_solid[] == 0 && (int)is_solid[-1] == 1)
//     {
//       if (point.j == pos_j)
//       {
//         thick = mt[];
//         pos = y;
//       }
//     }
//   }
//   if(pos > pos_old)
//   {
//     pos_old = pos;
//     thick_old = thick;
//   }
 
//   // if(t > 0.0)
//   // {
//   //   if(pid() == 0.0)
//   //   {
//   //     printf ("%g %g %g %g\n", t+tshift, pos_old / delta, thick_old, thick_old / delta);
//   //   }
//   //   dump();
//   //   exit(1);
//   // }
//   char name[80];
//   sprintf (name, "data/InitMicroLayer-%d", maxlevel);

//   if(pid() == 0)
//   {
//     static FILE * fp = fopen (name, "w");
//     fprintf (fp, "%g %g %g %g\n", t+tshift, pos_old, thick_old, thick_old / delta);
//     fflush (fp);
//   }
// }

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
  fprintf (fp, "%g %g\n",
      t+tshift, effective_radius);
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

  char filename[20];
  sprintf(filename, "data/dump-%05d", dump_num++);
  dump(filename);
}

event finalResult(t = final_time)
{
  //do nothing
}
