#define HYDRO 0

#include "grid/multigrid1D.h"

// redefine a_baro barotropic acceleration

scalar panom[];
// Adding in pressure anomaly in hectopascals (x100 pascals)
# define a_baro(eta, i) (gmetric(i)*(G*(eta[i-1] - eta[i])+sq(60.)*(100./1030.)*(panom[i-1]-panom[i]))/Delta)

#if GN
# include "green-naghdi.h"
# define MGD mgD.i
#else
# include "layered/hydro.h"
# if HYDRO
#   define MGD 0
# else
#   include "layered/nh.h"
#   define MGD mgp.i
# endif
# include "layered/perfs.h"
scalar h;
vector u;
#endif
#include "output.h"
#include "pressure_force1D.h"

#define MAXLEVEL 12
#define MINLEVEL 5
#define INITLEVEL 10
#define ETAE     2e-2 // error on free surface elevation (5 cm)
#define ETAMAXE  6e-2 // error on maximum free surface elevation (10 cm)
#define EXTENT 20000.
#define X_OFFSET -10000.
#define DEPTH 10.
#define PRDMN 1.0  // Specify the Froude number here
#define ATM_MAG 1.
#define ATM_R 100.

double t_end=20.;
double t_out=0.1;
double atm_spd;

int main()
{
  // the domain is 10 km 
  size (EXTENT);
  origin(X_OFFSET);
  // acceleration of gravity in m/min^2
  G = 9.81*sq(60.);
  atm_spd=PRDMN*sqrt(G*DEPTH);
  // 32^2 grid points to start with
#if TREE
  init_grid (1 << MAXLEVEL);
#else
  //N = 512; // 2^9
  //N = 1024; // 2^10
  //N = 2048; // 2^11
  N = 4096; // 2^12
  //N = 8192; // 2^13
#endif

  /**
  For the Green-Naghdi solver we do only one iteration of the
  multigrid solver to speed things up. This does not change the
  solution.

  In the non-hydrostatic case, we add a cut-off breaking parameter
  ([Popinet, 2020](/src/references.bib#popinet2020)). */
  
#if GN
  TOLERANCE = HUGE;
  NITERMIN = 1;
#else // !GN
#if NH
  CFL_H = 0.5; // this is necessary for stability at level 13
  //breaking = 0.07;
#endif
#endif // !GN
  DT = 0.1;
  run();
}

scalar etamax[];

//// ## Initial conditions


event init (i = 0)
{
  if (restore (file = "restart"))
    conserve_elevation();
  else {
    conserve_elevation();
    foreach() {
      foreach_layer() {
        h[] = (DEPTH) / nl;
      }
      panom[] = ATM_MAG*exp(-(x-atm_spd*t)*(x-atm_spd*t)/(ATM_R*ATM_R));
      etamax[]=0.;
    }


  }

  u.n[left]   = dirichlet(0.);
  u.n[right]  = + radiation(DEPTH-100./1030./G*ATM_MAG*exp(-(EXTENT-atm_spd*t)*(EXTENT-atm_spd*t)/(ATM_R*ATM_R)));
}

/**
## Quadratic bottom friction */

event etamax_calc (i++)
{
  foreach() {
    panom[] = ATM_MAG*exp(-(x-atm_spd*t)*(x-atm_spd*t)/(ATM_R*ATM_R));
    if (h[] > dry && h[] + zb[] > etamax[])
      etamax[] = h[] + zb[];
  }
}

event friction (i++)
{
  foreach() {
    double a_inv = h[] < dry ? 0. : h[]/(h[] + 1e-4*dt*norm(u));
    foreach_dimension()
      u.x[] *= a_inv;
  }
  if (i > 10 && dt < 1.e-7){
    fprintf(stderr,"Time-step becoming too small, t = %g, dt = %g\n",t,dt);
    char name[80];
    sprintf (name, "instability.txt", t);
    FILE *fp=fopen(name,"w");
    foreach(){
      fprintf(fp,"%g %g %g %g %g\n",x,eta[],etamax[],panom[],u.x[]);
    }
    exit(42);
  }
}

/**
## Outputs */

event logfile (i++) {
  stats s = statsf (h);
  norm n = normf (u.x);
  if (i == 0)
    fprintf (stderr,
	     "t i h.min h.max h.sum u.x.rms u.x.max dt mgD.i speed tn\n");
  fprintf (stderr, "%g %d %g %g %g %g %g %g %d %g %ld\n",
	   t, i, s.min, s.max, s.sum, n.rms, n.max, dt, MGD,
	   perf.speed, grid->tn);
}

event snapshots (t = t_out; t += t_out) {
  char name[80];
  sprintf (name, "dump-%g", t);
  dump (name);
}

event outputN (t=t_out; t+=t_out; t<=t_end) {
  char name[80];
  sprintf (name, "results-%g.txt", t);
  FILE *fp=fopen(name,"w");
  foreach(){
    fprintf(fp,"%g %g %g %g %g\n",x,eta[],etamax[],panom[],u.x[]);
  }
}