/**
# Test of gas injection in a box filled with liquid with tracer and oil.
*/
#define FCORRECT 0
/*
   We use the centered Navier--Stokes solver and log performance
   statistics. */
#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "three-phase.h"
#include "tension.h"
#include "view.h"

/**
  The density ratio and the dynamic viscosity ratio are defined */
#define RHOR 10.
#define MUR 54.

// Physical parameters of fluids 20°C
#define RHOL 1000.0
#define RHOG (RHOL/RHOR)
#define MUL 1.002e-3
#define MUG (MUL/RHOR)
#define RHOO (920.0)
#define MUO (0.079)
#define SIGMAGL (72e-3)
#define SIGMALO (25.5e-3)
#define SIGMAGO (31.7e-3)

/**
  We can control the maximum runtime. */

#define WIDTH 0.27
#define Ho 0.2
#define Ha 0.20667
#define VOFTHR 1e-10                                   
#define MAXTIME 0.5

//Geometric parameters of the jet
#define D 7.9e-3
#define RADIUS (D/2.)
#define LENGTH (2*RADIUS/3.3)
#define QG 1e-5
#define GRAV 9.81
#define AREAI (pi*sq(RADIUS))
#define VELOCITY (QG/AREAI)

scalar T[], f10[];
int maxlevel = 8;

//Boundary conditions
u.n[bottom]  = dirichlet(f10[]*VELOCITY);
u.t[bottom]  = dirichlet(0);
f1[bottom]   = f10[];
f2[bottom]   = 1.-f10[];
T[bottom]   = 1.-f10[];
f3[bottom]   = 0.;

u.n[top] = u.n[] > 0. ? neumann(0) : 0;
p[top]   = dirichlet(0);
f1[top]   = 1.;
f2[top]   = 0.;
f3[top]   = 0.;

/**
  The main function can take two optional parameters: the maximum level
  of adaptive refinement (as well as an optional maximum runtime). */

int main (int argc, char * argv[]) {
  if (argc > 1)
    maxlevel = atoi (argv[1]);

  /**
    We set the domain geometry and initial refinement. */
  size (WIDTH);
  origin (-L0/2, 0, -L0/2);
  init_grid (1 << maxlevel);
  /**
    We set the physical parameters: densities, viscosities and surface
    tension. */
  rho1 = RHOG;
  rho2 = RHOL;
  rho3 = RHOO;
  mu1 = MUG;
  mu2 = MUL;
  mu3 = MUO;
  f1.sigma = (SIGMAGO+ SIGMAGL- SIGMALO)/2.;
  f2.sigma = (-SIGMAGO+ SIGMAGL+ SIGMALO)/2.;
  f3.sigma = (SIGMAGO- SIGMAGL+ SIGMALO)/2.;
  f2.tracers = (scalar *){T};
  run();
}

/**
  For the initial conditions, we first try to restore the simulation
  from a previous "restart". */

event init (t = 0) {
  if (!restore (file = "restart")) {
    scalar fbis[];
    foreach()
      fbis[] = 1.*(y >= Ha);
    boundary ({fbis});
    fraction (f3, (y-Ho)*(Ha-y)); 
    fraction (f10, sq(RADIUS) - sq(x) - sq(z));
    foreach() {
      f1[] = ((y < LENGTH)? f10[] :fbis[]);
      u.y[] = f1[];
    }
    boundary ({f1,u.y});
    foreach() {
      f2[] = 1.-(f1[]+f3[]);
      T[] = 1.-(f1[]+f3[]);
    }
    boundary ({f2,T});
  }
}

#if FCORRECT
event fcorrect (i++) {
  static FILE * fp1 = fopen ("fluids.dat", "w");
  static FILE * fp3 = fopen ("injection.dat", "w");
  double summ = 0.;
  foreach() {
    summ = (f1[]+f2[]+f3[]);
    if (summ > 0.) {
      f1[] = clamp(f1[]/summ,0.,1.);
      f2[] = clamp(f2[]/summ,0.,1.);
      f3[] = clamp(f3[]/summ,0.,1.);
    }
  }
  boundary ({f1,f2,f3});
  stats r = statsf (f1);
  stats s = statsf (f2);
  stats o = statsf (f3);
  fprintf (fp1, "%g %g %g %g %g %g %g %g %g %g\n", t, r.min, r.max, r.sum, s.min, s.max, s.sum, o.min, o.max, o.sum);
  fflush (fp1);
}
#endif

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= GRAV;
  if (i <= 30) {
    TOLERANCE = 1e-6;
  }
  TOLERANCE = 1e-4;
}

event logfile (i++) {
  scalar c[],sum[];
  foreach() {
    sum[] = f1[]+f2[]+f3[];
    c[] = (f2[] > VOFTHR ? T[]/f2[] : 1.);
  }
  boundary({sum,c});
  stats s = statsf (c);
  stats p = statsf (T);
  stats q = statsf (sum);
  double Tquantf1 = 0., cquantf1 = 0., Tquantf2 = 0., cquantf2 = 0., Tquantf3 = 0., cquantf3 = 0., vol1 = 0., vol2 = 0., vol3 = 0.;
  foreach(reduction(+:Tquantf1) reduction(+:cquantf1) reduction(+:Tquantf2) reduction(+:cquantf2) reduction(+:Tquantf3) reduction(+:cquantf3) reduction(+:vol1) reduction(+:vol2) reduction(+:vol3)) {
    if (f1[] >= (1.-VOFTHR) && f2[] <= VOFTHR && f3[] <= VOFTHR) {
      Tquantf1 += T[]*dv();
      cquantf1 += c[]*dv();
      vol1 += dv();
    }
    else if (f2[] >= (1.-VOFTHR) && f1[] <= VOFTHR && f3[] <= VOFTHR) {
      Tquantf2 += T[]*dv();
      cquantf2 += c[]*dv();
      vol2 += dv();
    }
    else if (f3[] >= (1.-VOFTHR) && f1[] <= VOFTHR && f2[] <= VOFTHR) {
      Tquantf3 += T[]*dv();
      cquantf3 += c[]*dv();
      vol3 += dv();
    }
  }
  fprintf (ferr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
      t, vol1, vol2, vol3, Tquantf1/vol1, cquantf1/vol1, Tquantf2/vol2, cquantf2/vol2, Tquantf3/vol3, cquantf3/vol3, s.min-1., s.max-1., s.sum/s.volume-1., p.min, p.max, p.sum/p.volume, q.min, q.max, q.sum/q.volume);
}

event end (t = MAXTIME)
{
  scalar c[];
  foreach()
    c[] = (f2[] > VOFTHR ? T[]/f2[] : 1.);
  boundary({c});
view (fov = 24, quat = {0,0,0,1}, tx = 0.043206, ty = -0.439261, bg = {1,1,1}, width = 600, height = 600);
  box();
  draw_vof ("f1");
  draw_vof ("f3");
  squares ("c", spread = -1);
  save ("c.png");
}

/**
  Every time unit, we output a full snapshot of the simulation, to be
  able to restart and for visualisation. In three dimension. */

#if 0
event snapshot (t = 0; t +=0.1; t <= MAXTIME) {
scalar c[];
  foreach()
    c[] = (f2[] > VOFTHR ? T[]/f2[] : 1.);
  boundary({c});
  char name[80];
  sprintf (name, "snap-%g", t);
  dump (file = name);
}
#endif

/**
## Results

![c field (concentration field) at end of simulation ](bullage3f_mwe/cpatch.png)

![c field (concentration field) at end of simulation with normalization](bullage3f_mwe/cpatchfcor.png)

![maximum, average and minimum deviation of the concentration value in function of time](bullage3f_mwe/cdeviation.svg)

![maximum, average and minimum deviation of the concentration value in function of time with normalization](bullage3f_mwe/cdeviationfcor.svg)

![maximimum, average and minimum of the sum of fluid fraction value in function of time](bullage3f_mwe/sumdeviation.svg)

![maximum, average and minimum of the sum of fluid fraction in function of time with normalization](bullage3f_mwe/sumdeviationfcor.svg)

*/
