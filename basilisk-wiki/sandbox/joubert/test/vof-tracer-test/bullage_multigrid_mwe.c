/**
# Test of gas injection in a box filled with liquid wich contain a VOF tracer.
*/

/*
   We use the centered Navier--Stokes solver and log performance
   statistics. */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

/**
  The density ratio and the dynamic viscosity ratio are defined */
#define RHOR 10.
#define MUR 54.

// Physical parameters 
#define RHOL 1000.0
#define RHOG (RHOL/RHOR)
#define MUL 1.002e-3
#define MUG (MUL/RHOR)
#define SIGMAGL 72e-3

/**
  We can control the maximum runtime. */

#define WIDTH 0.27
#define Ha 0.21
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

scalar T[], f0[], fbis[];
int maxlevel = 8;

//Boundary conditions
u.n[bottom]  = dirichlet(f0[]*VELOCITY);
u.t[bottom]  = dirichlet(0.);
f[bottom]   = f0[];
T[bottom]   = (1.-f0[]);

u.n[top] = u.n[] > 0. ? neumann(0.) : 0.;
p[top]   = dirichlet(0.);
f[top]   = 1.;

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
  mu1 = MUG;
  mu2 = MUL;
  f.sigma = SIGMAGL;
  T.inverse = true;
  T.gradient = zero;//minmod2;
  f.tracers = (scalar *){T};
  run();
}

/**
  For the initial conditions, we first try to restore the simulation
  from a previous "restart". */

event init (t = 0) {
  if (!restore (file = "restart")) {
    foreach()
      fbis[] = 1.*(y >= Ha);
    boundary ({fbis});
    fraction (f0, sq(RADIUS) - sq(x) - sq(z));
    foreach() {
      f[] = (y < LENGTH)?f0[]:fbis[];
      u.y[] = f[];
      T[] = 1.-f[];
    }
    boundary ({f,u.y,T});
  }
}

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
  scalar c[];
  foreach()
    c[] = ((1.-f[]) > VOFTHR ? T[]/(1.-f[]) : 1.);
  stats s = statsf (c);
  double fquantf1 = 0., Tquantf1 = 0., cquantf1 = 0., fquantf2 = 0.,
	 Tquantf2 = 0., cquantf2 = 0., vol1 = 0., vol2 = 0.;
  foreach(reduction(+:fquantf1) reduction(+:Tquantf1) reduction(+:cquantf1)
      reduction(+:fquantf2) reduction(+:Tquantf2) reduction(+:cquantf2) 
      reduction(+:vol1) reduction(+:vol2)) {
    if (f[] >= (1.-VOFTHR)) {
      fquantf1 += f[]*dv();
      Tquantf1 += T[]*dv();
      cquantf1 += c[]*dv();
      vol1 += dv();
    }
    else if (f[] <= VOFTHR) {
      fquantf2 += f[]*dv();
      Tquantf2 += T[]*dv();
      cquantf2 += c[]*dv();
      vol2 += dv();
    }
  }
  fprintf (ferr, "%g %g %g %g %g %g %g %g %g %g %g %g\n", 
      t, vol1, vol2, fquantf1/vol1, Tquantf1/vol1, cquantf1/vol1-1.,
      fquantf2/vol2, Tquantf2/vol2, cquantf2/vol2-1., s.min-1.,
      s.max-1., s.sum/s.volume-1.);
}

event end (t = MAXTIME)
{
  scalar c[];
  foreach()
    c[] = ((1.-f[]) > VOFTHR ? T[]/(1.-f[]) : 1.);
view (fov = 24, quat = {0,0,0,1}, tx = 0.043206, ty = -0.439261, bg = {1,1,1}, width = 600, height = 600);
  box();
  draw_vof ("f");
  squares ("c", spread = -1);
  save ("c.png");
}

/**
  Every time unit, we output a full snapshot of the simulation, to be
  able to restart and for visualisation. In three dimension. */
#if 0 //debug
event snapshot (t = 0; t +=0.1; t <= MAXTIME) {
scalar c[];
  foreach()
    c[] = ((1.-f[]) > VOFTHR ? T[]/(1.-f[]) : 1.);
  boundary({c});
  char name[80];
  sprintf (name, "snap-%g", t);
  dump (file = name);
}
#endif

/**
## Results

![Concentration field at $t = 0.5$ without patch](bullage_multigrid_mwe/c.png)

~~~gnuplot maximimum, average and minimum deviation of the concentration value in function of time without the patch
reset
set xlabel 'time [s]'
set ylabel 'concentration-1 quantity []'
plot 'log' u 1:11 w l t 'maximum',\
'log' u 1:12 w l t 'average',\
'log' u 1:10 w l t 'minimum'
~~~

~~~gnuplot VOF fraction in function of time in the domain
reset
set xlabel 'time [s]'
set ylabel 'f quantity in gas [-]'
plot 'log' u 1:4 w l t '' lw 2
~~~

~~~gnuplot VOF fraction in function of time in the domain
reset
set xlabel 'time [s]'
set ylabel 'f quantity in liquid [-]'
plot 'log' u 1:7 w l t '' lw 2
~~~

All the results below are obtained with the vof-tracer patch.

![Concentration field at $t = 0.5$ with patch](bullage_multigrid_mwe/cpatch.png)

![maximimum, average and minimum deviation of the concentration value in function of time
with patch](bullage_multigrid_mwe/c_deviation.svg)

![VOF tracer quantity in function of time in the gas phase with patch](bullage_multigrid_mwe/tgas.svg)

![Deviation compare to 1 of concentration quantity in function of time in the gas phase with patch](bullage_multigrid_mwe/cgas.svg)

![VOF tracer quantity in function of time in the liquid phase with patch](bullage_multigrid_mwe/tliq.svg)

![Deviation compare to 1 of concentration quantity in function of time in the liquid phase with patch](bullage_multigrid_mwe/cliq.svg)

*/
