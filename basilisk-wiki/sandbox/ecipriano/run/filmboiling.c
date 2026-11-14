/**
# Film boiling

The film boiling configuration consists in a superheated solid wall, covered by
a thin vapor layer, which is perturbed by the phase change phenomena releasing
bubbles. The domain is initialized with a uniform constant temperature equal to
the saturation value. The bottom wall is maintained at a specific superheating
temperature, causing the vapor layer to expand and the perturbation triggers
Rayleigh-Taylor instability, leading to the formation, rising, and detachment of
the bubble. Notably, this test case combines the boiling model with embedded
boundaries and AXI metrics.

![Evolution of the temperature (left) and the velocity field (right)](filmboiling/movie.mp4)
*/

#include "embed.h"
#include "axi.h"
#include "navier-stokes/low-mach.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "boiling.h"
#include "view.h"

/**
Outflow boundary conditions at the top of the liquid column, while no-slip
condition on the surface of the solid wall. Symmetry elsewhere. */

u.n[left] = dirichlet (0.);
u.t[left] = dirichlet (0.);
p[left] = neumann (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

/**
We multiply by `cs` in order to limit the number of cells on the solid
(embed) phase. */

double Twall, Tsat;
T[left] = dirichlet (cs[]*Twall);

int maxlevel, minlevel = 6;
double lambdaw, gasvolume0;

int main (void) {
  rho1 = 200., rho2 = 5.;
  mu1 = 0.1, mu2 = 0.005;
  lambda1 = 40., lambda2 = 1.;
  cp1 = 400., cp2 = 200.;
  dhev = 1e+4;

  Tsat = 1., Twall = Tsat + 5.;
  TIntVal = TG0 = TL0 = Tsat;

  f.sigma = 0.1;
  G.x = -9.81;
  lambdaw = 2.*pi*sqrt (3.*f.sigma / (fabs (G.x)*(rho1 - rho2)));

  nv = 1;
  pcm.consistent = true;

  TOLERANCE = 1.e-6 [*];
  size (lambdaw);
  for (maxlevel = 6; maxlevel <= 8; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

#define sin(x,y) lambdaw/128.*(4. + cos(2.*pi*y/lambdaw))

event init (i = 0) {
  solid (cs, fs, -(y - 0.5*lambdaw));
  fractions_cleanup (cs, fs);
  fraction (f, x - sin(x,y));
  foreach()
    f[] = y <= lambdaw ? f[] : 0.;
#if defined(AXI) && defined(EMBED)
  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);
  restriction ({cs,fs,cm,fm});
#endif

  scalar TL = liq->T, TG = gas->T;
  foreach() {
    double xmax = sin(x,y);
    TG[] = clamp (Twall - x/xmax*(Twall - Tsat), Tsat, Twall)*(1. - f[]);
    TL[] = TL0*f[];
    TG[] *= (cs[] != 0.);
    TL[] *= (cs[] != 0.);
    T[] = TL[] + TG[];
  }
  copy_bcs ({TL,TG}, T);

  gasvolume0 = 0.;
  foreach (reduction(+:gasvolume0))
    gasvolume0 += (1. - f[])*dv();
}

/**
We adapt according to the interface position, the temperature, and velocity
field.  The simulation does not work if the interface is not maintained at the
maximum level of refinement. It would be nice to try and go deeper in this
problem. */

event adapt (i++) {
#if 1
  vector gf[];
  gradients ({f}, {gf});
  scalar mgf[];
  foreach()
    mgf[] = norm (gf);
  double mgfmax = statsf (mgf).max;
  foreach()
    mgf[] /= mgfmax;

  adapt_wavelet ({mgf,T,u.x,u.y}, {1e-2,1e-2,1e-1,1e-1}, maxlevel, minlevel);
#else
  adapt_wavelet_leave_interface ({T,u.x,u.y}, {f},
      (double[]){1e-2,1e-2,1e-2}, maxlevel, minlevel, 1);
#endif

#if defined(AXI) && defined(EMBED)
  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);
  restriction ({cs,fs,cm,fm});
#endif
}

/**
## Post-processing

Printing the interface at a specific time. */

event facets (t = 0.32) {
#if _MPI
  char name[80];
  sprintf (name, "facets-pid-%d", pid());
  FILE * fp = fopen (name, "w");

  scalar ff[];
  foreach()
    ff[] = cm[] ? f[] : 0.;
  output_facets (ff, fp);
  fclose (fp);

  MPI_Barrier (MPI_COMM_WORLD);

  if (pid() == 0) {
    char command[80];
    sprintf (command, "cat facets-pid-* > facets-%d && rm facets-pid-*",
        maxlevel);
    system (command);
  }
#else
  char name[80];
  sprintf (name, "facets-%d", maxlevel);
  FILE * fp = fopen (name, "w");

  scalar ff[];
  foreach()
    ff[] = cm[] ? f[] : 0.;
  output_facets (ff, fp);

  fclose (fp);
#endif
}

/**
Using MPI, we can verify the distribution of the cells among the different
processors. */

#if _MPI
scalar pid[];
event pids (i++) {
  foreach()
    pid[] = pid();
}
#endif

/**
We compute and write the space-averaged Nusselt number and the normalize gas
volume in time. */

event logfile (t += 0.001) {
  double gasvolume = 0.;
  foreach (reduction(+:gasvolume))
    gasvolume += (1. - f[])*dv();

  scalar TL = liq->T, TG = gas->T;
  foreach()
    T[] = TL[] + TG[];
  boundary ({T}); // foreach_boudary() does not trigger automatic bc

  double ld2 = sqrt (f.sigma/(rho1 - rho2)/fabs (G.x));
  double Nu = 0.;
  foreach_boundary (left, reduction(+:Nu))
    Nu += (T[] - T[-1,0])*cm[];
  Nu *= -ld2*2./(Twall - Tsat)/sq (0.5*lambdaw);

  fprintf (stderr, "level %d %g %g %g\n", maxlevel, t,
      gasvolume / gasvolume0, Nu), fflush (stderr);
}

/**
We write a movie with the evolution of the interface position, the temperature,
and the velocity fields. */

event movie (t += 0.002; t <= 0.35) {
  if (maxlevel == 8) {
    scalar ff[];
    foreach()
      ff[] = cm[] ? f[] : 0.;
    clear();
    view (ty = -0.5, psi=-pi/2.);
    draw_vof ("cs", "fs", filled = -1, fc = {1.,1.,1.});
    draw_vof ("ff", lw = 2.);
    squares ("T", min = Tsat, max = Twall);
    mirror ({0.,1.}) {
      draw_vof ("cs", "fs", filled = -1, fc = {1.,1.,1.});
      draw_vof ("ff", lw = 2.);
      squares ("(u.x^2 + u.y^2)^0.5", spread = -1);
    }
    save ("movie.mp4");
  }
}

/**
## Results

~~~gnuplot Convergence of the interface shape
set term push
set size ratio -1
unset xtics
unset ytics
unset border

plot "facets-6" u (-$2):($1) w l lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "facets-7" u (-$2):($1) w l lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "facets-8" u (-$2):($1) w l lw 1.2 lc -1 dt 1 t "LEVEL 8", \
     "facets-6" u ($2):($1)  w l lw 1.2 lc -1 dt 3 notitle, \
     "facets-7" u ($2):($1)  w l lw 1.2 lc -1 dt 2 notitle, \
     "facets-8" u ($2):($1)  w l lw 1.2 lc -1 dt 1 notitle
set term pop
~~~

~~~gnuplot Evolution of the space-averaged Nusselt number in time
reset
set xlabel "time [s]"
set ylabel "<Nu> [-]"
set grid
set xr[0:0.35]
#set yr[0:10]

plot "<grep 'level 6' log" u 3:5 w l lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' log" u 3:5 w l lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' log" u 3:5 w l lw 1.2 lc -1 dt 1 t "LEVEL 8"
~~~

~~~gnuplot Evolution of the normalized gas volume in time
reset
set xlabel "time [s]"
set ylabel "V_g / V_g^0 [-]"
set grid
set key left top

plot "<grep 'level 6' log" u 3:4 w l lw 1.2 lc -1 dt 3 t "LEVEL 6", \
     "<grep 'level 7' log" u 3:4 w l lw 1.2 lc -1 dt 2 t "LEVEL 7", \
     "<grep 'level 8' log" u 3:4 w l lw 1.2 lc -1 dt 1 t "LEVEL 8"
~~~
*/
