/**
# Vortex Ejection from a mode 3 instability

According to [Kizner et
al. (2013)](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/instabilities-of-the-flow-around-a-cylinder-and-emission-of-vortex-dipoles/2327C2CFA76059D27461A2BC6A09F146),
a flow around a no-slip cylinder with radius $R$ maybe unstable and eject three
dipolar vortex pairs. We study the flow using embedded boundaries and
the Navier-Stokes solver with the double-projection scheme. Furthermore, we will `view` our results.
*/
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "view.h"
/**
The maximum resolution is set to $\Delta_{min}=R/100$. This allows to run on the sandbox server (see `[Pushing the resolution](kizner.c#note-on-pushing-the-resolution.)' section).   
*/
int maxlevel = 12;
double Re = 30000.;
face vector muc[];

int main(){
  init_grid(64);
  L0 = 40.;
  mu = muc;
  X0 = Y0 = -L0/2.;
  /**
Rather than choosing stress-free outer-domain boundaries, periodic
boundaries are used. This markably increased the congergence
properties of the iterative Multigrid strategy applied to the
advection and viscous problems.
   */
  foreach_dimension()
    periodic(left);
  run();
}

event properties (i++){
  foreach_face()
    muc.x[] = fm.x[]/Re;
  boundary ((scalar*){muc});
}
/**
The cylinder is defined and the flow field is initialized c.f. Kizner
et al. (2013) with an $m=3$ perturbation.
 */
#define RAD (pow(sq(x) + sq(y), 0.5))
#define THETA(M) (M*asin(x/RAD))
#define RADP(P, M) ((1 + P*sin(THETA(M)))/(pow(1 + 0.5*sq(P), 0.5)))
double a1 = 1.5, b1 = 2.25;
double P = 0.005, m = 3;

event init(t = 0){
  double gamma = (sq(a1) - 1.)/(sq(b1) - sq(a1));
  refine (RAD < b1  && level < (maxlevel - 1));
  refine (RAD < 1.05 && RAD > 0.95  && level < maxlevel);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = RAD - 1;
  fractions(phi, cs, fs);
  foreach(){
    double r = RAD;
    double r1 = RADP(P,m)*r;
    double vr;
    if (r1 > 0.9 && r1 < a1)
      vr = r1 - 1./r1;
    else if (r1 >= a1 && r1 <= b1)
      vr = -gamma*r1 + ((1 + gamma)*sq(a1) - 1.)/r1;
    else // (0.9 > r || r > b)
      vr = 0;
    u.x[] = cm[]*0.5*vr*r*-y/(sq(x) + sq(y));
    u.y[] = cm[]*0.5*vr*r* x/(sq(x) + sq(y));
  }
  /**
The boundary conditions on the embedded boundary are set:
  */
  u.t[embed] = dirichlet (0.);
  u.n[embed] = dirichlet (0.);
  boundary (all);
  /**
Since the perturbed initialized solution is slightly inconsistent, the
Poisson solver is tuned to be extra robust for the first ten
timesteps.
   */
  CFL = 0.6;
  DT = 0.02;
  TOLERANCE = 1E-4;
  NITERMIN = 5;
}

event relax_a_little (i = 10)
  NITERMIN = 1;

/**
The grid is adaptedively refined and coarsened to properly represent
the boundary and the evolution of the flow field. We set
$\zeta_{u,v}\approx U_{max}/100$.
*/
event adapt (i++)
  adapt_wavelet ({cs, u.x, u.y}, (double[]){0.01, 0.004, 0.004}, maxlevel);

/**
## Ouput and Results

Movies are generated that display the vorticity dynamics and the
grid-cell structure. 

<video width="600" height="600" controls>
<source src="http://www.basilisk.fr/sandbox/Antoonvh/kizner12.mp4" type="video/mp4">
</video> 
<video width="600" height="600" controls>
<source src="http://www.basilisk.fr/sandbox/Antoonvh/kizner_cells12.mp4" type="video/mp4">
</video> 
 */
event movie (t += 0.4; t <= 100){
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  view (fov = 7, width = 600, height = 600, samples = 1);
  clear();
  draw_vof ("cs", filled = -1, fc = {1., 1., 1.});
  squares ("omega", min = -0.75, max = 0.75,
	   map = cool_warm, linear = true);
  draw_vof ("cs", "fs", lw = 2);
  save ("kizner12.mp4");
  clear();
  cells();
  save("kizner_cells12.mp4");
}
/**
Furthermore, we log the number of grid cells over time.
 */
event logger(t += 0.1){
  int cells = 0;
  foreach()
    cells++;
  printf("%g\t%d\t%d\t%g\t%g\t%d\t%d\t%d\n",
	 t, i, cells, dt, DT, mgp.i, mgpf.i, mgu.i);
}
/**
~~~gnuplot These numbers may be compared against the millions of cells that Kizner et al. (2013) employed.
set yr [0 : 22000]
set xlabel 'time'
set ylabel 'Cells'
set key off
plot 'out' u 1:3 w l lw 2
~~~

## Note on pushing the resolution.

A run with 13 levels of refinement was also performed offline. The
resulting animation can be viewed via [vimeo](https://vimeo.com/305325218). Note that the
refinement criterion was reduced ($\zeta_{u,v} = 0.002$), increasing
the maximum number of cells to 50000. Furthermore, the `CFL` number
was lowered (0.3) and the value of `NITERMIN` was set to `3`. It would
be nice if the methods were more robust.

##Reference

Kizner, Z., Makarov, V., Kamp, L., & Van Heijst, G. (2013). *Instabilities of the flow around a cylinder and emission of vortex dipoles*. Journal of Fluid Mechanics, 730, 419-441. doi:10.1017/jfm.2013.345
*/
