/**
# Flow around a cylinder in 2D

Following Stephane's setups for the [Von Karman-street
example](/src/examples/karman.c) and the [starting flow around a
cylinder test](/src/tests/starting.c), we also setup a similar
scenerio. The question we answer here is: What is a sufficient resolution
to determine the pressure drag, that is dominant in a configuration
where $Re = 500$.

![A simulation result](cyl2d/mov.mp4)

As such, we study the convergence of the drag coefficient ($C_d$) as a
function fo the used effort ($\approx cells \times
\mathrm{iteration}$). We vary the maximul level of refinement and the
refinement criterion ($c$, such that $\zeta_u = U/c$). The result is plotted below.

![Cells and $C_d$ for various runs. The legend lists the value of $c$ and the maximum level `ML` like: "c-ML" ](http://www.basilisk.fr/sandbox/Antoonvh/convergence_cyl.png)

We conclude that the simulation with $c = 15$ and `maxlevel` = 9 appears to capture the dynamics sufficiently, width a modest effort.
 */
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"

#define CYLINDER (sqrt(sq(x) + sq(y - 0.01*R)) - R)

double R = 1, U = 1, Re = 500, c = 15., ue;
int maxlevel = 9;

u.n[left]  = dirichlet (U);
uf.n[left] = U;
u.t[left]  = dirichlet (0);
p[left]    = dirichlet (0.);
pf[left]   = dirichlet (0.);

u.n[right] = neumann (0.);
p[right]   = neumann (0.);
pf[right]  = neumann (0.);

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

FILE * fp;
face vector muc[];
int main() {
  periodic (top);
  L0 = 30;
  X0 = -L0/3.5;
  Y0 = Z0 = -L0/2.;
  mu = muc;
  NITERMIN = 2;
  TOLERANCE = 1e-4;
  char logname[99];
  /**
The convergence test is conducted by uncommenting these lines:
   */
  //for (c = 5; c <= 20; c += 5)
  //  for (maxlevel = 7; maxlevel <= 10; maxlevel++) {
      ue = U/c;
      sprintf (logname, "log%d-%g", maxlevel, c);
      fp = fopen (logname, "w");
      run();
  //}
}

event init (t = 0) {
  refine (CYLINDER < 0.2*R && CYLINDER > -0.2*R && level < maxlevel);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = CYLINDER;
  boundary ({phi});
  fractions (phi, cs, fs);
  foreach()
    u.x[] = cs[] > 0;
}
/**
In order to omit issues with inconsitent boundary conditions, we use a
damping layer near the in and out flow boundaries.
 */
event damp (i++) {
  coord Uinf = {U, 0, 0};
  foreach() {
    if (fabs(x - (X0 + L0/2.)) > 4*L0/10.) 
      foreach_dimension()
	u.x[] += dt*(Uinf.x - u.x[])/2.;
  }
  boundary ((scalar*){u});
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*R*U/Re; 
  boundary ((scalar*){muc});
}

event adapt (i++) {
  adapt_wavelet ((scalar*){u}, (double[]){ue, ue, ue}, maxlevel, 5);
  unrefine (level > 5 && (x - X0) < L0/10. && (X0 + L0 - x) < L0/10.);
}

event logger (i += 5) {
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  fprintf (fp, "%d %g %g %g %g %g %ld\n",
	   i, t, Fp.x, Fp.y, Fmu.x, Fmu.y, grid->n);
}

#include "view.h"
event movies (t += 0.4) {
  view (fov = 6, width = 1200, height = 400, tx = -0.25);
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  translate (z = 0.05) {
    draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
    draw_vof ("cs", "fs", lw = 2);
  }
  squares ("omega", min = -2, max = 2, linear = true, map = cool_warm);
  cells();
  char str[99];
  sprintf (str, "Re = %g, C = %g, ML= %d", Re, c, maxlevel);
  draw_string (str, 1, lw = 3, lc = {1, 0, 1});
  save ("mov.mp4");
}

event stop (t = 150) 
  fclose (fp);
