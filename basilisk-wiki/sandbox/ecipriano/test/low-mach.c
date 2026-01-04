/**
# VOF advection with non-solenoidal velocity

We test the convergence of the VOF method in combination with a prescribed non
divergence-free velocity field. */

#include "navier-stokes/low-mach.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#define DIVERGENCE 0

scalar divu[];

u.n[top] = neumann (0.);
u.t[top] = neumann (0.);
p[top] = dirichlet (0.);

u.n[right] = neumann (0.);
u.t[right] = neumann (0.);
p[right] = dirichlet (0.);

int maxlevel;
double R0 = 0.05;

int main (void) {
  rho1 = rho2 = 1;
  mu1 = mu2 = 1;

  f.sigma = 0.03;

  for (maxlevel = 4; maxlevel <= 9; maxlevel++) {
    init_grid (1 << maxlevel);
    run();
  }
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {
  fraction (f, circle (x, y, R0));
}

double exact (double t) {
  return R0 + t;
}

event output (t += 0.01) {
  double R = sqrt (4./pi*statsf (f).sum);
  double err = fabs (R - exact (t)) / exact (t);
  fprintf (stderr, "maxres %d %g %g %g %g\n",
      maxlevel, t, R, exact (t), err);
}

void exact_velocity (vector u, face vector uf, scalar divu) {
  foreach_face() {
    coord o = {x,y,z};
    double mag = sqrt (sq(x) + sq(y));
    uf.x[] = o.x / (mag + 1e-10);
  }

  foreach() {
    coord o = {x,y,z};
    double mag = sqrt (sq(x) + sq(y));
    foreach_dimension()
      u.x[] = o.x / (mag + 1e-10);
  }

  foreach() {
    divu[] = 0.;
    foreach_dimension()
      divu[] += (uf.x[1] - uf.x[]);
    divu[] /= -Delta;
  }
}

void reset_velocity (vector u, face vector uf) {
  foreach()
    foreach_dimension()
      u.x[] = 0.;

  foreach_face()
    uf.x[] = 0.;
}

#if DIVERGENCE
event tracer_diffusion (i++) {
  exact_velocity (u, uf, drhodt);
  reset_velocity (u, uf);
}

event end_timestep (t = end) {
  vector ue[];
  face vector ufe[];
  exact_velocity (ue, ufe, drhodt);

  scalar e[];
  foreach()
    e[] = norm (u) - norm (ue);
  norm n = normf (e);

  fprintf (stderr, "error %d %g %g %g\n", maxlevel, n.avg, n.rms, n.max);
}
#else
event end_timestep (i++) {
  exact_velocity (u, uf, drhodt);
}
#endif

#if TREE
event adapt (i++) {
  scalar ff[];
  foreach()
    ff[] = f[];
  adapt_wavelet ({f}, (double[]){1e-3}, maxlevel);
}
#endif

event stop (t += 0.1; t <= 0.7) {
  scalar magu[];
  foreach()
    magu[] = norm (u);

#if TREE
  restriction ((scalar *){u});
#endif
  clear();
  view (tx = -0.5, ty = -0.5);
  draw_vof ("f", lw = 1.5);
  squares ("magu", spread = -1);
  vectors ("u", scale = 0.001, level = 5);
  save ("movie.mp4");
}

/**
# Results

~~~gnuplot Evolution of the droplet radius in time
set xlabel "time"
set ylabel "Droplet Radius"
set size square
set grid
set key top left

plot "<grep 'maxres 5' log" u 3:4 w l dt 1 lc 2  t "LEVEL 5", \
     "<grep 'maxres 6' log" u 3:4 w l dt 1 lc 3  t "LEVEL 6", \
     "<grep 'maxres 7' log" u 3:4 w l dt 1 lc 4  t "LEVEL 7", \
     "<grep 'maxres 8' log" u 3:4 w l dt 1 lc 5  t "LEVEL 8", \
     "<grep 'maxres 8' log" u 3:4 w l dt 2 lc -1 t "Exact"
~~~

~~~gnuplot Convergence rate
stats "<grep 'maxres 4' log | tail -n 1" u 6 nooutput name "LEVEL4"
stats "<grep 'maxres 5' log | tail -n 1" u 6 nooutput name "LEVEL5"
stats "<grep 'maxres 6' log | tail -n 1" u 6 nooutput name "LEVEL6"
stats "<grep 'maxres 7' log | tail -n 1" u 6 nooutput name "LEVEL7"
stats "<grep 'maxres 8' log | tail -n 1" u 6 nooutput name "LEVEL8"
stats "<grep 'maxres 9' log | tail -n 1" u 6 nooutput name "LEVEL9"

set print "errors"

print sprintf ("%d %.12f", 2**4, LEVEL4_mean)
print sprintf ("%d %.12f", 2**5, LEVEL5_mean)
print sprintf ("%d %.12f", 2**6, LEVEL6_mean)
print sprintf ("%d %.12f", 2**7, LEVEL7_mean)
print sprintf ("%d %.12f", 2**8, LEVEL8_mean)
print sprintf ("%d %.12f", 2**9, LEVEL9_mean)

unset print

reset
set xlabel "N"
set ylabel "Relative Error"

set xr[2**3:2**10]
set size square
set grid

set logscale x 2
set logscale y

f(x) = a*x**-b
fit f(x) "errors" u ($1 > 2**5 ? $1 : $1/0):2 via a,b
ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)

plot "errors" w p pt 8 title "Results", \
     f(x) w l t ftitle(a,b)
~~~
*/
