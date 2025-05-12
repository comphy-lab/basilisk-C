/**
# Vortex instability

A mode-two instability for the isolated circular vortex:

For $R < \pi$,  
$$v_{\theta} = \mathrm{sin}^3(r)$$  
 
and. 
$$v_{\theta} =  0, \text{Else where}$$

![Vorticity](ins/omg.mp4)

![Refinement](ins/lev.mp4)

~~~gnuplot MG-solver properties
set yr [0:5.5]
set ylabel 'Iterations'
set xlabel 'Time step number'
set size square
set grid
plot 'out' u 2:3 t 'MG cycles for p', '' u 2:4 t 'Relax. sweeps for p', \
   '' u 2:($5 + 0.1)  t 'MG cycles for p2 + 0.1', '' u 2:($6+0.1) t 'Relax. sweeps for p2 + 0.1'
 ~~~
*/
#define RKORDER (4)
#define ADV_3 1
#define PROJECT_4 1
#include "nsf2.h"

int maxlevel = 10;

int main() {
  foreach_dimension()
    periodic (left);
  L0 = 10*pi;
  X0 = Y0 = -L0/2.;
  N = (1 << (maxlevel - 2));
  run();
}

#define PERTURB (0.0005)
#define _X (x/(1 + PERTURB))
#define _Y (y*(1 + PERTURB))
#define R (sqrt(sq(_X) + sq(_Y)))

double u_x (double x, double y) {
  return R < pi ? -y/R*cube(sin(R)) : 0;
}

double u_y (double x, double y) {
  return R < pi ?  x/R*cube(sin(R)) : 0;
}

event init (t = 0) {
  TOLERANCE = 1e-4;
  refine (R < 2*pi && level < maxlevel - 1);
  refine (R < pi && level < maxlevel );
  scalar omg[], psi[];
  foreach_face()
    u.x[] = Gauss6_x (x, y, Delta, u_x);
  boundary ({u.x, u.y});
}

event adapt (i++) {
  vector uc[];
  foreach_dimension()
    uc.x.prolongation = refine_4th;
  boundary_flux ({u});
  foreach()
    foreach_dimension()
      uc.x[] = (u.x[1] + u.x[])/2.;
  boundary ((scalar*){uc});
  adapt_wavelet ((scalar*){uc}, (double[]){0.001, 0.001}, maxlevel);
}

event logger (i++) {
  if (i == 0)
    printf ("t i p.i p.nr p2.i p2.nr cells depth speed divmax\n");
  printf ("%g %d %d %d %d %d %ld %d %g %g\n",
	  t, i, mgp.i, mgp.nrelax, mgp2.i, mgp2.nrelax,
	  grid->tn, grid->maxdepth, perf.speed, mgp2.resa);
}

event mov (t += 0.5) {
  scalar omg[];
  boundary ((scalar*){u});
  vorticityf (u, omg);
  output_ppm (omg, file = "omg.mp4", n = 512, linear = true,
	      map = cool_warm, min = -1, max = 1);
  foreach()
    omg[] = level;
  output_ppm (omg, file = "lev.mp4", n = 512,
	      min = 0, max = maxlevel);
}

event stop (t = 50) {
  scalar omg[];
  vorticityf (u, omg);
  for (double xp = -2; xp <= 2; xp += L0/(1 << maxlevel))
      fprintf (stderr, "%g %g\n", xp, interpolate (omg, xp, 0));
}
