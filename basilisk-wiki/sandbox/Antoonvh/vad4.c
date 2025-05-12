/**
# 4th-order finite difference on an adaptive tree-grid

Using the vertices instead of cell centers

![Result](vad4/mov.mp4)

~~~gnuplot Convergence
set xr [30:110]
set xlabel 'Cells at t = 20'
set ylabel 'L1-error'
set grid
set size square
set logscale xy
plot 'out', 0.8e5*x**(-4), 1.5e6*x**(-5)
~~~
 */
#include "grid/bitree.h"
#define RKORDER 4
#include "lsrk.h"
#include "my_vertex.h"

const vector u0[] = {1};
(const) vector u;
double D = 0.01;
#define t0 (1./(4*D))
#define SOL(x) (sqrt(t0/(t + t0))*exp(-sq(x)/(4*D*(t + t0))))

foreach_dimension()
void pade_x (scalar s, scalar dsdx) {
  double alphaP = 1./4., aP = 1.5;
  scalar rhs[];
  foreach() {
    dsdx[] = ((s[1] - s[-1])/(2.*Delta));
    rhs[] = aP* (s[1] - s[-1])/(2.*Delta);
  }
  boundary ({dsdx});
  for (int rel = 0; rel < depth(); rel++) {
    foreach() 
      dsdx[] = rhs[] - alphaP*(dsdx[-1] + dsdx[1]);
    boundary ({dsdx});
  }
}

void adv_diff (scalar * sl, scalar * dsl) {
  boundary (sl);
  scalar dsdx[];
  dsdx.prolongation = prolongation_vert4;
  dsdx.restriction = dsdx.coarsen = restriction_vert;
  pade_x (sl[0], dsdx);
  foreach() {
    scalar s, ds;
    for (s, ds in sl, dsl) {
      ds[] = 0;
      foreach_dimension()
	ds[] += (-u.x[]*dsdx[] +
		 D*(-s[2] + 16*s[1] - 30*s[] + 16*s[-1] - s[-2])/(12*sq(Delta)));
    }
  }
}

vertex scalar s[];
double es;

int main() {
  periodic (left);
  L0 = 20;
  X0 = -L0/2;
  N = 512;
  DT = 0.0005;  
  u = u0;
  for (es = 0.01; es > 3e-5; es /= 2)
    run();
}

event init (t = 0) {
  s.prolongation = prolongation_vert4;
  s.refine = refine_vert4;
  s.coarsen = s.restriction = restriction_vert;
  foreach_vert()
    s[] = SOL(x);
  boundary ({s});
}

event advance (i++, last) {
  dt = dtnext (DT);
  A_Time_Step ({s}, dt, adv_diff);
}

#include "adapt_field.h"
event adapt (i++) {
  scalar w[];
  restriction ({s});
  for (int l = 1; l <= depth(); l++) {
    foreach_level(l)
      if (child.x == 1)
	w[] = fabs(s[] - interp5_x (point, s));
    foreach_level(l)
      if (child.x == -1)
	w[] = w[1];
  }
  adapt_field (w, es, es/1.5, 19);
}
    

event stop (t = 20) {
  double e = 0;
  foreach_vert()
    e += Delta*fabs(s[] - SOL(x));
  printf ("%ld %g\n", grid->tn, e);
}

/** 
## The movie
*/
#if 1
FILE * gnuplotPipe;
event init_pipe (t = 0) {
  if (es == 0.01/8) {
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf (gnuplotPipe,
           "set term pngcairo\n"
           "set xr [%g: %g]\n"
	   "set yr [-0.1: 1.1]\n"
           "set grid\n"
	   "set size square\n"
           "set title 'Advection \\& Diffusion'\n"
	   "set xlabel 'x'\n"
	   "set ylabel 's'\n"
	   "set key off\n",
	   X0, X0 + L0);
  }
}

int frame_nr = 0;
event frames (t += 0.1) {
  if (es == 0.01/8) {
     fprintf (gnuplotPipe, "set output 'plot%d.png'\n", frame_nr);
     fprintf (gnuplotPipe, "plot '-' w l lw 1 lc rgb \"gray\","
	      "'-' w l lw 3, '-' ps 1 pt 5\n");
     
  for (double xp = X0;  xp <= (X0 + L0); xp += L0/(300.0001))
    fprintf (gnuplotPipe, "%g %g\n", xp, SOL(xp));
  fprintf (gnuplotPipe, "e\n");
  double firstp = -1;
  foreach_vert() {
    if (firstp == -1)
      firstp = s[];
    fprintf (gnuplotPipe, "%g %g\n",x, s[]);;
  }
  fprintf (gnuplotPipe, "%g %g\n", X0 + L0, firstp);
  fprintf (gnuplotPipe, "e\n");
  foreach_vert() 
    fprintf (gnuplotPipe, "%g %g\n",x, s[]);
  fprintf (gnuplotPipe, "e\n");
  frame_nr++;
  }
}

event stop (t = 20) {
   if (es == 0.01/8) {
     pclose (gnuplotPipe);
     system ("rm mov.mp4");
     system ("ffmpeg -loglevel quiet -r 25 -f image2 -i plot%d.png \\"
	     "-c:v libx264 -vf format=yuv420p -y mov.mp4");
  system ("rm plot*");
  return 1;
   }
}
#endif
