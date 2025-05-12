/**
# Advection of the scalar front in a Poiseuille flow

We study the advection of a mixing front in a viscous tube flow.

![Advection (and diffusion) of the scalar front. The lines indicate
 different sections with constant scalar
 concentration](pois-front/mov.mp4)

Since the the mixing front gets longer over time, the iso-lines are
likely to move faster than the flow (they are not material lines). We track their
position and determine their speed.

~~~gnuplot
ftitle(a,b) = sprintf("z = %.02ft + %.02f", (a), -b)


f(x) = a*x + b
fit f(x) 'out' u 1:2 via a,b

g(x) = c*x + d
fit g(x) 'out' u 1:4 via c,d

h(x) = k*x + l
fit h(x) 'out' u 1:6 via k,l

set size square
set xlabel 'time'
set ylabel 'z'
set grid
set key top left
plot 'out' u 1:2 t 'c = 10%', f(x) t ftitle(a,b) lw 2,	\
'' u 1:4 t 'c = 50%', g(x) t ftitle(c,d) lw 2,\
'' u 1:6 t 'c = 90%', h(x) t ftitle(k,l) lw 2
~~~

Varying speeds are indeed diagnosed.
 */

#include "embed.h"
#include "axi.h"
#include "advection.h"
#include "diffusion.h" 
#include "run.h"
#include "view.h"

int maxlevel = 8;
scalar c[], * tracers = {c};

int main() {
  L0 = 20;
  periodic (left);
  DT = 0.1 ;
  run();
}

event init (t = 0) {
  double Um = 1;
  refine (y < 2 && level < maxlevel - 1);
  refine (y < 1.1 && y > 0.9 && level < maxlevel);
  solid (cs, fs, 1 - y);
  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);
  foreach_face(x)
    uf.x[] = y < 1 ? Um*(1 - sq(y)) : 0;
  foreach() 
    c[] = x < 1 ? cs[] : 0;
}

event init_c (i++) {
  foreach() 
    c[] = x < 1? cs[] : c[];
}

event traccer_diffusion (i++) {
  face vector kappa[];
  foreach_face()
    kappa.x[] = fs.x[]*0.1;
  diffusion (c, dt, kappa);
}

event mov (t += 0.5) {
  output_ppm (c, file = "c.mp4");
  view (fov = 7, tx = -0.5, width = 600, height = 200);
  draw_vof ("cs", "fs", filled = -1, fc = {0.4, .4, .4});

  squares ("c", min = 0, max = 1);
  isoline ("c", n = 5, min = 0.1, max = 0.9);   
  mirror ({0,1}) {
    draw_vof ("cs", "fs", filled = -1, fc = {0.4, .4, .4});   
    squares ("c", min = 0 , max = 1);
    isoline ("c", n = 5, min = 0.1, max = 0.9);
  }
  save ("mov.mp4");
}
/**
We log the position of the fonts.
*/

event logger (t = 5; t += 0.1) {
  int n = 5;
  double pos[5] = {0};
  double vol[5] = {0};
  foreach(reduction(+:pos[:n]) reduction(+:vol[:n])) {
    if (y < 0.9 && x < 9*L0/10) {
      for (int i = 0; i < n; i++) {
	double val = 0.1 + 0.2*i;
	if (fabs(c[] - val) < 0.05 ) {
	  pos[i] += x*cm[]*sq(Delta);
	  vol[i] += cm[]*sq(Delta);
	}
      }
    }
  }
  printf ("%g", t);
  for (int i = 0; i < n; i++)
    printf (" %g", vol[i] > 0 ? pos[i]/vol[i] : vol[i]);
  printf ("\n");
}

event stop (t = 15);
