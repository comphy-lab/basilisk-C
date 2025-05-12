/**
# Particle advection test

The order of accuracy with a time-dependent flow field is
diagnosed for two cases.

* The timestep is highly variable, to represent a worst-case
scenario ($\alpha \in \langle0.01, 0.99\rangle$).  
* The timestep is fixed to represent the best-case scenario

![An approximated particle path](tp18/mov.mp4)

~~~gnuplot Error convergence is at the intended order
set xr [100:3000]
  set logscale y
  set logscale x 2 
  set xlabel '#Steps'
  set ylabel 'Error' 
  set grid 
  set size square    
set key right outside
  plot 'rk2log' t 'Rk2 fixed step', 'rk2out' t 'RK2 variable step',\
        9e3*x**(-2) lw 2 t '2nd order', 'log' t 'RK3 fixed step',\
	'out' t 'RK3 variable step', 2e5*x**(-3) lw 2 t '3rd order'
~~~ 
*/

vector u[];
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

int its;
double Ub = 1, amp = 1.5, periods = 2.412;
bool fixed = false;

int roundUp (int num, int mult) {
  int rem = num % mult;
  if (rem == 0)
    return num;
  return num + mult - rem;
}

FILE * fp;
int main () {
  L0 = 3;
  X0 = Y0 = -L0/2;
  P_RK2 = true;
  fp = fopen ("rk2out", "w");
  for (DT = 0.1; DT > 0.005; DT /= 2) {
    its = roundUp (periods*2*pi/DT, 2);
    run();
  }
  fixed = true;
  fclose (fp);
  fp = fopen ("rk2log", "w");
  for (DT = 0.1; DT > 0.005; DT /= 2) {
    its = roundUp (periods*2*pi/DT, 2);
    run();
  }
  fclose (fp);
  fixed = false;
  P_RK2 = false;
  fp = stdout;
  for (DT = 0.1; DT > 0.005; DT /= 2) {
    its = roundUp (periods*2*pi/DT, 2);
    run();
  }
  fixed = true;
  fp = stderr;
  for (DT = 0.1; DT > 0.005; DT /= 2) {
    its = roundUp (periods*2*pi/DT, 2);
    run();
  }
}

event init (t = 0) {
  new_tracer_particles (1);
  foreach_particle() {
    p().x  = 1.;
    p().y  = 0.;
  }
  set_particle_attributes (0);
  dt = dtnext (DT);
}

event set_dtmax (i++) {
  if (fixed)
    dt = dtnext (DT);
  else
    dt = dtnext ((0.99*noise() + 1)*DT);
}

event velocity (i++) {
  foreach() {
    coord v = {y, -x};
    foreach_dimension()
      u.x[] = (amp*sin(t) + Ub)*v.x;
  }
  boundary (all);
}

event movie (i += 4) {
if (DT < 0.06 && DT > 0.03 && fixed && !P_RK2) {
    box();
    scatter (0);
    save ("mov.mp4");
  }
}

/**
## Diagnose the approximate particle locations

`loc_n` only represents the approximate particle positions
corresponding to $t=\mathtt{t}$ at even iterations. `loc` stores
intermediate RK-related guesses.
 */
event stop (i = its, last) {
  double theta = amp + t*Ub - amp*cos(t);
  coord sol = (coord){cos(theta) , -sin(theta)};
  foreach_particle()
    fprintf (fp, "%d %g\n",i,
	     sqrt(sq(p().locn.x - sol.x) + sq(p().locn.y - sol.y)));
  return 1;
}
