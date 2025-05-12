/**
# A sink in a stratified fluid

A circular sink wtih radius $R$ forces outflow at speed $U$. With
stratification strength defined by the BV frequency $N$, A parameter
can be defined,

$$\Pi = \frac{RN}{U}.$$

We are curious to see what fluid goes through the sink. 

## Results 

movies

![Simulation with $\Pi = 0.001$](sink/mov001.mp4)

![Simulation with $\Pi = 10$](sink/mov10.mp4)

We meassure the initial height of the particles that have sunk.

~~~gnuplot
set logscale x
set xr [0.0001 : 100]
set yr [0:1.5]
set grid
set xlabel 'Pi'
set ylabel 'Mean initial height [R]'
set size square
plot 'out'
~~~

## Conclusion:

Suction speeds $U > RN$ are required to suck in the warm air layers.

## Setup:
 */
#include "axi.h"
#include "navier-stokes/centered.h"
// particles are tagged with their initial height.
#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; double height;
#include "tracer-particles.h"
#include "tracer.h"
#include "view.h"
#include "scatter2.h"

double Nbv = 0.01, U = 1, R = 1;
scalar b[], *tracers = {b};

// Outflow with initial ramp
double tramp = 5;
u.n[left] = fabs(y) < R ? dirichlet (-(min(t/tramp, 1))) : dirichlet (0);
// Inflow conditions
u.n[right] = neumann (0);
p[right] = dirichlet (0);

int lev = 7;
face vector av[];

int main() {
  a = av;
  L0 = 15*R;
  X0 = -L0/2;
  DT = 0.1;
  for (Nbv = 0.001; Nbv <= 10; Nbv *= 10)
    run();
}

event init (t = 0) {
  init_tp_cells();
  foreach_particle() 
    p().height = p().x;
  refine (level < lev);
  foreach()
    b[] = sq(Nbv)*(x - X0);
  boundary ({b});
}

event acceleration (i++) {
  foreach_face(x)
    av.x[] = face_value(b,0);
}

event mov (t += 0.1) {
  if (Nbv == 0.001 || Nbv == 10) {
    view (psi = -pi/2.);
    squares ("u.x", min = -1.1, max = 1.1);
    scatter (0);
    draw_string ("Vertical velocity & tracers", size = 35, pos = 1);
    draw_string ("Buoyancy & Its isolines", size = 35, pos = 3);
    mirror ({0,1}) {
      squares ("b", 0, min = -0.1*sq(Nbv)*L0, max = 1.1*sq(Nbv)*L0, linear = true);
      for (double val = 0; val < sq(Nbv)*L0*1.1; val += sq(Nbv)*L0/20)
	isoline ("b", val);
    }
    if (Nbv == 0.001)
      save ("mov001.mp4");
    else
      save ("mov10.mp4");
  }
}

event stop (t = 4*tramp) {
  double mh = 0;
  int n = 0;
  foreach_particle() {
    if (x > R && p().height < -R) {
      mh += p().height;
      n++;
    }
  }
  printf ("%g %g\n", Nbv, mh/n - X0);
}
