/**
# Does potential flow yield optimal transport?

We consider the mapping of particle positions from one place to
another. We want to assign $n$ initial positions to $n$ final
positions, such that the "transport cost" ($D$) is minimized. Here we
take the squared distance between final and initial point as our cost
function, and sum it for all transported particles to find $D$.

~~~gnuplot Transport problem for approximately 140 particles.
set term svg size 900,200
set xr [-5: 5]
set yr [-1:1]
set size ratio -1
set key top outside
plot 'log' u 1:2 t 'Final positions', '' u 3:4 t 'initial positions'
~~~

Gues what! *A* solution to this problem can be found via potential flow:

![The proposed solution of the transport problem (indeed, the problem is derived from this solution)](optimal_transport/optimal.mp4)

Is this mapping, achieved by lagrangian tracing of a
potential flow field, the optimal map? A nesecarry condition is that if two
particles switch their mapped positions, $D$ will always increase.

~~~gnuplot Analysis of approx $10^4$ pairs, reveals some negative increases...
reset
set ylabel 'Increase in D'
set xlabel 'Particle-pair nr.'
plot 'out' u 3 t 'Increase'
~~~

It seems that the potential flow transport is not optimal. We plot the pair that would decrease $D$ the most if we switched their mapping.

~~~gnuplot
reset
set term svg size 900,200
set xr [-5: 5]
set yr [-1:1]
set size ratio -1
set key top outside
plot 'worst_score' u 1:2 t 'Final positions', '' u 3:4 t 'initial positions'
~~~

A clear case of particles "over taking" each other.

## Methods

We will extend the definition of a particle with its initial location.
*/
#include "run.h"
#define ADD_PART_MEM  coord u; coord u2; coord locn; long unsigned int tag; coord init_pos; 
#include "tracer-particles.h"
#include "embed.h"
#include "poisson.h"
#include "view.h"
#include "scatter2.h"

Particles ot;

scalar n[], a[], zeros[];
vector u[];

a[left] = dirichlet (cs[] > 0);
a[right] = dirichlet (0);

int main() {
  L0 = 10;
  X0 = Y0 = -L0/2.;
  DT = 0.01;
  N = 512;
  dt = dtnext(DT);
  run();
}

event dt_next (i++) {
  dt = dtnext (DT);
}

event init (t = 0) {
  vertex scalar phi[];
  scalar cs1[];
  face vector fs1[];
  do {
    // Geometry
    // Bottom boundary
    foreach_vertex()
      phi[] = 1 + y + 0.3*sin(2*x);
    fractions (phi, cs, fs);
    // top boundary
    foreach_vertex()
      phi[] = 1 - y + 0.5*cos(3*x);
    fractions (phi, cs1, fs1);
    foreach()
      cs[] = min (cs1[], cs[]);
    foreach_face()
      fs.x[] = min (fs.x[], fs1.x[]);
    // circle
    foreach_vertex()
      phi[] = sqrt(sq(x) + sq(y)) - 0.5;
    fractions (phi, cs1, fs1);
    foreach()
      cs[] = min (cs1[], cs[]);
    foreach_face()
      fs.x[] = min (fs.x[], fs1.x[]);
    // Potential 
    poisson (a, zeros, fs);
  } while (adapt_wavelet ({a}, (double[]){1e-6}, 10).nf > grid->tn/100);
  // Potential flow
  foreach_face() 
    fs1.x[] = fs.x[] * (a[] - a[-1])/Delta;
  foreach() {
    foreach_dimension()
      u.x[] = -(fs1.x[1] + fs1.x[0])/2;
  }
  boundary ({u.x, u.y});
  // Initialize particles
  ot = init_tp_circle (150, -4, 1e-3, l = 0.5);
  foreach_particle() 
    p().init_pos = (coord){x, y};
  tag_particles (ot);

  foreach()
    a[] = level;
  output_ppm (a, file = "grid.png", min = 5, max = 11);
}
/**
## Mesh

We inspect the mesh for its quality

![More red colors correspond to higher levels of refinement](optimal_transport/grid.png)
 */

event movie (t += 1) {
  scalar U[];
  foreach()
    U[] = sqrt(sq(u.x[]) + sq(u.y[]));
  squares ("U");
  translate (z = 1e-2) {
    draw_vof ("cs", "fs", filled = -1, fc = {0.6, 0.6, 0.6});
    draw_vof ("cs", "fs", lw = 2);
  }
  scatter (ot, s = 2);
  save ("optimal.mp4");
}
/**
At the boundaries, the potential flow field, and its interpolated
values are poorly estimated. As such, particles that get too close to
the boundary are removed.
 */
event to_remove (i += 2) {
  foreach_particle_in(ot) {
    if ((sq(p().x) + sq(p().y)) < sq(0.51)) {
      p().x = 999;
    }
  }
  remove_particles(ot, p().x == 999);
}

/**
## Analysis of particle-pair exchange

We may nest `foreach_particle()` iterators!
 */
event stop (t = 80) {
  double mind = 0;
  int taga = 0, tagb = 0;
  foreach_particle() {
    double xp1 = p().locn.x, yp1 = p().locn.y;
    double xpi1 = p().init_pos.x, ypi1 = p().init_pos.y;
    double d_1 = sq(xp1 - xpi1) + sq(yp1 - ypi1);
    long int tag1 = p().tag;
    fprintf (stderr, "%g %g %g %g\n", xp1, yp1, xpi1, ypi1);
    foreach_particle() {
      if (tag1 < p().tag) {
	double xp2 = p().locn.x, yp2 = p().locn.y;
	double xpi2 = p().init_pos.x, ypi2 = p().init_pos.y;
	double d1 = d_1 + sq(xp2 - xpi2) + sq(yp2 - ypi2);
	double d2 = sq(xp1 - xpi2) + sq(yp1 - ypi2) + sq(xp2 - xpi1) + sq(yp2 - ypi1);
	printf ("%g %g %g\n", d1, d2, d2 - d1);
	if (d2 - d1 < mind) {
	  mind = d2 - d1;
	  taga = tag1;
	  tagb = p().tag;
	}
      }
    }
  }
  FILE * fp = fopen ("worst_score", "w");
  foreach_particle() {
    if (p().tag == taga || p().tag == tagb)
      fprintf (fp, "%g %g %g %g\n", x, y, p().init_pos.x, p().init_pos.y);
  }
  fclose (fp);
}
