/**
# Continuously injecting and removing new particles. 

The [confetti.c]() example shows how particles can be added by
creating new lists. Here we show how particles can be added to an
existing lists.

![Particles and their pid()](inject_particles/parts.mp4)

~~~gnuplot Snapshot of particles and their spend in the domain 
set grid 
set size ratio -1
set xlabel 'x'
set ylabel 'y'
plot 'out' u 1:2:3 with points pt 4 palette t 'Residence time'
~~~

~~~gnuplot Particles and allocated space on pid 0
reset
set xr [0:100]
set xlabel 't'
plot 'log' u 1:2 t 'Particles', '' u 1:3 t 'Allocated space' 
~~~ 
*/
#include "navier-stokes/centered.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

u.n[left] = dirichlet (fabs(y) < 0.5);
u.t[left] = (fabs(y) < 0.5 ? dirichlet (0.3*sin(t)) : neumann (0));
p[left] = neumann (0);

u.n[right] = neumann (0);
p[right] = dirichlet (0);

Particles trackers;

int main() {
  L0 = 10;
  Y0 = X0 = -L0/2.;
  run();
}

event init (t = 0) {
  trackers = new_tracer_particles (0);
}

event adder (t = 1; t += 0.1) {
  for (double yp = -0.5 ; yp <= 0.5; yp += 0.1) {
    particle p = {.x = X0, .y = yp, .z = 0, .tag = t*10};
    set_a_particle_attributes (&p);
    add_particle (p, trackers);
  }
}

event remover (i++) {
  remove_particles (trackers, x > 3.2);
}

event damp (i++) {
  foreach()
    if (x > 3)
      u.y[] -= u.y[] * dt*(1. - exp(-(x - 3)))/2.;
  boundary ({u.y});
}

event adapt (i++) {
  adapt_wavelet ((scalar*){u}, (double[]){0.025, 0.025}, 7);
}

event trace_p (t += 0.1) {
  if (pid() == 0)
    fprintf (stderr, "%g %ld %ld\n", t, pn[trackers], pna[trackers]);
}

event mov (t += 0.5) {
  scatter (trackers, pc = {sin(pid()), cos(pid()), sin(2.5*pid())});
  box();
  save ("parts.mp4");
}

event stop (t = 100) {
  for (int i = 0; i < npe(); i++) {
    foreach_particle_in(trackers) {
      if (i == pid())
	printf ("%g %g %g\n", x, y, (t - p().tag/10));
    }
#if _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
}

