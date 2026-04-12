/**
# Transfer and residence times in an aneurysm

We are interested in some Lagrangian statistics of a pulsatile flow trough a (2D) aneurysm.

![Injected particles](aneurysm/injected.mp4)

We compute the mean age of the particles in the simulation and use that to estimate the residence time. We also compute the mean age of the particles that leave the domain, and use that to estimate the transfer time. 

~~~gnuplot Residence time scale does not converge.
set xlabel 'time'
set ylabel 'Time scale estimate'
set key top left
plot 'out' u 1:3 t 'Transfer time', '' u 1:5 t 'Residence time', '' u 1:($5/$1) t 'Residence time over total time ratio'
~~~

The statistics also include initially seeded particles. These are not flushed.

![Initially placed particles](aneurysm/initial.mp4)
 */
#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"
face vector muc[];
Particles injected, initial;
double Re = 500;
double U0 = 10;
double U_in (double t, double y) {
  return -U0*fabs(sin(t*pi))*(y - 1)*(y + 1);
}

u.n[left]  = dirichlet(U_in(t, y)*fs.x[]);

p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

int maxlevel = 8;

int main() {
  L0 = 5;
  X0 = Y0  = -L0/2;
  mu = muc;
  N = 1 << maxlevel;
  run();
}

event properties (i++) {
  foreach_face() {
    muc.x[] = U0*fs.x[]/Re;
  }
}

event init (t = 0) {
    vertex scalar phi[];
    do {
      foreach_vertex()
	phi[] = exp(-sq(x)) + 1 - fabs(y);
      fractions(phi, cs, fs);
      multigrid_restriction ({cs, fs});
    } while (adapt_wavelet({cs}, (double[]){1e-9}, maxlevel, 5).nf);
    initial = new_tracer_particles(0);
    int n = 0;
    while (n < 500) {
      particle p;
      coord pos = {X0 + L0/2 + L0/2*noise(), Y0 + L0/2 + L0/2*noise(), 0};
      if (interpolate (cs, pos.x, pos.y) > 0.9) {
	p.x = pos.x;
	p.y = pos.y;
	set_a_particle_attributes (&p);
	add_particle (p, initial);
	p.tag = 0;
	n++;
      }
    }
    injected = new_tracer_particles (0);
}

event injections (t += 0.2) {
  int np = 20;
  for (int j = 0; j < np; j++) {
    particle p;
    p.x = X0 + 0.1;
    p.y = Y0 + L0/2 - 1. + 2.*(j + 0.5)/np;
    p.tag = (long int)100*t;
    set_a_particle_attributes (&p);
    
    add_particle (p, injected);
  }
  particle_boundary(injected);
}

#define OUTFLOW_COND (x > X0 + 99*L0/100.)

double ta = 0;
int na = 0;
event particle_outflow (i++) {
  foreach_particle() {
    if (OUTFLOW_COND) {
      na++;
      ta += t - (double)p().tag/100.;
    }
  }
  remove_particles (injected, OUTFLOW_COND);
  remove_particles (initial, OUTFLOW_COND);
}

event log_data (t += 0.1) {
  double tt = 0;
  int nt = 0;
  foreach_particle() {
    nt++;
    tt += t - (double)p().tag/100.;
  }
  
  if (na)
    printf ("%g %d %g %d %g\n", t, na, ta/na, nt, tt/nt);
}

event movie (t += 0.05) {
  scatter (injected);
  
  draw_vof ("cs", "fs", filled = -1, fc = {0.5, 0.5, 0.5});
  squares ("u.x", min = -0.1, max = U0);
  save ("injected.mp4");

  scatter (initial);
  
  draw_vof ("cs", "fs", filled = -1, fc = {0.5, 0.5, 0.5});
  squares ("u.x", min = -0.1, max = U0);
  save ("initial.mp4");

}

event adapt (i++) {
  adapt_wavelet ({cs, u}, (double[]){1e-3, 1e-4, 1e-4}, maxlevel, 6);
}

event stop (t = 10);