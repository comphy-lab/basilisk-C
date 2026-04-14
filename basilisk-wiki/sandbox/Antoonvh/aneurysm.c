/**
# Transfer and residence times in an aneurysm

We are interested in some Lagrangian statistics of a pulsatile flow trough a (2D) aneurysm.

![Lagrangian particles and Eulerian fracer-field age](aneurysm/age.mp4)

We compute the mean age of the particles in the simulation and use
that to estimate the residence time. We also compute the mean age of
the particles that leave the domain, and use that to estimate the
transfer time.

~~~gnuplot Residence time scale does not converge.
set xlabel 'time'
set ylabel 'Time scale estimate'
set key top left
plot 'out' u 1:3 t 'Transfer time', '' u 1:5 t 'Residence time', '' u 1:($5/$1) t 'Residence time over total time ratio'
~~~

The statistics also include initially seeded particles. These are not flushed.

 */
#include "embed.h"
#include "navier-stokes/centered.h"
// Particle members are a union of the fpm and tracer members
#ifndef ADD_PART_MEM
#define ADD_PART_MEM double s; double * sl; coord u; coord u2; coord locn; long unsigned int tag; 
#endif
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"
#include "tracer.h"
#include "fpm.h"
double D = 0.1;
scalar age[], * tracers = {age};
face vector muc[];
Particles Lag;
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
age[left] = dirichlet(t);
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

int maxlevel = 8;
int ref_level = 5;
int main() {
  
  L0 = 5;
  D = L0/(1<<ref_level);
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
    foreach_vertex()
      phi[] = exp(-sq(x)) + 1 - fabs(y);
    fractions(phi, cs, fs);
    multigrid_restriction ({cs, fs});
    Lag = new_tracer_particles(0);
    foreach_level(ref_level) {
      if (cs[] > 0.51) {
	particle p;
	p.x = x;
	p.y = y;
	set_a_particle_attributes (&p);
	p.tag = 0;
	add_particle (p, Lag);
      }
    }
}

event injections (t += 0.01) {
  double yp = -1 +  + (1 + noise())*D;
  while (yp < 1) {
    scalar ref[];
    foreach_cell_all()
      ref[] = 0;
    //multigrid_restriction({ref});
    particle p;
    p.x = X0 + 0.1;
    p.y = yp;
    yp += (1 + noise())*D;
    p.tag = (long int)100*t;
    set_a_particle_attributes (&p);
    int index[2];
    double dist[2];
    assign_particles (Lag, ref);
    
    ref_outdated = false;
    find_nearest_particles ({p.x, p.y, 0}, 1, Lag, index, dist = dist, reference = ref, level = 3);
    
    if (dist[0] > sq(D)) {
    
      add_particle (p, Lag);
      	
    }
  }
  particle_boundary(Lag);
  
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
  remove_particles (Lag, OUTFLOW_COND);
  
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

event movie (t += 0.02) {
  scalar Li[];
  foreach()
    Li[] = t - age[];
  scatter (Lag);
  
  draw_vof ("cs", "fs", filled = -1, fc = {0.5, 0.5, 0.5});
  squares ("Li", min = -0.1, max = t + 0.1);
  
  save ("age.mp4");
  

}

event adapt (i++) {
  adapt_wavelet ({cs, age, u}, (double[]){1e-3, t/100., 1e-3, 1e-3}, maxlevel, 6);
}

 event stop (t = 10.1) {
   foreach_cell_all()
     reference[] = 0;
   return 0;
 }
