/**   
# Dust separator

We consider a stream of dusty air that enters a large box via an
inlet. It is sucked in by an outlet flow. Large particles may settle.

![The dust particles in a box. Red colors indicate high-speed regions](dust_separator/dust.mp4)

We compare the number of collected particles ($C$, bottom) against those that are sucked downstream ($D$, top).

~~~gnuplot About 70% of the particles are filtered.
set ylabel 'C/D'
set xlabel 'time'
set grid 
plot [4:20] 'out' u 1:($6/($8+1)) w l lw 2 t 'data'
~~~

We also compare the average radius of the particles that are collected
to those that are sucked downstream.

~~~gnuplot Large particles are collected
set ylabel 'Average radius'
set xlabel 'time'
set grid
plot 'out' u 1:7 w l lw 2 t 'Collected', '' u 1:9 w l lw 2 t 'Downstream'
~~~
 */
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "utils.h"
#include "stokes-particles.h"
#include "view.h"
#include "scatter2.h"
scalar fi[], fo[];

#define RAD (fabs(x))

double L = 0.5, Ri = 0.02;

double U0 = 0.5;

u.n[top] = dirichlet (U0*min(t , 1) * fo[]);
uf.n[top] = (U0*min(t , 1)*fo[]);
u.t[top] = dirichlet(0);
p[top]   = fo[]*neumann(0.) + (1 - fo[])*neumann (neumann_pressure(ghost));
pf[top]  = neumann(0.);

u.n[left] =  fi[] ? neumann (0) : dirichlet(0); 
uf.n[left] = fi[] ? neumann (0) : 0;

p[left]    = fi[]* dirichlet(0.)  + (1 - fi[]) *neumann (-neumann_pressure(0));
pf[left]   = fi[] * dirichlet(0.) + (1 - fi[])*neumann (0);
u.t[left] = dirichlet (0);
Particles dust;

int main() {
  L0 = L;
  X0 = Y0 = -L0/2.;
  N = 128;
  DT = 0.01;
  const face vector muc[] = {2e-3, 2e-3};
  mu = muc;
  G.y = -9.81;
  run();
}

event init (t = 0) {
  fraction (fi, (Ri) - sqrt(sq(y) + sq(z)));
  fraction (fo, (Ri) - sqrt(sq(x) + sq(z)));
  restriction ({fi, fo});
  output_ppm (fi, file = "fi.png");
  output_ppm (fo, file = "fo.png");
  dust = new_inertial_particles(0);
}

event add_dust_particles (t += 0.01) {
  for (double yp = -Ri/2. ; yp < Ri/2.; yp += Ri/10.) {
    particle p;
    p.x = X0 + 1e-3;
    p.y = yp;
    p.u2.y = (.6 + 0.5*noise())*5e-4;
    p.u2.x = 500.;
    p.u2.z = 0.;
    foreach_point (p.x, p.y, serial) 
      p.u.x = interpolate (u.x, p.y, p.y);
    p.u.y = 0.;
    add_particle (p, dust);
  }
}

#define DOWNSTREAM (y > Y0 + L0 - Ri/10 && fabs(x) < Ri)
#define COLLECTED (y < Y0 + Ri/10)

double rad_col = 0, rad_down = 0;
int col = 0, downstr = 0;

event remove_particles_event (i++) {
  foreach_particle_in(dust, reduction(+:rad_down) reduction(+:downstr) reduction(+:rad_col) reduction(+:col)) {
    if (DOWNSTREAM) {
      rad_down += p().u2.y;
      downstr++;
    }
    else if (COLLECTED) {
      rad_col += p().u2.y;
      col++;
    }
  }
  remove_particles (dust, (DOWNSTREAM || COLLECTED));
}


event mov (t += 0.02) {
  printf ("%g  %d %d %d %d %d %g %d %g\n", t, mgp.i, mgp.nrelax, mgpf.i, mgpf.nrelax,
	  col, col ? rad_col/col : 0, downstr, downstr? rad_down/downstr: 0);
  scalar U[];
  foreach()
    U[] = sqrt(sq(u.x[]) + sq(u.y[]));
  squares ("U", map = blue_white_red, min = -0.5, max = 0.5);
  scatter(dust, s = 2, pc = {0.3, 0.3, 0.3});
  save ("dust.mp4");
}

event stop (t = 20);