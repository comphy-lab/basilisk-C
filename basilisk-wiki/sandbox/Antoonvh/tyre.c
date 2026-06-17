/**
![Tyre inflation process](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcSzDacuwo35u4UZxSZzwwfPyZTTJNHLF21WHQ&s)
   
# Inflating a tyre

Inflating a tyre results in a convergent flow inside the tyre. We model it here:

![Using a potential flow solution](tyre/tyre.mp4)
 */
#include "run.h"
#include "embed.h"
#include "poisson.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"
vector u[], ub[];
scalar phi[];

double U = 1;
double W = 0.05, R = 1, tw = 0.1;

phi[embed] = (fabs(y) < W/2. && x > 0 && x < 1) ? neumann(-U) : neumann(0);
						     Particles pp;
							     
int main() {
  L0 = 2.5;
  X0 = Y0 = -L0/2;
  N = 256;
  run();
}

event init (t = 0) {
  vertex scalar V[];
  foreach_vertex()
    V[] = fabs(sqrt(sq(x) + sq(y)) - R);
  fractions (V, cs, fs, tw/2);
  foreach_face()
    fs.x[] = 1 - fs.x[];
  foreach()
    cs[] = 1 - cs[];
  scalar b[];
  foreach() {
    b[] = -1.043*(cs[] > 0)*W*U/(2*pi*R*tw);
  }
  TOLERANCE = HUGE;
  NITERMIN = 10;
  poisson (phi, b, alpha = fs);
  gradients ({phi}, {ub});
  
  
  pp = new_tracer_particles(0);
  int np = 0;
  while (np < 100) {
    double xp = L0 * noise();
    double yp = L0 * noise();
    
    foreach_point(xp, yp) {
      if (cs[] > 0.9) {
	particle p;
	p.x = xp, p.y = yp;
	add_particle (p, pp);
	np++;
      }
    }
  }
  DT = 0.001;
  dt = dtnext(DT);
}

event velocity (i++) {
  foreach()
    foreach_dimension() {
      u.x[] = ub.x[]*fabs(sin(t));
    }
}


event stability (i++) {
  dt = dtnext(DT);
}

event add_particles (t += 0.1) {
  int to_add = 2, np = 0;
  while (np < to_add) {
    double xp = (R - tw/2) + (tw/2.*noise());
    double yp = W*noise();
    foreach_point(xp, yp) {
      if (cs[] > 0.9) {
	particle p;
	p.x = xp, p.y = yp;
	set_a_particle_attributes(&p);
	add_particle (p, pp);
	np++;
      }
    }
  }
  
}

event movie (t += 0.1) {
  foreach_particle_in(pp) {
    foreach_point(p().x, p().y) {
      if (cs[] < 0.1)
	p().tag = 999;
    }
  }
  remove_particles (pp, p().tag == 999);
  
  scatter (pp);
  draw_vof ("cs", "fs", filled = -1, fc = {0.7,0.7, 0.7});
  save ("tyre.mp4");
}

event stop (t = 20);

