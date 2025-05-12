/**
# A test for for the Stokes-particles

We place nine Stokes Particles in a Taylor-green vortex, with gravity.
   
 ~~~gnuplot Particle trajectories
 set xr[0:3.1415]
 set yr[0:3.1415]
 set size ratio -1
 set key off
 plot 'out'
  ~~~

We may compare against the result of Dagan (2025)

![Looks very similar to Dagan's analytical trajectories. Mind that the x range extends past $x = \pi$,
 when comparing the trajectory with the sharp
 corner](https://www.antoonvanhooft.nl/media/dagan.png)
*/

vector u[];
face vector mu;
scalar rho[];// = {0};

#include "stokes-particles.h"

Particles nine;

int main() {
  N = 512;
  L0 = pi;
  run();
}

event init (t = 0) {
  G.y = -400;
  nine = new_inertial_particles(9);
  int j = 0;
  foreach_particle_in(nine) {
    p().x = (j++ + 1)*pi/18;
    p().y = pi/2.;
    p().u2.z = 0.001; //tau
    p().u2.x = 1;     //rho - rhof
  }
  foreach() {
    u.x[] = sin(x)*cos(y);
    u.y[] = -sin(y)*cos(x);
  }
  DT = 0.01;
}

event set_dtmax (i++) {
  dt = dtnext(DT);
}

event mov(i++) {
  remove_particles (nine, y > pi/1.99);
  foreach_particle_in(nine)
    printf ("%g %g\n", x, y);
}

event stop (t = 2*pi);

/**
## Reference

~~~bib
@article{dagan2025,
  title={Analytical solutions for particle dispersion in Taylor--Green vortex flows},
  author={Dagan, Yuval},
  journal={Theoretical and Computational Fluid Dynamics},
  volume={39},
  number={1},
  pages={14},
  year={2025},
  publisher={Springer}
}
~~~
 */
