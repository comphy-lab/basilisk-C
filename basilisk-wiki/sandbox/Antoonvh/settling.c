/**
# Settling of a sand granule

The settling velocity of a sand granule in water is computed.

~~~gnuplot Looks OK (mind that $\tau \approx 0.027$)
set ylabel 'velocity [m/s]'
set xlabel 't [s]'
set grid
plot 'out' u 1:3:4 w l palette lw 4 t 'dt', -0.13625*(1-exp(-x/0.02778)) w l lw 2 t 'Theory'
~~~

Theory:

$$v_p(t) = V_t(1-e^{\frac{t}{\tau}}),$$
with,
$$V_t = -\frac{gd^2}{18\mu}\left(\rho_p - \rho_f \right),$$
and
$$\tau = -\frac{\rho_p4r_p^2}{18\mu}.$$

*/
#include "stokes-particles.h"

#define MUZ 1e-3
const face vector mu[] = {MUZ, MUZ, MUZ}; //mu Water
const scalar rho[] = 1000;                //rho Water
vector u[];

int main() {
  G.y = -9.81;   // g at the earth's surface
  for (DT = 0.05; DT > 0.001 ; DT /= 2.)
    run();
}

event init (t = 0) {
  new_inertial_particles(1);
  foreach_particle() {
    p().y = 0.5;               
    p().u2.x = 2000;          // approx. rho Silicon
    p().u2.y = 0.25e-3;       // Medium sand-grain size
  }
}

event set_dtmax (i++) 
  dt = dtnext(DT);

event trac (i += 2) {
  foreach_particle()
    printf ("%g %g %g %g\n", t, y, p().u.y, dt);
}

event stop (t = 0.4);
