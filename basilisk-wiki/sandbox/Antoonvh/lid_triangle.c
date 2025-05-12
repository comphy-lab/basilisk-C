/**

# Lid-driven triangle

   We may extend "Popinet's trick" for a poor symmetry condtion.

   ![Vorticity field: The flow in bound to the upper triangle. The diagonal wall is free-slip. Ignore the lower triangle mess](lid_triangle/o.mp4)

   ![Velocity magnitude](lid_triangle/U.mp4)
*/

#include "navier-stokes/centered.h"

// driven lid
u.t[top] = dirichlet(x*(1 - x));

const face vector muc[] = {0.001, 0.001};

int main() {
  N = 128;
  mu = muc;
  DT = 0.1;
  run();
}

event wall (i++) {
  foreach() {
    if (y < x) {
      // mirror the flow vector
      coord U = {(u.x[-1,1] + u.y[-1,1])/sqrt(2), (-u.x[-1,1] + u.y[-1,1])/sqrt(2)};
      u.x[] = (U.x + U.y)/sqrt(2);
      u.y[] = (U.x - U.y)/sqrt(2);
    }
  }
}

event mov (t += 0.1) {
  scalar omg[];
  vorticity (u, omg);
  output_ppm (omg, file = "o.mp4", n = 300);
  scalar U[];
  foreach() {
    U[] = 0;
    foreach_dimension()
      U[] = sqrt(sq(u.x[]) + sq(u.y[]));
  }
  output_ppm (U, file = "U.mp4", n = 300, max = 0.2, min = -0.2);
}

event stop (t = 50);