/**
# Implicit interpolation from face to cell-centered values

If we estimate face values (using `f.x[] = (s[] + s[-1])/2`) from
scalar fields and back (using `s2[] = (f.x[] + f.x[1])/2`), the original
field values are lost. Further, if we continue going back and forth,
the field may be smoothened considerably. Each pass,

~~~literatec
u[] += (u[-1] - 2u[] + u[1])/4
~~~

Which looks a lot like the stencil for diffusion. In order to omit
this "problem" we could go from face values to the centered values,
using the same stencil with the implicit problem.

~~~literatec
(s[] + s[-1])/2 = f.x[]
~~~

Such a procedure could upgrade the approximate projection of a
non-solenoidal centered vector field to a discretely exact one,
without increasing the stencil size.

On this page we test, and perhaps illustrate this idea.
*/

#include "solve.h"
#include "utils.h"

/**
We define the implicit face-to-centered solver, which is prone to
even-odd decoupling on periodic grids.
 */

foreach_dimension()
void face_to_centered_x (scalar f, scalar s) { //f represents face values of each cell.
  solve (s, (s[-1] + s[0])/2., f[]);
}

int main() {
  L0 = 10;
  X0 = Y0 = -L0/2.;
  init_grid (N);
  TOLERANCE = 1e-4; // Gives a more concincing result
  /**
The test considers the projection of a centered vector field;
   */
  scalar p[];
  vector u[], up[];
  face vector uf[];
  foreach() {
    p[] = u.y[] = 0;
    u.x[] = exp (-sq(x) - sq(y));
  }
  foreach_face()
    uf.x[] = face_value (u.x, 0);
  project (uf, p);
  // `up` holds to implicit interpolation
  foreach_dimension()
    face_to_centered_x (uf.x, up.x);
  // `u` is projected using the approximate projection method
  foreach() {
    foreach_dimension()
      u.x[] -= (p[1] - p[-1])/(2*Delta);
  }
  // We look at the divergence fields
  scalar div_u[], div_up[];
  foreach() {
    div_u[] = div_up[] = 0;
    foreach_dimension() {
      div_u[]  += (u.x[1]  - u.x[-1]) /(2*Delta);
      div_up[] += (up.x[1] - up.x[-1])/(2*Delta);
    }
  }
  double scale = 1e-2;
  output_ppm (div_u,  file = "div_u.png",  n = 300, min = -scale, max = scale);
  output_ppm (div_up, file = "div_up.png", n = 300, min = -scale, max = scale);
}
/**
## Results

![It is well known that the divergence field of the approximately-projected vector field is not zero ...](exact_centered_projection/div_u.png)

![... but the divergence field of the proposed projection method is](exact_centered_projection/div_up.png)

Well done proposed-projection method, very discretely exact!
*/
