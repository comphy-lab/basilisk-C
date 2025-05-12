/**
# Face fractions from reconstructed interface
*/

#include "geometry.h"
#define VTOL 1.e-6

/**
The function below return the face fraction $sf.x$ at the selected cell face. As
it is done for the *sweep_x()* function of [vof.h](/src/vof.h), we use the
operator *foreach_dimension()* to automatize the derivation of the functions in
the other dimensions. Once the dimension is selected ($x$, $y$ or $z$), the
boolean variable *right* allows to select the face. If it is *TRUE* the face
selected is that separating the cells [] and [1]. If *FALSE* the one returned is
that between cells [] and [-1]. */

foreach_dimension()
static double interface_fraction_x (coord m, double alpha, bool right)
{
#if dimension == 2
  alpha += (m.x + m.y)/2;
  coord n = m;
  double xo = (right ? 1. : 0.);
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }
  if (n.y < 1e-4)
    return (n.x*(right ? 1 : -1) < 0. ? 1. : 0.);
  return clamp((alpha - n.x*xo)/n.y, 0., 1.);
#else

#endif  
}

/**
A unique value of the face fraction is calculated from the reconstructed
interfaces at the cells sharing the face by averaging as shown below,

![A unique value of the face fraction $s$ is calculated by averaging,
$s = \sqrt{vleft \times vright}$](frac.svg)

*/

void face_fraction (scalar f, face vector s)
{
  boundary({f});
  
  /**
  We compute the normal vector in each cell to apply *boundary* to the vector
  field in order to get consistent values in the ghost cells.*/
  
  vector normal_vector[];
  foreach() {
    coord m = mycs (point, f);
    foreach_dimension() 
      normal_vector.x[] = m.x;
  }
  boundary((scalar*){normal_vector});

  foreach_face() {
    if (f[-1] < VTOL || f[] < VTOL) // some cell is empty
      s.x[] = 0.;
    else if (f[-1] > 1.- VTOL && f[] > 1.- VTOL) // both cells are full
      s.x[] = 1.;
    else {
      double vleft = 1., vright = 1.;
      if (f[] < 1. - VTOL) {
        coord m;
        foreach_dimension()
          m.x = normal_vector.x[];
        double alpha = plane_alpha (f[], m);
        vleft = interface_fraction_x (m, alpha, false);
      }
      if (f[-1] < 1. - VTOL) {
        coord m;
        foreach_dimension()
          m.x = normal_vector.x[-1];
        double alpha = plane_alpha (f[-1], m);
        vright = interface_fraction_x (m, alpha, true);
      }
      s.x[] = sqrt(vleft*vright);
    }
  }
}
