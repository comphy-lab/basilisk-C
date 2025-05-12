/**
# Face fractions from reconstructed interface
*/

#include "geometry.h"
#define VTOL 1.e-6

/**
  This is a modified version of the function
  [three--phase_test](basilisk.fr/sandbox/lopez/fracface.h) of JM LOPEZ which intend
  to generalise it to multiple vof tracer.
  The function below return the face fraction $sf.x$ at the selected cell face. As
  it is done for the *sweep_x()* function of [vof.h](/src/vof.h), we use the
  operator *foreach_dimension()* to automatize the derivation of the functions in
  the other dimensions. Once the dimension is selected ($x$, $y$ or $z$), the
  boolean variable *right* allows to select the face. If it is *TRUE* the face
  selected is that separating the cells [] and [1]. If *FALSE* the one returned is
  that between cells [] and [-1]. */

static double interface_fraction (coord m, double alpha, bool right)
{
#if dimension == 2
  alpha += (m.x + m.y)/2;
  coord n = m;
  double xo = (right ? 1. : 0.);
  foreach_dimension()
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
  foreach_dimension()
    if (n.x < 1e-4)
      return (n.x*(right ? 1 : -1) < 0. ? 1. : 0.);
  return clamp((alpha - n.x*xo)/n.y, 0., 1.);

#else /* dimension ==3 */

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  double xo = (right ? 1. : 0.);
  foreach_dimension()
    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
  foreach_dimension()
    if (fabs (n.x) < 1e-4)
      return (n.x*(right ? 1 : -1) < 0. ? 1. : 0.);
  return clamp((alpha - (n.x*xo+n.z*x0))/n.y, 0., 1.);
#endif  
}

/**
  A unique value of the face fraction is calculated from the reconstructed
  interfaces at the cells sharing the face by averaging as shown below,

  ![A unique value of the face fraction $s$ is calculated by averaging,
  $s = \sqrt{vleft \times vright}$](frac.svg)

*/

trace
struct FaceFraction {
  scalar c;          // compulsory
  face vector s;     // compulsory 
};

void face_fraction (struct FaceFraction a,b)
{
  scalar c = a.c;
  face vector s = b.s;
  boundary({c});

  /**
    We compute the normal vector in each cell to apply *boundary* to the vector
    field in order to get consistent values in the ghost cells.*/

  vector normal_vector[];
  foreach() {
    coord m = mycs (point, c);
    foreach_dimension() 
      normal_vector.x[] = m.x;
  }
  boundary((scalar*){normal_vector});

  foreach_face() {
    if (c[-1] < VTOL || c[] < VTOL) // some cell is empty
      s.x[] = 0.;
    else if (c[-1] > 1.- VTOL && c[] > 1.- VTOL) // both cells are full
      s.x[] = 1.;
    else {
      double vleft = 1., vright = 1.;
      if (c[] < 1. - VTOL) {
	coord m;
	foreach_dimension()
	  m.x = normal_vector.x[];
	double alpha = plane_alpha (c[], m);
	vleft = interface_fraction (m, alpha, false);
      }
      if (c[-1] < 1. - VTOL) {
	coord m;
	foreach_dimension()
	  m.x = normal_vector.x[-1];
	double alpha = plane_alpha (c[-1], m);
	vright = interface_fraction (m, alpha, true);
      }
      s.x[] = sqrt(vleft*vright);
    }
  }
}
