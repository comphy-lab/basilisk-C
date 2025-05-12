/**
# The component-wise (Hadamard) vector product: $\star$

Given a vector $\vec{a} = a_x\vec{e}_x + a_y\vec{e}_y + a_z\vec{e}_z$,
and an other vector $\vec{b}$, we can define a component-wise product
$\vec{c} = \vec{a} \star \vec{b}$, with,

$$\vec{c} = \vec{a} \star \vec{b} =  a_xb_x\vec{e}_x + a_yb_y\vec{e}_y + a_zb_z\vec{e}_z.$$

Indeed,

$$\vec{e}_i\star\vec{e}_j = \begin{cases} \vec{e}_i, & \text{for } i = j \\ 0, & \text{for } i \neq j \end{cases} $$

Euclidian vectors equipped with such a binary operation would form
an Abelian group and should be less controversial than the
scalar/dot/inner product or vector/cross product. In order to
investigate why it is not popular, we look at its geometric properties:

1. The "length" of the product vector is equal to the inner product: $\|\vec{a}\star\vec{b}\| = \vec{a} \cdot \vec{b}$.
1. The identity element breaks the symmetry between negative and positive directions...
1. ... Hence, a product triplet does not retain its orientation under rotation.
1. Its direction ... ?

Lets see if we can get inspiration from some examples:

![It is not very intuitive. C.f. the right-hand rule of the cross product](component_wise/comp.mp4)
 */
#include "grid/multigrid3D.h"
#include "view.h"
#include "some_primitives.h"

coord component_product (coord a, coord b) {
  return (coord) {a.x*b.x, a.y*b.y, a.z*b.z};
}

double len (coord c) {
  return sqrt(sq(c.x) + sq(c.y) + sq(c.z));
}

void draw_ab (coord a, coord b) {
  coord c = component_product (a, b);
  draw_arrow (a, fc = {0.9, 0.2, 0.2}, len = 1, rad = 0.1);
  draw_arrow (b, fc = {0.2, 0.9, 0.2}, len = 1, rad = 0.1);
  draw_arrow (c, len = len(c), fc = {0.9, 0.9, 0.2}, rad = 0.1);
  int size = 20;
  float lw = 3;
  draw_string ("Red", 1, size, {0.9, 0.2, 0.2}, lw);
  draw_string ("   *", 1, size, lw = lw);
  draw_string ("    Green", 1, size, {0.2, 0.9, 0.2}, lw);
  draw_string ("          =", 1, size, lw = lw);
  draw_string ("            Yellow", 1, size, {0.9, 0.9, 0.2}, lw);
}

int main () {
  init_grid (1);
  L0 = 3;
  X0 = Y0 = Z0 = -L0/2;
  double tend = 9;
  for (double t = 0; t < tend; t += 0.01) {
    coord a = {cos(t), sin(pi*t), sin(t/2)};
     normalize (&a);
     coord b = {cos(2*t), cos(3*t), sin(sqrt(2)*t)};
     normalize (&b);
     draw_ab(a, b);
     save ("comp.mp4");
     clear();
  }
}