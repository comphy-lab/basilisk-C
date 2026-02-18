/**
# Interface gradients on trees

We test the compatibility between the interface gradients function and trees.
The setup is similar to the test case used for the convergence rate of interface
gradients ([gradients.c](gradients.c)). */

#include "fractions.h"
#include "gradients.h"
#include "run.h"

double fun (double r) {
  return r*r;
}

double exact (double r) {
  return 2.*r;
}

double radius (Point point, scalar f, face vector sf, bool vof) {
  coord m = vof ?
    interface_normal (point, f) :
    facet_normal (point, f, sf), p;
  double alpha = plane_alpha (f[], m);
  plane_area_center (m, alpha, &p);
  return sqrt (sq(x + p.x*Delta) + sq(y + p.y*Delta));
}

scalar f[], s[], gs[];
int vof = 1, maxlevel = 5;

int main (void) {
  size (1 [*]);
  init_grid (1 << maxlevel);
  run();
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {

  foreach() {
    double r = sqrt (sq(x) + sq(y));
    s[] = fun (r);
  }

  face vector sf[];
  solid (f, sf, circle (x, y, 0.5));
#if TREE
  f.prolongation = fraction_refine;
#endif

  while (adapt_wavelet ({f}, {1e-2}, maxlevel).nc)
    ;;
  unrefine (x >= 0.25);

  foreach()
    if (f[] > 0. && f[] < 1.) {
      double r = radius (point, f, sf, vof);
      gs[] = plic_gradient (point, s, f, sf, fun (r), vof, NULL);
      fprintf (stderr, "%g %g %.3g\n", x, y, gs[]);
    }
    else
      gs[] = nodata;
}

#include "view.h"

event picture (t = end) {
  clear();
  view (tx = -0.5, ty = -0.5, samples=12);
  draw_vof ("f", lw = 2.);
  squares ("gs");
  labels ("gs");
  cells();
  box();
  save ("picture.png");
}

/**
## Results

The values are uniform along the interface, confirming the correct interface
gradients calculation regardless of the level of refinement of the cells
involved in the normal probe.

![Interface gradients values](gradientstree/picture.png){width="70%"}
*/

