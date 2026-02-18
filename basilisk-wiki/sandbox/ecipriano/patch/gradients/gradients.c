/**
# Interface gradients

We test the convergence of the interface gradients calculation. Unlike the many
tests already available for embed, in this test we focus on the accuracy of the
vof-biased scheme and the classic embed-scheme, proving $2^{nd}$ order convergence
of the vof-biased scheme which turns out to be the best option when we don't have
accurate knowledge of the face fractions. */

#include "fractions.h"
#include "gradients.h"
#include "run.h"

/**
We define the analytical function, its exact gradient at the interface, and a
function which reconstructs the interface centroid posistion, where we have to
compute the gradients. */

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
int vof;

/**
We run the system at different levels of refinement, with and without vof-biased
scheme. */

int main (void) {
  size (1 [*]);
  for (vof = 0; vof <= 1; vof++) {
    for (int maxlevel = 5; maxlevel <= 12; maxlevel++) {
      init_grid (1 << maxlevel);
      run();
    }
  }
}

#define circle(x,y,R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {

  /**
  The functions for which we want to compute the interface gradients. */

  foreach() {
    double r = sqrt (sq(x) + sq(y));
    s[] = fun (r);
  }

  /**
  We initialize the face fractions using `solid()` in order to focus on the
  convergence of the interface gradients, neglecting errors in the
  reconstructions of the face fractions. Using `face_fractions()` the vof-biased
  scheme remains $2^{nd}$ order, while the embed-scheme gets worse. */

  face vector sf[];
  solid (f, sf, circle (x, y, 0.5));

  /**
  Compute numerical gradients. */

  foreach()
    if (f[] > 0. && f[] < 1.) {
      double r = radius (point, f, sf, vof);
      gs[] = plic_gradient (point, s, f, sf, fun (r), vof, NULL);
    }

  /**
  Compute the error and print the average and maximum values. */

  scalar e[];
  foreach()
    if (f[] > 0. && f[] < 1.) {
      double r = radius (point, f, sf, vof);
      e[] = gs[] - exact (r);
    }
    else
      e[] = 0.;
  norm n = normf (e);
  fprintf (stderr, "vof %d %d %.3g %.3g\n", vof, N, n.avg, n.max);
}

/**
## Results

We obserse $2^{nd}$ order convergence in $L^1$ for the vof-biased scheme and
$3^{rd}$ order convergence for the default embed-scheme. For the latter case, it
is pivotal to reconstruct the interface normals using `facet_normal()` instead
of `interface_normal()` if we want to get a correct convergence trend (to
verify, use `interface_normal()` either in the function `radius()` of this test
and in `plic_gradient()`). However, the embed calculation of normals require
accurate face fractions. Using vof, we reconstruct those using
`face_fractions()` which demonstrated [$2^{nd}$ order
accuracy](/src/test/facefraction.c). Therefore, using the embed scheme for
interface gradients in vof simulations is not beneficial.

~~~gnuplot Convergence rate of vof-biased gradients
set xlabel "Resolution"
set ylabel "Error"

set xr[2**4:2**13]
set size square
set grid

set logscale x 2
set logscale y

f(x) = a*x**-b
fit f(x) "<grep 'vof 1' log" u 3:4 via a,b

f1(x) = a1*x**-b1
fit f1(x) "<grep 'vof 1' log" u 3:5 via a1,b1

ftitle(a,b) = sprintf("%.3f/x^{%4.2f}", exp(a), -b)

plot "<grep 'vof 1' log" u 3:4 title "L^1", \
     "<grep 'vof 1' log" u 3:5 title "L^{inf}", \
     f(x) w l title ftitle(a, b), \
     f1(x) w l title ftitle(a1, b1)
~~~

~~~gnuplot Convergence rate of embed-scheme gradients
fit f(x) "<grep 'vof 0' log" u ($3 > 2**7 ? $3 : $3/0):4 via a,b
fit f1(x) "<grep 'vof 0' log" u ($3 > 2**7 ? $3 : $3/0):5 via a1,b1

plot "<grep 'vof 0' log" u 3:4 title "L^1", \
     "<grep 'vof 0' log" u 3:5 title "L^{inf}", \
     f(x) w l title ftitle(a, b), \
     f1(x) w l title ftitle(a1, b1)
~~~
*/

