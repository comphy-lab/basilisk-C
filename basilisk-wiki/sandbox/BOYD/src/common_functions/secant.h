#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FUNC(x) ((*func)(x))
double secant(double (*func)(double), double x1, double x2, double E) {
  int n = 0;
  int max_iters = 1e7;
  double xm = 0.0;
  double x0 = 0.0;
  double c = 0.0;
  if (FUNC(x1) * FUNC(x2) < 0.0) {
    while (n < max_iters) {
      // calculate the intermediate value

      x0 = ((x1 * FUNC(x2) - x2 * FUNC(x1)) /
            (FUNC(x2) - FUNC(x1)));

      // check if x0 is root of
      // equation or not
      c = FUNC(x1) * FUNC(x0);

      // update the value of interval
      x1 = x2;
      x2 = x0;

      // update number of iteration
      n += 1;

      // if x0 is the root of equation
      // then break the loop
      if (c == 0) {
        break;
      }
      xm = ((x1 * FUNC(x2) - x2 * FUNC(x1)) /
           (FUNC(x2) - FUNC(x1)));

      if (fabs(xm - x0) < E) {
        break;
      }
    }
    // print("Root of the given equation =",
    //                       round(x0, 6));
    // print("No. of iterations = ", n);
    return x0;
  } else {
      fprintf(stderr,"Can not find a root in the given interval");
  }
  // error is returns nan
  return NAN;
}
#undef FUNC

