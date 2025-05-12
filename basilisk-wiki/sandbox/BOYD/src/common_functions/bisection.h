#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FUNC(x) ((*func)(x))
double bisection(double (*func)(double), double a, double b, double E) {
  if (func(a) * func(b) >= 0) {
    fprintf(stderr, "You have not assumed right a and b\n");
    return NAN;
  }

  int n = 0;
  int max_iters = 1e7;
  double c = a;
  while ((b - a) >= E) {
    n++;
    // Find middle point
    c = (a + b) / 2;

    // Check if middle point is root
    if (func(c) == 0.0) break;

    // Decide the side to repeat the steps
    else if (func(c) * func(a) < 0)
      b = c;
    else
      a = c;
    if (n>max_iters){
      fprintf(stderr,"Can not find a root in the given interval");
    }
  }
  fprintf(stderr, "The value of root is : %g\n", c);
  return c;
}
#undef FUNC

