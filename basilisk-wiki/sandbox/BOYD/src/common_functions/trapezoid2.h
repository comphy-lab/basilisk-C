#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FUNC(x, y) ((*func)(x, y))
double trapzd(double (*func)(double, double), double a, double b, int n,
              double var) {
  double x, tnm, sum, del;
  static double s;
  int it, j;

  if (n == 1) {
    double fa = FUNC(a, var);
    double fb = FUNC(b, var);
    s = 0.5 * (b - a) * (fa + fb);
    return s;
  } else {
    for (it = 1, j = 1; j < n - 1; j++) it <<= 1;
    tnm = it;
    del = (b - a) / tnm;
    x = a + 0.5 * del;
    for (sum = 0.0, j = 1; j <= it; j++, x += del) sum += FUNC(x, var);
    s = 0.5 * (s + (b - a) * sum / tnm);
    return s;
  }
}
#undef FUNC

#define EPS 1.0e-5
#define JMAX 20

double qtrap(double (*func)(double, double), double a, double b, double var) {
    double trapzd(double (*func)(double, double), double a, double b, int n, double var);
    int j;
    double s, olds;

    olds = -1.0e30;
    for (j = 1; j <= JMAX; j++) {
      s = trapzd(func, a, b, j, var);
      if (j > 5)
        if (fabs(s - olds) < EPS * fabs(olds) || (s == 0.0 && olds == 0.0))
          return s;
      olds = s;
    }
    fprintf(stderr, "Too many steps in routine qtrap");
    exit(1);
}
#undef EPS
#undef JMAX
