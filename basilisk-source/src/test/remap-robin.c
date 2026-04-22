/**
# Test of Dirichlet, Robin and Neumann conditions when remapping 

We remap a parabola defined on 3 layers onto 12 layers. Since
remapping is third-order accurate, this should be exact. */

#include "layered/remapc.h"

double f (double x) {
  return 3.*sq(x) + 2.*x + 1.;
}

double df (double x) {
  return 6.*x + 2.;
}

double integral (double a, double b)
{
  return (cube(b) + sq(b) + b - cube(a) - sq(a) - a)/(b - a);
}

void tests (int nposinit)
{
  int nposend = 12;
  double zinit[nposinit], zend[nposend];
  for (int i = 0; i < nposinit; i++)
    zinit[i] = i/(nposinit - 1.);
  for (int i = 0; i < nposend; i++)
    zend[i] = i/(nposend - 1.)*(1. + 0.1*noise()); // fixme: not correct noise
  zend[nposend - 1] = 1.;
  
  double init[nposinit - 1][1], end[nposend - 1][1];
  for (int i = 0; i < nposinit - 1; i++)
    init[i][0] = integral (zinit[i], zinit[i+1]);

  /**
  We first test Dirichlet conditions. */

  remap_c (nposinit, nposend, zinit, zend, 1, init, end,
           f(0), 0, 0, f(1), 0, 0, true);
  for (int i = 0; i < nposend - 1; i++) {
    //    fprintf (stderr,"%d %g %g\n", i, (zend[i] + zend[i+1])/2., end[i] - integral (zend[i], zend[i+1]));
    assert (fabs (end[i][0] - integral (zend[i], zend[i+1])) < 1e-13);
  }

  /**
  We then use Robin conditions. */

  const double lambda_b = 1., lambda_t = 2.;
  remap_c (nposinit, nposend, zinit, zend, 1, init, end,
           f(0) - lambda_b*df(0), lambda_b, 0,
           f(1) - lambda_t*df(1), lambda_t, 0, true);
  for (int i = 0; i < nposend - 1; i++) {
    //    fprintf (stderr,"%d %g %g\n", i, (zend[i] + zend[i+1])/2., end[i] - integral (zend[i], zend[i+1]));
    assert (fabs (end[i][0] - integral (zend[i], zend[i+1])) < 1e-13);
  }

  /**
  And finally Neumann conditions. */

  remap_c (nposinit, nposend, zinit, zend, 1, init, end,
           0, 0, df(0),
           0, 0, df(1), true);
  for (int i = 0; i < nposend - 1; i++) {
    //    fprintf (stderr,"%d %g %g %g\n", i, (zend[i] + zend[i+1])/2., end[i] - integral (zend[i], zend[i+1]), end[i]);
    assert (fabs (end[i][0] - integral (zend[i], zend[i+1])) < 1e-13);
  }
}

int main()
{
  tests (4);
  tests (3);
}
