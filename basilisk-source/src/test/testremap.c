/**
# Comparison between Fortran and C remapping */

#include "ppr/ppr.h"
#include "layered/remapc.h"

/**
## Test functions

We validate our C implementation by using 3 test functions: a
polynomial of order 2, a piecewise constant function and a sum of
Gaussians. */

double polynomial (double z) {
  return sq(z);
}

double sharp (double z) {
  return (z >= -7. && z < -3. ? 4./10. :
	  z >= -3. && z < 1. ? 12./10 :
	  z >= 1. && z < 4. ? 8./10. :
	  exp(-1./2.*sq(z-9.)));
}

double sum_gaussian (double z) {
  return exp(-sq(z + 6.)) + 3./4.*exp(-sq(z + 3.)) +
    2./3.*exp(-sq(z)) + 1./2.*exp(-sq(z - 3.)) + 1./3.*exp(-sq(z - 6.));
}

/**
## Approximation of the integral */

double integral (double (* f)(double), double za, double zb)
{
  int N = 50;
  double integral = 0., zi[N];
  for (int n = 0; n < N; n++)
    zi[n] = za + n*(zb - za)/(N - 1.);
  for (int n = 0; n < N-1; n++)
    integral += (f(zi[n + 1]) + f(zi[n]))/2.*(zi[n + 1] - zi[n]);
  integral /= zb - za;
  return integral;
}

/**
## Comparison function */

void comparison (double (* f) (double),
		 int nposinit, int nposend, const double * zinit, const double * zend,
		 double f_b, double lambda_b, double f_t, double lambda_t,
		 int cell_lim)
{
  int ndof = 1, nvar = 1;
  // number of point = npos - 1 because those are mean values on each interval
  double init_pc[nposinit - 1], init_pf[nposinit - 1], init2_p[nposinit - 1];
  double end_pc[nposend - 1], end_pf[nposend - 1];

  for (int i = 0; i < nposinit - 1; i++) {
    init_pc[i] = integral (f, zinit[i], zinit[i+1]);
    init_pf[i] = integral (f, zinit[i], zinit[i+1]);
    init2_p[i] = init_pc[i];
  }

  int Nremap = 100; // number of remap
  for (int n = 0; n < Nremap; n++) {
    remap_c (nposinit, nposend, zinit, zend, init_pc, end_pc,
             f_b, lambda_b, 0, f_t, lambda_t, 0, true);
    remap_c (nposend, nposinit, zend, zinit, end_pc, init_pc,
             f_b, lambda_b, 0, f_t, lambda_t, 0, true);
    int edge_meth = p3e_method, cell_meth = ppm_method;
    my_remap (&nposinit, &nposend, &ndof, &nvar, zinit, zend, init_pf, end_pf,
	      &edge_meth, &cell_meth, &cell_lim);
    my_remap (&nposend, &nposinit, &ndof, &nvar, zend, zinit, end_pf, init_pf,
	      &edge_meth, &cell_meth, &cell_lim);
  }

  fputs ("\n\n", stderr);
  for (int i = 0; i < nposinit - 1; i++)
    fprintf (stderr,"%d %g %g %g %g\n", i, zinit[i], init2_p[i], init_pc[i], init_pf[i]);
}

/**
## main */

int main()
{

  /**
  We remap 60 initial bins to 54. */
  
  int nposinit = 60, nposend = 54;
  double zinit[nposinit], zend[nposend];

  for (int i = 0; i < nposinit; i++)
    zinit[i] = -10. + i*(10. - (-10.))/(nposinit - 1.);
  zend[0] = -10.;
  zend[nposend - 1] = 10.;
  for (int i = 1; i < nposend - 1; i++)
    zend[i] = -10. + i*(10 - (-10.))/(nposend - 1.) + 0.1*noise();
  
  comparison (polynomial,
	      nposinit, nposend, zinit, zend,
	      100., 0., 100., 0, null_limit);
  comparison (sharp,
	      nposinit, nposend, zinit, zend,
	      0., 0., 0., HUGE, mono_limit);
  comparison (sum_gaussian,
	      nposinit, nposend, zinit, zend,
	      0., 0., 0., HUGE, mono_limit);
  comparison (polynomial,
	      nposinit, nposend, zinit, zend,
	      40., - 3., 80., 1., null_limit);

  /**
  ~~~gnuplot Polynomial
  plot "log" index 0 using 2:3 title "init" with boxes fc 'gray', \
       "" index 0 using 2:4 title "remap C" lc 'red', \
       "" index 0 using 2:5 title "remap Fortran" lc 'black' pt 6
  ~~~

  ~~~gnuplot Piecewise constant
  plot "log" index 1 using 2:3 title "init" with boxes fc 'gray', \
       "" index 1 using 2:4 title "remap C" lc 'red', \
       "" index 1 using 2:5 title "remap Fortran" lc 'black' pt 6
  ~~~

  ~~~gnuplot Sum of Gaussians
  plot "log" index 2 using 2:3 title "init" with boxes fc 'gray', \
       "" index 2 using 2:4 title "remap C" lc 'red', \
       "" index 2 using 2:5 title "remap Fortran" lc 'black' pt 6
  ~~~
  
  ~~~gnuplot Polynomial Robin condition
  plot "log" index 3 using 2:3 title "init" with boxes fc 'gray', \
       "" index 3 using 2:4 title "remap C" lc 'red', \
       "" index 3 using 2:5 title "remap Fortran" lc 'black' pt 6
  ~~~
  */
}
