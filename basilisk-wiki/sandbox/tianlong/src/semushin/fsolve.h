/**
This file is copied from [Edorado's sandbox](http://basilisk.fr/sandbox/ecipriano/README). All credit to him!
*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#pragma autolink -lgsl -lgslcblas

typedef int (* nls_fun) (const gsl_vector * x, void * params, gsl_vector * f);

void fsolve_gsl (nls_fun fun,
    Array * arrUnk,
    void * params)
{
  const gsl_multiroot_fsolver_type * T;
  gsl_multiroot_fsolver * s;

  int status, iter = 0.;

  int size = arrUnk->len / sizeof(double);
  const size_t n = (size_t)(size);

  gsl_multiroot_function f = {fun, n, params};

  double * x_init = (double *)arrUnk->p;
  gsl_vector * x = gsl_vector_alloc (n);

  for (unsigned int i=0; i<size; i++)
    gsl_vector_set (x, i, x_init[i]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);

    if (status)   /* check if solver is stuck */
      break;

    status =
      gsl_multiroot_test_residual (s->f, 1.e-7);
  }
  while (status == GSL_CONTINUE && iter < 1000);

  double * res = (double *)arrUnk->p;
  for (unsigned int i=0; i<size; i++)
    res[i] = gsl_vector_get (s->x, i);

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
}

void fsolve (nls_fun fun,
    Array * arrUnk,
    void * params)
{
  fsolve_gsl (fun, arrUnk, params);
}