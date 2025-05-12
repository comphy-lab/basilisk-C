/**
# Test *fsolve()* using GSL

This is one of the multiroot test cases provided with the
gsl library. We want to solve the Rosenbrock system of
equations:
$$
  f_1(x,y) = a(1-x)
$$
$$
  f_2(x,y) = b(y-x^2)
$$
with $a=1,b=10$, knowing that the solution of this system
is $(x,y)=(1,1)$.
*/

/**
We have to define the use of gsl before including fsolve. */

#define USE_GSL
#include "run.h"
#include "fsolve-gsl.h"

/**
We define a struct with data to be passed to the non-linear
system function. */

typedef struct {
  double a, b;
} UserDataNls;

/**
We define the function with the non-linear system of
equations. For convience, when using different interfaces,
it is useful to define a common function that uses just
elementary C types. */

void function (const double * x, double * f, void * params) {
  UserDataNls * data = (UserDataNls *)params;
  double a = data->a;
  double b = data->b;

  f[0] = a*(1. - x[0]);
  f[1] = b*(x[1] - x[0]*x[0]);
}

/**
The gsl interface requires the use of a specific function
which can be implemented as follows. */

int function_gsl_interface (const gsl_vector * x, void * params, gsl_vector * f) {
  double * xdata = x->data;
  double * fdata = f->data;

  function (xdata, fdata, params);
  return GSL_SUCCESS;
}

int main (void) {
  /**
  We create an array with the first guess solution
  of the non-linear system. */

  Array * arrUnk = array_new();
  {
    double val1 = 10., val2 = 5.;
    array_append (arrUnk, &val1, sizeof(double));
    array_append (arrUnk, &val2, sizeof(double));

    double * unks = (double *)arrUnk->p;
    for (int i=0; i<2; i++)
      fprintf (stderr, "firstguess[%d] = %f\n", i, unks[i]);
  }

  /**
  We create the parameters to be passed to the nls
  function. */

  UserDataNls data;
  data.a = 1.;
  data.b = 10.;

  /**
  We solve the non-linear system, callig the function
  *fsolve()*. */

  fsolve (function_gsl_interface, arrUnk, &data);

  /**
  We recover the results of the system. */
  {
    double * unks = (double *)arrUnk->p;
    for (int i=0; i<2; i++)
      fprintf (stderr, "results[%d] = %f\n", i, unks[i]);
  }

  /**
  Cleanup operations. */

  array_free (arrUnk);
  run();
}

