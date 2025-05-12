/**
# Non-Linear System of Equations Solver

This module defines a function *fsolve()* which, in analogy with
the MATLAB function, provides a high-level interface for the solution
of non-linear systems of equations.
*/

/**
## GSL Interface

We use *gsl_multiroots* to solve the non-linear system of
equations. */

#ifdef USE_GSL

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#pragma autolink -lgsl -lgslcblas

#ifndef FSOLVE_RELTOL
# define FSOLVE_RELTOL 0.
#endif

#ifndef FSOLVE_ABSTOL
# define FSOLVE_ABSTOL 1.e-7
#endif

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

    if (status)   /* check if solver is stuck */ {
      fprintf (stderr, "WARNING: Non linear systems solver is stuck.\n");
      break;
    }

    if (FSOLVE_RELTOL > 0.)
      status =
        gsl_multiroot_test_delta (s->dx, s->x, FSOLVE_ABSTOL, FSOLVE_RELTOL);
    else
      status =
        gsl_multiroot_test_residual (s->f, FSOLVE_ABSTOL);

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

#endif // USE_GSL


/**
## KINSol Interface

if the [SUNDIALS library](https://github.com/LLNL/sundials) is used,
the function *fsolve()* relies on the KINSol solver. This
implementation works with Sundials 5.8 and it is not updated
for Sundials >= 6.0. */

#ifdef USE_SUNDIALS

#include <kinsol/kinsol.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sundials/sundials_types.h>

#define KIN_FTOL   1.e-6 // function tolerance 
#define KIN_STOL   1.e-6 // step tolerance

#pragma autolink -lsundials_kinsol

typedef int (* nls_fun)(N_Vector u, N_Vector fval, void *user_data);

void fsolve_sundials (nls_fun fun,
    Array * arrUnk,
    Point point)
{
  N_Vector u, s;
  SUNMatrix J;
  SUNLinearSolver LS;

  u = NULL;
  s = NULL;
  J = NULL;

  int size = arrUnk->len / sizeof(double);

  u = N_VNew_Serial (size);
  s = N_VNew_Serial (size);
  {
    realtype * udata = N_VGetArrayPointer (u);
    double * unk = (double *)arrUnk->p;

    for (unsigned int jj=0; jj<size; jj++)
      udata[jj] = unk[jj];
  }
  N_VConst (1.0, s);

  void * kmem;
  kmem = KINCreate();
  KINSetUserData (kmem, params);
  KINSetFuncNormTol (kmem, KIN_FTOL);
  KINSetScaledStepTol (kmem, KIN_STOL);
  KINInit (kmem, fun, u);
  J = SUNDenseMatrix (size, size);
  LS = SUNLinSol_Dense (u, J);
  KINSetLinearSolver (kmem, LS, J);
  KINSetMaxSetupCalls(kmem, 1);
  KINSetPrintLevel (kmem, 0);

  //TODO following two lines give problems
  //{
  //  char name[80];
  //  sprintf (name, "KINSolErr");
  //  const FILE * fKinErr = fopen (name, "a");
  //  KINSetErrFile (kmem, fKinErr);
  //}

  /**
  Solve non-linear system of equations. */

  KINSol (kmem, u, KIN_NONE, s, s);
  //KINSol (kmem, u, KIN_LINESEARCH, s, s);

  /**
  Recover Nls solution. */

  {
    realtype * udata = N_VGetArrayPointer (u);
    double * unk = (double *)arrUnk->p;
    for (unsigned jj=0; jj<size; jj++) {
      unk[jj] = udata[jj];
    }
  }

  /**
  Free memory. */

  N_VDestroy (u);
  N_VDestroy (s);
  KINFree (&kmem);
  SUNLinSolFree (LS);
  SUNMatDestroy (J);
  free (data);
}

void fsolve (nls_fun fun,
    Array * arrUnk,
    Point point)
{
  fsolve_sundials (fun, arrUnk, point);
}

#endif // USE_SUNDIALS

