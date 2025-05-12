#ifndef _POISSON_CONJUGATE_H
#define _POISSON_CONJUGATE_H

#if USE_MY_SOLID
extern scalar is_solid;
extern face vector is_solid_face;
extern void boundarySolidNeummanNoauto (scalar s);
#endif
#if USE_CONJUGATE_HEAT
void fixSmallCellsTemp(scalar TL, scalar TG);
void fixSmallCellsTempLevel(scalar TL, scalar TG, int l);
void imposeTempBoundaryConjugate(scalar TL, scalar TG, scalar TS, bool with_Rcc);
void imposeTempBoundaryConjugateLevel(scalar TL, scalar TG, scalar TS, int l);
#endif

void my_mg_cycle (scalar * a, scalar * res, scalar * da,
	       void (* relax) (scalar * da, scalar * res, 
			       int depth, void * data),
	       void * data,
	       int nrelax, int minlevel, int maxlevel)
{

  /**
  We first define the residual on all levels. */

  restriction (res);

  /**
  We then proceed from the coarsest grid (*minlevel*) down to the
  finest grid. */

  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {

    /**
    On the coarsest grid, we take zero as initial guess. */

    if (l == minlevel)
      foreach_level_or_leaf (l)
        for (scalar s in da)
          foreach_blockf (s)
            s[] = 0.;
    else
      foreach_level (l)
        for (scalar s in da)
          foreach_blockf (s)
            s[] = bilinear (point, s);

    boundary_level (da, l);

#if(USE_CONJUGATE_HEAT == 2)
    //fixSmallCellsTempLevel(a[0], a[1], l);
    imposeTempBoundaryConjugateLevel(da[0], da[1], da[2], l);
#endif

    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
#if(USE_CONJUGATE_HEAT == 2)
    //fixSmallCellsTempLevel(a[0], a[1], l);
    imposeTempBoundaryConjugateLevel(da[0], da[1], da[2], l);
#endif
    }
  }

  /**
  And finally we apply the resulting correction to *a*. */

  foreach() {
    scalar s, ds;
    for (s, ds in a, da)
      foreach_blockf (s)
	      s[] += ds[];
  }
#if(USE_CONJUGATE_HEAT == 2)
  //fixSmallCellsTemp(a[0], a[1]);
  imposeTempBoundaryConjugate(a[0], a[1], a[2], false);
#endif
}

struct myMGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
		       void * data);
  void (* relax) (scalar * da, scalar * res, int depth, 
		  void * data);
  void * data;
  
  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats my_mg_solve (struct myMGSolve p)
{

  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);

  /**
  The boundary conditions for the correction fields are the
  *homogeneous* equivalent of the boundary conditions applied to
  *a*. */

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];
  
  /**
  We initialise the structure storing convergence statistics. */

  mgstats s = {0};
  double sum = 0.;
  foreach (reduction(+:sum))
    for (scalar s in p.b)
      sum += s[];
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;
  
  /**
  Here we compute the initial residual field and its maximum. */

  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);

  /**
  We then iterate until convergence or until *NITERMAX* is reached. Note
  also that we force the solver to apply at least one cycle, even if the
  initial residual is lower than *TOLERANCE*. */

  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
        my_mg_cycle (p.a, res, da, p.relax, p.data,
            s.nrelax,
            p.minlevel,
            grid->maxdepth);
        s.resa = p.residual (p.a, p.b, res, p.data);

    /**
    We tune the number of relaxations so that the residual is reduced
    by between 2 and 20 for each cycle. This is particularly useful
    for stiff systems which may require a larger number of relaxations
    on the finest grid. */

#if 1
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
	      s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
	      s.nrelax--;
    }
#else
    if (s.resa == resb) /* convergence has stopped!! */
      break;
    if (s.resa > resb/1.1 && p.minlevel < grid->maxdepth)
      p.minlevel++;
#endif

    resb = s.resa;
  }
  s.minlevel = p.minlevel;
  
  /**
  If we have not satisfied the tolerance, we warn the user. */

  if (s.resa > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr, 
	     "WARNING: convergence for %s not reached after %d iterations\n"
	     "  res: %g sum: %g nrelax: %d\n", v.name,
	     s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }
    
  /**
  We deallocate the residual and correction fields and free the lists. */

  if (!p.res)
    delete (res), free (res);
  delete (da), free (da);

  return s;
}

struct myPoisson {
  scalar * a, * b;
  face vector * alpha;
  scalar * lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};

/**
We can now write the relaxation function. We first recover the extra
parameters from the data pointer. */

static void myrelax (scalar * al, scalar * bl, int l, void * data)
{
  struct myPoisson * p = (struct myPoisson *) data;

  scalar a, b;
  face vector alpha;
  scalar lambda;

  for(a, b, alpha, lambda in al, bl, p->alpha, p->lambda)
  {
#if JACOBI
    scalar c[];
#else
    scalar c = a;
#endif
    
    /**
    We use the face values of $\alpha$ to weight the gradients of the
    5-points Laplacian operator. We get the relaxation function. */

    foreach_level_or_leaf (l) {
      double n = - sq(Delta)*b[], d = - lambda[]*sq(Delta);
      foreach_dimension() {
        n += alpha.x[1]*a[1] + alpha.x[]*a[-1];
        d += alpha.x[1] + alpha.x[];
      }

#if USE_MY_SOLID
      if(!d)
        c[] = b[] = 0.0;
      else
#endif
        c[] = n/d;

      if(isnan(c[]) || isinf(c[]))
      {
        c[] = b[] = 0.0;
      }
    }

    /**
    For weighted Jacobi we under-relax with a weight of 2/3. */
    
#if JACOBI
    foreach_level_or_leaf (l)
      a[] = (a[] + 2.*c[])/3.;
#endif
    
#if TRASH
    scalar a1[];
    foreach_level_or_leaf (l)
      a1[] = a[];
    trash ({a});
    foreach_level_or_leaf (l)
      a[] = a1[];
#endif
  }
}

static double myresidual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  struct myPoisson * p = (struct myPoisson *) data;
  scalar a, b, res;
  face vector alpha;
  scalar lambda;
  double maxres = 0.;
#if(USE_CONJUGATE_HEAT == 2)
  extern double fix_small_thd;
  imposeTempBoundaryConjugate(al[0], al[1], al[2], false);
#endif
  int num_var = 0;
  for (a, b, res, alpha, lambda in al, bl, resl, p->alpha, p->lambda)
  {
    num_var++;
#if TREE
    /* conservative coarse/fine discretisation (2nd order) */
    face vector g[];
    foreach_face()
        g.x[] = alpha.x[] * face_gradient_x(a, 0);
    foreach (reduction(max:maxres),nowarning)
    {
      res[] = b[] - lambda[] * a[];
      foreach_dimension()
          res[] -= (g.x[1] - g.x[]) / Delta;
      // if(num_var == 1 && f[] < fix_small_thd)
      // {
      //   res[] = 0.0;
      // }
      // if(num_var == 2 && f[] > (1.0 - fix_small_thd))
      // {
      //   res[] = 0.0;
      // }
#if(USE_CONJUGATE_HEAT == 2)
      if(num_var < 3)
      {
        res[] *= (1.0 - is_solid[]);
      }
      else
      {
        res[] *= is_solid[];
      }
#endif
      if (fabs(res[]) > maxres)
        maxres = fabs(res[]);
    }
#else // !TREE
    /* "naive" discretisation (only 1st order on trees) */
    foreach (reduction(max:maxres), nowarning)
    {
      res[] = b[] - lambda[] * a[];
      foreach_dimension()
          res[] += (alpha.x[0] * face_gradient_x(a, 0) -
                    alpha.x[1] * face_gradient_x(a, 1)) /
                   Delta;

      if (fabs(res[]) > maxres)
        maxres = fabs(res[]);
    }
#endif // !TREE
  }

#if(USE_CONJUGATE_HEAT == 2)
  //imposeTempBoundaryConjugate(resl[0], resl[1], resl[2]);
#endif

  return maxres;
}

mgstats mypoisson (struct myPoisson p)
{

  scalar lambda;
  face vector alpha;
  for(lambda, alpha in p.lambda, p.alpha)
  {
    restriction({alpha, lambda});
  }

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar * a = p.a;
  scalar * b = p.b;

  mgstats s = my_mg_solve (a, b, myresidual, myrelax,
			&p, p.nrelax, p.res, minlevel = max(1, p.minlevel));

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}

#endif //_POISSON_CONJUGATE_H