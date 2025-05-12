/**
# Poisson-4th order solver. 

This is the implementation of the Poisson-4th Order
*/

#if LIMITED
  #include "Adapt_Limiters.h"
#else
  #include "Adapt_No_Limiters.h"
#endif

void mg_cycle (scalar * a, scalar * res, scalar * da, void (* relax) (scalar * da, scalar * res, int depth, void * data), void * data, int nrelax, int minlevel, int maxlevel)
{

  restriction (res);

  for (int l = minlevel; l <= maxlevel; l++) {

    if (l == minlevel)
      foreach_level_or_leaf (l)
	for (scalar s in da)
	  s[] = 0.;

    else
      foreach_level (l)
	for (scalar s in da)
          s[] = order5_refine ( point, s );

    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }

  foreach() {
    scalar s, ds;
    for (s, ds in a, da)
      s[] += ds[];
  }
  boundary (a);
}

int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-12;

typedef struct {
  int i;              // number of iterations
  double resb, resa;  // maximum residual before and after the iterations
  double sum;         // sum of r.h.s.
  int nrelax;         // number of relaxations
} mgstats;

struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
		       void * data);
  void (* relax) (scalar * da, scalar * res, int depth, 
		  void * data);
  void * data;
  
  int nrelax;
  scalar * res;
};

mgstats mg_solve (struct MGSolve p)
{

  /**
  We allocate a new correction and residual field for each of the scalars
  in *a*. */

  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    for (scalar s in p.a) {
      scalar r = new scalar;
      res = list_append (res, r);
    }

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

  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > TOLERANCE);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data, s.nrelax, 0, grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);

    /**
    We tune the number of relaxations so that the residual is reduced
    by between 2 and 20 for each cycle. This is particularly useful
    for stiff systems which may require a larger number of relaxations
    on the finest grid. */
    
    if (s.resa > TOLERANCE) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
	s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
	s.nrelax--;
    }
    resb = s.resa;
  }

  /**
  If we have not satisfied the tolerance, we warn the user. */

  if (s.resa > TOLERANCE)
    fprintf (ferr, 
	     "WARNING: convergence not reached after %d iterations\n"
	     "  res: %g sum: %g nrelax: %d\n", 
	     s.i, s.resa, s.sum, s.nrelax), fflush (ferr);

  /**
  We deallocate the residual and correction fields and free the lists. */

  if (!p.res)
    delete (res), free (res);
  delete (da), free (da);

  return s;
}

struct Poisson {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
  int nrelax;
  scalar * res;
};

void relax (scalar * al, scalar * bl, int l, void * data){
  
   scalar a = al[0], b = bl[0];
   struct Poisson * p = data;
   (const) face vector alpha = p->alpha;
   (const) scalar lambda = p->lambda;

   scalar c[];
   double n,d;

   foreach_level_or_leaf(l){ 
  
     n = -sq(Delta)*b[];
     d = -sq(Delta)*lambda[];

     foreach_dimension(){
        n += ( alpha.x[1]*(a[-1,0] + 15.*a[1,0] - a[2,0]) - alpha.x[]*(a[-2,0] - 15.*a[-1,0] - a[1,0]) )/12.;       
        d += (15.*(alpha.x[1]+alpha.x[]))/(12.);        
       }

     c[] = n/d;
    }
    
   foreach_level_or_leaf(l)
      a[] = (a[] + 2.*c[])/3.;
}

double residual (scalar * al, scalar * bl, scalar * resl, void * data){

   scalar a = al[0], b = bl[0], res = resl[0];
   struct Poisson * p = data;
   (const) face vector alpha = p->alpha;
   (const) scalar lambda = p->lambda;
   double maxres = 0;

#if TREE
 
   face vector g[];
   foreach_face()
      g.x[] = alpha.x[]*(a[-2,0]- 15.*a[-1,0] + 15.*a[] - a[1,0] )/(12.*Delta);
   boundary_flux({g});

   foreach(reduction(max:maxres)){
      res[] = b[]-lambda[]*a[];
      foreach_dimension()
          res[] += (g.x[]-g.x[1])/Delta;
       if (fabs(res[]) > maxres)
          maxres = fabs (res[]);
    }

#else

   foreach(){
      res[] = b[]-lambda[]*a[];
      foreach_dimension()        
         res[] -= ( alpha.x[1]*(a[-1,0] - 15.*a[] + 15.*a[1,0] - a[2,0]) - alpha.x[]*(a[-2,0] - 15.*a[-1,0] + 15.*a[0,0] - a[1,0]) )/(12.*sq(Delta)); 
      if(fabs(res[]) > maxres)
         maxres = fabs (res[]);
   }
#endif
   boundary({res});
   return maxres;
} 

mgstats poisson (struct Poisson p)
{

  if (!p.alpha.x.i) {
    const vector alpha[] = {1.,1.,1.};
    p.alpha = alpha;
  }
  if (!p.lambda.i) {
    const scalar lambda[] = 0.;
    p.lambda = lambda;
  }

  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  face vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction ({alpha,lambda});

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ({a}, {b}, residual, relax, &p, p.nrelax, p.res);

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}

struct Project {
  face vector u;
  scalar p;
  face vector alpha; // optional: default unityf
  double dt;         // optional: default one
  int nrelax;        // optional: default four
};

trace
mgstats project (struct Project q)
{
  face vector u = q.u;
  scalar p = q.p;
  (const) face vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;

  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += u.x[1] - u.x[];
    div[] /= dt*Delta;
  }

  mgstats mgp = poisson (p, div, alpha,
			 tolerance = TOLERANCE/sq(dt), nrelax = nrelax);

  foreach_face()
    u.x[] -= dt*alpha.x[]*(p[-2] - 15.*p[-1] + 15.*p[] - p[1])/(12.*Delta);
  boundary ((scalar *){u});

  return mgp;
}
