/**
# Poisson Helmholtz 4th order solver
*/

#include "Laplacian.h"

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
	  s[] = bicubic(point, s);

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


int NITERMAX = 100;
double TOLERANCE = 1e-3;

typedef struct {
  int i;     
  double resb, resa;
  double sum;         
} mgstats;

trace
mgstats mg_solve (scalar * a, scalar * b, double (* residual) (scalar * a, scalar * b, scalar * res, void * data), void (* relax) (scalar * da, scalar * res, int depth, void * data), void * data)
{

  scalar * da = list_clone (a), * res = NULL;
  for (scalar s in a) {
    scalar r = new scalar;
    res = list_append (res, r);
  }

  for (int b = 0; b < nboundary; b++)
    for (scalar s in da)
      s.boundary[b] = s.boundary_homogeneous[b];

  mgstats s = {0, 0., 0.};
  double sum = 0.;
  foreach()
    for (scalar s in b)
      sum += s[];
  s.sum = sum;

  s.resb = s.resa = residual (a, b, res, data);

  for (s.i = 0; s.i < NITERMAX && (s.i < 1 || s.resa > TOLERANCE); s.i++) {
    mg_cycle (a, res, da, relax, data, 4, 0, grid->maxdepth);
    s.resa = residual (a, b, res, data);
  }

  if (s.i == NITERMAX)
    fprintf (ferr, 
	     "WARNING: convergence not reached after %d iterations\n"
	     "  res: %g sum: %g\n", 
	     NITERMAX, s.resa, s.sum), fflush (ferr);

  delete (res); free (res);
  delete (da);  free (da);
  
  return s;
}

struct Poisson {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
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
        n += ( alpha.x[1]*(a[-1,2]-27*a[0,2]+27*a[1,2]-a[2,2]) - alpha.x[]*(a[-2,2]-27*a[-1,2]+27*a[0,2]-a[1,2]) )*(-17.)/(5760.*24.);
        n += ( alpha.x[1]*(a[-1,1]-27*a[0,1]+27*a[1,1]-a[2,1]) - alpha.x[]*(a[-2,1]-27*a[-1,1]+27*a[0,1]-a[1,1]) )*(308.)/(5760.*24.);
        n += ( alpha.x[1]*(a[-1,0]+27*a[1,0]-a[2,0]) - alpha.x[]*(a[-2,0]-27*a[-1,0]-a[1,0]) )*(5178.)/(5760.*24.);
        n += ( alpha.x[1]*(a[-1,-1]-27*a[0,-1]+27*a[1,-1]-a[2,-1]) - alpha.x[]*(a[-2,-1]-27*a[-1,-1]+27*a[0,-1]-a[1,-1]) )*(308.)/(5760.*24.);
        n += ( alpha.x[1]*(a[-1,-2]-27*a[0,-2]+27*a[1,-2]-a[2,-2]) - alpha.x[]*(a[-2,-2]-27*a[-1,-2]+27*a[0,-2]-a[1,-2]) )*(-17.)/(5760.*24.);
        d += (27.*(alpha.x[1]+alpha.x[])*5178.)/(5760.*24.);        
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
   foreach_face(){
      g.x[] = 0;
      g.x[] += alpha.x[]*(a[-2,2]-27*a[-1,2]+27*a[0,2]-a[1,2] )*(-17.)/(5760.*24.*Delta);
      g.x[] += alpha.x[]*(a[-2,1]-27*a[-1,1]+27*a[0,1]-a[1,1] )*(308.)/(5760.*24.*Delta);
      g.x[] += alpha.x[]*(a[-2,0]-27*a[-1,0]+27*a[0,0]-a[1,0] )*(5178.)/(5760.*24.*Delta);
      g.x[] += alpha.x[]*(a[-2,-1]-27*a[-1,-1]+27*a[0,-1]-a[1,-1] )*(308.)/(5760.*24.*Delta);
      g.x[] += alpha.x[]*(a[-2,-2]-27*a[-1,-2]+27*a[0,-2]-a[1,-2] )*(-17.)/(5760.*24.*Delta);
    }
   boundary_flux({g});
   foreach(){
      res[] = b[]-lambda[]*a[];
      foreach_dimension()
          res[] += (g.x[]-g.x[1])/Delta;
       if (fabs(res[]) > maxres)
          maxres = fabs (res[]);
    }

#else

   foreach(){
      res[] = b[]-lambda[]*a[];
      foreach_dimension(){
         res[] -= ( alpha.x[1]*(a[-1,2]-27*a[0,2]+27*a[1,2]-a[2,2]) - alpha.x[]*(a[-2,2]-27*a[-1,2]+27*a[0,2]-a[1,2]) )*(-17.)/(5760.*24.*sq(Delta));
         res[] -= ( alpha.x[1]*(a[-1,1]-27*a[0,1]+27*a[1,1]-a[2,1]) - alpha.x[]*(a[-2,1]-27*a[-1,1]+27*a[0,1]-a[1,1]) )*(308.)/(5760.*24.*sq(Delta));
         res[] -= ( alpha.x[1]*(a[-1,0]-27*a[0,0]+27*a[1,0]-a[2,0]) - alpha.x[]*(a[-2,0]-27*a[-1,0]+27*a[0,0]-a[1,0]) )*(5178.)/(5760.*24.*sq(Delta));
         res[] -= ( alpha.x[1]*(a[-1,-1]-27*a[0,-1]+27*a[1,-1]-a[2,-1]) - alpha.x[]*(a[-2,-1]-27*a[-1,-1]+27*a[0,-1]-a[1,-1]) )*(308.)/(5760.*24.*sq(Delta));
         res[] -= ( alpha.x[1]*(a[-1,-2]-27*a[0,-2]+27*a[1,-2]-a[2,-2]) - alpha.x[]*(a[-2,-2]-27*a[-1,-2]+27*a[0,-2]-a[1,-2]) )*(-17.)/(5760.*24.*sq(Delta));
       }
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
  
  face vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction({alpha,lambda});

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ({a}, {b}, residual, relax, &p);

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}