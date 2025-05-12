/**
## Computing a line integral

We want to compute $S$, defined as

$$ S(x, y_1,z) = \int_{y_1}^{\mathrm{top}} s \mathrm{d}y.$$

Altough this is easy on a cartesian grid (because of the particular
allignment), a tree grid may pose a challenge. Here we present a few
approaches.

Note that the equation implies that $S(x, \mathrm{top}, z) = 0$,
giving a boundary condition at the top boundary:

~~~literatec
S[top] = dirichlet (0.); // or just 0? 
~~~

On this page we present a few iterative strategies that make use of
multigrid acceleration.
*/

#include "poisson.h"
/**
## Option 1

Solve:
 
$$\frac{\partial a}{\partial y} = -b,$$

via,

$$\frac{\partial^2 a}{\partial y^2} = -\frac{\partial b}{\partial y},$$

Using the Poisson solver with $\alpha = \{0, 1, 0\}$). Mind the proper
boundary conditions. Note that this approach is not very practical
because the finite `TOLERANCE` has an indirect relation to the error in
$S$, and the BC on the derivative becomes very important.*/

struct Integrate_dy {
  scalar a, b;
  (const) scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};

mgstats integrate_dy (struct Integrate_dy p) {
  scalar a = p.a, b = p.b;
  boundary ({b});
  scalar dbdy[];
  foreach()
    dbdy[] = -(b[0,1] - b[0,-1])/(2.*Delta);
  const face vector alpha[] = {0., 1., 0.};
  return poisson (a, dbdy, alpha, p.lambda, p.tolerance,
		  p.nrelax, p.minlevel, p.res);
}

/**
## Option 2

Alternatively, we discretize

$$\frac{\partial a}{\partial y} = -b,$$

with 2nd-order accuracy as;

~~~literatec  
(a[0,1] - a[0,-1])/2./Delta = -b[]; 
~~~

The residual function is:
*/

static double residual_int (scalar * al, scalar * bl, scalar * resl, void * data) {
  scalar a = al[0], b = bl[0], res = resl[0];
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order)?? */
  face vector g[];
  foreach_face(y) {
    g.y[] = face_gradient_y (a, 0);
  }
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] + (g.y[0,1] + g.y[])/2.;
#else // !TREE
    foreach (reduction(max:maxres)) { 
      res[] = -b[] - (a[0,1] - a[0,-1])/(2.*Delta);
#endif // !TREE    
      if (fabs (res[]) > maxres) 
	maxres = fabs (res[]);
      //printf ("%g %g %g %g\n", x, y, res[], -b[]);
    }
    boundary (resl);
    printf ("%g \n", maxres);
    return maxres;
  }
/**
This centered derivative formulation does not provide a suggestion for
the relaxation of a[]. But perhaps it does for its neighbors. If they
are not ghosts, we modify both so that the discrete system is
statisfied.

~~~literatec
a[0,1]  <- (a[0,1]  - b[]*2*Delta + a[0,-1])/2.;
a[0,-1] <- (a[0,-1] + b[]*2*Delta + a[0,1]) /2.;
~~~

We use a `marker` field to indicate ghosts.
*/

scalar marker[];
marker[bottom] = nodata;
marker[top] = nodata;
 
static inline void refine_nodata (Point point, scalar s) {
  foreach_child()
    s[] = nodata;
}
 
 static void relax_int (scalar * al, scalar * bl, int l, void * data) {
   scalar a = al[0], b = bl[0];
#if JACOBI
  scalar c[];
#else
  scalar c = a;
#endif
  foreach_level_or_leaf (l) {
    if (marker[0,1] != nodata && marker[0,-1] != nodata) {
      double a01 = a[0,1]; 
      c[0,1]  = (a[0,1]  - (b[]*2*Delta + a[0,-1]))/2.;
      c[0,-1] = (a[0,-1] + (b[]*2*Delta + a01))    /2.;
      //c[0,1]  = (- b[]*2*Delta + a[0,-1]);
      //c[0,-1] = (  b[]*2*Delta + a01);
    
    }
    /**
If a neighbor is a ghosts, we only modify the other.

~~~literatec
a[0, +/- 1] = -/+ b[]*2*Delta + a[0, -/+ 1];
~~~
    */
    else if (marker[0,1] != nodata) {
      c[0,1]  = -b[]*2*Delta + a[0,-1];
    } else {
      c[0,-1] = b[]*2*Delta + a[0,1];
    }
  }
#if JACOBI
  foreach_level_or_leaf (l)
    a[] = (a[] + 2.*c[])/3.;
#endif
}

/**
   We also provinde a user-interface to this iterative
   Multigrid-accelerated method. Again, it does not work.
 */
    
struct Integrate_dn {
  scalar a, b;
  (const) face vector alpha;
  (const) scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
};

mgstats integrate_dn (struct Integrate_dn p) {
  scalar a = p.a, b = p.b;
  marker.prolongation = refine_nodata;
  boundary ({marker});
  return mg_solve({a}, {b}, residual_int, relax_int,
		  &p, p.nrelax, p.res,
		  tolerance = p.tolerance, minlevel = p.minlevel);
}
 
/**
## Option 3 

We may also discretize with first-order accuracy:

~~~literatec
(a[0,1] - a[])/Delta = -b[];
~~~

The residual is;
 */

  static double residual_1st (scalar * al, scalar * bl, scalar * resl, void * data) {
    scalar a = al[0], b = bl[0], res = resl[0];
    double maxres = 0.;
#if TREE
    /* conservative coarse/fine discretisation (2nd order)?? */
    face vector g[];
    foreach_face(y) {
      g.y[] = face_gradient_y (a, 0);
    }
    boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] + (g.y[0,1]);
#else // !TREE
    foreach (reduction(max:maxres)) { 
      res[] = b[] + (a[0,1] - a[0])/(Delta);
#endif // !TREE    
      if (fabs (res[]) > maxres) 
	maxres = fabs (res[]);
      //printf ("%g %g %g %g\n", x, y, res[], -b[]);
    }
    boundary (resl);
    return maxres;
  }
  
  /**
And the relaxation procedure:
   */
  
  static void relax_1st (scalar * al, scalar * bl, int l, void * data) {
    scalar a = al[0], b = bl[0];
#if JACOBI
    scalar c[];
#else
    scalar c = a;
#endif
    foreach_level_or_leaf(l) 
      c[]  = b[]*Delta + a[0,1];
    
    
#if JACOBI
  foreach_level_or_leaf (l)
    a[] = (a[] + 2.*c[])/3.;
#endif
  }

  /**
     Due to the biased nature of the stencil, the convergence
     properties become sensitive to the cell-iteration sequence. For
     cases with a top boundary condition for $a$, it appears we
     benefit from a top-down iterator for the relaxation.

     Using this solution, this method works quite OK. 
  */
  
#include "topdown-iterator.h"
  
  trace
    mgstats integrate_1st (struct Integrate_dn p) {
    update_cache_f2();
    scalar a = p.a, b = p.b;
    return mg_solve({a}, {b}, residual_1st, relax_1st,
		    &p, p.nrelax, p.res,
		    tolerance = p.tolerance, minlevel = p.minlevel);
    update_cache_f();
  }
  