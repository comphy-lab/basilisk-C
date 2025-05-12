

struct tr_Poisson {
  scalar a, b;
  (const) face vector alpha;
  (const) face vector beta;
  (const) scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;
#if EMBED
  bool (* embed_flux) (Point, scalar, vector, bool, double *);
#endif
};

/**
We can now write the relaxation function. We first recover the extra
parameters from the data pointer. */

static void tr_relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct tr_Poisson * p = (struct tr_Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) face vector beta = p->beta;
  (const) scalar lambda = p->lambda;

  /**
  We use either Jacobi (under)relaxation or we directly reuse values
  as soon as they are updated. For Jacobi, we need to allocate space
  for the new field *c*. Jacobi is useful mostly as it gives results
  which are independent of the order in which the cells are
  traversed. This is not the case for the simple traversal, which
  means for example that results will depend on whether a tree or
  a multigrid is used (because cells will be traversed in a different
  order). The same comment applies to OpenMP or MPI parallelism. In
  practice however Jacobi convergence tends to be slower than simple
  reuse. */
  
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
      n += alpha.x[1]*a[1] + alpha.x[]*a[-1]+0.5*Delta*( beta.x[1]*a[1] - beta.x[]*a[-1] );
      d += alpha.x[1] + alpha.x[] - 0.5*Delta*( beta.x[1] - beta.x[] );
	    //Subtract fluxes 
	     if (f[1] > 1. - 1e-6 && f[] < 1. - 1e-6){
	   n -= alpha.x[1]*a[1]*f[] +0.5*Delta*( beta.x[1]*a[1] )*f[];
	d -= alpha.x[1]*f[]  - 0.5*Delta*( beta.x[1]  )*f[];
	}
	    if (f[-1] > 1. - 1e-6 && f[] < 1. - 1e-6){
	   n -= alpha.x[]*a[-1]*f[] - 0.5*Delta*( beta.x[]*a[-1] )*f[];
	d -= alpha.x[]*f[]  + 0.5*Delta*( beta.x[]  )*f[];
	} 
	    
	     if (f[1] > 1e-6 && f[]<1e-6){
	   n -= alpha.x[1]*a[1] +0.5*Delta*( beta.x[1]*a[1] );
	d -= alpha.x[1]  - 0.5*Delta*( beta.x[1]  );
	}
	 if (f[-1] > 1e-6 && f[]<1e-6){
	   n -= alpha.x[]*a[-1] - 0.5*Delta*( beta.x[]*a[-1] );
	d -= alpha.x[]  + 0.5*Delta*( beta.x[]  );
	}
    }
    
#if EMBED
    if (p->embed_flux) {
      double c;
      if (p->embed_flux (point, a, alpha, true, &c))
	n -= c*sq(Delta);
      else
	n = 0., d = 1.;
    }
    if (!d)
      c[] = b[] = 0.;
    else
#endif // EMBED
     if (d!=0)
      c[] = n/d;
     else c[]=0;
  }

  /**
  For weighted Jacobi we under-relax by using a weight of 2/3. */
  
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

/**
The equivalent residual function is obtained in a similar way in the
case of a Cartesian grid, however the case of the tree mesh
requires more careful consideration... */

static double tr_residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct tr_Poisson * p = (struct tr_Poisson *) data;
  (const) face vector alpha = p->alpha;
  (const) face vector beta = p->beta;
  (const) scalar lambda = p->lambda;
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  face vector g[];
  foreach_face()
    g.x[] = alpha.x[]*face_gradient_x (a, 0) +
      beta.x[]*(a[] + a[-1])/2;
  boundary_flux ({g});
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension() {
      res[] -= (g.x[1] - g.x[])/Delta;
	//Subtract fluxes 
	if ( f[1] > 1. - 1e-6 && f[] < 1. - 1e-6){
	res[] += g.x[1]*f[]/Delta;
	}
	if ( f[-1] > 1. - 1e-6 && f[] < 1. - 1e-6){
	res[] += -g.x[]*f[]/Delta;
	}
	if ( f[1] > 1e-6 && f[]<1e-6){
	res[] += g.x[1]/Delta;
	}
	if ( f[-1] > 1e-6 && f[]<1e-6){
	res[] += -g.x[]/Delta;
	}
  }
#else // !TREE
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)) {
    res[] = b[] - lambda[]*a[];
    foreach_dimension()
      res[] -= (alpha.x[1]*face_gradient_x (a, 1) -
		alpha.x[0]*face_gradient_x (a, 0) +
		beta.x[1]*(a[1] + a[])/2.
		- beta.x[0]*(a[] + a[-1])/2.)/Delta;  	  
#endif // !TREE    
#if EMBED
    if (p->embed_flux) {
      double c;
      if (p->embed_flux (point, a, alpha, false, &c))
	res[] += c;
      else
	a[] = c, res[] = 0.;
    }
#endif // EMBED    
    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
  boundary (resl);
  return maxres;
}

/**
## User interface

Finally we provide a generic user interface for a Poisson--Helmholtz
equation of the form
$$
\nabla\cdot (\alpha\nabla a + \beta a) + \lambda a = b
$$ */


mgstats tr_poisson (struct tr_Poisson p)
{

  /**
  If $\alpha$ or $\lambda$ are not set, we replace them with constant
  unity vector (resp. zero scalar) fields. Note that the user is free to
  provide $\alpha$ and $\beta$ as constant fields. */

  if (!p.alpha.x.i)
    p.alpha = unityf;
   if (!p.beta.x.i)
    p.beta = zerof;
  if (!p.lambda.i)
    p.lambda = zeroc;
  

  /**
  We need $\alpha$ and $\lambda$ on all levels of the grid. */

  face vector alpha = p.alpha;
  face vector beta = p.beta;
  scalar lambda = p.lambda;
  restriction ({alpha,beta,lambda});

  /**
  If *tolerance* is set it supersedes the default of the multigrid
  solver. */

  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  //~ foreach(){
		 //~ double check1=a[];
		//~ if (check1>1e10){
	  //~ printf("ac");
	    //~ }
    //~ }
  mgstats s = mg_solve ({a}, {b}, tr_residual, tr_relax,
			&p, p.nrelax, p.res, minlevel = max(1, p.minlevel));

  /**
  We restore the default. */

  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}

struct tr_Diffusion {
  // mandatory
  scalar f;
  double dt;
  // optional
  face vector D;  // default 1
  face vector Dh;  // default 0
  scalar r, beta; // default 0
  scalar theta;   // default 1
};

trace
mgstats tr_diffusion (struct tr_Diffusion p)
{

  /**
  If *dt* is zero we don't do anything. */

  if (p.dt == 0.) {
    mgstats s = {0};
    return s;
  }

  /**
  We define $f$ and $r$ for convenience. */

  scalar f = p.f, r = automatic (p.r);

  /**
  We define a (possibly constant) field equal to $\theta/dt$. */
double dtp = p.dt;
  
   scalar theta_idt[];
   foreach(){
     theta_idt[] = - cm[]/dtp;
   }
  
   //~ scalar idt;
  //~ (const) scalar theta_idt = p.theta.i ? p.theta : idt;
 //~ scalar theta_idt;
  
  //~ if (p.theta.i) {
    //~ scalar theta_idt = p.theta;
    //~ foreach()
      //~ theta_idt[] *= idt[];
  //~ }
 //~ foreach(){
      //~ theta_idt[] =-1.0/p.dt;
 //~ }
  /**
  We use `r` to store the r.h.s. of the Poisson--Helmholtz solver. */

  if (p.r.i)
    foreach()
      r[] = theta_idt[]*f[] - r[];
  else // r was not passed by the user
    foreach()
      r[] = theta_idt[]*f[];

  /**
  If $\beta$ is provided, we use it to store the diagonal term $\lambda$. */

  scalar lambda = theta_idt;
  if (p.beta.i) {
    scalar beta = p.beta;
    foreach()
      beta[] += theta_idt[];
    lambda = beta;
  }

  boundary ({lambda});

  /**
  Finally we solve the system. */

  return tr_poisson (f, r, p.D,p.Dh, lambda);
}
