/**
# Augmented Lagrangian formulation for yield-stress fluids

Based on [Vinay et al, 2005](http://doi.org/10.1016/j.jnnfm.2005.04.005). 

Subcycling is turned off by default (see below). */

(const) face vector tau_0 = zerof;
tensor d[], lambda[];
double r = 1.;

vector up[];

event defaults (i = 0) {
  // each component of lambda must be a face vector to get consistent
  // boundary conditions
  foreach_dimension()
    init_face_vector (lambda.x, NULL);
}

static void add_to_viscosity (double a)
{
  if (is_constant(mu.x))
    foreach_dimension()
      _constant[mu.x.i - _NVARMAX] += a; // fixme -> mu.x[] += a;
  else {
    face vector muv = mu;
    foreach_face()
      muv.x[] += a;
  }  
}

event viscous_term (i++) {
  // save u
  foreach()
    foreach_dimension()
      up.x[] = u.x[];
  // add r/2 to viscosity coefficient
  add_to_viscosity (r/2.);
  
  // add Lagrangian constraint (eq. 3.35)
  foreach_dimension() {
    face vector lx = lambda.x, dx = d.x;
    foreach() {
      double d = 0.;
      foreach_dimension()
	d += lx.x[1] - lx.x[] - r*(dx.x[1] - dx.x[]);
      u.x[] += dt/rho[]*d/Delta; // fixme: metric
    }
  }
}

// compute the maximum difference between d and D
// only used for subcycling
static double delta_D()
{
  double max = 0.;
  foreach_dimension() {
    face vector dx = d.x;
    foreach_face(x) {
      double a = fabs (dx.x[] - (u.x[] - u.x[-1])/Delta);
      if (a > max)
	max = a;
    }
    #if dimension > 1
    foreach_face(y) {
      double a = fabs (dx.y[] - 
		       (u.x[] - u.x[0,-1] + 
			(u.y[1,-1] + u.y[1,0])/4. -
			(u.y[-1,-1] + u.y[-1,0])/4.)/(2.*Delta));
      if (a > max)
	max = a;
    }
    #endif
    #if dimension > 2
    foreach_face(z) {
      double a = fabs (dx.z[] -
		       (u.x[] - u.x[0,0,-1] + 
			(u.z[1,0,-1] + u.z[1,0,0])/4. -
			(u.z[-1,0,-1] + u.z[-1,0,0])/4.)/(2.*Delta));
      if (a > max)
	max = a;
    }
    #endif
  }
  return max;
}

static void update_lambda_d()
{
  // compute \lambda = \lambda + rD (eq. 3.37)
  foreach_dimension() {
    face vector lx = lambda.x;
    foreach_face(x)
      lx.x[] += r*(u.x[] - u.x[-1])/Delta;
    #if dimension > 1
    foreach_face(y)
      lx.y[] += r*(u.x[] - u.x[0,-1] + 
		   (u.y[1,-1] + u.y[1,0])/4. -
		   (u.y[-1,-1] + u.y[-1,0])/4.)/(2.*Delta);
    #endif
    #if dimension > 2
    foreach_face(z)
      lx.z[] += r*(u.x[] - u.x[0,0,-1] + 
		   (u.z[1,0,-1] + u.z[1,0,0])/4. -
		   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/(2.*Delta);
    #endif
  }
  
  // compute the Euclidean norm of \lambda
  face vector nlambda[];
  foreach_face()
    nlambda.x[] = sqrt(sq(lambda.x.x[])/2. +
		       sq(lambda.y.y[] + lambda.y.y[-1] +
			  lambda.y.y[0,1] + lambda.y.y[-1,1])/32. +
		       sq(lambda.y.x[]));
  
  // compute d_i^k (eq. 3.36) and
  // \lambda = \lambda - rd (eq. 3.37)
  foreach_dimension() {
    face vector dx = d.x, lx = lambda.x;
    foreach_face() {
      if (nlambda.x[] <= tau_0.x[])
	dx.x[] = 0.;
      else {
	double a = tau_0.x[]/nlambda.x[];
	dx.x[] = (1. - a)*lx.x[]/r;
	lx.x[] *= a;
      }
    }
  }
}

event end_timestep (i++) {
  update_lambda_d();

  /**
  This is subcycling, which is turned off by default. */

  int j = 0;
  while (delta_D() > 1e-6 && j++ < 0) {
    fprintf (stderr, "max: %g %d %g %d\n", t, i, delta_D(), j - 1);

    /**
    We can also choose to "sub-sub-cycle" on the projection. */
    
    int k = 0;
    do {
      // restore u
      foreach()
	foreach_dimension()
          u.x[] = up.x[];
      boundary ((scalar *){u});
      event ("viscous_term");
      event ("acceleration");
      event ("projection");
    } while (k++ < 0);
      
    update_lambda_d();
  }

  // remove r/2 from viscosity coefficient
  add_to_viscosity (- r/2.);
}
