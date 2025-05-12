/**
# Augmented Lagrangian formulation for a yield-stress fluid

The approach is based on the augmented Lagrangian approach presented
in [Saramito, 2016](#saramito2016). It is slightly different from the
the one presented in [Vinay et al,
2005](http://doi.org/10.1016/j.jnnfm.2005.04.005).

We define the yield stress *T*, the pseudo deformation tensor $d=2D$
($D$ is the deformation tensor), the augmention parameter *r* and the
Lagrangian multiplier $\lambda$, equivalent to a pseudo plastic stress
tensor. Both $d$ and $\lambda$ are evaluated on the faces of the cells
but are not multiplied by the metric.

Next we summarize the steps we follow to solve the Navier-Stokes
equations for a viscoplastic fluid:

1. Advection
$$
\frac{u^{\star} - u^n}{\Delta t} = (u \cdot \nabla u)^{n + 1 / 2}.
$$
We then loop onto the following steps until convergence:
2. Modified viscosity
$$
\frac{u^{\star\star,k} - u^{\star}}{\Delta t} = \alpha \nabla \cdot
(r\mu 2D)^{\star\star,k} + \alpha \nabla \cdot (\lambda - r\mu
d)^{\star\star,k-1}.
$$
3. Acceleration
$$
u^{\star\star,k}_f = \overline{u^{\star\star,h}} + \Delta t a_f.
$$
4. Projection
$$
u^{n + 1,k}_f = P (u^{\star\star,k}_f, p^{n + 1,k}).
$$
5. Centered pressure gradient correction
$$
g^{n + 1,k} = \overline{a_f - \alpha_f \nabla p^{n + 1,k}}
$$
$$
u^{n + 1,k} = u^{\star\star,k} + \Delta tg^{n + 1,k}.
$$
6. Lagrangian variables

$$
d^k = 
\begin{array}{cc}
0 & \mathrm{if} \quad ||\lambda^{k-1} + r\mu(2D^k)|| < T \\
\left(1 - \frac{T}{||\lambda^{k-1} +
r\mu(2D^k)||}\right)\frac{\lambda^{k-1} + r\mu(2D^k)}{(1 + r)\mu} &
\mathrm{if} \quad ||\lambda^{k-1} + r\mu(2D^k)|| \geq T
\end{array}
$$
$$
\lambda^k = \lambda_{k-1} + r\mu*(2D^k - d^k).
$$
*/

(const) face vector T = zerof;

tensor d[], lambda[];
double r = 1.;
#define REPS (1.e-30) // Avoid division by zero

/**
We perform sub-cycling on the velocity and Lagrangian fields up to
convergence between each time iteration. */

int iAL = 1000;
double TOLERANCE_AL = 1.e-6; // Tolerance for *d-2D*
scalar pk[]; // Previous sub-cycle pressure
vector up[], uk[]; // Previous time and sub-cycle velocities

/**
When dealing with multiphase flows, we can define the yield stress *T*
using an arithmetic average. The user can overload these
definitions. */

#ifndef T
#define T(f,T1,T2) (clamp(f,0.,1.)*(T1 - T2) + T2)
#endif

/**
We also define a scalar quantity to characterize the yielded and
unyielded regions of the domain: 0 for yielded, 1 for unyielded. */

scalar yuyz[];

static inline void refine_yuyz (Point point, scalar s)
{
  refine_bilinear (point, s);
  foreach_child()
    s[] = s[] > 0. ? 1. : 0.;
}

static inline void restriction_yuyz (Point point, scalar s)
{
  restriction_average (point, s);
  s[] = s[] > 0. ? 1. : 0.;
}

/**
## Initialization */

event defaults (i = 0)
{
  /**
  We declare each compontent of the tensors *lambda* and *d* as face
  vectors. */
  
  foreach_dimension() {
    init_face_vector (lambda.x, NULL);
  }
  
#if TREE
  yuyz.refine = yuyz.prolongation = refine_yuyz;
  yuyz.restriction = restriction_yuyz;
#endif // TREE
}

/**
We initialize *d=0* and *lambda=0*. */

event init (i = 0)
{
  foreach_dimension() {
    face vector dx = d.x, lx = lambda.x;
    foreach_face() {
      dx.x[] = 0.;
      lx.x[] = 0.;
    }
    boundary_flux ({dx, lx});
  }
}

/**
## Augmented Lagragian algorithm

#### Usefull functions

The following function multiplies the viscosity by *a*. */

static void multiply_viscosity (double a)
{
  if (is_constant(mu.x))
    foreach_dimension()
      _constant[mu.x.i - _NVARMAX] *= a; // fixme -> mu.x[] += a;
  else {
    face vector muv = mu;
    foreach_face()
      muv.x[] *= a;
    boundary ((scalar *) {muv});
  }
}

/**
The next function compute $\max (d - 2D(U))$, used only when
sub-cycling is turned on. Note that this function does not work in
parallel. */

static double delta_D()
{
  tensor twoD[];
  double max = 0.;
  foreach_dimension() {
    face vector twoDx = twoD.x, dx = d.x;
    foreach_face(x) {
      twoDx.x[] = 2.*(u.x[] - u.x[-1])/Delta;
      double a = fabs (dx.x[] - twoDx.x[]);
      if (a > max)
	max = a;
    }
#if dimension > 1
    foreach_face(y) {
      twoDx.y[] = 2.*((u.x[] - u.x[0,-1]) + // 1/2*duxdy, dy=Delta 
		      ((u.y[1,-1] + u.y[1,0])/4. - // 1/2duydx, dx=2*Delta
		       (u.y[-1,-1] + u.y[-1,0])/4.))/(2.*Delta);
      double a = fabs (dx.y[] - twoDx.y[]);
      if (a > max)
	max = a;
    }
#endif
#if dimension > 2
    foreach_face(z) {
      twoDx.z[] = 2.*((u.x[] - u.x[0,0,-1]) + // 1/2*duxdz 
		      ((u.z[1,0,-1] + u.z[1,0,0])/4. -
		       (u.z[-1,0,-1] + u.z[-1,0,0])/4.))/(2.*Delta);
      double a = fabs (dx.z[] - twoDx.z[]);
      if (a > max)
	max = a;
    }
#endif
  }
  return max;
}

/**
This function updates both *d* and *lambda*. */

static void update_lambda_d()
{ 
  /**
  We first compute $D$ and $\lambda = \lambda + r\mu 2 D$. */

  tensor twoD[];
  foreach_dimension() {
    face vector lx = lambda.x, twoDx = twoD.x;
    foreach_face(x) {
      twoDx.x[] = 2.*(u.x[] - u.x[-1])/Delta;
      lx.x[] += r*(mu.x[]/(fm.x[] + REPS))*twoDx.x[];
    }
#if dimension > 1
    foreach_face(y) {
      twoDx.y[] = 2.*((u.x[] - u.x[0,-1]) + // 1/2*duxdy, dy=Delta 
		      ((u.y[1,-1] + u.y[1,0])/4. - // 1/2duydx, dx=2*Delta
		       (u.y[-1,-1] + u.y[-1,0])/4.))/(2.*Delta);
      lx.y[] += r*(mu.y[]/(fm.y[] + REPS))*twoDx.y[];
    }
#endif
#if dimension > 2
    foreach_face(z) {
      twoDx.z[] = 2.*((u.x[] - u.x[0,0,-1]) + // 1/2*duxdz 
		      ((u.z[1,0,-1] + u.z[1,0,0])/4. -
		       (u.z[-1,0,-1] + u.z[-1,0,0])/4.))/(2.*Delta);
      lx.z[] += r*(mu.z[]/(fm.z[] + REPS))*twoDx.z[];
    }
#endif
    boundary_flux ({twoDx});
  }
  boundary ((scalar *) {lambda});
  
  /**
  We now compute the Euclidian norm of the updated value of $\lambda$,
  multiplied by the metric. */
  
  face vector nlambda[];
  foreach_face() {
    nlambda.x[] = sqrt (0.5*(
			     sq (lambda.x.x[])
#if dimension > 1
			     + sq (lambda.y.y[]    + lambda.y.y[-1] + 
				   lambda.y.y[0,1] + lambda.y.y[-1,1])/16.
			     + 2.*sq (lambda.y.x[])
#endif
#if dimension > 2
			     //TODO
#endif
			     ));
    nlambda.x[] *= fm.x[];
  }

  /**
  Finally, we udapte $d$ and $\lambda = \lambda - r*mu*d$. */
  
  foreach_dimension() {
    face vector dx = d.x, lx = lambda.x, twoDx = twoD.x;
    foreach_face() {
      // Unyielded
      if (T.x[] == 0.) {
	dx.x[] = twoDx.x[];
	lx.x[] = mu.x[]/(fm.x[] + REPS)*twoDx.x[];
      }
      // Yielded
      else if (nlambda.x[] <= T.x[]) {
	dx.x[] = 0.;
      }
      // Unyielded
      else {
	double a = T.x[]/nlambda.x[];
	dx.x[] = (1. - a)*lx.x[]/((1. + r)*mu.x[]/(fm.x[] + REPS));
	lx.x[] -= r*(mu.x[]/(fm.x[] + REPS))*dx.x[];
      }
    }
    boundary_flux ({dx, lx});
  }
}

/**
Finally, the following function updates the variable *yuyz* that
characterizes yielded and unyielded regions. */

void yielded_region()
{
  /**
  We use *d* to determine if a cell is yielded or unyielded. */

  foreach()
    yuyz[] = 1. - (fabs (d.x.x[])    == 0. &&
  		   fabs (d.x.x[1])   == 0. &&
  		   fabs (d.y.y[])    == 0. &&
  		   fabs (d.y.y[0,1]) == 0.);
  boundary ((scalar *) {yuyz});
}

/**
#### Algorithm 

In the following, we rely on the event inheritance features of
Basilisk to implement the augmented Lagrangian algorithm. 

We first solve the viscous problem using the enhanced viscosity. */

event viscous_term (i++)
{
  /**
  We save the velocity field to which has been added the nonlinear
  advection terms. The vector *up* is then used when sub-cycling is
  turned on. */
  
  foreach()
    foreach_dimension()
      up.x[] = u.x[];
  boundary ((scalar *) {up});

  /**
  We multiply the the viscosity by $r$. */
  
  multiply_viscosity ((r));

  /**
  We then add the Lagrangian constraint to the centered
  velocity field. */
  
  foreach_dimension() {
    face vector dx = d.x, lx = lambda.x;
    foreach() {
      double div = 0.;
      foreach_dimension()
	div += (
		alpha.x[1]*(lx.x[1] - mu.x[1]/(fm.x[1] + REPS)*dx.x[1]) -
		alpha.x[]* (lx.x[]  - mu.x[]/(fm.x[] + REPS)*dx.x[])
		)/(cm[]*Delta);
      u.x[] += dt*div;
    }
  }
}

/**
We remove $r$ from viscosity the coefficient. */

event acceleration (i++)
{
  multiply_viscosity (1./(r)); 
}

/**
We finally update *lambda* and *d* after the projection step and
eventually perform sub-cycling. */

event end_timestep (i++)
{
  /**
  We update $\lambda$ and $d$. */
  
  update_lambda_d();

  /**
  We then compute the maximum difference for the velocity between two
  sub-cycle iterations as well as the pressure. */
  
  double du_x = 0., du_y = 0., du_z = 0.;
  foreach_dimension()
    du_x = change (u.x, up.x);
  double du = sqrt (sq (du_x) + sq (du_y) + sq (du_z));
  double dp = change (p, pk);
  
  foreach() {
    pk[] = p[];
    foreach_dimension()
      uk.x[] = u.x[];
  }
  boundary ((scalar *) {pk, uk});

  /**
  This is sub-cycling, which is turned on by default. This is used to
  refine the computation of *lambda* and *d*. */

  int j = 0;
  while ((delta_D() > TOLERANCE_AL || du > TOLERANCE || dp > TOLERANCE) &&
	 j++ < iAL) {

    /**
    We update the velocity. */
    
    foreach() {
      foreach_dimension()
	u.x[] = up.x[];
    }
    boundary ((scalar *) {u});
    event ("viscous_term");
    event ("acceleration");
    event ("projection");

    /**
    We update the difference between the velocity for to sub-cycle
    iterations. */
      
    double du_x = 0., du_y = 0., du_z = 0.;
    foreach_dimension()
      du_x = change (u.x, uk.x);
    du = sqrt (sq (du_x) + sq (du_y) + sq (du_z));
    dp = change (p, pk);
    
    foreach() {
      pk[] = p[];
      foreach_dimension()
	uk.x[] = u.x[];
    }
    boundary ((scalar *) {pk, uk});

    /**
    We update $\lambda$ and $d$. */
    
    update_lambda_d();
  }

  fprintf (ferr, "#AL: %d %g %g %d %g %g %g\n",
	   i, t, dt, j-1,
	   delta_D(), du, dp);
  fflush (ferr);

  /**
  Finally, we determine if a cell is yielded or unyielded. */
    
  yielded_region ();
}
