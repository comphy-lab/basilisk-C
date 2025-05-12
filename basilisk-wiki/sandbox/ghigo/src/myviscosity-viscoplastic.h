/**
# Regularization method for a yield stress fluid

We define the yield stress *T* and the regularization parameters. */

(const) face vector T = zerof;
double eps = 1.e-4;

/**
When dealing with multiphase flows, we can define *T* using an
arithmetic average. The user can overload this definition. */

#ifndef T
#define T(f,T1,T2) (clamp(f,0.,1.)*(T1 - T2) + T2)
#endif

/**
We also define a scalar quantity to characterize the yielded and
unyielded regions of the domain: 0 for yielded, 1 for unyielded. Note
that when using a regularization method, this quantity is
post-processed and is not computed directly. */

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
#if TREE
  yuyz.refine = yuyz.prolongation = refine_yuyz;
  yuyz.restriction = restriction_yuyz;
#endif // TREE
}

event init (i = 0)
{
  foreach()
    yuyz[] = 0.; // 0 : yielded 
}

/**
## Regularization 

#### Helpful functions 

To compute the second invariant of the stress tensor $S_2$ on each
face of the grid, we define helpfull functions to compute gradients in
the directions orthogonal to the face normal. */

#if EMBED
foreach_dimension()
double yield_face_avg_gradient_t1_x (Point point, scalar a, int i) {
  double left = nodata, right = nodata;
  if (cs[i,0] && (emerged || csm1[i,0]))
    right = (fs.y[i,0] && fs.y[i,1] &&
	     cs[i,1] && cs[i,-1] &&
	     (emerged || (csm1[i,1] && csm1[i,-1])) ?  (a[i,1] - a[i,-1])/(2.*Delta) :
	     fs.y[i,1] && fs.y[i,2] &&
	     cs[i,1] && cs[i,2] &&
	     (emerged || (csm1[i,1] && csm1[i,2])) ?   (-a[i,2] + 4.*a[i,1] - 3.*a[i,0])/(2.*Delta) :
	     fs.y[i,0] && fs.y[i,-1] &&
	     cs[i,-1] && cs[i,-2] &&
	     (emerged || (csm1[i,-1] && csm1[i,-2])) ? (a[i,-2] - 4.*a[i,-1] + 3.*a[i,0])/(2.*Delta) :
	     fs.y[i,1] && cs[i,1] &&
	     (emerged || csm1[i,1]) ?                  (a[i,1] - a[i])/Delta :
	     fs.y[i,0] && cs[i,-1] &&
	     (emerged || csm1[i,-1]) ?                 (a[i] - a[i,-1])/Delta : nodata);
  if (cs[i-1,0] && (emerged || csm1[i-1,0]))
    left = (fs.y[i-1,0] && fs.y[i-1,1] &&
	     cs[i-1,1] && cs[i-1,-1] &&
	     (emerged || (csm1[i-1,1] && csm1[i-1,-1])) ?  (a[i-1,1] - a[i-1,-1])/(2.*Delta) :
	     fs.y[i-1,1] && fs.y[i-1,2] &&
	     cs[i-1,1] && cs[i-1,2] &&
	     (emerged || (csm1[i-1,1] && csm1[i-1,2])) ?   (-a[i-1,2] + 4.*a[i-1,1] - 3.*a[i-1,0])/(2.*Delta) :
	     fs.y[i-1,0] && fs.y[i-1,-1] &&
	     cs[i-1,-1] && cs[i-1,-2] &&
	     (emerged || (csm1[i-1,-1] && csm1[i-1,-2])) ? (a[i-1,-2] - 4.*a[i-1,-1] + 3.*a[i-1,0])/(2.*Delta) :
	     fs.y[i-1,1] && cs[i-1,1] &&
	     (emerged || csm1[i-1,1]) ?                    (a[i-1,1] - a[i-1])/Delta :
	     fs.y[i-1,0] && cs[i-1,-1] &&
	     (emerged || csm1[i-1,-1]) ?                   (a[i-1] - a[i-1,-1])/Delta : nodata);  
  return (right == nodata && left == nodata ? 0. :
	  right == nodata ? left :
	  left == nodata ? right :
	  fs.x[i] ? (left + right)/2. : 0.);
}

foreach_dimension()
double yield_face_avg_gradient_t2_x (Point point, scalar a, int i) {
  double left = nodata, right = nodata;
  if (cs[i] && (emerged || csm1[i]))
    right = (fs.z[i,0,0] && fs.z[i,0,1] &&
	     cs[i,0,1] && cs[i,0,-1] &&
	     (emerged || (csm1[i,0,1] && csm1[i,0,-1])) ?  (a[i,0,1] - a[i,0,-1])/(2.*Delta) :
	     fs.z[i,0,1] && fs.z[i,0,2] &&
	     cs[i,0,1] && cs[i,0,2] &&
	     (emerged || (csm1[i,0,1] && csm1[i,0,2])) ?   (-a[i,0,2] + 4.*a[i,0,1] - 3.*a[i,0,0])/(2.*Delta) :
	     fs.z[i,0,0] && fs.z[i,0,-1] &&
	     cs[i,0,-1] && cs[i,0,-2] &&
	     (emerged || (csm1[i,0,-1] && csm1[i,0,-2])) ? (a[i,0,-2] - 4.*a[i,0,-1] + 3.*a[i,0,0])/(2.*Delta) :
	     fs.z[i,0,1] && cs[i,0,1] &&
	     (emerged || csm1[i,0,1]) ?                    (a[i,0,1] - a[i])/Delta :
	     fs.z[i,0,0] && cs[i,0,-1] &&
	     (emerged || csm1[i,0,-1]) ?                   (a[i] - a[i,0,-1])/Delta : nodata);
  if (cs[i-1] && (emerged || csm1[i-1]))
    left = (fs.z[i-1,0,0] && fs.z[i-1,0,1] &&
	     cs[i-1,0,1] && cs[i-1,0,-1] &&
	     (emerged || (csm1[i-1,0,1] && csm1[i-1,0,-1])) ?  (a[i-1,0,1] - a[i-1,0,-1])/(2.*Delta) :
	     fs.z[i-1,0,1] && fs.z[i-1,0,2] &&
	     cs[i-1,0,1] && cs[i-1,0,2] &&
	     (emerged || (csm1[i-1,0,1] && csm1[i-1,0,2])) ?   (-a[i-1,0,2] + 4.*a[i-1,0,1] - 3.*a[i-1,0,0])/(2.*Delta) :
	     fs.z[i-1,0,0] && fs.z[i-1,0,-1] &&
	     cs[i-1,0,-1] && cs[i-1,0,-2] &&
	     (emerged || (csm1[i-1,0,-1] && csm1[i-1,0,-2])) ? (a[i-1,0,-2] - 4.*a[i-1,0,-1] + 3.*a[i-1,0,0])/(2.*Delta) :
	     fs.z[i-1,0,1] && cs[i-1,0,1] &&
	     (emerged || csm1[i-1,0,1]) ?                      (a[i-1,0,1] - a[i-1])/Delta :
	     fs.z[i-1,0,0] && cs[i-1,0,-1] &&
	     (emerged || csm1[i-1,0,-1]) ?                     (a[i-1] - a[i-1,0,-1])/Delta : nodata);  
  return (right == nodata && left == nodata ? 0. :
	  right == nodata ? left :
	  left == nodata ? right :
	  fs.x[i] ? (left + right)/2. : 0.);
}
#else
foreach_dimension()
double yield_face_avg_gradient_t1_x (Point point, scalar a, int i) { // avg dady, dy = 2*Delta
  double right = (a[i,1]   - a[i,-1])/(2.*Delta);
  double left  = (a[i-1,1] - a[i-1,-1])/(2.*Delta);
  return (left + right)/2.;
}

foreach_dimension()
double yield_face_avg_gradient_t2_x (Point point, scalar a, int i) { // avg dadz, dz = 2*Delta
  double right = (a[i,0,1]   - a[i,0,-1])/(2.*Delta);
  double left  = (a[i-1,0,1] - a[i-1,0,-1])/(2.*Delta);
  return (left + right)/2.;
}
#endif // EMBED


/**
The following function updates the variable *yuyz* that characterizes
yielded and unyielded regions. */

void yielded_region ()
{
  /**
  We first compute the second invariant of the stress tensor $S_2$. */

  face vector S2[];
  foreach_face() {
    double D_11 = face_gradient_x (u.x, 0); // uxx
    double D_22 = yield_face_avg_gradient_t1_x (point, u.y, 0); // avg of uyy
    double D_12 = 0.5*(face_gradient_x (u.y,0) +
		       yield_face_avg_gradient_t1_x (point, u.x, 0)); // 0.5*(uyx + avg uxy)
    S2.x[] = (sq (D_11) + sq (D_22) + 2.*sq (D_12));
#if dimension == 3 
    double D_33 = yield_face_avg_gradient_t2_x (u.z, 0); // avg uzz
    double D_13 = 0.5*(face_gradient_x (u.z,0) +
		       yield_face_avg_gradient_t2_x (point, u.x, 0)); // 0.5*(uzx + avg uxz)
    double D_23 = 0.5*(yield_face_avg_gradient_t1_x (point, u.z, 0) +
		       yield_face_avg_gradient_t2_x (point, u.y, 0)); // 0.5*(avg uzy + avg uyz)
    S2.x[] += (sq (D_33) + 2.*sq (D_13) + 2.*sq (D_23));
#endif // dimension == 3
  }
#if AXI
  foreach_face(x)
    S2.x[] += sq ((u.y[] + u.y[-1])/(2.*y));
  foreach_face(y)
    S2.y[] += sq ((u.y[] + u.y[0,-1])/(2.*y + 1.e-20)); // y=r, so avoid r=0
#endif //AXI

#if EMBED
  /**
  As $S_2$ is computed on each face of the domain, we do not account
  for the stress tensor on the embedded boundary. */
#endif // EMBED

  foreach_face()
    S2.x[] = mu.x[]*(2.*sqrt (S2.x[]/2.));

  boundary ((scalar *) {S2});

  /**
  Finally, we use *S2* to determine if a cell is yielded or
  unyielded. */
  
  foreach() {

    /**
    By default, the cell is unyielded. */

    yuyz[] = 1;

    /**
    We compute the surface-weighted average yield stress and invariant
    stress tensor in the cell. */

    double Ta = 0., Sa = 0., fa = 0.;
    foreach_dimension() {
      Ta += T.x[]  + T.x[1];
      Sa += S2.x[] + S2.x[1];
      fa += fm.x[] + fm.x[1];
    }
    Ta /= (fa + 1.e-20);
    Sa /= (fa + 1.e-20);
  
    if (cm[] <= 0. || (Ta && Sa <= Ta))
      yuyz[] = 0;
  }
  boundary ({yuyz});
}

/**
#### Modification of the viscosity */

event properties (i++)
{
  /**
  We compute here the square of second invariant *D2* (using the
  Frobenius norm) of the strain rate tensor (or deformation
  tensor). Note that *D2* is not multiplied by the metric. */

  face vector D2[];
  foreach_face() {
    double D_11 = face_gradient_x (u.x, 0); // uxx
    double D_22 = yield_face_avg_gradient_t1_x (point, u.y, 0); // avg of uyy
    double D_12 = 0.5*(face_gradient_x (u.y,0) +
		       yield_face_avg_gradient_t1_x (point, u.x, 0)); // 0.5*(uyx + avg uxy)
    D2.x[] = (sq (D_11) + sq (D_22) + 2.*sq (D_12));
#if dimension == 3 
    double D_33 = yield_face_avg_gradient_t2_x (u.z, 0); // avg uzz
    double D_13 = 0.5*(face_gradient_x (u.z,0) +
		       yield_face_avg_gradient_t2_x (point, u.x, 0)); // 0.5*(uzx + avg uxz)
    double D_23 = 0.5*(yield_face_avg_gradient_t1_x (point, u.z, 0) +
		       yield_face_avg_gradient_t2_x (point, u.y, 0)); // 0.5*(avg uzy + avg uyz)
    D2.x[] += (sq (D_33) + 2.*sq (D_13) + 2.*sq (D_23));
#endif // dimension == 3
  }

  /**
  When using cylindrical coordinates, we add the term
  $u_{\theta\theta} = u_r/r$. */

#if AXI
  foreach_face(x)
    D2.x[] += sq ((u.y[] + u.y[-1])/(2.*y));
  foreach_face(y)
    D2.y[] += sq ((u.y[] + u.y[0,-1])/(2.*y + 1.e-20)); // y=r, so avoid r=0
#endif //AXI

#if EMBED
  /**
  As $D_2$ is computed on each face of the domain, we do not account
  for the stress tensor on the embedded boundary. */
#endif // EMBED  

  boundary ((scalar *) {D2});

  /**
  We now modify the viscosity to account for yielded regions using a
  regularization method. */

  face vector muv = mu;

  foreach_face() {

    /**
    We introduce here a regularized viscosity, depending on the value
    of $||D_2||=\sqrt{D_2/2}$. Note here that the yield stress *T* is
    multiplied by the metric *fm*. */

    double DD = sqrt (D2.x[]/2.);

    muv.x[] += T.x[]*min (1./eps, 1./(2.*DD + 1.e-20)*(1. - exp (-2.*DD/eps)));
  }
  boundary ((scalar *) {muv});
}

/**
Finally, we determine if a cell is yielded or unyielded. */

event end_timestep (i++)
{
  yielded_region ();
}
