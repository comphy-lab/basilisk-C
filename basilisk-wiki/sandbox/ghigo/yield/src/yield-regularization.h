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

The following function updates the variable *yuyz* that characterizes
yielded and unyielded regions. */

#define center_gradient_x(a) ((a[1]     -     a[-1])/(2.*Delta))
#define center_gradient_y(a) ((a[0,1]   -   a[0,-1])/(2.*Delta))
#define center_gradient_z(a) ((a[0,0,1] - a[0,0,-1])/(2.*Delta))

#if EMBED
#undef center_gradient_x
#define center_gradient_x(a)						\
  (fs.x[] && fs.x[1] ? (a[1] - a[-1])/(2.*Delta) :			\
   fs.x[1] ? (a[1] - a[])/Delta :					\
   fs.x[]  ? (a[] - a[-1])/Delta : 0.)

#undef center_gradient_y
#define center_gradient_y(a)						\
  (fs.x[] && fs.x[1] ? (a[1] - a[-1])/(2.*Delta) :			\
   fs.x[1] ? (a[1] - a[])/Delta :					\
   fs.x[]  ? (a[] - a[-1])/Delta : 0.)

#undef center_gradient_z
#define center_gradient_z(a)						\
  (fs.x[] && fs.x[1] ? (a[1] - a[-1])/(2.*Delta) :			\
   fs.x[1] ? (a[1] - a[])/Delta :					\
   fs.x[]  ? (a[] - a[-1])/Delta : 0.)

#endif // EMBED

void yielded_region ()
{
  /**
  We first compute the second invariant of the stress tensor $S_2$. */

  face vector S2[];
  foreach_face() {
    double D_11 = face_gradient_x (u.x,0); //uxx
    double D_22 = 0.5*((fm.y[] && fm.y[0,1] ? (u.y[0,1] - u.y[0,-1])/(2.*Delta) :
			fm.y[0,1]           ? (u.y[0,1] - u.y[])/(Delta) :
			fm.y[]              ? (u.y[] - u.y[0,-1])/(Delta) : 0.)
		       +
		       (fm.y[-1] && fm.y[-1,1] ? (u.y[-1,1] - u.y[-1,-1])/(2.*Delta) :
			fm.y[-1,1]             ? (u.y[-1,1] - u.y[-1])/(Delta) :
			fm.y[-1]               ? (u.y[-1] - u.y[-1,-1])/(Delta) : 0.)); //avg of uyy, dy = 2*Delta
    double D_12 = 0.5*(face_gradient_x (u.y,0) +
		       0.5*((fm.y[] && fm.y[0,1] ? (u.x[0,1] - u.x[0,-1])/(2.*Delta) :
			     fm.y[0,1]           ? (u.x[0,1] - u.x[])/(Delta) :
			     fm.y[]              ? (u.x[] - u.x[0,-1])/(Delta) : 0.)
			    +
			    (fm.y[-1] && fm.y[-1,1] ? (u.x[-1,1] - u.x[-1,-1])/(2.*Delta) :
			     fm.y[-1,1]             ? (u.x[-1,1] - u.x[-1])/(Delta) :
			     fm.y[-1]               ? (u.x[-1] - u.x[-1,-1])/(Delta) : 0.))); //uyx + avg uxy, dy = 2*Delta
    S2.x[] = (sq (D_11) + sq (D_22) + 2.*sq (D_12));
#if dimension == 3 
    double D_33 = 0.5*((fm.z[] && fm.z[0,0,1] ? (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta) :
			fm.z[0,0,1]           ? (u.z[0,0,1] - u.z[])/(Delta) :
			fm.z[]                ? (u.z[] - u.z[0,0,-1])/(Delta) : 0.)
		       +
		       (fm.z[-1] && fm.z[-1,0,1] ? (u.z[-1,0,1] - u.z[-1,0,-1])/(2.*Delta) :
			fm.z[-1,0,1]             ? (u.z[-1,0,1] - u.z[-1])/(Delta) :
			fm.z[-1]                 ? (u.z[-1] - u.z[-1,0,-1])/(Delta) : 0.)); //avg uzz, dz = 2*Delta
    double D_13 = 0.5*(face_gradient_x (u.z,0) +
		       0.5*((fm.z[] && fm.z[0,0,1] ? (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta) :
			     fm.z[0,0,1]           ? (u.x[0,0,1] - u.x[])/(Delta) :
			     fm.z[]                ? (u.x[] - u.x[0,0,-1])/(Delta) : 0.)
			    +
			    (fm.z[-1] && fm.z[-1,0,1] ? (u.x[-1,0,1] - u.x[-1,0,-1])/(2.*Delta) :
			     fm.z[-1,0,1]             ? (u.x[-1,0,1] - u.x[-1])/(Delta) :
			     fm.z[-1]                 ? (u.x[-1] - u.x[-1,0,-1])/(Delta) : 0.))); //uzx + avg uxz, dz = 2*Delta
    double D_23 = 0.5*(0.5*((fm.y[] && fm.y[0,1] ? (u.z[0,1] - u.z[0,-1])/(2.*Delta) :
			     fm.y[0,1]           ? (u.z[0,1] - u.z[])/(Delta) :
			     fm.y[]              ? (u.z[] - u.z[0,-1])/(Delta) : 0.)
			    +
			    (fm.y[-1] && fm.y[-1,1] ? (u.z[-1,1] - u.z[-1,-1])/(2.*Delta) :
			     fm.y[-1,1]             ? (u.z[-1,1] - u.z[-1])/(Delta) :
			     fm.y[-1]               ? (u.z[-1] - u.z[-1,-1])/(Delta) : 0.)) +
		       0.5*((fm.z[] && fm.z[0,0,1] ? (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta) :
			     fm.z[0,0,1]           ? (u.y[0,0,1] - u.y[])/(Delta) :
			     fm.z[]                ? (u.y[] - u.y[0,0,-1])/(Delta) : 0.)
			    +
			    (fm.z[-1] && fm.z[-1,0,1] ? (u.y[-1,0,1] - u.y[-1,0,-1])/(2.*Delta) :
			     fm.z[-1,0,1]             ? (u.y[-1,0,1] - u.y[-1])/(Delta) :
			     fm.z[-1]                 ? (u.y[-1] - u.y[-1,0,-1])/(Delta) : 0.))); //avg uzy, dy = 2*Delta + avg uyz, dz = 2*Delta
    S2.x[] += (sq (D_33) + 2.*sq (D_13) + 2.*sq (D_23));
#endif // dimension == 3
  }
#if AXI
  foreach_face(x)
    S2.x[] += sq ((u.y[] + u.y[-1])/(2.*y));
  foreach_face(y)
    S2.y[] += sq ((u.y[] + u.y[0,-1])/(2.*y + 1.e-20)); // y=r, so avoid r=0
#endif //AXI

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

    if (constant(T.x) != 0.) {
      double Ta = 0., Sa = 0., fa = 0.;
      foreach_dimension() {
	Ta += T.x[]  + T.x[1];
	Sa += S2.x[] + S2.x[1];
	fa += fm.x[] + fm.x[1];
      }
      Ta /= (fa + 1.e-20);
      Sa /= (fa + 1.e-20);
  
      if (cs[] <= 0. || (Ta && Sa <= Ta))
	  yuyz[] = 0;
    }
  }
  boundary ((scalar *) {yuyz});
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
    double D_11 = face_gradient_x (u.x,0); //uxx
    double D_22 = 0.5*((fm.y[] && fm.y[0,1] ? (u.y[0,1] - u.y[0,-1])/(2.*Delta) :
			fm.y[0,1]           ? (u.y[0,1] - u.y[])/(Delta) :
			fm.y[]              ? (u.y[] - u.y[0,-1])/(Delta) : 0.)
		       +
		       (fm.y[-1] && fm.y[-1,1] ? (u.y[-1,1] - u.y[-1,-1])/(2.*Delta) :
			fm.y[-1,1]             ? (u.y[-1,1] - u.y[-1])/(Delta) :
			fm.y[-1]               ? (u.y[-1] - u.y[-1,-1])/(Delta) : 0.)); //avg of uyy, dy = 2*Delta
    double D_12 = 0.5*(face_gradient_x (u.y,0) +
		       0.5*((fm.y[] && fm.y[0,1] ? (u.x[0,1] - u.x[0,-1])/(2.*Delta) :
			     fm.y[0,1]           ? (u.x[0,1] - u.x[])/(Delta) :
			     fm.y[]              ? (u.x[] - u.x[0,-1])/(Delta) : 0.)
			    +
			    (fm.y[-1] && fm.y[-1,1] ? (u.x[-1,1] - u.x[-1,-1])/(2.*Delta) :
			     fm.y[-1,1]             ? (u.x[-1,1] - u.x[-1])/(Delta) :
			     fm.y[-1]               ? (u.x[-1] - u.x[-1,-1])/(Delta) : 0.))); //uyx + avg uxy, dy = 2*Delta
    D2.x[] = (sq (D_11) + sq (D_22) + 2.*sq (D_12));
#if dimension == 3 
    double D_33 = 0.5*((fm.z[] && fm.z[0,0,1] ? (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta) :
			fm.z[0,0,1]           ? (u.z[0,0,1] - u.z[])/(Delta) :
			fm.z[]                ? (u.z[] - u.z[0,0,-1])/(Delta) : 0.)
		       +
		       (fm.z[-1] && fm.z[-1,0,1] ? (u.z[-1,0,1] - u.z[-1,0,-1])/(2.*Delta) :
			fm.z[-1,0,1]             ? (u.z[-1,0,1] - u.z[-1])/(Delta) :
			fm.z[-1]                 ? (u.z[-1] - u.z[-1,0,-1])/(Delta) : 0.)); //avg uzz, dz = 2*Delta
    double D_13 = 0.5*(face_gradient_x (u.z,0) +
		       0.5*((fm.z[] && fm.z[0,0,1] ? (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta) :
			     fm.z[0,0,1]           ? (u.x[0,0,1] - u.x[])/(Delta) :
			     fm.z[]                ? (u.x[] - u.x[0,0,-1])/(Delta) : 0.)
			    +
			    (fm.z[-1] && fm.z[-1,0,1] ? (u.x[-1,0,1] - u.x[-1,0,-1])/(2.*Delta) :
			     fm.z[-1,0,1]             ? (u.x[-1,0,1] - u.x[-1])/(Delta) :
			     fm.z[-1]                 ? (u.x[-1] - u.x[-1,0,-1])/(Delta) : 0.))); //uzx + avg uxz, dz = 2*Delta
    double D_23 = 0.5*(0.5*((fm.y[] && fm.y[0,1] ? (u.z[0,1] - u.z[0,-1])/(2.*Delta) :
			     fm.y[0,1]           ? (u.z[0,1] - u.z[])/(Delta) :
			     fm.y[]              ? (u.z[] - u.z[0,-1])/(Delta) : 0.)
			    +
			    (fm.y[-1] && fm.y[-1,1] ? (u.z[-1,1] - u.z[-1,-1])/(2.*Delta) :
			     fm.y[-1,1]             ? (u.z[-1,1] - u.z[-1])/(Delta) :
			     fm.y[-1]               ? (u.z[-1] - u.z[-1,-1])/(Delta) : 0.)) +
		       0.5*((fm.z[] && fm.z[0,0,1] ? (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta) :
			     fm.z[0,0,1]           ? (u.y[0,0,1] - u.y[])/(Delta) :
			     fm.z[]                ? (u.y[] - u.y[0,0,-1])/(Delta) : 0.)
			    +
			    (fm.z[-1] && fm.z[-1,0,1] ? (u.y[-1,0,1] - u.y[-1,0,-1])/(2.*Delta) :
			     fm.z[-1,0,1]             ? (u.y[-1,0,1] - u.y[-1])/(Delta) :
			     fm.z[-1]                 ? (u.y[-1] - u.y[-1,0,-1])/(Delta) : 0.))); //avg uzy, dy = 2*Delta + avg uyz, dz = 2*Delta
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
    muv.x[] += min (T.x[]/eps, T.x[]/(2*DD + 1.e-20)*(1. - exp (-2.*DD/eps)));
  }
}

/**
Finally, we determine if a cell is yielded or unyielded. */

event end_timestep (i++)
{
  yielded_region ();
}
