/**
# Interface Gradients

This file is copied from [Edorado's sandbox](http://basilisk.fr/sandbox/ecipriano/README). All credit to him!

The calculation of the interface gradients can be performed
using the method developed in [embed.h](/src/embed.h). This
module was extended from [gradients.h](/sandbox/ggennari/phase_change/gradients.h)
and it reports the same interface gradients calculation
implemented in [embed.h](/src/embed.h), together with the
vof-averaged interface gradients proposed by [Fleckenstein and Bothe](#fleckenstein2015volume):

![VOF-averaged normal gradient scheme](/src/figures/dirichlet_gradient.svg)

$$
\left(\dfrac{\partial f}{\partial \mathbf{n}_\Gamma}\right)
=
\left(c\dfrac{f_\Gamma - f_0}{d_0} + \left(1 - c\right)\dfrac{f_\Gamma - f_1}{d_1}\right)
$$
*/

/**
## Copy of embed

Part of the embedded gradient functions are copied here, with
different names in order to eventually use the gradients also
with the [embed.h](/src/embed.h) module.
*/


#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double concentration_gradient_x (Point point, scalar s, scalar cs, face vector fs,
					   coord n, coord p, double bc,
					   bool third)
{
  foreach_dimension()
    n.x = - n.x;
  double d[2], v[2] = {nodata,nodata};
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	    !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
  if (v[0] == nodata) {

    /**
    This is a degenerate case, we set the gradient to zero. */
	
    return 0.;
  }

  /**
  For non-degenerate cases, the gradient can be  obtained using
  second-order, third-order or vof-averaged estimates. */
 
  if (v[1] != nodata && third) { // third-order gradient
    return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
  }
  else if (v[1] != nodata) {
    return (cs[]*(bc - v[0])/d[0] + (1. - cs[])*(bc - v[1])/d[1])/Delta; //vof-avg
  }
  return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double concentration_gradient (Point point, scalar s, scalar cs, face vector fs,
				coord n, coord p, double bc, bool third)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return concentration_gradient_x (point, s, cs, fs, n, p, bc, third);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return concentration_gradient_x (point, s, cs, fs, n, p, bc, third);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return concentration_gradient_y (point, s, cs, fs, n, p, bc, third);
  return concentration_gradient_z (point, s, cs, fs, n, p, bc, third);
#endif // dimension == 3
  return nodata;
}

coord my_facet_normal (Point point, scalar c, face vector s)
{
  if (s.x.i >= 0)
  { // compute normal from face fractions
    coord n;
    double nn = 0.;
    foreach_dimension()
    {
      n.x = s.x[] - s.x[1];
      nn += fabs(n.x);
    }
    if (nn > 0.)
    {
      foreach_dimension()
          n.x /= nn;
      return n;
    }
    else
      return interface_normal(point, c);
  }
  return interface_normal(point, c);
}

/**
## *ebmgrad()*: high-level interface for the calculation of the interface gradients:
* *tr*: scalar fields whose gradients must be computed
* *fL*: liquid phase volume fraction (*f*)
* *fG*: gas phase volume fraction (*1 - f*)
* *fsL*: face fraction for *fL*, computed using [fracface.h](fracface.h)
* *fsG*: face fraction for *fG*, computed using [fracface.h](fracface.h)
* *inverse*: true if tracer is in gas phase, false otherwise
* *trint*: interface value
* *success*: deprecated (fixme)
*/

double ebmgrad (Point point,
                  scalar tr,
                  scalar fL,
                  scalar fG,
                  face vector fsL,
                  face vector fsG,
                  bool inverse,
                  double trint,
                  bool* success)
{
  coord m = interface_normal(point, fL);
#if dimension == 2
  coord p = {0.,0.};
#else //dimension ==3
  coord p = {0.,0.,0.};
#endif

  double alpha = plane_alpha (fL[], m);
  plane_area_center (m, alpha, &p);
  normalize (&m);
  
#if dimension == 2
  coord n = {0., 0.};
#else
  coord n = {0., 0., 0.};
#endif

  if (inverse) {
    foreach_dimension()
      n.x = -m.x;
  }
  else {
    foreach_dimension()
      n.x = m.x;
  }

  bool third = false;
#ifdef INTGRAD_3rd
  third = true;
#endif

  double dirgrad = 0.;
  if (inverse)
    dirgrad = concentration_gradient (point, tr, fG, fsG, n, p, trint, third);
  else
    dirgrad = concentration_gradient (point, tr, fL, fsL, n, p, trint, third);

  return dirgrad;
}

/**
## References

~~~bib
@article{fleckenstein2015volume,
  title={A volume-of-fluid-based numerical method for multi-component mass transfer with local volume changes},
  author={Fleckenstein, Stefan and Bothe, Dieter},
  journal={Journal of Computational Physics},
  volume={301},
  pages={35--58},
  year={2015},
  publisher={Elsevier}
}
~~~
*/