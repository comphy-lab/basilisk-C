/**
This file aims to introduce, for the multilayer solver, basic characterisation functions of the interface : curvature, height metric, relative area and contact angle boundary condition. */

/**
# Curvature
To compute the curvature, we estimate the derivatives of the height
functions in a given direction (*x* or *y*). We then compute the curvature as
$$
\kappa = \frac{\eta_{xx}}{(1 + \eta_x^2)^{3/2}}
$$
in two dimensions, or
$$
\kappa = \frac{\eta_{xx}(1 +\eta _y^2) + \eta_{yy}(1 + h_x^2) - 2\eta_{xy}\eta_x\eta_y}
{(1 + \eta_x^2 + \eta_y^2)^{3/2}}
$$
in three dimensions.*/

static double kappa_c (Point point)
{
#if dimension == 1 // real dimension == 2
  double hx = (eta[1] - eta[-1])/(2.*Delta);
  double hxx = (eta[1] + eta[-1] - 2.*eta[])/sq(Delta);
  double kappa = hxx/pow(1. + sq(hx), 3/2.);

#else // real dimension == 3
  double hx = (eta[1] - eta[-1])/(2.*Delta);
  double hy = (eta[0,1] - eta[0,-1])/(2.*Delta);
  double hxy = (eta[1,1] + eta[-1,-1] - eta[1,-1] - eta[-1,1])/(4.*sq(Delta));

  /**
  We "filter" the curvature using a weighted sum of the three
  second-derivatives in the $x$ and $y$ directions. This is necessary
  to avoid a numerical mode when the curvature is used to compute
  surface tension. */
  
  double filter = 0.2; // Is this really necessary in our case ??
  double hxx = (filter*(eta[1,1] + eta[-1,1] - 2.*eta[0,1]) +
		(eta[1] + eta[-1] - 2.*eta[]) +
		filter*(eta[1,-1] + eta[-1,-1] - 2.*eta[0,-1]))/
    ((1. + 2.*filter)*sq(Delta));
  double hyy = (filter*(eta[1,1] + eta[1,-1] - 2.*eta[1]) +
		(eta[0,1] + eta[0,-1] - 2.*eta[]) +
		filter*(eta[-1,1] + eta[-1,-1] - 2.*eta[-1]))/
    ((1. + 2.*filter)*sq(Delta));
  
  double kappa = (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
#endif

    /**
     We limit the maximum curvature to $1/\Delta$. */
    if (fabs(kappa) > 1./Delta)
      kappa = (kappa > 0 ? 1 : -1)/Delta;
  return kappa;
}

/**
# Surface length and surface area
This function computes the surface distance compared to the cell width $\Delta$ and the surface ara.

$$
h_m = \frac{1}{\Delta} \cdot \sqrt{\Delta^2 + \delta \eta^2} = \sqrt{1 + \frac{\delta \eta^2}{\Delta^2}}
$$
*/

double hm (Point point, bool centered)
{
  int i = centered ? 1 : 0;
  double hx = (eta[i] - eta[-1])/((1. + i)*Delta);
  return sqrt(1. + sq(hx));
}

/**
$$
\mathcal{A} = h_{m,x} \cdot\, h_{m, y}
$$
*/


double area_m (Point point) {
  double area = 1;
  foreach_dimension () {
    double l = hm(point, true);
    double ki = kappa_c(point);
    
    if (ki == 0)
      area *= l;
    else 
      area *= asin(l*ki*Delta/2) * 2./(ki*Delta);
  }
  return area;
}

/**
#Contact Angle 
*/
#define contact_angle(theta, Delta) \
  val(_s) + Delta/tan(theta)