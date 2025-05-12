/**
# Integral formulation for surface tension

This is a beginning of implementation of the integral formulation for surface
tension depicted in
[Abu-Al-Saud, Popinet, Tchelepi, 2018](https://hal.archives-ouvertes.fr/hal-01706565/).
We adapt [integral.h](/sandbox/popinet/integral.h) (written by Stéphane in
Basilisk C but for a level set description of the interface) for a VOF
description, quasi 1D interface, with one set of height function.

In this version of the implementation of the integral formulation of the surface
tension, we extrapolate the interface in each cell as a circle arc. This removes
the spurious current of the contact mwe: [contact_circle.c](contact_circle.c)
is much better than [contact_line.c](contact_line.c). However, this
implementation using circles is just a first try, a lot of problem remains, 
notably it is not possible yet to have straight interface for the moment and the
prolongation and restriction fonction for the stress tensor are not written yet.
For a more versatile yet less accurate implementation see
[integral_vof_quasi1D.h](integral_vof_quasi1D.h).

We first copy some lines of [contact.h](/src/contact.h), but with an update of
the height funtions at the beginning of the tracer_advection event instead of
the end of vof, because at the end of VOF, with munal advection, interfaces is
still equal to NULL. Ideally, we should not have to do the advection manually
when we use f.tracers, and did not need to set interfaces to NULL and we could
use back contact.h.

////////////////////////////////////////////////////////////////////////////////

# Contact angles

This file is used to impose contact angles on boundaries for
interfaces described using a [VOF](vof.h) tracer and [height
functions](heights.h).

We first overload the default function used to compute the normal,
defined in [fractions.h](). */

coord interface_normal (Point point, scalar c);
#define interface_normal(point, c) interface_normal (point, c)

#include "fractions.h"
#include "curvature.h"

/**
We will compute the normal using height-functions instead. If this is
not possible (typically at low resolutions) we revert back to
the Mixed-Youngs-Centered approximation. */

coord interface_normal (Point point, scalar c)
{
  coord n;
  if (!c.height.x.i || (n = height_normal (point, c, c.height)).x == nodata)
    n = mycs (point, c);
  else {
    double nn = 0.;
    foreach_dimension()
      nn += fabs(n.x);
    foreach_dimension()
      n.x /= nn;
  }
  return n;
}

/**
The height functions are stored in the vector field associated with
each VOF tracer. They need to be updated every time the VOF field
changes. For the [centered Navier-Stokes
solver](navier-stokes/centered.h), this means after initialisation and
after VOF advection. 

Note that strictly speaking this should be done for each
[sweep](vof.h#sweep_x) of the direction-split VOF advection, which we
do not do here i.e. we use the normal at the beginning of the timestep
and assume it is constant during each sweep. This seems to work
fine. */

extern scalar * interfaces;

event init (i = 0) {
  for (scalar c in interfaces)
    if (c.height.x.i)
      heights (c, c.height);
}

event tracer_advection (i++) {
  for (scalar c in interfaces)
    if (c.height.x.i)
      heights (c, c.height);
}

/**
The macro below can be used to impose a contact angle on a boundary by
setting the corresponding tangential component of the height
function. 

Note that the equivalent function for the normal component of the
height function is not defined yet. This limits the range of
accessible contact angles, since values of the normal component of the
height function will be required to compute curvature at shallow
angles. */

double heights_to_curvature (Point point, scalar f) {
  if (interfacial (point, f)) {
    scalar hy = f.height.y;
    double yc = (2. + sq(hy[1]) + sq(hy[-1]) - 2.*sq(hy[]));
    yc /= (2.*(hy[-1] + hy[1] - 2.*hy[]));
    double xc = - (hy[] - hy[-1])*yc + 1./2.*(sq(hy[]) - sq(hy[-1]) -1.);
    return - sign(yc)/pow(sq(xc) + sq(hy[] - yc), 0.5)/Delta;
  }
  else
    return nodata;
}

void circle_curvature (scalar f, scalar kappa) {
  #if TREE
  kappa.refine = kappa.prolongation = curvature_prolongation;
  kappa.restriction = curvature_restriction;
  #endif
  foreach()
    kappa[] = heights_to_curvature (point, f);
  boundary({kappa});
}

double my_contact (Point point, double theta, bool left, scalar hy, scalar kappa) {
  int l = (left ? 1 : - 1);
  if (hy[l] == nodata)
    return nodata;
  else {
    double Dy = hy[] - hy[l];//y1 (x=0.5l) - y2 (x=1.5l)
    double Sy = hy[] + hy[l];
    double a = 1. - sq(cos(theta))*(1. + 1./sq(Dy));
    double c = 4. - (1. + sq(Dy))*sq(cos(theta));
    if (4. < a*c)
      return nodata;
    else {
      double alpha = -l*(2. + pow(4. - a*c, 0.5))/a;
      double x0 = l + alpha/2.;
      double R = -l*x0/cos(theta);
      double y0 = (alpha*l/Dy + Sy)/2.;
      if (fabs(R) < fabs(-l*0.5 - x0))
        return nodata;
      else
        return y0 - pow(sq(R) - sq(-l*0.5 - x0), 0.5);
    }
  }
}

#define contact_angle(theta, left) \
  (val(_s) == nodata ? nodata : my_contact(point, theta, left, _c.height.y, _c.kappa))

/**
////////////////////////////////////////////////////////////////////////////////

We associate to each VOF tracer a curvature field, a (potentially variable)
surface tension $\sigma$ and a contact angle $\theta_\mathrm{c}$ on the
different boundary (for the moment only left and right are available). This is
done easily by adding the following
[field attributes](/Basilisk C#field-attributes). */

attribute {
  scalar kappa;
  scalar sigma;
  double theta_left;
  double theta_right;
}

event defaults (i = 0) {
  
  /**
  Surface tension is a source term in the right-hand-side of the
  evolution equation for the velocity of the [centered Navier--Stokes
  solver](navier-stokes/centered.h) i.e. it is an acceleration. If
  necessary, we allocate a new vector field to store it. */

  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
}

/**
# Stability condition

The surface tension scheme is time-explicit so the maximum timestep is
the oscillation period of the smallest capillary wave.
$$
T = \sqrt{\frac{\rho_{m}\Delta_{min}^3}{\pi\sigma}}
$$
with $\rho_m=(\rho_1+\rho_2)/2$ and $\rho_1$, $\rho_2$ the densities
on either side of the interface. */

event stability (i++) {

  /**
  We first compute the minimum and maximum values of $\alpha/f_m =
  1/\rho$, as well as $\Delta_{min}$. */

  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin)) {
    if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
    if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
    if (Delta < dmin) dmin = Delta;
  }
  double rhom = (1./amin + 1./amax)/2.;

  /**
  We then consider each VOF interface with an associated value of
  $\sigma$ different from zero and set the maximum timestep. */
  
  for (scalar c in interfaces)
    if (c.sigma.i) {
      scalar sigma = c.sigma;
      if (is_constant(sigma)) {
	      double dt = sqrt (rhom*cube(dmin)/(pi*constant(sigma)));
	      if (dt < dtmax)
	        dtmax = dt;
      }
      else
	      foreach() { //foreach(reduction(max:dtmax)) {
	        double dt = sqrt (rhom*cube(dmin)/(pi*sigma[]));
	        if (dt < dtmax && sigma[] != nodata)
	          dtmax = dt;
	      }
    }
}

/**
# Functions

We define sereval functions to compute the tangents to the interface and one
to compute the weighted average of sereval variables possibliy equal to nodata.
*/

double do_average (double * values, double * weights, int length) {

  // int len = (sizeof values)/(sizeof (double));
  double weight_sum = 0.;
  double avg = 0.;
  for (int i = 0; i < length; i++) {
    weight_sum += (values[i] != nodata ? weights[i] : 0.);
    avg += (values[i] != nodata ? weights[i]*values[i] : 0.);
  }
  return (weight_sum != 0. ? avg/weight_sum : nodata);
}

coord normalize_Q (coord v) {
  double nv = 0.;
  foreach_dimension()
    nv += sq(v.x);
  nv = sqrt(nv);
  foreach_dimension()
    v.x = (nv != 0. ? v.x/nv : nodata);
  return v;
}

/**
R has to be normalized by Delta ! */

coord compute_circle (double R, coord A, coord B) {
  coord diff;
  foreach_dimension()
    diff.x = A.x - B.x;
  coord beta;
  beta.y = - sign(R)*fabs(diff.x)*pow(4.*sq(R)/(sq(diff.x) + sq(diff.y)) - 1., 0.5);
  beta.x = - beta.y*diff.y/diff.x;
  coord center;
  foreach_dimension()
    center.x = (beta.x + A.x + B.x)/2.;
  return center;
} 

double circle_x (double R, double hyl, double hyr, coord * center) {
  if (hyl == nodata || hyr == nodata)
    return nodata;
  else {
    double Dy = hyr - hyl;
    coord A = {0.5, hyr}; coord B = {-0.5, hyl};
    * center = compute_circle (R, A, B);
    if (fabs(R) < fabs(center->y))
      return nodata;
    else
      return (center->x - sign(R*Dy)*pow(sq(R) - sq(center->y), 0.5));
  }
}

double circle_y (double R, double hyl, double hyr, coord * center) {
  if (hyl == nodata || hyr == nodata)
    return nodata;
  else {
    coord A = {0.5, hyr}; coord B = {-0.5, hyl};
    * center = compute_circle (R, A, B);
    if (fabs(R) < fabs(center->x))
      return nodata;
    else
      return (center->y + sign(R)*pow(sq(R) - sq(center->x), 0.5));
  }
}

double Sxx_boundary (Point point, scalar hy, scalar kappa, scalar sigma,
                     bool left) {
      
  /**
  We check if the interface cross the vertical middle of the ghost cell with 
  fabs(hy[]) < 0.5. To take into account the case hy[] = 0.5 once we
  have to break the symetry : strict on one side and inclusive on the
  other. Since sign(x) = (x > 0 ? 1 : -1), the special hy == 0 will be
  handled as negative, then the logical choice is
  (hy[] < 0.5 && hy[] >= -0.5). */
  
  int ghost_Q = (left ? -1 : 1);
  
  if (hy[ghost_Q] < 0.5 && hy[ghost_Q] >= -0.5) {

    double xi = fabs(hy[ghost_Q]);
    
    /**
    Estimation of the curvature */
    
    int l = (hy[ghost_Q] > 0. ? 1 : -1);
    double second_value;
    second_value = (kappa[ghost_Q, l] != 0. && kappa[ghost_Q, l] != nodata ?
                    kappa[ghost_Q, l] :
                    coarse(kappa, ghost_Q, l) != 0. ? 
                    coarse(kappa, ghost_Q, l) : nodata);
    double kg = (kappa[ghost_Q] != nodata ? kappa[ghost_Q] :
                 second_value != nodata ? second_value : nodata);
    kg = (kg != nodata ? kg :
          kappa[ghost_Q, -l] != 0. && kappa[ghost_Q, -l] != nodata ?
          kappa[ghost_Q, -l] :
          coarse(kappa, ghost_Q, -l) != nodata ?
          coarse(kappa, ghost_Q, -l) : 0.);

    /**
    Computation of the tangents */
    
    coord A = {ghost_Q, hy[ghost_Q]}; coord B = {0., hy[0]};
    coord center = compute_circle (1./(kg*Delta), A, B);
    coord tangentg = {-(hy[ghost_Q] - center.y), ghost_Q-center.x};
    tangentg = normalize_Q (tangentg);
    
    /**
    Surface tension */
    
    second_value = (sigma[ghost_Q, l] != 0.
                    && sigma[ghost_Q, l] != nodata ? sigma[ghost_Q, l] :
                    coarse(sigma, ghost_Q, l) != 0. ?
                    coarse(sigma, ghost_Q, l) : nodata);
    double sigma_values[] = {sigma[ghost_Q], second_value};
    double weights[] = {1. - xi, xi};
    double sigmag = do_average (sigma_values, weights, 2);
    sigmag = (sigmag != nodata ? sigmag :
              sigma[ghost_Q, -l] != 0. && sigma[ghost_Q, -l] != nodata ?
              sigma[ghost_Q, -l] :
              coarse(sigma, ghost_Q, -l) != 0. ?
              coarse(sigma, ghost_Q, -l) : nodata);
    if (sigmag == nodata)
      return nodata;
    else {    
      return (sigmag*(fabs(tangentg.x)/Delta
                      + sign(hy[ghost_Q])*kg*(0.5 - xi)));
    }
  }
  return nodata;
}
/**
The prolongation function performs a similar averaging, but using the
same stencil as that used for bilinear interpolation, so that the
symmetries of the volume fraction field and curvature field are
preserved. */

scalar _c;

#if 0
static void stress_prolongation_x (Point point, scalar stress) {
  
  scalar hy = _c.height.y, kappa = _c.kappa, sigma = _c.sigma;
  coord tangent_parent_center = tangent_center (point, hy);
  /**
  coord tangent_parent_left = tangent_left (point, hy);  
  coord tangent_parent_right = tangent_left (neighborp(1), hy);
  */
  foreach_child() {
    //double hy = _c.height.y[];
    if (hy[] < 0.5 && hy[] >= -0.5) {
      double xi = fabs(hy[]);
      //int l = (hy[] > 0. ? 1 : -1);
      //double weights[] = {1. - xi, xi};
#if 1
      double ki = kappa[]; // to improve!
#elif LINEAR_CURVATURE // does not make much difference
      double values[] = {kappa[], kappa[0, l]};
      double ki = do_average (values, weights, 2);
#else
      double ki = (kappa[] != nodata ? kappa[] :
                   kappa[0, l] != nodata ? kappa[0, l] : 0.);
#endif
      coord tangenti = tangent_parent_center; // to improve!

#if 1
      double sigmai = sigma[]; // to improve!
#else
      double sigma_values[] = {sigma[], sigma[0, l]};
      double sigmai = do_average (sigma_values, weights, 2);
#endif
      stress[] = sigmai*(fabs(tangenti.x)/Delta + sign(hy[])*ki*(0.5 - xi));
    }
    else
      stress[] = nodata;
  }
}

static void stress_prolongation_y (Point point, scalar stress) {
  
  scalar hy = _c.height.y, kappa = _c.kappa, sigma = _c.sigma;
  coord tangent_parent_center = tangent_center (point, hy);
  /**
  coord tangent_parent_left = tangent_left (point, hy);  
  coord tangent_parent_right = tangent_left (neighborp(1), hy);
  */
  foreach_child() {
  
    stress[] = 0.;
    bool Syy_defined = false;
  
    for (int l = -1; l <= 1; l += 2) {
      
      /**
      We first prolongate hy to the neigbor (l, 0) if it not already defined
      (i.e. if (l, 0) is not a sibling cell). */ 
      double hy_l = nodata;
      if (l != child.x)
        hy_l = hy[l];
      else {
        hy_l = 0.75*coarse(hy, l) + 0.25*coarse(hy);
        hy_l = 2.*hy_l - 0.5*child.y;
      }
      if (hy[]*(hy[] + hy_l) <= 0.
          && (l < 0. || hy[] + hy_l != 0.)
          && hy[] != nodata && hy_l != nodata && hy[] != hy_l) {
        Syy_defined = true;
        double xi = hy[]/(hy[] - hy_l); // > 0
        //double weights[] = {1. - xi, xi};
  #if 1
        double ki = kappa[]; // to improve!
  #elif LINEAR_CURVATURE // does not make much difference
        double values[] = {kappa[], kappa[0, l]};
        double ki = do_average (values, weights, 2);
  #else
        double ki = (kappa[] != nodata ? kappa[] :
                     kappa[0, l] != nodata ? kappa[0, l] : 0.);
  #endif
        coord tangenti = tangent_parent_center; // to improve!

  #if 1
        double sigmai = sigma[]; // to improve!
  #else
        double sigma_values[] = {sigma[], sigma[0, l]};
        double sigmai = do_average (sigma_values, weights, 2);
  #endif
        stress[] = sigmai*(fabs(tangenti.y)/Delta + sign(hy[])*ki*(0.5 - xi));
      }
    }
    if (Syy_defined == false)
      stress[] = nodata;
  }
}
#endif // TREE

/**
# Surface tension term

The calculation of the acceleration is done by this event, overloaded
from [its definition](navier-stokes/centered.h#acceleration-term) in
the centered Navier--Stokes solver. */

/**
We set the boundary conditions on the height functions thanks to contact.h. */

event init(i = 0) {
  for (scalar c in interfaces) { 
    if (c.height.x.i) {
    
      /**
      c is not defined as a global variable. We reuse sucessively the global
      scalar _c to point toward each interface tracer. */ 

      _c = c;
      vector h = c.height;
      if (c.theta_left)
        h.t[left] = contact_angle (_c.theta_left, true);
      if (c.theta_right)
        h.t[right] = contact_angle (_c.theta_right, false);
    }
  }
}

event acceleration (i++) {
  
  /**
  We check whether surface tension is associated with any VOF tracer $c$. */

  for (scalar c in interfaces) {
    if (c.sigma.i) {

      /**
      c is not defined as a global variable. We reuse sucessively the global
      scalar _c to point toward each interface tracer. */ 

      _c = c;
      
      foreach()
        c[] = clamp(c[], 0., 1.);
      boundary({c});

      /**
      We define hy for convenience. */

      scalar hy;
      hy = c.height.y;
      boundary ({hy});
      
      /**
      We allocate kappa, compute the curvature and make c.kappa point to
      kappa (in order that _c.kappa points to kappa). */
      
      scalar kappa[];
      circle_curvature (c, kappa);
      c.kappa = kappa;

      /**
      We define sigma for convenience. */
      
      scalar sigma;
      sigma = c.sigma;
      boundary ({sigma});

      /**
      ## Diagonal compenents of the stress tensor
      
      We compute the diagonal components of the tensor. */

      tensor S[];
      
      for (scalar s in {S}) {
        foreach_dimension() {
          s[left] = neumann(0.);
          s[right] = neumann(0.);
        }
      }

      if (c.theta_left) {
        scalar s = S.x.x;
        s[left] = Sxx_boundary (point, _c.height.y, _c.kappa, _c.sigma, true);
      }
      if (c.theta_right) {
        scalar s = S.x.x;
        s[right] = Sxx_boundary (point, _c.height.y, _c.kappa, _c.sigma, false);
      }
      foreach() {
      
        /**
        S.x.x
        
        We check if the interface cross the vertical middle of the cell with 
        fabs(hy[]) < 0.5. To take into account the case hy[] = 0.5 once we
        have to break the symetry : strict on one side and inclusive on the
        other. Since sign(x) = (x > 0 ? 1 : -1), the special hy == 0 will be
        handled as negative, then the logical choice is
        (hy[] < 0.5 && hy[] >= -0.5). */
        
        if (hy[] < 0.5 && hy[] >= -0.5) {
          double xi = fabs(hy[]);
          int l = (hy[] > 0. ? 1 : -1);
          double ki = (kappa[] != nodata ? kappa[] :
                       kappa[0, l] != nodata ? kappa[0, l] : 0.);

          coord A = {-1., hy[-1]}; coord B = {1., hy[1]};
          coord center = compute_circle (1./(ki*Delta), A, B);
          coord tangenti = {-(hy[] - center.y), -center.x};
          tangenti = normalize_Q (tangenti);
 
          double weights[] = {1. - xi, xi};         
          double sigma_values[] = {sigma[], sigma[0, l]};
          double sigmai = do_average (sigma_values, weights, 2);
          
          S.x.x[] = sigmai*(fabs(tangenti.x)/Delta + sign(hy[])*ki*(0.5 - xi));
        }
        else
          S.x.x[] = nodata;

        /**
        S.y.y */
        
        S.y.y[] = 0.;
        
        /**
        We check if the interface cross the horizontal middle of the cell with
        (hy[]*(hy[] + hy[l]) < 0.). To take into account the special cases
        (hy[+/-1/2] = 0 and hy = 0), we include the cell if the interface is
        on the left, and exclude it if it is on the right. */
        bool Syy_defined = false;
        if (interfacial(point, f)) {
          for (int l = -1; l <= 1; l += 2) {
            double ki = (kappa[] != nodata ? kappa[] :
                       kappa[l] != nodata ?  kappa[l] : 0.);
            coord circle_center;
            double hx = circle_x (1./(ki*Delta), hy[(l - 1)/2], hy[(l + 1)/2],
                                  & circle_center);
            hx += l*0.5;
            if ((l == 1 && hx < 0.5 && hx >= 0.)
             || (l == -1 && hx >= -0.5 && hx < 0.)) {
              Syy_defined = true;
              double xi = fabs(hx); // > 0
              coord tangenti = {circle_center.y, (hx - l*0.5 - circle_center.x)};
              tangenti = normalize_Q (tangenti);

              double weights[] = {1. - xi, xi};              
              double sigma_values[] = {sigma[], sigma[l]};
              double sigmai = do_average (sigma_values, weights, 2);
              
              S.y.y[] += sigmai*(fabs(tangenti.y)/Delta + sign(hy[])*ki*(0.5 - xi));
            }
          }
        }
        if (Syy_defined == false)
          S.y.y[] = nodata;
      }

      /**
      We define the boundary conditions and, on trees, the prolongation and
      restriction functions for the diagonal terms. */

      /**
      If we want to apply a contact angle condition, we have to take care of the
      boundary on S.x.x and S.y.y. They cannot be handled simply with neumann()
      or dirichlet(), because, the cell where the stress is not null is not
      necessarly the ghost cell. */
      
      boundary ({S.x.x, S.y.y});

      /**
      ## Off-diagonal compenents of the stress tensor
      
      We compute the off-diagonal components of the tensor.  */
      
      foreach_vertex() {
        
        /**
        S.y.x

        We test if the interface cross the face i = -0.5, j in [-0.5, 0.5]
        thanks to (fabs(hyv) < 0.5) with hyv = (hy[] + hy[-1, 0])/2. + 0.5,
        the hight from the vertex.
        To take into account the case where hyv = 0.5, we rather test 
        (hyv <= 0.5 && hyv > -0.5). To not have allocate hyv if it is not
        useful, we use its expression:
        (hy[] + hy[-1, 0] <= 0. && hy[] + hy[-1, 0] > -2.). */

        double weights[] = {1., 1.};
        double values[] = {kappa[] != nodata ? kappa[] : kappa[0, -1],
                           kappa[-1] != nodata ? kappa[-1] : kappa[-1, -1]};
        double ki = do_average (values, weights, 2);
        if (ki != nodata) {
          coord circle_center;
          double hyv = circle_y (1./(ki*Delta), hy[-1] + 0.5, hy[] + 0.5,
                                 & circle_center);
          if (fabs(hyv) < 0.5 || hyv == 0.5) {
            coord tangenti = {- (hyv - circle_center.y), - circle_center.x};
            tangenti = normalize_Q (tangenti);

            double values[] = {sigma[-1], sigma[], sigma[-1, -1], sigma[0, -1]};
            double weights[] = {0.5 + hyv, 0.5 + hyv, 0.5 - hyv, 0.5 - hyv};
            double sigmai = do_average (values, weights, 4);

            S.y.x[] = sigmai*tangenti.y/Delta;
          }
          else
            S.y.x[] = nodata;

          /**
          S.x.y */
          
          double xi = circle_x (1./(ki*Delta), hy[-1] + 0.5, hy[] + 0.5, & circle_center);
          if (fabs(xi) < 0.5 || xi == 0.5) {
            coord tangenti = {circle_center.y, xi - circle_center.x};
            tangenti = normalize_Q (tangenti);

            double values[] = {sigma[0, -1], sigma[], sigma[-1, -1], sigma[-1]};
            double weights[] = {0.5 + xi, 0.5 + xi, 0.5 - xi, 0.5 - xi};
            double sigmai = do_average (values, weights, 4);

            S.x.y[] = sigmai*sign(hy[] - hy[-1])*tangenti.x/Delta;
          }
          else
            S.x.y[] = nodata;
        }
        else {
          S.y.x[] = nodata;
          S.x.y[] = nodata; 
        }
      }
      /**
      We do not apply boundary() to S.x.y and S.y.x because they are defined on
      vertex, which are already defined on the border of the domain. */
      
      /**
      ## Resulting force
      
      Finally, we add the divergence of the surface tension tensor to
      the acceleration. */
     
      face vector av = a;
      foreach_face() {
        double force = ((S.x.x[] != nodata ? S.x.x[] : 0.)
                      - (S.x.x[-1] != nodata ? S.x.x[-1] : 0.)
                      + (S.x.y[0,1] != nodata ? S.x.y[0,1] : 0.)
                      - (S.x.y[] != nodata ? S.x.y[] : 0.))/Delta;
        av.x[] += alpha.x[]/fm.x[]*force;
      }
    }
  }
}