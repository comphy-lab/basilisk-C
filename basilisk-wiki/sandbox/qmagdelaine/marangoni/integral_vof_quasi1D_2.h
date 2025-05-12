/**
# Integral formulation for surface tension

This is a beginning of implementation of the integral formulation for surface
tension depicted in
[Abu-Al-Saud, Popinet, Tchelepi, 2018](https://hal.archives-ouvertes.fr/hal-01706565/).
We adapt [integral.h](/sandbox/popinet/integral.h) (written by Stéphane in
Basilisk C but for a level set description of the interface) for a VOF
description, quasi 1D interface, with one set of height function. */

#include "curvature.h"
#include "contact.h"

/**
We associate to each VOF tracer a (potentially variable) surface tension
$\sigma$ and a contact angle $\theta_\mathrm{c}$ on the different boundary (for
the moment only left and right are available). This is done easily by adding
the following [field attributes](/Basilisk C#field-attributes). */

attribute {
  scalar sigma;
  double theta_c;
  bool left;
  bool right;
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
	      foreach() {
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

//  int len = (sizeof values)/(sizeof (double));

  double weight_sum = 0.;
  double avg = 0.;
  for (int i = 0; i < length; i++) {
    weight_sum += (values[i] != nodata ? weights[i] : 0.);
    avg += (values[i] != nodata ? weights[i]*values[i] : 0.);
  }
  return (weight_sum != 0. ? avg/weight_sum : 0.);
}

coord normalize_Q (coord v) {
  double nv = 0.;
  foreach_dimension()
    nv += sq(v.x);
  nv = sqrt(nv);
  foreach_dimension()
    v.x = (nv != 0. ? v.x/nv : 0.);
  return v;
}

coord tangent_left (Point point, scalar h) {
  coord T;
  if (h[] != nodata && h[-1] != nodata) {
    T.x = 1.;
    T.y = h[] - h[-1];
    T = normalize_Q(T);
  }
  else {
      T.x = nodata;
      T.y = nodata;
  }
  return T;
}

coord tangent_center (Point point, scalar h) {
  coord T0 = tangent_left (point, h);
  coord T1 = tangent_left (neighborp(1), h);
  coord T;
  foreach_dimension() {
    double values[] = {T0.x, T1.x}, weights[] = {1., 1.};
    T.x = do_average (values, weights, 2);
  }
  if (T.x != nodata)
    T = normalize_Q(T);
  return T;
}

#define LINEAR_CURVATURE 0 // set to 1 to use linear interpolation of curvature

/**
# Surface tension term

The calculation of the acceleration is done by this event, overloaded
from [its definition](navier-stokes/centered.h#acceleration-term) in
the centered Navier--Stokes solver. */

event acceleration (i++) {
  
  /**
  We check whether surface tension is associated with any VOF tracer $c$. */

  for (scalar c in interfaces) {
    if (c.sigma.i) {

      scalar dhy;
      dhy = c.height.y;

#if 0
      /**
      Boundary conditions on heigth functions. */

      dhy[left] = neumann(1./tan(f.theta_c)/Delta);
      dhy[right] = neumann(1./tan(f.theta_c)/Delta);
#endif

      foreach()
        c[] = clamp(c[], 0., 1.);
      boundary({c});

      /**
      We compute the curvature. */
      
      scalar kappa[];
      curvature (c, kappa);
      //height_y (c, dhy);
      boundary({dhy});
            
      /**
      We rename c.sigma sigma, for convenience. */
      
      scalar sigma = c.sigma;
      boundary({sigma});

      /**
      ## Diagonal components of the stress tensor
      
      We compute the diagonal components of the tensor. */

      tensor S[];
      for (scalar s in {S}) {
        foreach_dimension() {
          s[left] = neumann(0.);
          s[right] = neumann(0.);
        }
      }

      foreach() {      
        foreach_dimension()
          S.x.x[] = 0.;
        
        /**
        S.x.x
        
        We check if the interface cross the vertical middle of the cell with 
        fabs(dhy[]) < 0.5. To take into account the case dhy[] = 0.5 once we
        have to break the symetry : strict on one side and inclusive on the
        other. Since sign(x) = (x > 0 ? 1 : -1), the special dhy == 0 will be
        handled as negative, then the logical choice is
        (dhy[] < 0.5 && dhy[] >= -0.5). */
        
        if (dhy[] < 0.5 && dhy[] >= -0.5) {
          double xi = fabs(dhy[]);
          int l = (dhy[] > 0. ? 1 : -1);
          double weights[] = {1. - xi, xi};
#if LINEAR_CURVATURE // does not make much difference
          double values[] = {kappa[], kappa[0, l]};
          double ki = do_average (values, weights, 2);
#else
          double ki = (kappa[] != nodata ? kappa[] :
                       kappa[0, l] != nodata ? kappa[0, l] : 0.);
#endif
          coord tangenti = tangent_center (point, dhy);
          double sigma_values[] = {sigma[], sigma[0, l]};
          double sigmai = do_average (sigma_values, weights, 2);
          S.x.x[] += sigmai*(fabs(tangenti.x)/Delta + sign(dhy[])*ki*(0.5 - xi));
        }

        /**
        S.y.y
        
        We check if the interface cross the horizontal middle of the cell with
        (dhy[]*(dhy[] + dhy[l]) < 0.). To take into account the special cases
        (dhy[+/-1/2] = 0 and dhy = 0), we include the cell if the interface is
        on the left, and exclude it if it is on the right. */
        
        for (int l = -1; l <= 1; l += 2) {
          if (dhy[]*(dhy[] + dhy[l]) <= 0. && (l < 0. || dhy[] + dhy[l] != 0.) &&
              dhy[] != nodata && dhy[l] != nodata && dhy[] != dhy[l]) {
            double xi = dhy[]/(dhy[] - dhy[l]); // > 0
            double weights[] = {1. - xi, xi};
#if LINEAR_CURVATURE // does not make much difference
            double kappa_values[] = {kappa[], kappa[l]};
            double ki = do_average (kappa_values, weights, 2);
#else
            double ki = (kappa[] != nodata ? kappa[] :
                       kappa[l] != nodata ?  kappa[l] : 0.);
#endif
            coord principal_tangent = tangent_left(neighborp((1 + l)/2), dhy);
            coord other_tangent = tangent_left(neighborp((1 - l)/2), dhy);
            coord tangenti;
            double tangent_weights[] = {1. + 2.*xi, 1. - 2.*xi};
            foreach_dimension() {
              double tangent_values[] = {principal_tangent.x, other_tangent.x};
              tangenti.x = do_average (tangent_values, tangent_weights, 2);
            }
            tangenti = normalize_Q(tangenti);
            double sigma_values[] = {sigma[], sigma[l]};
            double sigmai = do_average (sigma_values, weights, 2);
            S.y.y[] += sigmai*(fabs(tangenti.y)/Delta + sign(dhy[])*ki*(0.5 - xi));
          }
        }
      }
 
     boundary ({S.x.x, S.y.y});

      /**
      ## Off-diagonal components of the stress tensor
      
      We compute the off-diagonal components of the tensor.  */

      foreach_vertex() {
        
        /**
        S.y.x

        We test if the interface cross the face i = -0.5, j in [-0.5, 0.5]
        thanks to (fabs(dhyv) < 0.5) with dhyv = (dhy[] + dhy[-1, 0])/2. + 0.5,
        the hight from the vertex.
        To take into account the case where dhyv = 0.5, we rather test 
        (dhyv <= 0.5 && dhyv > -0.5). To not have allocate dhyv if it is not
        useful, we use its expression:
        (dhy[] + dhy[-1, 0] <= 0. && dhy[] + dhy[-1, 0] > -2.). */
        
        if (dhy[] + dhy[-1, 0] <= 0. && dhy[] + dhy[-1, 0] > -2.) {
          double dhyv = (dhy[] + dhy[-1, 0])/2. + 0.5;
                         
          double values[] = {sigma[-1], sigma[], sigma[-1, -1], sigma[0, -1]};
          double weights[] = {(0.5 + dhyv)/2., (0.5 + dhyv)/2.,
                              (0.5 - dhyv)/2., (0.5 - dhyv)/2.};
          double sigmai = do_average (values, weights, 4);
          
          coord tangenti = tangent_left(point, dhy);
          S.y.x[] = sigmai*tangenti.y/Delta;
        }
        else
          S.y.x[] = 0.;
          
        /**
        S.x.y */
        
        if (((dhy[] + 0.5)*(dhy[-1] + 0.5) < 0. || dhy[] == -0.5)
             && dhy[] != nodata && dhy[-1] != nodata && dhy[] != dhy[-1]) {
          double dhyv = (dhy[] + dhy[-1])/2. + 0.5;
          double xi = - dhyv/(dhy[] - dhy[-1]);
          
          double values[] = {sigma[0, -1], sigma[], sigma[-1, -1], sigma[-1]};
          double weights[] = {(0.5 + xi)/2.,(0.5 + xi)/2.,
                              (0.5 - xi)/2., (0.5 - xi)/2.};
          double sigmai = do_average (values, weights, 4);
          
          int l = (xi > 0. ? 1 : -1);

          coord principal_tangent = tangent_left (point, dhy);
          coord other_tangent = tangent_left (neighborp(l), dhy);
          coord tangenti;
          double tangent_weights[] = {1. - fabs(xi), fabs(xi)};
          foreach_dimension() {
            double tangent_values[] = {principal_tangent.x, other_tangent.x};
            tangenti.x = do_average (tangent_values, tangent_weights, 2);
          }
          tangenti = normalize_Q (tangenti);
          S.x.y[] = sigmai*sign(dhy[] - dhy[-1])*tangenti.x/Delta;
        }
        else
          S.x.y[] = 0.;
      }
            
      /**
      ## Boundary conditions
      
      We do not apply boundary() to S.x.y and S.y.x because they are defined on
      vertex, which are already defined on the border of the domain. */
      
      /**
      If we want to apply a contact angle condition, we have to take care of the
      boundary on S.x.x and S.y.y. They cannot be handled simply with neumann()
      or dirichlet(), because, the cell where the stress is not null is not
      necessarly the ghost cell. */
      
      /**
      S.x.x on the left */
      
      if (c.left) {
        foreach_boundary(left) {
      
          S.x.x[ghost] = 0.;
          
          /**
          We check if the interface cross the vertical middle of the ghost cell with 
          fabs(dhy[]) < 0.5. To take into account the case dhy[] = 0.5 once we
          have to break the symetry : strict on one side and inclusive on the
          other. Since sign(x) = (x > 0 ? 1 : -1), the special dhy == 0 will be
          handled as negative, then the logical choice is 
          (dhy[] < 0.5 && dhy[] >= -0.5). */
          
          if (dhy[ghost] < 0.5 && dhy[ghost] >= -0.5) {
            int ghost_Q = -1;
            double xi = fabs(dhy[ghost]);
            int l = (dhy[ghost] > 0. ? 1 : -1);
            double weights[] = {1. - xi, xi};
#if LINEAR_CURVATURE // does not make much difference
            double values[] = {kappa[ghost], kappa[ghost_Q, l]};
            double kg = do_average (values, weights, 2);
#else
            double kg = (kappa[ghost] != nodata ? kappa[ghost] :
                         kappa[ghost_Q, l] != nodata ? kappa[ghost_Q, l] : 0.);
#endif
            /**
            We impose a dirichlet condition on the tangent:
            
            * "b" for on the boundary,
            * "i" for in the local cell,
            * "g" for in the gosth cell.
            
            $\mathbf{t}_\mathrm{b} = (\mathbf{t}_\mathrm{g} + \mathbf{t}_i)/2$
            */
            
            coord tangentb = {sin(c.theta_c), - cos(c.theta_c)};
            coord tangenti = tangent_center (point, dhy);
            coord tangentg;
            foreach_dimension()
              tangentg.x = 2.*tangentb.x - tangenti.x;
            tangentg = normalize_Q (tangentg);

            double sigma_values[] = {sigma[ghost], sigma[ghost_Q, l]};
            double sigmag = do_average (sigma_values, weights, 2);
            S.x.x[ghost] += (sigmag*(fabs(tangentg.x)/Delta
                             + sign(dhy[ghost])*kg*(0.5 - xi)));
          }
        }
      }
      
      /**
      S.x.x on the right
      The only difference is in the tangent on the boundary:
      
      * on the left it is $\{\sin(\theta), - \cos(\theta)\}$
      
      * on the right it is $\{\sin(\theta), \cos(\theta)\}$ */
       
      if (c.right) {
        foreach_boundary(right) {
          S.x.x[ghost] = 0.;
          
          /**
          We check if the interface cross the vertical middle of the ghost cell with 
          fabs(dhy[]) < 0.5. To take into account the case dhy[] = 0.5 once we
          have to break the symetry : strict on one side and inclusive on the
          other. Since sign(x) = (x > 0 ? 1 : -1), the special dhy == 0 will be
          handled as negative, then the logical choice is 
          (dhy[] < 0.5 && dhy[] >= -0.5). */
          
          if (dhy[ghost] < 0.5 && dhy[ghost] >= -0.5) {
            int ghost_Q = 1;
            double xi = fabs(dhy[ghost]);
            int l = (dhy[ghost] > 0. ? 1 : -1);
            double weights[] = {1. - xi, xi};
#if LINEAR_CURVATURE // does not make much difference
            double values[] = {kappa[ghost], kappa[ghost_Q, l]};
            double kg = do_average (values, weights, 2);
#else
            double kg = (kappa[ghost] != nodata ? kappa[ghost] :
                         kappa[ghost_Q, l] != nodata ? kappa[ghost_Q, l] : 0.);
#endif
            /**
            We impose a dirichlet condition on the tangent. */
            
            coord tangentb = {sin(c.theta_c), cos(c.theta_c)};
            coord tangenti = tangent_center (point, dhy);
            coord tangentg;
            foreach_dimension()
              tangentg.x = 2.*tangentb.x - tangenti.x;
            tangentg = normalize_Q (tangentg);

            double sigma_values[] = {sigma[ghost], sigma[ghost_Q, l]};
            double sigmag = do_average (sigma_values, weights, 2);
            S.x.x[ghost] += (sigmag*(fabs(tangentg.x)/Delta
                             + sign(dhy[ghost])*kg*(0.5 - xi)));
          }
        }
      }

      /**
      ## Resulting force
      
      Finally, we add the divergence of the surface tension tensor to
      the acceleration. */
      
      face vector av = a;
      foreach_face()
        av.x[] += (S.x.x[] - S.x.x[-1] + S.x.y[0,1] - S.x.y[])/Delta;
    }
  }
}
