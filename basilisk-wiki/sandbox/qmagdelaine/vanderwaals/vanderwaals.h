/**
# Reduced Van der Waals forces 

The implementation is due to [Mahady, Afkhami, Kondic, Phys. Fluids 28, 062002
(2016)](http://doi.org/10.1063/1.4949522). They show that Van der Waals-like
forces induced by a flat substrate on the liquid can be reduced to an
interfacial force by including them in the pressure, as it can be done for
the gravity (see [reduced.h](/src/reduced.h)):
$$
\mathbf{f} = - \xi\, F(y)\, \delta_S\, \mathbf{n} \quad \text{with} \quad
F(y) = \left[\left(\frac{h*}{y}\right)^m - \left(\frac{h*}{y}\right)^k\right]
$$
where $y$ is the distance to the substrate and $\xi =
\frac{1 - \cos(\theta_\mathrm{eq})}{h*}\,\frac{(m - 1)\, (k - 1)}{m - k}$.

This file has to be include after [two-phase.h](/src/two-phase.h). We define some
scalar attributes: the equilibrium contact angle $\theta$, the evanescent film
(used to regularise the equations) $h*$, the exponents $m$ and $k$, and the list
of strings to choose from which boundaries the forces are active. */

attribute {
  double theta, hs, m, k;
  char walls[6][10];
}

/**
We use [curvature.h](/src/curvature.h) to compure the position of the interface.
*/

#include "curvature.h"

/**
In the *defaults* event, we set all the string of *walls* to *no* in order that
the user has just to add the walls from which the force is active without taking
care of the others. */

event defaults (i = 0) {
  for (int i = 0; i < 6; i++)
    sprintf (f.walls[i], "no");
}

/**
In the *init* event, we set the VOF tracer to at the boundary from which the
forces are active because even if the liquid dewets, a evanescent film will
remain to regulatization purposes. */

event init (i = 0) {
  for (int i = 0; i < 6; i++) {
    if (!strcmp(f.walls[i], "left")) {f[left] = 1.;}
    else if (!strcmp(f.walls[i], "right")) {f[right] = 1.;}
    else if (!strcmp(f.walls[i], "bottom")) {f[bottom] = 1.;}
    else if (!strcmp(f.walls[i], "top")) {f[top] = 1.;}
#if dimension == 3
    else if (!strcmp(f.walls[i], "back")) {f[back] = 1.;}
    else if (!strcmp(f.walls[i], "front")) {f[front] = 1.;}
#endif
  }
}

/**
We do not use directly [iforce.h](/src/iforce.h), because it averages the
potential to get the face value, whereas here we rather average the distance to
get a face value and compute the face value of the potential from this
averaged value of the distance.

We overload the *acceleration()* event to add the contribution of Van der Waals
forces. */

event acceleration (i++) {
  foreach()
    f[] = clamp (f[], 0., 1.);
  boundary ({f});

  /**
  We compute the *position* of the interface with the function defined in
  [curvature.h](/src/curvature.h). This position corresponds to the distance to a
  boundary is the reference $Z$ is set correctly. */

  for (int i = 0; i < 6; i++) {
    coord direction = {0.,0.,0.}, Z = {X0,Y0,Z0};
    bool ok = true; 
    if (!strcmp(f.walls[i], "left"))
      direction.x = 1.;
    else if (!strcmp(f.walls[i], "right"))
      direction.x = 1., Z.x += L0;
    else if (!strcmp(f.walls[i], "bottom"))
      direction.y = 1.;
    else if (!strcmp(f.walls[i], "top"))
      direction.y = 1., Z.y += L0;
    else if (!strcmp(f.walls[i], "back"))
      direction.z = 1.;
    else if (!strcmp(f.walls[i], "front"))
      direction.z = 1., Z.z += L0;
    else
      ok = false;
    if (ok) {
    
      /**
      We compute $\xi$ with the attributes of the VOF tracer $f$: */
      
      double xi = (1. - cos(f.theta))/f.hs*(f.m - 1.)*(f.k - 1.)/(f.m - f.k);
      
      scalar pos[];
      position (f, pos, direction, Z, false);

      /**
      To get a well balanced scheme, all interfacial forces (surfance tension,
      reduced gravity or Van der Waals-like forces) have to be discretized
      consistently with the pressure. Therefore $\delta_S\, \mathbf{n}$ is
      evaluated with $\frac{f[] - f[-1]}{\Delta}$. To be really consistent, on 
      adaptative meshes, $f$ has to get the same prolongation function as the
      pressure. */
      
#if TREE
      f.prolongation = p.prolongation;
      boundary ({f});
#endif
      
      face vector av = a;
      foreach_face() {
        if (f[] != f[-1]) {
          double h = fabs(
                     (pos[] < nodata && pos[-1] < nodata) ?
                     (pos[] + pos[-1])/2. :
                     pos[] < nodata ? pos[] :
                     pos[-1] < nodata ? pos[-1] :
                     nodata);
          if (h != nodata)
            av.x[] -= alpha.x[]/fm.x[]*f.sigma*xi*(pow(f.hs/h, f.m)
                         - pow(f.hs/h, f.k))*(f[] - f[-1])/Delta;
        }
      }
  #if TREE
      f.prolongation = fraction_refine;
      boundary ({f});
  #endif
    }
  }
}

/**
# References

~~~bib
@article{mahady2016,
  title={A numerical approach for the direct computation of flows including
  fluid-solid interaction: modeling contact angle, film rupture, and dewetting},
  author={Mahady, Afkhami, Kondic},
  journal={Physics of Fluids},
  volume={28},
  issue={6},
  year={2016}
}
*/