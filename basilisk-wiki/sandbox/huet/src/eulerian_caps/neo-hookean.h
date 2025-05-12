/**
# Neo-hookean law for the Eulerian elasticity framework.

In this file we implement the neo-hookean law to compute the elastic surface stresses
of biological capsules. The computation of the stresses rely on the Eulerian
quantities defined in [elasticity.h](elasticity.h). The force are transferred
to the fluid in the whole membrane region computed in [capsule.h](capsule.h).
*/

/**
We define $G_s$, the elastic modulus of the material
*/
#ifndef E_S
  #define E_S 1.
#endif

event acceleration (i++) {
  foreach() {
    if (GAMMA) {
/**
Below, we compute the stress tensor for the neo-Hookean law:
$$
\mathbf{Ts} = \frac{G_s}{3} (\mathbf{G} - \frac{1}{J^3}\mathbf{P})
$$
Note that the elastic modulus $G_s$ will be multiplied to the stress tensor
during the acceleration step in [elasticity.h](http://basilisk.fr/_edit/sandbox/huet/src/elasticity.h),
and that the acceleration will be computed using the continuum surface force
(CSF) formulation, for which we need to multiply the divergence of the stress
tensor $\nabla_s \cdot \mathbf{Ts}$ by the norm of the gradient of the smoothed color
function $|\nabla \phi|$ (in our case, the smoothed color function is defined in [capsule.h](http://basilisk.fr/_edit/sandbox/huet/src/capsule.h) and named *caps* $\equiv \phi$).
It was shown in [Ii et al.](ii2012full) that since
$|\nabla \phi|$ is a smoothed 1-dimensional Dirac distribution, we have
$|\nabla \phi| \nabla_s \cdot \mathbf{Ts} = \nabla_s \cdot (|\nabla \phi| \mathbf{Ts})$.
As a result, we can multiply the stress tensor by $|\nabla \phi|$ at this
step.
*/
      pseudo_t P;
      foreach_dimension() {
        P.x.x = 1 - sq(extended_n.x[]);
        P.x.y = -extended_n.x[]*extended_n.y[];
        P.x.z = -extended_n.x[]*extended_n.z[];
      }
      foreach_dimension() {
        Ts.x.x[] = E_S*(G.x.x[] - P.x.x/cube(J[]))/3;
      }
      Ts.x.y[] = E_S*(G.x.y[] - P.x.y/cube(J[]))/3;
      Ts.x.z[] = E_S*(G.x.z[] - P.x.z/cube(J[]))/3;
      Ts.y.z[] = E_S*(G.y.z[] - P.y.z/cube(J[]))/3;

      #if !(DIRAC_IN_ACCELERATION)
        double dirac = 0.;
        #if SHARP_DIRAC
          foreach_dimension() {
            dirac += sq((f[1] - f[-1])/(2.*Delta));
          }
          dirac = sqrt(dirac);
        #else
          dirac = ngcaps[];
        #endif
        foreach_dimension() Ts.x.x[] *= dirac;
        Ts.x.y[] *= dirac;
        Ts.x.z[] *= dirac;
        Ts.y.z[] *= dirac;
      #endif
    }
    else {
      foreach_dimension() {
        Ts.x.x[] = 0.;
      }
      Ts.x.y[] = 0.;
      Ts.x.z[] = 0.;
      Ts.y.z[] = 0.;
    }
  }
  boundary((scalar *){Ts});
}

/**
# References
~~~bib
@Article{ii2012full,
  author  = {Ii, Satoshi and Gong, Xiaobo and Sugiyama, Kazuyasu and Wu, Jinbiao and Huang, Huaxiong and Takagi, Shu},
  journal = {Communications in Computational Physics},
  title   = {A full Eulerian fluid-membrane coupling method with a smoothed volume-of-fluid approach},
  year    = {2012},
  number  = {2},
  pages   = {544},
  volume  = {12},
  file    = {:ii2012full - A Full Eulerian Fluid Membrane Coupling Method with a Smoothed Volume of Fluid Approach.pdf:PDF},
  groups  = {Biological flows},
}
~~~
*/
