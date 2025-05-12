/**
# Explicit weak fluid-solid coupling

This module is an extension of the
[myembed-moving-multicolor.h](myembed-moving-multicolor.h) module.

We solve here the following equations for each particle:
$$
\begin{aligned}
&
\frac{\mathrm{d} \mathbf{x}_{\Gamma,i}}{\mathrm{d} t} =
\mathbf{u}_{\Gamma} \\
&
\frac{\mathrm{d} \mathbf{u}_{\Gamma,i}}{\mathrm{d} t} =
\frac{\mathbf{F}_{\Gamma,i}}{\rho_{\Gamma} V_{\Gamma,i}} + \left(1 -
\frac{\rho}{\rho_{\Gamma,i}}\right) \mathbf{g} \\
&
\frac{\mathrm{d}}{\mathrm{d} t} \left(\mathbf{I}_{\Gamma,i}
\mathbf{\omega}_{\Gamma,i}\right) = \mathbf{T}_{\Gamma,i},
\end{aligned}
$$
where $V_{\Gamma,i}$ and $\mathbf{I}_{\Gamma,i}$ are respectively the
volume and moment of inertia tensor of the rigid body $\Gamma_{i}$ and
$\mathbf{g}$ is the gravity acceleration vector. For the sake of
simplicity, we do not include here the time evolution equation for the
angular position $\mathbf{\theta}_{\Gamma,i}$ of the rigid body as we
consider only freely moving spherical particles in the following. The
vectors $\mathbf{F}_{\Gamma,i}$ and $\mathbf{T}_{\Gamma,i}$
respectively represent the hydrodynamic force and torque (about the
center of mass $\mathbf{x}_{\Gamma,i}$) exerted by the fluid on the
rigid body $\Gamma_{i}$:
$$
\begin{aligned}
&
\mathbf{F}_{\Gamma,i} = -\int_{\delta \Gamma_{i}} \left(-p \mathbb{I}
+ 2 \mu \mathbf{D}\right)\cdot \mathbf{n}_{\Gamma,i} \, \mathrm{d} S
\\
&
\mathbf{T}_{\Gamma,i} = -\int_{\delta \Gamma_{i}} \left(\mathbf{x} -
\mathbf{x}_{\Gamma,i}\right) \times \left(-p \mathbb{I} + 2 \mu
\mathbf{D}\right) \cdot \mathbf{n}_{\Gamma,i} \, \mathrm{d} S,
\end{aligned}
$$
where $\mathbf{n}_{\Gamma,i}$ is the inward (pointing from the fluid
towards the rigid body) unit normal vector to the rigid boundary
$\delta \Gamma_{i}$.

In practice, using this set of equations, we update here the
quantities defined in the module [myembed-moving-multicolor.h](): the
position *c*, velocities *u, w* and accelerations *au, aw* of the
discrete rigid body $\Gamma_{i,\Delta}$.

## Setup

The particles' density *r*, volume *v*, moment of inertia *i* are
defined by the user. */

#ifndef ADD_PARTICLE_PHY_PARAMS
#define ADD_PARTICLE_PHY_PARAMS double r; double v; coord i;
#endif // ADD_PARTICLE_PHY_PARAMS

#include "myembed-moving-multicolor.h"

/**
The gravity field *p_g* is also defined by the user. */

extern const coord  p_g; // Particles' gravity field

/**
## Help functions */

#define p_volume_cylinder(d) (pi*sq ((d)/2.))
#define p_moment_inertia_cylinder(d,r) (1./2.*((r)*(p_volume_cylinder ((d))))*sq ((d)/2.))

#define p_volume_sphere(d) (4./3.*pi*cube ((d)/2.))
#define p_moment_inertia_sphere(d,r) (2./5.*((r)*(p_volume_sphere ((d))))*sq ((d)/2.))

/**
## Prediction 

We compute the motion of the discrete rigid body $\Gamma_{i,\Delta}$,
from time $t^{n}$ to time $t^{n+1}$ using the following first-order
explicit time discretization:
$$
\begin{aligned}
&
\frac{\mathbf{u}_{\Gamma,i}^{n+1} - \mathbf{u}_{\Gamma,i}^n}{\Delta t}
= \frac{\mathbf{F}_{\Gamma,i}^{n}}{\rho_{\Gamma,i}V_{\Gamma,i}} +
\left(1 - \frac{\rho}{\rho_{\Gamma,i}}\right)\mathbf{g} \\
& 
\frac{\mathbf{I}_{\Gamma,i}^{n+1}\mathbf{\omega}_{\Gamma,i}^{n+1} -
\mathbf{I}_{\Gamma,i}^{n}\mathbf{\omega}_{\Gamma,i}^{n}}{\Delta t} =
\mathbf{T}_{\Gamma,i}^{n} \\
&
\frac{\mathbf{x}_{\Gamma,i}^{n+1} - \mathbf{x}_{\Gamma,i}^n}{\Delta t}
= \frac{\mathbf{u}_{\Gamma,i}^{n} + \mathbf{u}_{\Gamma,i}^{n+1}}{2} .
\end{aligned}
$$

In each cut-cell, we compute the pressure contribution to the force
and torque by linearly interpolating the pressure $p^{n}$ from the
center of the cell to the centroid $\mathbf{b}_{i}^{n}$ of the
discrete rigid boundary $\delta \Gamma_{i,\Delta}^{n}$ in the
cut-cell.

We then compute the viscous contribution to the force and torque,
assuming that the velocity $\mathbf{u}^{n}$ is constant along the
discrete rigid boundary $\delta \Gamma_{i,\Delta}^{n}$,
i.e. $\mathbf{\nabla} \mathbf{u} \rvert_{\delta \Gamma_{i,\Delta}^{n}}
\cdot \mathbf{\bar{t}}_{\Gamma,i}^{n} = 0$, where
$\mathbf{\bar{t}}_{\Gamma,i}^{n}$ is the tangential vector to the
discrete rigid boundary in the cut-cell.
*/

event advection_term (i++)
{
  /**
  We first compute the forces and torques on each particle. Note that
  we can only perform 1 force evaluation as the pressure is computed
  only at the projection step.

  We start by redefining the color fields to avoid interpolation
  errors on the geometry. */
  
  p_shape_col (pl);

  for (int i = 0; i < (p_n); i++) {

    /**
    If the particle has no density, it is considered as a wall. */
    
    if (pl[i].r) {

      /**
      We then compute the forces acting on the particle. */
      
      coord Fp, Fmu, Tp, Tmu;
      embed_color_force  (p, u, mu, pl[i].col, &Fp, &Fmu);
      embed_color_torque (p, u, mu, pl[i].col, pl[i].c, &Tp, &Tmu);

      /**
      We compute each particle's accelerations. */

      foreach_dimension() {
	pl[i].au.x = (Fp.x + Fmu.x)/(pl[i].r*pl[i].v) + (1. - 1./(pl[i].r))*p_g.x;
	pl[i].aw.x = (Tp.x + Tmu.x)/(pl[i].i.x);
      }

      /**
      We compute each particle's position. */

      foreach_dimension()
	pl[i].c.x += (dt)*(pl[i].u.x) + sq (dt)/2.*pl[i].au.x;

      /**
      We compute each particle's velocities. */

      foreach_dimension() {
	pl[i].u.x += (dt)*pl[i].au.x;
	pl[i].w.x += (dt)*pl[i].aw.x;
      }
    }
  }
}

/**
## References

~~~bib
@article{schneiders2016,
  title={An efficient conservative cut-cell method for rigid bodies interacting with viscous compressible flows},
  author={Schneiders, L. and Gunther, C. and Meinke, M. and Schroder, W.},
  journal={Journal of Comp_utational Physics},
  volume={311},
  pages={62--86},
  year={2016},
  publisher={Elsevier}
}
~~~
*/
