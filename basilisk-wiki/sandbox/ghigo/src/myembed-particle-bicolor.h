/**
# Explicit weak fluid-solid coupling

This module is an extension of the
[myembed-moving-bicolor.h](myembed-moving-bicolor.h) module. */

#include "myembed-moving-bicolor.h"

/**
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
quantities defined in the module [myembed-moving-bicolor.h](): the
position *p_p[i]*, velocities *p_u[i], p_w[i]* and accelerations
*p_au[i], p_aw[i]* of the discrete rigid body $\Gamma_{i,\Delta}$.

## Setup

The particles density *p_r[i]*, volume *p_v[i]*, moment of inertia
*p_i[i]* and the gravity field *p_g* are defined by the user. */

extern const double p_r[p_n]; // Particle's density
extern const double p_v[p_n]; // Particle's volume
extern const coord  p_i[p_n]; // Particle's moment of inertial
extern const coord  p_g; // Particle's gravity field

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
  
  p0_shape_col (p0_col);
  p1_shape_col (p1_col);
  restriction ({p0_col, p1_col}); // Since some BC might depend on p_col

  /**
  We then compute the forces acting on each particle. */

  coord Fp[p_n], Fmu[p_n], Tp[p_n], Tmu[p_n];
  // Particle 1
  embed_color_force  (p, u, mu, p0_col, &Fp[0], &Fmu[0]);
  embed_color_torque (p, u, mu, p0_col, p_p[0], &Tp[0], &Tmu[0]);
  // Particle 2
  embed_color_force  (p, u, mu, p1_col, &Fp[1], &Fmu[1]);
  embed_color_torque (p, u, mu, p1_col, p_p[1], &Tp[1], &Tmu[1]);

  for (int i = 0; i < (p_n); i++) {

    /**
    We compute each particle's accelerations. */

    foreach_dimension() {
      p_au[i].x = (Fp[i].x + Fmu[i].x)/(p_r[i]*p_v[i]) + (1. - 1./(p_r[i]))*p_g.x;
      p_aw[i].x = (Tp[i].x + Tmu[i].x)/(p_i[i].x);
    }

    /**
    We compute each particle's position. */

    foreach_dimension()
      p_p[i].x += (dt)*(p_u[i].x) + sq (dt)/2.*p_au[i].x;

    /**
    We compute each particle's velocities. */

    foreach_dimension() {
      p_u[i].x += (dt)*p_au[i].x;
      p_w[i].x += (dt)*p_aw[i].x;
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
