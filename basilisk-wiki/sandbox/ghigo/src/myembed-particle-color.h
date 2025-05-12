/**
# Explicit weak fluid-solid coupling

This module is an extension of the
[myembed-moving-color.h](myembed-moving-color.h) module. */

#include "myembed-moving-color.h"

/**
We solve here the following equations:
$$
\begin{aligned}
&
\frac{\mathrm{d} \mathbf{x}_{\Gamma}}{\mathrm{d} t} =
\mathbf{u}_{\Gamma} \\
&
\frac{\mathrm{d} \mathbf{u}_{\Gamma}}{\mathrm{d} t} =
\frac{\mathbf{F}_{\Gamma}}{\rho_{\Gamma} V_{\Gamma}} + \left(1 -
\frac{\rho}{\rho_{\Gamma}}\right) \mathbf{g} \\
&
\frac{\mathrm{d}}{\mathrm{d} t} \left(\mathbf{I}_{\Gamma}
\mathbf{\omega}_{\Gamma}\right) = \mathbf{T}_{\Gamma},
\end{aligned}
$$
where $V_{\Gamma}$ and $\mathbf{I}_{\Gamma}$ are respectively the
volume and moment of inertia tensor of the rigid body $\Gamma$ and
$\mathbf{g}$ is the gravity acceleration vector. For the sake of
simplicity, we do not include here the time evolution equation for the
angular position $\mathbf{\theta}_{\Gamma}$ of the rigid body as we
consider only freely moving spherical particles in the following. The
vectors $\mathbf{F}_{\Gamma}$ and $\mathbf{T}_{\Gamma}$ respectively
represent the hydrodynamic force and torque (about the center of mass
$\mathbf{x}_{\Gamma}$) exerted by the fluid on the rigid body
$\Gamma$:
$$
\begin{aligned}
&
\mathbf{F}_{\Gamma} = -\int_{\delta \Gamma} \left(-p \mathbb{I} + 2
\mu \mathbf{D}\right)\cdot \mathbf{n}_{\Gamma} \, \mathrm{d} S \\
&
\mathbf{T}_{\Gamma} = -\int_{\delta \Gamma} \left(\mathbf{x} -
\mathbf{x}_{\Gamma}\right) \times \left(-p \mathbb{I} + 2 \mu
\mathbf{D}\right) \cdot \mathbf{n}_{\Gamma} \, \mathrm{d} S,
\end{aligned}
$$
where $\mathbf{n}_{\Gamma}$ is the inward (pointing from the fluid
towards the rigid body) unit normal vector to the rigid boundary
$\delta \Gamma$.

In practice, using this set of equations, we update here the
quantities defined in the module
[myembed-moving-color](myembed-moving-color.h): the position *p_p*,
velocities *p_u, p_w* and acceleration *p_au, p_aw* of the discrete
rigid body $\Gamma_{\Delta}$.

## Setup

The particles density *p_r*, volume *p_v*, moment of inertia *p_i* and
the gravity field *p_g* are defined by the user. */

extern const double p_r; // Particle's density
extern const double p_v; // Particle's volume
extern const coord  p_i; // Particle's moment of inertial
extern const coord  p_g; // Particle's gravity field

/**
## Help functions */

#define p_volume_cylinder(d) (pi*sq ((d)/2.))
#define p_moment_inertia_cylinder(d,r) (1./2.*((r)*(p_volume_cylinder ((d))))*sq ((d)/2.))

#define p_volume_sphere(d) (4./3.*pi*cube ((d)/2.))
#define p_moment_inertia_sphere(d,r) (2./5.*((r)*(p_volume_sphere ((d))))*sq ((d)/2.))

/**
## Prediction 

We compute the motion of the discrete rigid body $\Gamma_{\Delta}$,
from time $t^{n}$ to time $t^{n+1}$ using the following first-order
explicit time discretization:
$$
\begin{aligned}
&
\frac{\mathbf{u}_{\Gamma}^{n+1} - \mathbf{u}_{\Gamma}^n}{\Delta t} =
\frac{\mathbf{F}_{\Gamma}^{n}}{\rho_{\Gamma}V_{\Gamma}} + \left(1 -
\frac{\rho}{\rho_{\Gamma}}\right)\mathbf{g} \\
& 
\frac{\mathbf{I}_{\Gamma}^{n+1}\mathbf{\omega}_{\Gamma}^{n+1} -
\mathbf{I}_{\Gamma}^{n}\mathbf{\omega}_{\Gamma}^{n}}{\Delta t} =
\mathbf{T}_{\Gamma}^{n} \\
&
\frac{\mathbf{x}_{\Gamma}^{n+1} - \mathbf{x}_{\Gamma}^n}{\Delta t} =
\frac{\mathbf{u}_{\Gamma}^{n} + \mathbf{u}_{\Gamma}^{n+1}}{2} .
\end{aligned}
$$

In each cut-cell, we compute the pressure contribution to the force
and torque by linearly interpolating the pressure $p^{n}$ from the
center of the cell to the centroid $\mathbf{b}^{n}$ of the discrete
rigid boundary $\delta \Gamma_{\Delta}^{n}$ in the cut-cell.

We then compute the viscous contribution to the force and torque,
assuming that the velocity $\mathbf{u}^{n}$ is constant along the
discrete rigid boundary $\delta \Gamma_{\Delta}^{n}$,
i.e. $\mathbf{\nabla} \mathbf{u} \rvert_{\delta \Gamma_{\Delta}^{n}}
\cdot \mathbf{\bar{t}}_{\Gamma}^{n} = 0$, where
$\mathbf{\bar{t}}_{\Gamma}^{n}$ is the tangential vector to the
discrete rigid boundary in the cut-cell.
*/

event advection_term (i++)
{
  /**
  We first compute the forces and torques. Note that we can only
  perform 1 force evaluation as the pressure is computed only at the
  projection step. 

  We start by redefining the color field to avoid interpolation errors
  on the geometry. */
  
  p_shape_col (p_col, (p_p));
  boundary ({p_col});

  coord Fp, Fmu, Tp, Tmu;
  embed_color_force (p, u, mu, p_col, &Fp, &Fmu);
  embed_color_torque (p, u, mu, p_col, (p_p), &Tp, &Tmu);

  /**
  We compute the particle's accelerations. */

  foreach_dimension() {
    p_au.x = (Fp.x + Fmu.x)/(p_r*p_v) + (1. - 1./(p_r))*p_g.x;
    p_aw.x = (Tp.x + Tmu.x)/(p_i.x);
  }

  /**
  We compute the particle's position. */

  foreach_dimension()
    p_p.x += (dt)*(p_u.x) + sq (dt)/2.*p_au.x;

  /**
  We compute the particle's velocities. */

  foreach_dimension() {
    p_u.x += (dt)*p_au.x;
    p_w.x += (dt)*p_aw.x;
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
