/**
Given a velocity field $\mathbf{u}$, this function returns the squared
magnitudes of both the strain-rate and the rotation-rate tensors at a point,
$$
S^2 = S_{ij}S_{ij}, \qquad \Omega^2 = \Omega_{ij}\Omega_{ij}
$$
with
$$
S_{ij} = \frac{1}{2}\left(\partial_i u_j + \partial_j u_i\right), \qquad
\Omega_{ij} = \frac{1}{2}\left(\partial_i u_j - \partial_j u_i\right)
$$
from one set of centred differences on the cell-centred velocity, instead of
differencing the velocity twice.

$\Omega_{ij}$ is antisymmetric ($\Omega_{xx} = \Omega_{yy} = \Omega_{zz} = 0$),
so only the off-diagonal terms contribute:
$\Omega_{ij}\Omega_{ij} = 2(\Omega_{xy}^2 + \Omega_{xz}^2 + \Omega_{yz}^2)$
in 3D, $2\Omega_{xy}^2$ in 2D. Equivalently
$\Omega_{ij}\Omega_{ij} = \frac{1}{2}|\boldsymbol\omega|^2$.

Together they give the Q-criterion for vortex identification,
$Q = \frac{1}{2}(\Omega_{ij}\Omega_{ij} - S_{ij}S_{ij})$.

The $S^2$ returned here is the same quantity as
[strain_rate_sq()](strain_rate.h)'s, which remains the right call when
$\Omega^2$ is not wanted. The two agree to roundoff rather than bit-for-bit:
the off-diagonal terms are summed as $2S_{xy}^2$ here and as
$S_{xy}^2 + S_{yx}^2$ there, which rounds differently. Same caveats as that
function: the metric factors
`fm` and `cm` are *not* taken into account, so this is valid for Cartesian
coordinates only, unlike [vorticity3d()](vorticity3d.h). Neither tensor is
trace-removed, which is consistent with $\nabla\cdot\mathbf{u} = 0$. */

#ifndef STRAIN_AND_VORTICITY_H
#define STRAIN_AND_VORTICITY_H

static inline void strain_and_vorticity_sq (Point point, vector u,
                                             double * S2, double * O2) {
  double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
  double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
  double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
  double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
  double Sxx = dudx;
  double Syy = dvdy;
  double Sxy = 0.5*(dudy + dvdx);
  double Oxy = 0.5*(dvdx - dudy);
  *S2 = sq(Sxx) + sq(Syy) + 2.*sq(Sxy);
  *O2 = 2.*sq(Oxy);
  #if dimension == 3
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double Szz = dwdz;
    double Sxz = 0.5*(dwdx + dudz);
    double Syz = 0.5*(dwdy + dvdz);
    double Oxz = 0.5*(dwdx - dudz);
    double Oyz = 0.5*(dwdy - dvdz);
    *S2 += sq(Szz) + 2.*sq(Sxz) + 2.*sq(Syz);
    *O2 += 2.*sq(Oxz) + 2.*sq(Oyz);
  #endif
}

#endif
