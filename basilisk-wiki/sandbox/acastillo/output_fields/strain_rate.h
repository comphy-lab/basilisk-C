/**
Given a velocity field $\mathbf{u}$, this function returns the squared
magnitude of the strain-rate tensor at a point,
$$
S^2 = S_{ij}S_{ij}, \qquad
S_{ij} = \frac{1}{2}\left(\partial_i u_j + \partial_j u_i\right)
$$
evaluated with centered differences on the cell-centered velocity. It is
mainly used for dissipation diagnostics, where the viscous dissipation rate
of an incompressible flow is $\varepsilon = 2\mu S_{ij}S_{ij}$.

Note that the metric factors `fm` and `cm` are *not* taken into account, so
this is valid for Cartesian coordinates only, unlike
[vorticity3d()](vorticity3d.h). The trace is not removed either, which is
consistent with $\nabla\cdot\mathbf{u} = 0$. */

#ifndef STRAIN_RATE_H
#define STRAIN_RATE_H

static inline double strain_rate_sq (Point point, vector u) {
  double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
  double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
  double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
  double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
  double Sxx = dudx;
  double Sxy = 0.5*(dudy + dvdx);
  double Syx = Sxy;
  double Syy = dvdy;
  double S2 = sq(Sxx) + sq(Sxy) + sq(Syx) + sq(Syy);
  #if dimension == 3
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    double Szz = dwdz;
    double Sxz = 0.5*(dwdx + dudz);
    double Syz = 0.5*(dwdy + dvdz);
    double Szx = Sxz;
    double Szy = Syz;
    S2 += sq(Szz) + sq(Sxz) + sq(Syz) + sq(Szx) + sq(Szy);
  #endif
  return S2;
}

#endif
