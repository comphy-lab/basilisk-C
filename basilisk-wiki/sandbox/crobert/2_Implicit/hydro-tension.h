/**
# Multilayer solver with surface tension

This file adds surface tension to the [multilayer solver](/src/layered/README)
 i.e. the Laplace pressure and the Marangoni effect that intervenes in the 
 continuity of the surface stress continuity. They are respectively  added 
 as a pressure term and a top boundary condition in viscosity.

$$
(\mathbf{T_{liquid}} - \mathbf{T_{gas}}) \cdot \mathbf{n} 
  = - \sigma \kappa\mathbf{n} + \mathbf{\nabla_S} \sigma
$$
with $T$, the stress tensors respectively in the gas and in the liquid, $\sigma$
 the surface tension coefficient, $\kappa$ the curvature of the free-surface
 and $\mathbf{\nabla_S}$ a surface gradient.

Note that this file is also compatible with the [implicit free-surface
extension](/src/layered/implicit.h) and with the [non-hydrostatic
extension](/src/layered/nh.h). In both cases the Laplace pressure term
is treated implicitly and does not restrict the timestep. 

The default surface tension coefficient $\sigma$ is constant and equal to
unity. */

(const) scalar sigma = unity;

/**
## Laplace pressure
The Laplace pressure corresponds to a barotropic pressure 
 $\phi(x) = - \sigma/\rho \kappa$ that is added to the hydrostatic
 equations (term in blue).
$$
\begin{aligned}
  \partial_t h_k + \mathbf{{\nabla}} \cdot \left( h \mathbf{u} \right)_k & =
  0,\\
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta)  
  \color{blue} + h_k \mathbf{{\nabla}} (\sigma \kappa)
\end{aligned}
$$

### Implementation

We will need auxilliary fields containing the "non-linear surface
tension coefficients", which are
$$
\sigma_n = \frac{\sigma}{n}
$$
in one dimension and
$$
\begin{aligned}
\sigma_x & = \sigma \frac{1 + (\partial_y\eta)^2}{n} \\
\sigma_y & = \sigma \frac{1 + (\partial_x\eta)^2}{n} \\
\sigma_n & = - 2 \sigma \frac{\partial_x\eta \partial_y\eta}{n}
\end{aligned}
$$
in two dimensions, with
$$
n = (1 + |\mathbf{\nabla}\eta|^2)^{3/2}
$$
*/

scalar sigma_n;
#if dimension > 1
vector sigma_d;
#endif

/**
We overload the default definition of the barotropic acceleration in
[hydro.h](/src/layered/hydro.h) which becomes
$$
\mathbf{a}_\text{baro} = - g \mathbf{{\nabla}} (\eta) 
                           \color{blue} + \mathbf{{\nabla}} (\sigma \kappa)
$$
where $\sigma \kappa$ is computed as
$$
\sigma \kappa = \sigma_n \frac{\partial^2\eta}{\partial x^2}
$$
in two dimensions and
$$
\sigma \kappa = \sigma_x \frac{\partial^2\eta}{\partial x^2} +
  \sigma_y \frac{\partial^2\eta}{\partial y^2} +
  \sigma_n \frac{\partial^2\eta}{\partial x y}
$$
in three dimensions.

In 3D, the weighing with the `0.2` coefficient is necessary to avoid odd-even
decoupling. */

#if dimension == 1
# define sigma_kappa(eta, i)					\
  (sigma_n[i]*(eta[i+1] + eta[i-1] - 2.*eta[i])/sq(Delta))
#else // dimension == 2
# define sigma_kappa(eta, i)						\
  ((sigma_d.x[i]*(0.2*(eta[i+1,1] + eta[i-1,1] - 2.*eta[i,1] +		\
		       eta[i+1,-1] + eta[i-1,-1] - 2.*eta[i,-1]) +	\
		  eta[i+1] + eta[i-1] - 2.*eta[i])/(1. + 2.*0.2) +	\
    sigma_d.y[i]*(0.2*(eta[i+1,1] + eta[i+1,-1] - 2.*eta[i+1] +		\
		       eta[i-1,1] + eta[i-1,-1] - 2.*eta[i-1]) +	\
		  eta[i,1] + eta[i,-1] - 2.*eta[i])/(1. + 2.*0.2) +	\
    sigma_n[i]*(eta[i+1,1] + eta[i-1,-1] -				\
		eta[i+1,-1] - eta[i-1,1]))/sq(Delta))
#endif // dimension == 2

#define p_baro(eta,i) (- G*eta[i] +	sigma_kappa(eta, i))
#define a_baro(eta, i)						\
  (gmetric(i)*(p_baro (eta, i) - p_baro (eta, i - 1))/Delta)

#include "layered/hydro.h"

/**
The non-linear surface tension coefficient(s) are defined using the
surface elevation at the beginning of the timestep. */

event face_fields (i++)
{
  sigma_n = new scalar;
#if dimension > 1
  scalar dx = new scalar, dy = new scalar;
  sigma_d.x = dx, sigma_d.y = dy;
#endif

  foreach() {
    double n = 1.;
    coord dh;
    foreach_dimension() {
      dh.x = (eta[1] - eta[-1])/(2.*Delta);
      n += sq(dh.x);
    }
    n = pow(n, 3./2.);
#if dimension == 1
    sigma_n[] = sigma[]/n;
#else
    sigma_n[] = - sigma[]*dh.x*dh.y/(2.*n);
    foreach_dimension()
      sigma_d.x[] = sigma[]*(1. + sq(dh.y))/n;
#endif
  }
#if dimension == 1
  boundary ({sigma_n});
  restriction ({sigma_n});
#else
  boundary ({sigma_n, sigma_d});
  restriction ({sigma_n, sigma_d});
#endif  
  
  /**
  In the case of a time-explicit integration (as controlled by CFL_H),
  we need to restrict the timestep based on the celerity of capillary
  waves (and not only gravity waves as done by the default solver). To
  do so, we imitate the code in
  [hydro.h](/src/layered/hydro.h#face_fields) which sets `dtmax`, but
  take into account the celerity of the shortest capillary waves (of
  wavelength $2\Delta$). */ 
  
  foreach_face (reduction (min:dtmax)) {
    double H = 0., um = 0.;
    foreach_layer() {
      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;
      if (fabs(uf) > um)
	um = fabs(uf);
      H += (h[] + h[-1])/2.;
    }
    if (H > dry) {
      double c = um/CFL + sqrt((G + max(sigma[], sigma[-1])*sq(pi/Delta))*
			       (hydrostatic ? H : Delta*tanh(H/Delta)))/CFL_H;
      double dt = min(cm[], cm[-1])*Delta/(c*fm.x[]);
      if (dt < dtmax)
	dtmax = dt;
    }
  }
}

/**
At the end of the timestep we delete the auxilliary field. */

event pressure (i++) {
  delete ({sigma_n});
#if dimension > 1
  delete ((scalar *){sigma_d});
#endif  
}

/**
# Marangoni stress 

The [Marangoni effect](https://en.wikipedia.org/wiki/Marangoni_effect)
 is a mass transfert along an interface due to surface tension gradients. 
 The Marangoni stress intervenes in the continuity of the tangential surface
 stress (term in blue):
$$ 
  \frac{\rho\,\nu}{1 + \eta_{x}^2}(1-\eta_{x}^2) (\partial_z u + \partial_x w)_{top} 
    - 4\,\eta_{x}\,\partial_x u|_{top} =
    \color{blue} \frac{\partial_x \sigma}{\sqrt{1 + \eta_{x}^2}}
$$

The Marangoni stress can be seen as a boundary condition for the vertical
 viscosity. In the thin film approximation, we get:
$$
\partial_z u|_{top} = \frac{l}{\rho\,\nu} \partial_x\sigma
$$

with $l = \frac{\sqrt{1+\eta_x^2}}{1-\eta_x^2} \approx 1$, the factor of
 non-linearity. This factor is kept in order to stay consistent with the
 complete equations added in 
[viscous_surface.h](/crobert/2_Implicit/viscous_srface.h).

### Implementation

The surface tension gradient is added as a boundary condition in the
 vertical viscosity solver, an auxilliary field $du_{Mrg}$ is needed to
 contain this condition.
*/

vector du_Mrg[];

event viscous_term (i++)
{
  if (nu > 0 && !is_constant(sigma)) {
    foreach()
      foreach_dimension(){
        double hx = (eta[1] - eta[-1])/(2.*Delta);
        double l = sqrt(1. + sq(hx))/(1. - sq(hx));
        du_Mrg.x[] = dut.x[] + (l/nu)*(sigma[1] - sigma[-1])/(2*Delta);
      }
    boundary ((scalar *){du_Mrg});
    dut = du_Mrg;
  }
}

event update_eta (i++) {
  dut = zerof;
}
