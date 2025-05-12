/**
# Surface stress continuity and viscosity

The complete surface stress continuity of a free surface (without surface
 tension) writes
 $$
(\mathbf{T_{liquid}} - \mathbf{T_{gas}}) \cdot \mathbf{n} = 0
$$

Projecting this equation along normal and tangential directions gives (in 2D)
$$
\begin{aligned}
  &\phi|_{top} = - \frac{2\,\nu}{1 + \eta_{x}^2} (\partial_x u|_{top} 
      (1 - \eta_{x}^2) - \eta_{x} (\partial_z u + \partial_x w))_{top} \\
  &\frac{\rho\,\nu}{1 + \eta_{x}^2}(1-\eta_{x}^2) (\partial_z u + \partial_x w)_{top} 
    = 4\,\eta_{x}\,\partial_x u|_{top}
\end{aligned}
$$
where &\phi_{nu}$ is the non-hydrostratique pressure. 

In the limits of small slopes, we get the simplified equations
$$
\begin{aligned}
  &\phi|_{top} = - 2\,\nu\,\partial_x u|_{top} \\
  &\partial_z u|_{top} = - \partial_x w|_{top} 
\end{aligned}
$$

These terms are added respecively in the momentum equation and in the
 vertical viscosity diffusion. */

bool visc_activate = true;
 
/**
## Normal stress continuity
The pressure contribution obtained in the normal stress continuty can be added
 to the momentum equation of the [multilayer solver](hydro.h) as
$$
\begin{aligned}
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta)
  {\color{blue} - \mathbf{{\nabla}} (h \phi_{\nu})_k + \left[ \phi_{\nu} 
  \mathbf{{\nabla}} z \right]_k}
\end{aligned}
$$
where the terms in blue have been added and $\phi_{\nu}$ is the viscous
pressure deviation due to velocity variations on the surface.

These terms are added to the acceleration of the [multilayer solver](hydro.h).*/

event acceleration (i++)
{
  /**
  The hydrostatic pressure deviation $\phi_{\nu}$ is stored on the interfaces
  between layers (consistently with the Keller box scheme
  discretisation, see [Popinet,
  2020](/Bibliography#popinet2020)). This gives the following vertical
  discrete integration scheme. */
  if (visc_activate){
    scalar phiNu = new scalar[nl];
    foreach(){
      double phiNu0 = 0.;
      foreach_dimension()
        phiNu0 -= 2.*nu*(u.x[1,0,nl-1] - u.x[-1,0,nl-1])/(2*Delta);
      foreach_layer()
        phiNu[] = phiNu0;
    }
    boundary ({phiNu});

    /**
    Once the pressure deviation is known, the terms in blue above are
    added to the face acceleration field `ha`, using the
    [pressure-gradient macro](hydro.h#horizontal-pressure-gradient). */
    
    foreach_face()
      hpg (pg, phiNu, 0,
	   ha.x[] += pg);
    boundary ((scalar *){ha});
  
    delete ({phiNu});
  }
}

/**
## Tangential stress continuity 
The viscous contribution is added as a boundary condition in the
 vertical viscosity solver, an auxilliary field $du_{\nu}$ is needed to
 contain this condition.
*/

#if NH
vector du_nu[];

event viscous_term (i++) 
{
  if (visc_activate){
    foreach()
      foreach_dimension ()
        du_nu.x[] = - (w[1,0,nl-1] - w[-1,0,nl-1])/(2.*Delta);
    boundary ((scalar *){du_nu});
    dut = du_nu;
  }
}
#endif
