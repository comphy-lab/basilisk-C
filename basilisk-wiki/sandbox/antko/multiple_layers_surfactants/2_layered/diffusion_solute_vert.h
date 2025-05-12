/**
# Solute diffusion in and between layers
This file is an adaptation of the file layered/diffusion.h for vertical viscosity friction. 
A neumann boudary condition is added at the bottom by default. 

Boundary conditions on the top and bottom layers need to be added to close the
system for the viscous stresses. We chose to impose a Neumann condition on the
free-surface i.e.
$$
\partial_z c |_t = \dot{c}_t
$$
and a Navier slip condition on the bottom i.e.
$$
c|_b = c_b + \lambda_cb \partial_z c|_b
$$
By default the diffusion coefficient is zero.*/
#include "diffusion.h"
double D = 0.;
(const) scalar dct = zeroc, dct_b = zeroc, c_b = zeroc, lambda_cb = zeroc;
double bottom_boundary_neumann = true;
/**
For stability, we discretise the diffusion term implicitly as
$$
\frac{(hc_l)^{n + 1} - (hc_l)^{\star}}{\Delta t} =
\nu \left( \frac{c_{l + 1} - c_l}{h_{l + 1 / 2}} -
\frac{c_l - c_{l - 1}}{h_{l - 1 / 2}} \right)^{n + 1}
$$
which can be expressed as the linear system
$$
\mathbf{Mc}^{n + 1} = \mathrm{rhs}
$$
where $\mathbf{M}$ is a 
[tridiagonal matrix](https://en.wikipedia.org/wiki/Tridiagonal_matrix). 
The lower, principal and upper diagonals are *a*, *b* and *c* respectively. */

void vertical_diffusion (Point point, scalar * hl, scalar * sl, double dt)
{
  double a[nl], b[nl], c[nl], rhs[nl];

  /**
  The *rhs* of the tridiagonal system is $h_lc_l$. */
      
  int l = 0;
  scalar s, h;
  for (s,h in sl,hl)
    rhs[l++] = s[]*h[];

  /**
  The lower, principal and upper diagonals $a$, $b$ and $c$ are given by
  $$
  a_{l > 0} = - \left( \frac{\nu \Delta t}{h_{l - 1 / 2}} \right)^{n + 1}
  $$
  $$
  c_{l < \mathrm{nl} - 1} = - \left( \frac{\nu \Delta t}{h_{l + 1 / 2}}
  \right)^{n + 1}
  $$
  $$
  b_{0 < l < \mathrm{nl} - 1} = h_l^{n + 1} - a_l - c_l
  $$
  */
  
  for (l = 1; l < nl - 1; l++) {
    scalar hm = hl[l-1], h = hl[l], hp = hl[l+1];
    a[l] = - 2.*D*dt/(hm[] + h[]);
    c[l] = - 2.*D*dt/(h[] + hp[]);
    b[l] = h[] - a[l] - c[l];
  }
    
  /**
  For the top layer the boundary conditions give the (ghost)
  boundary value
  $$
  c_{\mathrm{nl}} = c_{\mathrm{nl} - 1} + \dot{c}_t h_{\mathrm{nl} - 1},
  $$
  which gives the diagonal coefficient and right-hand-side
  $$
  b_{\mathrm{nl} - 1} = h_{\mathrm{nl} - 1}^{n + 1}
  - a_{\mathrm{nl} - 1}
  $$
  $$
  \mathrm{rhs}_{\mathrm{nl} - 1} = 
  (hc)_{\mathrm{nl} - 1}^{\star} + D \Delta t \dot{c}_t
  $$
  */

  scalar hm = hl[nl-2]; h = hl[nl-1];
  a[nl-1] = - 2.*D*dt/(hm[] + h[]);
  b[nl-1] = h[] - a[nl-1];
  rhs[nl-1] += D*dt*dct[];

  /**
  For the bottom layer a third-order discretisation of the Navier slip
  condition gives
  $$
  \begin{aligned}
  b_0 & = h_0 + 2 \Delta t D \left( \frac{1}{h_0 + h_1} + \frac{h^2_1 + 4
  h_0 h_1 + 4 h^2_0}{\det} \right),\\
  c_0 & = - 2 \Delta t D \left( \frac{1}{h_0 + h_1} + \frac{h^2_0}{\det}
  \right),\\
  \text{rhs}_0 & = (hu_0)^{\star} + 2 \Delta t D c_b  \frac{h^2_1 + 4 h_0
  h_1 + 3 h^2_0}{\det},\\
  \det & = h_0 h_1  (8 \lambda_cb + h_1) + h^2_0  (6 \lambda_cb + 3 h_1) + 2
  (h^2_1 \lambda_cb + h^3_0),
  \end{aligned}
  $$
  */
  scalar h1 = hl[1], h0 = hl[0];
  if (bottom_boundary_neumann) {
    c[0] = - 2.*D*dt/(h1[] + h0[]);
    b[0] = h0[] - c[0];
    rhs[0] += D*dt*dct_b[];
  }
  else {
    double den = h[0]*(8.*h[1]*lambda_cb[] + sq(h1[])) +
      sq(h0[])*(6.*lambda_cb[] + 3.*h[1]) +
      2.*(sq(h[1])*lambda_cb[] + cube(h0[]));
    b[0] = h0[] + 2.*dt*D*(1./(h0[] + h1[]) +
			    (sq(h1[]) + 4.*h[0]*h[1] + 4.*sq(h0[]))/den);
    c[0] = - 2.*dt*D*(1./(h0[] + h1[]) + sq(h0[])/den);
    rhs[0] += 2.*dt*D*c_b[]*(sq(h1[]) + 4.*h[0]*h[1] + 3.*sq(h[0]))/den;
  }
   
  /**
  We can now solve the tridiagonal system using the [Thomas
  algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
  
  for (l = 1; l < nl; l++) {
    b[l] -= a[l]*c[l-1]/b[l-1];
    rhs[l] -= a[l]*rhs[l-1]/b[l-1];
  }
  a[nl-1] = rhs[nl-1]/b[nl-1];
  s = sl[nl-1]; s[] = a[nl-1];
  for (l = nl - 2; l >= 0; l--) {
    s = sl[l];
    s[] = a[l] = (rhs[l] - c[l]*a[l+1])/b[l];
  }
}

/**
In the [layered solver](hydro.h), vertical viscosity is applied to the
velocity field just after advection, but before the pressure
gradient/acceleration term is applied. To take the pressure gradient
into account, we first apply the acceleration of the previous
timestep, apply vertical viscosity and then substract the previous
acceleration. */

event diffusive_term (i++,last)
{
  if (D > 0.) {
    struct { scalar * x, * y; } list = {NULL, NULL};
    foreach_dimension()
      for (scalar c in cl)
	list.x = list_append (list.x, c);
    foreach()
      foreach_dimension()
	      vertical_diffusion (point, hl, list.x, dt);
    boundary ((scalar *) cl);
    foreach_dimension()
      free (list.x);
  }
}

/**
## References

~~~bib
@Article{popinet2019,
  author = 	 {S. Popinet},
  title = 	 {A vertically-Lagrangian, non-hydrostatic, multilayer model
                  for multiscale free-surface flows},
  journal = 	 {Journal of Computational Physics},
  note =         {Submitted},
  year = 	 {2019}
}

@Article{devita2019,
  author =  { F. De Vita and P.Y. Lagrée and S. Chibbaro and S. Popinet},
  title =   {Beyond Shallow Water: Appraisal for a numerical approach to hydraulic jumps based upon the Boundary Layer theory},
  journal = {European Journal of Mechanics - B/Fluids},
  volume  = {79/C},
  pages =   {233-246},
  year =    {2019},
  url = {https://hal.archives-ouvertes.fr/hal-02295398}
}
~~~
*/
