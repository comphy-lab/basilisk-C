/**
# Tracer diffusion 
We consider the diffusion of a tracer $tr$ for the multilayer solver.
Diffusion is divided in two part : vertical diffusion using a Thomas algorithm 
and horizontal diffusion using the Poisson solver. */

/**
## Vertical tracer diffusion
This function is an adaptation of the file layered/diffusion.h for vertical viscosity friction. 
A neumann boudary condition is added at the bottom by default. 

Boundary conditions on the top and bottom layers need to be added to close the
system for diffusion. We chose to impose a Neumann condition on the
top i.e.
$$
\partial_z c |_t = \dot{c}_t
$$
and either a Neumann condition or a Navier slip condition on the bottom i.e.
$$
c|_b = c_b + \lambda_cb \partial_z c|_b
$$
The diffusion coefficient is D.*/
double bottom_boundary_neumann = true;
(const) scalar dct = zeroc, dct_b = zeroc, c_b = zeroc, lambda_cb = zeroc;

/**
For stability, we discretise the diffusion term implicitly as
$$
\frac{(hc_l)^{n + 1} - (hc_l)^{\star}}{\Delta t} =
D \left( \frac{c_{l + 1} - c_l}{h_{l + 1 / 2}} -
\frac{c_l - c_{l - 1}}{h_{l - 1 / 2}} \right)^{n + 1}
$$
which can be expressed as the linear system
$$
\mathbf{Mc}^{n + 1} = \mathrm{rhs}
$$
where $\mathbf{M}$ is a 
[tridiagonal matrix](https://en.wikipedia.org/wiki/Tridiagonal_matrix). 
The lower, principal and upper diagonals are *a*, *b* and *c* respectively. */

void vertical_diffusion2 (Point point, scalar h, scalar s, double D, double dt)
{
  double a[nl], b[nl], c[nl], rhs[nl];

  /**
  The *rhs* of the tridiagonal system is $h_lc_l$. */      
  foreach_layer()
    rhs[_layer] = s[]*h[];

  /**
  The lower, principal and upper diagonals $a$, $b$ and $c$ are given by
  $$
  a_{l > 0} = - \left( \frac{D \Delta t}{h_{l - 1 / 2}} \right)^{n + 1}
  $$
  $$
  c_{l < \mathrm{nl} - 1} = - \left( \frac{D \Delta t}{h_{l + 1 / 2}}
  \right)^{n + 1}
  $$
  $$
  b_{0 < l < \mathrm{nl} - 1} = h_l^{n + 1} - a_l - c_l
  $$
  */
  
  for (int l = 1; l < nl - 1; l++) {
    a[l] = - 2.*D*dt/(h[0,0,l-1] + h[0,0,l]);
    c[l] = - 2.*D*dt/(h[0,0,l] + h[0,0,l+1]);
    b[l] = h[0,0,l] - a[l] - c[l];
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
  a[nl-1] = - 2.*D*dt/(h[0,0,nl-2] + h[0,0,nl-1]);
  b[nl-1] = h[0,0,nl-1] - a[nl-1];
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
  if (bottom_boundary_neumann) {
    c[0] = - 2.*D*dt/(h[0,0,1] + h[]);
    b[0] = h[] - c[0];
    rhs[0] += D*dt*dct_b[];
  }
  else {
  double den = h[]*h[0,0,1]*(8.*lambda_b[] + h[0,0,1]) +
    sq(h[])*(6.*lambda_b[] + 3.*h[0,0,1]) +
    2.*(sq(h[0,0,1])*lambda_b[] + cube(h[]));
  b[0] = h[] + 2.*dt*D*(1./(h[] + h[0,0,1]) +
			  (sq(h[0,0,1]) + 4.*h[]*h[0,0,1] + 4.*sq(h[]))/den);
  c[0] = - 2.*dt*D*(1./(h[] + h[0,0,1]) + sq(h[])/den);
  rhs[0] += 2.*dt*D*c_b[]*(sq(h[0,0,1]) + 4.*h[]*h[0,0,1] + 3.*sq(h[0]))/den;
  }
   
  /**
  We can now solve the tridiagonal system using the [Thomas
  algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
  for (int l = 1; l < nl; l++) {
    b[l] -= a[l]*c[l-1]/b[l-1];
    rhs[l] -= a[l]*rhs[l-1]/b[l-1];
  }
  a[nl-1] = rhs[nl-1]/b[nl-1];
  s[0,0,nl-1] = a[nl-1];
  for (int l = nl - 2; l >= 0; l--)
    s[0,0,l] = a[l] = (rhs[l] - c[l]*a[l+1])/b[l];
}



/**
## Horizontal tracer diffusion
We consider the diffusion of a tracer $tr$ within a layer.
As the thickness $h_l$ is not constant, the diffusion equation within the layer becomes :

$$
h_l \frac{\partial tr_l}{\partial t} = - \frac{\partial}{\partial x} J\, h_l \quad \text{with} 
\quad J = - D\, \nabla_l tr_l \quad
\text{and} 
\quad \nabla_l \approx \frac{1}{\mathcal{A}}\, \nabla
$$
$$
\text{then, with a 1D gradient} \quad \, h_l\, \frac{\partial tr_l}{\partial t} =
\nabla(\frac{D\, h_l\,}{\mathcal{A}} \nabla tr_l)
$$
Therefore we see that to embody horizontal diffusion within the layer,
 we can used the classic diffusion solver if:

* the left term and the flux term are multiplied by $h_l$ to include the variable thickness;
* the diffusion coefficient is divided by the interface length to take into acount the layer curvature.

The inputs of the function are:
* $h$: height field;
* $tr$: diffusive tracer field;
* $D$: diffusion coefficient;
* $dt$: time step.
*/

#include "diffusion.h"
void horizontal_diffusion (scalar h, scalar tr, double D,  double dt)
{
  /**
  We allocate the tracer and parameters fields. */
  vector diffusion_coefficient = new face vector[nl];
  scalar theta_c = new scalar[nl];
  scalar tr_l = new scalar[nl];
  reset({diffusion_coefficient}, 0);
  reset({theta_c}, 0);
  reset({tr_l}, 0);

  /**
  The vector dz_f caracterises the layer curvature. */
  vector dz_f[];
  foreach_face()
    dz_f.x[] = (zb[] - zb[-1])/Delta;
  
  for (int l = 0; l<nl; l++) {
    /**
    In each layer, the tracer is copied and the theta field is computed as the cell height.*/
    foreach(){
      theta_c[] = h[0,0,l];
      tr_l[] = tr[0,0,l];
    }
    boundary({theta_c, tr_l});
  
    /**
    The diffusion coefficient is computed on each face.*/
    foreach_face() {
      dz_f.x[] += (h[0,0,l] - h[-1,0,l])/(2*Delta);
      double area_f = sqrt(1. + sq(dz_f.x[]));           
      diffusion_coefficient.x[] = D*hf.x[0,0,l]/area_f;
      dz_f.x[] += (h[0,0,l] - h[-1,0,l])/(2*Delta);
    }
    boundary((scalar*) {diffusion_coefficient});
  
    /**
    The implicit diffusion solver is called and the original tracer is updated.*/
    diffusion (tr_l, dt, D = diffusion_coefficient, theta = theta_c);
    foreach()
      tr[0,0,l] = tr_l[];
  }
  delete ({theta_c, tr_l, diffusion_coefficient});
}
