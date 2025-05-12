/**
#Vertical diffusion
We consider the vertical diffusion of a tracer $s$ for the multilayer solver.

For stability, we discretise the vertical diffusion equation implicitly as
$$
\frac{(hs_l)^{n + 1} - (hs_l)^{\star}}{\Delta t} =
D \left( \frac{s_{l + 1} - s_l}{h_{l + 1 / 2}} -
\frac{s_l - s_{l - 1}}{h_{l - 1 / 2}} \right)^{n + 1}
$$
which can be expressed as the linear system
$$
\mathbf{Ms}^{n + 1} = \mathrm{rhs}
$$
where $\mathbf{M}$ is a 
[tridiagonal matrix](https://en.wikipedia.org/wiki/Tridiagonal_matrix). 
The lower, principal and upper diagonals are *a*, *b* and *c* respectively. 

Boundary conditions on the top and bottom layers need to be added to close the
system. We chose to impose a Neumann condition on the free-surface i.e.
$$
\partial_z s |_t = \dot{s}_t
$$
On the bottom, either a Neumann condition or a Navier slip condition (for viscosity diffusion) is imposed i.e.
$$
\partial_z s |_b = \dot{s}_b
\textrm{ or }
s|_b = s_b + \lambda_b \partial_z s|_b
$$
*/

struct Vertical_Diffusion {
  // mandatory
  scalar h;
  scalar s;
  double dt;
  double D;
  // optional
  double dst;      //default 0
  double dsb;      //default 0
  double s_b;      //default 0
  double lambda_b; //default 0
};

bool Neum_b = false;
void vertical_diffusion (Point point, struct Vertical_Diffusion p)
{
  scalar h = p.h; scalar s = p.s; 
  double dt = p.dt; double D = p.D;
  double a[nl], b[nl], c[nl], rhs[nl];
  
  /**
  The *rhs* of the tridiagonal system is $h_l s_l$. */
      
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
  s_{\mathrm{nl}} = s_{\mathrm{nl} - 1} + \dot{s}_t h_{\mathrm{nl} - 1},
  $$
  which gives the diagonal coefficient and right-hand-side
  $$
  b_{\mathrm{nl} - 1} = h_{\mathrm{nl} - 1}^{n + 1}
  - a_{\mathrm{nl} - 1}
  $$
  $$
  \mathrm{rhs}_{\mathrm{nl} - 1} = 
  (hs)_{\mathrm{nl} - 1}^{\star} + D \Delta t \dot{s}_t
  $$
  */

  a[nl-1] = - 2.*D*dt/(h[0,0,nl-2] + h[0,0,nl-1]);
  b[nl-1] = h[0,0,nl-1] - a[nl-1];
  rhs[nl-1] += p.dst ? D*dt*p.dst : 0;

  /**
  For the bottom layer a third-order discretisation of the Navier slip
  condition gives
  $$
  \begin{aligned}
  b_0 & = h_0 + 2 \Delta t \nu \left( \frac{1}{h_0 + h_1} + \frac{h^2_1 + 3
  h_0 h_1 + 3 h^2_0}{\det} \right),\\
  c_0 & = - 2 \Delta t \nu \left( \frac{1}{h_0 + h_1} + \frac{h^2_0}{\det}
  \right),\\
  \text{rhs}_0 & = (hs_0)^{\star} + 2 \Delta t \nu u_b  \frac{h^2_1 + 3 h_0
  h_1 + 2 h^2_0}{\det},\\
  \det & = h_0 (h_0 + h_1)^2  + 2\lambda (3\,h_0 h_1 + 2\,h_0^2 + h_1^2),
  \end{aligned}
  $$
  */
  if (!Neum_b) {
    double lambda_b = p.lambda_b ? p.lambda_b : 0;
    double s_b = p.s_b ? p.s_b : 0;
    double den = h[]*sq(h[] + h[0,0,1]) 
    + 2.*lambda_b*(3.*h[]*h[0,0,1] + 2.*sq(h[]) + sq(h[0,0,1]));
    b[0] = h[] + 2.*dt*D*(1./(h[] + h[0,0,1]) +
  			  (sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 3.*sq(h[]))/den);
    c[0] = - 2.*dt*D*(1./(h[] + h[0,0,1]) + sq(h[])/den);
  rhs[0] += 2.*dt*D*s_b*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 2.*sq(h[0]))/den;

  if (nl == 1) {
    b[0] += c[0];
    rhs[0] += (- c[0]*h[] - D*dt) * (p.dst ? p.dst : 0);
  }
  
  }
  /**
  By default, a Neumann boundary condition is applied.*/
  else {
    c[0] = - 2.*D*dt/(h[0,0,1] + h[]);
    b[0] = h[] - c[0];
    rhs[0] += p.dsb ? D*dt*p.dsb : 0;
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
# Viscous friction between layers

By default the viscosity is zero and we impose free-slip on the
free-surface and no-slip on the bottom boundary i.e. $\dot{u}_t = 0$,
$\lambda_b = 0$, $u_b = 0$. */

double nu = 0.;
(const) scalar lambda_b = zeroc, u_b = zeroc;
(const) vector dut = zerof;

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

bool horvisco = false;
#include "diffusion.h"
void horizontal_diffusion_v (scalar h, scalar tr, double D,  double dt)
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


/**
In the [layered solver](hydro.h), vertical viscosity is applied to the
velocity field just after advection, but before the pressure
gradient/acceleration term is applied. To take the pressure gradient
into account, we first apply the acceleration of the previous
timestep, apply vertical viscosity and then substract the previous
acceleration. */

event viscous_term (i++,last)
{
  if (nu > 0.) {
    foreach() {
      foreach_layer()
	foreach_dimension()
	  u.x[] += dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
      foreach_dimension()
	vertical_diffusion (point, (struct Vertical_Diffusion){h, u.x, dt, nu, dut.x[], 0, u_b[], lambda_b[]});
      foreach_layer()
	foreach_dimension()
	  u.x[] -= dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
    }
    boundary ((scalar *) {u});
    if(horvisco)
      horizontal_diffusion_v(h, u.x, nu, dt);
    boundary ((scalar *) {u});
  }
}

/**
## References

~~~bib
@hal{popinet2020, hal-02365730}

@hal{devita2019, hal-02295398}
~~~
*/