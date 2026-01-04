/**
# Multilayer Saint-Venant system with mass exchanges

<div class="message">
Note that the [multilayer solver](layered/hydro.h) provides the same
functionality and should be prefered for most applications.</div>

The [Saint-Venant system](saint-venant.h) is extended to multiple
layers following [Audusse et al, 2011](references.bib#audusse2011) as
$$
\partial_th + \partial_x\sum_{l=0}^{nl-1}h_lu_l = 0
$$
with
$$
h_l = \mathrm{layer}_lh
$$
with $\mathrm{layer}_l$ the relative thickness of the layers satisfying
$$
\mathrm{layer}_l >= 0,\;\sum_{l=0}^{nl - 1}\mathrm{layer}_l = 1.
$$
The momentum equation in each layer is thus
$$
\partial_t(h\mathbf{u}_l) + \nabla\cdot\left(h\mathbf{u}_l\otimes\mathbf{u}_l + 
\frac{gh^2}{2}\mathbf{I}\right) = 
- gh\nabla z_b + \frac{1}{\mathrm{layer}_l}\left[\mathbf{u}_{l+1/2}G_{l+1/2} - 
\mathbf{u}_{l-1/2}G_{l-1/2}
+ \nu\left(\frac{u_{l+1} - u_l}{h_{l+1/2}} - 
\frac{u_{l} - u_{l-1}}{h_{l-1/2}}\right)\right]
$$
where $G_{l+1/2}$ is the relative vertical transport velocity between
layers and the second term corresponds to viscous friction between
layers.

These last two terms are the only difference with the [one layer
system](saint-venant.h). 

The horizontal velocity in each layer is stored in *ul* and the
vertical velocity between layers in *wl*. */

vector * ul = NULL;
scalar * wl = NULL;
double * layer;

/**
The index of the layer is set as a field attribute. */

attribute {
  int l;
}

/**
## Viscous friction between layers

Boundary conditions on the top and bottom layers need to be added to close the
system for the viscous stresses. We chose to impose a Neumann condition on the
top boundary i.e.
$$
\partial_z u |_t = \dot{u}_t
$$
and a Navier slip condition on the bottom i.e.
$$
u|_b = u_b + \lambda_b \partial_z u|_b
$$
By default the viscosity is zero and we impose free-slip on the top
boundary and no-slip on the bottom boundary i.e. $\dot{u}_t = 0$,
$\lambda_b = 0$, $u_b = 0$. */

double nu = 0.;
const scalar lambda0[] = 0, dut0[] = 0, u_b0[] = 0;
(const) scalar lambda_b = lambda0, dut = dut0, u_b = u_b0;
// fixme: should be able to replace the two lines above with:
// (const) scalar lambda_b[] = 0, dut[] = 0, u_b[] = 0;

/**
For stability, we discretise the viscous friction term implicitly as
$$
\frac{(hu_l)_{n + 1} - (hu_l)_{\star}}{\Delta t} =
\frac{\nu}{\mathrm{layer}_l}  \left( \frac{u_{l + 1} - u_l}{h_{l + 1 / 2}} -
\frac{u_l - u_{l - 1}}{h_{l - 1 / 2}} \right)_{n + 1}
$$
which can be expressed as the linear system
$$
\mathbf{Mu}_{n + 1} = \mathrm{rhs}
$$
where $\mathbf{M}$ is a 
[tridiagonal matrix](https://en.wikipedia.org/wiki/Tridiagonal_matrix). 
The lower, principal and upper diagonals are *a*, *b* and *c* respectively. */

void vertical_viscosity (Point point, double h, vector * ul, double dt)
{
  if (nu == 0.)
    return;
  
  double a[nl], b[nl], c[nl], rhs[nl];

  foreach_dimension() {

    /**
    The *rhs* of the tridiagonal system is $h_lu_l = h\mathrm{layer}_lu_l$. */
      
    int l = 0;
    for (vector u in ul)
      rhs[l] = h*layer[l]*u.x[], l++;

    /**
    The lower, principal and upper diagonals $a$, $b$ and $c$ are given by
    $$
    a_{l > 0} = - \left( \frac{\nu \Delta t}{h_{l - 1 / 2}} \right)_{n + 1}
    $$
    $$
    c_{l < \mathrm{nl} - 1} = - \left( \frac{\nu \Delta t}{h_{l + 1 / 2}}
    \right)_{n + 1}
    $$
    $$
    b_{0 < l < \mathrm{nl} - 1} = \mathrm{layer}_l h_{n + 1} - a_l - c_l
    $$
    */
    
    for (l = 1; l < nl - 1; l++) {
      a[l] = - 2.*nu*dt/(h*(layer[l-1] + layer[l]));
      c[l] = - 2.*nu*dt/(h*(layer[l] + layer[l+1]));
      b[l] = layer[l]*h - a[l] - c[l];
    }
    
    /**
    For the top layer the boundary conditions give the (ghost)
    boundary value
    $$
    u_{\mathrm{nl}} = u_{\mathrm{nl} - 1} + \dot{u}_t h_{\mathrm{nl} - 1},
    $$
    which gives the diagonal coefficient and right-hand-side
    $$
    b_{\mathrm{nl} - 1} = \mathrm{layer}_{\mathrm{nl} - 1} h_{n + 1}
    - a_{\mathrm{nl} - 1}
    $$
    $$
    \mathrm{rhs}_{\mathrm{nl} - 1} = \mathrm{layer}_{\mathrm{nl} - 1}  
    (hu_{\mathrm{nl} - 1})_{\star} + \nu \Delta t \dot{u}_t
    $$
    */

    a[nl-1] = - 2.*nu*dt/(h*(layer[max(0,nl-2)] + layer[nl-1]));
    b[nl-1] = layer[nl-1]*h - a[nl-1];
    rhs[nl-1] += nu*dt*dut[];

    /**
    For the bottom layer, the boundary conditions give the (ghost)
    boundary value $u_{- 1}$
    $$
    u_{- 1} = \frac{2 h_0}{2 \lambda_b + h_0} u_b + \frac{2 \lambda_b - h_0}{2
    \lambda_b + h_0} u_0,
    $$
    which gives the diagonal coefficient and right-hand-side
    $$
    b_0 = \mathrm{layer}_0 h_{n + 1} - c_0 + 
    \frac{2 \nu \Delta t}{2 \lambda_b + h_0}
    $$
    $$
    \mathrm{rhs}_0 = \mathrm{layer}_0  (hu_0)_{\star} + \frac{2 \nu \Delta t}{2
    \lambda_b + h_0} u_b
    $$
    */

    c[0] = - 2.*dt*nu/(h*(layer[0] + layer[min(1,nl-1)]));
    b[0] = layer[0]*h - c[0] + 2.*nu*dt/(2.*lambda_b[] + h*layer[0]);
    rhs[0] += 2.*nu*dt/(2.*lambda_b[] + h*layer[0])*u_b[];
    
    /**
    We can now solve the tridiagonal system using the [Thomas
    algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
    
    for (l = 1; l < nl; l++) {
      b[l] -= a[l]*c[l-1]/b[l-1];
      rhs[l] -= a[l]*rhs[l-1]/b[l-1];
    }
    vector u = ul[nl-1];
    u.x[] = a[nl-1] = rhs[nl-1]/b[nl-1];
    for (l = nl - 2; l >= 0; l--) {
      u = ul[l];
      u.x[] = a[l] = (rhs[l] - c[l]*a[l+1])/b[l];
    }
  }
}

/**
## Fluxes between layers

The relative vertical velocity between layers $l$ and $l+1$ is defined
as (eq. (2.22) of [Audusse et al, 2011](references.bib#audusse2011))
$$
G_{l+1/2} = \sum_{j=0}^{l}(\mathrm{div}_j + \mathrm{layer}_j\mathrm{dh})
$$
with
$$
\mathrm{div}_l = \nabla\cdot(h_l\mathbf{u}_l)
$$
$$
\mathrm{dh} = - \sum_{l=0}^{nl-1} \mathrm{div}_l
$$
*/

void vertical_fluxes (vector * evolving, vector * updates,
		      scalar * divl, scalar dh)
{
  foreach() {
    double Gi = 0., sumjl = 0.;
    for (int l = 0; l < nl - 1; l++) {
      scalar div = divl[l];
      Gi += div[] + layer[l]*dh[];
      sumjl += layer[l];
      scalar w = div;
      w[] = dh[]*sumjl - Gi;
      foreach_dimension() {

	/**
	To compute the vertical advection term, we need an estimate of
	the velocity at $l+1/2$. This is obtained using simple
	upwinding according to the sign of the interface velocity
	$\mathrm{Gi} = G_{l+1/2}$ and the values of the velocity in
	the $l$ and $l+1$ layers. Note that the inequality of
	upwinding is consistent with equs. (5.110) of [Audusse et al,
	2011](references.bib#audusse2011) and (77) of [Audusse et al,
	2011b](references.bib#audusse2011b) but not with eq. (2.23) of
	[Audusse et al, 2011](references.bib#audusse2011). */

	scalar ub = evolving[l].x, ut = evolving[l + 1].x;
	double ui = Gi < 0. ? ub[] : ut[];
	
	/**
	The flux at $l+1/2$ is then added to the updates of the bottom
	layer and substracted from the updates of the top layer. */
	
	double flux = Gi*ui;
	scalar du_b = updates[l].x, du_t = updates[l + 1].x;
	du_b[] += flux/layer[l];
	du_t[] -= flux/layer[l + 1];

	/**
	To compute the vertical velocity we use the definition of the
	mass flux term (eq. 2.13 of [Audusse et
	al, 2011](references.bib#audusse2011)):
	$$
	\mathrm{w}(\mathbf{x},z_{l+1/2}) = 
          \partial_t z_{l+1/2} - G_{l+1/2} + \mathbf{u}_{l+1/2}
	  \cdot \nabla z_{l+1/2}
	$$
	We can write the vertical position of the interface as:
	$$
	z_{l+1/2} = z_{b} + \sum_{j=0}^{l} h_{j}
	$$
        so that the vertical velocity is:
	$$
	\mathrm{w}(\mathbf{x},z_{l+1/2}) = 
          \mathrm{dh}\sum_{j=0}^{l}\mathrm{layer}_{j} - G_{l+1/2} + 
          \mathbf{u}_{l+1/2} \cdot \left[\nabla z_{b} + \nabla h 
				   \sum_{j=0}^{l}\mathrm{layer}_{j}\right]
	$$
	*/
	
	w[] += ui*((zb[1] - zb[-1]) + (h[1] - h[-1])*sumjl)/(2.*Delta);
      }
    }
  }
}
