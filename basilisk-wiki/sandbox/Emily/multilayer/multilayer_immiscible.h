/**
# Immiscible multilayer Saint-Venant system with variable $\rho$ and $\mu$...

The [Saint-Venant system](saint-venant.h) is extended to multiple
layers following [Audusse et al, 2011](references.bib#audusse2011) as
$$
\partial_t\sum_{l=0}^{nl-1}h_l + \partial_x\sum_{l=0}^{nl-1}h_lu_l = 0
$$
where
$$
h=\sum_{l=0}^{nl-1} h_l.
$$
The momentum equation in each layer is thus
$$
\partial_t(h\mathbf{u}_l) + \nabla\cdot\left(h\mathbf{u}_l\otimes\mathbf{u}_l + 
\frac{gh^2}{2}\mathbf{I}\right) 
= 
- g h_l \nabla \left(z_b + \sum_{j \ne l} \min \left( \frac{\rho_j}{\rho_l} , 1 \right)h_j\right) 
+ \frac{2\mu_l}{\rho_l}\left(\frac{u_{l+1} - u_l}{h_{l+1}+h_l} - 
\frac{u_{l} - u_{l-1}}{h_{l-1}+h_l}\right)
$$
where the final term corresponds to viscous friction between layers .


The thickness of each of the layers is stored in *hhl*, the horizontal velocity in each layer is stored in *ul* and the
vertical velocity between layers in *wl*. */

scalar * hhl = NULL;
vector * ul = NULL;
scalar * wl = NULL;
/**
## Viscous friction between layers

Boundary conditions on the top and bottom layers need to be added to close the
system for the viscous stresses. We chose to impose a Neumann condition on the
top boundary i.e.
$$
partial_z u |_t = \dot{u}_t
$$
and a Navier slip condition on the bottom i.e.
$$
u|_b = u_b + \lambda_b \partial_z u|_b
$$
By default the viscosity is zero and we impose free-slip on the top
boundary and no-slip on the bottom boundary i.e. $\dot{u}_t = 0$,
$\lambda_b = 0$, $u_b = 0$. */
double * mul;
double * rhol;
(const) scalar lambda_b = zeroc, dut = zeroc, u_b = zeroc;

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

void vertical_viscosity (Point point, scalar * hl, vector * ul, double dt)
{
  if (mul)   // This should be a vector... We'll start with it being the same for all layers
    return;
  
  double a[nl], b[nl], c[nl], rhs[nl], hcol[nl];

  foreach_dimension() {

    /**
    The *rhs* of the tridiagonal system is $h_lu_l = h\mathrm{layer}_lu_l$. */
    
    int l = 0;
    for (scalar h in hl)
      hcol[l]=h[], l++;
    l = 0;
    for (vector u in ul)
      rhs[l] = hcol[l]*u.x[], l++;

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
      a[l] = - 2.*mul[l]/rhol[l]*dt/(hcol[l-1] + hcol[l]);
      c[l] = - 2.*mul[l]/rhol[l]*dt/(hcol[l] + hcol[l+1]);
      b[l] = hcol[l] - a[l] - c[l];
    }
    
    /**
       CHECK STILL CORRECT!!!
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

    a[nl-1] = - 2.*mul[nl-1]/rhol[nl-1]*dt/(hcol[nl-2] + hcol[nl-1]);
    b[nl-1] = hcol[nl-1] - a[nl-1];
    rhs[nl-1] += mul[nl-1]/rhol[nl-1]*dt*dut[];

    /**
       CHECK STILL CORRECT!!!
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

    c[0] = - 2.*dt*mul[0]/rhol[0]/(hcol[0] + hcol[1]);
    b[0] = hcol[0] - c[0] + 2.*mul[0]/rhol[0]*dt/(2.*lambda_b[] + hcol[0]);
    rhs[0] += 2.*mul[0]/rhol[0]*dt/(2.*lambda_b[] + hcol[0])*u_b[];
    
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
REMOVED FROM THIS VERSION AS SET TO ZERO BY DEFAULT. WILL NEED TO REINSTATE IN MIXED VERSION
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

