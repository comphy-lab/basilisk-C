/** 
# Disclamer
We consider a copy of 
[Multilayer hydrostatic](http://basilisk.fr/src/layered/hydro.h)
The bottom boundary condition had been changed for a neumann BC


# Vertical diffusion

We consider the vertical diffusion of a tracer $s$ with a diffusion
coefficient $D$ for the multilayer solver.

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
and at the bottom i.e.
$$
\partial_z s |_b = \dot{s}_b
$$ */
void vertical_diffusion_neumann (Point point, scalar h, scalar s, double dt, double D,
                         double dst, double dsb)
{
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
  rhs[nl-1] += D*dt*dst;

  /**
  For the bottom layer Neumann condition gives:
  */
  c[0] = -2.*D*dt/(h[0,0,1]+h[]);
  b[0] = h[] - c[0];
  rhs[0] += dt*D*dsb;

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

void horizontal_diffusion_modif (scalar * list, double D, double dt)
{
  if (D > 0.) {
    scalar * d2sl = list_clone (list);
    foreach_layer() {
      foreach() {
        scalar s, d2s;
        for (s,d2s in list,d2sl) {
          double a = 0.;
          foreach_dimension(){
                  a += (hf.x[]*fm.x[]/(cm[-1] + cm[])*(s[-1] - s[]) +
                  hf.x[1]*fm.x[1]/(cm[1] + cm[])*(s[1] - s[]));
                 }
          d2s[] = 2.*a/(cm[]*sq(Delta));
        }
      }
      foreach()
        if (h[] > dry) {
          scalar s, d2s;
          for (s,d2s in list,d2sl)
            s[] += dt*D*d2s[]/h[];
        }
    }
    delete (d2sl);
    free (d2sl);
  }
}
