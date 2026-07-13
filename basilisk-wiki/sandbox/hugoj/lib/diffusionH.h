/**
# Vertical diffusion

This code builds upon [diffusion.h](https://basilisk.fr/src/layered/diffusion.h)
to implement a Neumann condition at the top and the bottom of the domain.

See the test case [here](https://basilisk.fr/sandbox/hugoj/test_diffusion_with_bott_neumann/diff_test.c)

*/


void vertical_diffusion2 (Point point, scalar h, scalar s, double dt, double D,
			 double dst, double dsb)
{
  double a[nl], b[nl], c[nl], rhs[nl];

  /**
  The *rhs* of the tridiagonal system is $h_l s_l$. */
      
  foreach_layer()
    rhs[point.l] = s[]*h[];

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
  For the bottom layer we use the same logic as for the top
  */

  b[0] = h[] + 2.*dt*D*(1./(h[] + h[0,0,1]));
  c[0] = - 2.*dt*D*(1./(h[] + h[0,0,1]));
  rhs[0] -= D*dt*dsb;

  if (nl == 1) {
    b[0] += c[0];
    rhs[0] += (- c[0]*h[] - D*dt) * dst; //Todo: add flux condition dst_b for case nl=1
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

/* # TODO

- fix the case nl=1

**/
