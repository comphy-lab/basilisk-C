/**
# Vertical diffusion

This is similar to [https://basilisk.fr/src/layered/diffusion.h](), but
adapted to a two-phase layered model. The layerwise kinematic viscosity
$D \equiv \nu_l$ is allowed to vary with the layer, and the implicit solve
is written for the layer momentum $\rho_l h_l s_l$.

The interface stress uses the dynamic viscosity
$$
\mu_l = \rho_l \nu_l
$$
with the code using the neighbouring-layer combination
$(\mu_l + \mu_{l+1})$ or $(\mu_l + \mu_{l-1})$ in the tridiagonal
coefficients.

We consider the vertical diffusion of a tracer/velocity $s$ for the
multilayer solver. The kinematic viscosity $D$ is prescribed per layer,
while the implicit operator uses the corresponding dynamic viscosity
$\mu_l = \rho_l D_l$.

For stability, we discretise the vertical diffusion equation implicitly as
$$
\frac{(\rho_l h_l s_l)^{n + 1} - (\rho_l h_l s_l)^{\star}}{\Delta t} =
\left( \mu_{l + 1 / 2} \frac{s_{l + 1} - s_l}{h_{l + 1 / 2}} -
\mu_{l - 1 / 2}\frac{s_l - s_{l - 1}}{h_{l - 1 / 2}} \right)^{n + 1}
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
and a Navier slip condition on the bottom i.e.
$$
s|_b = s_b + \lambda_b \partial_z s|_b
$$ */

void vertical_diffusion (Point point, scalar h, scalar s, double dt, double D,
			 double dst, double s_b, double lambda_b)
{
  double a[nl], b[nl], c[nl], rhs[nl];
  double nueq[nl];
  double mueq[nl];

  /**
  The *rhs* of the tridiagonal system is $\rho_l h_l s_l$. */
      
  foreach_layer()
    rhs[point.l] = rho(point.l)*s[]*h[];

  //fprintf (stdout, "diffusion check %g %g %g %g %g %g %g\n", t, x, y, dt, h[], s[], rhs[point.l] );

/** Layerwise kinematic viscosity and the corresponding
dynamic viscosity $\mu_l = \rho_l \nu_l$. */
  for (int l = 0; l < nl; l++) nueq[l] = D;
  //#if TWOPHASE
  foreach_layer(){
    if(point.l<NL) 
       nueq[point.l] = D;
    else
       nueq[point.l] = D*0.0;
    mueq[point.l] = rho(point.l)*nueq[point.l];
  }


  /**
  The lower, principal and upper diagonals $a$, $b$ and $c$ are given by
  $$
  a_{l > 0} = - \left( \frac{\mu_{l - 1 / 2} \Delta t}{h_{l - 1 / 2}} \right)^{n + 1}
  $$
  $$
  c_{l < \mathrm{nl} - 1} = - \left( \frac{\mu_{l + 1 / 2} \Delta t}{h_{l + 1 / 2}}
  \right)^{n + 1}
  $$
  $$
  b_{0 < l < \mathrm{nl} - 1} = \rho_l h_l^{n + 1} - a_l - c_l
  $$
  */

  for (int l = 1; l < nl - 1; l++) {
    a[l] = - (mueq[l] + mueq[l-1])*dt/(h[0,0,l-1] + h[0,0,l]);
    c[l] = - (mueq[l] + mueq[l+1])*dt/(h[0,0,l] + h[0,0,l+1]);
    b[l] = rho(l)*h[0,0,l] - a[l] - c[l];
    //fprintf (stdout, "diffusion check1 %g %g %g %d %g %g %g\n", t, x, y, l, a[l], b[l], c[l] );
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
  (\rho hs)_{\mathrm{nl} - 1}^{\star} + \mu \Delta t \dot{s}_t
  $$
  */

  a[nl-1] = - (mueq[nl-1] + mueq[nl-2])*dt/(h[0,0,nl-2] + h[0,0,nl-1]);
  b[nl-1] = rho(nl - 1)*h[0,0,nl-1] - a[nl-1];
  rhs[nl-1] += mueq[nl-1]*dt*dst;

  /**
  For the bottom layer a third-order discretisation of the Navier slip
  condition gives
  $$
  \begin{aligned}
  b_0 & = \rho_0 h_0 + 2 \Delta t \mu_0 \left( \frac{1}{h_0 + h_1} + \frac{h^2_1 + 3
  h_0 h_1 + 3 h^2_0}{\det} \right),\\
  c_0 & = - 2 \Delta t \mu_0 \left( \frac{1}{h_0 + h_1} + \frac{h^2_0}{\det}
  \right),\\
  \text{rhs}_0 & = (\rho hs_0)^{\star} + 2 \Delta t \mu_0 s_b  \frac{h^2_1 + 3 h_0
  h_1 + 2 h^2_0}{\det},\\
  \det & = h_0 (h_0 + h_1)^2  + 2\lambda (3\,h_0 h_1 + 2\,h_0^2 + h_1^2),
  \end{aligned}
  $$
  */

  double den = h[]*sq(h[] + h[0,0,1]) 
    + 2.*lambda_b*(3.*h[]*h[0,0,1] + 2.*sq(h[]) + sq(h[0,0,1]));
  b[0] = rho(0)*h[] + 2.*dt*mueq[0]*(1./(h[] + h[0,0,1]) +
			  (sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 3.*sq(h[]))/den);
  c[0] = - 2.*dt*mueq[0]*(1./(h[] + h[0,0,1]) + sq(h[])/den);

  rhs[0] += 2.*dt*mueq[0]*s_b*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 2.*sq(h[0]))/den;

  if (nl == 1) {
    b[0] += c[0];
    rhs[0] +=  (- c[0]*h[] - mueq[0]*dt) * dst;
  }
  //fprintf (stdout, "diffusion checkinit %g %g %g %g %g %g %g %g %g %g %g\n", t, x, y, dt, s[], a[nl-1], b[nl-1], rhs[nl-1], c[0], b[0], rhs[0] );

    
  /**
  We can now solve the tridiagonal system using the [Thomas
  algorithm](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm). */
  
  for (int l = 1; l < nl; l++) {
    b[l] -= a[l]*c[l-1]/b[l-1];
    rhs[l] -= a[l]*rhs[l-1]/b[l-1];
    //fprintf (stdout, "diffusion check %g %g %g %d %g %g %g %g\n", t, x, y, l, a[l], b[l], c[l], rhs[l] );
  }
  //fprintf (stdout, "out\n" );

  a[nl-1] = rhs[nl-1]/b[nl-1];
  s[0,0,nl-1] = a[nl-1];
  for (int l = nl - 2; l >= 0; l--)
    s[0,0,l] = a[l] = (rhs[l] - c[l]*a[l+1])/b[l];
  //fprintf(stdout, "u %g\n",s[]);
}

/**
# Viscous friction between layers

By default the viscosity is zero and we impose free-slip on the
free-surface and no-slip on the bottom boundary
i.e. $\dot{\mathbf{u}}_t = 0$, $\mathbf{\lambda}_b = 0$, $\mathbf{u}_b
= 0$. */

double nu = 0.;
(const) vector lambda_b[] = {0,0,0}, dut[] = {0,0,0}, u_b[] = {0,0,0};

/**
In the [layered solver](hydro.h), vertical viscosity is applied to the
velocity field just after advection, but before the pressure
gradient/acceleration term is applied. To take the pressure gradient
into account, we first apply the acceleration of the previous
timestep, apply vertical viscosity and then substract the previous
acceleration. */

event viscous_term (i++,last)
{
  //fprintf(stdout, "dtviscous t=%g dt=%g\n", t, dt);
  if (nu > 0.) {
    foreach() { 
      if (h[] > dry){
      foreach_layer()
	foreach_dimension()
	  u.x[] += dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
      foreach_dimension()
	vertical_diffusion (point, h, u.x, dt, nu,
			    dut.x[], u_b.x[], lambda_b.x[]);
      foreach_layer()
	foreach_dimension()
	  u.x[] -= dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
    }
    }
  }
}

/**
## Horizontal diffusion

This approximates
$$
h \partial_t s = D \nabla \cdot (h \nabla s)
$$
with $D$ the diffusion coefficient. Note that metric terms linked to
the slope of the layers are not taken into account. Note also that the
time discretisation is explicit so that the timestep must be limited
(manually) by $\min(\Delta^2/D)$. */

void horizontal_diffusion (scalar * list, double D, double dt)
{
  if (D > 0.) {
    scalar * d2sl = list_clone (list);
    foreach_layer() {
      foreach() {
	scalar s, d2s;
	for (s,d2s in list,d2sl) {
	  double a = 0.;
	  foreach_dimension()
	    a += (hf.x[]*fm.x[]/(cm[-1] + cm[])*(s[-1] - s[]) +
		  hf.x[1]*fm.x[1]/(cm[1] + cm[])*(s[1] - s[]));
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

/**
### Time-implicit integration

When the timestep constraint due to explicit time integration is too
strong, one can use the time-implicit version below. The diffusion equation
$$
h \partial_t s = D \nabla \cdot (h \nabla s)
$$
is discretised implicitly as
$$
h \frac{s^{n + 1} - s^n}{\Delta t} = \nabla \cdot (Dh \nabla s^{n + 1})
$$
which can be reordered as
$$
\nabla \cdot (\Delta tDh \nabla s^{n + 1}) - hs^{n + 1} = - hs^n
$$
i.e. a Poisson--Helmoltz equation of the form
$$
\nabla \cdot (\alpha \nabla s^{n + 1}) + \lambda s^{n + 1} = b
$$
which is solved using the [multigrid Poisson--Helmholz solver](/src/poisson.h#user-interface).

The function returns convergence statistics for the multigrid solver. */

#include "poisson.h"

trace
mgstats implicit_horizontal_diffusion (scalar * list, double D, double dt)
{
  mgstats mg = {0};
  if (D > 0.) {
    scalar b[], lambda[];
    face vector alpha[];
    foreach_layer() {
      foreach_face()
	alpha.x[] = 2.*dt*D*hf.x[]*fm.x[]/(cm[-1] + cm[]);
      foreach()
	lambda[] = - cm[]*h[];
      for (scalar s in list) {
	foreach()
	  b[] = lambda[]*s[];
	mg = poisson (s, b, alpha, lambda);
      }
    }
  }
  return mg;
}

/**
## References

~~~bib
@hal{popinet2020, hal-02365730}

@hal{devita2019, hal-02295398}
~~~
*/
