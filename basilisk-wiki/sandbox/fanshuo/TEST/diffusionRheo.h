/**
# Variable-Viscosity Momentum Diffusion
The vertical diffusion problem with a constant diffusivity $"D"$ was formulated in [http://basilisk.fr/src/layered/diffusion.h](). Now our objective is to extend the previous framework to cases where the $"D"$ varies in layers.

We recall that the vertical diffusion of a tracer $s$ in the layer $l$ is governed by
$$
\dfrac{(h_ls_l)^{n + 1} - (h_ls_l)^{\star}}{\Delta t} =
\left( D_{l + 1 / 2} \dfrac{s_{l + 1} - s_l}{h_{l + 1 / 2}} -
D_{l - 1 / 2}\dfrac{s_l - s_{l - 1}}{h_{l - 1 / 2}} \right)^{n + 1}.
$$
Here the second member represents the difference of intergated tracer flux passing through the upper and lower interface of layer $l$ (Fick's law). Here we take the generic form: the diffusivity may vary at layer interface. Of course, we can recover the expression in [http://basilisk.fr/src/layered/diffusion.h]() by setting $D = D_{l + 1 / 2} = D_{l - 1 / 2}$.


We recall that the dimension the diffusivity "$D$" is $L^2 T^{-1}$, being independent to the physcial problem. Here we study the transportation of quantity momentum, then $s = \rho_l u_l$, where $u_l$ is the averaged longitudinal velocity in layer $l$. Then the above equation becomes
$$
\dfrac{(\rho_l h_l u_l)^{n + 1} - (\rho_l h_l u_l)^{\star}}{\Delta t} =
\left( \tau_{l+1/2} - \tau_{l-1/2} \right)^{n + 1},
$$
where $\tau_{l+1/2}$ and $\tau_{l-1/2}$ are shear stress applying on the upper and lower interface of layer $l$, and 
$$
\tau_{l+1/2}=\rho_{1+1/2}D_{l + 1 / 2} \dfrac{u_{l + 1} - u_l}{h_{l + 1 / 2}},\quad
\tau_{l-1/2}=\rho_{1-1/2}D_{l - 1 / 2} \dfrac{u_{l} - u_{l-1}}{h_{l - 1 / 2}}.
$$
We may "immediately" notice that $D$ represents the cinematic viscosity. If we consider an imcompressible fluide ($\rho=1$), then the value of $D$ is the same as the dynamic viscosity.
An example of "two phase layered model" can be seen in [https://www.basilisk.fr/sandbox/ryang/multilayer_varydensity/my_diffusion.h](). As we consider here a single-phase probleme, the fluide density is constant and is set to one. 

In complex flows, the viscosity can be expressed in function of local quantities, such as the position, shear, pressure, and other rheological parameters, i.e.
$$
D=\nu_{eq}(y,\frac{\partial u}{\partial z},P,...)
$$
Idealy, the function $\nu_{eq}$ should take the values defined at the layer interface to reflect the flux going through at the interface. But here we take a different strategy: all the quantites are defined at the center of each layer, the avergerd layer velocity is also treated as pointwise velocity locating at the layer center. As a result, the $\nu_{eq}$ is defined at the layer center. When the value of $\nu_{eq}$ at an interface is required, e.g. $D_{l + 1 / 2}$ and $D_{l - 1 / 2}$, it is approximated by averaging the values defined in the two adjacent layers, such as
$$
D_{l + 1 / 2} = \dfrac{1}{2}(D_l + D_{l+1}),\quad D_{l - 1 / 2} = \dfrac{1}{2}(D_l + D_{l-1})
$$
With such convention, the values of $\nu_{eq}$ is storaged in a table `nu_eq[nl]`, whose dimension is same as the layer number $nl$*/



/** # The vertical diffusion function */
/**
For stability, we discretise the vertical diffusion equation implicitly as
$$
\dfrac{(h_l s_l)^{n + 1} - (h_l s_l)^{\star}}{\Delta t} =
\left( D_{l + 1 / 2} \dfrac{s_{l + 1} - s_l}{h_{l + 1 / 2}} -
D_{l - 1 / 2}\dfrac{s_l - s_{l - 1}}{h_{l - 1 / 2}} \right)^{n + 1}
$$
Reorganize it, we get
$$
\underbrace{\left( h_l^{n+1} + \dfrac{D_{l+1/2}\Delta t}{h_{l+1/2}}+\dfrac{D_{l-1/2}\Delta t}{h_{l-1/2}} \right)}_{b_l}s_l^{n+1} + \underbrace{\left( -\dfrac{D_{l+1/2}\Delta t}{h_{l+1/2}}\right)}_{c_l}s_{l+1}^{n+1}+\underbrace{\left( -\dfrac{D_{l-1/2}\Delta t}{h_{l-1/2}}\right)}_{a_l}s_{l-1}^{n+1}=(h_l s_l)^\star
$$


which can be expressed as the linear system
$$
\mathbf{Ms}^{n + 1} = \mathrm{rhs}
$$
where $\mathbf{M}$ is a 
[tridiagonal matrix](https://en.wikipedia.org/wiki/Tridiagonal_matrix). 
The lower, principal and upper diagonals are *a*, *b* and *c* respectively. 

Boundary conditions on the top and bottom layers need to be added to close the
system. We chose to impose a "generalized boundary condition"
$$
s\vert_{b,t} = s_{b,t} + \lambda_{b,t}(\partial_z s\vert_{b,t}-q_{b,t}).
$$
According to the value of $\lambda_{b,t}$, the boundary condition can be altered to a Dirichlet, Neumann, or Navier condition.
\begin{aligned}
\lambda_{b,t}=0 \to s\vert_{b,t} = s_{b,t}\quad \text{(Dirichlet condition)}\\
\lambda_{b,t}=\infty \to \partial_z s\vert_{b,t}=q_{b,t} \quad \text{(Neumann condition)}\\
\lambda_{b,t}\in[0, \infty], q_{b,t}=0 \to s\vert_{b,t} = s_{b,t} + \lambda_{b,t}\partial_z s\vert_{b,t} \quad \text{(Navier condition)}
\end{aligned}


*/

#if  RHEOLOGY
  #include "rheology.h"

#else
	 double nu = 0
#endif 


void vertical_diffusion (Point point, scalar h, scalar s, double dt, double D,
			 double dst, double dsb,double s_b, double lambda_b, double s_t, double lambda_t)
{
  double a[nl], b[nl], c[nl], rhs[nl];
  double nueq[nl];

  /**
  The *rhs* of the tridiagonal system is $h_l s_l$. */
      
  foreach_layer()
    rhs[point.l] = s[]*h[];

  //fprintf (stdout, "diffusion check %g %g %g %g %g %g %g\n", t, x, y, dt, h[], s[], rhs[point.l] );

  foreach_layer()
    rhs[point.l] = s[]*h[];


/** Regularisation for the nueq */
  for (int l = 1; l < nl; l++) nueq[l] = D;

#if RHEOLOGY
  regularization(point, s, h, nueq, D);
#endif

  /**
  The lower, principal and upper diagonals $a$, $b$ and $c$ are given by
  $$
  a_{l > 0} = - \left( \dfrac{D_{l - 1 / 2} \Delta t}{h_{l - 1 / 2}} \right)^{n + 1} = -\Delta t\dfrac{D_l + D_{l-1}}{h_l + h_{l-1}}
  $$
  $$
  c_{l < \mathrm{nl} - 1} = - \left( \dfrac{D_{l + 1 / 2} \Delta t}{h_{l + 1 / 2}}
  \right)^{n + 1}=-\Delta t\dfrac{D_l + D_{l+1}}{h_l + h_{l+1}}
  $$
  $$
  b_{0 < l < \mathrm{nl} - 1} = h_l^{n + 1} - a_l - c_l
  $$
  */

  for (int l = 1; l < nl - 1; l++) {
    a[l] = - (nueq[l]+nueq[l-1])*dt/(h[0,0,l-1] + h[0,0,l]);
    c[l] = - (nueq[l]+nueq[l+1])*dt/(h[0,0,l] + h[0,0,l+1]);
    b[l] = h[0,0,l] - a[l] - c[l];
    //fprintf (stdout, "diffusion check1 %g %g %g %d %g %g %g\n", t, x, y, l, a[l], b[l], c[l] );
  }
    
/** ## The boundary conditions*/
/**
Two additional points should be adressed. 

For bottom and top layers, the $D_{-1/2}$ and $D_{nl-1/2}$ can not be obtained by the "average" convention introduced before, as $D_{-1}$ and $D_{nl}$ do not exist. Nevertheless, the $D_{-1/2}$ and $D_{nl-1/2}$ have clear physcial sens, they represent the viscosities on the bottom and free surface boundaries, which may be complex to model.  We introduce two additional variables $D_b$ and $D_t$ to represent the viscosity at boundaries and we make the approximations :
$$
D_{-1/2} = D_b = \nu_{eq}[0], \quad D_{nl-1/2} = D_t = \nu_{eq}[nl-1]
$$
Another issue arises from the fact that the value of $s_{-1}$ and $s_{nl}$ are undefined. We suppose the ghost layers, where the values of $s$ satisfie the imposed boundary conditions. The detailed developments shoud be presented in forthcoming publications. Here, we outline the main ideas.

Take the bottom as example. We assume that the profile of $s$ is parabolic near the bottom, which is equivalent to a third-order Taylor expansion at the position y
$$
s(y) = P_0 + P_1 y +P_2 y^2 (+O(y^3))
$$
then the integrations of this profil over the layer give values of $s_{0}$ and $s_1$. meanwhile the profil satisfies bottom boudary condition $s\vert_b$. These yield
\begin{aligned}
s_0 &= \frac{1}{h_0}\int_0^{h_0} s(y)dy = P_0 + \frac{1}{2}h_0 P_1+ \frac{1}{3} h^2_0 P_2\\
s_1 &= \frac{1}{h_1}\int_{h_0}^{h_1+h_0} s(y)dy = P_0 + h_0 P_1 + h^2_0 P_2 +\frac{h_1}{2}P_1 + h_0h_1 P_2 + \frac{h^3_1}{3}P_2\\
s\vert_b &= s(0) = s_{b} + \lambda_{b}(\partial_y s\vert_{b}-q_{b}) \to P_0 - \lambda_b P_1 = s_b-\lambda_b q_b
\end{aligned}
write in matrice form
$$
\begin{bmatrix}
1 & \frac{h_0}{2} & \frac{h^2_0}{3}\\
1 & h_0+\frac{h_1}{2} & h^2_0+h_0 h_1+\frac{h^2_1}{3}\\
1 & -\lambda_b & 0\\
\end{bmatrix}
\begin{bmatrix}
P_0 \\
P_1 \\
P_2 \\
\end{bmatrix}
=
\begin{bmatrix}
s_0 \\
s_1 \\
s_b - \lambda_b q_b \\
\end{bmatrix}
$$
The coefficient $P_0$, $P_1$ and $P_2$ can be resolved, for example, with [Règle de Cramer](https://fr.wikipedia.org/wiki/Règle_de_Cramer). Then it can be readily noticed that the solutions for $P_0$, $P_1$, and $P_2$ share the same denominator, which is defined by the variable `det` (more or less..), meaning determinant of the $3\times 3$ matrice.

Since the profil of $s(y)$ is fully resolved with the knowledge of $P_0$, $P_1$ and $P_2$, we are able to compute the bottom ghost layer value 
$$
s_{-1} = \frac{1}{h_0}\int_{-h_0}^0 s(y)dy = P_0 - \frac{h_0}{2}P_1 + \frac{h_0}{3}P_2.
$$
Inserting the expression of $s_{-1}$ into the Equation. The coefficients in front of $s_0$ and $s_{1}$ give $b_0$ and $c_0$ and according $rhs_0$.
*/
  /**
  $$
  \begin{aligned}
  b_0 & = h_0 + \Delta t\left( \frac{D_0+D_1}{h_0 + h_1} + 2D_b\frac{h^2_1 + 3
  h_0 h_1 + 3 h^2_0}{\det} \right),\\
  c_0 & = - \Delta t \left( \frac{D_0+D_1}{h_0 + h_1} + 2D_b\frac{h^2_0}{\det}
  \right),\\
  \text{rhs}_0 & = (hs_0)^{\star} + 2 \Delta t D_b (s_b-\lambda_b q_b)  \frac{h^2_1 + 3 h_0
  h_1 + 2 h^2_0}{\det},\\
  \det & = h_0 (h_0 + h_1)^2  + 2\lambda _b(3\,h_0 h_1 + 2\,h_0^2 + h_1^2),
  \end{aligned}
  $$
  */
//
D_b = nueq[0];
   den = h[]*sq(h[] + h[0,0,1]) 
    + 2.*lambda_b*(3.*h[]*h[0,0,1] + 2.*sq(h[]) + sq(h[0,0,1]));
  b[0] = h[] + dt*((nueq[1]+nueq[0])/(h[] + h[0,0,1]) +
			 2*D_b*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 3.*sq(h[]))/den);
  c[0] = -dt*((nueq[0]+nueq[1])/(h[] + h[0,0,1]) + 2.*sq(h[])*D_b/den);

  rhs[0] += 2.*dt*D_b*(s_b-lambda_b*dsb)*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 2.*sq(h[0]))/den;

  if (nl == 1) {
    b[0] += c[0];
    rhs[0] +=  (- c[0]*h[] - nueq[0]*dt) * dst;
  }

/** 
For the top layer, we apply the same procedure. That yields
$$
  \begin{aligned}
  b_{nl-1} & = h_{nl-1} + \Delta t\left( \frac{D_{nl-1}+D_{nl-2}}{h_{nl-1} + h_{nl-2}} + 2D_t\frac{h^2_{nl-2} + 3
  h_{nl-1} h_{nl-2} + 3 h^2_{nl-1}}{\det} \right),\\
  a_{nl-1} & = - \Delta t \left( \frac{D_{nl-1}+D_{nl-2}}{h_{nl-1} + h_{nl-2}} + 2D_t\frac{h^2_{nl-1}}{\det}
  \right),\\
  \text{rhs}_{nl-1} & = (hs_{nl-1})^{\star} + 2 \Delta t D_t (s_t-\lambda_t q_t)  \frac{h^2_{nl-2} + 3 h_{nl-2}
  h_{nl-1} + 2 h^2_{nl-1}}{\det},\\
  \det & = h_{nl-1} (h_{nl-1} + h_{nl-2})^2  - 2\lambda _t(3\,h_{nl-1} h_{nl-2} + 2\,h_{nl-1}^2 + h_{nl-2}^2),
  \end{aligned}
  $$
*/
D_t = nueq[nl-1];
  double den = h[0,0,nl-1]*sq(h[0,0,nl-1] + h[0,0,nl-2]) 
    - 2.*lambda_t*(3.*h[0,0,nl-1]*h[0,0,nl-2] + 2.*sq(h[0,0,nl-1]) + sq(h[0,0,nl-2]));
  b[nl-1] = h[0,0,nl-1] + dt*((nueq[nl-1]+nueq[nl-2])/(h[0,0,nl-1] + h[0,0,nl-2]) +
			 2*D_t*(sq(h[0,0,nl-2]) + 3.*h[0,0,nl-1]*h[0,0,nl-2] + 3.*sq(h[0,0,nl-1]))/den);
  a[nl-1] = - dt*((nueq[nl-1]+nueq[nl-2])/(h[0,0,nl-1] + h[0,0,nl-2]) + 2*D_t*sq(h[0,0,nl-1])/den);

  rhs[nl-1] += 2.*dt*D_t*(s_t-lambda_t*dst)*(sq(h[0,0,nl-2]) + 3.*h[0,0,nl-2]*h[0,0,nl-1] + 2.*sq(h[0,0,nl-1]))/den;
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
 */



/**
In the [layered solver](hydro.h), vertical viscosity is applied to the
velocity field just after advection, but before the pressure
gradient/acceleration term is applied. To take the pressure gradient
into account, we first apply the acceleration of the previous
timestep, apply vertical viscosity and then substract the previous
acceleration. */
/*
event viscous_term (i++,last)
{
    if (nu > 0.) {
        foreach() {
            foreach_layer()
            foreach_dimension()
            u.x[] += dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
            foreach_dimension()
					vertical_diffusion (point, h, u.x, dt, nu, dut.x[], u_b.x[], lambda_b.x[]);
            foreach_layer()
            foreach_dimension()
            u.x[] -= dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
        }
        boundary ((scalar *) {u});
    }
}
*/

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
			    dut.x[],dub.x[],u_b.x[], lambda_b.x[],u_t.x[],lambda_t.x[]);
      foreach_layer()
	foreach_dimension()
	  u.x[] -= dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
    }
    }
  boundary ((scalar *) {u});
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
/*
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
*/
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

//#include "poisson.h"
/*
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
*/
/**
## References

~~~bib
@hal{popinet2020, hal-02365730}

@hal{devita2019, hal-02295398}
~~~
*/

