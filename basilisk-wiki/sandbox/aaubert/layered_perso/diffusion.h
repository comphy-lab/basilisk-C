/**
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
and a Navier slip condition on the bottom i.e.
$$
s|_b = s_b + \lambda_b \partial_z s|_b
$$ */

void vertical_diffusion (Point point,double x1, scalar h, scalar s, double dt,
			 face vector D, double dst, double s_b, double lambda_b,
			 bool horizontal,bool up, double s_t,double molvis2)
{
  
  double a[nl], b[nl], c[nl], rhs[nl];
  double d0=0.;
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
  double unnn=1.;   //fixme:needed for dimensions
  for (int l = 1; l < nl - 1; l++) {
    a[l] = - 2.*unnn*(D.x[0,0,l-1]+D.x[0,0,l])/2.*dt/(h[0,0,l-1] + h[0,0,l]);
    c[l] = - 2.*unnn*(D.x[0,0,l]+D.x[0,0,l+1])/2.*dt/(h[0,0,l] + h[0,0,l+1]);
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

  if (up) {
    a[nl-1] = - 2.*unnn*(D.x[0,0,nl-2]+D.x[0,0,nl-1])/2.*dt/(h[0,0,nl-2] + h[0,0,nl-1]);
    b[nl-1] = h[0,0,nl-1] - a[nl-1]+2.*unnn*D.x[0,0,nl-1]*dt/h[0,0,nl-1];
    rhs[nl-1]+=unnn*D.x[0,0,nl-1]*2.*s_t*dt/h[0,0,nl-1];
  }
  else {
    a[nl-1] =  - 2.*unnn*(D.x[0,0,nl-2]+D.x[0,0,nl-1])/2.*dt/(h[0,0,nl-2] + h[0,0,nl-1]);
    b[nl-1] = h[0,0,nl-1] - a[nl-1];
    rhs[nl-1] += unnn*D.x[0,0,nl-1]*dt*dst;
  }

  /**
  For the bottom layer a third-order discretisation of the Navier slip
  condition gives
  $$
  \begin{aligned}
  b_0 & = h_0 + 2 \Delta t D \left( \frac{1}{h_0 + h_1} + \frac{h^2_1 + 3
  h_0 h_1 + 3 h^2_0}{\det} \right),\\
  c_0 & = - 2 \Delta t D \left( \frac{1}{h_0 + h_1} + \frac{h^2_0}{\det}
  \right),\\
  \text{rhs}_0 & = (hs_0)^{\star} + 2 \Delta t D s_b  \frac{h^2_1 + 3 h_0
  h_1 + 2 h^2_0}{\det},\\
  \det & = h_0 (h_0 + h_1)^2  + 2\lambda (3\,h_0 h_1 + 2\,h_0^2 + h_1^2),
  \end{aligned}
  $$
  */
  double Dm1=x1<1.?D.x[0,0,0]:2.*molvis2-D.x[0,0,0];
  if ((horizontal)&&(x1<1.)) {   //neumann
    c[0] = - 2.*unnn*(D.x[0,0,0]+D.x[0,0,1])/2.*dt/(h[0,0,0] + h[0,0,1]);
    b[0] = h[0,0,0] - c[0];
    double zerooooo=0.;
    rhs[0] += D.x[0,0]*dt*zerooooo;
    //lambda_b=HUGE;
    //fourth order
    //double det=1/6.*pow(h[0,0,0]+h[0,0,1]+h[0,0,2]/2.,3.)-sq(h[0,0,0])/24.*(h[0,0,0]+h[0,0,1]+h[0,0,2]/2.)-1./3.*(h[0,0,0]+h[0,0,1]/2.)*(1./2.*sq(h[0,0,0]+h[0,0,1]+h[0,0,2]/2.)-sq(h[0,0,0])/8.);
    //double coef2=-sq(h[0,0,0])/24./det;
    //double coef1=-(1./2.*sq(h[0,0,0]+h[0,0,1]+h[0,0,2]/2.)-sq(h[0,0,0])/8.)*coef2/(1./2.*sq(h[0,0,0]+h[0,0,1]/2.)-sq(h[0,0,0])/8.);
    //double coefm1=(1-((h[0,0,0]/2.+h[0,0,1]/2.)*coef1+(h[0,0,0]/2.+h[0,0,1]+h[0,0,2]/2.)*coef2))/(-h[0,0,0]);
    //double coef0=-coefm1-coef1-coef2;
    //c[0]=-2.*unnn*(D.x[0,0,0]+D.x[0,0,1])/2.*dt/(h[0,0,0] + h[0,0,1])+dt*unnn*D.x[0,0,0]/h[0,0,0]*coef1/coefm1;
    //d0=unnn*D.x[0,0,0]*dt/h[0,0,0]*coef2/coefm1;
    //b[0]=h[0,0,0]+2.*unnn*(D.x[0,0,0]+D.x[0,0,1])/2.*dt/(h[0,0,0] + h[0,0,1])+unnn*D.x[0,0,0]*dt/h[0,0,0]*(1.+coef0/coefm1);
  }
  else {   //dirichlet
    //double den = h[]*sq(h[] + h[0,0,1]) 
    //  + 2.*lambda_b*(3.*h[]*h[0,0,1] + 2.*sq(h[]) + sq(h[0,0,1]));
    //b[0] = h[] + 2.*dt*unnn*(D.x[0,0,0]+D.x[0,0,1])/2.*(1./(h[] + h[0,0,1]) +
    //							(D.x[0,0,0]+Dm1)/2.*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 3.*sq(h[]))/den);
    //c[0] = - 2.*dt*unnn*((D.x[0,0,0]+D.x[0,0,1])/2./(h[] + h[0,0,1]) + (D.x[0,0,0]+Dm1)/2.*sq(h[])/den);
    //rhs[0] += 2.*dt*unnn*(D.x[0,0,0]+Dm1)/2.*s_b*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 2.*sq(h[]))/den;
    //my version
    double den=h[]*(h[]+h[0,0,1])*(h[]+h[0,0,1]/2.)+lambda_b*(3.*sq(h[])+4.*h[]*h[0,0,1]+sq(h[0,0,1]));
    b[0]=h[]+2*dt*unnn*(D.x[0,0,0]+D.x[0,0,1])/2./(h[0,0,0]+h[0,0,1])+dt*unnn*(D.x[0,0,0]+Dm1)/2.*(4.*sq(h[0,0,0])+4.*h[0,0,0]*h[0,0,1]+sq(h[0,0,1]))/den;
    c[0]=-2.*dt*unnn*(D.x[0,0,0]+D.x[0,0,1])/2.*(1./(h[] + h[0,0,1]))-dt*unnn*(D.x[0,0,0]+Dm1)/2.*sq(h[0,0,0])/den;
    rhs[0]+=dt*unnn*(D.x[0,0,0]+Dm1)/2.*s_b*(3.*sq(h[0,0,0])+4.*h[0,0,0]*h[0,0,1]+sq(h[0,0,1]))/den;
    //second order
    //b[0]=h[0,0]+2*dt*unnn*(D.x[0,0]+D.x[0,1])/2./(h[0,0]+h[0,1])+dt*molvis2/h[0,0]*2.;
    //c[0]=-2.*dt*unnn*(D.x[0,0]+D.x[0,1])/2.*(1./(h[0,0] + h[0,1]));
    //rhs[0]+=dt*molvis2*s_b*2./h[0,0];
  }

  if (nl == 1) {
    b[0] += c[0];
    rhs[0] += (- c[0]*h[] -unnn* D.x[0,0,0]*dt) * dst;
  }

#if 1
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
#else
  /**
     We will solve the pentadiagonal system using a generalization of the Thomas algorithm (see Numerical Metohds for Grid equations, Nauka, 1978) */

  double alpha[nl], gamma[nl+1];
  double beta_1,delta;

  for (int l=0;l<nl;l++) {
    if (l==0) {
      alpha[1]=-c[0]/b[0];
      beta_1=d0/b[0];
      gamma[1]=rhs[0]/b[0];
    }
    else if (l==1) {
      delta=b[1]+a[1]*alpha[1];
      alpha[2]=(-c[1]+beta_1*a[1])/delta;
      gamma[2]=(rhs[1]-a[1]*gamma[1])/delta;
    }
    else if (l<nl-1) {
      delta=b[l]+alpha[l]*a[l];
      alpha[l+1]=-c[l]/delta;
      gamma[l+1]=(rhs[l]-gamma[l]*a[l])/delta;
    }
    else if (l==nl-1) {
      delta=b[nl-1]+alpha[nl-1]*a[nl-1];
      gamma[nl]=(rhs[nl-1]-gamma[nl-1]*a[nl-1])/delta;
    }
  }
  s[0,0,nl-1]=gamma[nl];
  for (int l=nl-2;l>0;l--) {
    s[0,0,l]=alpha[l+1]*s[0,0,l+1]+gamma[l+1];
  }
  s[0,0,0]=alpha[1]*s[0,0,1]-beta_1*s[0,0,2]+gamma[1];
#endif  
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

void horizontal_diffusion (scalar * list, face vector D, double dt)
{
  scalar * d2sl = list_clone (list);
  foreach() {
    scalar s, d2s;
    for (s,d2s in list,d2sl) {
      foreach_layer() {
	double a = 0.;
	double unnn=1.;   //fixme:needed for dimensions
	foreach_dimension() {
	  a += (unnn*D.x[]*hf.x[]*fm.x[]/(cm[-1] + cm[])*(s[-1] - s[]) + 
		unnn*D.x[1]*hf.x[1]*fm.x[1]/(cm[1] + cm[])*(s[1] - s[]));
	}
	d2s[] = 2.*a/(cm[]*sq(Delta));
      }
    }
  }
  foreach() {
    if (h[] > dry) {
      scalar s, d2s;
      for (s,d2s in list,d2sl) {
	foreach_layer() {
	  s[] += dt*d2s[]/h[];
	}
      }
    }
  }
  delete (d2sl);
  free (d2sl);
}

/**
# Viscous friction between layers

By default the viscosity is zero and we impose free-slip on the
free-surface and no-slip on the bottom boundary
i.e. $\dot{\mathbf{u}}_t = 0$, $\mathbf{\lambda}_b = 0$, $\mathbf{u}_b
= 0$. */

face vector nu;
event defaults0(i=0) {
  nu=new face vector[nl];
  foreach_layer() {
    foreach_face() {
      dimensional (nu==Delta*Delta/t);
    }
  }
}

/*
fixme: should just be:
(const) vector lambda_b[] = {0,0,0}, dut[] = {0,0,0}, u_b[] = {0,0,0};
*/
const vector lambda0[] = {0,0,0}, dut0[] = {0,0,0}, u_b0[] = {0,0,0};
(const) vector lambda_b = lambda0, dut = dut0, u_b = u_b0;

/**
In the [layered solver](hydro.h), vertical viscosity is applied to the
velocity field just after advection, but before the pressure
gradient/acceleration term is applied. To take the pressure gradient
into account, we first apply the acceleration of the previous
timestep, apply vertical viscosity and then substract the previous
acceleration. */

event viscous_term (i++,last)
{
  
  foreach() {
    foreach_layer(){ 
      foreach_dimension() {
	u.x[] += dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
      }
    }
    foreach_dimension() {
      vertical_diffusion (point,x, h, u.x, dt, nu,
      			  dut.x[], u_b.x[], lambda_b.x[],true,false,1.,1/5e6);
    }
  }
  
  foreach() {	  
    foreach_layer() {
      foreach_dimension() {
	u.x[] -= dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
      }
    }
  }
}

event cleanup(t=end,last) {
  delete ((scalar *) {nu});
}


/**
## References

~~~bib
@hal{popinet2020, hal-02365730}

@hal{devita2019, hal-02295398}
~~~
*/