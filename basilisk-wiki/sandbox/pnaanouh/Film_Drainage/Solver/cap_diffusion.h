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
system. Any combination of Neumann and Navier slip conditons can be applied at the top any bottom boundaries. This is implemented using a seperate function for each combination but they should eventually be combined into one modular function.

The h_diffusion flag activates the horizontal diffusion
*/
bool h_diffusion = false;

/**
Uses a neumann condition at both the top and the bottom boundary
$$
\partial_z s |_t = \dot{s}_t \\
\partial_z s |_b = \dot{s}_b
$$
*/
void vertical_diffusion_NeumannNeumann (Point point, scalar h, scalar s, double dt, double D,
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
  For the bottom layer the boundary conditions give the (ghost)
  boundary value
  $$
  s_{-1} = s_{0} - \dot{s}_b h_{0},
  $$
  which gives the diagonal coefficient and right-hand-side
  $$
  b_{0} = h_{0}^{n + 1}
  - c_{0}
  $$
  $$
  \mathrm{rhs}_{0} = 
  (hs)_{0}^{\star} - D \Delta t \dot{s}_b
  $$
  */

  c[0]=- 2.*D*dt/(h[0,0,0] + h[0,0,1]);
  b[0]= h[] - c[0];
  rhs[0] -= D*dt*dsb;

  if (nl == 1) {
    b[0] = h[];
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
Uses a Neumann condition on the top and a Navier slip condition on the bottom
$$
\partial_z s |_t = \dot{s}_t \\
s|_b = s_b + \lambda_b \partial_z s|_b
$$
*/

void vertical_diffusion_NeumannNavier (Point point, scalar h, scalar s, double dt, double D,
				        double dst, double s_b, double lambda_b)
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
  \det & = h_0 (h_0 + h_1)^2  + 2\lambda_b (3\,h_0 h_1 + 2\,h_0^2 + h_1^2),
  \end{aligned}
  $$
  */

  double den = h[]*sq(h[] + h[0,0,1]) 
    + 2.*lambda_b*(3.*h[]*h[0,0,1] + 2.*sq(h[]) + sq(h[0,0,1]));
  b[0] = h[] + 2.*dt*D*(1./(h[] + h[0,0,1]) +
			  (sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 3.*sq(h[]))/den);
  c[0] = - 2.*dt*D*(1./(h[] + h[0,0,1]) + sq(h[])/den);
  rhs[0] += 2.*dt*D*s_b*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 2.*sq(h[0]))/den;
  

  if (nl == 1) {
    b[0] += c[0];
    rhs[0] += (- c[0]*h[] - D*dt) * dst;
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
Uses a Navier slip condition on the top and a Neumann condition on the bottom
$$
s|_t = s_t + \lambda_t \partial_z s|_t \\
\partial_z s |_b = \dot{s}_b
$$
*/

void vertical_diffusion_NavierNeumann (Point point, scalar h, scalar s, double dt, double D,
				       double s_t, double lambda_t, double dsb)
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
  For the top layer a third-order discretisation of the Navier slip
  condition gives
  $$
  \begin{aligned}
  a_{\mathrm{nl} - 1} & = - 2 \Delta t D \left( \frac{1}{h_{\mathrm{nl} - 1} + h_{\mathrm{nl} - 2}} + \frac{h^2_{\mathrm{nl} - 1}}{\det}
  \right),\\
  b_{\mathrm{nl} - 1}& = h_{\mathrm{nl} - 1} + 2 \Delta t D \left( \frac{1}{h_{\mathrm{nl} - 1} + h_{\mathrm{nl} - 2}} + \frac{h^2_{\mathrm{nl} - 2} + 3
  h_{\mathrm{nl} - 1} h_ {\mathrm{nl} - 2}+ 3 h^2_{\mathrm{nl} - 1}}{\det} \right),\\
  \text{rhs}_{\mathrm{nl} - 1} & = (hs_{\mathrm{nl} - 1})^{\star} + 2 \Delta t D s_t\frac{h^2_{\mathrm{nl} - 2} + 3 h_{\mathrm{nl} - 1}
  h_{\mathrm{nl} - 2} + 2 h^2_{\mathrm{nl} - 1}}{\det},\\
  \det & = h_{\mathrm{nl} - 1} (h_{\mathrm{nl} - 1} + h_{\mathrm{nl} - 2})^2  - 2\lambda_t (3\,h_ {\mathrm{nl} - 1}h_{\mathrm{nl} - 2} + 2\,h_{\mathrm{nl} - 1}^2 + h_{\mathrm{nl} - 2}^2),
  \end{aligned}
  $$
  */
  
  double den = h[0,0,nl-1]*sq(h[0,0,nl-1] + h[0,0,nl-2]) 
    - 2.*lambda_t*(3.*h[0,0,nl-1]*h[0,0,nl-2] + 2.*sq(h[0,0,nl-1]) + sq(h[0,0,nl-2]));
  
  a[nl-1] = - 2.*dt*D*(1./(h[0,0,nl-1] + h[0,0,nl-2]) + sq(h[0,0,nl-1])/den);
  b[nl-1] = h[0,0,nl-1] + 2.*dt*D*(1./(h[0,0,nl-1] + h[0,0,nl-2]) +
			  (sq(h[0,0,nl-2]) + 3.*h[0,0,nl-1]*h[0,0,nl-2] + 3.*sq(h[0,0,nl-1]))/den);
  rhs[nl-1] += 2.*dt*D*s_t*(sq(h[0,0,nl-2]) + 3.*h[0,0,nl-1]*h[0,0,nl-2] + 2.*sq(h[0,0,nl-1]))/den;
  
   /**
  For the bottom layer the boundary conditions give the (ghost)
  boundary value
  $$
  s_{-1} = s_{0} - \dot{s}_b h_{0},
  $$
  which gives the diagonal coefficient and right-hand-side
  $$
  b_{0} = h_{0}^{n + 1}
  - c_{0}
  $$
  $$
  \mathrm{rhs}_{0} = 
  (hs)_{0}^{\star} - D \Delta t \dot{s}_b
  $$
  */

  c[0]=- 2.*D*dt/(h[0,0,0] + h[0,0,1]);
  b[0]= h[] - c[0];
  rhs[0] -= D*dt*dsb;

  if (nl == 1) {
    b[0] = h[]+3.*dt*D/(h[]-3.*lambda_t);
    rhs[0] = s[]*h[]+3.*dt*D*(s_t/(h[]-3.*lambda_t)-1/2.*(h[]-2.*lambda_t)/(h[]-3.*lambda_t)*dsb); 
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
Uses a Navier slip condition on both the top and bottom
$$
s|_t = s_t + \lambda_t \partial_z s|_t \\
s|_b = s_b + \lambda_b \partial_z s|_b
$$
*/

void vertical_diffusion_NavierNavier (Point point, scalar h, scalar s, double dt, double D,
				      double s_t, double lambda_t, double s_b, double lambda_b)
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
  For the top layer a third-order discretisation of the Navier slip
  condition gives
  $$
  \begin{aligned}
  a_{\mathrm{nl} - 1} & = - 2 \Delta t D \left( \frac{1}{h_{\mathrm{nl} - 1} + h_{\mathrm{nl} - 2}} + \frac{h^2_{\mathrm{nl} - 1}}{\det}
  \right),\\
  b_{\mathrm{nl} - 1}& = h_{\mathrm{nl} - 1} + 2 \Delta t D \left( \frac{1}{h_{\mathrm{nl} - 1} + h_{\mathrm{nl} - 2}} + \frac{h^2_{\mathrm{nl} - 2} + 3
  h_{\mathrm{nl} - 1} h_ {\mathrm{nl} - 2}+ 3 h^2_{\mathrm{nl} - 1}}{\det} \right),\\
  \text{rhs}_{\mathrm{nl} - 1} & = (hs_{\mathrm{nl} - 1})^{\star} + 2 \Delta t D s_t\frac{h^2_{\mathrm{nl} - 2} + 3 h_{\mathrm{nl} - 1}
  h_{\mathrm{nl} - 2} + 2 h^2_{\mathrm{nl} - 1}}{\det},\\
  \det & = h_{\mathrm{nl} - 1} (h_{\mathrm{nl} - 1} + h_{\mathrm{nl} - 2})^2  - 2\lambda_t (3\,h_ {\mathrm{nl} - 1}h_{\mathrm{nl} - 2} + 2\,h_{\mathrm{nl} - 1}^2 + h_{\mathrm{nl} - 2}^2),
  \end{aligned}
  $$
  */

  double den_t = h[0,0,nl-1]*sq(h[0,0,nl-1] + h[0,0,nl-2]) 
    - 2.*lambda_t*(3.*h[0,0,nl-1]*h[0,0,nl-2] + 2.*sq(h[0,0,nl-1]) + sq(h[0,0,nl-2]));
  
  a[nl-1] = - 2.*dt*D*(1./(h[0,0,nl-1] + h[0,0,nl-2]) + sq(h[0,0,nl-1])/den_t);
  b[nl-1] = h[0,0,nl-1] + 2.*dt*D*(1./(h[0,0,nl-1] + h[0,0,nl-2]) +
			  (sq(h[0,0,nl-2]) + 3.*h[0,0,nl-1]*h[0,0,nl-2] + 3.*sq(h[0,0,nl-1]))/den_t);
  rhs[nl-1] += 2.*dt*D*s_t*(sq(h[0,0,nl-2]) + 3.*h[0,0,nl-1]*h[0,0,nl-2] + 2.*sq(h[0,0,nl-1]))/den_t;

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
  \det & = h_0 (h_0 + h_1)^2  + 2\lambda_b (3\,h_0 h_1 + 2\,h_0^2 + h_1^2),
  \end{aligned}
  $$
  */

  double den_b = h[]*sq(h[] + h[0,0,1]) 
    + 2.*lambda_b*(3.*h[]*h[0,0,1] + 2.*sq(h[]) + sq(h[0,0,1]));
  b[0] = h[] + 2.*dt*D*(1./(h[] + h[0,0,1]) +
			  (sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 3.*sq(h[]))/den_b);
  c[0] = - 2.*dt*D*(1./(h[] + h[0,0,1]) + sq(h[])/den_b);
  rhs[0] += 2.*dt*D*s_b*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 2.*sq(h[0]))/den_b;
  
  if (nl == 1) {
    b[0] = h[]+12.*D*dt*(h[]+lambda_b-lambda_t)/(sq(h[])+4.*h[]*(lambda_b-lambda_t)-12.*lambda_b*lambda_t);
    rhs[0] =  s[]*h[] + 6.*D*dt * ((h[]+2.*lambda_b)*s_t+(h[]-2.*lambda_t)*s_b)/(sq(h[])+4.*h[]*(lambda_b-lambda_t)-12.*lambda_b*lambda_t);
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
## Horizontal diffusion

This approximates
$$
h \partial_t s = D \nabla \cdot ( \nabla hs)
$$
with $D$ the diffusion coefficient. Note that the
time discretisation is explicit so that the timestep must be limited
(manually) by $\min(\Delta^2/D)$. */

void horizontal_diffusion_NeumannNeumann(scalar s, double D, double dt, (const) scalar dst, (const) scalar dsb)
{
  if (D > 0.) {
    
    scalar d2s= new scalar [nl];
    face vector dzl[];
    foreach_face(){
      dzl.x[] = (zb[]-zb[-1])/Delta;
    }

    // $\nabla^2(hu)_l$
    foreach(){
      foreach_layer(){
	foreach_dimension(){
	  d2s[] = (h[1]*s[1]-2.*h[]*s[]+h[-1]*s[-1])/sq(Delta);
	}
      }
    }

    // Inter_layer terms at first layer bottom
    foreach(){
      foreach_dimension(){
	if (nl>1){
	  d2s[0,0,0] += ((((3.*sq(h[1,0,0])+3.*h[1,0,0]*h[1,0,1]+sq(h[1,0,1]))*s[1,0,0]-sq(h[1,0,0])*s[1,0,1])/((h[1,0,0]+h[1,0,1])*(2.*h[1,0,0]+h[1,0,1]))-(1./2.)*dsb[1,0]*h[1,0,0]*(h[1,0,0]+h[1,0,1])/(2.*h[1,0,0]+h[1,0,1])) * dzl.x[1]  -  (((3.*sq(h[-1,0,0])+3.*h[-1,0,0]*h[-1,0,1]+sq(h[-1,0,1]))*s[-1,0,0]-sq(h[-1,0,0])*s[-1,0,1])/((h[-1,0,0]+h[-1,0,1])*(2.*h[-1,0,0]+h[-1,0,1]))-(1./2.)*dsb[-1,0]*h[-1,0,0]*(h[-1,0,0]+h[-1,0,1])/(2.*h[-1,0,0]+h[-1,0,1])) * dzl.x[0])/Delta;
	}
	else{
	  d2s[0,0,0] += ((s[1,0,0]-h[1,0,0]/6.*(dst[1,0]+2.*dsb[1,0])) * dzl.x[1]  -  (s[-1,0,0]-h[-1,0,0]/6.*(dst[-1,0]+2.*dsb[-1,0])) * dzl.x[0])/Delta;
	}

	d2s[0,0,0] -= dsb[0,0]*(sq((dzl.x[1] + dzl.x[])/2.));
      }
    }


    //Inter_layer terms on internal layer boundaries
    scalar b[]; //could just be a double
    for (int l = 1; l < nl; l++) {
      foreach_face(){
	dzl.x[] += (h[0,0,l-1]-h[-1,0,l-1])/Delta;
      }
      foreach(){
	foreach_dimension(){
	  b[] = (((h[1,0,l-1]*s[1,0,l]+h[1,0,l]*s[1,0,l-1])/(h[1,0,l]+h[1,0,l-1]))*dzl.x[1] - ((h[-1,0,l-1]*s[-1,0,l]+h[-1,0,l]*s[-1,0,l-1])/(h[-1,0,l]+h[-1,0,l-1]))*dzl.x[0])/Delta;
	  b[] -= 2.*(s[0,0,l]-s[0,0,l-1])/(h[0,0,l]+h[0,0,l-1]) * (sq((dzl.x[1] + dzl.x[])/2.));
  
	  d2s[0,0,l]   +=  b[];
	  d2s[0,0,l-1] -=  b[];
	}
      }
    }
    
    
    foreach_face(){
      dzl.x[] += (h[0,0,nl-1]-h[-1,0,nl-1])/Delta;
    }

    //Interlayer terms at last layer top
    foreach(){
      foreach_dimension(){
	if(nl>1){
	  d2s[0,0,nl-1] -= ((((3.*sq(h[1,0,nl-1])+3.*h[1,0,nl-1]*h[1,0,nl-2]+sq(h[1,0,nl-2]))*s[1,0,nl-1]-sq(h[1,0,nl-1])*s[1,0,nl-2])/((h[1,0,nl-1]+h[1,0,nl-2])*(2.*h[1,0,nl-1]+h[1,0,nl-2]))+(1./2.)*dst[1,0]*h[1,0,nl-1]*(h[1,0,nl-1]+h[1,0,nl-2])/(2.*h[1,0,nl-1]+h[1,0,nl-2])) * dzl.x[1]  -  (((3.*sq(h[-1,0,nl-1])+3.*h[-1,0,nl-1]*h[-1,0,nl-2]+sq(h[-1,0,nl-2]))*s[-1,0,nl-1]-sq(h[-1,0,nl-1])*s[-1,0,nl-2])/((h[-1,0,nl-1]+h[-1,0,nl-2])*(2.*h[-1,0,nl-1]+h[-1,0,nl-2]))+(1./2.)*dst[-1,0]*h[-1,0,nl-1]*(h[-1,0,nl-1]+h[-1,0,nl-2])/(2.*h[-1,0,nl-1]+h[-1,0,nl-2])) * dzl.x[0])/Delta;
	}
	else{
	  d2s[0,0,nl-1] -= ((s[1,0,nl-1]+h[1,0,nl-1]/6.*(2.*dst[1,0]+dsb[1,0])) * dzl.x[1]  -  (s[-1,0,nl-1]+h[-1,0,nl-1]/6.*(2.*dst[-1,0]+dsb[-1,0])) * dzl.x[0])/Delta;
	}
	
	d2s[0,0,nl-1] += dst[0,0]*(sq((dzl.x[1] + dzl.x[])/2.));
      }
    }

    //Diffusion term application
    foreach(){
      foreach_layer(){
	if (h[] > dry){
	  s[] += dt*D*d2s[]/(h[]);
	}
      }
    }
    delete ({d2s});
  }
}


void horizontal_diffusion_DirichletNeumann(scalar s, double D, double dt, (const) scalar st, (const) scalar dsb)
{
  if (D > 0.) {
    
    scalar d2s= new scalar [nl];
    face vector dzl[];
    foreach_face(){
      dzl.x[] = (zb[]-zb[-1])/Delta;
    }

    // $\nabla^2(hu)_l$
    foreach(){
      foreach_layer(){
	foreach_dimension(){
	  d2s[] = (h[1]*s[1]-2.*h[]*s[]+h[-1]*s[-1])/sq(Delta);
	}
      }
    }

    // Inter_layer terms at first layer bottom
    foreach(){
      foreach_dimension(){
	if (nl>1){
	  d2s[0,0,0] += ((((3.*sq(h[1,0,0])+3.*h[1,0,0]*h[1,0,1]+sq(h[1,0,1]))*s[1,0,0]-sq(h[1,0,0])*s[1,0,1])/((h[1,0,0]+h[1,0,1])*(2.*h[1,0,0]+h[1,0,1]))-(1./2.)*dsb[1,0]*h[1,0,0]*(h[1,0,0]+h[1,0,1])/(2.*h[1,0,0]+h[1,0,1])) * dzl.x[1]  -  (((3.*sq(h[-1,0,0])+3.*h[-1,0,0]*h[-1,0,1]+sq(h[-1,0,1]))*s[-1,0,0]-sq(h[-1,0,0])*s[-1,0,1])/((h[-1,0,0]+h[-1,0,1])*(2.*h[-1,0,0]+h[-1,0,1]))-(1./2.)*dsb[-1,0]*h[-1,0,0]*(h[-1,0,0]+h[-1,0,1])/(2.*h[-1,0,0]+h[-1,0,1])) * dzl.x[0])/Delta;
	}
	else{
	  d2s[0,0,0] += (((3.*s[1,0,0]-st[1,0])/2.-h[1,0,0]*dsb[1,0]/4.) * dzl.x[1]  -  ((3.*s[-1,0,0]-st[-1,0])/2.-h[-1,0,0] * dsb[-1,0]/4.)*dzl.x[0])/Delta;
	}

	d2s[0,0,0] -= dsb[0,0]*(sq((dzl.x[1] + dzl.x[])/2.));
      }
    }


    //Inter_layer terms on internal layer boundaries
    scalar b[]; //could just be a double
    for (int l = 1; l < nl; l++) {
      foreach_face(){
	dzl.x[] += (h[0,0,l-1]-h[-1,0,l-1])/Delta;
      }
      foreach(){
	foreach_dimension(){
	  b[] = (((h[1,0,l-1]*s[1,0,l]+h[1,0,l]*s[1,0,l-1])/(h[1,0,l]+h[1,0,l-1]))*dzl.x[1] - ((h[-1,0,l-1]*s[-1,0,l]+h[-1,0,l]*s[-1,0,l-1])/(h[-1,0,l]+h[-1,0,l-1]))*dzl.x[0])/Delta;
	  b[] -= 2.*(s[0,0,l]-s[0,0,l-1])/(h[0,0,l]+h[0,0,l-1]) * (sq((dzl.x[1] + dzl.x[])/2.));
  
	  d2s[0,0,l]   +=  b[];
	  d2s[0,0,l-1] -=  b[];
	}
      }
    }
    
    
    foreach_face(){
      dzl.x[] += (h[0,0,nl-1]-h[-1,0,nl-1])/Delta;
    }

    //Interlayer terms at last layer top
    foreach(){
      foreach_dimension(){
	d2s[0,0,0] -= (st[1]*dzl.x[1]-st[-1]*dzl.x[0])/Delta;
	
	if(nl>1){
	  d2s[0,0,nl-1] += (-2.*((3.*sq(h[0,0,nl-1])+3.*h[0,0,nl-1]*h[0,0,nl-2]+sq(h[0,0,nl-2]))*s[0,0,nl-1]-sq(h[0,0,nl-1])*s[0,0,nl-2]-(2.*h[0,0,nl-2]+h[0,0,nl-1])*(h[0,0,nl-1]+h[0,0,nl-2])*st[0,0])/(h[0,0,nl-1]*sq(h[0,0,nl-1]+h[0,0,nl-2]))) * (sq((dzl.x[1] + dzl.x[])/2.));
	}
	else{
	  d2s[0,0,nl-1] += (3.*(st[0,0]-s[0,0,nl-1])/h[0,0,nl-1]-(1./2.)*dsb[0,0]) * (sq((dzl.x[1] + dzl.x[])/2.));
	}
      }
    }

    //Diffusion term application
    foreach(){
      foreach_layer(){
	if (h[] > dry){
	  s[] += dt*D*d2s[]/(h[]);
	}
      }
    }
    delete ({d2s});
  }
}


void horizontal_diffusion_NeumannDirichlet(scalar s, double D, double dt, (const) scalar dst, (const) scalar sb)
{
  if (D > 0.) {
    
    scalar d2s= new scalar [nl];
    face vector dzl[];
    foreach_face(){
      dzl.x[] = (zb[]-zb[-1])/Delta;
    }

    // $\nabla^2(hu)_l$
    foreach(){
      foreach_layer(){
	foreach_dimension(){
	  d2s[] = (h[1]*s[1]-2.*h[]*s[]+h[-1]*s[-1])/sq(Delta);
	}
      }
    }

    // Inter_layer terms at first layer bottom
    foreach(){
      foreach_dimension(){
	d2s[0,0,0] += (sb[1]*dzl.x[1]-sb[-1]*dzl.x[0])/Delta;
	
	if (nl>1){
	  d2s[0,0,0] -= (2.*((3.*sq(h[0,0,0])+3.*h[0,0,0]*h[0,0,1]+sq(h[0,0,1]))*s[0,0,0]-sq(h[0,0,0])*s[0,0,1]-(2.*h[0,0,0]+h[0,0,1])*(h[0,0,0]+h[0,0,1])*sb[0,0])/(h[0,0,0]*sq(h[0,0,0]+h[0,0,1]))) * (sq((dzl.x[1] + dzl.x[])/2.));
	}
	else{
	  d2s[0,0,0] -= (3.*(s[0,0,0]-sb[0,0])/h[0,0,0]-(1./2.)*dst[0,0]) * (sq((dzl.x[1] + dzl.x[])/2.));
	}
      }
    }


    //Inter_layer terms on internal layer boundaries
    scalar b[]; //could just be a double
    for (int l = 1; l < nl; l++) {
      foreach_face(){
	dzl.x[] += (h[0,0,l-1]-h[-1,0,l-1])/Delta;
      }
      foreach(){
	foreach_dimension(){
	  b[] = (((h[1,0,l-1]*s[1,0,l]+h[1,0,l]*s[1,0,l-1])/(h[1,0,l]+h[1,0,l-1]))*dzl.x[1] - ((h[-1,0,l-1]*s[-1,0,l]+h[-1,0,l]*s[-1,0,l-1])/(h[-1,0,l]+h[-1,0,l-1]))*dzl.x[0])/Delta;
	  b[] -= 2.*(s[0,0,l]-s[0,0,l-1])/(h[0,0,l]+h[0,0,l-1]) * (sq((dzl.x[1] + dzl.x[])/2.));
  
	  d2s[0,0,l]   +=  b[];
	  d2s[0,0,l-1] -=  b[];
	}
      }
    }
    
    
    foreach_face(){
      dzl.x[] += (h[0,0,nl-1]-h[-1,0,nl-1])/Delta;
    }

    //Interlayer terms at last layer top
    foreach(){
      foreach_dimension(){
	if(nl>1){
	  d2s[0,0,nl-1] -= ((((3.*sq(h[1,0,nl-1])+3.*h[1,0,nl-1]*h[1,0,nl-2]+sq(h[1,0,nl-2]))*s[1,0,nl-1]-sq(h[1,0,nl-1])*s[1,0,nl-2])/((h[1,0,nl-1]+h[1,0,nl-2])*(2.*h[1,0,nl-1]+h[1,0,nl-2]))+(1./2.)*dst[1,0]*h[1,0,nl-1]*(h[1,0,nl-1]+h[1,0,nl-2])/(2.*h[1,0,nl-1]+h[1,0,nl-2])) * dzl.x[1]  -  (((3.*sq(h[-1,0,nl-1])+3.*h[-1,0,nl-1]*h[-1,0,nl-2]+sq(h[-1,0,nl-2]))*s[-1,0,nl-1]-sq(h[-1,0,nl-1])*s[-1,0,nl-2])/((h[-1,0,nl-1]+h[-1,0,nl-2])*(2.*h[-1,0,nl-1]+h[-1,0,nl-2]))+(1./2.)*dst[-1,0]*h[-1,0,nl-1]*(h[-1,0,nl-1]+h[-1,0,nl-2])/(2.*h[-1,0,nl-1]+h[-1,0,nl-2])) * dzl.x[0])/Delta;
	}
	else{
	  d2s[0,0,nl-1] -= (((3.*s[1,0,nl-1]-sb[1,0])/2.+h[1,0,nl-1]*dst[1,0]/4.) * dzl.x[1]  -  ((3.*s[-1,0,nl-1]-sb[-1,0])/2.+h[-1,0,nl-1]*dst[-1,0]/4.) * dzl.x[0])/Delta;
	}
	
	d2s[0,0,nl-1] += dst[0,0]*(sq((dzl.x[1] + dzl.x[])/2.));
      }
    }

    //Diffusion term application
    foreach(){
      foreach_layer(){
	if (h[] > dry){
	  s[] += dt*D*d2s[]/(h[]);
	}
      }
    }
    delete ({d2s});
  }
}


void horizontal_diffusion_DirichletDirichlet(scalar s, double D, double dt, (const) scalar st, (const) scalar sb)
{
  if (D > 0.) {
    
    scalar d2s= new scalar [nl];
    face vector dzl[];
    foreach_face(){
      dzl.x[] = (zb[]-zb[-1])/Delta;
    }

    // $\nabla^2(hu)_l$
    foreach(){
      foreach_layer(){
	foreach_dimension(){
	  d2s[] = (h[1]*s[1]-2.*h[]*s[]+h[-1]*s[-1])/sq(Delta);
	}
      }
    }

    // Inter_layer terms at first layer bottom
    foreach(){
      foreach_dimension(){
	d2s[0,0,0] += (sb[1]*dzl.x[1]-sb[-1]*dzl.x[0])/Delta;
	
	if (nl>1){
	  d2s[0,0,0] -= 2.*((3.*sq(h[0,0,0])+3.*h[0,0,0]*h[0,0,1]+sq(h[0,0,1]))*s[0,0,0]-sq(h[0,0,0])*s[0,0,1]-(2.*h[0,0,0]+h[0,0,1])*(h[0,0,0]+h[0,0,1])*sb[0,0])/(h[0,0,0]*sq(h[0,0,0]+h[0,0,1])) * (sq((dzl.x[1] + dzl.x[])/2.));
	}
	else{
	  d2s[0,0,0] -= ((6.*s[0,0,0]-2.*st[0,0]-4.*sb[0,0])/h[0,0,0]) * (sq((dzl.x[1] + dzl.x[])/2.));
	}
      }
    }


    //Inter_layer terms on internal layer boundaries
    scalar b[]; //could just be a double
    for (int l = 1; l < nl; l++) {
      foreach_face(){
	dzl.x[] += (h[0,0,l-1]-h[-1,0,l-1])/Delta;
      }
      foreach(){
	foreach_dimension(){
	  b[] = (((h[1,0,l-1]*s[1,0,l]+h[1,0,l]*s[1,0,l-1])/(h[1,0,l]+h[1,0,l-1]))*dzl.x[1] - ((h[-1,0,l-1]*s[-1,0,l]+h[-1,0,l]*s[-1,0,l-1])/(h[-1,0,l]+h[-1,0,l-1]))*dzl.x[0])/Delta;
	  b[] -= 2.*(s[0,0,l]-s[0,0,l-1])/(h[0,0,l]+h[0,0,l-1]) * (sq((dzl.x[1] + dzl.x[])/2.));
  
	  d2s[0,0,l]   +=  b[];
	  d2s[0,0,l-1] -=  b[];
	}
      }
    }
    
    
    foreach_face(){
      dzl.x[] += (h[0,0,nl-1]-h[-1,0,nl-1])/Delta;
    }

    //Interlayer terms at last layer top
    foreach(){
      foreach_dimension(){
	d2s[0,0,0] -= (st[1]*dzl.x[1]-st[-1]*dzl.x[0])/Delta;
	
	if(nl>1){
	  d2s[0,0,nl-1] += (-2.*((3.*sq(h[0,0,nl-1])+3.*h[0,0,nl-1]*h[0,0,nl-2]+sq(h[0,0,nl-2]))*s[0,0,nl-1]-sq(h[0,0,nl-1])*s[0,0,nl-2]-(2.*h[0,0,nl-2]+h[0,0,nl-1])*(h[0,0,nl-1]+h[0,0,nl-2])*st[0,0])/(h[0,0,nl-1]*sq(h[0,0,nl-1]+h[0,0,nl-2]))) * (sq((dzl.x[1] + dzl.x[])/2.));
	}
	else{
	  d2s[0,0,nl-1] += ((6.*s[0,0,nl-1]-4.*st[0,0]-2.*sb[0,0])/h[0,0,nl-1]) * (sq((dzl.x[1] + dzl.x[])/2.));
	}
      }
    }

    //Diffusion term application
    foreach(){
      foreach_layer(){
	if (h[] > dry){
	  s[] += dt*D*d2s[]/(h[]);
	}
      }
    }
    delete ({d2s});
  }
}

/**
# Viscous friction between layers

By default the viscosity is zero and we impose free-slip on the
free-surface and no-slip on the bottom boundary
i.e. $\dot{\mathbf{u}}_t = 0$, $\mathbf{\lambda}_b = 0$, $\mathbf{u}_b
= 0$. */

double nu = 0.;

(const) vector lambda_t[] = {0,0,0}, u_t[] = {0,0,0}, dut[] = {0,0,0}, lambda_b[] = {0,0,0}, u_b[] = {0,0,0}, dub[] = {0,0,0} ;



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
	if (symmetric_bathymetry){
	  if(navier_surface){
	    vertical_diffusion_NavierNeumann (point, h, u.x, dt, nu,
					      u_t.x[], lambda_t.x[], dub.x[]);
	  }
	  else{
	    vertical_diffusion_NeumannNeumann (point, h, u.x, dt, nu,
					       dut.x[], dub.x[]);
	  }
	}
	else{
	   if(navier_surface){
	     vertical_diffusion_NavierNavier (point, h, u.x, dt, nu, 
                                              u_t.x[], lambda_t.x[], u_b.x[], lambda_b.x[]);
	  }
	  else{
            vertical_diffusion_NeumannNavier (point, h, u.x, dt, nu,
					    dut.x[], u_b.x[], lambda_b.x[]);
	  }
	}
    }
    if (h_diffusion){
      foreach_dimension()
	if (symmetric_bathymetry){
	  if(navier_surface){
	    horizontal_diffusion_DirichletNeumann (u.x, nu, dt, u_t.x,dub.x);
	  }
	  else{
	    horizontal_diffusion_NeumannNeumann (u.x, nu, dt, dut.x,dub.x);
	  }
	}
	else{
	  if(navier_surface){
	    horizontal_diffusion_DirichletDirichlet (u.x, nu, dt, u_t.x,u_b.x);
	  }
	  else{
	    horizontal_diffusion_NeumannDirichlet (u.x, nu, dt, dut.x,u_b.x);
	  }
	}
    }
    foreach() {
      foreach_layer()
	foreach_dimension()
	  u.x[] -= dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
    }
  }
}

/**
## References

~~~bib
@hal{popinet2020, hal-02365730}

@hal{devita2019, hal-02295398}
~~~
*/
