/**
 
 This is `diffusionNN.h`   Non-Newtonian diffusion.
 This is exactly [http://basilisk.fr/src/layered/diffusion.h](), but we change the
 definition of viscosity which was a given constant by a table of values.
 $\nu$ changed in $(\nu_{l}+\nu_{l+1})/2$ or
 in $(\nu_{l}+\nu_{l-1})/2$ thereafter.
 
 
 These values are  function
 of shear
 $$\nu \rightarrow\nu_{eq}(\frac{\partial u}{\partial z},\nu)$$ the later `nu_eq(shear,nu)` is defined outside this file.
 
 
# Viscous friction between layers
 
 Boundary conditions on the top and bottom layers need to be added to close the
 system for the viscous stresses. We chose to impose a Neumann condition on the
 free-surface i.e.
 $$
 \partial_z u |_t = \dot{u}_t
 $$
 and a Navier slip condition on the bottom i.e.
 $$
 u|_b = u_b + \lambda_b \partial_z u|_b
 $$
 By default the viscosity is zero and we impose free-slip on the
 free-surface and no-slip on the bottom boundary i.e. $\dot{u}_t = 0$,
 $\lambda_b = 0$, $u_b = 0$. */

double nu = 0.;
(const) scalar lambda_b = zeroc, dut = zeroc, u_b = zeroc;

/**
 For stability, we discretise the viscous friction term implicitly as
 $$
 \frac{(hu_l)^{n + 1} - (hu_l)^{\star}}{\Delta t} =
 \nu \left( \frac{u_{l + 1} - u_l}{h_{l + 1 / 2}} -
 \frac{u_l - u_{l - 1}}{h_{l - 1 / 2}} \right)^{n + 1}
 $$
 which can be expressed as the linear system
 $$
 \mathbf{Mu}^{n + 1} = \mathrm{rhs}
 $$
 where $\mathbf{M}$ is a
 [tridiagonal matrix](https://en.wikipedia.org/wiki/Tridiagonal_matrix).
 The lower, principal and upper diagonals are *a*, *b* and *c* respectively. */

double nu_eq(double,double);


void vertical_viscosity (Point point, scalar h, scalar s, double dt)
{
    double nueq[nl];
    double a[nl], b[nl], c[nl], rhs[nl];
    double u[nl] ;
    
    /**
     The *rhs* of the tridiagonal system is $h_lu_l$. */
    
    foreach_layer()
    rhs[_layer] = s[]*h[];
    
    /**
     we get the value of velocity    */
    foreach_layer()
    u[_layer]= s[];
    
    /**
     definition of pressure,  `eta[] = zb[];  foreach_layer()  eta[] += h[]`;
     */
      double  press = max(G*(eta[]-zb[]),G*dry);
    /**
     definition of viscosity as a function of shear
     $$\nu \rightarrow \nu_{eq}(\frac{\partial u}{\partial z},\nu)$$
     */
    for (int l = 1; l < nl - 1; l++) {
        double h12 = (h[0,0,l-1] + h[0,0,l])/2;
      //  press -= G*h[0,0,l]/2;
        double shear = (h12>0? (u[l]-u[l-1])/h12: HUGE);
        nueq[l]= nu_eq(shear,press);
        press -= G*h[0,0,l] ;
    }
    nueq[0] = nueq[1];
    
    
    /**
     from here, $\nu$ changed in $(\nu_{l}+\nu_{l+1})/2$ or
     in $(\nu_{l}+\nu_{l-1})/2$ in `a`,`c` etc
     */
    
    
    /**
     The lower, principal and upper diagonals $a$, $b$ and $c$ are given by
     $$
     a_{l > 0} = - \left( \frac{\nu \Delta t}{h_{l - 1 / 2}} \right)^{n + 1}
     $$
     $$
     c_{l < \mathrm{nl} - 1} = - \left( \frac{\nu \Delta t}{h_{l + 1 / 2}}
     \right)^{n + 1}
     $$
     $$
     b_{0 < l < \mathrm{nl} - 1} = h_l^{n + 1} - a_l - c_l
     $$
     */
    
    for (int l = 1; l < nl - 1; l++) {
        a[l] = - (nueq[l]+nueq[l-1])*dt/(h[0,0,l-1] + h[0,0,l]);
        c[l] = - (nueq[l+1]+nueq[l])*dt/(h[0,0,l] + h[0,0,l+1]);
        b[l] = h[0,0,l] - a[l] - c[l];
    }
    
    /**
     For the top layer the boundary conditions give the (ghost)
     boundary value
     $$
     u_{\mathrm{nl}} = u_{\mathrm{nl} - 1} + \dot{u}_t h_{\mathrm{nl} - 1},
     $$
     which gives the diagonal coefficient and right-hand-side
     $$
     b_{\mathrm{nl} - 1} = h_{\mathrm{nl} - 1}^{n + 1}
     - a_{\mathrm{nl} - 1}
     $$
     $$
     \mathrm{rhs}_{\mathrm{nl} - 1} =
     (hu)_{\mathrm{nl} - 1}^{\star} + \nu \Delta t \dot{u}_t
     $$
     */
    
    a[nl-1] = - (nueq[nl-1]+nueq[nl-2])*dt/(h[0,0,nl-2] + h[0,0,nl-1]);
    b[nl-1] = h[0,0,nl-1] - a[nl-1];
    rhs[nl-1] += nueq[nl-1]*dt*dut[];
    
    /**
     For the bottom layer a third-order discretisation of the Navier slip
     condition gives
     $$
     \begin{aligned}
     b_0 & = h_0 + 2 \Delta t \nu \left( \frac{1}{h_0 + h_1} + \frac{h^2_1 + 3
     h_0 h_1 + 3 h^2_0}{\det} \right),\\
     c_0 & = - 2 \Delta t \nu \left( \frac{1}{h_0 + h_1} + \frac{h^2_0}{\det}
     \right),\\
     \text{rhs}_0 & = (hu_0)^{\star} + 2 \Delta t \nu u_b  \frac{h^2_1 + 3 h_0
     h_1 + 2 h^2_0}{\det},\\
     \det & = h_0 (h_0 + h_1)^2  + 2\lambda (3\,h_0 h_1 + 2\,h_0^2 + h_1^2),
     \end{aligned}
     $$
     */
    
    double den = h[]*sq(h[] + h[0,0,1])
    + 2.*lambda_b[]*(3.*h[]*h[0,0,1] + 2.*sq(h[]) + sq(h[0,0,1]));
    b[0] = h[] + dt*(nueq[0]+nueq[1])*(1./(h[] + h[0,0,1]) +
                           (sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 3.*sq(h[]))/den);
    c[0] = - dt*(nueq[0]+nueq[1])*(1./(h[] + h[0,0,1]) + sq(h[])/den);
    rhs[0] += dt*(nueq[0]+nueq[1])*u_b[]*(sq(h[0,0,1]) + 3.*h[]*h[0,0,1] + 2.*sq(h[0]))/den;
    
    if (nl == 1) {
        b[0] += c[0];
        rhs[0] += (- c[0]*h[] - nu*dt) * dut[];
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
            vertical_viscosity (point, h, u.x, dt);
            foreach_layer()
            foreach_dimension()
            u.x[] -= dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
        }
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
