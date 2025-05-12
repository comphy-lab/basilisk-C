/**
[Return to my homepage](http://basilisk.fr/sandbox/nlemoine/README)

# Friction term with detrended topography

## Rationale for topography detrending
The scalar continuity equation and the vector equation for momentum can be written as a single vector equation describing the evolution of the state vector of conservative variables:

$$
\displaystyle\frac{\partial}{\partial t} \left[\begin{array}{c} h \\ h\,u_x \\ h\,u_y \end{array}\right] = -\boldsymbol{\nabla}\cdot
\left[\begin{array}{cc} h\,u_x & h\,u_y \\ h\,u_x^2+ \frac{1}{2} g h^2 & h\,u_x\,u_y \\ h\,u_x\,u_y & h\,u_y^2 + \frac{1}{2} g h^2\end{array}\right]
\ -\ g\,h\left[\begin{array}{c} 0 \\ \partial_x z_b \\ \partial_y z_b \end{array}\right]
\ -\ \frac{1}{\rho}\left[\begin{array}{c} 0 \\ \tau_x \\ \tau_y \end{array}\right]$$

$$
\partial_t\mathbf{W}\qquad\quad = \quad\qquad-\boldsymbol{\nabla}\cdot\mathbf{F}(\mathbf{W})
\qquad\qquad +\qquad\qquad \mathbf{S_0}(\mathbf{W})
\qquad +\qquad \mathbf{S_f}(\mathbf{W})
$$
By default, the Saint-Venant solver in Basilisk takes the topographic source term $\mathbf{S_0}$ into account but the friction source term must be handled separately. We can thus change the decomposition between topographic and friction terms.

Assume that we can write bed elevation $z_b$ as the superposition of a signal $z^\ast_b$ with zero- or constant spatial mean, and a linear trend (plane) :

$$z_b(x,y) = z^\ast_b(x,y) - I_x\,x - I_y\,y$$

where $z^\ast_b(x,y)$ is the detrended topography (i.e. "flat in average"), and $I_x$ and $I_y$ the "*tilts*" in directions $x$ and $y$ respectively (put differently, if $z_b$ is a plane then the definitions are $I_x~\!=~\!-\partial_x z_b$ and $I_y=-\partial_y z_b$). In this case, the set of equations becomes:
$$\partial_t\mathbf{W}\quad = \quad-\boldsymbol{\nabla}\cdot\mathbf{F}(\mathbf{W})
\ -\ g\,h\left[\begin{array}{c} 0 \\ \partial_x z^\ast_b \\ \partial_y z^\ast_b \end{array}\right]
\qquad+\underbrace{g\,h\left[\begin{array}{c} 0 \\ I_x \\ I_y \end{array}\right]
\ -\ \frac{1}{\rho}\left[\begin{array}{c} 0 \\ \tau_x \\ \tau_y \end{array}\right]}_{}$$
$$
\partial_t\mathbf{W}\quad = \quad-\boldsymbol{\nabla}\cdot\mathbf{F}(\mathbf{W})
\quad+\qquad\mathbf{S^\ast_0}(\mathbf{W})
\quad+\qquad\qquad \mathbf{S'_f}(\mathbf{W})\qquad\qquad
$$

We can then give the detrended $z^\ast_b$ as input in Basilisk, and add the effect of regional tilt in the second source term $\mathbf{S'_f}$ that must be provided anyway to take friction into account. A steady-state solution, where the topographic source term exactly balances the friction term, will translate into a "lake-at-rest" water surface profile. This decomposition is also necessary for [periodic solutions](http://basilisk.fr/src/test/lake-tr.c).

##Implementation with Manning friction

In the Manning-Strickler model, bed shear stress is given by:
$$\boldsymbol{\tau}\quad=\quad\rho\,g\,n^2 h^{-\frac{1}{3}}\left|\mathbf{u} \right|\mathbf{u}\quad=\quad\rho\,g\,n^2 h^{-\frac{7}{3}}\left|\mathbf{q} \right|\mathbf{q}$$
The numerical scheme first solves the equation without the friction term, and then adds a correction. The corrector step amounts to solving the following equation:

$$
\displaystyle\frac{\partial}{\partial t} \left[\begin{array}{c}h \\ h\,u_x \\ h\,u_y \end{array}\right] =
\left[\begin{array}{c} 0 \\ +g\,h\,I_x - g\,n^2\,h^{-\frac{1}{3}}\left|\mathbf{u} \right|u_x \\
 +g\,h\,I_y - g\,n^2\,h^{-\frac{1}{3}}\left|\mathbf{u}\right|u_y
\end{array}\right]
$$

In this step, depth is assumed constant ($\partial_t h = 0$). We are left with:

$$
\displaystyle\frac{\partial}{\partial t} \left[\begin{array}{c}u_x \\ u_y \end{array}\right] =
\left[\begin{array}{c} +g\,I_x - g\,n^2\,h^{-\frac{4}{3}}\left|\mathbf{u} \right|u_x \\
 +g\,I_y - g\,n^2\,h^{-\frac{4}{3}}\left|\mathbf{u}\right|u_y
\end{array}\right]
$$

This equation can be rather stiff for high velocities so we cannot treat it in a fully explicit manner. Consider the evolution of component $u_x$. In order to compute its change from time $t_k$ to time $t_{k+1}=t_k+\Delta t$, we linearize the quadratic velocity term and we write it implicitely:

$$\left|\mathbf{u}\right|\,u_x \simeq |\mathbf{u}^{(k)}|\,u_x^{(k+1)}$$

It yields:

$$
u_x^{(k+1)}-u_x^{(k)} \simeq \Delta t \left[g\,I_x - g\,n^2\,|\mathbf{u}^{(k)}|h^{-\frac{4}{3}}u_x^{(k+1)}\right]$$
$$u_x^{(k+1)}\Big[1 + \underbrace{g\,n^2\,|\mathbf{u}^{(k)}|h^{-\frac{4}{3}}\Delta t}_{s}\Big] \simeq u_x^{(k)} + g\,I_x\,\Delta t
$$
We define the dimensionless quantity:
$$
s = g\,n^2\,|\mathbf{u}^{(k)}|h^{-\frac{4}{3}}\Delta t $$
And finally:
$$u_x^{(k+1)}=\frac{u_x^{(k)} + g\,I_x\,\Delta t}{1+s}$$
*/

// Manning coefficient
scalar nmanning[];
// Tilt
coord tilt = {.x = 0., .y = 0.};

event manningsourceterm(i++){
  
  foreach(){
    if(h[] > dry){      
      double s = dt*G*sq(nmanning[])*norm(u)/pow(h[],4./3.);
      foreach_dimension()
        u.x[] = (u.x[]+G*dt*tilt.x)/(1.+s);
    }
  }
  boundary ((scalar *){u});
}