/**
# Prescribed outgoing flow direction

This event updates ghost cell values on the right boundary so that:
$$\mathbf{u}\cdot\mathbf{t} = 0$$
$$\partial_{\mathbf{n}}\left(\mathbf{u}\cdot\mathbf{n}\right) = \mathbf{n}\cdot\boldsymbol{\nabla}\left(\mathbf{u}\cdot\mathbf{n}\right)= 0$$
$$\partial_{\mathbf{n}}h = \mathbf{n}\cdot\boldsymbol{\nabla}h= 0$$
$$\partial_{\mathbf{n}}z_b = \mathbf{n}\cdot\boldsymbol{\nabla}z_b = 0$$
where $\mathbf{n}$ is an outflow direction and $\mathbf{t}=[-n_y\,,n_x]^T$ the orthogonal direction, with $\mathbf{n}$ *not aligned* with the grid.
It amounts to setting a free outflow boundary condition but with a direction which is not normal to the square domain.

The condition for $h$ also reads:
$$n_x\partial_x h + n_y\partial_y h =0 $$
$$n_x\frac{h[\textrm{ghost}]-h[\ ]}{\Delta} +
n_y\frac{h[0,1]-h[0,-1]}{2\Delta} =0$$

Note that the derivative in $y$ is evaluated at domain cell center and not on the boundary in order to avoid involving more than one ghost cell (otherwise, we have to solve a linear system to update all ghost cell values consistently). Hence the pseudo-Neumann condition for $h$:

$$h[\textrm{ghost}]=h[\ ]-\left(\frac{n_y}{n_x}\right)
\frac{h[0,1]-h[0,-1]}{2}$$

The pseudo-Neumann condition for $z_b$ reads similarly (warning: the condition is *homogeneous*, so for a river reach [detrending](manning-tilt.h) is needed):

$$z_b[\textrm{ghost}]=z_b[\ ]-\left(\frac{n_y}{n_x}\right)
\frac{z_b[0,1]-z_b[0,-1]}{2}$$

The condition for $\mathbf{u}\cdot\mathbf{t}$ reads:

$$-u_x|_{\textrm{boundary}}\, n_y + u_y|_{\textrm{boundary}}\, n_x = 0$$
$$-\left(\frac{u.x[\ ]+u.x[\textrm{ghost}]}{2}\right)\frac{n_y}{n_x} + 
\left(\frac{u.y[\ ]+u.y[\textrm{ghost}]}{2}\right) = 0\qquad (1)
$$

The condition for $\mathbf{u}\cdot\mathbf{n}$ reads:

$$n_x^2\partial_x u_x + n_x n_y \partial_x u_y + n_x n_y \partial_y u_x + n_y^2 \partial_y u_y = 0$$
$$\partial_x u_x + \frac{n_y}{n_x}\partial_x u_x + \frac{n_y}{n_x}
\underbrace{\left[\partial_y u_x + \frac{n_y}{n_x}\partial_y u_y\right]}_{A/\Delta} = 0$$

$$\left(\frac{u.x[\textrm{ghost}]-u.x[\ ]}{\Delta}\right)+
\frac{n_y}{n_x}\left(\frac{u.y[\textrm{ghost}]-u.y[\ ]}{\Delta}\right)+
\frac{n_y}{n_x}\frac{A}{\Delta} = 0\qquad (2)
$$


From (1) we get:

$$u.y[\textrm{ghost}] = -u.y[\ ] + \textstyle\frac{n_y}{n_x}
\Big(u.x[\ ]+u.x[\textrm{ghost}]\Big)$$

Reinjecting in (2) yields:

$$\left(1+\textstyle\left(\frac{n_y}{n_x}\right)^2\right)u.x[\textrm{ghost}] = \left(1-\textstyle\left(\frac{n_y}{n_x}\right)^2\right)u.x[\ ]
+2\textstyle\frac{n_y}{n_x}u.y[\ ] - \textstyle\frac{n_y}{n_x}A
$$

TODO : adapt code for left, top, and bottom.
TODO 2 : use scalar field localangle[] so that outflow direction can vary along the border (useful when several streams leave the domain)
*/

event downstreamBC (i++)
{
  double tana = n_out.y / n_out.x ;

  foreach_boundary(right)	
  {
	zb[ghost] = zb[] - tana*(zb[0,1]-zb[0,-1]);
 	 h[ghost] =  h[] - tana*( h[0,1]- h[0,-1]);
    /**
Compute $$A = \Delta * \left[\partial_y u_x + \textstyle\frac{n_y}{n_x}\partial_y u_y\right]$$
    */	
    // Derivatives with respect to y are evaluated at domain cell center, not face center (avoids involving more than just 1 ghost cell)
        double A = 0.5*(u.x[0,1]-u.x[0,-1])+0.5*tana*(u.y[0,1]-u.y[0,-1]); 
        u.x[ghost] = (1.0-tana*tana)*u.x[] + 2.0*tana*u.y[] - tana*A;
        u.x[ghost] = u.x[ghost]/(1.0+tana*tana); 
        u.y[ghost] = -u.y[] + tana*(u.x[]+u.x[ghost]);
  }
}