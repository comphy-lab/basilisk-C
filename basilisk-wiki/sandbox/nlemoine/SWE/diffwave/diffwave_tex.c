/**
The diffusive wave approximation stems from the Saint-Venant Shallow Water Equations, which can be written in full form as:

$$
\displaystyle\frac{\partial}{\partial t} \left[\begin{array}{c} h \\ h\,u_x \\ h\,u_y \end{array}\right] = -\boldsymbol{\nabla}\cdot
\left[\begin{array}{cc} h\,u_x & h\,u_y \\ h\,u_x^2+ \frac{1}{2} g h^2 & h\,u_x\,u_y \\ h\,u_x\,u_y & h\,u_y^2 + \frac{1}{2} g h^2\end{array}\right]
\ -\ g\,h\left[\begin{array}{c} 0 \\ \partial_x z_b \\ \partial_y z_b \end{array}\right]
\ -\ \frac{1}{\rho}\left[\begin{array}{c} 0 \\ \tau_x \\ \tau_y \end{array}\right]$$

The first line is the continuity equation, while the second and third are the projections of the dynamic (momentum) equation. The diffusive wave approximation is a zero-inertia approximation, obtained by neglecting the inertial acceleration terms in the left-hand side and the convective acceleration terms in the right-hand side:

$$\displaystyle\left[\begin{array}{c} \frac{\partial h}{\partial t} \\
\cancel{\frac{\partial(h\,u_x)}{\partial t}} \\
\cancel{\frac{\partial(h\,u_y)}{\partial t}} 
\end{array}\right] = -\boldsymbol{\nabla}\cdot
\left[\begin{array}{cc} h\,u_x & h\,u_y \\ \frac{1}{2} g h^2 & 0 \\ 0 & \frac{1}{2} g h^2\end{array}\right]
-\quad\cancel{
\boldsymbol{\nabla}\cdot
\left[\begin{array}{cc} 0 & 0 \\ h\,u_x^2 & h\,u_x\,u_y \\ h\,u_x\,u_y & h\,u_y^2 \end{array}\right]
}
\ -\ g\,h\left[\begin{array}{c} 0 \\ \partial_x z_b \\ \partial_y z_b \end{array}\right]
\ -\ \frac{1}{\rho}\left[\begin{array}{c} 0 \\ \tau_x \\ \tau_y \end{array}\right]$$

Denoting $q_x = h\,u_x$ and $q_y = h\,u_y$ the components of the flux density vector, the system boils down to:

$$\begin{cases}
   \frac{\partial h}{\partial t} = -\boldsymbol{\nabla}\cdot\mathbf{q} & \\
   0 = -g h \frac{\partial h}{\partial x} -g h \frac{\partial z_b}{\partial x} - \frac{1}{\rho}\tau_x & \\
   0 = -g h \frac{\partial h}{\partial y} -g h \frac{\partial z_b}{\partial y} - \frac{1}{\rho}\tau_y & \\
   \end{cases}
$$

The dynamic equation (last two scalar equations) is limited to a quasi-equilibrium between friction and the resultant force of the pressure gradient. Indeed, denoting $\eta = z_b + h$ the water surface elevation, we have:

$$\begin{cases}
   \frac{\partial h}{\partial t} = -\boldsymbol{\nabla}\cdot\mathbf{q} & \\
   0 = -g h \frac{\partial \eta}{\partial x} -\frac{1}{\rho}\tau_x & \\
   0 = -g h \frac{\partial \eta}{\partial y} - \frac{1}{\rho}\tau_y & \\
   \end{cases}
$$

In the following we chose the Manning friction model, but any other model can be chosen. We then have $\quad\boldsymbol{\tau} = \rho\,g\,n^2\,h^{-\frac{7}{3}}\|\mathbf{q}\|\mathbf{q}\quad$ and then:

$$\begin{cases}
   q_x = -\displaystyle\frac{1}{n^2\|\mathbf{q}\|}h^{\frac{10}{3}}\frac{\partial\eta}{\partial x} & \\
   & \\
   q_y = -\displaystyle\frac{1}{n^2\|\mathbf{q}\|}h^{\frac{10}{3}}\frac{\partial\eta}{\partial y} & \\
   \end{cases}
$$
Hence
$$\|\mathbf{q}\| = \displaystyle\frac{1}{n^2\|\mathbf{q}\|}h^{\frac{10}{3}}\|\boldsymbol{\nabla}\eta\|$$

We see that we now have an explicit expression for the norm of the flux density vector, and hence for the vector itself:

$$\|\mathbf{q}\| = \displaystyle\frac{1}{n}h^{\frac{5}{3}}\|\boldsymbol{\nabla}\eta\|^{\frac{1}{2}}$$

$$\mathbf{q} = \displaystyle-\frac{1}{n}\big(\eta-z_b\big)^{\frac{5}{3}}\|\boldsymbol{\nabla}\eta\|^{-\frac{1}{2}}\boldsymbol{\nabla}\eta$$

Finally we are left with a single equation (the continuity equation) which can be written in the form of a diffusion equation:

$$\frac{\partial\eta}{\partial t} = +\boldsymbol{\nabla}\cdot\left(
\underbrace{\frac{1}{n}\big(\eta-z_b\big)^{\frac{5}{3}}\|\boldsymbol{\nabla}\eta\|^{-\frac{1}{2}}}_{
\normalsize D\left(\eta,\|\boldsymbol{\nabla}\eta\|\right)
}\boldsymbol{\nabla}\eta
\right)$$


The problem is nonlinear since the diffusivity $D$ depends both on $\eta$ and on the norm of the gradient of $\eta$. However it can be easily implemented in Basilisk thanks to the [Poisson solver](basilisk.fr/src/diffusion.h), using an implicit scheme with a staggered grid diffusivity.

\begin{figure}[h!]
 \centering\includegraphics[width=0.4\linewidth]{staggered.png}
 \caption{}
 \label{fig:staggered}
\end{figure}

It is necessary to desingularize the norm of the gradient of $\eta$, which has exponent $-\frac{1}{2}$ in the diffusivity. We use:

$$ D\left(\eta,\|\boldsymbol{\nabla}\eta\|\right) = \frac{1}{n}\big(\eta-z_b\big)^{\frac{5}{3}}\Big(\|\boldsymbol{\nabla}\eta\|+\epsilon\Big)^{-\frac{1}{2}} $$

*/