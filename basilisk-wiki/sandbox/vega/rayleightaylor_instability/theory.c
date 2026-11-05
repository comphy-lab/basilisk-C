/**

Initially, the fluid is at rest, so $\mathbf{U} = \mathbf{0}$. Using the incompressibilty condition on both side of the interface, the linearized Navier-Stokes equation can be written:

$$\partial_t \mathbf{u_i} = - \dfrac{1}{\rho_i} \nabla p_i + \nu_i \Delta \mathbf{u_i}$$

$$x \ in [0,L]$$

In the linearized settings, the curvature $\kappa = $ can be approximated by $\partial^2_{xx}\eta$, so the (dynamic) continuity equation at the interface is:

$$[P+p]_i = - \gamma \eta''$$

Where $P$ is the non-perturbed pressure, which is the hydrostatic one: $P = -\rho_i z$.

The last equation is the (cinematic) continuity equation at the interface:

$$\partial_t \eta = u_{zi}$$

Now, using the linearity and the x-invariance, we can focus on a single Fourier mode in $x$ : which will be called $k$. and using the linearity: $\partial_t = \sigma$.

Analysis

If we take the curl of the NS equation, it will leads to the Orr-Sommerfeld equation. But if we take its divergence, we will be the Poisson equation for the Pressure.

$$\Delta p = 0 \Longrightarrow (\partial^2_{zz} - k^2) p \Longrightarrow \partial^2_{zz} p = k^2 p \Longrightarrow p = p_{\pm} e^{\pm kz} $$

As the perturbation must vanishes as $\mid z \mid \longrightarrow + \infty $,  we have that $i=1$ (resp. $i=2$) is associated with the $+$ upperscript (resp. $-$)

Now, using this results in the linearized NS equation, we will get that $\mathbf{u}$ behave like $e^{ikx \pm kz}$, it implies that $\mathbf{u}$ is harmonic, hence the diffusive term vanishes.

Now, using the vertical component of the linearized equation we get $\sigma u_{zi} = \mp k/\rho_i p^{\pm}$, using the cinematic continuity equation, it becomes:

$$p_i^{\pm} = \mp \dfrac{\sigma^2 \rho_i}{k} \eta$$

The total difference of pressure at the interface is:

$$\gamma k^2 \eta = [P + p]_i = \Delta \rho g \eta - \dfrac{\sigma^2}{k} (\rho_1 + \rho_2) \eta$$


The growth rate of the perturbation can finally be written as follow : $$\sigma^2 = \dfrac{\Delta \rho g k}{2 \bar{\rho}} - \dfrac{\gamma k^3}{2 \bar\rho}$$


Note that $\sigma^2$ is real. Meaning that the growth rate is either purely real, or purely imaginary. In other words, the system is either unstable or neutrally stable. In this last case, the oscillations won't be linearly dumped.

## Dimensionless equations

We have two different choices for the lengthscale. We could either use $1/k$ or $L$. If the lengthscale is $1/k$, then, the interval will be $x \in [0, 2 \pi]$. On the other side, using $L$ will transform the interval to be $x \in [0,1]$, but the equations (dispersion relation) will make appear the dimensionless number $kL = 2\pi$.

We will use the following scale for the quantites (the scale $D$ is a typical lengthscale used for the coordinates and the interface) :
$$P = \dfrac{k}{\bar \rho g}, D = 1/k, U = \sqrt{\dfrac{g}{k}}, \tau = \dfrac{1}{\sqrt{gk}}$$

We have 2 dimensionless numbers :

$$r = \dfrac{\Delta \rho}{2 \bar\rho}, \text{Bo} = \dfrac{\bar\rho g}{\gamma k^2}$$

And we might need the Galilei numbers : $\text{Ga}_i = \dfrac{g}{k^3 \nu_i^2}$

The equations becomes (now, all the variables are non-dimensionnal):

$$\dfrac{Du}{Dt} = - \dfrac{1}{1 \pm r} \nabla p - \mathbf{e_z} + \dfrac{1}{\text{Ga}_i^{1/2}} \Delta \mathbf{u}$$

$$\Delta P = \dfrac{1}{Bo} \partial^2_{xx} \eta$$

$$x \in [0,1]$$

The dispersion relation becomes : $$\sigma^2 = 2\pi r - \dfrac{4 \pi^3}{\text{Bo}}$$

Now, we can simulate in a box of arbitrary size $L = 1$, focusing on the wavenumber ($k = 2 \pi$).

In the code, we have to set $\rho_i = 1 \pm r$, and the surface tension $\gamma = 1/\text{Bo}$
*/