/**
The growth rate of the perturbation can finally be written as follow : $$\sigma^2 = \dfrac{\Delta \rho g k}{2 \bar{\rho}} - \dfrac{\gamma k^3}{2 \bar\rho}$$

$$x \ in [0,L]$$

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