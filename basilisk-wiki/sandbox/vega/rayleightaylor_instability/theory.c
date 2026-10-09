/**
The growth rate of the perturbation can finally be written as follow : $$\sigma^2 = \dfrac{\Delta \rho g k}{2 \bar{\rho}} - \dfrac{\gamma k^3}{2 \bar\rho}$$

Note that $\sigma^2$ is real. Meaning that the growth rate is either purely real, or purely imaginary. In other words, the system is either unstable or neutrally stable. In this last case, the oscillations won't be linearly dumped.

We will use the following scale for the quantites (the scale $D$ is a typical lengthscale used for the coordinates and the interface) :
$$P = \dfrac{k}{\bar \rho g}, D = 1/k, U = \sqrt{\dfrac{g}{k}}, \tau = \dfrac{1}{\sqrt{gk}}$$

We have 2 dimensionless numbers :

$$r = \dfrac{\Delta \rho}{2 \bar\rho}, \text{Bo} = \dfrac{\bar\rho g}{\gamma k^2}$$

And we might need the Galilei numbers : $\text{Ga}_i = \dfrac{g}{k^3 \nu_i^2}$

The equations becomes (now, all the variables are non-dimensionnal):

$$\dfrac{Du}{Dt} = - \dfrac{1}{1 \pm r} \nabla p - \mathbf{e_z} + \dfrac{1}{\text{Ga}_i^{1/2}} \Delta \mathbf{u}$$

$$\Delta P = \dfrac{1}{Bo} \partial^2_{xx} \eta$$

The dispersion relation becomes : $$\sigma^2 = r - \dfrac{1}{2\text{Bo}}$$

Now, we can simulate in a box of arbitrary size $L = 1$, focusing on the wavenumber ($k = 2 \pi$).

In the code, we have to set $\rho_i = 1 \pm r$, and the surface tension $\gamma = 1/\text{Bo}$
*/