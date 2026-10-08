/**
The growth rate of the perturbation can finally be written as follow : $$\sigma^2 = \dfrac{\Delta \rho g k}{2 \bar{\rho}} - \dfrac{\gamma k^3}{2 \bar\rho}$$

Note that $\sigma^2$ is real. Meaning that the growth rate is either purely real, or purely imaginary. In other words, the system is either unstable or neutrally stable. In this last case, the oscillations won't be linearly dumped.

We will use the following scale for the quantites (the scale $D$ is a typical lengthscale used for the coordinates and the interface) :
$$P = \gamma k, D = 1/k, U = \left(\dfrac{\gamma k}{\bar\rho }\right)^{1/2}, \tau \sim \left(\dfrac{\bar\rho}{\gamma k^3}\right)^{1/2}$$

We have 2 dimensionless numbers :

$$r = \dfrac{\Delta \rho}{2 \bar\rho}, \text{Bo} = \dfrac{\bar\rho g}{\gamma k^2}$$

And we might need the Ohnesorge numbers : $Oh_i = \nu_i \left(\dfrac{\bar\rho}{\gamma L}\right)^{1/2}$

The equations becomes (now, all the variables are non-dimensionnal):

$$\dfrac{Du}{Dt} = - \dfrac{1}{1 \pm r} \nabla p - \text{Bo} \mathbf{e_z} + \text{Oh}_i \Delta \mathbf{u}$$

$$\Delta P = \partial^2_{xx} \eta$$

The dispersion relation becomes : $$\sigma^2 = r \text{Bo} - \dfrac{1}{2}$$

Now, we can simulate in a box of arbitrary size $L = 1$, focusing on the wavenumber ($k = 2 \pi$).

In the code, we have to set $\rho_i = 1 \pm r$, and the gravity $G = \text{Bo}$
*/