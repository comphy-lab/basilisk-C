/**
The growth rate of the perturbation can finally be written as follow : $$\sigma^2 = \dfrac{\Delta \rho g k}{2 \bar{\rho}} - \dfrac{\gamma k^3}{2 \bar\rho}$$

Note that $\sigma^2$ is real. Meaning that the growth rate is either purely real, or purely imaginary. In other words, the system is either unstable or neutrally stable. In this last case, the oscillations won't be linearly dumped.

We will use the following scale for the quantites :
$$P \sim \dfrac{\gamma}{L}, x,\eta \sim L \sim 1/k, U = \left(\dfrac{\gamma}{\bar\rho L}\right)^{1/2}, \tau \sim \left(\dfrac{\bar\rho L^3}{\gamma}\right)^{1/2}$$

We have 2 dimensionless numbers :

$$r = \dfrac{\Delta \rho}{2 \bar\rho}, Bo = \dfrac{\bar\rho g L^2}{\gamma}$$

And we might need the Ohnesorge numbers : $Oh_i = \nu_i \left(\dfrac{\bar\rho}{\gamma L}\right)^{1/2}$

The equations becomes (now, all the variables are non-dimensionnal):

$$\dfrac{Du}{Dt} = - \dfrac{1}{1 \pm r} \nabla p - Bo \mathbf{e_z} + Oh_i \Delta \mathbf{u}$$

$$\Delta P = \partial^2_{xx} \eta$$

The dispersion relation becomes : $$\sigma^2 = \dfrac{1}{Bo} - \dfrac{1}{2}$$
*/