/**
# Recoil Pressure

We include the recoil pressure contribution in the pressure jump as an
[interfacial force](/src/iforce.h):
$$
  \left[p\right]_\Gamma = -\dot{m}^2\left[\dfrac{1}{\rho}\right]_\Gamma
$$
According to [Mialhe et al. 2023](#mialhe2023extended), this term enables to
enforce momentum conservation at the interface, by counteracting the momentum
loss due to the phase change. */

/**
We need the interfacial force module as well as some
functions to compute the position of the interface. */

#include "iforce.h"
#include "curvature.h"

/**
We need external fields that are defined in [two-phase.h](/src/two-phase.h) and
in [evaporation.h](evaporation.h). */

extern scalar mEvapTot;

/**
The function *recoil_potential()* fills the interfacial potential $\phi$
according to the vaporization rate and the density jump at the interface. */

void recoil_potential (scalar f, scalar phi, bool add = false)
{
  foreach() {
    double hp = 0.;
    if (interfacial (point, f)) {
      hp = (rho2v[] > 0. && rho1v[] > 0.) ?
        -sq(mEvapTot[])*(1./rho2v[] - 1./rho1v[]) : 0.;
      if (add)
        phi[] += hp;
      else
        phi[] = hp;
    }
    else
      phi[] = nodata;
  }
}

/**
We overload the acceleration() event to add the contribution of
recoil pressure to the interfacial potential $\phi$.

If $\phi$ is already allocated, we add the contribution of recoil pressure,
otherwise we allocate a new field and set it to the same contribution. */

event acceleration (i++)
{
  scalar phi = f.phi;
  if (phi.i)
    recoil_potential (f, phi, add = true);
  else {
    phi = new scalar;
    recoil_potential (f, phi, add = false);
    f.phi = phi;
  }
}

/**
## References

~~~bib
@article{mialhe2023extended,
  title={An extended model for the direct numerical simulation of droplet evaporation. Influence of the Marangoni convection on Leidenfrost droplet},
  author={Mialhe, Guillaume and Tanguy, S{\'e}bastien and Tranier, L{\'e}o and Popescu, Elena-Roxana and Legendre, Dominique},
  journal={Journal of Computational Physics},
  volume={491},
  pages={112366},
  year={2023},
  publisher={Elsevier}
}
~~~
*/
