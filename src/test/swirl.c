/**
# Boundary layer on a rotating disk

[Von Kármán, 1921](#karman1921) showed that the steady flow of an
incompressible liquid of kinematic viscosity $\nu$ induced by an
infinite plane disk rotating at angular velocity $\Omega$ can be
described by a similarity solution. In effect, using
$\zeta=z\sqrt{\Omega/\nu}$ and setting the axial velocity $U$,
radial velocity $V$ and azimuthal velocity $W$ as
$$
U=\sqrt{\nu \Omega} F(\zeta) \quad V=\Omega r H(\zeta) 
\quad W = \Omega r G(\zeta)
$$
the Navier-Stokes equations reduce to a couple of ODEs:
$$
F'''-F\, F'' +F'^2/2 +2G^2  = 0 \quad \mathrm{and} \quad G''-F\,G'+G\,F' = 0
$$
with boundary conditions
$$
F(0)=F'(0)=0 \: G(0)=1. \quad \mathrm{and} \quad F'(\infty)=G(\infty)=0.
$$
where the prime denotes differentiation with respect to $\zeta$.

To reproduce this solution numerically, we use the axisymmetric
Navier--Stokes solver with azimuthal velocity (swirl). */

#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/swirl.h"

/**
The left boundary is the rotating disk with $\Omega = 1$ and a no-slip
condition for the tangential velocity i.e. */
  
u.t[left] = dirichlet(0);
w[left]   = dirichlet(y);

/**
We use an open (outflow) boundary condition for the right boundary. */

u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

/**
The top boundary condition is more tricky but the following seems to
work. */

u.n[top] = neumann(0);
p[top] = neumann(0);

/**
We use a constant viscosity but it needs to be weighted by the
(axisymmetric) metric. */

face vector muv[];

event properties (i++) {
  foreach_face()
    muv.x[] = 0.2*fm.x[];
}

/**
The computational domain is $12\times 12$ and we limit the
timestep. */

int main()
{
  size (12);
  N = 128;
  mu = muv;
  DT = 2e-2;
  run();
}

/**
We wait until the boundary layer is fully developed and
quasi-stationary. We only consider values close to the origin to
minimize the influence of boundaries (von Kármán's solution is valid
in an infinite domain). */

event end (t = 20)
{
  foreach()
    if (x*x + y*y < 8)
      fprintf (stderr, "%g %g %.4g %.4g\n", x, y, u.x[], w[]);
}

/**
~~~gnuplot Axial $F(\zeta)$ and azimuthal $G(\zeta)$ dimensionless velocity components
set xlabel '{/Symbol z}'
set key center right
nu = 0.2
plot [0:6]'analytical' u 1:2 w l t '-F({/Symbol z})',	     \
          'log' u ($1/sqrt(nu)):(-$3/sqrt(nu)) t 'Basilisk', \
          'analytical' u 1:3 w l t 'G({/Symbol z})',	     \
	  'log' u ($1/sqrt(nu)):($4/$2) t 'Basilisk'
~~~

## References

~~~bib
@article{karman1921,
  title={{\"U}ber laminare und turbulente Reibung},
  author={Karman, Th V},
  journal={ZAMM-Journal of Applied Mathematics and Mechanics/Zeitschrift 
           f{\"u}r Angewandte Mathematik und Mechanik},
  volume={1},
  number={4},
  pages={233--252},
  year={1921},
  publisher={Wiley Online Library}
}
~~~

## See also

* [Same case with 
   Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/swirl.html)
*/
