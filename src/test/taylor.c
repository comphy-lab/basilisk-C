/**
# Equilibrium of a droplet suspended in an electric field 

A droplet suspended in a fluid subjected to a uniform electric field
deforms due to the competing effect of electrical forces and surface
tension. If the electrification level is moderate an equilibrium shape
is reached. It was observed experimentally that droplets deforms as
prolate or oblate spheroids (i.e. larger elongation aligned with the
external electric field or viceversa). The analytical analysis of
[O'Konski \& Thacher (1953)](#konski1953) unexpectedly predicted
prolate forms whereas the experiments showed oblate ones (and
viceversa).

It was the genius of [Geoffrey Taylor (1966)](#taylor1966) who shed
light on the problem. The work of O'Konski \& Thacher assumed that
both fluids (inner and outer) were perfect dielectrics, given the low
conductivity of the fluids involved. Taylor realized that the
conductivity could be very low but it was not zero, so that the
charges could migrate through the "leaky" media, thus accumulating at
the fluid interface and altering radically the pattern of electrical
forces.

Moreover, approximating the electrostatic Maxwell equations with the
Taylor--Melcher leaky dielectric model and assuming small deformation
Taylor predicted recirculating velocities and provided analytical
expressions for the radial and azimuthal velocities as functions of
the dimensionless radius $r$ and the azimuthal coordinate $\theta$.
For $r < 1$ this gives,
$$
v_r =  Ar(1-r^2)(3\cos^2 \theta -1)  \quad \mathrm{and} \quad  
v_\theta = 3Ar/2 \left( 1-5r^2/3 \right)\cos 2\theta
$$
and for $r \ge 1$,
$$
v_r = A(r^{-4}-r^{-2})(3\cos^2 \theta -1) \quad \mathrm{and} \quad
v_\theta = -Ar^{-4} \sin 2\theta
$$
with
$$
A = -\frac{9}{10}  \frac{E_\infty^2 R_d \varepsilon_o}{
\mu_o} \frac{1}{1+\lambda}\frac{R-Q}{(R+2)^2},
$$
where $R$, $Q$ and $\lambda$ are the ratio of inner to outer
conductivity, permittivity and viscosity, respectively. $E_\infty$ is
the imposed electrid field, $R_d$ the droplet radius, $\varepsilon_o$
is the outer permittivity and $\mu_2$ is the outer viscosity. The
electrical forces induces recirculations in both (viscous) fluids.

This test case is also discussed in 
[Lopez-Herrera et al, 2011](/src/references.bib#lopez-herrera2011).

The problem is assumed to be axisymmetric. */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "ehd/implicit.h"
#include "ehd/stress.h"
#include "vof.h"
#include "tension.h"

/**
We need to track the interface with the volume fraction field *f*. The
viscosity is constant but the coefficients will vary due to the
axisymmetric metric terms. */

scalar f[], * interfaces = {f};
face vector muv[];

/**
The maximum level of resolution, LEVEL, will be varied from 8 to
10. */

int LEVEL = 8;

#define Ef 1.34 // External electric field
#define R0 0.1 // Radius of the droplet 
#define F 50. 
#define R 5.1 // Conductivity ratio
#define Q 10.0 // permittivity ratio
#define CMU 0.1 // Outer viscosity 
#define theta (M_PI/4.)
#define LAM 1. // Viscosity ratio
#define VC (sq(Ef)*R0/CMU) // characteristic velocity
#define A (-9./10.*(R - Q)/sq(R + 2.)/(1. + LAM))

/**
F is the conductivity of the outer medium. F has no influence on the
steady solution but decreases the characteristic electrical relaxation
time and consequently the electrical transient. */

#define cond(T) (F*((1. - (T)) + R*(T)))
#define perm(T) ((1. - (T)) + Q*(T))

/**
The electric potential is linear. */

phi[top]   = dirichlet(Ef*x);
phi[left]  = dirichlet(Ef*x);
phi[right] = dirichlet(Ef*x);

/**
We make sure there is no flow through the boundaries, otherwise the
compatibility condition for the Poisson equation can be violated. */

uf.n[left] = 0.;
uf.n[right] = 0.;
uf.n[top] = 0.;
uf.n[bottom] = 0.;

/**
The domain spans [0:2]. We will compute only a quarter of the droplet,
making use of axisymmetry and right-left symmetry. The surface tension
coefficient is unity. The viscosity coefficients are variable. */

int main() {

  /**
  Space and time are dimensionless. */

  L0 = 2[0];
  DT = HUGE[0];
  
  N = 1 << LEVEL;
  f.sigma = 1.;
  mu = muv;
  run(); 
}

event init (t = 0) {
  
  /**
  We initialize the volume fraction field corresponding to a circular
  interface of radius R0. */

  fraction (f, sq(R0) - sq(x) - sq(y));

  /**
  We initialize the electrical potential. */

  foreach()
    phi[] = Ef*x;
}


/**
Permittivity and electrical conductivity are face values and also
incorporate the metric factors. The viscosity is constant but the
viscosity coefficients need to incorporate the metric factors. */

event properties (i++)
{
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    epsilon.x[] = perm(ff)*fm.x[];
    K.x[] = cond(ff)*fm.x[];
    muv.x[] = CMU*fm.x[];
  }
}

/**
## Convergence

We store the horizontal component of the velocity to check its
convergence with time. */

scalar un[];

event init_un (i = 0) {
  foreach()
    un[] = u.x[];
}

event error (i += 20; t <= 10.) {

  /**
  We monitor the variation in the horizontal component of the velocity
  and the convergence of the multigrid solvers every 20 timesteps. */

  double du = change (u.x, un);
  fprintf (stdout, "%g %g %d %d %d %d %d %d %d %d\n", t, du,
	   mgp.i, mgp.nrelax, mgpf.i, mgpf.nrelax, mgu.i, mgu.nrelax,
	   mgphi.i, mgphi.nrelax);
  fflush (stdout);

  /**
  If the change is small enough (i.e. the solution has converged for
  this level of refinement), we increase the level of refinement. If
  the simulation has converged and the level of refinement is 10, we
  stop the simulation. */
  
  if (i > 0 && du < 1e-5)
    return (LEVEL++ == 10); /* stop */
}

/**
## Results

At the end of the simulation we create two files: *log* (standard
error) will contain the dimensionless radial and azimuthal velocities
and their theoretical values as functions of the dimensionless radial
coordinate along the line $\theta= 45^o$. *vector.svg* displays the
velocity field, interface position and isopotential lines, as
displayed by gfsview-batch. */

event result (t = end) {
  double h  = 0.35*L0/(2*99);
  for (int i = 1; i <= 100; i++) {
    double x = i*h, y = i*h, r = sqrt(sq(x) + sq(y))/R0;
    double ux = interpolate (u.x, x, y)/VC; // dimensionless velocities
    double uy = interpolate (u.y, x, y)/VC;
    double vrt, vtt; // theoretical radial and azimuthal velocities;
    if (r < 1.) {
      vrt = A*r*(1. - sq(r))*(3.*sq(sin(theta)) - 1.);
      vtt = 3*A/2*r*(1. - 5./3.*sq(r))*sin(2.*theta);
    }
    else {
      vrt = A/sq(r)*(1/sq(r) - 1.)*(3.*sq(sin(theta)) - 1.);
      vtt = - A*1./sq(sq(r))*sin(2.*theta);
    }
    fprintf (stderr, "%g %g %g %g %g\n", r, 
	     (ux*x + uy*y)/(R0*r), vrt, (-uy*x + ux*y)/(R0*r), vtt);
  }

  FILE * fp = popen ("gfsview-batch2D taylor.gfv", "w");
  output_gfs (fp);
  fprintf (fp, "Save vectors.svg { format = SVG }\n");
  pclose (fp);
}

/**
The mesh is adapted according to interpolation errors on the volume
fraction, charge density and velocity fields. */

event adapt (i += 20) {
  adapt_wavelet ({f, rhoe, u.x, u.y}, (double[]){1e-3, 1, 2e-4, 2e-4},
		 maxlevel = LEVEL);
}

/**
~~~gnuplot Radial profiles of radial azimuthal velocities compared to analytical results.
set xlabel 'r'
set ylabel 'v'
plot 'log' u 1:2 notitle, 'log' u 1:3 w l t "v_r",\
     'log' u 1:4 notitle, 'log' u 1:5 w l t "v_{/Symbol q}"
~~~

![Steady-state velocity vectors, interface position and equipotential lines.](taylor/vectors.svg)

## Bibliography

~~~bib
@article{konski1953,
  title={The distortion of aerosol droplets by an electric field},
  author={O'Konski, Chester T and Thacher Jr, Henry C},
  journal={The Journal of Physical Chemistry},
  volume={57},
  number={9},
  pages={955--958},
  year={1953},
  publisher={ACS Publications}
}

@article{taylor1966,
  title={Studies in electrohydrodynamics. I. The circulation produced 
         in a drop by an electric field},
  author={Taylor, Geoffrey Ingram},
  journal={Proc. R. Soc. Lond. A},
  volume={291},
  number={1425},
  pages={159--166},
  year={1966},
  publisher={The Royal Society}
}
~~~

## See also

* [Same test with 
   Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/electro.html)
*/
