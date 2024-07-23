/**
# Viscoelastic 2D drop in a Couette Newtonian shear flow

This problem has been used as benchmark (see for example [Chinyoka et
al. (2005)](#chinyoka2005)). Equations are usually made dimensionless 
with the outer density $\rho_2$, the droplet radius $a$ and the shear rate
$\dot{\gamma}$. The following dimensionless parameters appear
$$ 
We = \frac{\rho_2 a^3 \dot{\gamma}^2}{\sigma}, \quad
Re = \frac{\rho_2 a^2 \dot{\gamma}}{\mu_2}, \quad
\mu_r =\frac{\mu_1}{\mu_2}, \quad
m = \frac{\rho_1}{\rho_2}, \quad \beta = \frac{\mu_s}{\mu_1}, \quad
De = \dot{\gamma}\lambda
$$ 
where subscript "1" and "2" stand for the droplet and matrix fluid
respectively and $\mu_s$ is the solvent viscosity. Note that the
polymeric viscosity is $\mu_p=\mu_1-\mu_s$. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "log-conform.h"
#include "tension.h"

#define Ca 0.6     // Capillary number
#define Re 0.3     // Reynold number
#define We (Ca*Re) // Weber number
#define MUr 1.     // ratio of outer(matrix) to inner(drop) viscosity
#define M 1.       // ratio of outer to inner density
#define Deb 0.4    // Deborah number
#define Beta 0.5   // ratio of the solvent visc. to the total viscoelastic visc.

/**
We set a maximum level of 8 only so that the test case runs in less
than 30 minutes but note that the results presented below are for
MAXLEVEL = 9. */

int MAXLEVEL = 8;

scalar mupd[], lam[];

/**
The top and bottom boundary conditions are those of a Couette flow. */

u.t[top] = dirichlet (y);
u.t[bottom] = dirichlet (y);

/** 
The domain is a 16x16 box which will be masked later to a
bottom-left corner of the domain of coordinates ($x=-8, y=-4$). We set
a maximum timestep of 0.1 */

int main() {
  L0 = 16 ;
  origin (-L0/2, -L0/4.);
  periodic (right);
  DT = .1;

  /**
  We set the viscosities, densities, surface tension and visco-elastic
  parameters according to the analysis above. */
  
  mu1 = MUr*Beta/Re;
  mu2 = 1./Re;
  rho1 = M;
  rho2 = 1.;
  f.sigma = 1./We;
  lambda = lam;
  mup = mupd;

  init_grid (1 << MAXLEVEL);
  run();
}

event init (i = 0) {

  /**
  We mask above $y > 4$. The computational domain is now a 16x8
  rectangle with the origin in the center. */

  mask (y > 4 ? top : none);

  /**
  As initial conditions we set a viscoelastic droplet of radius 1
  and a linear velocity profile typical of the Couette flow. */

  fraction (f, 1. - (sq(x) + sq(y)));
  foreach()
    u.x[] = y;
}

/**
The event below set viscoelastic properties at each step. The
dimensionless relaxation parameter is the Deborah number, $De$, while
the dimensionless polymeric viscosity $\mu_p$ is
$$
\frac{\mu_p}{\rho_2 a^2 \dot{\gamma}} =
\frac{\mu_1 - \beta \mu_1}{\rho_2 a^2 \dot{\gamma}} =
\frac{1-\beta}{Re}\frac{\mu_1}{\mu_2}
$$
*/

event properties (i++) {
  foreach () {
    lam[] = Deb*f[];
    mupd[] = MUr*(1. - Beta)*f[]/Re;
  }
}

/**
The mesh is adapted according to the errors on volume fraction and
velocity. */

event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){1e-2, 1e-3, 1e-3}, MAXLEVEL, MAXLEVEL - 2);
}

/**
As outputs we plot the shape of the interface at instant t = 10 and the
time evolution of the deformation parameter. */

event interface (t = 10.) {
  output_facets (f);
}

event deformation (t += 0.1) {
  double rmax = -HUGE, rmin = HUGE ;
  foreach (reduction(max:rmax) reduction(min:rmin)) 
    if (f[] > 0 && f[] < 1) {
      coord p;
      coord n = mycs (point, f);
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      double rad  = sqrt(sq(x + Delta*p.x) + sq(y + Delta*p.y)); 
      if (rad > rmax)
	rmax = rad;
      if (rad < rmin)
	rmin = rad;
    }
  double D = (rmax - rmin)/(rmax + rmin);
  fprintf (stderr, "%g %g %g %g\n", t, rmin, rmax, D);
}

/**
We can optionally visualise the results with Basilisk View while we
run. */

#if 0
#include "view.h"

event viewing (i += 10) {
  static FILE * fp = popen ("bppm","w");
  view (width = 600, height = 300, fov = 10);
  clear();
  draw_vof ("f");
  // squares ("u.y", linear = true);
  squares ("level");
  save (fp = fp);
}
#endif

/**
We compare the results to those of [Figueiredo et al. (2016)](#figueiredo2016).

~~~gnuplot Time evolution of the deformation
set xlabel 't'
set ylabel 'D'
set key bottom Right
plot 'viscodrop.figueiredo' u 1:2 pt 7 t 'Figueiredo et al. (2016)', \
     'viscodrop.log-9' u 1:4 w l lw 2 t 'Basilisk'
~~~

~~~gnuplot Interface shape at t = 10
set xlabel 'x'
set ylabel 'y'
plot 'viscodrop.interface' u 1:2 pt 7 t 'Figueiredo et al. (2016)', \
     'viscodrop.out-9' u 1:2 w l lw 2 t 'Basilisk'
~~~

## References

~~~bib
@article{chinyoka2005,
  title={Two-dimensional study of drop deformation under 
  simple shear for {O}ldroyd-{B} liquids},
  author={Chinyoka, T and Renardy, YY and Renardy, M and Khismatullin, DB},
  journal={Journal of Non-Newtonian Fluid Mechanics},
  volume={130},
  number={1},
  pages={45--56},
  year={2005},
  publisher={Elsevier}
}

@article{figueiredo2016,
  title={A two-phase solver for complex fluids: Studies of 
  the {W}eissenberg effect},
  author={Figueiredo, RA and Oishi, CM and Afonso, AM and Tasso, 
  IVM and Cuminato, JA},
  journal={International Journal of Multiphase Flow},
  volume={84},
  pages={98--115},
  year={2016},
  publisher={Elsevier}
}
~~~
*/
