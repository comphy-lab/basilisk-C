/**
# Isopycnal gyres intersecting the bathymetry

This test verifies that isopycnals intersecting the bottom topography
do not lead to spurious imbalance. The setup is close, but not
identical, to the "two-layer model, high-cut topography" case of
section 3 of [Bleck & Smith, 1990](#bleck1990).

The velocity/height field can be compared to Fig.3 of [Bleck & Smith,
1990](#bleck1990).

![Velocity field and surface height](bleck/gyre.png)

At convergence the bottom layer must be stationary.

~~~gnuplot Evolution of the kinetic enery per unit area for both layers
set xlabel 'Time (years)'
set ylabel 'Kinetic energy per unit area (J/m^2)'
set logscale y
year = 86400.*365.
set grid
set key center right
plot [0:5]						\
   'clog' u ($1/year):3 w l t 'top layer',		\
   'clog' u ($1/year):2 w l t 'bottom layer'
~~~

## References

~~~bib
@article{bleck1990,
  title={A wind-driven isopycnic coordinate model of the north and 
  equatorial {A}tlantic {O}cean: 1. {M}odel development and supporting 
  experiments},
  author={Bleck, Rainer and Smith, Linda T},
  journal={Journal of Geophysical Research: Oceans},
  volume={95},
  number={C3},
  pages={3273--3285},
  year={1990},
  publisher={Wiley Online Library},
  doi={10.1029/JC095iC03p03273}
}
~~~
*/

#define hour 3600.
#define day (24.*hour)
#define year (365.*day)

/**
This flag changes the way face heights are computed for isopycnals
intersecting the bathymetry (in [hydro.h](/src/layered/hydro.h)). */

#include "layered/hydro.h"
#include "layered/implicit.h"

/**
## Coriolis

We include Coriolis acceleration on a $\beta$-plane. */

double Beta = 2e-11;
#define F0() (0.83e-4 + Beta*y)
#define alpha_H 0.5 // fixme: unstable for alpha_H 1
#include "layered/coriolis.h"

/**
## Vertical stratification

There are only two (isopycnal) layers. The bottom layer has a relative
density difference of 10^-2^. */

double * drho = (double []){ 1e-2, 0. };
#include "layered/isopycnal.h"

/**
## Surface "wind stress"

The "surface wind stress" is added directly as a (depth-weighted)
acceleration (along the x-direction) in the top layer. */

double tau0 = 1e-4;
  
event acceleration (i++)
{
  foreach_face(x)
    ha.x[0,0,1] += hf.x[] > dry ? tau0*cos(2.*pi*y/L0) : 0.;
}

/**
## Horizontal viscosity */

double nu_H = 1e5;

event acceleration (i++)
{
  if (nu_H > 0.) {
    vector d2u[];
    foreach_layer() {
      foreach()
	foreach_dimension()
	  d2u.x[] = (h[1]*u.x[1] + h[-1]*u.x[-1] +
		     h[0,1]*u.x[0,1] + h[0,-1]*u.x[0,-1] - 4.*h[]*u.x[])/sq(Delta);
      foreach()
	if (h[] > dry)
	  foreach_dimension()
	    u.x[] += dt*nu_H*d2u.x[]/max(h[],10.);
    }
  }
}

/**
## main */

int main()
{
  G = 9.8;
  nl = 2;
  size (6000e3 [1]);
  origin (0, - L0/2.);
  DT = 3000 [0,1];
#if TREE
  N = 16;
#else
  N = 32;
#endif
  theta_H = 0.51;
  run();
}

/**
## Initial conditions */

double sum0 = 0., sum1 = 0.;

event init (i = 0)
{
#if TREE
  refine (x < L0/4. && level < 5);
#endif
  /**
  No-slip (and dry) conditions on all boundaries. Note that this is
  only relevant for "wet" boundaries. For example, given that there is
  a (dry) coastline at x = 4400 km (see below), the `u.t[right]`
  boundary condition below is useless. */

  foreach_dimension() {
    u.t[left] = dirichlet(0);
    zb[left] = 1000.;
    h[left] = 0.;
    u.t[right] = dirichlet(0);
    zb[right] = 1000.;
    h[right] = 0.;
  }

  foreach() {

    /**
    The bathymetry. */
    
#if 0
    // low-cut case
    zb[] = x > 1400e3 ? - 4000. : - 1550. - 2450.*x/1400e3;
#elif 1
    // high-cut case
    zb[] = x > 2000e3 ? - 4000. : - 200. - 3800.*x/2000e3;
#else // flat bottom
    zb[] = -4000;
#endif

    /**
    A "coastline" at x = 4400 km. */
    
    if (x > 4400e3)
      zb[] = 1000.;

    /**
    The bottom layer is between the bottom and - 1000 m, the top layer
    is between 0 and -1000 m. */
    
    h[] = max (- 1000. - zb[], 0.);
    h[0,0,1] = max (- zb[] - h[], 0.);
  }

  sum0 = statsf(h).sum;
  sum1 = statsf(lookup_field("h1")).sum;
}

/**
## Outputs */

#if !_GPU // fixme: GPUs are not compatible with view.h yet
#include "view.h"
#endif

event end (t = 5*year)
{
#if !_GPU
  view (fov = 19.7885, tx = -0.35, ty = 0, width = 440, height = 600);
  squares ("zb < 100 ? eta : nodata", spread = -1);
  vectors ("u1", scale = 8e5);
#if TREE
  save ("gyre-tree.png");
#else
  save ("gyre.png");
#endif
#endif
}

event logfile (i += 300)
{
  double ke0 = 0., ke1 = 0., area = 0., rho0 = 1000 [0];
  scalar h1 = lookup_field ("h1");
  vector u1 = lookup_vector ("u1");
  foreach (reduction (+:ke0) reduction(+:ke1) reduction(+:area)) {
    ke0 += dv()*h[]*(sq(u.x[]) + sq(u.y[]))/2.;
    ke1 += dv()*h1[]*(sq(u1.x[]) + sq(u1.y[]))/2.;
    area += dv();
  }
  fprintf (stderr, "%g %g %g %g %g %g %d %d\n", t,
	   rho0*ke0/area, rho0*ke1/area,
	   (statsf (h).sum - sum0)/sq(L0), (statsf (h1).sum - sum1)/sq(L0),
	   dt, mgH.i, mgH.nrelax);
}
