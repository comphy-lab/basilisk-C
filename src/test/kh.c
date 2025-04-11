/**
# Lock exchange (Kelvin--Helmoltz shear instability)

A rectangular horizontal box is filled with a fluid with a non-uniform
density. The left half of the domain is "warmer" and thus lighter than
the right half. The fluid is viscous and the density field is
diffusive with a Péclet number of unity (i.e. the kinematic viscosity
of the fluid is identical to the diffusivity of the density).

This is a similar test case to that proposed by [Berntsen et
al. 2006](#berntsen2006) to illustrate the importance of
non-hydrostatic effects. The main difference is the symmetric slip
boundary conditions on the top and bottom boundaries.

![Animation of the density field](kh/T.mp4)

The results can be directly compared with the solution obtained using
the (completely different) [centered Navier--Stokes
solver](kh-ns.c). The agreement is excellent as illustrated below.

~~~gnuplot Ten equally-spaced contours of density at $t=$ 21.
set term svg enhanced size 1000,170 font ",10"
set size ratio -1
unset key
unset ytics
set xrange [-4:4]
plot 'log' u 1:2 w l lc rgbcolor "black"
~~~

~~~gnuplot Results with the [Navier--Stokes solver](kh-ns.c).
plot '../kh-ns/log' u 1:2 w l lc rgbcolor "black"
~~~

The hydrostatic solution is different.

~~~gnuplot Hydrostatic solution.
plot '../kh-hydro/log' u 1:2 w l lc rgbcolor "black"
~~~

## References

~~~bib
@article{berntsen2006,
  title={Assessment of non-hydrostatic ocean models using 
         laboratory scale problems},
  author={Berntsen, Jarle and Xing, Jiuxing and Alendal, Guttorm},
  journal={Continental Shelf Research},
  volume={26},
  number={12-13},
  pages={1433--1447},
  year={2006},
  publisher={Elsevier},
  doi={10.1016/j.csr.2006.02.014}
}
~~~

## Code

We use an x-z grid and the non-hydrostatic multilayer solver. */

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#if HYDRO
# include "layered/implicit.h"
#else
# include "layered/nh.h"
#endif

/**
The relative density difference is $10^{-2}$. */

#define drho(T) (- 0.01*(T))
#include "layered/dr.h"

#include "layered/remap.h"
#include "layered/perfs.h"

/**
The domain is $8 \times 1$. */

int main()
{
  L0 = 8.;
  X0 = -L0/2.;
  G = 9.81;
  nl = 64;
  N = nl*L0;
  nu = 2e-4;

  theta_H = 0.51;
  cell_lim = mono_limit;
  DT = 0.03;

  /**
  We use a slip boundary condition at the bottom, rather than the
  default no-slip (the top boundary is slip by default). */

  lambda_b[] = {HUGE};
  
  run();

  /**
  We generate movie and graphics at the end of the run. */
  
  system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
	  "ppm2mp4 T.mp4");
#if HYDRO
  system ("gnuplot -e 'set table' kh-hydro.plot | sed '/^#.*/d' > log");
#else
  system ("gnuplot -e 'set table' kh.plot | sed '/^#.*/d' > log");
#endif
}

/**
The initial density profile is a smooth hyperbolic tangent function. */

event init (i = 0)
{
  foreach() {
    double z = zb[];
    foreach_layer() {
      h[] = 1./nl;
      z += h[]/2.;
      T[] = - 0.5*tanh((x + 0.1*cos(pi*z/2.))/0.04);
      z += h[]/2.;
    }
  }
}

/**
Horizontal viscosity and diffusion of density (with the same
viscosity/diffusion coefficient). */

event viscous_term (i++)
{
#if NH
  horizontal_diffusion ({u,w,T}, nu, dt);
#else
  horizontal_diffusion ({u,T}, nu, dt);
#endif // NH
  
  foreach()
    vertical_diffusion (point, h, T, dt, nu, 0, 0, HUGE);
}

/**
Density field at $t=$ 21. */

event density (t = 21)
{
  foreach() {
    double z = zb[];
    foreach_layer() {
      z += h[]/2.;
      printf ("%g %g %g\n", x, z, T[]);
      z += h[]/2.;
    }
    printf ("\n");
  }
  fflush (stdout);
}

/**
We generate an animation of the density field. */

event gnuplot (t += 1; t <= 30)
{
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0) {
    system ("rm -f plot*.png");
    fprintf (fp,
	     "set term png font ',9' size 1000,250\n"
	     "set pm3d map\n" // interpolate 10,4\n"
	     "set size ratio -1\n"
	     "unset key\n"
	     "unset colorbox\n"
	     "unset ytics\n"
	     "# jet colormap\n"
	     "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	     " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, "
	     "0.625 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0,"
	     " 1 0.498 0 0 )\n"
	     "set cbrange [-0.5:0.5]\n");
  }
  fprintf (fp,
	   "set output 'plot-%02g.png'\n"
	   "sp [%g:%g][0:1.001]'-' u 1:2:3\n",
	   t, X0, X0 + L0);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g\n", x, z, T[]);
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g\n", x, z, T[]);
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 0.05\n");
  fflush (fp);
}
