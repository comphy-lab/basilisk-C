/**
# 3D breaking Stokes wave (multilayer solver)

A steep, 3D, third-order Stokes wave is unstable and breaks. This is
the 3D equivalent of this [test case](/src/test/stokes.c). The
bathymetry is given by
$$
z_b(y) = - 0.5 + \sin(\pi y)/4
$$
i.e. it is shallower toward the back of the domain which causes the
wave to break earlier there.

![Animation of the free-surface. The surface is coloured according to
 the $x$-component of the surface velocity.](breaking/movie.mp4)(
 width=100% )

The solution is obtained using the layered model and demonstrates its
robustness and a degree of realism even for this complex case. */

#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "view.h"

/**
The initial conditions are given by the wave steepness $ak$ and the
Reynolds number $Re=c\lambda/\nu$ with $c$ the phase speed of the
gravity wave and $\lambda$ its wavelength. */

double ak = 0.33;
double RE = 40000.;
double k_ = 2.*pi, h_ = 0.5, g_ = 1.;

#define T0  (k_*L0/sqrt(g_*k_))

/**
The domain is periodic in $x$ and resolved using 256$^2$
points and 30 layers. */

int main()
{
  origin (-L0/2., -L0/2.);
  periodic (right);
  N = 256;
  nl = 30;
  G = g_;
  nu = 1./RE;

  /**
  Some implicit damping is necessary to damp fast modes. This may be
  related to the slow/fast decoupling of the $\theta$-scheme mentioned
  by [Durran & Blossey, 2012](#durran2012). */
  
  theta_H = 0.51;

  run();
}

/**
The initial conditions for the free-surface and velocity are given by
the third-order Stokes solution. */

#include "../test/stokes.h"

event init (i = 0)
{

  /**
  We can use a larger CFL, in particular because we are not dealing
  with shallow-water/wetting/drying. */
  
  CFL = 0.8;

  /**
  The layer thicknesses follow a geometric progression, starting from
  a top layer with a thickness of 1/3 that of the regular
  distribution. */
  
  geometric_beta (1./3., true);

  double a = 0.25;
  foreach() {
    zb[] = - 0.5 + a*sin(k_/2.*y);
    double H = wave(x, 0) - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H*beta[point.l];
      z += h[]/2.;
      u.x[] = u_x(x, z);
      w[] = u_y(x, z);
      z += h[]/2.;
    }
  }
}

/**
We add (an approximation of) horizontal viscosity. */

event viscous_term (i++)
  horizontal_diffusion ((scalar *){u}, nu, dt);

/**
We log the evolution of the kinetic and potential energies.

~~~gnuplot Evolution of the kinetic, potential and total energy
set xlabel 't/T0'
plot [0:6]'log' u 1:2 w l t 'kinetic', '' u 1:3 w l t 'potential', \
     '' u 1:(($2+$3)/2.) w l t 'total/2'
~~~
*/

event logfile (i++; t <= 8.*T0)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
    }
    gpe += sq(eta[])*dv();
  }
  fprintf (stderr, "%g %g %g\n", t/T0, ke/2., g_*gpe/2.);
}

/**
And generate the movie of the free surface (this is quite
expensive). The movie is 45 seconds at 25 frames/second. */

event movie (t += 8.*T0/(45*25))
{
  view (fov = 17.3106, quat = {0.475152,0.161235,0.235565,0.832313},
	tx = -0.0221727, ty = -0.0140227, width = 1200, height = 768);
  char s[80];
  sprintf (s, "t = %.2f T0", t/T0);
  draw_string (s, size = 80);
  for (double x = -1; x <= 1; x++)
    translate (x) {
      squares ("u29.x", linear = true, z = "eta", min = -0.15, max = 0.6);
    }
  save ("movie.mp4");
}

/**
## Parallel run

The simulation was run in parallel on the [Occigen
machine](https://www.cines.fr/calcul/materiels/occigen/) on 64 cores,
using this script

~~~bash
local% qcc -source -D_MPI=1 breaking.c
local% scp _breaking.c user@occigen.cines.fr:
~~~

~~~bash
#!/bin/bash
#SBATCH -J breaking
#SBATCH --constraint=HSW24
#SBATCH --ntasks=64
#SBATCH --threads-per-core=1
#SBATCH --time=1:00:00
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=popinet@basilisk.fr

module purge
module load openmpi
module load intel gcc

NAME=breaking
mpicc -Wall -std=c99 -D_XOPEN_SOURCE=700 -O2 _$NAME.c -o $NAME \
    -L/home/popinet/local/gl -L/home/popinet/local/lib \
    -lglutils -lfb_tiny -lppr -lgfortran -lm

export LD_LIBRARY_PATH=/home/popinet/local/lib:$LD_LIBRARY_PATH
export PATH=$PATH:/home/popinet/local/bin
rm -f *.ppm
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS $NAME 2> log > out
~~~

The number of timesteps was 4159 and the runtime was 44 minutes with
movie generation.

## References

~~~bib
@article{durran2012,
  title={Implicit--explicit multistep methods for fast-wave--slow-wave problems},
  author={Durran, Dale R and Blossey, Peter N},
  journal={Monthly Weather Review},
  volume={140},
  number={4},
  pages={1307--1325},
  year={2012},
  publisher={American Meteorological Society},
  doi={10.1175/MWR-D-11-00088.1}
}
~~~
*/
