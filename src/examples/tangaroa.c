/**
# Two-phase flow around RV Tangaroa

This is an improved version of the [famous Gerris
example](http://gerris.dalembert.upmc.fr/gerris/examples/examples/tangaroa.html),
illustrating the combination of complex solid boundaries, air-water
turbulent flows and reduced gravity approach. 

We use the centered Navier--Stokes solver, two-phase flow and the
momentum-conserving option. Note that the momentum-conserving option
is crucial to obtain stable solutions for this air-water density ratio
configuration. */

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"

/**
We also need to compute distance functions (to describe the ship
geometry), use reduced gravity and visualisation functions. */

#include "distance.h"
#include "reduced.h"
#include "view.h"
#include "lambda2.h"

/**
On supercomputers we need to control the maximum runtime and we check
performances. */

#include "maxruntime.h"
#include "navier-stokes/perfs.h"

/**
## Importing the geometry 

This function computes the solid fraction given a pointer to an STL
file, a tolerance (maximum relative error on distance) and a 
maximum level. */

void fraction_from_stl (scalar f, FILE * fp, double eps, int maxlevel)
{

  /**
  We read the STL file and compute the bounding box of the model. */
  
  coord * p = input_stl (fp);
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  
  /**
  We initialize the distance field on the coarse initial mesh and
  refine it adaptively until the threshold error (on distance) is
  reached. */

  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){eps*maxl}, maxlevel, 5).nf);

  /**
  We also compute the volume fraction from the distance field. We
  first construct a vertex field interpolated from the centered field
  and then call the appropriate VOF functions. */

  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  fractions (phi, f);
}

/**
## Main function

We can change both the maximum level of refinement and the [Froude
number](https://en.wikipedia.org/wiki/Froude_number) at runtime. 

[RV Tangaroa](https://en.wikipedia.org/wiki/RV_Tangaroa) is 70 metres
long. If we assume that it moves at 20 knots (twice its actual cruise
speed), this gives a Froude number of approx 0.4. */

int LEVEL = 9;
double FROUDE = 0.4, Uship = 1., Lship = 1. [1];

/**
We need additional (fraction) fields for the ship geometry and for the
(inflow) boundary condition. */

scalar tangaroa[], f0[];

int main (int argc, char * argv[])
{
  maxruntime (&argc, argv);
  if (argc > 1)
    LEVEL = atoi(argv[1]);
  if (argc > 2)
    FROUDE = atof(argv[2]);
    
  init_grid (32);

  rho1 = 1. [0];    // water
  rho2 = rho1/815.; // air

  /**
  The domain is five times larger than the ship. We change the origin
  so that the ship is not too close to the inflow. */
	      
  size (5.*Lship);
  origin (-L0/2.,-L0/3.,-L0/2.);

  /**
  We need to tell the code that both `tangaroa` and `f0` are volume
  fraction fields. */
	      
  for (scalar s in {tangaroa,f0})
    s.refine = s.prolongation = fraction_refine;

  /**
  The acceleration of gravity is given by the Froude number i.e. */

  G.z = - sq(Uship/FROUDE)/Lship;
  run();
}

/**
## Boundary conditions

The inflow condition fixes the velocity (unity) and the water level
(using `f0`). */
	      
u.n[bottom] = dirichlet(Uship);
p[bottom]   = neumann(0.);
pf[bottom]  = neumann(0.);
f[bottom]   = f0[];

/**
Outflow uses standard Neumann/Dirichlet conditions.  */
	      
u.n[top]  = neumann(0.);
p[top]    = dirichlet(0.);
pf[top]   = dirichlet(0.);

/**
Boundary conditions for the solid and fraction tracers. */

tangaroa[back] = 0;
f[back] = 1;

/**
Not sure whether this is really useful. */
	      
uf.n[left] = 0.;
uf.n[right] = 0.;
	      
/**
## Initial conditions

We can optionally restart, otherwise we open the STL file and
initialize the corresponding fraction. We also initialize the `f0`
field used for the inflow condition and set the initial water level
and velocity field. */      

event init (t = 0) {
  if (!restore (file = "restart")) {
    FILE * fp = fopen ("tangaroa.stl", "r");
    fraction_from_stl (tangaroa, fp, 5e-4[0], LEVEL);
    fclose (fp);
    
    fraction (f0, - z);

    foreach() {
      f[] = f0[];
      u.y[] = Uship;
    }
  }

  /**
  We setup the default display. */

  display ("view (quat = {0.542,0.150,0.209,0.799},"
	   "      fov = 30, near = 0.01, far = 1000,"
	   "      tx = -0.0074, ty = -0.0138, tz = -0.3211);"
	   "draw_vof (c = 'f', color = 'z', min = -0.1, max = 0.1, "
	   "          linear = true, fc = {0.45,0.78,0.92});"
	   "draw_vof (c = 'tangaroa');", true);
}

/**
## Boundary conditions on the ship

We use a simple (but crude) imposition of $u=0$ in the solid. */
	      
event velocity (i++) {
  foreach()
    foreach_dimension()
      u.x[] = (1. - tangaroa[])*u.x[];
}

/**
## Animations

We generate animations of the ship surface (as represented by the
solid fraction) and of the air-water interface, colored with the
height.

Several classical features of ship wakes are recognisable: breaking
bow wave, breaking stern divergent wave, turbulent boundary layer,
Kelvin waves etc...

![Evolution of the air-water interface](tangaroa/movie.mp4)(width="800" height="600")

We also use the $\lambda_2$ criterion to display the turbulent
vortical structures generated by the airflow. The air recirculation at
the top of the steep primary Kelvin waves is particularly noticeable.

![Turbulent vortical structures](tangaroa/l2.mp4)(width="800" height="600")

The computations above were done on the Irene supercomputer using 12
levels of refinement. */
	      
event movie (t += 0.01; t <= 10)
{
  view (fov = 5.86528,
	quat = {0.515965,0.140691,0.245247,0.808605},
	tx = -0.07438, ty = -0.0612925,
	width = 1024, height = 768);

  clear();
  draw_vof ("tangaroa");
  scalar Z[];
  Z[back] = dirichlet (z);
  foreach()
    Z[] = z;
  draw_vof ("f", color = "Z", min = -0.1, max = 0.1, linear = true);
  save ("movie.mp4");

  draw_vof ("tangaroa", fc = {0.5,0.5,0.5});
  draw_vof ("f", color = "Z", min = -0.1, max = 0.1, linear = true);
  scalar l2[];
  lambda2 (u, l2);
  isosurface ("l2", -100);
  save ("l2.mp4");
}

#if DUMP
event snapshot (i += 100)
{
  char name[80];
  sprintf (name, "dump-%d", i);
  scalar l2[];
  lambda2 (u, l2);
  dump (file = name);
}
#endif

/**
## Mesh adaptation

This computation is only feasible thanks to mesh adaptation, based
both on volume fraction and velocity accuracy. */
	      
event adapt (i++) {
  double uemax = 0.1;
  adapt_wavelet ({f,tangaroa,u},
		 (double[]){0.01,0.01,uemax,uemax,uemax}, LEVEL, 5);
}

/**
## Running in parallel on Irene

Running with MPI-parallelism is a bit more complicated than usual
since the `distance()` function is not parallelised yet. A reasonably
simple workaround is to first generate a restart/dump file on the
local machine, without MPI, using something like:

~~~bash
CFLAGS=-DDUMP=1 make tangaroa.tst
~~~

(and also adjust the maximum level), then kill the running code (using
Ctrl-C) and do:

~~~bash
qcc -source -D_MPI=1 tangaroa.c
scp _tangaroa.c popinets@irene.ccc.cea.fr:/ccc/scratch/cont003/gen7325/popinets
scp tangaroa/dump-0 popinets@irene.ccc.cea.fr:/ccc/scratch/cont003/gen7325/popinets/restart
~~~

then on irene (to run on 480 cores for 10 hours, with 12 levels of refinement):

~~~bash
ssh popinets@irene.ccc.cea.fr
ccc_mpp -u popinets
sed 's/WALLTIME/36000/g' run.sh | ccc_msub -n 480
~~~

with the following `run.sh` script

~~~bash
#!/bin/bash
#MSUB -r tangaroa
#MSUB -T WALLTIME
#MSUB -@ popinet@basilisk.fr:begin,end
#MSUB -o basilisk_%I.out
#MSUB -e basilisk_%I.log
#MSUB -q skylake
#MSUB -A gen7760
#MSUB -m scratch
##MSUB -w

set -x
cd ${BRIDGE_MSUB_PWD}

mpicc -Wall -std=c99 -D_XOPEN_SOURCE=700 -O2 _tangaroa.c -o tangaroa \
      -L$HOME/gl -lglutils -lfb_tiny -lm
ccc_mprun -n ${BRIDGE_MSUB_NPROC} ./tangaroa -m WALLTIME 12 0.4 \
    2>> log-${BRIDGE_MSUB_NPROC} >> out-${BRIDGE_MSUB_NPROC}
~~~

### Generating MP4 movies

Note that when running on Irene the [ffmpeg](https://www.ffmpeg.org/)
MP4 encoder is not available and the logfile will contain the warning:

~~~bash
open_image(): cannot find 'ppm2mp4' or 'ffmpeg'/'avconv'
  falling back to raw PPM outputs
~~~

The MP4 files defined above will be renamed `movie.mp4.ppm` and
`l2.mp4.ppm`. As the extension indicates, these (large) files are now
raw (uncompressed) PPM images. To convert them to compressed MP4, you
will need to copy them to a machine where ffmpeg (and Basilisk) are
installed (i.e. your local machine) and do:

~~~bash
ppm2mp4 movie.mp4 < movie.mp4.ppm
ppm2mp4 l2.mp4 < l2.mp4.ppm
~~~
*/
