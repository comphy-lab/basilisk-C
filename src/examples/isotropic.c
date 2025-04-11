/**
# Forced isotropic turbulence in a triply-periodic box

We compute the evolution of forced isotropic turbulence (see [Rosales
& Meneveau, 2005](/src/references.bib#rosales2005)) and compare the
solution to that of the [hit3d](http://code.google.com/p/hit3d/)
pseudo-spectral code. The initial condition is an unstable solution to
the incompressible Euler equations. Numerical noise in the solution
eventually leads to the destabilisation of the base solution into a
fully turbulent flow where turbulent dissipation balances the linear
input of energy. */

#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"

/**
We use the $\lambda_2$ criterion and Basilisk View for visualisation
of vortices. */

#include "lambda2.h"
#include "view.h"

/**
We monitor performance statistics and control the maximum runtime. */

#include "navier-stokes/perfs.h"
#include "maxruntime.h"

const double MU = 0.01;

/**
We need to store the variable forcing term. */

face vector av[];

/**
The code takes the level of refinement as optional command-line
argument (as well as an optional maximum runtime). */

int maxlevel = 7;

int main (int argc, char * argv[])
{
  maxruntime (&argc, argv);
  if (argc > 1)
    maxlevel = atoi(argv[1]);

  /**
  The domain is $(2\pi)^3$ and triply-periodic. */
  
  L0 = 2.*pi [1];
  foreach_dimension()
    periodic (right);

  /**
  The viscosity is constant. The acceleration is defined below. The
  level of refinement is *maxlevel*. */

  mu[] = {MU,MU,MU};
  a = av;
  N = 1 << maxlevel;
  run();
}

/**
## Initial conditions

The initial condition is "ABC" flow. This is a laminar base flow that 
is easy to implement in both Basilisk and a spectral code. */

event init (i = 0) {
  double u0 = 1., k = 1.;
  if (!restore (file = "restart"))
    foreach() {
      u.x[] = u0*(cos(k*y) + sin(k*z));
      u.y[] = u0*(sin(k*x) + cos(k*z));
      u.z[] = u0*(cos(k*x) + sin(k*y));
    }
}

/**
## Linear forcing term

We compute the average velocity and add the corresponding linear
forcing term. */

event acceleration (i++) {
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  double a = 0.1;
  foreach_face()
    av.x[] += a*((u.x[] + u.x[-1])/2. - ubar.x);
}

/**
## Outputs

We log the evolution of viscous dissipation, kinetic energy, and
microscale Reynolds number. */

event logfile (i++; t <= 300) {
  coord ubar;
  foreach_dimension() {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
  
  double ke = 0., vd = 0., vol = 0.;
  foreach(reduction(+:ke) reduction(+:vd) reduction(+:vol)) {
    vol += dv();
    foreach_dimension() {
      // mean fluctuating kinetic energy
      ke += dv()*sq(u.x[] - ubar.x);
      // viscous dissipation
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
		  sq(u.x[0,1] - u.x[0,-1]) +
		  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  ke /= 2.*vol;
  vd *= MU/vol;

  if (i == 0)
    fprintf (stderr, "t dissipation energy Reynolds\n");
  fprintf (stderr, "%g %g %g %g\n",
	   t, vd, ke, 2./3.*ke/MU*sqrt(15.*MU/vd));
}

/**
We generate a movie of the vortices. 

![Animation of the $\lambda_2$ isosurface (a way to characterise
vortices) and cross-sections of velocity and vorticity.](isotropic/movie.mp4)
*/

event movie (t += 0.25; t <= 150)
{
  view (fov = 44, camera = "iso", ty = .2,
	width = 600, height = 600, bg = {1,1,1}, samples = 4);
  clear();
  squares ("u.y", linear = true);
  squares ("u.x", linear = true, n = {1,0,0});
  scalar omega[];
  vorticity (u, omega);
  squares ("omega", linear = true, n = {0,1,0});
  scalar l2[];
  lambda2 (u, l2);
  isosurface ("l2", -1);
  save ("movie.mp4");
}

/**
We can optionally try adaptivity. */

#if TREE
event adapt (i++) {
  double uemax = 0.2*normf(u.x).avg;
  adapt_wavelet ((scalar *){u}, (double[]){uemax,uemax,uemax}, maxlevel);
}
#endif

/**
## Running with MPI on occigen

On the local machine

~~~bash
%local qcc -source -D_MPI=1 isotropic.c
%local scp _isotropic.c occigen.cines.fr:
~~~

On occigen (using 512 cores)

~~~bash
module purge
module load openmpi
module load intel
mpicc -Wall -O2 -std=c99 -D_XOPEN_SOURCE=700 _isotropic.c -o isotropic \
        -I$HOME -L$HOME/gl -lglutils -lfb_tiny -lm
sed 's/WALLTIME/10:00/g' run.sh | sbatch
~~~

Note that this assumes that the Basilisk gl libraries have been
[installed](/src/gl/INSTALL#standalone-installation) in `$HOME/gl`.

With the `run.sh` script

~~~bash
#!/bin/bash
#SBATCH -J basilisk
#SBATCH --nodes=32
#SBATCH --constraint=HSW24
#SBATCH --ntasks-per-node=16
#SBATCH --threads-per-core=1
#SBATCH --time=WALLTIME
#SBATCH --output basilisk.output
#SBATCH --exclusive

LEVEL=7

module purge
module load openmpi
module load intel

srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS \
     ./isotropic -m WALLTIME $LEVEL \
     2> log-$LEVEL-$SLURM_NTASKS > out-$LEVEL-$SLURM_NTASKS
~~~

## Running with MPI on mesu

On the local machine

~~~bash
%local qcc -source -D_MPI=1 isotropic.c
%local scp _isotropic.c mesu.dsi.upmc.fr:
~~~

On mesu (using 512 cores)

~~~bash
module load mpt
mpicc -Wall -O2 -std=c99 -D_XOPEN_SOURCE=700 _isotropic.c -o isotropic \
      -L$HOME/gl -lglutils -lfb_tiny -lm
sed 's/WALLTIME/10:00/g' run.sh | qsub
~~~

with the `run.sh` script

~~~bash
#!/bin/bash 
#PBS -l select=22:ncpus=24:mpiprocs=24
#PBS -l walltime=WALLTIME
#PBS -N isotropic
#PBS -j oe  
# load modules 
module load mpt
# change to the directory where program job_script_file is located 
cd $PBS_O_WORKDIR 
# mpirun -np 64 !!!! does not work !!!!
NP=512
mpiexec_mpt -n $NP ./isotropic -m WALLTIME 2>> log.$NP >> out.$NP
~~~

## Results

The two codes agree at early time, or until the solution transitions
to a turbulent state. The statistics produced by the two codes agree
well after transition to turbulence.

~~~gnuplot Evolution of kinetic energy
set xlabel 'Time'
set ylabel 'Kinetic energy'
set logscale  y
plot 'isotropic.occigen' u 1:3 w l t 'Basilisk (occigen)', \
     'isotropic.mesu' u 1:3 w l t 'Basilisk (mesu)', \
     'isotropic.hit3d' u 1:($3*3./2.) w l t 'Spectral'
~~~

~~~gnuplot Evolution of microscale Reynolds number
set ylabel 'Microscale Reynolds number'
plot 'isotropic.occigen' u 1:4 w l t 'Basilisk (occigen)', \
     'isotropic.mesu' u 1:4 w l t 'Basilisk (mesu)', \
     'isotropic.hit3d' u 1:4 w l t 'Spectral'
~~~

~~~gnuplot Evolution of dissipation
set ylabel 'Dissipation function'
plot 'isotropic.occigen' u 1:2 w l t 'Basilisk (occigen)', \
     'isotropic.mesu' u 1:2 w l t 'Basilisk (mesu)', \
     'isotropic.hit3d' u 1:2 w l t 'Spectral'
~~~

The computational speed is respectable (for a relatively small 128^3^
problem on 512 cores). Note that these were obtained when switching
off movie outputs.

~~~gnuplot Computational speed in points.timesteps/sec/core
set ylabel 'Speed'
unset logscale
plot 'isotropic.occigen' u 1:($7/512) w l t 'occigen', \
     'isotropic.mesu' u 1:($7/512) w l t 'mesu'
~~~

## Scalability on irene

[Irene](http://www-hpc.cea.fr/en/complexe/tgcc-Irene.htm) is the
supercomputer at CEA.

On the local machine

~~~bash
%local qcc -source -D_MPI=1 isotropic.c
%local scp _isotropic.c irene.ccc.cea.fr:
~~~

On irene

~~~bash
mpicc -Wall -O2 -std=c99 -D_XOPEN_SOURCE=700 -xCORE-AVX512 \
      _isotropic.c -o isotropic \
      -L$HOME/gl -lglutils -lfb_tiny -lm
sed -e 's/WALLTIME/600/g' -e 's/LEVEL/7/g' run.sh | ccc_msub -n 512
~~~

with the `run.sh` script

~~~bash
#!/bin/bash
#MSUB -r isotropic
#MSUB -T WALLTIME
#MSUB -o basilisk_%I.out
#MSUB -e basilisk_%I.log
#MSUB -q skylake
#MSUB -A gen7760
#MSUB -m scratch

set -x
cd ${BRIDGE_MSUB_PWD}

ccc_mprun -n ${BRIDGE_MSUB_NPROC} ./isotropic -m WALLTIME LEVEL \
    2> log-${BRIDGE_MSUB_NPROC} > out-${BRIDGE_MSUB_NPROC}
~~~

~~~gnuplot Weak scaling on Irene (16^3^ per core)
set xlabel '# of cores'
set ylabel 'Speed (points x timestep/sec/core)'
set logscale x
plot '-' w lp t ''
8     366497
64    273009
512   196904
4096  167429
e
~~~

~~~gnuplot Weak scaling on Irene (32^3^ per core)
set xlabel '# of cores'
set ylabel 'Speed (points x timestep/sec/core)'
set logscale x
plot '-' w lp t ''
8     281250
64    189218
512   165312
4096  156683
e
~~~

## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/forcedturbulence.html)
*/
