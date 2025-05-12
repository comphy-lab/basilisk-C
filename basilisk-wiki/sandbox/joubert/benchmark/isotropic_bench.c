/**
# Forced isotropic turbulence in a triply-periodic box wtht view

This is a copy paste of the [isotropic code](http://basilisk.fr/src/examples/isotropic.c)
used to check the cpu scalability of Basilisk on supercomputers
using only 5 min of cp time. For weak scaling test in a cubic
domain number of processors should be a power of 8 ($8^1$,$8^2$,etc)
with multigrid in 3D.

We compute the evolution of forced isotropic turbulence (see [Rosales
& Meneveau, 2005](/src/references.bib#rosales2005)).
The initial condition is an unstable solution to the incompressible
Euler equations. */

#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"

/**
We monitor performance statistics and control the maximum runtime. */

#include "navier-stokes/perfs.h"

#define MU 0.01

/**
We need to store the variable forcing term. */

face vector av[];

/**
The code takes the level of refinement as optional command-line
argument (as well as an optional maximum runtime). */

int maxlevel = 4;

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi(argv[1]);

  /**
  The domain is $(2\pi)^3$ and triply-periodic. */
  
  L0 = 2.*pi;
  foreach_dimension()
    periodic (right);

  /**
  The viscosity is constant. The acceleration is defined below. The
  level of refinement is *maxlevel*. */

  const face vector muc[] = {MU,MU,MU};
  mu = muc;
  a = av;
  N = 1 << maxlevel;
  run();
}

/**
## Initial conditions

The initial condition is "ABC" flow. This is a laminar base flow that 
is easy to implement in both Basilisk and a spectral code. */

event init (i = 0) {
  if (!restore (file = "restart"))
    foreach() {
      u.x[] = cos(y) + sin(z);
      u.y[] = sin(x) + cos(z);
      u.z[] = cos(x) + sin(y);
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
  foreach_face()
    av.x[] += 0.1*((u.x[] + u.x[-1])/2. - ubar.x);
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
    fprintf (ferr, "t dissipation energy Reynolds\n");
  fprintf (ferr, "%g %g %g %g\n",
	   t, vd, ke, 2./3.*ke/MU*sqrt(15.*MU/vd));
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

Here is a list of the supercomputer with a description of the node
hardware and an example of compilation and submition of jobs I had
the oppurtunity to test:

## Running with MPI on Irene

[Irene](http://www-hpc.cea.fr/en/complexe/tgcc-JoliotCurie.htm) is the
supercomputer at CEA.

On the local machine

~~~bash
%local qcc -source -D_MPI=1 isotropic_bench.c
%local scp _isotropic_bench.c irene.ccc.cea.fr:
~~~

On irene

~~~bash
mpicc -Wall -O2 -std=c99 -xCORE-AVX512 _isotropic_bench.c -o isotropic_bench -lm
sed -e 's/WALLTIME/300/g' -e 's/LEVEL/7/g' run.sh | ccc_msub -n 512
~~~

with the `run.sh` script

~~~bash
#!/bin/bash
#MSUB -r isotropic_bench
#MSUB -T WALLTIME
#MSUB -o basilisk_%I.out
#MSUB -e basilisk_%I.log
#MSUB -q skylake
#MSUB -A gen7760
#MSUB -m scratch

set -x
cd ${BRIDGE_MSUB_PWD}

ccc_mprun -n ${BRIDGE_MSUB_NPROC} ./isotropic_bench -m WALLTIME LEVEL \
    2> log-${BRIDGE_MSUB_NPROC} > out-${BRIDGE_MSUB_NPROC}
~~~

## Running with MPI on Jean-Zay

[Jean-Zay](http://www.idris.fr/jean-zay/cpu/jean-zay-cpu-hw.html) is the
supercomputer at IDRIS.

On the local machine

~~~bash
%local qcc -source -D_MPI=1 isotropic_bench.c
%local scp _isotropic_bench.c irene.ccc.cea.fr:
~~~

On irene

~~~bash
module purge
module load intel-tools-19
mpicc -Wall -std=c99 -O2 _$NAME.c -o $NAME -lm
~~~

with the `run.sh` script

~~~bash
#!/bin/bash
#SBATCH --job-name=isotropic_bench      
#SBATCH --partition=cpu_gct3        
#SBATCH --ntasks=512                
#SBATCH --ntasks-per-node=40       
#SBATCH --hint=nomultithread       
#SBATCH --time=00:05:00            
#SBATCH --output=isotropic_bench%j.out  
#SBATCH --error=isotropic_bench%j.out    

module purge
module load intel-tools-19
 
LEVEL=7
NAME=isotropic_bench

srun --mpi=pmi2 ./$NAME $LEVEL \
     2> log-$LEVEL > out-$LEVEL
~~~

## Running with MPI on MeSu

[MeSu](https://hpcave.upmc.fr/index.php/resources/mesu-supercomputer/) is the
supercomputer at Sorbonne Université.

On the local machine

~~~bash
%local qcc -source -D_MPI=1 isotropic_bench.c
%local scp _isotropic_bench.c mesu.dsi.upmc.fr:
~~~

On mesu (using 512 cores)

~~~bash
module load mpt
mpicc -Wall -O2 -std=c99 _isotropic_bench.c -o isotropic_bench -lm
sed 's/WALLTIME/05:00/g' run.sh | qsub
~~~

with the `run.sh` script

~~~bash
#!/bin/bash 
#PBS -l select=22:ncpus=24:mpiprocs=24
#PBS -l walltime=WALLTIME
#PBS -N isotropic_bench
#PBS -j oe  
# load modules 
module load mpt
# change to the directory where program job_script_file is located 
cd $PBS_O_WORKDIR 
# mpirun -np 64 !!!! does not work !!!!
NP=512
mpiexec_mpt -n $NP ./isotropic_bench -m WALLTIME 2>> log.$NP >> out.$NP
~~~

## Running with MPI on Occigen

On the local machine

~~~bash
%local qcc -source -D_MPI=1 isotropic_bench.c
%local scp _isotropic_bench.c occigen.cines.fr:
~~~

On occigen (using 512 cores)

~~~bash
module purge
module load openmpi
module load intel
mpicc -Wall -O2 -std=c99 _isotropic_bench.c -o isotropic_bench -lm
sed 's/WALLTIME/05:00/g' run.sh | sbatch
~~~

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
     ./isotropic_bench -m WALLTIME $LEVEL \
     2> log-$LEVEL-$SLURM_NTASKS > out-$LEVEL-$SLURM_NTASKS
~~~

## Running with MPI on Piz-Daint

[Piz-daint](https://www.cscs.ch/computers/piz-daint/) is the
supercomputer at CSCS.

On the local machine

~~~bash
%local qcc -source -D_MPI=1 isotropic_bench.c
%local scp _isotropic_bench.c occigen.cines.fr:
~~~

On occigen (using 512 cores)

~~~bash
module load daint-mc
NAME=isotropic_bench                                                                                       
cc -Wall -std=c99 -O2 _$NAME.c -o $NAME -lm
~~~

With the `run.sh` script

~~~bash
#!/bin/bash -l
#SBATCH --job-name="isotropic"
#SBATCH --account="s1136"
#SBATCH --mail-type=ALL
#SBATCH --time=00:05:00
#SBATCH --nodes=15
#SBATCH --ntasks=512
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=mc
#SBATCH --hint=nomultithread

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
LEVEL=7
NAME=isotropic_bench                                                                                       

srun -n $SLURM_NTASKS ./$NAME $LEVEL \
2> log-$LEVEL-$SLURM_NTASKS > out-$LEVEL-$SLURM_NTASKS
~~~

## Scalability on supercomputers

Here we do a weak scaling (constant number of cells/core) test on
the previous list of supercomputer. We also test the different
partition available ( different hardware if there is several it is mentionned 
after the supercomputer name in the results) available for each supercomputer:

~~~gnuplot Weak scaling (16^3 per core)
set xlabel '# of cores'
set ylabel 'Speed (points x timestep/sec/core)'
set logscale x 2
plot '../weak_scaling/log_isotropic_mesu' i 0:0 u 12:($11/$12) w lp t 'MeSu beta',\
'../weak_scaling/log_isotropic_ireneSKL' i 0:0 u 1:2 w lp t 'Irene SKL',\
'../weak_scaling/log_isotropic_ireneROME' i 0:0 u 12:($11/$12) w lp t 'Irene ROME',\
'../weak_scaling/log_isotropic_occigenBDW' i 0:0 u 12:($11/$12) w lp t 'occigen BDW',\
'../weak_scaling/log_isotropic_occigenHSW' i 0:0 u 12:($11/$12) w lp t 'occigen HSW',\
'../weak_scaling/log_isotropic_jean-zay' i 0:0 u 12:($11/$12) w lp t 'Jean-Zay',\
'../weak_scaling/log_isotropic_pizdaint' i 0:0 u 12:($11/$12) w lp t 'Pizdaint'
~~~

~~~gnuplot Weak scaling (32^3 per core)
set xlabel '# of cores'
set ylabel 'Speed (points x timestep/sec/core)'
set logscale x 2
plot 'log_isotropic_mesu' i 1:1 u 12:($11/$12) w lp t 'MeSu beta',\
'../weak_scaling/log_isotropic_ireneSKL' i 1:1 u 1:2 w lp t 'Irene SKL',\
'../weak_scaling/log_isotropic_ireneROME' i 1:1 u 12:($11/$12) w lp t 'Irene ROME',\
'../weak_scaling/log_isotropic_occigenBDW' i 1:1 u 12:($11/$12) w lp t 'occigen BDW',\
'../weak_scaling/log_isotropic_occigenHSW' i 1:1 u 12:($11/$12) w lp t 'occigen HSW',\
'../weak_scaling/log_isotropic_jean-zay' i 1:1 u 12:($11/$12) w lp t 'Jean-Zay',\
'../weak_scaling/log_isotropic_pizdaint' i 1:1 u 12:($11/$12) w lp t 'Pizdaint'
~~~

~~~gnuplot Weak scaling (64^3 per core)
set xlabel '# of cores'
set ylabel 'Speed (points x timestep/sec/core)'
set logscale x 2
plot '../weak_scaling/log_isotropic_mesu' i 2:2 u 12:($11/$12) w lp t 'MeSu beta',\
'../weak_scaling/log_isotropic_ireneSKL' i 2:2 u 12:($11/$12) w lp t 'Irene SKL',\
'../weak_scaling/log_isotropic_ireneROME' i 2:2 u 12:($11/$12) w lp t 'Irene ROME',\
'../weak_scaling/log_isotropic_occigenBDW' i 2:2 u 12:($11/$12) w lp t 'occigen BDW',\
'../weak_scaling/log_isotropic_occigenHSW' i 2:2 u 12:($11/$12) w lp t 'occigen HSW',\
'../weak_scaling/log_isotropic_jean-zay' i 2:2 u 12:($11/$12) w lp t 'Jean-Zay',\
'../weak_scaling/log_isotropic_pizdaint' i 2:2 u 12:($11/$12) w lp t 'Pizdaint'
~~~

~~~gnuplot Weak scaling (128^3 per core)
set xlabel '# of cores'
set ylabel 'Speed (points x timestep/sec/core)'
set logscale x 2
plot '../weak_scaling/log_isotropic_mesu' i 3:3 u 12:($11/$12) w lp t 'MeSu beta',\
'../weak_scaling/log_isotropic_ireneSKL' i 3:3 u 12:($11/$12) w lp t 'Irene SKL',\
'../weak_scaling/log_isotropic_ireneROME' i 3:3 u 12:($11/$12) w lp t 'Irene ROME',\
'../weak_scaling/log_isotropic_occigenBDW' i 3:3 u 12:($11/$12) w lp t 'occigen BDW',\
'../weak_scaling/log_isotropic_occigenHSW' i 3:3 u 12:($11/$12) w lp t 'occigen HSW',\
'../weak_scaling/log_isotropic_jean-zay' i 3:3 u 12:($11/$12) w lp t 'Jean-Zay',\
'../weak_scaling/log_isotropic_pizdaint' i 3:3 u 12:($11/$12) w lp t 'Pizdaint'
~~~
*/
