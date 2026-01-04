/**
# Parallel scalability

This code can be used to test the scalability of common operations
(traversal, restriction, boundary conditions etc...) and their
combinations, in particular the multigrid Poisson solver. */

#include "poisson.h"
#include "utils.h"

scalar a[], b[];

static void mpi_print (timer t, int i, long tnc,
		       const char * name)
{
  double mpi[npe()];
  timing s = timer_timing (t, i, tnc, mpi);

#if 0
  scalar wt[];
  foreach()
    wt[] = mpi[pid()];
  char fname[80];
  sprintf (fname, "%s-%d.ppm", name, npe());
  FILE * fp = pid() ? NULL : fopen (fname, "w");
  output_ppm (wt, fp, n = 512);
#endif
  
  if (pid() == 0) {

    /**
    *s.min/i*, *s.avg/i*, *s.max/i* are the minimum, average and maximum
    *times spent in communication routines. */

    printf ("%d %g %g %g %s %g %g %g %.2g%% %ld %ld",
	    npe(), s.cpu/i, s.real/i, s.speed, name, s.min/i, s.avg/i, s.max/i,
	    100.*s.avg/s.real, s.mem, s.tnc/i);

    /**
    We also output the times spent in communication routines for each process. */
#if 0
    for (int j = 0; j < npe(); j++)
      printf (" %g", mpi[j]/i);
#endif

    /**
    If [memory profiling](/src/README.mtrace) is enabled we output the
    maximum allocated memory (in kilobytes). */

@if MTRACE
    printf (" %ld", pmtrace.max/1024);
@endif
    
    printf ("\n");
  }
  
  MPI_Barrier (MPI_COMM_WORLD);
}

int main (int argc, char * argv[])
{
  int maxlevel = argc > 1 ? atoi(argv[1]) : (dimension == 2 ? 8 : 5);
  timer t;

  double speed = 1e6; // approx speed in points/sec/core
  long tnc = ((long)1) << (dimension*maxlevel);
  int i, nloops = max(0.1*speed*npe()/(double)tnc, 1);
  if (!pid())
    fprintf (stderr, "grid: %s\nnloops = %d\n", GRIDNAME, nloops);

  if (tnc*nloops/(speed*npe()) > 100.) {
    fprintf (stderr, "this run would probably take more than 100 seconds."
	     " Aborting ...\n");
    exit(1);
  }

  size (1[0]);
  init_grid (N);
  foreach()
    a[] = 0., b[] = 0.;
  poisson (a, b); // to force allocation of extra fields
  
  MPI_Barrier (MPI_COMM_WORLD);
  t = timer_start();
  init_grid (1 << maxlevel);
  tnc = 0;
  foreach(reduction(+:tnc))
    tnc++;
  mpi_print (t, 1, tnc, "refine");

#if TREE  
  assert (tree_is_full());
#endif
  
  /**
  We fill `a` with a simple function. */

  MPI_Barrier (MPI_COMM_WORLD);
  i = 0;
  t = timer_start();
  while (i < nloops) {
    foreach()
      a[] = cos(pi*x)*cos(pi*y)*cos(pi*z);
#if 0
    boundary ({a});
#else
    prof_start ("barrier");
    MPI_Barrier (MPI_COMM_WORLD);
    prof_stop();
#endif
    i++;
  }
  mpi_print (t, i, tnc*i, "cos");
  boundary ({a});

  /**
  Here we compute
  $$
  b = \nabla^2 a
  $$
  using a 5-points Laplacian operator. */

  MPI_Barrier (MPI_COMM_WORLD);
  i = 0;
  t = timer_start();
  while (i < nloops) {
    foreach() {
      b[] = 0.;
      foreach_dimension()
        b[] += a[1] + a[-1];
      b[] = (b[] - 2.*dimension*a[])/sq(Delta);
    }
    boundary ({b});
    i++;
  }
  mpi_print (t, i, tnc*i, "laplacian");
  
  /**
  Something simpler: the sum of `a` over the entire mesh. */

  MPI_Barrier (MPI_COMM_WORLD);
  i = 0;
  t = timer_start();
  double sum = 0.;
  while (i < nloops) {
    sum = 0.;
    foreach(reduction(+:sum))
      sum += b[];
    i++;
  }
  mpi_print (t, i, tnc*i, "sum");
  if (pid() == 0)
    fprintf (stderr, "sum: %g\n", sum);

  /**
  The restriction operator. */

  MPI_Barrier (MPI_COMM_WORLD);
  i = 0;
  t = timer_start();
  while (i < nloops) {
    boundary ({b});
    restriction ({b});
    i++;
  }
  mpi_print (t, i, tnc*i, "restriction");

  /**
  And finally the Poisson solver. */

  MPI_Barrier (MPI_COMM_WORLD);
  i = 0;
  TOLERANCE = HUGE;
  t = timer_start();
  while (i < nloops) {
    poisson (a, b);
    i++;
  }
  mpi_print (t, i, tnc*i, "poisson");

  scalar e[];
  foreach()
    e[] = a[] - cos(pi*x)*cos(pi*y)*cos(pi*z);
  double max = normf(e).max;
  if (pid() == 0)
    fprintf (stderr, "error: %g\n", max);
  assert (normf(e).max < 1e-10);
  //  output_ppm (e, file = "error.png");

  sum = 0.;
  int n = 0;
  foreach (serial) {
    e[] = (long) &(a[]);
    sum += e[];
    n++;
  }
  foreach()
    e[] -= sum/n;
#if 0
  FILE * fp = pid() ? NULL : fopen ("map.ppm", "w");
  output_ppm (e, fp, n = 512);
#endif
  
  int nmin = n, nmax = n;
  mpi_all_reduce (nmin, MPI_INT, MPI_MIN);
  mpi_all_reduce (nmax, MPI_INT, MPI_MAX);
  if (pid() == 0)
    fprintf (stderr, "balance %d %d\n", nmin, nmax);
}

/**
## How to run on Occigen

This test is run on
[Occigen](https://www.cines.fr/calcul/materiels/occigen/configuration/). The
C99 source code is generated on a system with *qcc* installed and then
copied to Occigen using something like

~~~bash
qcc -grid=octree -source mpi-laplacian.c
scp _mpi-laplacian.c popinet@occigen.cines.fr:
~~~

On Occigen one can then compile the code using

~~~bash
module load bullxmpi
module load intel
mpicc -Wall -std=c99 -O2 -D_MPI=1 -D_GNU_SOURCE=1 _mpi-laplacian.c \
     -o mpi-laplacian -lm
~~~

The following script (*run.sh*) can then be used to run the code

~~~bash
#!/bin/bash
#SBATCH -J basilisk
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --threads-per-core=1
#SBATCH --time=00:10:00
#SBATCH --output basilisk.output
#SBATCH --exclusive

LEVEL=11

module purge
module load bullxmpi
module load intel

srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./mpi-laplacian $LEVEL \
     2> log-$LEVEL-$SLURM_NTASKS > out-$LEVEL-$SLURM_NTASKS
~~~

LEVEL is varied from 9 (~134 million grid points) to 11 (~8.6 billion
grid points) and the number of processes up to 16384 using something like

~~~bash
sbatch --ntasks=1 --ntasks-per-node=1 run.sh
sbatch --ntasks=2 --ntasks-per-node=2 run.sh
sbatch --ntasks=4 --ntasks-per-node=4 run.sh
sbatch --ntasks=8 --ntasks-per-node=8 run.sh
sbatch --ntasks=16 --ntasks-per-node=16 run.sh
sbatch --ntasks=32 --ntasks-per-node=16 --nodes=2 run.sh
sbatch --ntasks=64 --ntasks-per-node=16 --nodes=4 run.sh
sbatch --ntasks=128 --ntasks-per-node=16 --nodes=8 run.sh
sbatch --ntasks=256 --ntasks-per-node=16 --nodes=16 run.sh
sbatch --ntasks=512 --ntasks-per-node=16 --nodes=32 run.sh
sbatch --ntasks=1024 --ntasks-per-node=16 --nodes=64 run.sh
sbatch --ntasks=2048 --ntasks-per-node=16 --nodes=128 run.sh
sbatch --ntasks=4096 --ntasks-per-node=16 --nodes=256 run.sh
sbatch --ntasks=8192 --ntasks-per-node=16 --nodes=512 run.sh
sbatch --ntasks=16384 --ntasks-per-node=24 --nodes=683 run.sh
~~~

The results can then be collected using

~~~bash
tar czvf occigen.tgz out-*-*
~~~

## Results on Occigen

The memory usage per core is given below. The increase for a large number
of cores corresponds to the memory overhead of communication buffers
etc...

~~~gnuplot Memory usage on Occigen
# generate results for Occigen
cd 'occigen/3D'

# generate weak scaling curves
! bash weak.sh > weak

set logscale
set logscale x 2
set grid
set xrange [2:32768]
set format x '2^{%L}'
set xtics 2

set xlabel "# of cores"
set ylabel "Memory/core (GB)"
minlevel=9
maxlevel=11
plot [][0.1:] for [i=minlevel:maxlevel] \
     '< sh table.sh poisson '.i u 1:($2/$1) t ''.i.' levels' w lp, \
     18/x**0.9 t '18/x^{0.9}'
cd '../..'
~~~

The wall-clock time for one iteration of the multigrid Poisson solver
is given below.

The pink lines connect points corresponding with weak (or
*iso-granular*) scaling i.e. multiplying by eight both the computation
size and the number of cores. The ideal weak scaling would give
horizontal lines (i.e. constant computation time for
proportionally-increasing problem sizes and number of cores).

~~~gnuplot Wall-clock time on Occigen for the Poisson solver
cd 'occigen/3D'
set ylabel 'Time (sec)'
plot [][0.01:100] for [i=minlevel:maxlevel] \
     '< sh time.sh poisson '.i u 1:2 t ''.i.' levels' w lp, \
     'weak' u 1:2 w lp t 'weak scaling', \
     600/x**0.95 t '600/x^{0.95}'
cd '../..'
~~~

The time spent in communication routines is illustrated below.

~~~gnuplot Communication time on Occigen for the Poisson solver
cd 'occigen/3D'
plot [][1e-2:10] for [i=minlevel:maxlevel] \
     '< sh time.sh poisson '.i u 1:3 w lp t ''.i.' levels', \
     4.5/x**0.65 t '4.5/x^{0.65}'
cd '../..'
~~~

Similar results are obtained for a pure Laplacian.

~~~gnuplot Wall-clock time on Occigen for the Laplacian
cd 'occigen/3D'
plot [][1e-4:10] for [i=minlevel:maxlevel] \
     '< sh time.sh laplacian '.i u 1:2 w lp t ''.i.' levels', \
     50/x**0.93 t '50/x^{0.93}'
cd '../..'
~~~

~~~gnuplot Communication time on Occigen for the Laplacian
cd 'occigen/3D'
plot [][1e-4:1] for [i=minlevel:maxlevel] \
     '< sh time.sh laplacian '.i u 1:3 w lp t ''.i.' levels', \
     2./x**0.7 t '2/x^{0.7}'
cd '../..'
~~~

And for the restriction operator.

~~~gnuplot Wall-clock time on Occigen for the restriction operator
cd 'occigen/3D'
plot [][:1] for [i=minlevel:maxlevel] \
     '< sh time.sh restriction '.i u 1:2 w lp t ''.i.' levels', \
     18/x**0.85 t '18/x^{0.85}'
cd '../..'
~~~

~~~gnuplot Communication time on Occigen for the restriction operator
cd 'occigen/3D'
plot [][1e-3:1] for [i=minlevel:maxlevel] \
     '< sh time.sh restriction '.i u 1:3 w lp t ''.i.' levels', \
     2.8/x**0.66 t '2.8/x^{0.66}'
cd '../..'
~~~
*/
