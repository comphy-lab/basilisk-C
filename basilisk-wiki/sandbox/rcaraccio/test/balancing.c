/**
# A Preliminar test for a more general load balancing algorithm
The current load balancing between processors using MPI divides
the cells equally between each processor by cutting the space-filling
curve into equal pieces. This can cause issues when some cells require
significantly more computational power than others, since equal cell
counts no longer mean equal load. We are here to tackle this issue.

This first test case tries the new implementation on a static, uniformly
refined grid. We create a ficticius weight, representative of the cell computational
cost. The idea is to evenly split the cumulative weight function (CWF)
into equally loaded intervals. The CWF is computed according to the Z-order,
with the idea to still keep cells that are spatially close also on the same
processor.
*/


#include "grid/quadtree.h"
#include "test/refine_unbalanced.h"
#include "view.h"

const int maxlevel = 8;

/**
## New partition policy

Given a leaf's Z-order index, the global leaf count, the
number of processes and the cumulative-weight array, return the
process id
*/

static int new_balanced_pid (long index, long nt, int nproc, double* cf) {
  double il = cf[nt - 1]/nproc; // ideal load per process
  int idx = (int) floor ((cf[index] - cf[0])/il);
  return min (idx, nproc - 1);
}

/**
## Cumulative weight

`cf[i]` is the cumulative sum of `w` over leaves, in the same
order they are visited by the partitioning loop. */

void compute_cf (double* cf, long nt, scalar w) {
  long idx = 0;
  double acc = 0.;
  foreach_leaf() {
    acc += w[];
    cf[idx++] = acc;
  }
}

/**
## Per-process cell count and load
Performance monitor function that returns the number of cells
and the load in each processor*/

void balance_score (long* counter, double* load, scalar w) {
  for (int ne = 0; ne < npe(); ne++) {
    counter[ne] = 0;
    load[ne] = 0.;
  }

  foreach() {
    counter[pid()] += 1;
    load[pid()] += w[];
  }

  if (pid() == 0) {
    MPI_Reduce(MPI_IN_PLACE, counter, npe(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, load, npe(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(counter, NULL, npe(), MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(load, NULL, npe(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
}

/**
## Static partitioning using the new policy

Copy of `mpi_partitioning()` with `balanced_pid()` replaced by
`new_balanced_pid()`. */

trace
void new_mpi_partitioning (scalar w)
{
  prof_start ("new_mpi_partitioning");

  long nt = 0;
  foreach (serial) nt++;

  // We need to compute cf once before the loop
  double cf[nt];
  compute_cf (cf, nt, w);

  long i = 0;
  tree->dirty = true;
  foreach_cell_post (is_active (cell))
    if (is_active (cell)) {
      if (is_leaf (cell)) {
        cell.pid = new_balanced_pid (i++, nt, npe(), cf);
        if (cell.neighbors > 0) {
          int pid = cell.pid;
          foreach_child() cell.pid = pid;
        }
        if (!is_local (cell)) cell.flags &= ~active;
      }
      else {
        cell.pid = child(0).pid;
        bool inactive = true;
        foreach_child() if (is_active (cell)) { inactive = false; break; }
        if (inactive) cell.flags &= ~active;
      }
    }
  flag_border_cells();
  prof_stop();
  mpi_boundary_update_buffers();
}

/**
## Mesh setup

Single root cell, every cell active and locally owned, then refined
without MPI balancing, required by `*_mpi_partitioning()`.
See [this test](/src/test/balance5.c)*/

void prepare_grid (int maxlevel) {
  init_grid (1);
  foreach_cell() {
    cell.pid = pid();
    cell.flags |= active;
  }
  tree->dirty = true;
  refine_unbalanced (level < maxlevel, NULL);
}

/**
## Test weight field
A ficticius measure of the load in each processor*/

void set_weights (scalar w) {
  foreach()
    w[] = sqrt(sq(x) + sq(y));
}

int main() {

  /**
  ### New partitioning method */

  prepare_grid (maxlevel);
  scalar w[], mpi_idx[];
  set_weights (w);

  long nt = 0;
  foreach (serial) nt++;
  double cf[nt];
  compute_cf (cf, nt, w);

  FILE * fcf = fopen ("cum_weight", "w");
  fprintf (fcf, "#nt: %ld\n", nt);
  for (long i = 0; i < nt; i++)
    fprintf (fcf, "%ld %g\n", i, cf[i]);
  fclose (fcf);

  double ne = cf[nt - 1]/npe();
  FILE * fb = fopen ("cum_bounds", "w");
  for (int p = 1; p < npe(); p++) {
    double target = p*ne;
    long bidx = 0;
    while (bidx < nt && cf[bidx] < target) bidx++;
    fprintf (fb, "%ld %g\n", bidx, target);
  }
  fclose (fb);

  new_mpi_partitioning (w);

  long counter_new[npe()]; double load_new[npe()];
  balance_score (counter_new, load_new, w);

  foreach()
    mpi_idx[] = pid();

  view (tx = -0.5, ty = -0.5);
  box();
  cells();
  squares ("mpi_idx", min = 0, max = npe() - 1);
  save ("partition.png");
  squares ("w", spread = -1);
  save ("load.png");

  /**
  ### Baseline for comparison */

  prepare_grid (maxlevel);
  set_weights (w);

  mpi_partitioning();

  long counter_old[npe()]; double load_old[npe()];
  balance_score (counter_old, load_old, w);

  if (pid() == 0) {
    fprintf (stderr, "#pid(1) n_old(2) load_old(3) n_new(4) load_new(5)\n");
    for (int ne = 0; ne < npe(); ne++)
      fprintf (stderr, "%d %ld %g %ld %g\n", ne,
               counter_old[ne], load_old[ne],
               counter_new[ne], load_new[ne]);
  }
}

/**

## Results

![Weight field](balancing/load.png)
![Partitioning](aslamsharp/partitpartition)

~~~gnuplot Partitioning comparison: old vs new
set terminal png size 1000,500
set output "balance.png"
set multiplot layout 1,2 title "Partitioning comparison: old vs new"

set xlabel "PID"
set xtics 1

set ylabel "n cells"
set ytics nomirror
set yrange [0:*]

set y2label "load"
set y2tics nomirror
set y2range [0:*]

set grid xtics ytics
set key top right

set title "Old partitioning"
plot "log" using ($1-0.2):2 axes x1y1 with boxes lc rgb "yellow" title "n cells", \
     "log" using ($1+0.2):3 axes x1y2 with boxes lc rgb "blue" title "load"
set title "New partitioning"
plot "log" using ($1-0.2):4 axes x1y1 with boxes lc rgb "yellow" title "n cells", \
     "log" using ($1+0.2):5 axes x1y2 with boxes lc rgb "blue" title "load"
unset multiplot
unset output
~~~

~~~gnuplot Cumulative weight and equal-load partition boundaries
set terminal png size 800,600
set output "cum_weight.png"
set xlabel "cell index (Z-order)"
set ylabel "cumulative weight cf"
set grid
set key top left

stats "cum_bounds" using 1 nooutput
nb = STATS_records

plot "cum_weight" using 1:2 w l lw 2 lc rgb "blue" title "cf", \
     for [i=1:nb] "cum_bounds" every ::i-1::i-1 using (0):2:1:(0) \
     w vectors nohead lw 1 lc "red" notitle, \
     for [i=1:nb] "cum_bounds" every ::i-1::i-1 using 1:(0):(0):2 \
     w vectors nohead lw 1 lc "red" notitle
     
~~~
*/