/**
# Binning/Cell-Agglomeration Technique Test

We test the correct implementation of the binning algorithm. */

#include "common.h"
#include "run.h"
#include "binning.h"
#include "view.h"

int main (void) {
  init_grid (1 << 7);
  run();
}

#define radial(r,r1,r2,T1,T2) (-r1*r2/(r*(r1-r2))*(T1 - T2) + (r1*T1 - r2*T2)/(r1-r2))
#define gaussian(x, y, sigma) ( exp(-((sq(x - 0.) + sq(y - 0.)) / (2 * sq(sigma)))) )

/**
We assume two trends for the temperature T and mass fraction Y fields. The
binning tolerance is `eps`. The `binid` field is used to store bin id's for
visualization. TI and YI sore the initial values of the fields. */

scalar T[], Y[], TI[], YI[], * targets = {T,Y};
double eps = 1e-1;
scalar binid[];

event init (i = 0) {
  scalar mask[], * fields = targets;
  foreach() {
    double r = sqrt (sq(x) + sq(y));
    T[] = radial (r, 0.2, 0.8, 300., 800.);
    T[] = (r <= 0.2) ? 300. : (r >= 0.8) ? 800. : T[];
    Y[] = gaussian (x, y, 0.2);
    mask[] = (r <= 0.2) ? 0. : 1.;
    T[] *= mask[];
    Y[] *= mask[];

    TI[] = T[];
    YI[] = Y[];
  }

  /**
  We create the binning table, diving the domain in a number of bins. */

  BinTable * table = binning (fields, targets, (double[]){eps,eps}, mask = mask);

  /**
  We fill a scalar fields with the bin indeces, in order to visualize the
  domain partition. */

  binning_ids (table, binid);

  /**
  We compute and print binning statistics. */

  bstats bs = binning_stats (table);

  fprintf (stderr, "Total number of active bins         = %zu\n", bs.nactive);
  fprintf (stderr, "Average number of cells in the bins = %zu\n", bs.navg);
  fprintf (stderr, "Maximum number of cells in the bins = %zu\n", bs.nmax);
  fprintf (stderr, "Minimum number of cells in the bins = %zu\n", bs.nmin);
  fprintf (stderr, "Number of masked/excluded cells     = %zu\n", bs.nmask);
  fprintf (stderr, "\n");

  foreach_bin (table)
    fprintf (stderr, "BIN[%3zu] has %zu cells\n", bin->id, bin->ncells);
  fprintf (stderr, "\n");

  /**
  We print the values of the fields phi within every bin to ensure that empty
  values belonging to the mask are not considered bins. */

  foreach_bin (table) {
    foreach_bin_field (bin)
      fprintf (stderr, "BIN[%zu] has phi[%zu] = %g\n", i, j, phi);
    fprintf (stderr, "\n");
  }

  /**
  We perform a (fake) bin integration, which should change the field values. */

  foreach_bin (table) {
    foreach_bin_field (bin) {
      phi0 = phi;
      phi += 10.*dt;
    }
  }

  /**
  We map the bin values back to the cells. In the conditions of this test, the
  initial and final field values must be the same. */

  binning_remap (table, fields);

  /**
  We check that the fields changed as expected. */

  fprintf (stderr, "Total change of Y = %g\n", change (YI, Y));
  fprintf (stderr, "Total change of T = %g\n", change (TI, T));

  /**
  We clean the binning structures from the memory. */

  binning_cleanup (table);

  clear();
  view (tx = -0.5, ty = -0.5);
  squares ("T", spread = -1);
  save ("temperature.png");

  clear();
  view (tx = -0.5, ty = -0.5);
  squares ("Y", spread = -1);
  save ("massfrac.png");

  clear();
  view (tx = -0.5, ty = -0.5);
  squares ("mask", spread = -1);
  save ("mask.png");

  clear();
  view (tx = -0.5, ty = -0.5);
  squares ("binid", spread = -1);
  save ("ids.png");
}

/**
## Results

![Initial temperature field](binning/temperature.png)

![Initial mass fraction field](binning/massfrac.png)

![Mask field](binning/mask.png)

![Distribution of the bins](binning/ids.png)

*/

