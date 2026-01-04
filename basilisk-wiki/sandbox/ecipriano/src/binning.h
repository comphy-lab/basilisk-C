/**
# Binning/Cell-Agglomeration Technique

We wish to accelerate the direct integration of the reactive step by grouping
together cells with similar thermodynamic states in *bins*.  Finally, the bins
are integrated as if they were a single computational cell, and the results are
then re-mapped into the cell values in a conservative manner.  We follow the
algorithm proposed by [Goldin et al. 2009](#goldin2009cell).

## Usage

This algorithm can be written in a few lines of code using the functions and the
smart iterators defined in this module. The logic is:

* *binning*: divide the domain in a number of bins are return the table
* *integration*: integrate the ODE system each bin
* *remap*: map the bin values back to the cells

This procedure corresponds to the following code:

~~~literatec
BinTable * table = binning (fields, targets, eps, rho);

foreach_bin (table) {
  // ...
}
binning_remap (table, fields);
binning_cleanup (table);
~~~

### Missing features

1. ~~A different tolerance for each species~~
1. ~~Avoid empty bins to save memory~~
1. Energy-conserving remap
1. ~~Parallel simulations (each binning is local to the pid)~~

*/

/**
## Bin

A bin is an aggregate of cells. In practice, it is a structure containing a
global ID, the number of fields, the number of cells within the bin, and the
arrays `phi` and `phi0` which are the mass-averaged scalar field value within
the bin (with dimension equal to `nfields`).  The cells are stored using the
coordinate of the centroids, stored in the `cells` list of coordinates. */

typedef struct {
  size_t id;
  size_t nfields;
  size_t ncells;
  double * phi, * phi0;
  coord * cells;
  double rho, rho0;
  double cp, cp0;
} Bin;

Bin * new_bin (size_t id) {
  Bin * bin = malloc (sizeof (Bin));
  bin->id = id;
  bin->nfields = 0;
  bin->ncells = 0;
  bin->cells = NULL;
  bin->rho = 0.;
  bin->rho0 = 0.;
  bin->cp = 0.;
  bin->cp0 = 0.;
  return bin;
}

void free_bin (Bin * bin) {
  free (bin->phi);
  free (bin->phi0);
  free (bin->cells);
  free (bin);
}

void bin_append_cell (Point point, Bin * bin) {
  bin->ncells++;
  bin->cells = realloc (bin->cells, bin->ncells*sizeof (coord));

  coord o = {x,y,z};
  foreach_dimension()
    bin->cells[bin->ncells-1].x = o.x;
}

/**
## Bin Table

A bin table collects all the bins in the domain, and it facilitates loops and
operations over the bins. */

typedef struct {
  size_t ntargets;
  size_t nbins;
  Bin ** bins;
} BinTable;

BinTable * new_bintable_empty (size_t ntargets) {
  BinTable * table = malloc (sizeof (BinTable));
  table->ntargets = ntargets;
  table->nbins = 0;
  table->bins = NULL;
  return table;
}

BinTable * new_bintable (size_t ntargets, size_t nbins) {
  BinTable * table = new_bintable_empty (ntargets);
  table->nbins = nbins;
  table->bins = malloc (nbins*sizeof (Bin *));
  for (size_t i = 0; i < nbins; i++)
    table->bins[i] = new_bin (i);
  return table;
}

void free_bintable (BinTable * table) {
  for (size_t i = 0; i < table->nbins; i++)
    free_bin (table->bins[i]);
  free (table->bins);
  free (table);
}

size_t bintable_index (BinTable * table, size_t bin_id, bool * exist) {
  for (size_t i = 0; i < table->nbins; i++) {
    Bin * bin = table->bins[i];
    if (bin->id == bin_id) {
      *exist = true;
      return i;
    }
  }
  return table->nbins-1;
}

Bin * bintable_append (BinTable * table, size_t global_id) {
  bool exist = false;
  size_t local_id = bintable_index (table, global_id, &exist);
  if (exist) {
    Bin * bin = table->bins[local_id];
    bin->id = global_id;
    return bin;
  }
  else {
    table->nbins++;
    table->bins = realloc (table->bins, table->nbins*sizeof (Bin *));
    table->bins[table->nbins-1] = new_bin (table->nbins-1);
    Bin * bin = table->bins[table->nbins-1];
    bin->id = global_id;
    return bin;
  }
}

/**
## Iterators

We define useful iterators to simplify loops over the bins and their cells.  In
particular `foreach_bin` loops over each bin in the bin table;
`foreach_bin_field` iterates over the fields of a bin; `foreach_bin_cell`
iterates over the cells of a bin. */

macro foreach_bin (BinTable * table) {
  for (size_t i = 0; i < table->nbins; i++) {
    Bin * bin = table->bins[i]; NOT_UNUSED (bin);
    //if (bin->ncells > 0) // active bin
      {...}
  }
}

macro foreach_bin_field (Bin * bin) {
  for (size_t j = 0; j < bin->nfields; j++) {
    double phi = bin->phi[j], phi0 = bin->phi0[j];
    {...}
    bin->phi[j] = phi, bin->phi0[j] = phi0;
  }
}

macro foreach_bin_cell (Bin * bin) {
  for (size_t _k = 0; _k < bin->ncells; _k++) {
    coord o = bin->cells[_k];
    foreach_point (o.x,o.y,o.z, serial) {
      {...}
    }
  }
}

/**
## Program interface

The following functions can be called by the user to perform the binning,
remap the bins, and clean the memory. */

/**
Normalize a list of targets in order to obtain fields with fractional values
bounded between 0 and 1. */

void binning_normalize (scalar * targets) {
  int len = list_len (targets);
  double * tmax = malloc (len*sizeof (double));
  double * tmin = malloc (len*sizeof (double));

  for (int i = 0; i < len; i++) {
    scalar t = targets[i];
    tmax[i] = statsf (t).max;
    tmin[i] = statsf (t).min;
  }

  foreach() {
    for (int i = 0; i < len; i++) {
      scalar t = targets[i];
      t[] = (t[] - tmin[i])/(tmax[i] - tmin[i]);
    }
  }

  free (tmax);
  free (tmin);
}

/**
Create the bin table based on the list of (normalized) targets and the
user-defined tolerance.  The table creates the maximum potential bins, and
assigns the global index to each bin. Cells with null mask values are not added
to the bin, in order to exclude them from the averaging procedure. */

BinTable * binning_build_table (scalar * targets, double * eps,
    (const) scalar mask = unity)
{
  size_t len = (size_t)list_len (targets);

#if 0
  // All the possible combinations
  size_t nbins = 1;
  for (size_t i = 0; i < len; i++)
    nbins *= (L + 1);
#endif

  BinTable * table = new_bintable_empty (len);

  foreach() {
    size_t bin_j = 0;
    for (size_t i = 0; i < len; i++) {
      double L = 1./eps[i];
      scalar t = targets[i];
      size_t bin_i = (size_t)(t[]*L);
      if (bin_i > L) bin_i = L;
      bin_j += bin_i * pow (L + 1, i);
    }

    if (mask[]) {
      Bin * bin = bintable_append (table, bin_j);
      bin_append_cell (point, bin);
    }
  }
  return table;
}

/**
Populate the scalar fields of the bin based on the number of fields provided.
This is different from the list of targets, which does not necessarily have the
same dimensions of the list of fields. Fields are used for mass-averaging and
successive integration, while targets are used only for the binning. */

void binning_populate (BinTable * table, scalar * fields) {
  size_t nfields = list_len (fields);
  foreach_bin (table) {
    bin->nfields = nfields;
    bin->phi = malloc (nfields*sizeof (double));
    bin->phi0 = malloc (nfields*sizeof (double));
  }
}

/**
Calculate the value of the fields within the bins using a mass-averaged
approach. */

void binning_average (BinTable * table, scalar * fields,
    (const) scalar rho = unity)
{
  foreach_bin (table) {
    foreach_bin_field (bin) {
      scalar t = fields[j];
      double num = 0., den = 0.;
      foreach_bin_cell (bin) {
        num += t[]*rho[]*dv();
        den += rho[]*dv();
      }
      phi = (den > 0.) ? num / den : 0;
      phi0 = phi;
    }
  }
}

/**
Calculate the average of a given field in a bin. */

double bin_average (const Bin * bin, scalar field) {
  double num = 0.;
  foreach_bin_cell (bin)
    num += field[];
  return bin->ncells ? num / bin->ncells : 0;
}

double bin_volaverage (const Bin * bin, scalar field) {
  double num = 0., den = 0.;
  foreach_bin_cell (bin) {
    num += field[]*dv();
    den += dv();
  }
  return (den > 0.) ? num / den : 0.;
}

/**
This function automatically normalizes the list of fields, splits the domain in
bins with the correct mass-averaged field value, and returns the bin table. If
a density field `rho` is not provided, the averaging step is volume-averaged.
The mask field is an Heaviside function which excludes cells with null values.
This is useful in multiphase codes, when we want to integrate the systems just
in specific cells. */

BinTable * binning (scalar * fields, scalar * targets, double * eps,
    (const) scalar rho = unity, (const) scalar cp = unity,
    (const) scalar mask = unity)
{
  scalar * tnorm = list_clone (targets);
  foreach()
    for (size_t i = 0; i < list_len (targets); i++) {
      scalar tn = tnorm[i], t = targets[i];
      tn[] = t[];
    }
  binning_normalize (tnorm);
  BinTable * table = binning_build_table (tnorm, eps, mask);
  binning_populate (table, fields);
  binning_average (table, fields, rho);
  foreach_bin (table) {
    bin->rho = bin_volaverage (bin, rho);
    bin->cp = bin_average (bin, cp);
    bin->rho0 = bin->rho;
    bin->cp0 = bin->cp;
  }
  delete (tnorm), free (tnorm);
  return table;
}

/**
The bin field values are mapped-back to the original list of fields depending on
the variation of the bin field value. */

void binning_remap (BinTable * table, scalar * fields,
    (const) scalar rho = unity, (const) scalar cp = unity)
{
  foreach_bin (table) {
    foreach_bin_cell (bin) {
      foreach_bin_field (bin) {
        scalar t = fields[j];
        t[] += (phi - phi0);
      }
      if (!is_constant (rho)) rho[] += (bin->rho - bin->rho0);
      if (!is_constant (cp)) cp[] += (bin->cp - bin->cp0);
    }
  }
}

/**
Clean bin and bin table from memory. */

void binning_cleanup (BinTable * table) {
  free_bintable (table);
}

/**
### Post-Processing

The following functions allow the user to monitor the statistics of the binning
process. */

typedef struct {
  // total number of (active/non-empty) bins
  size_t nactive;
  // average number of cells per bin
  size_t navg;
  // maximum number of cells in a bin
  size_t nmax;
  // minimum number of cells in a bin
  size_t nmin;
  // number of masked cells
  size_t nmask;
} bstats;

bstats binning_stats (const BinTable * table) {
  bstats s = {0};
  s.nmax = 0, s.nmin = 1e9;

  size_t nreal = 0;
  foreach_bin (table) {
    foreach_bin_cell (bin) {
      nreal++;
    }
  }
  s.nmask = grid->n - nreal;

  foreach_bin (table) {
    s.nactive++;
    s.nmax = max (s.nmax, bin->ncells);
    s.nmin = min (s.nmin, bin->ncells);
  }
  s.navg = (s.nactive > 0) ? nreal / s.nactive : 0;
  return s;
}

void binning_ids (const BinTable * table, scalar ids) {
  foreach()
    ids[] = nodata;
  foreach_bin (table) {
    foreach_bin_cell (bin) {
      ids[] = bin->id;
    }
  }
}

/**
## References

~~~bib
@article{goldin2009cell,
  title={A cell agglomeration algorithm for accelerating detailed chemistry in CFD},
  author={Goldin, Graham M and Ren, Zhuyin and Zahirovic, Selma},
  journal={Combustion Theory and Modelling},
  volume={13},
  number={4},
  pages={721--739},
  year={2009},
  publisher={Taylor \& Francis}
}
~~~
*/
