/**
# Shared helpers for spectra sampling

`snap_to_cell()`, used by both [spectra_sample.h](spectra_sample.h)
(`locate()`-based) and
[spectra_sample_restricted.h](spectra_sample_restricted.h)
(`foreach_level()`-based). Both headers include this one instead of defining
it themselves, so a comparison test can include both without a
duplicate-symbol error.
*/

#ifndef SNAP_TO_CELL_DEFINED
#define SNAP_TO_CELL_DEFINED

/**
Snap a target height to the nearest cell centre. `foreach_region` returns the
cell containing the point, so a height on a face resolves to one side
arbitrarily -- $z = 0$ is a face whenever the grid size is even. */

double snap_to_cell (double h, int m)
{
  double del = L0/m;
  return Z0 + (floor ((h - Z0)/del) + 0.5)*del;
}

#endif
