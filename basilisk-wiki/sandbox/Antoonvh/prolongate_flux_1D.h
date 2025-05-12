/**
# Flux prolongation (in 1D) 

Cell-face fluxes may be approximated from a coarser-level
solution. This method may be attractive for the computation of fluxes
at resolution boundaries. It is the 'upside-down' version of
`halo_flux()`, reducing the importance of the accuracy of
halo-cell-centered values. For now it is only meant to work in 1D.

The function prolongates/injects coarse-face fluxes to finer
ones at resolution boundaries. Inspiration is taken from 
`halo_flux()` in `src/tree.h`.
 */

static void halo_flux_prolongate (vector * list) {
  vector * listv = NULL;
  for (vector v in list)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);
  if (listv) {
    for (int l = 0; l < depth(); l++) {
      foreach_halo (prolongation, l) {
	if (is_refined(neighbor(-1)) ){
	  int _i = 2*point.i - GHOSTS;
	  point.level++;
	  point.i = _i;
	  for (vector v in list)
	    v.x[] = coarse(v.x, 0);
	  point.i = (_i + GHOSTS)/2;
	  point.level--;
	}
	if (is_refined(neighbor(1))) {
	  int _i = 2*(point.i + 1) - GHOSTS;
	  point.level++;
	  point.i = _i;
	  for (vector v in list)
	    v.x[] = coarse(v.x, 0);
	  point.i = (_i + GHOSTS)/2 - 1;
	  point.level--;
	}
      }
    }
  }
}

/**
## Usage

It may be used in a future, ad-hoc, prove-of-the-point adaptive 1D advection-diffusion and [KdV](KdV.h) solver.
*/