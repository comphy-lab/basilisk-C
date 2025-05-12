/**
# Flux prolongation 

Cell-face fluxes may be approximated from a coarser-level
solution. This method may be attractive for the computation of fluxes
at resolution boundaries. It is the 'upside-down' version of
`halo_flux()`, reducing the importance of the accuracy of
halo-cell-centered values.

The function prolongates/injects coarse-face fluxes to finer ones at
resolution boundaries. Inspiration is taken from `src/tree.h`.
 */

void halo_flux_prolongate_level (vector * list, int l) {
  foreach_halo (prolongation, l) {
    foreach_dimension() {
#if dimension == 1
      /**
In 1D injectiontion is exact. 
       */
      if (is_refined(neighbor(-1)))
	for (vector v in list)
	  fine(v.x, 0, 0 ,0) = v.x[];
      if (is_refined(neighbor(1)))
	for (vector v in list)
	  fine(v.x, 2, 0 ,0) = v.x[1];
#elif dimension == 2
      /**
In 2D we use a 3rd-order-accurate conservative formulation.
       */
      if (is_refined(neighbor(-1))) {
	for (vector v in list) {
	  double dv = (v.x[0,1] - v.x[0,-1])/8.;
	  for (int j = 0; j <= 1; j++)
	    fine(v.x, 0, j, 0) = j == 0 ? v.x[] - dv : v.x[] + dv;
	}
      }
      if (is_refined(neighbor(1))) {
	for (vector v in list) {
	  double dv = (v.x[1,1] - v.x[1,-1])/8.;
	  for (int j = 0; j<=1; j++)
	    fine(v.x, 2, j, 0) =  j == 0 ? v.x[1] - dv : v.x[1] + dv;
	}
      }
#elif dimension == 3
      /**
In 3D we inject again. We will update with 3x3x3 linear
       */
      if (is_refined(neighbor(-1))) {
	for (vector v in list) {
	  fine(v.x, 0, 0, 0) = v.x[];
	  fine(v.x, 0, 1, 0) = v.x[];
	  fine(v.x, 0, 0, 1) = v.x[];
	  fine(v.x, 0, 1, 1) = v.x[];
	}
      }
      if (is_refined(neighbor(1))) { 
	for (vector v in list) {
	  fine(v.x, 2, 0, 0) = v.x[1];
	  fine(v.x, 2, 1, 0) = v.x[1];
	  fine(v.x, 2, 0, 1) = v.x[1];
	  fine(v.x, 2, 1, 1) = v.x[1];
	}
      }
#endif
    }
  }
}

static void halo_flux_prolongate (vector * list) {
  vector * listv = NULL;
  for (vector v in list)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);
  if (listv) 
    for (int l = 0; l < depth(); l++)
      halo_flux_prolongate_level (listv, l);
}

/**
## Usage

It may be used in a future, ad-hoc, prove-of-the-point adaptive advection-diffusion, [KdV](KdV.h) solver and a Poisson solver.
*/
