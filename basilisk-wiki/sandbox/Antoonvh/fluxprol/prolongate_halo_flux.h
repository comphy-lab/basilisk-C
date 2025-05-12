/**
# Flux prolongation 

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
In 3D we inject again. It could be 3rd-order as well on a 3x3x1 stencil
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

void restrict_face_level (vector * list, int l) {
  foreach_coarse_level(l) { 
    foreach_dimension()  
      if (!is_leaf(neighbor(-1)))
	for (vector v in list)
	  v.x[] = (fine(v.x, 0, 0, 0) + fine(v.x, 0, 1, 0))/2.; 
  }
}

void boundary_face_level (vector * list, int l) {
  for (vector v in list)
    boundary_iterate(level, (scalar*){v}, l);
}

static void restrict_face (vector * list) {
  for (int l = depth(); l > 0; l--) {
    if (l != depth())
      restrict_face_level (list, l);
    boundary_face_level (list, l);
  }
}

trace
static void flux_prolongate (vector * list) {
  vector * listv = NULL;
  for (vector v in list)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);
  if (listv) {
    restrict_face (listv);
    for (int l = 0; l < depth(); l++)
      halo_flux_prolongate_level (listv, l);
  }
  free (listv);
}

static void dummy (vectorl list) {
  ;
}

static void multigrid_restriction_level (scalar * list, int below_level) {
  for (int l = below_level -1; l >= 0; l--) {
    foreach_coarse_level(l, nowarning) {
      for (scalar s in list) {
	foreach_block()
	  s.restriction (point, s);
      }
    }
    boundary_iterate (level, list, l); 
  }
}

/**
## Usage

* [A Quadtree implementation](myQT.h)
*/
