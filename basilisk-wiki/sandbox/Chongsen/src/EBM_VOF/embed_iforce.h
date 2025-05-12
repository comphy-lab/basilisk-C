

attribute {
  scalar phi;
}


event defaults (i = 0) {  
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

event acceleration (i++)
{
  

  scalar * list = NULL;
  for (scalar f in interfaces)
  for (scalar tmp_c in tmp_interfaces)
    if (f.phi.i) {
      list = list_add (list, f);

      foreach(){
	  f[] = clamp (f[], 0., 1.);
	  tmp_c[] = clamp (tmp_c[], 0., 1.);
      }
    }

  
#if TREE
  for (scalar f in list) {
    f.prolongation = p.prolongation;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif


  face vector ia = a;
  for (scalar tmp_c in tmp_interfaces)
  foreach_face()
    for (scalar f in list)
      if (tmp_c[] != tmp_c[-1] && fm.x[] > 0.) {
	
	scalar phi = f.phi;
	double phif =
	  (phi[] < nodata && phi[-1] < nodata) ?
	  (phi[] + phi[-1])/2. :
	  phi[] < nodata ? phi[] :
	  phi[-1] < nodata ? phi[-1] :
	  0.;
	ia.x[] += alpha.x[]/(fm.x[] + SEPS)*phif*(tmp_c[] - tmp_c[-1])/Delta;
      }

  
#if TREE
  for (scalar f in list) {
    f.prolongation = fraction_refine;
    f.dirty = true; // boundary conditions need to be updated
  }
#endif

  for (scalar f in list) {
    scalar phi = f.phi;
    delete ({phi});
    f.phi.i = 0;
  }
  free (list);
}
