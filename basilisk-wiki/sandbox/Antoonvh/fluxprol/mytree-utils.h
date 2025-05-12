/**
# My tree utilities

Some functions that make life without halo-ghosts possible

## Centered-gradient computation

`my_gradients()` mimics `gradients()`, and functions like
[`centered_gradient()`](http://basilisk.fr/src/navier-stokes/centered.h#approximate-projection)
 */
void my_gradients (scalar * sl, vector * gl) {
  assert (list_len(sl) == vectors_len(gl));
  scalar s;
  vector g;
  for (s, g in sl, gl) {
    face vector gf[];
    foreach_face() 
      gf.x[] = (s[] - s[-1])/Delta;
    flux_prolongate ({gf});
    foreach() {
      foreach_dimension() 
	g.x[] = (gf.x[] + gf.x[1])/2.;
    }
  }
}

/**
# No prolongation

We can test if ghost-cells values are not inadvertently prolongated.
 */
static inline void no_prolongation (Point point, scalar s) {
  fputs ("# Error: No prolongation\n", stderr);
  assert (0);
}

/**
# My `tree_boundary_level()`

This is a re-implementation of `tree_boundary_level`, that does not
prolongate scalar values, it does do restriction!
 */
static void my_tree_boundary_level (scalar * list, int l) {
  
  int depth = l < 0 ? depth() : l;
  if (tree_is_full()) {
    boundary_iterate (level, list, depth);
    return;
  }
  /** 
Compose lists
  */
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL, * vlist = NULL;
  for (scalar s in list) 
    if (!is_constant (s)) {
      if (s.restriction == restriction_average) {
	listdef = list_add (listdef, s);
	list2 = list_add (list2, s);
      }
      else if (s.restriction != no_restriction) {
	listc = list_add (listc, s);
	if (s.face)
	  foreach_dimension()
	    list2 = list_add (list2, s.v.x);
	else {
	  list2 = list_add (list2, s);
	  if (s.restriction == restriction_vertex)
	    vlist = list_add (vlist, s);
	}
      }
    }

  if (vlist) // vertex scalars
#if dimension == 1
    foreach_vertex()
      if (is_refined(cell) || is_refined(neighbor(-1)))
	for (scalar s in vlist)
	  s[] = is_vertex (child(0)) ? fine(s) : nodata;
#elif dimension == 2
    foreach_vertex() {
      if (is_refined(cell) || is_refined(neighbor(-1)) ||
	  is_refined(neighbor(0,-1)) || is_refined(neighbor(-1,-1))) {
	// corner
	for (scalar s in vlist)
	  s[] = is_vertex (child(0)) ? fine(s) : nodata;
      }
      else
	foreach_dimension()
	  if (child.y == 1 &&
	      (is_prolongation(cell) || is_prolongation(neighbor(-1)))) {
	    // center of refined edge
	    for (scalar s in vlist)
	      s[] = is_vertex(neighbor(0,-1)) && is_vertex(neighbor(0,1)) ?
		(s[0,-1] + s[0,1])/2. : nodata;
	  }
    }
#else // dimension == 3
    foreach_vertex() {
      if (is_refined(cell) || is_refined(neighbor(-1)) ||
	  is_refined(neighbor(0,-1)) || is_refined(neighbor(-1,-1)) ||
	  is_refined(neighbor(0,0,-1)) || is_refined(neighbor(-1,0,-1)) ||
	  is_refined(neighbor(0,-1,-1)) || is_refined(neighbor(-1,-1,-1))) {
	// corner
	for (scalar s in vlist)
	  s[] = is_vertex (child(0)) ? fine(s) : nodata;
      }
      else
	foreach_dimension() {
	  if (child.y == 1 && child.z == 1 &&
	      (is_prolongation(cell) || is_prolongation(neighbor(-1)))) {
	    // center of refined face
	    for (scalar s in vlist)
	      s[] = is_vertex(neighbor(0,-1,-1)) && is_vertex(neighbor(0,1,-1))
		&& is_vertex(neighbor(0,-1,1)) && is_vertex(neighbor(0,1,1)) ?
		(s[0,-1,-1] + s[0,1,-1] + s[0,-1,1] + s[0,1,1])/4. : nodata;
	  }
	  else if (child.x == -1 && child.z == -1 && child.y == 1 &&
		   (is_prolongation(cell) || is_prolongation(neighbor(-1)) ||
		    is_prolongation(neighbor(0,0,-1)) ||
		    is_prolongation(neighbor(-1,0,-1)))) {
	    // center of refined edge
	    for (scalar s in vlist)
	      s[] = is_vertex(neighbor(0,-1)) && is_vertex(neighbor(0,1)) ?
		(s[0,-1] + s[0,1])/2. : nodata;
	  }
	}
    }
#endif // dimension == 3
  free (vlist);
  
  if (listdef || listc) {
    boundary_iterate (restriction, list2, depth);
    for (int l = depth - 1; l >= 0; l--) {
      foreach_coarse_level(l) {
	for (scalar s in listdef)
	  restriction_average (point, s);
	for (scalar s in listc)
	  s.restriction (point, s);
      }
      
      boundary_iterate (restriction, list2, l);
    }
    
    free (listdef);
    free (listc);
    free (list2);
  }
  /**
There used to be code here

~~~literatec
...
~~~
   */
  
}

/**
## Face prolongation

 */
double prol_face_linear_x(Point point, scalar s) {
  vector v = s.v;
  assert (child.x == -1);
  double g1 = (coarse(v.x,0,1,0) - coarse(v.x,0,-1,0))/8.;
  return coarse (v.x,0,0,0) + child.y*g1;
}

double prol_face_linear_y(Point point, scalar s) {
  vector v = s.v;
  assert (child.y == -1);
  double g1 = (coarse(v.y,1,0,0) - coarse(v.y,-1,0,0))/8.;
  return coarse(v.y,0,0,0) + child.x*g1;
}

void prolongate_faces_level (face vector f, int l) {
  foreach_face_level(l) {
    foreach_dimension() {
      is_face_x() {
	if (child.x == -1)
	  f.x[] = prol_face_linear_x(point, f.x);
      }
    }
  }
  foreach_face_level(l) {
    foreach_dimension() {
      if (child.x == 1)
	f.x[] = (f.x[-1] + f.x[1])/2.;
    }
  }
}

/**
## restriction on levels
 
 */


