/**
# My Vertex methods

For finite-difference applications on trees, the requirement of an
exact restriction operator warrants the use of vertex-points instead
of cell centers. Furtheromore, for 5-point-stencil schemes, not all
"Basilisk vertices" can be vertices.

## User-interface vertex iterator

Iterate leafs and vertices on levels
 */
@def foreach_vert()
  update_cache();
  foreach_cache(tree->leaves) {
    x -= Delta/2.;
#if dimension >= 2
    y -= Delta/2.;
#endif
#if dimension >= 3
    z -= Delta/2.;
#endif
@
@define end_foreach_vert() } end_foreach_cache()

@def foreach_vert_level(l) {
  if (l <= depth()) {
    update_cache();
    CacheLevel _active = tree->active[l];
    foreach_cache_level (_active,l) {
      x -= Delta/2.;
#if dimension >= 2
      y -= Delta/2.;
#endif
#if dimension >= 3
      z -= Delta/2.;
#endif
@
@define end_foreach_vert_level() } end_foreach_cache_level(); }}
/**
## Vertex restriction

Restrict only the local `point` instead of $2^{\mathtt{dimension}}$
 */
static inline void restriction_vert (Point point, scalar s) {
  s[] = fine(s,0,0,0);
}
/**
   Or a coarse estimate (not in 3D);
 */
static inline void restriction_coarsen_vert (Point point, scalar s) {
#if (dimension == 1)
  s[] = (fine(s,1,0,0) + 2*fine(s,0,0,0) + fine(s,-1,0,0))/4.;
#elif (dimension == 2)
  s[] = (fine(s,1,0,0) + 2*fine(s,0,0,0) + fine(s,-1,0,0) +
	 fine(s,0,1,0) + fine(s,0,-1,0))/6.;
#endif
}
/**
## High-Level Boundary treatment

2nd-order
 */
static inline void refine_vert (Point point, scalar s) { 
  // Injection
  fine (s,0,0,0) = s[];
  // Vertices with two nearest coarse neighbors
  fine (s,1,0,0) = (s[] + s[1])/2.;
#if dimension > 1
  fine (s,0,1,0) = (s[] + s[0,1])/2.;
#endif
#if dimension > 2
  fine (s,0,0,1) = (s[] + s[0,0,1])/2.;
#endif
  // Vertices with four nearest coarse neighbors
#if dimension > 1
  fine(s,1,1,0) = (s[0] + s[1] + s[0,1] + s[1,1])/4.;
#endif
#if dimension > 2 // dimension == 3
  fine(s,1,0,1) = (s[] + s[1] + s[0,0,1] + s[1,0,1])/4.;
  fine(s,0,1,1) = (s[] + s[0,1] + s[0,1] + s[0,1,1])/4.;
  // In 3D, there is a vertex with 8 nearest coarse neighbors
  fine(s,1,1,1) = (s[] + s[1,1,1] +
		   s[1] + s[0,1] + s[0,0,1] +
		   s[1,1] + s[0,1,1] + s[1,0,1])/8.;
#endif
}


static inline void prolongation_vert (Point point, scalar s) {
  // do not inject for error estimation
  refine_vert(point, s);
  //if (!is_leaf(cell)) 
      fine(s,0,0,0) = (s[-1] + s[1]
#if dimension > 1
		       + s[0,1] + s[0,-1]
#endif
		       )/(2*dimension);
}

/**
   4th-order
 */
static inline void refine_vert4 (Point point, scalar s) { 
  // Injection
  fine (s,0,0,0) = s[];
  // Vertices with two nearest coarse neighbors
  fine (s,1,0,0) = (9*(s[0] + s[1]) - (s[-1] + s[2]))/16.;
#if dimension > 1
  fine (s,0,1,0) = (9*(s[0] + s[0,1]) - (s[0,-1] + s[0,2]))/16.;
#endif
#if dimension > 2
  fine (s,0,0,1) = (9*(s[0] + s[0,0,1]) - (s[0,0,-1] + s[0,0,2]))/16.;
#endif
  // Vertices with four nearest coarse neighbors
#if dimension > 1
  fine(s,1,1,0) = (9.*(s[0] + s[1,1] + s[0,1] + s[1,0]) -
		   (s[-1,-1] + s[2,2] + s[-1,2] + s[2,-1]))/32.;
#endif
  // not 3D yet...
#if dimension > 2 // dimension == 3
  fine(s,1,0,1) = (s[] + s[1] + s[0,0,1] + s[1,0,1])/4.;
  fine(s,0,1,1) = (s[] + s[0,1] + s[0,1] + s[0,1,1])/4.;
  // In 3D, there is a vertex with 8 nearest coarse neighbors
  fine(s,1,1,1) = (s[] + s[1,1,1] +
		   s[1] + s[0,1] + s[0,0,1] +
		   s[1,1] + s[0,1,1] + s[1,0,1])/8.;
#endif
}

static inline void prolongation_vert4 (Point point, scalar s) { 
  refine_vert4(point, s);
  if (!is_leaf(cell)) {
    fine(s,0,0,0) = (4*(fine(s,-1,0,0) + fine(s,1,0,0)) -
		     (fine(s,-2,0,0) + fine(s,2,0,0))
#if dimension > 1
		     +4* (fine(s,0,-1,0) + fine(s,0,1,0)) -
		     (fine(s,0,-2,0) + fine(s,0,2,0))
#endif
		     )/(dimension*6);
  }
}

foreach_dimension()
double interp4_x (Point point, scalar s) {
  return (4*(s[-1] + s[1])- (s[-2] + s[2]))/6.;
}

/**
Biassed 5th order
*/
double interp5_x (Point point, scalar s) {
  return ((3./128)*coarse(s,-2,0,0) - (5./32)*coarse(s,-1,0,0) +
	  (45./64.)*coarse(s,0,0,0) + (15./32.)*coarse(s,1,0,0) -
	  (5./128.)*coarse(s,2,0,0));
}

double interp5_y (Point point, scalar s) {
  return ((3./128)*coarse(s,0,-2,0) - (5./32)*coarse(s,0,-1,0) +
	  (45./64.)*coarse(s,0,0,0) + (15./32.)*coarse(s,0,1,0) -
	  (5./128.)*coarse(s,0,2,0));
}

static inline void refine_vert5 (Point point, scalar s) { 
  // Injection
  fine (s,0,0,0) = s[];
  // Vertices with two nearest coarse neighbors
  fine (s,1,0,0) = ((3./128)*s[-2] - (5./32)*s[-1] + (45./64.)*s[]
		    + (15./32.)*s[1] - (5./128.)*s[2]);
#if dimension > 1
  fine (s,0,1,0) = ((3./128)*s[0,-2] - (5./32)*s[0,-1] + (45./64.)*s[]
		    + (15./32.)*s[0,1] - (5./128.)*s[0,2]);
#endif
#if dimension > 2
  fine (s,0,0,1) = ((3./128)*s[0,0,-2] - (5./32)*s[0,0,-1] + (45./64.)*s[]
		    + (15./32.)*s[0,0,1] - (5./128.)*s[0,0,2]);
#endif
  // Vertices with four nearest coarse neighbors
#if dimension > 1
  fine (s,1,1,0) = ((3./128)*s[-2,-2] - (5./32)*s[-1,-1] + (45./64.)*s[] +
		    (15./32.)*s[1,1] - (5./128.)*s[2,2]);
#endif
#if dimension > 2 // dimension == 3
  fine(s,1,0,1) = ((3./128)*s[-2,0,-2] - (5./32)*s[-1,0,-1] + (45./64.)*s[] +
		   (15./32.)*s[1,0,1] - (5./128.)*s[2,0,2]);
  fine(s,0,1,1) = ((3./128)*s[0,-2,-2] - (5./32)*s[0,-1,-1] + (45./64.)*s[] +
		   (15./32.)*s[0,1,1] - (5./128.)*s[0,2,2]);
  // In 3D, there is a vertex with 8 nearest coarse neighbors
  fine(s,1,1,1) = ((3./128)*s[-2,-2,-2] - (5./32)*s[-1,-1,-1] + (45./64.)*s[] +
		   (15./32.)*s[1,1,1] - (5./128.)*s[2,2,2]);
#endif
}



/**
## To do

* Other than periodic box boundaries
 */

    
