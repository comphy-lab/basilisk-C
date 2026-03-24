

macro2 foreach_vertex_level (int l, char flags = 0, Reduce reductions = None) {
  OMP_PARALLEL(reductions) {
    int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
    Point point = {0};
    point.level = l;
    SET_DIMENSIONS();
    int _k;
    OMP(omp for schedule(static))
      for (_k = GHOSTS; _k <= point.n.x + GHOSTS; _k++) {
        point.i = _k;
#if dimension > 1
        for (point.j = GHOSTS; point.j <= point.n.y + GHOSTS; point.j++)
#if dimension > 2
          for (point.k = GHOSTS; point.k <= point.n.z + GHOSTS; point.k++)
#endif
            {
#endif
              POINT_VARIABLES();
              x -= Delta/2.;
#if dimension > 1
              y -= Delta/2.;
#endif
#if dimension > 2
              z -= Delta/2.;
#endif
              {...}

#if dimension > 1
            }
#endif
      }
  }
}


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
