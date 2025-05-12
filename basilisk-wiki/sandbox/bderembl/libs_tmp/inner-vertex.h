
// TODO: wil not generalize to MPI

@def foreach_inner_face_generic()
  OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = depth(); point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static))
  for (_k = GHOSTS + 1  ; _k <= point.n + GHOSTS - 1; _k++) {
    point.i = _k;
#if dimension > 1
    for (point.j = GHOSTS + 1; point.j <= point.n + GHOSTS - 1; point.j++)
#if dimension > 2
      for (point.k = GHOSTS + 1; point.k <= point.n + GHOSTS - 1; point.k++)
#endif
        {
#endif
	  POINT_VARIABLES
@
@def end_foreach_inner_face_generic()
#if dimension > 1
	}
#endif
  }
}
@

@def foreach_inner_vertex()
foreach_inner_face_generic() {  
  x -= Delta/2.;
#if dimension > 1  
  y -= Delta/2.;  
#endif
#if dimension > 2
  z -= Delta/2.;  
#endif
@
@define end_foreach_inner_vertex() } end_foreach_inner_face_generic()



@def foreach_inner_vertex_level(l)
OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = l; point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static))
  for (_k = GHOSTS + 1; _k < point.n + GHOSTS; _k++) {
    point.i = _k;
#if dimension > 1
    for (point.j = GHOSTS + 1; point.j < point.n + GHOSTS; point.j++)
#if dimension > 2
      for (point.k = GHOSTS + 1; point.k < point.n + GHOSTS; point.k++)
#endif
	{
#endif
          POINT_VARIABLES
  x -= Delta/2.;
#if dimension > 1  
  y -= Delta/2.;  
#endif
#if dimension > 2
  z -= Delta/2.;  
#endif
@
@def end_foreach_inner_vertex_level()
#if dimension > 1
	}
#endif
  }
}
@

@def foreach_vertex_level(l)
OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = l; point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static))
  for (_k = GHOSTS; _k < point.n + GHOSTS + 1; _k++) {
    point.i = _k;
#if dimension > 1
    for (point.j = GHOSTS; point.j < point.n + GHOSTS + 1; point.j++)
#if dimension > 2
      for (point.k = GHOSTS; point.k < point.n + GHOSTS + 1; point.k++)
#endif
	{
#endif
          POINT_VARIABLES
  x -= Delta/2.;
#if dimension > 1  
  y -= Delta/2.;  
#endif
#if dimension > 2
  z -= Delta/2.;  
#endif
@
@def end_foreach_vertex_level()
#if dimension > 1
	}
#endif
  }
}
@
