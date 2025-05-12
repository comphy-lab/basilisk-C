
// TODO: wil not generalize to MPI

@def foreach_inner_face_generic()
  OMP_PARALLEL() {

  int index_left = 1, index_right = 1;
  int index_top = 1, index_bottom = 1;
  int index_front = 1, index_back = 1;
#if _MPI
  if (mpi_coords[0] != 0) index_left = 0;
  if (mpi_coords[0] != mpi_dims[0] - 1) index_right = 0;
  if (mpi_coords[1] != 0) index_bottom = 0;
  if (mpi_coords[1] != mpi_dims[1] - 1) index_top = 0;
#if dimension > 2
  if (mpi_coords[2] != 0) index_front = 0;
  if (mpi_coords[2] != mpi_dims[2] - 1) index_back = 0;
#endif
#endif

  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = depth(); point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static))
  for (_k = GHOSTS + index_left  ; _k <= point.n + GHOSTS - index_right; _k++) {
    point.i = _k;
#if dimension > 1
    for (point.j = GHOSTS + index_bottom; point.j <= point.n + GHOSTS - index_top; point.j++)
#if dimension > 2
      for (point.k = GHOSTS + index_front; point.k <= point.n + GHOSTS - index_back; point.k++)
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

  int index_left = 1, index_right = 1;
  int index_top = 1, index_bottom = 1;
  int index_front = 1, index_back = 1;
#if _MPI
  if (mpi_coords[0] != 0) index_left = 0;
  if (mpi_coords[0] != mpi_dims[0] - 1) index_right = 0;
  if (mpi_coords[1] != 0) index_bottom = 0;
  if (mpi_coords[1] != mpi_dims[1] - 1) index_top = 0;
#if dimension > 2
  if (mpi_coords[2] != 0) index_front = 0;
  if (mpi_coords[2] != mpi_dims[2] - 1) index_back = 0;
#endif
#endif

  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0};
  point.level = l; point.n = 1 << point.level;
  int _k;
  OMP(omp for schedule(static))
  for (_k = GHOSTS + index_left; _k <= point.n + GHOSTS - index_right; _k++) {
    point.i = _k;
#if dimension > 1
    for (point.j = GHOSTS + index_bottom; point.j <= point.n + GHOSTS - index_top; point.j++)
#if dimension > 2
      for (point.k = GHOSTS + index_front; point.k <= point.n + GHOSTS - index_back; point.k++)
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


@ifndef foreach_coarse_vertex_level
@ define foreach_coarse_vertex_level      foreach_level
@ define end_foreach_coarse_vertex_level  end_foreach_level
@endif
