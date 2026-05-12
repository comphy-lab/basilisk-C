

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
  /* if (point.i < GHOSTS +1 || point.i > point.n.x + GHOSTS - 1 || */
  /*     point.j < GHOSTS +1 || point.j > point.n.y + GHOSTS - 1) */
  if (x <= X0 + 0.5*Delta || x >= X0 + L0 - 0.5*Delta ||
    y <= Y0 + 0.5*Delta || y >= Y0 + L0 - 0.5*Delta)
    s[] = fine(s,0,0,0);  
  else  
    s[] = (fine(s,1,0,0) + 2*fine(s,0,0,0) + fine(s,-1,0,0) +  
           fine(s,0,1,0) + fine(s,0,-1,0))/6.;  
#endif  
}  

static inline double bilinear_vertex (Point point, scalar s)
{
#if dimension == 1
  bool on_x = ((point.i + GHOSTS) % 2 == 0);
  if (on_x)
    return coarse(s, 0);
  else
    return (coarse(s, 0) + coarse(s, 1))/2.;

#elif dimension == 2
  bool on_x = ((point.i + GHOSTS) % 2 == 0);
  bool on_y = ((point.j + GHOSTS) % 2 == 0);
  if (on_x && on_y)
    return coarse(s, 0, 0);
  else if (on_x)
    return (coarse(s, 0, 0) + coarse(s, 0, 1))/2.;
  else if (on_y)
    return (coarse(s, 0, 0) + coarse(s, 1, 0))/2.;
  else
    return (coarse(s, 0, 0) + coarse(s, 1, 0) +
            coarse(s, 0, 1) + coarse(s, 1, 1))/4.;

#elif dimension == 3
  bool on_x = ((point.i + GHOSTS) % 2 == 0);
  bool on_y = ((point.j + GHOSTS) % 2 == 0);
  bool on_z = ((point.k + GHOSTS) % 2 == 0);
  if (on_x && on_y && on_z)
    return coarse(s, 0, 0, 0);
  else if (!on_x && on_y && on_z)
    return (coarse(s, 0, 0, 0) + coarse(s, 1, 0, 0))/2.;
  else if (on_x && !on_y && on_z)
    return (coarse(s, 0, 0, 0) + coarse(s, 0, 1, 0))/2.;
  else if (on_x && on_y && !on_z)
    return (coarse(s, 0, 0, 0) + coarse(s, 0, 0, 1))/2.;
  else if (!on_x && !on_y && on_z)
    return (coarse(s, 0, 0, 0) + coarse(s, 1, 0, 0) +
            coarse(s, 0, 1, 0) + coarse(s, 1, 1, 0))/4.;
  else if (!on_x && on_y && !on_z)
    return (coarse(s, 0, 0, 0) + coarse(s, 1, 0, 0) +
            coarse(s, 0, 0, 1) + coarse(s, 1, 0, 1))/4.;
  else if (on_x && !on_y && !on_z)
    return (coarse(s, 0, 0, 0) + coarse(s, 0, 1, 0) +
            coarse(s, 0, 0, 1) + coarse(s, 0, 1, 1))/4.;
  else
    return (coarse(s, 0, 0, 0) + coarse(s, 1, 0, 0) +
            coarse(s, 0, 1, 0) + coarse(s, 1, 1, 0) +
            coarse(s, 0, 0, 1) + coarse(s, 1, 0, 1) +
            coarse(s, 0, 1, 1) + coarse(s, 1, 1, 1))/8.;
#endif
}
