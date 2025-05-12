/** Injection and restriction functions for markers and color vertex.
  For the variables without interdependence. 
*/

static void myface_injection (Point point, vector v) {
  foreach_dimension(){
    // 0.5 for corner case where the interface passes through the cell vertex, should be test
    // fine(v.x, 0, 0) == 0. && fine(v.x, 0, 1) == 0. but with a markers at fine(v.x, 0, 1)
    // why it shows "stack smashing detected" when v.x[] = 0.5 ??? check this
    v.x[] = 0.;
    v.x[1] = 0.;
    if (fine(v.x, 0, 0) > 0.)
      v.x[] = fine(v.x, 0, 0) / 2.;

    if (fine(v.x, 0, 1) > 0.)
      v.x[] = 0.5 + fine(v.x, 0, 1) / 2.;
    
    if (fine(v.x, 2, 0) > 0.)
      v.x[1] = fine(v.x, 2, 0) / 2.;

    if (fine(v.x, 2, 1) > 0.)
      v.x[1] = 0.5 + fine(v.x, 2, 1) / 2.;
  }    
}

static inline void myrestriction_face (Point point, scalar s){
  myface_injection (point, s.v);
}

/** Functions for the scalar used in adapt_wavelet, which make sure that the cells near the interface are refined to maximum level. */

scalar mask_intf[]; // mask used for adapt_wavelet AMR
static inline void restriction_intf (Point point, scalar s) {
  s[] = 0.;
}

static inline void prolongation_intf (Point point, scalar s) {
  foreach_child()
    s[] = 0.;
}

// we need this to set the mask_intf to the correct value 
// at the resolution boundary between coarse and fine level leaf cells.
static inline void restriction_conf (Point point, scalar s){
  double conf_max = -HUGE;
  foreach_child()
    conf_max = max(conf_max, s[]);
  s[] = conf_max;
}

// Strange, I get a correct result after copying out this function from
// multigrid-common.h
// The function restrition_vertex is very special. The scalar bound with
// restriction_vertex will be treated specially by "tree_boundary_level" function
// in the tree-common.h. The behavior of automatic boundary condition is not what we want.
static inline void my_restriction_vertex (Point point, scalar s)
{
  for (int i = 0; i <= 1; i++) {
    s[i] = fine(s,2*i);
#if dimension >= 2  
    s[i,1] = fine(s,2*i,2);
#endif
#if dimension >= 3
    for (int j = 0; j <= 1; j++)
      s[i,j,1] = fine(s,2*i,2*j,2);
#endif
  }
}

// static inline void my_restriction_vertex_m (Point point, scalar s)
// {
//   for (int i = 0; i <= 1; i++) {
//     double c1 = fine(s, 2*i + 1), c2 = fine(s,2*i);
//     s[i] = c1 > 0. ? 0.5 * (1. + c1) : 0.5 * c2;
// #if dimension >= 2
//     c1 = fine(s, 2*i + 1, 2), c2 = fine(s,2*i, 2);
//     s[i,1] = c1 > 0. ? 0.5 * (1. + c1) : 0.5 * c2;
// #endif
// #if dimension >= 3
//     for (int j = 0; j <= 1; j++) {
//       c1 = fine(s, 2*i + 1, 2), c2 = fine(s,2*i, 2);
//       s[i,j,1] = fine(s,2*i,2*j,2);
//     }
// #endif
//   }
// }

static inline void my_restriction_vertex_zero (Point point, scalar s)
{
  for (int i = 0; i <= 1; i++) {
    s[i] = 0.;
#if dimension >= 2
    s[i, 1] = 0.;
#endif
#if dimension >= 3
    for (int j = 0; j <= 1; j++)
      s[i,j,1] = 0.;
#endif
  }
}

event defaults (i = 0) {
  mask_intf.restriction = restriction_intf;
  mask_intf.prolongation = prolongation_intf;
  mask_intf.refine = refine_injection;
}

void set_mask(scalar intf) {
  intf.restriction = restriction_conf;
  boundary({intf});
  foreach() {
    int with_itf = intf[] > 0.5;
    foreach_neighbor(1){
      with_itf = with_itf || (intf[] > 0.5);
    }

    mask_intf[] = with_itf ? 1. : 0.;
  }
  boundary({mask_intf});
}

foreach_dimension()
static void refine_face_injection_x (Point point, scalar s)
{
  vector v = s.v;
#if dimension <= 2
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {
    for (int j = 0; j <= 1; j++)
      fine(v.x,0,j) = v.x[];
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors &&
      (is_local(cell) || is_local(neighbor(1)))) {
    for (int j = 0; j <= 1; j++)
      fine(v.x,2,j) = v.x[1];
  }
  if (is_local(cell)) {
    for (int j = 0; j <= 1; j++)
      fine(v.x,1,j) = (v.x[] + v.x[1])/2.;
  }
#else // dimension > 2
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) {
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	      fine(v.x,0,j,k) = v.x[];
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors &&
      (is_local(cell) || is_local(neighbor(1)))) {
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	      fine(v.x,2,j,k) = v.x[1];
  }
  if (is_local(cell)) {
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	      fine(v.x,1,j,k) = (v.x[] + v.x[1])/2.;
  }
#endif // dimension > 2
}

static void my_prolongation_vertex (Point point, scalar s)
{
  fine(s,1,1,1) = (s[] + s[1] + s[0,1] + s[1,1] +
		   s[0,0,1] + s[1,0,1] + s[0,1,1] + s[1,1,1])/8.;

  for (int i = 0; i <= 1; i++) {
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
        if (allocated_child(2*i,2*j,2*k))
          fine(s,2*i,2*j,2*k) = s[i,j,k];

    foreach_dimension()
      if (neighbor(i).neighbors) {
        fine(s,2*i,1,1) = (s[i,0,0] + s[i,1,0] + s[i,0,1] + s[i,1,1])/4.;
        fine(s,2*i,1,0) = (s[i,0,0] + s[i,1,0])/2.;
        fine(s,2*i,0,1) = (s[i,0,0] + s[i,0,1])/2.;
        if (allocated_child(2*i,1,2))
          fine(s,2*i,1,2) = (s[i,0,1] + s[i,1,1])/2.;
        if (allocated_child(2*i,2,1))
          fine(s,2*i,2,1) = (s[i,1,0] + s[i,1,1])/2.;
      }
  }
}

/** All these are used for debugging*/
static void my_no_restriction_r (Point point, scalar s)
{
  s.dirty = false;
  // double ii, jj, kk;
  // ii = x / Delta;
  // jj = y / Delta;
  // kk = z / Delta;
  // if (x < 0.60 && x > 0.53 && y < 0.78 && y > 0.68 ) {
  //   printf("Refine: pid: %d, x:%g, y:%g, z:%g\n",\
  //     pid(), x, y, z);
  //   printf("pid: %d, x:%g, y:%g, z:%g,\
  //    s|x 1:%g 2:%g 3:%g 4:%g 5:%g 6:%g 7:%g 8:%g\n",\
  //     pid(), x / Delta, y / Delta, z / Delta, s[], s[1], s[0,1], s[1,1],
	// 	   s[0,0,1], s[1,0,1], s[0,1,1], s[1,1,1]);
  //   printf("\n");
  // }
}

static void my_no_restriction_p (Point point, scalar s)
{
  s.dirty = false;
  // double ii, jj, kk;
  // ii = x / Delta;
  // jj = y / Delta;
  // kk = z / Delta;
  // if (x < 0.6 && x > 0.53 && y < 0.78 && y > 0.68 && z < 0.78 && z > 0.68 ) {
  //   printf("Prolong: pid: %d, x:%g, y:%g, z:%g\n",\
  //     pid(), x, y, z);
  //   printf("pid: %d, x:%g, y:%g, z:%g,\
  //    s|x 1:%g 2:%g 3:%g 4:%g 5:%g 6:%g 7:%g 8:%g\n",\
  //     pid(), x / Delta, y / Delta, z / Delta, s[], s[1], s[0,1], s[1,1],
	// 	   s[0,0,1], s[1,0,1], s[0,1,1], s[1,1,1]);
  //   printf("\n");
  // }
}

static void my_no_restriction_c (Point point, scalar s)
{
  s.dirty = false;
  // double ii, jj, kk;
  // ii = x / Delta;
  // jj = y / Delta;
  // kk = z / Delta;
  // if (x < 0.60 && x > 0.53 && y < 0.78 && y > 0.68 ) {
  //   printf("Coarsen: pid: %d, x:%g, y:%g, z:%g\n",\
  //     pid(), x, y, z);
  //   printf("pid: %d, x:%g, y:%g, z:%g,\
  //    s|x 1:%g 2:%g 3:%g 4:%g 5:%g 6:%g 7:%g 8:%g\n",\
  //     pid(), x / Delta, y / Delta, z / Delta, s[], s[1], s[0,1], s[1,1],
	// 	   s[0,0,1], s[1,0,1], s[0,1,1], s[1,1,1]);
  //   printf("\n");
  // }
}


static void my_no_restriction_res (Point point, scalar s)
{
  s.dirty = false;
  // double ii, jj, kk;
  // ii = x / Delta;
  // jj = y / Delta;
  // kk = z / Delta;
  // if (x < 0.60 && x > 0.53 && y < 0.78 && y > 0.68 ) {
  //   printf("Restriction: pid: %d, x:%g, y:%g, z:%g\n",\
  //     pid(), x, y, z);
  //   printf("pid: %d, x:%g, y:%g, z:%g,\
  //    s|x 1:%g 2:%g 3:%g 4:%g 5:%g 6:%g 7:%g 8:%g\n",\
  //     pid(), x / Delta, y / Delta, z / Delta, s[], s[1], s[0,1], s[1,1],
	// 	   s[0,0,1], s[1,0,1], s[0,1,1], s[1,1,1]);
  //   printf("\n");
  // }
}
