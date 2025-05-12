/** Injection and restriction functions for markers and color vertex.
  For the variables without interdependence. 
*/

static void myface_injection (Point point, vector v) {
  foreach_dimension() {
    // 0.5 for corner case where the interface passes through the cell vertex, should be test
    // fine(v.x, 0, 0) == 0. && fine(v.x, 0, 1) == 0. but with a markers at fine(v.x, 0, 1)
    // why it shows "stack smashing detected" when v.x[] = 0.5 ??? check this
    v.x[] = 0.;
    v.x[1] = 0.;
    if (fine(v.x, 0, 0) > 0.)
      v.x[] = fine(v.x, 0, 0)/2.;

    if (fine(v.x, 0, 1) > 0.)
      v.x[] = 0.5 + fine(v.x, 0, 1)/2.;
    
    if (fine(v.x, 2, 0) > 0.)
      v.x[1] = fine(v.x, 2, 0)/2.;

    if (fine(v.x, 2, 1) > 0.)
      v.x[1] = 0.5 + fine(v.x, 2, 1)/2.;
  }    
}

static inline void myrestriction_face (Point point, scalar s) {
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
static inline void restriction_conf (Point point, scalar s) {
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
static inline void my_restriction_vertex (Point point, scalar s) {
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

/** For color vertex*/
static void refine_vertex_ebit (Point point, scalar ss) {
  fine(ss, 1, 1) = ss[];

  // should check the meaning of allocated_child function
  for (int i = 0; i <= 1; i++) {
    for (int j = 0; j <= 1; j++)
      if (allocated_child(2*i, 2*j))
        fine(ss,2*i,2*j) = ss[i,j];
    
      foreach_dimension()
        if (neighbor(i).neighbors) {
          fine(ss,2*i,1) = ss[i];
        }
  }
}

static inline void my_restriction_vertex_zero (Point point, scalar s) {
  for (int i = 0; i <= 1; i++) {
    s[i] = 0.;
#if dimension >= 2
    s[i,1] = 0.;
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

  mask_intf.nodump = true;
}

void set_mask (scalar intf, int n_iter = 1) {
  intf.restriction = restriction_conf;
  boundary ({intf});
  scalar mask_tmp[];
  mask_intf.restriction = restriction_conf;
  foreach()
    mask_intf[] = intf[];

  for (int iter = 0; iter < n_iter; iter++) {
    foreach() {
      bool is_mask = false;
      foreach_neighbor(1)
        if (mask_intf[] > 0.5)
          is_mask = true;

      mask_tmp[] = is_mask ? 1. : 0.;
    }

    foreach()
      mask_intf[] = mask_tmp[];

    boundary ({mask_intf});
  }
  mask_intf.restriction = restriction_intf;
}


/** Warper function of adapt_wavelet. Add the criteria of the EBIT method and 
  the embedded buondary. 
  1. For the EBIT method, the cells within the 3*3 stencil 
  of the interfacial cell are refined to maximum level.
  2. For embedded boundary, the cell in the solid region is not allowed to merge
  with the mix cell.
*/
trace
astats adapt_wavelet_ebit (scalar * slist,       // list of scalars
		      double * max,         // tolerance for each scalar
		      int maxlevel,         // maximum level of refinement
		      int minlevel = 1,     // minimum level of refinement
		      scalar * list = all)  // list of fields to update
{
  int ns = 0;
  scalar * ilist = list_copy(slist);

  for (scalar s in slist)
    ns++;

  ilist = list_add(ilist, mask_intf);
  #if EMBED
  double smax[ns + 2];
  smax[ns] = 0.02;
  smax[ns + 1] = 1.e-6;
  scalar mask_embed[];
  mask_embed.restriction = mask_embed.coarsen = restriction_conf;
  mask_embed.prolongation = mask_embed.refine = refine_injection;

  foreach() {
    mask_embed[] = cs[] > 0. ? 1. : 0.;
  }
  boundary ({mask_embed});
  ilist = list_add(ilist, mask_embed);
  #else
  double smax[ns + 1];
  smax[ns] = 0.02;
  #endif

  ns = 0;
  for (scalar s in ilist) {
    smax[ns] = max[ns];
    ns++;
  }
  
  astats st = {0, 0};
  st = adapt_wavelet (ilist, smax, maxlevel, minlevel, list);
  
  return st;
}