/**
# A flexible adaptation algorithm

A special purpose adaptation algorithm. consider
using `adapt_wavelet()` instead.
*/
struct Adapt_Field {
  scalar f;       // The field with leaf-filled values
  double max;     // refinement criterion
  double min;     // coarseing criterion 
  int maxlevel;   // maximum level of refinement
  int minlevel;   // minimum level of refinement (default 1)
  scalar * list;  // list of fields to update (default all)
};

astats adapt_field (struct Adapt_Field p);

extern vector u;

struct Adapt_Flow {
  double ue;
  int maxlevel;
  double p;
  double cfac;
};

attribute {
  void (* interpolant) (Point, scalar);
}

void wavelet_interpolant (scalar s, scalar w)
{
  restriction ({s});
  for (int l = grid->maxdepth - 1; l >= 0; l--) {
    foreach_coarse_level (l) {
      foreach_child()
        w[] = s[];
      s.interpolant (point, s);
      foreach_child() {
        double sp = s[];
        s[] = w[];
        /* difference between fine value and its prolongation */
        w[] -= sp;
      }
    }
    boundary_level ({w}, l + 1);
  }
  /* root cell */
  foreach_level(0) 
    w[] = s[];
  boundary_level ({w}, 0);
}

astats adapt_flow (struct Adapt_Flow pa) {
  double ue = pa.ue;
  int maxlevel = pa.maxlevel;
  double p = pa.p;
  double uec = ue/1.5;
  if (pa.cfac)
    uec = ue/pa.cfac; 
  scalar w[];
  vector wv[];
  foreach_dimension() 
    wavelet_interpolant (u.x, wv.x);
    
  for (int l = 1; l < grid->maxdepth; l++)
    foreach_coarse_level (l) { 
      double maxw = 0.;
      foreach_child() {
	double a = 0;
	foreach_dimension()
	  a += fabs(wv.x[]);
	maxw = max(maxw, a);
      }
      w[] = maxw;
    }
#if _MPI
  w.restriction = no_restriction;
  w.prolongation = no_restriction;
  boundary ({w});
#endif
  for (int l = grid->maxdepth; l > 1; l--)  
    foreach_level(l)  {// 3x3(x3) Gaussian kernel
#if (dimension == 2)
      w[] = pow(Delta, p)*(1.  * coarse(w,0,0,0)   +
			   .5  *(coarse(w,1,0,0)   + coarse(w,0,1,0)   +
				 coarse(w,-1,0,0)  + coarse(w,0,-1,0)) +
			   0.25*(coarse(w,1,1,0)   + coarse(w,1,-1,0)  +
				 coarse(w,-1,-1,0) + coarse(w,-1,1,0)))/4.;
#else //(dimension == 3) 
      w[] = pow(Delta, p)*(1.   * coarse(w,0,0,0)   +
			   .5   *(coarse(w,1,0,0)   + coarse(w,-1,0,0)  +
				  coarse(w,0,1,0)   + coarse(w,0,-1,0)  +
				  coarse(w,0,0,1)   + coarse(w,0,0,-1)) +
			   0.25 *(coarse(w,1,1,0)   + coarse(w,1,-1,0)  +
				  coarse(w,-1,-1,0) + coarse(w,-1,1,0)  +
				  coarse(w,0,1,1)   + coarse(w,0,-1,1)  +
				  coarse(w,0,-1,-1) + coarse(w,0,1,-1)  +
				  coarse(w,1,0,1)   + coarse(w,1,0,-1)  +
				  coarse(w,-1,0,-1) + coarse(w,-1,0,1)) +
			   0.125*(coarse(w,1,1,1)   + coarse(w,-1,1,1)  +
				  coarse(w,1,-1,1)  + coarse(w,1,1,-1)  +
				  coarse(w,-1,-1,1) + coarse(w,-1,1,-1) +
				  coarse(w,1,-1,-1) + coarse(w,-1,-1,-1)))/8.;
#endif
    }
  return adapt_field (w, ue, uec, maxlevel);
}



struct Adapt_List {
  scalar * list;
  double * ues;
  int maxlevel;
  double p;
  double fac;
};

astats adapt_list (struct Adapt_List pa) {
  scalar * list = pa.list;
  double * ues = pa.ues;
  int maxlevel = pa.maxlevel;
  double p = pa.p;
  if (!pa.fac)
    pa.fac = 1.5;
  scalar w[], * wl = list_clone (list);
 scalar s, ws;
  for (s,ws in list,wl)
    wavelet (s, ws);
  for (int l = 1; l < depth() - 1; l++)
    foreach_coarse_level (l) { 
      double maxw = 0.;
      foreach_child() {
	double a = 0;
	scalar s, ws;
	NOT_UNUSED(s);
	int is = 0;
	for (s, ws in list, wl) {
	  a += fabs(ws[])/ues[is]; //scale w with ue
	  is++;
	}
	maxw = max(maxw, a);
      }
      w[] = maxw;
    }
#if _MPI
  w.restriction = no_restriction;
  w.prolongation = no_restriction;
  boundary ({w});
#endif
  for (int l = depth(); l > 1; l--)  
    foreach_level(l)  {// 3x3(x3) Gaussian kernel
#if (dimension == 2)
      w[] = pow(Delta, p)*(1.  * coarse(w,0,0,0)   +
			   .5  *(coarse(w,1,0,0)   + coarse(w,0,1,0)   +
				 coarse(w,-1,0,0)  + coarse(w,0,-1,0)) +
			   0.25*(coarse(w,1,1,0)   + coarse(w,1,-1,0)  +
				 coarse(w,-1,-1,0) + coarse(w,-1,1,0)))/4.;
#else //(dimension == 3) 
      w[] = pow(Delta, p)*(1.   * coarse(w,0,0,0)   +
			   .5   *(coarse(w,1,0,0)   + coarse(w,-1,0,0)  +
				  coarse(w,0,1,0)   + coarse(w,0,-1,0)  +
				  coarse(w,0,0,1)   + coarse(w,0,0,-1)) +
			   0.25 *(coarse(w,1,1,0)   + coarse(w,1,-1,0)  +
				  coarse(w,-1,-1,0) + coarse(w,-1,1,0)  +
				  coarse(w,0,1,1)   + coarse(w,0,-1,1)  +
				  coarse(w,0,-1,-1) + coarse(w,0,1,-1)  +
				  coarse(w,1,0,1)   + coarse(w,1,0,-1)  +
				  coarse(w,-1,0,-1) + coarse(w,-1,0,1)) +
			   0.125*(coarse(w,1,1,1)   + coarse(w,-1,1,1)  +
				  coarse(w,1,-1,1)  + coarse(w,1,1,-1)  +
				  coarse(w,-1,-1,1) + coarse(w,-1,1,-1) +
				  coarse(w,1,-1,-1) + coarse(w,-1,-1,-1)))/8.;
#endif
    }
  return adapt_field (w, 1., 1./pa.fac, maxlevel); //w is already scaled
}



trace
astats adapt_field (struct Adapt_Field p) {
  scalar f = p.f;
  if (p.list == NULL)
    p.list = all;
  scalar * listc = NULL;
  for (scalar s in p.list)
    if (!is_constant(s) && s.restriction != no_restriction)
      listc = list_add (listc, s);
  astats st = {0, 0};
  // refinement
  if (p.minlevel < 1)
    p.minlevel = 1;
  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  foreach_cell() {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
  	if (cell.flags & too_coarse) {
  	  cell.flags &= ~too_coarse;
  	  refine_cell (point, listc, refined, &tree->refined);
  	  st.nf++;
  	}
  	continue;
      }
      else { // !is_leaf (cell)
  	if (cell.flags & refined) {
  	  // cell has already been refined, skip its children
  	  cell.flags &= ~too_coarse;
  	  continue;
  	}
  	// check whether the cell or any of its children is local
  	bool local = is_local(cell);
  	if (!local)
  	  foreach_child()
  	    if (is_local(cell)) {
  	      local = true;
	      break;
	    }
  	if (local) {
  	  static const int just_fine = 1 << (user + 3);
	  double max = p.max, min = p.min;
	  foreach_child() {
	    double e = fabs(f[]);
	    if (e > max && level < p.maxlevel) {
	      cell.flags &= ~too_fine;
	      cell.flags |= too_coarse;
	    }
	    else if ((e <= min || level > p.maxlevel) &&
		     !(cell.flags & (too_coarse|just_fine))) {
	      if (level >= p.minlevel)
		cell.flags |= too_fine;
	    }
	    else if (!(cell.flags & too_coarse)) {
	      cell.flags &= ~too_fine;
	      cell.flags |= just_fine;
	    }
	  }
	  foreach_child() {
  	    cell.flags &= ~just_fine;
  	    if (!is_leaf(cell)) {
  	      cell.flags &= ~too_coarse;
  	      if (level >= p.maxlevel)
  		cell.flags |= too_fine;
  	    }
  	    else if (!is_active(cell))
  	      cell.flags &= ~too_coarse;
  	  }
  	}
      }
    }
    else // inactive cell
      continue;
  }
  mpi_boundary_refine (listc);
  
  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth(); l >= 0; l--) {
    foreach_cell()
      if (!is_boundary(cell)) {
  	if (level == l) {
  	  if (!is_leaf(cell)) {
  	    if (cell.flags & refined)
  	      // cell was refined previously, unset the flag
  	      cell.flags &= ~(refined|too_fine);
  	    else if (cell.flags & too_fine) {
  	      if (is_local(cell) && coarsen_cell (point, listc))
  		st.nc++;
  	      cell.flags &= ~too_fine; // do not coarsen parent
  	    }
  	  }
  	  if (cell.flags & too_fine)
  	    cell.flags &= ~too_fine;
  	  else if (level > 0 && (aparent(0).flags & too_fine))
  	    aparent(0).flags &= ~too_fine;
  	  continue;
  	}
  	else if (is_leaf(cell))
  	  continue;
      }
    mpi_boundary_coarsen (l, too_fine);
    }
  free (listc);
  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (p.list);
  return st;
}
