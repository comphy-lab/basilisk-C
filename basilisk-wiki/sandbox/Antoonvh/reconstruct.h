/**
# Dynamic reconstruction of a geometry defined by `GEOM(x,y)`.

A (dis)proof-of-concept in 2D.
 */
#if TREE

/**
The idea is to compute the 2D face fractions `fs` and then compute
the volume fraction...
 */

static inline void recompute_face_fraction_refine_2D_x (Point point, scalar s) {
  vector v = s.v;
  if (!is_refined(neighbor(-1)) &&
      (is_local(cell) || is_local(neighbor(-1)))) { //lhs
    double ll = GEOM (x - Delta/2, y - Delta/2);
    double ml = GEOM (x - Delta/2, y);
    double ul = GEOM (x - Delta/2, y + Delta/2);
    double fl = 0.5, fu = 0.5;
    if (ll != ml)
      fl = ll / (ll - ml);
    if (ml != ul)
      fu = ml / (ml - ul);
    if (ll * ml < 0) 
      fine(v.x,0,0,0) = ll < 0 ? 1. - fl : fl ;
    else
      fine(v.x,0,0,0) = (ll > 0 || ml > 0);
    if (ml * ul < 0) 
      fine(v.x,0,1,0) = ml < 0 ? 1. - fu : fu ;
    else
      fine(v.x,0,1,0) = (ml > 0 || ul > 0);
  }
  
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors &&
      (is_local(cell) || is_local(neighbor(1)))) { //rhs
    double lr = GEOM (x + Delta/2., y - Delta/2);
    double mr = GEOM (x + Delta/2., y);
    double ur = GEOM (x + Delta/2., y + Delta/2);
    double fl = 0.5, fu = 0.5;
    if (lr != mr)
      fl = lr / (lr - mr);
    if (mr != ur)
      fu = mr / (mr - ur);
    if (lr * mr <= 0) 
      fine(v.x,2,0,0) = lr < 0 ? 1 - fl : fl ;
    else
      fine(v.x,2,0,0) = (lr > 0 || mr > 0);
    if (mr * ur <= 0) 
      fine(v.x,2,1,1) = mr < 0 ? 1. - fu : fu ;
    else
      fine(v.x,2,1,1) = (mr > 0 || ur > 0);
  }
  if (is_local(cell)) {
    double lm = GEOM (x, y - Delta/2);
    double mm = GEOM (x, y);
    double um = GEOM (x, y + Delta/2);
    double fl = 0.5, fu = 0.5;
    if (lm != mm)
      fl = lm / (lm - mm);
    if (mm != um)
      fu = mm / (mm - um);
    if (lm * mm < 0) 
      fine(v.x,1,0,0) = lm < 0 ? 1. - fl : fl ;
    else
      fine(v.x,1,0,0) = (lm > 0 || mm > 0);
    if (mm * um < 0) 
      fine(v.x,1,1,0) = mm < 0 ? 1. - fu : fu ;
    else
      fine(v.x,1,1,0) = (mm > 0 || um > 0);
    }
}

static inline void recompute_face_fraction_refine_2D_y (Point point, scalar s) {
  fflush (stdout);
  fflush (stdout);
  vector v = s.v;
  if (!is_refined(neighbor(0,-1)) &&
      (is_local(cell) || is_local(neighbor(0,-1)))) { //lower
    double ll = GEOM (x - Delta/2., y - Delta/2.);
    double ml = GEOM (x           , y - Delta/2.);
    double ul = GEOM (x + Delta/2., y - Delta/2.);
    double fl = 0.5, fu = 0.5;
    if (ll != ml)
      fl = ll / (ll - ml);
    if (ml != ul)
      fu = ml / (ml - ul);
    if (ll * ml < 0) 
      fine(v.y,0,0,0) = ll < 0 ? 1. - fl : fl ;
    else
      fine(v.y,0,0,0) = (ll > 0 || ml > 0);
    if (ml * ul < 0) 
      fine(v.y,1,0,0) = ml < 0 ? 1. - fu : fu ;
    else
      fine(v.y,1,0,0) = (ml > 0 || ul > 0);
  }
  if (!is_refined(neighbor(0,1)) && neighbor(0,1).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1)))) { //Upper
    double lr = GEOM (x - Delta/2., y + Delta/2.);
    double mr = GEOM (x           , y + Delta/2.);
    double ur = GEOM (x + Delta/2., y + Delta/2.);
    double fl = 0.5, fu = 0.5;
    if (lr != mr)
      fl = lr / (lr - mr);
    if (mr != ur)
      fu = mr / (mr - ur);
    if (lr * mr < 0) 
      fine(v.y,0,2,0) = lr < 0 ? 1. - fl : fl ;
    else
      fine(v.y,0,2,0) = (lr > 0 || mr > 0);
    if (mr * ur < 0) 
      fine(v.y,1,2,0) = mr < 0 ? 1. - fu : fu ;
    else
      fine(v.y,1,2,0) = (mr > 0 || ur > 0);
  }
  if (is_local(cell)) {
    double lm = GEOM (x - Delta/2., y);
    double mm = GEOM (x           , y);
    double um = GEOM (x + Delta/2., y);
    double fl = 0.5, fu = 0.5;
    if (lm != mm)
      fl = lm / (lm - mm);
    if (mm != um)
      fu = mm / (mm - um);
    if (lm * mm < 0) 
      fine(v.y,0,1,0) = lm < 0 ? 1. - fl : fl ;
    else
      fine(v.y,0,1,0) = (lm > 0 || mm > 0);
    if (mm * um < 0) 
      fine(v.y,1,1,0) = mm < 0 ? 1. - fu : fu ;
    else
      fine(v.y,1,1,0) = (mm > 0 || um > 0);
  }  
}

double intercept (coord n, Point point, face vector fs) {
  double alpha = 0;
  int ni = 0;
  if (fs.x[] > 0. && fs.x[] < 1) {
    alpha += sign(GEOM(x - Delta/2., y - Delta/2.))*n.y*(fs.x[] - 0.5) + n.x * -0.5;
    ni++;
  }
  if (fs.x[1] > 0. && fs.x[1] < 1.) {
    alpha += sign(GEOM(x + Delta/2., y - Delta/2.))*n.y*(fs.x[1] - 0.5) + n.x * 0.5;
    ni++;
  }
  if (fs.y[] > 0. && fs.y[] < 1.) {
    alpha += sign(GEOM(x - Delta/2., y - Delta/2.))*n.x*(fs.y[0] - 0.5) + n.y * -0.5;
    ni++;
  }
  if (fs.y[0,1] > 0. && fs.y[0,1] < 1.) {
    alpha += sign(GEOM(x - Delta/2., y + Delta/2.))*n.x*(fs.y[0,1] - 0.5) + n.y * 0.5;
    ni++;
  }
  if (ni != 0)
    return (alpha/ni);
  else 
    return (max (fs.x[], 0));
}

static void refine_embed_fraction_recompute (Point point, scalar c) {
  foreach_child() {
    if (fs.x[] == 0. && fs.x[1] == 0. && fs.y[] == 0. && fs.y[0,1] == 0.) {
	c[] = 0.;
    } else if (fs.x[] == 1. && fs.x[1] == 1. && fs.y[] == 1. && fs.y[0,1] == 1.) {
      c[] = 1.;
    } else {
      coord n = facet_normal (point, c, fs); 
      double alpha = intercept (n, point, fs);
      c[] = line_area (n.x, n.y, alpha);
    }
  }
}

#if EMBED
event set_cs_and_fs (i = 0) {
  /**
`fs` should appear before `cs` in lists.
   */
  scalar * temp = NULL;
  temp = list_append (temp, fs.x);
  temp = list_append (temp, fs.y);
  temp = list_append (temp, cs);
  for (scalar s in all) {
    if (s.i != cs.i && !(s.i >= fs.x.i && s.i <= fs.x.i + dimension - 1))
      temp = list_append (temp, s);
  }
  all = temp;
  /**
oef...
  */
  foreach_dimension() {
    fs.x.refine = recompute_face_fraction_refine_2D_x;
    fs.x.restriction = restriction_face; // :S
  }
  cs.refine = refine_embed_fraction_recompute;
}
#endif
/**
   `refine_embed_linear` is useful but it does not like fraction soups
   that emerge from degenerate representation of `GEOM`etries when
   they are not consistent between levels. Thus we may simply ignore
   the assertions and hope for good results.
 */

static inline void refine_embed_linear_not_assert (Point point, scalar s)
{
  foreach_child() {
    if (!cs[])
      s[] = 0.;
    else {
//assert (coarse(cs));
      int i = (child.x + 1)/2, j = (child.y + 1)/2;
#if dimension == 2
      if (coarse(fs.x,i) && coarse(fs.y,0,j) &&
	  (coarse(cs) == 1. || coarse(cs,child.x) == 1. ||
	   coarse(cs,0,child.y) == 1. || coarse(cs,child.x,child.y) == 1.)) {
//assert (coarse(cs,child.x) && coarse(cs,0,child.y));
	if (coarse(fs.x,i,child.y) && coarse(fs.y,child.x,j)) {
          // bilinear interpolation
//assert (coarse(cs,child.x,child.y));
	  s[] = (9.*coarse(s) + 
		 3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
		 coarse(s,child.x,child.y))/16.;
	}
	else
	  // triangular interpolation	  
	  s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      }
      else if (coarse(cs,child.x,child.y) &&
	       ((coarse(fs.x,i) && coarse(fs.y,child.x,j)) ||
		(coarse(fs.y,0,j) && coarse(fs.x,i,child.y)))) {
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
      }
#else // dimension == 3
      int k = (child.z + 1)/2;
      if (coarse(fs.x,i) > 0.25 && coarse(fs.y,0,j) > 0.25 &&
	  coarse(fs.z,0,0,k) > 0.25 &&
	  (coarse(cs) == 1. || coarse(cs,child.x) == 1. ||
	   coarse(cs,0,child.y) == 1. || coarse(cs,child.x,child.y) == 1. ||
	   coarse(cs,0,0,child.z) == 1. || coarse(cs,child.x,0,child.z) == 1. ||
	   coarse(cs,0,child.y,child.z) == 1. ||
	   coarse(cs,child.x,child.y,child.z) == 1.)) {
//assert (coarse(cs,child.x) && coarse(cs,0,child.y) &&
         coarse(cs,0,0,child.z));
      if (coarse(fs.x,i,child.y) && coarse(fs.y,child.x,j) &&
	  coarse(fs.z,child.x,child.y,k) &&
	  coarse(fs.z,child.x,0,k) && coarse(fs.z,0,child.y,k)) {
//assert (coarse(cs,child.x,child.y) && coarse(cs,child.x,0,child.z) &&
		  coarse(cs,0,child.y,child.z) &&
		  coarse(cs,child.x,child.y,child.z));
// bilinear interpolation
	  s[] = (27.*coarse(s) + 
		 9.*(coarse(s,child.x) + coarse(s,0,child.y) +
		     coarse(s,0,0,child.z)) + 
		 3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
		     coarse(s,0,child.y,child.z)) + 
		 coarse(s,child.x,child.y,child.z))/64.;
	}
	else
	  // tetrahedral interpolation
	  s[] = (coarse(s) + coarse(s,child.x) + coarse(s,0,child.y) +
		 coarse(s,0,0,child.z))/4.;
      }
      else if (coarse(cs,child.x,child.y,child.z) &&
	       ((coarse(fs.z,child.x,child.y,k) &&
		 ((coarse(fs.x,i) && coarse(fs.y,child.x,j)) ||
		  (coarse(fs.y,0,j) && coarse(fs.x,i,child.y))))
		||
		(coarse(fs.z,0,0,k) &&
		 ((coarse(fs.x,i,0,child.z) && coarse(fs.y,child.x,j,child.z)) ||
		  (coarse(fs.y,0,j,child.z) && coarse(fs.x,i,child.y,child.z))))
		||
		(coarse(fs.z,child.x,0,k) &&
		 coarse(fs.x,i) && coarse(fs.y,child.x,j,child.z))
		||
		(coarse(fs.z,0,child.y,k) &&
		 coarse(fs.y,0,j) && coarse(fs.x,i,child.y,child.z))
		))
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
	// Pathological cases, use 1D gradients.
	s[] = coarse(s);
	foreach_dimension() {
	  if (coarse(fs.x,(child.x + 1)/2) && coarse(cs,child.x))
	    s[] += (coarse(s,child.x) - coarse(s))/4.;
	  else if (coarse(fs.x,(- child.x + 1)/2) && coarse(cs,- child.x))
	    s[] -= (coarse(s,- child.x) - coarse(s))/4.;
	}
      }
    }
  }
}
# if EMBED

event properties (i = 0) {
  for (scalar s in all) {
    if (s.prolongation == refine_embed_linear)
      s.prolongation = refine_embed_linear_not_assert;
    if (s.refine == refine_embed_linear)
      s.refine = refine_embed_linear_not_assert;
  }
}
#endif //embed

#endif //tree
