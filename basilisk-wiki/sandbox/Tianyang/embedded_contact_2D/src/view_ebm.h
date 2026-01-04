#include "view.h"

trace
bool draw_vof_ebm (char * cs, char * c, char * cet, char * s = NULL, bool edges = false,
	       double larger = 0., int filled = 0,
	       char * color = NULL,
	       double min = 0, double max = 0, double spread = 0,
	       bool linear = false,
	       Colormap map = jet,
	       float fc[3] = {0}, float lc[3] = {0}, float lw = 1.,
	       bool expr = false,
	       COLORBAR_PARAMS)
{
  scalar d = lookup_field (c);
  scalar dcs = lookup_field (cs);
  scalar dcet = lookup_field (cet);
  if (d.i < 0) {
    fprintf (stdout, "draw_vof(): no field named '%s'\n", c);
    return false;
  }
  face vector fs = lookup_vector (s);
  
  colorize_args();
  
  double cmin = 1e-3; // do not reconstruct fragments smaller than this

#if TREE
  // make sure we prolongate properly
  void (* prolongation) (Point, scalar) = d.prolongation;
  if (prolongation != fraction_refine) {
    d.prolongation = fraction_refine;
    d.dirty = true;
  }
#endif // TREE
    
  bview * view = draw();
#if dimension == 2
  if (filled) {
    glColor3f (fc[0], fc[1], fc[2]);
    glNormal3d (0, 0, view->reversed ? -1 : 1);
    foreach_visible (view) {
      if ((filled > 0 && d[] >= 1.) || (filled < 0 && d[] <= 0.)) {
	glBegin (GL_QUADS);
	glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
	glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
	glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
	glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
	glEnd();
	view->ni++;
      }
      else if (d[] > 0. && d[] < 1.) {
	coord n = facet_normal (point, d, fs), r = {1.,1.};
	if (filled < 0)
	  foreach_dimension()
	    n.x = - n.x;
	double alpha = plane_alpha (filled < 0. ? 1. - d[] : d[], n);
	alpha += (n.x + n.y)/2.;
	foreach_dimension()
	  if (n.x < 0.) alpha -= n.x, n.x = - n.x, r.x = - 1.;
	coord v[5];
	int nv = 0;
	if (alpha >= 0. && alpha <= n.x) {
	  v[nv].x = alpha/n.x, v[nv++].y = 0.;
	  if (alpha <= n.y)
	    v[nv].x = 0., v[nv++].y = alpha/n.y;
	  else if (alpha >= n.y && alpha - n.y <= n.x) {
	    v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;
	    v[nv].x = 0., v[nv++].y = 1.;
	  }
	  v[nv].x = 0., v[nv++].y = 0.;
	}
	else if (alpha >= n.x && alpha - n.x <= n.y) {
	  v[nv].x = 1., v[nv++].y = (alpha - n.x)/n.y;
	  if (alpha >= n.y && alpha - n.y <= n.x) {
	    v[nv].x = (alpha - n.y)/n.x, v[nv++].y = 1.;
	    v[nv].x = 0., v[nv++].y = 1.;
	  }
	  else if (alpha <= n.y)
	    v[nv].x = 0., v[nv++].y = alpha/n.y;
	  v[nv].x = 0., v[nv++].y = 0.;
	  v[nv].x = 1., v[nv++].y = 0.;
	}
	glBegin (GL_POLYGON);
	if (r.x*r.y < 0.)
	  for (int i = nv - 1; i >= 0; i--)
	    glvertex2d (view, x + r.x*(v[i].x - 0.5)*Delta,
			y + r.y*(v[i].y - 0.5)*Delta);
	else
	  for (int i = 0; i < nv; i++)
	    glvertex2d (view, x + r.x*(v[i].x - 0.5)*Delta,
			y + r.y*(v[i].y - 0.5)*Delta);
	glEnd ();
	view->ni++;
      }
    }
  }
  else // !filled
    draw_lines (view, lc, lw) {
      glBegin (GL_LINES);
      foreach_visible (view)
        if (d[] > 0. && d[] < dcs[] && cfilter (point, dcet, cmin)) {
          coord n = interface_normal (point, d);
          double alpha = plane_alpha (dcet[], n);
          coord segment[2];
          if (facets_ebm (n, alpha, segment) == 2) {
            glvertex2d (view, x + segment[0].x*Delta, y + segment[0].y*Delta);
            glvertex2d (view, x + segment[1].x*Delta, y + segment[1].y*Delta);
            view->ni++;
          }
        }
      glEnd ();
    }
#else // dimension == 3
  if (!larger)
    larger = edges || (color && !linear) ? 1. : 1.1;
  if (edges)
    draw_lines (view, lc, lw) {
      foreach_visible (view)
	if (cfilter (point, d, cmin)) {
	  coord n = facet_normal (point, d, fs);
	  double alpha = plane_alpha (d[], n);
	  coord v[12];
	  int m = facets_ebm (n, alpha, v, larger);
	  if (m > 2) {
	    glBegin (GL_LINE_LOOP);
	    for (int i = 0; i < m; i++)
	      glvertex3d (view,
			  x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
	    glEnd ();
	    view->ni++;
	  }
	}
    }
  else // !edges
    colorize() {
      foreach_visible (view)
	if (cfilter (point, d, cmin)) {
	  coord n = facet_normal (point, d, fs);
	  double alpha = plane_alpha (d[], n);
	  coord v[12];
	  int m = facets_ebm (n, alpha, v, larger);
	  if (m > 2) {
	    glBegin (GL_POLYGON);
	    for (int i = 0; i < m; i++) {
	      if (linear) {
		color_vertex (interp (point, v[i], col));
	      }
	      else {
		color_facet();
	      }
	      glnormal3d (view, n.x, n.y, n.z);
	      glvertex3d (view,
			  x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
	    }
	    glEnd ();
	    view->ni++;
	  }
	}
    }
#endif // dimension == 3

#if TREE
  // revert prolongation
  if (prolongation != fraction_refine) {
    d.prolongation = prolongation;
    d.dirty = true;
  }
#endif // TREE

  if (expr) delete({col});
  return true;
}