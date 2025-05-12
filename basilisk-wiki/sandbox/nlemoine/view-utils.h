#include "view.h"

//bool surf_cells ((const) scalar s, float lc[3] = {0}, float lw = 1.)
bool surf_cells ((const) scalar Z, scalar mask = {-1}, float lc[3] = {0}, float lw = 1.)
{
  bview * view = draw();
  draw_lines (view, lc, lw) {
//    foreach_visible (view) {
      foreach_leaf(){
      bool display_cell = (mask.i < 0) || (mask[]>0.);
      if(display_cell)
      { 
        glBegin (GL_LINE_LOOP);
/*        glvertex3d (view, x - Delta_x/2., y - Delta_y/2., interpolate (s, x - Delta_x/2., y - Delta_y/2.));
        glvertex3d (view, x + Delta_x/2., y - Delta_y/2., interpolate (s, x + Delta_x/2., y - Delta_y/2.));
        glvertex3d (view, x + Delta_x/2., y + Delta_y/2., interpolate (s, x + Delta_x/2., y + Delta_y/2.));
        glvertex3d (view, x - Delta_x/2., y + Delta_y/2., interpolate (s, x - Delta_x/2., y + Delta_y/2.));
*/
	glvertex3d (view, x - Delta_x/2., y - Delta_y/2.,
			   (Z[] + Z[-1] + Z[-1,-1] + Z[0,-1])/4.);
	glvertex3d (view, x + Delta_x/2., y - Delta_y/2.,
			   (Z[] + Z[1] + Z[1,-1] + Z[0,-1])/4.);
	glvertex3d (view, x + Delta_x/2., y + Delta_y/2.,
			   (Z[] + Z[1] + Z[1,1] + Z[0,1])/4.);
	glvertex3d (view, x - Delta_x/2., y + Delta_y/2.,
			   (Z[] + Z[-1] + Z[-1,1] + Z[0,1])/4.);
        glEnd();
        view->ni++;
      }
    }
  }
  return true;
}

trace
bool masked_squares (char * color,
	      char * z = NULL,
	      double min = 0, double max = 0, double spread = 0,
	      bool linear = false,
	      Colormap map = jet,
	      float fc[3] = {0}, float lc[3] = {0},
	      bool expr = false,
	      scalar mask = {-1},
	      coord n = {0,0,1},
	      double alpha = 0,
	      float lw = 1,
	      COLORBAR_PARAMS)
{
#if dimension == 2
  scalar Z = {-1};
  vector fn;
  bool zexpr = false;
  if (z) {
    Z = compile_expression (z, &zexpr);
    if (Z.i < 0)
      return false;
    fn = new vector;
    foreach()
      foreach_dimension()
        fn.x[] = (Z[1] - Z[-1])/(2.*Delta_x);
    boundary ({fn});
  }
#endif
  colorize_args();
  scalar f = col;
  
  bview * view = draw();
  glShadeModel (GL_SMOOTH);
  if (linear) {
    colorize() {
#if dimension == 2
      if (Z.i < 0) {
	glNormal3d (0, 0, view->reversed ? -1 : 1);
	foreach_visible (view)
        {
          bool display_cell = (f[]!= nodata) && ((mask.i < 0) || (mask[]>0.)) ;
	  if (display_cell) {
	    glBegin (GL_TRIANGLE_FAN);
	    color_vertex ((4.*f[] +
			   2.*(f[1] + f[-1] + f[0,1] + f[0,-1]) +
			   f[-1,-1] + f[1,1] + f[-1,1] + f[1,-1])/16.);
	    glvertex2d (view, x, y);
	    color_vertex ((f[] + f[-1] + f[-1,-1] + f[0,-1])/4.);
	    glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
	    color_vertex ((f[] + f[1] + f[1,-1] + f[0,-1])/4.);
	    glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
	    color_vertex ((f[] + f[1] + f[1,1] + f[0,1])/4.);
	    glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
	    color_vertex ((f[] + f[-1] + f[-1,1] + f[0,1])/4.);
	    glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
	    color_vertex ((f[] + f[-1] + f[-1,-1] + f[0,-1])/4.);
	    glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
	    glEnd();
	    view->ni++;
	  }
        }
      }
      else // Z.i > 0
	foreach_leaf() // fixme: foreach_visible() would be better
        {
          bool display_cell = (f[]!= nodata) && ((mask.i < 0) || (mask[]>0.)) ;
	  if (display_cell) {
	    glBegin (GL_TRIANGLE_FAN);
	    color_vertex ((4.*f[] +
			   2.*(f[1] + f[-1] + f[0,1] + f[0,-1]) +
			   f[-1,-1] + f[1,1] + f[-1,1] + f[1,-1])/16.);
	    glvertex3d (view, x, y, Z[]);
	    color_vertex ((f[] + f[-1] + f[-1,-1] + f[0,-1])/4.);
	    glvertex3d (view, x - Delta_x/2., y - Delta_y/2.,
			       (Z[] + Z[-1] + Z[-1,-1] + Z[0,-1])/4.);
	    color_vertex ((f[] + f[1] + f[1,-1] + f[0,-1])/4.);
	    glvertex3d (view, x + Delta_x/2., y - Delta_y/2.,
			       (Z[] + Z[1] + Z[1,-1] + Z[0,-1])/4.);
	    color_vertex ((f[] + f[1] + f[1,1] + f[0,1])/4.);
	    glvertex3d (view, x + Delta_x/2., y + Delta_y/2.,
			       (Z[] + Z[1] + Z[1,1] + Z[0,1])/4.);
	    color_vertex ((f[] + f[-1] + f[-1,1] + f[0,1])/4.);
	    glvertex3d (view, x - Delta_x/2., y + Delta_y/2.,
			       (Z[] + Z[-1] + Z[-1,1] + Z[0,1])/4.);
	    color_vertex ((f[] + f[-1] + f[-1,-1] + f[0,-1])/4.);
	    glvertex3d (view, x - Delta_x/2., y - Delta_y/2.,
			       (Z[] + Z[-1] + Z[-1,-1] + Z[0,-1])/4.);
	    glEnd();
	    view->ni++;	    
	  }
        }
#else // dimension == 3
      foreach_visible_plane (view, n, alpha)
	if (f[] != nodata) {
	  coord v[12];
	  int m = facets (n, alpha, v, 1.);
	  if (m > 2) {
	    coord c = {0,0,0};
	    for (int i = 0; i < m; i++)
	      foreach_dimension()
		c.x += v[i].x/m;
	    glBegin (GL_TRIANGLE_FAN);
	    color_vertex (interp (point, c, f));
	    glvertex3d (view, x + c.x*Delta, y + c.y*Delta, z + c.z*Delta);
	    for (int i = 0; i < m; i++) {
	      color_vertex (interp (point, v[i], f));
	      glvertex3d (view,
			  x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
	    }
	    color_vertex (interp (point, v[0], f));
	    glvertex3d (view,
			x + v[0].x*Delta, y + v[0].y*Delta, z + v[0].z*Delta);
	    glEnd ();
	    view->ni++;
	  }
	}
#endif // dimension == 3
    }
  }
  else { // !linear
#if dimension == 2
    glNormal3d (0, 0, view->reversed ? -1 : 1);
    glBegin (GL_QUADS);
    foreach_visible (view)
      if (f[] != nodata) {
	color_facet();
	glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
	color_facet();
	glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
	color_facet();
	glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
	color_facet();
	glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
	view->ni++;
      }
    glEnd();
#else // dimension == 3
    foreach_visible_plane (view, n, alpha)
      if (f[] != nodata) {
	coord v[12];
	int m = facets (n, alpha, v, 1.);
	if (m > 2) {
	  glBegin (GL_POLYGON);
	  for (int i = 0; i < m; i++) {
	    color_facet();
	    glvertex3d (view,
			x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
	  }
	  glEnd ();
	  view->ni++;
	}
      }
#endif // dimension == 3
  }
  if (expr) delete ({col});
#if dimension == 2
  if (zexpr) delete ({Z});
  if (z) delete ((scalar *){fn});
#endif
  return true;
}

bool masked_squares_checkerboard (char * color,
	      char * z = NULL,
	      double min = 0, double max = 0, double spread = 0,
	      bool linear = false,
	      Colormap map = jet,
	      float fc[3] = {0}, float lc[3] = {0},
	      bool expr = false,
	      scalar mask = {-1},
	      coord n = {0,0,1},
	      double alpha = 0,
	      float lw = 1,
	      COLORBAR_PARAMS)
{
#if dimension == 2
  scalar Z = {-1};
  vector fn;
  bool zexpr = false;
  if (z) {
    Z = compile_expression (z, &zexpr);
    if (Z.i < 0)
      return false;
    fn = new vector;
    foreach()
      foreach_dimension()
        fn.x[] = (Z[1] - Z[-1])/(2.*Delta_x);
    boundary ({fn});
  }
#endif
  colorize_args();
  scalar f = col;
  
  bview * view = draw();
  glShadeModel (GL_SMOOTH);
  if (linear) {
    colorize() {
#if dimension == 2
      if (Z.i < 0) {
	glNormal3d (0, 0, view->reversed ? -1 : 1);
	foreach_visible (view)
        {
          bool display_cell = (f[]!= nodata) && ((mask.i < 0) || (mask[]>0.)) ;
	  if (display_cell) {
	    glBegin (GL_TRIANGLE_FAN);
	    color_vertex ((4.*f[] +
			   2.*(f[1] + f[-1] + f[0,1] + f[0,-1]) +
			   f[-1,-1] + f[1,1] + f[-1,1] + f[1,-1])/16.);
	    glvertex2d (view, x, y);
	    color_vertex ((f[] + f[-1] + f[-1,-1] + f[0,-1])/4.);
	    glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
	    color_vertex ((f[] + f[1] + f[1,-1] + f[0,-1])/4.);
	    glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
	    color_vertex ((f[] + f[1] + f[1,1] + f[0,1])/4.);
	    glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
	    color_vertex ((f[] + f[-1] + f[-1,1] + f[0,1])/4.);
	    glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
	    color_vertex ((f[] + f[-1] + f[-1,-1] + f[0,-1])/4.);
	    glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
	    glEnd();
	    view->ni++;
	  }
        }
      }
      else // Z.i > 0
	foreach_leaf() // fixme: foreach_visible() would be better
        {
          bool display_cell = (f[]!= nodata) && ((mask.i < 0) || (mask[]>0.)) ;
	  if (display_cell) {
	    glBegin (GL_TRIANGLE_FAN);
            color_vertex (f[]);
	    glvertex3d (view, x, y, Z[]);
	    color_vertex (f[]);
	    glvertex3d (view, x - Delta_x/2., y - Delta_y/2.,
			       (Z[] + Z[-1] + Z[-1,-1] + Z[0,-1])/4.);
	    color_vertex (f[]);
	    glvertex3d (view, x + Delta_x/2., y - Delta_y/2.,
			       (Z[] + Z[1] + Z[1,-1] + Z[0,-1])/4.);
	    color_vertex (f[]);
	    glvertex3d (view, x + Delta_x/2., y + Delta_y/2.,
			       (Z[] + Z[1] + Z[1,1] + Z[0,1])/4.);
	    color_vertex (f[]);
	    glvertex3d (view, x - Delta_x/2., y + Delta_y/2.,
			       (Z[] + Z[-1] + Z[-1,1] + Z[0,1])/4.);
	    color_vertex (f[]);
	    glvertex3d (view, x - Delta_x/2., y - Delta_y/2.,
			       (Z[] + Z[-1] + Z[-1,-1] + Z[0,-1])/4.);
	    glEnd();
	    view->ni++;	    
	  }
        }
#else // dimension == 3
      foreach_visible_plane (view, n, alpha)
	if (f[] != nodata) {
	  coord v[12];
	  int m = facets (n, alpha, v, 1.);
	  if (m > 2) {
	    coord c = {0,0,0};
	    for (int i = 0; i < m; i++)
	      foreach_dimension()
		c.x += v[i].x/m;
	    glBegin (GL_TRIANGLE_FAN);
	    color_vertex (interp (point, c, f));
	    glvertex3d (view, x + c.x*Delta, y + c.y*Delta, z + c.z*Delta);
	    for (int i = 0; i < m; i++) {
	      color_vertex (interp (point, v[i], f));
	      glvertex3d (view,
			  x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
	    }
	    color_vertex (interp (point, v[0], f));
	    glvertex3d (view,
			x + v[0].x*Delta, y + v[0].y*Delta, z + v[0].z*Delta);
	    glEnd ();
	    view->ni++;
	  }
	}
#endif // dimension == 3
    }
  }
  else { // !linear
#if dimension == 2
    glNormal3d (0, 0, view->reversed ? -1 : 1);
    glBegin (GL_QUADS);
    foreach_visible (view)
      if (f[] != nodata) {
	color_facet();
	glvertex2d (view, x - Delta_x/2., y - Delta_y/2.);
	color_facet();
	glvertex2d (view, x + Delta_x/2., y - Delta_y/2.);
	color_facet();
	glvertex2d (view, x + Delta_x/2., y + Delta_y/2.);
	color_facet();
	glvertex2d (view, x - Delta_x/2., y + Delta_y/2.);
	view->ni++;
      }
    glEnd();
#else // dimension == 3
    foreach_visible_plane (view, n, alpha)
      if (f[] != nodata) {
	coord v[12];
	int m = facets (n, alpha, v, 1.);
	if (m > 2) {
	  glBegin (GL_POLYGON);
	  for (int i = 0; i < m; i++) {
	    color_facet();
	    glvertex3d (view,
			x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
	  }
	  glEnd ();
	  view->ni++;
	}
      }
#endif // dimension == 3
  }
  if (expr) delete ({col});
#if dimension == 2
  if (zexpr) delete ({Z});
  if (z) delete ((scalar *){fn});
#endif
  return true;
}