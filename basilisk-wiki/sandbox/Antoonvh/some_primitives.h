/**
# Some primitives for bview

One may draw
* A cone
* A cylinder
* An arrow
* A curve

*/

coord coord_diff (coord a, coord b) {
  return (coord){a.x - b.x, a.y - b.y, a.z - b.z};
}

coord cross (coord a, coord b) {
  return (coord){a.y*b.z - b.y*a.z,
		 a.z*b.x - a.x*b.z,
		 a.x*b.y - a.y*b.x};
}

coord rejection (coord a, coord b) {
  double s = a.x*b.x + a.y*b.y + a.z*b.z;
  return (coord){b.x - s*a.x, b.y - s*a.y, b.z - s*a.z};
}

bool draw_arrow_head (coord dir = {1, 0, 0},
		 coord base = {0, 0, 0},
		 double len = 1,
		 double rad = 0.5,
		 float fc[3] = {.9, 0, 0.9},
		 int sections = 20,
		 bool draw_base = true){
  bview * view = draw();
  normalize (&dir);
  double dtheta = 2*pi/sections;
  glShadeModel (GL_SMOOTH);
  glColor3f(fc[0], fc[1], fc[2]);
  coord top = {base.x + len*dir.x , base.y + len*dir.y , base.z + len*dir.z};
  coord random = {noise(), noise(), noise()};
  coord n1 = rejection (dir, random);
  normalize (&n1);
  coord n2 = cross (dir, n1);
  // cone
  glBegin (GL_TRIANGLE_FAN);
  glnormal3d (view, 0, 0, 0); // Absence of a normal
  glvertex3d (view, top.x, top.y, top.z);
  for (double the = 0; the < 2*pi + 1e-6; the += dtheta) {
    double then = the + dtheta;
    coord b1 = {base.x + rad*(cos(the)*n1.x + sin(the)*n2.x),
		base.y + rad*(cos(the)*n1.y + sin(the)*n2.y),
		base.z + rad*(cos(the)*n1.z + sin(the)*n2.z)};
    coord b2 = {base.x + rad*(cos(then)*n1.x + sin(then)*n2.x),
		base.y + rad*(cos(then)*n1.y + sin(then)*n2.y),
		base.z + rad*(cos(then)*n1.z + sin(then)*n2.z)};
    coord n = cross (coord_diff (top, b1), coord_diff(top, b2));
    glnormal3d (view, n.x, n.y, n.z);
    glvertex3d (view, b1.x, b1.y, b1.z);
  }
  glEnd();
  if (draw_base) {
  // base circle
    glBegin (GL_TRIANGLE_FAN);
    glnormal3d (view, 0,0,0); 
    glvertex3d (view, base.x, base.y, base.z);
    for (double the = 0; the > -2*pi - 1e-6; the -= dtheta) {
      double then = the - dtheta;
      coord b1 = {base.x + rad*(cos(the)*n1.x + sin(the)*n2.x),
		  base.y + rad*(cos(the)*n1.y + sin(the)*n2.y),
		  base.z + rad*(cos(the)*n1.z + sin(the)*n2.z)};
      coord b2 = {base.x + rad*(cos(then)*n1.x + sin(then)*n2.x),
		  base.y + rad*(cos(then)*n1.y + sin(then)*n2.y),
		  base.z + rad*(cos(then)*n1.z + sin(then)*n2.z)};
      coord n = cross (coord_diff (base, b1), coord_diff(base, b2));
      glnormal3d (view, n.x, n.y, n.z);
      glvertex3d (view, b1.x, b1.y, b1.z);
    }
  }
  glEnd();
  return true;
}

bool draw_cylinder (coord dir = {1, 0, 0},
	       coord base = {0, 0, 0},
	       double len = 2,
	       double rad = 0.5,
	       float fc[3] = {.9, 0, 0.9},
	       int sections = 20,
	       bool base_cap = true,
	       bool end_cap = true) {
  bview * view = draw();
  normalize (&dir);
  double dtheta = 2*pi/sections;
  glShadeModel (GL_SMOOTH);
  glColor3f(fc[0], fc[1], fc[2]);
  coord random = {noise(), noise(), noise()};
  coord n1 = rejection (dir, random);
  normalize (&n1);
  coord n2 = cross (dir, n1);
  coord base2 = {base.x + len*dir.x, base.y + len*dir.y, base.z + len*dir.z};
  
  // Mantle
  for (double the = 0; the < 2*pi + 1e-6; the += dtheta) {
    double then = the + dtheta;
    glBegin (GL_POLYGON);
    coord b1 = {base.x + rad*(cos(the)*n1.x + sin(the)*n2.x),
		base.y + rad*(cos(the)*n1.y + sin(the)*n2.y),
		base.z + rad*(cos(the)*n1.z + sin(the)*n2.z)};
    coord b2 = {base.x + rad*(cos(then)*n1.x + sin(then)*n2.x),
		base.y + rad*(cos(then)*n1.y + sin(then)*n2.y),
		base.z + rad*(cos(then)*n1.z + sin(then)*n2.z)};
    coord b4 = {base2.x + rad*(cos(the)*n1.x + sin(the)*n2.x),
		base2.y + rad*(cos(the)*n1.y + sin(the)*n2.y),
		base2.z + rad*(cos(the)*n1.z + sin(the)*n2.z)};
    coord b3 = {base2.x + rad*(cos(then)*n1.x + sin(then)*n2.x),
		base2.y + rad*(cos(then)*n1.y + sin(then)*n2.y),
		base2.z + rad*(cos(then)*n1.z + sin(then)*n2.z)};
    coord n = cross (coord_diff (b1, b2), coord_diff(b1, b4));
    coord coords[4] = {b1, b2, b3, b4};
    for (int i = 0; i < 4; i++) {
      glnormal3d (view, n.x, n.y, n.z);
      glvertex3d (view, coords[i].x, coords[i].y, coords[i].z);
    }
    view->ni++;
    glEnd();
  }
  // caps
  if (base_cap || end_cap) {
    for (double a = len*(!base_cap); a <= len*end_cap; a += len) {
      coord basen = (coord){base.x + a*dir.x, base.y + a*dir.y, base.z + a*dir.z};
      glBegin (GL_TRIANGLE_FAN);
      glnormal3d (view, 0, 0, 0); 
      glvertex3d (view, base.x, base.y, base.z);
      for (double the = 0; the > -2*pi - 1e-6; the -= dtheta) {
	double then = the - dtheta;
	coord b1 = {basen.x + rad*(cos(the)*n1.x + sin(the)*n2.x),
		    basen.y + rad*(cos(the)*n1.y + sin(the)*n2.y),
		    basen.z + rad*(cos(the)*n1.z + sin(the)*n2.z)};
	coord b2 = {basen.x + rad*(cos(then)*n1.x + sin(then)*n2.x),
		    basen.y + rad*(cos(then)*n1.y + sin(then)*n2.y),
		    basen.z + rad*(cos(then)*n1.z + sin(then)*n2.z)};
	coord n = cross (coord_diff (basen, b1), coord_diff(basen, b2));
	glnormal3d (view, n.x, n.y, n.z);
	glvertex3d (view, b1.x, b1.y, b1.z);
      }
      glEnd();
    }
  }
  return true;
}

bool draw_arrow (coord dir = {1, 0, 0},
		 coord base = {0, 0, 0},
		 double len = 2,
		 double rad = 0.5,
		 double head_rad_f = 2, //twice the radius of the cylinder
		 double head_len_rad = 2, // one dimater "length"
		 float fc[3] = {.9, 0, .9},
		 int sections = 20) {
  normalize (&dir);
  draw_cylinder (dir, base, len, rad, fc = fc, sections = sections);
  coord arrow_base = {base.x + len*dir.x, base.y + len*dir.y, base.z + len*dir.z};
  double arrow_rad = rad*head_rad_f;
  double arrow_len = arrow_rad*head_len_rad;
  draw_arrow_head (dir, arrow_base, arrow_len, arrow_rad, fc, sections);
  return true;
}

bool draw_axis_cross (coord base = {0, 0, 0},
		double len = 1,
		double rad = 0.1,
		int sections = 20) {
  int dim = 0;
  foreach_dimension() {
    coord dir = {0, 0, 0};
    coord new_base = base;
    float col[3] = {.2, .2, .2};
    if (dim == 0) 
      col[0] = .8;
    else if (dim == 1)
      col[0] = 0.8, col[1] = .8;
    else if (dim == 2)
      col[1] = .8;
    dim++;
    dir.x = 1;
    new_base.x -= rad;
    draw_arrow (dir, new_base, len, rad, fc = col);
  }
  return true;
}

bool draw_paramterization (double rad, double ts, double te,  coord (*my_fun)(double), int sections, float fc[3]) {
  double dta = (te - ts) / sections;
  double overlap = dta/10.;
  for (int i = 0; i < sections; i++) {
    coord base = my_fun(ts + dta*i - overlap);
    coord end = my_fun(ts + dta*(i + 1) + overlap);
    coord dir = coord_diff (end, base);
    double len = sqrt (sq(dir.x) + sq(dir.y) + sq(dir.z));
    if (i == 0)
      draw_cylinder (dir, base, len, rad, fc = fc, sections = 25, end_cap = false);
    else if (i == (sections - 1))
      draw_cylinder (dir, base, len, rad, fc = fc, sections = 25, base_cap = false);
    else
      draw_cylinder (dir, base, len, rad, fc = fc, sections = 25, base_cap = false, end_cap = false);
    // Undo translation
    glTranslatef (0, 0, 1e-4); 
  }
  return true;
}

/**
   
## Test
* [A test](test_primitives.c)
*/
