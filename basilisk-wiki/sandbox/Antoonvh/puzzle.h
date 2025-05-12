/**
# Display the grid as a puzzle

You can set some parameters for the lock and use it in 2D
 */

double puzzle_base = 0.12, pb_height = -0.03, pb_radius = 0.12, side_radius = 2.0;

/**
A random number for each cell
 */

union double_int {
  double d;
  int i;
};

int randomizer (int i1, int i2, int i3) {
  union double_int d1;
  d1.d =(double)i1*i2/757. + i2*i3/787. + i1*i3/797.;
  return d1.i;
}
 /**
 Drawing:
 */

void puzzle_side_x (Point point, coord c, bview * view) {
  int nrs = 20;
  int rando = randomizer (point.i,point.j, level)% 2 == 0 ? 1 : -1 ;
  
  double theta0 = asin (puzzle_base/(2.*pb_radius));
  double theta1 = asin (0.5/side_radius);
  double theta2 = asin (puzzle_base/(2.*side_radius));
  glBegin (GL_LINE_STRIP);
  glvertex2d (view, c.x, c.y - Delta_y/2.);
  if (fabs(c.x - X0) > 1e-6 && fabs(c.x - X0 - L0) > 1e-6) {
    for (double theta = theta1; theta >= theta2 ; theta -= pi/(8.*nrs)) {
      glvertex2d (view, c.x + rando*Delta_x*(side_radius*cos(theta1) - side_radius*cos(theta)),
		  c.y - Delta_y*side_radius*sin(theta));
    }
    for (double theta = theta0; theta <= 2*pi-theta0; theta += pi/nrs) {
      glvertex2d (view, c.x + rando*Delta_x*(pb_height + pb_radius - pb_radius*cos(theta)),
		  c.y - Delta_y*pb_radius*sin(theta));
    }
    for (double theta = theta2; theta <= theta1 ; theta += pi/(8.*nrs)) {
      glvertex2d (view, c.x + rando*Delta_x*(side_radius*cos(theta1) - side_radius*cos(theta)),
		  c.y + Delta_y*side_radius*sin(theta));
    }
    
    //glvertex2d (view, c.x, c.y + Delta_y*puzzle_base/2);
  }
  glvertex2d (view, c.x, c.y + Delta_y/2.);
  glEnd();
}

void puzzle_side_y (Point point, coord c, bview * view) {
  int nrs = 20;
  int rando = randomizer(point.i, level, point.j)% 2 == 0 ? 1 : -1 ;
  double theta0 = asin (puzzle_base/(2.*pb_radius));
  double theta1 = asin (0.5/side_radius);
  double theta2 = asin (puzzle_base/(2.*side_radius));
  glBegin (GL_LINE_STRIP);
  glvertex2d (view, c.x - Delta_x/2., c.y);
  if (fabs(c.y - Y0) > 1e-6 && fabs(c.y - Y0 - L0) > 1e-6) {
    for (double theta = theta1; theta >= theta2 ; theta -= pi/(8.*nrs)) {
      glvertex2d (view, c.x - Delta_x*side_radius*sin(theta),
		  c.y + rando*Delta_y*(side_radius*cos(theta1) - side_radius*cos(theta)));
    }
    for (double theta = theta0; theta <= 2*pi-theta0; theta += pi/nrs) {
      glvertex2d (view, c.x - Delta_y*pb_radius*sin(theta),
		  c.y + rando*Delta_y*(pb_height + pb_radius - pb_radius*cos(theta)) );
    }
    for (double theta = theta2; theta <= theta1 ; theta += pi/(8.*nrs)) {
      glvertex2d (view, c.x + Delta_x*side_radius*sin(theta),
		  c.y + rando*Delta_y*(side_radius*cos(theta1) - side_radius*cos(theta)));
    }
  }
    glvertex2d (view, c.x  + Delta_x/2., c.y);
  glEnd();
}

trace
bool puzzle_cells (struct _cells p)
{
  bview * view = draw();
  draw_lines (view, p.lc, p.lw) {
#if dimension == 2
    foreach_face (x) {
      coord c = {x, y};
      puzzle_side_x (point, c, view);
      view->ni++;
    }
    foreach_face (y) {
      coord c = {x,y};
      puzzle_side_y (point, c, view);
      view->ni++;
    }
#else // dimension == 3
    assert (dimension == 3);
#endif // dimension == 3
  }
  return true;
}
/**
## Test
* [A simple test](test_puzzle.c)
*/
