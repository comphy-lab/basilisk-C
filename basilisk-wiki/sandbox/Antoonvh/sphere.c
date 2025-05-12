/**
# Mapping a single Quadtree to a sphere in bview. 

It requires cheating....

![Mapping in view](sphere/circle.mp4)


*/

#include "view.h"

void map_circle (coord * p) {
  double x = p->x, y = p->y;
  p->x = x*sqrt(1 - sq(y)/2.);
  p->y = y*sqrt(1 - sq(x)/2.);
}

struct _cells2 {
  char * z;
  float lc[3], lw; // the line color and width
};
bool cells2 (struct _cells2 p);

int main() {
  L0 = 2 - 1e-5;
  X0 = Y0 = -L0/2;
  init_grid (N);
  scalar zp[], zp2[], m[];
  foreach_dimension() {
    for (scalar s in {zp, zp2}) {
      s[left] =  dirichlet(0);
      s[right]= dirichlet(0);
    }
  }
  foreach() {
    m[] = x;
    double xp = x*sqrt(1 - sq(y)/2.);
    double yp = y*sqrt(1 - sq(x)/2.);
    zp[] = sqrt(1. - sq(xp) - sq(yp));
    zp2[] = -zp[];
  }
  for (double thet = 0; thet <= 2*pi; thet += 0.05*pi) {
    view (map = map_circle, theta = thet, psi = sin(thet));
    squares("m", z = "zp", linear = true);
    squares("m", z = "zp2", linear = true);
    cells2("zp");
    cells2("zp2");
    save ("circle.mp4");
    clear();
  }
}



trace
bool cells2 (struct _cells2 p)
{
  bview * view = draw();
  draw_lines (view, p.lc, p.lw) {
    scalar Z = {-1};
    bool zexpr = false;
    Z = compile_expression (p.z, &zexpr);
    if (Z.i < 0)
      return false;
    foreach_visible (view) {
      glBegin (GL_LINE_LOOP);
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
  return true;
}
