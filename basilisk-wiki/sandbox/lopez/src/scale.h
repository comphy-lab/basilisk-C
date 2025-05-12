/**
# A Segment scale bview function

The `scale` is a very simple bview function that allows you to add the
classical segment scale. The origin of the segment is determined with
`pos` and its end position with the coordinate variable `length`. As
usual the width and the color of the segment can also be set in a
similar way to [draw_vof()](/src/draw.h#draw_vof).

*/

void scale (coord pos = {-0.95, -0.95, 0}, coord length = {0.05, 0., 0.},
	   float lw = 1, float lc[3] = {0}){
  bview * view = draw();
  draw_lines (view, lc, lw) {
    foreach_visible (view) {
      glBegin (GL_LINES);
      glvertex2d (view, pos.x, pos.y);
      glvertex2d (view, pos.x + length.x, pos.y + length.y);
      glEnd();
      view->ni++;
    }
  }
}
/**
##test

*[A simple test / example](../scale.c)
*/