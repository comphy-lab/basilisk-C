/**
![David Hilbert on displayed on one of his favortice Curves [(by
 Stansy via
 Wikimedia)](https://commons.wikimedia.org/wiki/File:The_seventh_order_of_a_Hilbert_curve.jpg)](https://upload.wikimedia.org/wikipedia/commons/5/59/The_seventh_order_of_a_Hilbert_curve.jpg){
 width="300" }

# Hilbert curve grid iterator

![It fills the space](hilbert/hilbert_iterator.mp4)
 */
#include "utils.h"

void Hilbert (Point point, int lp, Cache * c) {
  if (point.level == depth()) {
    cache_append (c, point, 0);
  } else {
    Point p1 = point;
    p1.level++;
    switch (lp) {
    case 1:
      p1.i = _LEFT;  p1.j = _BOTTOM; Hilbert (p1,4,c);
      p1.i = _LEFT;  p1.j = _TOP;    Hilbert (p1,1,c);
      p1.i = _RIGHT; p1.j = _TOP;    Hilbert (p1,1,c);
      p1.i = _RIGHT; p1.j = _BOTTOM; Hilbert (p1,2,c);
      break;
    case 2:
      p1.i = _RIGHT; p1.j = _TOP;    Hilbert (p1,3,c);
      p1.i = _LEFT;  p1.j = _TOP;    Hilbert (p1,2,c);
      p1.i = _LEFT;  p1.j = _BOTTOM; Hilbert (p1,2,c);
      p1.i = _RIGHT; p1.j = _BOTTOM; Hilbert (p1,1,c);
      break;
    case 3:
      p1.i = _RIGHT; p1.j = _TOP;    Hilbert (p1,2,c);
      p1.i = _RIGHT; p1.j = _BOTTOM; Hilbert (p1,3,c);
      p1.i = _LEFT;  p1.j = _BOTTOM; Hilbert (p1,3,c);
      p1.i = _LEFT;  p1.j = _TOP;    Hilbert (p1,4,c);
      break;
    case 4:
      p1.i = _LEFT;   p1.j = _BOTTOM; Hilbert (p1,1,c);
      p1.i = _RIGHT;  p1.j = _BOTTOM; Hilbert (p1,4,c);
      p1.i = _RIGHT;  p1.j = _TOP;    Hilbert (p1,4,c);
      p1.i = _LEFT;   p1.j = _TOP;    Hilbert (p1,3,c);
      break;
    }
  }
}

int main() {
  init_grid (32);
  Cache H = {0};
  foreach_level(0) 
    Hilbert (point, 1, &H);
  cache_shrink (&H);
  scalar s[];
  foreach()
    s[] = nodata;
  double a = 0;
  foreach_cache(H) {
    s[] = a++;
    output_ppm (s, file = "hilbert_iterator.mp4", n = 500);
  }
  free(H.p);
}
