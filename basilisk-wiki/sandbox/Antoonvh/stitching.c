/**
# Stitching domains; a poor man's method

![A snake is left to roam outside its box!](stitching/ss.mp4)

time to cube a spehere?
 */
#include "utils.h"
#include "run.h"
#include "view.h"
scalar sl[], sr[], sb[];

int dirx = -1, diry = 0;

int main() {
  init_grid (32);

  sl[left] = val(sr,0,0);
  sr[left] = val(sl,0,0);
  
  sr[bottom] = val(sb,0,0);
  sb[bottom] = val(sr,0,0);
  
  run();
}

event init (t = 0) {
  foreach() {
    for (scalar s in {sl, sr, sb})
      s[] = nodata;
    if (fabs(x - 0.4) < Delta/2 &&
	fabs(y - 0.3) < Delta/2) {
      sr[] = 1;
    }
  }
}


event advance (i++) {
  foreach() {
    for (scalar s in {sl, sr, sb})
      if (s[] != nodata)
	s[]++;
  }
  foreach() {
    for (scalar s in {sl, sr, sb})
      if (s[] > 15)
	s[] = nodata;
  }
  foreach() {
    // Directions change between the fields
    if (sr[-dirx, -diry] == 2.)
      sr[] = 1;
    if (sl[dirx, -diry] == 2.)
      sl[] = 1;
    if (sb[-dirx, diry] == 2.)
      sb[] = 1;
  }
}

event go_up (i = 35) {
  dirx = 0;
  diry = 1;
}


event go_right (i = 50) {
  dirx = 1;
  diry = 0;
}


event go_down (i = 100) {
  dirx = 0;
  diry = -1;
}


event go_left (i = 140) {
  dirx = -1;
  diry = 0;
}


event go_up_again (i = 150) {
  dirx = 0;
  diry = 1;
}

event mov (i++) {
  view (fov = 40, width = 750, height = 750);
  cells();
  squares ("sr", min = 1, max = 2);
  mirror({-1}) {
    cells();
    squares ("sl", min = 1, max = 2);
  }
  mirror ({0,-1}) {
    cells();
    squares ("sb", min = 1, max = 2);
  }
  save ("ss.mp4");
}

event stop (i = 170) {;}