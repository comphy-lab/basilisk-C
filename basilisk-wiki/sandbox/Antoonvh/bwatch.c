/**
# `bwatch`; a Basilisk-based ray-casting rendering system

We visualize a twisted and pink version of [Alexis Berny's `csg`
geometry](/sandbox/aberny/csgBool.c) against some squares, the grid
and a mirror.

![](bwatch/nja.mp4)(loop)
*/
#include "grid/octree.h"
#include "bwatch.h"

double cubeF (double x, double y, double z, coord center, double size) {
  double P1_Plus = x - size/2. + center.x;
  double P1_Minus = x + size/2. + center.x;
  double P1 = max (P1_Plus, -P1_Minus);
  double P2_Plus = y - size/2. + center.y;
  double P2_Minus = y + size/2. + center.y;
  double P2 = max (P2_Plus, -P2_Minus);
  double P3_Plus = z - size/2. + center.z;
  double P3_Minus = z + size/2. + center.z;
  double P3 = max (P3_Plus, -P3_Minus);
  double c = max (P1, max (P2, P3));
  return c;
}

double spherez (double x, double y, double z, coord center, double radius)  {
  return (sq(x - center.x) + sq (y - center.y) + sq (z - center.z)
          - sq (radius));
}

double geometry (double x, double y, double z) {
  coord center = {0,0,0};
  double s = spherez (x, y, z, center, 0.25);
  double c = cubeF (x, y, z, center, 0.38);
  double sIc = max (s, c);
  double cylinder1 = sq(x) + sq(y) - sq(0.12);
  double cylinder2 = sq(z) + sq(y) - sq(0.12);
  double cylinder3 = sq(x) + sq(z) - sq(0.12);
  double cylinderInter = min (cylinder1, cylinder2);
  double cylinderUnion = min (cylinderInter, cylinder3);
  return max(sIc, -cylinderUnion);
}

#define Xcc (x*cos(twist*y) + z*sin(twist*y))
#define Zcc (z*cos(twist*y) - x*sin(twist*y))
double twist = 4;

int main() {
  X0 = Y0 = Z0 = -L0/2 ;
  N = 64;
  init_grid (N);
  scalar f[];
#if TREE
  f.prolongation = f.refine = fraction_refine;
  do {
    fraction (f, geometry (Xcc, y, Zcc));
    boundary ({f});
  } while (adapt_wavelet ({f}, (double[]){0.005}, 9).nc > grid->tn/2000);
#else
  fraction (f, geometry (Xcc, y, Zcc));
#endif
  foreach()
    f[] = 1 - f[];
  boundary    ({f});
  restriction ({f});
  scalar s[];
  foreach()
    s[] = (point.i + point.j + point.k) % 2;
  
  // Draw a background, some squares, a slice of the grid, the vof object and a reflector:
  sphere (11, mat = {.dull = true});
  quadriangles (s, -0.3, mat = {.min = -0.25, .max = 1.25});
  lattice (-0.3);
  sketch_vof (f, mat = {.col = {250, 1, 150}});
  disk (0.3, {-.3, 0.2, -0.2}, {1, -0.7, 1.5});
  
  FILE * fp = popen ("ppm2mp4 nja.mp4", "w");
  int frames = 60;
  
  // Rotate and jump around object
  double R = 5*L0;
  for (int frame = 0; frame <= frames; frame++) {
    cam.O   = (coord){X0 + L0/2 + R*sin(2*pi*frame/frames),
		      Y0 + L0/2 + sin(4*pi*frame/frames) ,
		      Z0 + L0/2 + R*cos(2*pi*frame/frames)};
    store (fp);
    printf ("frame: %d\n", frame);
  }
  // Rotate around center-ray
  for (int frame = 0; frame <= frames; frame++) {
    cam.up   = (coord) {sin(2*pi*frame/frames),
			cos(2*pi*frame/frames) ,
			0.};
    store (fp);
    printf ("frame: %d\n", frame);
  }
  // Zoom using fov
  for (int frame = 0; frame <= frames; frame++) {
    cam.fov   = 1.1 + 0.5*sin(2*pi*frame/frames);
    store (fp);
    printf ("frame: %d\n", frame);
  }
  // Dolly zoom
  double Oz = cam.O.z;
  for (int frame = 0; frame <= frames; frame++) {
    cam.O.z   = Oz*(1 + 0.95*cos(2*pi*frame/frames));
    store (fp);
    printf ("frame: %d\n", frame);
  }
  // Look about
  for (int frame = 0; frame <= frames; frame++) {
    cam.poi.x = L0*sin(2*pi*frame/frames);
    store (fp);
    printf ("frame: %d\n", frame);
  }
  pclose (fp);
  plain();
}
