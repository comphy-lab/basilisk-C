/**
# Volumetric visualization

A proof of a inefficient concept.

![It works a bit, perhaps a bit of shading helps to see
 depth.](plume/mov.mp4)
*/
#include "grid/multigrid3D.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "integrator2.h"

scalar s[], * tracers = {s};

double sigma = 0.1; // size
double att = 0.1;   // Attenuation coefficient
double prf = 4e3;    // Light intensity of some sort

trace
double get_c (ray r, scalar Sl, scalar Sc, scalar s) {
  double c = 0;
  foreach_ray_cell_intersection(r) {
    double l = 0;
    foreach_dimension()
      l += sq(_a[1].x - _a[0].x);
    l = l > 0 ? sqrt(l) : 0;
    c += l*
      interpolate(s, (_a[0].x +_a[1].x)/2, (_a[0].y +_a[1].y)/2, (_a[0].z +_a[1].z)/2)*
      interpolate(Sl, (_a[0].x +_a[1].x)/2, (_a[0].y +_a[1].y)/2, (_a[0].z +_a[1].z)/2)*
      interpolate(Sc, (_a[0].x +_a[1].x)/2, (_a[0].y +_a[1].y)/2, (_a[0].z +_a[1].z)/2);
  }
  return fabs(c);
}

trace
void make_a_frame (scalar s, FILE * fp) {
  scalar Sl[], Sc[];
  coord ildir = {1, 1, 3};
  normalize (&ildir);
  integrate_dn (Sl, s, ildir);
  foreach()
    Sl[] = exp(-Sl[]/att);
  boundary ({Sl});
  coord icdir = {0.01, 0.01, 1};
  normalize (&icdir);
  integrate_dn (Sc, s, icdir);
  foreach()
    Sc[] = exp(-Sc[]/att);
  boundary ({Sc});
  fprintf (fp, "P6\n%d %d\n%d\n", cam.nx, cam.ny, 255); // .ppm header
  foreach_ray() {
    unsigned char px[3] = {10.*(ii + jj)/cam.nx + 20,
			   25.*(ii + jj)/cam.nx + 120,
			   50.*jj/cam.nx + 120};
    double c = get_c (_r, Sl, Sc, s);
    double occ = fabs(BC_int(_r.O, _r.dir, s));
    for (int i = 0; i < 3; i++) 
      px[i] =  min(max(px[i]*exp(-occ/att) + prf*c, 1), 255);
    fwrite (px, 1, 3, fp);
  }
  fflush (fp);
}

FILE * fp;
int main() {
  X0 = Y0 = Z0 = -0.5;
  N = 64;
  fp = popen ("ppm2mp4 -r 10 mov.mp4", "w");
  run();
}

event init (t = 0) {
  fraction (s, sq(sigma) - (sq(x) + sq(y + 0.25) + sq(z)));
  foreach() 
    u.y[] = 0.4*s[];
  boundary ({s, u});
}

event frames (t += 0.2) {
  make_a_frame (s, fp);
}

event stop (t = 5) {
  pclose (fp);
  return 1;
}
