/**
![[Nicolas Tesla](https://www.nemokennislink.nl/publicaties/nikola-tesla-de-excentriekeling/)](https://assets.kennislink.nl/system/files/000/178/381/medium/Nicola-Tesla-in-39-True-photo-without-photoshop-nikola-tesla-26914162-940-1260.jpg?1495013772)



# Testing Tesla's valve

We test a simple Tesla valve design:

![Interesting history effects](tesla_valve/tesla.mp4)

The flow rate is direction dependend:

~~~gnuplot
set xlabel 'Time'
set ylabel 'Velocity'
set grid
plot 'out' u 1:3 t "->", 'log' u 1:3 t "<-", 0 lc 'black'
~~~

Perhaps other parameters work even better
 */
#include "embed.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "PointTriangle.h"
#include "view.h"

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);

double muv = 0.0001;

face vector muc[];

double angle = 20; //channel Angle
double L = 2.8;    //Length of a channel
double R = 0.3;    //Radius of cirle
double sect = 25;  //sections in circle
 

int get_coords_Tesla (coord * p) {
  double dydx = tan(pi*angle/180.);
  double O = R/(dydx); 
  double y1 = dydx*L/2, y2 = dydx * (L/2 - O);
  int c = 0;
  double start = -2;
  for (int j = 0; j < 5; j++) {
    for (int l = 1; l > -2; l -= 2) {
      p[c++] = (coord){start, 0.};
      p[c++] = (coord){start + L/2, l*y1};
      start += L/2 - O;
      p[c++] = (coord){start, l*y2};
      p[c++] = (coord){start + y2/dydx, 0.};
      for (double theta = 0; theta < pi; theta += pi/sect + 1e-3) {
	double xr = start + O, yr = l*(y1 - R);
	p[c++] = (coord){xr + R*sin(theta), yr + R*cos(theta)};
	p[c++] = (coord){xr + R*sin(theta + pi/sect), yr + R*cos(theta + pi/sect)};
      }
      start += L/2 - O;
    }
  }
  p[c].x = nodata;
  return c;
}

int main() {
  periodic (left);
  double dydx = tan(pi*angle/180.);
  double O = R/(dydx); 
  L0 = 4*(L - 2*O);
  printf ("%g %g\n", L0, O);
  N = 256;
  Y0 = -2;
  const face vector av[] = {0.01, 0};
  a = av;
  mu = muc;
  run();
  const face vector av2[] = {-0.01, 0};
  a = av2;
  run();
}

event properties (i++) {
  foreach_face()
    muc.x[] = muv*fs.x[];
  boundary((scalar*){muc});
}

event init (t = 0) {
  coord p[999];
  get_coords_Tesla(p);
  vertex scalar d[];
  foreach_vertex() {
    double min = HUGE;
    coord * a = p;
    coord vert = {x, y};
    while (a->x != nodata) {
      coord r;
      double d, d2 = PointSegmentDistance (&vert, a, a + 1, &r, &d);
      if (d2 < min)
	min = d2;
      a += 2;
    }
    d[] = min > 0 ? -sqrt (min) : -100;
  }  fractions (d, cs, fs, -0.1);
}

event adapt (i++) {
  adapt_wavelet({cs, u}, (double[]){1e-3, 1e-2, 1e-2}, 8);
}

event movie (t += 0.1) {
  scalar U[];
  foreach() 
    U[] = sqrt(sq(u.x[]) + sq(u.y[]));
  boundary ({U});
  view (fov = 10, tx = -0.5,
	height = 300, width = 600);
  draw_vof ("cs", "fs", filled = -1, fc = {0.6, 0.6, 0.6});
  draw_vof ("cs", "fs");
  squares ("U", min = 0, max = 3e-1);
  if (constant(a.x) > 0)
    draw_string ("  Flow ->", 1, lw = 4, size = 30);
  else
    draw_string ("  <- Flow", 1, lw = 4, size = 30);
  save ("tesla.mp4");
  if (fabs(t - 49.9) < 1e-5 && constant(a.x) > 0) 
    save ("tesla.png");
}

event flow_rate (i++) {
  double flx = 0;
  foreach_face(x) {
    if (fabs(x - X0) < 1e-5) 
      flx += uf.x[]*fm.x[]*Delta;
  }
  if (constant(a.x) > 0)
    printf ("%g %g %g\n", t, flx, statsf(u.x).sum);
  else
    fprintf (stderr, "%g %g %g\n", t, flx, statsf(u.x).sum);
}

event stop (t = 50);
