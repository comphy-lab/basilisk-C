/**
# Flux partitioning at an interface

A proof-of-concept for an embedded flux condtion where the scalar
spreads into either side of the interface. E.g. the heat from
irradiation is partitioned between the solid and the transparent
environment, whilst the temperature does not jump over the interface.

![Evolution of temperature and grid structure](test_flux/mov.mp4)

 */
#include "embed.h"
#include "diffusion.h"
#include "utils.h"
#include "run.h"
#include "view.h"
#define GEOM (0.1*sin(2*pi*x) - 0.4 + y)
scalar s1[], s2[], grad1[], grad2[];
(const) scalar FLUX;

s1[top] = dirichlet (0);
s1[embed] = neumann (grad1[]);
s2[embed] = neumann (grad2[]);
s2[bottom] = dirichlet (0);
foreach_dimension()
double interpolate_nearby_x  (Point point, scalar s, scalar cs, face vector fs,
			      coord n, coord p, double v[2], double d[2]) {
  foreach_dimension()
    n.x = - n.x;
  for (int i = 0; i < 2; i++) {
    d[i] = nodata;
    v[i] = nodata;
  }
  bool defined = true;
  foreach_dimension()
    if (defined && !fs.x[(n.x > 0.)])
      defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
      if (fs.x[i + (i < 0),j] && fs.y[i,j] && fs.y[i,j+1] &&
	  cs[i,j-1] && cs[i,j] && cs[i,j+1])
	v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = fs.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!fs.y[i,j,k+m] || !fs.y[i,j+1,k+m] ||
	    !fs.z[i,j+m,k] || !fs.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else
	break;
    }
  assert(v[1] != nodata);
  return v[0];
}

void set_interface_val (scalar T1, scalar T2,
			    face vector kappa1, face vector kappa2,
			scalar grad1, scalar grad2) {
  double v1[2], d1[2];//!!
  double v2[2], d2[2];//!!
  scalar cs2[];
  face vector fs2[];
  
  foreach() 
    cs2[] = 1. - cs[];
  foreach_face()
      fs2.x[] = 1. - fs.x[];
  
  boundary ({cs2, fs2});
  foreach() {
    if (cs[] > 0. && cs[] < 1.) {
      coord p, n;
      embed_geometry (point, &p, &n);
      double mua1 = 0., fa1 = 0.;
      double mua2 = 0., fa2 = 0.;

#if dimension == 2
      foreach_dimension()
	if (fabs(n.x) >= fabs(n.y)) {
	  interpolate_nearby_x (point, s1, cs, fs, n, p, v1, d1);
	  foreach_dimension()
	    n.x = -n.x;
	  interpolate_nearby_x (point, s2, cs2, fs2, n, p, v2, d2);
	}
      
#else // dimension == 3
      assert (dimension == 2);
#endif // dimension == 3
      
      foreach_dimension() {
	mua1 += kappa1.x[] + kappa1.x[1];
	fa1  += fs.x[] + fs.x[1];
	mua2 += kappa2.x[] + kappa2.x[1];
	fa2  += 2. - fs.x[] - fs.x[1];
      }
      
      double Ts0 = (v1[0]*(mua1/fa1)/(Delta*d1[0]) + v2[0]*(mua2/fa2)/(Delta*d2[0]) + FLUX[])/
	((mua1/fa1)/(Delta*d1[0]) + (mua2/fa2)/(Delta*d2[0]));
      grad1[]  = -(v1[0] - Ts0)/(d1[0]*Delta);
      grad2[]  = -(v2[0] - Ts0)/(d2[0]*Delta);
    }
  }
}

int main() {
  DT = 5;
  periodic (left);
  const scalar flx[] = 0.01;
  FLUX = flx;
  N = 64;
  run();
  
}

face vector kappa1[], kappa2[];

event init (t = 0) {
  vertex scalar phi[];
  foreach_vertex()
    phi[] = GEOM;
  fractions (phi, cs, fs);
  
  
  foreach() {
    s1[] = 0;
    s2[] = 0;
    grad1[] = grad2[] =  nodata;
  }
}

event properties (i++) {
  foreach_face() {
    kappa1.x[] = fs.x[]*0.00001846; // Air
    kappa2.x[] = (1. - fs.x[])*0.000001; // soil
  }
}

event set_dt(i++) {
  dt = dtnext(DT);
}

event set_interface_temp(i++) {
  set_interface_val (s1, s2, kappa1, kappa2, grad1, grad2);
}

void invert_embed (scalar s, face vector f) {
  foreach()
    s[] = 1. - s[];
  foreach_face()
    f.x[] = 1. - f.x[];
}

event tracer_diffusion (i++) {
  diffusion (s1, dt, kappa1);
  invert_embed (cs, fs);
  diffusion (s2, dt, kappa2);
  invert_embed (cs, fs);
}

event outputer (t += 20) {
  scalar s[];
  foreach()
    s[] = cs[]*s1[] + (1. - cs[])*s2[];

  view (fov = 20, width = 800, height = 400, ty = -0.5);
  translate (x = -1) {
    squares ("s", min = -1, max = 100, linear = true);
    draw_vof ("cs", "fs");
  }
  cells();
  save ("mov.mp4");
}

event adapt (i++) {
  adapt_wavelet ({cs, s1, s2}, (double[]){1e-6, 5, 5}, 8);
}
event stop (t = 3600);


