/**
# Flow past a cylinder in 3D

<center> <video controls> <source
src="https://surfdrive.surf.nl/files/index.php/s/N79po3rbsweYyxi/download"
type="video/mp4"> <caption><p align="center">Via Surfdrive</caption>
</video></center>
 
 */
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "navier-stokes/perfs.h"
#define CYLINDER (sqrt(sq(x) + sq(y - 0.01*R)) - R)

double R = 1, U = 1, Re = 1000., c = 15., ue, nu;
int maxlevel = 9;

u.n[left]  = dirichlet (U);
uf.n[left] = U;
u.t[left]  = dirichlet (0);
p[left]    = dirichlet (0.);
pf[left]   = dirichlet (0.);

u.n[right] = neumann (0.);
p[right]   = neumann (0.);
pf[right]  = neumann (0.);

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
#if (dimension == 3)
u.r[embed] = dirichlet (0.);
u.r[left]  = dirichlet (0);
#endif
FILE * fp;
face vector muc[];
int main() {
  nu = U*R/Re;
  periodic (top);
#if (dimension == 3)
  periodic (back);
#endif
  L0 = 30;
  X0 = -L0/3.5;
  Y0 = Z0 = -L0/2.;
  mu = muc;
  char logname[99];
  ue = U/c;
  sprintf (logname, "log3rd%g3D%d-%g",Re, maxlevel, c);
  fp = fopen (logname, "w");
  N = 64;
  run();
}

event init (t = 0) {
  refine (CYLINDER < 0.2*R && CYLINDER > -0.2*R && level < maxlevel);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = CYLINDER;
  boundary ({phi});
  fractions (phi, cs, fs);
  foreach() {
    u.x[] = cs[] > 0;
    u.y[] = noise()*U/200.;
#if (dimension == 3)
    u.z[] = noise()*U/200.;
#endif
  }
}

event damp (i++) {
  coord Uinf = {U, 0, 0};
  foreach() {
    if (fabs(x - (X0 + L0/2.)) > 4*L0/10.) 
      foreach_dimension()
	u.x[] += dt*(Uinf.x - u.x[])/2.;
  }
  boundary ((scalar*){u});
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*nu; 
  boundary ((scalar*){muc});
}

void prolongate_ratio (Point point, scalar s) {
  foreach_child() {
    if (s[] != nodata)
      s[] += s[]*Delta;
  }
}

event adapt (i++) {
  scalar res[];
  foreach() {
    res[] = nodata;
    if (cs[] > 0 && cs[] < 1) {
      res[] = U/sqrt(R*nu);
    }
  }
  res.prolongation = prolongate_ratio;
  adapt_wavelet ((scalar*){res, u}, (double[]){1., ue, ue, ue}, maxlevel, 4);
  unrefine (level > 4 && (x - X0) < L0/10. && (X0 + L0 - x) < L0/10.);
}

double embed_interpolate2 (Point point, scalar s, coord p) {
  int i = sign(p.x), j = sign(p.y);
#if dimension == 2
  if (cs[i] && cs[0,j] && cs[i,j])
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) + 
	    (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
#else //dimension == 3, see cartesian-common.h
  int k = sign(p.z);
  x = fabs(p.x); y = fabs(p.y); z = fabs(p.z);
  /* trilinear interpolation */
  if (cs[i] && cs[0,j] && cs[i,j] && cs[0,0,k] &&
      cs[i,0,k] && cs[0,j,k] && cs[i,j,k]) {
    return (((s[]*(1. - x) + s[i]*x)*(1. - y) + 
	     (s[0,j]*(1. - x) + s[i,j]*x)*y)*(1. - z) +
	    ((s[0,0,k]*(1. - x) + s[i,0,k]*x)*(1. - y) + 
	     (s[0,j,k]*(1. - x) + s[i,j,k]*x)*y)*z);
  }
#endif
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (cs[i])
	val += fabs(p.x)*(s[i] - s[]);
      else if (cs[-i])
	val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}

event logger (i += 5) {
#if (dimension == 2)
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  fprintf (fp, "%d %g %g %g %g %g %ld\n",
	   i, t, Fp.x, Fp.y, Fmu.x, Fmu.y, grid->n);
#elif (dimension == 3)
  double Fpx = 0., Fpy = 0, Fpz = 0;
  foreach (reduction(+:Fpx) reduction(+:Fpy) reduction(+:Fpz))
    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);
      double Fn = area*embed_interpolate2 (point, p, b);
      Fpx += Fn*n.x;
      Fpy += Fn*n.y;
      Fpz += Fn*n.z;
    }
  fprintf (fp, "%d %g %g %g %g %ld\n",
	   i, t, Fpx, Fpy, Fpz, grid->n);
  fflush (fp);
  printf ( "%d %g %g %g %ld %d %d %d\n",
	   i, t, Fpx, Fpy, grid->n, mgpf.i, mgp.i, mgu.i);
#endif
}

#include "view.h"
#include "lambda2.h"
event movies (t += 0.2) {
#if (dimension == 2)
  view (fov = 7, width = 1500, height = 400, tx = -0.25);
  scalar omega[];
  vorticity (u, omega);
  boundary ({omega});
  translate (z = 0.05) {
    draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
    draw_vof ("cs", "fs", lw = 2);
  }
  squares ("omega", min = -2, max = 2, linear = true, map = cool_warm);
  cells();
#elif (dimension == 3)
  scalar l2[];
  lambda2 (u, l2);
  view (fov = 18.1854, quat = {0.431384,-0.216693,-0.317091,0.816338},
	tx = 0, ty = 0, bg = {0.3,0.4,0.6}, width = 1080,
	height = 1080, samples = 3);
  isosurface ("l2", -0.01);
  cells (alpha = -L0/2);
  draw_vof ("cs", "fs", fc = {0.5,0.1,0.2});
#endif
  char str[99];
  sprintf (str, "Re = %g, C = %g, ML= %d", Re, c, maxlevel);
  draw_string (str, 1, lw = 3, lc = {1, 0, 1});
  save ("mov20003D.mp4");
}

event dumper (t += 10) {
  
  p.nodump = false;
  char str[99];
  sprintf (str, "dump3D%g", t);
  dump(str);
}

 event stop (t = 150) {
  fclose (fp);
}
