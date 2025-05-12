/**
# Flow past a fractal tree

<div class="figure">
<video controls="" preload="metadata" width="700px">
<source src="https://surfdrive.surf.nl/files/index.php/s/5eWRNDyDKHLDNY7/download" type="video/mp4">
Your browser does not support the video tag. </video>
<p class="caption">
Flow past a fractal tree with one generation of branches (via surfdrive)
</p>
</div>

<div class="figure">
<video controls="" preload="metadata" width="700px">
<source src="https://surfdrive.surf.nl/files/index.php/s/6JmnOyn4DsOjoPe/download" type="video/mp4">
Your browser does not support the video tag. </video>
<p class="caption">
Flow past a fractal tree with two generations of branches (via surfdrive)
</p>
</div>


See also this [page](endmovie.c), for the generation of the last
section of the movie.
 */
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "navier-stokes/perfs.h"
#include "treegen.h"
#include "view.h"
#include "lambda2.h"

double U = 1, Re = 1000., c = 10., ue, nu;
int maxlevel = 11;
int minlevel = 5;
double FILLVAL = 1e5;

scalar J[], res[];

u.n[left]  = dirichlet (U);
uf.n[left] = U;
u.t[left]  = dirichlet (0.);
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

void prolongate_ratio (Point point, scalar s) {
  foreach_child() {
    if (s[] != FILLVAL)
      s[] += s[]*Delta;
  }
}

FILE * fp;
face vector muc[];
int main() {
  nu = U*stem_rad/Re;
#if (dimension == 3)
  periodic (back);
#endif
  L0 = 100;
  X0 = -L0/3.;
  Z0 = -L0/2.;
  mu = muc;
  char logname[99];
  ue = U/c;
  sprintf (logname, "log3rd%g3D%d-%g",Re, maxlevel, c);
  fp = fopen (logname, "w");
  N = 128;
  run();
}

Branch * trees;
event init (t = 0) {
  levels = 1;
  srand (0); //Reproducible and identical among threads
  trees = tree_skeleton();
  res.prolongation = prolongate_ratio;
  do {
    tree_interface (trees, cs, fs, J);
    boundary ({J});
    fractions_cleanup (cs, fs);
    foreach() {
      res[] = FILLVAL;
      if (cs[] > 0 && cs[] < 1) {
	res[] = U/sqrt(trees[(int)(J[] + 0.5)].R*nu);
      }
    }
    boundary ({res});
  } while (adapt_wavelet ({res}, (double[]){2.}, maxlevel, minlevel).nf > grid->tn/1000);
  tree_interface (trees, cs, fs, J);
  fractions_cleanup (cs, fs);
  
  dump();
  view (height = 1080, width = 1080, fov = 9, ty = -0.2);
  draw_vof ("cs", "fs", color = "J");
  cells();
  save ("tree.ppm");
  foreach() {
    u.x[] = cs[] > 0;
    u.y[] = noise()*U/200.;
#if (dimension == 3)
    u.z[] = noise()*U/200.;
#endif
  }
}

event properties (i++) {
  foreach_face()
    muc.x[] = fm.x[]*nu; 
  boundary ((scalar*){muc});
}
	  
event damp (i++) {
  coord Uinf = {U, 0, 0};
  double Lay = L0/10;
  double tau = 0.5;
  foreach() {
    if (x - X0  < Lay || (X0 + L0) - x < Lay) {
      foreach_dimension()
	u.x[] += dt*(Uinf.x - u.x[])/tau;
    }
  }
  boundary ((scalar*){u});
}

event adapt (i++) {
  foreach() {
    res[] = FILLVAL;
    if (cs[] > 0 && cs[] < 1) {
      res[] = U/sqrt(trees[(int)(J[] + 0.5)].R*nu);
    }
  }
  boundary ({res});
  adapt_wavelet ((scalar*){res, u},
		 (double[]){2., ue, ue, ue}, maxlevel, minlevel);
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
  if (pid() == 0) {
    fprintf (fp, "%d %g %g %g %g %ld\n",
	     i, t, Fpx, Fpy, Fpz, grid->n);
    fflush (fp);
    printf ( "%d %g %g %g %ld %d %d %d\n",
	     i, t, Fpx, Fpy, grid->n, mgpf.i, mgp.i, mgu.i);
  }
#endif
}

event movies (t += 0.2) {
  scalar l2[];
  lambda2 (u, l2);
  view (fov = 12, theta = -0.8, phi = 0.4, 
	tx = 0, ty = -0.15, bg = {65./256,157./256,217./256},
	width = 1080, height = 1080);
  isosurface ("l2", -0.1);
  if (t < 50 || t > 100)
    cells();
  draw_vof ("cs", "fs", fc = {98./256,78./256,44./256});
  char str[99];
  sprintf (str, "Re = %g, C = %g, ML= %d", Re, c, maxlevel);
  draw_string (str, 1, lw = 3, lc = {1, 0, 1});
  double py = 20;
  translate (y = -py)
    squares ("u.x", n = {0,1,0}, alpha = py,
	     min = -1.1, max = 1.1, map = cool_warm);
  save ("movtree.mp4");
}

event dumper (t += 10) {
  p.nodump = false;
  char str[99];
  sprintf (str, "dump3D%g", t);
  dump(str);
}

event stop (t = 150) {
  fclose (fp);
  free (trees);
}
