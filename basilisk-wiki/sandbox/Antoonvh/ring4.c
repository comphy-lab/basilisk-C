/**
# Two opposing Gaussian vortex rings  

Two Gaussian vortex rings in a triply-periodic domain.

<div class="figure">
<video controls="" preload="metadata" width="900">
<source src="https://surfdrive.surf.nl/files/index.php/s/HGUFYkO3gcHxCiK/download" type="video/mp4">
Your browser does not support the video tag. </video>
<p class="caption">
Volumetric rendering of the negative $\lambda_2$ field (via surfdrive)
</p>
</div>
*/
#include "grid/octree.h"
#include "nsf4t.h"
scalar * tracers = NULL;
#include "lambda2.h"
#include "view.h"
#include "filaments.h"
#include "bwatch.h"

#define R      (sqrt(sq(x) + sq(y)))
#define radi(zp)   (R <= 0.01 ? major_rt				\
		    : (sqrt(sq(x - major_rt*(x/R)) +			\
			    sq(y - major_rt*(y/R)) + sq(z - zp))))

double ue = 5e-3;
int maxlevel = 9;
int segs = 654;
double major_rt = 3.0; // Ring' major radius (minor = 1)
double d1 = 10;    // Initial distance
double muv = 1./3000.;   // Fluid's viscosity

double amp = 0.1;
int ampn = 7;

int main(int argc, char ** argv) {
  if (argc > 1)
    ue = atof (argv[1]);
  if (argc > 2)
    maxlevel = atoi (argv[2]);
  if (argc > 3)
    segs = atoi (argv[3]);
  if (argc > 4)
     muv = atof (argv[4]);
  if (pid() == 0)
    printf ("# ue = %g, ilev = %d, segs = %d, Re = %g \n",
	    ue, maxlevel, segs, 1./muv);
  foreach_dimension()
    periodic(left);
  L0 = 40*major_rt;
  origin (-L0/2 + pi, -L0/2 + exp(1), -L0/2 - sqrt(2));
  const scalar muc[] = muv;
  nu = muc;
  run();
}

coord ring1 (double t) {
  coord C;
  C.x = major_rt*cos(t);
  C.y = major_rt*sin(t);
  C.z = d1 + amp*cos(ampn*t); ;
  return C;
}

coord ring2 (double t) {
  coord C;
  C.x =  major_rt*cos(t);
  C.y = -major_rt*sin(t);
  C.z = -d1 + amp*sin(ampn*t);
  return C;
}

event init (t = 0) {
  refine ((radi(-d1) < 12 || radi(d1) < 12) && level < maxlevel - 2);
  refine ((radi(-d1) < 6  || radi(d1) < 6)  && level < maxlevel - 1);
  refine ((radi(-d1) < 3  || radi(d1) < 3)  && level < maxlevel);
 vector A[], omg[], omg2[];
 TOLERANCE = 1e-5;
 get_vor_vector (omg,  ring1, 0, 2*pi, segs, Lamb_Oseen);
 get_vor_vector (omg2, ring2, 0, 2*pi, segs, Lamb_Oseen);
 foreach_dimension()
   A.x.prolongation = refine_4th;
 foreach() {
   foreach_dimension() {
     A.x[] = 0.;
     omg.x[] += omg2.x[];
   }
 }
 foreach_dimension() {
   stats o = statsf (omg.x);
   foreach() 
     omg.x[] -= o.sum/o.volume;
 }
 boundary ((scalar*){A});
 foreach_dimension()
   poisson (A.x, omg.x);
 vector uc[];
 foreach_dimension()
   uc.x.prolongation = refine_4th;
 foreach() { 
   foreach_dimension() {
     uc.x[] = (8.*(A.z[0,1] - A.z[0,-1]) + A.z[0,-2] - A.z[0,2] -
	       8.*(A.y[0,0,1] - A.y[0,0,-1]) - A.y[0,0,-2] + A.y[0,0,2])/(12*Delta);
   }
 }
 boundary ((scalar*){uc});
 vector_to_face (uc);
 project (u, p);
}

event dumping (t = {5, 100, 200, 300, 350, 375, 400, 450, 500}) {
  vector uc[];
  face_to_vector (uc);
  char fname[99];
  sprintf (fname, "dump%g", t);
  dump (fname, (scalar*){uc});
}

event logger (i += 5) {
  double ke2 = 0., vd2 = 0.;
  double ke4 = 0., vd4 = 0.;
  vector uv[];
  scalar dissv[], Ev[];
  face_to_vector(uv);
  foreach(reduction (+:ke2) reduction (+:vd2)) {
    foreach_dimension() {
      ke2 += dv()*sq(uv.x[]);
      vd2 += dv()*(sq(uv.x[1]     - uv.x[-1]) +
		   sq(uv.x[0,1]   - uv.x[0,-1]) +
		   sq(uv.x[0,0,1] - uv.x[0,0,-1]))/sq(2.*Delta);
    }
  }
  foreach_dimension() {
    uv.x.prolongation = refine_vert5;
    uv.x.restriction = restriction_vert;
  }
  dissv.prolongation = refine_vert5;
  dissv.restriction = restriction_vert;
  Ev.prolongation = refine_vert5;
  Ev.restriction = restriction_vert;
  foreach() {
    foreach_dimension()
      uv.x[] = FACE_TO_VERTEX_4 (u.x);
  }
  boundary ((scalar*){uv});
  foreach() {
    Ev[] = dissv[] = 0;
    foreach_dimension() {
      Ev[] += sq(u.x[]);
      dissv[] += (sq(8*(u.x[1] - u.x[-1]) + u.x[-2] - u.x[2]) +
		 sq(8*(u.x[0,1] - u.x[0,-1]) + u.x[0,-2] - u.x[0,2]) +
		 sq(8*(u.x[0,0,1] - u.x[0,0,-1]) + u.x[0,0,-2] - u.x[0,0,2]))/sq(12*Delta);
    }
  }
  boundary ({dissv, Ev});
  foreach(reduction (+:ke4) reduction (+:vd4)) {
    coord cc = {x, y , z};
    foreach_child() {
      coord ccc = {0};
      foreach_dimension()
	ccc.x = cc.x + child.x*Delta/sqrt(3);
      ke4 += dv()*interpolate_vertex_4 (Ev,    ccc.x, ccc.y, ccc.z);
      vd4 += dv()*interpolate_vertex_4 (dissv, ccc.x, ccc.y, ccc.z);
    }
  }
  ke2 /= 2.;
  ke4 /= 2.;
  vd2 *= muv;
  vd4 *= muv;
  fprintf (stderr, "%d %g %d %d %d %d %ld %d %g %g %g %g\n",
	   i, t, mgp.i, mgp.nrelax, mgp2.i, mgp2.nrelax,
	   grid->tn, grid->maxdepth, ke2, vd2, ke4, vd4);
}

event adapt (i++) 
  adapt_flow (ue, 99, 1);

event mov (i += 5) {
  scalar l2[];
  vector uc[];
  face_to_vector (uc);
  lambda2 (uc, l2);
  view (fov = 12, phi = 0.5, theta = 0.4);
  isosurface ("l2", -0.001);
  cells (n = {0,1,0});
  char str[99];
  sprintf (str,"t = %3.3g", t);
  draw_string (str, size = 25, pos = 2);
  save ("l2.mp4");
}

#if BWATCH
event movb (t += 1) {
  scalar l2[];
  vector uc[];
  face_to_vector (uc);
  lambda2 (uc, l2);
  foreach()
    l2[] = l2[] > 0 ? 0 : -l2[];
  boundary ({l2});
  static FILE * fp = popen ("ppm2mp4 l2b.mp4", "w");
  watch (O = {30, 30., 100} ,
	 fov = 55, nx = 900, ny = 900); 
  volume (l2, cols = true, sc = 0.005, mval = 0.001,
	  min = -0.1, max = 0.1, shading = 1);
  store (fp);
  plain();
}
#endif

event stop (t = 1000) {
  return 1;
}
