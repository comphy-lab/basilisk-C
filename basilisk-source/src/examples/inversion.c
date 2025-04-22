// #include "grid/multigrid3D.h"
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"

scalar f[], * interfaces = {f};

#define rho1 900.
#define rho2 1000.
#define mu1  0.1
#define mu2  0.001
#define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
#define mu(f) (1./(clamp(f,0,1)*(1./(mu1) - 1./(mu2)) + 1./(mu2)))

#define H 0.1
#define G 9.81
#define Ug ((rho2 - rho1)/rho1*sqrt(H*G/2.))
#define tc (H/(2.*Ug))

face vector alphav[], muv[], av[];
scalar rhov[];

int maxlevel = 6;

#if 0
uf.n[left] = 0.;
uf.n[right] = 0.;
p[left] = neumann(0);
p[right] = neumann(0);

f[top] = 0.;
f[right] = 0.;
f[left] = 0.;
f[bottom] = 0.;

#if dimension == 3
f[front] = 0.;
f[back] = 0.;
#endif
#endif

timer tt;

int main (int argc, char * argv[]) {
  maxlevel = argc > 1 ? atoi(argv[1]) : 7;
  size (H);
  origin (-H/2., -H/2., -H/2.);
#if !TREE
  N = 1 << maxlevel;
#endif
  a = av;
  mu = muv;
  alpha = alphav;
  rho = rhov;
  f.sigma = 0.045;
  DT = 2e-2;
  tt = timer_start();
  run();
}

event init (i = 0) {
#if TREE
  scalar f1[];
  foreach()
    f1[] = (x <= 0 && y <= 0 && z <= 0);
  astats s;
  do {
    s = adapt_wavelet ({f1}, (double[]){0.0}, maxlevel, list = NULL);
    foreach()
      f1[] = (x <= 0 && y <= 0 && z <= 0);
  } while (s.nf);
  foreach()
    f[] = (x <= 0 && y <= 0 && z <= 0);
#else
  foreach()
    f[] = (x <= 0 && y <= 0 && z <= 0);
#endif
}

event acceleration (i++) {
  foreach_face(y)
    av.y[] -= G;
}

event properties (i++) {
#if TREE
  f.prolongation = refine_bilinear;
  f.dirty = true;
#endif

  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    muv.x[] = fm.x[]*mu(ff);
  }
  foreach()
    rhov[] = cm[]*rho(f[]);

#if TREE
  f.prolongation = fraction_refine;
  f.dirty = true;
#endif
}

event logfile (i++; t <= 20) {
  double ke1 = 0., ke2 = 0., vd = 0., vol1 = 0.;
  double ep1 = 0., ep2 = 0.;
  double er1 = 0., er2 = 0.;
  double area = 0.;
  int nc = 0;
  static long tnc = 0;
  foreach(reduction(+:ke1) reduction(+:ke2) reduction(+:vd)
	  reduction(+:vol1) reduction(+:ep1) reduction(+:ep2)
	  reduction(+:er1) reduction(+:er2) reduction(+:area)
	  reduction(+:nc)) {
    if (y > H/2. - H/8.)
      vol1 += f[]*dv();
    ep1 += rho1*f[]*G*(y + H/2.)*dv();
    ep2 += rho2*(1. - f[])*G*(y + H/2.)*dv();
    // interfacial area
    if (f[] > 1e-4 && f[] < 1. - 1e-4) {
      coord m = mycs (point, f);
      double alpha = plane_alpha (f[], m);
      coord p;
      area += sq(Delta)*plane_area_center (m, alpha, &p);
    }
    double w2 = 0.;
    foreach_dimension() {
      // kinetic energy
      ke1 += dv()*f[]*rho1*sq(u.x[]);
      ke2 += dv()*(1. - f[])*rho2*sq(u.x[]);
      // viscous dissipation
      vd += dv()*(sq(u.x[1] - u.x[-1]) +
		  sq(u.x[0,1] - u.x[0,-1]) +
		  sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
      // enstrophy
      w2 += sq(u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0]);
    }
    w2 /= sq(2.*Delta);
    er1 += dv()*f[]*w2;
    er2 += dv()*(1. - f[])*w2;
    nc++;
  }
  ke1 /= 2.;
  ke2 /= 2.;
  er1 /= 2.;
  er2 /= 2.;
  //  vd *= MU/vol;

  if (i == 0)
    fprintf (stderr,
	     "t ke1 ke2 ep1 ep2 er1 er2 R2 area mgp.i mgu.i nc time speed\n");
  double elapsed = timer_elapsed (tt);
  tnc += nc;
  fprintf (stderr, "%g %g %g %g %g %g %g %g %g %d %d %d %g %g\n",
	   t/tc,
	   ke1/(1./16.*rho1*sq(Ug)*cube(H)),
	   ke2/(1./16.*rho2*sq(Ug)*cube(H)),
	   ep1/(rho1*G*15.*sq(H)*sq(H)/128.),
	   ep2/(rho2*G*49.*sq(H)*sq(H)/128.),
	   er1/0.0733,
	   er2/1.3759,
	   8.*vol1/cube(H),
	   area,
	   mgp.i, mgu.i, nc, elapsed, tnc/elapsed);
#if 0
  nc = 0;
  foreach()
    nc++;
  fprintf (stderr, "nc: %d\n", nc);
  fflush (stderr);
#endif
}

#if 0
event movies (t += 0.1*tc) {
  char name[80];
  sprintf (name, "f-%d.ppm", maxlevel);
  static FILE * fp = fopen (name, "w");
  output_ppm (f, fp, min = 0, max = 1, n = 256);
}
#endif

#if !_MPI
event gfsview (i += 10) {
  scalar pid[];
  foreach()
    pid[] = tid();
#if dimension == 3
  static FILE * fp = popen ("gfsview3D inversion.gfv", "w");
#else
  static FILE * fp = popen ("gfsview2D inversion2D.gfv", "w");
#endif
  output_gfs (fp);
}
#endif

#if 0 //TREE && !_MPI
event gfsview (t += 0.1*tc) {
@if _MPI
  char name[80];
  sprintf (name, "output-%g.gfs", t);
  FILE * fp = fopen (name, "w");
@else
  static FILE * fp =
    popen ("gfsview3D ../inversion.gfv", "w");
@endif
  output_gfs (fp, translate = true);
@if _MPI
  fclose (fp);
@endif
  //  fprintf (fp, "Save stdout { format = PPM width = 512 height = 512 }\n");
}
#endif

#if 0
event snapshot (i = 100; i += 100) {
  dump (file = "snapshot", t = t);
  char name[80];
  sprintf (name, "snapshot-%d.gfs", i);
  scalar pid[];
  foreach()
    pid[] = tid();
  output_gfs (file = name);
}
#endif

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.005,0.005,0.005,0.005}, maxlevel);
}
#endif
