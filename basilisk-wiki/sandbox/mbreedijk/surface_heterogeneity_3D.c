#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "tracer.h"
#include "diffusion.h"
#include "profile5c.h"
#include "view.h"
#include <stdlib.h>
#include "lambda2.h"

double  heaviside ( double x, double minx1, double maxx1, double minx2, double maxx2 ) {
  double a = 0;
  if ((x > minx1 && x < maxx1) || (x > minx2 && x < maxx2))
    a = 1;
  return a;
}
    
scalar b[], * tracers = {b};
double b0, Nbv, k, nuv, kappav, scaling_e, L, L0;
double Pi1 = 2, Pi2, Pi3 = 1, Pi4 = 20000, Pi5 = 60;
double minx1, maxx1, minx2, maxx2, w, d;
double t_end = 150;

int maxlevel = 11, minlevel = 8;

b[bottom] = dirichlet (b0*heaviside(x,minx1,maxx1,minx2,maxx2));
b[top]    = neumann (sq(Nbv));
u.t[bottom] = dirichlet (0); //no slip
#if (dimension == 3)
u.r[bottom] = dirichlet (0); //no slip
#endif

face vector av[];

int main() { 
  periodic (left);
  periodic (back);
  Nbv = 1;
  d = 4;
  w = 0.5;
  Pi2 = d/w;
  
  // variables dimensionless groups
  L = Pi1*w;
  L0 = (2*w*Pi5)+d;
  nuv = sq(b0)/(sq(Nbv)*Nbv)/Pi4;
  kappav = nuv*Pi3;
  b0 = L*sq(Nbv);
  scaling_e = sq(Nbv/(b0))/w/L/L0;

  printf("w = %f\n",w);
  printf("L = %f\n",L);
  printf("L0 = %f\n", L0);
  printf("b0 = %f\n", b0);
  printf("Pi2 = %f\n", Pi2);

  // define source
  minx1 = (L0/2)-(w/2)-(d/2);
  maxx1 = minx1 + w;
  minx2 = minx1 + d;
  maxx2 = minx2 + w;
  
  const face vector muc[] = {nuv, nuv, nuv};
  // Link fields
  mu = muc;
  a  = av;

  p.prolongation = p.refine = refine_linear; //3rd-order interpolation (vertical)
  foreach_dimension()
    u.x.refine = refine_linear;
  N = 1 << minlevel; // = bitwise 2^minlevel
  run();
}

event init (t = 0) {
  refine (level < maxlevel && heaviside(x,minx1,maxx1,minx2,maxx2) && y < L); //Refine near heat source
  DT = 0.1/Nbv;
  foreach() 
    b[] = sq(Nbv)*y;
  boundary ({b});
}

event tracer_diffusion (i++) {
  const face vector kappa[] = {kappav, kappav, kappav};
  diffusion (b, dt, kappa);
}

event acceleration (i++) {
  coord g = {0, 1, 0};
  foreach_face()
    av.x[] = g.x*face_value(b, 0);
}

event adapt (i++) {
  adapt_wavelet ({b, u}, (double[]){b0/150., 0.01, 0.01, 0.01},maxlevel,minlevel); // refinement
}

event damp (i++) {
  foreach() {
    if (y > L0/2) {
      foreach_dimension()
	u.x[] *= exp(-(y - L0/2)/5.);
    }
  }
  boundary ((scalar*){u});
}


event time_series (i += 25) {
  double e = 0, diss = 0, sbflx = 0, e_dim = 0;
  foreach(reduction(+:diss) reduction(+:e) reduction(+:sbflx)) {
    foreach_dimension() {
      e += sq(u.x[])*dv(); 
      diss += dv()*(sq(u.x[1] - u.x[-1]) +
		    sq(u.x[0,1] - u.x[0,-1]) +
		    sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
    if (Delta > y)
      sbflx += (b[0,-1] - b[])*dv()/sq(Delta);
  }
  diss *= -nuv;
  sbflx *= nuv;
  e /= 2.;
  e_dim += e*scaling_e;
  
  static FILE * fp = fopen ("timeseries_3d", "w_two3d");
  if (i == 0)
    fprintf (fp, "t\ti\tn\twct\tspeed\te\te_dim\tdiss\tsbflx\n");
  fprintf (fp, "%g\t%d\t%ld\t%g\t%g\t%g\t%g\t%g\t%g\n",
	   t, i, grid->tn, perf.t, perf.speed, e, e_dim, diss, sbflx);
  fflush (fp);
}


event mov (t += 1) {
  output_ppm (b, file = "b_3d.mp4", n = 300, box = {{0,0},{L0,L0}},min = -1, max = 1.5, linear = true);
  scalar lev[];
  foreach()
    lev[] = level;
  output_ppm (lev, file = "l_3d.mp4", n = 300, box = {{0,0},{L0,L0}}, min = 1, max = maxlevel + 0.5);
  foreach() {
    lev[] = 0;
    foreach_dimension()
    lev[] += sq(u.x[]);
  }
  output_ppm (lev, file = "e_3d.mp4", n = 300, box = {{0,0},{(L0),L0}}, min = 0);
  vector db[];
  scalar m[];
  boundary ({b});
  gradients ({b}, {db});
  foreach() {
    m[] = 0;
    foreach_dimension()
      m[] += sq(db.x[]);
    if (m[] > 0)
      m[] = log(sqrt(m[])+1.);
  }
  boundary ({m});
  output_ppm (m, file = "db_3d.mp4", n = (1 << maxlevel), min = 0, max = 2, linear = true);
}

event mp4 (t += 1) {
  scalar l2[];
  lambda2 (u,l2);  
  scalar e_movie[];
  foreach() {
    e_movie[] = 0;
    foreach_dimension()
      e_movie[] += sq(u.x[]);
  }
  view (fov = 30, tx = -0.5, ty = -0.1, phi = 0.3);
  squares("b", linear = true,  n = {0,1,0}, alpha = 0, min = 0, max = b0);
  cells();
  isosurface("l2", -0.001, color = "u.y");
  save ("movie_3d.mp4");
}

event stop (t = t_end);
