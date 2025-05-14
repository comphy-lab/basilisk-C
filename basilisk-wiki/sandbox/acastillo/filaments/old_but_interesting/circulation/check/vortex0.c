#include "navier-stokes/centered.h"
#include "view.h"
#include "fractions.h"

static double omega_mean;

double lamb_oseen (double r, double a){
  return 1.0/(pi*sq(a)) * exp( -sq(r/a) );
}

int maxlevel = 10;
int minlevel = 5;

int main(){
  L0 = 16 ;
  X0 = Y0 = -L0/2 ;
  omega_mean = -0.5*1.0/sq(L0);
  //TOLERANCE = 1e-12 ;
  // double reynolds= 100 ;
  // const face vector muc[] = {1./reynolds,1./reynolds,1./reynolds} ;
  // mu = muc;

  periodic(top);
  periodic(left);
  init_grid (1 << maxlevel);
  display_control (maxlevel, 6, 12);
  run();
}

event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){1e-5,1e-5}, maxlevel, minlevel);
}


event init (t = 0){
  scalar c[], b[];

  double bsum = 0.;
  foreach(reduction(+:bsum)) {
    double r = sqrt(sq(x) + sq(y));
    b[] = lamb_oseen(r, 1.0);
    bsum += b[]*dv();
  }
  fprintf (stderr, "bsum (before): %g\n", bsum);

  omega_mean = -0.5*bsum/sq(L0);
  fprintf (stderr, "omega_mean: %g\n", omega_mean);

  bsum = 0.;
  foreach(reduction(+:bsum)) {
    double r = sqrt(sq(x) + sq(y));
    c[] = -sq(r)*omega_mean/2.0;
    b[] += 2*omega_mean;
    bsum += b[]*dv();
  }
  boundary ({c,b});
  fprintf (stderr, "bsum: %g\n", bsum);

  /** The Poisson equation is solved. */
  scalar res[];
  scalar * list = {res};

  timer t = timer_start();
  mgstats s = poisson (c, b, tolerance = 1e-15, minlevel = 4, res=list);
  double dt = timer_elapsed (t);
  printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

  coord f = {-1.,1.};
  foreach()
    foreach_dimension()
      u.y[] = f.y * center_gradient(c);
  boundary ((scalar *){u});

  scalar omega[];
  vorticity (u, omega);

  scalar psi[];
  foreach()
    psi[] = c[];
  boundary ({psi, omega});

  scalar res2[];
  scalar * list2 = {res2};
  t = timer_start();
  s = poisson (psi, omega, tolerance = 1e-15, minlevel = 4, res=list2);
  dt = timer_elapsed (t);
  printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);

  scalar e1[], e2[];
  foreach(){
    e1[] = c[]-psi[];
    e2[] = b[]-omega[];
  }

  /**
  The solution is displayed using bview.
  */

  squares ("c", spread = -1);
  isoline ("c", n = 21);
  save ("a.png");
  clear();

  squares ("b", spread = -1);
  isoline ("b", n = 21);
  save ("b.png");
  clear();

  squares ("res", spread = -1);
  save ("res.png");
  clear();

  squares ("psi", spread = -1);
  isoline ("psi", n = 21);
  save ("psi.png");
  clear();

  squares ("omega", spread = -1);
  isoline ("omega", n = 21);
  save ("omega.png");
  clear();

  squares ("res2", spread = -1);
  save ("res2.png");
  clear();

  squares ("u.x", spread = -1);
  isoline ("u.x", n = 21);
  save ("u.png");
  clear();

  squares ("u.y", spread = -1);
  isoline ("u.y", n = 21);
  save ("v.png");
  clear();

  squares ("e1", spread = -1);
  isoline ("e1", n = 21);
  save ("e1.png");
  clear();

  squares ("e2", spread = -1);
  isoline ("e2", n = 21);
  save ("e2.png");
  clear();
  dump ("dump0");

  stats s0 ;
  s0 = statsf (c); fprintf (stderr, "a: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (b); fprintf (stderr, "b: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (psi); fprintf (stderr, "psi: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (omega); fprintf (stderr, "omega: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (e1); fprintf (stderr, "e(psi): %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (e2); fprintf (stderr, "e(omega): %g %g %g \n", s0.min, s0.max, s0.sum);

  fprintf (stderr, "\n");

  event("logfile");
}

#define K0() (0.)
#define F0() (2.0 * omega_mean)
#define alpha_H 0.5
event end_timestep (i++){
  foreach(){
   coord b0 = { - K0(), - K0() }, b1 = { F0(), -F0() };
	  coord m0 = { 1. - alpha_H*dt*b0.x, 1. - alpha_H*dt*b0.y };
	  coord m1 = { - alpha_H*dt*b1.x, - alpha_H*dt*b1.y };
	  double det = m0.x*m0.y - m1.x*m1.y;
    coord r;
	  foreach_dimension() {
	     r.x = u.x[] + (1. - alpha_H)*dt*(b0.x*u.x[] + b1.x*u.y[]) + dt*g.x[];
	  }
	  foreach_dimension(){
      u.x[] = (m0.y*r.x - m1.x*r.y)/det - dt*g.x[];
    }
  }
  boundary ((scalar *){u});
}










event movie (i++) {
  squares ("u.x", spread = -1);
  isoline ("u.x", n = 21);
  save ("u.mp4");

  squares ("u.y", spread = -1);
  isoline ("u.y", n = 21);
  save ("v.mp4");

  scalar omega[];
  vorticity (u, omega);
  squares ("omega", spread = -1);
  isoline ("omega", n = 21);
  save ("omega.mp4");

  scalar psi[];
  poisson (psi, omega, tolerance = 1e-10, minlevel = 4);
  squares ("psi", spread = -1);
  isoline ("psi", n = 21);
  save ("psi.mp4");
}

event output (t = 2.0)
  dump();






















void global_quantities(FILE * fp, scalar omega, vector u){
  double ekin = 0.0, circ = 0.0, enst = 0.0, area_tot=0.0;
  foreach(reduction(+:ekin), reduction(+:circ), reduction(+:enst), reduction(+:area_tot)){
    ekin += dv() * (sq(u.x[]) + sq(u.y[]));
    circ += dv() * omega[];
    enst += dv() * sq(omega[]);
    area_tot += dv();
  }

  stats s0 = statsf (omega);

  double circ2[1] = {0.}, mu_x[1] = {0.}, mu_y[1] = {0.};
  foreach(reduction(+:circ2[:1]), reduction(+:mu_x[:1]), reduction(+:mu_y[:1])){
    coord p = {x,y,z};
    int nj = 0;
    if (sqrt(sq(x)+sq(y)) < L0*3/8.){
      circ2[nj] += dv() * omega[] ;
      mu_x[nj] += dv() * omega[]  * p.x;
      mu_y[nj] += dv() * omega[]  * p.y;
    }
  }

  for (int nj = 0; nj < 1; nj++){
    mu_x[nj] /= circ2[nj];
    mu_y[nj] /= circ2[nj];
  }

  double M20[1] = {0.}, M02[1] = {0.}, M11[1] = {0.};
  foreach(reduction(+:M20[:1]), reduction(+:M02[:1]), reduction(+:M11[:1])){
    coord p = {x,y,z};
    int nj = 0;
    if (sqrt(sq(x)+sq(y)) < L0*3/8.){
      M20[nj] += dv() * omega[] * sq(p.x - mu_x[nj]);
      M02[nj] += dv() * omega[] * sq(p.y - mu_y[nj]);
      M11[nj] += dv() * omega[] * (p.x - mu_x[nj]) * (p.y - mu_y[nj]);
    }
  }

  for (int nj = 0; nj < 1; nj++){
    M20[nj] /= circ2[nj];
    M02[nj] /= circ2[nj];
    M11[nj] /= circ2[nj];
  }

  double va[1] = {0.}, vb[1] = {0.}, vc[1] = {0.}, ve[1] = {0.};
  for (int nj = 0; nj < 1; nj++){
    va[nj] = sqrt(M20[nj] + M02[nj]);
    vb[nj] = sqrt((M20[nj] + M02[nj] + sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj]))));
    vc[nj] = sqrt((M20[nj] + M02[nj] - sqrt(4*sq(M11[nj]) + sq(M20[nj] - M02[nj]))));
    ve[nj] = sqrt(1 - sq(vc[nj])/sq(vb[nj]));
  }

  int nj = 0;
  fprintf (fp, "%.5f %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, ekin, circ, enst, area_tot, s0.max, circ2[nj], mu_x[nj], mu_y[nj], interpolate(omega, mu_x[nj], mu_y[nj]), va[nj], vb[nj], vc[nj], ve[nj]);
}

void global_quantities2(FILE * fp, scalar omega, vector u){

  vertex scalar phi[];
  foreach_vertex()
    phi[] = sqrt(sq(x)+sq(y)) - L0*3/8.;

  boundary ({phi});

  scalar c[];
  face vector s[];
  fractions (phi, c, s);

  squares ("c", spread = -1);
  save ("c2.png");

  double area = 0., utan=0.;
  FILE * fp2 = fopen("interface.asc", "w");
  foreach (reduction(+:area), reduction(+:utan)){
    if (c[] > 0 && c[] < 1.) {
      coord n = interface_normal (point, c), p;
      double alpha = plane_alpha (c[], n);
      double dl = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
      normalize (&n);

      double uvel = interpolate (u.x, x+p.x*Delta, y+p.y*Delta, z + p.z*Delta );
      double vvel = interpolate (u.y, x+p.x*Delta, y+p.y*Delta, z + p.z*Delta );
      double nvel = ( uvel * n.x + vvel * n.y);
      double tvel = ( uvel * n.y - vvel * n.x);
      utan += tvel * dl;
      area += dl;

      fprintf (fp2, "%g %g %g %g %g %g %g %g %g \n", x+p.x*Delta, y+p.y*Delta, n.x, n.y, uvel, vvel, nvel, tvel, dl);
    }
  }
  fclose(fp2);
  fprintf (fp, "%g %.15g Circulation: %.15g \n", t, area, utan);
}


event logfile (t+=0.01) {
  scalar omega[];
  vector u_inert[];
  vorticity (u, omega);

  foreach(){
    omega[] -= 2*omega_mean;
    u_inert.x[] = u.x[] + y*omega_mean;
    u_inert.y[] = u.y[] - x*omega_mean;
  }

  global_quantities( stderr, omega, u_inert);
  global_quantities2(stdout, omega, u_inert);
}
