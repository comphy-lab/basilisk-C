#include "navier-stokes/centered.h"
#include "view.h"
#include "fractions.h"
#include "tag.h"

static double omega_mean;

double lamb_oseen (double x, double y, double a, double b){
  return 1.0/(pi*sq(a)) * exp( -(sq(x-b)+sq(y))/sq(a) );
}

int maxlevel = 10;
int minlevel = 5;

int main(){
  L0 = 16 ;
  X0 = Y0 = -L0/2 ;
  omega_mean = -0.5*2.0/sq(L0);
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

event init (t = 0){
  scalar c[], b[];

  double bsum = 0.;
  foreach(reduction(+:bsum)) {
    b[] = lamb_oseen(x, y, 0.2, 1.0) + lamb_oseen(x, y, 0.2, -1.0);
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

  FILE * fp = fopen("globals.asc", "w"); fclose(fp);
  fp = fopen("closed_paths.asc", "w"); fclose(fp);

  event("logfile");
}

#define K0() (0.)
#define F0() (2.0 * omega_mean)
#define alpha_H 1.0
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

event output (t = 5)
  dump();




struct VorticityCentroids {
  int nm;
  double * circulation;
  double * centroids_x;
  double * centroids_y;
  double * centroids_z;
  double * area_total;
  double * omax;
};

struct VorticityMoments {
  double * M20;
  double * M02;
  double * M11;
  double * va;
  double * vb;
  double * vc;
  double * ve;
};

#include "PointTriangle.h"
#define sortme(x) a = x[i]; x[i] = x[j]; x[j] = a;
void moments(scalar omega, scalar m, struct VorticityCentroids pcen, struct VorticityMoments pmom){
  int nm = pcen.nm;
  double var0[nm], var1[nm], var2[nm], var3[nm], var4[nm];
  for (int nj = 0; nj < nm; nj++){
    pcen.centroids_x[nj] = 0.0;
    pcen.centroids_y[nj] = 0.0;
    var0[nj] = 0.0;
    var1[nj] = 0.0;
    var2[nj] = 0.0;
    var3[nj] = 0.0;
    var4[nj] = 0.0;
  }

  // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
  // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
  foreach(reduction(+:var0[:nm]), reduction(+:var1[:nm]), reduction(+:var2[:nm]), reduction(+:var3[:nm]), reduction(+:var4[:nm])){
    coord p = {x,y,z};
    for (int nj = 0; nj < nm; nj++){
      if (m[] == nj){
        var0[nj] += dv() ;
        var1[nj] += dv() * omega[];
        var2[nj] += dv() * omega[] * p.x;
        var3[nj] += dv() * omega[] * p.y;
        var4[nj] += dv() * omega[] * p.z;
      }
    }
  }

  for (int nj = 0; nj < nm; nj++){
    pcen.area_total[nj]  = var0[nj];
    pcen.circulation[nj] = var1[nj];
    if (pcen.circulation[nj] != 0) {
      pcen.centroids_x[nj] = var2[nj]/var1[nj];
      pcen.centroids_y[nj] = var3[nj]/var1[nj];
      pcen.centroids_z[nj] = var4[nj]/var1[nj];
      pcen.omax[nj] = interpolate(omega, pcen.centroids_x[nj], pcen.centroids_y[nj]);
    }
  }

  for (int nj = 0; nj < nm; nj++){
    var0[nj] = 0.0;
    var1[nj] = 0.0;
    var2[nj] = 0.0;
  }
  foreach(reduction(+:var0[:nm]), reduction(+:var1[:nm]), reduction(+:var2[:nm])){
    coord p = {x,y,z};
    for (int nj = 0; nj < nm; nj++){
      if (m[] == nj){
        var0[nj] += dv() * omega[] * sq(p.x - pcen.centroids_x[nj]);
        var1[nj] += dv() * omega[] * sq(p.y - pcen.centroids_y[nj]);
        var2[nj] += dv() * omega[] * (p.x - pcen.centroids_x[nj]) * (p.y - pcen.centroids_y[nj]);
      }
    }
  }

  for (int nj = 0; nj < nm; nj++){
    pmom.M20[nj] = var0[nj]/pcen.circulation[nj];
    pmom.M02[nj] = var1[nj]/pcen.circulation[nj];
    pmom.M11[nj] = var2[nj]/pcen.circulation[nj];

    double mplus  = pmom.M20[nj] + pmom.M02[nj];
    double mminus = pmom.M20[nj] - pmom.M02[nj];

    pmom.va[nj] = sqrt(mplus);
    pmom.vb[nj] = sqrt((mplus + sqrt(4*sq(pmom.M11[nj]) + sq(mminus))));
    pmom.vc[nj] = sqrt((mplus - sqrt(4*sq(pmom.M11[nj]) + sq(mminus))));
    pmom.ve[nj] = sqrt(1 - sq(pmom.vc[nj])/sq(pmom.vb[nj]));
  }

  for (int i = 0; i < nm; i++){
    for (int j = i + 1; j < nm; ++j){
      if (abs(pcen.circulation[i]) < abs(pcen.circulation[j])){
        double a;
        sortme(pcen.circulation);
        sortme(pcen.centroids_x);
        sortme(pcen.centroids_y);
        sortme(pcen.centroids_z);
        sortme(pcen.area_total);
        sortme(pcen.omax);
        sortme(pmom.M20);
        sortme(pmom.M02);
        sortme(pmom.M11);
        sortme(pmom.va);
        sortme(pmom.vb);
        sortme(pmom.vc);
        sortme(pmom.ve);
      }
    }
  }
}

void global_quantities(FILE * fp, scalar omega, vector u){

  double ekin = 0.0, circ = 0.0, enst = 0.0, area_tot=0.0, xmax=1e30, ymax=1e30;
  stats s0 = statsf (omega);

  foreach(reduction(+:ekin), reduction(+:circ), reduction(+:enst), reduction(+:area_tot)){
    ekin += dv() * (sq(u.x[]) + sq(u.y[]));
    circ += dv() * omega[];
    enst += dv() * sq(omega[]);
    area_tot += dv();

    if (omega[] == s0.max) {
      xmax = x;
      ymax = y;
    }
  }
  fprintf (fp, "%.5f %.15f %.15f %.15f %.15f %.15f %.15f %.15f \n", t, ekin, circ, enst, area_tot, s0.max, xmax, ymax);
}

void global_quantities2(FILE * fp, scalar omega, vector u, double xc, double yc, double radius, int nj){

  vertex scalar phi[];
  foreach_vertex()
    phi[] = sqrt(sq(x-xc)+sq(y-yc)) - radius;
  boundary ({phi});

  scalar c[];
  face vector s[];
  fractions (phi, c, s);
  squares ("c", spread = -1);
  save ("c2.png");

  double area = 0., utang=0., unorm=0.;
  foreach (reduction(+:area), reduction(+:utang), reduction(+:unorm)){
    if (c[] > 0 && c[] < 1.) {
      coord n = interface_normal (point, c), p;
      double alpha = plane_alpha (c[], n);
      double dl = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
      normalize (&n);

      double uvel = interpolate (u.x, x+p.x*Delta, y+p.y*Delta, z + p.z*Delta );
      double vvel = interpolate (u.y, x+p.x*Delta, y+p.y*Delta, z + p.z*Delta );
      double tvel = ( uvel * n.y - vvel * n.x);
      double nvel = ( uvel * n.x + vvel * n.y);
      utang += tvel * dl;
      unorm += nvel * dl;
      area += dl;
    }
  }
  fprintf (fp, "%.5f %d %.15f %.15f %.15f \n", t, nj, area, utang, unorm);
}


event logfile (t+=0.025) {
  scalar omega[];
  vector u_inert[];
  vorticity (u, omega);

  foreach(){
    omega[] -= 2*omega_mean;
    u_inert.x[] = u.x[] + y*omega_mean;
    u_inert.y[] = u.y[] - x*omega_mean;
  }
  stats s = statsf(omega);

  FILE * fp = fopen("globals.asc", "a");
  global_quantities( fp, omega, u_inert);
  fclose(fp);

  scalar m[];
  foreach()
    m[] = omega[] > s.stddev/10;
  int nm = tag (m)+1;
  squares ("m", spread = -1);
  save ("m.mp4");

  {
    double area[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], M20[nm], M02[nm], M11[nm], va[nm], vb[nm], vc[nm], ve[nm], omax[nm];
    moments(omega, m, (struct VorticityCentroids){nm, circulation, centroids_x, centroids_y, centroids_z, area, omax}, (struct VorticityMoments) {M20, M02, M11, va, vb, vc, ve});
    for (int nj = 0; nj < nm; nj++){
      fprintf (stdout, "%.5f %d %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f \n", t, nj, area[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va[nj], vb[nj], vc[nj], ve[nj], omax[nj]);
    }

    nm = 2;
    coord c1 = {centroids_x[1],centroids_y[1],0};
    coord c2 = {centroids_x[2],centroids_y[2],0};
    foreach(){
      coord p = {x,y,z};
      m[] = vecdist2(p, c1) > vecdist2(p, c2) ? 0 : 1;
    }
    squares ("m", spread = -1);
    save ("m2.mp4");
  }
  {
    double area[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], M20[nm], M02[nm], M11[nm], va[nm], vb[nm], vc[nm], ve[nm], omax[nm];
    moments(omega, m, (struct VorticityCentroids){nm, circulation, centroids_x, centroids_y, centroids_z, area, omax}, (struct VorticityMoments) {M20, M02, M11, va, vb, vc, ve});
    for (int nj = 0; nj < nm; nj++){
      fprintf (stderr, "%.5f %d %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f \n", t, nj, area[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va[nj], vb[nj], vc[nj], ve[nj], omax[nj]);

      FILE * fp = fopen("closed_paths.asc", "a");
      global_quantities2(fp, omega, u_inert, centroids_x[nj], centroids_y[nj], 1.0, nj);
      fclose(fp);
    }
  }
}
