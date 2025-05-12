#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "fractions.h"

static double omega_mean;

double lamb_oseen (double x, double y, double a, double b){
  return 1.0/(pi*sq(a)) * exp( -(sq(x-b)+sq(y))/sq(a) );
}

#include "moments3d.h"

int maxlevel = 10;
int minlevel = 5;

int main(){
  L0 = 16 ;
  X0 = Y0 = Z0 = -L0/2 ;
  omega_mean = -0.5*2.0/sq(L0);
  //TOLERANCE = 1e-12 ;
  // double reynolds= 100 ;
  // const face vector muc[] = {1./reynolds,1./reynolds,1./reynolds} ;
  // mu = muc;

  periodic(top);
  periodic(left);
  periodic(front);
  init_grid (1 << minlevel);
  display_control (maxlevel, 6, 12);
  run();
}

event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){1e-5,1e-5,1e-5}, maxlevel, minlevel);
}


event init (t = 0){

  refine  ((sqrt(sq(x) + sq(y))) < 16 && level < maxlevel-3);
  refine  ((sqrt(sq(x) + sq(y))) < 8  && level < maxlevel-2);
  refine  ((sqrt(sq(x) + sq(y))) < 4  && level < maxlevel-1);
  refine  ((sqrt(sq(x) + sq(y))) < 2  && level < maxlevel);

  vector c[], b[];

  double bsum = 0.;
  foreach(reduction(+:bsum)) {
    foreach_dimension(){
      b.x[] = 0.;
      c.x[] = 0.;
    }
    b.z[] = lamb_oseen(x, y, 0.2, 1.0) + lamb_oseen(x, y, 0.2, -1.0);
    bsum += b.z[]*dv();
  }
  fprintf (stderr, "bsum (before): %g\n", bsum);

  omega_mean = -0.5*bsum/pow(L0,3);
  fprintf (stderr, "omega_mean: %g\n", omega_mean);

  bsum = 0.;
  foreach(reduction(+:bsum)) {
    double r = sqrt(sq(x) + sq(y));
    c.z[] = -sq(r)*omega_mean/2.0;
    b.z[] += 2*omega_mean;
    bsum += b.z[]*dv();
  }
  boundary ((scalar *){c,b});
  fprintf (stderr, "bsum: %g\n", bsum);

  /** The Poisson equation is solved. */
  foreach_dimension(){
    timer t = timer_start();
    mgstats s = poisson (c.x, b.x, tolerance = 1e-10, minlevel = 4);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);
  }

  foreach(){
    u.x[] = -((c.z[0,1,0] - c.z[0,-1,0]) - (c.y[0,0,1] - c.y[0,0,-1]))/(2.*Delta);
    u.y[] = -((c.x[0,0,1] - c.x[0,0,-1]) - (c.z[1,0,0] - c.z[-1,0,0]))/(2.*Delta);
    u.z[] = -((c.y[1,0,0] - c.y[-1,0,0]) - (c.x[0,1,0] - c.x[0,-1,0]))/(2.*Delta);
  }
  boundary ((scalar *){u});

  vector omega[];
  vorticity3d (u, omega);

  vector psi[];
  foreach()
    foreach_dimension()
      psi.x[] = c.x[];
  boundary ((scalar *){psi, omega});

  foreach_dimension(){
    timer t = timer_start();
    mgstats s = poisson (psi.x, omega.x, tolerance = 1e-10, minlevel = 4);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);
  }

  scalar e1[], e2[];
  foreach(){
    e1[] = c.z[]-psi.z[];
    e2[] = b.z[]-omega.z[];
  }

  /**
  The solution is displayed using bview.
  */

  view(camera="iso");
  squares ("c.z", spread = -1);
  save ("a.png");
  clear();

  squares ("b.z", spread = -1);
  save ("b.png");
  clear();

  squares ("psi.z", spread = -1);
  save ("psi.png");
  clear();

  squares ("omega.z", spread = -1);
  save ("omega.png");
  clear();

  squares ("u.x", spread = -1);
  save ("u.png");
  clear();

  squares ("u.y", spread = -1);
  save ("v.png");
  clear();

  squares ("u.z", spread = -1);
  save ("w.png");
  clear();

  squares ("e1", spread = -1);
  save ("e1.png");
  clear();

  squares ("e2", spread = -1);
  save ("e2.png");
  clear();

  cells();
  save ("cells.png");
  clear();

  dump ("dump0");

  stats s0 ;
  s0 = statsf (c.z); fprintf (stderr, "a: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (b.z); fprintf (stderr, "b: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (psi.z); fprintf (stderr, "psi: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (omega.z); fprintf (stderr, "omega: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (e1); fprintf (stderr, "e(psi): %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (e2); fprintf (stderr, "e(omega): %g %g %g \n", s0.min, s0.max, s0.sum);

  fprintf (stderr, "\n");

  event("logfile");
  event("slices");

  FILE * fp;
  char header[180];
  sprintf(header, "[1]t\t [2]tag\t [3]Atot\t [4]Gamma\t [5]mu_x\t [6]mu_y\t [7]mu_z\t [8]M20\t [9]M02\t [10]M11\t [11]a\t [12]b\t [13]c\t [14]e\t [15]maxvor \n");
  fp = fopen("vortex_z0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_z0.asc", "w"); fputs(header, fp); fclose (fp);
}


#define K0() (0.)
#define F0() (2.0 * omega_mean)
#define alpha_H 0.5
#include "coriolis3d.h"


#include "lambda2.h"
event movie (i++) {
  squares ("u.x", alpha=-L0/2);
  save ("u.mp4");

  squares ("u.y", alpha=-L0/2);
  save ("v.mp4");

  squares ("u.z", alpha=-L0/2);
  save ("w.mp4");

  vector omega[];
  vorticity3d (u, omega);

  squares ("omega.z", alpha=-L0/2);
  save ("omega.mp4");

  scalar l2[];
  lambda2 (u, l2);
  squares ("l2", alpha=-L0/2);
  isosurface ("l2", -1e-3);
  save ("l2.mp4");
}

event output (t = 2.0)
  dump();

event logfile (t += 0.01) {

  vector omega[], u_inert[];
  vorticity3d(u, omega);

  foreach(){
    u_inert.x[] = u.x[] + y*omega_mean;
    u_inert.y[] = u.y[] - x*omega_mean;
    u_inert.z[] = u.z[];
    omega.z[] -= 2.0*omega_mean;
  }

  // Write the global quantities for verification
  global_quantities(stderr, omega, u_inert);
  global_quantities2(stdout, u_inert, 0.0, 0.0, 3.0, 0);

  scalar l2[];
  lambda2 (u, l2);

  // Compute the vorticity moments on the plane x0, once using the l2 criterion
  scalar m[];
  foreach()
    m[] = l2[]< -1e-3 ;

  int nm = tag (m) + 1;
  {
    double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], omax[nm];
    double M20[nm], M02[nm], M11[nm], va[nm], vb[nm], vc[nm], ve[nm];
    vorticity_moments_plane(  omega.z, m, (struct VorticityCentroids){nm, (coord){0,0,1}, 0.0, 0.0, circulation, centroids_x, centroids_y, centroids_z, Atot}, (struct VorticityMoments) {M20, M02, M11, va, vb, vc, ve, omax});
    for (int nj = 0; nj < nm; nj++){
      if ((pid()==0) & (circulation[nj] != 0)) {
        FILE * fp = fopen("vortex_z0_l2.asc", "a");
        fprintf (fp, "%.5g %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, Atot[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va[nj], vb[nj], vc[nj], ve[nj], omax[nj]);
        fclose(fp);
      }
    }

    view(camera="iso");
    squares ("m", spread = -1, n = {0,0,1} );
    for (int nj = 1; nj < nm; nj++){
      scatter (centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    }
    save ("m.png");

    foreach(){
      m[] = 0;
      coord p = {x,y,z};
      for (int nj = 1; nj < nm; nj++){
        coord c1 = {centroids_x[nj],centroids_y[nj],z};
        if (vecdist2(p, c1) < sq(4.0))
          m[] = nj;
      }
    }
    squares ("m", spread = -1,  n = {0,0,1} );
    for (int nj = 1; nj < nm; nj++){
      scatter (centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    }
    save ("m2.png");
  }
  // then, using a region of interest around the centroids obtained from l2
  // this avoids cutting the low intensity vorticity
  {
    double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], omax[nm];
    double M20[nm], M02[nm], M11[nm], va[nm], vb[nm], vc[nm], ve[nm];

    vorticity_moments_plane(  omega.z, m, (struct VorticityCentroids){nm, (coord){0,0,1}, 0.0, 0.0, circulation, centroids_x, centroids_y, centroids_z, Atot}, (struct VorticityMoments) {M20, M02, M11, va, vb, vc, ve, omax});

    for (int nj = 0; nj < nm; nj++){
      if ((pid()==0) & (circulation[nj] != 0)) {
        FILE * fp = fopen("vortex_z0.asc", "a");
        fprintf (fp, "%.5g %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, Atot[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va[nj], vb[nj], vc[nj], ve[nj], omax[nj]);
        fclose(fp);
      }
    }
  }
}







#include "../../output_fields/output_matrix_normal.h"
event slices (t += 0.5) {
  static int nslices = 0;

  vector omega[], u_inert[];
  vorticity3d(u, omega);

  foreach(){
    u_inert.x[] = u.x[] + y*omega_mean;
    u_inert.y[] = u.y[] - x*omega_mean;
    u_inert.z[] = u.z[];
    omega.z[] -= 2.0*omega_mean;
  }
  boundary ((scalar *){u_inert, omega});

  coord tvec = {0,0,1};
  coord bvec = {1,0,0};
  coord nvec = {0,1,0};
  coord P = {0,0,0};

  char name[80];
  int nres = 256;
  double len_window = 4.0;

  sprintf(name, "local_z0a_%3.3d.vts", nslices);
  output_vts_normal((vector *) {u_inert, omega}, name, nres, len_window, P, tvec, nvec, bvec);

  sprintf(name, "local_x0b_%3.3d.vts", nslices);
  output_vts_normal2((vector *) {u_inert, omega}, name, nres, len_window, P, tvec, nvec, bvec);

  nslices++;
}
