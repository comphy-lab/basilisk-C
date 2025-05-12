#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "lambda2.h"
#include "view.h"
#include "PointTriangle.h"

#define sortme(x) a = x[i]; x[i] = x[j]; x[j] = a;
#define RAD (sqrt(sq(x) + sq(y)))
#define MINLEVEL 5

int n_seg;
double Gamma= 1.0, as;

#include "../filaments.h"
#include "../filaments.c"
#include "../../output_fields/output_matrix_normal.h"

void helical_pair_filaments(vector omega, char * name){

  double Rs, Hs, as, circulation, axialvel;

  // Step 1. Define the helix as the spacecurve c
  FILE * fp = fopen(name, "r");
  fscanf(fp, "%d",  &n_seg);
  fscanf(fp, "%lf", &Rs);
  fscanf(fp, "%lf", &as);
  fscanf(fp, "%lf", &Hs);
  fscanf(fp, "%lf", &circulation);
  fscanf(fp, "%lf", &axialvel);

  double t0[n_seg];
  coord c[n_seg], dc[n_seg], d2c[n_seg];
  for (int i = 0; i < n_seg; i++){
    fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &t0[i], &c[i].x, &c[i].y, &c[i].z, &dc[i].x, &dc[i].y, &dc[i].z, &d2c[i].x, &d2c[i].y, &d2c[i].z);
    c[i].z -= L0/2;
  }
  fclose (fp);

  view (camera="iso", fov=40);
  box();
  cells( alpha=-L0/2 );
  cells( n = {1,0,0} );
  draw_space_curve(n_seg, c);
  char subname[80];
  sprintf(subname, "%s.png", name);
  save (subname);

  // Step 2. Compute a Frenet-Serret Frame
  coord tvec[n_seg], nvec[n_seg], bvec[n_seg];
  for (int i = 0; i < n_seg; i++){
    foreach_dimension(){
      tvec[i].x =  dc[i].x/sqrt(vecdot( dc[i],  dc[i]));
      nvec[i].x = d2c[i].x/sqrt(vecdot(d2c[i], d2c[i]));
    }
    bvec[i] = vecdotproduct(tvec[i], nvec[i]);
  }

  // Steps 3. Initialize the vorticity field
  foreach(){
    coord pcar = {x,y,z};
    coord val_omega = get_vorticity_filament(pcar, n_seg, as, t0, c, tvec, nvec, bvec, 0);
    foreach_dimension()
      omega.x[] += val_omega.x;
  }
}

coord omega_mean = {0.,0.,0.};
#include "moments3d.h"

int main(){
  L0 = 8*4.5;
  X0 = Y0 = Z0 = -L0/2;
  omega_mean = (coord) {0.,0.,(-0.5*2.0/sq(L0))};

  init_grid (1 << MINLEVEL);

  //TOLERANCE = 1e-6;
  //double reynolds= 16100;
  double reynolds= 5000;
  const face vector muc[] = {1./reynolds,1./reynolds,1./reynolds};
  mu = muc;

  periodic(top);
  periodic(left);
  periodic(front);
  run();
}

event adapt (i++){
  adapt_wavelet ((scalar*){u}, (double[]){5e-3, 5e-3, 5e-3}, MAXLEVEL, MINLEVEL);
}

#ifdef _INIT
event output (i = 2){
  printf( "Done!... \n" );
}
#else
event output (t = _TMAX){
  printf( "Done!... \n" );

  char name[180];
  sprintf (name, "./dump-%g", t);
  dump (name, list = (scalar *){u});
  printf( "Saving!... \n" );
}
#endif

event init (t = 0){

  FILE * fp;
  char header[180];
  sprintf(header, "[1]t\t [2]tag\t [3]Atot\t [4]Gamma\t [5]mu_x\t [6]mu_y\t [7]mu_z\t [8]M20\t [9]M02\t [10]M11\t [11]a\t [12]b\t [13]c\t [14]e\t [15]maxvor \n");
  fp = fopen("vortex_x0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_y0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_x0.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("vortex_y0.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("globals.asc", "w"); fclose (fp);
  fp = fopen("closed_paths.asc", "w"); fclose (fp);

#ifdef _INIT
  printf( "Creating initial conditions... \n" );

  printf( "Initial conditions from filamentary solutions \n" );

  // Refine around the region of interest
  {
    char name[180];
    double Hs1, Rs1;
    double Hs2, Rs2;

    sprintf (name, "../base_hext.asc");
    fp = fopen(name, "r");
    fscanf(fp, "%d",  &n_seg);
    fscanf(fp, "%lf", &Rs1);
    fscanf(fp, "%lf", &as);
    fscanf(fp, "%lf", &Hs1);
    fclose (fp);

    sprintf (name, "../base_hint.asc");
    fp = fopen(name, "r");
    fscanf(fp, "%d",  &n_seg);
    fscanf(fp, "%lf", &Rs2);
    fscanf(fp, "%lf", &as);
    fscanf(fp, "%lf", &Hs2);
    fclose (fp);

    if (MAXLEVEL < 9){
      refine  (RAD < Rs1 + 32*as && RAD > Rs2 - 32*as && level < MAXLEVEL - 3);
      refine  (RAD < Rs1 + 16*as && RAD > Rs2 - 16*as && level < MAXLEVEL - 2);
      refine  (RAD < Rs1 +  8*as && RAD > Rs2 -  8*as && level < MAXLEVEL - 1);
      refine  (RAD < Rs1 +  4*as && RAD > Rs2 -  4*as && level < MAXLEVEL);
    }
    else {
      sprintf (name, "./level%d/dump0", MAXLEVEL-1);
      fprintf (stderr, "%s", name);
      restore(name);
      adapt_wavelet ((scalar*){u}, (double[]){5e-3, 5e-3, 5e-3}, MAXLEVEL, MINLEVEL);
    }
  }

  vector c[], b[];
  foreach(){
    foreach_dimension(){
      b.x[] = 0.;
      c.x[] = 0.;
    }
  }
  helical_pair_filaments(b, "../base_hext.asc");
  //helical_pair_filaments(b, "../base_hint.asc");

  double bsum_x = 0., bsum_y = 0., bsum_z = 0.;
  foreach(reduction(+:bsum_x), reduction(+:bsum_y), reduction(+:bsum_z))
    foreach_dimension()
      bsum_x += b.x[]*dv();


  fprintf (stderr, "bsum (before): %g %g %g\n", bsum_x, bsum_y, bsum_z);
  fprintf (stderr, "omega_mean (before): %g %g %g\n", omega_mean.x, omega_mean.y, omega_mean.z);

  omega_mean = (coord) {-0.5*bsum_x/pow(L0,3), -0.5*bsum_y/pow(L0,3), -0.5*bsum_z/pow(L0,3)};

  fprintf (stderr, "omega_mean (after): %g %g %g\n", omega_mean.x, omega_mean.y, omega_mean.z);

  bsum_x = 0., bsum_y = 0., bsum_z = 0.;
  foreach(reduction(+:bsum_x), reduction(+:bsum_y), reduction(+:bsum_z)) {
    double r = sqrt(sq(x) + sq(y));
    foreach_dimension(){
      c.x[] -= sq(r)*omega_mean.x/2.0;
      b.x[] += 2.0*omega_mean.x;
      bsum_x += b.x[]*dv();
    }
  }
  boundary ((scalar *){c,b});

  fprintf (stderr, "bsum (after): %g %g %g\n", bsum_x, bsum_y, bsum_z);


  /** The Poisson equation is solved. */
  foreach_dimension(){
    timer t = timer_start();
    mgstats s = poisson (c.x, b.x, tolerance=1e-15);
    double dt = timer_elapsed (t);
    printf ("%d %g %d %d\n", N, dt, s.i, s.nrelax);
  }
  boundary ((scalar *){c});

  foreach(){
    u.x[] = -((c.z[0,1,0] - c.z[0,-1,0]) - (c.y[0,0,1] - c.y[0,0,-1]))/(2.*Delta);
    u.y[] = -((c.x[0,0,1] - c.x[0,0,-1]) - (c.z[1,0,0] - c.z[-1,0,0]))/(2.*Delta);
    u.z[] = -((c.y[1,0,0] - c.y[-1,0,0]) - (c.x[0,1,0] - c.x[0,-1,0]))/(2.*Delta);
  }
  boundary ((scalar *){u});

  vector omega[];
  vorticity3d (u, omega);

  dump ("dump0", list = (scalar *){u});

  /**
  The solution is displayed using bview.
  */


  squares ("c.z", spread = -1);
  save ("a.png");
  clear();

  squares ("b.z", spread = -1);
  save ("b.png");
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

  cells();
  save ("cells.png");
  clear();
  printf( "Done! \n" );

  stats s0 ;
  s0 = statsf (c.z); fprintf (stderr, "a: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (b.z); fprintf (stderr, "b: %g %g %g \n", s0.min, s0.max, s0.sum);
  s0 = statsf (omega.z); fprintf (stderr, "omega: %g %g %g \n", s0.min, s0.max, s0.sum);
  fprintf (stderr, "\n");

#else
  printf( "Loading initial conditions... \n" );
  restore("dump0");
  printf( "Done! \n" );

#endif
  event ("logfile");
  event ("slices");
  event ("movie");

  scalar l2[];
  lambda2 (u, l2);
  view(camera="iso", fov=15, width=1080, height=1080);
  squares ("l2", linear = false, alpha=-L0/2);
  isosurface ("l2", -1e-3);
  cells(alpha=-L0/2);
  save ("lambda2.png");
}

#define K0() (0.)
#define F0() omega_mean
#define alpha_H 0.5
#include "coriolis3d.h"

event movie (t += 0.01) {
  view(camera="iso", fov=15, width=1080, height=1080);
  squares ("u.x", alpha=-L0/2);
  squares ("u.x", n = {1,0,0} );
  save ("u.png");

  squares ("u.y", alpha=-L0/2);
  squares ("u.y", n = {1,0,0} );
  save ("v.png");

  squares ("u.z", alpha=-L0/2);
  squares ("u.z", n = {1,0,0} );
  save ("w.png");

  vector omega[];
  vorticity3d (u, omega);

  squares ("omega.z", alpha=-L0/2);
  squares ("omega.x", n = {1,0,0} );
  save ("omega.png");

  scalar l2[];
  lambda2 (u, l2);
  squares ("l2", alpha=-L0/2);
  isosurface ("l2", -1e-2);
  save ("l2.png");
}

event logfile (t += 0.01) {

  vector omega[], u_inert[];
  vorticity3d(u, omega);

  foreach(){
    u_inert.x[] = u.x[] + y*omega_mean.z - z*omega_mean.y;
    u_inert.y[] = u.y[] - x*omega_mean.z + z*omega_mean.x;
    u_inert.z[] = u.z[] - y*omega_mean.x + x*omega_mean.y;
    foreach_dimension()
      omega.x[] -= 2.0*omega_mean.x;
  }
  boundary((scalar*){u_inert, omega});

  // Write the global quantities for verification
  global_quantities("globals.asc", omega, u_inert);

  //global_quantities2("closed_paths.asc", u_inert, 0.0, 0.0, L0*3/8, 0);

  scalar l2[];
  lambda2 (u, l2);

  // Compute the vorticity moments on the plane x0, once using the l2 criterion
  scalar m[];
  foreach()
    m[] = ( ( ( abs(x) <= 4*Delta ) & ( y >= 0 ) ) & ( l2[]< -1e-2 ) );

  FILE * fp;
  int nm = tag (m) + 1;
  {
    double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], omax[nm];
    double M20[nm], M02[nm], M11[nm], va[nm], vb[nm], vc[nm], ve[nm];
    vorticity_moments_plane(  omega.x, m, (struct VorticityCentroids){nm, (coord){1,0,0}, 0.0, 0.0, circulation, centroids_x, centroids_y, centroids_z, Atot}, (struct VorticityMoments) {M20, M02, M11, va, vb, vc, ve, omax});
    for (int nj = 0; nj < nm; nj++){
      if (pid()==0) {
        fp = fopen("vortex_x0_l2.asc", "a");
        fprintf (fp, "%.5g %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, Atot[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va[nj], vb[nj], vc[nj], ve[nj], omax[nj]);
        fclose(fp);
      }
    }

    // view(camera="iso", fov=15, width=1080, height=1080);
    // squares ("m", spread = -1, n = {1,0,0} );
    // for (int nj = 1; nj < nm; nj++){
    //   scatter (centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    // }
    // save ("m.png");

    foreach(){
      m[] = 0;
      coord p = {x,y,z};
      for (int nj = 1; nj < nm; nj++){
        coord c1 = {centroids_x[nj],centroids_y[nj],centroids_z[nj]};
        if (vecdist2(p, c1) < sq(0.5))
          m[] = nj;
      }
    }
    // squares ("m", spread = -1,  n = {1,0,0} );
    // for (int nj = 1; nj < nm; nj++){
    //   scatter (centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    // }
    // save ("m2.png");
  }
  // then, using a region of interest around the centroids obtained from l2
  // this avoids cutting the low intensity vorticity
  {
    double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], omax[nm];
    double M20[nm], M02[nm], M11[nm], va[nm], vb[nm], vc[nm], ve[nm];

    vorticity_moments_plane(  omega.x, m, (struct VorticityCentroids){nm, (coord){1,0,0}, 0.0, 0.0, circulation, centroids_x, centroids_y, centroids_z, Atot}, (struct VorticityMoments) {M20, M02, M11, va, vb, vc, ve, omax});

    for (int nj = 0; nj < nm; nj++){
      if (pid()==0) {
        fp = fopen("vortex_x0.asc", "a");
        fprintf (fp, "%.5g %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, Atot[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va[nj], vb[nj], vc[nj], ve[nj], omax[nj]);
        fclose(fp);
      }
    }
  }

  // Compute the vorticity moments on the plane x0, once using the l2 criterion
  foreach()
    m[] = ( ( ( abs(y) <= 4*Delta ) & ( x <= 0 ) ) & ( l2[]< -1e-2 ) );

  nm = tag (m) + 1;
  {
    double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], omax[nm];
    double M20[nm], M02[nm], M11[nm], va[nm], vb[nm], vc[nm], ve[nm];
    vorticity_moments_plane(  omega.y, m, (struct VorticityCentroids){nm, (coord){0,1,0}, 0.0, 0.0, circulation, centroids_x, centroids_y, centroids_z, Atot}, (struct VorticityMoments) {M20, M02, M11, va, vb, vc, ve, omax});
    for (int nj = 0; nj < nm; nj++){
      if (pid()==0) {
        fp = fopen("vortex_y0_l2.asc", "a");
        fprintf (fp, "%.5g %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, Atot[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va[nj], vb[nj], vc[nj], ve[nj], omax[nj]);
        fclose(fp);
      }
    }


    // squares ("m", spread = -1, n = {0,1,0} );
    // for (int nj = 1; nj < nm; nj++){
    //   scatter (centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    // }
    // save ("m3.png");

    foreach(){
      m[] = 0;
      coord p = {x,y,z};
      for (int nj = 1; nj < nm; nj++){
        coord c1 = {centroids_x[nj],centroids_y[nj],centroids_z[nj]};
        if (vecdist2(p, c1) < sq(0.5))
          m[] = nj;
      }
    }
    // squares ("m", spread = -1,  n = {0,1,0} );
    // for (int nj = 1; nj < nm; nj++){
    //   scatter (centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    // }
    // save ("m4.png");
  }
  // then, using a region of interest around the centroids obtained from l2
  // this avoids cutting the low intensity vorticity
  {
    double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm], omax[nm];
    double M20[nm], M02[nm], M11[nm], va[nm], vb[nm], vc[nm], ve[nm];

    vorticity_moments_plane(  omega.y, m, (struct VorticityCentroids){nm, (coord){0,1,0}, 0.0, 0.0, circulation, centroids_x, centroids_y, centroids_z, Atot}, (struct VorticityMoments) {M20, M02, M11, va, vb, vc, ve, omax});

    for (int nj = 0; nj < nm; nj++){
      if (pid()==0) {
        fp = fopen("vortex_y0.asc", "a");
        fprintf (fp, "%.5g %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, Atot[nj], circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj], M20[nj], M02[nj], M11[nj], va[nj], vb[nj], vc[nj], ve[nj], omax[nj]);
        fclose(fp);
      }
    }
  }
}












event slices (t += 0.1) {
  static int nslices = 0;

  vector omega[], u_inert[];
  vorticity3d(u, omega);

  foreach(){
    u_inert.x[] = u.x[] + y*omega_mean.z - z*omega_mean.y;
    u_inert.y[] = u.y[] - x*omega_mean.z + z*omega_mean.x;
    u_inert.z[] = u.z[] - y*omega_mean.x + x*omega_mean.y;
    foreach_dimension()
      omega.x[] -= 2.0*omega_mean.x;
  }
  boundary ((scalar *){u_inert, omega});

  scalar l2[];
  lambda2 (u, l2);

  // Compute the vorticity centroids on the plane x0
  scalar m[];
  foreach()
    m[] = ( ( ( abs(x) <= 4*Delta ) & ( y >= 0 ) ) & ( l2[]< -1e-2 ) );

  int nm = tag (m) + 1;

  double Atot[nm], circulation[nm], centroids_x[nm], centroids_y[nm], centroids_z[nm];
  vorticity_centroids_plane(  omega.x, m, (struct VorticityCentroids){nm, (coord){1,0,0}, 0.0, 0.0, circulation, centroids_x, centroids_y, centroids_z, Atot});

  double cen_omega_x=0., cen_omega_y=0., cen_omega_z=0.;
  for (int nj = 1; nj < nm; nj++){
    coord P = {centroids_x[nj], centroids_y[nj], centroids_z[nj]};

    cen_omega_x = interpolate (omega.x, centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    cen_omega_y = interpolate (omega.y, centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    cen_omega_y = interpolate (omega.z, centroids_x[nj], centroids_y[nj], centroids_z[nj]);

    @if _MPI
      MPI_Allreduce (MPI_IN_PLACE, &cen_omega_x, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce (MPI_IN_PLACE, &cen_omega_y, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce (MPI_IN_PLACE, &cen_omega_z, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    @endif
    coord tvec={cen_omega_x, cen_omega_y, cen_omega_z};
    normalize(&tvec);

    double _alpha = vecdot(tvec, P);
    coord nvec = {(_alpha/tvec.x)-P.x, 0-P.y, 0-P.z};
    normalize(&nvec);
    coord bvec = vecdotproduct(tvec, nvec);

    char name[80];
    int nres = 256;
    double len_window = 2.0;
    if ((sqrt(sq(P.x)+sq(P.y)+sq(P.z)) + 1.1*len_window) < (L0/2)) {
      fprintf (stderr, "%.5g %d %.5g %.5g %.5g %.5g \n", t, nj, circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj]);

      // Slice in cartesian coordinates
      sprintf(name, "single_x0c_%3.3d_%3.3d.vts", nslices, nj);
      output_vts_normal((vector *) {u_inert, omega}, name, nres, len_window, P, tvec, nvec, bvec);

      // Slice in local coordinates
      sprintf(name, "single_x0h_%3.3d_%3.3d.vts", nslices, nj);
      output_vts_normal2((vector *) {u_inert, omega}, name, nres, len_window, P, tvec, nvec, bvec);

    }
    if (nj >= 10) break;
  }

  int imin[nm];
  double mindist[nm];
  for (int nj = 0; nj < nm; nj++){
    imin[nj]=0;
    mindist[nj] = 1e30;
  }

  for (int nj = 1; nj < nm; nj++){
    coord Pj = {centroids_x[nj], centroids_y[nj], centroids_z[nj]};
    for (int ni = 1; ni < nm; ni++){
      coord Pi = {centroids_x[ni], centroids_y[ni], centroids_z[ni]};
      if ((ni != nj) & (mindist[nj] > vecdist2(Pj, Pi))) {
        mindist[nj] = vecdist2(Pj, Pi);
        imin[nj] = ni;
      }
    }
  }

  for (int nj = 1; nj < nm; nj++){
    fprintf (stderr, "%.5g %d %d %.5g \n", t, nj, imin[nj], mindist[nj]);
  }

  for (int nj = 1; nj < nm; nj++){
    coord P = {0.5*(centroids_x[nj]+centroids_x[imin[nj]]), 0.5*(centroids_y[nj]+centroids_y[imin[nj]]), 0.5*(centroids_z[nj]+centroids_z[imin[nj]])};
    cen_omega_x = interpolate (omega.x, centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    cen_omega_y = interpolate (omega.y, centroids_x[nj], centroids_y[nj], centroids_z[nj]);
    cen_omega_y = interpolate (omega.z, centroids_x[nj], centroids_y[nj], centroids_z[nj]);

    @if _MPI
      MPI_Allreduce (MPI_IN_PLACE, &cen_omega_x, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce (MPI_IN_PLACE, &cen_omega_y, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce (MPI_IN_PLACE, &cen_omega_z, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    @endif
    coord tvec={cen_omega_x, cen_omega_y, cen_omega_z};
    normalize(&tvec);

    double _alpha = vecdot(tvec, P);
    coord nvec = {(_alpha/tvec.x)-P.x, 0-P.y, 0-P.z};
    normalize(&nvec);
    coord bvec = vecdotproduct(tvec, nvec);

    char name[80];
    int nres = 256;
    double len_window = 4.0;
    if ((sqrt(sq(P.x)+sq(P.y)+sq(P.z)) + 1.05*len_window) < (L0/2)) {
      fprintf (stderr, "%.5g %d %.5g %.5g %.5g %.5g \n", t, nj, circulation[nj], centroids_x[nj], centroids_y[nj], centroids_z[nj]);

      // Slice in cartesian coordinates
      sprintf(name, "paired_x0c_%3.3d_%3.3d.vts", nslices, nj);
      output_vts_normal((vector *) {u_inert, omega}, name, nres, len_window, P, tvec, nvec, bvec);

      // Slice in local coordinates
      sprintf(name, "paired_x0h_%3.3d_%3.3d.vts", nslices, nj);
      output_vts_normal2((vector *) {u_inert, omega}, name, nres, len_window, P, tvec, nvec, bvec);

    }
    if (nj >= 10) break;
  }
  nslices++;
}
