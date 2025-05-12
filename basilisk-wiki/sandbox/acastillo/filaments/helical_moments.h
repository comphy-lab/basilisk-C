
scalar m[];

event init (t = 0){
  char header[180];
  sprintf(header, "[1]t\t [2]tag\t [3]Gamma\t [4]mu_x\t [5]mu_y\t [6]mu_z\t [7]M20\t [8]M02\t [9]M11\t [10]a\t [11]b\t [12]c\t [13]e\t [14]maxvor \n");

  FILE * fp;
  // fp = fopen("moments0_x0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  // fp = fopen("moments3_x0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  // fp = fopen("moments2_z0_l2.asc", "w"); fputs(header, fp); fclose (fp);
  // fp = fopen("moments3_z0_l2.asc", "w"); fputs(header, fp); fclose (fp);

  fp = fopen("moments0_x0_om.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("moments2_z0_om.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("momentsn_x0_om.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("momentsn_z0_om.asc", "w"); fputs(header, fp); fclose (fp);

  // fp = fopen("moments3_x0_om.asc", "w"); fputs(header, fp); fclose (fp);
  // fp = fopen("moments3_z0_om.asc", "w"); fputs(header, fp); fclose (fp);
}

event moments (t += tsample1) {

  scalar omega_mag[];
  foreach()
    omega_mag[] = sqrt(sq(omega.x[]) + sq(omega.y[]) + sq(omega.z[]));
  boundary ((scalar*){omega_mag});
  stats f = statsf (omega_mag);

  //////////////////////////////////////////////////////////////////////////////
  // Vorticity moments using the l2 contours
  // lambda2 (u, l2);
  // scalar m1[], m2[];
  // foreach(){
  //   m1[] = (((y>0) & (abs(x) <= 2*as)) & (l2[] < 0));
  //   m2[] = ((abs(z) <= 2*as) & (l2[] < 0));
  // }
  //
  // int nm1 = tag (m1);
  // int nm2 = tag (m2);
  //
  // vorticity_moments_plane(omega.x,   m1, nm1, "moments0_x0_l2.asc", (coord){1,0,0}, 0.);
  // vorticity_moments_plane(omega_mag, m1, nm1, "moments3_x0_l2.asc", (coord){1,0,0}, 0.);
  //
  // vorticity_moments_plane(omega.z,   m2, nm2, "moments2_z0_l2.asc", (coord){0,0,1}, 0.);
  // vorticity_moments_plane(omega_mag, m2, nm2, "moments3_z0_l2.asc", (coord){0,0,1}, 0.);

  //////////////////////////////////////////////////////////////////////////////
  // Vorticity moments using the vorticity amplitude

  foreach()
    m[] = (((y>0) & (abs(x) <= 2*as)) & (omega_mag[] > f.stddev/10));

  int nm = tag (m);

  vorticity_moments_plane(omega.x,   m, nm, "moments0_x0_om.asc", (coord){1,0,0}, 0.);
  vorticity_moments2_plane(omega,    m, nm, "momentsn_x0_om.asc", (coord){1,0,0}, 0.);
  //vorticity_moments_plane(omega_mag, m3, nm3, "moments3_x0_om.asc", (coord){1,0,0}, 0.);

  foreach()
    m[] = ((abs(z) <= 2*as) & (omega_mag[] > f.stddev/10));

  nm = tag (m);

  vorticity_moments_plane(omega.z,   m, nm, "moments2_z0_om.asc", (coord){0,0,1}, 0.);
  vorticity_moments2_plane(omega,    m, nm, "momentsn_z0_om.asc", (coord){0,0,1}, 0.);
  //vorticity_moments_plane(omega_mag, m4, nm4, "moments3_z0_om.asc", (coord){0,0,1}, 0.);

}








#include "../output_fields/output_matrix_normal.h"
event slices (t += tsample2) {
  static int nslices = 0;

  // Create a mask around the intersection of the vortices and the plane x=0

  scalar omega_mag[];
  foreach()
    omega_mag[] = sqrt(sq(omega.x[]) + sq(omega.y[]) + sq(omega.z[]));
  boundary ((scalar*){omega_mag});
  stats f = statsf (omega_mag);

  foreach()
    m[] = ((((y>=0) & (abs(x) <= 4*as)) & (abs(z) <= 2*Hs)) & (omega_mag[] > f.stddev/10));

  int nm = tag (m);

  /** For the first vortex, we compute the first moment of the vorticity field
    C.x = int(omega.x da)
    C.y = int(omega.y da)
    C.z = int(omega.z da)

    and the corresponding centroid using the vorticity amplitude
    P.x = int( x ||omega|| da)/int( ||omega|| da)
    P.y = int( y ||omega|| da)/int( ||omega|| da)
    P.z = int( z ||omega|| da)/int( ||omega|| da)
  */

  coord P={0,0,0}, tvec={1,0,0}, nvec={0,1,0}, bvec={0,0,1};


  double tol=1e-3;
  double _alpha=0;
  coord n={1,0,0};
  double var1, var2, var3, var4, var5, var6, var7;
  for (int nj = 0; nj < nm; nj++) {

    var1 = 0., var2 = 0., var3 = 0., var4 = 0., var5 = 0., var6 = 0., var7 = 0.;
    foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4), reduction(+:var5), reduction(+:var6), reduction(+:var7)){
      double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
      if (fabs(alpha) > 0.87)
        continue;
      if (is_active(cell) && is_local(cell)) {
        coord p = {x,y,z};

        if (m[] == nj+1){
          var1 += dA * 0.5*(val(omega.x) + val(omega.x,1,0,0));
          var2 += dA * 0.5*(val(omega.y) + val(omega.y,1,0,0));
          var3 += dA * 0.5*(val(omega.z) + val(omega.z,1,0,0));
          var4 += dA * 0.5*(val(omega_mag) + val(omega_mag,1,0,0));
          var5 += dA * 0.5*(val(omega_mag) + val(omega_mag,1,0,0)) * p.x;
          var6 += dA * 0.5*(val(omega_mag) + val(omega_mag,1,0,0)) * p.y;
          var7 += dA * 0.5*(val(omega_mag) + val(omega_mag,1,0,0)) * p.z;
        }
      }
    }
    if (abs(var4) > tol) {
      var5 /= var4;
      var6 /= var4;
      var7 /= var4;


  /**
  We assume that the vortex core passes through P and is aligned with C
  This way, we may write the plane perpendicular to the vortex core using only
  the unit vector tangent to the filament (normal to the plane) and the point
  P. Any point X along such plane satisfies

  (X - P) . tvec = 0

  alternatively written as

  Ax + By + Cz = alpha

  where:

  A = tvec.x
  B = tvec.y
  C = tvec.z
  alpha = (tvec . P)

  */
      P =(coord){var5, var6, var7};
      tvec=(coord){var1, var2, var3};
      normalize(&tvec);
      break;
    }
  }

  /**
  We may also define the plane in terms of two vectors orthogonal to tvec.
  Here, nvec, corresponds to the intersection of the plane and the x-axis,
  i.e., y=z=0, while bvec = tvec x nvec
  */

  _alpha = vecdot(tvec, P);
  nvec = (coord){(_alpha/tvec.x)-P.x, 0-P.y, 0-P.z};
  normalize(&nvec);
  bvec = vecdotproduct(tvec, nvec);

  fprintf (stdout, "%.15g %.15g %.15g %.15g \n", t, P.x, P.y, P.z);

  /**
  Finally, we project the velocity and vorticity fields into the local frame
  */
  vector u_loc[], omega_loc[];
  foreach(){
    u_loc.x[] = u.x[]*tvec.x + u.y[]*tvec.y + u.z[]*tvec.z;
    u_loc.y[] = u.x[]*nvec.x + u.y[]*nvec.y + u.z[]*nvec.z;
    u_loc.z[] = u.x[]*bvec.x + u.y[]*bvec.y + u.z[]*bvec.z;

    omega_loc.x[] = omega.x[]*tvec.x + omega.y[]*tvec.y + omega.z[]*tvec.z;
    omega_loc.y[] = omega.x[]*nvec.x + omega.y[]*nvec.y + omega.z[]*nvec.z;
    omega_loc.z[] = omega.x[]*bvec.x + omega.y[]*bvec.y + omega.z[]*bvec.z;
  }
  boundary ((scalar*){u_loc, omega_loc});

  /**
  and store everything
  */

  view (camera="iso", fov=4*L0);
  box();
  isosurface ("l2", 0);
  squares("omega_loc.x", n = tvec, alpha=_alpha);
  _alpha = vecdot(nvec, P);
  squares("omega_loc.x", n = nvec, alpha=_alpha);
  save ("helical5.mp4");

  FILE * fp ;
  char name[80];
  int nres = 512; //(1<<MAXLEVEL)*(int)((16*as)/L0);
  sprintf(name, "local_slice_l2%3.3d.bin", nslices); fp = fopen(name, "w"); output_matrix_normal(l2,           fp, nres, 8*as, P, nvec, bvec); fclose (fp);
  sprintf(name, "local_slice_ut%3.3d.bin", nslices); fp = fopen(name, "w"); output_matrix_normal(u_loc.x,      fp, nres, 8*as, P, nvec, bvec); fclose (fp);
  sprintf(name, "local_slice_un%3.3d.bin", nslices); fp = fopen(name, "w"); output_matrix_normal(u_loc.y,      fp, nres, 8*as, P, nvec, bvec); fclose (fp);
  sprintf(name, "local_slice_ub%3.3d.bin", nslices); fp = fopen(name, "w"); output_matrix_normal(u_loc.z,      fp, nres, 8*as, P, nvec, bvec); fclose (fp);
  sprintf(name, "local_slice_ot%3.3d.bin", nslices); fp = fopen(name, "w"); output_matrix_normal(omega_loc.x,  fp, nres, 8*as, P, nvec, bvec); fclose (fp);
  sprintf(name, "local_slice_on%3.3d.bin", nslices); fp = fopen(name, "w"); output_matrix_normal(omega_loc.y,  fp, nres, 8*as, P, nvec, bvec); fclose (fp);
  sprintf(name, "local_slice_ob%3.3d.bin", nslices); fp = fopen(name, "w"); output_matrix_normal(omega_loc.z,  fp, nres, 8*as, P, nvec, bvec); fclose (fp);
  nslices++;
}
