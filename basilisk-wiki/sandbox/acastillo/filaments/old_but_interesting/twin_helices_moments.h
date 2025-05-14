event init (t = 0){
  char header[180];
  sprintf(header, "[1]t\t [2]tag\t [3]Gamma\t [4]mu_x\t [5]mu_y\t [6]mu_z\t [7]M20\t [8]M02\t [9]M11\t [10]a\t [11]b\t [12]c\t [13]e\t [14]maxvor \n");

  FILE * fp;
  fp = fopen("moments0_x0_om.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("moments2_z0_om.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("momentsn_x0_om.asc", "w"); fputs(header, fp); fclose (fp);
  fp = fopen("momentsn_z0_om.asc", "w"); fputs(header, fp); fclose (fp);
}

event moments (t += tsample1) {

  vector omega[];
  vorticity3d(u, omega);

  scalar l2[];
  lambda2 (u, l2);

  scalar m[];
  foreach()
    m[] = (((y>0) & (abs(x) <= 0.2)) & (l2[]<=0.0));

  int nm = tag (m);

  vorticity_moments_plane(omega.x,   m, nm, "moments0_x0_om.asc", (coord){1,0,0}, 0.);
  vorticity_moments2_plane(omega,    m, nm, "momentsn_x0_om.asc", (coord){1,0,0}, 0.);

  foreach()
    m[] = ((abs(z) <= 0.2) & (l2[]<=0.0));

  nm = tag (m);

  vorticity_moments_plane(omega.z,   m, nm, "moments2_z0_om.asc", (coord){0,0,1}, 0.);
  vorticity_moments2_plane(omega,    m, nm, "momentsn_z0_om.asc", (coord){0,0,1}, 0.);
}


#include "../output_fields/output_matrix_normal.h"

event slices (t += tsample2) {

  static int nslices = 0;

  // Create a mask around the intersection of the vortices and the plane x=0
  vector omega[];
  vorticity3d(u, omega);

  scalar l2[];
  lambda2 (u, l2);

  scalar m[];
  foreach()
    m[] = (((y>0) & (abs(x) <= 0.2)) & (l2[]<=0.0));

  int nm = tag (m);

  /** For the first vortex, we compute the first moment of the vorticity field
    \vec{C} = int \vec{omega} da

    and the corresponding centroid using the vorticity amplitude
    \vec{P} = int \vec{x} \vert\vec{omega}\vert da / int \vert\vec{omega}\vert da
  */

  coord P={0,0,0}, tvec={1,0,0}, nvec={0,1,0}, bvec={0,0,1}, n={1,0,0};

  double tol=1e-14, _alpha=0, var0, var1, var2, var3, var4, var5, var6, var7;
  for (int nj = 0; nj < nm; nj++) {
    var1 = 0., var2 = 0., var3 = 0., var4 = 0., var5 = 0., var6 = 0., var7 = 0.;
    foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4), reduction(+:var5), reduction(+:var6), reduction(+:var7)){
      double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
      if (fabs(alpha) > 0.87)
        continue;
      if (is_active(cell) && is_local(cell)) {
        coord p = {x,y,z};
        if (m[] == nj+1){
          var0 = sqrt(sq(val(omega.x)) + sq(val(omega.y)) + sq(val(omega.z))) + sqrt(sq(val(omega.x,1,0,0)) + sq(val(omega.y,1,0,0)) + sq(val(omega.z,1,0,0)));

          var1 += dA * 0.5 * (val(omega.x) + val(omega.x,1,0,0));
          var2 += dA * 0.5 * (val(omega.y) + val(omega.y,1,0,0));
          var3 += dA * 0.5 * (val(omega.z) + val(omega.z,1,0,0));
          var4 += dA * 0.5 * var0 ;
          var5 += dA * 0.5 * var0 * p.x;
          var6 += dA * 0.5 * var0 * p.y;
          var7 += dA * 0.5 * var0 * p.z;
        }
      }
    }
    if (var4 != 0) {
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

  A = tvec.x,  B = tvec.y, C = tvec.z, alpha = (tvec . P)

  */
      P =(coord){var5, var6, var7};
      tvec=(coord){var1, var2, var3};
      normalize(&tvec);

      /**
      We may also define the plane in terms of two vectors orthogonal to tvec.
      Here, nvec, corresponds to the intersection of the plane and the x-axis,
      i.e., y=z=0, while bvec = tvec x nvec
      */

      _alpha = vecdot(tvec, P);
      nvec = (coord){(_alpha/tvec.x)-P.x, 0-P.y, 0-P.z};
      normalize(&nvec);
      bvec = vecdotproduct(tvec, nvec);
      fprintf (stdout, "%.15g %.15g %.15g %.15g %.15g %.15g \n", t, P.x, P.y, P.z, sqrt(sq(P.x)+sq(P.y)+sq(P.z)), _alpha);

      /**
      Finally, we project the velocity and vorticity fields into the local frame
      and store everything (we re-use l2 to reduce the memory use)
      */

      FILE * fp ;
      char name[80];
      int nres = 256;
      double len_window = 0.5;
      if ((sqrt(sq(P.x)+sq(P.y)+sq(P.z)) + 1.1*len_window) < (L0/2)){
        vector * list = {u, omega};
        for (vector v in list){
          foreach()
            l2[] = v.x[]*tvec.x + v.y[]*tvec.y + v.z[]*tvec.z;
          boundary ((scalar *){l2});

          sprintf(name, "local_slice_t%s_%3.3d_%3.3d.bin", v.x.name, nslices, nj); fp = fopen(name, "w"); output_matrix_normal(l2, fp, nres, len_window, P, nvec, bvec); fclose (fp);

          foreach()
            l2[] = v.x[]*nvec.x + v.y[]*nvec.y + v.z[]*nvec.z;
          boundary ((scalar *){l2});

          sprintf(name, "local_slice_n%s_%3.3d_%3.3d.bin", v.x.name, nslices, nj); fp = fopen(name, "w"); output_matrix_normal(l2, fp, nres, len_window, P, nvec, bvec); fclose (fp);

          foreach()
            l2[] = v.x[]*bvec.x + v.y[]*bvec.y + v.z[]*bvec.z;
          boundary ((scalar *){l2});

          sprintf(name, "local_slice_b%s_%3.3d_%3.3d.bin", v.x.name, nslices, nj); fp = fopen(name, "w"); output_matrix_normal(l2, fp, nres, len_window, P, nvec, bvec); fclose (fp);
        }
        if (nj >= 3)
          break;
      }
    }
  }
  nslices++;
}
