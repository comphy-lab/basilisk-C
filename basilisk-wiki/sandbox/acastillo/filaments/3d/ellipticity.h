#include "tag.h"
#define dA (0.5 * sq(Delta))

trace
void vorticity3d(vector u, vector omega){
  vector du[], dv[], dw[]; // Would be nice to use a tensor here.
  foreach(){
    du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
    dv.x[] = (u.y[1] - u.y[-1])/(2.*Delta);
    dw.x[] = (u.z[1] - u.z[-1])/(2.*Delta);

    du.y[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    dv.y[] = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    dw.y[] = (u.z[0,1] - u.z[0,-1])/(2.*Delta);

    du.z[] = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    dv.z[] = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    dw.z[] = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
  }

  foreach(){
    omega.x[] = dw.y[] - dv.z[];
    omega.y[] = du.z[] - dw.x[];
    omega.z[] = dv.x[] - du.y[];
  }
  boundary ((scalar*){omega});
}

trace
void vorticity_moments_plane(scalar omega, scalar m, int nm, char * name, coord n, double _alpha){

  double circulation, M20, M02, M11, va, vb, vc, ve, omax;
  coord centroids={0.,0.,0.};

  for (int nj = 0; nj < nm; nj++){
    // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
    // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
    double var1 = 0., var2 = 0., var3 = 0., var4 = 0.;
    foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4)){
      double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
      if (fabs(alpha) > 0.87)
        continue;
      if (is_active(cell) && is_local(cell)) {
        coord p = {x,y,z};

        double sval;
        if (n.x == 1)
          sval = 0.5*(val(omega) + val(omega,1,0,0));
        else if (n.y == 1)
          sval = 0.5*(val(omega) + val(omega,0,1,0));
        else
          sval = 0.5*(val(omega) + val(omega,0,0,1));

        if (m[] == nj+1){
          var1 += dA * sval;
          var2 += dA * sval * p.x;
          var3 += dA * sval * p.y;
          var4 += dA * sval * p.z;
        }
      }
    }
    circulation = var1;


    if (circulation != 0) {
      centroids = (coord) {var2/var1, var3/var1, var4/var1};

      // For the ellipticity we require the centered vorticity moments
      // $ M_{p,q} = \iint (x - \mu_x)^p (y - \mu_y)^q \omega dA $

      var1 = 0., var2 = 0., var3 = 0., var4 = 0.;
      double var5 = 0., var6 = 0.;

      foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4), reduction(+:var5), reduction(+:var6)){
        double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
        if (fabs(alpha) > 0.87)
          continue;
        if (is_active(cell) && is_local(cell)) {
          coord p = {x,y,z};

          double sval;
          if (n.x == 1)
            sval = 0.5*(val(omega) + val(omega,1,0,0));
          else if (n.y == 1)
            sval = 0.5*(val(omega) + val(omega,0,1,0));
          else
            sval = 0.5*(val(omega) + val(omega,0,0,1));

          if (m[] == nj+1){
            var1 += dA * sval * sq(p.x - centroids.x);
            var2 += dA * sval * sq(p.y - centroids.y);
            var3 += dA * sval * sq(p.z - centroids.z);
            var4 += dA * sval * (p.x - centroids.x) * (p.y - centroids.y);
            var5 += dA * sval * (p.y - centroids.y) * (p.z - centroids.z);
            var6 += dA * sval * (p.z - centroids.z) * (p.x - centroids.x);
          }
        }
      }

      if (n.x == 1){
        M20 = var2/circulation;
        M02 = var3/circulation;
        M11 = var5/circulation;
      }
      else if (n.y == 1){
        M20 = var1/circulation;
        M02 = var3/circulation;
        M11 = var6/circulation;
      }
      else{
        M20 = var1/circulation;
        M02 = var2/circulation;
        M11 = var4/circulation;
      }

      /** The mean radius is obtained from $a=\sqrt{M_{2,0} + M_{0,2}}$, whereas
      the major and minos axes of the equivalent ellipse are given by
      $b = \sqrt{M_{2,0} + M_{0,2} + \sqrt{4M_{1,1}^2 + (M_{2,0} - M_{0,2})^2}}$
      and
      $c = \sqrt{M_{2,0} + M_{0,2} - \sqrt{4M_{1,1}^2 + (M_{2,0} - M_{0,2})^2}}$,
      respectively, while the ellipticity reads as $e = \sqrt{1 - c^2/b^2}$.
      */

      va = sqrt(M20 + M02);
      vb = sqrt((M20 + M02 + sqrt(4*sq(M11) + sq(M20 - M02))));
      vc = sqrt((M20 + M02 - sqrt(4*sq(M11) + sq(M20 - M02))));
      ve = sqrt(1 - sq(vc)/sq(vb));
    }
    //fprintf (stdout, "%.15g %d %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, pid(), circulation, centroids.x, centroids.y, centroids.z, M20, M02, M11, va, vb, vc, ve);

    omax = interpolate (omega, centroids.x, centroids.y, centroids.z);
    if (pid() == 0){
      @if _MPI
        MPI_Reduce (MPI_IN_PLACE, &omax, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      @endif
    }
    @if _MPI
    else
      MPI_Reduce (&omax, NULL, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    @endif

    FILE * fp;
    if ((pid()==0) & (circulation != 0)) {
      fp = fopen(name, "a");
      fprintf (fp, "%.4f %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, circulation, centroids.x, centroids.y, centroids.z, M20, M02, M11, va, vb, vc, ve, omax);
      fclose(fp);
    }
  }
}























trace
void vorticity_moments2_plane(vector omega, scalar m, int nm, char * name, coord n, double _alpha){

  double circulation, M20, M02, M11, va, vb, vc, ve, omax;
  coord centroids={0.,0.,0.}, tvec={0.,0.,0.};

  for (int nj = 0; nj < nm; nj++){
    // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
    // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
    double var1 = 0., var2 = 0., var3 = 0., var4 = 0.;
    foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4)){
      double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
      if (fabs(alpha) > 0.87)
        continue;
      if (is_active(cell) && is_local(cell)) {
        coord sval;
        if (n.x == 1){
          sval.x = val(omega.x) + val(omega.x,1,0,0);
          sval.y = val(omega.y) + val(omega.y,1,0,0);
          sval.z = val(omega.z) + val(omega.z,1,0,0);
        }
        else if (n.y == 1){
          sval.x = val(omega.x) + val(omega.x,0,1,0);
          sval.y = val(omega.y) + val(omega.y,0,1,0);
          sval.z = val(omega.z) + val(omega.z,0,1,0);
        }
        else{
          sval.x = val(omega.x) + val(omega.x,0,0,1);
          sval.y = val(omega.y) + val(omega.y,0,0,1);
          sval.z = val(omega.z) + val(omega.z,0,0,1);
        }

        if (m[] == nj+1){
          var1 += dA * 0.5 * sval.x;
          var2 += dA * 0.5 * sval.y;
          var3 += dA * 0.5 * sval.z;
          var4 += dA * 0.25 * sqrt(sq(sval.x) + sq(sval.y) + sq(sval.z));
        }
      }
    }

    circulation = var4;

    double tol = 1e-14;
    if (circulation != 0) {
      tvec=(coord){var1, var2, var3};
      normalize(&tvec);

      double var1 = 0., var2 = 0., var3 = 0., var4 = 0.;
      foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4)){
        double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
        if (fabs(alpha) > 0.87)
          continue;
        if (is_active(cell) && is_local(cell)) {
          coord p = {x,y,z};

          double sval;
          if (n.x == 1)
            sval = 0.5*(val(omega.x) + val(omega.x,1,0,0)) * tvec.x + 0.5*(val(omega.y) + val(omega.y,1,0,0)) * tvec.y + 0.5*(val(omega.z) + val(omega.z,1,0,0)) * tvec.z;
          else if (n.y == 1)
            sval = 0.5*(val(omega.x) + val(omega.x,0,1,0)) * tvec.x + 0.5*(val(omega.y) + val(omega.y,0,1,0)) * tvec.y + 0.5*(val(omega.z) + val(omega.z,0,1,0)) * tvec.z;
          else
            sval = 0.5*(val(omega.x) + val(omega.x,0,0,1)) * tvec.x + 0.5*(val(omega.y) + val(omega.y,0,0,1)) * tvec.y + 0.5*(val(omega.z) + val(omega.z,0,0,1)) * tvec.z;

          if (m[] == nj+1){
            var1 += dA * sval;
            var2 += dA * sval * p.x;
            var3 += dA * sval * p.y;
            var4 += dA * sval * p.z;
          }
        }
      }

      centroids = (coord) {var2/var1, var3/var1, var4/var1};

      // For the ellipticity we require the centered vorticity moments
      // $ M_{p,q} = \iint (x - \mu_x)^p (y - \mu_y)^q \omega dA $

      var1 = 0., var2 = 0., var3 = 0., var4 = 0.;
      double var5 = 0., var6 = 0.;

      foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3), reduction(+:var4), reduction(+:var5), reduction(+:var6)){
        double alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;
        if (fabs(alpha) > 0.87)
          continue;
        if (is_active(cell) && is_local(cell)) {
          coord p = {x,y,z};

          double sval;
          if (n.x == 1)
            sval = 0.5*(val(omega.x) + val(omega.x,1,0,0)) * tvec.x + 0.5*(val(omega.y) + val(omega.y,1,0,0)) * tvec.y + 0.5*(val(omega.z) + val(omega.z,1,0,0)) * tvec.z;
          else if (n.y == 1)
            sval = 0.5*(val(omega.x) + val(omega.x,0,1,0)) * tvec.x + 0.5*(val(omega.y) + val(omega.y,0,1,0)) * tvec.y + 0.5*(val(omega.z) + val(omega.z,0,1,0)) * tvec.z;
          else
            sval = 0.5*(val(omega.x) + val(omega.x,0,0,1)) * tvec.x + 0.5*(val(omega.y) + val(omega.y,0,0,1)) * tvec.y + 0.5*(val(omega.z) + val(omega.z,0,0,1)) * tvec.z;

          if (m[] == nj+1){
            var1 += dA * sval * sq(p.x - centroids.x);
            var2 += dA * sval * sq(p.y - centroids.y);
            var3 += dA * sval * sq(p.z - centroids.z);
            var4 += dA * sval * (p.x - centroids.x) * (p.y - centroids.y);
            var5 += dA * sval * (p.y - centroids.y) * (p.z - centroids.z);
            var6 += dA * sval * (p.z - centroids.z) * (p.x - centroids.x);
          }
        }
      }

      if (n.x == 1){
        M20 = var2/circulation;
        M02 = var3/circulation;
        M11 = var5/circulation;
      }
      else if (n.y == 1){
        M20 = var1/circulation;
        M02 = var3/circulation;
        M11 = var6/circulation;
      }
      else{
        M20 = var1/circulation;
        M02 = var2/circulation;
        M11 = var4/circulation;
      }

      /** The mean radius is obtained from $a=\sqrt{M_{2,0} + M_{0,2}}$, whereas
      the major and minos axes of the equivalent ellipse are given by
      $b = \sqrt{M_{2,0} + M_{0,2} + \sqrt{4M_{1,1}^2 + (M_{2,0} - M_{0,2})^2}}$
      and
      $c = \sqrt{M_{2,0} + M_{0,2} - \sqrt{4M_{1,1}^2 + (M_{2,0} - M_{0,2})^2}}$,
      respectively, while the ellipticity reads as $e = \sqrt{1 - c^2/b^2}$.
      */

      va = sqrt(M20 + M02);
      vb = sqrt((M20 + M02 + sqrt(4*sq(M11) + sq(M20 - M02))));
      vc = sqrt((M20 + M02 - sqrt(4*sq(M11) + sq(M20 - M02))));
      ve = sqrt(1 - sq(vc)/sq(vb));
      //fprintf (stdout, "%.15g %d %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, pid(), circulation, centroids.x, centroids.y, centroids.z, M20, M02, M11, va, vb, vc, ve);
    }

    var1 = interpolate (omega.x, centroids.x, centroids.y, centroids.z);
    var2 = interpolate (omega.y, centroids.x, centroids.y, centroids.z);
    var3 = interpolate (omega.z, centroids.x, centroids.y, centroids.z);
    omax = (var1 * tvec.x) + (var2 * tvec.y) + (var3 * tvec.z);
    if (pid() == 0){
      @if _MPI
        MPI_Reduce (MPI_IN_PLACE, &omax, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      @endif
    }
    @if _MPI
    else
      MPI_Reduce (&omax, NULL, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    @endif

    FILE * fp;
    if ((pid()==0) & (circulation != 0)) {
      fp = fopen(name, "a");
      fprintf (fp, "%.4f %d %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g \n", t, nj, circulation, centroids.x, centroids.y, centroids.z, M20, M02, M11, va, vb, vc, ve, omax);
      fclose(fp);
    }
  }
}
