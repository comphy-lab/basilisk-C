#include "tag.h"
void vorticity_moments(scalar omega, scalar m, FILE * fp){
  // We use [tag.h]() to identify each vortex
  int n = tag (m);

  double circulation[n], M20[n], M02[n], M11[n], va[n], vb[n], vc[n], ve[n], omax;
  coord centroids[n];

  for (int j = 0; j < n; j++){
    // Compute the circulation $\Gamma = \iint \omega dA$ and position of the
    // vortex centroids $\vec{\mu} = \iint \vec{x}\omega dA$
    double var1 = 0., var2 = 0., var3 = 0.;
    foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3)){
      coord p = {x,y,z};

      if (m[] == j+1){
        var1 += dv() * omega[];
        var2 += dv() * omega[] * p.x;
        var3 += dv() * omega[] * p.y;
      }
    }
    circulation[j] = var1;
    centroids[j].x = var2/var1;
    centroids[j].y = var3/var1;


    // For the ellipticity we require the centered vorticity moments
    // $ M_{p,q} = \iint (x - \mu_x)^p (y - \mu_y)^q \omega dA $
    var1 = 0., var2 = 0., var3 = 0.;
    foreach(reduction(+:var1), reduction(+:var2), reduction(+:var3)){
      coord p = {x,y,z};
      if (m[] == j+1){
        var1 += dv() * omega[] * sq(p.x - centroids[j].x);
        var2 += dv() * omega[] * sq(p.y - centroids[j].y);
        var3 += dv() * omega[] * (p.x - centroids[j].x) * (p.y - centroids[j].y);
      }
    }

    M20[j] = var1/circulation[j];
    M02[j] = var2/circulation[j];
    M11[j] = var3/circulation[j];

    /** The mean radius is obtained from $a=\sqrt{M_{2,0} + M_{0,2}}$, whereas
    the major and minos axes of the equivalent ellipse are given by
    $b = \sqrt{M_{2,0} + M_{0,2} + \sqrt{4M_{1,1}^2 + (M_{2,0} - M_{0,2})^2}}$
    and
    $c = \sqrt{M_{2,0} + M_{0,2} - \sqrt{4M_{1,1}^2 + (M_{2,0} - M_{0,2})^2}}$,
    respectively, while the ellipticity reads as $e = \sqrt{1 - c^2/b^2}$.
    */
    va[j] = sqrt(M20[j] + M02[j]);
    vb[j] = sqrt((M20[j] + M02[j] + sqrt(4*sq(M11[j]) + sq(M20[j] - M02[j]))));
    vc[j] = sqrt((M20[j] + M02[j] - sqrt(4*sq(M11[j]) + sq(M20[j] - M02[j]))));
    ve[j] = sqrt(1 - sq(vc[j])/sq(vb[j]));

    omax = interpolate (omega, centroids[j].x, centroids[j].y);

    if (pid() == 0) {
      @if _MPI
        MPI_Reduce (MPI_IN_PLACE, &omax, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      @endif

      fprintf (fp, "%g %d %g %g %g %g %g %g %g %g %g %g %g\n", t, j, circulation[j], centroids[j].x, centroids[j].y, M20[j], M02[j], M11[j], va[j], vb[j], vc[j], ve[j], omax);
    }
    @if _MPI
    else
      MPI_Reduce (&omax, NULL, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    @endif
  }
}
