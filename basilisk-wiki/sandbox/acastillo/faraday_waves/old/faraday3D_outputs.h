
#include "Blues.h"
event movie(t += 0.05){  
  for (scalar s in{f}){
    view(camera = "top");
    box();
    draw_vof("f");
    save("top.mp4");
  }
}

event time_series(t += 0.1){
  
  double ekin = 0., epot[3] = {0.};
  foreach (reduction(+ : ekin), reduction(+ : epot[:3])){
    epot[0] += dv() * rhov[] * x;
    epot[1] += dv() * rhov[] * y;
    epot[2] += dv() * rhov[] * z;
    foreach_dimension(){
      ekin += dv() * 0.5 * rhov[] * sq(u.x[]);
    }
  }

  if (pid() == 0){
    FILE *fp = fopen("globals.asc", "a");
    fprintf(fp, "%f %.9g %.9g %.9g %.9g %.9g %.9g \n", t, force.Gn, ekin, epot[0], epot[1], epot[2], force.ramp);
    fclose(fp);
  }
}

event save_profiles(t += 0.2){
  
  scalar mfrac[];
  foreach (){
    mfrac[] = rho1*f[]/(f[]*(rho1 - rho2) + rho2);
  }

  scalar *list = {u, p, mfrac};
  int n = (1 << MAXLEVEL);
  int m1 = (1 << MAXLEVEL) / (L0 / _H0);
  int m2 = (1 << MAXLEVEL) / (L0 / _D0);  
  profile_foreach_region2(list, rhov, filename="profiles_favre.asc", xmin = -L0/2., xmax = L0/2., ymin = -_D0 / 2., ymax = _D0 / 2., hmin = -_H0 / 2. + _mindel, hmax = _H0 / 2. - _mindel, n = n, m1 = m1, m2 = m2);
}

#include "curvature.h"
event save_facets(t += 0.2){
  scalar kappa[];
  curvature(f, kappa);

  char fname[99];
  sprintf(fname, "interface_t%.5f", t);
  output_facets_vtu(f, kappa, fname);
}

event save_slices(t += 0.7854){
  char filename[99];
  sprintf(filename, "slice_x_t%09.5f.vtkhdf", t);
  output_vtkhdf_slice({f, p}, {u}, filename, (coord){1, 0, 0}, 0);

  sprintf(filename, "slice_y_t%09.5f.vtkhdf", t);
  output_vtkhdf_slice({f, p}, {u}, filename, (coord){0, 1, 0}, 0);

  sprintf(filename, "slice_z_t%09.5f.vtkhdf", t);
  output_vtkhdf_slice({f, p}, {u}, filename, (coord){0, 0, 1}, 0);
}

event save_snapshots(t += 6.2832){
  char filename[99];
  sprintf(filename, "snapshot_t%09.5f.vtkhdf", t);
  output_vtkhdf({f,p}, {u}, filename);
}
