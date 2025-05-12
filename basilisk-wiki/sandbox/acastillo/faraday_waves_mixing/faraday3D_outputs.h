
#include "Blues.h"
event movie(t += 0.1){  
  for (scalar s in {f}){
    view(camera = "top");
    box();
    draw_vof("f");
    save("top.mp4");
  }
}

int dissipation_rate (double* rates){
  double rateWater = 0.0;
  double rateAir = 0.0;
  double rateWater2 = 0.0;
  foreach (reduction (+:rateWater), reduction (+:rateAir), reduction (+:rateWater2)) {
    if (cs[] == 1) {
      double dcdx = (c[1]     - c[-1]    )/(2.*Delta);
      double dcdy = (c[0,1]   - c[0,-1]  )/(2.*Delta);
      double dcdz = (c[0,0,1] - c[0,0,-1])/(2.*Delta);
      double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
      double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
      double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
      double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
      double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
      double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
      double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
      double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
      double SDeformxx = dudx;
      double SDeformxy = 0.5*(dudy + dvdx);
      double SDeformxz = 0.5*(dudz + dwdx);
      double SDeformyx = SDeformxy;
      double SDeformyy = dvdy;
      double SDeformyz = 0.5*(dvdz + dwdy);
      double SDeformzx = SDeformxz;
      double SDeformzy = SDeformyz;
      double SDeformzz = dwdz;
      double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
            sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
            sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)) ;
      rateWater += mu1/rho[]*f[]*sqterm; //water
      rateAir   += mu2/rho[]*(1. - f[])*sqterm; //air
      rateWater2 += c.D1*f[]*dv()*(sq(dcdx) + sq(dcdy) + sq(dcdz)); //water
    }
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  rates[2] = rateWater2;
  return 0;
}

#include "../output_fields/available_potential.h"
event time_series(t += 0.1){
  
  scalar zref[];
  double zmix = reference_height(zref, f, c, 0.0, 1.0, true, _H0/2.);

  double ekin = 0., epot[3] = {0.}, bpot = 0.;
  foreach (reduction(+ : ekin), reduction(+ : epot[:3]), reduction(+: bpot)){
    epot[0] += dv() * rhov[] * x;
    epot[1] += dv() * rhov[] * y;
    epot[2] += dv() * rhov[] * z;
    bpot += dv() * f[] * rhov[] * zref[];
    foreach_dimension(){
      ekin += dv() * 0.5 * rhov[] * sq(u.x[]);
    }
  }

  double rates[3];
  dissipation_rate(rates);
  double disvWater = rates[0];
  double disvAir   = rates[1];
  double discWater = rates[2];

  if (pid() == 0){
    FILE *fp = fopen("globals.asc", "a");
    fprintf(fp, "%f %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g \n", t, force.Gn, ekin, epot[0], epot[1], epot[2], bpot, zmix, disvWater, disvAir, discWater, force.ramp);
    fclose(fp);
  }
}

event save_profiles(t += 0.2){
  
  scalar mfrac[];
  foreach (){
    mfrac[] = rho1*f[]/(f[]*(rho1 - rho2) + rho2);
  }

  scalar *list = {u, p, mfrac, c, zref};
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

  reference_height(yref, f, c, 0.0, 1.0, true, _H0/2.);
  foreach()
    yref[] *= cs[]*f[];
    
  char filename[99];
  sprintf(filename, "slice_x_t%09.5f.vtkhdf", t);
  output_vtkhdf_slice({f, p, c, rhov, yref}, {u}, filename, (coord){1, 0, 0}, 0);

  sprintf(filename, "slice_y_t%09.5f.vtkhdf", t);
  output_vtkhdf_slice({f, p, c, rhov, yref}, {u}, filename, (coord){0, 1, 0}, 0);

  sprintf(filename, "slice_z_t%09.5f.vtkhdf", t);
  output_vtkhdf_slice({f, p, c, rhov, yref}, {u}, filename, (coord){0, 0, 1}, 0);
}

event save_snapshots(t += 6.2832){
  char filename[99];
  sprintf(filename, "snapshot_t%09.5f.vtkhdf", t);
  output_vtkhdf({f, p, c, rhov, yref}, {u}, filename);
}
