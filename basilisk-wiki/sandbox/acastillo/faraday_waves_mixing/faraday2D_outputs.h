#include "Blues.h"
event movie(t += 0.1){

  char str1[99];
  for (scalar s in {rhov}){
    view(camera = "front");    
    draw_vof("f", lc = {1, 0, 0}, lw = 2.);
    sprintf(str1, "cs[] > 0. ? %s[] : nodata", s.name);
    squares(str1, linear = false, spread=-1, map = BluesBlues_r);
    save("front.mp4");
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
      double dudx = (u.x[1]   - u.x[-1]    )/(2.*Delta);
      double dudy = (u.x[0,1] - u.x[0,-1]  )/(2.*Delta);
      double dvdx = (u.y[1]   - u.y[-1]    )/(2.*Delta);
      double dvdy = (u.y[0,1] - u.y[0,-1]  )/(2.*Delta);
      double SDeformxx = dudx;
      double SDeformxy = 0.5*(dudy + dvdx);
      double SDeformyx = SDeformxy;
      double SDeformyy = dvdy;
      double sqterm = 2.*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformyx) + sq(SDeformyy)) ;
      rateWater += mu1/rho[]*f[]*sqterm;        //water
      rateAir   += mu2/rho[]*(1. - f[])*sqterm; //air
      rateWater2 += c.D1*f[]*dv()*(sq(dcdx) + sq(dcdy)); //water
    }
  }

  rates[0] = rateWater;
  rates[1] = rateAir;
  rates[2] = rateWater2;
  return 0;
}


#include "../output_fields/available_potential.h"
event time_series(t += 0.1){

  double ymix = reference_height(yref, f, c, 0.0, 1.0, true, _H0/2.);
  foreach()
    yref[] *= cs[]*f[];

  double ekin = 0., epot[2] = {0.}, bpot = 0.;
  foreach (reduction(+ : ekin), reduction(+ : epot[:2]), reduction(+: bpot)){
    epot[0] += dv() * rhov[] * x;
    epot[1] += dv() * rhov[] * y;
    bpot += dv() * f[] * rhov[] * yref[];
    foreach_dimension()	{
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
    fprintf(fp, "%f %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g %.9g \n", t, force.Gn, ekin, epot[0], epot[1], bpot, ymix, disvWater, disvAir, discWater, force.ramp);
    fclose(fp);
  }
}

event save_profiles(t += 0.2){

  scalar mfrac[];
  foreach (){
    mfrac[] = rho1*f[]/(f[]*(rho1 - rho2) + rho2);
  }

  scalar *list = {u, p, mfrac, c, yref};
  int n = (1 << MAXLEVEL);
  int m1 = (1 << MAXLEVEL)*(_H0/L0);
  profile_foreach_region2(list, rhov, filename="profiles_favre.asc", xmin=-L0/2., xmax=L0/2., hmin=-_H0/2.+_mindel, hmax=_H0/2.-_mindel, n=n, m1=m1);
  
}

#include "curvature.h"
event save_facets(t += 0.2){
  scalar kappa[];
  curvature(f, kappa);

  char fname[99];
  sprintf(fname, "interface_t%09.5f", t);
  output_facets_vtu(f, kappa, fname);
}

event save_snapshots(t += 0.7854){
  
  reference_height(yref, f, c, 0.0, 1.0, true, _H0/2.);
  foreach()
    yref[] *= cs[]*f[];

  char filename[99];	
  sprintf(filename, "snapshot_t%09.5f.vtkhdf", t);
  output_vtkhdf({f,p,c,rhov,yref}, {u}, filename);
}



