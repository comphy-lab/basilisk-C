/**
We want to calculate the energy in the liquid and in the air

We also want to compute the viscous dissipation rate

Pay attention, this is done for an axi-symetric case (cylindrical coordinate)

We assume that:
 * the y coordinate correspond to the radial coordinate
 * the x coordinate correspond to the z coordinate
 * the velocity does not depend of $\theta$
 * there is no component of the velocity allong $\theta$
*/


int dissipation_rate (double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {

    /**
    Gradient of the velocity*/
    
    double dudr = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    double dudt = 0;
    double dudz = (u.y[1]   - u.x[-1]  )/(2.*Delta);

    double dvdr = 0;
    double dvdt = u.y[]/y;
    double dvdz = 0;

    double dwdr = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    double dwdt = 0;
    double dwdz = (u.x[1]   - u.x[-1]  )/(2.*Delta);

    /**
    Matrix of deformation*/

    double SDeformrr = dudr;
    double SDeformrt = 0.5*(dudt + dvdr);
    double SDeformrz = 0.5*(dudz + dwdr);

    double SDeformtr = SDeformrt;
    double SDeformtt = dvdt;
    double SDeformtz = 0.5*(dvdz + dwdt);

    double SDeformzr = SDeformrz;
    double SDeformzt = SDeformtz;
    double SDeformzz = dwdz; 

    double sqterm = 2.*dv()*(sq(SDeformrr) + sq(SDeformrt) + sq(SDeformrz) +
           sq(SDeformtr) + sq(SDeformtt) + sq(SDeformtz) +
           sq(SDeformzr) + sq(SDeformzt) + sq(SDeformzz)) ;
    rateWater += mu1*f[]*sqterm*dt; //water
    rateAir   += mu2*(1. - f[])*sqterm*dt; //air
  }
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}

/**
We compute the surface of the interface, in a fluid cell, in an axisymetric case*/

// double surfaceCyl(coord p1, coord p2) {
//   double s = M_PI*fabs(p2.y+p1.y)*fabs(p2.x-p1.x);
//   s+= M_PI*fabs(sq(p2.y)-sq(p1.y));
//   return s;
// }

double surfaceCyl(coord p1, coord p2) {
  double s = 0;
  if(fabs(p2.y-p1.y)>0.){
    double tanPhi = (p2.x-p1.x)/(p2.y-p1.y);
    s = M_PI*fabs(sq(p2.y)-sq(p1.y))/cos(atan(tanPhi));
  }
  else
    s =2*M_PI*fabs(p2.x-p1.x)*p1.y; // 2*pi*r*Delta(z)
  // double s = 2./3.*M_PI*fabs(p2.y*p2.x-p1.y*p1.x+p1.y*p2.x-p2.y*p1.x);
  // s+= M_PI*fabs(sq(p2.y)-sq(p1.y));
  return s;
}

double surfaceEnergy (){
  double surf = 0;
  face vector s;
  s.x.i = -1;
  foreach(reduction(+:surf)) {
    if (f[]> 1e-6 && f[]<1. - 1e-6){
      coord n = facet_normal(point, f, s);
      double alpha = plane_alpha(f[], n);
      coord segment[2];
      if (facets(n, alpha, segment) == 2) {
        coord p1 = {x + segment[0].x*Delta, y + segment[0].y*Delta};
        coord p2 = {x + segment[1].x*Delta, y + segment[1].y*Delta};
        surf += surfaceCyl(p1, p2);
      }
    }
  }
  return surf;
}


double dissWaterTot = 0;
double dissAirTot = 0;
double dissTot = 0;


event energy(i++){
  static FILE * fpwater = fopen("budgetWater.dat", "w");
  static FILE * fpair = fopen("budgetAir.dat", "w");
  static FILE * fpglobal = fopen("budgetGlobal.dat", "w");
  double ke = 0., gpe = 0.;
  double keAir = 0., gpeAir = 0.;
  foreach(reduction(+:ke) reduction(+:gpe) 
    reduction(+:keAir) reduction(+:gpeAir)) {
    double norm2 = 0.;
    foreach_dimension()
      norm2 += sq(u.x[]);
    ke += rho[]*norm2*f[]*dv();
    keAir += rho[]*norm2*(1.0-f[])*dv();
    gpe += -rho1*G.x*x*f[]*dv();
    gpeAir += -rho2*G.x*x*(1.0-f[])*dv();
  }
  double rates[2];
  dissipation_rate(rates);
  double dissWater = rates[0];
  double dissAir   = rates[1];
  double dissGlobal = dissAir+dissWater;

  double kSig = Sigma*surfaceEnergy();

  dissWaterTot += dissWater;
  dissAirTot += dissAir;
  dissTot = dissAirTot + dissWaterTot;
  double keTot = ke + keAir;
  double gpeTot = gpe + gpeAir;
    if (i == 0) {
    fprintf (fpwater, "i t ke gpe dissipation dissipation_Tot\n");
    fprintf (fpair, "i t ke gpe dissipation dissipation_Tot\n");
    fprintf(fpglobal , "i t ke gpe dissipation dissipation_Tot kSig\n");
    }
  fprintf (fpwater, "%d %g %g %g %g %g\n",
     i, t, ke, gpe, dissWater, dissWaterTot);
  fprintf (fpair, "%d %g %g %g %g %g\n",
     i, t, keAir, gpeAir, dissAir, dissAirTot);
    fprintf (fpglobal, "%d %g %g %g %g %g %g\n",
     i, t, keTot, gpeTot, dissGlobal, dissTot, kSig);
  fflush(fpwater);
  fflush(fpair);
  fflush(fpglobal);
}