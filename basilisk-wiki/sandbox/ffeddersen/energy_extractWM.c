
/**
 Iterate through dump files and extract the KE, PE, and dissipation rates in air & water
 */

#include "embed.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"

#include "tracer.h"
#include "navier-stokes/perfs.h"

#define FLOOR_BOUNDARY_EMBED y+max(0.37109375,1-max(0,0.0693*(x-30)))

#define RHOratio 0.000980392156862745
#define MUratio 0.018099999999999998

const double h0 = 1;

scalar dye[];
scalar * tracers = {dye};
scalar umag[];
scalar dissrate[];
scalar mask[];

//face vector airforce[];
int dissipation_rate (vector u, double* rates, scalar drate);

double rho1, rho2;
double mu1, mu2;

double maxv_water (scalar a) {

  double maxi = -1.0e100;

  foreach (reduction(max:maxi)){
    if (fabs(a[]*f[]) > maxi)
      maxi = fabs(a[]*f[]);
  }
  return maxi;
}

double maxv_air (scalar a) {

  double maxi = -1.0e100;

  foreach (reduction(max:maxi)){
    if (fabs(a[]*(1-f[])) > maxi)
      maxi = fabs(a[]*(1-f[]));
  }
  return maxi;
}

#define RE 40000

int main(int argc, char * argv[])
{

  
  double norm2;
  double keWater,keAir,gpeWater,gpeAir;
  double dissWater, dissAir, maxs;
  double rates[2];
  double grav = 1;  // this replaed  "g_"
  
  FILE * fpenergy = fopen("EnergyBudget.dat", "w");
  
  char  dumpname[1000];



  rho1 = 1.;
  rho2 = RHOratio;

  mu1 = 1.0/RE;
  mu2 = (1.0/RE) * MUratio;

  fprintf (fpenergy, "%% t ke_w gpe_w dissipation_w  ke_a  gpe_a  diss_a  max(|u|) \n");

  int dumpnum = 0;
  sprintf(dumpname,"shoal_dump/d%04d",dumpnum);  
  //sprintf(dumpname,"../fshoal_RE2400_BO2000/shoal_dump/d%04d",dumpnum);
  //sprintf(dumpname,"/home/kghanson/basilisk/shoaling/runs/test_2022_11_09/shoal_RE2400_BO2000/shoal_dump/d%04d",dumpnum);
  //sprintf(dumpname,"/home/kghanson/basilisk/shoaling/runs/test_2022_11_09/shoal_RE-1800_BO1000/shoal_dump/d%04d",dumpnum);
  //  sprintf(dumpname,"/home/kghanson/basilisk/shoaling/runs/test_2022_11_09/shoal_RE1200_BO1000/shoal_dump/d%04d",dumpnum);
  
  printf("First file to read: %s\n",dumpname);
  while(restore(file=dumpname))  {
    printf("\r %s: t=%f        ",dumpname,t);
    fflush(stdout);

    /*    if (t>1022.29) {
      break;
      }*/
    fraction(mask, FLOOR_BOUNDARY_EMBED);
/*    foreach() {
      if (FLOOR_BOUNDARY_EMBED > 0)
	mask[] = 1;
      else
	mask[] = 0;
    }*/
    double massAir = 0.;
    double massWater = 0.; 
    keWater = 0.0; gpeWater = 0.0;
    keAir = 0.0; gpeAir = 0.0;
    foreach(reduction(+:keWater) reduction(+:gpeWater) 
	  reduction(+:keAir) reduction(+:gpeAir) 
	  reduction(+:massAir) reduction(+:massWater) ) {
      norm2 = 0.0;
      foreach_dimension()
        norm2 += sq(u.x[]);
      keWater += rho1*norm2*f[]*dv() * mask[]; //(1.0-beach[]);
      keAir += rho2*norm2*(1.0-f[])*dv() * mask[]; //*(1.0-beach[]);
      gpeWater += rho1*grav*y*f[]*dv() * mask[]; //*(1.0-beach[]);   
      gpeAir += rho2*grav*y*(1.0-f[])*dv() * mask[]; //*(1.0-beach[]);
      massWater += rho1*f[]*dv()*mask[];
      massAir += rho2*(1.-f[])*dv()*mask[];
    umag[] = sqrt(sq(u.x[]) + sq(u.y[]));
    } 
  keWater *=0.5;
  keAir *=0.5;
//  gpeWater += 75000.;
  dissipation_rate(u, rates, dissrate);
  double dissWater = rates[0];
  double dissAir   = rates[1];

  double maxs_water = maxv_water(umag);
  double maxs_air = maxv_air(umag);   
    fprintf (fpenergy, "%g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",
	     t, keWater, gpeWater, dissWater, maxs_water, keAir, gpeAir, dissAir, maxs_air, massWater, massAir);
    fprintf (stdout, "%g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n",    
	     t, keWater, gpeWater, dissWater, maxs_water, keAir, gpeAir, dissAir, maxs_air, massWater, massAir);
  
    dumpnum++;
    sprintf(dumpname,"shoal_dump/d%04d",dumpnum);    
  //  sprintf(dumpname,"../fshoal_RE2400_BO2000/shoal_dump/d%04d",dumpnum);
    //sprintf(dumpname,"/home/kghanson/basilisk/shoaling/runs/test_2022_11_09/shoal_RE2400_BO2000/shoal_dump/d%04d",dumpnum);
    //sprintf(dumpname,"/home/kghanson/basilisk/shoaling/runs/test_2022_11_09/shoal_RE-1800_BO1000/shoal_dump/d%04d",dumpnum);
    //    sprintf(dumpname,"/home/kghanson/basilisk/shoaling/runs/test_2022_11_09/shoal_RE1200_BO1000/shoal_dump/d%04d",dumpnum);

  } 
  fclose(fpenergy);
  printf("Done.                                          \n");

}





//------------------DIAGNOSTICS---------------------//
/*Define functions for determining kinetic and potential energy*/

int dissipation_rate (vector u, double* rates, scalar drate)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1,0] - u.x[-1,0])/(2.0*Delta);
    double dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    //    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1,0]   - u.y[-1,0]  )/(2.*Delta);
    double dvdy = (u.y[0,1] - u.y[0,-1])/(2.0*Delta);
    // double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    //double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    //double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    //double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.0*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    //double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    //double SDeformyz = 0.5*(dvdz + dwdy);
    //double SDeformzx = SDeformxz;
    //double SDeformzy = SDeformyz;
    //double SDeformzz = dwdz; 
    /*    double sqterm = 2.0*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			      sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			      sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz)); */
    double sqterm = 2.0*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformyx) + sq(SDeformyy));
    double localRate = (mu1/rho1)*f[]*sqterm * mask[]; //*(1.0-beach[]);
    drate[] = localRate;
    rateWater += localRate; //water
    rateAir += (mu2/rho2)*(1.0-f[])*sqterm * mask[];  //*(1.0-beach[]); //air
  }
  //fprintf (fix, "\n");
  //fprintf (fiy, "\n");
  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}



