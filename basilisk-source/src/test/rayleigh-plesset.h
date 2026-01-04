#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct RPdata {
  double pliq, rhol, p0, gamma, sigma, R0, visc, cson;
} RPdata;

struct RK4data {
 double t0, h;
 double * y;
 int nvars; 
 void * data; 
};

static void RK4( void (*f)(double ti, double * y, void * data, double * dydt), double t0, double h, double * y, int nvars, void * data) 
{

/* INPUTS 
*/
   double k1[nvars], k2[nvars], k3[nvars], k4[nvars], ytmp[nvars];

   f(t0,y,data,k1);

   for (int j=0; j<nvars; j++) 
     ytmp[j] = y[j] + 0.5*k1[j]*h;
   f(t0+0.5*h,ytmp,data,k2);
   
   for (int j=0; j<nvars; j++) 
     ytmp[j] = y[j] + 0.5*k2[j]*h;
   f(t0+0.5*h,ytmp,data, k3);

   for (int j=0; j<nvars; j++) 
     ytmp[j] = y[j] + k3[j]*h;
   f(t0 + h,ytmp,data, k4);

   for (int j=0; j<nvars; j++)
      y[j] += (k1[j]+2*k2[j]+2*k3[j]+k4[j])/6 * h;

}

void RP (double ti, double * y, void * data, double * dydt){

//y[1] the radius R double dote
//y[2] velocity at the interface R dote
//
  struct RPdata * RPd = data;

  double pbub = RPd->p0*pow(RPd->R0/y[0],(3.*RPd->gamma));
  double pliqInt = pbub - 2.*RPd->sigma/y[0] - 4.*RPd->visc*y[1]/y[0];
  dydt[0] = y[1];
  dydt[1] = ((pliqInt - RPd->pliq) / RPd->rhol - 3.*pow(y[1],2.)/2.)/y[0];  

}

void KM (double ti, double * y, void * data, double * dydt){

//y[1] the radius R double dote
//y[2] velocity at the interface R dote
//
  struct RPdata * RPd = data;

  double pbub = RPd->p0*pow(RPd->R0/y[0],(3.*RPd->gamma));
  double pliqInt = pbub - 2.*RPd->sigma/y[0] - 4.*RPd->visc*y[1]/y[0];
  
  double k1 = 1.-y[1]/RPd->cson + 4*RPd->visc/RPd->rhol/RPd->cson/y[0];
  double Cfactor = y[1]/RPd->cson*(pow(0.5*y[1],2) - (pliqInt - RPd->pliq)/RPd->rhol - 3*RPd->gamma*pbub/RPd->rhol + 2*RPd->sigma/y[0]/RPd->rhol + 4*RPd->visc*y[1]/y[0]/RPd->rhol);
  dydt[0] = y[1];
  dydt[1] = ((pliqInt - RPd->pliq) / RPd->rhol - 3.*pow(y[1],2.)/2. + Cfactor )/y[0]/k1;  


}



void Integrate_RP (FILE * fp, double ti, double tend, void * data){

  struct RPdata * RPd = data;

  double y[2], yKM[2];
  y[0] = yKM[0] = RPd->R0;
  y[1] = yKM[1] =  0.;
 
  double dt = RPd->R0*sqrt(RPd->rhol/RPd->pliq)*1.e-3;

  for (double ti=0.; ti<=tend; ti+=dt) {
    RK4( RP,ti,dt,y,2, RPd );
    RK4( KM,ti,dt,yKM,2, RPd );
    fprintf(fp, "%g %g %g \n",ti, y[0], yKM[0]);
  }

}
