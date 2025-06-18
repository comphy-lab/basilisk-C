
#define ML 1
#define HYDRO 1



#include "grid/multigrid1D.h"
#if !ML
# include "saint-venant.h"
#else // ML
# include "./hydroF.h"
# define phi q
# if !HYDRO
#   include "layered/nh.h"
# endif
# include "layered/remap.h"
# include "layered/perfs.h"
#endif // ML


const double NU = 0.1, T0 = 10000;
const double HR = 1. [1];

/*
Functions reserved for Bagnold flow
*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - */
double Kb(){

  return I0*(tan(slope)-mu0)*pow(G*cos(slope),0.5)/(deltamu*dg);

}
double tc()
{
  double ans;
  ans = I0*tauc*pow(G/cos(slope),0.5)/(deltamu*rho*G*dg);
  return ans;

}
double Ubc(Point point)
{
  double yc = 0. [1];
  yc = HR-tc()/Kb();
  if (yc<=0){
    printf("No flow, change the rheology parameters.\n");
	return 1;
  }

  double zc = 0.;
  for (int l = - point.l; l < nl - point.l; l++) {
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;

  if(zc<yc) {

    return 2./3.*Kb()*(1.-pow(1-zc/HR,1.5))*pow(HR,1.5)-2.*tc()*pow(HR,0.5)*(1-pow(1-zc/HR,0.5));
  }
  else{

    return 2./3.*Kb()*pow(HR,1.5)-2.*tc()*pow(HR,0.5)+(4./3.)*pow(tc(),1.5)/pow(Kb(),0.5);
	 
  }
}


double Ubc2(double zc)
{
  double yc;
  yc = HR-tc()/Kb();
  if(zc<yc) {

    return 2./3.*Kb()*(1.-pow(1-zc/HR,1.5))*pow(HR,1.5)-2.*tc()*pow(HR,0.5)*(1-pow(1-zc/HR,0.5));
  }
  else{

    return 2./3.*Kb()*pow(HR,1.5)-2.*tc()*pow(HR,0.5)+(4./3.)*pow(tc(),1.5)/pow(Kb(),0.5);
	 
  }
}


double Sbc(double zc)
{
  double yc;
  yc = HR-tc(I0,tauc,deltamu,rho,dg)/Kb(I0,mu0,deltamu,dg);
  if(zc<yc) {

    return Kb()*pow(1-zc/HR,0.5)*pow(HR,.5)-tc()*pow(HR,-0.5)*pow(1-zc/HR,-0.5);
  }
  else{

    return 0;
	 
  }
}


double pressionHydro(Point point, scalar h,int layer){

	double H = 0.;
	double zc = 0.;
    for (int l = 0; l < layer; l++) {
		H+=h[0,0,l];
	}
	zc = H + 0.5*h[0,0,layer];
	return rho*G*cos(slope)*(eta[]-zb[]-zc);
}


double Nueq(Point point, scalar s, scalar h,int layer){
	double shear;
    double term1;
	double term2;
	if (layer>0&&layer<(nl-1)){ 
      shear = (s[0,0,layer+1]-s[0,0,layer-1])/(h[0,0,layer]+0.5*(h[0,0,layer+1]+h[0,0,layer-1]));}
	else if(layer==0){
	  shear = (s[0,0,layer+1]-s[0,0,layer])/(0.5*(h[0,0,layer]+h[0,0,layer+1]));
	}
	else if(layer==nl-1){
	  shear = (s[0,0,layer]-s[0,0,layer-1])/(0.5*(h[0,0,layer]+h[0,0,layer-1]));
	}

	double H = 0.;
    for (int l = 0; l < layer; l++) {
		H+=h[0,0,l];
	}
   term1 = deltamu*dg*pow(rho*pressionHydro(point,h,layer),0.5)/I0;
   term2 = (tauc+mu0*pressionHydro(point,h,layer))/(shear);
   return term1+term2;
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - */
FILE *file1;
int main()
{
  periodic(right);
 
  slope = 0.6[0];  
  L0 = 1.[1];
  G = 1. [1,-2];
  N = 16; 
  
  nu = NU;
  nl = 512; 
#if ML
#if NOMETRIC
  max_slope = 0.;
#endif
#if !HYDRO  
  NITERMIN = 2;
#endif
#endif

  I0 = 0.3;
  mu0 = 0.38;
  deltamu = 0.26;
  dg = 0.04;
  rho = 1.0;
  tauc = 0.125;

  file1 = fopen("profil.dat", "w");
  fclose(file1);

  run();
}



/**
We initialise the topography and the initial thickness of each layer *h*. */

event init (i = 0)
{
  foreach() {
    zb[] = 0.;
#if !ML
    h[] = HR - zb[];
#else
    foreach_layer()
      h[] = (HR)/nl;
#endif
 	eta[] = HR;  
  }



//Mettre une vitesse initiale dans le domaine 
 foreach (serial) {
	 foreach_layer(){
          u.x[] = Ubc(point);
	}
    }



 }

event acc(i++){
  foreach (serial) {
	 foreach_layer(){
          u.x[] = u.x[] + G*sin(slope)*dt;
     }
  }
}




/**
We check for convergence. */
double uold = 0 [1,-1];
double unew;

event logfile (t=1;t += 0.1 ) {

   unew = 0. [1,-1];
  foreach(){

    unew = unew + u.x[0,0,nl-1];

  }

  file1 = fopen ("conv.dat", "a");

  if(fabs(unew-uold)<0.1) fprintf(file1,"%g %g\n", t,fabs(unew-uold));
  printf("%lf\n",fabs(unew-uold));
  if (i > 1 && fabs(unew-uold) < 1e-5)
    return 1;
  else uold = unew;
  fclose(file1);
}

/**
## Outputs

We save profiles at regular intervals. */

event profiles (t += 5;t<=T0)
{
  foreach (serial) {
#if !ML
    double H = h[];
#else
    double H = 0.;
    foreach_layer()
      H += h[];
#endif
    fprintf (stderr, "%g %g %g\n", x, zb[] + H, zb[]);
  }
  fprintf (stderr, "\n\n");
}

/**
For the hydrostatic case, we compute a diagnostic vertical velocity
field `w`. Note that this needs to be done within this event because
it relies on the fluxes `hu` and face heights `hf`, which are only
defined temporarily in the [multilayer solver](hydro.h#update_eta). */

#if HYDRO
scalar w = {-1};

event update_eta (i++)
{
  if (w.i < 0)
    w = new scalar[nl];
  vertical_velocity (w, hu, hf);

  /**
  The layer interface values are averaged at the center of each
  layer. */
  
  foreach() {
    double wm = 0.;
    foreach_layer() {
      double w1 = w[];
      w[] = (w1 + wm)/2.;
      wm = w1;
    }
  }
}
#endif // HYDRO

/**
We also save the velocity and non-hydrostatic pressure fields. */

event output (t++)
{
  double shear = 0.0 [0,-1];
  file1 = fopen ("profil.dat", "w");
  foreach (serial) {
        double z = zb[];
#if HYDRO
    printf ("%g %g %g \n", x, z, u.x[]);
    foreach_layer() {
      //Compute the shear---
      shear = 0.0;
      if(point.l>=0&&point.l<(nl-1)){
	    shear = (u.x[0,0,1]-u.x[0,0,0])/(0.5*(h[0,0,0]+h[0,0,1]));}
	  else if(point.l==nl-1){
        shear = (u.x[0,0,0]-u.x[0,0,-1])/(0.5*(h[0,0,0]+h[0,0,-1]));}
    //- - - - - - - - - - - - - - -

      z += h[];
      fprintf (file1,"%g %g %g %g %g %g\n", x, z, u.x[],Ubc2(z),shear,Sbc(z));
    }
#elif // ML
    printf ("%g %g %g %g %g\n", x, z, u.x[], w[], phi[]);
    foreach_layer() {
      z += h[];
      printf ("%g %g %g %g %g\n", x, z, u.x[], w[], phi[]);
    }

#endif // ML
    printf ("\n");

  fprintf (file1,"\n \n");   
  }
//#if HYDRO
//  delete ({w});
//#endif
  printf ("# end = %g\n", t);
  fclose(file1);
}

/**

~~~gnuplot Velocity and shear profiles for Bingham flow
 set xlabel "y"
 set ylabel "u, shear"
 p "profil.dat" u 2:3 w p t'U computed', "profil.dat" u 2:4 w l t'Uexact', "" u 2:5 w p t"shear computed", "" u 2:6 w l t'shear exact'
~~~
*/ 
