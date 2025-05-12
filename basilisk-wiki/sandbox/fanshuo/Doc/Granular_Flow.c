/**

$\mu(I)$ rheology is implemented in a "multilayer solver" (hydro.h)
velocity profils of kind "Bingham" and "Bagnold" are recovered.

More details about this solver (hydro.h) are given in [Popinet (2020)](/Bibliography#popinet2020). 

Following codes are tested on 15/05/2024 
*/


//Flow configurations
#define Poiseuille 0
#define Bingham 0
#define Bagnold 0  //dry granular flow
#define BagnoldC 1 //For cohesive case

//Remapping parameters
#define Geometric 0
#define HALF 0
#define Top  0

#include "grid/multigrid1D.h"
#include "./src/functionF.h"
#include "./src/hydroF.h"
#include "layered/implicit.h"
#include "./src/remapF.h"
#include "layered/perfs.h"




double  HR = 1., T0=5000;
FILE * file1, * file2;
int main()
{

  periodic(right); //periodic condition for the flow

  G = 1.;

  //Geometrie parameters
  L0 = 1;
  slope = 0.2;
  N = 16 ;
  nl = 512; 

  //Rheology parameters
  #if Poiseuille
    nu = 0.1; //In the case where viscosity is a constant

  #elif Bingham
    tauy = 0.1;
    mu = 0.1;

  #elif Bagnold

    I0 = 0.3; 
    mu0 = 0.1;
    deltamu = 0.26;
    dg = 0.04;
    rho = 1.0;

  #elif BagnoldC
    I0 = 0.3;
    mu0 = 0.1;
    deltamu = 0.26;
    dg = 0.04; 
    rho = 1.0; 
    tauc = 0.01;
  #endif

  NITERMIN = 2;
  TOLERANCE = 1e-5;

  run();
}

event init (i = 0)
{
	
  foreach() {
   	zb[] = 0.; //!!
	eta[] = HR;  

	//By default, same relative thickness for each layer 
	foreach_layer()  h[] = HR/nl;
   }

  #if Geometric 
  geometric_beta (1.01,Top);
  #endif

  //Analytical solution is set as inital velocity
  foreach (serial) {
    foreach_layer(){
      #if Poiseuille
        u.x[] = Up(point,HR);
      #elif Bingham
        u.x[] = Ubgm(point,HR,tauy,mu);
      #elif Bagnold
	    u.x[] = Ub(point,HR,I0,mu0,deltamu,rho,dg);
      #elif BagnoldC
        u.x[] = Ubc(point,HR,I0,tauc,deltamu,rho,dg);
      #endif
	}
  }
    
}


event acc(i++)
{
  foreach (serial) {
    foreach_layer(){
      u.x[] = u.x[] + G*sin(slope)*dt;
	}
  }
}

// Check for the convergence
double uold,unew;

event convergence (t = 1; t += 5; t <= T0){

  unew = 0.;
  foreach() 
    foreach_layer()
      unew = unew + u.x[];
  printf("%lf\n",fabs(unew-uold)/(nl*N));
  if (i>1 && fabs(unew-uold)/(nl*N)<TOLERANCE ) return 1;
  else uold = unew;

}


event profiles (t += 5;t<=T0)
{
  double H ;
  double shear;
  double z ;
  double Hf;
  file1 = fopen ("profile.dat", "w");


  if ((int)t%1==0){
       printf("time =  %g/%g \n",t,T0);
     }

  foreach (serial){
    Hf = 0.0;
	z = 0.0;
	for (int l = 0; l<nl;l++){
      Hf = Hf + h[0,0,l];
      H = 0.;
	  for (int ll=0;ll<l;ll++){
	    H += h[0,0,ll];
	  }
      z = H+ 0.5*h[0,0,l];
      //Compute the shear---
      shear = 0.0;
      if(l>=0&&l<(nl-1)){
        shear = (u.x[0,0,l+1]-u.x[0,0,l])/(0.5*(h[0,0,l]+h[0,0,l+1]));}
      else if(l==0){
        shear = (u.x[0,0,l]-0.0)/(0.5*h[0,0,l]);}
      else if(l==nl-1){
        shear = (u.x[0,0,l]-u.x[0,0,l-1])/(0.5*(h[0,0,l]+h[0,0,l-1]));}
      //- - - - - - - - - - - - - - -

      #if Poiseuille
        fprintf (file1,"%g %g %g %g %g %g\n", x,z,u.x[0,0,l],Up2(z,HR),shear,Sp(z,HR));
      #elif Bingham
        fprintf (file1,"%g %g %g %g %g %g\n", x,z,u.x[0,0,l],Ubgm2(z,HR,tauy,mu),shear,Sbgm(z,HR,tauy,mu));
      #elif Bagnold
        fprintf (file1,"%g %g %g %g %g %g\n", x,z,u.x[0,0,l],Ub2(z,HR,I0,mu0,deltamu,rho,dg), shear, Sb(z,HR,I0,mu0,deltamu,rho,dg));
      #elif BagnoldC
        fprintf (file1,"%g %g %g %g %g %g\n", x,z,u.x[0,0,l],Ubc2(z,HR,I0,tauc,deltamu,rho,dg), shear,Sbc(z,HR,I0,tauc,deltamu,rho,dg));
      #endif
  }

    fprintf (file1,"\n \n");   
  } 

  fclose (file1);
}

/**

~~~gnuplot Velocity and stress profiles for flow
 set xlabel "y"
 set ylabel "u, shear"
p 'profile.dat' u 2:3 w lp t'U comp', "" u 2:4 w l t'U exact', '' u 2:5 w lp t'shear comp', '' u 2:6 w l t'shear exact'

~~~
*/ 