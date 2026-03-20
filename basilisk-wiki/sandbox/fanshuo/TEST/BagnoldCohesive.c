
#define ML 1
#define HYDRO 1
#define RHEOLOGY 1
#define BAGNOLDCOH 1


#include "grid/multigrid1D.h"
#if !ML
# include "saint-venant.h"
#else // ML
# include "hydroMT.h"
# define phi q
# if !HYDRO
#   include "layered/nh.h"
# endif
# include "layered/remap.h"
//# include "layered/perfs.h"
#endif // ML


const double NU = 0.1, T0 = 10000;
const double HR = 1. [1];
scalar uold[];
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - */
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
 foreach () {
	 foreach_layer(){
          u.x[] = Ubc(point);
           uold[] = 0.;
	}
    }

 }

event acc(i++){
  foreach () {
	 foreach_layer(){
          u.x[] = u.x[] + G*sin(slope)*dt;
     }
  }
}




/**
We check for convergence. */

event logfile (t=1;t += 0.1 ) {

  double du = change(u.x, uold);
  if(i>0 && du < 1e-6)
    return 1;
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


event output (t = end)
{
    foreach () {
      double z = zb[];
    foreach_layer() {
      z += h[];
      fprintf (stdout, "%g %g %g \n", x, z, u.x[]);
    }
      fprintf (stdout,"\n \n");   
  }
#if HYDRO
  delete ({w});
#endif
 
}

/**

~~~gnuplot Velocity and shear profiles for Bingham flow
 set xlabel "y"
 set ylabel "u, shear"
 p "profil.dat" u 2:3 w p t'U computed', "profil.dat" u 2:4 w l t'Uexact', "" u 2:5 w p t"shear computed", "" u 2:6 w l t'shear exact'
~~~
*/ 
