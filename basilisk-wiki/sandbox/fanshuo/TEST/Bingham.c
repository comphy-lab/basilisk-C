/**
##   Free surface flow of a Bingham fluid with layered solver
An example of 2D complex flow  over a plate with a free surface is presented here.
The configuration is periodic. 
*/

#define ML 1
#define HYDRO 1
#define RHEOLOGY 1
#define BINGHAM 1

#include "grid/multigrid1D.h"
#if !ML
# include "saint-venant.h"
#else // ML
# include "hydroF.h"

# define phi q
# if !HYDRO
#   include "layered/nh.h"
# endif
# include "layered/remap.h"
//# include "layered/perfs.h"
#endif // ML


const double HR = 1., NU = 0.1, T0 = 30;
//double slope = 0.25[0];
scalar uold[];

/*
Functions reserved for Bingham flow
*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - */
double Ubgm(Point point)
{
  double Yc = 0.0;
  double zc = 0.0;
  double Umax;
  for (int l = - point.l; l < nl - point.l; l++) {
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;

  Yc = HR - tauy/(G*sin(slope));
  if (Yc>0) Umax = pow(tauy-HR*G*sin(slope),2)/(2*mu*G*sin(slope));
  else {
    printf("No flow, change the rheology parameters.\n");
	return 1;
  }

  if(zc>=Yc) return Umax;
  else return Umax*(1-pow(zc/Yc-1,2));

}

double Ubgm2(double z,double HR)
{
  double Yc = 0.0;
  double Umax;

  Yc = HR - tauy/(G*sin(slope));
  if (Yc>0) Umax = pow(tauy-HR*G*sin(slope),2)/(2*mu*G*sin(slope));
  else {
    printf("No flow, change the rheology parameters.\n");
	return 1;
  }

  if(z>=Yc) return Umax;
  else return Umax*(1-pow(z/Yc-1,2));

}

double Sbgm(double z,double HR)
{
  double Yc = 0.0;
  double Umax;

  Yc = HR - tauy/(G*sin(slope));
  if (Yc>0) Umax = pow(tauy-HR*G*sin(slope),2)/(2*mu*G*sin(slope));
  else {
    printf("No flow, change the rheology parameters.\n");
	return 1;
  }

  if(z>=Yc) return 0;
  else return -2*Umax/Yc*(z/Yc-1);
}


/*- - - - - - - - - - - - - - - - - - - - - - - - - */
int main()
{
  periodic(right);
  L0 = 1.;
  G = 1.;
  N = 16; 
  nu = NU;
  nl = 64; 
	lambda_b[] = {0.,0.,0.};
	u_b[] = {0.,0.,0.};
	lambda_t[] = {-HUGE,0.,0.};
	u_t[] = {0.0,0.,0.};
#if !HYDRO  
  NITERMIN = 2;
#endif

//Bingham parameters
  tauy = 0.1;
  mu = 0.1;
  slope = 0.25;
 for (N = 8; N<= 32; N*= 2){
    for (nl = 4; nl <= 256; nl *= 2)
      run();
    fprintf (stderr,"\n\n");
  } 
 
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

//We initialize the velocity 
  foreach (){ 
    foreach_layer()
     // u.x[] = Ubgm(point);
    uold[] = 0; 
  }
/**
In the non-hydrostatic case, a boundary condition is required for
the non-hydrostatic pressure $\phi$ of each layer. */
  
#if !HYDRO && ML
  phi[right] = dirichlet(0.);
#endif

}


event acceleration (i++, t<=T0)
{
  foreach_face()
    foreach_layer()
      ha.x[] += G*sin(slope)*hf.x[]; 
}


/**
## Outputs


event error (t = end ) {
  scalar e[];
  foreach(){
    e[] = u.x[] - Ubgm2(z,HR);
  }
  norm n = normf (e);
fprintf (stderr, "%d %d %g %g %g\n", N, nl, n.avg, n.rms, n.max);
}


We save profiles at regular intervals. */

#if 0
event profiles (t += 5;t<=T0)
{
  foreach() {
#if !ML
    double H = h[];
#else
    double H = 0.;
    foreach_layer()
      H += h[];
#endif
    fprintf (stderr, "%g %g %g %g\n", x, (zb[] + H), H, t);
  }
  fprintf (stderr, "\n\n");
}
#endif

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
event error (t = end)
{
  int i = 0;
  foreach() {
    if (i++ == N/2) {
      double z = zb[];
      double norm = 0, norm2 = 0, normax = 0; 
 
      foreach_layer() {
        double e = fabs(u.x[] - Ubgm2(z+0.5*h[],HR) );
        norm += e*h[];
        norm2 += sq(e)*h[];
        normax = max( normax, e );
        z += h[];
      }
      norm = norm/z;
      norm2 = sqrt(norm2/z);
      fprintf (stderr, "%d %d %g %g %g %g %g\n", N, nl, norm, norm2, normax, dt, t);
    }
  }
}

event output (t=end )
{
  double shearnum = 0.;
  
  if ( N == 16 && nl == 128 ){
    foreach () {
      double z = zb[];
      foreach_layer() { 
      shearnum = shear(point, u.x, h, 0, point.l);

      fprintf (stdout,"%g %g %g %g %g %g %g\n", x, z, u.x[], w[],Ubgm2(z,HR),Sbgm(z,HR),shearnum);
         z += h[];
    }
  fprintf (stdout,"\n \n");
  }
  
	fclose(stdout);
  }
	
}

/**
## Result: Velocity and shear profile

~~~gnuplot Velocity and shear profiles for Bingham flow (N = 32, nl = 32)
 set xlabel "y"
 set ylabel "u, shear"
 p "out" u 2:3 w p t'U computed',"" u 2:5 w l t'Uexact', "" u 2:7 w p t"shear computed", "" u 2:6 w l t'shear exact'
~~~

# Convergence
~~~gnuplot Variation of the relative error with number of layers, for different grid $N$.
set key bottom left
set xlabel 'nl'
set ylabel 'Max |e|'
set logscale
#set format x "%.0e"
set format y "%.1e"

set xrange [2:512]
set cbrange [1:2]
set xtics 2,2,512

set yrange [1e-4:1]
#set cbrange [1:2]
#set xtics 1e-5,10,1e-1
set label 1 "N^{-2}" at 8,0.4*16**-2 font ',12' textcolor rgb 'purple' offset 0,1
set label 2 "N^{-1}" at 8,1.4*16**-1 font ',12' textcolor rgb 'red' offset 0,1


set grid ytics
plot for [i=0:3] 'log' index i u 2:5 t "N=".columnhead(1) with lp lw 2 ps 1.5 lt i+2,\
         [4:8<<2] 0.5*x**-2 t '' w l lw 2 lc rgb 'purple',\
[4:8<<2] 0.8*x**-1 t '' w l lw 2 lc rgb 'red'
~~~
*/






















