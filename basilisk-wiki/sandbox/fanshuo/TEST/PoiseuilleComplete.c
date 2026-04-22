

#define ML 1
#define HYDRO 1
#define RHEOLOGY 1
#define BINGHAM 0

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


const double HR = 1., NU = 0.1, T0 = 100;
//double slope = 0.25[0];
scalar uold[];

/*
Analytic	velocity for Poiseuille flow
*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - */
double Up(Point point)
{
  double omega = 0.0;
  double m = 0.0;
  double n = 0.;
	double zc;
  for (int l = - point.l; l < nl - point.l; l++) {
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;
  
	omega = G*sin(slope)/nu;
	m = (-2.*u_b.x[]+2*u_t.x[]+HR*(HR-2*lambda_t.x[])*omega)/(2*(HR+lambda_b.x[]-lambda_t.x[]));
	n = (2*u_t.x[]*lambda_b.x[]+2.*u_b.x[]*(HR-lambda_t.x[])+HR*lambda_b.x[]*(HR-2*lambda_t.x[])*omega)/(2*(HR+lambda_b.x[]-lambda_t.x[]));

return -1./2.*omega*pow(z,2)+m*z+n;

}


/*- - - - - - - - - - - - - - - - - - - - - - - - - */
int main()
{
  periodic(right);
  L0 = 1.;
  G = 1.;
  N = 16; 
  nu = NU;
  nl = 256; 
	lambda_b[] = {0.,0.,0.};
	u_b[] = {0.,0.,0.};
	lambda_t[] = {0.,0.,0.};
	u_t[] = {0.0,0.,0.};
#if !HYDRO  
  NITERMIN = 2;
#endif

//Bingham parameters
  slope = 0.3;
  
 for (N = 16; N<= 64; N*= 2){
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
    //  u.x[] = Up(point);
    uold[] = 0; 
  }
/**
In the non-hydrostatic case, a boundary condition is required for
the non-hydrostatic pressure $\phi$ of each layer. */
  
#if !HYDRO && ML
  phi[right] = dirichlet(0.);
#endif

}
#if !ML
event acc (i++, t<=T0){
  foreach () 
    for (vector u in ul) {
      u.x[] = u.x[] + G*sin(slope)*dt;
  }
}
#else
event acceleration (i++, t<=T0)
{
  foreach_face()
    foreach_layer()
      ha.x[] += G*sin(slope)*hf.x[]; 
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

event output (t=end )
{
	double omega, m,n;	
	
	omega = G*sin(slope)/nu;
  if ( N == 16 && nl == 128 ){
    foreach () {
		m = (-2.*u_b.x[]+2*u_t.x[]+HR*(HR-2*lambda_t.x[])*omega)/(2*(HR+lambda_b.x[]-lambda_t.x[]));
		n = (2*u_t.x[]*lambda_b.x[]+2.*u_b.x[]*(HR-lambda_t.x[])+HR*lambda_b.x[]*(HR-2*lambda_t.x[])*omega)/(2*(HR+lambda_b.x[]-lambda_t.x[]));
      double z = zb[];
      foreach_layer() { 
        z += h[];
      fprintf (stdout,"%g %g %g %g\n", x, z, u.x[], -1./2.*omega*pow(z,2)+m*z+n);
    }
  fprintf (stdout,"\n \n");
  }
  
	fclose(stdout);
  }
	
}

event error (t = end)
{
  double omega, m,n;
  int i = 0;
  foreach() {
    omega = G*sin(slope)/nu;
    m = (-2.*u_b.x[]+2*u_t.x[]+HR*(HR-2*lambda_t.x[])*omega)/(2*(HR+lambda_b.x[]-lambda_t.x[]));
		n = (2*u_t.x[]*lambda_b.x[]+2.*u_b.x[]*(HR-lambda_t.x[])+HR*lambda_b.x[]*(HR-2*lambda_t.x[])*omega)/(2*(HR+lambda_b.x[]-lambda_t.x[]));
    if (i++ == N/2) {
      double z = zb[];
      double norm = 0, norm2 = 0, normax = 0; 
    #if ML
      foreach_layer() {
        double e = fabs(u.x[] - (-1./2.*omega*pow((z+h[]/2.),2)+m*(z+h[]/2.)+n) );
        norm += e*h[];
        norm2 += sq(e)*h[];
        normax = max( normax, e );
        z += h[];
      }
    #else
        for (int l = 0; l <= nl - 1; l++) {
          vector u = ul[l]; 
          double e = fabs(u.x[] -(-1./2.*omega*pow(z,2)+m*z+n) );
          norm += e*h[]*layer[l];
          norm2 += sq(e)*h[]*layer[l];
          normax = max( normax,e );
          z += h[]*layer[l];
        }
    #endif
      norm = norm/z;
      norm2 = sqrt(norm2/z);
      fprintf (stderr, "%d %d %g %g %g %g %g\n", N, nl, norm, norm2, normax, dt, t);
    }
  }
}


/**
# Velocity  profile

~~~gnuplot 
 set xlabel "y"
 set ylabel "u"
 p "out" u 3:2 w p t'U computed',"" u 4:2 w l t'Uexact'
~~~

# Convergence
We track the value of the relative error on the velocity for various
number of layers $\text{nl}$. The error decrease with 2nd ordrer precision.

~~~gnuplot Variation of the relative error with number of layers, for different grid $N$.
set key bottom left
set xlabel 'nl'
#set ylabel 'Max |e|'
set logscale
#set format x "%.0e"
set format y "%.2e"

set xrange [2:512]
set cbrange [1:2]
set xtics 2,2,512

#set yrange [1e-6:1e-1]
#set cbrange [1:2]
#set xtics 1e-5,10,1e-1

set grid ytics
plot for [i=0:3] 'log' index i u 2:5 t "N=".columnhead(1) with lp lw 2 ps 1.5 lt i+2,\
         [4:8<<2] 0.02*x**-2 t '' w l lw 2 lc rgb 'purple'

~~~


*/
























