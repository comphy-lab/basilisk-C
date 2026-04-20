/**
##   Free surface flow of a Bingham fluid with layered solver
An example of 2D complex flow  over a plate with a free surface is presented here.
The configuration is periodic. 
*/

#define ML 1
#define HYDRO 1
#define RHEOLOGY 1
#define BINGHAM 0

#include "grid/multigrid1D.h"
#if !ML
# include "saint-venant.h"
#else // ML
# include "./hydroMT.h"

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


FILE *file1;
/*- - - - - - - - - - - - - - - - - - - - - - - - - */
int main()
{
  periodic(right);
  L0 = 1.;
  G = 1.;
  N = 64; 
  nu = NU;
  nl = 256; 
	lambda_b[] = {0.1,0.,0.};
	u_b[] = {0.,0.,0.};
	lambda_t[] = {-0.2,0.,0.};
	u_t[] = {0.0,0.,0.};
#if !HYDRO  
  NITERMIN = 2;
#endif

//Bingham parameters
  slope = 0.3;
  
	file1 = fopen("profil.dat","w");
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

//We initialize the velocity 
  foreach (){ 
    foreach_layer()
      u.x[] = Up(point);
    //uold[] = 0; 
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

/** We stop when the difference between two time steps is small enough */
/*
event logfile (t += 0.1; t<=T0) {

  double du = change (u.x,uold);
    if (i > 0 && du < 1e-6)
    return 1;  stop 
}
*/
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

event output (t+=0.1;t<=T0 )
{
	double omega, m,n;	
	
	omega = G*sin(slope)/nu;
	file1 = fopen("profil.dat","w");
		printf("%g \n\n",t);
    foreach () {
		m = (-2.*u_b.x[]+2*u_t.x[]+HR*(HR-2*lambda_t.x[])*omega)/(2*(HR+lambda_b.x[]-lambda_t.x[]));
		n = (2*u_t.x[]*lambda_b.x[]+2.*u_b.x[]*(HR-lambda_t.x[])+HR*lambda_b.x[]*(HR-2*lambda_t.x[])*omega)/(2*(HR+lambda_b.x[]-lambda_t.x[]));
      double z = zb[];
      foreach_layer() { 
        z += h[];
      fprintf (file1,"%g %g %g %g\n", x, z, u.x[], -1./2.*omega*pow(z,2)+m*z+n);
    }
  fprintf (file1,"\n \n");
  }
  
	fclose(file1);
	
}

/**
##Velocity ans shear profile

~~~gnuplot Velocity and shear profiles for Bingham flow (N = 32, nl = 32)
 set xlabel "y"
 set ylabel "u, shear"
 p "out" u 2:3 w p t'U computed',"" u 2:5 w l t'Uexact', "" u 2:7 w p t"shear computed", "" u 2:6 w l t'shear exact'
~~~
*/

























