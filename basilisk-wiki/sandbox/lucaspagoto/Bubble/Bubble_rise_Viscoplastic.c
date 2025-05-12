/**
# Bubble rise in a Bingham surrounding 

This code is for the case of a bubble rising n a Bingham (inelastic viscoplastic) surrounding

*/
  
#include "axi.h"         // Axisymmetric simulation
#include "navier-stokes/centered.h"
#include "vof.h"         // The header two-phase could be used instead
#include "tension.h"
#include "view.h"

/**Here, we define the dimensionless numbers governing the flow:*/

#define Fr 25.0          // Froude number
#define Bo 5.0           // Bond Number
#define Pl 0.02          // Plastic NUmber
#define MUR 100.         // Viscosity Ratio
#define RHOR 100.        // Density Ratio

/** Some flow geometries*/
/** $Ar = \sqrt{\rho_1|\Delta \rho|g D^3} / \mu_1$; $Bo = |\Delta \rho|g D^2 / \sigma$; $\mu_r = \mu_1 / \mu_2$; $\rho_r = \rho_1 / \rho_2$*/
/** The indexes 1 and 2 are for the surrounding and bubble phases, respectively.*/

#define R0 0.5           // Drop Radius
#define HEIGHT 50*R0     // Domain Size
#define Xo 20.*R0        // Drop Initial Position

/** The variables of the flow */
scalar rhov[];
face vector alphav[], muv[], muv2[];
scalar f[], * interfaces = {f};

double xdd = 0.;         // Drop center of mass position
double vdd = 0.;         // Drop center o mass velocity
double uddx = 0.;        // Drop volume

double rho1, rho2, mu1, mu2, U, tau_y, mu_p;

/** By balacing the bouyant ($\tau_b = |\Delta\rho| g D$) and viscous ($\tau_i = \mu_1 U/D = \tau_y + \mu_p U/D$) stress, the characteristic velocity $U = \cfrac{|\Delta \rho|gD^2 - \tau_y D}{\mu_p}$ is obtained.
*/

/* Simulations time */
double tEnd = 12.01;

/* Mesh refinement levels */
int LEVEL = 10;
int MINLEVEL = 6;


int main ()
{
  size (HEIGHT);
  init_grid (1 << MINLEVEL);

  rho1 = 1.;
  rho2 = rho1/RHOR;
  U = sqrt(Fr*(rho1 - rho2));
  tau_y = Pl*(rho1 - rho2);
  mu_p = U/Fr - tau_y/U;
  mu1 = U/Fr;
  mu2 = mu1/MUR;  
  f.sigma = (rho1 - rho2)/Bo;

  alpha = alphav;
  rho = rhov;
  mu = muv;

  TOLERANCE = 1e-4;  // Tolerance
  NITERMAX = 100;    // Maximum number of iterations
  run();
}


/** Gravity acceleration */
event acceleration (i++)
{
  face vector av = a;
  foreach_face(x)
    av.x[] = -1.;
}


/** Initialization */
event init (t = 0)
{
  if (!restore (file = "restart"))
  {
    refine (sq(x - Xo) + sq(y) - sq(1.2*R0) < 0 && level < LEVEL);
    fraction (f, -sq(x - Xo) - sq(y) + sq(R0));
  }
}


/** Calculation of the vicosity and density in each phase */
scalar yielded[];    // Yielded region
event properties (i++) 
{
  scalar eta[];        // Viscosity function
  scalar fa[];         // For a smoothing of the viscosity jump
  scalar sr_norm[];    // Strain rate norm

/** Here, we calculate the strain rate magnitude and the viscosity field*/

  double Nreg = 1e4;                    // Dimensionless regularization parameters
  double epsilon = tau_y/(Nreg*mu1);    // Regularization parameter

  double shear_xx, shear_xy, shear_yy, shear_oo;
  shear_xx = 0.;
  shear_xy = 0.;
  shear_yy = 0.; 
  shear_oo = 0.;

  foreach(serial)
  {
    shear_xx = 2*(u.x[1,0] - u.x[-1,0])/(2*Delta);
    shear_xy = 2*((u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0]))/(4.*Delta);
    shear_yy = 2*(u.y[0,1] - u.y[0,-1])/(2*Delta);
    shear_oo = 2*u.y[]/y;

    sr_norm[] = sqrt(0.5)*sqrt( sq(shear_xx) + 2*sq(shear_xy) + sq(shear_yy) + sq(shear_oo) );

    eta[] = (tau_y/(sr_norm[] + epsilon) + mu_p);


    fa[] =  (4.*f[] + 2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) + f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;    
    
    rhov[] = cm[]*((1 - f[])*rho1 + f[]*rho2);

    yielded[] = (f[] < 0.5) ? (eta[]*sr_norm[] < tau_y ? 0 : 1) : 1;
  }
//  boundary ({fa, eta, rhov});


// Here we pass the viscosity, specific volume and density.
  foreach_face()
  {
    double fm1 = (fa[] + fa[-1])/2.;
    double etam = (eta[] + eta[-1])/2.;   // Viscosity is passed at the faces of the cells
    
    muv.x[] = fm.x[] / ((1. - fm1)/etam + fm1/mu2);
    alphav.x[] = fm.x[]/ ((1 - fm1)*rho1 + fm1*rho2);
  }
//  boundary ({muv,alphav,rhov});
}


/** Mesh adaptation  */
event adapt (i++)
{
    adapt_wavelet ({f, u, yielded}, (double[]){1e-3, 1e-3, 1e-3, 1e-2}, LEVEL, MINLEVEL);
}


/** Calculating drop position and velocity */

event thin_film (i++)
{  
  scalar drox[], droy[], intx[], inty[];
  position (f, drox, {1,0});
  position (f, droy, {0,1});
  double positmax, radiusmax, positmin;
  positmax = statsf(drox).max;   // Position of the bubble front
  radiusmax = statsf(droy).max;  // Maximum position of the bubble side
  positmin = statsf(drox).min;   // Position of the bubble back
  
  
  double xd1 = 0.;
  double vd1 = 0.;
  double udx1 = 0.;
  foreach(reduction(+:vd1) reduction(+:xd1) reduction(+:udx1))
  {
     vd1 += dv()*f[];
     xd1 += dv()*f[]*x;
     udx1 += u.x[]*dv()*f[];
  } 

  xdd = xd1;
  vdd = vd1;
  uddx = udx1;

  char nameu[50];
  sprintf (nameu, "pos.txt");
  static FILE * fname = fopen (nameu, "w");
  fprintf (fname, "%.8g %.8g %.8g %.8g %.8g %.8g %.8g\n", t, vd1, xd1/vd1, udx1/vd1, positmax, radiusmax, positmin);
  fflush (fname);
}


// Print
event logfile(i+=1)
{
  scalar le[];
  int j = 0;
  foreach(reduction(+:j))
  {
    le[] = level;
    j++;
  }
  stats s = statsf(le);

  printf ("i = %d t = %g, # = %d, l.mi = %g, l.ma = %g\n", i, t, j, s.min, s.max);
  fflush(stdout);
}



/** Generating the bubble shape and variable field */
/*
event outputs (t += 1.0)
{

  char namei1[80];
  sprintf (namei1, "interface-bubble-%.2f.txt", t);
  FILE* fp1 = fopen (namei1, "w");
  output_facets (f, fp1);
  fclose (fp1);

  char field2[80];
  sprintf (field2, "fields-output-%.2f.txt", t);
  FILE* fld2 = fopen (field2, "w");
  output_field ((scalar *){u, ,f}, fld2, n = 500, linear = true,  box = {{xdd/vdd - 5*R0, 0.},{xdd/vdd + 5*R0, 5*R0}});
  fclose (fld2);
}
*/


/** Backup */
/*
event snapshot (t = 0; t <= tEnd; t += 1.0)
{
  char named[80];
  sprintf (named, "dump-%.2f", t);
  dump (file = named);
}
*/


event movie_mesh (t+=0.04)
{
  clear(); 
  view(fov = 5, tx = 0., ty = -0.4, psi = -pi/2);
  
  translate (x =  Xo - xdd/vdd)
  {
  draw_vof ("f", lw = 4, lc = {1,0,0});
  cells();
  mirror (n = {0,-1})
        {
          draw_vof ("f", lw = 4, lc = {1,0,0});
          cells();
        }
  }
  save("movie_mesh.mp4");
}

event movie_fields (t+=0.04)
{
  clear(); 
  view(fov = 10.0, tx = 0., ty = -0.6, psi = -pi/2);

  draw_vof ("f");
  squares("f", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f");
          squares("f", linear = true);
        }
  save("movie_vf.mp4");

  draw_vof ("f");
  squares("yielded", linear = true);
  mirror (n = {0,-1})
        {
          draw_vof ("f");
          squares("yielded", linear = true);
        }
  save("movie_yielded.mp4");
}



event end (t = tEnd) {}


/**
## Results
### Volume fraction field
![Volume fraction field](Bubble_rise_Viscoplastic/movie_vf.mp4)

### Yielded region
![Yielded region](Bubble_rise_Viscoplastic/movie_yielded.mp4)

### Mesh
![Mesh](Bubble_rise_Viscoplastic/movie_mesh.mp4)
*/

/**
~~~gnuplot Velocity profile: comparison between numerical and analytical solutions
reset
set title "Center of mass velocity"
set xlabel 'Time'
set ylabel 'Velocity'
set xrange [0:12]
set yrange [0:0.2]
plot 'pos.txt' u 1:4 w l lw 3 lc 7 notitle
~~~
*/