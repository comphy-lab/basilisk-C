#include "grid/multigrid1D.h"
#include "../../Emily/green-naghdi_two-layer.h"
//#include "../../Emily/two_layer.h"

#define MAXLEVEL 12 // 0.006 m resolution for 6.27 m tank 


double t_end=6; 
double mul[2];

/**
![Full tank view](gn/movie.mp4)(autoplay loop width="100%")

*/

//define Geometry

#define theta0 0.2638 //angle
#define width 0.338/0.265 //width of flow
#define gradient_m -1.*tan(theta0) //gradient of slope
#define Hi 1. //initial water depth
#define slope_exposed 1.
#define Hel (slope_exposed/0.265)*sin(theta0)
#define c_intercept tan(theta0)*((width+ (Hi+Hel)/tan(theta0)))
#define height_cm 0.18/0.265
#define frac 3./4.
#define c_intercept2 height_cm + ((frac*height_cm)/(1- frac))
//For testing, rho ratio is default

#define H_MFRIC 0.005
#define HS_FRIC 0.03
#define MANNINGS 0  // If set to 1 will run Manning's friction - uses H_MFRIC and HS_MFRIC otherwise quadratric and uses H_FRIC and HS_FRIC
#define H_FRIC 0.01
#define EPSILON 0.0

int main() {
  L0=32;     //size(TANK_SIZE);
  G = 1.;
  N = 1 << MAXLEVEL;
  system ("rm -r gnuplot");
  system ("mkdir gnuplot");
  //init_grid(1 << MAXLEVEL);
  run();
 }

double slide_depth,slide_vol;

event init (i=0) {
  
  foreach(){
    zb[] =  x < width ? Hel + Hi : (x  < (width + (Hi+Hel)/tan(theta0)) ?
      			   gradient_m*x + c_intercept  : 0); 
    eta[] = 1.;
    if ( x <= width*frac) {
      hs[]=height_cm;
    }
    else if ( x < width){
      hs[] = (-1*x*(height_cm/(width-(width*frac))) + c_intercept2 );
    }
    else {
      hs[] = 0.;
    }
    h[]=max(0.,eta[]-zb[]-hs[]);
    eta[] = zb[]+hs[]+h[];
  }

   
  boundary(all);
  RHOratio = 1.;
  mul[0]=0.001/1000; //Slide 
  mul[1]=0.001/1000; // water 
  DT=0.005;

  FILE * fp = fopen ("init", "w");
  foreach() {
    fprintf (fp, "%g %g %g %g %g %g\n", x, eta[], zb[], hs[], h[], tan(theta0) );
  }
  fprintf (fp, "\n");
  fclose (fp);
}

// Quadratic/Manning's friction implementation


event friction (i++) {
  foreach() {
    if (hs[] < dry) {
	#if MANNINGS // Manning's Friction implementation
      double a1_inv = h[] < dry ? 0. : h[]/(h[] + G*H_MFRIC*H_MFRIC*dt*norm(u)/pow(h[],0.3333333333333)+2*mul[1]/h[]*dt);
	#else // Quadratic Friction implementation
      double a1_inv = h[] < dry ? 0. : h[]/(h[] + G*H_FRIC*dt*norm(u)+2*mul[1]/h[]*dt);
	#endif
      foreach_dimension(){
	us.x[]=0.;
	u.x[] *= a1_inv;
      }
    }
    else {
      if (h[] < dry){
	#if MANNINGS // Manning's Friction implementation
	double a2_inv = hs[]/(hs[] + G*HS_MFRIC*HS_MFRIC*dt*norm(us)/pow(hs[],0.3333333333333)+2*mul[0]/hs[]*dt); // for landslide 
	#else // Quadratic Friction implementation
	double a2_inv = hs[]/(hs[] + G*HS_FRIC*dt*norm(us)+2*mul[0]/hs[]*dt); // for landslide
	#endif
	foreach_dimension(){
	  us.x[] *= a2_inv;
	  u.x[] *= 0.;
	}
      }
      else {
	#if MANNINGS // Manning's Friction implementation
	double a = 1 + G*HS_MFRIC*HS_MFRIC*dt*norm(us)/pow(hs[],1.3333333333333)+ 2*mul[0]*(2*h[]+hs[])/hs[]/hs[]/(hs[]+h[])*dt;
	#else // Quadratic Friction implementation
	double a = 1 + G*HS_FRIC*dt*norm(us)/hs[] + 2*mul[0]*(2*h[]+hs[])/hs[]/hs[]/(hs[]+h[])*dt;
	#endif
	double b =- 2*mul[0]/hs[]/(hs[]+h[])*dt;
	double c = -  2*mul[1]/hs[]/(hs[]+h[])*dt;
	#if MANNINGS // Manning's Friction implementation
	double d = 1 + G*H_MFRIC*H_MFRIC*dt*norm(u)/pow(h[],1.3333333333333)+2*mul[1]/h[]/(hs[]+h[])*dt;
	#else // Quadratic Friction implementation
	double d = 1 + G * H_FRIC*dt*norm(u)/h[] + 2*mul[1]/h[]/(hs[]+h[])*dt;
	#endif
	double invdet=1./(a*d-b*c);
	foreach_dimension(){
	  double us_old=us.x[];
	  double u_old=u.x[];
	  us.x[] = invdet*(d*us_old - b*u_old);
	  u.x[] = invdet*(-c*us_old + a*u_old);
	}
      }
    }
  }
  boundary ((scalar *){u,us});
}


void plot_density (FILE * fp, double t)
{
  fprintf (fp,
	   "set title 't = %.2f'\n" 
           "set term pngcairo font \",10\" size 1024,400\n"
           "set size ratio -1\n"
	   "p [0.:32][0.:3.]'-' u 1:3:($3+$4) w filledcu lc 4 t '',"
	   " '' u 1:($3+$4):2 w filledcu lc 3 t '',"
	   " '' u 1:(-1):3 w filledcu lc -1 t ''\n", t);
  foreach () {
    fprintf (fp, "%g %g %g %g %g\n", x, eta[], zb[], hs[], h[]);
  }

  fprintf (fp, "e\n\n");
  fflush (fp);
  fprintf (stderr, "%.3f %.3f\n", t, statsf(u.x).max);
}


event gnuplot (t += 0.01)
{
  //static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  static FILE * fp = popen ("gnuplot", "w");
  fprintf (fp,
	   "set output 'gnuplot/plot-%06d.png'\n", i);
  plot_density (fp, t);
 
}


/**We output profiles in separate files indexed by the time. */
event output (t += 0.1) {
  //event output (i+=100) {
  char name[80];
  sprintf (name, "out-%g", t);
  //sprintf (name, "grilli2016_out-%d", i);
  FILE * fp = fopen (name, "w");
  foreach() {
    fprintf (fp, "%g %g %g %g %g\n", x, eta[], zb[], hs[], h[]);
  }
  fprintf (fp, "\n");
  fclose (fp);
}

Gauge gauges[] = {
  // file x y description
  {"WG_3m.txt",3,0.,"% Wave Gauge 1"},
  {"WG_4m.txt",4,0.,"% Wave Gauge 2"},
  {"WG_5m.txt",5,0.,"% Wave Gauge 3"},
  {"WG_6m.txt",6,0.,"% Wave Gauge 4"},
  {NULL}
};

event gauges0 (i++) output_gauges (gauges ,{eta,h,hs});

event end (t = t_end)
{ 
  system ("for f in gnuplot/plot-*.png; do"
	  " convert $f ppm:- && rm -f $f; done | "
	  "ppm2mp4 movie.mp4");
  fprintf (stderr, "\n\nDone\n");
}
