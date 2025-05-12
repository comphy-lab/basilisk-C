/** 
# Stability of the meniscus

The stability of a meniscus is not an easy test in numerical simulation. It is 
therefore a strong test case to check Laplace pressure calculation. In this 
example, we check the ability of the multilayer solver to capture the shape of a 
stable meniscus imposed by an homogeneous surface tension, zero gravity and contact 
angle boundary conditions.
![Meniscus](meniscus/movie.mp4)
*/ 

/**
## Includes and parameters
This test case uses the multilayer solver described in 
[Popinet, 2020](/Bibliography#popinet2020) and the surface tension additive.*/
#include "grid/multigrid1D.h"
#include "../hydro.h"
#include "layered/remap.h"
#include "../nh_PhiS.h"
#include "../tension.h"

/**
## Geometry */
#define L 4.                 // size of the box
#define T_END 20.            //end of simulation
#define H_ERR 1e-12
#define LEVEL 6              //default horizontal resolution
#define layers 4             //default vertical resolution
#define liquid_theta (60.*pi/180.)  // contact angle of the liquid on the walls

/**
We define the shape of the liquid between the two walls. */
#define Rf (L/(2.*cos(liquid_theta)))
#define Yc (1. + L/pow(2.*cos(liquid_theta), 2)*(pi/2. - liquid_theta + cos(liquid_theta)*sin(liquid_theta)))
#define shape(x) (Yc - pow((sq(Rf) - sq(x - L/2.)), 0.5))

/**
## Main function 
The test is done with viscosity but without gravity. Errors are very few 
dependant on the number of layers, except for very high horizontal which need
a sufficiant number of layers to be stable*/
FILE * fp = NULL;
int main() {
  corr_dux = true;
  L0 = L;
  origin (0.);
  G = 0;
  nl = layers;
  gradient = NULL ;
  nu = 1.;
  CFL = CFL_H = 0.5;
  
  for (N = 4; N <= 128; N *= 2)
    for (nl = 2; nl <= 4; nl += 2)
      run();
  
  nl = layers;
  N = 1 << LEVEL;
  blue = 0;
  run();
}

/**
## Initialisation */
scalar Hi[];
event init (i = 0) {
  /** 
  We set the boundary conditions on the height functions thanks to curvature.h.
  Surface tension is homogeneous. */
  eta[left] = contact_angle (liquid_theta, L0/N);
  eta[right] = contact_angle (liquid_theta, L0/N);
  foreach(){
    double H = shape(x);
    Hi[] = H;
    foreach_layer() {
      h[] = H/nl;
      sigma[] = 1.;
    }
  }
  boundary({eta, sigma, Hi}); 

  /** 
  We store the initial shape of the meniscus. */
  if (N == 1 << LEVEL && nl == layers) {
    char name[80];
    sprintf (name, "velocity-blue%d", blue);
    fp = fopen (name, "w");
    
    FILE * fps = fopen ("init_shape", "w");
    foreach() {
      fprintf(fps, "%g %g\n", x, Hi[]);
      fflush(fps);
    }
  }
}

/**
## Outputs */

event output (i+=20; t <= T_END)
{
  /**
  We output the evolution of the height convergence and the maximum velocity. */
  double dH = change (eta, Hi)*N/L0;
  if (N == 1 << LEVEL && nl == layers) {
    double un = 0;
    foreach()
      foreach_layer()
        un = max(un, norm(u));
    fprintf (fp, "%g %g %g\n", t, un, dH);
  }
}

event logfile (t = end) {
  /**
  At the end of the simulation, we recompute the reference solution. */
  scalar Href[];
  Href[left] = contact_angle (liquid_theta, L0/N);
  Href[right] = contact_angle (liquid_theta, L0/N);
  foreach()
    Href[] = shape(x);
  boundary({Href});
  
  /**
  And output the norm of the velocity *un* and the shape error *ec*. */
  scalar ec[];
  double un = 0;
  foreach() {
    foreach_layer() 
        un = max(un, norm(u));
    ec[] = (eta[] - Href[])/Delta;
  }

  norm ne = normf (ec);
  fprintf (stderr, "%d %d %d %g %g %g %g %g\n", 
	         N, nl, blue, liquid_theta, un, 
	         ne.avg, ne.rms, ne.max);
  fflush(stderr);

  /** We save the final shape of the interface.*/ 
  if (N == 1 << LEVEL && nl == layers) {
    char name[80];
    sprintf (name, "final_shape-blue%d", blue);
    FILE * fps = fopen (name, "w");
    foreach() {
      fprintf(fps, "%g %g\n", x, eta[]);
      fflush(fps);
    }
  }
}

/**# Movie
Plot the evolution in time with the max speed*/
#if 1
void setup (FILE * fp)
{
  fprintf (fp,
	   "set pm3d map interpolate 2,2\n"
	   "# jet colormap\n"
	   "set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25"
	   " 0 0.5647 1, 0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625"
	   " 1 0.9333 0, 0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )\n"
	   "unset key\n"
	   "set cbrange [0:5e-6]\n"
	   "set xlabel 'x'\n"
	   "set ylabel 'height'\n"
	   "set xrange [0:4]\n"
	   "set yrange [0:1.5]\n"
     "set rmargin at screen 0.8\n"
	   );
}

void plot (FILE * fp)
{
  fprintf (fp,
	   "set title 't = %.2f'\n"
	   "sp '-' u 1:2:3\n", t);
  foreach (serial) {
    double z = zb[];
    fprintf (fp, "%g %g %g\n", x, z, fabs(u.x[]));
    foreach_layer() {
      z += h[];
      fprintf (fp, "%g %g %g\n", x, z, fabs(u.x[]));
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 1\n");
  fflush (fp);  
}

event gnuplot (t += 0.2; t <= T_END/2.)
{
  if (N == 1 << LEVEL && nl == layers && blue == 1) {
    static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
    fprintf (fp,
  	   "set term pngcairo font \",10\" size 900,300\n"
  	   "set output 'plot-%04d.png'\n", (int)(t/0.2));
    if (i == 0)
      setup (fp);
    plot (fp);
  }
}

event moviemaker (t = end)
{

  if (N == 1 << LEVEL && nl == layers && blue == 1)
    system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
	    "ppm2mp4 movie.mp4");
}
#endif

/**
# Results

~~~gnuplot Final shapes using multilayers

darkgray="#666666"
purple="#9E5DB6"
blue="#5082DC"
turquoise="#008C7D"
forest="#149632"
orange="#FF780F"
raspberry="#FA0F50"
set style line 1 pt 7 ps 0.7

set terminal svg enhanced size 640,640 font ",8"
set output '_final_shape.svg'
set key 
set border
set tics
set xlabel "x"
set ylabel "y"
set size ratio 0.2

h0 = 1.
L0 = 4.
theta = pi/3.
Yc = h0 + L0/(4.*cos(theta))*((pi/2. - theta)/cos(theta) + sin(theta))
r = L0/(2.*cos(theta))

plot \
  '../contact_tension/final_shape' u 1:2 w l lc rgb orange t 'VOF' , \
  'final_shape-blue1' u 1:2 w l lc rgb blue t 'multilayer with blue term', \
  'final_shape-blue0' u 1:2 w l lc rgb raspberry t 'multilayer without blue term', \
  'final_shape-blue1' u 1:(Yc - (r**2 - ($1 - L0/2.)**2)**(0.5)) \
    w l lc rgb forest t 'model', \
  'init_shape' u 1:2 w l lc rgb darkgray t 'init'
~~~
  
The interface looks good for this implementation of the surface tension (and 
contact angle).

~~~gnuplot Error on the final shapes with the height functions
set output '_delta_final_shape.svg'
set xlabel "x"
set ylabel "y"
set size ratio 1

plot \
  '../contact_tension/init_shape' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb darkgray t 'VOF : init - model', \
  '../contact_tension/final_shape' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb orange t 'VOF - model' , \
  'final_shape-blue1' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb blue t 'multilayer with blue term - model', \
  'final_shape-blue0' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb raspberry t 'multilayer without blue term - model', \
  'init_shape' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb darkgray t 'multilayer : init - model'
~~~

~~~gnuplot Error on the final shapes with the height functions
set output '_delta_final_shape2.svg'

plot 'final_shape-blue1' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb blue t 'multilayer with blue term - model', \
  'init_shape' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb darkgray t 'init - model'
~~~

Looking more precisely at the error (the difference between the obtained interface 
and the theoretical position of the interface), we observe firstly that multilayer 
describes more accurately the circle shape (initial shape is more precise) than 
the VOF description. Then, the implementation of surface tension is very stable, 
event better than surface tension in VOF : the final error is smaller. The 
so-called "blue term" appears to be essential in this formulation of surface 
tension. 
  
~~~gnuplot max velo
set terminal svg font ",8"
set output '_velocity.svg'
set xlabel "t"
set ylabel "u"
set logscale y
set format y "%.2e"
set size ratio 1
plot [0:20]'../contact_tension/velocity' u 1:2 w p ls 1 lc rgb orange t 'VOF', \
     'velocity-blue1' u 1:2 w p ls 1 lc rgb blue t 'with blue term', \
     'velocity-blue0' u 1:2 w p ls 1 lc rgb raspberry t 'without blue term'
~~~

~~~gnuplot multi changes
set terminal svg enhanced size 640,640 font ",8"
set output '_multi_changes.svg'
set xlabel "t"
set ylabel "df"
set size ratio 1
plot [0:20]'../contact_tension/velocity' u 1:3 w p ls 1 lc rgb orange t 'VOF', \
     'velocity-blue1' u 1:3 w p ls 1 lc rgb blue t 'with blue term', \
     'velocity-blue0' u 1:3 w p ls 1 lc rgb raspberry t 'without blue term'
~~~

~~~gnuplot Final error with horizontal resolution
set output '_errorN.svg'
set xlabel "N"
set ylabel "error"
plot 'log' every 2::0::10 u 1:6 w l ls 1 t 'nl = 2', \
  'log' every 2::1::11 u 1:6 w l t 'nl = 4'
~~~

Precision do not depend of the number of layer (which is coherent as this problem 
is only geometric). However, the error reduces slightly with the horizontal 
resolution.*/