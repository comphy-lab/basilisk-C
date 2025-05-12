/**
# Stability of a meniscus

The stability of a meniscus is a difficult case for Laplace pressure
balance.  In this example we check the ability of the multilayer
solver to capture the shape of a stable meniscus due to an homogeneous
surface tension, zero gravity and contact angle boundary conditions.

![Meniscus](meniscus/movie.mp4)

This test case uses the multilayer solver described in [Popinet,
2020](/Bibliography#popinet2020), with its implicit implementation and
the surface tension additive. */

#include "grid/multigrid1D.h"
#include "crobert/2_Implicit/hydro-tension.h"
#include "layered/implicit.h"
#include "layered/remap.h"
// #include "layered/perfs.h"

/**
## Geometry */

#define L 4.                 // size of the box
#define T_END 20.            //end of simulation
#define H_ERR 1e-12
#define LEVEL 6              //default horizontal resolution
#define layers 4             //default vertical resolution
#define contact_angle (60.*pi/180.)  // contact angle of the liquid on the walls

/**
We define the shape of the liquid between the two walls. */

#define radius (L/(2.*cos(contact_angle)))
#define zc (1. + radius*sin(contact_angle))
#define meniscus_shape(x) (zc - sqrt(sq(radius) - sq(x - L/2.)))
  
/**
## Main function 

The test is done with viscosity but without gravity. The error does
not depend on the number of layers. */

int main()
{
  L0 = L;
  origin (0.);
  G = 0;
  nu = 1;
  
  N = 1 << LEVEL;
  nl = layers;
  TOLERANCE = 3e-12;

  /**
  We use a large timestep and the backward-Euler scheme to converge as
  quickly as possible toward the stationary solution. */
  
  CFL_H = 40;
  theta_H = 1;
  
  system ("rm -f plot-*.png");
#if 0
  run();
#else
  for (N = 4; N <= 128; N *= 2)
    for (nl = 2; nl <= 4; nl += 2)
      run();
#endif
}

/**
## Initialisation */

scalar Hi[];

event init (i = 0)
{

  /**
  This could be improved using a second-order discretisation. */
  
  eta[left]  = neumann (1./tan(contact_angle));
  eta[right] = neumann (1./tan(contact_angle));

  foreach() {
    double H = meniscus_shape(x);
    Hi[] = H;
    foreach_layer() {
      h[] = H/nl;
    }
  }
  boundary({Hi}); 

  /** 
  We store the initial shape of the meniscus. */
  if (N == 1 << LEVEL && nl == layers) {
    FILE * fps = fopen ("init_shape", "w");
    foreach()
      fprintf(fps, "%g %g\n", x, Hi[]);
    fclose (fps);
  }
}

/**
## Left and right walls 

This is necessary to balance the surface tension acceleration on the
side walls. */

event half_advection (i++)
{
  ha.n[left] = 0.;
  ha.n[right] = 0.;
  hu.n[left] = 0.;
  hu.n[right] = 0.;
  boundary ((scalar *){ha, hu});
}

/**
## Outputs */

event output (i++; t <= T_END)
{
  /**
  We output the evolution of the height convergence and the maximum velocity. */

  double dH = change (eta, Hi)*N/L0;
  if (N == 1 << LEVEL && nl == layers) {
    static FILE * fp = fopen ("velocity", "w");
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
  Href[left] = neumann (1./tan(contact_angle));
  Href[right] = neumann (1./tan(contact_angle));
  foreach()
    Href[] = meniscus_shape(x);
  boundary ({Href});
  
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
  fprintf (stderr, "%d %d %g %g %g %g %g %d\n", 
	         N, nl, contact_angle, un, 
	   ne.avg, ne.rms, ne.max, mgH.i);
  fflush (stderr);

  /**
  We save the final shape of the interface. */
  
  if (N == 1 << LEVEL && nl == layers) {
    FILE * fps = fopen ("final_shape", "w");
    foreach()
      fprintf (fps, "%g %g\n", x, eta[]);
    fclose (fps);
  }
}

/**# Movie
Plot the evolution in time with the max speed*/
#if 1
void setup (FILE * fp)
{
  fprintf (fp,
	   "set term pngcairo font \",10\" size 1024,400\n"
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
	   "set yrange [0:1.2]\n"
	   "set size ratio -1\n"
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

event gnuplot (i += 20)
{
  if (N == 1 << LEVEL && nl == layers) {
    static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
    if (i == 0)
      setup (fp);
    fprintf (fp, "set output 'plot-%04d.png'\n", (int)(t/0.2));
    plot (fp);
  }
}

event moviemaker (t = end)
{

  if (N == 1 << LEVEL && nl == layers)
    system ("for f in plot-*.png; do convert $f ppm:- && rm -f $f; done | "
	    "ppm2mp4 movie.mp4");
}
#endif

/**
# Results

~~~gnuplot Final shapes using multilayers
set term @SVG size 800,200
darkgray="#666666"
purple="#9E5DB6"
blue="#5082DC"
turquoise="#008C7D"
forest="#149632"
orange="#FF780F"
raspberry="#FA0F50"
set style line 1 pt 7 ps 0.7

set xlabel "x"
set ylabel "y"

h0 = 1.
L0 = 4.
theta = pi/3.
r = L0/(2.*cos(theta))
zc = h0 + (r*sin(theta))

set size ratio -1
set key below

plot \
  'final_shape' u 1:2 w l lc rgb blue t 'multilayer', \
  'final_shape' u 1:(zc - (r**2 - ($1 - L0/2.)**2)**(0.5)) \
    w l lc rgb forest t 'model', \
  'init_shape' u 1:2 w l lc rgb darkgray t 'init'
~~~
  
The interface looks good for this implementation of the surface tension (and 
contact angle boundary conditions).

~~~gnuplot Error on the final shapes - comparison with VOF
reset
set xlabel "x"
set ylabel "y"

plot \
  '../meniscus_VOF/init_shape' u 1:($2 - (zc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb darkgray t 'VOF : init - model', \
  '../meniscus_VOF/final_shape' u 1:($2 - (zc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb orange t 'VOF - model' , \
  'final_shape' u 1:($2 - (zc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb blue t 'multilayer - model', \
  'init_shape' u 1:($2 - (zc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb darkgray t 'multilayer : init - model'
~~~

~~~gnuplot Error on the final shapes
plot 'final_shape' u 1:($2 - (zc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb blue t 'Multilayer', \
  'init_shape' u 1:($2 - (zc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb darkgray t 'init - model'
~~~

Looking more precisely at the error (the difference between the final
interface and the theoretical position of the interface), we observe
that multilayer model describes more accurately the meniscus shape
(initial shape is more precise) than the VOF description. Then, the
implementation of surface tension is very stable, slightly better than
surface tension in VOF : the final error is smaller.
  
~~~gnuplot Maximum velocity as a function of time.
set xlabel "t"
set ylabel "u"
set logscale y
set format y "%.2e"
plot [0:20]'../meniscus_VOF/velocity' u 1:2 w p ls 1 lc rgb orange t 'VOF', \
     'velocity' u 1:2 w p ls 1 lc rgb blue t 'Multilayer'
~~~

~~~gnuplot Change in the free-surface height with time.
set xlabel "t"
set ylabel "df"
plot [0:20]'../meniscus_VOF/velocity' u 1:3 w p ls 1 lc rgb orange t 'VOF', \
     'velocity' u 1:3 w p ls 1 lc rgb blue t 'Multilayer'
~~~

~~~gnuplot Final error dependence on the horizontal resolution.
set xlabel "N"
set ylabel "error"
set logscale
set xtics 4,2,128

ftitle(a,b) = sprintf("order %4.2f", -b)
f1(x) = a1 + b1*x
fit [log(10):] f1(x) 'log' every 2::0::10 u (log($1)):(log($6)) via a1,b1
set xrange [3:192]
plot \
  exp (f1(log(x))) t ftitle(a1,b1), \
  'log' every 2::0::10 u 1:6 w p ls 1 t 'nl = 2', \
  'log' every 2::1::11 u 1:6 w p t 'nl = 4'
~~~

The accuracy does not depend of the number of layers (which is
coherent as this problem is only geometric), however convergence is
only first-order in the horizontal resolution. This may be due to the
first-order boundary condition on the contact angle (see comment
above). */
