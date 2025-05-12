/**
# Example for the integral formulation of the surface tension with contact angles with a multilayer solver.

We show here an example using an integral formulation of the surface tension
force, as depicted in
[Abu-Al-Saud, Popinet, Tchelepi, 2018](https://hal.archives-ouvertes.fr/hal-01706565/),
adapted to a multilayer solver.

For the moment, limitations remain. Nevertheless, this example already
shows that this implementation of the surface tension allow to describe
contact angles. 
*/

/**
## Includes */
#include "grid/multigrid1D.h"
#include "../hydro_c.h"
#include "../nh_c.h"
#include "layered/remap.h"
#include "../solute.h"
#include "../surface.h"

/**
## Paremeters 
The horizontal resolution, the number of layers and the contact angles must be wisely chosen so as to keep the dH in a cell below the cell size Delta. */

#define L 4. // size of the box
#define LEVEL 6
#define T_END 20. //4.
#define DELTA_T (T_END/100.)
#define H_ERR 1e-10
#define layers 2
#define _g 0.

/**
We set the value of the contact angle of the liquid on the walls. */
#define liquid_theta (60.*pi/180.)

/**
We define the geometry of the liquid between the two walls. */

#define Rf (L/(2.*cos(liquid_theta)))
#define Yc (1. + L/pow(2.*cos(liquid_theta), 2)*(pi/2. - liquid_theta + cos(liquid_theta)*sin(liquid_theta)))
#define shape(x) (Yc - pow((sq(Rf) - sq(x - L/2.)), 0.5))


/**
##Main function */

FILE * fp = NULL;

int main() {
  L0 = L;
  origin (0.);
  N = 1 << LEVEL;
  nl = layers;
  G = _g;
  run();
}

/**
## Boundary conditions

We set the boundary conditions on the height functions thanks to contact.h. */
H[left] = contact_angle (liquid_theta, L0/N, 1);
H[right] = contact_angle (liquid_theta, L0/N, 1);

/**
## Initialization
*/
scalar Hi[];
event init (i = 0) {
  scalar c, h;
  for(c, h in cl, hl)
    foreach(){
      double H = shape(x);
      c[] = 0;
      h[] = H/nl;
      gam[] = 1;
    }
  heights(hl, H);
  boundary({gam}); 
  foreach()
    Hi[] = H[];
  fp = fopen ("velocity", "w");
}

event advection_term (i++) {
  if (i == 0) {
    static FILE * fps = fopen ("init_shape", "w");
    heights(hl, H);
    foreach() {
      fprintf(fps, "%g %g\n", x, H[]);
      fflush(fps);
    }
  }
}

/**
## Movie and outputs */

event logfile (i++; t <= T_END)
{
  /**
  At every timestep, we check whether the volume fraction field has
  converged. */
  heights(hl,H);
  double dH = change (H, Hi)*N/L0;
  if (i > 1 && dH < H_ERR)
    return 1; /* stop */

  /**
  And we output the evolution of the maximum velocity. */

  scalar un[];
  foreach() {
    un[] = 0.;
    for (vector u in ul) 
      un[] = max(un[],norm(u));
  }  
  fprintf (fp, "%g %g %g\n", t, normf(un).max, dH);
}

event error (t = end) {
  
  /**
  We recompute the reference solution. */
  scalar Href[];
  Href[left] = contact_angle (liquid_theta, L0/N, 1);
  Href[right] = contact_angle (liquid_theta, L0/N, 1);
  foreach()
    Href[] = shape(x);
  boundary({Href});
  
  /**
  And compute the norm of the velocity *un* and the shape error *ec*. */
  scalar un[], ef[];
  foreach() {
    un[] = 0.;
    for (vector u in ul) 
      un[] = max(un[],norm(u));
    ef[] = (H[] - Href[])/Delta;
  }
  
  /**
  We output these on standard error (i.e. the *log* file). */
  norm ne = normf (ef);
  fprintf (stderr, "%d %g %g %g %g %g\n", 
	         LEVEL, liquid_theta, normf(un).max, 
	         ne.avg, ne.rms, ne.max);
}


/**
At the end, we save the final shape of the interface :*/ 
event final_shape (t = end) {
  static FILE * fps = fopen ("final_shape", "w");
  heights(hl, H);
  foreach() {
    fprintf(fps, "%g %g\n", x, H[]);
    fflush(fps);
  }
}

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

set terminal @PNG enhanced size 640,640 font ",8"
set output '_final_shape.png'
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
  '../contact_tension/final_shape' u 1:2 w l lc rgb blue t 'tension.h' , \
  'final_shape' u 1:2 w l lc rgb raspberry t 'multilayer', \
  'final_shape' u 1:(Yc - (r**2 - ($1 - L0/2.)**2)**(0.5)) \
    w l lc rgb forest t 'model', \
  'init_shape' u 1:2 w l lc rgb orange t 'init'
~~~
  
The interface deforms in time and tends to a triangle shape.

~~~gnuplot Error on the final shapes with the height functions
set output '_delta_final_shape.png'
set xlabel "x"
set ylabel "y"
set size ratio 1

plot \
  '../contact_tension/final_shape' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb blue t 'tension.h - model' , \
  'final_shape' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb raspberry t 'multilayer - model', \
  'init_shape' u 1:($2 - (Yc - (r**2 - ($1 - L0/2.)**2)**(1./2.))) \
    w l lc rgb orange t 'init - model'
~~~

Looking at the error (the difference between the obtained interface and the
theoretical position of the interface), we have, in the multilayer formulation, strong errors in the evolution of this surface, which is supposed to be stable at the original place.

~~~gnuplot max velo
set terminal @PNG font ",8"
set output '_velocity.png'
set xlabel "t"
set ylabel "u"
set logscale y
set format y "%.2e"
set size ratio 1
plot [0:20]'../contact_tension/velocity' u 1:2 w p ls 1 lc rgb blue t 'tension.h', \
     'velocity' u 1:2 w p ls 1 lc rgb raspberry t 'integral'
~~~

~~~gnuplot multi changes
set terminal @PNG enhanced size 640,640 font ",8"
set output '_multi_changes.png'
set xlabel "t"
set ylabel "df"
set size ratio 1
plot [0:20]'../contact_tension/velocity' u 1:3 w p ls 1 lc rgb blue t 'tension.h', \
     'velocity' u 1:3 w p ls 1 lc rgb raspberry t 'integral'
~~~

Those errors strongly reduce for increasing reolution or lower stifness of the interface.
How can we supress those errors ? Improve area interpolation ? Adding boundary conditions on the potential Phi ?
*/