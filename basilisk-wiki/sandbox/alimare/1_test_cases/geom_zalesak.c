/** We want to reproduce the geometry from the Zalesak's test case. We use the
[Constructive Geometry Paradigm](https://en.wikipedia.org/wiki/Constructive_solid_geometry). For a more general
example see the sandbox of [A. Berny](http://basilisk.fr/sandbox/aberny/csgBool.c)

#Constructing the Zalesak's geometry

The Zalesak's disk (Z) is a notched disk, which can be seen as the difference
between a circle (C) and a rectangle (R):
$$ Z = C - R $$
*/

#include "utils.h"
#include "fractions.h"
#include "distance.h"
#include "../alex_functions.h"
#define Pi 3.141592653589793

double geometry(double x, double y) {

  coord center_circle, center_rectangle, size_rectangle;
  center_circle.x = center_rectangle.x = 0.5;
  center_circle.y = 0.75;

  center_rectangle.y = 0.725;

  size_rectangle.x = 0.05;
  size_rectangle.y = 0.25; 

  double s = -circle (x, y, center_circle, 0.15);
  double r = -rectangle (x, y, center_rectangle, size_rectangle);

  double zalesak = difference(s,r) ;

  return zalesak;
}

int main() {

  /**
  We initialise the grid with 7 levels. */
  init_grid(1<<7);

  scalar f[];

  int iteration = 0;
  do {
    iteration++;
    fprintf(stderr,"#iteration %d \n", iteration);
    fraction(f, geometry (x, y));
  }while (adapt_wavelet({f}, (double []){0.2},
    maxlevel = 9, 2).nf != 0 && iteration <= 10);

/**
The obtained geometry can then be displayed:

![Reconstructed VOF surface.](geom_zalesak/vof.png)  

Here too:

![Levels of refinement.](geom_zalesak/levels.png)  

We recall that the facets obtained by the fraction() function do not define a
closed surface. As such they are not suitable for the initialization of the
level set function.
 */
  
  scalar l[];
  foreach()
    l[] = f[];
  output_ppm (l, file = "vof.png", n = 800, min = 0, max = 1);
  output_facets(f);

  foreach()
    l[] = level;
  output_ppm (l, file = "levels.png", n = 800, min = 0, max = 9);
  output_facets(f);

  FILE * fp1 = fopen ("zalesak.gnu","w");


/**
#Definition of the analytic curve of the 0-level set function

We trace a polyline with 2 points on the circular part and then add the last two
points. Of course we start and end with the same point. Thus, we have properly
defined a closed contour which can be used for the distance() function.

*/
  double theta_i = acos(sqrt(1-sq(0.025/0.15)));
  double theta_f = 2.*Pi-theta_i;
  double delta_theta = 0.5*2.*Pi/360.;
  double theta;
  fprintf(fp1,"0.525 0.625\n" ); // first point

  for(theta=theta_i; theta<theta_f; theta+=delta_theta){
    double yy = 0.75-0.15*cos(theta);
    double xx = 0.5+0.15*sin(theta);
    fprintf(fp1,"%g %g\n", xx,yy ); 
    }

  // and last but not least
  fprintf(fp1,"0.475 0.625\n" );
  fprintf(fp1,"0.475 0.85\n");
  fprintf(fp1,"0.525 0.85\n");

  fprintf(fp1,"0.525 0.625" ); // first point

  fclose(fp1);

  coord * p = input_xy (fopen ("zalesak.gnu", "r"));
  scalar d[];
  distance (d, p);
  boundary({d});

  foreach()
    l[] = d[];
  output_ppm (l, file = "dist_init.png", n = 800);
}


/**
~~~gnuplot Shapes of the interface
reset
set size ratio -1
plot [0:1][0:1]'zalesak.gnu' w l t "Zalesak's disk interfaces"
~~~
*/
