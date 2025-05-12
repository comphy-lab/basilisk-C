/**
# The curvature of a circle with a normalized radius ($R$)  
*/

#include "grid/quadtree.h"
#include "fractions.h"
#include "curvature.h"
#include "utils.h"
#include "interface_iterator.h"

#define func (sq((x-0.01))+sq(y+(M_PI/100.))-1.) 

scalar f[], curv[];

int main(){
  /**
  First we try an equidistant grid with $\frac{\Delta}{R} = 3/256$
  */
  init_grid (256);
  L0 = 3.;
  X0 = Y0 = -L0/2.;
  f.prolongation = f.refine = fraction_refine;
  fraction (f,func);
  boundary ({f});
  FILE * fp1 = fopen ("facets1","w");
  output_facets (f,fp1);
  fclose (fp1);
  curvature (f, curv);
  FILE * fp2 = fopen ("curvature1","w");
  double xyn[2];
  xyn[0] = X0;
  xyn[1] = Y0 + L0;
  loop_interfacial_cells (fp2, f, curv, xyn);
  fclose (fp2);
 /**
 ~~~gnuplot The interface
 set size square
 plot 'facets1' w lines lw 3
 ~~~
 
  ~~~gnuplot The interface
  set xlabel 'length along circle'
 plot 'curvature1' u 4:5 w lines lw 2 t 'curvature',\
   'curvature1' u 4:2 w lines lw 2 t 'x-coordinate' ,\
   'curvature1' u 4:3 w lines lw 2 t 'y-coordinate'
 ~~~
 
 Second, the left half of the domain is coarsend. 
 */
  unrefine (x<0);
  fraction (f,func);
  boundary ({f});
  FILE * fp3 = fopen ("facets2","w");
  output_facets (f,fp3);
  fclose (fp3);
  curvature (f, curv);
  FILE * fp4 = fopen ("curvature2","w");
  loop_interfacial_cells (fp4, f, curv, xyn);
  fclose (fp4);
  /**
   ~~~gnuplot The interface
   unset xlabel
 plot 'facets2' w lines lw 3
 ~~~
 
  ~~~gnuplot Curvature and location
  set xlabel 'length along circle'
 plot 'curvature2' u 4:5 w lines lw 2 t 'curvature',\
   'curvature2' u 4:2 w lines lw 2 t 'x-coordinate' ,\
   'curvature2' u 4:3 w lines lw 2 t 'y-coordinate'
 ~~~
 
 Third, the bubble is initialized consistent with the wavelet-based adaptation algorithm.  
 */
  astats s;
  do {
    fraction (f, func);
    boundary ({f});
    s = adapt_wavelet ({f}, (double []){0.005}, maxlevel = 9);
  } while (s.nf > 15 || s.nc > 15);
  fraction (f, func);
  boundary ({f});
  FILE * fp5 = fopen ("facets3","w");
  output_facets (f,fp5);
  fclose (fp5);
  curvature (f, curv);
  FILE * fp6 = fopen ("curvature3","w");
  loop_interfacial_cells (fp6, f, curv, xyn);
  fclose (fp6);
  scalar lev[];
  foreach()
    lev[]=level;
  FILE * fp7 = fopen ("level3","w");
  loop_interfacial_cells(fp7, f, lev, xyn);
  fclose(fp7);
    /**
   ~~~gnuplot The interface
   unset xlabel
 plot 'facets3' w lines lw 3
 ~~~
 
  ~~~gnuplot Curvature, location and level of refinement
  set yr [-2 : 11]
  set xlabel 'length along circle'
  set key box
 plot 'curvature3' u 4:5 w lines lw 2 t 'curvature',\
   'curvature3' u 4:2 w lines lw 2 t 'x-coordinate' ,\
   'curvature3' u 4:3 w lines lw 2 t 'y-coordinate' ,\
      'level3' u 4:5 w lines lw 2 t 'Refinement level' 
 ~~~
 
 Finally, we compare the fidelity of the curvature estimation for the three different grids directly 
  
  ~~~gnuplot Curvature, location and level of refinement
  set yr [-1.6 : -0.4]
  set ylabel 'Curvature'
  set xlabel 'length along circle (not corresponding in locations)'
  set key box
 plot 'curvature1' u 4:5 w lines lw 3 t 'Level = 8 equidistant',\
   'curvature2' u 4:5 w lines lw 2 t 'Partially Coarsened' ,\
   'curvature3' u 4:5 w lines lw 1 t 'Max Level = 9 Adaptive grid'
 ~~~
 
 Adaptation seems to behave a bit strange?
*/
}
