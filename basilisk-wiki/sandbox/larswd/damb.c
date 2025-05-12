/**
Importing relevant multilayer packages
*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"

/** 
Initializing a domain with 128 cells and length 100 m
*/

int main(){
  L0 = 50;
  X0 = -10;
  N = 128;
  nl = 16;
  G = 1.0;
  CFL_H = 1.;
  run();
}

/**
Similar setup as the dam break example in [this dam break example](http://basilisk.fr/sandbox/M1EMN/Exemples/damb.c) 
Flat bottom and initially fluid at rest. The leftmost 25 meters are subjet to a 
greater layer thickness $h$. The initial $h$-field is given by

$$ h(x) = \frac{1}{n_l}\begin{cases}
                       1.1 & x < -1.0, \\
                       0.1 & x > 1.0, \\
                       1.1 - \tfrac{1}{2}(x+1.0) & -1.0 \leq x \leq 1.0.
                      \end{cases}$$
Here $n_l$ denotes the number of layers. 
*/

event init(i=0){
  u.x[left] = neumann(0);
  u.x[right] = neumann(0);
  foreach(){
    zb[] = 0.0;
    foreach_layer(){
      u.x[] = 0;
      if (x < -1){
        h[] = 1.1/nl;
      } else if ( x > 1.0){
        h[] = 0.1/nl;
      } else {
        h[] = (1.1 - 1.0*(x +1.0)/2.0)/nl;
      }
        
    }
  }
}

/**
We plot the evolution of the field at every 0.01 seconds and let the simulation run for three seconds

*/
void plot(FILE * fp){
  foreach(){
    double zc = zb[];
    foreach_layer(){
      zc += h[]/2.;
      fprintf(fp, "%g %g %g\n", x, zc, u.x[]);
      zc += h[]/2.; 
    }
    fprintf(fp, "\n");
  }
}

event plot (t <= 10; t += 0.1)
{
  char name[80];
  sprintf(name, "log%g", t);
  FILE * fp = fopen(name, "w");
  plot(fp);
  fclose(fp);
}


/**
And here we can see that the water has taken ajelly-like conistency with
no change in either surface elevation or  water velocity. 

~~~gnuplot  Initial velocity profile
reset
set size square
set xr [-10: 40]
set yr [0: 2]
set pm3d map
set grid
set xlabel 'x [m]'
set ylabel 'z [m]'
splot "log0"
~~~

~~~gnuplot  Velocity profile after 1 second
reset
set size square
set xr [-10: 40]
set yr [0: 2]
set pm3d map
set grid
set xlabel 'x [m]'
set ylabel 'z [m]'
splot "log1"
~~~

~~~gnuplot  Velocity profile after 10 seconds
reset
set size square
set xr [-10: 40]
set yr [0: 2]
set pm3d map
set grid
set xlabel 'x [m]'
set ylabel 'z [m]'
splot "log10"
~~~
*/
