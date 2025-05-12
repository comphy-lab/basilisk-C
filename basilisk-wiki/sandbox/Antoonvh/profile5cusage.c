/**
# An example usage of profile function 5b.
*/
#include "grid/octree.h"    // Octree is 3D
//#include "grid/quadtree.h //<- Using a 2D grid gives identical results for this example
#include "profile5c.h"

scalar wave[], height[], lev[];

int main(){
  init_grid (1 << 5);
  L0 = 2.*M_PI;
  refine (level < 7 && ((y + x/3.) < 3.));
  foreach(){
    wave[] = 2.*sin(y) + 3.;
    height[] = y;
    lev[] = level;
  }
  profile ();
  profile (fname = "equi", n = 10);
}
  /**
  ## Result
  The code above prints to `stdout`, whos content are plotted below using an automated script:
  
  ~~~gnuplot Adapive
  datafile = 'out'
  firstrow = system('head -1 '.datafile)
  set yr [-0.5:8]
  set xr [-0.5:7]
  set ylabel word(firstrow, 1)
  set key autotitle columnheader box on
  set size square
  plot datafile u 2:1 t 'wave' ,\
      datafile u 3:1  t 'y',\
      datafile u 4:1  t 'level'
  ~~~
  
  The content of the file "equi" is also plotted:
  
    ~~~gnuplot Equidistant profile for 10 heights
  datafile = 'equi'
  firstrow = system('head -1 '.datafile)
  set yr [-0.5:8]
  set xr [-0.5:7]
  set ylabel word(firstrow, 1)
  set key autotitle columnheader box on
  plot datafile u 2:1 t 'wave' ,\
      datafile u 3:1 t 'y' ,\
      datafile u 4:1 t 'level' 
  ~~~
  
  We are happy with the result as it shows some minimal consistency.
  */