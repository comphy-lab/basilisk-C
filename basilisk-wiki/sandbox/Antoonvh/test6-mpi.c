/**
# An example of the functions under `profile6.h`
 */
#include "profile6.h"
scalar gaussian[], wave[], lev[];

/**
We define a domain-centered radial coordinate:
 */
#define RADIUS sqrt(sq(x - L0/2.) + sq(y - L0/2.))

int main(){
  init_grid (1 << 5);
  L0 = 2.*M_PI;
  refine (level < 7 && ((y + x/3.) < 3.));
  foreach(){
    wave[] = 2.*sin(y) + 3.;
    lev[] = level;
    gaussian[] = 2*exp(-sq(RADIUS)/4.) + 3.;
  }
  /**
   The `profiles()` function replaces the trusty `profile()` function
   from [`profile5c.c`](profile5c.c) in the default mode.
   */
  profiles();
  /**
The results:

~~~gnuplot Adapive
  datafile = 'out'
  firstrow = system('head -1 '.datafile)
  
  set ylabel word(firstrow, 1)
  set key autotitle columnheader box on
  set size square
  plot datafile u 2:1 t 'Gaussian' ,\
      datafile u 3:1  t 'wave' ,\
      datafile u 4:1  t 'level'
  ~~~
    
  We also average two fields, 40 times along the radial coordinate
  (`RADIUS`) and write it to a file named ``equi`'.
  */
  scalar * list = {gaussian, lev};
  profile_equi (list, RADIUS, 40, "equi");

  /**
The content of the file is also plotted:
  
~~~gnuplot Equidistant profile for 10 heights
  datafile = 'equi'
  firstrow = system('head -1 '.datafile)
  
  unset ylabel
  set xlabel word(firstrow, 1)
  set key autotitle columnheader box on
  plot datafile u 1:2 w l ,\
      datafile u 1:3   
     
  ~~~

Well done functions under [`profile6.h`](profile6.h)!
  */
}