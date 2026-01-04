

/** I want to understand how the domain is traversed in the foreach statement. In the case of quadtree, a Z-curve is used (see https://basilisk.fr/sandbox/twitkamp/foreach.c) . I want to check how it is done with the multigrid solver*/

#include "grid/multigrid.h"
#include "run.h"

int main(){
  for (int size = 4; size <= 16; size *= 2){ 
    N = size;
    L0 = N;
    run();
  }
}

/** Output the position of the cell */

event init(t = 0){
  char file[80];
  sprintf(file, "foreach_N_%d", N);
  FILE * fp =  fopen(file, "w");

  foreach(){
    fprintf(fp, "x, y, %g, %g\n", x, y);
  }
}

/**
## Results

Using the multigrid solver, the domain is traversed in (reversed ?) column-major order

~~~gnuplot 
set xrange [0:4]
set yrange [0:4]
set xlabel 'x'
set ylabel 'y'
set style line 11 lc rgb 'black' lt 1 lw 2
set xtics 1
set ytics 1
unset mxtics
unset mytics
stats 'foreach_N_4' using 3:4 nooutput
first_x = STATS_min_x
first_y = STATS_min_y
set label "start" at 1.1*first_x, first_y left

set grid xtics ytics ls 11

plot 'foreach_N_4' u 3:4 w l  notitle, \
     '' using 3:4 every ::0::0 with points pt 7 ps 1.5 lc rgb 'red' notitle
~~~

~~~gnuplot 
set xrange [0:8]
set yrange [0:8]
set xlabel 'x'
set ylabel 'y'
set style line 11 lc rgb 'black' lt 1 lw 2
set xtics 1
set ytics 1
unset mxtics
unset mytics
stats 'foreach_N_8' using 3:4 nooutput
first_x = STATS_min_x
first_y = STATS_min_y
set label "start" at 1.1*first_x, first_y left

set grid xtics ytics ls 11

plot 'foreach_N_8' u 3:4 w l  notitle, \
     '' using 3:4 every ::0::0 with points pt 7 ps 1.5 lc rgb 'red' notitle
~~~

~~~gnuplot 
set xrange [0:16]
set yrange [0:16]
set xlabel 'x'
set ylabel 'y'
set style line 11 lc rgb 'black' lt 1 lw 2
set xtics 1
set ytics 1
unset mxtics
unset mytics
stats 'foreach_N_16' using 3:4 nooutput
first_x = STATS_min_x
first_y = STATS_min_y
set label "start" at 1.1*first_x, first_y left

set grid xtics ytics ls 11

plot 'foreach_N_16' u 3:4 w l  notitle, \
     '' using 3:4 every ::0::0 with points pt 7 ps 1.5 lc rgb 'red' notitle
~~~
*/