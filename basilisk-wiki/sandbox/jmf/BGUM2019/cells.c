/**
#  Cells : static refine/unrefine
      
*/
#include "run.h"

int LEVEL = 4;

int main()
  {
    init_grid (1<<LEVEL);
    origin (-0.5, -0.5);
    run();
  }


event init (i=0)
{

output_cells (stderr);
  
  /**
  Static refine/unrefine 
  */

refine (x > 0.25 && y > 0.25 && level < (LEVEL+2));
unrefine (x < 0.25 && y < 0.25 && level >= (LEVEL-1));

  /**
  Output refined cells
  */
  
output_cells (stdout);
}

/**
## Results
 ~~~gnuplot
 plot "out" u 1:2 w lp t "Ref/unref"
 ~~~
 ~~~gnuplot
 plot "log" u 1:2 w l t "original"
 ~~~
 */
   

