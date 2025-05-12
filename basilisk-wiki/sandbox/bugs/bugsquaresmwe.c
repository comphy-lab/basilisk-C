/**
# Squares() bug

Unexpected syntax error when there is an underscore in the squares() function of view. */

#include "fractions.h"
#include "view.h"

int main() 
{
  origin (-0.5,-0.5,-0.5);
  init_grid (128);
  scalar LSview[], LS_view[];

  foreach(){
    LS_view[] = sq(x) + sq(y) + sq(z) - sq(0.3);  
    LSview[] = LS_view[];
  }
  scalar f[];
  fraction (f, sq(x) + sq(y) + sq(z) - sq(0.3));
  draw_vof ("f");
  squares ("LSview");
  save ("LSview.png");

  draw_vof ("f");
  squares ("LS_view");
  save ("LS_view.png");
}

/**
![Working](bugsquaresmwe/LSview.png)

![Not working](bugsquaresmwe/LS_view.png)
*/
