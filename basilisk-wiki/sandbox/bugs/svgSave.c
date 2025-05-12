/**
The save function is not working for svg file format.
*/

#include "utils.h"
#include "view.h"

int main(int argc, char const *argv[]){
  /**
  Initialisation of the grid
  */
  init_grid(1<<6);

  /**
  Definition of a scalar f to have something to plot*/
  scalar f[]; 

  foreach()
    f[] = sq(x)+sq(y)-1;

  /**
  We can dump the results, to see what the domain looks like.*/

  // dump("initial");

  clear();
  /**
  We will represent the scalar f with the squares command*/
  squares("f", min = -1, max = 1);

  /**
  We save the results in a svg file. This line is not working and give the following error:

  basilisk/src/view.h:438: redraw_feedback: Assertion `p->history' failed.

  We did get an svg file, but it is empty
  */
  save(file = "test.svg", format = "svg");

  return 0;
}