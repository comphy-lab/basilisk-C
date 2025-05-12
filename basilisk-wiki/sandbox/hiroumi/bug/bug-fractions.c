/**
This bug could cause a problem when the fraction is calculated from a distance function (level-set) by using "fractions()". 
In this example, an interface defined by the distance function crosses the cell center of some cells. 
The fraction of such cells is supposed to be 0.5; however, "fractions()" gives 0 of the fraction to them. 
*/

#include "grid/quadtree.h"
#include "run.h"
#include "fractions.h"

scalar levelset[];
scalar f[];

int main ()
{
  init_grid (16);
  origin (0, 0, 0);
  L0 = 0.01;
  run();
}

event init (i = 0) {
  f.refine = f.prolongation = fraction_refine;
  boundary ({levelset,f});
  foreach() levelset[] = x+y-5e-3;

  vertex scalar phi[];
  foreach_vertex()
      phi[]=0.25*(levelset[]+levelset[-1]+levelset[0,-1]+levelset[-1,-1]);
  fractions(phi,f);
}

/**
  ![Profiles of the fraction and the interface](sandbox/hiroumi/bug/Bug-in-fractions.png)
*/
