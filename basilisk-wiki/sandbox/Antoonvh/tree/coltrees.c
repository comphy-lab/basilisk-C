/**
# Space-colonized trees. 

On this page, a quick(er) method for volume fraction computation is
tested: Only the distance to new vertices is computed. 

![Looks good](coltrees/all.png)

*/
int maxlevel = 9;
int ja = 10000; //Number of Attraction points
double R1 = 3, H = 1.5, sy = 1.4; //Size, height of crown bottom, prolateness
double rstem = 0.125, smt = 0.1;

#include "grid/octree.h"
#include "embed.h"
#include "colonization.h"
#include "view.h"
// Seed points.
void atrp (coord * atr) {
  int j = 0;
  while (j < ja) {
    double R = R1*pow(fabs(noise()), 1./3.);
    coord vec = {noise(), noise(), noise()};
    normalize (&vec);
    if (vec.y > 0) {
      foreach_dimension()
	atr[j].x = R*vec.x;
      atr[j].y *= sy;
      atr[j].y += H;
      j++;
    }
  }
}
  
int main() {
  L0 = 1.1*(sy*R1 + H);
  X0 = Z0 = -L0/2;
  coord atr[ja];
  // Loop over density parameter
  for (double dk = 0.7; dk >= 0.25; dk -= 0.1) {
    printf ("\n#dk = %g\n", dk);
    init_grid (N);
    srand (0);
    atrp (atr); // Modified by `colonize()'
    int tn = 1;
    Tnode * nodes = calloc (tn, sizeof(Tnode));
    foreach_dimension() //Root
      nodes[0].x = 0;
    tn = colonize (atr, ja, &nodes, 1, 2, dk/1.5, dk);
    Branch * trees = calloc (499, sizeof(Branch)); //?
    nb = nodestotree (nodes, tn, &trees, rstem);
    // Construct
    tree_interface_adapt (trees, cs, fs, smooth = smt, alist = {cs},
			  crit = (double[]){0.01}, maxlevel = maxlevel,
			  ulist = {cs}, stopc = 1e4);
    view (fov = 20, ty = -0.45, width = 700, height = 700);
    draw_vof ("cs", "fs");
    char fname[99];
    sprintf (fname, "tree%g.png", 10*dk);
    save (fname);
    free (nodes); nodes = NULL;
    free (trees); trees = NULL;
  }
  system ("convert *.png -append all.png");
}
