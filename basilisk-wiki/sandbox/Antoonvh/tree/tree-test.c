/**
# A tree on a tree

![A fractal tree; The interface in color coded with the local branch
 radius.](tree-test/tree.mp4)
 */

#include "grid/octree.h"
#include "treegen.h"
#include "view.h"

Branch * trees;

int maxlevel = 9;

int main() {
  levels = 4;
  srand (time(NULL));
  trees = tree_skeleton();
  L0 = 55.;
  X0 = Z0 = -L0/2;
  N = 64;
  init_grid (N);
  scalar cs[];
  scalar J[]; //Index of closes branch
  cs.prolongation = cs.refine = fraction_refine;
  face vector fs[];
  int it = maxlevel - depth(), j = 0;
  while (j < it + 2) {
    printf ("%d %ld\n", j++, grid->n);
    tree_interface (trees, cs, fs);
    boundary ({cs});
    adapt_wavelet ({cs}, (double[]){0.01}, maxlevel);
  }
  tree_interface (trees, cs, fs, J); //for `fs` and `J`
  foreach() 
    J[] = trees[(int)(J[] + 0.5)].R; //replace index with radius
  boundary ({J});
  for (double theta = 0; theta <= 2*pi; theta += 0.025) {
    view (ty = -0.5, theta = theta);
    draw_vof ("cs", "fs", color = "J", min = 0, max = trees[0].R);
    cells();
    save ("tree.mp4");
  }
  free (trees);
}

