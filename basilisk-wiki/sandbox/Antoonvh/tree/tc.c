/**
# Example of the space colonization algorithm

![Growing branches kill attraction points](tc/mov.mp4)

## howto
First the relevant header files are included
 */
#include "grid/octree.h"
#include "view.h"
#define BVIEW 1
int n_part;
#include "../scatter.h"
#include "colonization.h"
/**
A function is defined that plots the growing root and the remaining attraction points.
 */
trace
void movie_maker (Tnode * tnodes, int tn, coord * atr, int na) {
  Branch * trees = calloc (499, sizeof(Branch)); //?
  nb = nodestotree (tnodes, tn, &trees);
  scalar cs[], R[];
  cs.prolongation = cs.refine = fraction_refine;
  face vector fs[];
  tree_interface (trees, cs, fs, R); //This maybe expensive
  boundary({cs});
  
  view (fov = 31, psi = -pi/4, bg = {0.6, 0.45, 0.33});
#if dimension == 2
  draw_vof ("cs", "fs", filled = -1, fc = {0.8,0.8,0.8});
#else
  view (quat = {0.433278,0.701835,0.0938225,0.557577}, ty = -.3);
  draw_vof ("cs", "fs", fc = {0.8,0.8,0.8});
#endif
  n_part = na;
  scatter (atr);
  save ("mov.mp4");
  adapt_wavelet ({cs}, (double[]){0.01}, 9, 5);
  free (trees);
  dump();
}
/**
We setup a scenario with a single seed node at the edge of the domain
that is attracted to 500 attraction points, randomly placed in a
piramid.
*/
int main() {
  init_grid (128);
  L0 = 2.1;
  X0 = Y0 = Z0 = -L0/2;
  int tn = 1;
  Tnode * nodes = calloc (tn, sizeof(Tnode));
  foreach_dimension()
    nodes[0].x = 1;
  int na = 250;
  coord atr[na];
  for (int k = 0; k < na; k++) {
    coord ap;
    do {
      ap = (coord){noise(), noise(), noise()};
    } while (ap.x + ap.y + ap.z < 0);
    atr[k] = ap;
  }
      
#if _MPI
  MPI_Bcast (&atr[0], na*sizeof(coord), MPI_BYTE, 0, MPI_COMM_WORLD); 
#endif
  /**
     The nodes are computed and the `movie_maker` function is
     called each iteration.
   */
  tn =  colonize (atr, na, &nodes, tn, 0.5, .025, 0.1,
	    myfun = movie_maker);
  free (nodes);
}
