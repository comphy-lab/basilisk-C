/**
# Plantenna

This page is part of the Plantenna project. 

![](plantenna/letters.mp4)(autoplay)
 */
#include "distance.h"
#include "view.h"
#define BVIEW 1
int n_part;
#include "../scatter.h"
#include "colonization.h"

/**
A movie-making function is defined below:.
 */
void movie_maker (Tnode * tnodes, int tn, coord * atr, int na) {
  coord * nlocs = (coord*)malloc (sizeof(coord)*tn);
  for (int j = 0; j < tn; j++) 
    foreach_dimension() 
      nlocs[j].x = tnodes[j].x;
  n_part = tn;
  scatter (nlocs, s = 10);
  save ("letters.mp4");
  free (nlocs);
}

int main() {
  init_grid (8);
  size (350);
  origin (50, 350);
  /**
Letters are read following [this test
case](/src/test/basilisk.c). 
   */
  coord * p = input_xy (fopen ("plantenna.gnu", "r"));
  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){1e-2}, 9).nf);
  boundary ({d});
  /**
The atraction points are seeded inside the letters.
   */
  int nam = 3000, na = 0;
  coord * atr = malloc (sizeof(coord)*nam);
  while (na < nam) {
    coord a = {X0 + L0*(noise() + 1)/2., Y0 + L0/2 + 50*noise()};
    if (interpolate (d, a.x, a.y) > 0)
      atr[na++] = a;
  }
  /**
Nine seed nodes are also initialized. 
   */
  int tn = 9;
  Tnode * nodes = calloc (tn, sizeof(Tnode));
  nodes[0].x = 100;
  nodes[0].y = 530;
  nodes[1].x = 130;
  nodes[1].y = 530;
  nodes[2].x = 160;
  nodes[2].y = 530;
  nodes[3].x = 200;
  nodes[3].y = 540;
  nodes[4].x = 220;
  nodes[4].y = 530;
  nodes[5].x = 240;
  nodes[5].y = 530;
  nodes[6].x = 270;
  nodes[6].y = 550;
  nodes[7].x = 300;
  nodes[7].y = 540;
  nodes[8].x = 340;
  nodes[8].y = 530;
  /**
The camera is set and the nodes are computed with a link to the
`movie_maker` function.
   */
  view (fov = 3., width = 900, height = 180,
	tx = -0.65, ty = -1.56);
  colonize (atr, na, &nodes, tn, 4.2, 0.2, 0.5, myfun = movie_maker);
  free (atr);
  free (nodes);
}
