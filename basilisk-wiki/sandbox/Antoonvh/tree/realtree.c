/**
# A "realistic" tree in Basilisk

Based on
[in-house](https://www.tudelft.nl/en/ceg/about-faculty/departments/geoscience-remote-sensing/)
laser-scanner data (in collaboration with [Adriaan van
Natijne](https://www.tudelft.nl/citg/over-faculteit/afdelingen/geoscience-remote-sensing/staff/phd-students/ir-al-adriaan-van-natijne/))
we aim to model a tree using the [QSM tree model
code](https://github.com/InverseTampere/TreeQSM/). We check how a tree
gets represented as an embedded boundary in Basilisk. 

A Tree resonstruction with branches color-coded according to their
radius. I did not realize TU Delft's Mekel Park is so ugly:

<img src="realtree/tree.png" alt="drawing" width=70%/>

There is room for improvement.
 */

#include "grid/octree.h"
#include "embed.h"
#include "treegen.h"
#include "view.h"

scalar J[];
int maxlevel = 10;
#define SIZEOFBRANCH (60) //!= sizeof(branch)
int main () {
  Branch * branches;
  FILE * fp = fopen ("realtree.branches", "rb");
  fseek (fp, 0, SEEK_END); 
  int size = ftell (fp); // get file size
  fseek (fp, 0, SEEK_SET);
  assert (!(size % SIZEOFBRANCH)); //check filesize
  nb = size/SIZEOFBRANCH;
  printf ("The tree consists of %d segments\n", nb);
  branches = (Branch*) malloc (sizeof(Branch)*nb); 
  for (int i = 0; i < nb; i++){
    Branch branch;
    fread (&branch, SIZEOFBRANCH, 1, fp);
    coord trans = {50, 19, 0};
    foreach_dimension() {
      branch.start.x -= trans.x;
      branch.end.x -= trans.x;
    }
    branches[i] = branch;
  }
  printf ("Read the branch data: 100%%\n");
  L0 = 20;
  X0 = -5.;
  Y0 = -10;
  Z0 = 0;
  init_grid (128);
 printf ("Reconstruction..\nIt.\tCells\tMax. level\n");
 for (int i = 0; i <= 5 ; i++) {
   tree_interface (branches, cs, fs);
   boundary ({cs});
   adapt_wavelet ({cs}, (double[]){0.001}, maxlevel);
   printf ("%d\t%ld\t%d\n", i , grid->n, grid->maxdepth);
 }
 tree_interface (branches, cs, fs, J);
 fractions_cleanup (cs, fs);
 
 printf ("Generating image\n");
 foreach()
   J[] = branches[(int)(J[] + 0.5)].R;
 view (fov = 24, quat = {0.632164,0.374264,0.346885,0.583067},
       tx = 0.0840115, ty = -0.372052, bg = {0.3,0.4,0.6},
       width = 1600, height = 1600, samples = 2);
 draw_vof ("cs", "fs", color = "J");
 box();
 save ("tree.png");
 free (branches);
}
