/**
# Multiple branches with automatically scaling refinement.

There are 3 branches with a factor of four difference in radius. With
$\Delta_{min}R^{-1} \propto Re^{-0.5}$, this entails that the resolution
should be doubled between each branch. We accomplish this:

![Each branch is at a different level](branchrefine/split.png)

![With Smoothening](branchrefine/splitsmt.png)

 */
#include "grid/octree.h"
#include "embed.h"
#include "treegen.h"
#include "view.h"

double U = 1., nu = 1./100., fillval = 1e5;
void prolongate_ratio (Point point, scalar s) {
  foreach_child() {
    if (s[] != fillval)
      s[] += s[]*Delta;
  }
}

int main () {
  /**
     We have three adhoc branches;
  */
  nb = 3;
  Branch branches[nb];
  branches[0].start = (coord){0,-10,0};
  branches[0].end = (coord){0,10,0};
  branches[0].R = 1;
  
  branches[1].start = (coord){1,0,0};
  branches[1].end = (coord){10,0,0};
  branches[1].R = 0.25 ;
  
  branches[2].start = (coord){3,0.25,0};
  branches[2].end = (coord){5,2.5,0};
  branches[2].R = 0.06125;

  
  L0 = 10;
  X0 = Y0 = Z0 = -L0/2.;
  init_grid (32);

  /**
     The `adapt_wavelet`-compatible algorithm is defined below.
   */
  scalar J[], res[];
  res.prolongation = prolongate_ratio;
  
  double smt = 0.5;
  do {
    tree_interface (branches, cs, fs, J, smt);
    foreach() {
      res[] = fillval;
      if (cs[] > 0 && cs[] < 1) {
	res[] = U/sqrt(branches[(int)(J[] + 0.5)].R*nu);
      }
    }
    boundary ({res});
  } while (adapt_wavelet ({res}, (double[]){2.0}, maxlevel = 99).nf > 0);
  
  tree_interface (branches, cs, fs, J, smt);
  fractions_cleanup (cs, fs);
  view (fov = 13, width = 800, height = 800, tx = -0.2);
  draw_vof("cs", "fs", edges = true, color = "J", min = -1, max = 4);
  cells();
  save ("splitsmt.png");
  double fillval = 1e5;
  //Again, no smoothening
  
  do {
    tree_interface (branches, cs, fs, J);
    foreach() {
      res[] = fillval;
      if (cs[] > 0 && cs[] < 1) {
	res[] = U/sqrt(branches[(int)(J[] + 0.5)].R*nu);
      }
    }
    boundary ({res});
  } while (adapt_wavelet ({res}, (double[]){2.0}, maxlevel = 99).nf > 0);
  
  tree_interface (branches, cs, fs, J);
  fractions_cleanup (cs, fs);
  draw_vof("cs", "fs", edges = true, color = "J", min = -1, max = 4);
  cells();
  save ("split.png");

}
  
  
