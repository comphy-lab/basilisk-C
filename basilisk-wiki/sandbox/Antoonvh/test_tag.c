/**
## Volume fractions at the bottom boundary

We aim to find the bubbles that are attached to the bottom
boundary. Before we start, we generate a bunch of random bubbles.
 */
#include "fractions.h"
#include "utils.h"
#include "tag.h"
#include "view.h"

double xp, yp, rx, ry;
#define BUBBLE (sq(1.) - sq((x - xp)/rx) - sq((y - yp)/ry))
/**
We need a volume fraction field (`f`) and a `temp` field for various
purposes.
 */
scalar f[], temp[];

int main() {
  srand (2);
  L0 = 2.;
  X0 = Y0 = -L0/2.;
  init_grid (256);
  for (int j = 0; j < 10; j++) {
    xp = noise() - 0.2, yp = noise();
    rx = (noise() + 1.1)/4, ry = (noise() + 1.1)/4.;
    fraction (temp, BUBBLE);
    foreach()
      f[] = min(temp[] + f[], 1);
  }
  output_ppm (f, file = "f.png", n = 256);
  /**
     The bubbles look like this:  
     
     ![](test_tag/f.png)
     
Now we tag the connected bubbles;
  */
  foreach()
    temp[] = f[];
  int n = tag (temp);
  output_ppm (temp, file = "tag.png", n = 256);
  /**
The `tag()` function identified regions:

![tagged bubbles](test_tag/tag.png)

Next, lets see what regions appear at the border of interest.
   */
  int * bttm = calloc (n + 1, sizeof(int));;
  foreach_boundary(bottom)
    bttm[(int)temp[]] = 1;
# if _MPI
  MPI_Allreduce (MPI_IN_PLACE, bttm, n + 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
#endif
  /**
The cells with the relevant tag values obtain the original volume
fraction.
   */
  foreach() 
    temp[] = bttm[(int)temp[]] == 1 ? f[] : 0;
  free (bttm);
  boundary ({temp});
  
  scalar p[];
  foreach()
    p[] = pid();
  box();
  squares("p", min = 0, max = npe() - 1);
  draw_vof ("temp");
  save ("final.png");

  /**
![We can inspect the result. The colors indicate the `pid()` of the MPI domain decomposition](test_tag/final.png)
  */
}
