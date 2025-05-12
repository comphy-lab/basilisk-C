/**
# Gravity Potential of a 3D Moosive star

Inspired by the seminal work of [Anders et
al. (2022)](https://arxiv.org/pdf/2204.00002.pdf), we study a moosive star.
Before one can have convection, the details of the gravity potential
should be identified. As such we load a 3D model for a moose:

![Moose model color-coded by the level of grid refinement](moose/moose.png)

And we solve the associated Poisson problem for the gravity
potential.

$$\nabla^2 \phi = \begin{cases}
1& \text{inside the star} \\
0& \text{outside the star}
\end{cases}$$

![Gravity potential on the surface of the Moosive star](moose/vof.png)

We also make a volumetric rendering of the potential

![The red blob is the gravitational centre](moose/blob.png)

The research is only a week late...
 */
#include "grid/octree.h"
#include "poisson.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"
#include "bwatch.h"

int main()
{

  /**
     The 3D model is fetched from [LordKylo's
     site](https://www.thingiverse.com/lordkylo/designs) via
     [Thingiverse](www.thingiverse.com). Kindly shared with the
     [CC-BY-SA
     lisence](https://creativecommons.org/licenses/by-sa/3.0/).
*/
  
  system ("test -f moose.stl  ||"
	  "wget https://cdn.thingiverse.com/assets/28/c8/04/94/e9/moose.stl");
  
  /**
  We read the STL file, compute the bounding box of the model and set
  the domain center and size using this bounding box. */
  
  coord * p = input_stl (fopen ("moose.stl", "r"));
  coord min, max;
  bounding_box (p, &min, &max);  
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;
  
  init_grid (8);
  size (2.*maxl);
  origin ((max.x + min.x)/2. - L0/2,
	  (max.y + min.y)/2. - L0/2,
	  (max.z + min.z)/2. - L0/2);
  
  /**
  We initialize the distance field on the coarse initial mesh and
  refine it adaptively until the threshold error (on distance) is
  reached. */

  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){5e-4*L0}, 8).nf);

	 /**
  We display an isosurface of the distance function coloured with the
  level of refinement. */
  view (fov = 20, ty = 0.55);
  box();
  cells(alpha = X0 + L0/2.5);
  isosurface ("d", 0, color = "level", min = 5, max = 10);
  save ("moose.png");

  /**
  We also compute the volume and surface fractions from the distance
  field. */

  scalar f[];
  face vector s[];
  solid (f, s, (d[] + d[-1] + d[0,-1] + d[-1,-1] +
		d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.);

  
  scalar phi[];
  foreach_dimension() {
    phi[left] = dirichlet (0);
    phi[right] = dirichlet (0);
  }
  /** 
Compute the gravity potential
*/
  foreach()
    phi[] = 0;
  poisson (phi, f);
  
  
  clear();
  view (fov = 15, ty = 0.45, bg = {0.8, 0.8, 0.8});
  draw_vof ("f", "s", color = "phi");
  save ("vof.png");
  
  scalar phi2[];
  double mval = -5;
  foreach() {
    phi2[] = phi[] < mval ? -phi[] : -1;
  }

  watch (fov = 1.5*L0, O = {-L0*cos(t), L0*sin(t), 0},
	 poi = {X0 + L0/2, Y0 + L0/2, Z0 + L0/2});
  volume (phi2, cols = true, sc = 5000, mval = mval, max = 200, shading = 1);
  store (fopen ("moose.ppm", "w"));
  system ("convert moose.ppm blob.png");
}