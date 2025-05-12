
/**

In this script, we present the core concept of AMR (Adaptive Mesh Refinement) through an analogy with pictures where we consider that an element of the mesh corresponds to one pixel of a picture. 

*/


#include "utils.h"
#include "input.h"
#include "run.h"

/**
Below: for .vtk output. The .vtk files have no interest to be computed on the Basilisk website, so I commented the concerned lines. I didn't erase them, as they can be usefull as an example for whose who want to use these output in their simulations
*/

/* #include "./vtknew_cell.h"  */


int main(){
    init_grid(1<<8);
    run();
}

event tend(i=10){}


/**

AMR
---

We import a gray-scale version of the Basilisk logo and, in each cell, we put one value corresponding to the gray value.

*/

event movie(i++){
      
    scalar gris[];
    FILE * fp = fopen("./logo_ascii.h","r");  // import the (modified and grayed) basilisk logo
    input_pgm(gris,fp);
    fclose(fp);

    foreach()
      if (y>0.82 && gris[]==0.)
        gris[]=1.;

/**
Once the values are stored in the mesh elements, we perform AMR (i.e. we change the element size repartition in the domain).
Each new cell can be finer or coarser than the original cells, and still contains one value of gray.
*/

    
    // AMR
    for (int il=0;il<=12;il++){
      int k = 100-10*i;
      double eps=0.005+(0.5-0.005)*cube(1.*k/100.);
      adapt_wavelet({gris},(double[]){eps},maxlevel=8);
    }

/**

Commented lines below: .vtk outputs.

*/
    
    /* scalar lev[]; */
    /* foreach() */
    /* 	 lev[]=level; */
    /* boundary({lev}); */

    /* char fname[80]; */
    /* sprintf(fname,"fields_%i.vtk",k); */
    /* fp = fopen(fname,"w"); */
    /* output_vtk((scalar*){lev,gris},fp); */
    /* fclose(fp); */


/**
We output two video: one with the gray value in each element, and one with the level of the elements (level of element directly linked to element size).
*/
    
    output_ppm (gris, file = "gris.mp4", box = {{0.,0.},{1.,1.}}, linear = true, spread = 2, n = 256); // it would be nicer with a gray scale and a reduced framerate

    scalar l[];
    foreach()
      l[] = level;

    output_ppm (l, file = "level.mp4", box = {{0.,0.},{1.,1.}}, linear = false, min = 3, max = 8, n = 256); // it would be nicer a reduced framerate 

}


/**

We show the two video.
We see on the starting coarse mesh (mesh containing big elements), that the corresponding picture is underresolved and fuzzy.
As soon as the mesh element are refined, the the logo becomes sufficiently resolved.

   ![Evolution of the gray level](AMRlogo/gris.mp4)

   ![Evolution of the level of refinement](AMRlogo/level.mp4)

If we imagine that the logo is the solution of an EDP, we can say that the solution is not correct on the coarse mesh (the picture is underrsolved, which corresponds to a solution containing huge errors), whereas on the finer meshes, the solution seems correct (small errors).
On the opposite, finer meshes implies more elements, and thus more expensive computations.
In that perspective, doing AMR means searching an optimized mesh (optimized cell size repartition), which minimizes the error for a given number of element, or which minimizes the number of element for a given error tolerance. 

*/








