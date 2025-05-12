/**# a simple test of reading stl files aortic_arch.stl
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"

/**
##Importing the geometry */
void fraction_from_stl (scalar f, FILE * fp)
{
  /**
We read the STL file and compute the bounding box of the model. */
  coord * p = input_stl (fp);
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
    maxl = max.x - min.x;
    scalar d[];
  distance (d, p);
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
             d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  fractions (phi, f);
}

/**
## volume fraction field declaire */
scalar f0[];
int main()
{ 
  init_grid (16);
  size (8.);
  origin (-1*L0/4 ,-2*L0/5, 0.);

  FILE * fp = fopen ("aortic_arch.stl", "r");
  char name[10];
  char name2[10];
  // a trick loop in order to reduce the mesh number(there will be something smarter, just I can't find it, do not take it seriously)
  for (int ii=5; ii<=10; ii++){	
    fraction_from_stl (f0, fp);
    adapt_wavelet ({f0}, (double[]){0.00001}, ii);	
    // a rough check view
    sprintf (name, "lv_%d.png",ii); 
    view (fov = 19.9677, quat = {0.614306,0.390179,0.155159,0.66807}, 
          tx = -0.120804, ty = -0.340923, width = 640, height = 480);
    draw_vof ("f0");
    squares ("f0", min = 0, max = 1);
    cells();
    save (name);
    //cell number
    int cellnumber= 0;
    foreach()
      cellnumber +=1;
    fprintf(stderr,"The mesh cell number of lv %d is %d \n", ii,cellnumber);
    //dump file
    sprintf (name2, "dump_%d",ii);  	
    dump(file = name2);
  }

  /**
ALL the view value can be done by choose the view angle you want in bview3D and then save it by save("view.bv")

input and output view*/
  view (fov = 20, quat = {0,0,0,1}, tx = -0.19684, ty = 0.00981146,  bg = {0.3,0.4,0.6}, width = 600, height = 592);
  box();
  draw_vof("f0"); 	
  squares ("f0", alpha = 0.4,min = 0, max = 1);
  save ("1.png");
  /**
main chanel view*/
  clear();
  view (fov = 24, quat = {0.055611,-0.381449,-0.626467,0.677454}, tx = -0.100814, ty = -0.309643, bg = {0.3,0.4,0.6}, width = 600, height = 600, samples = 1);
  box();
  draw_vof("f0");
  squares("f0",n = {1,-1,0},alpha =2.4, min = 0, max = 1);
  squares("f0",n = {-1,-1,1.5},alpha = -3.5, min = 0, max = 1);
  save ("2.png");
  /**
cross section view*/
  clear();
  view (fov = 11.1239, quat = {0.658889,-0.274942,-0.26525,0.648006}, tx = -0.207846, ty = -0.306235, bg = {0.3,0.4,0.6}, width = 600, height = 600, samples = 1);
  box();
  draw_vof("f0");
  squares("f0",n = {-1,-1,0},alpha =-2.4, min = 0, max = 1);
  save ("3.png");
  /**
a save dumpfile*/
  dump(file="initdump");
  run();
}
