/**
#Inertial matrix calculation test case

*/

#include "grid/octree.h"
#define QUADRATIC 1
#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))
#include "../alex_functions.h"
#include "../basic_geom.h"
#include "distance.h"
#include "../inertial.h"
#include "utils.h"
int main(){

  /**
  We read the STL file, compute the bounding box of the model and set
  the domain center and size using this bounding box. */
  
  // coord * p = input_stl (fopen ("windTurbineGeom.stl", "r"));


  // coord min, max;
  // bounding_box (p, &min, &max);  
  // double maxl = -HUGE;
  // foreach_dimension()
  //   if (max.x - min.x > maxl)
  //     maxl = max.x - min.x;
  
  
  // size (1.2*maxl);
  // origin ((max.x + min.x)/2. - L0/2,
  //   (max.y + min.y)/2. - L0/2,
  //   (max.z + min.z)/2. - L0/2);

  // /**
  // We initialize the distance field on the coarse initial mesh and
  // refine it adaptively until the threshold error (on distance) is
  // reached. */

  // scalar d[];
  // distance (d, p);
  // while (adapt_wavelet ({d}, (double[]){5e-4*L0}, 8).nf);

  init_grid (32);
  origin(-L0/2,-L0/2,-L0/2);

  scalar dist[];
  scalar cs[];
  face vector fs[];
  // for (int i = 0; i < 6; ++i)
  // {
    double R_init = 0.2;
    coord center = {0.1,0.00,0.};
    foreach(){
        dist[] = sphere(x,y,z,center,R_init);
    }
    boundary ({dist});
    restriction({dist});


    vertex scalar dist_n[];
    cell2node(dist,dist_n);

    fractions (dist_n, cs, fs);
    boundary({cs,fs});
    restriction({cs,fs});
    // adapt_wavelet({cs},(double[]){1.e-2},8);
    // fprintf(stderr, "%ld\n",grid->tn );
  // }
  dump();
  // exit(1);
  
  coord cogP,refVertex;
  double volMom2[6];

  refVertex.x = 0.;
  refVertex.y = 0.;
  refVertex.z = 0.;
  double vol = 0.;
  polyhedronMom(&cogP,volMom2,refVertex,&vol,cs,fs);
  fprintf(stderr, "CALCULATED\n");
  fprintf(stderr, "%g %g %g %g\n", vol, cogP.x,cogP.y,cogP.z);
  fprintf(stderr, "THEORY\n");
  double voltheo = 4./3*M_PI*R_init*R_init*R_init;
  fprintf(stderr, "%g %g %g %g\n", voltheo,
    center.x,center.y,center.z);
  exit(1);
}