/**
This simulation is the bursting of a bubble on a water-air interface, in 3D
*/

#define dimension 3

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"

// #include "distance.h"
#include "view.h"

#define TETRA 1

/**
## Numerical Parameters

LEVEL stand for the maximum refinement of the grid

L0 is for the size of the simulation domain*/

int LEVEL = 10;
#define L0 4.


int  main(int argc, char const *argv[]) {

  size (L0);
  // origin (-L0/2.,-L0/2.,-L0/2.);
  init_grid (1 << (4));

  run();
}

#if TETRA
/**
The coordinate of my tetrahedron*/

coord A = {0.1, 0.1, 0.1};
coord B = {0.9, 0.1, 0.1};
coord C = {0.1, 0.9, 0.1};
coord D = {0.1, 0.1, 0.9};
#else
/**
The coordinate of 4 points of my plane*/
coord A = {0, 0, 0.5};
coord B = {0.8, 0, 0};
coord C = {0.8, 1, 0};
coord D = {0, 1, 0.5};
#endif

/** Generation of the array with the triangles coordinate. This function is
based on the input_stl function*/

coord * inputTest(void) {
  Array * h = array_new();

  /**
  First triangle*/
  {  
    coord toto = {A.x, A.y, A.z};
    array_append(h, &toto, sizeof(coord));
  }
  {  
    coord toto = {B.x, B.y, B.z};
    array_append(h, &toto, sizeof(coord));
  }
  {
    coord toto = {C.x, C.y, C.z};
    array_append(h, &toto, sizeof(coord));
  }

  /**
  Second triangle*/
  {
    coord toto = {A.x, A.y, A.z};
    array_append(h, &toto, sizeof(coord));
  }
  {
    coord toto = {C.x, C.y, C.z};
    array_append(h, &toto, sizeof(coord));
  }
  {
    coord toto = {D.x, D.y, D.z};
    array_append(h, &toto, sizeof(coord));
  }

#if TETRA

  /**
  Instead of a plane, we can generate a tetrahedron

  Third triangle*/
  { 
    coord toto = {A.x, A.y, A.z};
    array_append(h, &toto, sizeof(coord));
  }
  {
    coord toto = {D.x, D.y, D.z};
    array_append(h, &toto, sizeof(coord));
  }
  {
    coord toto = {B.x, B.y, B.z};
    array_append(h, &toto, sizeof(coord));
  }

  /**
  Fourth triangle*/

  {
    coord toto = {B.x, B.y, B.z};
    array_append(h, &toto, sizeof(coord));
  }
  {
    coord toto = {D.x, D.y, D.z};
    array_append(h, &toto, sizeof(coord));
  }
  {
    coord toto = {C.x, C.y, C.z};
    array_append(h, &toto, sizeof(coord));
  }
#endif
  coord p = {nodata};
  array_append(h, &p, sizeof(coord));

  return (coord *) array_shrink(h);
}

/** We can play with Antoon tric to display an isosurface. We base it on the
distance field obtain with its new function distance_to_surface.*/

#include "../Antoonvh/frac-dist.h"
#include "output_stl.h"


event init (t = 0) {

  coord * h = inputTest();

  scalar d[];
  distance (d, h);
  while (adapt_wavelet ({d}, (double[]){0.01}, LEVEL).nf);
  vertex scalar phi[];

  foreach_vertex()
    phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1] +
       d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  face vector s[];
  fractions (phi, f, s);

  scalar dist[];
  distance_to_surface(c = f, d = dist);

  dump("initial");
  clear();
  view (fov = 33.2344, quat = {0.206574,0.0297392,-0.0399921,0.977161}, tx = -0.370845, ty = -0.692088, bg = {0.3,0.4,0.6}, width = 600, height = 600, samples = 1);
  draw_vof("f");
  box();
  save(file = "vof.png");

  stl_output_binary(f, "3D.stl");
  // output_stl_binary(f);

  clear();
  view (fov = 28.5561, quat = {0.0297249,-0.905678,-0.289561,0.308253},
        tx = 0.148153, ty = -0.393888, bg = {0.3,0.4,0.6}, 
        width = 600, height = 600, samples = 4);
  isosurface("dist", 0);   
  save("isosurface.png");
  
  clear();
  view (fov = 40.2904, quat = {-0.260473,0.358927,0.116761,0.888649}, 
        tx = -0.0447443, ty = 0.0754991, bg = {0.3,0.4,0.6},
        width = 600, height = 600, samples = 4);
  squares("dist", min = -0.1, max = 0.1);
  box();
  squares("dist", min = -0.1, max = 0.1, n = {1,0,0});
  squares("dist", min = -0.1, max = 0.1, n = {0,1,0});
  cells();
  // cells(n={1,0,0});
  // cells(n={0,1,0});
  save("dist.png");
}

/**
Our interface looks like:

  ![Reconstructed VOF surface.](test3D/vof.png)

We can observe the isosurface:

  ![isosurface](test3D/isosurface.png)

We can also observe the distance field to reconstruct the isosurface on the box

  ![distance field](test3D/dist.png)
*/