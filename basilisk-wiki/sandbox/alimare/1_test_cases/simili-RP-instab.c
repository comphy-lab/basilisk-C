/**
# Skeleton of a Rayleigh-Plateau instability

The goal is to build the skeleton of a simili Rayleigh-Plateau instability. In
particular we are interested in the ability of retaining a skeleton even
though the VOF variable wants to pinch the interface.
*/



#include "grid/octree.h"
#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"
#include "../thinning.h"
int n_part = 0;
#include "../../Antoonvh/scatter.h"
#define QUADRATIC 1
#include "../alex_functions.h"
#include "../basic_geom.h"

#define Pi 3.141592653589793
double mini_snake(double x, double y, double z, double NB_width,
  double variation){
  coord center    = {0,0,0};
  double size     = L0/48.*(1.+variation*cos(8*Pi*z));
  return clamp(-circle(x,y,center, size) , -NB_width, NB_width);
}

int main()
{

  init_grid (64);
  origin(-0.5,-0.5,-0.5);

  /**
  We initialize the distance field on the coarse initial mesh and
  refine it adaptively until the threshold error (on distance) is
  reached. */
  scalar dist[];
  scalar f[];
  face vector s[];
  int MAXLEVEL    = 8;
  
  for (int j = 0; j < 8; j++){
    double NB_width = 2*L0/(1 << MAXLEVEL);
    double variation = 0.2+j*0.1;
    char name[80];
    foreach()
      dist[] = mini_snake(x,y,z, NB_width, variation);
    boundary({dist});
    restriction({dist});

    for (int i = 0; i < 2; i++){
        /* code */
      adapt_wavelet ({dist}, (double[]){1.e-4}, MAXLEVEL, MAXLEVEL-4);
      foreach()
        dist[] = mini_snake(x,y,z, NB_width, variation);
      boundary({dist});
      restriction({dist});
      fprintf(stderr, "loop %ld\n", grid->n );
    }
    
    NB_width = 8*L0/(1 << MAXLEVEL); // for vizualisation we extend the width
    foreach()
        dist[] = mini_snake(x,y,z, NB_width, variation);
    boundary({dist});
    restriction({dist});

    vertex scalar distn[];
    cell2node(dist,distn);


    fractions (distn, f, s);
    
    // clear();
    /**
    We display an isosurface of the distance function coloured with the
    level of refinement and the surface reconstructed from volume fractions. */

    double offset = 0.005;
    view (quat = {0.159, 0.838, 0.231, -0.468},
          fov = 10, 
          tx = -0.009, ty = 0.028,
          width = 600, height = 600);
    sprintf (name, "vof-%g.png", variation);
    draw_vof ("f");
    save (name);
    view(fov = 5);
    draw_vof ("f");
    cells (n = {1,0,0}, alpha = offset);
    sprintf (name, "vofzoom-%g.png", variation);
    save (name);

    squares ("fabs(dist[]) < 0.9*0.03125 ? dist[] : nodata", min = -0.03125, max = 0.03125, n = 
    {1,0,0}, alpha = offset);
    cells (n = {1,0,0}, alpha = offset);
    sprintf (name, "dist-%g.png", variation);
    save(name);


    scalar c[];
    foreach(){
      if(f[])c[] =1;
      else c[] = 0;
    }
    boundary({c});
    thinning3D(c);

  /**
  To output the skeleton, we will use the Lagrangian particles of Antoonvh, the
  remaining cells will be shown as spheres.
  */

    Cache skeleton = {0}; 

    foreach(){
      if(c[]){
        cache_append (&skeleton, point, 0);
        n_part++;
      }
    }
    if(n_part == 0){
      fprintf(stderr, "pas de squelette\n");
      exit(1);
    }
    coord * loc = malloc (n_part*sizeof(coord));
    int n = 0;
    foreach_cache(skeleton) {
      coord cc = {x, y, z};
      foreach_dimension()
        loc[n].x = cc.x;
      n++;
    }
    squares ("fabs(dist[]) < 0.9*0.03125 ? dist[] : nodata", min = -0.03125, max = 0.03125, n = 
    {1,0,0}, alpha = offset);
    cells (n = {1,0,0}, alpha = offset);
    scatter(loc, s = 40);
    sprintf (name, "medial_axis-%g.png", variation);
    save(name);

    free (loc);
    free (skeleton.p);
  }
}

/**
We correctly detect the medial axis of our simili-Rayleigh-Plateau instability
even when the VOF interface 

|  Fullview | VOF + mesh   |      Skeleton + distance + mesh      | 
|:-------------:|:-------------:|:-------------:|
| ![](simili-RP-instab/vof-0.2.png) | ![](simili-RP-instab/vofzoom-0.2.png) |![](simili-RP-instab/medial_axis-0.2.png)
| ![](simili-RP-instab/vof-0.3.png) | ![](simili-RP-instab/vofzoom-0.3.png) |![](simili-RP-instab/medial_axis-0.3.png)
| ![](simili-RP-instab/vof-0.4.png) | ![](simili-RP-instab/vofzoom-0.4.png) |![](simili-RP-instab/medial_axis-0.4.png)
| ![](simili-RP-instab/vof-0.5.png) | ![](simili-RP-instab/vofzoom-0.5.png) |![](simili-RP-instab/medial_axis-0.5.png)
| ![](simili-RP-instab/vof-0.6.png) | ![](simili-RP-instab/vofzoom-0.6.png) |![](simili-RP-instab/medial_axis-0.6.png)
| ![](simili-RP-instab/vof-0.7.png) | ![](simili-RP-instab/vofzoom-0.7.png) |![](simili-RP-instab/medial_axis-0.7.png)
| ![](simili-RP-instab/vof-0.8.png) | ![](simili-RP-instab/vofzoom-0.8.png) |![](simili-RP-instab/medial_axis-0.8.png)
| ![](simili-RP-instab/vof-0.9.png) | ![](simili-RP-instab/vofzoom-0.9.png) |![](simili-RP-instab/medial_axis-0.9.png)


*/