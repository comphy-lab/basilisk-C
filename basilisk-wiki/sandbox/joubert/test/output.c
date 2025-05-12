/**
# Test to compute variables from previous dump*/
//Seem functionnal in 2D.   

#include "utils.h"
#include "view.h"

#define WIDTH 0.27
#define Ho 0.2
#define Ha 0.21
#define MAXTIME 10.

/**
  The default maximum level of refinement is 8 and the error threshold
  on velocity is 0.01. */

int maxlevel = 8;
double uemax = 0.01;

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    uemax = atof (argv[2]);

  size (WIDTH);
  origin (-L0/2, 0, -L0/2);
  init_grid (128);

  vector sd_u[], cum_u[], cumt_u[], cum2_u[];
  scalar f1_cum[], f3_cum[], f1_cumt[], f3_cumt[];
  restore (file = "restart");
  char name[80];
  sprintf (name, "snapshot2-%g", t);
  foreach(){
    foreach_dimension(){
      cum_u.x[] = cumt_u.x[]/t;
      sd_u.x[] = sqrt(fabs((cum2_u.x[]/t - sq(cumt_u.x[]/t))));
    }
    f1_cum[] = f1_cumt[]/t;
    f3_cum[] = f3_cumt[]/t;
  }
  boundary ((scalar *) {cum_u,sd_u,f1_cum,f3_cum});
  dump (name);

  view (fov = 25.3263, quat = {0,0,0,1}, tx = 0.00531804, ty = -0.452867, bg = {1,1,1}, width = 960, height = 720);
  squares("cum_u.y",alpha=1e-6, linear=true);
  box();
  save ("cum_u.y.png");

  clear();
  box();
  squares("f1_cum",alpha=1e-6, linear=true);
  save ("f1_cum.png");

  clear();
  box();
  squares("f3_cum",alpha=1e-6, linear=true);
  save ("f3_cum.png");
}

/**

## Results

![Average y velocity.](output/cum_u.y.png)(width="800" height="600")

![Average air fraction.](output/f1_cum.png)(width="800" height="600")

![Average oil fraction.](output/f3_cum.png)(width="800" height="600")

*/
