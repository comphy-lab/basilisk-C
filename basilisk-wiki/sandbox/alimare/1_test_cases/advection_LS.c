#define BICUBIC 1
#define LevelSet     1
#define Pi 3.14159265358979323846

#include "embed.h"
#include "curvature.h"
#include "view.h"
#include "../advection_A.h"

#include "../level_set.h"
#include "../LS_advection.h"

scalar dist[];
scalar * tracers = NULL;
scalar * level_set = {dist};

#define plane(x, y, n) (y  + 0.2*sin(2.*n*Pi*x))

int main(){
  periodic(right);
  origin(-0.5*L0, -0.5*L0);
  init_grid(128);
  DT = 0.2*L0/(1 << (grid->maxdepth));

  foreach(){
     dist[] = clamp(plane(x,y,2),-0.25,0.25);
  }
  boundary ({dist});
  restriction({dist});

  vertex scalar distn[];
  cell2node(dist,distn);

  fractions (distn, cs, fs);
  fractions_cleanup(cs,fs);
  boundary({cs,fs});
  restriction({cs,fs});

/**
Loop

*/
  for(int i = 0; i<200; i++){

    double deltat  = 0.45*L0 / (1 << (grid->maxdepth));  // Delta
    int err = 0;
    int klimit = 0;

    vector vpcr[];
    foreach(){
      foreach_dimension(){
        vpcr.x[] = 0.;
      }
      if(interfacial(point,cs))vpcr.y[] = 1.;
    }
    boundary((scalar * ){vpcr});
    restriction((scalar * ){vpcr});

    scalar * speedrecons  = {vpcr.x,vpcr.y};
    recons_speed(dist, deltat, speedrecons,
     klimit, 1.e-5, &err, 
     300, overshoot = 0.3,
     cs, fs);
    draw_vof("cs");
    squares("dist", min = -0.2, max = 0.2);
    save("vpcry.mp4");
    stats s = statsf(vpcr.y);
    fprintf(stderr, "%g %g\n", s.min, s.max);
    double maxd = 0.20;
    double dtLS = timestep_LS (vpcr, DT, dist, maxd);
    RK2(dist, vpcr, dtLS, maxd);
    
    boundary ({dist});
    restriction({dist});

/**
After the advection, we need to redistance the level-set function.
*/
    LS_reinit(dist, it_max = 3, RK2 = 1);

    scalar distn[];
    cell2node(dist,distn);
    fractions (distn, cs, fs);
    fractions_cleanup(cs,fs);

    boundary({cs,fs});
    restriction({cs,fs});
  }

}

/**
Visualization of the results

![Level-set function](advection_LS/vpcry.mp4)
*/