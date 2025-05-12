
/**
# Solidification with a plane interface

We investigate the solidification of a block of undercooledwater, between four plates, two at -20°C and two at +20°C. It serves at minimum working example to show how the functions defined in [elementary_body.h](/sandbox/alimare/elementary_body.h) which is dervied from : [elementary_body.h](/sandbox/qmagdelaine/phase_change/elementary_body.h). I simply added the possbility for tr_{eq} to be a field in the dirichlet_diffusion definition.

![Solification of a block of water with the temperature
field](solidification_mwe_corner/video_solidification.mp4)

We define the geometrical, temporal and resolution parameters: */

#define L 10. // size of the box

#define MIN_LEVEL 3
#define LEVEL     5
#define MAX_LEVEL 7
#define H0 L0/(1 << MIN_LEVEL)


#define T_END   10.
#define DT_MAX  L0/(1 << MAX_LEVEL)*0.8
#define DELTA_T 0.1 // for videos and measurements
#define Pi 3.141592653589793
#define NB_width L0/(1 << (LEVEL-2))// NarrowBand_width for level_set

/**
Solvers used :

Navier-Stokes
Surface tension is also modeled.
*/
#define LevelSet 1

#include "navier-stokes/centered.h"
#include "alimare/elementary_body.h"

#if LevelSet
# include "alimare/level_set.h"
#endif
/**
Level set function is used
*/

#define Ray_min               10.*L0
#define Precoeff              5.*(T_eq-TS_inf)/(SIGMA*Ray_min)

/**
We allocate several scalar fields to describe both the
interface and temperature fields. */

scalar f[], dist[];
scalar * interfaces = {f};
scalar * level_set = {dist};

/**
The main function of the program, where we set the domain geometry to
be ten times larger than the initial thickness of ice: */

int main()
{
  size (L);
  // origin (0., -L0/2.);
  origin (0., 0.);
  N  = 1 << LEVEL;
  init_grid (N);
  DT = DT_MAX;
  run();
}

/**
The initial position of the interface is defined with this function: */

// #define plane(x, y, H) (sqrt((x-L0/3.)*(x-L0/3.)/45.+(y)*(y)/15.) - H)
// #define plane(x, y, H) (sqrt((x-L0/2.)*(x-L0/2.)/350.+(y-L0/2.)*(y-L0/2.)/5.) - H )
// double hexagon (double x, double y, double h)
// {
//   double theta = atan2(x,-y);
//   double threshold1 = Pi/3.;
//   double threshold2 = 2.*Pi/3.;
//   return (theta < threshold1 ? -x+sqrt(3.*h)*(y+1.)  :
//         (theta < threshold2 ? -x+sqrt(3.*h)/2. : -x+sqrt(3.*h)*(1.-y)));
//   // return (-x+sqrt(3.)/2. );
//   // return ( x+sqrt(3.)*(y-1.) );
          
// }

double plane (double x, double y, double h)
{
  // double theta = atan2(x,-y);
  // double threshold1 = Pi/3.;
  // double threshold2 = 2.*Pi/3.;

  // return (y-fabs(sin(3.*Pi*x/L0))-h);
  return (sqrt(powf(y-L0/4.,2.)+powf(x-L0/4.,2.))-h);
  // return (sqrt(powf(y-L0/2.,2.)+powf(x-L0/2.,2.))-h);
  // return (x-h);          
}

double plane2 (double x, double y, double h)
{
  double x0 = L0/4., y0= L0/4.;
  return ((0.1+powf(x-x0,2.)+powf(y-y0,2.))*
            sqrt(powf(y-y0,2.)+powf(x-x0,2.))-h);
}


/**
Before the first step, we initialize the temperature fields (after having
refined the grid around the future interface). */

event init (i = 0) {
  scalar curve[];
#if LevelSet
  vertex scalar LS_vert[];
#endif

#if TREE
    refine (level < MAX_LEVEL && plane(x, y, (H0 - NB_width)) > 0.
            && plane(x, y, (H0 + NB_width)) < 0.);
#endif

  fraction (f, plane(x, y, H0));

  boundary({f});
  curvature (f,curve);
  boundary({curve});
#if LevelSet
  foreach_vertex(){
    LS_vert[] = plane2(x,y,H0);
  }
  boundary({LS_vert});
#endif

  foreach() {


#if LevelSet
    
// bilinear interpolation of the LS_vert field

    double delta1 = LS_vert[1,0] - LS_vert[];
    double delta2 = LS_vert[0,1] - LS_vert[];
    double delta3 = LS_vert[1,1] + LS_vert[] - (LS_vert[1,0]+LS_vert[0,1]);

    dist[] = 1./2.*(delta1 + delta2 +delta3/2.) + LS_vert[];
#endif
  }
    

  output_ppm (f, n=600, file = "init.png"\
  , min = 0., max = 1.); 
  output_ppm (dist, n=600, file = "level_set_init.png",
       min = -NB_width, max = NB_width); 

  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, n=600, min = 0, max = MAX_LEVEL, file = "level_init.png");

}

event stability (i++) {
  DT = DT_MAX;
  face vector tv[];
  foreach_face(x){
    uf.x[]    += 0.25;
  }
  foreach_face(y){
    uf.y[]    += 0.5;
  }
  double dtmax2 = DT_MAX;
  timestep (uf, dtmax2);
  boundary((scalar*){uf});
}

/**
After the *vof()* event, the evaporation velocity has to be erased. */

event tracer_advection (i++) {
  foreach_face()
    uf.x[] = 0.;
  boundary((scalar*){uf});
}

/**
## Diffusion with immersed dirichlet condition

The temperature diffuses at each timestep. We need for that the maximal
level in the simulation. */


#if TREE
event adapt (i++) {
  adapt_wavelet ({dist},\
                 (double[]){ 1e-2}, minlevel = MIN_LEVEL, \
                 maxlevel = MAX_LEVEL);
}
#endif


// Reinitialization of the LS fucntion

event LS_reinitialization(i+=50,last){
  if(i>15){
    LS_reinit2(dist,L0/(1 << MAX_LEVEL), NB_width/2.);
  }
}


/**
## Video

We now juste have to save the video.
*/



event image_finale(i = 401)
{
  boundary ({dist});
  output_ppm (dist, n=600, file = "level_set_final_reinit2.png", 
    min = -NB_width, max = NB_width); 
}


event movie (i+=5  ; i<=401)
{

  printf ( "%d %f \n", i, t);
  boundary ({f, dist});

  scalar l[];
  foreach()
    l[] = level;
  // static FILE * fp = fopen ("grid.gif", "w");
  output_ppm (l, n = 512, file="grid.gif", min = MIN_LEVEL, max = MAX_LEVEL);
  #if LevelSet
  output_ppm (dist, n = 512, linear = true, file = "LS_reinit.gif", 
    opt = "--delay 1",min = -NB_width, max = NB_width);
  #endif
}
