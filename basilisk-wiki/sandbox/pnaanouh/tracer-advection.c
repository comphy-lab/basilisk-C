/**
# Showcase of a bug in the Tracer advection of two-phase flows

 We're using a multigrid mesh with the centred Navier-stokes solver. no-coalescence.h prevents the numerical coalescence of the droplets. The bug only occurs when the surface tension of the droplets is defined*/

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

/** we define the mesh size and simulation length */

#define TMAX 2.
#define LEVEL 7 

/** Define the tag field globaly*/
scalar tag_level[];
scalar f1[], f2[], * interfaces = {f1,f2};

int main(){
  f1.tracers={tag_level};
  f2.tracers={tag_level};
  size(3);
  init_grid(1 << (LEVEL));
  origin(-L0/2.,-L0/2.);
  const face vector muc[] = {0.1,0.1};
  mu=muc;
  f1.sigma=f2.sigma=1e4;
  foreach_dimension()
    periodic(right);
  run();
}

/** tag field is initialized to 1 for one of the drops and 2 for the other. The entire computational domain is moving up*/
event init (t=0){
  fraction(f1, -sq(x)-sq(y-0.75)+sq(0.5));
  fraction(f2, -sq(x)-sq(y+0.75)+sq(0.5));
  foreach(){
    tag_level[]=2*f1[]+f2[];
    u.x[]=0;
    u.y[]=5.;
  }
}

/** outputs a movie of the tag_level field */
event movies(t=0; t<=TMAX; t+=0.04){
  box();
  draw_vof("f1", lw = 2.);
  draw_vof("f2", lw = 2.);
  squares(color="tag_level",map=cool_warm,min=0,max=2);
  save("tag.mp4");
}
