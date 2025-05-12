/**
# Bubble rising 3D 
(to get a precise results, it would caused a long simulation time, since we need MAX_LEVEL = 10+)
it's a simple version edit from the exemples on basilisk, intend to compared the results with axi.h:
[Bubble.c](http://basilisk.fr/src/examples/bubble.c) &
[rising.c](http://mail.basilisk.fr/src/test/rising.c)
*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "view.h"
#define Z0 0.5
#define tmax 3.
#define bubble(x,y,z)  (sq(x) + sq(y - Z0) + sq(z)) - sq(0.25)
/**
##The main function . */

int LEVEL = 7;

int MAX_LEVEL = 10;

int main() {
  //We set the domain geometry and initial refinement.
  size (2);
  origin (-L0/2., 0, -L0/2.);
  init_grid (64);
  /**
We set the physical parameters: densities, viscosities and surface tension. */
  rho1 = 1000., mu1 = 10.;
  rho2 = 1., mu2 = 0.1, f.sigma = 1.96;
  TOLERANCE = 1e-4;
  run();
}


//cylinder bc, slip laterial wall
bid tube;
uf.n[tube] =0.;
pf[top]=dirichlet(0.);
p[top]=dirichlet(0.);


//no slip flow at two side wall
u.t[top] =  dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.r[top] =  dirichlet(0.);
u.r[bottom] = dirichlet(0.);



event init (t = 0) {
  //we refine the mesh around the bubble,in order to get a better fraction field
  refine (sq(x) + sq(y - Z0) + sq(z) - sq(0.3) < 0 && level < LEVEL);
  fraction (f,  bubble(x, y, z));
  //we also set a cylinder field
  mask (sq(x)+sq(z) > sq(0.5) ? tube :none);
  view(fov = 30.65,ty= - 0.1, width = 640, height = 640);	
  clear();
  box();
  draw_vof("f", edges = true, lw = 0.5);
  save("f.png");

}

/**
We add the acceleration of gravity (unity) in the downward (-y)
direction. */
event acceleration (t < tmax;i++) {
  face vector av = a;
  foreach_face(y){
    av.y[] -= 1.;
  }
}


/**
## Outputs
Every 100 iteration, we output a fraction field. */
int j=0;
event bviewer (i += 5; t <= tmax){
  j += 1;
  char legend[1000];
  char photo[40];
  sprintf(legend, "t = %0.2g", t);
  sprintf (photo, "pic-%d", j);
  clear();
  view (fov = 20, ty = -0.5, width = 720, height = 720);
  draw_string(legend, 1, size = 30., lw = 2.);	
  squares ("f", min = 0, max = 1);
  if ((j % 100) == 0)
    save (file = photo);
  save ("field.mp4");

  clear();
  view (fov = 20, ty = -0.5, width = 720, height = 720);
  squares ("u.y", min = -0.5, max = 0.5);
  save ("axivelo.mp4");

  clear();
  view (fov = 20, ty = -0.5, width = 720, height = 720);
  squares ("u.x", min = -0.5, max = 0.5);
  save ("planvelo.mp4");
}


int k=0;
event snapshot (t += 0.5)
{
  k += 1;
  char name[40];
  char velofield[40];
  sprintf (velofield, "ve-%d", j);
  sprintf (name, "dump-%d", i);
  dump (file = name);
  output_field((scalar *){u}, fopen (velofield, "w"));
}

/**
Every ten timesteps, we output the time, volume, position, and
velocity of the bubble.*/
event loggfile (i++) {
  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0.;
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb)
          reduction(+:vbx) reduction(+:vby) reduction(+:vbz)
          reduction(+:sb)) {
    double dv = (1. - f[])*dv();
    xb += x*dv;
    yb += y*dv;
    zb += z*dv;
    vbx += u.x[]*dv;
    vby += u.y[]*dv;
    vbz += u.z[]*dv;
    sb += dv;
  }
  fprintf (ferr,
           "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
           t, sb, xb/sb, yb/sb, zb/sb,
           vbx/sb, vby/sb, vbz/sb);
  fflush (ferr);
}

/**
We adapt the mesh by controlling the error on the volume fraction and
velocity field. */
event adapt (i++) {
  double uemax = 1e-2;
  //  adapt_wavelet ({f,u}, (double[]){1e-3,uemax/10.,uemax/10.,uemax/10.}, MAX_LEVEL, 1 );
  adapt_wavelet ({f,u}, (double[]){5e-3,uemax,uemax,uemax}, LEVEL, 3);
}
 
/** 
##output


![Bubble Interface](bubble3d/field.mp4)

~~~gnuplot Rise velocity as a function of time 
set grid
set xlabel 'Time'
set key bottom right
plot  'log' u 1:7 w l t'real 3D velocity'
~~~
*/



