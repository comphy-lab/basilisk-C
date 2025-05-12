#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "reduced.h"
#include "curvature.h" // for position()
#include "view.h"
#include "adapt_wavelet_leave_interface.h"


#define JETTIME 5.0
#define theta 0.34 // angle of the slope in radians


u.t[left] = neumann(0);
// How to set a neumann boundary condition at an angle? 
p[left] = dirichlet(0);
// value of pressure goes to zero at left boundary (why?)3aa

u.n[right] = dirichlet (-65.*(y < 2. && y > 0.1)*(t < JETTIME));
u.n[right] = dirichlet (0.);
p[right] = dirichlet(0);
// value of pressure goes to zero at right boundary (needed? why?)


#define TRIANGLE (y < -tan(theta)*x- tan(theta)*L0/2.)
#define MAXLEVEL 8
#define WATERDEPTH 1

int main()
{
  size (BOXWIDTH);
  origin(-L0); // bottom right origin
  rho2 = 0.01; // increased the density
  mu1 = 0.02;
  mu2 = 0.0004;
  G.y = -1. * cos(theta);
  G.x = -1* sin(theta);
  N = 1 << MAXLEVEL;
  run();
}
 
event init (i = 0) {
  foreach()
    f[] = y < -1.*tan(theta)*x - 0.1 ;
  boundary ({f}); // put water in the bottom half

 // refine ( TRIANGLE && level < 8.);
 // mask ( TRIANGLE ? slope : none);
}

/*
scalar triangle_slope[];
event make_slope (i++) {
 coord vc = {0.,0.}; 
  fraction (triangle_slope, TRIANGLE);
  foreach()
    foreach_dimension()
      u.x[] = triangle_slope[]*u.x[];
  boundary ((scalar *){u});
}
*/
/*
event logfile (i++) { // This log file is applicable 
 /* //fprintf (stderr, "%g\n", t);
  double ke_gas = 0., pe_gas = 0.,  ke_water = 0.,  pe_water = 0., c = 0., d = 0., e = 0., XX1 = 0., XX2 = 5., XX3 = 10.;
  foreach(reduction(+:ke_gas) reduction(+:pe_gas) reduction(+:pe_water) reduction(+:ke_water)) { // why do we need these reduction operators here?
    double r = rho(f[]);
    if (r == 0.01) {
      ke_gas += r*sq(norm(u))/2.*dv();
      pe_gas += r*G.y*y*dv();
    }
    else if (r == 1.) {
      ke_water += r*sq(norm(u))/2.*dv();
      pe_water += r*G.y*y*dv();
    }
  }



  scalar X[]; // Field of x values (strictly speaking y values)
  position (f, X, (coord){1.,0.});
  scalar Y[]; // Field of y values (strictly speaking x values)
  position (f, Y, (coord){0.,1.});
  double difference = HUGE;
  foreach(){
    if (Y[]!=nodata && fabs(X[]-XX1) <= difference){
      difference = fabs(X[]-XX1);
      c = Y[]; // Outputs the free surface height at XX1 = 0, i.e. x = 0, at the centre of the cavity
    }
    if (Y[]!=nodata && fabs(X[]-XX2) <= difference){
      difference = fabs(X[]-XX2);
      d= Y[]; // Outputs the free surface height at XX2 = 10, i.e. x = 10, away from the cavity
    }
  if (Y[]!=nodata && fabs(X[]-XX3) <= difference){
      difference = fabs(X[]-XX3);
      e= Y[]; // Outputs the free surface height at XX2 = 10, i.e. x = 10, away from the cavity
    }
  }

  fprintf (stderr, "%g %g %g %g %g %g %g %g %g\n", t, ke_gas, pe_gas, ke_water, pe_water, c, d, e, statsf(Y).min);



}
*/

 
event outputs (i++) {

 dump ();
}
 
event movie (t += 0.03; t <= 20.) {
  view (fov = 21.0715,
  ty =  (WATERDEPTH-(L0/2.))/BOXWIDTH, // centre the box regardless of water depth
  bg = {1,1,1},
  width = 1208, height = 666, samples = 4);
  box (notics = true);
  draw_vof ("f", filled = 0, fc = {1,1,1});
  squares ("u.x"); // colour according to velocity in
  mirror (n = {1,0}) {
    draw_vof ("f", filled = -1, fc = {1,1,1});
    squares ("u.y");
    box (notics = true);
  }
  save ("movie.mp4");
}
 
//Adapting by velocity (keeping refined about the interface) 
#if TREE
event adapt (i++) {
  //scalar omega[];
//vorticity (u, omega);
  adapt_wavelet_leave_interface ((scalar *){u},{f}, (double[]){0.001,0.001,0.001}, 8, 3, 1);
}
#endif



