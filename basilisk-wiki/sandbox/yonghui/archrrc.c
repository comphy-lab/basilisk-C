/**
#real 3d arch with KV model ADD THE RRC LATER*/
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "utils.h"
#include "distance.h"
#include "fractions.h"
#include "view.h"
#include "lambda2.h"

#define MIN_LEVEL 5
#define LEVEL 7
#define MAX_LEVEL 10
#define eq_r 1.
#define ADAPT 0


#define tmax   1*2.*M_PI
#define alpha_w  10.


/**
##Importing the geometry */
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
  /**
We read the STL file, compute the bounding box of the model and set the domain center and size using this bounding box. */
  init_grid (8);
  size (8.);
  origin (-1*L0/4 ,-2*L0/5, 0.);
  run();
}

/**
We then set an inflow velocity on the top and free outflow on the bottum.
axi- Y
Z - front z =L0 back z=0 
Y - top bottom 
X - left right
*/

u.n[back] = ( y >= 0)? dirichlet(10.*f0[]*cos(t)):neumann(0);
u.t[back] = ( y >= 0)? dirichlet(0):neumann(0);
u.r[back] = ( y >= 0)? dirichlet(0):neumann(0);


p[back] = ( y >= 0)? dirichlet(0):neumann(0) ;
pf[back] = ( y >= 0)? dirichlet(0):neumann(0) ;


p[front] = neumann(0) ;
pf[front] = neumann(0) ;


event init (t = 0) {
  if (!restore (file = "restart")) {
    FILE * fp = fopen ("aortic_arch.stl", "r");
    char name2[10];
    for (int ii=MIN_LEVEL; ii<=LEVEL; ii++){	
      fraction_from_stl (f0, fp);
      adapt_wavelet ({f0}, (double[]){0.00001}, ii);	
      //dump file
      sprintf (name2, "dump_%d",ii);  	
      dump(file = name2);
    }
    fraction_from_stl (f0, fp);
    view (fov = 19.9677, quat = {0.614306,0.390179,0.155159,0.66807}, 
          tx = -0.120804, ty = -0.340923, width = 1920, height = 1080);
    draw_vof ("f0");
    squares ("f0", min = -10, max =10);
    cells();
    save ("init.png");
    double viscosity =  sq(eq_r) / sq(alpha_w);
    const face vector muc[] = {viscosity, viscosity, viscosity};
    mu = muc;
  }
}


event bc (t <= tmax; i += 1) {
  //BC : lateral velocity = 0
  foreach()
    foreach_dimension()
    u.x[] = f0[]*u.x[];

    boundary ((scalar *){u , p});
  //position du test point
  double px = 3.5;
  double py = 0.8;
  double pz = 2.;
  fprintf(stderr, "%d %g %g %g \n", i, t, dt,interpolate(u.x , px, py, pz))
}


event snapshot (i += 1000)
{
  char name[80];
  sprintf (name, "dump-%d", i);
  scalar l2[];
  lambda2 (u,l2);
  dump (file = name);
}

int j = 0;
event profill (t += 0.4*M_PI)
{
  if (t > tmax - 2.*M_PI){
    j += 1;
    char name2[80];
    sprintf (name2, "test-%d", j);	
    FILE *fp2;
    fp2 = fopen (name2, "w");
    double uyy = -0.5 ;
    double uxx = 2.5;
    double velou,velov,velo;
    for (double uzz = 0.5; uzz <= 4.; uzz += 4./100. ){
      velou = interpolate(u.x , uxx, uyy, uzz);
      velov = interpolate(u.y , uxx, uyy, uzz);
      velo = velou * cos(M_PI/4.) + velov * sin(M_PI/4.);
      fprintf(fp2,"%g %g %g %g %g\n", t, uzz-2.5, velo, velou, velov);
    }
  }
}



event video (t += tmax/300.){
  scalar l2[];
  lambda2 (u,l2);
  /**
input and output view*/
  view (fov = 20, quat = {0,0,0,1}, tx = -0.19684, ty = 0.00981146,  bg = {0.3,0.4,0.6}, width = 640, height = 640);
  // draw_vof("f0"); 	
  squares ("u.z", alpha = 0.4);
  cells( alpha = 0.4);
  save ("input_uz.mp4");
  /**
main chanel view*/
  clear();
  view (fov = 18.0787, quat = {0.126228,-0.561414,-0.620277,0.533052}, tx = 0.00499754, ty = -0.373772, bg = {0.3,0.4,0.6}, width = 640, height = 640, samples = 1);
  //draw_vof("f0");
  squares("u.z",n = {1,-1,0},alpha =2.4 );
  squares("u.z",n = {-1,-1,1.5},alpha = -3.5 );
  cells(n = {1,-1,0},alpha =2.4 );
  cells(n = {-1,-1,1.5},alpha = -3.5 );
  save ("mainchannel_uz.mp4");

  clear();
  view (fov = 18.0787, quat = {0.126228,-0.561414,-0.620277,0.533052}, tx = 0.00499754, ty = -0.373772, bg = {0.3,0.4,0.6}, width = 640, height = 640, samples = 1);
  //draw_vof("f0");
  squares("l2",n = {1,-1,0},alpha =2.4 );
  squares("l2",n = {-1,-1,1.5},alpha = -3.5 );
  cells(n = {1,-1,0},alpha =2.4 );
  cells(n = {-1,-1,1.5},alpha = -3.5 );
  save ("mainchannel_l2.mp4");

  clear();
  view (fov = 18.0787, quat = {0.126228,-0.561414,-0.620277,0.533052}, tx = 0.00499754, ty = -0.373772, bg = {0.3,0.4,0.6}, width = 640, height = 640, samples = 1);
  //draw_vof("f0");
  squares("u.x",n = {1,-1,0},alpha =2.4 );
  squares("u.x",n = {-1,-1,1.5},alpha = -3.5 );
  cells(n = {1,-1,0},alpha =2.4 );
  cells(n = {-1,-1,1.5},alpha = -3.5 );
  save ("mainchannel_ux.mp4");

  /**
cross section view*/
  view (fov = 11.1239, quat = {0.658889,-0.274942,-0.26525,0.648006}, tx = -0.207846, ty = -0.306235, bg = {0.3,0.4,0.6}, width = 640, height = 640, samples = 1);
  //  draw_vof("f0");
  squares("u.x",n = {-1,-1,0},alpha =-2.4);
  cells(n = {-1,-1,0},alpha =-2.4);
  save ("section_ux.mp4");

  squares("l2",n = {-1,-1,0},alpha =-2.4);
  cells(n = {-1,-1,0},alpha =-2.4);
  save ("section_l2.mp4");

}

#if ADAPT
event adapt (i++) {
  double uemax = 0.01;
  adapt_wavelet ({u}, (double[]){uemax,uemax,uemax}, MAX_LEVEL, 4);
}
#endif











