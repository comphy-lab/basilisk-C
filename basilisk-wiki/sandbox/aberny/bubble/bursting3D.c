/**
This simulation is the bursting of a bubble on a water-air interface, in 3D
*/

// #define dimension 3

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"

#include "distance.h"

#include "drop_stat.h"

#include "navier-stokes/perfs.h"
#include "profiling.h"

#include "output_stl.h"

/** 
## Physical Parameters

There are 2 parameters for this simulation. The first one is the Laplace
number. The second one is the Bond number. The default value used for the Bond
number is 0.07. It's used to compute the initial shape of the bubble.*/

// double La = 25100;
double La = 50000;
// double Bond = 0.017;
double Bond = 5;

/**
There is 2 geometrical parameters*/
double thick = 0.01;

/**
The surface tension is Sigma. We define it equal to 1.*/

#define Sigma 1.

/**
The simulation should end at t = END*/

#define END 1.

/**
If we want to have the gravitational effect, we define Gravity to 1. If not,
it should be defined at 0.*/

#define Gravity 0

#if Gravity
#include "reduced.h"
#endif

#define FILM 1

#define Bview 1

#if Bview
  #include "view.h"
#endif

/**
## Numerical Parameters

LEVEL stand for the maximum refinement of the grid

L0 is for the size of the simulation domain*/

int LEVEL = 9;
#define L0 6.

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

double intemax = 0.005;
double uemax = 0.2;

int  main(int argc, char const *argv[]) {

  size (L0);
  origin (-L0/2., -L0/2., -L0/2.);
  init_grid (1 << (6));

  /**
  By default, the value of the Laplace number is 10 000. But we can
  change that with an initial input argument. */

  if (argc >= 2){
    La = atof (argv[1]);
  }

  /**
  By default, the value of the Bond number is 0.07. But we can change that 
  with an initial input argument. */

  if (argc >=3){
    Bond = atof (argv[2]);
  }

  /**
  By default, the value of the maximum level refinement is 10. But we can change
  that with an initial input argument. */
  
  if (argc >=4){
    LEVEL = atoi (argv[3]);
  } 

  /**
  We divide by $La$ to obtain the Oh number, so $La$ must be
  positive. Concerning the $Bo$ number, we prefer to avoid numerical problems, 
  so we forced him to be positive. */

  assert(La > 0.);
  assert(Bond > 0.);

  double Oh = sqrt(1./La);

  /**
  We then define the properties of the fluids. In this case, the fluid
  1 corresponds to water, and the fluid 2 corresponds to air. */

  rho1 = 1., mu1 = Oh;
  rho2 = 1./998., mu2 = Oh/55.;
  
  /**
  The surface tension is defined with the Ohnesorge number. */

  f.sigma = Sigma;

  /**
  Last, if the gravity is defined, we set the gravity action to the Bond
  number value. */

  #if Gravity
  G.x = -Bond;
  #endif

  /**
  We print the characteristics of the fluid in the log file. */
  fprintf (stderr, "props %f %f %f %f %f %i\n \n",
     rho1, mu1, La, Bond, Sigma, LEVEL);
  run();
}

#include "findBond.h"
#include "view.h"


double cylinder(double x, double y, double z, double y0, double z0, double r) {
  double circle =  -sq(y-y0) - sq(z-z0) + sq(r);
  double plan = x+0.1;
  return min(circle, plan);
}


double pertCylinder(double x, double y, double z, 
                    double y0, double z0, double r0,
                    double eps, int omega) {
  double ri = sqrt(sq(y)+sq(z));
  double theta;
  // if (z>0)
  //   theta = atan(y/z);
  // else if (z<0)
  //   theta = atan(y/z)+pi;
  // else{
  //   if (y>0)
  //     theta = pi/2.;
  //   else
  //     theta = -pi/2.;
  //   theta = 0;
  // }

  if (y>0)
    theta = acos(z/ri);
  else
    theta = -acos(z/ri);

  double circle = ri-r0*(1+eps*cos(omega*theta));
  // double circle = r;
  double plan = x;
  return min(-circle, plan);
}

event init (t = 0) {
  if (!restore (file = "restart") ){


    /**
    Here, dataShape will contain the coordinate of the water/air interface. */

    coord* bubble;

    int* sizeBubble = NULL;
    int tailleBubble = 0;
    sizeBubble = &tailleBubble;

    fprintf(stderr, "bubble shape 2D generation\n");


    #if FILM

    bubble = bubbleOnly(Bond, thick, sizeBubble);

    fprintf(stderr, "bubble shape 3D generation\n");
    coord * shapeVolume = shape3D(bubble, tailleBubble);

    int* sizeFreeSurface = NULL;
    int  tailleSurface = 0;
    sizeFreeSurface = &tailleSurface;

    coord * topSurface;

    fprintf(stderr, "freeSurface shape 2D generation\n");
    topSurface = freeSurface(Bond, sizeFreeSurface);
    fprintf(stderr, "freeSurface shape 3D generation\n");
    coord * topVolume = shape3D(topSurface, tailleSurface);

    #else

    bubble = bubbleFree(Bond, sizeBubble);

    fprintf(stderr, "bubble shape 3D generation\n");
    coord * shapeVolume = shape3D(bubble, tailleBubble);
    #endif


    /**
    We are using the distance function, applied on the list of coordinate
    dataShape. This will generate a distance function. */

    scalar d[];
    scalar c[];
    scalar e[];

    fprintf(stderr, "distance function\n");

    #if FILM

    distance (d, shapeVolume);
    distance (c, topVolume);
    double r = thick/2.;

    // fprintf(stderr, "adapt\n");
    // {
    //   foreach(){
    //     if(fabs(c[])<5*thick && fabs(d[])<5*thick)
    //       e[] = 4000*d[];
    //     else
    //       e[] = d[];
    //   }
    //   fprintf(stderr, "toto\n");
    // }while(adapt_wavelet ({c,d,e}, (double[]){0.005, 0.005, 0.005}, LEVEL).nf);
    while(adapt_wavelet ({c,d}, (double[]){0.005, 0.0000001}, LEVEL).nf);
    fprintf(stderr, "cells %ld\n",grid->tn);   

    #else
    distance (d, shapeVolume);
    while(adapt_wavelet ({d}, (double[]){0.001}, LEVEL).nf); 
    fprintf(stderr, "cells %ld\n",grid->tn); 
    #endif  


    // dump("beforeF");

    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */

    #if FILM

    {

      vertex scalar phi[];

      foreach_vertex(){
        // phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
                // d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
        double a = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
                d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
        double b = (c[] + c[-1] + c[0,-1] + c[-1,-1] +
                c[0,0,-1] + c[-1,0,-1] + c[0,-1,-1] + c[-1,-1,-1])/8.;
        // double c = -sq (topCap.r+0.05) + sq (x - topCap.y) + sq (y - topCap.x) + sq( z - topCap.x);
        double g = min(b,a);
        double l = pertCylinder(x, y, z, 0., 0., 1.2*r, 0.2, 16);

        double geom = min(g,-l);

        phi[] = geom;
        // phi[] = g;
      }
      fprintf(stderr, "fluid phase partition\n");
      face vector s[];

      /**
      We can now initialize the volume fraction of the domain. */

      fractions (phi, f, s);

    } while (adapt_wavelet ({f}, (double[]){0.0005}, LEVEL).nf);

    #else
    {

      vertex scalar phi[];

      foreach_vertex(){
        double a = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
                d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
        phi[] = a;
      }
      fprintf(stderr, "fluid phase partition\n");
      face vector s[];

      /**
      We can now initialize the volume fraction of the domain. */

      fractions (phi, f, s);

    } while (adapt_wavelet ({f}, (double[]){0.05}, LEVEL).nf);

    #endif

    // fraction(f, h);
    
    scalar l[];
    foreach()
      l[] = level;

    dump("initial");
    stl_output_binary(f, "initial.stl");

    FILE * fp1 = fopen("interfaceTop.ppm", "w");
    clear();
    view (fov = 9.6042, quat = {0.31096,0.612435,-0.331667,0.646703}, tx = 0.00257234, ty = 0.0230582, bg = {0.3,0.4,0.6}, width = 1280, height = 720, samples = 4);

    draw_vof("f");
    save(fp = fp1);
  }
  else{
    stl_output_binary(f, "restart.stl");
    static FILE * fp1 = fopen("interfaceInitTop.ppm", "w");
      clear();
      view (fov = 9.6042, quat = {0.31096,0.612435,-0.331667,0.646703}, tx = 0.00257234, ty = 0.0230582, bg = {0.3,0.4,0.6}, width = 1280, height = 720, samples = 4);

      draw_vof("f");
      save(fp = fp1);

    static FILE * fp2 = fopen("interfaceInitBottom.ppm", "w");
      clear();
      view (fov = 9.6042, quat = {0.585635,0.457361,-0.524115,0.416122}, tx = 0.00257234, ty = 0.0230582, bg = {0.3,0.4,0.6}, width = 1280, height = 720, samples = 4);
      draw_vof("f");
      save(fp = fp2);
  }
}

/**

## Adaptation

At each step, we will adapt the mesh. We make sure that the mesh size will
be at least $1/2^7$ */

event adapt (i++) {
  adapt_wavelet ({f,u}, (double []){intemax, uemax, uemax, uemax},
     maxlevel = LEVEL);
}


int j = 0;
int jj = 257;


event outputInterface (i+=10) {
  scalar omega[];
  vorticity (u, omega);

  scalar normU[];
  foreach()
    normU[] = sqrt(sq(u.x[])+sq(u.y[]));

  char dumpFile[80];
  sprintf (dumpFile, "dump-%d", j);
  dump (file = dumpFile);
  j++;
  if (j == 5)
    j = 0;

  #if Bview
    #if _MPI
    {
      static FILE * fp1 = fopen("interfaceBottom.ppm", "w");
      clear();
      view (fov = 9.6042, quat = {0.585635,0.457361,-0.524115,0.416122}, tx = 0.00257234, ty = 0.0230582, bg = {0.3,0.4,0.6}, width = 1280, height = 720, samples = 4);
      draw_vof("f");
      save(fp = fp1);
    }

    {
      static FILE * fp1 = fopen("interfaceTop.ppm", "w");
      clear();
      view (fov = 9.6042, quat = {0.31096,0.612435,-0.331667,0.646703}, tx = 0.00257234, ty = 0.0230582, bg = {0.3,0.4,0.6}, width = 1280, height = 720, samples = 4);

      draw_vof("f");
      save(fp = fp1);
    }

    {
      static FILE * fp1 = fopen("f.ppm", "w");
      clear();
      view (fov = 12.621, quat = {0,0,-0.707107,0.707107}, tx = 0.00566754, ty = -0.0425876, bg = {0.3,0.4,0.6}, width = 1028, height = 805, samples = 1);
      squares("f", min =0, max = 1);
      save(fp = fp1);
    }
    #if FILM  
    {
      static FILE * fp1 = fopen("fZoom.ppm", "w");
      clear();
      view (fov = 3.17719, quat = {0,0,-0.707107,0.707107}, ty = -0.187812, bg = {0.3,0.4,0.6}, width = 1028, height = 805, samples = 1);
squares("f", min = 0, max = 1);
      save(fp = fp1);
    }

    {
      static FILE * fp1 = fopen("u.ppm", "w");
      clear();
      view (fov = 9.85764, quat = {0,0,-0.707107,0.707107}, tx = -0.28242, ty = -0.0396485, bg = {0.3,0.4,0.6}, width = 945, height = 862, samples = 1);
      squares("normU", min = 0, max = 60);
      save(fp = fp1);
    }
      #endif
    #else
    {
      static FILE * fp1 = popen("ppm2mp4 interfaceBottom.mp4", "w");
      clear();
      view (fov = 9.6042, quat = {0.585635,0.457361,-0.524115,0.416122}, tx = 0.00257234, ty = 0.0230582, bg = {0.3,0.4,0.6}, width = 1280, height = 720, samples = 4);
      draw_vof("f");
      save(fp = fp1);
    }

    {
      static FILE * fp1 = popen("ppm2mp4 interfaceTop.mp4", "w");
      clear();
      view (fov = 9.6042, quat = {0.31096,0.612435,-0.331667,0.646703}, tx = 0.00257234, ty = 0.0230582, bg = {0.3,0.4,0.6}, width = 1280, height = 720, samples = 4);

      draw_vof("f");
      save(fp = fp1);
    }
    #endif
  #endif

  // {
    // static FILE * fp = fopen ("interface", "w");
    // char interfaceFile[80];
    // sprintf(interfaceFile, "interface-%03d",jj);
    // static FILE * fp = fopen (interfaceFile, "w");
    // FILE * fp = fopen (interfaceFile, "w");
    // output_facets (f, fp);
    // fprintf(fp, "\n");
  // }
}

/**
The simulation evolution looks like:

<p><center>
<video width="1280" height="720" controls>
  <source src="burstingBubble/interface.mp4" type "video/mp4">
  Your browser does not support the video tag
</video>
<br>
The bubble bursting
</center></p>


We output a GFS file in order to monitor the simulation (evolution of the
mesh, local value of the velocity...), to be sure that the simulation ran
correctly. We also dump a file called "runningSimu" if we need to restart the
simulation.

The time step corresponds to approximately 10 dump during a simulation (and a
maximum of 21 outputs since we stop the simulation at $t=1$)*/

event simuInfo (i += 100; t<= END) {
  
  /**   
  It could be usefull to observe the vorticity of the flow. That's why we
  compute it.*/

  scalar omega[];
  vorticity (u, omega);

  char dumpFile[80];
  sprintf (dumpFile, "simulation-%07d", i);
  dump (file = dumpFile);

  char dumpSTL[80];
  sprintf (dumpSTL, "simulation-%07d.stl", i);
  stl_output_binary(f, dumpSTL);
}


/**
A general evolution pictures can be:

~~~gnuplot Evolution of the interface
unset key
set size ratio -1
plot [-2:2]'interface' every :10 u 2:1 w l, 'interface' u (-$2):1 w l lt 1
~~~

## Measure

measure.h includes all we need for the droplet measurement. It provides a new
structure to store the droplet data, it's initialisation and a comparison
function for the qsort function (comming from stdlib.h).

This library also contains an event for the physical measure. If we just want
the video output, we just have to set MEASURE to 0*/

#define MEASURE 1

#if MEASURE
#include "measure.h"
#endif


event logfile (i++) {
  fprintf(stderr, "%g %i %i %i\n", t, i, mgpf.i, mgpf.nrelax);
}


/**
## End

As verification, we output the final state. */

event finalState (t = end) {

  scalar m[];
  foreach()
    m[] = f[] > 1e-3;
  tag(m);

  dump(file = "final");
}
