/**
This simulation is the bursting of a bubble on a water-air interface, in plane 2D
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"

#include "distance.h"

#include "drop_stat.h"

#include "navier-stokes/perfs.h"
#include "profiling.h"

/** 
## Physical Parameters

There are 2 parameters for this simulation. The first one is the Laplace
number. The second one is the Bond number. The default value used for the Bond
number is 0.07. It's used to compute the initial shape of the bubble.*/

double La = 20000;
double Bond = 1;


/**
There is 2 geometrical parameters*/
double thick = 0.01;

/**
The surface tension is Sigma. We define it equal to 1.*/

#define Sigma 1.

/**
The simulation should end at t = END*/

#define END 0.1

/**
If we want to have the gravitational effect, we define Gravity to 1. If not,
it should be defined at 0.*/

#define Gravity 0

#define FILM 1

#if Gravity
#include "reduced.h"
#endif

/**
## Numerical Parameters

LEVEL stand for the maximum refinement of the grid

L0 is for the size of the simulation domain*/

int LEVEL = 11;
#define L0 4.

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

double intemax = 0.005;
double uemax = 0.2;

int  main(int argc, char const *argv[]) {

  size (L0);
  origin (-L0/2.,0);
  init_grid (1 << (5));

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


double hole(double x, double y, double r) {
  double line1 = y - r;
  double line2 = x;
  return -min(-line1, line2);
}

event init (t = 0) {
  if (!restore (file = "restart")){


    /**
    Here, dataShape will contain the coordinate of the water/air interface. */

    coord* bubble;

    int* sizeBubble = NULL;
    int tailleBubble = 0;
    sizeBubble = &tailleBubble;

    bubble = bubbleOnly(Bond, thick, sizeBubble);

    // coord * shapeVolume = input_stl(fopen("interface.stl", "r"));

    // fprintf(stderr, "x= %f y= %f r= %f\n",fillet.x, fillet.y, fillet.r );

    int* sizeFreeSurface = NULL;
    int  tailleSurface = 0;
    sizeFreeSurface = &tailleSurface;

    coord * topSurface;
    topSurface = freeSurface(Bond, sizeFreeSurface);

    coord * capAllone;
    capAllone = capOnly(Bond, thick);

    /**
    We are using the distance function, applied on the list of coordinate
    dataShape. This will generate a distance function. */

    scalar c[];
    scalar d[];
    scalar e[];
    scalar g[];

    distance (d, bubble);
    distance (c, topSurface);
    // distance (toto, capAllone);

    

    double r = thick/2.;
    // {
    //   foreach(){
    //     e[] = c[]*d[];
    //     if (fabs(e[])<0.01 && x>0.1 && y<1.5) {
    //       g[] = 10*fabs(e[]);
    //     }
    //     else
    //         g[] = 0.1;
    //   }
    //   boundary({e});
    //   boundary({g});
    // }while(adapt_wavelet ({c,d,g}, (double[]){0.005, 0.005, 0.0001}, LEVEL).nf);
    
    while(adapt_wavelet ({c,d}, (double[]){0.01, 0.0000001}, LEVEL).nf); 
    fprintf(stderr, "cells %ld\n",grid->tn);  
    dump(file = "cells"); 

    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */

    {

      vertex scalar phi[];

      foreach_vertex(){
        // phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
                // d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
        double a = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
        double b = (c[] + c[-1] + c[0,-1] + c[-1,-1])/4.;
        // double c = -sq (topCap.r+0.05) + sq (x - topCap.y) + sq (y - topCap.x) + sq( z - topCap.x);
        double g = min(-b,-a);
        double l = hole(x, y, r);

        double geom = min(g,l);

        phi[] = geom;
        // phi[] = l;
      }

      /**
      We can now initialize the volume fraction of the domain. */

      fractions (phi, f);

    } while (adapt_wavelet ({f}, (double[]){0.0005}, LEVEL).nf);

    // fraction(f, h);
    
    scalar l[];
    foreach()
      l[] = level;

    dump("initial");

    // clear();
    // view (fov = 4.69417, quat = {0,0,-0.707108,0.707106}, tx = 0., ty = 0.0855609, bg = {1,1,1}, width = 1280, height = 720, samples = 4);
    // draw_vof("f", filled = 1);
    // mirror({0,1,0}, 0)
    //   draw_vof("f", filled = 1);
    // save("initial.ppm");
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

/**
## General output

We output 2 files. One will display the interface to observe the evolution of 
the simulation. The second is the evolution of the mesh, in order to check the 
non-saturation of the mesh in the simulation.*/

int jj = 0;

event outputInterface (t+=0.0005) {
  
  scalar omega[];
  vorticity (u, omega);
  
  {
    // static FILE * fp = fopen ("interface", "w");
    jj++;
    char interfaceFile[80];
    sprintf(interfaceFile, "interface-%04d",jj);
    // static FILE * fp = fopen (interfaceFile, "w");
    FILE * fp = fopen (interfaceFile, "w");
    output_facets (f, fp);
    fclose(fp);
    // fprintf(fp, "\n");
  }

  // {
  //   static FILE * fp = fopen ("mesh", "w");
  //   output_cells (fp);
  //   fprintf(fp, "\n");
  // }

  {
    static FILE * fp = popen ("ppm2mp4 f.mp4", "w");
    output_ppm (f, fp, min = 0, max = 1, box = {{-2,0},{2,4}});
  }

  // // {
  //   static FILE * fp = popen ("ppm2mp4 level.mp4", "w");
  //   scalar l[];
  //   foreach()
  //     l[] = level;
  //   output_ppm (l, fp, min = 0, max = LEVEL, box = {{-2,0},{2,4}});
  // }

  // {
  //   static FILE * fp = popen ("ppm2mp4 omega.mp4", "w");
  //   scalar omega[];
  //   vorticity (u, omega);
  //   output_ppm (omega, fp, linear = true, box = {{-2,0},{2,4}}, spread = 100);
  // // }

  /** 
  If we want to use bview to observe the evolution of the simulation, then we
  define Bview to 1.*/

  #define Bview 1


  #if Bview

  {
    static FILE * fp1 = popen ("ppm2mp4 interface.mp4", "w");
    clear();
    view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
    draw_vof("f", filled = 1);
    mirror({0,1,0}, 0)
      draw_vof("f", filled = 1);
    save(fp = fp1);
  }

  {
    static FILE * fp2 = popen ("ppm2mp4 omega.mp4", "w");

    scalar l[];
    foreach()
      l[] = level;
    clear();
    view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
    draw_vof("f");
    squares("omega", map = cool_warm);
    mirror({0,1,0}, 0) {
      draw_vof("f");
      squares("omega", map = cool_warm);
    }
    save(fp = fp2);
  }

  {
    static FILE * fp3 = popen ("ppm2mp4 level.mp4", "w");
    scalar l[];
    foreach()
      l[] = level;
    clear();
    view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
    draw_vof("f");
    squares("l",  min = 1, max = LEVEL);
    mirror({0,1,0}, 0) {
      draw_vof("f");
      squares("l", min = 1, max = LEVEL);
    }
    save(fp = fp3);
  }

  #endif

  dump ("dump");
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

event simuInfo (i += 5000; t<= END) {
  
  /**   
  It could be usefull to observe the vorticity of the flow. That's why we
  compute it.*/

  scalar omega[];
  vorticity (u, omega);

  char dumpFile[80];
  sprintf (dumpFile, "simulation-%07d", i);
  dump (file = dumpFile);
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

#define MEASURE 0

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