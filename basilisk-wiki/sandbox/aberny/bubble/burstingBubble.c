/**
In this simulation we observe the bursting of a bubble on a water-air interface.
There is several option in this simulation we can activate to study different
things arround the bursting Bubble process

First, this simulation is long enough to allow us to observe all the jet drops
coming from a bursting bubble event

## Available option
 * ACOUSTIC
 * Gravity
 * BRUIT
 * AdaptOn
 * ROMEO
 * DIV
 * IntTracking
 * Bview
 * MEASURE
 * ENERGY


### ACOUSTIC

This option will increase the domain size and it will add a pressure measurement
far away from the bubble. This was a request from an acoustician friends, in
order to compare the pressure he measured experimentaly, and the one observe by
Basilisk, under the incompressible assumption 
*/

/**
### Gravity

This option will set the gravity forces to true. They are defined by the Bond
number (Bo), with $Bo = \frac{\rho g R^2}{\gamma}$, with $\rho$ being the
liquid density, g the acceleration of gravity, R the radius of the bubble, and
$\gamma$ the surface tension.

### BRUIT

This option (BRUIT is noise in french) will add an initial noise in the
simulation. This noise is a random velocity in the liquid field. This velocity
field is centerd on 0 and has an amplitude based on the velocity of the first
drop.

### AdaptOn

When set to true, this option will allow the use of the adaptiv mesh of Basilisk

### ROMEO

Juliette Pierre asked me if it was possible to track some bubble taht could be
trapped below the simulation. This option allow the simulation to track those
kind of bubble in the liquid phase

### DIV

This option will allow us to observe the divergence field in the simulation

### IntTracking

This option will allow the simulation to output the interface of the simulation
for the first steps of the simulation

### Bview

This option allow the use of bview to output images and video

### MEASURE

This option allow the use of the measure module. This module is a
post-processing tool to measure the position of the interface allong the axis
of symmetry, but also the size, the position and the velocity of all the
droplets coming from the simulation

### ENERGY

This option allow the use of the energy module (which is not working). This
module is a post-processing tool to measure the energy in the simulation

## Publication linked to this simulation

The results on all the drops were published [here](https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.5.033605) in Physical Review of Fluids.

A detail study on the statistic is available [here](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2021GL092919) in Geophysical Research Letters.

## Basilisk option

We want to perform an axi-symetric simulation, with the centered Navier-Stokes
solver, on a 2 fluid simulations.
*/

#include <time.h>

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

/** 
tag.h is a module use to tag the different phase of bubble or liquid
*/
#include "tag.h"

/** 
distance.h is a module use to get the height function
*/

#include "distance.h"

/**
We add profiling option to track the perfs of the simulation. Those module are
not mandatory for the simulation to work properly, but they provide usefull
information*/
#include "navier-stokes/perfs.h"
#include "profiling.h"

/**  
Last option: the gravity*/
#define Gravity 1

#if Gravity
#include "reduced.h"
#endif

/** 
## Physical Parameters

There are 2 parameters for this simulation. The first one is the Laplace
number. The second one is the Bond number. The default value used for the (James) Bond
number is 0.07. It's used to compute the initial shape of the bubble.*/

double La = 10000;
double Bond = 0.07;

/**
The surface tension $\sigma$ is defined by the parameter Sigma. We define it
equal to 1.*/

#define Sigma 1.

/**
The simulation should end at t = END. Iend is used in case of noise to mark the
end of an output*/

#define END 1.2

#define Iend 8500

/**
## Options

Those option have been describe at the begining of this file
*/
#define ACOUSTIC 0


#define BRUIT 0
double amplitude = 0.01; // the noise amplitude

/**
Alternatively, we can add some noise on the interface (only if there is no noise
in the liquid phase). This option is disabled and may not work anymore*/

// #if !BRUIT
//   #define InterfaceNoise 1
// #endif



#define AdaptOn 1

#define ROMEO 0

#define DIV 0

#define IntTracking 0

/**
## Numerical Parameters

LEVEL stand for the maximum refinement of the grid

L0 is for the size of the simulation domain. We can see here the impact of the
ACOUSTIC option*/

int LEVEL = 10;

#if ACOUSTIC
#define L0 25.
#else
#define L0 10.
#endif

/**
The MinLevel option define the minimum level of refinement for the adaptive mesh 
*/
#define MinLevel 4

/**
If we want to output a simulation when there is creation of new liquid
phase, we should define DumpSimulation to 1. Otherwise, it should be 0.*/

#define DumpSimulation 1

/**
We add the boundary condition to allow the liquid drops to escape the simulation
domain once the reach to top of the domain 
*/
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

/**  
intemax, uemax1 and uemax2 are numerical parameters that controll the condition
when the mesh will increase or decrease its refinement level*/
double intemax = 0.005;
double uemax1 = 0.2;
double uemax2 = 0.15;

/**
## The simulation setup*/

int  main(int argc, char const *argv[]) {

  size (L0); // size of the domain
  #if ACOUSTIC
  origin (-20, 0.);
  init_grid (1 << (11)); // initial grid in ACOUSTIC: $2^11$ cells in each direction
  #else
  origin (-L0/2., 0.);
  init_grid (1 << (9)); // initial grid: $2^9$ cells in each direction
  #endif
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
  By default the adaptiv criterion on the velocity is at 0.15, but we can change
  that with an initial input argument.*/
  if (argc >=5){
    uemax2 = atof (argv[4]);
  }

  /**
  By default the initial noise amplitude is at 0.01, but we can change that with
  an initial input argument.*/
  if (argc >=6){
    amplitude = atof (argv[5]);
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
  We have defined the density, the viscosity and the shape with dimensionless
  number. Therefore, we have to fix the surface tension to 1. */

  f.sigma = Sigma;

  /**
  Last, if the gravity is defined, we set the gravity action to the Bond
  number value. */

  #if Gravity
  G.x = -Bond;
  #endif

  /**
  In case of simulation with noise, we initialize the random seed to make sure
  that all the simulation with noise will be different*/
  srand(time(NULL));

  /**
  We print the characteristics of the fluid in the log file. */
  fprintf (stderr, "props %f %f %f %f %f %i\n \n",
	   rho1, mu1, La, Bond, Sigma, LEVEL);
  run();
}

/**
## Some functions

We compute the velocity of Ganan-Calvo (PRL 2017) (with the correction from
Deike et al PRF 2018). This will be the velocity of reference for the noise.*/

double refVelo(double La, double Bo){
  /*
  */
  return sqrt(La)*19*pow(( (1+2.2*Bo) * La * (pow(500,-1./2.)-pow(La,-1./2.)) ),-3./4.);
}

/**
If we want to add noise in the simulation, there is a function for that.*/

double velocityNoise(double amplitude) {
  return( (2*(rand()/(double)RAND_MAX)-1)*amplitude );
}

/**

## Initialisation

We initialize the simulation with the bubble shape. We include the library
findBond in the code. It will solve the bubble shape equation and return a liste
of coordinate.

Then, we use a second function, geometry that will take the hollow circle and
the bubble shape, to proceed to the union of the 2 distance fields.*/

#include "findBond.h"

#include "view.h" //include here because we didn't need to perform output before


event init (t = 0) {
  if (!restore (file = "restart") ){

    /**
    Here, dataShape will contain the coordinate of the water/air interface. */

    coord* dataShape;

    Circle* hollow = NULL;
    Circle fillet;
    hollow = &fillet;

    Circle* cap = NULL;
    Circle topCap;
    cap = &topCap;

    dataShape = shapeBond(Bond, hollow, cap);


    // fprintf(stderr, "x= %f y= %f r= %f\n",fillet.x, fillet.y, fillet.r );
    // for debug purpose

    /**
    We are using the distance function, applied on the list of coordinate
    dataShape. This will generate a distance function. */

    scalar d[];
    distance (d, dataShape);

    #if BRUIT
    double Ca = refVelo(La, Bond); 
    fprintf(stderr, "ref Velocity = %f, velo ref adim = %f \n", 
      Ca, Ca*sqrt(1/La));
    fprintf(stderr, "amplitude = %f\n", amplitude);
    foreach()
      foreach_dimension(){
        if (d[]<0.){
          u.x[] = velocityNoise(amplitude)*Ca;
        }
      }

    /** 
    Since we changed the velocity field, we need to redifine the boundary
    condition*/

    boundary ({u.x});
    boundary ({u.y});
    dump("noiseInit");
    FILE * fp2 = fopen("noiseField", "w");
    // output_field({u.x, u.y}, fp2,  n = 1<<9);
    #endif

    #if AdaptOn
    while (adapt_wavelet ({d}, (double[]){intemax}, LEVEL,9).nf);
    #endif
    

    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */

    #if AdaptOn
    do {
    #endif
      vertex scalar phi[];

      foreach_vertex(){
        double a = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
        double b = -sq (fillet.r) + sq (x - fillet.y) + sq (y - fillet.x);
        double c = min(a,b);

        phi[] = -c;
      }

      /**
      We can now initialize the volume fraction of the domain. */

      fractions (phi, f);
    #if AdaptOn
    } while(adapt_wavelet({f}, (double[]){intemax}, LEVEL, 9).nf);
    #endif

    /**
    We output the initial situation */
    static FILE * fp = fopen("initial.png", "w");
    output_ppm(f, fp, 400, min = 0, max = 1, box = {{-2.5,0},{2.5,3}}); // a png  file, without bview
    dump("initial");
    clear();
    view (fov = 4.69417, quat = {0,0,-0.707108,0.707106}, tx = 0., ty = 0.0855609, bg = {1,1,1}, width = 1280, height = 720, samples = 4);
    draw_vof("f", filled = 1);
    mirror({0,1,0}, 0)
      draw_vof("f", filled = 1);
    save("initial.ppm");

    #if BRUIT
      clear();
      view (fov = 12, quat = {0,0,-0.707108,0.707106}, tx = 0., ty = 0.0855609, bg = {1,1,1}, width = 1280, height = 720, samples = 4);
      draw_vof("f");
      squares("u.x", min = -0.1, max = 0.1, map = cool_warm);
      mirror({0,1,0}, 0){
        draw_vof("f");
        squares("u.y", min = -0.1, max = 0.1, map = cool_warm);}
      save("initialNoise.png");
    #endif //BRUIT
  }
  else { // restart Case
    FILE * fp = fopen ("initialRestore.png", "w");
    output_ppm(cm, fp, 400);
  }
}

/**
## Outputing noise in the simulation

When we have noise, we output in mp4 file the first Iend steps.

We show the velocity field and the vorticity field*/
#if (BRUIT || InterfaceNoise)

event noiseMovie (i++; i < Iend) {
  output_field();
  {
    static FILE * fpNoise = popen ("ppm2mp4 noise.mp4", "w");
    clear();
    view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
    draw_vof("f");
    squares("u.x", map = cool_warm, min = -0.5, max = 0.5);
    mirror({0,1,0}, 0) {
      draw_vof("f");
      squares("u.y",map = cool_warm, min = -0.5, max = 0.5);
    }
    save(fp = fpNoise);
  }
  {
    scalar omega[];
    vorticity (u, omega);
    static FILE * fpNoise = popen ("ppm2mp4 noiseOmega.mp4", "w");
    clear();
    view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
    draw_vof("f");
    squares("omega", min = -5, max = 5, linear = true);
    mirror({0,1,0}, 0) {
      draw_vof("f");
      squares("omega", min = -5, max = 5, linear = true);
    }
    save(fp = fpNoise);
  }
}

#endif

/**

## Adaptation

At each step, we will adapt the mesh. We make sure that the mesh size will
be at least $1/2^MinLevel$ */

#if AdaptOn
event adapt (i++) {
  if (i<20)
    adapt_wavelet ({f,u}, (double []){intemax, uemax1, uemax1, uemax1},
		 maxlevel = LEVEL, minlevel = MinLevel);
  else 
    adapt_wavelet ({f,u}, (double []){intemax, uemax2, uemax2, uemax2},
     maxlevel = LEVEL, minlevel = MinLevel);
}
#endif

/**
## General output

We output a few files. One will display the interface to observe the evolution of 
the simulation. The other are video files, output when Bview is define to 1*/

int jj = 0;

event outputInterface (i += 500) {
// event outputFilm (t += 0.001) { // it's better to have an output based on the steps rather than the time
  
  /**
  We can output the interface in a txt file 
  */
  // {
    // static FILE * fp = fopen ("interface", "w");
    // jj++;
    // char interfaceFile[80];
    // sprintf(interfaceFile, "interface-%03d",jj);
    // static FILE * fp = fopen (interfaceFile, "w");
    // FILE * fp = fopen (interfaceFile, "w");
    // output_facets (f, fp);
    // fprintf(fp, "\n");
  // }

  /** 
  If we want to use bview to observe the evolution of the simulation, we
  define Bview to 1.*/

  #define Bview 1


  #if Bview
  scalar omega[];
  vorticity (u, omega);
  /**
  We display the evolution of the interface. The fluid will be in black, the gas
  in white*/
  {
    static FILE * fp1 = popen ("ppm2mp4 interface.mp4", "w");
    clear();
    view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
    draw_vof("f", filled = 1);
    mirror({0,1,0}, 0)
      draw_vof("f", filled = 1);
    save(fp = fp1);
  }
  /* 
  We show the vorticity field*/
  {
    static FILE * fp2 = popen ("ppm2mp4 omega.mp4", "w");
    // scalar omega[];
    // vorticity (u, omega);
    // scalar l[];
    // foreach()
    //   l[] = level;
    clear();
    view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
    draw_vof("f");
    squares("omega", map = jet, min = -50, max = 50, linear = true);
    // squares("omega", map = bwr, min = -50, max = 50, linear = true);
    mirror({0,1,0}, 0) {
      draw_vof("f");
      squares("-omega", map = jet, min = -50, max = 50, linear = true);
      // squares("-omega", map = bwr, min = -50, max = 50, linear = true);
    }
    save(fp = fp2);
  }
  // The colormap bwr is not yet in the general version of basilisk.
  
  /**
  The velocity field.*/

  // {
  //   static FILE * fp2 = popen ("ppm2mp4 ux.mp4", "w");
  //   scalar omega[];
  //   vorticity (u, omega);
  //   scalar l[];
  //   foreach()
  //     l[] = level;
  //   clear();
  //   view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
  //   draw_vof("f");
  //   squares("u.x", map = jet, min = -30, max = 30);
  //   mirror({0,1,0}, 0) {
  //     draw_vof("f");
  //     squares("u.x", map = jet, min = -30, max = 30);
  //   }
  //   save(fp = fp2);
  // }
  /**
  The vorticity field, with the cool_warm colormap. The liquid is in black, the
  vorticity is only showed in the gas*/

  // { static FILE * fp2 = popen
  //   ("ppm2mp4 omegaFilled.mp4", "w");

  //   scalar l[];
  //   foreach()
  //     l[] = level;
  //   clear();
  //   view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
  //   draw_vof("f", filled = 1);
  //   squares("omega", map = cool_warm, min = -30, max = 30);
  //   mirror({0,1,0}, 0) {
  //     draw_vof("f", filled = 1);
  //     squares("omega", map = cool_warm, min = -30, max = 30);
  //   }
  //   save(fp = fp2);
  // }

  /**
  The "level" field. The color is based on the level of resolution*/
  {
    static FILE * fp3 = popen ("ppm2mp4 level.mp4", "w");
    scalar l[];
    foreach()
      l[] = level;
    clear();
    view (fov = 12.699, quat = {0,0,0.707107,-0.707107}, ty = 0.00253494, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
    draw_vof("f", filled= 1);
    squares("l",  min = 1, max = LEVEL);
    mirror({0,1,0}, 0) {
      draw_vof("f");
      squares("l", min = 1, max = LEVEL);
    }
    save(fp = fp3);
  }
  /**
  The "level" field again, but with a zoom.*/
  {
    static FILE * fplz = popen ("ppm2mp4 levelZoom.mp4", "w");
    scalar l[];
    foreach()
      l[] = level;
    clear();
    view (fov = 4.42681, quat = {0,0,0,1}, tx = 0.0579169, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
    draw_vof("f");
    squares("l",  min = 1, max = LEVEL);
    cells();
    mirror({0,1,0}, 0) {
      draw_vof("f");
      squares("l", min = 1, max = LEVEL);
      cells();
    }
    save(fp = fplz);
  }

  /**
  The vorticity field, with a zoom this time*/
  {
    static FILE * fpz = popen("ppm2mp4 omegaZoom.mp4", "w");
    clear();
    view (fov = 6.69986, quat = {0,0,0,1}, tx = -0.05, ty = 0, bg = {1,1,1}, width = 1920, height = 1080, samples = 4);
  
    squares("omega", min = -200, max = 200, map = jet, linear = true);
    // squares("omega", min = -200, max = 200, map = bwr, linear = true);
    draw_vof("f");
    mirror(n = {0,1,0}, alpha = 0){
      squares("-omega", min = -200, max = 200, map = jet, linear = true);
      // squares("omega", min = -200, max = 200, map = bwr, linear = true);
      draw_vof("f");
    }
    save(fp = fpz);
  }

  #if ROMEO
  /** 
  The liquid/gas interface, below the bubble to observe a potential tiny bubble.*/
  {
  static FILE * romeo = popen("ppm2mp4 zoom.mp4", "w");
  scalar l[];
  foreach()
    l[] = level;
  clear;
  view (fov = 0.647604, quat = {0,0,-0.707107,0.707107}, tx = -0.000410862, ty = 0.160506, bg = {1,1,1}, width = 1476, height = 870, samples = 4);
  draw_vof("f", filled = 1);
    squares("l",  min = 1, max = LEVEL);
    mirror({0,1,0}, 0) {
      draw_vof("f");
      squares("l", min = 1, max = LEVEL);
    }
  save(fp = romeo);
  }
  #endif


  #endif
  {
    char inter[80];
    sprintf(inter, "interface-t-%g", t);
    FILE * fp = fopen(inter, "w");
    output_facets(f, fp);
    fclose(fp);
  }

  dump ("dump"); // for a potential restart
}
/**
This small function is used to detect a time step. This allow to have an output
at a wanted time steps, without forcing the solver to reach it precisely*/
int timeDetect (double t, double tPrevious, double target) {
  if (tPrevious == -1)
    return -1;
  if ((tPrevious - target <0.) && (t - target > 0.))
    return 1;
  else
    return -1;
}


double tAvant = -1;
/**
We output a movie that travel allong the jet

The simulation evolution looks like:

<p><center>
<video width="1280" height="720" controls>
  <source src="burstingBubble/interface.mp4" type "video/mp4">
  Your browser does not support the video tag
</video>
<br>
The bubble bursting
</center></p>


We output a file at regular time step until the end of the simulation
*/


event simuInfo (i += 500; t<= END) {
// event simuInfo (i += 100; t<= END) {
// event simuInfo(i+=2500; i<=Iend){
  /**   
  It could be usefull to observe the vorticity of the flow. That's why we
  compute it.*/

  vector h[];
  heights (f, h);

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

#define MEASURE 1

#if MEASURE
#include "measure.h"
#endif


#define ENERGY 0

#if ENERGY
#include "energy.h"
#endif

#if ACOUSTIC
#include "acoustic.h"
#endif

event logfile (i++) {
  fprintf(stderr, "%g %i %i %i\n", t, i, mgpf.i, mgpf.nrelax);
}


/**
## End

As verification, we output the final state. */


// event finalState (i = Iend) {
//   dump(file = "final");
// }

event finalState (t = end) {
  dump(file = "final");
}


/**
## Post-Processing
We are interested in the velocity of the interface along the axis of 
symmetry.

~~~gnuplot Velocity along the axis of symmetry
set ratio -1
set xlabel 'Time'
set ylabel 'Velocity'
set yrange [0:0.4]
set grid
plot 'out' u 2:6 
~~~

We can also follow the evolution of the position of the different elements
(drops and interface).

~~~gnuplot Position of the elements along the axis of symmetry
set xlabel 'Time'
set ylabel 'Position on the axis'
plot 'out' u 2:4 w l t "jet position", 'drop' u 2:4 every ::0::0 t "first drop", 'drop' u 2:4 every ::1::1 t "second drop", 'drop' u 2:4 every ::2::2 t "third drop"
~~~

Max number of cells:

~~~bash
max=`awk '{max_val=($6<max_val)?max_val:$6;} END{print max_val;}' log`
~~~

Evolution of the jet position:

~~~gnuplot Position of the jet and of the droplets
set xlabel 'Times'
set ylabel 'x'
set xrange [0:1]
set yrange [-2:]

plot 'out' u 2:4 every :::0::0 w l t "jet" lt rgb "#FF0000",\
     'dropCorr' u 2:4 every ::0::0 t "highest drop" ps 0.1,\
     'dropCorr' u 2:4 every ::1::1 t "second highest drop" ps 0.1,\
     'dropCorr' u 2:4 every ::2::2 t "third highest drop" ps 0.1,\
     'dropCorr' u 2:4 every ::3::3 t "fourth highest drop" ps 0.1
~~~
*/
