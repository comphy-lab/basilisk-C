/** 
This small code is a test case for the initialisation of the code
burstingBubble.c

Here we will compute the shape of a bubble of an equivalent radius 1. The 
volume of this bubble should then be equal to $4/3\pi$


Since we will use the basilisk structure to measure the volume of our bubble, 
we need to include several library.

We will work with 2 phases. The bubble will be one of them. Then, the rest of
the domain will be the second phase. We wont pay attention to the free
surface.*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

/**
We will need the distance library for the interface computation.*/
#include "distance.h"


/** bubbleShape and findBond are our custome library we are testing with this
code.*/

#define WRONG 0

#include "bubbleShape.h"
#include "findBond.h"

/**
To display what we are observing, we will use bview*/

#include "view.h"

/** 
We will work with a Bond number at 1. This is indeed a case where the bubble
emerge a lot above the free surface. It means it's a case where we can make the
biggest mistake.*/

double Bond = 1;

/** We set sigma, mu and rho to 1 since we wont run the simulation after the
initialisation.*/

#define Sigma 1.
#define Mu 1.
#define Rho 1.

/**
Numerical parameters*/

#define LEVEL 12
#define L0 5.

int main(int argc, char const *argv[]){
  
  size (L0);
  origin (-L0/2., 0.);
  init_grid (1 << (8));

  if (argc >= 2){
      Bond = atof (argv[1]);
  }


  f.sigma = Sigma;
  rho1 = 1., mu1 = Mu;
  rho2 = 1., mu2 = Mu;

  run();
}

event init (t = 0) {

  /**
  Here, dataShape will contain the coordinate of the water/air interface. */

  coord* dataShape;

  /**   
  The function bubbleOnly will gave us the shape of the bubble and the
  bubble only. We can eventually translate this bubble by adding a second 
  argument (which is 0 by default   */

  int* size = NULL;
  int taille = 0;
  size = &taille;
  
  dataShape = bubbleOnly(Bond, 0, size);

  /**
  We are using the distance function, applied on the list of coordinate
  dataShape. This will generate a distance function. */

  scalar d[];
  distance (d, dataShape);

  while (adapt_wavelet ({d}, (double[]){0.005}, LEVEL).nf);


  /**
  The distance function is defined at the center of each cell, we have
  to calculate the value of this function at each vertex. */

  vertex scalar phi[];

  foreach_vertex(){
    double a = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    phi[] = a;
  }

  /**
  We can now initialize the volume fraction of the domain. */

  fractions (phi, f);

  adapt_wavelet ({f}, (double []){0.001},
     maxlevel = LEVEL);
  
  clear();
  view (fov = 13, quat = {0,0,-0.707108,0.707106}, tx = 0., ty = 0.0855609, bg = {1,1,1}, width = 1280, height = 720, samples = 4);
  draw_vof("f", filled = 1);
  mirror({0,1,0}, 0)
    draw_vof("f", filled = 1);
  save("initial.png");

  /**
  The air bubble we are computing looks like that:
  ![a bubble](bubbleRef/initial.png)
  */

  /**
  We tag the liquid and we measure the volume of the surface*/
  double v = 0;

  foreach(reduction(+:v)){
    v += dv()*f[];
  }

  /**   In axi, the volume of a sphere with radius equal to one is $2/3$. We
compare our volume with $2/3$*/

  double C = v/(2./3.);
  C = pow(C,1./3.);

  fprintf(stderr, "volume calc: %f\n", v);
  fprintf(stderr, "volume th: %f\n", (2./3.));
  dump("initial");
  #if WRONG
  fprintf(stderr, "C(Bo*) = C(%f) = %f\n",Bond,C );
  #endif
}