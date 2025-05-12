/**
This example illustrates the slow diffusion of particles from a ighly heterogenous to a fairly homogenous state. 
*/
#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "particle.h"

/**
Initializing a rectangular domain 200 x 60 m with 30 layers and 16 cells
*/
Particles particle_collection;
double depth = 60;
double sigma = 5.0;
int k = 0;
int main() {
  L0 = 200;
  X0 = -100;
  nl = 20;
  N = 16;
  run();
}

/**
Uniform depth and 50000 particles placed uniformly in the top left corner. 
*/
event init (t = 0) {
  foreach(){
    zb[] = -depth;
    foreach_layer(){
        u.x[] = 0;
        h[]   = -zb[]/(double)nl;// + 1.0*exp(-(x-40)*(x-40)/5);
      }
  }
  printf("Particles made!\n");
  particle_collection = new_particles (50000);
  foreach_particle () {
    double a = rand()/(double)(RAND_MAX);
    double b = rand()/(double)(RAND_MAX);
    p().x = -100*a;
    p().z = -5*b;
  }
  periodic(left);
  periodic(right);
}

/**
Simple diffusive event at each time step
*/
event diffuse(i++){
  random_step(particle_collection, 0.5);
  particle_boundary(particle_collection);
}

/**
Plot pdf every three hours
*/

void plot(FILE * fp){
  scalar s = new scalar[nl];
  particle_pdf(particle_collection, s);
  foreach(){
    double zc = zb[];
    foreach_layer(){
      zc += h[]/2.;
      fprintf(fp, "%g %g %g\n", x, zc, s[]);
      zc += h[]/2.; 
    }
    fprintf(fp, "\n");
  }
}

event plot (t <= 28800; t += 7200)
{
  char name[80];
  sprintf(name, "log%d", k);
  FILE * fp = fopen(name, "w");
  plot(fp);
  fclose(fp);
  k++;
}

/** 
We plot the initial state
~~~gnuplot Initial state
reset
set size square
set xr [-100: 100]
set yr [-60: 5]
set pm3d map
set grid
set xlabel 'x [m]'
set ylabel 'z [m]'
splot "log0"
~~~

A middle state

~~~gnuplot Middle state
reset
set size square
set xr [-100: 100]
set yr [-60: 5]
set pm3d map
set grid
set xlabel 'x [m]'
set ylabel 'z [m]'
splot "log1" 
~~~

And end state 

~~~gnuplot End state
reset
set size square
set xr [-100: 100]
set yr [-60: 5]
set pm3d map
set grid
set xlabel 'x [m]'
set ylabel 'z [m]'
splot "log2"
~~~
*/